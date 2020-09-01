# standard error

#setwd('/Users/rebeccakahn/Dropbox/Git_repo/ADAGIO/')

#Parameters --------
setwd("/n/holyscratch01/lipsitch_lab/rkahn/ADAGIO")

source('./simfixoutbreak.R')
source('./make_network.R')
source('./network_epidemic.R')
source('./analyze_data.R')
source('./clustering algorithm.R')


require(NetSurv)
require(Matrix)
require(Rlab)
require(igraph)
require(deSolve)
require(reshape2)
require(ggplot2)
require(NetSurv)
require(Matrix)
require(Rlab)
require(igraph)
require(deSolve)
require(caTools)
library(survival)
library(coxme)
library(frailtypack)
require(tidyverse)


args=(commandArgs(TRUE))


# Number of simulated trials
nsim<-1

# Population structure parameters:
# Average size of one community
ave_community_size<-5000
# Range of community sizes (sizes are uniformly distributed on this range)
community_size_range<-0
# Number of communities
num_communities<-2
# Probability of an edge between two nodes in the same community
rate_within<-0.02
# Probability of an edge between two nodes in different communities
rate_between<-0.001

# Disease characteristics:
# Per-time-step hazard of infection for a susceptible nodes from an infectious neighbour
beta<-0.0029
# Expected number of importations to the population over two years
num_introductions<-20

# All or nothing efficacy of vaccine
direct_VE<-0.6
infect_VE <-0.3

# Gamma-distribution parameters of incubation and infectious period
incperiod_shape<-3.11
incperiod_rate<-0.32
infperiod_shape<-1.13
infperiod_rate<-0.226
ave_inc_period <- ceiling(incperiod_shape/incperiod_rate)

# First day of trial enrollment, relative to start of epidemic
trial_startday<-1
# Days of follow-up
trial_length<-300
# Number of days over which subjects are enrolled
enrollment_period<-1
# Target community enrollment proportion
num_clusters_enrolled_per_day<-2
if (num_clusters_enrolled_per_day > num_communities) {stop("Enrolling too many communities!")}
cluster_coverage<-1 # not enrolling whole population so for cluster randomization will get overall (not total)
# if imperfect sampling, percent that are sampled:
sample_percent <- 1

gen_len <- 18958  #Campbell 2018
mut.rate <- 0.012 #0.000025
bn <- 10

sim <- 1

nbootstrap <- 100

# serial interval distribution
N <- 10000
inf.period.start <- rgamma(N,shape=incperiod_shape,rate=incperiod_rate)
inf.period <- rgamma(N,shape=infperiod_shape,rate=infperiod_rate)
inf.period.end <- inf.period.start + inf.period
distribution <- rep(NA,30)
for (t in 1:trial_length){
  distribution[t] <- sum(ifelse(c(inf.period.start<=t & inf.period.end>=t),1,0))/N
}
distribution_norm <- distribution/sum(distribution)

# Calculate R0
R0 <- (1 - (infperiod_rate/(infperiod_rate+beta))^infperiod_shape) *
  (((ave_community_size-1)*(1-rate_within)*rate_within + 
      (num_communities-1)*ave_community_size*(1-rate_between)*rate_between + 
      ((ave_community_size-1)*rate_within + (num_communities-1)*ave_community_size*rate_between)^2)/
     ((ave_community_size-1)*rate_within + (num_communities-1)*ave_community_size*rate_between)-1)
R0

for (i in 1:length(args)) {
  eval (parse (text = args[[i]] ))
}

j<-j

list <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
  args <- as.list(match.call())
  args <- args[-c(1:2,length(args))]
  length(value) <- length(args)
  for(i in seq(along=args)) {
    a <- args[[i]]
    if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
  }
  x
}

js: 1:500
Bootstrap_master <- data.frame()

for (x in js){
  cat(x,"\n")

  load(paste0(x,"_g",R0,"_",infect_VE,"_",direct_VE,"_",mut.rate,"_",bn,"_",infperiod_shape,"_",sample_percent,"_","sample.Rdata"))
  results_edge_list_analysis <- read.csv(paste0(x,"_results_edge_list_analysis",R0,"_",infect_VE,"_",direct_VE,"_",mut.rate,"_",bn,"_",infperiod_shape,"_",sample_percent,"_","sample_new.csv"))
  results_edge_list_analysis$Infector[is.na(results_edge_list_analysis$Infector) & results_edge_list_analysis$eventstatus==1] <- 0
  results_infected <- results_edge_list_analysis[results_edge_list_analysis$eventstatus==1,]
  results_infected$DayRecovered[is.na(results_infected$DayRecovered)] <- trial_startday + trial_length + 1
  #print(head(results_infected))
  
  results_contacts2 <- read.csv(paste0(x,"_results_contacts2",R0,"_",infect_VE,"_",direct_VE,"_",mut.rate,"_",bn,"_",infperiod_shape,"_",sample_percent,"_","sample_new.csv"))
  results_edge_list_analysis_scoremat <- read.csv(paste0(x,"_results_edge_list_analysis_DS_scoremat",R0,"_",infect_VE,"_",direct_VE,"_",mut.rate,"_",bn,"_",infperiod_shape,"_",sample_percent,"_","sample_new.csv"))
  results_edge_list_analysis_scoremat <- results_edge_list_analysis_scoremat[!(is.na(results_edge_list_analysis_scoremat$prob_SI_norm)),]
  
  results_edge_list_analysis_hybrid <- read.csv(paste0(x,"_results_edge_list_analysis_DS_hybrid",R0,"_",infect_VE,"_",direct_VE,"_",mut.rate,"_",bn,"_",infperiod_shape,"_",sample_percent,"_","sample_new.csv"))
  results_edge_list_analysis_dist <- read.csv(paste0(x,"_results_edge_list_analysis_DS_dist",R0,"_",infect_VE,"_",direct_VE,"_",mut.rate,"_",bn,"_",infperiod_shape,"_",sample_percent,"_","sample_new.csv"))
  results_edge_list_analysis_W <- read.csv(paste0(x,"_results_edge_list_analysis_DS_W",R0,"_",infect_VE,"_",direct_VE,"_",mut.rate,"_",bn,"_",infperiod_shape,"_",sample_percent,"_","sample_new.csv"))
  
  coxmodel <- function(data) {
    
    survmodel<-try(coxph(Surv(DayInfected,eventstatus)~TrialStatus+strata(Community),data),silent=T)
    usesurvmod <- !inherits(survmodel, 'try-error')
    
    if (usesurvmod && vcov(survmodel)>=0){
      # If no error was thrown and the variance is positive, use the results of the model
      
      vaccEffEst <- 1-exp(survmodel$coefficient + c(0, 1.96, -1.96)*sqrt(survmodel$var))
      zval <- survmodel$coefficient/sqrt(survmodel$var)
      pval <- pnorm(zval, lower.tail = vaccEffEst[1]>0)*2
      
    } else {
      
      vaccEffEst<-c(NA,NA,NA)
      pval <- NA
      
    }
    
    return(vaccEffEst)
    
  }
  
  VEs <- coxmodel(results_edge_list_analysis)[1]
  
  # Initialize vectors ----
  VEIs_bootstrap <- rep(NA,nbootstrap)
  VEIs_contact_bootstrap_prob <- rep(NA,nbootstrap)
  VEIs_contact_bootstrap_cluster2 <- rep(NA,nbootstrap)
  VEIs_SV_cluster_bootstrap <- rep(NA,nbootstrap)
  VEIs_SV_contacts_restrict_bootstrap <- rep(NA,nbootstrap)
  VEIs_hybrid_cluster_bootstrap <- rep(NA,nbootstrap)
  VEIs_hybrid_contacts_restrict_bootstrap <- rep(NA,nbootstrap)
  
  # VEI bootstrap ----
  results_infected2 <- results_infected[,c("InfectedNode","TrialStatus","Infector","infector_state")]
  head(results_infected2)
  for (i in 1:nbootstrap){
    ID <- unique(results_infected2$InfectedNode[!is.na(results_infected2$InfectedNode)])
    ID_bootstrap <- sample(ID,length(ID),replace=TRUE)
    results_edge_list_analysis_bootstrap <- data.frame("InfectedNode"=NA,"TrialStatus"=NA,"Infector"=NA,"infector_state"=NA)
    for (k in 1:length(ID_bootstrap)){
      results_edge_list_analysis_bootstrap <- 
        rbind(results_edge_list_analysis_bootstrap,results_infected2[results_infected2$InfectedNode==ID_bootstrap[k],])
    }
    summary_infectors <- aggregate(results_edge_list_analysis_bootstrap$InfectedNode,by=list(results_edge_list_analysis_bootstrap$infector_stat),length)
    names(summary_infectors)<-c("infector_stat","number")
    vacc_inf <- summary_infectors$number[summary_infectors$infector_stat==1]
    cont_inf<- summary_infectors$number[summary_infectors$infector_stat==0]
    if (length(vacc_inf)==0 & length(cont_inf)==0){
      RRs_bootstrap <- NA
    } else if (length(vacc_inf)==0 & length(cont_inf)>0){
      RRs_bootstrap <- 0
    } else if (length(vacc_inf)>0 & length(cont_inf)==0){
      RRs_bootstrap <- 1
    } else{
      RRs_bootstrap <- vacc_inf/cont_inf
    } 
    if (length(RRs_bootstrap)>0){
      VEIs_bootstrap[i] <- 1-RRs_bootstrap/(1-VEs[1])
    }
    #print(i)
  }
  ###
  
  VEIs_bootstrap %>%
    as.data.frame() %>%
    setNames("Estimate") %>%
    add_column(j=x) %>%
    add_column(method="VEIs") -> VEIs_bootstrap
  
  
  # VEI contacts bootstraps prob ----
  results_contacts2 <- results_contacts2[,c("InfectedNode","node_stat","Contact","infector_stat","SI_prob_norm")]
  names(results_contacts2) <- c("InfectedNode","TrialStatus","Infector","infector_state","SI_prob_norm")
  head(results_contacts2)
  for (i in 1:nbootstrap){
    ID <- unique(results_contacts2$InfectedNode[!is.na(results_contacts2$InfectedNode)])
    ID_bootstrap <- sample(ID,length(ID),replace=TRUE)
    results_edge_list_analysis_bootstrap <- data.frame("InfectedNode"=NA,"TrialStatus"=NA,"Infector"=NA,"infector_state"=NA,"SI_prob_norm"=NA)
    for (k in 1:length(ID_bootstrap)){
      results_edge_list_analysis_bootstrap <- 
        bind_rows(results_edge_list_analysis_bootstrap,results_contacts2[results_contacts2$InfectedNode==ID_bootstrap[k],])
    }
    summary_infectors <- aggregate(results_edge_list_analysis_bootstrap$SI_prob_norm,by=list(results_edge_list_analysis_bootstrap$infector_stat),sum)
    names(summary_infectors)<-c("infector_stat","number")
    vacc_inf <- summary_infectors$number[summary_infectors$infector_stat==1]
    cont_inf<- summary_infectors$number[summary_infectors$infector_stat==0]
    if (length(vacc_inf)==0 & length(cont_inf)==0){
      RRs_bootstrap <- NA
    } else if (length(vacc_inf)==0 & length(cont_inf)>0){
      RRs_bootstrap <- 0
    } else if (length(vacc_inf)>0 & length(cont_inf)==0){
      RRs_bootstrap <- 1
    } else{
      RRs_bootstrap <- vacc_inf/cont_inf
    } 
  if (length(RRs_bootstrap)>0){
    VEIs_contact_bootstrap_prob[i] <- 1-RRs_bootstrap/(1-VEs[1])
  }
    #print(i)
  }
  
  VEIs_contact_bootstrap_prob %>%
    as.data.frame() %>%
    setNames("Estimate") %>%
    add_column(j=x) %>%
    add_column(method="Contacts_prob") -> VEIs_contact_bootstrap_prob
  
  
  # VEI contacts bootstraps cluster ----
  head(results_contacts2)
  if (nrow(results_contacts2)>1){
    for (i in 1:nbootstrap){
      ID <- unique(results_contacts2$InfectedNode[!is.na(results_contacts2$InfectedNode)])
      ID_bootstrap <- sample(ID,length(ID),replace=TRUE)
      results_edge_list_analysis_bootstrap <- data.frame("InfectedNode"=NA,"TrialStatus"=NA,"Infector"=NA,"infector_state"=NA,"SI_prob_norm"=NA)
      for (k in 1:length(ID_bootstrap)){
        results_edge_list_analysis_bootstrap <- 
          bind_rows(results_edge_list_analysis_bootstrap,results_contacts2[results_contacts2$InfectedNode==ID_bootstrap[k],])
      }
      cluster <- results_edge_list_analysis_bootstrap
      names(cluster)[1] <- "indIDVar.2"
      names(cluster)[3] <- "indIDVar.1"
      cluster <- cluster[!is.na(cluster$indIDVar.2),]
      cluster$SI_prob_norm[is.na(cluster$SI_prob_norm)] <- 0
      cluster_sum <- aggregate(cluster$indIDVar.2,by=list(cluster$indIDVar.2),length) 
      if (mean(cluster_sum$x)>1){
        cluster_test_hc <- clusterInfectors(cluster, "indIDVar", "SI_prob_norm",
                                            clustMethod = "hc_absolute",
                                            0.2)
        cluster_test_hc %>%
          subset(cluster==1) %>%
          group_by(indIDVar.2) %>%
          mutate(prob_SI_norm_c = SI_prob_norm/sum(SI_prob_norm)) -> cluster_test_hc_norm
        
        vacc_inf_contacts_cluster1 <- sum(cluster_test_hc_norm$prob_SI_norm_c[cluster_test_hc_norm$infector_state==1 & cluster_test_hc_norm$cluster==1],na.rm=TRUE)
        cont_inf_contacts_cluster1<- sum(cluster_test_hc_norm$prob_SI_norm_c[cluster_test_hc_norm$infector_state==0  & cluster_test_hc_norm$cluster==1],na.rm=TRUE)
        if (length(vacc_inf_contacts_cluster1)==0 & length(cont_inf_contacts_cluster1)==0){
          RRs_bootstrap <- NA
        } else if (length(vacc_inf_contacts_cluster1)==0 & length(cont_inf_contacts_cluster1)>0){
          RRs_bootstrap <- 0
        } else if (length(vacc_inf_contacts_cluster1)>0 & length(cont_inf_contacts_cluster1)==0){
          RRs_bootstrap <- 1
        } else{
          RRs_bootstrap <- vacc_inf_contacts_cluster1/cont_inf_contacts_cluster1
        }   
      } else{
        RRs_bootstrap <- NA
      }
      if (length(RRs_bootstrap)>0){
        VEIs_contact_bootstrap_cluster2[i] <- 1-RRs_bootstrap/(1-VEs[1])
      }
      #print(i)
    }
  }
  
  VEIs_contact_bootstrap_cluster2 %>%
    as.data.frame() %>%
    setNames("Estimate") %>%
    add_column(j=x) %>%
    add_column(method="Contacts_cluster") -> VEIs_contact_bootstrap_cluster2
  
  if (length(results_edge_list_analysis_scoremat)>2){
    
    # Shared variants  ----
    head(results_edge_list_analysis_scoremat)
    results_edge_list_analysis_scoremat2 <- results_edge_list_analysis_scoremat[,c("InfectedNode","node_stat","Infector","infector_stat","prob_SI_norm")]
    
    for (i in 1:nbootstrap){
      ID <- unique(results_edge_list_analysis_scoremat2$InfectedNode[!is.na(results_edge_list_analysis_scoremat2$InfectedNode)])
      ID_bootstrap <- sample(ID,length(ID),replace=TRUE)
      results_edge_list_analysis_scoremat2_bootstrap <- data.frame("InfectedNode"=NA,"node_stat"=NA,"Infector"=NA,"infector_stat"=NA,"prob_SI_norm"=NA)
      for (k in 1:length(ID_bootstrap)){
        results_edge_list_analysis_scoremat2_bootstrap <- 
          rbind(results_edge_list_analysis_scoremat2_bootstrap,results_edge_list_analysis_scoremat2[results_edge_list_analysis_scoremat2$InfectedNode==ID_bootstrap[k],])
      }  
      cluster <- results_edge_list_analysis_scoremat2_bootstrap
      names(cluster)[1] <- "indIDVar.2"
      names(cluster)[3] <- "indIDVar.1"
      cluster <- cluster[!is.na(cluster$indIDVar.2),]
      cluster$prob_SI_norm[is.na(cluster$prob_SI_norm)] <- 0
      cluster_sum <- aggregate(cluster$indIDVar.2,by=list(cluster$indIDVar.2),length) 
      if (mean(cluster_sum$x)>1){
        cluster_test_hc <- clusterInfectors(cluster, "indIDVar", "prob_SI_norm",
                                            clustMethod = "hc_absolute",
                                            0.2)
        cluster_test_hc %>%
          subset(cluster==1) %>%
          group_by(indIDVar.2) %>%
          mutate(prob_SI_norm_c = prob_SI_norm/sum(prob_SI_norm)) -> cluster_test_hc_norm
        
        vacc_inf_contacts_cluster1 <- sum(cluster_test_hc_norm$prob_SI_norm_c[cluster_test_hc_norm$infector_stat==1 & cluster_test_hc_norm$cluster==1],na.rm=TRUE)
        cont_inf_contacts_cluster1<- sum(cluster_test_hc_norm$prob_SI_norm_c[cluster_test_hc_norm$infector_stat==0  & cluster_test_hc_norm$cluster==1],na.rm=TRUE)
        if (length(vacc_inf_contacts_cluster1)==0 & length(cont_inf_contacts_cluster1)==0){
          RRs_bootstrap <- NA
        } else if (length(vacc_inf_contacts_cluster1)==0 & length(cont_inf_contacts_cluster1)>0){
          RRs_bootstrap <- 0
        } else if (length(vacc_inf_contacts_cluster1)>0 & length(cont_inf_contacts_cluster1)==0){
          RRs_bootstrap <- 1
        } else{
          RRs_bootstrap <- vacc_inf_contacts_cluster1/cont_inf_contacts_cluster1
        }   
      } else{
        RRs_bootstrap <- NA
      }
      if (length(RRs_bootstrap)>0){
        VEIs_SV_cluster_bootstrap[i] <- 1-RRs_bootstrap/(1-VEs[1])
      }
      #print(i)
    }

    ###
    
    # Shared variants restrict contacts -----
    results_edge_list_analysis_scoremat_contacts <- results_edge_list_analysis_scoremat[results_edge_list_analysis_scoremat$contacts!=0 & !is.na(results_edge_list_analysis_scoremat$contacts),]
    results_edge_list_analysis_scoremat_contacts2 <- results_edge_list_analysis_scoremat_contacts[,c("InfectedNode","node_stat","Infector","infector_stat","prob_SI_norm")]
    
    if (length(unique(results_edge_list_analysis_scoremat_contacts2$InfectedNode))>1){
    
      for (i in 1:nbootstrap){
        ID <- unique(results_edge_list_analysis_scoremat_contacts2$InfectedNode[!is.na(results_edge_list_analysis_scoremat_contacts2$InfectedNode)])
        ID_bootstrap <- sample(ID,length(ID),replace=TRUE)
        results_edge_list_analysis_scoremat_contacts2_bootstrap <- data.frame("InfectedNode"=NA,"node_stat"=NA,"Infector"=NA,"infector_stat"=NA,"prob_SI_norm"=NA)
        for (k in 1:length(ID_bootstrap)){
          results_edge_list_analysis_scoremat_contacts2_bootstrap <- 
            rbind(results_edge_list_analysis_scoremat_contacts2_bootstrap,results_edge_list_analysis_scoremat_contacts2[results_edge_list_analysis_scoremat_contacts2$InfectedNode==ID_bootstrap[k],])
        }  
        cluster <- results_edge_list_analysis_scoremat_contacts2_bootstrap
        names(cluster)[1] <- "indIDVar.2"
        names(cluster)[3] <- "indIDVar.1"
        cluster <- cluster[!is.na(cluster$indIDVar.2),]
        cluster$prob_SI_norm[is.na(cluster$prob_SI_norm)] <- 0
        cluster_sum <- aggregate(cluster$indIDVar.2,by=list(cluster$indIDVar.2),length) 
        if (mean(cluster_sum$x)>1){
          cluster_test_hc <- clusterInfectors(cluster, "indIDVar", "prob_SI_norm",
                                              clustMethod = "hc_absolute",
                                              0.2)
          cluster_test_hc %>%
            subset(cluster==1) %>%
            group_by(indIDVar.2) %>%
            mutate(prob_SI_norm_c = prob_SI_norm/sum(prob_SI_norm)) -> cluster_test_hc_norm
          
          vacc_inf_contacts_cluster1 <- sum(cluster_test_hc_norm$prob_SI_norm_c[cluster_test_hc_norm$infector_stat==1 & cluster_test_hc_norm$cluster==1],na.rm=TRUE)
          cont_inf_contacts_cluster1<- sum(cluster_test_hc_norm$prob_SI_norm_c[cluster_test_hc_norm$infector_stat==0  & cluster_test_hc_norm$cluster==1],na.rm=TRUE)
          if (length(vacc_inf_contacts_cluster1)==0 & length(cont_inf_contacts_cluster1)==0){
            RRs_bootstrap <- NA
          } else if (length(vacc_inf_contacts_cluster1)==0 & length(cont_inf_contacts_cluster1)>0){
            RRs_bootstrap <- 0
          } else if (length(vacc_inf_contacts_cluster1)>0 & length(cont_inf_contacts_cluster1)==0){
            RRs_bootstrap <- 1
          } else{
            RRs_bootstrap <- vacc_inf_contacts_cluster1/cont_inf_contacts_cluster1
          }   
        } else{
          RRs_bootstrap <- NA
        }
        if (length(RRs_bootstrap)>0){
          VEIs_SV_contacts_restrict_bootstrap[i] <- 1-RRs_bootstrap/(1-VEs[1])
        }
        #print(i)
      }
  } 
}   
  
    VEIs_SV_cluster_bootstrap %>%
      as.data.frame() %>%
      setNames("Estimate") %>%
      add_column(j=x) %>%
      add_column(method="SV_cluster") -> VEIs_SV_cluster_bootstrap
    
    VEIs_SV_contacts_restrict_bootstrap %>%
      as.data.frame() %>%
      setNames("Estimate") %>%
      add_column(j=x) %>%
      add_column(method="SV_cluster_contacts") -> VEIs_SV_contacts_restrict_bootstrap

if (length(results_edge_list_analysis_hybrid)>2){
  # hybrid ----
  head(results_edge_list_analysis_hybrid)
  results_edge_list_analysis_hybrid2 <- results_edge_list_analysis_hybrid[,c("InfectedNode","node_stat","Infector","infector_stat","prob_SI_norm")]
  
  for (i in 1:nbootstrap){
    ID <- unique(results_edge_list_analysis_hybrid2$InfectedNode[!is.na(results_edge_list_analysis_hybrid2$InfectedNode)])
    ID_bootstrap <- sample(ID,length(ID),replace=TRUE)
    results_edge_list_analysis_hybrid2_bootstrap <- data.frame("InfectedNode"=NA,"node_stat"=NA,"Infector"=NA,"infector_stat"=NA,"prob_SI_norm"=NA)
    for (k in 1:length(ID_bootstrap)){
      results_edge_list_analysis_hybrid2_bootstrap <- 
        rbind(results_edge_list_analysis_hybrid2_bootstrap,results_edge_list_analysis_hybrid2[results_edge_list_analysis_hybrid2$InfectedNode==ID_bootstrap[k],])
    }  
    cluster <- results_edge_list_analysis_hybrid2_bootstrap
    names(cluster)[1] <- "indIDVar.2"
    names(cluster)[3] <- "indIDVar.1"
    cluster <- cluster[!is.na(cluster$indIDVar.2),]
    cluster$prob_SI_norm[is.na(cluster$prob_SI_norm)] <- 0
    cluster_sum <- aggregate(cluster$indIDVar.2,by=list(cluster$indIDVar.2),length) 
    if (mean(cluster_sum$x)>1){
      cluster_test_hc <- clusterInfectors(cluster, "indIDVar", "prob_SI_norm",
                                          clustMethod = "hc_absolute",
                                          0.2)
      cluster_test_hc %>%
        subset(cluster==1) %>%
        group_by(indIDVar.2) %>%
        mutate(prob_SI_norm_c = prob_SI_norm/sum(prob_SI_norm)) -> cluster_test_hc_norm
      
      vacc_inf_contacts_cluster1 <- sum(cluster_test_hc_norm$prob_SI_norm_c[cluster_test_hc_norm$infector_stat==1 & cluster_test_hc_norm$cluster==1],na.rm=TRUE)
      cont_inf_contacts_cluster1<- sum(cluster_test_hc_norm$prob_SI_norm_c[cluster_test_hc_norm$infector_stat==0  & cluster_test_hc_norm$cluster==1],na.rm=TRUE)
      if (length(vacc_inf_contacts_cluster1)==0 & length(cont_inf_contacts_cluster1)==0){
        RRs_bootstrap <- NA
      } else if (length(vacc_inf_contacts_cluster1)==0 & length(cont_inf_contacts_cluster1)>0){
        RRs_bootstrap <- 0
      } else if (length(vacc_inf_contacts_cluster1)>0 & length(cont_inf_contacts_cluster1)==0){
        RRs_bootstrap <- 1
      } else{
        RRs_bootstrap <- vacc_inf_contacts_cluster1/cont_inf_contacts_cluster1
      }   
    } else{
      RRs_bootstrap <- NA
    }
    if (length(RRs_bootstrap)>0){
      VEIs_hybrid_cluster_bootstrap[i] <- 1-RRs_bootstrap/(1-VEs[1])
    }
    #print(i)
  }
  
  
  # hybrid restrict contacts -----
  results_edge_list_analysis_hybrid_contacts <- results_edge_list_analysis_hybrid[results_edge_list_analysis_hybrid$contacts!=0 & !is.na(results_edge_list_analysis_hybrid$contacts),]
  results_edge_list_analysis_hybrid_contacts2 <- results_edge_list_analysis_hybrid_contacts[,c("InfectedNode","node_stat","Infector","infector_stat","prob_SI_norm")]
  
  for (i in 1:nbootstrap){
    ID <- unique(results_edge_list_analysis_hybrid_contacts2$InfectedNode[!is.na(results_edge_list_analysis_hybrid_contacts2$InfectedNode)])
    ID_bootstrap <- sample(ID,length(ID),replace=TRUE)
    results_edge_list_analysis_hybrid_contacts2_bootstrap <- data.frame("InfectedNode"=NA,"node_stat"=NA,"Infector"=NA,"infector_stat"=NA,"prob_SI_norm"=NA)
    for (k in 1:length(ID_bootstrap)){
      results_edge_list_analysis_hybrid_contacts2_bootstrap <- 
        rbind(results_edge_list_analysis_hybrid_contacts2_bootstrap,results_edge_list_analysis_hybrid_contacts2[results_edge_list_analysis_hybrid_contacts2$InfectedNode==ID_bootstrap[k],])
    }  
    cluster <- results_edge_list_analysis_hybrid_contacts2_bootstrap
    names(cluster)[1] <- "indIDVar.2"
    names(cluster)[3] <- "indIDVar.1"
    cluster <- cluster[!is.na(cluster$indIDVar.2),]
    cluster$prob_SI_norm[is.na(cluster$prob_SI_norm)] <- 0
    cluster_sum <- aggregate(cluster$indIDVar.2,by=list(cluster$indIDVar.2),length) 
    if (mean(cluster_sum$x)>1){
      cluster_test_hc <- clusterInfectors(cluster, "indIDVar", "prob_SI_norm",
                                          clustMethod = "hc_absolute",
                                          0.2)
      cluster_test_hc %>%
        subset(cluster==1) %>%
        group_by(indIDVar.2) %>%
        mutate(prob_SI_norm_c = prob_SI_norm/sum(prob_SI_norm)) -> cluster_test_hc_norm
      
      vacc_inf_contacts_cluster1 <- sum(cluster_test_hc_norm$prob_SI_norm_c[cluster_test_hc_norm$infector_stat==1 & cluster_test_hc_norm$cluster==1],na.rm=TRUE)
      cont_inf_contacts_cluster1<- sum(cluster_test_hc_norm$prob_SI_norm_c[cluster_test_hc_norm$infector_stat==0  & cluster_test_hc_norm$cluster==1],na.rm=TRUE)
      if (length(vacc_inf_contacts_cluster1)==0 & length(cont_inf_contacts_cluster1)==0){
        RRs_bootstrap <- NA
      } else if (length(vacc_inf_contacts_cluster1)==0 & length(cont_inf_contacts_cluster1)>0){
        RRs_bootstrap <- 0
      } else if (length(vacc_inf_contacts_cluster1)>0 & length(cont_inf_contacts_cluster1)==0){
        RRs_bootstrap <- 1
      } else{
        RRs_bootstrap <- vacc_inf_contacts_cluster1/cont_inf_contacts_cluster1
      }   
    } else{
      RRs_bootstrap <- NA
    }
    if (length(RRs_bootstrap)>0){
      VEIs_hybrid_contacts_restrict_bootstrap[i] <- 1-RRs_bootstrap/(1-VEs[1])
    }
    #print(i)
  }
}  
  VEIs_hybrid_cluster_bootstrap %>%
    as.data.frame() %>%
    setNames("Estimate") %>%
    add_column(j=x) %>%
    add_column(method="hybrid_cluster") -> VEIs_hybrid_cluster_bootstrap
  
  VEIs_hybrid_contacts_restrict_bootstrap %>%
    as.data.frame() %>%
    setNames("Estimate") %>%
    add_column(j=x) %>%
    add_column(method="hybrid_cluster_contacts") -> VEIs_hybrid_contacts_restrict_bootstrap
  

  
  
  Bootstrap_master <- bind_rows(Bootstrap_master,VEIs_bootstrap,
                                VEIs_contact_bootstrap_prob,VEIs_contact_bootstrap_cluster2,
                                VEIs_SV_cluster_bootstrap,VEIs_SV_contacts_restrict_bootstrap,
                                VEIs_hybrid_cluster_bootstrap,VEIs_hybrid_contacts_restrict_bootstrap)
  
  write.csv(Bootstrap_master,"Bootstrap_master.csv")

}

