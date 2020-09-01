# Code for estimating VEI from sequence and contact tracing data

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
#### load files -----
js <- 1:500
for (x in js){
  cat(x,"\n")
  
  load(paste0(x,"_g",R0,"_",infect_VE,"_",direct_VE,"_",mut.rate,"_",bn,"_",infperiod_shape,"_",sample_percent,"_","sample.Rdata"))
  results_edge_list_analysis <- read.csv(paste0(x,"_results_edge_list_analysis",R0,"_",infect_VE,"_",direct_VE,"_",mut.rate,"_",bn,"_",infperiod_shape,"_",sample_percent,"_","sample_new.csv"))

  results_contacts <- read.csv(paste0(x,"_results_contacts",R0,"_",infect_VE,"_",direct_VE,"_",mut.rate,"_",bn,"_",infperiod_shape,"_",sample_percent,"_","sample_new.csv"))
  results_edge_list_analysis_DS_scoremat <- read.csv(paste0(x,"_results_edge_list_analysis_DS_scoremat",R0,"_",infect_VE,"_",direct_VE,"_",mut.rate,"_",bn,"_",infperiod_shape,"_",sample_percent,"_","sample_new.csv"))
  results_edge_list_analysis_DS_hybrid <- read.csv(paste0(x,"_results_edge_list_analysis_DS_hybrid",R0,"_",infect_VE,"_",direct_VE,"_",mut.rate,"_",bn,"_",infperiod_shape,"_",sample_percent,"_","sample_new.csv"))
  results_edge_list_analysis_DS_dist <- read.csv(paste0(x,"_results_edge_list_analysis_DS_dist",R0,"_",infect_VE,"_",direct_VE,"_",mut.rate,"_",bn,"_",infperiod_shape,"_",sample_percent,"_","sample_new.csv"))
  results_edge_list_analysis_W <- read.csv(paste0(x,"_results_edge_list_analysis_DS_W",R0,"_",infect_VE,"_",direct_VE,"_",mut.rate,"_",bn,"_",infperiod_shape,"_",sample_percent,"_","sample_new.csv"))
  
  for (df in 1:4){
    print(df)
    if (df==1){
      results_edge_list_analysis_DS <- results_edge_list_analysis_DS_scoremat
    } else if (df==2){
      results_edge_list_analysis_DS <- results_edge_list_analysis_DS_hybrid
    } else if (df==3){
      results_edge_list_analysis_DS <- results_edge_list_analysis_DS_dist
    } else if (df==4){
      results_edge_list_analysis_DS <- results_edge_list_analysis_W
    } 
     
    
  if (nrow(results_edge_list_analysis_DS)>1){
    
    #1# look only at most likely
    summary_infectors_DS_max <- aggregate(results_edge_list_analysis_DS$max_prob_SI_norm,by=list(results_edge_list_analysis_DS$infector_stat),sum,na.rm=TRUE)
    names(summary_infectors_DS_max)<-c("infector_stat","number")
    vacc_inf_DS_max <- summary_infectors_DS_max$number[summary_infectors_DS_max$infector_stat==1]
    cont_inf_DS_max<- summary_infectors_DS_max$number[summary_infectors_DS_max$infector_stat==0]
    if (length(vacc_inf_DS_max)==0 & length(cont_inf_DS_max)==0){
      RR_inf_DS_max <- NA
    } else if (length(vacc_inf_DS_max)==0 & length(cont_inf_DS_max)>0){
      RR_inf_DS_max <- 0
    } else if (length(vacc_inf_DS_max)>0 & length(cont_inf_DS_max)==0){
      RR_inf_DS_max <- 1
    } else{
      RR_inf_DS_max <- vacc_inf_DS_max/cont_inf_DS_max
    }    
    
    #2# aggregate by probability instead of number of infectors
    summary_infectors_DS_prob <- aggregate(results_edge_list_analysis_DS$prob_SI_norm,by=list(results_edge_list_analysis_DS$infector_stat),sum,na.rm=TRUE)
    names(summary_infectors_DS_prob)<-c("infector_stat","prob_SI_norm")
    vacc_inf_DS_prob <- summary_infectors_DS_prob$prob_SI_norm[summary_infectors_DS_prob$infector_stat==1]
    cont_inf_DS_prob<- summary_infectors_DS_prob$prob_SI_norm[summary_infectors_DS_prob$infector_stat==0]
    if (length(vacc_inf_DS_prob)==0 & length(cont_inf_DS_prob)==0){
      RR_inf_DS_prob <- NA
    } else if (length(vacc_inf_DS_prob)==0 & length(cont_inf_DS_prob)>0){
      RR_inf_DS_prob <- 0
    } else if (length(vacc_inf_DS_prob)>0 & length(cont_inf_DS_prob)==0){
      RR_inf_DS_prob <- 1
    } else{
      RR_inf_DS_prob <- vacc_inf_DS_prob/cont_inf_DS_prob
    }   
    
    
    #3# clustering algorithm 0.1
    cluster <- results_edge_list_analysis_DS[,2:ncol(results_edge_list_analysis_DS)]
    names(cluster)[1] <- "indIDVar.2"
    names(cluster)[3] <- "indIDVar.1"
    cluster$prob_SI_norm[is.na(cluster$prob_SI_norm)] <- 0
    cluster_sum <- aggregate(cluster$indIDVar.2,by=list(cluster$indIDVar.2),length) 
    if (mean(cluster_sum$x)>1){
      cluster_test_hc <- clusterInfectors(cluster, "indIDVar", "prob_SI_norm",
                                          clustMethod = "hc_absolute",
                                          0.1)
      cluster_test_hc %>%
        subset(cluster==1) %>%
        group_by(indIDVar.2) %>%
        mutate(prob_SI_norm_c = prob_SI_norm/sum(prob_SI_norm)) -> cluster_test_hc_norm
      
      vacc_inf_DS_cluster1 <- sum(cluster_test_hc_norm$prob_SI_norm_c[cluster_test_hc_norm$infector_stat==1 & cluster_test_hc_norm$cluster==1],na.rm=TRUE)
      cont_inf_DS_cluster1<- sum(cluster_test_hc_norm$prob_SI_norm_c[cluster_test_hc_norm$infector_stat==0  & cluster_test_hc_norm$cluster==1],na.rm=TRUE)
      if (length(vacc_inf_DS_cluster1)==0 & length(cont_inf_DS_cluster1)==0){
        RR_inf_DS_cluster1 <- NA
      } else if (length(vacc_inf_DS_cluster1)==0 & length(cont_inf_DS_cluster1)>0){
        RR_inf_DS_cluster1 <- 0
      } else if (length(vacc_inf_DS_cluster1)>0 & length(cont_inf_DS_cluster1)==0){
        RR_inf_DS_cluster1 <- 1
      } else{
        RR_inf_DS_cluster1 <- vacc_inf_DS_cluster1/cont_inf_DS_cluster1
      }   
    } else{
      RR_inf_DS_cluster1 <- NA
    }
    
    #4# clustering algorithm 0.2
    cluster <- results_edge_list_analysis_DS[,2:ncol(results_edge_list_analysis_DS)]
    names(cluster)[1] <- "indIDVar.2"
    names(cluster)[3] <- "indIDVar.1"
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
      
      vacc_inf_DS_cluster2 <- sum(cluster_test_hc_norm$prob_SI_norm_c[cluster_test_hc_norm$infector_stat==1 & cluster_test_hc_norm$cluster==1],na.rm=TRUE)
      cont_inf_DS_cluster2<- sum(cluster_test_hc_norm$prob_SI_norm_c[cluster_test_hc_norm$infector_stat==0  & cluster_test_hc_norm$cluster==1],na.rm=TRUE)
    if (length(vacc_inf_DS_cluster2)==0 & length(cont_inf_DS_cluster2)==0){
        RR_inf_DS_cluster2 <- NA
      } else if (length(vacc_inf_DS_cluster2)==0 & length(cont_inf_DS_cluster2)>0){
        RR_inf_DS_cluster2 <- 0
      } else if (length(vacc_inf_DS_cluster2)>0 & length(cont_inf_DS_cluster2)==0){
        RR_inf_DS_cluster2 <- 1
      } else{
        RR_inf_DS_cluster2 <- vacc_inf_DS_cluster2/cont_inf_DS_cluster2
      }   
    } else{
      RR_inf_DS_cluster2 <- NA
    }
    
    
    #5 -6 contacts
    results_edge_list_analysis_DS_contacts <- results_edge_list_analysis_DS[results_edge_list_analysis_DS$contacts!=0 & !is.na(results_edge_list_analysis_DS$contacts),]
    if (nrow(results_edge_list_analysis_DS_contacts)>0){
      #7# look at maximum only among contacts -- max
      ID_contacts <- unique(results_edge_list_analysis_DS_contacts$InfectedNode)
      for (i in 1:length(ID_contacts)){
        num_max_contacts <- length(which(results_edge_list_analysis_DS_contacts$prob_SI_norm[results_edge_list_analysis_DS_contacts$InfectedNode==ID_contacts[i]]==max(results_edge_list_analysis_DS_contacts$prob_SI_norm[results_edge_list_analysis_DS_contacts$InfectedNode==ID_contacts[i]],na.rm=TRUE)))
        results_edge_list_analysis_DS_contacts$max_prob_SI_norm[results_edge_list_analysis_DS_contacts$InfectedNode==ID_contacts[i] & results_edge_list_analysis_DS_contacts$prob_SI_norm == max(results_edge_list_analysis_DS_contacts$prob_SI_norm[results_edge_list_analysis_DS_contacts$InfectedNode==ID_contacts[i]],na.rm=TRUE)] <-1/num_max_contacts
      }
      summary_infectors_DS_contacts <- aggregate(results_edge_list_analysis_DS_contacts$max_prob_SI_norm,by=list(results_edge_list_analysis_DS_contacts$infector_stat),sum,na.rm=TRUE)
      names(summary_infectors_DS_contacts)<-c("infector_stat","number")
      vacc_inf_DS_contacts<- summary_infectors_DS_contacts$number[summary_infectors_DS_contacts$infector_stat==1]
      cont_inf_DS_contacts <- summary_infectors_DS_contacts$number[summary_infectors_DS_contacts$infector_stat==0]
      if (length(vacc_inf_DS_contacts)==0 & length(cont_inf_DS_contacts)==0){
        RR_inf_DS_contacts <- NA
      } else if (length(vacc_inf_DS_contacts)==0 & length(cont_inf_DS_contacts)>0){
        RR_inf_DS_contacts <- 0
      } else if (length(vacc_inf_DS_contacts)>0 & length(cont_inf_DS_contacts)==0){
        RR_inf_DS_contacts <- 1
      } else{
        RR_inf_DS_contacts <- vacc_inf_DS_contacts/cont_inf_DS_contacts
      } 
      
      #8# look only among contacts -- weighted
      ID_contacts <- unique(results_edge_list_analysis_DS_contacts$InfectedNode)
      # normalize
      for (i in 1:length(ID_contacts)){
        results_edge_list_analysis_DS_contacts$prob_SI_norm[results_edge_list_analysis_DS_contacts$InfectedNode==ID_contacts[i]] <- results_edge_list_analysis_DS_contacts$prob_SI_norm[results_edge_list_analysis_DS_contacts$InfectedNode==ID_contacts[i]] / sum(results_edge_list_analysis_DS_contacts$prob_SI_norm[results_edge_list_analysis_DS_contacts$InfectedNode==ID_contacts[i]])
      }
      summary_infectors_DS_contacts_prob <- aggregate(results_edge_list_analysis_DS_contacts$prob_SI_norm,by=list(results_edge_list_analysis_DS_contacts$infector_stat),sum,na.rm=TRUE)
      names(summary_infectors_DS_contacts_prob)<-c("infector_stat","number")
      vacc_inf_DS_contacts_prob<- summary_infectors_DS_contacts_prob$number[summary_infectors_DS_contacts_prob$infector_stat==1]
      cont_inf_DS_contacts_prob <- summary_infectors_DS_contacts_prob$number[summary_infectors_DS_contacts_prob$infector_stat==0]
      if (length(vacc_inf_DS_contacts_prob)==0 & length(cont_inf_DS_contacts_prob)==0){
        RR_inf_DS_contacts_prob <- NA
      } else if (length(vacc_inf_DS_contacts_prob)==0 & length(cont_inf_DS_contacts_prob)>0){
        RR_inf_DS_contacts_prob <- 0
      } else if (length(vacc_inf_DS_contacts_prob)>0 & length(cont_inf_DS_contacts_prob)==0){
        RR_inf_DS_contacts_prob <- 1
      } else{
        RR_inf_DS_contacts_prob <- vacc_inf_DS_contacts_prob/cont_inf_DS_contacts_prob
      } 
      
      #9# clustering algorithm 0.1
      cluster <- results_edge_list_analysis_DS_contacts[,2:ncol(results_edge_list_analysis_DS_contacts)]
      names(cluster)[1] <- "indIDVar.2"
      names(cluster)[3] <- "indIDVar.1"
      cluster$prob_SI_norm[is.na(cluster$prob_SI_norm)] <- 0
      cluster_sum <- aggregate(cluster$indIDVar.2,by=list(cluster$indIDVar.2),length) 
      if (mean(cluster_sum$x)>1){
        cluster_test_hc <- clusterInfectors(cluster, "indIDVar", "prob_SI_norm",
                                            clustMethod = "hc_absolute",
                                            0.1)
        cluster_test_hc %>%
          subset(cluster==1) %>%
          group_by(indIDVar.2) %>%
          mutate(prob_SI_norm_c = prob_SI_norm/sum(prob_SI_norm)) -> cluster_test_hc_norm
        
        vacc_inf_DS_cluster1 <- sum(cluster_test_hc_norm$prob_SI_norm_c[cluster_test_hc_norm$infector_stat==1 & cluster_test_hc_norm$cluster==1],na.rm=TRUE)
        cont_inf_DS_cluster1<- sum(cluster_test_hc_norm$prob_SI_norm_c[cluster_test_hc_norm$infector_stat==0  & cluster_test_hc_norm$cluster==1],na.rm=TRUE)
        if (length(vacc_inf_DS_cluster1)==0 & length(cont_inf_DS_cluster1)==0){
          RR_inf_DS_contacts_cluster1 <- NA
        } else if (length(vacc_inf_DS_cluster1)==0 & length(cont_inf_DS_cluster1)>0){
          RR_inf_DS_contacts_cluster1 <- 0
        } else if (length(vacc_inf_DS_cluster1)>0 & length(cont_inf_DS_cluster1)==0){
          RR_inf_DS_contacts_cluster1 <- 1
        } else{
          RR_inf_DS_contacts_cluster1 <- vacc_inf_DS_cluster1/cont_inf_DS_cluster1
        }   
      } else{
        RR_inf_DS_contacts_cluster1 <- NA
      }
      
      #10# clustering algorithm 0.2
      cluster <- results_edge_list_analysis_DS_contacts[,2:ncol(results_edge_list_analysis_DS_contacts)]
      names(cluster)[1] <- "indIDVar.2"
      names(cluster)[3] <- "indIDVar.1"
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
        
        vacc_inf_DS_cluster1 <- sum(cluster_test_hc_norm$prob_SI_norm_c[cluster_test_hc_norm$infector_stat==1 & cluster_test_hc_norm$cluster==1],na.rm=TRUE)
        cont_inf_DS_cluster1<- sum(cluster_test_hc_norm$prob_SI_norm_c[cluster_test_hc_norm$infector_stat==0  & cluster_test_hc_norm$cluster==1],na.rm=TRUE)
        if (length(vacc_inf_DS_cluster1)==0 & length(cont_inf_DS_cluster1)==0){
          RR_inf_DS_contacts_cluster2 <- NA
        } else if (length(vacc_inf_DS_cluster1)==0 & length(cont_inf_DS_cluster1)>0){
          RR_inf_DS_contacts_cluster2 <- 0
        } else if (length(vacc_inf_DS_cluster1)>0 & length(cont_inf_DS_cluster1)==0){
          RR_inf_DS_contacts_cluster2 <- 1
        } else{
          RR_inf_DS_contacts_cluster2 <- vacc_inf_DS_cluster1/cont_inf_DS_cluster1
        }   
      } else{
        RR_inf_DS_contacts_cluster2 <- NA
      }
      
      
      
    } else{
      RR_inf_DS_contacts <- NA
      RR_inf_DS_contacts_prob <- NA
      RR_inf_DS_contacts_cluster1 <- NA
      RR_inf_DS_contacts_cluster2 <- NA
    }
  } else{
    RR_inf_DS_max <- NA
    RR_inf_DS_prob <- NA
    RR_inf_DS_contacts <- NA
    RR_inf_DS_contacts_prob <- NA
    RR_inf_DS_cluster1 <- NA
    RR_inf_DS_cluster2 <- NA
    RR_inf_DS_contacts_cluster1 <- NA
    RR_inf_DS_contacts_cluster2 <- NA
    results_edge_list_analysis_DS <- NA
  }
  
  
  
  if (df==1){
    RR_inf_DS_scoremat_max <- RR_inf_DS_max
    RR_inf_DS_scoremat_prob <- RR_inf_DS_prob
    RR_inf_DS_scoremat_contacts <- RR_inf_DS_contacts
    RR_inf_DS_scoremat_contacts_prob <- RR_inf_DS_contacts_prob
    RR_inf_DS_scoremat_cluster1 <- RR_inf_DS_cluster1
    RR_inf_DS_scoremat_cluster2 <- RR_inf_DS_cluster2
    RR_inf_DS_scoremat_contacts_cluster1 <- RR_inf_DS_contacts_cluster1
    RR_inf_DS_scoremat_contacts_cluster2 <- RR_inf_DS_contacts_cluster2
    results_edge_list_analysis_DS_scoremat <- results_edge_list_analysis_DS
  } else if (df==2){
    RR_inf_DS_hybrid_max <- RR_inf_DS_max
    RR_inf_DS_hybrid_prob <- RR_inf_DS_prob
    RR_inf_DS_hybrid_contacts <- RR_inf_DS_contacts
    RR_inf_DS_hybrid_contacts_prob <- RR_inf_DS_contacts_prob
    RR_inf_DS_hybrid_cluster1 <- RR_inf_DS_cluster1
    RR_inf_DS_hybrid_cluster2 <- RR_inf_DS_cluster2
    RR_inf_DS_hybrid_contacts_cluster1 <- RR_inf_DS_contacts_cluster1
    RR_inf_DS_hybrid_contacts_cluster2 <- RR_inf_DS_contacts_cluster2
    results_edge_list_analysis_DS_hybrid <- results_edge_list_analysis_DS
  } else if (df==3){
    RR_inf_DS_dist_max <- RR_inf_DS_max
    RR_inf_DS_dist_prob <- RR_inf_DS_prob
    RR_inf_DS_dist_contacts <- RR_inf_DS_contacts
    RR_inf_DS_dist_contacts_prob <- RR_inf_DS_contacts_prob
    RR_inf_DS_dist_cluster1 <- RR_inf_DS_cluster1
    RR_inf_DS_dist_cluster2 <- RR_inf_DS_cluster2
    RR_inf_DS_dist_contacts_cluster1 <- RR_inf_DS_contacts_cluster1
    RR_inf_DS_dist_contacts_cluster2 <- RR_inf_DS_contacts_cluster2
    results_edge_list_analysis_DS_dist <- results_edge_list_analysis_DS
  } else if (df==4){
    RR_inf_W_max <- RR_inf_DS_max
    RR_inf_W_prob <- RR_inf_DS_prob
    RR_inf_W_contacts <- RR_inf_DS_contacts
    RR_inf_W_contacts_prob <- RR_inf_DS_contacts_prob
    RR_inf_W_cluster1 <- RR_inf_DS_cluster1
    RR_inf_W_cluster2 <- RR_inf_DS_cluster2
    RR_inf_W_contacts_cluster1 <- RR_inf_DS_contacts_cluster1
    RR_inf_W_contacts_cluster2 <- RR_inf_DS_contacts_cluster2
    results_edge_list_analysis_W <- results_edge_list_analysis_DS
  } 
}  
  # # VEs
  VEs <- coxmodel(results_edge_list_analysis)[1]
  
  results_infected <- results_edge_list_analysis[results_edge_list_analysis$eventstatus==1,]
  # for those who remain infected after end of the trial, make recovered day after trial ends
  results_infected$DayRecovered[is.na(results_infected$DayRecovered)] <- trial_startday + trial_length
  
  # VEI
  summary_infectors <- aggregate(results_infected$InfectedNode,by=list(results_infected$infector_stat),length)
  names(summary_infectors)<-c("infector_stat","number")
  vacc_inf <- summary_infectors$number[summary_infectors$infector_stat==1]
  cont_inf<- summary_infectors$number[summary_infectors$infector_stat==0]
  if (length(vacc_inf)==0 & length(cont_inf)==0){
    RR_inf <- NA
  } else if (length(vacc_inf)==0 & length(cont_inf)>0){
    RR_inf <- 0
  } else if (length(vacc_inf)>0 & length(cont_inf)==0){
    RR_inf <- 1
  } else{
    RR_inf <- vacc_inf/cont_inf
  } 
  
  summary_infectors_contacts <- aggregate(results_contacts$InfectedNode,by=list(results_contacts$infector_stat),length)
  names(summary_infectors_contacts)<-c("infector_stat","number")
  vacc_inf_contacts<- summary_infectors_contacts$number[summary_infectors_contacts$infector_stat==1]
  cont_inf_contacts <- summary_infectors_contacts$number[summary_infectors_contacts$infector_stat==0]
  if (length(vacc_inf_contacts)==0 & length(cont_inf_contacts)==0){
    RR_inf_contacts <- NA
  } else if (length(vacc_inf_contacts)==0 & length(cont_inf_contacts)>0){
    RR_inf_contacts <- 0
  } else if (length(vacc_inf_contacts)>0 & length(cont_inf_contacts)==0){
    RR_inf_contacts <- 1
  } else{
    RR_inf_contacts <- vacc_inf_contacts/cont_inf_contacts
  } 
  
  results <- merge(results_edge_list_analysis[,c("InfectedNode","DayInfected","node_stat")],results_contacts,by="InfectedNode")
  results_contacts2 <- merge(results,results_edge_list_analysis[,c("InfectedNode","DayInfected")],by.x="Contact",by.y="InfectedNode",)
  
  for (i in 1:nrow(results_contacts2)){
    #cat(i,"\n")
    if(results_contacts2$DayInfected.x[i] > results_contacts2$DayInfected.y[i]){
      results_contacts2$SI_prob[i] <- distribution_norm[results_contacts2$DayInfected.x[i] - results_contacts2$DayInfected.y[i]]
    }
  }

  results_contacts2 %>%
    group_by(InfectedNode) %>%
    mutate(SI_prob_norm = case_when(sum(SI_prob)>0~SI_prob/sum(SI_prob),
                                    TRUE~0),
           max = case_when(SI_prob_norm == max(SI_prob_norm) ~ 1,
                           TRUE~0),
           max2 = max/sum(max)) -> results_contacts2
  
  summary_infectors_contacts <- aggregate(results_contacts2$SI_prob_norm,by=list(results_contacts2$infector_stat),sum)
  names(summary_infectors_contacts)<-c("infector_stat","number")
  vacc_inf_contacts<- summary_infectors_contacts$number[summary_infectors_contacts$infector_stat==1]
  cont_inf_contacts <- summary_infectors_contacts$number[summary_infectors_contacts$infector_stat==0]
  if (length(vacc_inf_contacts)==0 & length(cont_inf_contacts)==0){
    RR_inf_contacts_prob <- NA
  } else if (length(vacc_inf_contacts)==0 & length(cont_inf_contacts)>0){
    RR_inf_contacts_prob <- 0
  } else if (length(vacc_inf_contacts)>0 & length(cont_inf_contacts)==0){
    RR_inf_contacts_prob <- 1
  } else{
    RR_inf_contacts_prob <- vacc_inf_contacts/cont_inf_contacts
  } 
  
  summary_infectors_contacts_max <- aggregate(results_contacts2$max2,by=list(results_contacts2$infector_stat),sum)
  names(summary_infectors_contacts_max)<-c("infector_stat","number")
  vacc_inf_contacts_max<- summary_infectors_contacts_max$number[summary_infectors_contacts_max$infector_stat==1]
  cont_inf_contacts_max <- summary_infectors_contacts_max$number[summary_infectors_contacts_max$infector_stat==0]
  if (length(vacc_inf_contacts_max)==0 & length(cont_inf_contacts_max)==0){
    RR_inf_contacts_max <- NA
  } else if (length(vacc_inf_contacts_max)==0 & length(cont_inf_contacts_max)>0){
    RR_inf_contacts_max <- 0
  } else if (length(vacc_inf_contacts_max)>0 & length(cont_inf_contacts_max)==0){
    RR_inf_contacts_max <- 1
  } else{
    RR_inf_contacts_max <- vacc_inf_contacts_max/cont_inf_contacts_max
  } 
  
  #9# clustering algorithm 0.1
  cluster <- results_contacts2
  names(cluster)[2] <- "indIDVar.2"
  names(cluster)[1] <- "indIDVar.1"
  cluster$SI_prob_norm[is.na(cluster$SI_prob_norm)] <- 0
  cluster_sum <- aggregate(cluster$indIDVar.2,by=list(cluster$indIDVar.2),length) 
  if (mean(cluster_sum$x)>1){
    cluster_test_hc <- clusterInfectors(cluster, "indIDVar", "SI_prob_norm",
                                        clustMethod = "hc_absolute",
                                        0.1)
    cluster_test_hc %>%
      subset(cluster==1) %>%
      group_by(indIDVar.2) %>%
      mutate(prob_SI_norm_c = SI_prob_norm/sum(SI_prob_norm)) -> cluster_test_hc_norm
    
    vacc_inf_contacts_cluster1 <- sum(cluster_test_hc_norm$prob_SI_norm_c[cluster_test_hc_norm$infector_stat==1 & cluster_test_hc_norm$cluster==1],na.rm=TRUE)
    cont_inf_contacts_cluster1<- sum(cluster_test_hc_norm$prob_SI_norm_c[cluster_test_hc_norm$infector_stat==0  & cluster_test_hc_norm$cluster==1],na.rm=TRUE)
    if (length(vacc_inf_contacts_cluster1)==0 & length(cont_inf_contacts_cluster1)==0){
      RR_inf_contacts_cluster1 <- NA
    } else if (length(vacc_inf_contacts_cluster1)==0 & length(cont_inf_contacts_cluster1)>0){
      RR_inf_contacts_cluster1 <- 0
    } else if (length(vacc_inf_contacts_cluster1)>0 & length(cont_inf_contacts_cluster1)==0){
      RR_inf_contacts_cluster1 <- 1
    } else{
      RR_inf_contacts_cluster1 <- vacc_inf_contacts_cluster1/cont_inf_contacts_cluster1
    }   
  } else{
    RR_inf_contacts_cluster1 <- NA
  }
  
  #10# clustering algorithm 0.2
  cluster <- results_contacts2
  names(cluster)[2] <- "indIDVar.2"
  names(cluster)[1] <- "indIDVar.1"
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
    
    vacc_inf_contacts_cluster1 <- sum(cluster_test_hc_norm$prob_SI_norm_c[cluster_test_hc_norm$infector_stat==1 & cluster_test_hc_norm$cluster==1],na.rm=TRUE)
    cont_inf_contacts_cluster1<- sum(cluster_test_hc_norm$prob_SI_norm_c[cluster_test_hc_norm$infector_stat==0  & cluster_test_hc_norm$cluster==1],na.rm=TRUE)
    if (length(vacc_inf_contacts_cluster1)==0 & length(cont_inf_contacts_cluster1)==0){
      RR_inf_contacts_cluster2 <- NA
    } else if (length(vacc_inf_contacts_cluster1)==0 & length(cont_inf_contacts_cluster1)>0){
      RR_inf_contacts_cluster2 <- 0
    } else if (length(vacc_inf_contacts_cluster1)>0 & length(cont_inf_contacts_cluster1)==0){
      RR_inf_contacts_cluster2 <- 1
    } else{
      RR_inf_contacts_cluster2 <- vacc_inf_contacts_cluster1/cont_inf_contacts_cluster1
    }   
  } else{
    RR_inf_contacts_cluster2 <- NA
  }
  
  VEIs <- 1-RR_inf/(1-VEs) 
  VEIs_contacts <- 1-RR_inf_contacts/(1-VEs)
  VEIs_contacts_prob <- 1-RR_inf_contacts_prob/(1-VEs)
  VEIs_contacts_max <- 1-RR_inf_contacts_max/(1-VEs)
  VEIs_contacts_cluster1 <- 1-RR_inf_contacts_cluster1/(1-VEs)
  VEIs_contacts_cluster2 <- 1-RR_inf_contacts_cluster2/(1-VEs)
  
  VEIs_W_max <-1-RR_inf_W_max/(1-VEs)
  VEIs_W_prob <-1-RR_inf_W_prob/(1-VEs)
  VEIs_W_contacts<-1-RR_inf_W_contacts/(1-VEs)
  VEIs_W_contacts_prob<-1-RR_inf_W_contacts_prob/(1-VEs)
  VEIs_W_cluster1 <- 1-RR_inf_W_cluster1/(1-VEs)
  VEIs_W_cluster2 <- 1-RR_inf_W_cluster2/(1-VEs)
  VEIs_W_contacts_cluster1 <- 1-RR_inf_W_contacts_cluster1/(1-VEs)
  VEIs_W_contacts_cluster2 <- 1-RR_inf_W_contacts_cluster2/(1-VEs)
  
  VEIs_scoremat_max <-1-RR_inf_DS_scoremat_max/(1-VEs)
  VEIs_scoremat_prob <-1-RR_inf_DS_scoremat_prob/(1-VEs)
  VEIs_scoremat_contacts <- 1-RR_inf_DS_scoremat_contacts/(1-VEs)
  VEIs_scoremat_contacts_prob <- 1-RR_inf_DS_scoremat_contacts_prob/(1-VEs)
  VEIs_scoremat_cluster1 <- 1-RR_inf_DS_scoremat_cluster1/(1-VEs)
  VEIs_scoremat_cluster2 <- 1- RR_inf_DS_scoremat_cluster2/(1-VEs)
  VEIs_scoremat_contacts_cluster1 <- 1-RR_inf_DS_scoremat_contacts_cluster1/(1-VEs)
  VEIs_scoremat_contacts_cluster2 <- 1- RR_inf_DS_scoremat_contacts_cluster2/(1-VEs)
  
  VEIs_hybrid_max <-1-RR_inf_DS_hybrid_max/(1-VEs)
  VEIs_hybrid_prob <-1-RR_inf_DS_hybrid_prob/(1-VEs)
  VEIs_hybrid_contacts <- 1-RR_inf_DS_hybrid_contacts/(1-VEs)
  VEIs_hybrid_contacts_prob <- 1-RR_inf_DS_hybrid_contacts_prob/(1-VEs)
  VEIs_hybrid_cluster1 <- 1-RR_inf_DS_hybrid_cluster1/(1-VEs)
  VEIs_hybrid_cluster2 <- 1- RR_inf_DS_hybrid_cluster2/(1-VEs)
  VEIs_hybrid_contacts_cluster1 <- 1-RR_inf_DS_hybrid_contacts_cluster1/(1-VEs)
  VEIs_hybrid_contacts_cluster2 <- 1- RR_inf_DS_hybrid_contacts_cluster2/(1-VEs)
  
  VEIs_dist_max <-1-RR_inf_DS_dist_max/(1-VEs)
  VEIs_dist_prob <-1-RR_inf_DS_dist_prob/(1-VEs)
  VEIs_dist_contacts <- 1-RR_inf_DS_dist_contacts/(1-VEs)
  VEIs_dist_contacts_prob <- 1-RR_inf_DS_dist_contacts_prob/(1-VEs)
  VEIs_dist_cluster1 <- 1-RR_inf_DS_dist_cluster1/(1-VEs)
  VEIs_dist_cluster2 <- 1- RR_inf_DS_dist_cluster2/(1-VEs)
  VEIs_dist_contacts_cluster1 <- 1-RR_inf_DS_dist_contacts_cluster1/(1-VEs)
  VEIs_dist_contacts_cluster2 <- 1- RR_inf_DS_dist_contacts_cluster2/(1-VEs)
  
  
  args=(commandArgs(TRUE))
  for (i in 1:length(args)) {
    eval (parse (text = args[[i]] ))
  }
  
  genome_results <- as.data.frame(cbind(x,infect_VE,direct_VE,R0,mut.rate,bn,numevents_cont,numevents_vacc,
                                        VEs,VEIs,VEIs_contacts,VEIs_contacts_prob,VEIs_contacts_max,VEIs_contacts_cluster1,VEIs_contacts_cluster2,
                                        VEIs_W_max,VEIs_W_prob,VEIs_W_contacts,VEIs_W_contacts_prob,VEIs_W_cluster1,VEIs_W_cluster2,VEIs_W_contacts_cluster1,VEIs_W_contacts_cluster2,
                                        VEIs_scoremat_max,VEIs_scoremat_prob,VEIs_scoremat_contacts,VEIs_scoremat_contacts_prob,VEIs_scoremat_cluster1,VEIs_scoremat_cluster2,VEIs_scoremat_contacts_cluster1,VEIs_scoremat_contacts_cluster2,
                                        VEIs_hybrid_max,VEIs_hybrid_prob,VEIs_hybrid_contacts,VEIs_hybrid_contacts_prob,VEIs_hybrid_cluster1,VEIs_hybrid_cluster2,VEIs_hybrid_contacts_cluster1,VEIs_hybrid_contacts_cluster2,
                                        VEIs_dist_max,VEIs_dist_prob,VEIs_dist_contacts,VEIs_dist_contacts_prob,VEIs_dist_cluster1,VEIs_dist_cluster2,VEIs_dist_contacts_cluster1,VEIs_dist_contacts_cluster2)) 
  
write.csv(genome_results,paste0(x,'_',R0,'_',infect_VE,"_",direct_VE,"_",mut.rate,"_",bn,"_",infperiod_shape,"_",sample_percent,'_genome_results_new2.csv'))
write.csv(results_contacts2, paste0(x,"_results_contacts2",R0,"_",infect_VE,"_",direct_VE,"_",mut.rate,"_",bn,"_",infperiod_shape,"_",sample_percent,"_","sample_new.csv"))
}
