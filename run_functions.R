# run functions and produce sequence and contact tracing data
# network code written by Matt Hitchings and adapted by Rebecca Kahn
# simulated sequence data (using adpated version of seedy by Colin Worby)

#Parameters --------
setwd("/n/holyscratch01/lipsitch_lab/rkahn/ADAGIO")

source('./simfixoutbreak.R')
source('./make_network.R')
source('./network_epidemic.R')
source('./analyze_data.R')
source('./clustering algorithm.R')


require(ggplot2)
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
mut.rate <- 0.012 #0.003
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
####

# end functions and beginning of runs #-----

g<-make_network(ave_community_size, community_size_range, num_communities,rate_within, rate_between)

list[results,results_edge_list,trial_nodes,edge_list,e_nodes_master,g]<-
  network_epidemic(g,beta,num_introductions,direct_VE,infect_VE,edge_list,
                   incperiod_shape,incperiod_rate,infperiod_shape,infperiod_rate,1,0,
                   trial_startday,trial_length,num_clusters_enrolled_per_day,
                   enrollment_period,cluster_coverage,sim)

save(g,file=paste0(j,"_g",R0,"_",infect_VE,"_",direct_VE,"_",mut.rate,"_",bn,"_",infperiod_shape,"_",sample_percent,"_","sample.Rdata"))

list[VE,pval,events_vacc,events_cont,analysed_trialsize,
     RR_inf,mean_R,med_R,vacc_vacc_inf,vacc_cont_inf,cont_vacc_inf,cont_cont_inf,importation,results_edge_list_analysis,results_infected]<-
  analyse_data(results_edge_list,trial_nodes,trial_startday,trial_length,ave_inc_period,num_clusters_perarm,0)

write.csv(results_edge_list_analysis, paste0(j,"_results_edge_list_analysis",R0,"_",infect_VE,"_",direct_VE,"_",mut.rate,"_",bn,"_",infperiod_shape,"_",sample_percent,"_","sample_new.csv"))


############

# simulate pathogen evolution (consensus sequence)

ID <- results_infected$InfectedNode
ID
exp.times <- results_infected$DayExposed
inf.times <- results_infected$DayInfected
rec.times <-results_infected$DayRecovered + 1
truesource <- results_infected$Infector

# consensus sequences
W <- simfixoutbreak2(ID, exp.times, inf.times, rec.times,truesource, mut.rate=mut.rate,equi.pop=5000, shape=flat,samples.per.time = 1,
                     samp.schedule="random",samp.freq=10,inoc.size=bn,imp.var=25,
                     full=FALSE, feedback=1, glen=gen_len,
                     ref.strain=NULL)

# Assess transmission routes
GD <- gd(W$sampledata[,3], W$libr, W$nuc, W$librstrains)

ID <- W$sampledata[,1]
sample.times <- numeric(length(ID))
inf.times <- numeric(length(ID))
rec.times <- numeric(length(ID))
truesource <- numeric(length(ID))
for (i in 1:length(ID)) {
  sample.times[i] <- W$sampledata[which(W$sampledata[,1]==ID[i]),2]
  inf.times[i] <- W$epidata[which(W$epidata[,1]==ID[i]),3]
  rec.times[i] <- W$epidata[which(W$epidata[,1]==ID[i]),4]
  truesource[i] <- W$epidata[which(W$epidata[,1]==ID[i]),5]
}

# calculate likelihood of each potential transmission route using geometric-Poisson approximation of SNP distance
K <- transroutes(ID=ID, GD=GD, sample.times=sample.times, inf.times=inf.times,
                 rec.times=rec.times, mut.rate=mut.rate, eq.size=5000,
                 bottle.size=bn, p.level=0.95, summary=TRUE)

# create a matrix with most likely infector for each infectee
n <- length(ID)
Wmat_max <- matrix(0,n,n)
true <- matrix(0,n,n)
for (i in 1:n){
  Wmat_max[i,which(W$sampledata[,1]==K$maxpostsource[i])] <- 1
  true[i,which(W$sampledata[,1]==truesource[i])] <- 1
}

# normalize the matrix
for (i in 1:n){
  if (sum(Wmat_max[i,])>0) {
    Wmat_max[i,] <- Wmat_max[i,]/sum(Wmat_max[i,])
  }
}

# pick the ones that won't be sampled if imperfect sampling
missing <- sample(c(1:nrow(Wmat_max)),(1-sample_percent)*nrow(Wmat_max))
if (sample_percent <1){
  Wmat_max <- Wmat_max[-missing,-missing]
  true <- true[-missing,-missing]
}

# make matrix with posterior probabilites for each potential infector of each infectee
Wmat <- K$posterior

# normalize
n <- nrow(Wmat)
Wmat_norm <- matrix(0,n,n) # normalised score matrix
for (i in 1:n){
  if (sum(Wmat[i,])>0) {
    Wmat_norm[i,] <- Wmat[i,]/sum(Wmat[i,])
  }
}

# remove the ones that won't be sampled
if (sample_percent <1){
  Wmat_norm <- Wmat_norm[-missing,-missing]
}


# simulate pathogen evolution with deep sequencing
ID <- results_infected$InfectedNode
exp.times <- results_infected$DayExposed
inf.times <- results_infected$DayInfected 
rec.times <-results_infected$DayRecovered + 1
truesource <- results_infected$Infector
# by changing full to TRUE generate deep sequences
DS <- simfixoutbreak2(ID, exp.times, inf.times, rec.times,truesource, mut.rate=mut.rate,equi.pop=5000, shape=flat,samples.per.time = 1,
                      samp.schedule="random",samp.freq=10,inoc.size=bn,imp.var=25,
                      full=TRUE, feedback=1, glen=gen_len,
                      ref.strain=NULL)

observed.strains <- unique(unlist(DS$obs.strain)) # unique observed strains
polymorphic.loci <- NULL
for (i in 1:length(observed.strains)) {
  polymorphic.loci <- c(polymorphic.loci, DS$libr[[observed.strains[i]]])
}
# polymorphisms observed over all samples
polymorphic.loci <- unique(polymorphic.loci) 

n <- nrow(DS$sampledata)
host.polys <- list() # polymorphisms by host
host.poly.freq <- list() # frequencies for each polymorphism
for (i in 1:n) {
  hoststrains <- unlist(DS$obs.strain[[i]])
  hostfreq <- unlist(DS$obs.freq[[i]])
  hp <- NULL
  hf <- NULL
  for (q in 1:length(hoststrains)) {
    hp <- c(hp, DS$libr[[hoststrains[q]]]) # all locations of mutations
    hf <- c(hf, rep(hostfreq[q],length(DS$libr[[hoststrains[q]]]))) # frequency
  }
  if (sum(is.na(hp)>0)) {
    hp <- hp[-which(is.na(hp))] #takes out NAs
  }
  host.polys[[i]] <- unique(hp)
  hf1 <- NULL
  for (q in 1:length(host.polys[[i]])) {
    hf1 <- c(hf1, sum(hf[which(hp==host.polys[[i]][q])])) ## (# strains in person i with that polymorphism * frequency of strain)
  }
  host.poly.freq[[i]] <- hf1/sum(DS$obs.freq[[i]])
}

polytable <- matrix(0, length(polymorphic.loci), n) ## make matrix with each polymorphism and frequency in each person

for (i in 1:n) {
  for (q in 1:length(host.polys[[i]])) {
    polytable[which(polymorphic.loci==host.polys[[i]][q]),i] <- host.poly.freq[[i]][q]
  }
}

sharers <- apply(polytable,1,function(x){sum(!x%in%c(0,1))}) ## ignore polymorphisms that are in every strain in someone or not in anyone
sharedpolylocs <- which(sharers>1)
sharedpolys <- polymorphic.loci[sharedpolylocs]

if (sum(sharers==0)>0) {
  sharers <- sharers[-which(sharers==0)]
}
sharedpolylocs <- which(sharers>1)
sharedpolys <- polymorphic.loci[sharedpolylocs]

cand <- list()


#DS_data <- merge(DS$sampledata,DS$epidata,by.x="pID",by.y="ID")

scoremat <- matrix(0,n,n) # number of shared variants between hosts
for (i in 1:n) {
  if (sum(polymorphic.loci[which(!polytable[,i]%in%c(0,1))]%in%sharedpolys)>1) { 
    # shared polys for this host
    B <- polymorphic.loci[which(!polytable[,i]%in%c(0,1))]
    B <- which(polymorphic.loci%in%B[which(B%in%sharedpolys)]) # for each person, which polymorphisms are shared across
    candidates <- NULL
    for (q in 1:length(B)) {
      candidates <- c(candidates, which(!polytable[B[q],]%in%c(0,1))) # which people have shared polymorphism
    }
    candidates <- unique(candidates)
    candidates <- candidates[-which(candidates==i)]
    
    # only have candidates for people infected after
    candidates_ID <- DS$sampledata[candidates]
    cutoff_time_high <- DS$epidata[which(DS$epidata[,1] %in% DS$sampledata[i]),3]
    possible_infectors <- DS$epidata[DS$epidata[,3]<cutoff_time_high,1]
    candidates_time <- candidates_ID[which(candidates_ID %in% possible_infectors)]
    candidates2 <- c()
    for (r in 1:length(candidates_time)){
      candidates2 <- c(candidates2,which(DS$sampledata[,1]==candidates_time[r]))
    }
    
    if (length(candidates2)>0) {
      score <- numeric(length(candidates2))
      for (q in 1:length(score)) {
        score[q] <- length(intersect(which(!polytable[,candidates2[q]]%in%c(0,1)),
                                     which(!polytable[,i]%in%c(0,1))))
      }
      cand[[i]] <- candidates2
      scoremat[i,candidates2] <- score
    }
  }
}

# which are the (joint) likeliest sources? (share most vars)

transmat <- matrix(0,n,n) 
transmat2 <- matrix(0,n,n)
scoremattrans <- matrix(0,n,n) # normalised score matrix
truemat <- matrix(0,n,n) # true transmission matrix
parsimony <- apply(scoremat,1,function(x){which(x==max(x))}) #one shares most with
for (i in 1:n) {
  if (!i%in%parsimony[[i]]) {
    transmat[i,parsimony[[i]]] <- 1 #transmat is max: which one(s) share most with (scoremat includes more with amount shared, scoremattrans is standardized)
  }
  if (sum(transmat[i,])>0) {
    transmat2[i,] <- transmat[i,]/sum(transmat[i,])
  }
  truemat[i, which(DS$sampledata[,1]==DS$epidata[
    which(DS$epidata[,1]==DS$sampledata[i,1]),5])] <- 1
  if (sum(scoremat[i,])>0) {
    scoremattrans[i,] <- scoremat[i,]/sum(scoremat[i,])
  }
}

scorematflat <- scoremat
scorematflat[which(scorematflat>1)] <- 1


# Create matrices for consensus sequence purely distance
ROC_close <- list()
ROC_weight <- list()
for (j in 1:1) { #50
  singlesamp <- numeric(n)
  for (i in 1:n) {
    cat(j)
    if (length(DS$obs.strain[[i]])>1) {
      singlesamp[i] <- sample(DS$obs.strain[[i]], 1, prob=DS$obs.freq[[i]])
    } else {
      singlesamp[i] <- DS$obs.strain[[i]]
    }
  }
  
  # use consensus sequence matrix from above
  distmat <- GD 
  
  distmat2 <- matrix(1,nrow=n,ncol=n)
  # only allow candidates for people infected after
  for (row in 1:n){
    cutoff_time_high <- W$epidata[which(W$epidata[,1] %in% W$sampledata[row]),3]
    possible_infectors <- W$epidata[W$epidata[,3]<cutoff_time_high,1]
    distmat2[row,which(W$sampledata[,1] %in% possible_infectors)] <- distmat[row,which(W$sampledata[,1] %in% possible_infectors)]
    not_infector <- setdiff(W$sampledata[,1],possible_infectors)
    distmat2[row,which(W$sampledata[,1] %in% not_infector)] <- 0
  }
  
  genclose <- matrix(0,n,n)
  genclosenorm <- matrix(0,n,n)
  gen_weight <- matrix(10,n,n)
  gen_weight[which(distmat2>0)] <- 1/distmat2[which(distmat2>0)]
  diag(gen_weight) <- 0
  gen_weight[which(gen_weight==10)] <- 0
  
  # make matrix with potential infector with minimum genetic distance (closest)
  for (i in 1:n) {
    genclose[i,which(distmat2[i,]==min(distmat2[i,-i]))] <- 1
    genclosenorm[i,] <- genclose[i,]/sum(genclose[i,])
    # normalize weighted
    gen_weight[i,] <- gen_weight[i,]/sum(gen_weight[i,])
    
  }
  # normalize
  for (i in 1:n) {
    genclose[i,which(distmat2[i,]==min(distmat2[i,-i]))] <- 1
    genclosenorm[i,] <- genclose[i,]/sum(genclose[i,])
  }
  diag(genclose) <- 0
  gen_weight[which(is.nan(gen_weight))] <- 0
  
  # remove ones that aren't sampled
  if (sample_percent <1){
    genclosenorm <- genclosenorm[-missing,-missing]
    gen_weight <- gen_weight[-missing,-missing]
  }
  
  ROC_close[[j]] <- ROC_calc(genclosenorm, true, problevels=c(0,1/(20:1)))
  ROC_weight[[j]] <- ROC_calc(gen_weight, true, problevels=c(0,1/(100:1)))
  #lines(ROC_close[[j]]$FPR, ROC_close[[j]]$TPR, col=rgb(0,0.5,0,0.1))
}



# create hybrid -- if scoremat doesn't provide any info, use geometric Poisson SNP distance weights
hybrid_weight <- scoremattrans
for (i in 1:nrow(hybrid_weight)){
  if (sum(hybrid_weight[i,])==0){
    Wmat_row <- which(W$sampledata[,1]==DS$sampledata[i,1])
    for (j in 1:ncol(hybrid_weight)){
      Wmat_col <- which(W$sampledata[,1]==DS$sampledata[j,1])
      hybrid_weight[i,j] <- Wmat_norm[Wmat_row,Wmat_col]
      #cat(i,Wmat_row,j,Wmat_col,"\n")
    }
  }
}

# create hybrid -- if transmat doesn't provide any info, use geometirc Poisson SNP distance 
hybrid_max <- transmat2
for (i in 1:nrow(hybrid_max)){
  if (sum(hybrid_max[i,])==0){
    Wmat_row <- which(W$sampledata[,1]==DS$sampledata[i,1])
    for (j in 1:ncol(hybrid_max)){
      Wmat_col <- which(W$sampledata[,1]==DS$sampledata[j,1])
      hybrid_max[i,j] <- Wmat_max[Wmat_row,Wmat_col]
      #cat(i,Wmat_row,j,Wmat_col,"\n")
    }
  }
}

# recreate edge list from above and calculate RR for deep sequence data
for (df in 1:4){
  print(df)
  if (df==1){
    data <- scoremattrans
  } else if (df==2){
    data <- hybrid_weight
  } else if (df==3){
    data <- gen_weight
  } else if (df==4){
    data <- Wmat_norm
  } 
  results_edge_list_analysis_DS <- data.frame("InfectedNode"=NA,"node_stat"=NA,"Infector"=NA,"infector_stat"=NA,"prob"=NA,"max"=NA,"contacts"=NA,"SIprob"=NA)
  if (sum(data)>0){
    for (i in 1:nrow(data)){
      if (df==3 | df==4){
        cat(i,"\n")
        closest_sources <- which(data[i,]!=0)
        for (inf_source in 1:length(closest_sources)){
          if (sum(closest_sources)==0){
            results_edge_list_analysis_DS <- 
              rbind(results_edge_list_analysis_DS,c(W$sampledata[i],V(g)$trialstatus[W$sampledata[i]],
                                                    0,NA,
                                                    0,0,0,NA))
          } else{
            closest_sources_table <- as.data.frame(cbind(closest_sources,c(data[i,closest_sources])))
            names(closest_sources_table) <- c("closest_sources","prob")
            max_sources <- which(data[i,] == max(data[i,]))
            max <- ifelse(closest_sources[inf_source] %in% max_sources,1/length(max_sources),0)
            contact <- ifelse(W$sampledata[closest_sources[inf_source]] %in% neighbors(g,W$sampledata[i,1]),1,0)
            SI_prob <- distribution_norm[as.numeric(W$epidata[which(W$epidata[,1]==W$sampledata[i]),3] - W$epidata[which(W$epidata[,1]==W$sampledata[closest_sources[inf_source]]),3])]
            results_edge_list_analysis_DS <- 
              rbind(results_edge_list_analysis_DS,c(W$sampledata[i],V(g)$trialstatus[W$sampledata[i]],
                                                    W$sampledata[closest_sources[inf_source]],V(g)$trialstatus[W$sampledata[closest_sources[inf_source]]],
                                                    data[i,closest_sources[inf_source]],max,contact,SI_prob))
          }
        }
      } else {
        # track all sources
        cat(i,"\n")
        closest_sources <- which(data[i,]!=0)
        for (inf_source in 1:length(closest_sources)){
          if (sum(closest_sources)==0){
            results_edge_list_analysis_DS <- 
              rbind(results_edge_list_analysis_DS,c(DS$sampledata[i],V(g)$trialstatus[DS$sampledata[i]],
                                                    0,NA,
                                                    0,0,0,NA))
          } else{
            closest_sources_table <- as.data.frame(cbind(closest_sources,c(data[i,closest_sources])))
            names(closest_sources_table) <- c("closest_sources","prob")
            max_sources <- which(data[i,] == max(data[i,]))
            max <- ifelse(closest_sources[inf_source] %in% max_sources,1/length(max_sources),0)
            contact <- ifelse(DS$sampledata[closest_sources[inf_source]] %in% neighbors(g,DS$sampledata[i,1]),1,0)
            SI_prob <- distribution_norm[as.numeric(DS$epidata[which(DS$epidata[,1]==DS$sampledata[i]),3] - DS$epidata[which(DS$epidata[,1]==DS$sampledata[closest_sources[inf_source]]),3])]
            results_edge_list_analysis_DS <- 
              rbind(results_edge_list_analysis_DS,c(DS$sampledata[i],V(g)$trialstatus[DS$sampledata[i]],
                                                    DS$sampledata[closest_sources[inf_source]],V(g)$trialstatus[DS$sampledata[closest_sources[inf_source]]],
                                                    data[i,closest_sources[inf_source]],max,contact,SI_prob))
          }
        }
      }
    }
    results_edge_list_analysis_DS <- unique(results_edge_list_analysis_DS)
    
    # multiply probability by SI
    results_edge_list_analysis_DS$prob_SI <- results_edge_list_analysis_DS$prob*results_edge_list_analysis_DS$SIprob
    
    # normalize prob SI
    prob_SI_sums <- aggregate.data.frame(results_edge_list_analysis_DS$prob_SI,by=list(results_edge_list_analysis_DS$InfectedNode),sum,na.rm=TRUE)
    names(prob_SI_sums) <- c("InfectedNode","prob_SI_sums")
    results_edge_list_analysis_DS <- merge(results_edge_list_analysis_DS,prob_SI_sums,by="InfectedNode")
    results_edge_list_analysis_DS$prob_SI_norm <- results_edge_list_analysis_DS$prob_SI/results_edge_list_analysis_DS$prob_SI_sums
    
    SI_norm_max <- aggregate(results_edge_list_analysis_DS$prob_SI_norm,by=list(results_edge_list_analysis_DS$InfectedNode),FUN=function(x) max(x,na.rm=TRUE))
    names(SI_norm_max) <- c("InfectedNode","max_prob_SI")
    results_edge_list_analysis_DS <- merge(results_edge_list_analysis_DS, SI_norm_max, by = "InfectedNode", all.y = T)
    results_edge_list_analysis_DS$max_prob_SI[is.na(results_edge_list_analysis_DS$prob_SI_norm)] <- NA
    results_edge_list_analysis_DS$max_prob_SI<- ifelse(results_edge_list_analysis_DS$max_prob_SI == results_edge_list_analysis_DS$prob_SI_norm, 1, 0)
    SI_max <- aggregate(results_edge_list_analysis_DS$max_prob_SI,by=list(results_edge_list_analysis_DS$InfectedNode),sum,na.rm=TRUE)
    names(SI_max) <- c("InfectedNode","num_max")
    results_edge_list_analysis_DS <- merge(results_edge_list_analysis_DS, SI_max, by = "InfectedNode", all.y = T)
    results_edge_list_analysis_DS$max_prob_SI_norm <- results_edge_list_analysis_DS$max_prob_SI/results_edge_list_analysis_DS$num_max
    
    
    if (df==1){
      results_edge_list_analysis_DS_scoremat <- results_edge_list_analysis_DS
    } else if (df==2){
      results_edge_list_analysis_DS_hybrid <- results_edge_list_analysis_DS
    } else if (df==3){
      results_edge_list_analysis_DS_dist <- results_edge_list_analysis_DS
    } else if (df==4){
      results_edge_list_analysis_W <- results_edge_list_analysis_DS
    } 
  }
}

# # VEs
numevents_vacc <- nrow(results_edge_list_analysis[(results_edge_list_analysis$eventstatus==1) & (results_edge_list_analysis$TrialStatus==1),])
numevents_cont <- nrow(results_edge_list_analysis[(results_edge_list_analysis$eventstatus==1) & (results_edge_list_analysis$TrialStatus==0),])
num_vacc <- nrow(results_edge_list_analysis[(results_edge_list_analysis$TrialStatus==1),])
num_cont <- nrow(results_edge_list_analysis[(results_edge_list_analysis$TrialStatus==0),])
VEs <- VE[1]

args=(commandArgs(TRUE))
for (i in 1:length(args)) {
  eval (parse (text = args[[i]] ))
}
j<-j

genome_results <- as.data.frame(cbind(j,infect_VE,direct_VE,R0,mut.rate,bn,numevents_cont,numevents_vacc,VE))

# Write files ----
write.csv(genome_results,paste0(j,'_',R0,'_',mut.rate,'_genome_results_new.csv'))
write.csv(results_contacts, paste0(j,"_results_contacts",R0,"_",infect_VE,"_",direct_VE,"_",mut.rate,"_",bn,"_",infperiod_shape,"_",sample_percent,"_","sample_new.csv"))
write.csv(results_edge_list_analysis_DS_scoremat, paste0(j,"_results_edge_list_analysis_DS_scoremat",R0,"_",infect_VE,"_",direct_VE,"_",mut.rate,"_",bn,"_",infperiod_shape,"_",sample_percent,"_","sample_new.csv"))
write.csv(results_edge_list_analysis_DS_hybrid, paste0(j,"_results_edge_list_analysis_DS_hybrid",R0,"_",infect_VE,"_",direct_VE,"_",mut.rate,"_",bn,"_",infperiod_shape,"_",sample_percent,"_","sample_new.csv"))
write.csv(results_edge_list_analysis_DS_dist, paste0(j,"_results_edge_list_analysis_DS_dist",R0,"_",infect_VE,"_",direct_VE,"_",mut.rate,"_",bn,"_",infperiod_shape,"_",sample_percent,"_","sample_new.csv"))
write.csv(results_edge_list_analysis_W, paste0(j,"_results_edge_list_analysis_DS_W",R0,"_",infect_VE,"_",direct_VE,"_",mut.rate,"_",bn,"_",infperiod_shape,"_",sample_percent,"_","sample_new.csv"))




  
