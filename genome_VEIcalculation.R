# simulate pathogen evolution using seedy (written by Colin Worby); adapted by Rebecca Kahn
# estimate VEI 

# ROC_calc function ----
ROC_calc <- function(mat, truemat, problevels=NULL) {
  row_sums <- apply(mat,1,sum)
  if(sum(row_sums) >0){
    mat <- mat[which(row_sums!=0),]
    truemat2 <- truemat[which(row_sums!=0),]
  } else {
    truemat2 <- truemat 
  }
  n <- nrow(mat)
  
  if (is.null(problevels)) {
    problevels <- sort(unique(c(0,1,mat)))
  }
  
  TPR <- numeric(length(problevels))
  FPR <- numeric(length(problevels))
  
  for (i in 1:length(problevels)) {
    tot_s <- sum(mat>problevels[i])
    for (j in 1:n) {
      if (sum(truemat2[j,])>0) {
        if (mat[j,which(truemat2[j,]==1)]>problevels[i]) {
          TPR[i] <- TPR[i]+1
        }
      }
    }
    FPR[i] <- tot_s-TPR[i]
  }
  FPR <- rev(FPR)/(n*(n-1))
  TPR <- rev(TPR)/n
  
  FPR <- c(FPR, 1)
  TPR <- c(TPR, TPR[length(TPR)])
  
  AUC <- sum((FPR[-1]-FPR[-length(FPR)])*(TPR[-1]+TPR[-length(FPR)])/2)
  
  return(invisible(list(TPR=TPR, FPR=FPR, AUC=AUC, problevels=problevels)))
}

# Parameters  -------

betas <- c(0.0026,0.0029,0.003,0.0032,0.0035,0.0045,0.006)
direct_VEs <- c(0.6,0.8)
infect_VEs <- c(0,0.3,0.5,0.7)

m<-2
n<-1
p<-2

j<-1
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
rate_between<-.001

# Disease characteristics:
# Per-time-step hazard of infection for a susceptible nodes from an infectious neighbour
beta<-betas[m]
# Expected number of importations to the population over two years
num_introductions<-20
# Leaky multiplicative efficacy of vaccine
direct_VE<-direct_VEs[n]
infect_VE <- infect_VEs[p]
# Gamma-distribution parameters of infectious period
infperiod_shape<-500 #average 5 days, 4 is recovery average 6.67 days
infperiod_rate<-100
# First day of trial enrollment, relative to start of epidemic
trial_startday<-1
# Days of follow-up
trial_length<-300
# Number of clusters targeted for enrollment
# Must be less than or equal to the number of communities
num_clusters_enrolled_per_day<-num_communities
if (num_clusters_enrolled_per_day > num_communities) {stop("Enrolling too many communities!")}
# Number of days over which subjects are enrolled
enrollment_period<-1
# Target community enrollment proportion
cluster_coverage<-1 # not enrolling whole population so for cluster randomization will get overall (not total)
# if imperfect sampling, percent that are sampled:
sample_percent <- 1

# mutation rate source: http://science.sciencemag.org/content/348/6230/117.long - correct interpretation?

gen_len <- 13155  #Campbell 2018
mut.rate <- 0.0000119*gen_len/4 #Campbell 2018 gives 0.0000119 per site(?) per day
#but need per generation time (same article says 3 days but leaving as 1 for now)
bn <- 10

# Calculate R0
R0 <- (1 - (infperiod_rate/(infperiod_rate+beta))^infperiod_shape) *
  (((ave_community_size-1)*(1-rate_within)*rate_within + 
      (num_communities-1)*ave_community_size*(1-rate_between)*rate_between + 
      ((ave_community_size-1)*rate_within + (num_communities-1)*ave_community_size*rate_between)^2)/
     ((ave_community_size-1)*rate_within + (num_communities-1)*ave_community_size*rate_between)-1)
R0

args=(commandArgs(TRUE))
for (i in 1:length(args)) {
  eval (parse (text = args[[i]] ))
}
j<-j

# Load files -----

require(seedy)
require(igraph)

# read in results list and network
filename <- paste0('/n/home00/rkahn/seq/',j,'_results_edge_list_analysis1.51_0.3_0.6_0.039136125_10_500_1_sample.csv')
results_edge_list_analysis <- read.csv(filename)
results_edge_list_analysis$Infector[is.na(results_edge_list_analysis$Infector) & results_edge_list_analysis$eventstatus==1] <- 0
results_infected <- results_edge_list_analysis[results_edge_list_analysis$eventstatus==1,]
# for those who remain infected after end of the trial, make recovered day after trial ends
results_infected$DayRecovered[is.na(results_infected$DayRecovered)] <- trial_startday + trial_length

g <- paste0('/n/home00/rkahn/seq/',j,'_g1.51_0.3_0.6_0.039136125_10_500_1_sample.Rdata')
load(g)

# Genome ------
# simulate pathogen evolution
ID <- results_infected$InfectedNode
ID
inf.times <- results_infected$DayInfected # assuming no incubation period and day infected = day infectious
rec.times <-results_infected$DayRecovered
truesource <- results_infected$Infector

# consensus sequences
W <- simfixoutbreak(ID, inf.times, rec.times,truesource, mut.rate=mut.rate, shape=flat,
                    samp.schedule="random",equi.pop=5000,samp.freq=10,inoc.size=bn,imp.var=25,
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
  inf.times[i] <- W$epidata[which(W$epidata[,1]==ID[i]),2]
  rec.times[i] <- W$epidata[which(W$epidata[,1]==ID[i]),3]
  truesource[i] <- W$epidata[which(W$epidata[,1]==ID[i]),4]
}

# calculate likelihood of each potential transmission route using geometric-Poisson approximation of SNP distance
K <- transroutes(ID=ID, GD=GD, sample.times=sample.times, inf.times=inf.times,
                 rec.times=rec.times, mut.rate=mut.rate, eq.size=5000,
                 bottle.size=bn, p.level=0.95, summary=TRUE)

# track true source and most likely source
sources <- data.frame(cbind(W$sampledata[,1],truesource,K$maxpostsource))

# make NA any source who was externally infected
replace <- which(sources$V3==0)
sources$V3[replace] <- NA

# create a matrix with most likely infector for each infectee
n <- length(ID)
Wmat_max <- matrix(0,n,n)
true <- matrix(0,n,n)
for (i in 1:n){
  source <- sources[i,3]
  true_source <- sources[i,2]
  Wmat_max[i,which(W$sampledata[,1]==K$maxpostsource[i])] <- 1
  true[i,which(W$sampledata[,1]==true_source)] <- 1
}

# normalize the matrix
for (i in 1:n){
  if (sum(Wmat_max[i,])>0) {
    Wmat_max[i,] <- Wmat_max[i,]/sum(Wmat_max[i,])
  }
}

# set a probability cutoff for including potential infectors
#Wmat_max[Wmat_max<0.25] <- 0

# pick the ones that won't be sampled if imperfect sampling
missing <- sample(c(1:nrow(Wmat_max)),(1-sample_percent)*nrow(Wmat_max))
if (sample_percent <1){
  Wmat_max <- Wmat_max[-missing,-missing]
  true <- true[-missing,-missing]
}

ROC_W_max <-ROC_calc(Wmat_max,true)
pairwise_AUC_max <- ROC_W_max$AUC

# make matrix with posterior probabilites for each potential infector of each infectee
Wmat <- K$posterior

# set a cutoff for probability
#Wmat[Wmat<0.25] <- 0

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

# calculate AUC
ROC_W <- ROC_calc(Wmat_norm,true)
pairwise_AUC <- ROC_W$AUC

# recreate edge list results table from above but with all potential infectors for each infectee
data <- Wmat_norm
results_edge_list_analysis_W <- data.frame("InfectedNode"=NA,"node_stat"=NA,"Infector"=NA,"infector_stat"=NA,"prob"=NA,"max"=NA,"contacts"=NA)
for (i in 1:nrow(Wmat_norm)){
  # track all sources
  cat(i,"\n")
  closest_sources <- which(data[i,]!=0)
  for (inf_source in 1:length(closest_sources)){
    if (sum(closest_sources)==0){
      results_edge_list_analysis_W <- 
        rbind(results_edge_list_analysis_W,c(W$sampledata[i],V(g)$trialstatus[W$sampledata[i]],
                                             0,NA,
                                             0,0,0))
    } else{
      closest_sources_table <- as.data.frame(cbind(closest_sources,c(data[i,closest_sources])))
      names(closest_sources_table) <- c("closest_sources","prob")
      # top2_source_probs <- closest_sources_table[order(closest_sources_table$prob,decreasing=TRUE),][1:2,]
      # top2_source <- top2_source_probs[!is.na(top2_source_probs[,2]),]
      # top2_source <- top2_source[,1]
      # norm_prob <- sum(top2_source_probs[,2],na.rm=TRUE)
      # top2 <- ifelse(closest_sources[inf_source] %in% top2_source,top2_source_probs[which(top2_source_probs[,1]==closest_sources[inf_source]),2],0)
      max_sources <- which(data[i,] == max(data[i,]))
      max <- ifelse(closest_sources[inf_source] %in% max_sources,1/length(max_sources),0)
      contact <- ifelse(W$sampledata[closest_sources[inf_source]] %in% neighbors(g,W$sampledata[i,1]),1,0)
      results_edge_list_analysis_W <- 
        rbind(results_edge_list_analysis_W,c(W$sampledata[i],V(g)$trialstatus[W$sampledata[i]],
                                             W$sampledata[closest_sources[inf_source]],V(g)$trialstatus[W$sampledata[closest_sources[inf_source]]],
                                             data[i,closest_sources[inf_source]],max,contact))
    }
  }
}
results_edge_list_analysis_W <- unique(results_edge_list_analysis_W)

# get average number of infectors
results_edge_list_analysis_W_prob <- results_edge_list_analysis_W[results_edge_list_analysis_W$prob!=0,]
av_num_infectors_table <- aggregate(results_edge_list_analysis_W_prob$prob,by=list(results_edge_list_analysis_W_prob$InfectedNode),length)
av_num_infectors_W <- mean(av_num_infectors_table$x,na.rm=TRUE)

# note the true infector
results_edge_list_analysis_true <- results_edge_list_analysis[1:nrow(results_infected),c(2,10)]
names(results_edge_list_analysis_true) <- c("InfectedNode","TrueInfector")
results_edge_list_analysis_W_true <- merge(results_edge_list_analysis_W_prob,results_edge_list_analysis_true,by="InfectedNode")
results_edge_list_analysis_W_true$truth <- ifelse(results_edge_list_analysis_W_true$Infector==results_edge_list_analysis_W_true$TrueInfector,1,0)
mean_prob_truth_W <- mean(results_edge_list_analysis_W_true$prob[results_edge_list_analysis_W_true$truth==1],na.rm=TRUE)
sd_prob_truth_W <- sd(results_edge_list_analysis_W_true$prob[results_edge_list_analysis_W_true$truth==1],na.rm=TRUE)

# note true infector stat
results_edge_list_analysis_W_true$true_stat <- ifelse(results_edge_list_analysis_W_true$TrueInfector==0,NA,V(g)$trialstatus[results_edge_list_analysis_W_true$TrueInfector])

# note if potential infector has same trial status as truth
results_edge_list_analysis_W_true$correctstat <- ifelse(results_edge_list_analysis_W_true$infector_stat==results_edge_list_analysis_W_true$true_stat,1,0)

# sum of correct probabilities
results_edge_list_analysis_W_true$correctstatprob <- results_edge_list_analysis_W_true$prob*results_edge_list_analysis_W_true$correctstat
correct_stat_probs <- aggregate(results_edge_list_analysis_W_true$correctstatprob,by=list(results_edge_list_analysis_W_true$InfectedNode),sum)
names(correct_stat_probs) <- c("InfectedNode","CorrectProb")
mean_prob_correct_stat_W <- mean(correct_stat_probs$CorrectProb,na.rm=TRUE)
sd_prob_correct_stat_W <- sd(correct_stat_probs$CorrectProb,na.rm=TRUE)

# Average number of maximum likely infectors with the correct trial status
results_edge_list_analysis_W_true$correctstatprobmax <- results_edge_list_analysis_W_true$max*results_edge_list_analysis_W_true$correctstat
correct_stat_probs_max <- aggregate(results_edge_list_analysis_W_true$correctstatprobmax,by=list(results_edge_list_analysis_W_true$InfectedNode),sum)
names(correct_stat_probs_max) <- c("InfectedNode","CorrectProb")
mean_prob_correct_stat_max_W <- mean(correct_stat_probs_max$CorrectProb,na.rm=TRUE)
sd_prob_correct_stat_max_W <- sd(correct_stat_probs_max$CorrectProb,na.rm=TRUE)

#1# summarize by trial status of infector
summary_infectors_W_stat <- aggregate(results_edge_list_analysis_W$InfectedNode,by=c(list(results_edge_list_analysis_W$InfectedNode),list(results_edge_list_analysis_W$infector_stat)),length)
names(summary_infectors_W_stat)<-c("InfectedNode","infector_stat","number")
summary_infectors_W <- aggregate(summary_infectors_W_stat$number,by=list(summary_infectors_W_stat$InfectedNode),sum)
names(summary_infectors_W)<-c("InfectedNode","total")
summary_infectors_W_merged <- merge(summary_infectors_W_stat,summary_infectors_W)

summary_infectors_W_merged$perc <- summary_infectors_W_merged$number/summary_infectors_W_merged$total
summary_infectors_W_final <- aggregate(summary_infectors_W_merged$perc,by=list(summary_infectors_W_merged$infector_stat),sum)
names(summary_infectors_W_final)<-c("infector_stat","number")
RR_inf_W <- summary_infectors_W_final$number[summary_infectors_W_final$infector_stat==1]/summary_infectors_W_final$number[summary_infectors_W_final$infector_stat==0]

#2# look only at most likely
summary_infectors_W_max <- aggregate(results_edge_list_analysis_W$max,by=list(results_edge_list_analysis_W$infector_stat),sum)
names(summary_infectors_W_max)<-c("infector_stat","number")
vacc_inf_W_max <- summary_infectors_W_max$number[summary_infectors_W_max$infector_stat==1]
cont_inf_W_max<- summary_infectors_W_max$number[summary_infectors_W_max$infector_stat==0]
RR_inf_W_max <- vacc_inf_W_max/cont_inf_W_max

#2c# use contacts to break ties
results_edge_list_analysis_W_max_contacts <- results_edge_list_analysis_W[results_edge_list_analysis_W$max!=0 & !is.na(results_edge_list_analysis_W$max),]
max_ties <- unique(results_edge_list_analysis_W_max_contacts$InfectedNode[(results_edge_list_analysis_W_max_contacts$max!=1) 
                                                                          & !is.na(results_edge_list_analysis_W_max_contacts$max)])
for (i in 1:length(max_ties)){
  max_contacts <- results_edge_list_analysis_W_max_contacts$Infector[results_edge_list_analysis_W_max_contacts$InfectedNode==max_ties[i] & results_edge_list_analysis_W_max_contacts$contacts==1]
  not_contacts <- results_edge_list_analysis_W_max_contacts$Infector[results_edge_list_analysis_W_max_contacts$InfectedNode==max_ties[i] & results_edge_list_analysis_W_max_contacts$contacts==0]
  if (length(max_contacts) > 0){
    results_edge_list_analysis_W_max_contacts$max[results_edge_list_analysis_W_max_contacts$InfectedNode==max_ties[i] & results_edge_list_analysis_W_max_contacts$Infector %in% max_contacts] <- 1/length(max_contacts)
    results_edge_list_analysis_W_max_contacts$max[results_edge_list_analysis_W_max_contacts$InfectedNode==max_ties[i] & results_edge_list_analysis_W_max_contacts$Infector %in% not_contacts] <- 0
  }
}
summary_infectors_W_max_contacts <- aggregate(results_edge_list_analysis_W_max_contacts$max,by=list(results_edge_list_analysis_W_max_contacts$infector_stat),sum)
names(summary_infectors_W_max_contacts)<-c("infector_stat","number")
vacc_inf_W_max_contacts <- summary_infectors_W_max_contacts$number[summary_infectors_W_max_contacts$infector_stat==1]
cont_inf_W_max_contacts<- summary_infectors_W_max_contacts$number[summary_infectors_W_max_contacts$infector_stat==0]
RR_inf_W_max_contacts <- vacc_inf_W_max_contacts/cont_inf_W_max_contacts

#3# aggregate by probability instead of number of infectors
summary_infectors_W_prob <- aggregate(results_edge_list_analysis_W$prob,by=list(results_edge_list_analysis_W$infector_stat),sum)
names(summary_infectors_W_prob)<-c("infector_stat","number")
vacc_inf_W_prob <- summary_infectors_W_prob$number[summary_infectors_W_prob$infector_stat==1]
cont_inf_W_prob<- summary_infectors_W_prob$number[summary_infectors_W_prob$infector_stat==0]
RR_inf_W_prob <- vacc_inf_W_prob/cont_inf_W_prob

#4# aggregate by probability instead of number of infectors (max)
results_edge_list_analysis_W_max <- results_edge_list_analysis_W[results_edge_list_analysis_W$max!=0 & !is.na(results_edge_list_analysis_W$max),]
summary_infectors_W_max_prob <- aggregate(results_edge_list_analysis_W_max$prob,by=list(results_edge_list_analysis_W_max$infector_stat),sum)
names(summary_infectors_W_max_prob)<-c("infector_stat","number")
vacc_inf_W_max_prob <- summary_infectors_W_max_prob$number[summary_infectors_W_max_prob$infector_stat==1]
cont_inf_W_max_prob<- summary_infectors_W_max_prob$number[summary_infectors_W_max_prob$infector_stat==0]
RR_inf_W_max_prob <- vacc_inf_W_max_prob/cont_inf_W_max_prob

#5# look at maximum of total probabilities for infector status
# aggregate by infected node
results_infector_stat_W <- aggregate(results_edge_list_analysis_W$prob,by=c(list(results_edge_list_analysis_W$InfectedNode),list(results_edge_list_analysis_W$infector_stat)),sum)
names(results_infector_stat_W) <- c("Infected Node","infector_stat","prob")
# note which trial status is the maximum for each infected node
results_infector_stat_W$max <- ifelse(results_infector_stat_W$prob>=0.5,1,0)
results_infector_stat_W_max <- results_infector_stat_W[results_infector_stat_W$max==1,]
# assign 1 to trial status with most probability
RR_inf_W_max_class <- nrow(results_infector_stat_W_max[results_infector_stat_W_max$infector_stat==1,]) / nrow(results_infector_stat_W_max[results_infector_stat_W_max$infector_stat==0,])

#6# assign total probability to trial status with most probability
RR_inf_W_prob_max_class <- sum(results_infector_stat_W_max$prob[results_infector_stat_W_max$infector_stat==1]) / sum(results_infector_stat_W_max$prob[results_infector_stat_W_max$infector_stat==0])

#7# look at maximum only among contacts
#restrict to only contacts
results_edge_list_analysis_W_contacts <- results_edge_list_analysis_W[results_edge_list_analysis_W$contacts!=0 & !is.na(results_edge_list_analysis_W$contacts),]
#look at infected nodes with identified contacts
ID_contacts <- unique(results_edge_list_analysis_W_contacts$InfectedNode)
# loop through each infected node
for (i in 1:length(ID_contacts)){
  # identify how many most likely infectors there are among the contacts
  num_max_contacts <- length(which(results_edge_list_analysis_W_contacts$prob[results_edge_list_analysis_W_contacts$InfectedNode==ID_contacts[i]]==max(results_edge_list_analysis_W_contacts$prob[results_edge_list_analysis_W_contacts$InfectedNode==ID_contacts[i]],na.rm=TRUE)))
  # make max = 1 for most likely infector (or 0.5 if multiple)
  results_edge_list_analysis_W_contacts$max[results_edge_list_analysis_W_contacts$InfectedNode==ID_contacts[i] & results_edge_list_analysis_W_contacts$prob == max(results_edge_list_analysis_W_contacts$prob[results_edge_list_analysis_W_contacts$InfectedNode==ID_contacts[i]],na.rm=TRUE)] <-1/num_max_contacts
}
summary_infectors_W_contacts <- aggregate(results_edge_list_analysis_W_contacts$max,by=list(results_edge_list_analysis_W_contacts$infector_stat),sum)
names(summary_infectors_W_contacts)<-c("infector_stat","number")
vacc_inf_W_contacts<- summary_infectors_W_contacts$number[summary_infectors_W_contacts$infector_stat==1]
cont_inf_W_contacts <- summary_infectors_W_contacts$number[summary_infectors_W_contacts$infector_stat==0]
RR_inf_W_contacts <- vacc_inf_W_contacts/cont_inf_W_contacts

#8# look only among contacts -- weighted
ID_contacts <- unique(results_edge_list_analysis_W_contacts$InfectedNode)
for (i in 1:length(ID_contacts)){
  results_edge_list_analysis_W_contacts$prob[results_edge_list_analysis_W_contacts$InfectedNode==ID_contacts[i]] <- results_edge_list_analysis_W_contacts$prob[results_edge_list_analysis_W_contacts$InfectedNode==ID_contacts[i]] / sum(results_edge_list_analysis_W_contacts$prob[results_edge_list_analysis_W_contacts$InfectedNode==ID_contacts[i]])
}
summary_infectors_W_contacts_prob <- aggregate(results_edge_list_analysis_W_contacts$prob,by=list(results_edge_list_analysis_W_contacts$infector_stat),sum)
names(summary_infectors_W_contacts_prob)<-c("infector_stat","number")
vacc_inf_W_contacts_prob<- summary_infectors_W_contacts_prob$number[summary_infectors_W_contacts_prob$infector_stat==1]
cont_inf_W_contacts_prob <- summary_infectors_W_contacts_prob$number[summary_infectors_W_contacts_prob$infector_stat==0]
RR_inf_W_contacts_prob <- vacc_inf_W_contacts_prob/cont_inf_W_contacts_prob

# simulate pathogen evolution with deep sequencing
ID <- results_infected$InfectedNode
inf.times <- results_infected$DayInfected # assume no incubation period
rec.times <-results_infected$DayRecovered
truesource <- results_infected$Infector
# by changing full to TRUE generate deep sequences?
DS <- simfixoutbreak(ID, inf.times, rec.times,truesource, mut.rate=mut.rate,  shape=flat,
                     samp.schedule="random",equi.pop=5000,samp.freq=10,inoc.size=bn,
                     full=TRUE, feedback=1, glen=gen_len,
                     ref.strain=NULL)

observed.strains <- unique(unlist(DS$obs.strain)) # unique observed strains
polymorphic.loci <- NULL
for (i in 1:length(observed.strains)) {
  polymorphic.loci <- c(polymorphic.loci, DS$libr[[observed.strains[i]]])
}
# polymorphisms observed over all samples
polymorphic.loci <- unique(polymorphic.loci) 

n <- nrow(DS$sampledata) #should be same as nrow(DS$epidata)
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

sharers <- apply(polytable,1,function(x){sum(!x%in%c(0,1))}) ## ignore polymorphisms that are in every strain in someone or not any anyone
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
      candidates <- c(candidates, which(!polytable[B[q],]%in%c(0,1))) # which people have share polymorphism
    }
    candidates <- unique(candidates)
    candidates <- candidates[-which(candidates==i)]
    
    # only have candidates for people infected after
    candidates_ID <- DS$sampledata[candidates]
    cutoff_time_high <- DS$epidata[which(DS$epidata[,1] %in% DS$sampledata[i]),2]
    #restrict to 2 serial intervals before
    cutoff_time_low <- cutoff_time_high - 2*(infperiod_shape/infperiod_rate)
    possible_infectors <- DS$epidata[DS$epidata[,2]>=cutoff_time_low & DS$epidata[,2]<cutoff_time_high,1]
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
    which(DS$epidata[,1]==DS$sampledata[i,1]),4])] <- 1
  if (sum(scoremat[i,])>0) {
    scoremattrans[i,] <- scoremat[i,]/sum(scoremat[i,])
  }
}

scorematflat <- scoremat
scorematflat[which(scorematflat>1)] <- 1

# remove ones that aren't sampled
if (sample_percent<1){
  scoremattrans <- scoremattrans[-missing,-missing]
  transmat2 <- transmat2[-missing,-missing]
  truemat <- truemat[-missing,-missing]
}

# look only at ones for which there is only 1 max
transmat2_1 <- transmat2
transmat2_1[transmat2_1<1] <-0

if (sum(scoremat)>0){
  ROC_scoremattrans <- ROC_calc(scoremattrans, truemat)
  ROC_transmat2 <- ROC_calc(transmat2, truemat)
  ROC_transmat2_1 <- ROC_calc(transmat2_1, truemat)
  scoremattrans_AUC <- ROC_scoremattrans$AUC
  transmat2_AUC <- ROC_transmat2$AUC
  transmat2_1_AUC <- ROC_transmat2_1$AUC
} else {
  scoremattrans_AUC <-NA
  transmat2_AUC <- NA
  transmat2_1_AUC <- NA
}

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
  
  distmat2 <- matrix(0,nrow=n,ncol=n)
  # only allow candidates for people infected after
  for (row in 1:n){
    cutoff_time_high <- DS$epidata[which(DS$epidata[,1] %in% DS$sampledata[row]),2]
    #restrict to 2 serial intervals before
    cutoff_time_low <- cutoff_time_high - 2*(infperiod_shape/infperiod_rate)
    possible_infectors <- DS$epidata[DS$epidata[,2]>=cutoff_time_low & DS$epidata[,2]<cutoff_time_high,1]
    distmat2[row,which(DS$sampledata[,1] %in% possible_infectors)] <- distmat[row,which(DS$sampledata[,1] %in% possible_infectors)]
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
  
  ROC_close[[j]] <- ROC_calc(genclosenorm, truemat, problevels=c(0,1/(20:1)))
  ROC_weight[[j]] <- ROC_calc(gen_weight, truemat, problevels=c(0,1/(100:1)))
  #lines(ROC_close[[j]]$FPR, ROC_close[[j]]$TPR, col=rgb(0,0.5,0,0.1))
}

ROCmeanTPR <- numeric(22)
ROCmeanFPR <- numeric(22)
ROCqt <- matrix(0,22,2)
for (i in 1:22) {
  ROCmeanTPR[i] <- mean(unlist(lapply(ROC_close, function(x){x$TPR[i]})))
  ROCmeanFPR[i] <- mean(unlist(lapply(ROC_close, function(x){x$FPR[i]})))
  ROCqt[i,] <- quantile(unlist(lapply(ROC_close, function(x){x$TPR[i]})), 
                        c(0.025,0.975))
}
ROCweightmeanTPR <- numeric(102)
ROCweightmeanFPR <- numeric(102)
ROCweightqt <- matrix(0,102,2)
for (i in 1:102) {
  ROCweightmeanTPR[i] <- mean(unlist(lapply(ROC_weight, function(x){x$TPR[i]})))
  ROCweightmeanFPR[i] <- mean(unlist(lapply(ROC_weight, function(x){x$FPR[i]})))
  ROCweightqt[i,] <- quantile(unlist(lapply(ROC_weight, function(x){x$TPR[i]})), 
                              c(0.025,0.975))
}


ROC_close_AUC <-numeric(j)
ROC_weight_AUC <- numeric(j)
for (i in 1:j){
  ROC_close_AUC[i] <- ROC_close[[i]]$AUC
  ROC_weight_AUC[i] <- ROC_weight[[i]]$AUC
}
close_AUC <- mean(ROC_close_AUC,na.rm=TRUE)
weight_AUC <- mean(ROC_weight_AUC,na.rm=TRUE)


# create hybrid -- if scoremat doesn't provide any info, use geometric Poisson SNP distance weights
hybrid_weight <- scoremattrans
for (i in 1:nrow(hybrid_weight)){
  if (sum(hybrid_weight[i,])==0){
    hybrid_weight[i,] <- Wmat_norm[i,]
  }
}

ROC_hw <-ROC_calc(hybrid_weight,truemat)
hybrid_weight_AUC <- ROC_hw$AUC

# create hybrid -- if transmat doesn't provide any info, use geometirc Poisson SNP distance 
hybrid_max <- transmat2
for (i in 1:nrow(hybrid_max)){
  if (sum(hybrid_max[i,])==0){
    hybrid_max[i,] <- Wmat_max[i,]
  }
}

ROC_hm <-ROC_calc(hybrid_max,truemat)
hybrid_max_AUC <- ROC_hm$AUC

# use shared variants to break geom Poisson SNP distance ties (only keep ones that are in the max shared variants matrix)
scoremat_max_tie <- Wmat_max
for (i in 1:nrow(scoremat_max_tie)){
  if (sum(scoremat_max_tie[i,]!=0)>1){
    infectors <- which(scoremat_max_tie[i,]!=0)
    max_SM_W <- which(transmat2[i,] == max(transmat2[i,]))
    scoremat_max_tie[i,max_SM_W] <- 1/length(max_SM_W)
    scoremat_max_tie[i,infectors[which(!(infectors %in% max_SM_W))]] <- 0
  }
  # if none identified by geom Poisson SNP, use shared variants
  if (sum(scoremat_max_tie[i,])==0){
    scoremat_max_tie[i,] <- transmat2[i,]
  }
}

# incorporate epi data into scoremat (restrict to contacts)
scoremattrans_epi <- scoremattrans
for (i in 1:nrow(scoremattrans_epi)){
  infectors <- which(scoremattrans_epi[i,]!=0)
  contacts <- which(DS$sampledata[,1] %in% neighbors(g,DS$sampledata[i,1]))
  overlap <- contacts[which(contacts %in% infectors)]
  # for those not contacts, make probability of being infector 0
  if (length(infectors)>0){
    scoremattrans_epi[i,infectors[which(!(infectors %in% contacts))]] <- 0
  }
  # renormalize
  for (i in 1:nrow(scoremattrans_epi)){
    if (sum(scoremattrans_epi[i,])>0) {
      scoremattrans_epi[i,] <- scoremattrans_epi[i,]/sum(scoremattrans_epi[i,])
    }
  }
}

ROC_scoremattrans_epi <- ROC_calc(scoremattrans_epi, truemat)
scoremattrans_epi_AUC_restrict <- ROC_scoremattrans_epi$AUC

# take only max of contacts (restrict to contacts)
scoremattrans_epi_max <- scoremattrans_epi
max_list <- apply(scoremattrans_epi_max,1,function(x){which(x==max(x))})
for (i in 1:nrow(scoremattrans_epi_max)){
  infectors <- which(scoremattrans_epi_max[i,]!=0)
  scoremattrans_epi_max[i,infectors[which(!(infectors == max_list[[i]]))]] <- 0
  infectors_max <- which(scoremattrans_epi_max[i,]!=0)
  scoremattrans_epi_max[i,infectors[which((infectors == max_list[[i]]))]] <- 1/length(infectors_max)
}

ROC_scoremattrans_epi_max <- ROC_calc(scoremattrans_epi_max, truemat)
scoremattrans_epi_AUC_max_restrict <- ROC_scoremattrans_epi_max$AUC

# take only max of contacts -- only 1 maximum (restrict to contacts)
scoremattrans_epi_max2 <- scoremattrans_epi_max
for (i in 1:nrow(scoremattrans_epi_max2)){
  infectors <- which(scoremattrans_epi_max2[i,]!=0)
  trial_status <- V(g)$trialstatus[DS$sampledata[infectors]]
  if (!(length(unique(trial_status) %in% (c(0,1))))){
    scoremattrans_epi_max2[i,] <- 0
  }
}

ROC_scoremattrans_epi_max2 <- ROC_calc(scoremattrans_epi_max2, truemat)
scoremattrans_epi_AUC_max_restrict2 <- ROC_scoremattrans_epi_max$AUC


# incorporate epi data into Wmat_norm (restrict)
Wmat_norm_epi <- Wmat_norm
for (i in 1:nrow(Wmat_norm_epi)){
  infectors <- which(Wmat_norm_epi[i,]!=0)
  contacts <- which(W$sampledata[,1] %in% neighbors(g,W$sampledata[i,1]))
  overlap <- contacts[which(contacts %in% infectors)]
  # for those not contacts, make probability of being infector 0
  if (length(infectors)>0){
    Wmat_norm_epi[i,infectors[which(!(infectors %in% contacts))]] <- 0
  }
  # renormalize
  for (i in 1:nrow(Wmat_norm_epi)){
    if (sum(Wmat_norm_epi[i,])>0) {
      Wmat_norm_epi[i,] <- Wmat_norm_epi[i,]/sum(Wmat_norm_epi[i,])
    }
  }
}

ROC_Wmat_norm_epi <- ROC_calc(Wmat_norm_epi, true)
Wmat_norm_epi_AUC_restrict <- ROC_Wmat_norm_epi$AUC

# take only max of contacts
Wmat_norm_epi_max <- Wmat_norm_epi
max_list <- apply(Wmat_norm_epi_max,1,function(x){which(x==max(x))})
for (i in 1:nrow(Wmat_norm_epi_max)){
  infectors <- which(Wmat_norm_epi_max[i,]!=0)
  Wmat_norm_epi_max[i,infectors[which(!(infectors == max_list[[i]]))]] <- 0
  infectors_max <- which(Wmat_norm_epi_max[i,]!=0)
  Wmat_norm_epi_max[i,infectors[which((infectors == max_list[[i]]))]] <- 1/length(infectors_max)
  
}

ROC_Wmat_norm_epi_max <- ROC_calc(Wmat_norm_epi_max, true)
Wmat_norm_epi_AUC_max_restrict <- ROC_Wmat_norm_epi_max$AUC

# incorporate epi data into hybrid (restrict to contacts)
hybrid_epi <- hybrid_weight
for (i in 1:nrow(hybrid_epi)){
  infectors <- which(hybrid_epi[i,]!=0)
  contacts <- which(DS$sampledata[,1] %in% neighbors(g,DS$sampledata[i,1]))
  overlap <- contacts[which(contacts %in% infectors)]
  # for those not contacts, make probability of being infector 0
  if (length(infectors)>0){
    hybrid_epi[i,infectors[which(!(infectors %in% contacts))]] <- 0
  }
  # renormalize
  for (i in 1:nrow(hybrid_epi)){
    if (sum(hybrid_epi[i,])>0) {
      hybrid_epi[i,] <- hybrid_epi[i,]/sum(hybrid_epi[i,])
    }
  }
}

ROC_hybrid_epi <- ROC_calc(hybrid_epi, truemat)
hybrid_epi_AUC_restrict <- ROC_hybrid_epi$AUC

# take only max of contacts (restrict to contacts)
hybrid_epi_max <- hybrid_epi
max_list <- apply(hybrid_epi_max,1,function(x){which(x==max(x))})
for (i in 1:nrow(hybrid_epi_max)){
  infectors <- which(hybrid_epi_max[i,]!=0)
  hybrid_epi_max[i,infectors[which(!(infectors == max_list[[i]]))]] <- 0
  infectors_max <- which(hybrid_epi_max[i,]!=0)
  hybrid_epi_max[i,infectors[which((infectors == max_list[[i]]))]] <- 1/length(infectors_max)
}

ROC_hybrid_epi_max <- ROC_calc(hybrid_epi_max, truemat)
hybrid_epi_AUC_max_restrict <- ROC_hybrid_epi_max$AUC

# incorporate epi data into dist (restrict to contacts)
dist_epi <- gen_weight
for (i in 1:nrow(dist_epi)){
  infectors <- which(dist_epi[i,]!=0)
  contacts <- which(DS$sampledata[,1] %in% neighbors(g,DS$sampledata[i,1]))
  overlap <- contacts[which(contacts %in% infectors)]
  # for those not contacts, make probability of being infector 0
  if (length(infectors)>0){
    dist_epi[i,infectors[which(!(infectors %in% contacts))]] <- 0
  }
  # renormalize
  for (i in 1:nrow(dist_epi)){
    if (sum(dist_epi[i,])>0) {
      dist_epi[i,] <- dist_epi[i,]/sum(dist_epi[i,])
    }
  }
}

ROC_dist_epi <- ROC_calc(dist_epi, truemat)
dist_epi_AUC_restrict <- ROC_dist_epi$AUC

# take only max of contacts (restrict to contacts)
dist_epi_max <- dist_epi
max_list <- apply(dist_epi_max,1,function(x){which(x==max(x))})
for (i in 1:nrow(dist_epi_max)){
  infectors <- which(dist_epi_max[i,]!=0)
  dist_epi_max[i,infectors[which(!(infectors == max_list[[i]]))]] <- 0
  infectors_max <- which(dist_epi_max[i,]!=0)
  dist_epi_max[i,infectors[which((infectors == max_list[[i]]))]] <- 1/length(infectors_max)
}

ROC_dist_epi_max <- ROC_calc(dist_epi_max, truemat)
dist_epi_AUC_max_restrict <- ROC_dist_epi_max$AUC

# incorporate epi data into transmat2 (tie breaker)
transmat2_epi <- transmat2
for (i in 1:nrow(transmat2_epi)){
  infectors <- which(transmat2_epi[i,]!=0)
  contacts <- which(DS$sampledata[,1] %in% neighbors(g,DS$sampledata[i,1]))
  overlap <- contacts[which(contacts %in% infectors)]
  if (length(overlap)>0){
    transmat2_epi[i,overlap] <- 1/length(overlap)
    transmat2_epi[i,infectors[which(!(infectors %in% contacts))]] <- 0
  }
}

ROC_transmat2_epi <- ROC_calc(transmat2_epi, truemat)
transmat2_epi_AUC_tiebreaker <- ROC_transmat2_epi$AUC

# incorporate epi data into hybrid max (tie breaker)
hybrid_epi <- hybrid_max
for (i in 1:nrow(hybrid_epi)){
  infectors <- which(hybrid_epi[i,]!=0)
  contacts <- which(DS$sampledata[,1] %in% neighbors(g,DS$sampledata[i,1]))
  overlap <- contacts[which(contacts %in% infectors)]
  if (length(overlap)>0){
    hybrid_epi[i,overlap] <- 1/length(overlap)
    hybrid_epi[i,infectors[which(!(infectors %in% contacts))]] <- 0
  }
}

ROC_hybrid_epi <- ROC_calc(hybrid_epi, truemat)
hybrid_epi_AUC_max_tiebreaker <- ROC_hybrid_epi$AUC

# incorporate epi data into dist (tie breaker)
dist_epi <- genclosenorm
for (i in 1:nrow(dist_epi)){
  infectors <- which(dist_epi[i,]!=0)
  contacts <- which(DS$sampledata[,1] %in% neighbors(g,DS$sampledata[i,1]))
  overlap <- contacts[which(contacts %in% infectors)]
  if (length(overlap)>0){
    dist_epi[i,overlap] <- 1/length(overlap)
    dist_epi[i,infectors[which(!(infectors %in% contacts))]] <- 0
  }
}

ROC_dist_epi <- ROC_calc(dist_epi, truemat)
dist_epi_AUC_max_tiebreaker <- ROC_dist_epi$AUC

# incorporate epi data into Wmat_max (tie breaker)
W_epi <- Wmat_max
for (i in 1:nrow(W_epi)){
  infectors <- which(W_epi[i,]!=0)
  contacts <- which(W$sampledata[,1] %in% neighbors(g,W$sampledata[i,1]))
  overlap <- contacts[which(contacts %in% infectors)]
  if (length(overlap)>0){
    W_epi[i,overlap] <- 1/length(overlap)
    W_epi[i,infectors[which(!(infectors %in% contacts))]] <- 0
  }
}

ROC_W_epi <- ROC_calc(W_epi, true)
W_epi_AUC_max_tiebreaker <- ROC_W_epi$AUC

# set cut offs
# scoremattrans[scoremattrans<0.25] <- 0
# hybrid_weight[hybrid_weight<0.25] <- 0
# gen_weight[gen_weight<0.25] <- 0
# scoremat_max_tie[scoremat_max_tie<0.25] <- 0
# scoremattrans_epi[scoremattrans_epi<0.25] <- 0


# recreate edge list from above and calculate RR for deep sequence data
for (df in 1:6){
  print(df)
  if (df==1){
    data <- scoremattrans
  } else if (df==2){
    data <- hybrid_weight
  } else if (df==3){
    data <- gen_weight
  } else if (df==4){
    data <- hybrid_epi
  } else if (df==5){
    data <- transmat2_epi
  } else {
    data <- dist_epi
  }
  results_edge_list_analysis_DS <- data.frame("InfectedNode"=NA,"node_stat"=NA,"Infector"=NA,"infector_stat"=NA,"prob"=NA,"max"=NA,"contacts"=NA)
  if (sum(data)>0){
    for (i in 1:nrow(scoremattrans)){
      # track all sources
      cat(i,"\n")
      closest_sources <- which(data[i,]!=0)
      for (inf_source in 1:length(closest_sources)){
        if (sum(closest_sources)==0){
          results_edge_list_analysis_DS <- 
            rbind(results_edge_list_analysis_DS,c(DS$sampledata[i],V(g)$trialstatus[DS$sampledata[i]],
                                                  0,NA,
                                                  0,0,0))
        } else{
          closest_sources_table <- as.data.frame(cbind(closest_sources,c(data[i,closest_sources])))
          names(closest_sources_table) <- c("closest_sources","prob")
          # top5_source_probs <- closest_sources_table[order(closest_sources_table$prob,decreasing=TRUE),][1:5,]
          # top5_source <- top5_source_probs[!is.na(top5_source_probs[,2]),]
          # top5_source <- top5_source[,1]
          # norm_prob <- sum(top5_source_probs[,2],na.rm=TRUE)
          # top5 <- ifelse(closest_sources[inf_source] %in% top5_source,top5_source_probs[which(top5_source_probs[,1]==closest_sources[inf_source]),2],0)
          max_sources <- which(data[i,] == max(data[i,]))
          max <- ifelse(closest_sources[inf_source] %in% max_sources,1/length(max_sources),0)
          contact <- ifelse(DS$sampledata[closest_sources[inf_source]] %in% neighbors(g,DS$sampledata[i,1]),1,0)
          results_edge_list_analysis_DS <- 
            rbind(results_edge_list_analysis_DS,c(DS$sampledata[i],V(g)$trialstatus[DS$sampledata[i]],
                                                  DS$sampledata[closest_sources[inf_source]],V(g)$trialstatus[DS$sampledata[closest_sources[inf_source]]],
                                                  data[i,closest_sources[inf_source]],max,contact))
        }
      }
    }
    results_edge_list_analysis_DS <- unique(results_edge_list_analysis_DS)
    
    # get average number of infectors
    results_edge_list_analysis_DS_prob <- results_edge_list_analysis_DS[results_edge_list_analysis_DS$prob!=0,]
    av_num_infectors_table <- aggregate(results_edge_list_analysis_DS_prob$prob,by=list(results_edge_list_analysis_DS_prob$InfectedNode),length)
    av_num_infectors <- mean(av_num_infectors_table$x,na.rm=TRUE)
    
    # note the true infector
    results_edge_list_analysis_true <- results_edge_list_analysis[1:nrow(results_infected),c(2,10)]
    names(results_edge_list_analysis_true) <- c("InfectedNode","TrueInfector")
    results_edge_list_analysis_DS_true <- merge(results_edge_list_analysis_DS_prob,results_edge_list_analysis_true,by="InfectedNode")
    results_edge_list_analysis_DS_true$truth <- ifelse(results_edge_list_analysis_DS_true$Infector==results_edge_list_analysis_DS_true$TrueInfector,1,0)
    mean_prob_truth <- mean(results_edge_list_analysis_DS_true$prob[results_edge_list_analysis_DS_true$truth==1],na.rm=TRUE)
    sd_prob_truth <- sd(results_edge_list_analysis_DS_true$prob[results_edge_list_analysis_DS_true$truth==1])
    
    # note true infector stat
    results_edge_list_analysis_DS_true$true_stat <- ifelse(results_edge_list_analysis_DS_true$TrueInfector==0,NA,V(g)$trialstatus[results_edge_list_analysis_DS_true$TrueInfector])
    
    # note if potential infector has same trial status as truth
    results_edge_list_analysis_DS_true$correctstat <- ifelse(results_edge_list_analysis_DS_true$infector_stat==results_edge_list_analysis_DS_true$true_stat,1,0)
    
    # sum of correct probabilities
    results_edge_list_analysis_DS_true$correctstatprob <- results_edge_list_analysis_DS_true$prob*results_edge_list_analysis_DS_true$correctstat
    correct_stat_probs <- aggregate(results_edge_list_analysis_DS_true$correctstatprob,by=list(results_edge_list_analysis_DS_true$InfectedNode),sum)
    names(correct_stat_probs) <- c("InfectedNode","CorrectProb")
    mean_prob_correct_stat <- mean(correct_stat_probs$CorrectProb,na.rm=TRUE)
    sd_prob_correct_stat <- sd(correct_stat_probs$CorrectProb,na.rm=TRUE)
    
    # sum of correct probabilities for maximum
    results_edge_list_analysis_DS_true$correctstatprobmax <- results_edge_list_analysis_DS_true$max*results_edge_list_analysis_DS_true$correctstat
    correct_stat_probs_max <- aggregate(results_edge_list_analysis_DS_true$correctstatprobmax,by=list(results_edge_list_analysis_DS_true$InfectedNode),sum)
    names(correct_stat_probs_max) <- c("InfectedNode","CorrectProb")
    mean_prob_correct_stat_max <- mean(correct_stat_probs_max$CorrectProb,na.rm=TRUE)
    sd_prob_correct_stat_max <- sd(correct_stat_probs_max$CorrectProb,na.rm=TRUE)
    
    #1# summarize by trial status of infector
    summary_infectors_DS_stat <- aggregate(results_edge_list_analysis_DS$InfectedNode,by=c(list(results_edge_list_analysis_DS$InfectedNode),list(results_edge_list_analysis_DS$infector_stat)),length)
    names(summary_infectors_DS_stat)<-c("InfectedNode","infector_stat","number")
    summary_infectors_DS <- aggregate(summary_infectors_DS_stat$number,by=list(summary_infectors_DS_stat$InfectedNode),sum)
    names(summary_infectors_DS)<-c("InfectedNode","total")
    summary_infectors_DS_merged <- merge(summary_infectors_DS_stat,summary_infectors_DS)
    summary_infectors_DS_merged$perc <- summary_infectors_DS_merged$number/summary_infectors_DS_merged$total
    summary_infectors_DS_final <- aggregate(summary_infectors_DS_merged$perc,by=list(summary_infectors_DS_merged$infector_stat),sum)
    names(summary_infectors_DS_final)<-c("infector_stat","number")
    RR_inf_DS <- summary_infectors_DS_final$number[summary_infectors_DS_final$infector_stat==1]/summary_infectors_DS_final$number[summary_infectors_DS_final$infector_stat==0]
    
    #2# look only at most likely
    summary_infectors_DS_max <- aggregate(results_edge_list_analysis_DS$max,by=list(results_edge_list_analysis_DS$infector_stat),sum)
    names(summary_infectors_DS_max)<-c("infector_stat","number")
    vacc_inf_DS_max <- summary_infectors_DS_max$number[summary_infectors_DS_max$infector_stat==1]
    cont_inf_DS_max<- summary_infectors_DS_max$number[summary_infectors_DS_max$infector_stat==0]
    RR_inf_DS_max <- vacc_inf_DS_max/cont_inf_DS_max
    
    #3# aggregate by probability instead of number of infectors
    summary_infectors_DS_prob <- aggregate(results_edge_list_analysis_DS$prob,by=list(results_edge_list_analysis_DS$infector_stat),sum)
    names(summary_infectors_DS_prob)<-c("infector_stat","prob")
    vacc_inf_DS_prob <- summary_infectors_DS_prob$prob[summary_infectors_DS_prob$infector_stat==1]
    cont_inf_DS_prob<- summary_infectors_DS_prob$prob[summary_infectors_DS_prob$infector_stat==0]
    RR_inf_DS_prob <- vacc_inf_DS_prob/cont_inf_DS_prob
    
    #4# add up probabilities for max likely
    results_edge_list_analysis_DS_max <- results_edge_list_analysis_DS[results_edge_list_analysis_DS$max!=0,]
    summary_infectors_DS_max_prob <- aggregate(results_edge_list_analysis_DS_max$prob,by=list(results_edge_list_analysis_DS_max$infector_stat),sum)
    names(summary_infectors_DS_max_prob)<-c("infector_stat","number")
    vacc_inf_DS_max_prob <- summary_infectors_DS_max_prob$number[summary_infectors_DS_max_prob$infector_stat==1]
    cont_inf_DS_max_prob<- summary_infectors_DS_max_prob$number[summary_infectors_DS_max_prob$infector_stat==0]
    RR_inf_DS_max_prob <- vacc_inf_DS_max_prob/cont_inf_DS_max_prob
    
    #5# look at maximum of total probabilities for infector status 
    # aggregate by infected node
    results_infector_stat_DS <- aggregate(results_edge_list_analysis_DS$prob,by=c(list(results_edge_list_analysis_DS$InfectedNode),list(results_edge_list_analysis_DS$infector_stat)),sum)
    names(results_infector_stat_DS) <- c("Infected Node","infector_stat","prob")
    # note which trial status is the maximum for each infected node
    results_infector_stat_DS$max <- ifelse(results_infector_stat_DS$prob>=0.5,1,0)
    results_infector_stat_DS_max <- results_infector_stat_DS[results_infector_stat_DS$max==1,]
    # assign 1 to trial status with most probability
    RR_inf_DS_max_class <- nrow(results_infector_stat_DS_max[results_infector_stat_DS_max$infector_stat==1,]) / nrow(results_infector_stat_DS_max[results_infector_stat_DS_max$infector_stat==0,])
    
    #6# assign total probability to trial status with most probability
    RR_inf_DS_prob_max_class <- sum(results_infector_stat_DS_max$prob[results_infector_stat_DS_max$infector_stat==1]) / sum(results_infector_stat_DS_max$prob[results_infector_stat_DS_max$infector_stat==0])
    #aggregate(results_edge_list_analysis_DS_max$prob,by=list(results_edge_list_analysis_DS_max$InfectedNode),sum,na.rm=TRUE)
    
    #7# look at maximum only among contacts -- max
    results_edge_list_analysis_DS_contacts <- results_edge_list_analysis_DS[results_edge_list_analysis_DS$contacts!=0 & !is.na(results_edge_list_analysis_DS$contacts),]
    ID_contacts <- unique(results_edge_list_analysis_DS_contacts$InfectedNode)
    for (i in 1:length(ID_contacts)){
      num_max_contacts <- length(which(results_edge_list_analysis_DS_contacts$prob[results_edge_list_analysis_DS_contacts$InfectedNode==ID_contacts[i]]==max(results_edge_list_analysis_DS_contacts$prob[results_edge_list_analysis_DS_contacts$InfectedNode==ID_contacts[i]],na.rm=TRUE)))
      results_edge_list_analysis_DS_contacts$max[results_edge_list_analysis_DS_contacts$InfectedNode==ID_contacts[i] & results_edge_list_analysis_DS_contacts$prob == max(results_edge_list_analysis_DS_contacts$prob[results_edge_list_analysis_DS_contacts$InfectedNode==ID_contacts[i]],na.rm=TRUE)] <-1/num_max_contacts
    }
    summary_infectors_DS_contacts <- aggregate(results_edge_list_analysis_DS_contacts$max,by=list(results_edge_list_analysis_DS_contacts$infector_stat),sum)
    names(summary_infectors_DS_contacts)<-c("infector_stat","number")
    vacc_inf_DS_contacts<- summary_infectors_DS_contacts$number[summary_infectors_DS_contacts$infector_stat==1]
    cont_inf_DS_contacts <- summary_infectors_DS_contacts$number[summary_infectors_DS_contacts$infector_stat==0]
    RR_inf_DS_contacts <- vacc_inf_DS_contacts/cont_inf_DS_contacts
    
    #8# look only among contacts -- weighted
    ID_contacts <- unique(results_edge_list_analysis_DS_contacts$InfectedNode)
    for (i in 1:length(ID_contacts)){
      results_edge_list_analysis_DS_contacts$prob[results_edge_list_analysis_DS_contacts$InfectedNode==ID_contacts[i]] <- results_edge_list_analysis_DS_contacts$prob[results_edge_list_analysis_DS_contacts$InfectedNode==ID_contacts[i]] / sum(results_edge_list_analysis_DS_contacts$prob[results_edge_list_analysis_DS_contacts$InfectedNode==ID_contacts[i]])
    }
    summary_infectors_DS_contacts_prob <- aggregate(results_edge_list_analysis_DS_contacts$prob,by=list(results_edge_list_analysis_DS_contacts$infector_stat),sum)
    names(summary_infectors_DS_contacts_prob)<-c("infector_stat","number")
    vacc_inf_DS_contacts_prob<- summary_infectors_DS_contacts_prob$number[summary_infectors_DS_contacts_prob$infector_stat==1]
    cont_inf_DS_contacts_prob <- summary_infectors_DS_contacts_prob$number[summary_infectors_DS_contacts_prob$infector_stat==0]
    RR_inf_DS_contacts_prob <- vacc_inf_DS_contacts_prob/cont_inf_DS_contacts_prob
  } else{
    RR_inf_DS <-NA 
    RR_inf_DS_max <-NA 
    RR_inf_DS_prob <-NA 
    RR_inf_DS_max_prob  <-NA 
    RR_inf_DS_contacts <-NA 
    RR_inf_DS_contacts_prob <-NA 
    RR_inf_DS_max_class <-NA 
    RR_inf_DS_prob_max_class <-NA 
    results_edge_list_analysis_DS <-NA 
    av_num_infectors <-NA 
    mean_prob_truth <-NA 
    sd_prob_truth <-NA 
    mean_prob_correct_stat <-NA 
    sd_prob_correct_stat <-NA 
    mean_prob_correct_stat_max <-NA 
    sd_prob_correct_stat_max <-NA 
  }
  
  if (df==1){
    RR_inf_DS_scoremat <- RR_inf_DS
    RR_inf_DS_scoremat_max <- RR_inf_DS_max
    RR_inf_DS_scoremat_prob <- RR_inf_DS_prob
    RR_inf_DS_scoremat_max_prob  <- RR_inf_DS_max_prob
    RR_inf_DS_scoremat_contacts <- RR_inf_DS_contacts
    RR_inf_DS_scoremat_contacts_prob <- RR_inf_DS_contacts_prob
    RR_inf_DS_scoremat_max_class <- RR_inf_DS_max_class
    RR_inf_DS_scoremat_prob_max_class <- RR_inf_DS_prob_max_class
    results_edge_list_analysis_DS_scoremat <- results_edge_list_analysis_DS
    av_num_infectors_scoremat <- av_num_infectors
    mean_prob_truth_scoremat <- mean_prob_truth
    sd_prob_truth_scoremat <- sd_prob_truth
    mean_prob_correct_stat_scoremat <- mean_prob_correct_stat
    sd_prob_correct_stat_scoremat <- sd_prob_correct_stat
    mean_prob_correct_stat_max_scoremat <- mean_prob_correct_stat_max
    sd_prob_correct_stat_max_scoremat <- sd_prob_correct_stat_max
  } else if (df==2){
    RR_inf_DS_hybrid <- RR_inf_DS
    RR_inf_DS_hybrid_max <- RR_inf_DS_max
    RR_inf_DS_hybrid_prob <- RR_inf_DS_prob
    RR_inf_DS_hybrid_max_prob  <- RR_inf_DS_max_prob
    RR_inf_DS_hybrid_contacts <- RR_inf_DS_contacts
    RR_inf_DS_hybrid_contacts_prob <- RR_inf_DS_contacts_prob
    RR_inf_DS_hybrid_max_class <- RR_inf_DS_max_class
    RR_inf_DS_hybrid_prob_max_class <- RR_inf_DS_prob_max_class
    results_edge_list_analysis_DS_hybrid <- results_edge_list_analysis_DS
    av_num_infectors_hybrid <- av_num_infectors
    mean_prob_truth_hybrid <- mean_prob_truth
    sd_prob_truth_hybrid <- sd_prob_truth
    mean_prob_correct_stat_hybrid <- mean_prob_correct_stat
    sd_prob_correct_stat_hybrid <- sd_prob_correct_stat
    mean_prob_correct_stat_max_hybrid <- mean_prob_correct_stat_max
    sd_prob_correct_stat_max_hybrid <- sd_prob_correct_stat_max
  } else if (df==3){
    RR_inf_DS_dist <- RR_inf_DS
    RR_inf_DS_dist_max <- RR_inf_DS_max
    RR_inf_DS_dist_prob <- RR_inf_DS_prob
    RR_inf_DS_dist_max_prob  <- RR_inf_DS_max_prob
    RR_inf_DS_dist_contacts <- RR_inf_DS_contacts
    RR_inf_DS_dist_contacts_prob <- RR_inf_DS_contacts_prob
    RR_inf_DS_dist_max_class <- RR_inf_DS_max_class
    RR_inf_DS_dist_prob_max_class <- RR_inf_DS_prob_max_class
    results_edge_list_analysis_DS_dist <- results_edge_list_analysis_DS
    av_num_infectors_dist <- av_num_infectors
    mean_prob_truth_dist <- mean_prob_truth
    sd_prob_truth_dist <- sd_prob_truth
    mean_prob_correct_stat_dist <- mean_prob_correct_stat
    sd_prob_correct_stat_dist <- sd_prob_correct_stat
    mean_prob_correct_stat_max_dist <- mean_prob_correct_stat_max
    sd_prob_correct_stat_max_dist <- sd_prob_correct_stat_max
  } else if (df==4){
    RR_inf_DS_hybrid_epi_data <- RR_inf_DS
    RR_inf_DS_hybrid_epi_data_max <- RR_inf_DS_max
    RR_inf_DS_hybrid_epi_data_prob <- RR_inf_DS_prob
    RR_inf_DS_hybrid_epi_data_max_prob  <- RR_inf_DS_max_prob
    RR_inf_DS_hybrid_epi_data_contacts <- RR_inf_DS_contacts
    RR_inf_DS_hybrid_epi_data_contacts_prob <- RR_inf_DS_contacts_prob
    RR_inf_DS_hybrid_epi_data_max_class <- RR_inf_DS_max_class
    RR_inf_DS_hybrid_epi_data_prob_max_class <- RR_inf_DS_prob_max_class
    results_edge_list_analysis_DS_hybrid_epi_data <- results_edge_list_analysis_DS
    av_num_infectors_hybrid_epi_data <- av_num_infectors
    mean_prob_truth_hybrid_epi_data <- mean_prob_truth
    sd_prob_truth_hybrid_epi_data<- sd_prob_truth
    mean_prob_correct_stat_hybrid_epi_data <- mean_prob_correct_stat
    sd_prob_correct_stat_hybrid_epi_data <- sd_prob_correct_stat
    mean_prob_correct_stat_max_hybrid_epi_data <- mean_prob_correct_stat_max
    sd_prob_correct_stat_max_hybrid_epi_data <- sd_prob_correct_stat_max
  } else if (df==5){
    RR_inf_scoremat_epi_data <- RR_inf_DS
    RR_inf_scoremat_epi_data_max <- RR_inf_DS_max
    RR_inf_scoremat_epi_data_prob <- RR_inf_DS_prob
    RR_inf_scoremat_epi_data_max_prob  <- RR_inf_DS_max_prob
    RR_inf_scoremat_epi_data_contacts <- RR_inf_DS_contacts
    RR_inf_scoremat_epi_data_contacts_prob <- RR_inf_DS_contacts_prob
    RR_inf_scoremat_epi_data_max_class <- RR_inf_DS_max_class
    RR_inf_scoremat_epi_data_prob_max_class <- RR_inf_DS_prob_max_class
    results_edge_list_analysis_scoremat_epi_data <- results_edge_list_analysis_DS
    av_num_infectors_epi_data <- av_num_infectors
    mean_prob_truth_epi_data <- mean_prob_truth
    sd_prob_truth_epi_data<- sd_prob_truth
    mean_prob_correct_stat_epi_data <- mean_prob_correct_stat
    sd_prob_correct_stat_epi_data <- sd_prob_correct_stat
    mean_prob_correct_stat_max_epi_data <- mean_prob_correct_stat_max
    sd_prob_correct_stat_max_epi_data <- sd_prob_correct_stat_max
  } else {
    RR_inf_dist_epi_data <- RR_inf_DS
    RR_inf_dist_epi_data_max <- RR_inf_DS_max
    RR_inf_dist_epi_data_prob <- RR_inf_DS_prob
    RR_inf_dist_epi_data_max_prob  <- RR_inf_DS_max_prob
    RR_inf_dist_epi_data_contacts <- RR_inf_DS_contacts
    RR_inf_dist_epi_data_contacts_prob <- RR_inf_DS_contacts_prob
    RR_inf_dist_epi_data_max_class <- RR_inf_DS_max_class
    RR_inf_dist_epi_data_prob_max_class <- RR_inf_DS_prob_max_class
    results_edge_list_analysis_dist_epi_data <- results_edge_list_analysis_DS
    av_num_infectors_dist_epi_data <- av_num_infectors
    mean_prob_truth_dist_epi_data <- mean_prob_truth
    sd_prob_truth_dist_epi_data<- sd_prob_truth
    mean_prob_correct_stat_dist_epi_data <- mean_prob_correct_stat
    sd_prob_correct_stat_dist_epi_data <- sd_prob_correct_stat
    mean_prob_correct_stat_max_dist_epi_data <- mean_prob_correct_stat_max
    sd_prob_correct_stat_max_dist_epi_data <- sd_prob_correct_stat_max
  }
  
  if (!is.na(sum(results_edge_list_analysis_DS_scoremat))){
    # restrict scoremat to only 1 max
    results_edge_list_analysis_DS_max_summary <- aggregate(results_edge_list_analysis_DS_scoremat$max,by=c(list(results_edge_list_analysis_DS_scoremat$InfectedNode),list(results_edge_list_analysis_DS_scoremat$infector_stat)),sum)
    names(results_edge_list_analysis_DS_scoremat_max_summary)<-c("InfectedNode","infector_stat","max")
    results_edge_list_analysis_DS_scoremat_max <- results_edge_list_analysis_DS_scoremat_max_summary[results_edge_list_analysis_DS_scoremat_max_summary$max==1,]
    summary_infectors_DS_max <- aggregate(results_edge_list_analysis_DS_scoremat_max$max,by=list(results_edge_list_analysis_DS_scoremat_max$infector_stat),sum)
    names(summary_infectors_DS_max)<-c("infector_stat","number")
    vacc_inf_DS_max <- summary_infectors_DS_max$number[summary_infectors_DS_max$infector_stat==1]
    cont_inf_DS_max<- summary_infectors_DS_max$number[summary_infectors_DS_max$infector_stat==0]
    if (length(vacc_inf_DS_max/cont_inf_DS_max)>0){
      RR_inf_SM_max <- vacc_inf_DS_max/cont_inf_DS_max
    } else {
      RR_inf_SM_max <- NA
    }
    
    # restrict scoremat to only 1 max among contacts
    results_edge_list_analysis_scoremat_contacts <- results_edge_list_analysis_DS_scoremat[results_edge_list_analysis_DS_scoremat$contacts!=0 & !is.na(results_edge_list_analysis_DS_scoremat$contacts),]
    ID_contacts <- unique(results_edge_list_analysis_scoremat_contacts$InfectedNode)
    for (i in 1:length(ID_contacts)){
      num_max_contacts <- length(which(results_edge_list_analysis_scoremat_contacts$prob[results_edge_list_analysis_scoremat_contacts$InfectedNode==ID_contacts[i]]==max(results_edge_list_analysis_scoremat_contacts$prob[results_edge_list_analysis_scoremat_contacts$InfectedNode==ID_contacts[i]],na.rm=TRUE)))
      results_edge_list_analysis_scoremat_contacts$max[results_edge_list_analysis_scoremat_contacts$InfectedNode==ID_contacts[i] & results_edge_list_analysis_scoremat_contacts$prob == max(results_edge_list_analysis_scoremat_contacts$prob[results_edge_list_analysis_scoremat_contacts$InfectedNode==ID_contacts[i]],na.rm=TRUE)] <-1/num_max_contacts
    }
    
    results_edge_list_analysis_scoremat_max_summary <- aggregate(results_edge_list_analysis_scoremat_contacts$max,by=c(list(results_edge_list_analysis_scoremat_contacts$InfectedNode),list(results_edge_list_analysis_scoremat_contacts$infector_stat)),sum)
    names(results_edge_list_analysis_scoremat_max_summary)<-c("InfectedNode","infector_stat","max")
    results_edge_list_analysis_scoremat_max <- results_edge_list_analysis_scoremat_max_summary[results_edge_list_analysis_scoremat_max_summary$max==1,]
    summary_infectors_scoremat_max <- aggregate(results_edge_list_analysis_scoremat_max$max,by=list(results_edge_list_analysis_scoremat_max$infector_stat),sum)
    names(summary_infectors_scoremat_max)<-c("infector_stat","number")
    vacc_inf_scoremat_max_contacts <- summary_infectors_scoremat_max$number[summary_infectors_scoremat_max$infector_stat==1]
    cont_inf_scoremat_max_contacts<- summary_infectors_scoremat_max$number[summary_infectors_scoremat_max$infector_stat==0]
    if (length(vacc_inf_scoremat_max_contacts/cont_inf_scoremat_max_contacts)>0){
      RR_inf_SM_max_contacts <- vacc_inf_scoremat_max_contacts/cont_inf_scoremat_max_contacts
    } else {
      RR_inf_SM_max_contacts <- NA
    }
  } else{
    RR_inf_SM_max <- NA
    RR_inf_SM_max_contacts <- NA
  }
}

# VEs
numevents_vacc <- nrow(results_edge_list_analysis[(results_edge_list_analysis$eventstatus==1) & (results_edge_list_analysis$TrialStatus==1),])
numevents_cont <- nrow(results_edge_list_analysis[(results_edge_list_analysis$eventstatus==1) & (results_edge_list_analysis$TrialStatus==0),])
num_vacc <- nrow(results_edge_list_analysis[(results_edge_list_analysis$TrialStatus==1),])
num_cont <- nrow(results_edge_list_analysis[(results_edge_list_analysis$TrialStatus==0),])  

total_vacc_pt <- sum(results_edge_list_analysis$DayInfected[results_edge_list_analysis$TrialStatus==1])
total_cont_pt <- sum(results_edge_list_analysis$DayInfected[results_edge_list_analysis$TrialStatus==0])
VEs <- 1 - (numevents_vacc/num_vacc)/(numevents_cont/num_cont)

# VEI
summary_infectors <- aggregate(results_infected$InfectedNode,by=list(results_infected$infector_stat),length)
names(summary_infectors)<-c("infector_stat","number")
vacc_inf <- summary_infectors$number[summary_infectors$infector_stat==1]
cont_inf<- summary_infectors$number[summary_infectors$infector_stat==0]
RR_inf <- vacc_inf/cont_inf

# VEI_contacts
results_contacts <- as.data.frame(cbind("InfectedNode"=NA,"Contact"=NA,"infector_stat"=NA))
infected_nodes <- results_infected$InfectedNode
for (i in 1:nrow(results_infected)){
  infected_contacts <- infected_nodes[which(infected_nodes %in% neighbors(g,results_infected[i,2]))]
  cutoff_time_high <- DS$epidata[which(DS$epidata[,1] == infected_nodes[i]),2]
  #restrict to 2 serial intervals before
  cutoff_time_low <- cutoff_time_high - 2*(infperiod_shape/infperiod_rate)
  possible_infectors <- DS$epidata[DS$epidata[,2]>=cutoff_time_low & DS$epidata[,2]<cutoff_time_high,1]
  infected_contacts_cutofftime <- intersect(infected_contacts,possible_infectors)
    if (length(infected_contacts_cutofftime >0)){
      for (m in 1:length(infected_contacts_cutofftime)){
        contact_info <- as.data.frame(cbind(results_infected[i,2],infected_contacts_cutofftime[m],V(g)$trialstatus[infected_contacts_cutofftime[m]]))
        names(contact_info) <- c("InfectedNode","Contact","infector_stat")
        results_contacts <- rbind(results_contacts,contact_info)
      }
    }
}
summary_infectors_contacts <- aggregate(results_contacts$InfectedNode,by=list(results_contacts$infector_stat),length)
names(summary_infectors_contacts)<-c("infector_stat","number")
vacc_inf_contacts<- summary_infectors_contacts$number[summary_infectors_contacts$infector_stat==1]
cont_inf_contacts <- summary_infectors_contacts$number[summary_infectors_contacts$infector_stat==0]
RR_inf_contacts <- vacc_inf_contacts/cont_inf_contacts


# create a matrix with contacts
n <- nrow(results_infected)
contacts_mat <- matrix(0,n,n)
for (i in 1:n){
  infected_contacts <- DS$sampledata[which(DS$sampledata[,1] %in% neighbors(g,DS$sampledata[i,1])),1]
  cutoff_time_high <- DS$epidata[which(DS$epidata[,1] %in% DS$sampledata[i]),2]
  #restrict to 2 serial intervals before
  cutoff_time_low <- cutoff_time_high - 2*(infperiod_shape/infperiod_rate)
  possible_infectors <- DS$epidata[DS$epidata[,2]>=cutoff_time_low & DS$epidata[,2]<cutoff_time_high,1]
  infected_contacts_cutofftime <- intersect(infected_contacts,possible_infectors)
  contacts_mat[i,which(DS$sampledata[,1] %in% infected_contacts_cutofftime)] <- 1
}

# normalize
contacts_norm <- matrix(0,n,n) # normalised score matrix
for (i in 1:n){
  if (sum(contacts_mat[i,])>0) {
    contacts_norm[i,] <- contacts_mat[i,]/sum(contacts_mat[i,])
  }
}

ROC_contacts <-ROC_calc(contacts_norm,truemat)
contacts_AUC <- ROC_contacts$AUC

# calculate VEIs from ratios of infector statuses
VEIs <- 1-RR_inf/(1-VEs) 
VEIs_contacts <- 1-RR_inf_contacts/(1-VEs)
VEIs_W <-1-RR_inf_W/(1-VEs)
VEIs_W_max <-1-RR_inf_W_max/(1-VEs)
VEIs_W_prob <-1-RR_inf_W_prob/(1-VEs)
VEIs_W_max_prob <-1-RR_inf_W_max_prob/(1-VEs)
VEIs_W_contacts<-1-RR_inf_W_contacts/(1-VEs)
VEIs_W_contacts_prob<-1-RR_inf_W_contacts_prob/(1-VEs)
VEIs_W_max_class <- 1-RR_inf_W_max_class/(1-VEs)
VEIs_W_prob_max_class <- 1-RR_inf_W_prob_max_class/(1-VEs)
VEIs_W_epi_data <- 1-RR_inf_W_max_contacts/(1-VEs)

VEIs_scoremat <-1-RR_inf_DS_scoremat/(1-VEs)
VEIs_scoremat_max <-1-RR_inf_DS_scoremat_max/(1-VEs)
VEIs_scoremat_prob <-1-RR_inf_DS_scoremat_prob/(1-VEs)
VEIs_scoremat_max_prob <-1-RR_inf_DS_scoremat_max_prob/(1-VEs)
VEIs_scoremat_contacts <- 1-RR_inf_DS_scoremat_contacts/(1-VEs)
VEIs_scoremat_contacts_prob <- 1-RR_inf_DS_scoremat_contacts_prob/(1-VEs)
VEIs_scoremat_max_class <- 1-RR_inf_DS_scoremat_max_class/(1-VEs)
VEIs_scoremat_prob_max_class <- 1- RR_inf_DS_scoremat_prob_max_class/(1-VEs)
VEIs_scoremat_max1only <- 1- RR_inf_SM_max/(1-VEs)
VEIs_scoremat_max1only_contacts <- 1-RR_inf_SM_max_contacts/(1-VEs)

VEIs_hybrid <-1-RR_inf_DS_hybrid/(1-VEs)
VEIs_hybrid_max <-1-RR_inf_DS_hybrid_max/(1-VEs)
VEIs_hybrid_prob <-1-RR_inf_DS_hybrid_prob/(1-VEs)
VEIs_hybrid_max_prob <-1-RR_inf_DS_hybrid_max_prob/(1-VEs)
VEIs_hybrid_contacts <- 1-RR_inf_DS_hybrid_contacts/(1-VEs)
VEIs_hybrid_contacts_prob <- 1-RR_inf_DS_hybrid_contacts_prob/(1-VEs)
VEIs_hybrid_max_class <- 1-RR_inf_DS_hybrid_max_class/(1-VEs)
VEIs_hybrid_prob_max_class <- 1- RR_inf_DS_hybrid_prob_max_class/(1-VEs)

VEIs_dist <-1-RR_inf_DS_dist/(1-VEs)
VEIs_dist_max <-1-RR_inf_DS_dist_max/(1-VEs)
VEIs_dist_prob <-1-RR_inf_DS_dist_prob/(1-VEs)
VEIs_dist_max_prob <-1-RR_inf_DS_dist_max_prob/(1-VEs)
VEIs_dist_contacts <- 1-RR_inf_DS_dist_contacts/(1-VEs)
VEIs_dist_contacts_prob <- 1-RR_inf_DS_dist_contacts_prob/(1-VEs)
VEIs_dist_max_class <- 1-RR_inf_DS_dist_max_class/(1-VEs)
VEIs_dist_prob_max_class <- 1- RR_inf_DS_dist_prob_max_class/(1-VEs)

VEIs_hybrid_epi_data <-1-RR_inf_DS_hybrid_epi_data/(1-VEs)
VEIs_hybrid_epi_data_max <-1-RR_inf_DS_hybrid_epi_data_max/(1-VEs)
VEIs_hybrid_epi_data_prob <-1-RR_inf_DS_hybrid_epi_data_prob/(1-VEs)
VEIs_hybrid_epi_data_max_prob <-1-RR_inf_DS_hybrid_epi_data_max_prob/(1-VEs)
VEIs_hybrid_epi_data_contacts <- 1-RR_inf_DS_hybrid_epi_data_contacts/(1-VEs)
VEIs_hybrid_epi_data_contacts_prob <- 1-RR_inf_DS_hybrid_epi_data_contacts_prob/(1-VEs)
VEIs_hybrid_epi_data_max_class <- 1-RR_inf_DS_hybrid_epi_data_max_class/(1-VEs)
VEIs_hybrid_epi_data_prob_max_class <- 1- RR_inf_DS_hybrid_epi_data_prob_max_class/(1-VEs)

VEIs_scoremat_epi_data <-1-RR_inf_scoremat_epi_data/(1-VEs)
VEIs_scoremat_epi_data_max <-1-RR_inf_scoremat_epi_data_max/(1-VEs)
VEIs_scoremat_epi_data_prob <-1-RR_inf_scoremat_epi_data_prob/(1-VEs)
VEIs_scoremat_epi_data_max_prob <-1-RR_inf_scoremat_epi_data_max_prob/(1-VEs)
VEIs_scoremat_epi_data_contacts <- 1-RR_inf_scoremat_epi_data_contacts/(1-VEs)
VEIs_scoremat_epi_data_contacts_prob <- 1-RR_inf_scoremat_epi_data_contacts_prob/(1-VEs)
VEIs_scoremat_epi_data_max_class <- 1-RR_inf_scoremat_epi_data_max_class/(1-VEs)
VEIs_scoremat_epi_data_prob_max_class <- 1- RR_inf_scoremat_epi_data_prob_max_class/(1-VEs)

VEIs_dist_epi_data <-1-RR_inf_dist_epi_data/(1-VEs)
VEIs_dist_epi_data_max <-1-RR_inf_dist_epi_data_max/(1-VEs)
VEIs_dist_epi_data_prob <-1-RR_inf_dist_epi_data_prob/(1-VEs)
VEIs_dist_epi_data_max_prob <-1-RR_inf_dist_epi_data_max_prob/(1-VEs)
VEIs_dist_epi_data_contacts <- 1-RR_inf_dist_epi_data_contacts/(1-VEs)
VEIs_dist_epi_data_contacts_prob <- 1-RR_inf_dist_epi_data_contacts_prob/(1-VEs)
VEIs_dist_epi_data_max_class <- 1-RR_inf_dist_epi_data_max_class/(1-VEs)
VEIs_dist_epi_data_prob_max_class <- 1- RR_inf_dist_epi_data_prob_max_class/(1-VEs)

numevents_vacc <- nrow(results_edge_list_analysis[(results_edge_list_analysis$eventstatus==1) & (results_edge_list_analysis$TrialStatus==1),])
numevents_cont <- nrow(results_edge_list_analysis[(results_edge_list_analysis$eventstatus==1) & (results_edge_list_analysis$TrialStatus==0),])

args=(commandArgs(TRUE))
for (i in 1:length(args)) {
  eval (parse (text = args[[i]] ))
}
j<-j

genome_results <- as.data.frame(cbind(j,infect_VE,direct_VE,R0,mut.rate,bn,numevents_cont,numevents_vacc,
                 pairwise_AUC,pairwise_AUC_max,scoremattrans_AUC,transmat2_AUC,transmat2_1_AUC,
                 close_AUC,weight_AUC,hybrid_weight_AUC,hybrid_max_AUC,
                 contacts_AUC,
                 scoremattrans_epi_AUC_restrict,scoremattrans_epi_AUC_max_restrict,scoremattrans_epi_AUC_max_restrict2,
                 Wmat_norm_epi_AUC_restrict,Wmat_norm_epi_AUC_max_restrict,dist_epi_AUC_restrict,dist_epi_AUC_max_restrict,hybrid_epi_AUC_restrict,hybrid_epi_AUC_max_restrict,
                 W_epi_AUC_max_tiebreaker,hybrid_epi_AUC_max_tiebreaker,dist_epi_AUC_max_tiebreaker,transmat2_epi_AUC_tiebreaker,
                 VEs,VEIs,VEIs_contacts,
                 VEIs_W,VEIs_W_max,VEIs_W_prob,VEIs_W_max_prob,VEIs_W_contacts,VEIs_W_contacts_prob,VEIs_W_max_class,VEIs_W_prob_max_class,VEIs_W_epi_data,
                 VEIs_scoremat,VEIs_scoremat_max,VEIs_scoremat_prob,VEIs_scoremat_max_prob,VEIs_scoremat_max_class,VEIs_scoremat_prob_max_class,VEIs_scoremat_max1only,VEIs_scoremat_max1only_contacts,
                 VEIs_hybrid,VEIs_hybrid_max,VEIs_hybrid_prob,VEIs_hybrid_max_prob,VEIs_hybrid_max_class,VEIs_hybrid_prob_max_class,
                 VEIs_dist,VEIs_dist_max,VEIs_dist_prob,VEIs_dist_max_prob,VEIs_dist_max_class,VEIs_dist_prob_max_class,
                 VEIs_scoremat_epi_data,VEIs_scoremat_epi_data_max,VEIs_scoremat_epi_data_prob,VEIs_scoremat_epi_data_max_prob,VEIs_scoremat_epi_data_max_class,VEIs_scoremat_epi_data_prob_max_class,
                 VEIs_hybrid_epi_data,VEIs_hybrid_epi_data_max,VEIs_hybrid_epi_data_prob,VEIs_hybrid_epi_data_max_prob,VEIs_hybrid_epi_data_max_class,VEIs_hybrid_epi_data_prob_max_class,
                 VEIs_dist_epi_data,VEIs_dist_epi_data_max,VEIs_dist_epi_data_prob,VEIs_dist_epi_data_max_prob,VEIs_dist_epi_data_max_class,VEIs_dist_epi_data_prob_max_class,
                 VEIs_scoremat_contacts,VEIs_hybrid_contacts,VEIs_dist_contacts,VEIs_scoremat_epi_data_contacts,VEIs_hybrid_epi_data_contacts,VEIs_dist_epi_data_contacts,
                 VEIs_scoremat_contacts_prob,VEIs_hybrid_contacts_prob,VEIs_dist_contacts_prob,VEIs_scoremat_epi_data_contacts_prob,VEIs_hybrid_epi_data_contacts_prob,VEIs_dist_epi_data_contacts_prob,
                 av_num_infectors_W,mean_prob_truth_W,sd_prob_truth_W,mean_prob_correct_stat_W,sd_prob_correct_stat_W,
                 av_num_infectors_scoremat,mean_prob_truth_scoremat,sd_prob_truth_scoremat,mean_prob_correct_stat_scoremat,sd_prob_correct_stat_scoremat,
                 av_num_infectors_hybrid,mean_prob_truth_hybrid,sd_prob_truth_hybrid,mean_prob_correct_stat_hybrid,sd_prob_correct_stat_hybrid,
                 av_num_infectors_dist,mean_prob_truth_dist,sd_prob_truth_dist,mean_prob_correct_stat_dist,sd_prob_correct_stat_dist,
                 mean_prob_correct_stat_max_scoremat,mean_prob_correct_stat_max_hybrid,mean_prob_correct_stat_max_dist,
                 sd_prob_correct_stat_max_scoremat,sd_prob_correct_stat_max_hybrid,sd_prob_correct_stat_max_dist,
                 mean_prob_correct_stat_max_W,sd_prob_correct_stat_max_W,
                 av_num_infectors_hybrid_epi_data,mean_prob_truth_hybrid_epi_data,sd_prob_truth_hybrid_epi_data,mean_prob_correct_stat_hybrid_epi_data,
                 sd_prob_correct_stat_hybrid_epi_data,mean_prob_correct_stat_max_hybrid_epi_data,sd_prob_correct_stat_max_hybrid_epi_data,
                 av_num_infectors_epi_data,mean_prob_truth_epi_data,sd_prob_truth_epi_data,mean_prob_correct_stat_epi_data,
                 sd_prob_correct_stat_epi_data,mean_prob_correct_stat_max_epi_data,sd_prob_correct_stat_max_epi_data,
                 av_num_infectors_dist_epi_data,mean_prob_truth_dist_epi_data,sd_prob_truth_dist_epi_data,mean_prob_correct_stat_dist_epi_data,
                 sd_prob_correct_stat_dist_epi_data,mean_prob_correct_stat_max_dist_epi_data,sd_prob_correct_stat_max_dist_epi_data))

write.csv(genome_results,paste0(j,'_genome_results.csv'))
# write.csv(results_edge_list_analysis_DS_scoremat, paste0(j,"_results_edge_list_analysis_DS_scoremat",R0,"_",infect_VE,"_",direct_VE,"_",mut.rate,"_",bn,"_",infperiod_shape,"_",sample_percent,"_","sample.csv"))
# write.csv(results_edge_list_analysis_W, paste0(j,"_results_edge_list_analysis_W",R0,"_",infect_VE,"_",direct_VE,"_",mut.rate,"_",bn,"_",infperiod_shape,"_",sample_percent,"_","sample.csv"))
# write.csv(results_edge_list_analysis, paste0(j,"_results_edge_list_analysis",R0,"_",infect_VE,"_",direct_VE,"_",mut.rate,"_",bn,"_",infperiod_shape,"_",sample_percent,"_","sample.csv"))
# write.csv(results_edge_list_analysis_DS_hybrid, paste0(j,"_results_edge_list_analysis_hybrid",R0,"_",infect_VE,"_",direct_VE,"_",mut.rate,"_",bn,"_",infperiod_shape,"_",sample_percent,"_","sample.csv"))
# write.csv(results_edge_list_analysis_DS_dist, paste0(j,"_results_edge_list_analysis_dist",R0,"_",infect_VE,"_",direct_VE,"_",mut.rate,"_",bn,"_",infperiod_shape,"_",sample_percent,"_","sample.csv"))
