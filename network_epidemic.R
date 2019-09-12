# make a network and simulate an outbreak
# iRCT
# all or nothing vaccine
# written by Matt Hitchings and adapted by Rebecca Kahn

setwd("/n/home00/rkahn/seq")
require(seedy)
require(ggplot2)

betas <- c(0.0027,0.0029,0.003,0.0032,0.0035,0.0045,0.006)
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


make_network <- function(ave_community_size, community_size_range, num_communities, rate_within, rate_between) {
  # Function to make a network of closely-connected communities that are more sparsely
  # connected with each other. I use a stochastic block model.
  # Inputs:
  # Ave_community_size is the average number of people in a community
  # Community_size_range is the range of community sizes. Currently size of a community
  # is being chosen from a uniform distribution on (ave-range/2, ave+range/2)
  # Num_communities is the number of communities in the study population
  # Note that as the code stands the study population size is not fixed. 
  # rate_within is the probability of an edge between any two nodes within the same community
  # rate_between is the probability of an edge between any two nodes in different communities
  
  require(NetSurv)
  require(Matrix)
  require(Rlab)
  require(igraph)
  require(deSolve)
  require(reshape2)
  require(ggplot2)
  
  # Create the network, and assign all members a community number
  community_sizes <- ave_community_size + round(runif(num_communities,-community_size_range/2,community_size_range/2))
  studypop_size <- sum(community_sizes)
  # Currently all communities have the same connectedness, and all communities are equally
  # connected to each other
  within_rates <- diag(nrow=num_communities,ncol=num_communities,x=rate_within)
  between_rates <- matrix(rate_between,nrow=num_communities,ncol=num_communities) -
    diag(nrow=num_communities,ncol=num_communities,x=rate_between)
  rates<-within_rates+between_rates
  
  g <- sample_sbm(studypop_size,rates,community_sizes)
  # Give the nodes a name so that igraph remembers them
  V(g)$name<-1:studypop_size
  V(g)$community<-rep(1:num_communities,community_sizes)
  # Trial status will track whether a node is not in the trial (NA), in the control arm (0) or
  # in the vaccine arm (1)
  ## Enrollment day is the day a node is enrolled into the trial
  ## Symptomatic will track whether they are symptomatic or asymptomatically infected
  V(g)$trialstatus<-NA
  V(g)$enrollmentday<-NA
  V(g)$eventstatus<-0
  
  return(g)
  
}


network_epidemic<-function(g,beta,num_introductions,VE,VEI,edge_list,
                           infperiod_shape,infperiod_rate,
                           bTrial,bCluster,
                           trial_startday,trial_length,
                           num_enrolled_per_day,enrollment_period,cluster_coverage,simnum) {
  # Inputs:
  # g - the graph to run the epidemic on
  # beta - Every infectious individual contacts all their neighbours in a time step
  # and infects each susceptible with hazard beta. So beta represents the hazard of infection from one
  # contact between an infectious and a susceptible.
  # num_introductions - how many separate introductions we expect on average from the main epidemic. This is used
  # to calibrate the external force of infection
  # VE - direct leaky efficacy of the vaccine
  # bTrial - whether we are running a trial or not
  # bCluster - indicator of whether we are running the cRCT (1) or the iRCT (0)
  # trial_startday - first day of trial enrollment
  # trial_length - end of follow-up of trial partcipants, counting from the first day of enrollment
  # num_enrolled_per_day - number of individuals/clusters enrolled per day
  # enrollment_period - length of enrollment period
  # cluster_coverage - The proportion of each cluster we expect to enroll
  
  require(NetSurv)
  require(Matrix)
  require(Rlab)
  require(igraph)
  require(deSolve)
  require(reshape2)
  require(ggplot2)
  require(caTools)
  
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
  
  # Recover and spread functions
  recover<-function(i_nodes,r_nodes,infperiod_shape,infperiod_rate,t,results) {
    # Input is a list of the infected nodes, 
    # with number of days since infection and total infectious
    # period.
    # For each of these nodes, we will add it to recovered if the number of days infected has
    # reached the total length of the infectious period.
    
    # Advance infectious nodes 
    indices_to_remove <- i_nodes[2,]>=i_nodes[3,]
    newremoved<-as.vector(i_nodes[1,indices_to_remove])
    
    # Add one day to the length of each infected individual's time infected
    i_nodes[2,] <- i_nodes[2,]+1
    
    # Remove any recovered from i_nodes and add to r_nodes
    i_nodes <- i_nodes[,!(i_nodes[1,] %in% newremoved),drop=FALSE]
    r_nodes <- c(r_nodes,newremoved)
    results$DayRecovered[results$InfectedNode %in% newremoved] <-t
    
    list(i_nodes, r_nodes,results)
  }
  
  spread<-function(g, s_nodes, v_nodes,i_nodes, edge_list,t,results,
                   beta, VE, VEI,
                   connected_nodes,external_inf_F,source_num_inf){
    # Spread will create new infected nodes from two sources: infectious nodes within the the study
    # population, and external pressure from the source population
    # Inputs:
    # g is the graph, used to find neighbours of infected nodes
    # s_nodes, e_nodes and i_nodes are susceptible, infected nodes
    # beta is the hazard of infection for one contact
    # connected_nodes is a list of nodes that are connected the the source population
    # external_inf_F is a constant of proportionality that defines infectious pressure from source population 
    # to an individual
    # source_num_inf is the number of infectious individuals in the source population
    
    # Process: go through list of i_nodes, and choose a random number of its susceptible neighbours to
    # be infected, according to beta and choose a random number of its susceptible vaccinated neighbours to
    # be infected, according to beta and VE
    # Then go through list of nodes that are connected to source population and infect each susceptible
    # one with probability 1-exp(-FI), where I is the number/proportion of infectious, and F is a constant
    
    if (ncol(i_nodes)>0) {
      
      # make betas
      beta_unvacc <- beta
      beta_vacc <- beta_unvacc*(1-VEI)
      
      vacc <- V(g)[trialstatus==1 & !is.na(trialstatus)]$name
      
      i_nodes_unvacc <- setdiff(i_nodes[1,],vacc)
      i_nodes_vacc <- intersect(i_nodes[1,],vacc)
      #cat(t,length(i_nodes_unvacc),length(i_nodes_vacc),"\n")
      
      # start with vaccinated infected nodes
      if (length(i_nodes_vacc)>0){
        # Get a list of all unvaccinated neighbours of all infected vaccinated nodes
        potential_contacts_vacc<-lapply(i_nodes_vacc,function(x) neighbors(g,x))
        susc_contacts_vacc<-lapply(potential_contacts_vacc,function(x,susceptibles) intersect(x,susceptibles),susceptibles=s_nodes)
        num_neighbours_susc_vacc<-rapply(susc_contacts_vacc,length)
        # Sample from each group of neighbours in turn
        # First choose how many neighbours each node infects
        num_contacts_susc_vacc<-rbinom(length(num_neighbours_susc_vacc),num_neighbours_susc_vacc,1-exp(-beta_vacc))
        # Then sample from the neighbours
        # If one node gets picked twice by different nodes, just discard the duplicate.
        # In the very rare case that each i_nodes makes a number of new infectees equal to the number
        # of infectious nodes, mapply will make a matrix and unlist won't work. Therefore, the c() around
        # it ensures we turn the matrix into a vector. Unique then removes duplicates.
        infectees_susc_vacc<-c(unlist(mapply(function(x,y) x[sample.int(length(x),y)],x=susc_contacts_vacc,y=num_contacts_susc_vacc)))
        infectees_susc_vacc<-unique(infectees_susc_vacc)
        infector <- c()
        for (node in infectees_susc_vacc){
          # get list of neighbors for each infected
          neighbors <- neighbors(g,node)
          infector_node_vacc <- c()
          # figure out which neighbors are in vaccinated i_nodes and save those as possible infectors
          for (n in neighbors){
            if (n %in% i_nodes_vacc){
              infector_node_vacc <- c(infector_node_vacc,n)
              #print(n)
              #print(i_nodes[1,])
            }
          }
          if (length(infector_node_vacc)>0){
            # randomly pick from possible infectors
            infector <- infector_node_vacc[1]
          } else{
            # if no neighbors are in i_nodes then infector is unknown
            infector <- NA
          }
          #print(infector)
          # save to edge list]
          edge_list <- rbind(edge_list,cbind(rep(sim,length(infector)),infector,node))
        }
        
        if (length(v_nodes)>0) {
          # Same as above but with vaccinated susceptible nodes who are neighbors of vaccinated infected
          vacc_contacts_vacc<-lapply(potential_contacts_vacc,function(x,vacc) intersect(x,vacc),vacc=v_nodes)
          num_neighbours_vacc_vacc<-rapply(vacc_contacts_vacc,length)
          num_contacts_vacc_vacc<-rbinom(length(num_neighbours_vacc_vacc),num_neighbours_vacc_vacc,1-exp(-beta_vacc))
          infectees_vacc_vacc<-c(unlist(mapply(function(x,y) x[sample.int(length(x),y)],x=vacc_contacts_vacc,y=num_contacts_vacc_vacc)))
          infectees_vacc_vacc<-unique(infectees_vacc_vacc)
        } else {
          infectees_vacc_vacc<-c()
        }
        infector <- c()
        for (node in infectees_vacc_vacc){
          # get list of neighbors for each infected
          neighbors <- neighbors(g,node)
          infector_node_vacc <- c()
          # figure out which neighbors are in i_nodes and save those as possible infectors
          for (n in neighbors){
            if (n %in% i_nodes_vacc){
              infector_node_vacc <- c(infector_node_vacc,n)
              #print(n)
              #print(i_nodes[1,])
            }
          }
          if (length(infector_node_vacc)>0){
            # randomly pick from possible infectors
            infector <- infector_node_vacc[1]
          } else{
            # if no neighbors are in i_nodes then infector is unknown
            infector <- NA
          }
          #print(infector)
          # save to edge list]
          edge_list <- rbind(edge_list,cbind(rep(sim,length(infector)),infector,node))
        }
      } else {
        infectees_susc_vacc <- c()
        infectees_vacc_vacc <- c()
      } 
      
      # same for unvaccinated infected nodes
      if (length(i_nodes_unvacc)>0){
        # Get a list of all unvaccinated neighbours of all infected unvaccinated nodes
        potential_contacts_unvacc<-lapply(i_nodes_unvacc,function(x) neighbors(g,x))
        susc_contacts_unvacc<-lapply(potential_contacts_unvacc,function(x,susceptibles) intersect(x,susceptibles),susceptibles=s_nodes)
        num_neighbours_susc_unvacc<-rapply(susc_contacts_unvacc,length)
        # Sample from each group of neighbours in turn
        # First choose how many neighbours each node infects
        num_contacts_susc_unvacc<-rbinom(length(num_neighbours_susc_unvacc),num_neighbours_susc_unvacc,1-exp(-beta_unvacc))
        # Then sample from the neighbours
        # If one node gets picked twice by different nodes, just discard the duplicate.
        # In the very rare case that each i_nodes makes a number of new infectees equal to the number
        # of infectious nodes, mapply will make a matrix and unlist won't work. Therefore, the c() around
        # it ensures we turn the matrix into a vector. Unique then removes duplicates.
        infectees_susc_unvacc<-c(unlist(mapply(function(x,y) x[sample.int(length(x),y)],x=susc_contacts_unvacc,y=num_contacts_susc_unvacc)))
        infectees_susc_unvacc<-unique(infectees_susc_unvacc)
        infector <- c()
        for (node in infectees_susc_unvacc){
          # get list of neighbors for each infected
          neighbors <- neighbors(g,node)
          infector_node_unvacc <- c()
          # figure out which neighbors are in i_nodes and save those as possible infectors
          for (n in neighbors){
            if (n %in% i_nodes_unvacc){
              infector_node_unvacc <- c(infector_node_unvacc,n)
              #print(n)
              #print(i_nodes[1,])
            }
          }
          if (length(infector_node_unvacc)>0){
            # randomly pick from possible infectors
            infector <- infector_node_unvacc[1]
          } else{
            # if no neighbors are in i_nodes then infector is unknown
            infector <- NA
          }
          #print(infector)
          # save to edge list]
          edge_list <- rbind(edge_list,cbind(rep(sim,length(infector)),infector,node))
        }
        
        if (length(v_nodes)>0) {
          # Same as above but with vaccinated susceptible nodes who are neighbors of unvaccinated infected
          vacc_contacts_unvacc<-lapply(potential_contacts_unvacc,function(x,vacc) intersect(x,vacc),vacc=v_nodes)
          num_neighbours_vacc_unvacc<-rapply(vacc_contacts_unvacc,length)
          num_contacts_vacc_unvacc<-rbinom(length(num_neighbours_vacc_unvacc),num_neighbours_vacc_unvacc,1-exp(-beta_unvacc))
          infectees_vacc_unvacc<-c(unlist(mapply(function(x,y) x[sample.int(length(x),y)],x=vacc_contacts_unvacc,y=num_contacts_vacc_unvacc)))
          infectees_vacc_unvacc<-unique(infectees_vacc_unvacc)
        } else {
          infectees_vacc_unvacc<-c()
        }
        infector <- c()
        for (node in infectees_vacc_unvacc){
          # get list of neighbors for each infected
          neighbors <- neighbors(g,node)
          infector_node_unvacc <- c()
          # figure out which neighbors are in i_nodes and save those as possible infectors
          for (n in neighbors){
            if (n %in% i_nodes_unvacc){
              infector_node_unvacc <- c(infector_node_unvacc,n)
              #print(n)
              #print(i_nodes[1,])
            }
          }
          if (length(infector_node_unvacc)>0){
            # randomly pick from possible infectors
            infector <- infector_node_unvacc[1]
          } else{
            # if no neighbors are in i_nodes then infector is unknown
            infector <- NA
          }
          #print(infector)
          # save to edge list]
          edge_list <- rbind(edge_list,cbind(rep(sim,length(infector)),infector,node))
        }
      } else {
        infectees_susc_unvacc <- c()
        infectees_vacc_unvacc <- c()
      } 
      
      infectees_susc <- c(infectees_susc_unvacc,infectees_susc_vacc)
      infectees_vacc <- c(infectees_vacc_unvacc,infectees_vacc_vacc)
      
    } else {
      infectees_susc <- c()
      infectees_vacc <- c()
    } 
    
    infectees <- c(infectees_susc, infectees_vacc)
    
    # Pick out the nodes connected to the source that are still susceptible 
    # and haven't just been infected
    target_cnodes_susc <- setdiff(intersect(connected_nodes,s_nodes),infectees_susc)
    target_cnodes_vacc <- setdiff(intersect(connected_nodes,v_nodes),infectees_vacc)
    
    # Make a vector to represent external infection hazard for each individual
    communities_s <- V(g)[target_cnodes_susc]$community
    communities_v <- V(g)[target_cnodes_vacc]$community
    comm_sizes_s <- sapply(1:num_communities,function(x) sum(communities_s==x))
    comm_sizes_v <- sapply(1:num_communities,function(x) sum(communities_v==x))
    
    # Hazard of infection
    extFs_s<-rep(extF,comm_sizes_s)
    extFs_v<-rep(extF,comm_sizes_v)
    
    # Probability of infection
    prob_inf_fromsource <- 1 - exp(-mean(extFs_s)*source_num_inf)
    prob_inf_fromsource_v <- 1 - exp(-(1-VE)*mean(extFs_v)*source_num_inf)
    
    # Choose a number of individuals to be infected, then sample those individuals
    if (length(target_cnodes_susc)>0) {
      num_conn_inf_susc <- rbinom(1,length(target_cnodes_susc),prob_inf_fromsource)
      conn_inf_susc <- target_cnodes_susc[sample.int(length(target_cnodes_susc),num_conn_inf_susc,prob=extFs_s)]
    } else {
      conn_inf_susc <- c()
    }
    
    # Same as above, but for vaccinated individuals
    if (length(target_cnodes_vacc)>0) {
      num_conn_inf_vacc <- rbinom(1,length(target_cnodes_vacc),prob_inf_fromsource_v)
      conn_inf_vacc <- target_cnodes_vacc[sample.int(length(target_cnodes_vacc),num_conn_inf_vacc,prob=extFs_v)]
    } else {
      conn_inf_vacc <- c()
    }
    
    newinfected_susc <- c(infectees_susc,conn_inf_susc)
    newinfected_vacc <- c(infectees_vacc,conn_inf_vacc)
    newinfected <- c(newinfected_susc, newinfected_vacc)
    newinfected <- unique(newinfected)
    newimports <- c(conn_inf_susc,conn_inf_vacc)
    imports <- as.data.frame(cbind(rep(sim,length(newimports)),rep(NA,length(newimports)),newimports))
    names(imports) <- c("V1","infector","node")
    edge_list <- rbind(edge_list,imports)
    
    
    if (length(newinfected)>0) {
      
      # Add them to i_nodes and remove from s_nodes and v_nodes
      # Remove any progressing from e_nodes and add to i_nodes
      inf_periods <- rgamma(length(newinfected),infperiod_shape,infperiod_rate)
      i_nodes <- cbind(i_nodes,rbind(newinfected,rep(0,length(newinfected)),inf_periods))
      s_nodes<-setdiff(s_nodes,newinfected_susc)
      v_nodes <- setdiff(v_nodes,newinfected_vacc)
    }
    
    list(s_nodes, v_nodes, i_nodes,edge_list,results,g,sort(newinfected))
  }
  
  #### RUN THE EPIDEMIC IN THE SOURCE POPULATION ####
  # This is to define external infectious pressure to the network
  
  model <- function(t, y, parms) {
    with(as.list(c(y,parms)), {
      
      beta <- betahat * (1 - a2/(1 + exp(-a1 * (t - atau))))
      
      dS <- -beta * S * (I1+I2+I3) / (S+I1+I2+I3+R)
      dI1 <- beta * S * (I1+I2+I3) / (S+I1+I2+I3+R) - gamma * 3 * I1
      dI2 <- gamma * 3 * I1 - gamma * 3 * I2
      dI3 <- gamma * 3 * I2 - gamma * 3 * I3
      dR <- gamma * 3 * I3
      list(c(dS,dI1,dI2,dI3,dR))
    })
  }
  N <- 50000
  y<- c(S=N-1,I1=1,I2=0,I3=0,R=0)
  times<-seq(0,730,1)
  parms<-c(betahat=0.94,a1=0.19,a2=0.6,atau=27.79,gamma=0.33)
  out<-as.data.frame(lsoda(y,times,model,parms))
  
  #### RUN THE EPIDEMIC IN THE STUDY POPULATION ####
  
  # Define how the study population is linked to the source population
  # Connect all individuals to source population at same hazard
  # Constant of proportionality varies by community
  studypop_size<-length(V(g))
  connected_to_source <- V(g)$name
  
  # Calibrate extF to the number of introductions, given the progression of the epidemic in the source population
  num_communities <- max(V(g)$community)
  comm_sizes <- sapply(1:num_communities,function(x) length(V(g)[community==x]))
  sumsqrt <- sum(sqrt(comm_sizes))
  extF <- -log(1-num_introductions/(sqrt(comm_sizes)*sumsqrt))/trapz(times,out$I1+out$I2+out$I3)
  
  # Number of timesteps to run the epidemic - only need to go until the end of the trial
  num_timesteps <- trial_startday + trial_length + enrollment_period - 1
  
  if (bTrial) {
    
    # Parameters to do with trial recruitment
    # Enrollment per day is number of clusters enrolled per day
    enrollment_schedule <- rep(num_enrolled_per_day,enrollment_period)
    enroll_endday <- trial_startday+enrollment_period-1
    
    non_trial_clusters <- 1:max(V(g)$community)
    
  }
  
  # Initialize the S, E, I, and R nodes. I seed the epidemic from an SIR curve in a source population,
  # so initially all nodes in the study population are susceptible
  # i_nodes is a matrix. The first row is the  identity of the node. The second row
  # is the number of days since infection/infectiousness. The third row is the total nfectious period, 
  # drawn from a distribution when it becomes infected/infectious.
  # We are only going to consider a vaccine with effects on susceptibility, so only need one
  # vaccinated class
  i_nodes<-matrix(nrow=3,ncol=0)
  v_nodes<-c()
  s_nodes <- as.vector(V(g))
  r_nodes <- c()
  edge_list <- data.frame(NULL)
  
  # Initialize results.
  # Results will be, for each newly-infected node, the identity of the node, the day it was infected,
  # the community of which it is a member and its trial status at the time of infection. 
  # This is enough information to run a Cox PH with gamma frailty.
  # Make a data frame the size of the study pop and fill it in, then trim at the end
  results<-data.frame("SimulationNumber"=rep(NA,studypop_size),
                      "InfectedNode"=rep(NA,studypop_size),
                      "DayInfected"=rep(NA,studypop_size),
                      "Community"=rep(NA,studypop_size),
                      "TrialStatus"=rep(NA,studypop_size),
                      "DayEnrolled"=rep(NA,studypop_size),
                      "DayRecovered"=rep(NA,studypop_size))
  numinfected<-0
  
  for (t in 1:num_timesteps) {
    
    # I'm recovering first, so I need to ensure that everyone has at least one chance to infect.
    # I do this by initializing an infectious node with 0 days since infection, seeing whether they
    # recover, then advancing them one day along their infectious period.
    
    if (bTrial) {
      
      # Recruit and randomize if during the enrollment period
      if ((t>=trial_startday) && (t<=enroll_endday)) {
        
        num_to_enroll <- enrollment_schedule[t-trial_startday+1]
        
        if (bCluster == 0) {
          # Individually-randomized trial, stratifying on community
          
          # Need to choose from the clusters not already enrolled
          new_clusters <- sample(non_trial_clusters,num_to_enroll)
          
          # From the chosen clusters, choose a fraction of the non-infectious individual. That fraction is defined in the inputs
          # I will then vaccinate half of each chosen sample
          # For each new cluster, I sample from that cluster a proportion of the whole cluster,
          # but only from the susceptible individuals. If I'm trying to sample more than are available (because
          # there are lots of infectious/recovered individuals), just sample all of them.
          # These are the people who are recruited - I then assign half and half to vaccine or control
          new_recruits <- lapply(new_clusters,
                                 function(x) sample(intersect(V(g)[community==x]$name,c(s_nodes)),
                                                    min(round(cluster_coverage*length(V(g)[community==x]$name)),
                                                        length(intersect(V(g)[community==x]$name,c(s_nodes))))))
          new_vacc <- unlist(lapply(new_recruits,
                                    function(x) sample(x,round(length(x)/2))))
          new_controls <- setdiff(unlist(new_recruits),new_vacc)
          
          
          non_trial_clusters<-setdiff(non_trial_clusters,new_clusters)
          
        } else {
          
          # We try and enroll as many from the cluster as you can. I have set an
          # enrollment rate rather than cluster size, e.g. 70% enrolled in each cluster.
          # It means that every simulated trial would have slightly different numbers enrolled 
          # (=coverage*ave_community_size)
          
          # Need to choose from the clusters not already enrolled
          new_clusters <- sample(non_trial_clusters,num_to_enroll)
          new_clusters_v <- sample(new_clusters,num_to_enroll/2)
          new_clusters_c <- setdiff(new_clusters,new_clusters_v)
          # From the chosen clusters, a fraction of the non-infectious individual. That fraction is defined in the inputs
          # This looks complicated: for each new cluster, I sample from that cluster a proportion of the whole cluster,
          # but only from the susceptible  individuals. If I'm trying to sample more than are available (because
          # there are lots of infectious/recovered individuals), just sample all of them.
          new_vacc <- unlist(lapply(new_clusters_v,
                                    function(x) sample(intersect(V(g)[community==x]$name,c(s_nodes)),
                                                       min(round(cluster_coverage*length(V(g)[community==x]$name)),
                                                           length(intersect(V(g)[community==x]$name,c(s_nodes)))))))
          new_controls <- unlist(lapply(new_clusters_c,
                                        function(x) sample(intersect(V(g)[community==x]$name,c(s_nodes)),
                                                           min(round(cluster_coverage*length(V(g)[community==x]$name)),
                                                               length(intersect(V(g)[community==x]$name,c(s_nodes)))))))
          
          non_trial_clusters<-setdiff(non_trial_clusters,new_clusters)
        }
        
        V(g)[name %in% new_controls]$trialstatus<-0
        V(g)[name %in% new_vacc]$trialstatus<-1
        V(g)[name %in% new_controls]$enrollmentday<-t
        V(g)[name %in% new_vacc]$enrollmentday<-t
        
        # Move the vaccinated susceptibles to from s_nodes to v_nodes
        vacc_susc <- intersect(s_nodes,new_vacc)
        s_nodes <- setdiff(s_nodes,vacc_susc)
        v_nodes <- c(v_nodes,vacc_susc)
        
        # all or nothing vaccine
        num_immune <- rbinom(1,length(vacc_susc),VE)
        new_immune <- sample(vacc_susc,num_immune)
        v_nodes <- setdiff(v_nodes,new_immune)
        r_nodes <- c(r_nodes,new_immune)
        
      }
      
    }
    
    # Only need to recover if there are any infected
    if ((ncol(i_nodes)>0)) {
      list[i_nodes,r_nodes,results]<-
        recover(i_nodes,r_nodes,infperiod_shape,infperiod_rate,t,results)
    }
    list[s_nodes,v_nodes,i_nodes,edge_list,results,g,newinfected]<-
      spread(g,s_nodes,v_nodes,i_nodes,edge_list,t,results,
             beta,VE,VEI,
             connected_to_source,extF,out$I1[t]+out$I2[t]+out$I3[t])
    
    numnewinfected<-length(newinfected)
    if (numnewinfected>0) {
      
      newcommunities <- V(g)[name %in% newinfected]$community
      
      # Update results
      results$SimulationNumber[(numinfected+1):(numinfected+numnewinfected)]<-rep(simnum,numnewinfected)
      results$InfectedNode[(numinfected+1):(numinfected+numnewinfected)]<-newinfected
      results$DayInfected[(numinfected+1):(numinfected+numnewinfected)]<-rep(t,numnewinfected)
      results$Community[(numinfected+1):(numinfected+numnewinfected)]<-newcommunities
      results$TrialStatus[(numinfected+1):(numinfected+numnewinfected)]<-V(g)[name %in% newinfected]$trialstatus
      results$DayEnrolled[(numinfected+1):(numinfected+numnewinfected)]<-V(g)[name %in% newinfected]$enrollmentday
      
      numinfected <- numinfected+numnewinfected
      
      V(g)[name %in% newinfected]$eventstatus<-1
      
    }
    
  }
  
  trial_nodes <- V(g)[!is.na(V(g)$trialstatus)]$name
  trial_nodes_info<-data.frame("SimulationNumber"=rep(simnum,length(trial_nodes)),
                               "Node"=trial_nodes,
                               "Community"=V(g)[trial_nodes]$community,
                               "TrialStatus"=V(g)[trial_nodes]$trialstatus,
                               "DayEnrolled"=V(g)[trial_nodes]$enrollmentday)
  
  # Tidy up results
  if (numinfected>0) {
    results<-results[1:numinfected,]
  } else {
    results<-results[1,]
    results$SimulationNumber[1]<-simnum
  }
  
  # summarize edge list
  as.data.frame(edge_list)
  edge_list$infector_stat <- V(g)$trialstatus[edge_list$infector]
  edge_list$node_stat <- V(g)$trialstatus[edge_list$node]
  
  # merge edge list and results
  names(edge_list) <- c("SimNum","Infector","InfectedNode","infector_stat","node_stat")
  results_edge_list <- merge(results,edge_list)
  
  
  list(results,results_edge_list,trial_nodes_info,edge_list,g,newinfected)
  
}

analyse_data <- function(results_edge_list,trial_nodes,trial_startday,trial_length,
                         numclusters_perarm,bCluster) {  
  #results_edge_list$DayInfected <- results_edge_list$DayInfected - results_edge_list$DayEnrolled
  # Get a list of nodes that were enrolled in the trial but never infected
  noninf<-setdiff(trial_nodes$Node,results_edge_list$InfectedNode)
  # Get list of nodes that became infectious while they were in the trial
  # This is the step that excludes those who were R at time of enrollment
  
  # for now not excluding those not in trial from genomics analysis
  results_edge_list_analysis<-results_edge_list
  #results_edge_list_analysis<-results_edge_list[!is.na(results_edge_list$TrialStatus),]
  
  # Get a list of nodes who were infected after their follow-up time was over
  # (i.e. those enrolled at the beginning but infected right at the end)
  results_edge_list_analysis$eventstatus<-rep(1,nrow(results_edge_list_analysis))
  censored <- results_edge_list_analysis[results_edge_list_analysis$DayInfected>trial_length,]
  results_edge_list_analysis<-results_edge_list_analysis[results_edge_list_analysis$DayInfected<=trial_length,]
  # Assign them eventstatus=1 for the Cox analysis
  # Make data frame for those who were never infected (i.e. censored by end of study)
  noninfdf<-data.frame(InfectedNode=noninf,SimulationNumber=rep(sim,length(noninf)),
                       DayInfected=rep(trial_length,length(noninf)),
                       Community=trial_nodes$Community[trial_nodes$Node %in% noninf],
                       TrialStatus=trial_nodes$TrialStatus[trial_nodes$Node %in% noninf],
                       DayEnrolled=trial_nodes$DayEnrolled[trial_nodes$Node %in% noninf],
                       DayRecovered=rep(0,length(noninf)),
                       SimNum=rep(sim,length(noninf)),Infector=rep(NA,length(noninf)),
                       infector_stat=rep(NA,length(noninf)),
                       node_state=trial_nodes$TrialStatus[trial_nodes$Node %in% noninf],eventstatus=rep(0,length(noninf)))
  
  if (nrow(censored)>0) {
    censored$DayInfected<-trial_length
    censored$eventstatus<-0
  }
  
  names(results_edge_list_analysis)<-c("InfectedNode","SimulationNumber","DayInfected","Community","TrialStatus","DayEnrolled","DayRecovered","SimNum",
                                       "Infector","infector_stat","node_stat","eventstatus")
  names(censored)<-c("InfectedNode","SimulationNumber","DayInfected","Community","TrialStatus","DayEnrolled","DayRecovered","SimNum",
                     "Infector","infector_stat","node_stat","eventstatus")
  names(noninfdf)<-c("InfectedNode","SimulationNumber","DayInfected","Community","TrialStatus","DayEnrolled","DayRecovered","SimNum",
                     "Infector","infector_stat","node_stat","eventstatus")
  results_edge_list_analysis<-rbind(results_edge_list_analysis,noninfdf,censored)

  # Finally, exclude any cases who were infected during the first n days of follow-up
  # This tries to rid of those who were already latently infected when enrolled
  
  numevents_vacc <- nrow(results_edge_list_analysis[(results_edge_list_analysis$eventstatus==1) & (results_edge_list_analysis$TrialStatus==1),])
  numevents_cont <- nrow(results_edge_list_analysis[(results_edge_list_analysis$eventstatus==1) & (results_edge_list_analysis$TrialStatus==0),])
  num_vacc <- nrow(results_edge_list_analysis[(results_edge_list_analysis$TrialStatus==1),])
  num_cont <- nrow(results_edge_list_analysis[(results_edge_list_analysis$TrialStatus==0),])  
  
  total_vacc_pt <- sum(results_edge_list_analysis$DayInfected[results_edge_list_analysis$TrialStatus==1])
  total_cont_pt <- sum(results_edge_list_analysis$DayInfected[results_edge_list_analysis$TrialStatus==0])
  VE_pointest <- 1 - (numevents_vacc/num_vacc)/(numevents_cont/num_cont)
  pval <- NA
  
  sample_size <- nrow(results_edge_list_analysis)
  
  # get list of infectors and infectees status
  summary_infectors <- aggregate(results_edge_list_analysis$InfectedNode,by=c(list(results_edge_list_analysis$infector_stat),list(results_edge_list_analysis$node_stat)),length)
  names(summary_infectors)<-c("infector_stat","node_state","number")
  vacc_vacc_inf <- summary_infectors[which(summary_infectors[,1]==1 & summary_infectors[,2]==1),3]
  vacc_cont_inf <- summary_infectors[which(summary_infectors[,1]==1 & summary_infectors[,2]==0),3]
  cont_vacc_inf <- summary_infectors[which(summary_infectors[,1]==0 & summary_infectors[,2]==1),3]
  cont_cont_inf <- summary_infectors[which(summary_infectors[,1]==0 & summary_infectors[,2]==0),3]
  
  importation <- nrow(results_edge_list_analysis[is.na(results_edge_list_analysis$infector_stat) & results_edge_list_analysis$eventstatus==1,])
  
  RR_inf <- (vacc_vacc_inf+vacc_cont_inf)/(cont_vacc_inf+cont_cont_inf)
  num_inf_by_vacc <- vacc_vacc_inf+vacc_cont_inf
  num_inf_by_cont <- cont_vacc_inf+cont_cont_inf

  results_edge_list_analysis$Infector[is.na(results_edge_list_analysis$Infector) & results_edge_list_analysis$eventstatus==1] <- 0
  results_infected <- results_edge_list_analysis[results_edge_list_analysis$eventstatus==1,]
  
  # for those who remain infected after end of the trial, make recovered day after trial ends
  results_infected$DayRecovered[is.na(results_infected$DayRecovered)] <- trial_startday + trial_length
  
  list(VE_pointest,pval,numevents_vacc,numevents_cont,sample_size,
       RR_inf,vacc_vacc_inf,vacc_cont_inf,cont_vacc_inf,cont_cont_inf,importation,results_edge_list_analysis,
       results_infected)
  
}

# end functions and beginning of runs #
RRs_inf <- rep(NA,nsim)
pvals<-rep(NA,nsim)
VEs<-rep(NA,nsim)
numevents_cont<-rep(NA,nsim)
numevents_vacc<-rep(NA,nsim)
ss <- rep(NA,nsim)
vacc_vacc_infs <-rep(NA,nsim)
vacc_cont_infs <-rep(NA,nsim)
cont_vacc_infs <-rep(NA,nsim)
cont_cont_infs <-rep(NA,nsim)
importations <-rep(NA,nsim)

sim <- 1

for (sim in 1:nsim) {

    g<-make_network(ave_community_size, community_size_range, num_communities,rate_within, rate_between)
    
    list[results,results_edge_list,trial_nodes,edge_list,g]<-
      network_epidemic(g,beta,num_introductions,direct_VE,infect_VE,edge_list,
                       infperiod_shape,infperiod_rate,1,0,
                       trial_startday,trial_length,num_clusters_enrolled_per_day,
                       enrollment_period,cluster_coverage,sim)
    
    list[VE,pval,events_vacc,events_cont,analysed_trialsize,
         RR_inf,vacc_vacc_inf,vacc_cont_inf,cont_vacc_inf,cont_cont_inf,importation,results_edge_list_analysis,
         results_infected]<-
      analyse_data(results_edge_list,trial_nodes,trial_startday,trial_length,num_clusters_perarm,0)
    
write.csv(results_edge_list_analysis, paste0(sim,"_results_edge_list_analysis",R0,"_",infect_VE,"_",direct_VE,"_",mut.rate,"_",bn,"_",infperiod_shape,"_",sample_percent,"_","sample.csv"))
save(g,file=paste0(sim,"_g",R0,"_",infect_VE,"_",direct_VE,"_",mut.rate,"_",bn,"_",infperiod_shape,"_",sample_percent,"_","sample.Rdata"))
}



