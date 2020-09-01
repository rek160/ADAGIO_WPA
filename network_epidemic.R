# run epidemic and trial
# network code written by Matt Hitchings and adapted by Rebecca Kahn

network_epidemic<-function(g,beta,num_introductions,VE,VEI,edge_list,
                           incperiod_shape,incperiod_rate,infperiod_shape,infperiod_rate,
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
  recover<-function(e_nodes,i_nodes,r_nodes,infperiod_shape,infperiod_rate,t,results) {
    # Input is a list of the exposed nodes, 
    # with number of days since infection and total incubation/latent
    # period, and equivalently for the infectious nodes.
    # For each of these nodes, we will add it to newinfectious if the number of days exposed has
    # reached the total length of the incubation period, and equivalently for the infectious nodes.
    
    # Advance infectious nodes first, otherwise you will doubly advance any nodes switching from
    # exposed to infectious at this time step
    indices_to_remove <- i_nodes[2,]>=i_nodes[3,]
    newremoved<-as.vector(i_nodes[1,indices_to_remove])
    
    # Add one day to the length of each infected individual's time infected
    i_nodes[2,] <- i_nodes[2,]+1
    
    # Remove any recovered from i_nodes and add to r_nodes
    i_nodes <- i_nodes[,!(i_nodes[1,] %in% newremoved),drop=FALSE]
    r_nodes <- c(r_nodes,newremoved)
    results$DayRecovered[results$InfectedNode %in% newremoved] <-t
    
    # Now advance exposed nodes
    indices_to_remove <- e_nodes[2,]>=e_nodes[3,]
    newinfectious<-as.vector(e_nodes[1,indices_to_remove])
    
    # Add one day to the length of each infected individual's time infected
    e_nodes[2,] <- e_nodes[2,]+1
    
    # Remove any progressing from e_nodes and add to i_nodes
    e_nodes <- e_nodes[,!(e_nodes[1,] %in% newinfectious),drop=FALSE]
    inf_periods <- rgamma(length(newinfectious),infperiod_shape,infperiod_rate)
    i_nodes <- cbind(i_nodes,rbind(newinfectious,rep(0,length(newinfectious)),inf_periods))
    list(results,e_nodes, i_nodes, r_nodes, sort(newinfectious))
  }
  
  spread<-function(g, s_nodes, v_nodes, e_nodes, i_nodes, edge_list, 
                   beta, VE, VEI,
                   incperiod_shape, incperiod_rate,
                   connected_nodes,external_inf_F,source_num_inf){
    # Spread will create new infected nodes from two sources: infectious nodes within the the study
    # population, and external pressure from the source population
    # Inputs:
    # g is the graph, used to find neighbours of infected nodes
    # s_nodes, e_nodes and i_nodes are susceptible, exposed and infected nodes
    # beta is the hazard of infection for one contact
    # incperiod_shape and rate are used to assign each newly exposed node a latent/incubation period
    # length, currently drawn from a gamma distribution
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
      
      
      # Make a beta vector
      # first is vaccination status of infector; second is vaccination status of infectee
      beta_unvacc_unvacc <- rep(beta,ncol(i_nodes))
      beta_unvacc_vacc <- beta_unvacc_unvacc*(1-VE)
      beta_vacc_unvacc <- beta_unvacc_unvacc*(1-VEI)
      beta_vacc_vacc <- beta_unvacc_unvacc*(1-VE)*(1-VEI)
      
      vacc <- V(g)[trialstatus==1 & !is.na(trialstatus)]$name
      
      i_nodes_unvacc <- setdiff(i_nodes[1,],vacc)
      i_nodes_vacc <- intersect(i_nodes[1,],vacc)
      #cat(t,length(i_nodes_unvacc),length(i_nodes_vacc),"\n")
      
      # start with unvaccinated infected nodes
      if (length(i_nodes_unvacc)>0){
        # Get a list of all unvaccinated neighbours of all infected unvaccinated nodes
        potential_contacts_unvacc<-lapply(i_nodes_unvacc,function(x) neighbors(g,x))
        susc_contacts_unvacc<-lapply(potential_contacts_unvacc,function(x,susceptibles) intersect(x,susceptibles),susceptibles=s_nodes)
        num_neighbours_susc_unvacc<-rapply(susc_contacts_unvacc,length)
        # Sample from each group of neighbours in turn
        # First choose how many neighbours each node infects
        num_contacts_susc_unvacc<-rbinom(length(num_neighbours_susc_unvacc),num_neighbours_susc_unvacc,1-exp(-beta_unvacc_unvacc))
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
          num_contacts_vacc_unvacc<-rbinom(length(num_neighbours_vacc_unvacc),num_neighbours_vacc_unvacc,1-exp(-beta_unvacc_vacc))
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
      
      # remove those infected from susceptible and v_nodes so not doubly infected next (rare occurrence)
      s_nodes2 <- setdiff(s_nodes,infectees_susc_unvacc)
      v_nodes2 <- setdiff(v_nodes,infectees_vacc_unvacc)
      
      # same for vaccinated infected nodes
      if (length(i_nodes_vacc)>0){
        # Get a list of all unvaccinated neighbours of all infected vaccinated nodes
        potential_contacts_vacc<-lapply(i_nodes_vacc,function(x) neighbors(g,x))
        susc_contacts_vacc<-lapply(potential_contacts_vacc,function(x,susceptibles) intersect(x,susceptibles),susceptibles=s_nodes2)
        num_neighbours_susc_vacc<-rapply(susc_contacts_vacc,length)
        # Sample from each group of neighbours in turn
        # First choose how many neighbours each node infects
        num_contacts_susc_vacc<-rbinom(length(num_neighbours_susc_vacc),num_neighbours_susc_vacc,1-exp(-beta_vacc_unvacc))
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
        
        if (length(v_nodes)>0) {
          # Same as above but with vaccinated susceptible nodes who are neighbors of vaccinated infected
          vacc_contacts_vacc<-lapply(potential_contacts_vacc,function(x,vacc) intersect(x,vacc),vacc=v_nodes2)
          num_neighbours_vacc_vacc<-rapply(vacc_contacts_vacc,length)
          num_contacts_vacc_vacc<-rbinom(length(num_neighbours_vacc_vacc),num_neighbours_vacc_vacc,1-exp(-beta_vacc_vacc))
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
      
      # Give each newly exposed node an incubation/latent period
      inc_periods <- rgamma(length(newinfected),incperiod_shape,incperiod_rate)
      # Add them to e_nodes and remove from s_nodes and v_nodes
      e_nodes <- cbind(e_nodes,rbind(newinfected,rep(0,length(newinfected)),inc_periods))
      e_nodes_master <- cbind(e_nodes_master,e_nodes)
      s_nodes<-setdiff(s_nodes,newinfected_susc)
      v_nodes <- setdiff(v_nodes,newinfected_vacc)
      
    }
    
    list(s_nodes, v_nodes, e_nodes, edge_list,e_nodes_master)
  }
  
  #### RUN THE EPIDEMIC IN THE SOURCE POPULATION ####
  # This is to define external infectious pressure to the network
  
  model <- function(t, y, parms) {
    with(as.list(c(y,parms)), {
      
      beta <- betahat * (1 - a2/(1 + exp(-a1 * (t - atau))))
      
      dS <- -beta * S * (I1+I2+I3) / (S+E1+E2+E3+I1+I2+I3+R)
      dE1 <- beta * S * (I1+I2+I3) / (S+E1+E2+E3+I1+I2+I3+R) - sigma * 3 * E1
      dE2 <- sigma * 3 * E1 - sigma * 3 * E2
      dE3 <- sigma * 3 * E2 - sigma * 3 * E3
      dI1 <- sigma * 3 * E3 - gamma * 3 * I1
      dI2 <- gamma * 3 * I1 - gamma * 3 * I2
      dI3 <- gamma * 3 * I2 - gamma * 3 * I3
      dR <- gamma * 3 * I3
      list(c(dS,dE1,dE2,dE3,dI1,dI2,dI3,dR))
    })
  }
  N <- 50000
  y<- c(S=N-1,E1=0,E2=0,E3=0,I1=1,I2=0,I3=0,R=0)
  times<-seq(0,300,1)
  parms<-c(betahat=0.94,a1=0.19,a2=0.6,atau=27.79,sigma=0.14,gamma=0.33)
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
  # i_nodes and e_nodes are matrices. The first row is the  identity of the node. The second row
  # is the number of days since infection/infectiousness. The third row is the total incubation/infectious period, 
  # drawn from a distribution when it becomes infected/infectious.
  # We are only going to consider a vaccine with effects on susceptibility, so only need one
  # vaccinated class
  e_nodes<-matrix(nrow=3,ncol=0)
  e_nodes_master <-matrix(nrow=3,ncol=0)
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
                      "DayExposed"=rep(NA,studypop_size),
                      "DayRecovered"=rep(NA,studypop_size))
  numinfectious<-0
  
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
          # but only from the susceptible or exposed individuals. If I'm trying to sample more than are available (because
          # there are lots of infectious/recovered individuals), just sample all of them.
          # These are the people who are recruited - I then assign half and half to vaccine or control
          new_recruits <- lapply(new_clusters,
                                 function(x) sample(intersect(V(g)[community==x]$name,c(e_nodes[1,],s_nodes)),
                                                    min(round(cluster_coverage*length(V(g)[community==x]$name)),
                                                        length(intersect(V(g)[community==x]$name,c(e_nodes[1,],s_nodes))))))
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
          # but only from the susceptible or exposed individuals. If I'm trying to sample more than are available (because
          # there are lots of infectious/recovered individuals), just sample all of them.
          new_vacc <- unlist(lapply(new_clusters_v,
                                    function(x) sample(intersect(V(g)[community==x]$name,c(e_nodes[1,],s_nodes)),
                                                       min(round(cluster_coverage*length(V(g)[community==x]$name)),
                                                           length(intersect(V(g)[community==x]$name,c(e_nodes[1,],s_nodes)))))))
          new_controls <- unlist(lapply(new_clusters_c,
                                        function(x) sample(intersect(V(g)[community==x]$name,c(e_nodes[1,],s_nodes)),
                                                           min(round(cluster_coverage*length(V(g)[community==x]$name)),
                                                               length(intersect(V(g)[community==x]$name,c(e_nodes[1,],s_nodes)))))))
          
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
        
      }
      
    }
    
    # Only need to recover if there are any infected or exposed
    if ((ncol(i_nodes)>0)||(ncol(e_nodes)>0)) {
      list[results,e_nodes,i_nodes,r_nodes,newinfectious]<-
        recover(e_nodes,i_nodes,r_nodes,infperiod_shape,infperiod_rate,t,results)
      
    } else {
      newinfectious <- c()
    }
    
    list[s_nodes,v_nodes,e_nodes,edge_list,e_nodes_master]<-
      spread(g,s_nodes,v_nodes,e_nodes,i_nodes,edge_list,
             beta,VE,VEI,incperiod_shape,incperiod_rate,
             connected_to_source,extF,out$I1[t]+out$I2[t]+out$I3[t])
    
    
    numnewinfectious<-length(newinfectious)
    
    if (numnewinfectious>0) {
      
      newcommunities <- V(g)[name %in% newinfectious]$community
      
      # Update results
      results$SimulationNumber[(numinfectious+1):(numinfectious+numnewinfectious)]<-rep(simnum,numnewinfectious)
      results$InfectedNode[(numinfectious+1):(numinfectious+numnewinfectious)]<-newinfectious
      results$DayInfected[(numinfectious+1):(numinfectious+numnewinfectious)]<-rep(t,numnewinfectious)
      results$Community[(numinfectious+1):(numinfectious+numnewinfectious)]<-newcommunities
      results$TrialStatus[(numinfectious+1):(numinfectious+numnewinfectious)]<-V(g)[name %in% newinfectious]$trialstatus
      results$DayEnrolled[(numinfectious+1):(numinfectious+numnewinfectious)]<-V(g)[name %in% newinfectious]$enrollmentday
      
      numinfectious <- numinfectious+numnewinfectious
      
    }
    
  }
  
  trial_nodes <- V(g)[!is.na(V(g)$trialstatus)]$name
  trial_nodes_info<-data.frame("SimulationNumber"=rep(simnum,length(trial_nodes)),
                               "Node"=trial_nodes,
                               "Community"=V(g)[trial_nodes]$community,
                               "TrialStatus"=V(g)[trial_nodes]$trialstatus,
                               "DayEnrolled"=V(g)[trial_nodes]$enrollmentday)
  
  # Tidy up results
  if (numinfectious>0) {
    results<-results[1:numinfectious,]
  } else {
    results<-results[1,]
    results$SimulationNumber[1]<-simnum
  }
  
  #add exposed date
  for (i in 1:nrow(results)){
    results[i,7] <- results[i,3]-ceiling(e_nodes_master[3,which(e_nodes_master[1,]==results[i,2])[1]])
    print(i)
  }
  
  
  # summarize edge list
  as.data.frame(edge_list)
  edge_list$infector_stat <- V(g)$trialstatus[edge_list$infector]
  edge_list$node_stat <- V(g)$trialstatus[edge_list$node]
  
  # merge edge list and results
  names(edge_list) <- c("SimNum","Infector","InfectedNode","infector_state","node_stat")
  results_edge_list <- merge(results,edge_list)
  
  
  list(results,results_edge_list,trial_nodes_info,edge_list,e_nodes_master,g)
  
}