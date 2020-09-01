# written by Matt Hitchings
analyse_data <- function(results_edge_list,trial_nodes,trial_startday,trial_length,ave_inc_period,
                         numclusters_perarm,bCluster) {
  
  library(survival)
  library(coxme)
  library(frailtypack)
  #library(lme4)
  
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
  
  coxmodel <- function(data,VEpointest) {
    
    survmodel<-try(coxph(Surv(DayInfected,eventstatus)~TrialStatus+strata(Community),data),silent=T)
    usesurvmod <- !inherits(survmodel, 'try-error')
    
    if (usesurvmod && vcov(survmodel)>=0){
      # If no error was thrown and the variance is positive, use the results of the model
      
      vaccEffEst <- 1-exp(survmodel$coefficient + c(0, 1.96, -1.96)*sqrt(survmodel$var))
      zval <- survmodel$coefficient/sqrt(survmodel$var)
      pval <- pnorm(zval, lower.tail = vaccEffEst[1]>0)*2
      
    } else {
      
      vaccEffEst<-c(VEpointest,NA,NA)
      pval <- NA
      
    }
    
    list(vaccEffEst,pval)
    
  }
  
  #results_edge_list$DayInfected <- results_edge_list$DayInfected - results_edge_list$DayEnrolled
  
  
  # Get a list of nodes that were enrolled in the trial but never infected
  noninf<-setdiff(trial_nodes$Node,results_edge_list$InfectedNode)
  # Get list of nodes that became infectious while they were in the trial
  # This is the step that excludes those who were R at time of enrollment
  results_edge_list_analysis<-results_edge_list[!is.na(results_edge_list$TrialStatus),]
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
                       DayExposed=rep(NA,length(noninf)),
                       DayRecovered=rep(NA,length(noninf)),
                       SimNum=rep(sim,length(noninf)),Infector=rep(NA,length(noninf)),
                       infector_state=rep(NA,length(noninf)),
                       node_state=trial_nodes$TrialStatus[trial_nodes$Node %in% noninf],eventstatus=rep(0,length(noninf)))
  if (nrow(censored)>0) {
    censored$DayInfected<-trial_length
    censored$eventstatus<-0
  }
  
  names(results_edge_list_analysis)<-c("InfectedNode","SimulationNumber","DayInfected","Community","TrialStatus","DayEnrolled","DayExposed","DayRecovered","SimNum",
                                       "Infector","infector_state","node_stat","eventstatus")
  names(censored)<-c("InfectedNode","SimulationNumber","DayInfected","Community","TrialStatus","DayEnrolled","DayExposed","DayRecovered","SimNum",
                     "Infector","infector_state","node_stat","eventstatus")
  names(noninfdf)<-c("InfectedNode","SimulationNumber","DayInfected","Community","TrialStatus","DayEnrolled","DayExposed","DayRecovered","SimNum",
                     "Infector","infector_state","node_stat","eventstatus")
  results_edge_list_analysis<-rbind(results_edge_list_analysis,noninfdf,censored)
  #print(results_edge_list_analysis)
  
  # Finally, exclude any cases who were infected during the first n days of follow-up
  # This tries to rid of those who were already latently infected when enrolled
  results_edge_list_analysis<-results_edge_list_analysis[(results_edge_list_analysis$DayInfected-results_edge_list_analysis$DayEnrolled)>ave_inc_period,]
  
  numevents_vacc <- nrow(results_edge_list_analysis[(results_edge_list_analysis$eventstatus==1) & (results_edge_list_analysis$TrialStatus==1),])
  numevents_cont <- nrow(results_edge_list_analysis[(results_edge_list_analysis$eventstatus==1) & (results_edge_list_analysis$TrialStatus==0),])
  
  total_vacc_pt <- sum(results_edge_list_analysis$DayInfected[results_edge_list_analysis$TrialStatus==1])
  total_cont_pt <- sum(results_edge_list_analysis$DayInfected[results_edge_list_analysis$TrialStatus==0])
  VE_pointest <- 1 - (numevents_vacc/total_vacc_pt)/(numevents_cont/total_cont_pt)
  
  sample_size <- nrow(results_edge_list_analysis)
  
  results_edge_list_analysis$Infector[is.na(results_edge_list_analysis$Infector) & results_edge_list_analysis$eventstatus==1] <- 0
  results_infected <- results_edge_list_analysis[results_edge_list_analysis$eventstatus==1,]
  
  # for those who remain infected after end of the trial, make recovered day after trial ends
  results_infected$DayRecovered[is.na(results_infected$DayRecovered)] <- trial_startday + trial_length + 1
  
  # for infectors not in the trial, make infector external source
  results_infected$Infector[!(results_infected$Infector %in% results_infected$InfectedNode)] <- 0
  
  
  # Analysis for iRCT
  if ((numevents_vacc>0)&&(numevents_cont>0)) {
    # If we have events in both arms, can try a Cox PH. It can still be singular, so if it
    # throws an error, the trial has failed (not enough events)
    list[vaccEffEst,pval] <- coxmodel(results_edge_list_analysis,VE_pointest)
    
  } else if ((numevents_vacc>0)&&(numevents_cont==0)) {
    # If there are no events in the control arm but
    # events in the vaccine arm, VE estimate is -1 and p-value is 1
    vaccEffEst<-c(-1,-1,-1)
    pval <- 1
    
  } else if ((numevents_vacc==0)&&(numevents_cont>0)) {
    # If there are no events in the vaccine arm and events in the control arm, VE is 1, but
    # for the p-value we add one event to both arms and do a Cox regression on that data
    # I give both events the median time among control events.
    newevent_v_rownum <- min(which((results_edge_list_analysis$eventstatus==0)&(results_edge_list_analysis$TrialStatus==1)))
    newevent_c_rownum <- min(which((results_edge_list_analysis$eventstatus==0)&(results_edge_list_analysis$TrialStatus==0)))
    
    eventtime <- median(results_edge_list_analysis$DayInfected[results_edge_list_analysis$eventstatus==1])
    
    results_edge_list_analysis$DayInfected[newevent_v_rownum] <- eventtime
    results_edge_list_analysis$eventstatus[newevent_v_rownum] <- 1
    
    results_edge_list_analysis$DayInfected[newevent_c_rownum] <- eventtime
    results_edge_list_analysis$eventstatus[newevent_c_rownum] <- 1
    
    list[,pval] <- coxmodel(results_edge_list_analysis,VE_pointest)
    vaccEffEst<-1
    
  } else {
    # If no events are observed in either arm, the trial has failed and no result can be obtained
    vaccEffEst<-c(NA,NA,NA)
    pval <- NA
    
  }
  
  # create data frames of 4 combinations of infector/infectee vaccination status; 
  # first is infector, second is infectee; ex.: vacc_cont; vaccinated infecting a control
  results_edge_list_analysis_vacc_vacc <- results_edge_list_analysis[(results_edge_list_analysis$infector_stat=="1") & !is.na(results_edge_list_analysis$infector_stat) & (results_edge_list_analysis$node_stat=="1") & !is.na(results_edge_list_analysis$node_stat) ,]
  results_edge_list_analysis_vacc_cont <- results_edge_list_analysis[(results_edge_list_analysis$infector_stat=="1") & !is.na(results_edge_list_analysis$infector_stat) & (results_edge_list_analysis$node_stat=="0") & !is.na(results_edge_list_analysis$node_stat) ,]
  results_edge_list_analysis_cont_vacc <- results_edge_list_analysis[(results_edge_list_analysis$infector_stat=="0") & !is.na(results_edge_list_analysis$infector_stat) & (results_edge_list_analysis$node_stat=="1") & !is.na(results_edge_list_analysis$node_stat) ,]
  results_edge_list_analysis_cont_cont <- results_edge_list_analysis[(results_edge_list_analysis$infector_stat=="0") & !is.na(results_edge_list_analysis$infector_stat) & (results_edge_list_analysis$node_stat=="0") & !is.na(results_edge_list_analysis$node_stat) ,]
  
  # get summary of vaccinated infecting vaccinated
  if (nrow(results_edge_list_analysis_vacc_vacc)>0){
    infectors_summary_vacc_vacc <- aggregate(x = results_edge_list_analysis_vacc_vacc, 
                                             by = list(unique.values = results_edge_list_analysis_vacc_vacc$Infector), 
                                             FUN = length)
    infectors_summary_vacc_vacc <- infectors_summary_vacc_vacc[,1:2]
  } else { 
    infectors_summary_vacc_vacc <- as.data.frame(cbind(NA,NA))
  }
  names(infectors_summary_vacc_vacc) <- c("infector","# infectees")
  #print(infectors_summary_vacc_vacc)
  
  # get summary of vaccinated infecting control
  if (nrow(results_edge_list_analysis_vacc_cont)>0){
    infectors_summary_vacc_cont <- aggregate(x = results_edge_list_analysis_vacc_cont, 
                                             by = list(unique.values = results_edge_list_analysis_vacc_cont$Infector), 
                                             FUN = length)
    infectors_summary_vacc_cont <- infectors_summary_vacc_cont[,1:2]
  } else { 
    infectors_summary_vacc_cont <- as.data.frame(cbind(NA,NA))
  }
  names(infectors_summary_vacc_cont) <- c("infector","# infectees")
  
  
  # get summary of control infecting vaccinated
  if (nrow(results_edge_list_analysis_cont_vacc)>0){
    infectors_summary_cont_vacc <- aggregate(x = results_edge_list_analysis_cont_vacc, 
                                             by = list(unique.values = results_edge_list_analysis_cont_vacc$Infector), 
                                             FUN = length)
    infectors_summary_cont_vacc <- infectors_summary_cont_vacc[,1:2]
  } else { 
    infectors_summary_cont_vacc <- as.data.frame(cbind(NA,NA))
  }
  names(infectors_summary_cont_vacc) <- c("infector","# infectees")
  
  
  # get summary of control infecting control
  if (nrow(results_edge_list_analysis_cont_cont)>0){
    infectors_summary_cont_cont <- aggregate(x = results_edge_list_analysis_cont_cont, 
                                             by = list(unique.values = results_edge_list_analysis_cont_cont$Infector), 
                                             FUN = length)
    infectors_summary_cont_cont <- infectors_summary_cont_cont[,1:2]
  } else { 
    infectors_summary_cont_cont <- as.data.frame(cbind(NA,NA))
  }
  names(infectors_summary_cont_cont) <- c("infector","# infectees")
  
  vacc_vacc_inf <- sum(infectors_summary_vacc_vacc$`# infectees`)
  if (is.na(vacc_vacc_inf)){
    vacc_vacc_inf<-0
  }
  
  cont_vacc_inf <- sum(infectors_summary_cont_vacc$`# infectees`)
  if (is.na(cont_vacc_inf)){
    cont_vacc_inf<-0
  }
  
  vacc_cont_inf <- sum(infectors_summary_vacc_cont$`# infectees`)
  if (is.na(vacc_cont_inf)){
    vacc_cont_inf<-0
  }
  
  cont_cont_inf <- sum(infectors_summary_cont_cont$`# infectees`)
  if (is.na(cont_cont_inf)){
    cont_cont_inf<-0
  }
  
  importation <- nrow(results_edge_list_analysis[is.na(results_edge_list_analysis$infector_stat) & results_edge_list_analysis$eventstatus==1,])
  
  # define RR by # secondary infections by vaccinated/by controls stratified by vaccination status of infected
  # RR_vacc <- vacc_vacc_inf/cont_vacc_inf
  # RR_cont <- vacc_cont_inf/cont_cont_inf
  # prop_events_vacc <- numevents_vacc/(numevents_vacc+numevents_cont)
  # 
  # # combine weighting by proportion of events among vaccinated and controls
  # RR_inf <- RR_vacc*prop_events_vacc + RR_cont*(1-prop_events_vacc)
  #RR_inf_var
  
  RR_inf <- (vacc_vacc_inf+vacc_cont_inf)/(cont_vacc_inf+cont_cont_inf)
  
  infectors_summary <-  aggregate(x = results_edge_list_analysis,by = list(unique.values = results_edge_list_analysis$Infector),FUN = length)
  infectors_summary <- infectors_summary[,1:2]
  names(infectors_summary) <- c("infector","# infectees")
  
  # summmary statistics of Re
  mean_R <- mean(infectors_summary$`# infectees`)
  #cat("Average R: ",mean_R,"\n")
  med_R <- median(infectors_summary$`# infectees`)
  #cat("Median R: ",med_R,"\n")
  
  list(vaccEffEst,pval,numevents_vacc,numevents_cont,sample_size,
       RR_inf,mean_R,med_R,vacc_vacc_inf,vacc_cont_inf,cont_vacc_inf,cont_cont_inf,importation,results_edge_list_analysis,
       results_infected)
  
}