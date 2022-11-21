
DiD_getSEs_EventTime_noOLS <- function(data_cohort,varnames){

  # check that there are treated and control units available
  data_cohort[,available_treated := sum(treated), list(Cohort,EventTime)]
  data_cohort[,available_control := sum(1-treated), list(Cohort,EventTime)]
  data_cohort = data_cohort[available_treated > 1 & available_control > 1]
  data_cohort[,available_treated := NULL]
  data_cohort[,available_control := NULL]

  # set up variable names
  time_name = varnames$time_name
  outcome_name = paste0(varnames$outcome_name,"_diff")
  cohort_name = varnames$cohort_name
  id_name = varnames$id_name

  # get the SE
  all_events = data_cohort[,sort(unique(EventTime))]
  covmat_list = list()
  weight_list = list()
  for(event in all_events){
    data_event = data_cohort[EventTime==event]
    cohorts = data_event[,sort(unique(Cohort))]
    numcc = length(cohorts)
    covmat_treated_treated = matrix(0,nrow=numcc,ncol=numcc)
    covmat_control_control = matrix(0,nrow=numcc,ncol=numcc)
    covmat_control_treated = matrix(0,nrow=numcc,ncol=numcc)
    covmat_treated_control = matrix(0,nrow=numcc,ncol=numcc)
    weight_treated = rep(0,numcc)
    for(cc_row_iter in 1:numcc){
      cc_row_val = cohorts[cc_row_iter]
      for(cc_col_iter in 1:numcc){
        cc_col_val = cohorts[cc_col_iter]
        if(cc_col_val == cc_row_val){
          # diagonals
          covmat_treated_treated[cc_row_iter, cc_col_iter] = data_event[Cohort==cc_row_val & treated==1, var(get(outcome_name))/(.N)]
          covmat_control_control[cc_row_iter, cc_col_iter] = data_event[Cohort==cc_row_val & treated==0, var(get(outcome_name))/(.N)]
          weight_treated[cc_row_iter] = data_event[Cohort==cc_row_val & treated==1, .N]
        }
        if(cc_col_val < cc_row_val){ # we will make this symmetric at the end
          # control-control
          cc_row_control = data_event[Cohort==cc_row_val & treated==0]
          cc_col_control = data_event[Cohort==cc_col_val & treated==0]
          cc_row_control_col_control = merge(cc_row_control, cc_col_control, by=id_name)
          if(nrow(cc_row_control_col_control)){
            covmat_control_control[cc_row_iter, cc_col_iter] = (cc_row_control_col_control[, cov(get(paste0(outcome_name,".x")), get(paste0(outcome_name,".y")))]/sqrt(cc_row_control[,as.numeric(.N)]))* 1/sqrt(cc_col_control[,as.numeric(.N)])
          }
        }
        if(cc_col_val > cc_row_val){ # no past treated group can appear in the current or any future control group
          # treated-control
          cc_row_control = data_event[Cohort==cc_row_val & treated==0]
          cc_col_treated = data_event[Cohort==cc_col_val & treated==1]
          cc_row_control_col_treated = merge(cc_row_control, cc_col_treated, by=id_name)
          if(nrow(cc_row_control_col_treated)){
            covmat_control_treated[cc_row_iter, cc_col_iter] = (cc_row_control_col_treated[, cov(get(paste0(outcome_name,".x")), get(paste0(outcome_name,".y")))]/sqrt(cc_row_control[,as.numeric(.N)])) * 1/sqrt(cc_col_treated[,as.numeric(.N)])
          }
        }
      }
    }
    makeSymm <- function(m) {
      m[upper.tri(m)] <- t(m)[upper.tri(m)]
      return(m)
    }
    covmat_control_control = makeSymm(covmat_control_control)
    covmat_cohort = rbind(
      cbind(covmat_treated_treated, t(covmat_control_treated)),
      cbind(covmat_control_treated, covmat_control_control)
    )
    covmat_list[[as.character(event)]] = covmat_cohort
    weight_list[[as.character(event)]] = c(weight_treated/sum(weight_treated), -weight_treated/sum(weight_treated))
  }

  # collect the event-specific DiDe SEs
  ATTe_SEs = data.table()
  for(event in all_events){
    covmat_cohort = covmat_list[[as.character(event)]]
    weight_treated = as.matrix(weight_list[[as.character(event)]])
    ATTe_SE = sqrt(as.numeric((t(weight_treated) %*% covmat_cohort) %*% weight_treated))
    ATTe_SEs = rbindlist(list(ATTe_SEs, data.table(EventTime=event, ATTe_SE)))
  }

  return(ATTe_SEs)
}

getSEs_multipleEventTimes_noOLS <- function(data_cohort,varnames,Eset){

  # check that there are treated and control units available
  data_cohort[,available_treated := sum(treated), list(Cohort,EventTime)]
  data_cohort[,available_control := sum(1-treated), list(Cohort,EventTime)]
  data_cohort = data_cohort[available_treated > 1 & available_control > 1]
  data_cohort[,available_treated := NULL]
  data_cohort[,available_control := NULL]

  # set up variable names
  time_name = varnames$time_name
  outcome_name = paste0(varnames$outcome_name,"_diff")
  cohort_name = varnames$cohort_name
  id_name = varnames$id_name

  # get the SE
  data_cohortevents = data_cohort[EventTime %in% Eset]
  cohortevents = data.table()
  for(ee in Eset){
    cohortevents = rbindlist(list(cohortevents, data.table(cohort=data_cohortevents[EventTime==ee,sort(unique(Cohort))],event=ee)))
  }
  cohortevents = cohortevents[order(cohort,event)]
  cohortevents1 = copy(cohortevents)[, treated := 1]
  cohortevents0 = copy(cohortevents)[, treated := 0]
  cohortevents = rbindlist(list(cohortevents1,cohortevents0))
  numce = nrow(cohortevents)
  cohortevents_covmat = matrix(0, numce, numce)
  cohortevents_weights = rep(NA,numce/2)
  cohortevents_means = rep(NA,numce)

  for(cc_row_iter in 1:numce){
    # row
    cohort_row = cohortevents[cc_row_iter]$cohort
    event_row = cohortevents[cc_row_iter]$event
    treated_row = cohortevents[cc_row_iter]$treated
    data_row = data_cohortevents[Cohort == cohort_row & EventTime==event_row & treated==treated_row]
    # row weight
    cohortevents_means[cc_row_iter] = data_row[,mean(get(outcome_name))]
    if(treated_row==1){
      cohortevents_weights[cc_row_iter] = data_row[, .N]
    }
    for(cc_col_iter in 1:numce){
      # column
      cohort_col = cohortevents[cc_col_iter]$cohort
      event_col = cohortevents[cc_col_iter]$event
      treated_col = cohortevents[cc_col_iter]$treated
      data_col = data_cohortevents[Cohort == cohort_col & EventTime==event_col & treated==treated_col]
      # covariance
      data_rowcol_merge = merge(data_row, data_col, by=id_name)
      if(nrow(data_rowcol_merge)>0){
        cohortevents_covmat[cc_row_iter,cc_col_iter] = (data_rowcol_merge[, cov(get(paste0(outcome_name,".x")), get(paste0(outcome_name,".y")))]/sqrt(data_row[,as.numeric(.N)])) * 1/sqrt(data_col[,as.numeric(.N)])
      }
    }
  }
  cohortevents_weights = cohortevents_weights/sum(cohortevents_weights)
  cohortevents_weights = c(cohortevents_weights,-1*cohortevents_weights)
  cohortevents_weights = as.matrix(cohortevents_weights)
  ATT_Eset = as.numeric(t(cohortevents_weights) %*% cohortevents_means)
  ATT_Eset_SE = sqrt(as.numeric((t(cohortevents_weights) %*% cohortevents_covmat) %*% (cohortevents_weights)))
  return(data.table(Eset = paste0(Eset,collapse=","), ATT_Eset=ATT_Eset, ATT_Eset_SE=ATT_Eset_SE))
}

