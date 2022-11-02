
DiDe <- function(inputdata, varnames, control_group = "all", baseperiod=-1, min_event=NULL, max_event=NULL, return_data=FALSE){

  # set up variable names
  time_name = varnames$time_name
  outcome_name = varnames$outcome_name
  cohort_name = varnames$cohort_name
  id_name = varnames$id_name

  # set up and checks
  cohorts = sort(inputdata[, unique(get(cohort_name))])
  nevertreated_exist = sum(is.infinite(cohorts))>0
  cohorts = cohorts[!is.infinite(cohorts)]

  upperbound_outcomeyear = Inf

  if(control_group=="never-treated"){
    if(!nevertreated_exist){
      stop("You specified control_group='never-treated', but there are no never-treated observations. Note: you must code never-treated observations as infinity (Inf).")
    }
  }

  if(control_group=="future-treated"){
    if(length(cohorts) <= 1){
      stop("You specified control_group='future-treated', but there is only one treatment cohort with finite cohort time, so no control units are available based on future-treated observations.")
    }
    upperbound_outcomeyear = max(cohorts)
    cohorts = cohorts[!(cohorts == max(cohorts))]
  }

  # drop cohorts that are useless since they were treated before the baseperiod
  min_time_actual = inputdata[,min(get(time_name))]
  too_early_cohorts = cohorts - baseperiod < min_time_actual
  if(sum(too_early_cohorts) > 0){
    warning(sprintf("We cannot provide ATT estimates for cohort %s due to the absence of a base_preperiod in the data.", paste0(cohorts[too_early_cohorts], collapse=",")))
    cohorts_todrop = cohorts[too_early_cohorts]
    cohorts = cohorts[!too_early_cohorts]
    inputdata = inputdata[!(get(cohort_name) %in% cohorts_todrop)]
  }

  # cohort-specific estimation
  results_cohort = data.table()
  data_cohort = data.table()
  for(cc in cohorts){
    # set up cohort-specific event times
    times_for_cohort = sort(inputdata[get(cohort_name) == cc, unique(get(time_name))])
    times_for_cohort = times_for_cohort[times_for_cohort < upperbound_outcomeyear]
    event_periods = times_for_cohort - cc
    if(!is.null(min_event)){
      event_periods = event_periods[event_periods >= min_event]
    }
    if(!is.null(max_event)){
      event_periods = event_periods[event_periods <= max_event]
    }
    # loop DiD over the event times for this cohort
    for(event_postperiod in event_periods){
      res = DiDge(inputdata, cohort_time = cc, event_postperiod = event_postperiod, baseperiod = baseperiod,
                  varnames=varnames, control_group = control_group, return_data=TRUE)
      results_cohort = rbindlist(list(results_cohort,res$results))
      data_cohort = rbindlist(list(data_cohort, res$data))
    }
  }
  results_cohort=results_cohort[order(Cohort,EventTime)]
  data_cohort=data_cohort[order(Cohort,EventTime)]
  ATTe_SEs = DiD_getSEs_EventTime(data_cohort=data_cohort,varnames=varnames)

  # take the average across cohorts
  results_cohort[, Ntreated_event := sum(Ntreated), by="EventTime"]
  results_cohort[, cohort_weights := Ntreated/Ntreated_event]
  results_average = results_cohort[, list(ATTe=sum(cohort_weights * ATTge),
                                          ATTe_SE_nocorr=sqrt(sum(cohort_weights^2 * ATTge_SE^2)),
                                          Etreated_post=sum(cohort_weights*Etreated_post),
                                          Etreated_pre=sum(cohort_weights*Etreated_pre),
                                          Etreated_SE=sqrt(sum(cohort_weights^2*Etreated_SE^2)),
                                          Econtrol_post=sum(cohort_weights*Econtrol_post),
                                          Econtrol_pre=sum(cohort_weights*Econtrol_pre),
                                          Econtrol_SE=sqrt(sum(cohort_weights^2*Econtrol_SE^2)),
                                          Ntreated=sum(Ntreated),
                                          Ncontrol=sum(Ncontrol)
  ), list(EventTime,Baseperiod)][order(EventTime,Baseperiod)]
  results_average = merge(results_average, ATTe_SEs, by="EventTime")[order(EventTime,Baseperiod)]

  results_average = results_average[,.SD,.SDcols=c("EventTime","Baseperiod","ATTe","ATTe_SE","ATTe_SE_nocorr","Etreated_post","Etreated_pre","Etreated_SE","Econtrol_post","Econtrol_pre","Econtrol_SE","Ntreated","Ncontrol")]

  if(return_data){
    return(list(results_cohort=results_cohort, results_average=results_average, data_cohort=data_cohort))
  }
  return(list(results_cohort=results_cohort, results_average=results_average))
}

DiD_getSEs_EventTime <- function(data_cohort,varnames){

  # set up variable names
  time_name = varnames$time_name
  outcome_name = varnames$outcome_name
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
          covmat_treated_treated[cc_row_iter, cc_col_iter] = data_event[Cohort==cc_row_val & treated==1, var(Y_diff)/.N]
          covmat_control_control[cc_row_iter, cc_col_iter] = data_event[Cohort==cc_row_val & treated==0, var(Y_diff)/.N]
          weight_treated[cc_row_iter] = data_event[Cohort==cc_row_val & treated==1, .N]
        }
        if(cc_col_val < cc_row_val){ # we will make this symmetric at the end
          # control-control
          cc_row_control = data_event[Cohort==cc_row_val & treated==0]
          cc_col_control = data_event[Cohort==cc_col_val & treated==0]
          cc_row_control_col_control = merge(cc_row_control, cc_col_control, by=id_name)
          if(nrow(cc_row_control_col_control)){
            covmat_control_control[cc_row_iter, cc_col_iter] = (cc_row_control_col_control[, cov(Y_diff.x, Y_diff.y)]/sqrt(cc_row_control[,as.numeric(.N)]))* 1/sqrt(cc_col_control[,as.numeric(.N)])
          }
        }
        if(cc_col_val > cc_row_val){ # no past treated group can appear in the current or any future control group
          # treated-control
          cc_row_control = data_event[Cohort==cc_row_val & treated==0]
          cc_col_treated = data_event[Cohort==cc_col_val & treated==1]
          cc_row_control_col_treated = merge(cc_row_control, cc_col_treated, by=id_name)
          if(nrow(cc_row_control_col_treated)){
            covmat_control_treated[cc_row_iter, cc_col_iter] = (cc_row_control_col_treated[, cov(Y_diff.x, Y_diff.y)]/sqrt(cc_row_control[,as.numeric(.N)])) * 1/sqrt(cc_col_treated[,as.numeric(.N)])
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

  ATTe_SEs = data.table()
  for(event in all_events){
    covmat_cohort = covmat_list[[as.character(event)]]
    weight_treated = as.matrix(weight_list[[as.character(event)]])
    ATTe_SE = sqrt(as.numeric((t(weight_treated) %*% covmat_cohort) %*% weight_treated))
    ATTe_SEs = rbindlist(list(ATTe_SEs, data.table(EventTime=event, ATTe_SE)))
  }

  return(ATTe_SEs)
}

getSEs_multipleEventTimes <- function(data_cohort,varnames,Eset){

  # set up variable names
  time_name = varnames$time_name
  outcome_name = varnames$outcome_name
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
    cohortevents_means[cc_row_iter] = data_row[,mean(Y_diff)]
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
        cohortevents_covmat[cc_row_iter,cc_col_iter] = (data_rowcol_merge[, cov(Y_diff.x, Y_diff.y)]/sqrt(data_row[,as.numeric(.N)])) * 1/sqrt(data_col[,as.numeric(.N)])
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

getSEs_covariates_multipleEventTimes <- function(data_cohort,varnames,Eset){

  # set up variable names
  time_name = varnames$time_name
  outcome_name = paste0(varnames$outcome_name,"_diff")
  cohort_name = varnames$cohort_name
  id_name = varnames$id_name
  covariate_names = paste0(varnames$covariate_names,"_diff")

  # get the SE
  data_cohortevents = data_cohort[EventTime %in% Eset]
  cohortevents = data.table()
  for(ee in Eset){
    cohortevents = rbindlist(list(cohortevents, data.table(cohort=data_cohortevents[EventTime==ee,sort(unique(Cohort))],event=ee)))
  }
  cohortevents = cohortevents[order(cohort,event)]
  numce = nrow(cohortevents)
  cohortevents_covmat = matrix(0, numce, numce)
  cohortevents_weights = rep(NA,numce/2)
  cohortevents_means = rep(NA,numce)

  intercepts = c()
  treateds = c()
  covariates = c()
  treated_weights = c()
  for(rowiter in 1:numce){
    cc = cohortevents[rowiter][,cohort]
    ee = cohortevents[rowiter][,event]
    intercepts = c(intercepts, sprintf("intercept_%s_%s",cc,ee))
    treateds = c(treateds, sprintf("treated_%s_%s",cc,ee))
    data_cohortevents[, (sprintf("intercept_%s_%s",cc,ee)) := (Cohort==cc)*(EventTime==ee)]
    data_cohortevents[, (sprintf("treated_%s_%s",cc,ee)) := (treated==1)*(Cohort==cc)*(EventTime==ee)]
    treated_weights = c(treated_weights, data_cohortevents[, sum(get(sprintf("treated_%s_%s",cc,ee))) ] )
    for(vv in covariate_names){
      covariates = c(covariates, sprintf("%s_%s_%s",vv,cc,ee))
      data_cohortevents[, (sprintf("%s_%s_%s",vv,cc,ee)) := (Cohort==cc)*(EventTime==ee)*get(vv)]
    }
  }

  OLSformula = paste0(outcome_name," ~ -1 + ",
                      paste0(treateds, collapse=" + "),
                      " + ",
                      paste0(intercepts, collapse=" + "),
                      " + ",
                      paste0(covariates, collapse=" + "))

  OLSres = lm(OLSformula, data=data_cohortevents)
  OLSmeans = summary(OLSres)$coefficients
  OLSmeans = as.numeric(OLSmeans[treateds,][,"Estimate"])
  OLSvcov = vcov(OLSres)
  OLSvcov = OLSvcov[treateds, treateds]
  treated_weights = treated_weights/sum(treated_weights)
  ATT_Eset = as.numeric(t(treated_weights) %*% OLSmeans)
  ATT_Eset_SE = sqrt(as.numeric((t(treated_weights) %*% OLSvcov) %*% treated_weights))
  return(data.table(Eset = paste0(Eset,collapse=","), ATT_Eset=ATT_Eset, ATT_Eset_SE=ATT_Eset_SE))
}

#' Combine DiD estimates across cohorts and event times.
#'
#' @description
#' Estimate DiD for all possible cohorts and event time pairs (g,e), as well as the average across cohorts for each event time (e).
#'
#' @param inputdata A data.table.
#' @param varnames A list of the form varnames = list(id_name, time_name, outcome_name, cohort_name), where all four arguments of the list must be a character that corresponds to a variable name in inputdata.
#' @param control_group There are three possibilities: control_group="never-treated" uses the never-treated control group only; control_group="future-treated" uses those units that will receive treatment in the future as the control group; and control_group="all" uses both the never-treated and the future-treated in the control group. Default is control_group="all".
#' @param baseperiod This is the base pre-period that is normalized to zero in the DiD estimation. Default is baseperiod=-1.
#' @param min_event This is the minimum event time (e) to estimate. Default is NULL, in which case, no minimum is imposed.
#' @param max_event This is the maximum event time (e) to estimate. Default is NULL, in which case, no maximum is imposed.
#' @param Esets If a list of sets of event times is provided, it will loop over those sets, computing the average ATT_e across event times e. Default is NULL.
#' @returns A list with two components: results_cohort is a data.table with the DiDge estimates (by event e and cohort g), and results_average is a data.table with the DiDe estimates (by event e, average across cohorts g).
#' @examples
#' # simulate some data
#' simdata = SimDiD(sample_size=1000, ATTcohortdiff = 2)$simdata
#'
#' # define the variable names as a list()
#' varnames = list()
#' varnames$time_name = "year"
#' varnames$outcome_name = "Y"
#' varnames$cohort_name = "cohort"
#' varnames$id_name = "id"
#'
#' # estimate the ATT for all cohorts at event time 1 only
#' DiD(simdata, varnames, min_event=1, max_event=1)
#'
#' # change the base period to -3
#' DiD(simdata, varnames, baseperiod=-3, min_event=1, max_event=1)
#'
#' # check the pre-periods -4 through -2
#' DiD(simdata, varnames, control_group = "never-treated", min_event=-4, max_event=-2)
#'
#' # use only the never-treated control group, estimate events -4 through 1
#' DiD(simdata, varnames, control_group = "never-treated", min_event=-4, max_event=1)
#'
#' # use only the future treated control group, estimate events -4 through 1
#' DiD(simdata, varnames, control_group = "future-treated", min_event=-4, max_event=1)
#'
#' # estimate average ATTe across sets of events
#' DiD(simdata, varnames, min_event=-4, max_event=6, Esets=list(c(-3,-2,-1),c(1,2,3)))
#'
#' # simulate data with time-varying covariates
#' sim = SimDiD(sample_size=2000,time_covars=TRUE)
#' simdata = sim$simdata
#'
#' # add the time-varying covariate names to the varnames list
#' varnames$covariate_names = c("X1","X2")
#'
#' # run estimation that controls for time-varying covariates
#' DiD(simdata, varnames, min_event=1, max_event=2)
#'
#' # run estimation that controls for time-varying covariates
#' DiD(simdata, varnames, min_event=1, max_event=2, Esets=list(c(1,2)))
#' @export
DiD <- function(inputdata, varnames, control_group = "all", baseperiod=-1, min_event=NULL, max_event=NULL, Esets=NULL){
  if(is.null(Esets)){
    results = DiDe(inputdata=inputdata, varnames=varnames, control_group=control_group, baseperiod=baseperiod, min_event=min_event, max_event=max_event, return_data=FALSE)
    return(results)
  }
  if(!is.null(Esets)){
    results = DiDe(inputdata=inputdata, varnames=varnames, control_group=control_group, baseperiod=baseperiod, min_event=min_event, max_event=max_event, return_data=TRUE)
    data_cohort = results$data_cohort
    results_Esets = data.table()
    for(Eset in Esets){
      if(is.null(varnames$covariate_names)){
        results_Esets = rbindlist(list(results_Esets,getSEs_multipleEventTimes(data_cohort,varnames,Eset)))
      }
      if(!is.null(varnames$covariate_names)){
        results_Esets = rbindlist(list(results_Esets,getSEs_covariates_multipleEventTimes(data_cohort,varnames,Eset)))
      }
    }
    results$results_Esets = results_Esets
    results$data_cohort = NULL
    return(results)
  }
}

