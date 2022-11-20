
DiD_getSEs_EventTime_OLS <- function(data_cohort,varnames,robust=FALSE){

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
  cluster_names = varnames$cluster_names
  covariate_names = NULL
  if(!is.null(varnames$covariate_names)){
    covariate_names = paste0(varnames$covariate_names,"_diff")
  }

  # get the SE
  all_events = data_cohort[,sort(unique(EventTime))]
  ATTe_SEs = data.table()
  for(ee in all_events){
    intercepts = NULL
    treateds = NULL
    treated_weights = NULL
    covariates = NULL
    data_event = data_cohort[EventTime==ee]
    all_cohorts = data_event[,unique(Cohort)]
    for(cc in all_cohorts){
      intercepts = c(intercepts, sprintf("intercept_%s",cc))
      treateds = c(treateds, sprintf("treated_%s",cc))
      data_event[, (sprintf("intercept_%s",cc)) := as.numeric((Cohort==cc))]
      data_event[, (sprintf("treated_%s",cc)) := as.numeric((treated==1)*(Cohort==cc))]
      treated_weights = c(treated_weights, data_event[, sum(get(sprintf("treated_%s",cc))) ] )
      if(!is.null(covariate_names)){
        for(vv in covariate_names){
          covariates = c(covariates, sprintf("%s_%s",vv,cc))
          data_event[, (sprintf("%s_%s",vv,cc)) := as.numeric((Cohort==cc)*get(vv))]
        }
      }
    }

    OLSformula = paste0(outcome_name," ~ -1 + ",
                        paste0(treateds, collapse=" + "),
                        " + ",
                        paste0(intercepts, collapse=" + "))
    if(!is.null(covariate_names)){
      OLSformula = paste0(OLSformula,
                          " + ",
                          paste0(covariates, collapse=" + "))
    }

    # regression
    OLSvcov = NULL
    OLSlm = lm(OLSformula, data=data_event)
    OLSmeans = as.numeric(OLSlm$coefficients[treateds])

    # covariance matrix for errors
    if(is.null(cluster_names) & !robust){
      OLSvcov = vcov(OLSlm)
    }
    if(!is.null(cluster_names)){
      library(sandwich, warn.conflicts = F, quietly = T)
      data_event <<- copy(data_event) # due to a well-known scoping bug in R's base lm.predict that no one will fix despite years of requests, this redundancy is necessary!
      CLformula = as.formula(paste0(" ~ ", paste0(cluster_names, collapse=" + ")))
      OLSvcov = vcovCL(OLSlm, cluster = CLformula)
    }
    if(robust & is.null(cluster_names)){
      library(sandwich, warn.conflicts = F, quietly = T)
      OLSvcov = vcovHC(OLSlm, "HC1")
    }

    # finish
    OLSvcov = OLSvcov[treateds, treateds]
    OLSvcov = as.matrix(OLSvcov)
    treated_weights = treated_weights/sum(treated_weights)
    ATTe = as.numeric(t(treated_weights) %*% OLSmeans)
    ATTe_SE = sqrt(as.numeric((t(treated_weights) %*% OLSvcov) %*% treated_weights))
    ATTe_SEs = rbindlist(list(ATTe_SEs, data.table(EventTime=ee, ATTe_OLS=ATTe, ATTe_SE=ATTe_SE)))
  }

  return(ATTe_SEs)
}


getSEs_multipleEventTimes_OLS <- function(data_cohort,varnames,robust=FALSE,Eset){

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
  cluster_names = varnames$cluster_names
  covariate_names = NULL
  if(!is.null(varnames$covariate_names)){
    covariate_names = paste0(varnames$covariate_names,"_diff")
  }

  # get the SE
  data_event = data_cohort[EventTime %in% Eset]
  cohortevents = data.table()
  for(ee in Eset){
    cohortevents = rbindlist(list(cohortevents, data.table(cohort=data_event[EventTime==ee,sort(unique(Cohort))],event=ee)))
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
    data_event[, (sprintf("intercept_%s_%s",cc,ee)) := (Cohort==cc)*(EventTime==ee)]
    data_event[, (sprintf("treated_%s_%s",cc,ee)) := (treated==1)*(Cohort==cc)*(EventTime==ee)]
    treated_weights = c(treated_weights, data_event[, sum(get(sprintf("treated_%s_%s",cc,ee))) ] )
    if(!is.null(covariate_names)){
      for(vv in covariate_names){
        covariates = c(covariates, sprintf("%s_%s_%s",vv,cc,ee))
        data_event[, (sprintf("%s_%s_%s",vv,cc,ee)) := (Cohort==cc)*(EventTime==ee)*get(vv)]
      }
    }
  }

  OLSformula = paste0(outcome_name," ~ -1 + ",
                      paste0(treateds, collapse=" + "),
                      " + ",
                      paste0(intercepts, collapse=" + "))
  if(!is.null(covariate_names)){
    OLSformula = paste0(OLSformula,
                        " + ",
                        paste0(covariates, collapse=" + "))
  }

  # regression
  OLSvcov = NULL
  OLSlm = lm(OLSformula, data=data_event)
  OLSmeans = as.numeric(OLSlm$coefficients[treateds])

  # covariance matrix for errors
  if(is.null(cluster_names) & !robust){
    OLSvcov = vcov(OLSlm)
  }
  if(!is.null(cluster_names)){
    library(sandwich, warn.conflicts = F, quietly = T)
    data_event <<- copy(data_event) # due to a well-known scoping bug in R's base lm.predict that no one will fix despite years of requests, this redundancy is necessary!
    CLformula = as.formula(paste0(" ~ ", paste0(cluster_names, collapse=" + ")))
    OLSvcov = vcovCL(OLSlm, cluster = CLformula)
  }
  if(robust & is.null(cluster_names)){
    library(sandwich, warn.conflicts = F, quietly = T)
    OLSvcov = vcovHC(OLSlm, "HC1")
  }

  OLSvcov = OLSvcov[treateds, treateds]
  treated_weights = treated_weights/sum(treated_weights)
  ATT_Eset = as.numeric(t(treated_weights) %*% OLSmeans)
  ATT_Eset_SE = sqrt(as.numeric((t(treated_weights) %*% OLSvcov) %*% treated_weights))
  return(data.table(Eset = paste0(Eset,collapse=","), ATT_Eset=ATT_Eset, ATT_Eset_SE=ATT_Eset_SE))
}

