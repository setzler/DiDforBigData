
DiDge_main <- function(inputdata, varnames, cohort_time, event_postperiod, baseperiod = -1, control_group = "all", return_data=FALSE, forceOLS=TRUE, robust=FALSE){

  # set up variable names
  time_name = varnames$time_name
  outcome_name = varnames$outcome_name
  cohort_name = varnames$cohort_name
  id_name = varnames$id_name
  covariate_names = varnames$covariate_names
  cluster_names = varnames$cluster_names
  keep_vars = c(outcome_name,covariate_names)
  all_keep_vars = c(outcome_name,covariate_names,cluster_names)

  # prepare time periods
  pre_time = cohort_time + baseperiod
  post_time = cohort_time + event_postperiod
  time_set = sort(inputdata[get(cohort_name) == cohort_time, unique(get(time_name))])
  if(!(pre_time %in% time_set)){
    stop(print(sprintf("error: for cohort %s, preperiod %s is unavailable.",cohort_time,baseperiod)))
  }
  if(!(post_time %in% time_set)){
    stop(print(sprintf("error: for cohort %s, postperiod %s is unavailable.",cohort_time,event_postperiod)))
  }

  # restrict to pre and post time periods of interest, then restrict to observations present in both time periods
  inputdata = inputdata[get(time_name)==pre_time | get(time_name)==post_time]
  inputdata = inputdata[,.SD,.SDcols=c(id_name, time_name, cohort_name, all_keep_vars)]
  nn0 = nrow(inputdata)
  inputdata = na.omit(inputdata) # remove any rows with missing values
  nn1 = nrow(inputdata)
  if(nn1 < nn0){
    warning(sprintf("%s out of %s observations dropped due to missing values",(nn0-nn1),nn0))
  }
  inputdata[, nobs := .N, by=id_name]
  inputdata = inputdata[nobs==2 | (get(time_name)==pre_time & get(time_name)==post_time)] # the second condition keeps the base period
  gc()

  # define treated data and get means
  treated_data_prepost = merge(
    inputdata[(get(cohort_name) == cohort_time) & (get(time_name) == pre_time), .SD, .SDcols=c(id_name,all_keep_vars)],
    inputdata[(get(cohort_name) == cohort_time) & (get(time_name) == post_time), .SD, .SDcols=c(id_name,keep_vars)],
    by=c(id_name)
  )
  names(treated_data_prepost) = gsub("\\.x","_pre",names(treated_data_prepost))
  names(treated_data_prepost) = gsub("\\.y","_post",names(treated_data_prepost))
  setnames(treated_data_prepost, paste0(keep_vars,"_pre"), paste0("treated_",keep_vars,"_pre"))
  setnames(treated_data_prepost, paste0(keep_vars,"_post"), paste0("treated_",keep_vars,"_post"))
  treated_data_prepost[, treated := 1.0]

  # define control data and get means
  control_data_prepost = NULL
  if(control_group == "all"){
    control_data_prepost = merge(
      inputdata[(get(cohort_name) > max(post_time,cohort_time)) & (get(time_name) == pre_time), .SD, .SDcols=c(id_name,all_keep_vars)],
      inputdata[(get(cohort_name) > max(post_time,cohort_time)) & (get(time_name) == post_time), .SD, .SDcols=c(id_name,keep_vars)],
      by=c(id_name)
    )
  }
  if(control_group == "never-treated"){
    control_data_prepost = merge(
      inputdata[(get(cohort_name) > max(post_time,cohort_time)) & is.infinite(get(cohort_name)) & (get(time_name) == pre_time), .SD, .SDcols=c(id_name,all_keep_vars)],
      inputdata[(get(cohort_name) > max(post_time,cohort_time)) & is.infinite(get(cohort_name)) & (get(time_name) == post_time), .SD, .SDcols=c(id_name,keep_vars)],
      by=c(id_name)
    )
  }
  if(control_group == "future-treated"){
    control_data_prepost = merge(
      inputdata[(get(cohort_name) > max(post_time,cohort_time)) & is.finite(get(cohort_name)) & (get(time_name) == pre_time), .SD, .SDcols=c(id_name,all_keep_vars)],
      inputdata[(get(cohort_name) > max(post_time,cohort_time)) & is.finite(get(cohort_name)) & (get(time_name) == post_time), .SD, .SDcols=c(id_name,keep_vars)],
      by=c(id_name)
    )
  }
  names(control_data_prepost) = gsub("\\.x","_pre",names(control_data_prepost))
  names(control_data_prepost) = gsub("\\.y","_post",names(control_data_prepost))
  setnames(control_data_prepost, paste0(keep_vars,"_pre"), paste0("control_",keep_vars,"_pre"))
  setnames(control_data_prepost, paste0(keep_vars,"_post"), paste0("control_",keep_vars,"_post"))
  control_data_prepost[, treated := 0.0]

  # levels
  treated_pre_outcome = paste0("treated_",outcome_name,"_pre")
  treated_post_outcome = paste0("treated_",outcome_name,"_post")
  treated_pre = treated_data_prepost[, list(Ntreated_pre=sum(!is.na(get(treated_pre_outcome))),Etreated_pre=mean(get(treated_pre_outcome)))]
  treated_post = treated_data_prepost[, list(Ntreated_post=sum(!is.na(get(treated_post_outcome))),Etreated_post=mean(get(treated_post_outcome)),Etreated_var=var(get(treated_post_outcome)))]
  treated_results = cbind(treated_pre,treated_post)
  treated_results$treated_diff_var = treated_data_prepost[, var(get(treated_post_outcome) - get(treated_pre_outcome))]
  treated_results[, Ntreated_pre := NULL]
  setnames(treated_results,"Ntreated_post","Ntreated")

  control_pre_outcome = paste0("control_",outcome_name,"_pre")
  control_post_outcome = paste0("control_",outcome_name,"_post")
  control_pre = control_data_prepost[, list(Ncontrol_pre=sum(!is.na(get(control_pre_outcome))),Econtrol_pre=mean(get(control_pre_outcome)))]
  control_post = control_data_prepost[, list(Ncontrol_post=sum(!is.na(get(control_post_outcome))),Econtrol_post=mean(get(control_post_outcome)),Econtrol_var=var(get(control_post_outcome)))]
  control_results = cbind(control_pre,control_post)
  control_results$control_diff_var = control_data_prepost[, var(get(control_post_outcome) - get(control_pre_outcome))]
  control_results[, Ncontrol_pre := NULL]
  setnames(control_results,"Ncontrol_post","Ncontrol")

  # set up results
  results = cbind(control_results,treated_results)
  results[, Cohort := cohort_time]
  results[, Preperiod := pre_time]
  results[, CalendarTime := post_time]
  results[, Baseperiod := baseperiod]
  results[, EventTime := event_postperiod]
  results[, Econtrol_SE := sqrt(Econtrol_var/(Ncontrol-1))]
  results[, Etreated_SE := sqrt(Etreated_var/(Ntreated-1))]
  results[, pred_Etreated_post := Etreated_pre + (Econtrol_post - Econtrol_pre)]
  results[, ATTge := (Etreated_post - pred_Etreated_post)]
  results[, ATTge_SE := sqrt(treated_diff_var/(Ntreated-1) + control_diff_var/(Ncontrol-1))]

  results_variables_order = c("Cohort","EventTime","Baseperiod","CalendarTime","ATTge","ATTge_SE",
                              "Econtrol_pre","Econtrol_post","Econtrol_SE",
                              "Etreated_pre","Etreated_post","Etreated_SE","pred_Etreated_post",
                              "Ncontrol","Ntreated")

  # OLS
  data_prepost = NULL
  if(!is.null(covariate_names) | !is.null(cluster_names) | return_data | forceOLS){
    names(treated_data_prepost) = gsub("treated_","",names(treated_data_prepost))
    names(control_data_prepost) = gsub("control_","",names(control_data_prepost))
    data_prepost = rbindlist(list(treated_data_prepost,control_data_prepost))
    for(ii in keep_vars){
      data_prepost[, (paste0(ii,"_diff")) := get(paste0(ii,"_post")) - get(paste0(ii,"_pre"))]
    }
    data_prepost = data_prepost[,.SD,.SDcols=c(id_name,"treated",paste0(keep_vars,"_diff"),cluster_names)]
  }

  if(!is.null(covariate_names) | !is.null(cluster_names) | forceOLS){
    # reg formula, no covariates
    OLSformula = paste0(paste0(outcome_name,"_diff"), paste0(" ~ treated"))
    if(!is.null(covariate_names)){
      # reg formula, with covariates
      OLSformula = paste0(OLSformula, " + ", paste0(paste0(covariate_names,"_diff"),collapse=" + "))
      # since covariates change the estimates, keep a copy of the old estimate
      results[, ATTge_nocovars := copy(ATTge)]
      results[, ATTge_SE_nocovars := copy(ATTge_SE)]
      results_variables_order = c(results_variables_order,"ATTge_nocovars","ATTge_SE_nocovars")
    }
    # execute the regression
    OLSlm = lm(as.formula(OLSformula),data=data_prepost)
    newATT = as.numeric(OLSlm$coefficients["treated"])
    # check if the treated coefficient is missing
    if(!is.na(newATT)){
      results[, ATTge := newATT]
      if(is.null(cluster_names) & !robust){
        OLSvcov = vcov(OLSlm)
      }
      if(!is.null(cluster_names)){
        library(sandwich, warn.conflicts = F, quietly = T)
        CLformula = as.formula(paste0(" ~ ", paste0(cluster_names, collapse=" + ")))
        OLSvcov = vcovCL(OLSlm, cluster = CLformula)
      }
      if(is.null(cluster_names) & robust){
        library(sandwich, warn.conflicts = F, quietly = T)
        OLSvcov = vcovHC(OLSlm, "HC1")
      }
      OLSvcov = OLSvcov["treated", "treated"]
      newATTSE = sqrt(as.numeric(OLSvcov))
      results[, ATTge_SE := newATTSE]
    }
    if(is.na(newATT)){
      results[, ATTge := NA]
      results[, ATTge_SE := NA]
    }
  }

  # combine means into an output table
  results = results[,.SD,.SDcols=results_variables_order]
  results = results[order(Cohort,EventTime)]
  if(return_data){
    data_prepost[, Cohort := cohort_time]
    data_prepost[, EventTime := event_postperiod]
    export_vars = c(id_name,"Cohort","EventTime","treated",paste0(keep_vars,"_diff"),cluster_names)
    data_prepost = data_prepost[,.SD,.SDcols=export_vars]
    return(list(results=results,data_prepost=data_prepost))
  }
  return(results)
}

DiDge_bins <- function(inputdata, varnames, cohort_time, event_postperiod, baseperiod = -1, control_group = "all"){

  # set up variable names
  time_name = varnames$time_name
  outcome_name = varnames$outcome_name
  cohort_name = varnames$cohort_name
  id_name = varnames$id_name
  bin_name = varnames$bin_name

  # get bins
  pre_time = cohort_time + baseperiod
  post_time = cohort_time + event_postperiod
  bin_set_treated_pre = inputdata[(get(cohort_name) == cohort_time) & (get(time_name) == pre_time)][,.N,bin][N>1][,sort(bin)]
  bin_set_treated_post = inputdata[(get(cohort_name) == cohort_time) & (get(time_name) == pre_time)][,.N,bin][N>1][,sort(bin)]
  bin_set_control_pre = inputdata[(get(cohort_name) > max(post_time,cohort_time)) & (get(time_name) == pre_time)][,.N,bin][N>1][,sort(bin)]
  bin_set_control_post = inputdata[(get(cohort_name) > max(post_time,cohort_time)) & (get(time_name) == post_time)][,.N,bin][N>1][,sort(bin)]
  bin_set = Reduce(intersect, list(bin_set_treated_pre,bin_set_treated_post,bin_set_control_pre,bin_set_control_post))

  # loop DiD over bins
  results_bins = NULL
  for(binval in bin_set){
    res = DiDge_main(inputdata[get(bin_name)==binval], cohort_time = cohort_time, event_postperiod = event_postperiod, baseperiod = baseperiod,
                     varnames=varnames, control_group = control_group)
    res[, (bin_name) := binval]
    results_bins = rbindlist(list(results_bins,res))
  }

  # take the average across bins
  original_names = names(results_bins)
  original_names = original_names[original_names != bin_name]
  results_bins[, Ntreated_bin := sum(Ntreated), list(Cohort,EventTime,Baseperiod,CalendarTime)]
  results_bins[, bin_weights := Ntreated/Ntreated_bin]
  results_average = results_bins[, list(ATTge=sum(bin_weights * ATTge),
                                        ATTge_SE=sqrt(sum(bin_weights^2 * ATTge_SE^2)),
                                        Etreated_post=sum(bin_weights*Etreated_post),
                                        Etreated_pre=sum(bin_weights*Etreated_pre),
                                        Etreated_SE=sqrt(sum(bin_weights^2*Etreated_SE^2)),
                                        Econtrol_post=sum(bin_weights*Econtrol_post),
                                        Econtrol_pre=sum(bin_weights*Econtrol_pre),
                                        Econtrol_SE=sqrt(sum(bin_weights^2*Econtrol_SE^2)),
                                        pred_Etreated_post=sum(bin_weights*pred_Etreated_post),
                                        Ntreated=sum(Ntreated),
                                        Ncontrol=sum(Ncontrol)
  ), list(Cohort,EventTime,Baseperiod,CalendarTime)][order(Cohort,EventTime)]
  results_average = results_average[,.SD,.SDcols=original_names]
  return(results_average)
}


#' Estimate DiD for a single cohort (g) and a single event time (e).
#'
#' @param inputdata A data.table.
#' @param varnames A list of the form varnames = list(id_name, time_name, outcome_name, cohort_name), where all four arguments of the list must be a character that corresponds to a variable name in inputdata.
#' @param control_group There are three possibilities: control_group="never-treated" uses the never-treated control group only; control_group="future-treated" uses those units that will receive treatment in the future as the control group; and control_group="all" uses both the never-treated and the future-treated in the control group. Default is control_group="all".
#' @param baseperiod This is the base pre-period that is normalized to zero in the DiD estimation. Default is baseperiod=-1.
#' @param min_event This is the minimum event time (e) to estimate. Default is NULL, in which case, no minimum is imposed.
#' @param max_event This is the maximum event time (e) to estimate. Default is NULL, in which case, no maximum is imposed.
#' @param return_data If true, this returns the treated and control differenced data. Default is FALSE.
#' @param forceOLS Compute standard errors using OLS, even if analytic expression is available.
#' @param robust Compute robust standard errors, using the same "HC1" approach used by the Stata ", robust" option.
#' @returns A list with two components: results_cohort is a data.table with the DiDge estimates (by event e and cohort g), and results_average is a data.table with the DiDe estimates (by event e, average across cohorts g).
#' @examples
#' # simulate some data
#' simdata = SimDiD(sample_size=200)$simdata
#'
#' # define the variable names as a list()
#' varnames = list()
#' varnames$time_name = "year"
#' varnames$outcome_name = "Y"
#' varnames$cohort_name = "cohort"
#' varnames$id_name = "id"
#'
#' # estimate the ATT for cohort 2007 at event time 1
#' DiDge(simdata, varnames, cohort_time=2007, event_postperiod=1)
#'
#' # change the base period to -3
#' DiDge(simdata, varnames, baseperiod=-3, cohort_time=2007, event_postperiod=1)
#'
#' # use only the never-treated control group
#' DiDge(simdata, varnames, control_group = "never-treated", cohort_time=2007, event_postperiod=1)
#'
#' # use only the never-treated control group
#' DiDge(simdata, varnames, control_group = "future-treated", cohort_time=2007, event_postperiod=1)
#'
#' # simulate some data with covariates, add the covariates to the varnames, update the estimates
#' sim = SimDiD(sample_size=200,time_covars=TRUE)
#' varnames$covariate_names = c("X1","X2")
#' DiDge(inputdata=copy(sim$simdata), varnames, cohort_time=2007, event_postperiod = 3)
#'
#' varnames$covariate_names = NULL # we are done with this example
#'
#' # simulate some data with bins, add the bins to the varnames, update the estimates
#' sim = SimDiD(sample_size=3000, bin_covars=TRUE)
#'
#' # if we do not add the bin_name to varnames, it will ignore bins and give biased estimates:
#' DiDge(inputdata=copy(sim$simdata), varnames, cohort_time=2007, event_postperiod = 3)
#'
#' # now we account for bins:
#' varnames$bin_name = c("bin")
#' DiDge(inputdata=copy(sim$simdata), varnames, cohort_time=2007, event_postperiod = 3)
#'
#' # now we cluster on bins:
#' varnames$bin_name = NULL
#' varnames$cluster_names = c("bin")
#' DiDge(inputdata=copy(sim$simdata), varnames, cohort_time=2007, event_postperiod = 3)
#'
#' @export
DiDge <- function(inputdata, varnames, cohort_time, event_postperiod, baseperiod = -1, control_group = "all", return_data=FALSE, forceOLS=TRUE, robust=FALSE){
  if(is.null(varnames$bin_name)){
    return(DiDge_main(inputdata=inputdata, varnames=varnames, cohort_time=cohort_time, event_postperiod=event_postperiod, baseperiod=baseperiod, control_group = control_group, return_data=return_data, forceOLS=forceOLS, robust=robust))
  }
  if(!is.null(varnames$bin_name)){
    return(DiDge_bins(inputdata=inputdata, varnames=varnames, cohort_time=cohort_time, event_postperiod=event_postperiod, baseperiod=baseperiod, control_group = control_group))
  }
}

