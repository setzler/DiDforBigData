
#' Estimate DiD for a single cohort (g) and a single event time (e).
#'
#' @param inputdata A data.table.
#' @param varnames A list of the form varnames = list(id_name, time_name, outcome_name, cohort_name), where all four arguments of the list must be a character that corresponds to a variable name in inputdata.
#' @param control_group There are three possibilities: control_group="never-treated" uses the never-treated control group only; control_group="future-treated" uses those units that will receive treatment in the future as the control group; and control_group="all" uses both the never-treated and the future-treated in the control group. Default is control_group="all".
#' @param baseperiod This is the base pre-period that is normalized to zero in the DiD estimation. Default is baseperiod=-1.
#' @param min_event This is the minimum event time (e) to estimate. Default is NULL, in which case, no minimum is imposed.
#' @param max_event This is the maximum event time (e) to estimate. Default is NULL, in which case, no maximum is imposed.
#' @returns A list with two components: results_cohort is a data.table with the DiDge estimates (by event e and cohort g), and results_average is a data.table with the DiDe estimates (by event e, average across cohorts g).
#' @examples
#' # simulate some data
#' simdata = SimDiD(sample_size=200)
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
#' @export
DiDge <- function(inputdata, varnames, cohort_time, event_postperiod, baseperiod = -1, control_group = "all"){

  # set up variable names
  time_name = varnames$time_name
  outcome_name = varnames$outcome_name
  cohort_name = varnames$cohort_name
  id_name = varnames$id_name

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
  inputdata[, nobs := .N, by=id_name]
  inputdata = inputdata[nobs==2 | (get(time_name)==pre_time & get(time_name)==post_time)]
  gc()

  # define treated data and get means
  treated_data_prepost = merge(
    inputdata[(get(cohort_name) == cohort_time) & (get(time_name) == pre_time), list(tmp_idvar=get(id_name),treated_pre_outcome=get(outcome_name))],
    inputdata[(get(cohort_name) == cohort_time) & (get(time_name) == post_time), list(tmp_idvar=get(id_name),treated_post_outcome=get(outcome_name))],
    by=c("tmp_idvar")
  )
  treated_pre = treated_data_prepost[, list(Ntreated_pre=sum(!is.na(treated_pre_outcome)),Etreated_pre=mean(treated_pre_outcome))]
  treated_post = treated_data_prepost[, list(Ntreated_post=sum(!is.na(treated_post_outcome)),Etreated_post=mean(treated_post_outcome),Etreated_var=var(treated_post_outcome))]
  treated_results = cbind(treated_pre,treated_post)
  treated_results$treated_diff_var = treated_data_prepost[, var(treated_post_outcome - treated_pre_outcome)]

  # define control data and get means
  control_data_prepost = NULL
  if(control_group == "all"){
    control_data_prepost = merge(
      inputdata[(get(cohort_name) > max(post_time,cohort_time)) & (get(time_name) == pre_time), list(tmp_idvar=get(id_name),control_pre_outcome=get(outcome_name))],
      inputdata[(get(cohort_name) > max(post_time,cohort_time)) & (get(time_name) == post_time), list(tmp_idvar=get(id_name),control_post_outcome=get(outcome_name))],
      by=c("tmp_idvar")
    )
  }
  if(control_group == "never-treated"){
    control_data_prepost = merge(
      inputdata[(get(cohort_name) > max(post_time,cohort_time)) & is.infinite(get(cohort_name)) & (get(time_name) == pre_time), list(tmp_idvar=get(id_name),control_pre_outcome=get(outcome_name))],
      inputdata[(get(cohort_name) > max(post_time,cohort_time)) & is.infinite(get(cohort_name)) & (get(time_name) == post_time), list(tmp_idvar=get(id_name),control_post_outcome=get(outcome_name))],
      by=c("tmp_idvar")
    )
  }
  if(control_group == "future-treated"){
    control_data_prepost = merge(
      inputdata[(get(cohort_name) > max(post_time,cohort_time)) & is.finite(get(cohort_name)) & (get(time_name) == pre_time), list(tmp_idvar=get(id_name),control_pre_outcome=get(outcome_name))],
      inputdata[(get(cohort_name) > max(post_time,cohort_time)) & is.finite(get(cohort_name)) & (get(time_name) == post_time), list(tmp_idvar=get(id_name),control_post_outcome=get(outcome_name))],
      by=c("tmp_idvar")
    )
  }
  control_pre = control_data_prepost[, list(Ncontrol_pre=sum(!is.na(control_pre_outcome)),Econtrol_pre=mean(control_pre_outcome))]
  control_post = control_data_prepost[, list(Ncontrol_post=sum(!is.na(control_post_outcome)),Econtrol_post=mean(control_post_outcome),Econtrol_var=var(control_post_outcome))]
  control_results = cbind(control_pre,control_post)
  control_results$control_diff_var = control_data_prepost[, var(control_post_outcome - control_pre_outcome)]

  # combine means into an output table
  results = cbind(control_results,treated_results)
  results[, Cohort := cohort_time]
  results[, Preperiod := pre_time]
  results[, CalendarTime := post_time]
  results[, Baseperiod := baseperiod]
  results[, EventTime := event_postperiod]
  results[, Econtrol_SE := sqrt(Econtrol_var/Ncontrol_post)]
  results[, Etreated_SE := sqrt(Etreated_var/Ntreated_post)]
  results[, pred_Etreated_post := Etreated_pre + (Econtrol_post - Econtrol_pre)]
  results[, ATTge := (Etreated_post - pred_Etreated_post)]
  results[, ATTge_SE := sqrt(treated_diff_var/Ntreated_post + control_diff_var/Ncontrol_post)]
  results = results[,.(EventTime,Baseperiod,Cohort,CalendarTime,ATTge,ATTge_SE,
                       Econtrol_pre,Econtrol_post,Econtrol_SE,Etreated_pre,Etreated_post,Etreated_SE,pred_Etreated_post,
                       Ncontrol_pre,Ncontrol_post,Ntreated_pre,Ntreated_post)]
  return(results)
}



DiDe <- function(inputdata, varnames, control_group = "all", baseperiod=-1, min_event=NULL, max_event=NULL){

  # set up variable names
  time_name = varnames$time_name
  outcome_name = varnames$outcome_name
  cohort_name = varnames$cohort_name
  id_name = varnames$id_name

  # set up and checks
  cohorts = sort(inputdata[, unique(get(cohort_name))])
  nevertreated_exist = sum(is.infinite(cohorts))>0
  cohorts = cohorts[!is.infinite(cohorts)]

  if(control_group=="never-treated"){
    if(!nevertreated_exist){
      stop("You specified control_group='never-treated', but there are no never-treated observations. Note: you must code never-treated observations as infinity (Inf).")
    }
  }

  if(control_group=="future-treated"){
    if(length(cohorts) <= 1){
      stop("You specified control_group='future-treated', but there is only one treatment cohort with finite cohort time, so no control units are available based on will-be-treated observations.")
    }
    cohorts = cohorts[!(cohorts == max(cohorts))]
  }

  min_time_actual = inputdata[,min(get(time_name))]
  invalid_cohorts = cohorts - baseperiod < min_time_actual
  if(sum(invalid_cohorts) > 0){
    print(sprintf("We cannot provide ATT estimates for cohort %s due to the absence of a base_preperiod in the data.", paste0(cohorts[invalid_cohorts], collapse=",")))
    cohorts_todrop = cohorts[invalid_cohorts]
    cohorts = cohorts[!invalid_cohorts]
    inputdata = inputdata[!(get(cohort_name) %in% cohorts_todrop)]
  }

  # cohort-specific estimation
  results_cohort = data.table()
  for(cc in cohorts){
    # set up cohort-specific event times
    times_for_cohort = sort(inputdata[get(cohort_name) == cc, unique(get(time_name))])
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
                            varnames=varnames, control_group = control_group)
      results_cohort = rbindlist(list(results_cohort,res))
    }
  }
  results_cohort=results_cohort[order(Cohort,EventTime)]

  # take the average across cohorts
  results_cohort[, Ntreated_event := sum(Ntreated_post), by="EventTime"]
  results_cohort[, cohort_weights := Ntreated_post/Ntreated_event]
  results_average = results_cohort[, list(ATTe=sum(cohort_weights * ATTge),
                                          ATTe_SE=sqrt(sum(cohort_weights^2 * ATTge_SE^2)),
                                          Etreated_post=sum(cohort_weights*Etreated_post),
                                          Etreated_pre=sum(cohort_weights*Etreated_pre),
                                          Etreated_SE=sqrt(sum(cohort_weights^2*Etreated_SE^2)),
                                          Econtrol_post=sum(cohort_weights*Econtrol_post),
                                          Econtrol_pre=sum(cohort_weights*Econtrol_pre),
                                          Econtrol_SE=sqrt(sum(cohort_weights^2*Econtrol_SE^2)),
                                          Ntreated_post=sum(Ntreated_post),
                                          Ntreated_pre=sum(Ntreated_pre),
                                          Ncontrol_post=sum(Ncontrol_post),
                                          Ncontrol_pre=sum(Ncontrol_pre)
                                          ), list(EventTime,Baseperiod)][order(EventTime,Baseperiod)]

  return(list(results_cohort=results_cohort, results_average=results_average))
}

#' Estimate DiD for all possible cohorts and event time pairs (g,e), as well as the average across cohorts for each event time (e).
#'
#' @param inputdata A data.table.
#' @param varnames A list of the form varnames = list(id_name, time_name, outcome_name, cohort_name), where all four arguments of the list must be a character that corresponds to a variable name in inputdata.
#' @param control_group There are three possibilities: control_group="never-treated" uses the never-treated control group only; control_group="future-treated" uses those units that will receive treatment in the future as the control group; and control_group="all" uses both the never-treated and the future-treated in the control group. Default is control_group="all".
#' @param baseperiod This is the base pre-period that is normalized to zero in the DiD estimation. Default is baseperiod=-1.
#' @param min_event This is the minimum event time (e) to estimate. Default is NULL, in which case, no minimum is imposed.
#' @param max_event This is the maximum event time (e) to estimate. Default is NULL, in which case, no maximum is imposed.
#' @returns A list with two components: results_cohort is a data.table with the DiDge estimates (by event e and cohort g), and results_average is a data.table with the DiDe estimates (by event e, average across cohorts g).
#' @examples
#' # simulate some data
#' simdata = SimDiD(sample_size=1000, ATTcohortdiff = 2)
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
#' @export
DiD <- function(inputdata, varnames, control_group = "all", baseperiod=-1, min_event=NULL, max_event=NULL){

  results = DiDe(inputdata=inputdata, varnames=varnames, control_group=control_group, baseperiod=baseperiod, min_event=min_event, max_event=max_event)

  return(results)
}


