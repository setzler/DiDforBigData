
PlotDiDgt <- function(DiD_results, ci_factor = 1.96, max_jitter=0.05, lower_event=NULL, upper_event=NULL, colors=NULL){
  # ci_factor = 1.96; max_jitter=0.05; lower_event=-5; upper_event=5

  # read the data
  if(is.list(DiD_results)){
    figdata_temp = copy(DiD_results$results_cohort)
  }
  if(is.data.table(DiD_results)){
    figdata_temp = copy(DiD_results)
  }

  # set up the time variables
  all_calendartimes = figdata_temp[, sort(unique(CalendarTime))]
  figdata_temp[, CalendarTime := as.numeric(CalendarTime)]
  figdata_temp[, EventTime := CalendarTime - Cohort]
  if(!is.null(lower_event)){
    figdata_temp = figdata_temp[EventTime >= lower_event]
  }
  if(!is.null(upper_event)){
    figdata_temp = figdata_temp[EventTime <= upper_event]
  }
  figdata_temp = figdata_temp[order(CalendarTime,Cohort)]

  # set up jitter
  cohorts = figdata_temp[, unique(Cohort) ]
  jitter_vals = seq(-max_jitter, max_jitter, length.out=length(cohorts))
  for(ii in 1:length(jitter_vals)){
    figdata_temp[Cohort==cohorts[ii], jitter_val := jitter_vals[ii]]
  }
  figdata_temp[, CalendarTime_jitter := CalendarTime + jitter_val ]

  # make ggplot
  gg <- ggplot(aes( x = CalendarTime_jitter, y = ATTge, colour = factor(Cohort)), data = figdata_temp) +
    labs(x = "Calendar Time", color = "Cohort", y ="") +
    geom_point() + theme_bw(base_size = 20) +
    geom_errorbar(aes(ymin = ATTge - ci_factor * ATTge_SE, ymax = ATTge + ci_factor * ATTge_SE), width=0.1) +
    scale_x_continuous(breaks = all_calendartimes) +
    scale_y_continuous(breaks = pretty_breaks()) +
    geom_hline(yintercept=0, color="black", linetype="dashed") +
    theme(legend.position = "bottom")

  # set colors manually
  if(!is.null(colors)){
    if(length(cohorts) == length(colors)){
      gg = gg + scale_color_manual(values=colors)
    }
    if(length(cohorts) == length(colors)){
      warning(sprintf("The provided colors could not be used, as there are %s colors and %s cohorts in the results.", length(cohorts), length(colors)))
    }
  }

  return(gg)

}

PlotDiDge <- function(DiD_results, ci_factor = 1.96, max_jitter=0.05, lower_event=NULL, upper_event=NULL, colors=NULL){
  # ci_factor = 1.96; max_jitter=0.05; lower_event=-5; upper_event=5

  # read the data
  if(is.list(DiD_results)){
    figdata_temp = copy(DiD_results$results_cohort)
  }
  if(is.data.table(DiD_results)){
    figdata_temp = copy(DiD_results)
  }

  # set up the time variables
  all_events = figdata_temp[, sort(unique(EventTime))]
  omitted_event_time = figdata_temp[,sort(unique(Baseperiod))]
  figdata_temp = figdata_temp[EventTime != Baseperiod]
  figdata_temp[, EventTime := as.numeric(EventTime)]
  if(!is.null(lower_event)){
    figdata_temp = figdata_temp[EventTime >= lower_event]
  }
  if(!is.null(upper_event)){
    figdata_temp = figdata_temp[EventTime <= upper_event]
  }
  figdata_temp = figdata_temp[order(EventTime,Cohort)]

  # set up jitter
  cohorts = figdata_temp[, unique(Cohort) ]
  jitter_vals = seq(-max_jitter, max_jitter, length.out=length(cohorts))
  for(ii in 1:length(jitter_vals)){
    figdata_temp[Cohort==cohorts[ii], jitter_val := jitter_vals[ii]]
  }
  figdata_temp[, EventTime_jitter := EventTime + jitter_val ]

  # make ggplot
  gg <- ggplot(aes( x = EventTime_jitter, y = ATTge, colour = factor(Cohort)), data = figdata_temp) +
    labs(x = "Event Time", color = "Cohort", y ="") +
    geom_point() + theme_bw(base_size = 20) +
    geom_errorbar(aes(ymin = ATTge - ci_factor * ATTge_SE, ymax = ATTge + ci_factor * ATTge_SE), width=0.1) +
    scale_x_continuous(breaks = all_events) +
    scale_y_continuous(breaks = pretty_breaks()) +
    geom_hline(yintercept=0, color="black", linetype="dashed") +
    annotate(geom = "text", x = omitted_event_time, y = 0, label = "Normalized\nto Zero", size = 4) +
    theme(legend.position = "bottom")

  # set colors manually
  if(!is.null(colors)){
    if(length(cohorts) == length(colors)){
      gg = gg + scale_color_manual(values=colors)
    }
    if(length(cohorts) == length(colors)){
      warning(sprintf("The provided colors could not be used, as there are %s colors and %s cohorts in the results.", length(cohorts), length(colors)))
    }
  }

  return(gg)

}


PlotDiDe <- function(DiD_results, ci_factor = 1.96, max_jitter=0.05, lower_event=NULL, upper_event=NULL, colors=NULL){
  # ci_factor = 1.96; max_jitter=0.05; lower_event=-5; upper_event=5

  # prepare the data
  if(is.list(DiD_results)){
    figdata_temp = copy(DiD_results$results_average)
  }
  if(is.data.table(DiD_results)){
    figdata_temp = copy(DiD_results)
  }

  # set up the time variables
  all_events = figdata_temp[, sort(unique(EventTime))]
  omitted_event_time = figdata_temp[,sort(unique(Baseperiod))]
  figdata_temp = figdata_temp[EventTime != Baseperiod]
  figdata_temp[, EventTime := as.numeric(EventTime)]
  if(!is.null(lower_event)){
    figdata_temp = figdata_temp[EventTime >= lower_event]
  }
  if(!is.null(upper_event)){
    figdata_temp = figdata_temp[EventTime <= upper_event]
  }
  figdata_temp = figdata_temp[order(EventTime)]

  # make ggplot
  gg <- ggplot(aes( x = EventTime, y = ATTe), data = figdata_temp) +
    labs(x = "Event Time", y ="") +
    geom_point() + theme_bw(base_size = 20) +
    geom_errorbar(aes(ymin = ATTe - ci_factor * ATTe_SE, ymax = ATTe + ci_factor * ATTe_SE), width=0.1) +
    scale_x_continuous(breaks = all_events) +
    scale_y_continuous(breaks = pretty_breaks()) +
    geom_hline(yintercept=0, color="black", linetype="dashed") +
    annotate(geom = "text", x = omitted_event_time, y = 0, label = "Normalized\nto Zero", size = 4) +
    theme(legend.position = "bottom")

  return(gg)

}

#' Plot the results from the DiD estimation.
#'
#' @param DiD_results Output from the DiD function.
#' @param plot_type There are 3 options: default option is plot_type="ge" plots separate points for each cohort (g) and event time (e); option plot_type="gt" plots separate points for each cohort (g) and calendar time (t); and option plot_type="e" plots the average across cohorts for each event time (e).
#' @param ci_factor This is the multiplier that converts a standard errors (SE) into a confidence interval (CI). That is, CI = estimate +/- ci_factor*SE. Default option is ci_factor = 1.96, which corresponds to a 95% CI. A natural alternative is ci_factor = 1.645, which corresponds to a 90% CI.
#' @param max_jitter This jitters the points so that they do not overlap. If plot_type="ge" or plot_type="gt", there will be multiple points overlapping due to multiple cohorts, and max_jitter>0 will dodge the points by at most max_jitter. Default option is max_jitter=0.05.
#' @param lower_event Minimum event time (e) to show in the plot. Default option is NULL.
#' @param upper_event Maximum event time (e) to show in the plot. Default option is NULL.
#' @param colors Manually set the colors for each cohort. Only relevant if plot_type="ge" or plot_type="gt". Default option is NULL.
#' @returns A ggplot object.
#' @export
PlotDiD <- function(DiD_results, plot_type="ge", ci_factor = 1.96, max_jitter=0.05, lower_event=NULL, upper_event=NULL, colors=NULL){
  gg = NULL
  if(plot_type=="gt"){
    gg = PlotDiDgt(DiD_results=DiD_results, lower_event=lower_event, upper_event=upper_event, ci_factor = ci_factor, max_jitter=max_jitter, colors=colors)
  }
  if(plot_type=="ge"){
    gg = PlotDiDge(DiD_results=DiD_results, lower_event=lower_event, upper_event=upper_event, ci_factor = ci_factor, max_jitter=max_jitter, colors=colors)
  }
  if(plot_type=="e"){
    gg = PlotDiDe(DiD_results=DiD_results, lower_event=lower_event, upper_event=upper_event, ci_factor = ci_factor, max_jitter=max_jitter, colors=colors)
  }
  return(gg)
}





PlotESge <- function(DiD_results, cohort=2007, ci_factor = 1.96, max_jitter=0.05, lower_event=-3, upper_event=5, colors=NULL){
  # cohort=2007; ci_factor = 1.96; max_jitter=0.05; lower_event=-5; upper_event=5

  # prepare the data
  if(is.list(DiD_results)){
    figdata_temp = copy(DiD_results$results_cohort)
  }
  if(is.data.table(DiD_results)){
    figdata_temp = copy(DiD_results)
  }
  if(!(cohort %in% figdata_temp[,unique(Cohort)])){
    stop(sprintf("You provided cohort=%s, which is not a valid cohort in the DiD_results.", cohort))
  }

  # set up the time variables
  all_calendartimes = figdata_temp[, sort(unique(CalendarTime))]
  figdata_temp[, CalendarTime := as.numeric(CalendarTime)]
  figdata_temp[, EventTime := CalendarTime - Cohort]
  figdata_temp = figdata_temp[Cohort==cohort]
  figdata_temp[, EventTime := as.numeric(EventTime)]
  if(!is.null(lower_event)){
    figdata_temp = figdata_temp[EventTime >= lower_event]
  }
  if(!is.null(upper_event)){
    figdata_temp = figdata_temp[EventTime <= upper_event]
  }

  # reshape data for treatment and control groups
  figdata_temp = figdata_temp[,.(CalendarTime, Econtrol_post, Econtrol_SE, Etreated_post, Etreated_SE)]
  figdata_temp1 = melt(figdata_temp[,.(CalendarTime, Econtrol_post, Etreated_post)], id.vars="CalendarTime")
  figdata_temp1[, Group := "Control"]
  figdata_temp1[variable=="Etreated_post", Group := "Treated"]
  figdata_temp2 = melt(figdata_temp[,.(CalendarTime, Econtrol_SE, Etreated_SE)], id.vars="CalendarTime")
  figdata_temp2[, Group := "Control"]
  figdata_temp2[variable=="Etreated_SE", Group := "Treated"]
  figdata_temp1 = figdata_temp1[,list(CalendarTime, Group, Mean=value)]
  figdata_temp2 = figdata_temp2[,list(CalendarTime, Group, SE=value)]
  figdata_temp3 = merge(figdata_temp1, figdata_temp2, by=c("CalendarTime","Group"))

  # set up jitter
  groups = figdata_temp3[, unique(Group) ]
  jitter_vals = seq(-max_jitter, max_jitter, length.out=length(groups))
  for(ii in 1:length(jitter_vals)){
    figdata_temp3[Group==groups[ii], jitter_val := jitter_vals[ii]]
  }
  figdata_temp3[, CalendarTime_jitter := CalendarTime + jitter_val ]

  # make ggplot
  gg <- ggplot(aes( x = CalendarTime_jitter, y = Mean, colour = factor(Group)), data = figdata_temp3) +
    labs(x = "Calendar Time", color = "", y ="") +
    geom_point() + theme_bw(base_size = 20) +
    geom_errorbar(aes(ymin = Mean - ci_factor * SE, ymax = Mean + ci_factor * SE), width=0.1) +
    scale_x_continuous(breaks = all_calendartimes) +
    scale_y_continuous(breaks = pretty_breaks()) +
    geom_vline(xintercept=cohort, color="black", linetype="dashed") +
    theme(legend.position = "bottom")

  if(!is.null(colors)){
    gg = gg + scale_color_manual(values=colors)
  }

  return(gg)

}



PlotESe <- function(DiD_results, ci_factor = 1.96, max_jitter=0.05, lower_event=-3, upper_event=5, colors=NULL){
  # ci_factor = 1.96; max_jitter=0.05; lower_event=-3; upper_event=5

  # prepare the data
  if(is.list(DiD_results)){
    figdata_temp = copy(DiD_results$results_average)
  }
  if(is.data.table(DiD_results)){
    figdata_temp = copy(DiD_results)
  }

  # set up the time variables
  all_events = figdata_temp[, sort(unique(EventTime))]
  figdata_temp[, EventTime := as.numeric(EventTime)]
  if(!is.null(lower_event)){
    figdata_temp = figdata_temp[EventTime >= lower_event]
  }
  if(!is.null(upper_event)){
    figdata_temp = figdata_temp[EventTime <= upper_event]
  }
  figdata_temp = figdata_temp[order(EventTime)]

  # the treatment group composition changes across event times, so we renormalize the levels to the omitted event time
  omitted_event_time = figdata_temp[,sort(unique(Baseperiod))]
  figdata_temp[, Econtrol_pre_init := Econtrol_pre[EventTime==omitted_event_time]]
  figdata_temp[, Etreated_pre_init := Etreated_pre[EventTime==omitted_event_time]]
  figdata_temp[, Econtrol_pre_diff := Econtrol_pre - Econtrol_pre_init ]
  figdata_temp[, Etreated_pre_diff := Etreated_pre - Etreated_pre_init ]
  figdata_temp[, Econtrol_pre := Econtrol_pre - Econtrol_pre_diff]
  figdata_temp[, Etreated_pre := Etreated_pre - Etreated_pre_diff]
  figdata_temp[, Econtrol_post := Econtrol_post - Econtrol_pre_diff]
  figdata_temp[, Etreated_post := Etreated_post - Etreated_pre_diff]

  # reshape data for treatment and control groups
  figdata_temp = figdata_temp[,.(EventTime, Econtrol_post, Econtrol_pre, Econtrol_SE, Etreated_post, Etreated_pre, Etreated_SE)]
  figdata_temp1 = melt(figdata_temp[,.(EventTime, Econtrol_post, Etreated_post)], id.vars="EventTime")
  figdata_temp1[, Group := "Control"]
  figdata_temp1[variable=="Etreated_post", Group := "Treated"]
  figdata_temp2 = melt(figdata_temp[,.(EventTime, Econtrol_SE, Etreated_SE)], id.vars="EventTime")
  figdata_temp2[, Group := "Control"]
  figdata_temp2[variable=="Etreated_SE", Group := "Treated"]
  figdata_temp1 = figdata_temp1[,list(EventTime, Group, Mean=value)]
  figdata_temp2 = figdata_temp2[,list(EventTime, Group, SE=value)]
  figdata_temp3 = merge(figdata_temp1, figdata_temp2, by=c("EventTime","Group"))

  # set up jitter
  groups = figdata_temp3[, unique(Group) ]
  jitter_vals = seq(-max_jitter, max_jitter, length.out=length(groups))
  for(ii in 1:length(jitter_vals)){
    figdata_temp3[Group==groups[ii], jitter_val := jitter_vals[ii]]
  }
  figdata_temp3[, EventTime_jitter := EventTime + jitter_val ]

  # make ggplot
  gg <- ggplot(aes( x = EventTime_jitter, y = Mean, colour = factor(Group)), data = figdata_temp3) +
    labs(x = "Event Time", color = "", y ="") +
    geom_point() + theme_bw(base_size = 20) +
    geom_errorbar(aes(ymin = Mean - ci_factor * SE, ymax = Mean + ci_factor * SE), width=0.1) +
    scale_x_continuous(breaks = all_events) +
    scale_y_continuous(breaks = pretty_breaks()) +
    geom_vline(xintercept=0, color="black", linetype="dashed") +
    theme(legend.position = "bottom")

  if(!is.null(colors)){
    gg = gg + scale_color_manual(values=colors)
  }

  return(gg)

}


#' Plot an event study in the original units.
#'
#' @param DiD_results Output from the DiD function.
#' @param plot_type There are 2 options: default option is plot_type="ge" plots all the event time (e) for a given cohort (g); and option plot_type="e" plots the average across cohorts for each event time (e).
#' @param cohort Select one cohort to plot if plot_type="ge".
#' @param ci_factor This is the multiplier that converts a standard errors (SE) into a confidence interval (CI). That is, CI = estimate +/- ci_factor*SE. Default option is ci_factor = 1.96, which corresponds to a 95% CI. A natural alternative is ci_factor = 1.645, which corresponds to a 90% CI.
#' @param max_jitter This jitters the points so that they do not overlap. If plot_type="ge" or plot_type="gt", there will be multiple points overlapping due to multiple cohorts, and max_jitter>0 will dodge the points by at most max_jitter. Default option is max_jitter=0.05.
#' @param lower_event Minimum event time (e) to show in the plot. Default option is NULL.
#' @param upper_event Maximum event time (e) to show in the plot. Default option is NULL.
#' @param colors Manually set the colors for each cohort. Only relevant if plot_type="ge" or plot_type="gt". Default option is NULL.
#' @returns A ggplot object.
#' @export
PlotES <- function(DiD_results, plot_type, cohort=NULL, ci_factor = 1.96, max_jitter=0.05, lower_event=-3, upper_event=5, colors=NULL){
  gg = NULL
  if(plot_type=="ge"){
    gg = PlotESge(DiD_results=DiD_results, cohort=cohort, lower_event=lower_event, upper_event=upper_event, ci_factor = ci_factor, max_jitter=max_jitter, colors=colors)
  }
  if(plot_type=="e"){
    gg = PlotESe(DiD_results=DiD_results, lower_event=lower_event, upper_event=upper_event, ci_factor = ci_factor, max_jitter=max_jitter, colors=colors)
  }
  return(gg)
}

