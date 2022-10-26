
#' Compare speed and memory efficiency of various DiD estimators available in R.
#' @description
#' This function compares speed and memory efficiency of various DiD estimators available in R:
#' - `DiDforBigData`: This package's DiDforBigData::DiD()
#' - `saturatedOLS`: This package's fully-interacted OLS (internal)
#' - `ES`: Novgorodsky and Setzler's Stacked Event-study eventStudy::ES()
#' - `CSreg`: Callaway and Sant'Anna's reg estimator did::att_gt(..., est_method = "reg")
#' - `CSdr`: Callaway and Sant'Anna's dr estimator did::att_gt(..., est_method = "dr")
#' - `CH`: DIDmultiplegt::did_multiplegt()
#'
#' @export
EfficiencyTestDiD <- function(estimators = c("DiDforBigData","CSreg","CSdr","CH"), sample_sizes=c(1e3,1e4,1e5), reps=3){
  library(profmem)
  library(did)
  library(eventStudy)
  library(DIDmultiplegt)

  deChaisemartin <- function(inputdata,varname_options,dynamic=5,placebo=3,brep=2){
    # set up variable names
    time_name = varname_options$time_name
    outcome_name = varname_options$outcome_name
    cohort_name = varname_options$cohort_name
    id_name = varname_options$id_name
    #
    inputdata[, D := as.numeric(year >= get(cohort_name))]
    CH_res = did_multiplegt(df=copy(inputdata), Y=outcome_name, T = time_name, G = id_name, D = "D", dynamic=dynamic, placebo=3, brep=brep, cluster=id_name)
    if(dynamic>0){
      CH_res_dyn = data.table( EventTime=0, ATTe=CH_res$effect, ATTe_SE=CH_res$se_effect)
      for(ee in 1:dynamic){
        CH_res_dyn = rbindlist(list(CH_res_dyn, data.table( EventTime=ee, ATTe=CH_res[[paste0("dynamic_",ee)]], ATTe_SE=CH_res[[paste0("se_dynamic_",ee)]])))
      }
      if(placebo>0){
        for(ee in 1:placebo){
          CH_res_dyn = rbindlist(list(CH_res_dyn, data.table( EventTime=(-1*ee), ATTe=CH_res[[paste0("placebo_",ee)]], ATTe_SE=CH_res[[paste0("se_placebo_",ee)]])))
        }
      }
      CH_res_dyn = CH_res_dyn[order(EventTime)]
      return(CH_res_dyn)
    }
  }

  stackedES <- function(inputdata, varname_options){
    # set up variable names
    time_name = varname_options$time_name
    outcome_name = varname_options$outcome_name
    cohort_name = varname_options$cohort_name
    id_name = varname_options$id_name
    # run the event study
    inputdata[, (cohort_name) := as.integer(get(cohort_name))]
    inputdata[is.na(get(cohort_name)), (cohort_name) := 1e9] # ES() doesn't like Inf
    ES_res <- ES(long_data=inputdata, outcomevar=outcome_name,
                  unit_var=id_name, cal_time_var=time_name,
                  onset_time_var=cohort_name, cluster_vars=id_name)
    return(ES_res)
  }

  # compare to Callaway & Sant'Anna
  CallawaySantanna <- function(inputdata, varname_options, est_method="reg"){
    # set up variable names
    time_name = varname_options$time_name
    outcome_name = varname_options$outcome_name
    cohort_name = varname_options$cohort_name
    id_name = varname_options$id_name
    # run Callaway-Sant'Anna estimator
    callawaysantanna =  att_gt(yname = outcome_name,
                               gname = cohort_name,
                               idname = id_name,
                               tname = time_name,
                               xformla = ~1,
                               data = inputdata,
                               est_method = est_method,
                               control_group = "notyettreated",
                               base_period="universal",
                               anticipation=0)
    CS = data.table(CalendarTime=callawaysantanna$t, Cohort=callawaysantanna$group, ATTge=callawaysantanna$att, ATTge_SE=callawaysantanna$se)
    CS[, EventTime := CalendarTime - Cohort]
    CS_ES = aggte(callawaysantanna,type="dynamic")
    CS_ES = data.table(EventTime=CS_ES$egt, ATTe=CS_ES$att.egt, ATTe_SE=CS_ES$se.egt)
    return(list(results_cohort=CS, results_average=CS_ES))
  }


  # Compare to OLS
  saturatedOLS <- function(inputdata){
    # set up variable names
    time_name = varname_options$time_name
    outcome_name = varname_options$outcome_name
    cohort_name = varname_options$cohort_name
    id_name = varname_options$id_name
    # rename
    setnames(inputdata,cohort_name,"cohort")
    setnames(inputdata,time_name,"year")
    setnames(inputdata,outcome_name,"Y")
    setnames(inputdata,id_name,"id")
    # cohort and year sets
    Gset = inputdata[,unique(cohort)]
    Gset = Gset[is.finite(Gset)]
    Gset = sort(Gset)
    Tset = inputdata[,unique(year)]
    Tset = sort(Tset)
    # construct the cohort dummies and cohort-time interactions
    cohort_names = c()
    int_names = c()
    for(GG in Gset){
      for(TT in Tset){
        if(!(TT - GG == -2)){
          int_names = c(int_names, paste0("intTT",TT,"GG",GG))
          inputdata[, (paste0("intTT",TT,"GG",GG)) := as.numeric(year==TT & cohort==GG)]
        }
      }
    }
    for(GG in Gset){
      cohort_names = c(cohort_names, paste0("cohortGG",GG))
      inputdata[, (paste0("cohortGG",GG)) := as.numeric(cohort==GG)]
    }
    # run OLS
    DiD_formula = paste0("Y ~ -1 + factor(year) + ",paste0(cohort_names, collapse=" + ")," + ",paste0(int_names, collapse=" + "))
    DiD_model = lm(DiD_formula, data=inputdata)
    # extract coefficients
    OLS_mat = data.table(summary(DiD_model)$coefficients)
    OLS_mat$rownames = rownames(summary(DiD_model)$coefficients)
    OLS_mat = OLS_mat[14:nrow(OLS_mat)]
    OLS_mat[, CalendarTime := substr(rownames, 6,9)]
    OLS_mat[, Cohort := substr(rownames, 12,15)]
    setnames(OLS_mat,c("Estimate","Std. Error"),c("OLS","OLS_SE"))
    OLS_mat = OLS_mat[,list(CalendarTime=as.integer(CalendarTime),Cohort=as.integer(Cohort),OLS, OLS_SE)]
    return(OLS_mat)

  }


  # prepare variable names
  varname_options = list()
  varname_options$time_name = "year"
  varname_options$outcome_name = "Y"
  varname_options$cohort_name = "cohort"
  varname_options$id_name = "id"
  # define speed tester
  speed_tester <- function(sample_size,reps){
    all_res = data.table()
    for(seed in 1:reps){
      this_res = data.table()
      # simulate
      inputdata = SimDiD(sample_size = sample_size, seed = seed, minyear=2004, maxyear=2013)
      # my package
      if("DiDforBigData" %in% estimators){
        time0 = proc.time()[3]
        My_p <- profmem({ My_res <- DiD(copy(inputdata), varname_options = varname_options, min_event = -3, max_event = 5) })
        time1 = proc.time()[3]
        My_time = (time1 - time0)/60
        My_mem = sum(My_p$bytes,na.rm=T)/1e9
        My_ATT = My_res$results_average[EventTime==1]$ATTe
        My_ATTse = My_res$results_average[EventTime==1]$ATTe_SE
        this_res = rbindlist(list(this_res, data.table(method="DiDforBigData", ATTe=My_ATT,ATTse=My_ATTse)))
        gc()
        print(sprintf("finished DiDforBigData in %s",My_time))
      }
      # OLS
      if("saturatedOLS" %in% estimators){
        time0 = proc.time()[3]
        OLS_p <- profmem({ OLS_res <- saturatedOLS(copy(inputdata)) })
        time1 = proc.time()[3]
        OLS_time = (time1 - time0)/60
        OLS_mem = sum(OLS_p$bytes,na.rm=T)/1e9
        OLS_ATT = OLS_res[CalendarTime==2007 & Cohort==2007]$OLS
        OLS_ATTse = OLS_res[CalendarTime==2007 & Cohort==2007]$OLS_SE
        this_res = rbindlist(list(this_res, data.table(method="saturedOLS", ATTe=OLS_ATT,ATTse=OLS_ATTse)))
        gc()
        print(sprintf("finished saturedOLS in %s",OLS_time))
      }
      # Novgorodsky Setzler eventStudy
      if("ES" %in% estimators){
        time0 = proc.time()[3]
        ES_p <- profmem({ ES_res <- stackedES(copy(inputdata), varname_options) })
        time1 = proc.time()[3]
        ES_time = (time1 - time0)/60
        ES_mem = sum(ES_p$bytes,na.rm=T)/1e9
        ES_ATT = ES_res[ref_onset_time==2007 & ref_event_time==0 & rn=="catt"]$estimate
        ES_ATTse = ES_res[ref_onset_time==2007 & ref_event_time==0 & rn=="catt"]$cluster_se
        this_res = rbindlist(list(this_res, data.table(method="ES", ATTe=ES_ATT,ATTse=ES_ATTse)))
        gc()
        print(sprintf("finished ES in %s",ES_time))
      }
      # Callaway Sant'Anna, reg
      if("CSreg" %in% estimators){
        time0 = proc.time()[3]
        CSreg_p <- profmem({ CSreg_res <- CallawaySantanna(copy(inputdata), varname_options, est_method = "reg") })
        time1 = proc.time()[3]
        CSreg_time = (time1 - time0)/60
        CSreg_mem = sum(CSreg_p$bytes,na.rm=T)/1e9
        CSreg_ATT = CSreg_res$results_average[EventTime==1]$ATTe
        CSreg_ATTse = CSreg_res$results_average[EventTime==1]$ATTe_SE
        this_res = rbindlist(list(this_res, data.table(method="CSreg", ATTe=CSreg_ATT, ATTse=CSreg_ATTse)))
        gc()
        print(sprintf("finished CSreg in %s",CSreg_time))
      }
      # Callaway Sant'Anna, dr
      if("CSdr" %in% estimators){
        time0 = proc.time()[3]
        CSdr_p <- profmem({ CSdr_res <- CallawaySantanna(copy(inputdata), varname_options, est_method = "dr") })
        time1 = proc.time()[3]
        CSdr_time = (time1 - time0)/60
        CSdr_mem = sum(CSdr_p$bytes,na.rm=T)/1e9
        CSdr_ATT = CSdr_res$results_average[EventTime==1]$ATTe
        CSdr_ATTse = CSdr_res$results_average[EventTime==1]$ATTe_SE
        this_res = rbindlist(list(this_res, data.table(method="CSdr", ATTe=CSdr_ATT, ATTse=CSdr_ATTse)))
        gc()
        print(sprintf("finished CSdr in %s",CSdr_time))
      }
      # de Chaisemartin & D'Haultfoeuille
      if("CH" %in% estimators){
        time0 = proc.time()[3]
        CH_p <- profmem({ CH_res <- deChaisemartin(copy(inputdata),varname_options,dynamic=5,placebo=3) })
        time1 = proc.time()[3]
        CH_time = (time1 - time0)/60
        CH_mem = sum(CH_p$bytes,na.rm=T)/1e9
        CH_ATT = CH_res[EventTime==1]$ATTe
        CH_ATTse = CH_res[EventTime==1]$ATTe_SE
        this_res = rbindlist(list(this_res, data.table(method="CH", ATTe=CH_ATT, ATTse=CH_ATTse)))
        gc()
        print(sprintf("finished CH in %s",CH_time))
      }
      # collect
      all_res = rbindlist(list(all_res,this_res))
    }
    if(reps>1){
      all_res = all_res[,list(ATTe=median(ATTe),ATTse=median(ATTse)),method]
    }
    all_res$sample_size = sample_size
    return(all_res)
  }

  all_speeds = data.table()
  for(nn in sample_sizes){
    speeds = speed_tester(sample_size=nn,reps=reps)
    all_speeds = rbindlist(list(all_speeds, speeds))
    gc()
    print(all_speeds[])
    print(sprintf("sample %s done",nn))
  }

  write.csv(all_speeds, file="inst/speed_test.csv", row.names = FALSE)

  return(all_speeds)

#
#   all_speeds_long = melt(all_speeds, id.vars = "sample_size")
#   all_speeds_long[variable=="OLS_time", variable := "OLS"]
#   all_speeds_long[variable=="CS_time", variable := "Callaway & Sant'Anna"]
#   all_speeds_long[variable=="My_time", variable := "My package"]
#   all_speeds_long[, value := value/60]
#
#   gg = ggplot(aes(x=factor(sample_size), y=value, fill=variable),data=all_speeds_long)+
#     theme_bw(base_size=22) + theme(legend.position = "bottom") +
#     labs(x="Sample Size",y="Speed (Minutes)",fill="") +
#     scale_y_continuous(breaks=pretty_breaks()) +
#     geom_bar(stat='identity', position='dodge') +
#     scale_fill_manual(values=c('blue','black','green'))
#   ggsave(gg,filename="inst/speed_levels.pdf",width=10,height=5)
#
#   gg = ggplot(aes(x=factor(sample_size), y=log(value), fill=variable),data=all_speeds_long)+
#     theme_bw(base_size=18) + theme(legend.position = "bottom") +
#     labs(x="Sample Size",y="Speed (Log Minutes)",fill="") +
#     scale_y_continuous(breaks=pretty_breaks()) +
#     geom_bar(stat='identity', position='dodge') +
#     scale_fill_manual(values=c('blue','black','green'))
#   ggsave(gg,filename="inst/speed_logs.pdf",width=8,height=5)



}


