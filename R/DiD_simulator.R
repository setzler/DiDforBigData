#' DiD data simulator with staggered treatment.
#' @description
#' Simulate data from the model Y_it =  alpha_i + mu_t + ATT*(t >= G_i) + epsilon_it,
#' where i is individual, t is year, and G_i is the cohort.
#' The ATT formula is ATTat0 + EventTime*ATTgrowth + \*cohort_counter\*ATTcohortdiff,
#' where cohort_counter is the order of treated cohort (first, second, etc.).
#' @param seed Set the random seed. Default is seed=1.
#' @param sample_size Number of individuals. Default is sample_size=100.
#' @param cohorts Vector of years at which treatment onset occurs. Default is cohorts=c(2007,2010,2012).
#' @param ATTat0 Treatment effect at event time 0. Default is 1.
#' @param ATTgrowth Increment in the ATT for each event time after 0. Default is 1.
#' @param ATTcohortdiff Incrememnt in the ATT for each cohort. Default is 0.5.
#' @param anticipation Number of years prior to cohort to allow 50% treatment effects. Default is anticipation=0.
#' @param minyear Minimum calendar year to include in the data. Default is minyear=2003.
#' @param maxyear Maximum calendar year to include in the data. Default is maxyear=2013.
#' @param idvar Variance of individual fixed effects (alpha_i). Default is idvar=1.
#' @param yearvar Variance of year effects (mu_i). Default is yearvar=1.
#' @param shockvar Variance of idiosyncratic shocks (epsilon_it). Default is shockvar=1.
#' @param add_covariates Add 2 covariates, called "X1" and "X2". Deault is FALSE.
#' @returns A list with two data.tables.
#' The first data.table is simulated data with variables (id, year, cohort, Y), where Y is the outcome variable.
#' The second data.table contains the true ATT values, both at the (event,cohort) level and by event averaging across cohorts.
#' @examples
#' # simulate data with default options
#' SimDiD()
#'
#' # change the random seed
#' SimDiD(seed=123)
#'
#' # increase sample size
#' SimDiD(sample_size=200)
#'
#' # increase the initial ATT for the first treated cohort
#' SimDiD(ATTat0=2)
#'
#' # change cross-time and cross-cohort ATT variation, where the formula is ATTat0 + EventTime*ATTgrowth
#' SimDiD(ATTgrowth=2,ATTcohortdiff=3)
#'
#' # change the range of years considered and treatment cohorts
#' SimDiD(minyear=1994, maxyear=1999, cohorts=c(1995,1996))
#'
#' # add two periods of anticipation prior to treatment. The effect will be 0.5*ATTat0.
#' SimDiD(anticipation=2)
#' @export
SimDiD <- function(seed=1,sample_size=100, cohorts=c(2007,2010,2012), ATTat0=1, ATTgrowth=1, ATTcohortdiff=0.5, anticipation=0, minyear=2003, maxyear=2013, idvar=1, yearvar=1, shockvar=1, add_covariates=FALSE){
  # seed=1; sample_size=1000; cohorts=c(2007,2010,2012); ATTat0=1; ATTgrowth=1; ATTcohortdiff=0.5; anticipation=0; minyear=2003; maxyear=2013; idvar=1; yearvar=1; shockvar=1; add_covariates=FALSE
  set.seed(seed)

  # create id-by-year data
  simdata = setDT(expand.grid(id=1:sample_size,year=minyear:maxyear))

  # simulate unobservables
  simdata[, shock := rnorm(nrow(simdata),sd=shockvar)]
  simdata[, individualFE := rnorm(1,sd=idvar), id]
  simdata[, yearFE := rnorm(1,sd=yearvar), year]
  simdata[, individualFE_ecdf := ecdf(individualFE)(individualFE)]

  # simulate treatment cohort
  simcohorts = c(cohorts,Inf)
  cutoffs = rep(1,length(simcohorts))/length(simcohorts)
  cutoffs = c(0,cumsum(cutoffs))
  for(ii in 1:(length(cutoffs)-1)){
    simdata[individualFE_ecdf >= cutoffs[ii] & individualFE_ecdf <= cutoffs[ii+1], cohort := simcohorts[ii]]
  }
  simdata[, event := year - cohort]

  # simulate ATT
  # simdata[event < -1*anticipation, event := 0]
  simdata[, cohort_counter := 0.0]
  for(ii in 1:length(cohorts)){
    simdata[cohort==cohorts[ii], cohort_counter := (ii-1)]
  }
  simdata[, ATT := 0.0]
  simdata[year >= cohort, ATT := (ATTat0 + event*ATTgrowth + cohort_counter*ATTcohortdiff)]
  simdata[year < cohort & year >= (cohort-anticipation), ATT := 0.5*ATTat0]
  true_ATTge = simdata[(year >= (cohort-anticipation)),list(ATTge = mean(ATT)), list(cohort, event)][order(cohort,event)][, cohort := as.character(cohort)]
  true_ATTe = simdata[(year >= (cohort-anticipation)),list(cohort="Average",ATTge = mean(ATT)), list(event)][order(event)]
  true_ATT = rbindlist(list(true_ATTge,true_ATTe),use.names=T)

  # simulate outcome
  keep_vars = c("id","year","cohort","Y")
  simdata[, Y := 10 + individualFE + yearFE + ATT + shock]
  if(add_covariates){
    simdata[, X1 := year - mean(year) + rnorm(nrow(simdata))] # this causes no issues, since it differs out on average
    simdata[, X2 := 0.0]
    simdata[event > 0, X2 := event]
    simdata[, X2 := log(1 + cohort_counter + X2) + rnorm(nrow(simdata))] # this causes problems because there is an event-cohort interaction
    simdata[, covariate_term := -1.0*X1 + 1.0*X2]
    simdata[, Y := Y + covariate_term]
    keep_vars = c(keep_vars, "X1", "X2")
  }
  simdata = simdata[order(id,year)]
  simdata = simdata[,.SD,.SDcols=keep_vars]
  return(list(simdata=simdata,true_ATT=true_ATT))
}
