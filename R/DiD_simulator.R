#' DiD data simulator with staggered treatment.
#' @description
#' Simulate data from the model Y_it =  alpha_i + mu_t + ATT*(t >= G_i) + epsilon_it, where i is individual, t is year, and G_i is the cohort.
#' @param seed Set the random seed. Default is seed=1.
#' @param sample_size Number of individuals. Default is sample_size=100.
#' @param cohorts Vector of years at which treatment onset occurs. Default is cohorts=c(2007,2010,2012).
#' @param ATTs Vector of treatment effects. Must be same length as cohorts. Default is ATTs=c(1,1,1).
#' @param anticipation Number of years prior to cohort to allow 50% treatment effects. Default is anticipation=0.
#' @param minyear Minimum calendar year to include in the data. Default is minyear=2003.
#' @param maxyear Maximum calendar year to include in the data. Default is maxyear=2013.
#' @param idvar Variance of individual fixed effects (alpha_i). Default is idvar=1.
#' @param yearvar Variance of year effects (mu_i). Default is yearvar=1.
#' @param shockvar Variance of idiosyncratic shocks (epsilon_it). Default is shockvar=1.
#' @returns A data.table with variables (id, year, cohort, Y), where Y is the outcome variable.
#' @examples
#' SimDiD()
#' @export
SimDiD <- function(seed=1,sample_size=100, cohorts=c(2007,2010,2012), ATTs=c(1,1,1), anticipation=0, minyear=2003, maxyear=2013, idvar=1, yearvar=1, shockvar=1){
  # seed=1; sample_size=1000; cohorts=c(2007,2010,2012); ATTs=c(1,1,1); anticipation=0; minyear=2003; maxyear=2013; idvar=1; yearvar=1; shockvar=1
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
  ATTdata = data.table(cohort=c(cohorts,Inf), ATT=c(ATTs,0))
  simdata = merge(simdata,ATTdata,by="cohort")

  # simulate outcome
  simdata[, Y := 10 + individualFE + yearFE + ATT*(year >= cohort) + ATT*0.5*(year < cohort & year >= (cohort-anticipation)) + shock]
  simdata = simdata[order(id,year)]
  simdata = simdata[,.(id,year,cohort,Y)]
  return(simdata)
}
