
Errorhandle.mclapply <- function (X, FUN, mc.cores=1, rbind_output=FALSE, print_errors = TRUE,
                                  try_again_if_failure = TRUE, stop_if_failure = FALSE, show_progress = TRUE){
  # define the group split
  suppressWarnings({
    define_groups = split(X, 1:(length(X)/mc.cores))
  })

  # define the progress bar
  if(show_progress){
    # print(sprintf("Beginning parallelization at time %s",proc.time()[3]))
    pb <- progress_bar$new(
      format = "Progress [:bar] Completed tasks: :current/:total (:percent) Elapsed (:elapsed) ETA (:eta)",
      total = length(define_groups), clear = FALSE, width= 120 )
  }

  # define error check
  check_for_errors <- function(res_inner){
    there_was_failure = 0
    for(mmm in 1:length(res_inner)){
      check_val = res_inner[[mmm]]
      if('try-error' %in% class(check_val)){
        there_was_failure = 1
        if(print_errors){
          error_output = list("failed_group"=ggg, "which_item"=mmm, "error_message"=check_val)
          print(error_output)
        }
      }
    }
    return(there_was_failure)
  }

  # execute the loop
  res_outer = NULL
  for(jjj in 1:length(define_groups)){
    # first attempt
    ggg = define_groups[[jjj]]
    res_inner = parallel::mclapply(ggg,FUN,mc.cores=mc.cores)
    there_was_failure = check_for_errors(res_inner)
    # second attempt (if specified)
    if(there_was_failure & try_again_if_failure){
      print("There was a failure. Trying again.")
      res_inner = parallel::mclapply(ggg,FUN,mc.cores=mc.cores)
      there_was_failure = check_for_errors(res_inner)
    }
    # stop if failure (if specified)
    if(stop_if_failure & there_was_failure){
      return(list("failed_group"=ggg, "which_item"=mmm, "error_message"=check_val) )
    }
    # finish
    res_outer = c(res_outer, res_inner)
    if(show_progress){
      pb$tick()
    }
  }


  # try to rbind, with error handling
  if(rbind_output){
    tryCatch(
      {
        res_outer = rbindlist(res_outer)
        return(res_outer)
      },
      error = function(e) {
        return(res_outer)
      }
    )
    return(res_outer)
  }
  return(res_outer)
}
