#-------------------------------------------------------------------------------
#
#   Null model standardized effect size and diagnostic function
#
#-------------------------------------------------------------------------------

# Author: William K. Annis

# Created: 03/16/2026


#' Null model standardized effect size and diagnostic function
#'
#' @description Estimates null model standardized effects sizes and returns
#' diagnostics statistics for null distributions.
#'
#' @param obs vector, dist, or matrix containing observed diversity values
#' @param null_list list including null iterations of diversity metric. Must be
#' in same format as obs
#' @param emp_padding value used to pad empirical p values to prevent inf
#' effect sizes. If p value = zero, effect is infinite. This can be corrected
#' by adding 1 to the numerator and denominator of p value
#' @param diag Logical. If true, returns measures of skew, and kurtosis
#' for each null distribution. Used to assess normality of null distribution.
#' @param mc.cores numerical. Number of cores for parallel processing
#'
#' @details
#' his function is flexible and takes a range of input formats such as
#' data.frames, vectors, matrices, and distance objects, and maintains this
#' format in the exported values. The function estimates standardize effect
#' sizes in the traditional z score method (SES). Additionally empirical
#' p-values based effect sizes (ES) are calculated. Empirical p-values describe
#' the proportion of the null model that is more extreme than the observed
#' values and p-value based effect sizes (ES) are calculated. Finally the
#' function reports optional diagnostic metric to assess if null distributions
#' are symmetrical and normal. asymmetrical null distributions should be
#' assessed using empirical p-value based effect sizes rather than z-score
#' based SES. See Botta-Dukát (2018) for more information on selecting SES or
#'  p-value based ES.
#'
#' @returns List containing:means of null distributions (null_mean), standard
#' deviation of null distribution (null_sd), z-score based standardized effect
#' sizes (ses), empirical_pvalues, empirical_pvalue-based effect size (ES),
#' skewness of null distribution (skew), and kurtosis of null
#' distribution (kurt)


ses_fun <- function(obs,null_list,emp_padding=0,diag = T, mc.cores = 1) {

  # Format all input data into vectors (obs) and matices (null) for
  # streamlined calculations
  if(is(obs, "dist")) {
    stopifnot("please supply null and obs values in same format" =
                all(dim(obs) == dim(null_list[[1]])))
    obs_format <- as.vector(obs)
    null <- do.call(rbind,lapply(null_list,as.vector))
  }

  if (is.matrix(obs) | is.data.frame(obs)) {
    stopifnot("please supply null and obs values in same format" =
                all(dim(obs) == dim(null_list[[1]])))
    stopifnot("If matrix, input has to be less than 3 dimensions"=length(dim(obs)) < 3)
    obs <- as.matrix(obs)
    obs_format <- unlist(as.vector(obs))
    null <- do.call(rbind,lapply(null_list,function(x) unlist(as.vector(as.matrix(x)))))
  }

  if(is.null(dim(obs))) {
    stopifnot("please supply null and obs values in same format" =
                length(obs) == length(null_list[[1]]))
    obs_format <- obs
    null <- do.call(rbind,null_list)
  }

  ## Summarize mean and sd of null model  ##
  null_mean <- matrixStats::colMeans2(null)
  null_sd <- matrixStats::colSds(null)

  ## Calculate ses  ##
  ses <- (obs_format-null_mean)/null_sd

  ##  Calaculate empirical p values ##

  # Transpose for better memory utlilization
  null_t <- t(null)

  # p value function
  emp_p <- sapply(1:length(obs_format), function(i){
    o <- obs_format[i]
    n <- null_t[i,]
    p <- ((sum(o<n)+sum(o==n)/2)+emp_padding)/(length(n)+emp_padding)
    p[p==1] <- 1-1/length(n)  # if p value is 1, apply the padding to upperside to prevent -inf es
    p
  })

  ## Estimate effect size based on emperical p  ##
  es <-VGAM::probitlink(1-emp_p)

  ## Prepare output  ##

  # Compile all results
  out <-list(
    obs = obs_format,
    null_mean = null_mean,
    null_sd = null_sd,
    ses = ses,
    empirical_pvalue = emp_p,
    empirical_es = es
  )

  ## normality diagnostics  ##
  if(diag){
    out$skew <-apply(null,2,DescTools::Skew)
    out$kurt <-apply(null,2,DescTools::Kurt)
  }

  # Revert output back to input structure
  lapply(
    out,
    structure,
    class = attr(obs,"class"),
    Size = attr(obs, "Size"),
    Labels =attr(obs, "Labels"),
    dim = attr(obs, "dim"),
    dimnames =attr(obs, "dimnames"),
    names = attr(obs, "names")
  )
}
