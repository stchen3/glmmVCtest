#' Approximate RLRT for variance components in GLMMs
#'
#' Conducts the approximate restricted likelihood ratio test (aRLRT) for a single
#' random effect (variance component) in a generalized linear mixed model (GLMM).
#'
#' This is our recommended method for testing variance components in GLMMs.
#' It assumes that the first specified random effect is the factor to be tested.
#' This random effect must be independent of all other random effects, but can have any
#' covariance structure.
#' Please use the \code{glmmPQL.mod} function for model estimation,
#' it has minor modifications from \code{MASS::glmmPQL} to allow for testing.
#'
#' NOTE: The glmmPQL model cannot use the ~(1|group) notation
#' for any random effects (fixed are ok). Please specify a random interept directly
#' by adding an intercept column (column of 1s) and using ~(0+'name'|group). Variable name
#' does not matter. See example below.
#'
#' The function modifies and re-estimates the wLMM from PQL estimation to have iid
#' residuals, then re-estimates and formats model estimates for \code{RLRsim::exactRLRT}
#' to conduct testing. That function uses a finite-sample null distribution
#' that improves testing performance compared to the typical asymptotic result
#' from Self and Liang (1987).
#'
#' The function does not currently support simultaneous testing of fixed and random effects.
#'
#' @param fit.glmmPQL Alternative model estimated using the \code{glmmPQL.mod} function.
#'
#' @return Returns a list of testing results:
#' \itemize{
#' \item \code{nLRT}: outcome of the aRLRT test (statistic and p-value)
#' \item \code{std.data}: data.frame of standardized normalized responses (Ytilde),
#' fixed effects design, and random effects design
#' \item \code{fit.alt}: lme under the alternative hypothesis
#' \item \code{fit.null}: lme under the null hypothesis
#' \item \code{fit.test}: lme for just the random effect to be tested
#' }
#'
#' @examples
#' \dontrun{
#' library(glmm)
#' data(salamander) # load data
#' salamander$int <- 1 # add column of 1s for random intercept
#' # Fit alternative hypothesis
#' glmm.fit <- glmmPQL.mod(Mate~0+Cross, random=list(~0+int|Male),
#' data=salamander, family=binomial, weights=rep(1,nrow(salamander)))
#' # Test significance of random intercept for Male salamanders
#' aRLRT <- test.aRLRT(glmm.fit)
#' }
#'
#' @seealso
#' Chen, S. T., Xiao, L., Staicu, A. M. (in prep).
#' Restricted Likelihood Ratio Tests for Variance Components in Generalized Linear Models.
#'
test.aRLRT<-function(fit.glmmPQL){
  # Extract and standardize to Ytilde
  mcall.orig<-fit.glmmPQL$mcall # lme call (X and Z haven't been adjusted to iid)
  mcall.std<-std.glmmPQL(fit.glmmPQL) # updated mcall with standardized Ytilde, Xtilde, Ztilde
  std.data<-mcall.std$data

  # refit under null and alt hypothesis
  fit.alt<-try(eval(mcall.std),silent=T)
  if('try-error' %in% class(fit.alt)){
    stop('Error in lme model estimation under alternative. Consider simplifying
         or rescaling variables.')
  }

  if(length(mcall.std$random)==1){ # no nuisance random effects
    fit.null<-try(lm(mcall.std$fixed,data=std.data)) # fixed effects only
    if('try-error' %in% class(fit.null)){
      stop('Error in lm model estimation under null. Consider simplifying
           or rescaling variables.')
    }
    fit.test<-fit.alt # same as alternative model
  } else { # nuisance r.effect
    mcall.null<-mcall.test<-mcall.std # update from alt model fit
    mcall.null$random<-mcall.std$random[2:length(mcall.std$random)] #only nuis
    mcall.test$random<-mcall.std$random[1] # only test
    fit.null<-try(eval(mcall.null),silent=T)
    if('try-error' %in% class(fit.null)){
      stop('Error in lme model estimation under null. Consider simplifying
           or rescaling variables.')
    }
    fit.test<-try(eval(mcall.test),silent=T)
    if('try-error' %in% class(fit.test)){
      stop('Error in lme model estimation for testing variable. Consider rescaling.')
    }
  }

  # testing
  aRLRT<-exactRLRT(m=fit.test,mA=fit.alt,m0=fit.null)

  return(list(aRLRT=aRLRT,std.data=std.data,fit.alt=fit.alt,fit.null=fit.null,fit.test=fit.test))
}
