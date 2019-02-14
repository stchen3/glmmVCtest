#' Asymptotic-approximate RLRT for variance components in GLMMs
#'
#' Conducts an asymptotic-approximate restricted-likelihood ratio test (as-aRLRT)
#' for random effects (variance components) in generalized linear mixed models (GLMMs).
#' The method applies RLRT to the working linear mixed model (LMM) approximation
#' for the GLMM used in PQL estimation. The test statistic is compared to the
#' asymptotic null distribution from Self and Liang (1987).
#'
#' Compared to the recommended aRLRT, this method instead uses an asymptotic null
#' distribution (Self and Liang 1987), which typically has conservative type I error
#' rates and lower power. The function assumes that the first specified random effect is
#' the factor of interest.
#' This random effect must be independent of all other random effects, but can have any
#' covariance structure.
#' Please use the \code{glmmPQL.mod} function in this package; it has minor modifications
#' from \code{MASS::glmmPQL}.
#' NOTE: The glmmPQL model cannot use the ~(1|group) notation
#' for any random effects (fixed are ok). Please specify a random interept directly
#' by adding an intercept column (column of 1s) and using ~(0+'name'|group). Variable name
#' does not matter.
#'
#' The function does not currently support simultaneous testing of fixed and random effects.
#' That may be added in the future.
#'
#' @param fit.glmmPQL Alternative model estimated using the \code{glmmPQL.mod} function.
#'
#' @return Returns a list of testing results:
#' \itemize{
#' \item \code{asaRLRT}: outcome of the asaRLRT test (statistic and p-value)
#' \item \code{std.data}: data.frame of standardized normalized responses (Ytilde),
#' fixed effects design, and random effects design
#' \item \code{fit.alt}: lme under the alternative hypothesis
#' \item \code{fit.null}: lme under the null hypothesis
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
#' # Test significance of random subject-specific intercept for Male salamanders
#' asaRLRT <- test.asaRLRT(glmm.fit)
#' }
#'
#' @seealso
#' \code{\link{test.aRLRT}}

test.asaRLRT<-function(fit.glmmPQL,nsim=10000){
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
  n.fix.alt<-length(fit.alt$coefficients$fixed)

  if(length(mcall.std$random)==1){ # no nuisance random effects
    fit.null<-try(lm(mcall.std$fixed,data=std.data)) # fixed effects only
    if('try-error' %in% class(fit.null)){
      stop('Error in lm model estimation under null. Consider simplifying
           or rescaling variables.')
    }
    n.fix.null<-length(fit.null$coefficients)
  } else { # nuisance r.effect
    mcall.null<-mcall.test<-mcall.std # update from alt model fit
    mcall.null$random<-mcall.std$random[2:length(mcall.std$random)] #only nuis
    mcall.test$random<-mcall.std$random[1] # only test
    fit.null<-try(eval(mcall.null),silent=T)
    if('try-error' %in% class(fit.null)){
      stop('Error in lme model estimation under null. Consider simplifying
           or rescaling variables.')
    }
    n.fix.null<-length(fit.null$coefficients$fixed)
  }

  # testing (only works for r.effects for now)
  asaRLRT.stat <- max(0,as.numeric(-2*(logLik(fit.null,REML=T)-logLik(fit.alt,REML=T))))
  n.fixed.diff <- n.fix.alt - n.fix.alt
  null.dist <- asymptotic.null(n.fixed.diff)
  asaRLRT.p <- mean(asaRLRT.stat<=null.dist)
  asaRLRT <- list(statistic=asaRLRT.stat,p=asaRLRT.p)

  return(list(asaRLRT=asaRLRT,std.data=std.data,fit.alt=fit.alt,fit.null=fit.null))
}
