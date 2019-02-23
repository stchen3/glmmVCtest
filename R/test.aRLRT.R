#' Approximate RLRT for variance components in GLMMs
#'
#' Conducts the approximate restricted likelihood ratio test (aRLRT) for a single
#' random effect (variance component) in a generalized linear mixed model (GLMM).
#'
#' This is our recommended method for testing variance components in GLMMs.
#' It assumes that the first specified random effect is the factor to be tested.
#' This random effect must be independent of all other random effects, but can have any
#' covariance structure.
#' For generalized responses, please use the \code{glmmPQL.mod} function for estimation,
#' it has minor modifications from \code{MASS::glmmPQL} to allow for testing. Please use
#' the \code{nlme::lme} function for estimating models with normal responses
#' (linear mixed models).
#'
#' NOTE: The model cannot use the ~(1|group) notation
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
#' @param fit Alternative model estimated using the \code{glmmPQL.mod} function (generalized responses)
#' or \code{\link[nlme:lme]{lme}} function (normal responses).
#'
#' @return Returns a list of testing results:
#' \itemize{
#' \item \code{aRLRT}: outcome of the (a)RLRT test (statistic and p-value)
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
test.aRLRT<-function(fit){
  if(class(fit)=='glmmPQL'){
    stop('Please use glmmPQL.mod to estimate the alternative hypothesis')
  }
  if(class(fit)=='lme'){ # normal responses
    if(length(fit$call$random)==2){ # only 1 random effect
      RLRT <- exactRLRT(fit)
      return(list(aRLRT=RLRT,std.data=fit$data,fit.alt=fit,fit.null=NULL,fit.test=NULL))
    } else { # multiple randome effects
      null.call <- test.call <- fit$call
      null.call$random<- fit$call$random[-2] # remove effect being tested
      fit.null<-eval(null.call)
      test.call$random <- fit$call$random[1:2] # effect being tested only
      fit.test <- eval(test.call)
      RLRT <- exactRLRT(m=fit.test,mA=fit,m0=fit.null)
      return(list(aRLRT=RLRT,std.data=fit$data,fit.alt=fit,fit.null=fit.null,fit.test=fit.test))
    }
  } else if(class(fit)=='list'){ # glmmPQL.mod output
    if(class(fit$fit)[1]=='glmmPQL'){ # generalized responses
      # Extract and standardize to Ytilde
      mcall.orig<-fit$mcall # lme call (X and Z haven't been adjusted to iid)
      mcall.std<-std.glmmPQL(fit) # updated mcall with standardized Ytilde, Xtilde, Ztilde
      std.data<-mcall.std$data

      # refit under null and alt hypothesis
      fit.alt.std<-try(eval(mcall.std),silent=T)
      if('try-error' %in% class(fit.alt.std)){
        stop('Error in lme model estimation under alternative. Consider simplifying
             or rescaling variables.')
      }

      if(length(mcall.std$random)==1){ # no nuisance random effects
        fit.null.std<-try(lm(mcall.std$fixed,data=std.data)) # fixed effects only
        if('try-error' %in% class(fit.null.std)){
          stop('Error in lm model estimation under null. Consider simplifying
               or rescaling variables.')
        }
        fit.test.std<-fit.alt.std # same as alternative model
        } else { # nuisance r.effect
          mcall.null<-mcall.test<-mcall.std # update from alt model fit
          mcall.null$random<-mcall.std$random[2:length(mcall.std$random)] #only nuis
          mcall.test$random<-mcall.std$random[1] # only test
          fit.null.std<-try(eval(mcall.null),silent=T)
          if('try-error' %in% class(fit.null.std)){
            stop('Error in lme model estimation under null. Consider simplifying
                 or rescaling variables.')
          }
          fit.test.std<-try(eval(mcall.test),silent=T)
          if('try-error' %in% class(fit.test.std)){
            stop('Error in lme model estimation for testing variable. Consider rescaling.')
          }
        }
      # testing
      aRLRT<-exactRLRT(m=fit.test.std,mA=fit.alt.std,m0=fit.null.std)
      return(list(aRLRT=aRLRT,std.data=std.data,fit.alt=fit.alt.std,fit.null=fit.null.std,fit.test=fit.test.std))
      } else {
        stop('Only lme and glmmPQL.mod class models are supported!')
      }
  } else {
    stop('Only lme and glmmPQL.mod class models are supported!')
  }
}

