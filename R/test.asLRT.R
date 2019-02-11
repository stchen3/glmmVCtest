#' Asymptotic LRT for variance components in GLMMs
#'
#' Conducts the asymptotic likelihood ratio test (asLRT) from Molenberghs and Verbeke
#' (2007) for a single random effect (variance component) in a GLMM.
#' This method uses Laplace approximation (via lme4::glmer) to directly calculate
#' the likelihood.
#'
#' Compared to the null model, the alternative model should contain one additional
#' random effect (variance component) to be tested and (optionally) additional
#' fixed effects to be tested. Note that simultaneous testing of fixed and random effects
#' has not be extensively studied. Use results with caution.
#'
#' This method is known to have conservative
#' type I error rates due to the asymptotic null distribution.
#' Models without grouping variables for random effects (such as basis approximations)
#' cannot be estimated with \code{lme4::glmer}.
#'
#' @param fit.alt alternative model estimated using \code{MASS::glmer} function.
#' @param fit.null null model estimated using \code{MASS::glmer} function.
#'
#' @return Returns a list of components:
#' \itemize{
#' \item \code{statistic}: estimated LRT statistic.
#' \item \code{p}: p-value.
#' }
#'
#' @examples
#' \dontrun{
#' library(glmm)
#' data(salamander)
#' nobs <- nrow(salamander) # number of obs
#' glmer.alt<-glmer(Mate~0+Cross+(1|Male),
#'                 data=salamander,family=binomial,weights=rep(1,nobs))
#' glm.null<-glm(Mate~0+Cross,
#'                  data=salamander,family=binomial,weights=rep(1,nobs))
#' # Test significance of random intercept for Male salamanders
#' test.asympLRT(glmer.alt,glm.null)
#' }
#'
#' @seealso
#' \code{\link{test.nRLRT}}
#'
#' Chen, S T., Xiao, L., and Staicu, A. M. (in prep).
#' Restricted likelihood ratio tests for variance components in generalized linear mixed models.
#'
#' Molenberghs, G. and Verbeke, G. (2007).
#' Likehood ratio, score, and wald tests in a constrained parameter space.
#' \emph{The American Statistician} \strong{61}, 22--27.

test.asLRT<-function(fit.alt,fit.null){
  n.fix.alt<- length(fixef(fit.alt))
  n.fix.null<- length(fixef(fit.null))
  n.fixed.diff <- n.fix.alt - n.fix.alt
  asympGLRT.stat<- max(0,as.numeric(-2*(logLik(fit.null)-logLik(fit.alt))))
  null.dist <- asymptotic.null(n.fixed.diff)

  asympGLRT.p <- mean(asympGLRT.stat<=null.dist)

  return(list(statistic=asympGLRT.stat,p=asympGLRT.p))
}
