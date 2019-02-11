#' Sample from the asymptotic null distribution.
#'
#' Sample from the asymptotic null distribution for zero-value
#' variance component from Self and Liang (1987).
#'
#' @param n.fix.test Number of fixed effects being tested.
#' @param nsim Number of samples to draw.
#' @examples
#' # To test 1 random effect and 0 fixed effects
#' asymptotic.null(0)
#' @seealso Self, G. S. and Liang K. Y. (1987).
#' Asymptotic properties of maximum likelihood estimators and likelihood ratio tests under nonstandard conditions.
#' \emph{Journal of the American Statistical Association} \strong{82}, 605--610.
#'
#' Stram, D. O. and Lee, J. W. (1994).
#' Variance component testing in the longitudinal mixed effects model.
#' \emph{Biometrics} \strong{50}, 1171--1177.
asymptotic.null<-function(n.fix.test,nsim=10000){
  c(rchisq(nsim/2,df=n.fix.test),rchisq(nsim/2,df=(n.fix.test+1)))
}
