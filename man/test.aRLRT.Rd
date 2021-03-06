% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test.aRLRT.R
\name{test.aRLRT}
\alias{test.aRLRT}
\title{Approximate RLRT for variance components in GLMMs}
\usage{
test.aRLRT(fit)
}
\arguments{
\item{fit}{Alternative model estimated using the \code{glmmPQL.mod} function (generalized responses)
or \code{\link[nlme:lme]{lme}} function (normal responses).}
}
\value{
Returns a list of testing results:
\itemize{
\item \code{aRLRT}: outcome of the (a)RLRT test (statistic and p-value)
\item \code{std.data}: data.frame of standardized normalized responses (Ytilde),
fixed effects design, and random effects design
\item \code{fit.alt}: lme under the alternative hypothesis
\item \code{fit.null}: lme under the null hypothesis
\item \code{fit.test}: lme for just the random effect to be tested
}
}
\description{
Conducts the approximate restricted likelihood ratio test (aRLRT) for a single
random effect (variance component) in a generalized linear mixed model (GLMM).
}
\details{
This is our recommended method for testing variance components in GLMMs.
It assumes that the first specified random effect is the factor to be tested.
This random effect must be independent of all other random effects, but can have any
covariance structure.
Please use the \code{glmmPQL.mod} function for model estimation,
it has minor modifications from \code{MASS::glmmPQL} to allow for testing.

NOTE: The glmmPQL model cannot use the ~(1|group) notation
for any random effects (fixed are ok). Please specify a random interept directly
by adding an intercept column (column of 1s) and using ~(0+'name'|group). Variable name
does not matter. See example below.

The function modifies and re-estimates the wLMM from PQL estimation to have iid
residuals, then re-estimates and formats model estimates for \code{RLRsim::exactRLRT}
to conduct testing. That function uses a finite-sample null distribution
that improves testing performance compared to the typical asymptotic result
from Self and Liang (1987).

The function does not currently support simultaneous testing of fixed and random effects.
}
\examples{
\dontrun{
library(glmm)
data(salamander) # load data
salamander$int <- 1 # add column of 1s for random intercept
# Fit alternative hypothesis
glmm.fit <- glmmPQL.mod(Mate~0+Cross, random=list(~0+int|Male),
data=salamander, family=binomial, weights=rep(1,nrow(salamander)))
# Test significance of random intercept for Male salamanders
aRLRT <- test.aRLRT(glmm.fit)
}

}
\seealso{
Chen, S. T., Xiao, L., Staicu, A. M. (in prep).
Restricted Likelihood Ratio Tests for Variance Components in Generalized Linear Models.
}
