% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aScore.R
\name{aScore}
\alias{aScore}
\title{Approximate Score test for variance components in GLMMs.}
\usage{
aScore(null.resids, X, Z.test, V = NULL, Vinv = NULL, dVdv)
}
\arguments{
\item{null.resids}{vector of fixed-effects residuals from null model}

\item{X}{design matrix for fixed effects}

\item{Z.test}{design matrix for the random effect to be tested, with identity correlation structure.}

\item{V}{conditional variance of Y (optional if Vinv is given)}

\item{Vinv}{inverse of conditional variance of Y (optional if V is given)}

\item{dVdv}{column-bind of derivative of V in terms of all variance components (see Details)}
}
\value{
Returns a list of test components:
\itemize{
\item \code{U}: bias-corrected score statistic
\item \code{p}: p-value
\item \code{kappa}: scaling factor for null distribution (weighted chi-square)
\item \code{v}: df for null distribution (weighted chi-square)
}
}
\description{
This function conducts the approximate (bias-corrected) Score (aScore) test from Zhang and Lin (2003)
for testing zero-value variance components in generalized linear mixed models (GLMMs).
}
\details{
This function assumes that all random effects have identity correlation structures.
For more general correlations, consider other methods or modifying this function.
Model estimation are under the NULL hypothesis at convergence, and should use functions
glm (stats) or glmmPQL.mod (glmmVCtest) for estimation.
This method requires the design matrices for the fixed effects,
random effect being tested, conditional variance structure, and first derivative of the
conditional variance.
Thus, it is trickier to apply than the other methods in this package.

V is the variance of Y, conditional on all random effects under the null hypothesis.
The function will invert V using solve(V), and Vinv can be given instead if preferred.
dVdv is a column-bind of the derivative of V in terms of all variance components.
Order does matter, and the derivative of the residual variance should be last (right-most).

Example 1: If there are no random effects under the null hypothesis
(no nuisance effects), then the null model is a generalized linear model to be
estimated with the glm function in (stats). V and dVdv are identical,
and are simply the inverse of the estimated weights from glm.
X is the fixed effects matrix and Z.test is the design matrix of the random effect being
tested.

Example 2: If nuisance random effects are present, the null model should be
estimated using the glmmPQL.mod function. Normalized responses do not need to be
standardized. Then V = sigsq*I + sigsq_1*t(Z1)Z1 + sigsq_2*t(Z2)Z2 + ...,
and dVdv = cbind(t(Z1)Z1, ... , I).
Z.test is the design matrix of the random effect being tested.
}
\examples{
\dontrun{
library(glmm)
data(salamander) # load data
# Fit null hypothesis
null.fit <- glm(Mate~0+Cross, data=salamander, family=binomial)

# Test significance of random intercept for Male salamanders
w<- null.fit$weights
W <- diag(w)
Winv <- diag(1/w)
V <- Winv
dVdv <- Winv
null.resids <- null.fit$residuals

# Construct design matrices
Xmat <- as.matrix(model.matrix(null.fit))
Zmat.male<-matrix(0,nrow=n,ncol=length(unique(salamander$male)))
males<-unique(salamander$Male)
for(i in c(1:n)){
 col.id<-match(salamander$Male[i],males)
 Zmat.male[i,col.id]<-1
}

aScore.result<-aScore(null.resids,X=Xmat,V=V,dVdv=dVdv,Z.test=Zmat.male)
}
}
\seealso{
Lin, X. (1997).
Variance component testing in generalised linear models with random effects.
\emph{Biometrika} \strong{84}, 309--326.

Zhang, D. and Lin, X. (2003).
Hypothesis testing in semiparametric additive mixed models.
\emph{Biostatistics} \strong{4}, 57--74.
}
