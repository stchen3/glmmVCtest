#' Update function call with standardized data
#'
#' Updates the function call from glmmPQL.mod so the working LMM has
#' standardized residuals, for testing. The normalized response and design matrices
#' are all standardized accordingly.
#'
#' @param fit.glmmPQL Model estimate from glmmPQL.mod.
#'
#' @return Returns an updated function call for PQL estimation.
#'
std.glmmPQL<-function(fit.glmmPQL){
  mcall.unstd<-fit.glmmPQL$mcall # not standardized
  fit.unstd<-fit.glmmPQL$fit
  unstd.data<-mcall.unstd$data

  W.half<-diag(sqrt(1/unstd.data$invwt)) # sqaure-root weight matrix
  Ytilde<- W.half %*% unstd.data$zz # std responses

  # fixed effects
  X.design<-data.frame(model.matrix(fit.unstd,data=unstd.data)) # fixed effects
  Xtilde<- data.frame(W.half %*% as.matrix(X.design)) # std fixed effects
  Xtilde.names<- c('0',names(Xtilde)) # what about auto intercept?
  Xtilde.formula<-as.formula(paste("Ytilde~",paste(Xtilde.names,collapse="+"))) # new fixed formula
  std.fixed<-data.frame(Ytilde=Ytilde,Xtilde)

  # random effects
  group.names<-names(fit.unstd$groups) # grouping variables
  extract.groups <- unstd.data[,(names(unstd.data) %in% c(group.names,'invwt'))] #doesn't keep names
  n.grouping<-length(mcall.unstd$random)
  std.rand.all<-list()
  for(i in 1:n.grouping){ # std all r.effects
    r.effects<-data.frame(model.matrix(formula(fit.unstd$modelStruct$reStr)[[i]],data=unstd.data))
    std.rand<-data.frame(as.matrix(r.effects)*diag(W.half)) #
    std.rand.all[[i]]<-std.rand
  }
  std.rand<-do.call(cbind,std.rand.all)

  std.data.full<-cbind(std.fixed,std.rand,extract.groups) # standardized data

  mcall.std<-mcall.unstd # update mcall
  mcall.std$fixed<-Xtilde.formula # new fixed effects formula
  mcall.std$data<-std.data.full # std data
  mcall.std$weights<-NULL # iid errors

  return(mcall.std)
}
