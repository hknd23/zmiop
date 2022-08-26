##################################################################
################Write a robust IOP model Program##################
##################################################################


#' @useDynLib IOP
#' @importFrom stats dgamma runif as.formula model.frame model.matrix model.response na.omit na.pass
#' @import alr3
#' @import RcppArmadillo
#' @importFrom Rcpp sourceCpp
#' @importFrom mvtnorm rmvnorm dmvnorm
#' @export

#' @title ziop
#' @description Likelihood function for ZIOP model.
#'
#' @param est starting values for the estimation. Vector of length of the number of parameters.
#' @param Y the ordinal dependent variable.
#' @param X covariates for the ordered probit stage.
#' @param Z covariates for the inflation (split population) stage.
#' @param data dataset that contains the dependent and independent variables.
#' @param weights an optional vector of weights to be used in the fitting process.
#' @param offsetx offest value for the ordered probit stage covariates (X). This can be used to specify an a priori known component to be included in the linear predictor during fitting. For more information, see \code{\link{offset}}.
#' @param offsetz offest value for the inflation (split population) stage covariates (Z). This can be used to specify an a priori known component to be included in the linear predictor during fitting. For more information, see \code{\link{offset}}.
#'
#' @return Likelihood of the ZIOP model specification.
#'
#' @export
#'
##### ZIOP MODEL #####
ziop<- function(est,Y,X,Z,data, weights, offsetx, offsetz) {
  n=nrow(data)
  llik <- matrix(0, nrow=n, ncol = 1)
  y.cat<-nlevels(as.factor(Y))
  y0<-sort(unique(Y))
  V<-matrix(NA,nrow(data),y.cat)
  for(k in 1:y.cat){
    V[,k]<-Y==y0[k]
  }
  tau<-rep(1,y.cat)
  for (i in 1:(y.cat-1)){
    tau[i]<-ifelse(i==1,est[i],tau[i-1]+exp(est[i]))
  }
  gamma<-est[(y.cat):(y.cat+ncol(Z)-1)]
  beta<-est[(y.cat+ncol(Z)):length(est)]
  ZG<-Z%*%gamma + offsetz
  XB<-X%*%beta + offsetx


  cprobs<-matrix(nrow=length(XB),ncol=y.cat)
  probs<-matrix(nrow=n,ncol=y.cat)

  for(j in 1:(y.cat-1)){

    cprobs[,j]<-pnorm(tau[j]-XB)

  }


  probs[,1]<-(cprobs[,1])*pnorm(ZG)+(1-pnorm(ZG))	#first category

  for(j in 2:(y.cat-1)){

    probs[,j]<-(cprobs[,j]-cprobs[,(j-1)])*(pnorm(ZG)) #middle categories

  }

  probs[,y.cat]<-(1-cprobs[,(y.cat-1)])*pnorm(ZG) #last category


  llik<--1*sum( log(probs[V])*weights)

  return(llik)

}


#' @title miop
#' @description Likelihood function for MIOP model.
#'
#' @param est starting values for the estimation. Vector of length of the number of parameters.
#' @param Y the ordinal dependent variable.
#' @param X covariates for the ordered probit stage.
#' @param Z covariates for the inflation (split population) stage.
#' @param data dataset that contains the dependent and independent variables.
#' @param weights an optional vector of weights to be used in the fitting process.
#' @param offsetx offest value for the ordered probit stage covariates (X). This can be used to specify an a priori known component to be included in the linear predictor during fitting. For more information, see \code{\link{offset}}.
#' @param offsetz offest value for the inflation (split population) stage covariates (Z). This can be used to specify an a priori known component to be included in the linear predictor during fitting. For more information, see \code{\link{offset}}.
#'
#' @return Likelihood of the MIOP model specification.
#'
#' @export
#'
#'
##### MIOP MODEL #####
miop <- function(est,Y,X,Z,data, weights, offsetx, offsetz) {
  n=nrow(data)
  llik <- matrix(0, nrow=n, ncol = 1)
  y.cat<-nlevels(as.factor(Y))
  y0<-sort(unique(Y))
  V<-matrix(NA,nrow(data),y.cat)
  for(k in 1:y.cat){
    V[,k]<-Y==y0[k]
  }
  tau<-rep(1,y.cat)
  for (i in 1:(y.cat-1)){
    tau[i]<-ifelse(i==1,est[i],tau[i-1]+exp(est[i]))
  }
  gamma<-est[(y.cat):(y.cat+ncol(Z)-1)]
  beta<-est[(y.cat+ncol(Z)):length(est)]

  ZG<-Z%*%gamma + offsetz
  XB<-X%*%beta + offsetx

  cprobs<-matrix(nrow=length(XB),ncol=y.cat)
  probs<-matrix(nrow=n,ncol=y.cat)

  for(j in 1:(y.cat-1)){

    cprobs[,j]<-pnorm(tau[j]-XB)

  }

  probs[,1]<- pnorm(ZG)*(cprobs[,1]) #first category

  for(j in 2:(y.cat-1)){

    probs[,j]<- (1-pnorm(ZG)) + pnorm(ZG)*(cprobs[,j]-cprobs[,(j-1)]) #middle category

  }

  probs[,y.cat]<-(1-cprobs[,(y.cat-1)])*pnorm(ZG) #last category

  llik<--1*sum( log(probs[V])*weights)

  return(llik)

}


#' @title iop.mod.Est
#' @description Raw form of the \code{ \link{iop}} function. For user-friendly formula-oriented command, use \code{ \link{iop}}.
#'
#' @param x covariates for the ordered probit stage.
#' @param z covariates for the inflation (split population) stage.
#' @param y the ordinal dependent variable.
#' @param weights an optional vector of weights to be used in the fitting process.
#' @param offsetx offest value for the ordered probit stage covariates (X). This can be used to specify an a priori known component to be included in the linear predictor during fitting. For more information, see \code{\link{offset}}.
#' @param offsetz offest value for the inflation (split population) stage covariates (Z). This can be used to specify an a priori known component to be included in the linear predictor during fitting. For more information, see \code{\link{offset}}.
#' @param na.action a function indicating what should happen when NAs are included in the data. Options are "na.omit" or "na.fail". The default is "na.omit".
#' @param type type of inflation ordered probit model to be used. Options are "ziop" or "miop". The type of the inflation model must be specified.
#'
#'
#' @export
#'
#A function for the iop.mod model
iop.mod.Est<- function(x, z, y, weights, offsetx, offsetz, na.action, type) {

	## Define data set
	dataset<-cbind(y,z[,1],scale(z[,2:ncol(z)],center=F),scale(x,center=F))
	xnew<-scale(x,center=F)
	znew<-cbind(z[,1],scale(z[,2:ncol(z)],center=F))

	y.category<-nlevels(as.factor(y))



	## weights
	weights <- as.vector(weights)


	##offsetx
	offsetx <- as.matrix(offsetx)
	offsetz <- as.matrix(offsetz)

	##na.action
	na.action <- na.action

	## Optimize by model

	type <- type

	## ZIOP
	if(type=="ziop"){

	## Define starting parameters
	est<-rep(.01, (ncol(x)+ncol(z)+(y.category-1)))

	output.iop.mod<-optim(fn=ziop, par=est, Y=y, X=xnew, Z=znew, weights=weights, offset=offsetx,
	             offsetz= offsetz, method="BFGS", control=list(maxit=500), data=dataset,hessian=TRUE)

	#df
	n<-nrow(dataset)
	df.null<- nrow(dataset)-((y.category-1)+1)
	df.fitted <- nrow(dataset)-(ncol(x)+ncol(z)+(y.category-1))

	## coefficients
	coef<-NULL
	cutnames<-NULL
	covariates<-cbind(z,x)
	covariates.only<-covariates
	for(i in 1:(y.category-1)){
	  descaledest<-output.iop.mod$par[i]
	  coef<-c(coef,descaledest)
	  name<-paste("cut",i,sep=" ")
	  cutnames<-c(cutnames,name)
	  covariates<-cbind(1,covariates)
	}
	for(i in 1:(ncol(covariates.only))){
	  descaledest<-output.iop.mod$par[i+(y.category-1)]/attributes(scale(covariates.only[,i],center=F))$"scaled:scale"
	  coef<-c(coef,descaledest)
	}
	vcv<-solve(output.iop.mod$hessian)
	vcov<-matrix(NA,nrow(vcv),ncol(vcv))
	for(i in 1:nrow(vcov)){
	  for(j in 1:ncol(vcov)){
	    vcov[i,j]<-vcv[i,j]/((attributes(scale(covariates[,i],center=F))$"scaled:scale")*(attributes(scale(covariates[,j],center=F))$"scaled:scale"))
	  }
	}

	names(coef)<-c(cutnames,paste("inflation",ifelse(colnames(z)=="(Intercept)","Intercept",colnames(z)),sep="."),colnames(x))
	coef<-list(cutpoints=coef[1:(y.category-1)],inflation=coef[y.category:(y.category+ncol(z)-1)],ordered=coef[(y.category+ncol(z)):length(est)])
	colnames(znew)<-paste("inflation",ifelse(colnames(z)=="(Intercept)","Intercept",colnames(z)),sep=".")

	list(coefficients =  coef, vcov = vcov,n=n,df.null=df.null,df.residual = df.fitted,
	     ll=output.iop.mod$value*-1,iterations=output.iop.mod$counts[[1]], convergence=output.iop.mod$convergence==0,optim=output.iop.mod,
	     weights= if(all(weights==1)) NULL else weights, na.action=na.action, offset = list(count = if(all(offsetx==0)) NULL else offsetx,
	                                                                                        inflation = if(all(offsetz==0)) NULL else offsetz), y=y,x=as.data.frame(xnew),z=as.data.frame(znew),y.categories=y.category)
	} else{


	## MIOP
	if(type=="miop"){

	  ## Define starting parameters
	  est<-rep(.01, (ncol(x)+ncol(z)+(y.category-1)))

	  output.iop.mod<-optim(fn=miop, par=est, Y=y, X=xnew, Z=znew, weights=weights, offset=offsetx,
              offsetz= offsetz, method="BFGS", control=list(maxit=500), data=dataset,hessian=TRUE)

	  #df
	  n<-nrow(dataset)
	  df.null<- nrow(dataset)-((y.category-1)+1)
	  df.fitted <- nrow(dataset)-(ncol(x)+ncol(z)+(y.category-1))

	  ## coefficients
	  coef<-NULL
	  cutnames<-NULL
	  covariates<-cbind(z,x)
	  covariates.only<-covariates
	  for(i in 1:(y.category-1)){
	    descaledest<-output.iop.mod$par[i]
	    coef<-c(coef,descaledest)
	    name<-paste("cut",i,sep=" ")
	    cutnames<-c(cutnames,name)
	    covariates<-cbind(1,covariates)
	  }
	  for(i in 1:(ncol(covariates.only))){
	    descaledest<-output.iop.mod$par[i+(y.category-1)]/attributes(scale(covariates.only[,i],center=F))$"scaled:scale"
	    coef<-c(coef,descaledest)
	  }
	  vcv<-solve(output.iop.mod$hessian)
	  vcov<-matrix(NA,nrow(vcv),ncol(vcv))
	  for(i in 1:nrow(vcov)){
	    for(j in 1:ncol(vcov)){
	      vcov[i,j]<-vcv[i,j]/((attributes(scale(covariates[,i],center=F))$"scaled:scale")*(attributes(scale(covariates[,j],center=F))$"scaled:scale"))
	    }
	  }

	  names(coef)<-c(cutnames,paste("inflation",ifelse(colnames(z)=="(Intercept)","Intercept",colnames(z)),sep="."),colnames(x))
	  coef<-list(cutpoints=coef[1:(y.category-1)],inflation=coef[y.category:(y.category+ncol(z)-1)],ordered=coef[(y.category+ncol(z)):length(est)])
	  colnames(znew)<-paste("inflation",ifelse(colnames(z)=="(Intercept)","Intercept",colnames(z)),sep=".")

	  list(coefficients =  coef, vcov = vcov,n=n,df.null=df.null,df.residual = df.fitted,
	       ll=output.iop.mod$value*-1,iterations=output.iop.mod$counts[[1]], convergence=output.iop.mod$convergence==0,optim=output.iop.mod,
	       weights= if(all(weights==1)) NULL else weights, na.action=na.action, offset = list(count = if(all(offsetx==0)) NULL else offsetx,
	                                                                                          inflation = if(all(offsetz==0)) NULL else offsetz), y=y,x=as.data.frame(xnew),z=as.data.frame(znew),y.categories=y.category)
	}
}
}


#' @title iop.mod.default
#' @description Default method for a \code{ \link{iop}}.
#'
#' @param object an object of class \code{iop.mod} (output of \code{\link{iop}}).
#'
#' @export
#'
#To add a default method, we define a function called iop.mod.default():
#This function tries to convert its first argument x into a matrix (and will throw an error if not successful)
#When then call our function for parameter estimation, and add fitted values, residuals and the function call to the results.
#Finally we set the class of the return object to "iop".
iop.mod.default<-function(x, z, y, weights, offsetx, offsetz, na.action, type,...)
{
	x <- as.matrix(x)
	z <- as.matrix(z)
	y <- as.numeric(y)
	weights <- as.vector(weights)
	offsetx <- as.vector(offsetx)
	offsetz <- as.vector(offsetz)
	na.action <- na.action
	type <- type
	est <-iop.mod.Est(x, z, y, weights, offsetx, offsetz, na.action, type)
	est$call <- match.call()
	class(est) <- "iop.mod"
	est
}

#' @title print.iop.mod
#' @description Print method for a \code{ \link{iop}} object.
#'
#' @param object an object of class \code{iop.mod} (output of \code{\link{iop}}).
#'
#' @export
#'
#Defining the print() method for our class as
print.iop.mod <- function(object, ...)
{
	cat("Call:\n")
	print(object$call)
	cat("\nCoefficients:\n")
	print(object$coefficients)
}


#' @title summary.iop.mod
#' @description Summary method for a \code{ \link{iop}} object.
#'
#' @param object an object of class \code{iop.mod} (output of \code{\link{iop}}.
#'
#' @export
#'
#summary method
summary.iop.mod <- function(object, ...)
{
	#build coefficient vector for cutpoints
	cutpoints<-object$coefficients$cutpoints
	#extract correct standard errors
	cut.se<-sqrt(object$vcov[1,1])
	cut.coeff<-cutpoints[1]
	for(k in 2:length(cutpoints)){
		coeff<-cutpoints[1:k]
		varcov<-matrix(NA,k,k)
		for(i in 1:k){
			for(j in 1:k){
			varcov[i,j]<-object$vcov[i,j]
			}
		}

		cut.coeff<-c(cut.coeff,cut.coeff[k-1]+exp(cutpoints[k]))
	}
	newcoeff<-c(object$coefficients$inflation,object$coefficients$ordered,cut.coeff)
	names(newcoeff)<-c(names(object$coefficients$inflation),names(object$coefficients$ordered),names(object$coefficients$cutpoints))
	se <- c(sqrt(diag(object$vcov[(length(cutpoints)+1):ncol(object$vcov),(length(cutpoints)+1):ncol(object$vcov)])),
	        sqrt(diag(object$vcov[(1:length(cutpoints)),1:(length(cutpoints))])))
	tval <- newcoeff/ se

	TAB <- cbind(Estimate = newcoeff,
		StdErr = se,
		t.value = tval,
		p.value = 2*pt(-abs(tval), df=object$df.residual))
  	res<- list(call=object$call,
		coefficients=TAB,
		parameters=length(coef(object)),
		iterations=object$iterations,
		ll=round(object$ll,2))

	class(res) <-"summary.iop.mod"
	res
}

#' @title print.iop.mod
#' @description Print method for a \code{ \link{iop}} object.
#'
#' @param object an object of class \code{iop.mod} (output of \code{\link{iop}}).
#'
#' @export
#'
#The utility function printCoefmat() can be used to print the matrix with appropriate rounding and some decoration
print.summary.iop.mod <- function(object, ...)
{
	cat("Call:\n")
	print(object$call)
	cat("\n")

	printCoefmat(object$coefficients, P.value=TRUE, has.Pvalue=TRUE)

	cat("\n")
	cat(noquote(paste("Number of iterations in BFGS optimization: ", object$iterations,sep="")))
	cat("\n")
	cat(noquote(paste("Log-likelihood: ",object$ll, " on ", object$parameters, " Df",sep="")))
	cat("\n")
}



#' @title iop
#' @description \code{iop} fits an ordered probit model with inflations in either the "zero (bottom)" or "middle" categories.
#' @param formula a formula in the form Y ~ X1 + X2... | Z1 + Z2 ... where Y is the ordered probit dependent variable; Xi are the ordered probit stage covariates; and the Zi are the inflation (split population) stage covariates. See \code{link{formula}}.
#' @param data list object of data.
#' @param weights an optional vector of weights to be used in the fitting process. Default is NULL.
#' @param offset This can be used to specify an a priori known component to be included in the linear predictor during fitting. The same offset is applied to both stages. See \code{\link{offset}}.
#' @param na.action a function indicating what should happen when NAs are included in the data. Options are "na.omit" or "na.fail". The default is "na.omit".
#' @param type type of inflation ordered probit model to be used. Options are "ziop" or "miop". The type of the inflation model must be specified.
#'
#' @examples
#'  model1 <- iop(Y ~ X1 + X2|Z1 + Z2, data=data, type=c('ziop'))
#'
#' @export
#'
iop<- function(formula, data=list(), weights=NULL, offset=NULL, na.action=c("na.omit","na.fail"),
                            type=c("ziop", "miop"), ...)
{ if (missing(na.action)){na.action <- "na.omit"}

  n <- nrow(data)

  if(is.null(weights)) {weights <- rep.int(1, n)}
  if(length(weights) == 1) {weights <- rep.int(weights, n)}
  weights <- as.matrix(weights)

  if(is.null(offset)) {offset <- rep.int(0, n)}
  if(length(offset) == 1) {offset <- rep.int(offset, n)}
  offset <- as.matrix(offset)

	equations<-as.character(formula)
	formula1<-paste(equations[2],equations[1],strsplit(equations[3], "|", fixed = TRUE)[[1]][1],"-1",sep="")
	formula2<-paste(equations[2],equations[1],strsplit(equations[3], "|", fixed = TRUE)[[1]][2],sep="")
	mf1 <- model.frame(formula=as.formula(formula1), data=data, weights=weights, offset=offset, na.action = na.pass)
	mf2 <- model.frame(formula=as.formula(formula2), data=data, weights=weights, offset=offset, na.action = na.pass)
	X <- model.matrix(attr(mf1, "terms"), data=mf1)
	Z <- model.matrix(attr(mf2, "terms"), data=mf2)
	Y <- model.response(mf2)

	weights <- model.weights(mf2)

	offsetx <- model.offset(mf1)
	if (length(offsetx)==1){offsetx <- rep.int(offsetx, n)}

	offsetz <- model.offset(mf2)
	if (length(offsetz)==1){offsetz <- rep.int(offsetx, n)}

	offset <- offset

	na.action <- na.action

	dataset <- data.frame(cbind(Y, X, Z, weights, offsetx, offsetz))


	if(na.action=="na.omit"){
	  dataset <- na.omit(dataset)
	  col <- ncol(dataset)
	  y <- as.matrix(dataset[,1], ncol=1)
	  colnames(y) <- colnames(Y)
	  x <- data.frame(dataset[,2:(ncol(X)+1)])
	  names(x) <- colnames(X)
	  z <- dataset[,(ncol(X)+2):(ncol(X)+ncol(Z)+1)]
	  names(z) <- colnames(Z)
	  weights <- dataset[,(col-2)]
	  offsetx <- dataset[,(col-1)]
	  offsetz <- dataset[,col]

	est<-iop.mod.default(x, z, y, weights, offsetx, offsetz, na.action, type, ...)
	est$call <- match.call()
	est$formula <- formula
	est
	}
	else{
	  if(na.action=="na.fail"){
	    if(all(is.numeric(dataset))){
	      x <- X
	      z <- Z
	      y <- Y
	      est<-iop.mod.default(x, z, y, weights, offsetx, offsetz, na.action, type, ...)
	      est$call <- match.call()
	      est$formula <- formula
	      est
	    }
	    else{
	      x <- X
	      z <- Z
	      y <- Y
	      if(any(is.na(x))) warning("missing values in object x")
	      if(any(is.na(z))) warning("missing values in object z")
	      if(any(is.na(y))) warning("missing values in object y")
	      if(any(is.na(weights))) warning("missing values in object weights")
	      if(any(is.na(offsetx))) warning("missing values in object offset")
	      if(any(is.na(offsetz))) warning("missing values in object offset")
	      if(any(is.na(offset))) warning("missing values in object offset")
	      if(any(is.na(type))) warning("missing values in object type")

	      }
	  }
	  else{		dataset <- na.omit(dataset)
	  col <- ncol(dataset)
	  y <- as.matrix(dataset[,1], ncol=1)
	  colnames(y) <- colnames(Y)
	  x <- data.frame(dataset[,2:(ncol(X)+1)])
	  names(x) <- colnames(X)
	  z <- dataset[,(ncol(X)+2):(ncol(X)+ncol(Z)+1)]
	  names(z) <- colnames(Z)
	  weights <- dataset[,(col-2)]
	  offsetx <- dataset[,(col-1)]
	  offsetz <- dataset[,col]

	  est<-iop.mod.default(x, z, y, weights, offsetx, offsetz, na.action, type, ...)
	  est$call <- match.call()
	  est$formula <- formula
	  est

	  }
	}
}


#' @title fitted
#' @description A function that extracts fitted values from an object of class \code{iop.mod}.
#'
#' @param object an object of class \code{iop.mod} (output of \code{\link{iop}}).
#' @param newdata An optional data frame in which to look for variables to use when model fitting.
#' @param type the tye of equation to be fitted. Options include "response.full" (both inflation and ordered probit stages), "response.ordered" (ordered probit stage only) and "response.inflation" (inflation stage only.
#'
#'
#' @examples
#'  model1 <- iop(Y ~ X1 + X2|Z1 + Z2, data=data, type=c('ziop'))
#'  fitted(model1, type=c("response.full"))
#'
#' @export
#'
#'
#'
fitted <- function(object, newdata=NULL, type = c("response.full","response.ordered", "response.inflation", "linear"), ...)
{
	if(missing(type)){
		if(is.null(newdata)){
			ZG<-as.matrix(object$z) %*% object$coefficients$inflation
			XB<-as.matrix(object$x) %*% object$coefficients$ordered
			cprobs<-matrix(nrow=object$y.categories-1,ncol=1)
			n <- nrow(object$z)
			probs<-matrix(nrow=n,ncol=object$y.categories)
			cprobs[1,1]<-object$coefficients$cutpoints[1]
			for(j in 2:(object$y.categories-1)){

				cprobs[j,1]<-cprobs[j-1,1]+exp(object$coefficients$cutpoints[j])
			}

			probs[,object$y.categories]<-(pnorm(ZG))*((1-pnorm(cprobs[object$y.categories-1,1]-XB)))
			probs[,1]<-((1-pnorm(ZG))+(pnorm(ZG))*(pnorm(cprobs[1,1]-XB)))

			for(j in 2:(object$y.categories-1)){

			probs[,j]<-(pnorm(ZG))*((pnorm(cprobs[j,1]-XB)) - (pnorm(cprobs[j-1,1]-XB)))
			}

		}
 		else{
			equations<-as.character(object$formula)
			formula1<-paste(equations[2],equations[1],strsplit(equations[3], "|", fixed = TRUE)[[1]][1],"-1",sep="")
			formula2<-paste(equations[2],equations[1],strsplit(equations[3], "|", fixed = TRUE)[[1]][2],sep="")
			mf1 <- model.frame(formula=as.formula(formula1), data=newdata, na.action = na.omit)
			mf2 <- model.frame(formula=as.formula(formula2), data=newdata, na.action = na.omit)
			x <- model.matrix(attr(mf1, "terms"), data=mf1)
			z <- model.matrix(attr(mf2, "terms"), data=mf2)
			ZG<-as.matrix(z) %*% object$coefficients$inflation
			XB<-as.matrix(x) %*% object$coefficients$ordered
			cprobs<-matrix(nrow=object$y.categories-1,ncol=1)
			n <- nrow(object$z)
			probs<-matrix(nrow=n,ncol=object$y.categories)

			cprobs[1,1]<-object$coefficients$cutpoints[1]
			for(j in 2:(object$y.categories-1)){

				cprobs[j,1]<-cprobs[j-1,1]+exp(object$coefficients$cutpoints[j])
			}

			probs[,object$y.categories]<-(pnorm(ZG))*((1-pnorm(cprobs[object$y.categories-1,1]-XB)))
			probs[,1]<-((1-pnorm(ZG))+(pnorm(ZG))*(pnorm(cprobs[1,1]-XB)))

			for(j in 2:(object$y.categories-1)){

				probs[,j]<-(pnorm(ZG))*((pnorm(cprobs[j,1]-XB)) - (pnorm(cprobs[j-1,1]-XB)))
			}

		}
		probs<-cbind(sort(unique(object$y))[apply(probs, 1, which.max)])
		colnames(probs)<-c("Fitted.Y")
	}
  else{
	if(type=="response.full"){
		if(is.null(newdata)){
			ZG<-as.matrix(object$z) %*% object$coefficients$inflation
			XB<-as.matrix(object$x) %*% object$coefficients$ordered
			cprobs<-matrix(nrow=object$y.categories-1,ncol=1)
			n <- nrow(object$z)
			probs<-matrix(nrow=n,ncol=object$y.categories)

			cprobs[1,1]<-object$coefficients$cutpoints[1]
			for(j in 2:(object$y.categories-1)){

				cprobs[j,1]<-cprobs[j-1,1]+exp(object$coefficients$cutpoints[j])
			}

			probs[,object$y.categories]<-(pnorm(ZG))*((1-pnorm(cprobs[object$y.categories-1,1]-XB)))
			probs[,1]<-((1-pnorm(ZG))+(pnorm(ZG))*(pnorm(cprobs[1,1]-XB)))

			for(j in 2:(object$y.categories-1)){

			probs[,j]<-(pnorm(ZG))*((pnorm(cprobs[j,1]-XB)) - (pnorm(cprobs[j-1,1]-XB)))
			}

		}
 		else{
			equations<-as.character(object$formula)
			formula1<-paste(equations[2],equations[1],strsplit(equations[3], "|", fixed = TRUE)[[1]][1],"-1",sep="")
			formula2<-paste(equations[2],equations[1],strsplit(equations[3], "|", fixed = TRUE)[[1]][2],sep="")
			mf1 <- model.frame(formula=as.formula(formula1), data=newdata, na.action = na.omit)
			mf2 <- model.frame(formula=as.formula(formula2), data=newdata, na.action = na.omit)
			x <- model.matrix(attr(mf1, "terms"), data=mf1)
			z <- model.matrix(attr(mf2, "terms"), data=mf2)
			ZG<-as.matrix(z) %*% object$coefficients$inflation
			XB<-as.matrix(x) %*% object$coefficients$ordered
			cprobs<-matrix(nrow=object$y.categories-1,ncol=1)
			n <- nrow(object$z)
			probs<-matrix(nrow=n,ncol=object$y.categories)

			cprobs[1,1]<-object$coefficients$cutpoints[1]
			for(j in 2:(object$y.categories-1)){

				cprobs[j,1]<-cprobs[j-1,1]+exp(object$coefficients$cutpoints[j])
			}

			probs[,object$y.categories]<-(pnorm(ZG))*((1-pnorm(cprobs[object$y.categories-1,1]-XB)))
			probs[,1]<-((1-pnorm(ZG))+(pnorm(ZG))*(pnorm(cprobs[1,1]-XB)))

			for(j in 2:(object$y.categories-1)){

				probs[,j]<-(pnorm(ZG))*((pnorm(cprobs[j,1]-XB)) - (pnorm(cprobs[j-1,1]-XB)))
			}

		}
		probs<-cbind(sort(unique(object$y))[apply(probs, 1, which.max)])
		colnames(probs)<-c("Fitted.Y")
	}
	if(type=="response.ordered"){
		if(is.null(newdata)){
			XB<-as.matrix(object$x) %*% object$coefficients$ordered
			cprobs<-matrix(nrow=object$y.categories-1,ncol=1)
			n <- nrow(object$z)
			probs<-matrix(nrow=n,ncol=object$y.categories)

			cprobs[1,1]<-object$coefficients$cutpoints[1]
			for(j in 2:(object$y.categories-1)){

				cprobs[j,1]<-cprobs[j-1,1]+exp(object$coefficients$cutpoints[j])
			}

			probs[,object$y.categories]<-((1-pnorm(cprobs[object$y.categories-1,1]-XB)))
			probs[,1]<-((pnorm(cprobs[1,1]-XB)))

			for(j in 2:(object$y.categories-1)){

			probs[,j]<-((pnorm(cprobs[j,1]-XB)) - (pnorm(cprobs[j-1,1]-XB)))
			}

		}
 		else{
			equations<-as.character(object$formula)
			formula1<-paste(equations[2],equations[1],strsplit(equations[3], "|", fixed = TRUE)[[1]][1],"-1",sep="")
			formula2<-paste(equations[2],equations[1],strsplit(equations[3], "|", fixed = TRUE)[[1]][2],sep="")
			mf1 <- model.frame(formula=as.formula(formula1), data=newdata, na.action = na.omit)
			mf2 <- model.frame(formula=as.formula(formula2), data=newdata, na.action = na.omit)
			x <- model.matrix(attr(mf1, "terms"), data=mf1)
			z <- model.matrix(attr(mf2, "terms"), data=mf2)
			XB<-as.matrix(x) %*% object$coefficients$ordered
			cprobs<-matrix(nrow=object$y.categories-1,ncol=1)
			n <- nrow(object$z)
			probs<-matrix(nrow=n,ncol=object$y.categories)

			cprobs[1,1]<-object$coefficients$cutpoints[1]
			for(j in 2:(object$y.categories-1)){

				cprobs[j,1]<-cprobs[j-1,1]+exp(object$coefficients$cutpoints[j])
			}

			probs[,object$y.categories]<-((1-pnorm(cprobs[object$y.categories-1,1]-XB)))
			probs[,1]<-((pnorm(cprobs[1,1]-XB)))

			for(j in 2:(object$y.categories-1)){

				probs[,j]<-((pnorm(cprobs[j,1]-XB)) - (pnorm(cprobs[j-1,1]-XB)))
			}

		}
		probs<-cbind(sort(unique(object$y))[apply(probs, 1, which.max)])
		colnames(probs)<-c("Fitted.Y")
	}
	if(type=="response.inflation"){
		if(is.null(newdata)){
			ZG<-as.matrix(object$z) %*% object$coefficients$inflation
			cprobs<-matrix(nrow=object$y.categories-1,ncol=1)
			n <- nrow(object$z)
			probs<-matrix(nrow=n,ncol=1)

			probs[,1]<-1-pnorm(ZG)

			colnames(probs)<-c("Pr(Inflation)")
		}
 		else{
			equations<-as.character(object$formula)
			formula1<-paste(equations[2],equations[1],strsplit(equations[3], "|", fixed = TRUE)[[1]][1],"-1",sep="")
			formula2<-paste(equations[2],equations[1],strsplit(equations[3], "|", fixed = TRUE)[[1]][2],sep="")
			mf1 <- model.frame(formula=as.formula(formula1), data=newdata, na.action = na.omit)
			mf2 <- model.frame(formula=as.formula(formula2), data=newdata, na.action = na.omit)
			x <- model.matrix(attr(mf1, "terms"), data=mf1)
			z <- model.matrix(attr(mf2, "terms"), data=mf2)
			ZG<-as.matrix(z) %*% object$coefficients$inflation
			n <- nrow(object$z)
			probs<-matrix(nrow=n,ncol=1)

			probs[,1]<-1-pnorm(ZG)
			colnames(probs)<-c("Pr(Inflation)")
		}
	probs<-ifelse(probs>.5,1,0)
	colnames(probs)<-c("Y=Inflated")
	}
	if(type=="linear"){
		if(is.null(newdata)){
			ZG<-as.matrix(object$z) %*% object$coefficients$inflation
			XB<-as.matrix(object$x) %*% object$coefficients$ordered
			probs<-cbind(ZG,XB)

			colnames(probs)<-c("ZG","XB")
		}
 		else{
			equations<-as.character(object$formula)
			formula1<-paste(equations[2],equations[1],strsplit(equations[3], "|", fixed = TRUE)[[1]][1],"-1",sep="")
			formula2<-paste(equations[2],equations[1],strsplit(equations[3], "|", fixed = TRUE)[[1]][2],sep="")
			mf1 <- model.frame(formula=as.formula(formula1), data=newdata, na.action = na.omit)
			mf2 <- model.frame(formula=as.formula(formula2), data=newdata, na.action = na.omit)
			x <- model.matrix(attr(mf1, "terms"), data=mf1)
			z <- model.matrix(attr(mf2, "terms"), data=mf2)
			ZG<-as.matrix(z) %*% object$coefficients$inflation
			XB<-as.matrix(x) %*% object$coefficients$ordered
			probs<-cbind(ZG,XB)

			colnames(probs)<-c("ZG","XB")
		}
	}}
	probs
}



#' @title predict
#' @description A function to extract predicted values from the IOP model object.
#'
#' @param object an object of class \code{iop.mod} (output of \code{\link{iop}}).
#' @param newdata An optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used.
#' @param type
#'
#'
#' @examples
#'  model1 <- iop(Y ~ X1 + X2|Z1 + Z2, data=data, type=c('ziop'))
#'  predict(model1, type=c("response.full"))
#'
#' @export
#'
#'
#'
predict <- function(object, newdata=NULL, type = c("prob.full", "prob.ordered", "prob.inflation","response.full","response.ordered", "response.inflation","linear"), ...)
	{
	if(missing(type)){
		if(is.null(newdata)){
			ZG<-as.matrix(object$z) %*% object$coefficients$inflation
			XB<-as.matrix(object$x) %*% object$coefficients$ordered
			cprobs<-matrix(nrow=object$y.categories-1,ncol=1)
			n <- nrow(object$z)
			probs<-matrix(nrow=n,ncol=object$y.categories)

			cprobs[1,1]<-object$coefficients$cutpoints[1]
			for(j in 2:(object$y.categories-1)){

				cprobs[j,1]<-cprobs[j-1,1]+exp(object$coefficients$cutpoints[j])
			}

			probs[,object$y.categories]<-(pnorm(ZG))*((1-pnorm(cprobs[object$y.categories-1,1]-XB)))
			probs[,1]<-((1-pnorm(ZG))+(pnorm(ZG))*(pnorm(cprobs[1,1]-XB)))

			for(j in 2:(object$y.categories-1)){

			probs[,j]<-(pnorm(ZG))*((pnorm(cprobs[j,1]-XB)) - (pnorm(cprobs[j-1,1]-XB)))
			}

			colnames(probs)<-paste("Pr(Y=",sort(unique(object$y)),")",sep="")
		}
 		else{
			equations<-as.character(object$formula)
			formula1<-paste(equations[2],equations[1],strsplit(equations[3], "|", fixed = TRUE)[[1]][1],"-1",sep="")
			formula2<-paste(equations[2],equations[1],strsplit(equations[3], "|", fixed = TRUE)[[1]][2],sep="")
			mf1 <- model.frame(formula=as.formula(formula1), data=newdata, na.action = na.omit)
			mf2 <- model.frame(formula=as.formula(formula2), data=newdata, na.action = na.omit)
			x <- model.matrix(attr(mf1, "terms"), data=mf1)
			z <- model.matrix(attr(mf2, "terms"), data=mf2)
			ZG<-as.matrix(z) %*% object$coefficients$inflation
			XB<-as.matrix(x) %*% object$coefficients$ordered
			cprobs<-matrix(nrow=object$y.categories-1,ncol=1)
			n <- nrow(object$z)
			probs<-matrix(nrow=n,ncol=object$y.categories)

			cprobs[1,1]<-object$coefficients$cutpoints[1]
			for(j in 2:(object$y.categories-1)){

				cprobs[j,1]<-cprobs[j-1,1]+exp(object$coefficients$cutpoints[j])
			}

			probs[,object$y.categories]<-(pnorm(ZG))*((1-pnorm(cprobs[object$y.categories-1,1]-XB)))
			probs[,1]<-((1-pnorm(ZG))+(pnorm(ZG))*(pnorm(cprobs[1,1]-XB)))

			for(j in 2:(object$y.categories-1)){

				probs[,j]<-(pnorm(ZG))*((pnorm(cprobs[j,1]-XB)) - (pnorm(cprobs[j-1,1]-XB)))
			}

			colnames(probs)<-paste("Pr(Y=",sort(unique(object$y)),")",sep="")
		}
	}else{
	if(type=="prob.full"){
		if(is.null(newdata)){
			ZG<-as.matrix(object$z) %*% object$coefficients$inflation
			XB<-as.matrix(object$x) %*% object$coefficients$ordered
			cprobs<-matrix(nrow=object$y.categories-1,ncol=1)
			n <- nrow(object$z)
			probs<-matrix(nrow=n,ncol=object$y.categories)

			cprobs[1,1]<-object$coefficients$cutpoints[1]
			for(j in 2:(object$y.categories-1)){

				cprobs[j,1]<-cprobs[j-1,1]+exp(object$coefficients$cutpoints[j])
			}

			probs[,object$y.categories]<-(pnorm(ZG))*((1-pnorm(cprobs[object$y.categories-1,1]-XB)))
			probs[,1]<-((1-pnorm(ZG))+(pnorm(ZG))*(pnorm(cprobs[1,1]-XB)))

			for(j in 2:(object$y.categories-1)){

			probs[,j]<-(pnorm(ZG))*((pnorm(cprobs[j,1]-XB)) - (pnorm(cprobs[j-1,1]-XB)))
			}

			colnames(probs)<-paste("Pr(Y=",sort(unique(object$y)),")",sep="")
		}
 		else{
			equations<-as.character(object$formula)
			formula1<-paste(equations[2],equations[1],strsplit(equations[3], "|", fixed = TRUE)[[1]][1],"-1",sep="")
			formula2<-paste(equations[2],equations[1],strsplit(equations[3], "|", fixed = TRUE)[[1]][2],sep="")
			mf1 <- model.frame(formula=as.formula(formula1), data=newdata, na.action = na.omit)
			mf2 <- model.frame(formula=as.formula(formula2), data=newdata, na.action = na.omit)
			x <- model.matrix(attr(mf1, "terms"), data=mf1)
			z <- model.matrix(attr(mf2, "terms"), data=mf2)
			ZG<-as.matrix(z) %*% object$coefficients$inflation
			XB<-as.matrix(x) %*% object$coefficients$ordered
			cprobs<-matrix(nrow=object$y.categories-1,ncol=1)
			n <- nrow(object$z)
			probs<-matrix(nrow=n,ncol=object$y.categories)

			cprobs[1,1]<-object$coefficients$cutpoints[1]
			for(j in 2:(object$y.categories-1)){

				cprobs[j,1]<-cprobs[j-1,1]+exp(object$coefficients$cutpoints[j])
			}

			probs[,object$y.categories]<-(pnorm(ZG))*((1-pnorm(cprobs[object$y.categories-1,1]-XB)))
			probs[,1]<-((1-pnorm(ZG))+(pnorm(ZG))*(pnorm(cprobs[1,1]-XB)))

			for(j in 2:(object$y.categories-1)){

				probs[,j]<-(pnorm(ZG))*((pnorm(cprobs[j,1]-XB)) - (pnorm(cprobs[j-1,1]-XB)))
			}

			colnames(probs)<-paste("Pr(Y=",sort(unique(object$y)),")",sep="")
		}
	}
	if(type=="prob.ordered"){
		if(is.null(newdata)){
			XB<-as.matrix(object$x) %*% object$coefficients$ordered
			cprobs<-matrix(nrow=object$y.categories-1,ncol=1)
			n <- nrow(object$z)
			probs<-matrix(nrow=n,ncol=object$y.categories)

			cprobs[1,1]<-object$coefficients$cutpoints[1]
			for(j in 2:(object$y.categories-1)){

				cprobs[j,1]<-cprobs[j-1,1]+exp(object$coefficients$cutpoints[j])
			}

			probs[,object$y.categories]<-((1-pnorm(cprobs[object$y.categories-1,1]-XB)))
			probs[,1]<-((pnorm(cprobs[1,1]-XB)))

			for(j in 2:(object$y.categories-1)){

			probs[,j]<-((pnorm(cprobs[j,1]-XB)) - (pnorm(cprobs[j-1,1]-XB)))
			}

			colnames(probs)<-paste("Pr(Y=",sort(unique(object$y)),")",sep="")
		}
 		else{
			equations<-as.character(object$formula)
			formula1<-paste(equations[2],equations[1],strsplit(equations[3], "|", fixed = TRUE)[[1]][1],"-1",sep="")
			formula2<-paste(equations[2],equations[1],strsplit(equations[3], "|", fixed = TRUE)[[1]][2],sep="")
			mf1 <- model.frame(formula=as.formula(formula1), data=newdata, na.action = na.omit)
			mf2 <- model.frame(formula=as.formula(formula2), data=newdata, na.action = na.omit)
			x <- model.matrix(attr(mf1, "terms"), data=mf1)
			z <- model.matrix(attr(mf2, "terms"), data=mf2)
			XB<-as.matrix(x) %*% object$coefficients$ordered
			cprobs<-matrix(nrow=object$y.categories-1,ncol=1)
			n <- nrow(object$z)
			probs<-matrix(nrow=n,ncol=object$y.categories)

			cprobs[1,1]<-object$coefficients$cutpoints[1]
			for(j in 2:(object$y.categories-1)){

				cprobs[j,1]<-cprobs[j-1,1]+exp(object$coefficients$cutpoints[j])
			}

			probs[,object$y.categories]<-((1-pnorm(cprobs[object$y.categories-1,1]-XB)))
			probs[,1]<-((pnorm(cprobs[1,1]-XB)))

			for(j in 2:(object$y.categories-1)){

				probs[,j]<-((pnorm(cprobs[j,1]-XB)) - (pnorm(cprobs[j-1,1]-XB)))
			}

			colnames(probs)<-paste("Pr(Y=",sort(unique(object$y)),")",sep="")
		}
	}
	if(type=="prob.inflation"){
		if(is.null(newdata)){
			ZG<-as.matrix(object$z) %*% object$coefficients$inflation
			cprobs<-matrix(nrow=object$y.categories-1,ncol=1)
			n <- nrow(object$z)
			probs<-matrix(nrow=n,ncol=1)

			probs[,1]<-1-pnorm(ZG)

			colnames(probs)<-c("Pr(Inflation)")
		}
 		else{
			equations<-as.character(object$formula)
			formula1<-paste(equations[2],equations[1],strsplit(equations[3], "|", fixed = TRUE)[[1]][1],"-1",sep="")
			formula2<-paste(equations[2],equations[1],strsplit(equations[3], "|", fixed = TRUE)[[1]][2],sep="")
			mf1 <- model.frame(formula=as.formula(formula1), data=newdata, na.action = na.omit)
			mf2 <- model.frame(formula=as.formula(formula2), data=newdata, na.action = na.omit)
			x <- model.matrix(attr(mf1, "terms"), data=mf1)
			z <- model.matrix(attr(mf2, "terms"), data=mf2)
			ZG<-as.matrix(z) %*% object$coefficients$inflation
			n <- nrow(object$z)
			probs<-matrix(nrow=n,ncol=1)

			probs[,1]<-1-pnorm(ZG)
			colnames(probs)<-c("Pr(Inflation)")
		}
	}
	if(type=="response.full"){
		if(is.null(newdata)){
			ZG<-as.matrix(object$z) %*% object$coefficients$inflation
			XB<-as.matrix(object$x) %*% object$coefficients$ordered
			cprobs<-matrix(nrow=object$y.categories-1,ncol=1)
			n <- nrow(object$z)
			probs<-matrix(nrow=n,ncol=object$y.categories)

			cprobs[1,1]<-object$coefficients$cutpoints[1]
			for(j in 2:(object$y.categories-1)){

				cprobs[j,1]<-cprobs[j-1,1]+exp(object$coefficients$cutpoints[j])
			}

			probs[,object$y.categories]<-(pnorm(ZG))*((1-pnorm(cprobs[object$y.categories-1,1]-XB)))
			probs[,1]<-((1-pnorm(ZG))+(pnorm(ZG))*(pnorm(cprobs[1,1]-XB)))

			for(j in 2:(object$y.categories-1)){

			probs[,j]<-(pnorm(ZG))*((pnorm(cprobs[j,1]-XB)) - (pnorm(cprobs[j-1,1]-XB)))
			}

		}
 		else{
			equations<-as.character(object$formula)
			formula1<-paste(equations[2],equations[1],strsplit(equations[3], "|", fixed = TRUE)[[1]][1],"-1",sep="")
			formula2<-paste(equations[2],equations[1],strsplit(equations[3], "|", fixed = TRUE)[[1]][2],sep="")
			mf1 <- model.frame(formula=as.formula(formula1), data=newdata, na.action = na.omit)
			mf2 <- model.frame(formula=as.formula(formula2), data=newdata, na.action = na.omit)
			x <- model.matrix(attr(mf1, "terms"), data=mf1)
			z <- model.matrix(attr(mf2, "terms"), data=mf2)
			ZG<-as.matrix(z) %*% object$coefficients$inflation
			XB<-as.matrix(x) %*% object$coefficients$ordered
			cprobs<-matrix(nrow=object$y.categories-1,ncol=1)
			n <- nrow(object$z)
			probs<-matrix(nrow=n,ncol=object$y.categories)

			cprobs[1,1]<-object$coefficients$cutpoints[1]
			for(j in 2:(object$y.categories-1)){

				cprobs[j,1]<-cprobs[j-1,1]+exp(object$coefficients$cutpoints[j])
			}

			probs[,object$y.categories]<-(pnorm(ZG))*((1-pnorm(cprobs[object$y.categories-1,1]-XB)))
			probs[,1]<-((1-pnorm(ZG))+(pnorm(ZG))*(pnorm(cprobs[1,1]-XB)))

			for(j in 2:(object$y.categories-1)){

				probs[,j]<-(pnorm(ZG))*((pnorm(cprobs[j,1]-XB)) - (pnorm(cprobs[j-1,1]-XB)))
			}

		}
		probs<-cbind(sort(unique(object$y))[apply(probs, 1, which.max)])
		colnames(probs)<-c("Fitted.Y")
	}
	if(type=="response.ordered"){
		if(is.null(newdata)){
			XB<-as.matrix(object$x) %*% object$coefficients$ordered
			cprobs<-matrix(nrow=object$y.categories-1,ncol=1)
			n <- nrow(object$z)
			probs<-matrix(nrow=n,ncol=object$y.categories)

			cprobs[1,1]<-object$coefficients$cutpoints[1]
			for(j in 2:(object$y.categories-1)){

				cprobs[j,1]<-cprobs[j-1,1]+exp(object$coefficients$cutpoints[j])
			}

			probs[,object$y.categories]<-((1-pnorm(cprobs[object$y.categories-1,1]-XB)))
			probs[,1]<-((pnorm(cprobs[1,1]-XB)))

			for(j in 2:(object$y.categories-1)){

			probs[,j]<-((pnorm(cprobs[j,1]-XB)) - (pnorm(cprobs[j-1,1]-XB)))
			}

		}
 		else{
			equations<-as.character(object$formula)
			formula1<-paste(equations[2],equations[1],strsplit(equations[3], "|", fixed = TRUE)[[1]][1],"-1",sep="")
			formula2<-paste(equations[2],equations[1],strsplit(equations[3], "|", fixed = TRUE)[[1]][2],sep="")
			mf1 <- model.frame(formula=as.formula(formula1), data=newdata, na.action = na.omit)
			mf2 <- model.frame(formula=as.formula(formula2), data=newdata, na.action = na.omit)
			x <- model.matrix(attr(mf1, "terms"), data=mf1)
			z <- model.matrix(attr(mf2, "terms"), data=mf2)
			XB<-as.matrix(x) %*% object$coefficients$ordered
			cprobs<-matrix(nrow=object$y.categories-1,ncol=1)
			n <- nrow(object$z)
			probs<-matrix(nrow=n,ncol=object$y.categories)

			cprobs[1,1]<-object$coefficients$cutpoints[1]
			for(j in 2:(object$y.categories-1)){

				cprobs[j,1]<-cprobs[j-1,1]+exp(object$coefficients$cutpoints[j])
			}

			probs[,object$y.categories]<-((1-pnorm(cprobs[object$y.categories-1,1]-XB)))
			probs[,1]<-((pnorm(cprobs[1,1]-XB)))

			for(j in 2:(object$y.categories-1)){

				probs[,j]<-((pnorm(cprobs[j,1]-XB)) - (pnorm(cprobs[j-1,1]-XB)))
			}

		}
		probs<-cbind(sort(unique(object$y))[apply(probs, 1, which.max)])
		colnames(probs)<-c("Fitted.Y")
	}
	if(type=="response.inflation"){
		if(is.null(newdata)){
			ZG<-as.matrix(object$z) %*% object$coefficients$inflation
			cprobs<-matrix(nrow=object$y.categories-1,ncol=1)
			n <- nrow(object$z)
			probs<-matrix(nrow=n,ncol=1)

			probs[,1]<-1-pnorm(ZG)

			colnames(probs)<-c("Pr(Inflation)")
		}
 		else{
			equations<-as.character(object$formula)
			formula1<-paste(equations[2],equations[1],strsplit(equations[3], "|", fixed = TRUE)[[1]][1],"-1",sep="")
			formula2<-paste(equations[2],equations[1],strsplit(equations[3], "|", fixed = TRUE)[[1]][2],sep="")
			mf1 <- model.frame(formula=as.formula(formula1), data=newdata, na.action = na.omit)
			mf2 <- model.frame(formula=as.formula(formula2), data=newdata, na.action = na.omit)
			x <- model.matrix(attr(mf1, "terms"), data=mf1)
			z <- model.matrix(attr(mf2, "terms"), data=mf2)
			ZG<-as.matrix(z) %*% object$coefficients$inflation
			n <- nrow(object$z)
			probs<-matrix(nrow=n,ncol=1)

			probs[,1]<-1-pnorm(ZG)
			colnames(probs)<-c("Pr(Inflation)")
		}
	probs<-ifelse(probs>.5,1,0)
	colnames(probs)<-c("Y=Inflated")
	}
	if(type=="linear"){
		if(is.null(newdata)){
			ZG<-as.matrix(object$z) %*% object$coefficients$inflation
			XB<-as.matrix(object$x) %*% object$coefficients$ordered
			probs<-cbind(ZG,XB)

			colnames(probs)<-c("ZG","XB")
		}
 		else{
			equations<-as.character(object$formula)
			formula1<-paste(equations[2],equations[1],strsplit(equations[3], "|", fixed = TRUE)[[1]][1],"-1",sep="")
			formula2<-paste(equations[2],equations[1],strsplit(equations[3], "|", fixed = TRUE)[[1]][2],sep="")
			mf1 <- model.frame(formula=as.formula(formula1), data=newdata, na.action = na.omit)
			mf2 <- model.frame(formula=as.formula(formula2), data=newdata, na.action = na.omit)
			x <- model.matrix(attr(mf1, "terms"), data=mf1)
			z <- model.matrix(attr(mf2, "terms"), data=mf2)
			ZG<-as.matrix(z) %*% object$coefficients$inflation
			XB<-as.matrix(x) %*% object$coefficients$ordered
			probs<-cbind(ZG,XB)

			colnames(probs)<-c("ZG","XB")
		}
	}}
	probs
}


#' @title coef
#' @description A function to extract coefficients from \code{iop} model results.
#'
#' @param object an object of class \code{iop.mod} (output of \code{\link{iop}}).
#'
#'
#' @examples
#'  model1 <- iop(Y ~ X1 + X2|Z1 + Z2, data=data, type=c('ziop'))
#'  coef(model1)
#'
#' @export
#'
#'
#expanded object-oriented function for extracting coefficients...
coef<-function (object, model = c("full", "cutpoints", "inflation","ordered"), ...){
		model <- match.arg(model)
		rval <- object$coefficients
		rval <- switch(model, full = structure(c(rval$cutpoints, rval$inflation, rval$ordered),
		.Names = c(names(rval$cutpoints),names(rval$inflation),paste("ordered", names(rval$ordered), sep = "."))),
		cutpoints= rval$cutpoints,inflation = rval$inflation, ordered = rval$ordered)
    	rval
}


#' @title residuals
#' @description A function to extract residuals from \code{iop} model results.
#'
#' @param object an object of class \code{iop.mod} (output of \code{\link{iop}}).
#'
#'
#' @examples
#'  model1 <- iop(Y ~ X1 + X2|Z1 + Z2, data=data, type=c('ziop'))
#'  residuals(model1)
#'
#' @export
#'
#'
#expanded function for residuals
residuals<-function (object, type = c("response"), ...)
	{
		type <- match.arg(type)
		res <- as.data.frame(object$y)-predict(object,type="response.full")
		switch(type, response = {
	return(res)
	})
}


#' @title vcov
#' @description A function to extract variance-covariance matrix from \code{iop} model results.
#'
#' @param object an object of class \code{iop.mod} (output of \code{\link{iop}}).
#'
#'
#' @examples
#'  model1 <- iop(Y ~ X1 + X2|Z1 + Z2, data=data, type=c('ziop'))
#'  vcov(model1)
#'
#' @export
#'
#'
#expanded object-oriented function for vcov
vcov<-function (object, model = c("full", "cutpoints","inflation", "ordered"), ...)
	{
		model <- match.arg(model)
		rval <- object$vcov
		if(model == "full"){
 			return(rval)
			}
		else{
			if(model=="cutpoints"){
				cf <- object$coefficients[[model]]
				wi <- seq(along = object$coefficients$cutpoints)
				rval <- if (model == "cutpoints")
				rval[wi, wi]
				colnames(rval) <- rownames(rval) <- names(cf)
			}
			if(model=="inflation"){
				cf <- object$coefficients[[model]]
				wi <- seq(along = object$coefficients$inflation)
				rval <- if (model == "inflation")
				rval[wi, wi]
				colnames(rval) <- rownames(rval) <- names(cf)
			}
			if(model=="ordered"){
				cf <- object$coefficients[[model]]
				wi <- seq(along = object$coefficients$ordered)
				rval <- if (model == "ordered")
				rval[wi, wi]
				colnames(rval) <- rownames(rval) <- names(cf)
			}
		}
    return(rval)
}

