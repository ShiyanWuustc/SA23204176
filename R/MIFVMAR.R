#' @title MIFVMA Method
#' @name MIFVMA
#' @description Using forward-validation model averaging method when given continuous mixed data
#' @param y the continuous response variable
#' @param x the independent variable with mixed frequency
#' @param per the initial time point in FVMA
#' @param k the largest lag 
#' @return a list contain model averaging weight and one-step prediction
#' @examples
#' \dontrun{
#'    theta1=7*10^(-4)
#'    theta2=-6*10^(-3)
#'    n=200
#'    x.sim<-matrix(arima.sim(model=list(ar=0.5), n=k*(n+3)), ncol=k,byrow = TRUE)
#'    x<-cbind(x.sim[1:(n+2),],x.sim[2:(n+3),])
#'    beta<-exp(theta1*seq(2*k)+theta2*seq(2*k)^2)/sum(exp(theta1*seq(2*k)+theta2*seq(2*k)^2))
#'    y=rep(0,n+2)
#'    for( i in 2:(n+2)) {
#'      y[i]=0.5+-0.2*y[i-1] + 1.7*x[i,]%*%beta
#'    }
#'    MIFVMA(y=y,x=x,per=0.7,k=5)
#' }
#' @export
MIFVMA<-function(y,x,per,k){   
  n<-length(y)-2
  t2<-0.8*n
  t3<-0.6*n
  t4<-0.4*n    
  tmin<-min(t2,t3,t4)
  v=n-per*tmin
  e<-list(n-v)
  for(j in (v+1):n){
    x_train1<-lapply(1:k, function(i) cbind(y[1:(j-1)],       x[1:(j-1)       ,1:i]))
    x_train2<-lapply(1:k, function(i) cbind(y[(n-t2+1):(j-1)],x[(n-t2+1):(j-1),1:i]))
    x_train3<-lapply(1:k, function(i) cbind(y[(n-t3+1):(j-1)],x[(n-t3+1):(j-1),1:i]))
    x_train4<-lapply(1:k, function(i) cbind(y[(n-t4+1):(j-1)],x[(n-t4+1):(j-1),1:i]))
    x_test<-lapply(1:k, function(i) cbind(y[j],t(x[j,1:i])))
    coef1<-lapply(1:k, function(i) glm(y[2:j]       ~x_train1[[i]],family = gaussian,control = list(maxit=200))$coef)
    coef2<-lapply(1:k, function(i) glm(y[(n-t2+2):j]~x_train2[[i]],family = gaussian,control = list(maxit=200))$coef)
    coef3<-lapply(1:k, function(i) glm(y[(n-t3+2):j]~x_train3[[i]],family = gaussian,control = list(maxit=200))$coef)
    coef4<-lapply(1:k, function(i) glm(y[(n-t4+2):j]~x_train4[[i]],family = gaussian,control = list(maxit=200))$coef)
    e[[j-v]]<-c(unlist(lapply(1:k, function(i) y[(j+1)]-c(1,x_test[[i]])%*%coef1[[i]]))
                ,unlist(lapply(1:k, function(i) y[(j+1)]-c(1,x_test[[i]])%*%coef2[[i]]))
                ,unlist(lapply(1:k, function(i) y[(j+1)]-c(1,x_test[[i]])%*%coef3[[i]]))
                ,unlist(lapply(1:k, function(i) y[(j+1)]-c(1,x_test[[i]])%*%coef4[[i]])))
  } 
  
  e_tlide<-do.call(rbind,e)
  s<-t(e_tlide)%*%e_tlide/(n-v)
  Dmat<-2*s
  l=ncol(s)
  dvec<-rep(0,l)
  Amat<-t(rbind(rep(1,l),diag(x=1,l,l),diag(x=-1,l,l)))  
  bvec<-c(1,rep(0,l),rep(-1,l))
  
  w_hat<-solve.QP(Dmat, dvec, Amat, bvec, meq=1)$solution   
  
  m.candi1<-lapply(1:k, function(i) glm(y[2:(n+1)]       ~cbind(y[1:n],x[1:n,1:i]),family = gaussian,control = list(maxit=200)))
  m.candi2<-lapply(1:k, function(i) glm(y[(n-t2+2):(n+1)]~cbind(y[(n-t2+1):n],x[(n-t2+1):n,1:i]),family = gaussian,control = list(maxit=200)))
  m.candi3<-lapply(1:k, function(i) glm(y[(n-t3+2):(n+1)]~cbind(y[(n-t3+1):n],x[(n-t3+1):n,1:i]),family = gaussian,control = list(maxit=200)))
  m.candi4<-lapply(1:k, function(i) glm(y[(n-t4+2):(n+1)]~cbind(y[(n-t4+1):n],x[(n-t4+1):n,1:i]),family = gaussian,control = list(maxit=200)))
  
  x_pred<-lapply(1:k, function(i) cbind(y[n+1],t(x[n+1,1:i])))  
  y_hat_i<-c(unlist(lapply(1:k,  function(i) c(1,x_pred[[i]])%*%m.candi1[[i]]$coef))  
             ,unlist(lapply(1:k, function(i) c(1,x_pred[[i]])%*%m.candi2[[i]]$coef))
             ,unlist(lapply(1:k, function(i) c(1,x_pred[[i]])%*%m.candi3[[i]]$coef))
             ,unlist(lapply(1:k, function(i) c(1,x_pred[[i]])%*%m.candi4[[i]]$coef)))
  y_hat_w_fvma<-t(y_hat_i)%*%w_hat  
  return(list(prediction=y_hat_w_fvma,weight=w_hat))
}

#' @title Bootstrap in R
#' @name BootstrapR
#' @description Bootstrap in R used for comparing the computation time with Rcpp.
#' @param x data
#' @param B sampling time
#' @return a list of Bootstrap estimator
#' @examples
#' \dontrun{
#' data=rnorm(100)
#' bootsR(data,1000)
#' }
#' @export
bootsR <- function(x,B){
  thetastar <- numeric(B)
  theta <- mean(x)
  for(b in 1:B){
    xstar <- sample(x,replace=TRUE)
    thetastar[b] <- mean(xstar)
  }
  return(c(bias=mean(thetastar)-theta,se.boot=sd(thetastar),se.samp=sd(x)/sqrt(length(x))))
}


f <- function(x, sigma) {
  if (any(x < 0)) return (0)
  stopifnot(sigma > 0)
  return((x / sigma^2) * exp(-x^2 / (2*sigma^2)))
}

#' @title Metropolis-Hastings sampler generating random number of Rayleigh distribution in R
#' @name MCMCMHR
#' @description Metropolis-Hastings sampler generates random number of Rayleigh distribution in R used for comparing the computation time with Rcpp.
#' @param m sample size
#' @param sigma the parameter of Rayleigh distribution 
#' @param b burn-in size
#' @return a list of Bootstrap estimator
#' @examples
#' \dontrun{
#' RayleighR(m=5e3,sigma=1,b=2e3)
#' }
#' @export
RayleighR <- function(m,sigma,b){
  x <- numeric(m)
  x[1] <- rchisq(1, df=1)
  k <- 0
  u <- runif(m)
  
  for (i in 2:m) {
    xt <- x[i-1]
    y <- rchisq(1, df = xt)
    num <- f(y, sigma) * dchisq(xt, df = y)
    den <- f(xt, sigma) * dchisq(y, df = xt)
    if (u[i] <= num/den){
      x[i] <- y
    } else {
      x[i] <- xt
      k <- k+1     #y is rejected
    }
  }
  chain <- x[(b+1):m]
  return(chain)
}

#' @title Benchmark R and Rcpp Bootstraps.
#' @name benchmarks
#' @description Use R package \code{microbenchmark} to compare the performance of Bootstrap and Rayleigh Random number generating in R and Rcpp.
#' @examples
#' \dontrun{
#' data <- rnorm(100)
#' tm1 <- microbenchmark::microbenchmark(
#'   R = bootsR(data, 1000),
#'   C = bootsC(data, 1000)
#' )
#' print(summary(tm1)[,c(1,3,5,6)])
#' 
#' tm2 <- microbenchmark::microbenchmark(
#'   R = RayleighR(5e3,1,2e3),
#'   C = RayleighC(5e3,1,2e3)
#' )
#' print(summary(tm2)[,c(1,3,5,6)])
#' }
#' @import quadprog
#' @import microbenchmark
#' @import gcookbook
#' @import ggplot2
#' @import bootstrap
#' @import boot
#' @import DAAG
#' @import Deriv
#' @import coda
#' @importFrom Rcpp evalCpp
#' @importFrom stats rnorm rgamma rexp runif glm gaussian dchisq rchisq sd
#' @useDynLib SA23204176
NULL