## ----eval=FALSE---------------------------------------------------------------
#  function(y,x,per,k){
#    #y为变量，x为混频数据,k是选取得最大滞后阶数
#    n<-length(y)-2
#    t2<-0.8*n
#    t3<-0.6*n
#    t4<-0.4*n    #选取一般时间间隔，这里可以设置一个参数
#    tmin<-min(t2,t3,t4)
#    v=n-per*tmin
#    e<-list(n-v)
#    for(j in (v+1):n){
#      x_train1<-lapply(1:k, function(i) cbind(y[1:(j-1)],       x[1:(j-1)       ,1:i]))
#      x_train2<-lapply(1:k, function(i) cbind(y[(n-t2+1):(j-1)],x[(n-t2+1):(j-1),1:i]))
#      x_train3<-lapply(1:k, function(i) cbind(y[(n-t3+1):(j-1)],x[(n-t3+1):(j-1),1:i]))
#      x_train4<-lapply(1:k, function(i) cbind(y[(n-t4+1):(j-1)],x[(n-t4+1):(j-1),1:i]))
#      x_test<-lapply(1:k, function(i) cbind(y[j],t(x[j,1:i])))
#      coef1<-lapply(1:k, function(i) glm(y[2:j]       ~x_train1[[i]],family = gaussian,control = list(maxit=200))$coef)
#      coef2<-lapply(1:k, function(i) glm(y[(n-t2+2):j]~x_train2[[i]],family = gaussian,control = list(maxit=200))$coef)
#      coef3<-lapply(1:k, function(i) glm(y[(n-t3+2):j]~x_train3[[i]],family = gaussian,control = list(maxit=200))$coef)
#      coef4<-lapply(1:k, function(i) glm(y[(n-t4+2):j]~x_train4[[i]],family = gaussian,control = list(maxit=200))$coef)
#      e[[j-v]]<-c(unlist(lapply(1:k, function(i) y[(j+1)]-c(1,x_test[[i]])%*%coef1[[i]]))
#                  ,unlist(lapply(1:k, function(i) y[(j+1)]-c(1,x_test[[i]])%*%coef2[[i]]))#使用j时刻的test数据拟合j+1时刻真实数据集
#                  ,unlist(lapply(1:k, function(i) y[(j+1)]-c(1,x_test[[i]])%*%coef3[[i]]))
#                  ,unlist(lapply(1:k, function(i) y[(j+1)]-c(1,x_test[[i]])%*%coef4[[i]])))
#    }
#  
#    e_tlide<-do.call(rbind,e)
#    s<-t(e_tlide)%*%e_tlide/(n-v)
#    Dmat<-2*s
#    l=ncol(s)
#    dvec<-rep(0,l)
#    Amat<-t(rbind(rep(1,l),diag(x=1,l,l),diag(x=-1,l,l)))
#    bvec<-c(1,rep(0,l),rep(-1,l))
#  
#    w_hat<-solve.QP(Dmat, dvec, Amat, bvec, meq=1)$solution   #FVMA权重
#  
#    m.candi1<-lapply(1:k, function(i) glm(y[2:(n+1)]       ~cbind(y[1:n],x[1:n,1:i]),family = gaussian,control = list(maxit=200)))
#    m.candi2<-lapply(1:k, function(i) glm(y[(n-t2+2):(n+1)]~cbind(y[(n-t2+1):n],x[(n-t2+1):n,1:i]),family = gaussian,control = list(maxit=200)))
#    m.candi3<-lapply(1:k, function(i) glm(y[(n-t3+2):(n+1)]~cbind(y[(n-t3+1):n],x[(n-t3+1):n,1:i]),family = gaussian,control = list(maxit=200)))
#    m.candi4<-lapply(1:k, function(i) glm(y[(n-t4+2):(n+1)]~cbind(y[(n-t4+1):n],x[(n-t4+1):n,1:i]),family = gaussian,control = list(maxit=200)))
#  
#    x_pred<-lapply(1:k, function(i) cbind(y[n+1],t(x[n+1,1:i])))  #候选模型list
#    y_hat_i<-c(unlist(lapply(1:k,  function(i) c(1,x_pred[[i]])%*%m.candi1[[i]]$coef))  #使用1:n的数据拟合,得到y_hat_n+1
#               ,unlist(lapply(1:k, function(i) c(1,x_pred[[i]])%*%m.candi2[[i]]$coef))
#               ,unlist(lapply(1:k, function(i) c(1,x_pred[[i]])%*%m.candi3[[i]]$coef))
#               ,unlist(lapply(1:k, function(i) c(1,x_pred[[i]])%*%m.candi4[[i]]$coef)))
#    y_hat_w_fvma<-t(y_hat_i)%*%w_hat  #FVMA加权平均后的y_hat
#    return(list(prediction=y_hat_w_fvma,weight=w_hat))
#  }

## -----------------------------------------------------------------------------
library(SA23204176)
#数据生成：考虑Almon双指数MIDAS模型
theta1=7*10^(-4)
theta2=-6*10^(-3)
n=200
l=4
x.sim<-matrix(arima.sim(model=list(ar=0.5), n=l*(n+3)), ncol=l,byrow = TRUE)
x<-cbind(x.sim[1:(n+2),],x.sim[2:(n+3),])
beta<-exp(theta1*seq(2*l)+theta2*seq(2*l)^2)/sum(exp(theta1*seq(2*l)+theta2*seq(2*l)^2))
y=rep(0,n+2)
for( i in 2:(n+2)) {
  y[i]=0.5+-0.2*y[i-1] + 1.7*x[i,]%*%beta
}
#MIFVMA函数使用
MIFVMA(y=y,x=x,per=0.7,k=5)

## ----eval=FALSE---------------------------------------------------------------
#  NumericVector boots(NumericVector x,int B) {
#    NumericVector thetastar(B);
#    double theta = mean(x);
#    int n = x.size();
#  
#    for(int b = 0; b < B; b++) {
#      NumericVector xstar = Rcpp::sample(x, n, true);
#      thetastar[b] = mean(xstar);
#    }
#  
#    double bias = mean(thetastar) - theta;
#    double se_boot = sd(thetastar);
#    double se_samp = sd(x) / sqrt(n);
#  
#    return NumericVector::create(Named("bias") = bias,
#                                 Named("se.boot") = se_boot,
#                                 Named("se.samp") = se_samp);
#  }

## ----eval=TRUE----------------------------------------------------------------
library(Rcpp)
library(SA23204176)
data <- rnorm(100)
bootstrap_results <-boots(data, 1000)
bootstrap_results

## ----eval=FALSE---------------------------------------------------------------
#  #R函数
#  bootsR <- function(x,B){
#    thetastar <- numeric(B)
#    theta <- mean(x)
#    for(b in 1:B){
#      xstar <- sample(x,replace=TRUE)
#      thetastar[b] <- mean(xstar)
#    }
#    return(c(bias=mean(thetastar)-theta,se.boot=sd(thetastar),se.samp=sd(x)/sqrt(length(x))))
#  }

## -----------------------------------------------------------------------------
library(microbenchmark)
tm1 <- microbenchmark(
  R = bootsR(data, 1000),
  C = boots(data, 1000)
)
data <- rnorm(100)
knitr::kable(summary(tm1)[,c(1,3,5,6)])

## ----eval=FALSE---------------------------------------------------------------
#  #非用户水平函数
#  double f(double x, double sigma) {
#    if (x < 0) return 0;
#    return (x / pow(sigma, 2)) * exp(-pow(x, 2) / (2 * pow(sigma, 2)));
#  }
#  
#  #用户水平函数
#   NumericVector Rayleigh(int m, double sigma, int b) {
#     sigma = 4;
#     NumericVector x(m);
#     x[0] = R::rchisq(1);
#     int k = 0;
#     NumericVector u = runif(m);
#  
#     Environment stats("package:stats");
#     Function dchisq = stats["dchisq"];
#     Function rchisq = stats["rchisq"];
#  
#     for (int i = 1; i < m; ++i) {
#       double xt = x[i - 1];
#       double y = R::rchisq(1);
#       double num = f(y, sigma) * as<double>(dchisq(xt, y));
#       double den = f(xt, sigma) * as<double>(dchisq(y, xt));
#  
#       if (u[i] <= num / den) {
#         x[i] = y;
#       } else {
#         x[i] = xt;
#         k++; // y is rejected
#       }
#     }
#     return x[Range(b, m - 1)];
#   }

## -----------------------------------------------------------------------------
RayRandom=Rayleigh(m=2000,sigma=1,b=100)
head(RayRandom)

## ----eval=FALSE---------------------------------------------------------------
#  f <- function(x, sigma) {
#    if (any(x < 0)) return (0)
#    stopifnot(sigma > 0)
#    return((x / sigma^2) * exp(-x^2 / (2*sigma^2)))
#  }
#  RayleighR <- function(m,sigma,b){
#      sigma <- 4
#      x <- numeric(m)
#      x[1] <- rchisq(1, df=1)
#      k <- 0
#      u <- runif(m)
#  
#      for (i in 2:m) {
#          xt <- x[i-1]
#          y <- rchisq(1, df = xt)
#          num <- f(y, sigma) * dchisq(xt, df = y)
#          den <- f(xt, sigma) * dchisq(y, df = xt)
#          if (u[i] <= num/den){
#            x[i] <- y
#          } else {
#            x[i] <- xt
#            k <- k+1     #y is rejected
#          }
#      }
#      chain <- x[(b+1):m]
#      return(chain)
#  }

## -----------------------------------------------------------------------------
tm2 <- microbenchmark(
  R = RayleighR(5e3,1,2e3),
  C = Rayleigh(5e3,1,2e3)
)
knitr::kable(summary(tm2)[,c(1,3,5,6)])

