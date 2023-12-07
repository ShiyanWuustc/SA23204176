## ----fig.align="center"-------------------------------------------------------
set.seed(8)       #随机种子
x <- rnorm(100)
y <- rexp(100)     
layout(matrix(c(1,2,3,4,5,5,6,7,8),3,3,byrow=T),widths = c(1:0.5), heights = c(0.5:0.5))
par(mai=c(0.1,0.3,0.15,0.45), cex.main=0.9)
barplot(runif(8,1,8),col=2:7,main='条形图')
pie(1:12, col=rainbow(6), labels='', border=NA, main='饼图')
plot(x, y, pch=19, col=c(1,2,4), xlab='', ylab='', main='散点图')
plot(rnorm(25), rnorm(25), cex=(y+2), col=2:4, lwd=2, xlab='', ylab='', main='气泡图')
plot(density(y), col=4, lwd=1, xlab='', ylab='', main='核密度图')
hist(rnorm(500), col=3, xlab='', ylab='', main='直方图')
boxplot(x, col=2, main='箱线图')
qqnorm(y, col=1:7, pch=19, xlab='', ylab='',main='Q-Q图')

## ----results='hold'-----------------------------------------------------------
A <- c(1, 2, 4, 5)  
B <- c(2, 3)  
union(A, B)  # 求集合A与B的并集
intersect(A, B)  # 求集合A与B的交集
setdiff(A, B)  # 求集合A与B的差(A-B)

## ----results='hold'-----------------------------------------------------------
(f <- expression(asin(x)))
D(f, "x")

## ----results='hold'-----------------------------------------------------------
(f <- expression(log(cos(exp(x)))))
D(f, "x")

## ----results='hold'-----------------------------------------------------------
library(Deriv)
f <- function(x) log(x + sqrt(x ^ 2 + 1))
(dy <- Deriv(f, "x"))

## ----results='hold'-----------------------------------------------------------
f <- function(x) sin(x)
result <- integrate(f, lower=0, upper=pi)
result[["value"]]
f <- function(x) log(x)
result <- integrate(f, lower=0, upper=5)
result[["value"]]

## ----results='hold'-----------------------------------------------------------
library(gcookbook)                    
sum <- summary(heightweight)
knitr::kable(sum)#数据可视化

## ----results='hold'-----------------------------------------------------------
model <- lm(heightIn~ageYear,heightweight)     
summary(model) 

## ----collapse=TRUE,fig.align="center"-----------------------------------------
library(ggplot2)
ggplot(heightweight, aes(x=ageYear, y=heightIn,color="#5e616d"))+      
  geom_point()+                                                                     
  stat_smooth(method = lm,color="black")+                      
  annotate("text", label = "R^2==0.4249",parse=T,x=16,y=52)+    
  annotate("text", label = "y_hat=37.4356+1.7483x",x=16,y=50)  

## ----table.align="center"-----------------------------------------------------
scores <- c(63, 72, 74, 79, 82, 82, 87, 89, 90, 90)
# R的自定义函数
cdf_table <- function(x){
x <- sort(x) 
n <- length(x) 
tab <- unname(c(table(x))) #频数
pct = tab / n #频率
d <- data.frame( 
x = sort(unique(x)), #删掉重复值并排序
freq = tab, #x中数值出现的频数(出现的次数)
pct = pct, #x中数值对应的频率
cumfreq = cumsum(tab), #频数的累加和
cumpct = cumsum(pct)) #频率的累加和
d #返回d
}
knitr::kable(cdf_table(scores))

## ----fig.align="center"-------------------------------------------------------
d <- data.frame(scores) 
dfreq <- cdf_table(d$scores) #返回上述函数结果
p <- ggplot(data = dfreq, mapping = aes(x = x, y = cumpct )) 
#对上图做进一步细化设定
p + stat_ecdf(geom = "step") +scale_y_continuous( limits = c(-.01, 1.01), expand = c(0, 0),name = "累计百分比") 

## -----------------------------------------------------------------------------
my.sample <- function(x,size,prob){
  if(is.null(prob)){#判断是否指定概率分布
    p0 <- 1/length(x)
    prob <- numeric(length(x)) 
    prob <- rep(p0,length(x)) #若未指定概率分布，默认设置为等概率抽样  
  }
  cp <- cumsum(prob)
  U <- runif(size)
  r <- x[findInterval(U,cp)+1]#抽样
  ct <- as.vector(table(r))#计算频数
  print(ct/sum(ct)/prob)#检验抽样频率与概率的关系：结果越接近1越好
  return(r)
}

## -----------------------------------------------------------------------------
test1 <- my.sample(x=c(1:10),size = 1000,prob = NULL)#抽取对象：数字；无指定概率
test2 <- my.sample(x=letters,size = 1000,prob = NULL)#抽取对象：字母；无指定概率
test3 <- my.sample(x=c(1:10),size = 1000,prob = c(1:10/55))#抽取对象：字母；有指定概率
test3

## -----------------------------------------------------------------------------
n <- 1000
u <- runif(n)
x <- ifelse(u <= 1/2 , log(2*u) ,-log(2-2*u))#逆函数
hist(x, prob = TRUE, main = expression(f(x)==0.5*e^-abs(x)),ylim = c(0, 0.5))#绘画直方图
y <- seq(min(x) ,max(x), .01)
z <- ifelse(y <= 0 , 1/2*exp(y) ,1/2*exp(-y))
lines(y, z)#绘制目标分布曲线

## -----------------------------------------------------------------------------
BETA <- function(n, a, b) {
j <- k <- 0
y <- numeric(n)
while (k < n) {
  u <- runif(1)
  j <- j + 1
  x <- runif(1) #~g(.)
  if (x^(a - 1) * (1 - x)^(b - 1) > u) 
    {
     k <- k + 1
     y[k] <- x    #接受x
    }
  }
print(j) #可观察c的值
y
  }

## -----------------------------------------------------------------------------
a <- 3
b <- 2
y <- BETA(1000, a, b)
hist(y , prob = TRUE)
z <- seq(0, 1, .01)
fz <- 12 * z^(a-1) * (1 - z)^(b-1) #常数c=gama(5)/gama(3)gama(2)=12
lines(z, fz)

## -----------------------------------------------------------------------------
Epan <- function(n) {
u1 <- runif(n, -1, 1)
u2 <- runif(n, -1, 1)
u3 <- runif(n, -1, 1)
x <- numeric(n)
  for (i in 1:n) {
  if ((abs(u3[i]) >= abs(u2[i])) && (abs(u3[i]) >=abs(u1[i])))
  x[i] <- u2[i]
  else x[i] <- u3[i]
  }
  x
}

## -----------------------------------------------------------------------------
z <- Epan(1000)
hist(z, prob = TRUE)
z <- seq(-1, 1, .01)
lines(z, 3/4* (1 - z^2))

## -----------------------------------------------------------------------------
k <- 100#重复次数
n <- 1e6
l1 <- 1
l2 <- 0.8
l3 <- 0.5
d <- 1
Vpihat1 <- 0
Vpihat2 <- 0
Vpihat3 <- 0
for(i in 1:k){
X <- runif(n,0,d/2)
Y <- runif(n,0,pi/2)
p1 <- mean(l1/2*sin(Y)>X)
p2 <- mean(l2/2*sin(Y)>X)
p3 <- mean(l3/2*sin(Y)>X)
Vpihat1 <- Vpihat1+4*(l1/d)^2*(1-p1)/p1^3#计算渐进方差，rou=1
Vpihat2 <- Vpihat2+4*(l2/d)^2*(1-p2)/p2^3#rou=0.8
Vpihat3 <- Vpihat3+4*(l3/d)^2*(1-p3)/p3^3#rou=0.5
}
Vpihat1/k;Vpihat2/k;Vpihat3/k#比较渐进方差

## -----------------------------------------------------------------------------
m <- 1e4
n <- 1e4
mc <- vector()
anti <- vector()
i=1
while (i <= n ){
  U <- runif(m)
  U1 <- runif(m/2)
  U2 <- 1-U1
  mc[i]<-mean(exp(U))#simple monte carlo 
  anti[i]<-mean((exp(U1)+exp(U2))/2)#antithetic variate
  i=i+1
}
V1 <- var(mc)
V2 <- var(anti)
redu <- (V1-V2)/V1
print(redu)

## ----fig.align="center"-------------------------------------------------------
x <- seq(1, 10,0.01)
g <- x^2/sqrt(2*pi)*exp(-x^2/2)
f1 <- 2*dnorm(x-1, mean = 0, sd = 1)
f2 <- dgamma(x-1, 1.5, 1)
plot(x, g, type='l',col=12, ylim = c(0, 1))
lines(x, f1, col=54)
lines(x, f2, col=61)
legend("topright", inset = 0.02, legend = c("g(x)", "f1","f2"), lty = 1,col=c(12,54,61))

## ----fig.align="center"-------------------------------------------------------
plot(x, g/f1, type='l',col=12, ylim = c(0,3))
lines(x, g/f2,col=54)
legend("topright", inset = 0.02, legend = c("g/f1","g/f2"), lty = 1,col=c(12,54))
V1 <- var(g/f1);V2 <- var(g/f2)
c(V1,V2)#g/f1、g/f2的方差

## -----------------------------------------------------------------------------
m <- 1e4
n <- 1e3
theta1 <- vector()
theta2 <- vector()
i <- 1 

for (i in 1:n){
x <- sqrt(rchisq(m, 1)) + 1
g <- x^2/sqrt(2*pi)*exp(-x^2/2)
f1 <- 2*dnorm(x-1, mean = 0, sd = 1)
theta1[i] <- mean(g/f1)

x <- rgamma(m, 1.5, 1) + 1
g <- x^2/sqrt(2*pi)*exp(-x^2/2)
f2 <- dgamma(x-1, 1.5, 1)
theta2[i] <- mean(g/f2)
}

M1 <- mean(theta1);M2 <-mean(theta2)
c(M1,M2)#基于f1和f2的估计值
V1 <- var(theta1);V2 <- var(theta2)
c(V1,V2)#基于f1和f2的估计量方差

## -----------------------------------------------------------------------------
m <- 1e4 
k1 <- 1 ; k2 <- 5 #使用分层抽样
r2 <- m/k2
i <- 1 ; j<- 1
V2 <- vector()

#Example 5.10：不使用分层抽样法
T2 <- vector()
u1 <- runif(m, 0, 1) 
x <- -log(1-(1-exp(-1))*u1)
g <- exp(-x)/(1+x^2)
f1 <- (k1/(1-exp(-1)))*exp(-x)
est1 <- mean(g/f1)
s1 <- sd(g/f1)

#使用分层抽样法
for (j in 1:k2){
  u2 <- runif(m/k2, (j-1)/k2, j/k2)
  x <- -log(1-(1-exp(-1))*u2)
  g <- exp(-x)/(1+x^2)
  f2 <- (k2/(1-exp(-1)))*exp(-x)
  T2[j] <- mean(g/f2)
  V2[j] <- var(g/f2)
  } 
est2 <- sum(T2)
s2 <- sqrt(mean(V2))

#结果展示
c(est1,est2)#无分层 vs 有分层 估计值结果
c(s1,s2)#无分层 vs 有分层 估计量标准差

## -----------------------------------------------------------------------------
#t-interval
m <- 1e5 #repeat times
n <- 20 
tL <- qt(0.025, df=n-1) ; tU <- qt(0.975, df=n-1)
LCL <- UCL <- vector() 
for(i in 1:m){
  x <- rchisq(n, df=2)
  LCL[i] <- mean(x)+tL*sd(x)/sqrt(n)
  UCL[i] <- mean(x)+tU*sd(x)/sqrt(n)
}
mean(LCL<2 & UCL>2) #t-interval CP

#Example 6.4
n <- 20
alpha <- .05
UCL <- replicate(1000, expr = {
x <- rchisq(n, df = 2)
(n-1) * var(x) / qchisq(alpha, df = n-1)
} )
mean(UCL > 4) #CP of interval for variance

## -----------------------------------------------------------------------------
m <- 1e5 #repeat times
n <- 1e3 #sample size
alpha <- .05 #nominal significance level
x<- matrix(0,3,n)
pval <- matrix(0,3,m)
for(i in 1:m){
  set.seed(i)
  x[1,] <- rchisq(n, df=1)
  x[2,] <- runif(n,0,2)
  x[3,] <- rexp(n,rate=1)
  for (j in 1:3){
    pval[j,i] <- 2*(1-pt(abs( (mean(x[j,])-1) / (sd(x[j,]/ sqrt(n) ) )),n-1) )
  }
}
c(mean(pval[1,]<=alpha),mean(pval[2,]<=alpha),mean(pval[3,]<=alpha))

## ----results='hold'-----------------------------------------------------------
m <- 1000 
M <- 1000 #repeat times
alpha <- 0.1
FWER <- FDR <- TPR <- matrix(nrow=2,ncol=M)

for (i in 1:M){
set.seed(i)
p1 <- runif(m*0.95,0,1) #原假设
p2 <- rbeta(m*0.05,0.1,1) #对立假设
p <- c(p1,p2)#合并

p.adj1 <- p.adjust(p,method='bonferroni')#bonferroni
p.adj2 <- p.adjust(p,method='BH') #B-H   

#FWER：原假设为真情况下，故仅使用前950个原假设
FWER[1,i] <- ifelse( sum(p.adj1[1:950]<alpha)>0 ,1,0)#bonferroni
FWER[2,i] <- ifelse( sum(p.adj2[1:950]<alpha)>0 ,1,0)#B-H

#FDR
FDR[1,i] <- sum(p.adj1[1:950]<alpha)/sum(p.adj1<alpha)
FDR[2,i] <- sum(p.adj2[1:950]<alpha)/sum(p.adj2<alpha)

#TPR
TPR[1,i] <- sum(p.adj1[950:m]<alpha)/50
TPR[2,i] <- sum(p.adj2[950:m]<alpha)/50
}
c(mean(FWER[1,]),mean(FWER[2,]))#FWER
c(mean(FDR[1,]),mean(FDR[2,]))#FDR
c(mean(TPR[1,]),mean(TPR[2,]))#TPR



## ----results='hold'-----------------------------------------------------------
lamda <- 2 #true value 
B <- 1000 #replicates
m <- 1000 #repeat times
n <- c(5,10,20) #sample size

#theoretical
the_bia <- lamda/(n-1)
the_sd <- lamda*n/((n-1)*sqrt(n-2))

#bootstrap
thetastar <- vector()
boo_bia <- boo_sd <- matrix(nrow=3,ncol=m)

for(i in 1:m ){
  for (j in 1:3){
    set.seed(i+j)
    x <- rexp(n[j],lamda)
    for (b in 1:B){
    xstar <- sample(x, replace = TRUE )
    thetastar[b] <- 1/mean(xstar)     
    }
    boo_bia[j,i] <- mean(thetastar)-1/mean(x)
    boo_sd[j,i] <- sd(thetastar)
  }
}
the_bia #theoretical
apply(boo_bia,MARGIN = 1,mean)# bootstrap
the_sd #theoretical
apply(boo_sd,MARGIN = 1,mean) #bootstrap


## -----------------------------------------------------------------------------
library(bootstrap)
B1 <- c(100,200,500,1000)   #replicates numbers
k <- 100                    #repeat times for generating se_hat
n <- nrow(law)              #sample size
S <- numeric(k)  
alpha <- 0.05

for (q in 1:4){
  B <- B1[q]  #令replicates numbers取值不同
  R <- se <- numeric(B)  
  for (b in 1:B) {
    set.seed(b)
    i <- sample(1:n, size = n, replace = TRUE) #抽序号
    LSAT <- law$LSAT[i]    #一列向量
    GPA <- law$GPA[i]
    R[b] <- cor(LSAT, GPA)
    for (j in 1:k){
      set.seed(b+j)
      l <- sample(i,replace=TRUE)
      S[j] <- cor(law$LSAT[l],law$GPA[l]) 
    }
    se[b] <- sd(S)        #Se^(theta^(b))
  }
  theta <- mean(R)
  SD <- sd(R)             #Se^(theta^)
  t <- sort((R-theta)/se) #t^*_a
  print(c(theta-t[length(t)*(1-alpha/2)]*SD,theta-t[length(t)*alpha/2]*SD)) 
  #B=100,200,500,1000下置信区间
}

## -----------------------------------------------------------------------------
library(boot)
x <- c(3,5,7,18,43,85,91,98,100,130,230,487)
lamda <- function(y,i) mean(y[i])
de <- boot(data=x,statistic=lamda, R = 5e3)
de  #lambda
ci <- boot.ci(de,type=c("norm","basic","perc","bca"))
ci.norm<-ci$norm[2:3]
ci.basic<-ci$basic[4:5]
ci.perc<-ci$percent[4:5]
ci.bca<-ci$bca[4:5]
CI <- cbind(ci.norm, ci.basic, ci.perc, ci.bca)   
CI              #bootstrap confidence intervals
CI[2,]-CI[1,]   #length of ci

## -----------------------------------------------------------------------------
data <- as.matrix(bootstrap::scor)
n <- nrow(data)
theta.jack <- numeric(n)
for (i in 1:n) {
  sigma <- cov(data[-i,])
  lambda <- sort(eigen(sigma)$values)
  theta.jack[i] <- lambda[length(lambda)]/sum(lambda) #jackknife process
}
lambda <- sort(eigen(cov(data))$values)               #sort eigen 
theta.hat <- lambda[length(lambda)]/sum(lambda)       #\hat\theta
bias.jack <- (n-1)*(mean(theta.jack)-theta.hat)       #bias
se.jack <- sqrt((n-1)*mean((theta.jack-theta.hat)^2)) #standard error
round(c(original=theta.hat,bias.jack=bias.jack,se.jack=se.jack),5) #display bias & se

## -----------------------------------------------------------------------------
library(DAAG); attach(ironslag)
n <- length(magnetic) 
e1 <- e2 <- e3 <- e4 <- matrix(nrow=n,ncol=n)
for (i in 1:n-1) {
  for (j in (i+1):n){       #i \ne j
    k <- c(i,j)
    y <- magnetic[-k]
    x <- chemical[-k]       #leave-two-out samples
    
    J1 <- lm(y ~ x)         #model 1
    yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k]
    e1[i,j] <- sum((magnetic[k] - yhat1)^2)
    
    J2 <- lm(y ~ x + I(x^2))#model 2
    yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] +
    J2$coef[3] * chemical[k]^2
    e2[i,j] <- sum((magnetic[k] - yhat2)^2)
    
    J3 <- lm(log(y) ~ x)    #model 3
    logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k]
    yhat3 <- exp(logyhat3)
    e3[i,j] <- sum((magnetic[k] - yhat3)^2)
    
    J4 <- lm(log(y) ~ log(x))#model 4
    logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[k])
    yhat4 <- exp(logyhat4)
    e4[i,j] <- sum((magnetic[k] - yhat4)^2) 
  }
}
c(mean(e1,na.rm=T), mean(e2,na.rm=T), mean(e3,na.rm=T), mean(e4,na.rm=T))  #average squared prediction error
J2  #model 2

## -----------------------------------------------------------------------------
attach(chickwts)            #load data
boxplot(formula(chickwts))  #graphical summary
x1 <- as.vector(weight[feed == "soybean"])
x2 <- as.vector(weight[feed == "linseed"])
x3 <- as.vector(weight[feed == "sunflower"])
detach(chickwts)
R <- 999      #repeat times

CVM <-  function(x1,x2){
  n <- length(x1)
  m <- length(x2)
  N <- n+m
  k0 <- 1:N #index
  z <- c(x1, x2)  #pooled sample
  Fn <- Gm <- numeric(N)
  
  #W2 statistic
  for (i in 1:N) {
    Fn[i] <- mean(as.integer(z[i] <= x1))
    Gm[i] <- mean(as.integer(z[i] <= x2))
  }
  W2 <- m*n/N^2*sum((Fn-Gm)^2)  
  
  #W2* statitic from permutation
  W2star <- numeric(R)
  for(i in 1:R){
    k <- sample(k0)
    zstar <- z[k]
    x1star <- zstar[1:n]
    x2star <- zstar[(n+1):N]  
    for (j in 1:N) {
      Fn[j] <- mean(as.integer(zstar[j] <= x1star))
      Gm[j] <- mean(as.integer(zstar[j] <= x2star))
    }
    W2star[i] <- m*n/N^2*sum((Fn-Gm)^2)
  }
  W2star <- c(W2star,W2)
  p <- mean(W2star>=W2)
  return(p) #p.value
}
CVM(x1,x2)  #Example 8.1
CVM(x2,x3)  #Example 8.2

## -----------------------------------------------------------------------------
R <- 999  #repeat times
maxout <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  return(max(c(outx, outy)))
}
Count5 <- function(x, y) {
  z <- c(x, y)
  n <- length(x)
  m <- length(y)
  N <- n+m
  stats <- replicate(R, expr = {
    k <- sample(1:N)
    zstar <- z[k]
    xstar <- zstar[1:n]
    ystar <- zstar[(n+1):N]
    maxout(xstar, ystar)
    })
  stat <- maxout(x, y)
  statstar <- c(stats, stat)
  return(list(estimate = stat, p = mean(statstar>=stat)))
 }

## -----------------------------------------------------------------------------
#unequal sample sizes 
n1 <- 20
n2 <- 40
mu1 <- mu2 <- 0

#1.equal variance
sigma1 <- sigma2 <- 1
x <- rnorm(n1, mu1, sigma1)
y <- rnorm(n2, mu2, sigma2)
Count5(x, y)

#2. unequal variance 
sigma1 <- 1
sigma2 <- 2
x <- rnorm(n1, mu1, sigma1)
y <- rnorm(n2, mu2, sigma2)
Count5(x, y)

## -----------------------------------------------------------------------------
Model <- function(N,b1,b2,b3,f0){
  g <- function(alpha){  
    x1 <- rpois(N, 1)
    x2 <- rexp(N, 1)
    x3 <- rbinom(N, 1, 0.5)
    tmp <- exp(-alpha-b1*x1-b2*x2-b3*x3)
    p <- 1/(1+tmp)
    mean(p) - f0
    }
  solution <- uniroot(g,c(-20,20))
  return(solution)
}

## -----------------------------------------------------------------------------
N = 10e6
b1 = 0
b2 = 1
b3 = -1
f0 = c(0.1,0.01,0.001,0.0001)
x <- matrix(nrow =4,ncol = 5)
k <- 1
f0 <- c(0.1,0.01,0.001,0.0001)
for(i in f0){
  x[k,] <- unlist(Model(N,b1,b2,b3,i))
  k <- k+1
}
x <- x[, -c(4:5), drop = FALSE]
colnames(x) <- c("root", "f.root", "iter")
rownames(x) <- c("f0=0.1", "f0=0.01", "f0=0.001","f0=0.0001")
x

## -----------------------------------------------------------------------------
plot(x[,1],-log(f0), type = "l")

## -----------------------------------------------------------------------------
rw.Laplace <- function(N, x0, sigma) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)

  update_value <- function(i, xt, sigma, u) {
    y <- rnorm(1, xt, sigma)
    ifelse(u[i] <= exp(abs(xt) - abs(y)), y, xt)
  }

  for (i in 2:N) {
    x[i] <- update_value(i, x[i - 1], sigma, u)
  }

  k <- sum(x[-1] == x[-N])
  return(list(x = x, k = k))
}

run_simulation <- function(N, x0, sigmas) {
  simulations <- lapply(sigmas, function(sigma) rw.Laplace(N, x0, sigma))
  list(
    k_values = sapply(simulations, function(sim) sim$k),
    x_values = lapply(simulations, function(sim) sim$x)
  )
}

N <- 5e3
sigma <- c(0.5, 1, 4, 16)
x0 <- rnorm(1)
results <- run_simulation(N, x0, sigma)
print(results$k_values/N)  

## -----------------------------------------------------------------------------
y1 <- results$x_values[[1]]
y2 <- results$x_values[[2]]
y3 <- results$x_values[[3]]
y4 <- results$x_values[[4]]

b <- 2e3
plot(y1, type = "l")
plot(y2, type = "l")
plot(y3, type = "l")
plot(y4, type = "l")

#throw buin-in sample
p <- ppoints(100)#QQ plots
y <- qexp(p, 1)
z <- c(-rev(y), y)

Q1 <- quantile(y1[(b + 1):N], p)
qqplot(z, Q1, cex = 0.2)
abline(0, 1)
Q2 <- quantile(y2[(b + 1):N], p)
qqplot(z, Q2, cex = 0.2)
abline(0, 1)
Q3 <- quantile(y3[(b + 1):N], p)
qqplot(z, Q3, cex = 0.2)
abline(0, 1)
Q4 <- quantile(y4[(b + 1):N], p)
qqplot(z, Q4, cex = 0.2)
abline(0, 1)

fx <- exp(-abs(z))/2
hist(y1[(b + 1):N], breaks = "Scott", freq = FALSE, ylim = c(0,0.5))
lines(z, fx)
hist(y2[(b + 1):N], breaks = "Scott", freq = FALSE, ylim = c(0,0.5))
lines(z, fx)
hist(y3[(b + 1):N], breaks = "Scott", freq = FALSE, ylim = c(0,0.5))
lines(z, fx)
hist(y4[(b + 1):N], breaks = "Scott", freq = FALSE, ylim = c(0,0.5))
lines(z, fx)

## -----------------------------------------------------------------------------
N <- 5e3
b <- 1e3
rho <- 0.9
mu1 <- mu2 <- 0
sigma1 <- sigma2 <- 1
s1 <- sqrt(1 - rho^2) * sigma1
s2 <- sqrt(1 - rho^2) * sigma2

generate_observation <- function(x_prev, mu1, mu2, rho, sigma1, sigma2, s1, s2) {
  m1 <- mu1 + rho * (x_prev[2] - mu2) * sigma1 / sigma2 #conditional mean
  x1 <- rnorm(1, m1, s1)
  m2 <- mu2 + rho * (x1 - mu1) * sigma2 / sigma1
  x2 <- rnorm(1, m2, s2)
  return(c(x1, x2))
}

X <- matrix(0, N, 2)
X[1, ] <- c(mu1, mu2)

for (i in 2:N) {
  X[i, ] <- generate_observation(X[i - 1, ], mu1, mu2, rho, sigma1, sigma2, s1, s2)
}

x <- X[(b+1):N, ]
X <- x[, 1]
Y <- x[, 2]
plot(X, Y,cex = 0.5)

## -----------------------------------------------------------------------------
M <- lm(Y ~ X)
summary(M)#model

qqnorm(M$res)#Normality
qqline(M$res)

plot(M$fit, M$res, cex = 0.5) #constant

## -----------------------------------------------------------------------------
library(coda)
Gelman.Rubin <- function(psi) {
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  psi.means <- rowMeans(psi)
  B <- n * var(psi.means)
  psi.w <- apply(psi, 1, "var")
  W <- mean(psi.w)
  v.hat <- W * (n - 1)/n + (B/n)
  r.hat <- v.hat/W
  return(r.hat)
 }

f <- function(x, sigma) {
  if (sigma <= 0) {
    stop("error set")
  }
  if (x < 0) {
    return(0)
  }
  result <- (x / sigma^2) * exp(-x^2 / (2 * sigma^2))
  result
}

Rayleigh.MH.chain1 <- function(sigma, m, x0) {
  x <- numeric(m)
  x[1] <- x0
  u <- runif(m)

  #generate chain
  for (i in 2:m) {
    xt <- x[i - 1]
    y <- rchisq(1, df = xt)
    num <- f(y, sigma) * dchisq(xt, df = y)
    den <- f(xt, sigma) * dchisq(y, df = xt)
    accept_ratio <- num / den
    if (u[i] <= accept_ratio) {
      x[i] <- y
    } else {
      x[i] <- xt
    }
  }
  x
}

#set
sigma = 4
x0 = c(1/sigma^2, 1/sigma, sigma^2, sigma^3)
k = 4
m = 2e3
X = matrix(0, nrow = k, ncol = m)
for (i in 1:k) X[i, ] = Rayleigh.MH.chain1(sigma, m,x0[i])
psi = t(apply(X, 1, cumsum))
for (i in 1:nrow(psi)) psi[i, ] = psi[i, ]/(1:ncol(psi))
rhat=Gelman.Rubin(psi)
rhat

# use the coda package to check for convergence of the chain
library(coda)
mcmc_chains <- lapply(1:4, function(i) as.mcmc(X[i, ]))
Y <- mcmc.list(mcmc_chains)
print(gelman.diag(Y))
gelman.plot(Y, col = c(1, 1))


## -----------------------------------------------------------------------------
u=c(11,8,27,13,16,0,23,10,24,2)
v=c(12,9,28,14,17,1,24,11,25,3)

eps=10e-1
lam=0.7
flag=1
while (abs(flag)>eps)
{
  lam1=lam-(sum(v-u*exp(lam)))/sum(-u*exp(lam))
  flag=lam1-lam
  lam=lam1
}
lam #MLE

## -----------------------------------------------------------------------------
eps=10e-1
lamb=0.7
flag=1
n=10
while (abs(flag)>eps)
{
  lam1=n*lam/(n-lam*sum((v-u*exp(lam))/(exp(lam)-1)))
  flag=lam1-lam
  lam=lam1
}
lam #EM

## -----------------------------------------------------------------------------
solve.game <- function(A) {
min.A <- min(A)
A <- A - min.A
max.A <- max(A)
A <- A/max(A)
m <- nrow(A)
n <- ncol(A)
it <- n^3
a <- c(rep(0, m), 1)
A1 <- -cbind(t(A), rep(-1, n))
b1 <- rep(0, n)
A3 <- t(as.matrix(c(rep(1, m), 0)))
b3 <- 1
sx <- simplex(a = a, A1 = A1, b1 = b1, A3 = A3, b3 = b3,
maxi = TRUE, n.iter = it)
a <- c(rep(0, n), 1)
A1 <- cbind(A, rep(-1, m))
b1 <- rep(0, m)
A3 <- t(as.matrix(c(rep(1, n), 0)))
b3 <- 1
sy <- simplex(a = a, A1 = A1, b1 = b1, A3 = A3, b3 = b3,
maxi = FALSE, n.iter = it)
soln <- list(A = A * max.A + min.A, x = sx$soln[1:m],
y = sy$soln[1:n], v = sx$soln[m + 1] * max.A +min.A)
soln
}

## -----------------------------------------------------------------------------
A <- matrix(c(0, -2, -2, 3, 0, 0, 4, 0, 0, 2, 0, 0, 0,
-3, -3, 4, 0, 0, 2, 0, 0, 3, 0, 0, 0, -4, -4, -3,
0, -3, 0, 4, 0, 0, 5, 0, 0, 3, 0, -4, 0, -4, 0, 5,
0, 0, 3, 0, 0, 4, 0, -5, 0, -5, -4, -4, 0, 0, 0,
5, 0, 0, 6, 0, 0, 4, -5, -5, 0, 0, 0, 6, 0, 0, 4,
0, 0, 5, -6, -6, 0), 9, 9)
library(boot)
B <- A + 2
s <- solve.game(B)
s$v
round(cbind(s$x, s$y), 7)
round(s$x * 61, 7)

## -----------------------------------------------------------------------------
a <- list('a',1,1.5)
str(unlist(a))
str(as.vector(a))
b <- c('1','2',3)
str(as.vector(b))

## -----------------------------------------------------------------------------
a <- c(1,2)
dim(a)

## -----------------------------------------------------------------------------
x <- matrix()
is.matrix(x)
is.array(x)

## -----------------------------------------------------------------------------
df <- data.frame(
  A = c(1, 2, 3),      
  B = c('a', 'b', 'c'), 
  C = c(1.5, 2.5, 3.5)
)
m <- as.matrix(df)
print(m)
str(m)

## -----------------------------------------------------------------------------
#data frame with 0 rows
df1<- data.frame(A = numeric(), B = character(), C = integer())

#data frame with 0 columns
df2<- data.frame()

print(str(df1))
print(str(df2))


## -----------------------------------------------------------------------------
scale01 = function(x) {
         rng = range(x, na.rm = TRUE)
         (x - rng[1]) / (rng[2] - rng[1])
}

## -----------------------------------------------------------------------------
df1 <- data.frame(matrix(data = 1:12,nrow=4,ncol=3))
df2 <- data.frame(
  A = c(1, 2, 3),      
  B = c('a', 'b', 'c'), 
  C = c(1.5, 2.5, 3.5),
  D = c(TRUE,FALSE,TRUE)
)

#1.apply it to every column of a data frame
data.frame(lapply(df1, function(x) scale01(x)))  

#2.apply it to every numeric column in a data frame (i.e. exclude the non-numeric columns)
data.frame(lapply(df2, function(x) if (is.numeric(x)) scale01(x) else x))

## -----------------------------------------------------------------------------
df1 <- data.frame(matrix(data = 1:12,nrow=4,ncol=3))
df2 <- data.frame(
  A = c(1, 2, 3),      
  B = c('a', 'b', 'c'), 
  C = c(1.5, 2.5, 3.5),
  D = c(TRUE,FALSE,TRUE)
)

#a)
data.frame(vapply(df1, function(x) sd(x),numeric(1))) #sd

#b)
df2b <- df2[vapply(df2, function(x) is.numeric(x), logical(1))]  #choose numeric column
data.frame(vapply(df2b, function(x) sd(x), numeric(1))) #sd

## -----------------------------------------------------------------------------
R <- function( a,b,n,N,burn ){
  x <- y <- rep(0, N)
  x[1] <- 7 #initial value
  y[1] <- rbeta(1, x[1] + a, n - x[1] + b) 
  for (i in 2:N){  #generate chain
    x[i] <- rbinom(1, prob = y[i - 1], size = n) #conditional distribution
    y[i] <- rbeta(1, x[i] + a, n - x[i] + b)
  }
  x <- x[(burn+1):N]
  return(x)
}

## -----------------------------------------------------------------------------
library(Rcpp)
cppFunction('std::vector<int> rcpp(int a, int b, int n, int N, int burn) {
  std::vector<int> x(N);
  std::vector<double> y(N);
  std::vector<int> result(N - burn);
  x[0] = 7;
  y[0] = R::rbeta(x[0] + a, n - x[0] + b);

  for (int i = 1; i < N; i++) {
      x[i] = R::rbinom(n, y[i - 1]);
      y[i] = R::rbeta(x[i] + a, n - x[i] + b);
    }
  for (int i = burn; i < N; i++) {
      result[i - burn] = x[i];
  }
  return result;
}'
)

## -----------------------------------------------------------------------------
a <- 2
b <- 3
n <- 10 #fixed para
N <- 1e4 
burn <- 5e3 #burn-in

library(microbenchmark)
rbind(microbenchmark(R(a,b,n,N,burn)), microbenchmark(rcpp(a,b,n,N,burn))) #computation time

