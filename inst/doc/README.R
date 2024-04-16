## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----example------------------------------------------------------------------
library("newIMVC")
library("mvtnorm")
###The new IMVC measure###
n=200
x=rnorm(n)
y=x^2+rt(n,2)
IMVC(y,x,K=10,type="nonlinear")
###IMVC based feature screening###
n=200
p=300
pho1=0.8
mean_x=rep(0,p)
sigma_x=matrix(NA,nrow = p,ncol = p)
for (i in 1:p) {
  for (j in 1:p) {
    sigma_x[i,j]=pho1^(abs(i-j))
  }
}
x=rmvnorm(n, mean = mean_x, sigma = sigma_x,method = "chol")
x1=x[,1]
x2=x[,2]
x3=x[,12]
x4=x[,22]
y=2*x1+0.5*x2+3*x3*ifelse(x3<0,1,0)+2*x4+rnorm(n)
IMVCS(y,x,K=5,d=round(n/log(n)),type="nonlinear")
###IMVC based hypothesis test###
n=100
x=rnorm(n)
y=2*x+rt(n,2)
IMVCT(x,y,K=5,type = "linear")
y=2*cos(x)+rt(n,2)
IMVCT(x,y,K=5,type = "nonlinear",num_per = 100)
###IMVC based FDR control###
n=200
p=100
pho1=0.5
mean_x=rep(0,p)
sigma_x=matrix(NA,nrow = p,ncol = p)
for (i in 1:p) {
  for (j in 1:p) {
    sigma_x[i,j]=pho1^(abs(i-j))
  }
}
x=rmvnorm(n, mean = mean_x, sigma = sigma_x,method = "chol")
x1=x[,1]
x2=x[,2]
x3=x[,3]
x4=x[,4]
x5=x[,5]
y=x1+x2+x3+x4+x5+rnorm(n)
IMVCFDR(y,x,K=5,numboot=100,timeboot=50,true_signal=c(1,2,3,4,5),null_method="hist",alpha=0.2)


