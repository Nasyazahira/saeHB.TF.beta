## code to prepare `DATASET` dataset goes here

set.seed(123)
narea=10
nsub=3
n=narea*nsub
codearea=c()
for(i in 1:narea){
  codearea=c(codearea,rep(i,nsub))
}
X1=runif(n,0,1)
X2=runif(n,0,1)
b0=b1=b2=0.5
v=rnorm(narea,0,1) #Area Random Effect
u=rnorm(n,0,1) #Subarea Random Effect
vim=c()
for(i in 1:n){vim[i]=v[codearea[i]]}
Mu <- exp(b0+b1*X1+b2*X2+u+vim)/(1+exp(b0+b1*X1+b2*X2+u+vim))
pi=rgamma(n,1,0.5)
A=Mu*pi
B=(1-Mu) * pi
y=rbeta(n,A,B)
y <- ifelse(y<1,y,0.99)
y <- ifelse(y>0,y,0.01)
MU=A/(A+B)
vardir = (A*B)/((A+B)^2*(A+B+1))
w <- runif(n,0.2,0.7)
dataBeta <- data.frame(y, X1, X2, vardir,codearea,w)

usethis::use_data(dataBeta, overwrite = TRUE)
