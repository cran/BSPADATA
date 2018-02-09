
hetero_sar=function(y,X,Z,W,nsim,burn,step,b_pri,B_pri,g_pri,G_pri,beta_0,gammas_0,rho_0,kernel=NULL,plot=TRUE)
{
########Lectura de la informaci?n
rowst=function(x){
  x1=c()
  x1=(x)/sum(x)
}
y=as.matrix(y)
  if (is.null(X) | is.null(y) ){
    stop("No data")
  }
  if(burn>nsim | burn<0){
  stop("Burn must be between 0 and nsim")
  }
  if(nsim<=0){
    stop("There must be more than 0 simulations")
  }
  if(step<0 | step > nsim){
    stop("Jump length must not be lesser than 0 or greater than nsim")
  }
  if(class(W)=="nb"){
   matstand=nb2mat(W)
   mat0=nb2listw(W,style="B")
   mat=listw2mat(mat0)
  }
else{
if(class(W)=="listw"){
mat=listw2mat(W)
matstand=apply(mat,2,rowst)
matstand=t(matstand)
}
else{
  if(sum(rowSums(W))==nrow(X))
  {
  matstand=W
  mat=matrix(nrow=nrow(X),ncol=nrow(X))
  for(i in 1:nrow(mat)){
  for(j in 1:ncol(mat)){
if(matstand[i,j]==0){mat[i,j]=0}
else{mat[i,j]=1/matstand[i,j]}
  }
  }
  }
  else{
  mat=W
  matstand=apply(mat,2,rowst)
  matstand=t(matstand)
  }
  }
  }
dpost <- function(betas,gammas,rho) {
A=diag(nrow(X))-rho*matstand
Sigma=diag(c(exp(Z%*%gammas)))
solvSigma=diag(1/c(exp(Z%*%gammas)))
k=t(A%*%y-X%*%(betas))%*%solvSigma%*%(A%*%y-X%*%(betas))
fc.y=k
fc.beta=t(betas - b_pri)%*%solve(B_pri)%*%(betas-b_pri)
fc.gamma=t(gammas - g_pri)%*%solve(G_pri)%*%(gammas-g_pri)
dp <- (det(Sigma)^(-1/2))*det(A)*exp(-0.5*fc.y)*exp(-0.5*fc.beta)*exp(-0.5*fc.gamma)
dp
}

#Generaci?n de valores para la distribuci?n propuesta
r.proposal=function(Gammas){
a.now=Z%*%Gammas
A=diag(nrow(X))-Rho*matstand
b.now=A%*%y-X%*%betas.now
y.now=a.now+(b.now^2/exp(a.now))-1
G_pos=solve(solve(G_pri)+0.5*t(Z)%*%Z)
g_pos=G_pos%*%(solve(G_pri)%*%g_pri+0.5*(t(Z)%*%y.now))
gammas.pro=rmvnorm(1,g_pos,G_pos)
gammas.pro
}

dproposal<-function(gammas.now, gammas.old){
a.now=Z%*%gammas.old
A=diag(nrow(X))-Rho*matstand
b.now=A%*%y-X%*%betas.now
y.now=a.now+(b.now^2/exp(a.now))-1
G_pos=solve(solve(G_pri)+0.5*t(Z)%*%Z)
g_pos=G_pos%*%(solve(G_pri)%*%g_pri+0.5*(t(Z)%*%y.now))
dmvnorm(gammas.now,g_pos,G_pos)
}

dproposal_rho<-function(rho){
Sigma=diag(c(exp(Z%*%Gammas)))
a=t(y)%*%t(matstand)%*%solve(Sigma)%*%matstand%*%y
b=t(y)%*%t(matstand)%*%solve(Sigma)%*%(y-X%*%betas.now)
dnorm(rho,b/a,1/sqrt(a))
}


#Algoritmo Metropolis Hastings
beta.mcmc=matrix(NA,nrow=nsim,ncol(X))
gamma.mcmc=matrix(NA,nrow=nsim,ncol(Z))
rho.mcmc=c()
ind1=rep(0,nsim)
ind2=rep(0,nsim)
logV_DIC=c()
pb <- txtProgressBar(min = 0, max = nsim, style = 3)


if(kernel=="uniform"){
for(i in 1:nsim){
#Valores a posteriori condicional
if(i==1){
Gammas=gammas_0
Sigma=diag(c(exp(Z%*%Gammas)))
Rho=rho_0
}
else{
Sigma=diag(c(exp(Z%*%Gammas)))
}
B_pos=solve(solve(B_pri)+t(X)%*%solve(Sigma)%*%X)
b_pos=B_pos%*%(solve(B_pri)%*%b_pri+t(X)%*%solve(Sigma)%*%y-Rho*t(X)%*%solve(Sigma)%*%matstand%*%y)
#Beta a posteriori condicional
betas.now=c(rmvnorm(1,b_pos,B_pos))

#A posteriori condicional completa para Sigma2
gammas.now=c(r.proposal(Gammas))
q1.1=dproposal(gammas.now,Gammas)
q2.1=dproposal(Gammas,gammas.now)
p1.1=dpost(betas.now,gammas.now,Rho)
p2.1=dpost(betas.now,Gammas,Rho)
T.val=min(1,(p1.1/p2.1)*(q1.1/q2.1))
u<-runif(1)
if(p2.1==0){T.val=0}
if(q2.1==0){T.val=0}
if (u <=T.val) {
Gammas= gammas.now
ind1[i] = 1
}
#A posteriori condicional completa para Rho
rho.now=runif(1,1/abs(min(eigen(mat)$values)),1)
p1.2=dpost(betas.now,Gammas,rho.now)
p2.2=dpost(betas.now,Gammas,Rho)
T.val2=min(1,p1.2/p2.2)
u<-runif(1)
if(p2.2==0){T.val2=0}
if (u <=T.val2) {
Rho <- rho.now
ind2[i] = 1
}
beta.mcmc[i,]<-betas.now
gamma.mcmc[i,]<-gammas.now
rho.mcmc[i]<-rho.now
Sigma=diag(c(exp(Z%*%gamma.mcmc[i,])))
detS=det(Sigma)
detB=det(diag(nrow(X))-rho.mcmc[i]*matstand)
Yg=(diag(nrow(X))-rho.mcmc[i]*matstand)%*%y-X%*%beta.mcmc[i,]
logV_DIC[i]=(-(nrow(X)/2)*log(pi))+log(detB)-0.5*log(detS)-0.5*t(Yg)%*%solve(Sigma)%*%Yg
Sys.sleep(0.000000001)
# update progress bar
setTxtProgressBar(pb, i)
      }
}

if(kernel=="normal")
{
for(i in 1:nsim){
#Valores a posteriori condicional
if(i==1){
Gammas=gammas_0
Sigma=diag(c(exp(Z%*%Gammas)))
Rho=rho_0
}
else{
Sigma=diag(c(exp(Z%*%Gammas)))
}
B_pos=solve(solve(B_pri)+t(X)%*%solve(Sigma)%*%X)
b_pos=B_pos%*%(solve(B_pri)%*%b_pri+t(X)%*%solve(Sigma)%*%y-Rho*t(X)%*%solve(Sigma)%*%matstand%*%y)
#Beta a posteriori condicional
betas.now=c(rmvnorm(1,b_pos,B_pos))
#A posteriori condicional completa para Sigma2
gammas.now=c(r.proposal(Gammas))
q1.1=dproposal(gammas.now,Gammas)
q2.1=dproposal(Gammas,gammas.now)
p1.1=dpost(betas.now,gammas.now,Rho)
p2.1=dpost(betas.now,Gammas,Rho)
T.val=min(1,(p1.1/p2.1)*(q1.1/q2.1))
u<-runif(1)
if(p2.1==0){T.val=0}
if(q2.1==0){T.val=0}
if (u <=T.val) {
Gammas= gammas.now
ind1[i] = 1
}
###Propuesta de Rho
Sigma=diag(c(exp(Z%*%Gammas)))
solvSigma=diag(1/c(exp(Z%*%Gammas)))
a=t(y)%*%t(matstand)%*%solvSigma%*%matstand%*%y
b=t(y)%*%t(matstand)%*%solvSigma%*%(y-X%*%betas.now)
rho.now=rnorm(1,b/a,1/sqrt(a))
p1.2=dpost(betas.now,Gammas,rho.now)
p2.2=dpost(betas.now,Gammas,Rho)
q1.2=dproposal_rho(rho.now)
q2.2=dproposal_rho(Rho)
T.val2=min(1,(p1.2/p2.2)*(q1.2/q2.2))
u<-runif(1)
if(p2.2==0){T.val2=0}
if (u <=T.val2) {
Rho <- rho.now
ind2[i] = 1
}
beta.mcmc[i,]<-betas.now
gamma.mcmc[i,]<-gammas.now
rho.mcmc[i]<-rho.now
Sigma=diag(c(exp(Z%*%gamma.mcmc[i,])))
detS=det(Sigma)
detB=det(diag(nrow(X))-rho.mcmc[i]*matstand)
Yg=(diag(nrow(X))-rho.mcmc[i]*matstand)%*%y-X%*%beta.mcmc[i,]
logV_DIC[i]=(-(nrow(X)/2)*log(pi))+log(detB)-0.5*log(detS)-0.5*t(Yg)%*%solve(Sigma)%*%Yg
Sys.sleep(0.000000001)
# update progress bar
setTxtProgressBar(pb, i)
}
}
beta.mcmc_1=beta.mcmc[(burn+1):nsim,]
  gamma.mcmc_1=gamma.mcmc[(burn+1):nsim,]
  rho.mcmc_1=rho.mcmc[(burn+1):nsim]
  beta.mcmc_2=matrix(NA,nrow=(nsim-burn+1)/step,ncol(X))
  gamma.mcmc_2=matrix(NA,nrow=(nsim-burn+1)/step,ncol(Z))
  rho.mcmc_2=c()
  for (i in 1:(nsim-burn+1))
  {
    if(i%%step==0)
    {
      beta.mcmc_2[i/step,]=beta.mcmc_1[i,]
      gamma.mcmc_2[i/step,]=gamma.mcmc_1[i,]
      rho.mcmc_2[i/step]=rho.mcmc_1[i]
    }
  }

if(plot=="TRUE"){
 for (i in 1:ncol(X)) {
            dev.new()
            ts.plot(beta.mcmc[, i], main = paste("Complete chain for beta",
                i), xlab = "number of iterations", ylab = paste("parameter beta",
                i))}
for (i in 1:ncol(Z)) {
            dev.new()
            ts.plot(gamma.mcmc[, i], main = paste("Complete chain for gamma",
                i), xlab = "number of iterations", ylab = paste("parameter gamma",
                i))}
dev.new()
ts.plot(rho.mcmc,main ="Complete chain for Rho", xlab = "number of iterations", ylab ="parameter rho")
}
 rho.mcmc_3=rho.mcmc_2[rho.mcmc_2<=1]
  Bestimado = colMeans(beta.mcmc_2)
  Gammaest = colMeans(gamma.mcmc_2)
  Rhoest=mean(rho.mcmc_3)
DesvBeta <- apply(beta.mcmc_2,2,sd)
DesvGamma <- apply(gamma.mcmc_2,2,sd)
DesvRho<-sd(rho.mcmc_3)
AccRate1<-sum(ind1)/nsim
AccRate2<-sum(ind2)/nsim

Sigma1=diag(c(exp(Z%*%Gammaest)))
solvSigma1=diag(1/c(exp(Z%*%Gammaest)))
detS=det(Sigma1)
detA=det(diag(nrow(X))-Rho*matstand)
Yg=(diag(nrow(X))-Rho*matstand)%*%y-X%*%Bestimado
Veros=detA*((detS)^(-0.5))*exp(-0.5*t(Yg)%*%solvSigma1%*%Yg)
p=ncol(X)+ncol(Z)+1
BIC=-2*log(Veros)+p*log(nrow(X))
logV_DIC=logV_DIC[is.nan(logV_DIC)==FALSE]
Dbar=mean(-2*logV_DIC)
logV1_DIC=(-(nrow(X)/2)*log(pi))+log(detB)-0.5*log(detS)-0.5*t(Yg)%*%solve(Sigma)%*%Yg
Dev=-2*logV1_DIC
DIC=2*Dbar+Dev
yestimado=solve(diag(nrow(X))-Rhoest*matstand)%*%X%*%Bestimado
residuals=(y)-yestimado
  list(Bestimado = Bestimado, Gammaest = Gammaest, Rhoest= Rhoest,DesvBeta=DesvBeta,DesvGamma=DesvGamma,DesvRho=DesvRho,AccRate1=AccRate1,AccRate2=AccRate2,BIC=BIC,DIC=DIC)
}
