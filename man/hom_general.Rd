\name{hom_general}
\alias{hom_general}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Bayesian fitting of Spatial General Model with homoscedastic normal error term.
%%  ~~function to do ... ~~
}
\description{
Performs the Bayesian fitting of Homoscedastic General Model with normal error term
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hom_general(y, X, W1, W2=NULL, nsim, burn, step, b_pri, B_pri, r_pri, lambda_pri, beta_0,
sigma2_0, rho_0, lambda_0, kernel = NULL,
plot = TRUE, mateq=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ Object of class matrix, with the dependent variable
%%     ~~Describe \code{y} here~~
}
  \item{X}{ Object of class matrix, with covariates of model
%%     ~~Describe \code{X} here~~
}
  \item{W1}{ Object of class matrix, nb or listw related to Spatial Contiguity Matrix for response variable, Anselin(1988)
%%     ~~Describe \code{W} here~~
}
  \item{W2}{ Object of class matrix, nb or listw related to Spatial Contiguity Matrix for disturbance terms, Anselin(1988)
%%     ~~Describe \code{W} here~~
}
  \item{nsim}{
  A number that indicates the amount of iterations
%%     ~~Describe \code{nsim} here~~
}
  \item{burn}{
  A number that indicates the amount of iterations to be burn at the beginning of the chain
%%     ~~Describe \code{burn} here~~
}
  \item{step}{
%%     ~~Describe \code{step} here~~
A number that indicates the length between samples in chain that generate the point estimates for each parameter.
}
  \item{b_pri}{
%%     ~~Describe \code{b_pri} here~~
A vector with the prior mean of beta
}
  \item{B_pri}{
%%     ~~Describe \code{B_pri} here~~
A matrix with the prior variance of beta
}
  \item{r_pri}{
%%     ~~Describe \code{r_pri} here~~
A number with the prior shape parameter of sigma^2
}
  \item{lambda_pri}{
%%     ~~Describe \code{lambda_pri} here~~
A number with the prior rate parameter of sigma^2
}
  \item{beta_0}{
%%     ~~Describe \code{beta_0} here~~
A vector with start values for beta chain
}
  \item{sigma2_0}{
%%     ~~Describe \code{sigma2_0} here~~
A number with start value for sigma^2 chain
}

  \item{rho_0}{
%%     ~~Describe \code{sigma2_0} here~~
A number with start value for rho chain
}

  \item{lambda_0}{
%%     ~~Describe \code{rho_0} here~~
A number with start value for lambda chain
}
  \item{kernel}{
%%     ~~Describe \code{kernel} here~~
Distribution used in transition kernel to get samples of rho and lambda, it can be  "uniform" or "normal"
}
  \item{plot}{
%%     ~~Describe \code{plot} here~~
If it is TRUE present the graph of the chains
}
\item{mateq}{
Logical variable indicating whether W1=w2 or not.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
hom_general is a function made in order to fit Spatial General Model with a normal homoscedatic disturbance term through MCMC methods as Metropolis-Hastings algorithm, under two proposals for trasition kernel to get samples of spatial responde and error lag parameters, rho and lambda, respectively.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
List with the following:
\item{Bestimado}{Estimated coefficients of beta}
\item{Sigma2est}{Estimated coefficient of sigma^2}
\item{Rhoest}{Estimated coefficient of rho}
\item{Lambdaest}{Estimated coefficient 1of lambda}
\item{DesvBeta}{Estimated standard deviations of beta}
\item{DesvGamma}{Estimated standard deviation of gamma}
\item{DesvRho}{Estimated standard deviation of rho}
\item{DesvLambda}{Estimated standard deviation of lambda}
\item{AccRate1}{Acceptance Rate for samples of rho}
\item{AccRate2}{Acceptance Rate for samples of lambda}
\item{BIC}{Value of Bayesian Information Criterion}
\item{DIC}{Value of Deviance Information Criterion}
}
\references{
%% ~put references to the literature/web site here ~
1. Cepeda C. E. (2001). Modelagem da variabilidade em modelos lineares generalizados. Unpublished Ph.D. tesis. Instituto de Matematicas.
Universidade Federal do Rio do Janeiro.

2.Cepeda, E. and Gamerman D. (2005). Bayesian Methodology for modeling parameters in the two-parameter exponential family. Estadistica 57, 93 105.

3.Cepeda C., E. and Gamerman D. (2001). Bayesian Modeling of Variance Heterogeneity in Normal Regression Models. Brazilian Journal of Probability and Statistics. 14, 207-221.

4.Luc Anselin, Spatial Econometrics: Methods and Models, Kluwer Academic, Boston, 1988.

5. D. Gamerman, Markov Chains Monte Carlo: Stochastic Simulation for bayesian Inference, Chapman and Hall, 1997.

6. James Le Sage and Kelley Pace, Introduction to Spatial Econometrics, Chapman & Hall/CRC, Boca Raton, 2009.

}
\author{
Jorge  Sicacha-Parada <jasicachap@unal.edu.co>, Edilberto Cepeda-Cuervo <ecepedac@unal.edu.co>
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

library(spdep)
library(mvtnorm)
library(pscl)
n=49
x0=rep(1,n)
x1=runif(n,0,400)
x2=runif(n,10,23)
X=cbind(x0,x1,x2)
sigma2=rep(45,n)
Sigma=diag(sigma2)
data(oldcol)
W1=COL.nb
matstand=nb2mat(W1)
A=diag(n)-0.75*matstand
B=diag(n)-0.20*matstand
miu=solve(A)\%*\%(18+0.026*x1-0.4*x2)
Sigma2=t(solve(A)\%*\%solve(B))\%*\%Sigma\%*\%solve(A)\%*\%solve(B)
y=rmvnorm(1,miu,Sigma2)
y_1=t(y)
y=y_1
data(oldcol)
hom_general(y,X,W1=COL.nb,nsim=500,burn=25,step=5,b_pri=rep(0,3),B_pri=diag(rep(1000,3)),
r_pri=0.01,lambda_pri=0.01,beta_0=rep(0,3),
sigma2_0=90,rho_0=0.5,lambda_0=0.5,kernel="normal",
plot=FALSE,mateq=TRUE)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Bayesian}
\keyword{Metropolis Hastings}
\keyword{Spatial Econometric Models}
\keyword{General Model}
