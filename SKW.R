#-------------------------------------------------------------------------------
#Seno KUMARASWAMY(mu,sigma)
#SKW(mu,sigma)
#-------------------------------------------------------------------------------
#LOG LIKELIHOOD ----
LSKW<-expression(log(pi)+log(mu)+log(sigma)+(mu - 1)*log(y)+sigma*log(1-y^mu)+
                      log((2*(1-y^mu)^sigma-3))+log(cos((pi*(1-(1-y^mu)^sigma)*(2-(1-y^mu)^sigma))/4))-
                      log(4)-log((y^mu-1)))
#-------------------------------------------------------------------------------
#SCORE VECTOR ----
dfdm<-D(LSKW,"mu")
d2fd2m<-D(dfdm,"mu")
dfds<-D(LSKW,"sigma")
d2fd2s<-D(dfds,"sigma")
d2fdmds<-D(dfdm,"sigma")

#-------------------------------------------------------------------------------
#GAMLSS FUNCTION ----
SKW<-function(mu.link = "log", sigma.link = "log") 
{
  mstats <- checklink("mu.link", "SKW", substitute(mu.link), 
                      c("logit","probit","cloglog","identity","log","own"))
  dstats <- checklink("sigma.link", "SKW", substitute(sigma.link), 
                      c("inverse","log","identity","own"))
  structure(list(family = c("SKW"), 
                 parameters = list(mu = TRUE, sigma = TRUE), 
                 nopar = 2, 
                 type = "Continuous", 
                 mu.link = as.character(substitute(mu.link)), 
                 sigma.link = as.character(substitute(sigma.link)), 
                 mu.linkfun = mstats$linkfun, 
                 sigma.linkfun = dstats$linkfun, 
                 mu.linkinv = mstats$linkinv, 
                 sigma.linkinv = dstats$linkinv, 
                 mu.dr = mstats$mu.eta, 
                 sigma.dr = dstats$mu.eta, 
                 dldm = function(y, mu, sigma) {
                   dldm <- eval(dfdm)
                   dldm
                 }, 
                 d2ldm2 = function(y,mu, sigma) {
                   d2ldm2 <- eval(d2fd2m)
                   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15) 
                   d2ldm2
                 }, 
                 dldd = function(y, mu, sigma) {
                   dldd <- eval(dfds)
                   dldd
                 },
                 d2ldd2 = function(y,mu, sigma) {
                   d2ldd2 = eval(d2fd2s)
                   d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)  
                   d2ldd2
                 },
                 d2ldmdd = function(y,mu, sigma) {
                   d2ldmdd = eval(d2fdmds)
                   d2ldmdd<-ifelse(is.na(d2ldmdd)==TRUE,0,d2ldmdd)
                   d2ldmdd  
                 }, 
                 #--------------------------------------------------------------
                 G.dev.incr = function(y, mu, sigma, w, ...) -2 * log(dSKW(y=y, mu=mu, sigma=sigma)), 
                 rqres = expression(
                   rqres(pfun = "pSKW",  type = "Continuous", y = y, mu = mu, sigma = sigma)
                 ),
                 mu.initial = expression(mu <- rep(mean(y),length(y))),   
                 sigma.initial = expression(sigma<- rep(0.5, length(y))),
                 mu.valid = function(mu) all(mu > 0), 
                 sigma.valid = function(sigma)  all(sigma > 0),
                 y.valid = function(y) all(y > 0 & y < 1)
  ), 
  class = c("gamlss.family", "family"))
}

#-------------------------------------------------
# MATH FUNCTIONS - Unit Sine Kumaraswamy
#-------------------------------------------------
#Probability Distribution Function----
dSKW<-function(y,mu,sigma,log=FALSE)
{
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", "")) 
  if (any(y <= 0) | any(y >= 1)) stop(paste("y must be between 0 and 1", "\n", ""))
  fy1<-(pi*mu*sigma*y^(mu-1)*(1-y^mu)^sigma*(2*(1-y^mu)^sigma-3)*cos((pi*(1-(1-y^mu)^sigma)*(2-(1-y^mu)^sigma))/4))/(4*(y^mu-1))
  if(log==FALSE) fy<-fy1 else fy<-log(fy1)
  return(fy)
}
#integrate(dSKW,lower=0,upper=1,mu=0.5,sigma=1.2)
#-------------------------------------------------------------------------------
#Cumulative Distribution Function----
pSKW<-function(q, mu, sigma, lower.tail = TRUE, log.p = FALSE){
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", ""))
  if (any(q <= 0) | any(q >= 1)) stop(paste("q must be between 0 and 1", "\n", ""))
  cdf1<-sin((pi*(1-(1-q^mu)^sigma)*(2-(1-q^mu)^sigma))/4)
  if(lower.tail==TRUE) cdf<-cdf1 else cdf<- 1-cdf1
  if(log.p==FALSE) cdf<- cdf else cdf<- log(cdf)
  return(cdf)
}
#-------------------------------------------------------------------------------
#Quantile Function ----
qSKW<-function(p,mu,sigma)
{
  y<-(1-(3-sqrt(9-4*(2-(4*asin(p))/pi)))^(1/sigma)/2^(1/sigma))^(1/mu)
  return(y)
} 
#-------------------------------------------------------------------------------
#Random Generation Function ----
rSKW<-function(n,mu,sigma)
{
  u<-runif(n)
  y<-(1-(3-sqrt(9-4*(2-(4*asin(u))/pi)))^(1/sigma)/2^(1/sigma))^(1/mu)
  return(y)
}
#-------------------------------------------------------------------------------

               