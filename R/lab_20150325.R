#### Functions for Fay lab meeting 3/25/15
## Model fitting using Maximum Likelihood


#' Schaefer biomass dynamics function
#' @author Gavin Fay
#' @description This function calculates biomass based on the Schaefer production function
#' @param bio
#' @param r
#' @param k
#' @param catch
schaefer.bio <- function(bio=100,r=0.2,k=200,catch=0)
{
  #calculate the biomass
  new.bio <- bio + bio*r*(1-bio/k) - catch
  #prevent negative biomass (we might want to add a penalty for this later)
  new.bio <- max(c(1,new.bio),na.rm=TRUE)
  #function returns the value for the new biomass
  return(new.bio)
}

#' Predict the biomass time series
#' @author Gavin Fay
#' @description calculates predicted biomass time series based on the Schaefer
#' production model given values for starting biomass, r, and k. 
predict.bio <- function(start.bio=100,catch=rep(0,10),r=0.2,k=200)
{
  #get length of time series for predictions
  nyears <- length(catch)+1
  #set up our vector for the predicted biomass
  pred.bio <- rep(NA,nyears)
  #initialize first year biomass
  pred.bio[1] <- start.bio
  #predict the biomass for remaining years
  for (year in 2:nyears)
    pred.bio[year] <- schaefer.bio(pred.bio[year-1],r,k,catch[year-1])
  #function returns the biomass time series
  return(pred.bio)
}
  
nll.function.conc <- function(predict.bio,obs) 
{
  nobs <- 0
  q.resid <- rep(NA,length(obs))
  #first evaluate MLE for q
  for (iobs in 1:length(obs))
  {
    #check if there was an observation this year. NB Change this!
    if (obs[iobs]!=-99)
    {
      nobs <- nobs + 1
      q.resid[iobs] <- log(obs[iobs]/
                             mean(predict.bio[iobs:(iobs+1)],na.rm=TRUE))
    }
  }
  # equation for q, note that na.rm=TRUE in sum() means the 
  # function will ignore the missing values
  q <- exp((1/nobs)*sum(q.resid,na.rm=TRUE))
  
  #now evaluate sigma
  resid <- rep(NA,length(obs))
  for (iobs in 1:length(obs))
  {
    #check if there was an observation this year. NB Change this!
    if (obs[iobs]!=-99)
    {
      resid[iobs] <- (log(obs[iobs])-log(q)-
                        log(mean(predict.bio[iobs:(iobs+1)],na.rm=TRUE)))^2
    }
  }
  sigma <- sqrt(sum(resid,na.rm=TRUE)/nobs)

  neg.log.like <- nobs*log(sigma) + nobs/2
  
  return(neg.log.like)
}  
  
nll.function <- function(predict.bio,obs,q,sigma) 
{
  log.like <- 0
  for (iobs in 1:length(obs))
   {
    #check if there was an observation this year. NB Change this!
    if (obs[iobs]!=-99)
    {
     predicted <- q*mean(predict.bio[iobs:(iobs+1)],na.rm=TRUE)
     log.like <- log.like + loglike(obs[iobs],predicted,sigma) 
    }
   }
  return(-1.*log.like)
}
  
loglike <- function(observed,predicted,sigma)
{
  log.like <- (-1./pi)*log(sigma) -
    (1./(2*sigma^2))*(log(observed)-log(predicted))^2
  return(log.like)
}
  
objfun <- function(params=c(log.startbio=1,log.r=-0.2,log.k=10,log.q=-1,
                          log.sigma=-0.2),catch=rep(0,10),obs=rep(10,10),flag=1)
{
  #transform parameters
  start.bio <- exp(params['log.startbio'])
  r <- exp(params['log.r'])
  k <- exp(params['log.k'])
  q <- exp(params['log.q'])
  sigma <- exp(params['log.sigma'])
  
  #predict the biomass time series
  predict.bio <- predict.bio(start.bio,catch,r,k)

  #compute the (negative log) likelihood function
  nll <- nll.function(predict.bio,obs,q,sigma)
  
  #return the value for the negative log-likelihood
  if (flag==1)
    {
    results <- NULL
    results$nll <- nll  
    results$predict.bio <- predict.bio
    return(results)
   }
  if (flag==0) return(nll)
}

#This version uses the concentrated likelihood where q and sigma are 
#solved analytically
objfun.conc <- function(params=c(log.startbio=1,log.r=-0.2,log.k=10),
                        catch=rep(0,10),obs=rep(10,10),flag=1)
{
  #transform parameters
  r <- exp(params[2])
  k <- exp(params[3])
  start.bio <- exp(params[1])*k
  
  
  #predict the biomass time series
  predict.bio <- predict.bio(start.bio,catch,r,k)
  
  #compute the (negative log) likelihood function
  nll <- nll.function.conc(predict.bio,obs)
  
  #return the value for the negative log-likelihood
  if (flag==1)
  {
    results <- NULL
    results$nll <- nll  
    results$predict.bio <- predict.bio
    return(results)
  }
  if (flag==0) return(nll)
}


TheData <- read.csv("Yellowtail.csv")
head(TheData)
init.params <- c(log.startbio=0.25,
                 log.r=log(0.5),
                 log.k=log(200))
#log.q=log(0.9),
#log.sigma=-0.2)
obs <- TheData$Survey..kg.tow.
catch <- TheData$Yield..kt.


#objfun(init.params,catch,obs,flag=1)
objfun.conc(init.params,catch,obs,flag=1)


#a <- optim(par=init.params,fn=objfun,catch=catch,obs=obs,flag=0)
a <- optim(par=init.params,fn=objfun.conc,catch=catch,obs=obs,flag=0)

#output the function value, estimated parameters, and predicted biomass
objfun.conc(a$par,catch,obs,flag=1)
exp(a$par)


#likelihood profile - for next meeting?
#
lprof.par <- 2
lprof.par.vals <- log(seq(0.1,0.8,0.05))
a <- NULL
for (iprof in 1:length(lprof.par.vals))
  a[[iprof]] <- optim(par=init.params,fn=objfun.lprof,catch=catch,obs=obs,
                      flag=0,lprof.par,lprof.par.vals[iprof])
#We need to turn objfun.conc into objfun.lprof!

