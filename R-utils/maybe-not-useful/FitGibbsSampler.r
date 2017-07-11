#==========================================================================================#
#==========================================================================================#
#     The following functions are small adaptations of the supporting online material      #
# available at:                                                                            #
#                                                                                          #
# Condit, R., et al., 2006: The importance of demographic niches to tree diversity.        #
#      Science, 313, 98-101.                                                               #
#                                                                                          #
#     The main modification is to allow mortality and growth to be based on either number  #
# of individuals or biomass.                                                               #
#                                                                                          #
# (Below is the original text from Condit et al. (2006).                                   #
# Functions for fitting a hyperdistribution of mortality or growth rates with a log-normal #
# distribution using the Gibbs sampler. These functions will (should) source in R2.1.0 as  #
# is (no add-on packages necessary), and given a table of data in the correct format,      #
# will produce the posterior distribution of all parameters.                               #
#                                                                                          #
# Included functions                                                                       #
#                                                                                          #
# (1)  metrop.lnormMort.Gibbs                                                              #
# (2)  mu.mortGibbs                                                                        #
# (3)  sd.mortGibbs                                                                        #
# (4)  spmean.mortGibbs                                                                    #
# (5)  metrop1step                                                                         #
# (6)  lnorm.growth.Gibbs                                                                  #
# (7)  GibbsCommVar                                                                        #
# (8)  GibbsSppVar                                                                         #
# (9)  GibbsCommMean                                                                       #
# (10) GibbsTheta                                                                          #
#                                                                                          #
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
# FUNCTION 1 -- The main function for mortality, accepting data and producing the Gibbs    #
#               sampler.  Functions 2-5 are its subroutines. Run this function alone for   #
#               mortality data.                                                            #
#                                                                                          #
#    The table of data has one row per species, and must have the following columns:       #
#    -> species names (or codes) as rownames;                                              #
#    -> a column N with the number of individuals in the initial census;                   #
#    -> a column S with the number of survivors in the later census;                       #
#    -> a column time with a time interval for each species;                               #
#    -> a column rate with a calculated mortality rate, (log(N)-log(S))/time.              #
#    -> There should be no missing values, with all N > 0.                                 #
#                                                                                          #
#    The two hyperparameters are the mean (logMu) and SD (logSD) of log(rate).  The        #
# argument start.param sets their starting values.                                         #
#                                                                                          #
#    Argument unity sets a typical value for one individual.  This should be one when      #
# using mortality rate based on individual count, but it should be the biomass in case     #
# mortality rate is found as biomass loss.                                                 #
#                                                                                          #
#    The argument div sets the step size used by the Metropolis algorithm for choosing new #
# values of each parameter.  The default values worked for the Barro Colorado plot, but    #
# will have to be adjusted for other datasets.                                             #
#                                                                                          #
#    The number of cycles run by the Gibbs sampler is the argument steps, and the current  #
# state will be printed to the screen every showstep steps unless verbose is set to FALSE. #
# The argument burn.in.phase is used to discard the first steps (burn-in-phase); estimates #
# based on the Gibbs sampler are drawn from the steps after this.  The acceptance rate is  #
# printed as well, allowing easy adjustment of the step size (argument div).               #
#                                                                                          #
#    Data are output as a list, including all post-burn-in values of every parameter plus  #
# their means, medians, and quantiles.                                                     #
#------------------------------------------------------------------------------------------#
metrop.lnormMort.Gibbs <<- function( data
                                   , outrate
                                   , start.param   = c(-3,.8)
                                   , div           = c(.35,.35,.15)
                                   , steps         = 1000
                                   , burn.in.phase = 100
                                   , verbose       = TRUE
                                   , showstep      = 100
                                   ){

   #---- Find the number of elements for PDF and CDF. -------------------------------------#
   nrate = length(outrate)

   #---- Find the indices of the post burn-in phase. --------------------------------------#
   whichuse = seq(from=1,to=steps,by=1) + burn.in.phase
   steps    = steps + burn.in.phase
   nuse     = length(whichuse)

   #---- Set matrices. --------------------------------------------------------------------#
   logMu=logSD = matrix(nrow=steps,ncol=2)
   logMu[1,1]  = start.param[1]
   logSD[1,1]  = start.param[2]

   nospp   = dim(data)[1]
   spmean  = matrix(ncol=steps,nrow=nospp)
   draws   = matrix(ncol=steps,nrow=nospp)
   rownames(spmean) = rownames(draws) = rownames(data)
   spmean[,1] = data$rate

   #----- Tweak mortality for the cases where there are no survivors or no deaths. --------#
   nomort           = which(data$S == 0. | data$N == data$S)
   spmean[nomort,1] = ( ( log(data$N[nomort]+ 2)
                        - log(data$S[nomort]+ 1) )
                      / data$time[nomort] )
   raresp = which(data$N < 5)[1]
   maxsp  = which.max(data$N)


   for (i in 2:steps){
      logMu[i,] = metrop1step( func       = mu.mortGibbs
                             ,start.param = logMu[i-1,1]
                             ,scale.param = div[1]
                             ,spmean      = spmean[,i-1]
                             ,logSD       = logSD[i-1,1]
                             )#end function metrop1step
      logSD[i,] = metrop1step( func        = sd.mortGibbs
                             , start.param = logSD[i-1,1]
                             , scale.param = div[2]
                             , spmean      = spmean[,i-1]
                             , logMu       = logMu[i,1]
                             )#end function metrop1step

      for (j in 1:nospp){
         nextmean = metrop1step( func        = spmean.mort.Gibbs
                               , start.param = spmean[j,i-1]
                               , scale.param = div[3]/data$N[j]^(.33)
                               , N           = data$N[j]
                               , S           = data$S[j]
                               , time        = data$time[j]
                               , logMu       = logMu[i,1]
                               , logSD       = logSD[i,1]
                               )#end function metrop1step

         spmean[j,i] = nextmean[1]
         draws[j,i]  = nextmean[2]
      }#end for

      if (verbose && (i == 2 | i %% showstep == 0)){
         cat("            -> Step ", i, ": ", round(logMu[i,1],2)
                                            , round(logSD[i,1],2)
                                            , round(spmean[raresp,i],4)
                                            , round(spmean[maxsp,i],4) )
         cat(   " .. Accept: "              , round(1-sum(logMu[,2],na.rm=T)/i,2)
                                            , round(1-sum(logSD[,2],na.rm=T)/i,2)
                                            , round(1-sum(draws[raresp,],na.rm=T)/i,2)
                                            , round(1-sum(draws[maxsp,],na.rm=T)/i,2)
                                            , "\n")
      }#end if
   }#end for
  
 spmean = spmean[,whichuse]
 logMu  = logMu[whichuse,1]
 logSD  = logSD[whichuse,1]


 
 mediantheta = apply(spmean,1,median)
 meantheta   = apply(spmean,1,mean)
 uppertheta  = apply(spmean,1,quantile,probs=.975)
 lowertheta  = apply(spmean,1,quantile,probs=.025)
 
 meanMu  = mean(logMu)
 meanRt  = exp(meanMu)
 meanSD  = mean(logSD)
 upperMu = quantile(logMu,probs=.975)
 upperRt = exp(upperMu)
 lowerMu = quantile(logMu,probs=.025)
 lowerRt = exp(lowerMu)
 upperSD = quantile(logSD,probs=.975)
 lowerSD = quantile(logSD,probs=.025)

 #----- Find mean CDF and PDF using the expected statistics. ------------------------------#
 cdf.mean = plnorm(q=outrate,meanlog=meanMu,sdlog=meanSD)
 pdf.mean = diff(cdf.mean)
 pdf.mean = c(pdf.mean,pdf.mean[nrate-1])
 
 #----- Find the cumulative distribution function for each member. ------------------------#
 logMu.mat = matrix(logMu,nrow=nrate,ncol=nuse,byrow=TRUE)
 logSD.mat = matrix(logSD,nrow=nrate,ncol=nuse,byrow=TRUE)
 rate.mat  = matrix(outrate,nrow=nrate,ncol=nuse,byrow=FALSE)
 cdf.mat   = matrix(plnorm(q=rate.mat,meanlog=logMu.mat,sdlog=logSD.mat)
                   ,nrow=nrate,ncol=nuse)
 pdf.mat   = apply(cdf.mat,MARGIN=2,FUN=diff)
 pdf.mat   = rbind(pdf.mat,pdf.mat[nrate-1,])
 
 cdf.q025  = tapply(X=cdf.mat,INDEX=row(cdf.mat),FUN=quantile,probs=0.025)
 cdf.q975  = tapply(X=cdf.mat,INDEX=row(cdf.mat),FUN=quantile,probs=0.975)
 pdf.q025  = tapply(X=pdf.mat,INDEX=row(pdf.mat),FUN=quantile,probs=0.025)
 pdf.q975  = tapply(X=pdf.mat,INDEX=row(pdf.mat),FUN=quantile,probs=0.975)
 
 distrib           = cbind(pdf.mean,pdf.q025,pdf.q975,cdf.mean,cdf.q025,cdf.q975)
 pcdf.names        = paste("pcr",sprintf("%.3f",100*outrate),sep=".")
 dimnames(distrib) = list (pcdf.names,c("pdf.median","pdf.lower","pdf.upper"
                                       ,"cdf.median","cdf.lower","cdf.upper") )

 ans     = list( fullMu    = logMu
               , fullSD    = logSD
               , draws     = draws
               , fulltheta = spmean
               , means     = c( Mu      = meanMu
                              , lowerMu = lowerMu
                              , upperMu = upperMu
                              , Rate    = meanRt
                              , lowerRt = lowerRt
                              , upperRt = upperRt
                              , SD      = meanSD
                              , lowerSD = lowerSD
                              , upperSD = upperSD
                              )#end c
               , table     = data.frame( N         = data$N
                                       , S         = data$S
                                       , time      = data$time
                                       , obsmean   = data$rate
                                       , fitmort   = meantheta
                                       , lowermean = lowertheta
                                       , uppermean = uppertheta
                                       )#end data.frame
               , distrib   = distrib
               )#end list
   return(ans)
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#   FUNCTION 2 -- This is the conditional likelihood function for the probability of       #
#                 observing the parameter logMu given the data, species means, and logSD.  #
#                 It is called each step of the Gibbs sampler.                             #
#------------------------------------------------------------------------------------------#
mu.mortGibbs <<- function(logMu,spmean,logSD){
  if (logSD <= 0) return(-Inf)
 
  llike = dlnorm(spmean,meanlog=logMu,sdlog=logSD,log=TRUE)
  return(sum(llike))
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#   FUNCTION 3 -- This is the conditional likelihood function for the probability of       #
#                 observing the parameter logSD given the data, species means, and logMu.  #
#                 It is used each step of the Gibbs sampler.                               #
#------------------------------------------------------------------------------------------#
sd.mortGibbs <<- function(logSD,spmean,logMu){
 if(logSD <= 0) return(-Inf)

 llike = dlnorm(spmean,meanlog=logMu,sdlog=logSD,log=TRUE)
 return(sum(llike))
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#   FUNCTION 4 -- This is the conditional likelihood function for the probability of       #
#                 observing the parameters spmean mortality parameter of each species)     #
#                 given the data, logMu, and logSD. It is used each step of the Gibbs      #
#                 sampler, once per each species.                                          #
#------------------------------------------------------------------------------------------#
spmean.mort.Gibbs <<- function(spmean,N,S,time,logMu,logSD){
   if(is.na(spmean)) browser()
   if(spmean <= 0) return(-Inf)
   
   theta = exp(-spmean*time)
   
   llike = ( dlnorm(x=spmean,meanlog=logMu,sdlog=logSD,log=TRUE)
           + dbinom(x=S,size=N,prob=theta,log=TRUE) )
   return(llike)
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
# FUNCTION 5 -- A generic routine for taking a single Metropolis step given any            #
#               probability function, func, an initial parameter, start.param, and a       #
#               scaling parameter for the step size, scale.param.  It returns the next     #
#               parameter value as well as an acceptance indicator (0 if the value was     #
#               accepted, 1 if rejected).                                                  #
#------------------------------------------------------------------------------------------#
metrop1step <<- function(func,start.param,scale.param,...){
   origlike = func(start.param,...)
   newval   = rnorm(1,mean=start.param,sd=scale.param)
   newlike  = func(newval,...)
   
   if (newlike >= origlike){
      return(c(newval,0))
   }else{
      likeratio=exp(newlike-origlike)
   }#end if

   if (runif(1) < likeratio){
      return(c(newval,0))
   }else{
      return(c(start.param,1))
   }#end if
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
# FUNCTION 6 -- The main function for growth, accepting data and producing the Gibbs       #
#               sampler.  Functions 7-10 are its subroutines. Run this function alone for  #
#               growth data.                                                               #
#                                                                                          #
#    The table of data has one row per individual, and must have the following columns:    #
#    -> a column sp with species names (or codes);                                         #
#    -> a column growth with the growth rate to be used.  If log-transformed data are      #
#       needed, this column should be the log-tranformation; it will not be transformed in #
#       this program.  There should be no missing values in either column.                 #
#                                                                                          #
#    The two hyperparameters are the mean (mu) and variance (comm.var) of the distribution #
# of mean growth rates across the community (the hyperdistribution).  The third parameter  #
# is the within-species variance (spp.var). The argument start.param sets the 3 starting   #
# values, or it can be NULL and the starting values are set automatically.                 #
#                                                                                          #
#    The argument prior sets the parameters defining the prior distributions.              #
#                                                                                          #
#    The number of cycles run by the Gibbs sampler is the argument steps, and the current  #
# state will be printed to the screen every showstep steps, unless verbose is FALSE.  The  #
# argument burn.in.phase defines the number of initial time steps (burn-in-phase) that     #
# will be discarded for the final analysis, so the Gibbs sampler are drawn from those      #
# steps after burn.in.phase. The acceptance rate is printed as well, allowing easy         #
# adjustment of the step size (argument div).                                              #
#                                                                                          #
#    Data are output as a list, including all post-burn-in values of every parameter plus  #
# their means, medians, and quantiles. A table of mean growth rates per species, plus the  #
# estimated mean from the Gibbs sampler, is also returned.                                 #
#------------------------------------------------------------------------------------------#
lnorm.growth.Gibbs <<- function( growth
                               , outrate
                               , start.param   = c(0,1,1)
                               , prior         = c(0,0,0,0,1e5,1e5)
                               , steps         = 1000
                               , burn.in.phase = 100
                               , verbose       = TRUE
                               , showstep      = 100
                               ){

   #---- Find the number of elements for PDF and CDF. -------------------------------------#
   nrate = length(outrate)

   #----- Define the post-burn-in phase. --------------------------------------------------#
   whichuse = burn.in.phase + seq(from=1,to=steps,by=1)
   steps    = steps + burn.in.phase
   nuse     = length(whichuse)

   gcol             = which(colnames(growth)=="growth")
   mu               = comm.var = spp.var = numeric()
   abund            = table(growth$sp)
   nospp            = length(abund)
   spmean           = matrix(nrow=nospp,ncol=steps)
   rownames(spmean) = names(abund)

   spmean[,1]       = obsmean = tapply(growth[,gcol],growth$sp,mean,na.rm=TRUE)
   obssd            = tapply(growth[,gcol],growth$sp,sd,na.rm=TRUE)
   m                = match(growth$sp,names(obsmean))
   growth$spmean    = obsmean[m]

   if (is.null(start.param)){
      mu      [1] = 0
      comm.var[1] = var(spmean[,1])
      spp.var [1] = 1
   }else{
    mu        [1] = start.param[1]
    comm.var  [1] = start.param[2]
    spp.var   [1] = start.param[3]
   }#end if

   a1     = prior[1]
   b1     = prior[2]
   a2     = prior[3]
   b2     = prior[4]
   sigma0 = prior[5]

   raresp = which(abund < 5)[1]
   maxsp  = which.max(abund)

   for (i in 2:steps){
      comm.var [i] = GibbsCommVar(a1,b1,spmean[,i-1],mu[i-1],abund)
      spp.var  [i] = GibbsSppVar(a2,b2,growth,spmean[,i-1])
      mu       [i] = GibbsCommMean(comm.var[i],sigma0,mu[1],spmean[,i-1],abund)
      spmean  [,i] = GibbsTheta(comm.var[i],spp.var[i],spmean[,1],mu[i],abund)
     
      if (verbose && (i == 2 | i %% showstep == 0)){
         cat("            -> Step ", i, ": ", round(mu[i],2)
                                            , round(comm.var[i]^.5,2)
                                            , round(spp.var[i]^.5,2)
                                            , round(spmean[raresp,i],4)
                                            , round(spmean[maxsp,i],4)
                                            , "\n" 
                                            )#end cat
      }#end if
   }#end for

   spmean  = spmean[,whichuse]
   comm.sd = comm.var[whichuse]^.5
   spp.sd  = spp.var[whichuse]^.5
   mu      = mu[whichuse]

   mediantheta = apply(spmean,1,median)
   meantheta   = apply(spmean,1,mean)
   uppertheta  = apply(spmean,1,quantile,prob=.975)
   lowertheta  = apply(spmean,1,quantile,prob=.025)

   meanCommSD  = mean(comm.sd)
   meanSppSD   = mean(spp.sd)

   meanMu       = mean(mu)
   upperMu      = quantile(mu,prob=.975)
   lowerMu      = quantile(mu,prob=.025)

   meanRt       = mean(mu)
   upperRt      = quantile(mu,prob=.975)
   lowerRt      = quantile(mu,prob=.025)

   upperCommSD  = quantile(comm.sd,prob=.975)
   lowerCommSD  = quantile(comm.sd,prob=.025)
   upperSppSD   = quantile(spp.sd,prob=.975)
   lowerSppSD   = quantile(spp.sd,prob=.025)


   #----- Find mean CDF and PDF using the expected statistics. ----------------------------#
   cdf.mean = pnorm(q=outrate,mean=meanMu,sd=meanCommSD)
   pdf.mean = diff(cdf.mean)
   pdf.mean = c(pdf.mean,pdf.mean[nrate-1])
 
   #----- Find the cumulative distribution function for each member. ----------------------#
   Mu.mat    = matrix(mu     ,nrow=nrate,ncol=nuse,byrow=TRUE)
   SD.mat    = matrix(comm.sd,nrow=nrate,ncol=nuse,byrow=TRUE)
   rate.mat  = matrix(outrate,nrow=nrate,ncol=nuse,byrow=FALSE)
   cdf.mat   = matrix(pnorm(q=rate.mat,mean=Mu.mat,sd=SD.mat)
                     ,nrow=nrate,ncol=nuse)
   pdf.mat   = apply(cdf.mat,MARGIN=2,FUN=diff)
   pdf.mat   = rbind(pdf.mat,pdf.mat[nrate-1,])
 
   cdf.q025  = tapply(X=cdf.mat,INDEX=row(cdf.mat),FUN=quantile,probs=0.025)
   cdf.q975  = tapply(X=cdf.mat,INDEX=row(cdf.mat),FUN=quantile,probs=0.975)
   pdf.q025  = tapply(X=pdf.mat,INDEX=row(pdf.mat),FUN=quantile,probs=0.025)
   pdf.q975  = tapply(X=pdf.mat,INDEX=row(pdf.mat),FUN=quantile,probs=0.975)
 
   distrib           = cbind(pdf.mean,pdf.q025,pdf.q975,cdf.mean,cdf.q025,cdf.q975)
   pcdf.names        = paste("pcr",sprintf("%.3f",100*outrate),sep=".")
   dimnames(distrib) = list (pcdf.names,c("pdf.median","pdf.lower","pdf.upper"
                                         ,"cdf.median","cdf.lower","cdf.upper") )


   result.table = data.frame( N         = abund
                            , obsmean   = obsmean
                            , obssd     = obssd
                            , fitgrow   = meantheta
                            , lowermean = lowertheta
                            , uppermean = uppertheta
                            )#end data.frame
   colnames(result.table)[1:2] = c("sp","N")
   rownames(result.table)      = NULL
 
   ans = list( fullCommSD = comm.sd
             , fullSppSD  = spp.sd
             , fullMu     = mu
             , fulltheta  = spmean
             , means      = c( Mu          = meanMu
                             , lowerMu     = lowerMu
                             , upperMu     = upperMu
                             , Rt          = meanRt
                             , lowerRt     = lowerRt
                             , upperRt     = upperRt
                             , commSD      = meanCommSD
                             , lowerCommSD = lowerCommSD
                             , upperCommSD = upperCommSD
                             , sppSD       = meanSppSD
                             , lowerSppSD  = lowerSppSD
                             , upperSppSD  = upperSppSD
                             )#end c
             , table      = result.table
             , distrib    = distrib
             )#end list
   return(ans)
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
# FUNCTION 7 --- Function that draws the next comm.var parameter in the Gibbs sampler.     #
#------------------------------------------------------------------------------------------#
GibbsCommVar <<- function(a,b,theta,mu,ab){
   shape = a + 0.5 * length(theta)
   scale = b + 0.5 * sum((theta - mu)^2)
 
   result = rinvgamma(1,shape,scale)
   if (is.na(result) | result < 1e-34) result = 1e-34
 
   return(result)
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
# FUNCTION 8 -- Function that draws the next spp.var parameter in the Gibbs sampler.       #
#------------------------------------------------------------------------------------------#
GibbsSppVar <<- function(a,b,growth,theta,J,measure="rgr"){
   shape = a + 0.5*dim(growth)[1]
   gcol  = which(colnames(growth)=="growth") 
 
   devsq = (growth[,gcol]-growth$spmean)^2
   scale = b + 0.5*sum(devsq)
   return(rinvgamma(1,shape,scale))
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
# FUNCTION 9 -- Function that draws the next comm.mean parameter in the Gibbs sampler.     #
#------------------------------------------------------------------------------------------#
GibbsCommMean <<- function(cvar,sigma0,mu,theta,ab){
   Nmean = ( cvar * mu + sigma0 * sum(theta)) / ( cvar + sigma0 * length(theta))
   Nvar  = cvar * sigma0 / ( cvar + sigma0 * length(theta))

   return(rnorm(1,mean=Nmean,sd=Nvar^.5))
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
# FUNCTION 10 -- Function that draws the next species mean in the Gibbs sampler, for a     #
#                single species.                                                           #
#------------------------------------------------------------------------------------------#
GibbsTheta <<- function(cvar,svar,yhat,mu,ab){
   Nmean = ab * cvar * yhat / (ab * cvar + svar) + svar * mu / ( ab * cvar + svar )
   Nvar  = cvar * svar / ( ab * cvar + svar )
 
   return(rnorm(length(yhat),mean=Nmean,sd=Nvar^.5))
}#end function
#==========================================================================================#
#==========================================================================================#


