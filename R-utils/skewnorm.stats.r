n.min.skew.norm <<- 15
sn.version      <<- as.numeric(substring(packageVersion(pkg="sn"),1,1))


#==========================================================================================#
#==========================================================================================#
#    This function finds all three parameters for skew normal distribution.                #
#------------------------------------------------------------------------------------------#
sn.stats <<- function(x,na.rm=FALSE,maxit=9999){
   #----- Stop if package fGarch isn't loaded. --------------------------------------------#
   if (! "package:sn" %in% search()){
      stop("Function sn.stats requires package sn!")
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- If na.rm is TRUE, we select only the data that is finite. -----------------------#
   if (na.rm){
      sel = is.finite(x)
      xsel = x[sel]
      rm(sel)
   }else{
      xsel = x
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Package sn is very different depending on version.                                 #
   #---------------------------------------------------------------------------------------#
   if (sn.version == 0){


      #------------------------------------------------------------------------------------#
      #     The default is to try the maximum likelihood estimator.  It is fast, but it    #
      # doesn't always converge.  In case it fails, we fall back to EM algorithm, which is #
      # more robust but a lot slower.  If none of them work, then falls back to Gaussian,  #
      # as this usually happens when there are very few points, or all points being the    #
      # same.                                                                              #
      #------------------------------------------------------------------------------------#
      myfit.mle = try(sn.mle(y=xsel,plot.it=FALSE,control=list(maxit=maxit)),silent=TRUE)
      if ("try-error" %in% is(myfit.mle)){
         warning(" - MLE failed, falling back to EM")
         converged = FALSE
      }else{
         dp         = cp.to.dp(myfit.mle$cp)
         ans        = c(dp,0)
         rm(dp)
         converged  = myfit.mle$optim$convergence == 0
      }#end if
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Run EM if it failed.                                                           #
      #------------------------------------------------------------------------------------#
      if (! converged){
         myfit.em = try(sn.em(y=xsel),silent=TRUE)
         if ("try-error" %in% is(myfit.em)){
            #------------------------------------------------------------------------------#
            #      Not enough point, assume Gaussian and warn the user...                  #
            #------------------------------------------------------------------------------#
            warning("EM failed too, using mean instead.")
            ans = c(mean(xsel,na.rm=TRUE),sd(xsel,na.rm=TRUE),0.,1.)
            #------------------------------------------------------------------------------#
         }else{
            ans = c(myfit.em$dp,0)
         }#end if
         rm(myfit.em)
      }#end if
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Delete everything else.                                                        #
      #------------------------------------------------------------------------------------#
      rm(xsel,myfit.mle,converged)
      #------------------------------------------------------------------------------------#
   }else{
      #------------------------------------------------------------------------------------#
      #     Version 1 uses the skew elliptical error distribution to estimate parameters.  #
      #------------------------------------------------------------------------------------#
      myfit.mle = try( selm( formula = xsel~1
                           , family  = "SN"
                           , data    = data.frame(xsel=xsel)
                           )#end selm
                     , silent = TRUE
                     )#end try
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      if ("try-error" %in% is(myfit.mle)){
         warning(" - MLE failed, using mean instead.")
         ans = c(mean(xsel,na.rm=TRUE),sd(xsel,na.rm=TRUE),0.,1.)
      }else{
         dp         = myfit.mle@param$dp
         ans        = c(dp,0)
         rm(dp)
      }#end if
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Delete everything else.                                                        #
      #------------------------------------------------------------------------------------#
      rm(xsel,myfit.mle)
      #------------------------------------------------------------------------------------#
   }#end if (sn.version == 0)
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Return the answer.                                                                #
   #---------------------------------------------------------------------------------------#
   names(ans) = c("location","scale","shape","gaussian")
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function sn.stats
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#    This function finds the location parameter for a skew normal distribution.  In case   #
# the distribution is not skewed, location parameter becomes the mean.                     #
#------------------------------------------------------------------------------------------#
sn.location <<- function(x,na.rm=FALSE,maxit=9999){
   ans = sn.stats(x,na.rm,maxit)["location"]
   #---------------------------------------------------------------------------------------#
   #     Return the answer.                                                                #
   #---------------------------------------------------------------------------------------#
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#    This function finds the scale parameter for a skew normal distribution.  In case      #
# the distribution is not skewed, scale parameter becomes the standard deviation.          #
#------------------------------------------------------------------------------------------#
sn.scale <<- function(x,na.rm=FALSE,maxit=9999){
   ans = sn.stats(x,na.rm,maxit)["scale"]


   #---------------------------------------------------------------------------------------#
   #     Return the answer.                                                                #
   #---------------------------------------------------------------------------------------#
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#    This function finds the shape parameter for a skew normal distribution.  In case      #
# the distribution is not skewed, shape parameter becomes 0.                               #
#------------------------------------------------------------------------------------------#
sn.shape <<- function(x,na.rm=FALSE,maxit=9999){
   ans = sn.stats(x,na.rm,maxit)["shape"]


   #---------------------------------------------------------------------------------------#
   #     Return the answer.                                                                #
   #---------------------------------------------------------------------------------------#
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#    This function finds the shape parameter for a skew normal distribution.  In case      #
# the distribution is not skewed, location parameter becomes the standard deviation.       #
#------------------------------------------------------------------------------------------#
sn.converged <<- function(x,na.rm=FALSE){
   ans = sn.stats(x,na.rm,maxit)["gaussian"] == 0


   #---------------------------------------------------------------------------------------#
   #     Return the answer.                                                                #
   #---------------------------------------------------------------------------------------#
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     This function normalises the dataset, using the skewed normal distribution.  The     #
# output values is the normalised value equivalent to the normal distribution.             #
#------------------------------------------------------------------------------------------#
skew2normal <<- function(x,location,scale,shape,idx=rep(1,times=length(x))){
   #----- Stop if package fGarch isn't loaded. --------------------------------------------#
   if (! "package:sn" %in% search()){
      stop("Function skew2normal requires package sn!")
   }#end if
   #---------------------------------------------------------------------------------------#

   #------ Find the size of x. ------------------------------------------------------------#
   nx = length(x)
   #---------------------------------------------------------------------------------------#


   #------ Find the number of possible indices. -------------------------------------------#
   unique.idx = unique(idx)
   nidx       = length(unique.idx)
   #---------------------------------------------------------------------------------------#



   #----- Stop if number of shapes . ------------------------------------------------------#
   if ((length(location) != nidx) | (length(scale) != nidx) | (length(shape) != nidx)){
      stop(paste(" Length of statistics must match the number of unique indices: ("
                ,nidx," in this case.",sep=""))
   }#end if
   #---------------------------------------------------------------------------------------#



   #------ Find the cumulative density function. ------------------------------------------#
   cdf.skew      = x * NA
   for (n in sequence(nidx)){
      sel          = is.finite(x) & idx == unique.idx[n]
      stats.ok     = is.finite(location[n]) && is.finite(scale[n]) && is.finite(shape[n])
      if (any(sel) && stats.ok && R.major <= 2){
         cdf.skew[sel] = psn(x[sel],location=location[n],scale=scale[n],shape=shape[n])
      }else if (any(sel) && stats.ok){
         cdf.skew[sel] = psn(x[sel],xi=location[n],omega=scale[n],alpha=shape[n])
      }#end if
   }#end for
   #---------------------------------------------------------------------------------------#



   #----- Find the equivalent quantile for the Gaussian distribution. ---------------------#
   xnorm    = qnorm(p=cdf.skew,mean=0,sd=1)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Return the answer.                                                                #
   #---------------------------------------------------------------------------------------#
   rm(nx,unique.idx,nidx,cdf.skew,sel,stats.ok)
   return(xnorm)
   #---------------------------------------------------------------------------------------#
}#end function
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     This function normalises the dataset, using the skewed normal distribution.  The     #
# output values is the normalised value equivalent to the normal distribution.             #
#------------------------------------------------------------------------------------------#
normal2skew <<- function(xnorm,location,scale,shape,idx=rep(1,times=length(xnorm))){
   #----- Stop if package fGarch isn't loaded. --------------------------------------------#
   if (! "package:sn" %in% search()){
      stop("Function normal2skew requires package sn!")
   }#end if
   #---------------------------------------------------------------------------------------#


   #------ Find the number of possible indices. -------------------------------------------#
   unique.idx = unique(idx)
   nidx       = length(unique.idx)
   #---------------------------------------------------------------------------------------#



   #----- Stop if number of shapes . ------------------------------------------------------#
   if ((length(location) != nidx) | (length(scale) != nidx) | (length(shape) != nidx)){
      stop(paste(" Length of statistics must match the number of unique indices: (",nidx
                ," in this case.",sep=""))
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Find the cumulative distribution function of the normalised quantiles. ----------#
   cdf.skew = pnorm(q=xnorm,mean=0,sd=1)
   x        = NA * cdf.skew
   for (n in 1:nidx){
      sel      = is.finite(cdf.skew) & idx == unique.idx[n]
      stats.ok = is.finite(location[n]) && is.finite(scale[n]) && is.finite(shape[n])
      if (any(sel) && stats.ok && R.major <= 2){
         x[sel] = qsn(p=cdf.skew[sel],location=location[n],scale=scale[n],shape=shape[n])
      }else if (any(sel) && stats.ok){
         x[sel] = qsn(p=cdf.skew[sel],xi=location[n],omega=scale[n],alpha=shape[n])
      }#end if
   }#end for
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Return the answer.                                                                #
   #---------------------------------------------------------------------------------------#
   rm(unique.idx,nidx,cdf.skew,sel,stats.ok)
   return(x)
   #---------------------------------------------------------------------------------------#
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#       Estimate the probability of multiple independent events given the distribution.    #
#------------------------------------------------------------------------------------------#
psn.mult <<- function(z,n,nsample=10000,location=0,scale=1,shape=0,lower.tail=TRUE
                     ,log.p=FALSE,force.sample=FALSE){


   if ( n == 1 && ! force.sample){
      #----- Use the standard probability distribution function (no sampling). ------------#
      if (R.major <= 2){
         prob = psn(x=z,location=location,scale=scale,shape=shape,lower.tail=lower.tail
                   ,log.p=log.p)
      }else{
         prob = psn(x=z,xi=location,omega=scale,alpha=shape,lower.tail=lower.tail
                   ,log.p=log.p)
      }#end if (R.major <= 2)
      #------------------------------------------------------------------------------------#
   }else{
      #----- Create an array with multiple data sets. -------------------------------------#
      if (R.major <= 2){
         x    = rsn(n=n*nsample,location=location,scale=scale,shape=shape)
      }else{
         x    = rsn(n=n*nsample,xi=location,omega=scale,alpha=shape)
      }#end if
      i    = rep(sequence(nsample),each=n)
      ztry = tapply(X=x,INDEX=i,FUN=sum)
      #------------------------------------------------------------------------------------#


      #----- Estimate probability based on the number of successes. -----------------------#
      if (lower.tail){
         prob  = sum(ztry <= z) / nsample
      }else{
         prob  = sum(ztry >= z) / nsample
      }#end if (lower.tail
      prob.min = 0 + 0.5 / nsample
      prob.max = 1 - 0.5 / nsample
      prob     = pmax(prob.min,pmin(prob.max,prob))
      if (log.p) prob = log(prob)
      #------------------------------------------------------------------------------------#



      #----- Free memory. -----------------------------------------------------------------#
      rm(x,i,ztry,prob.min,prob.max)
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Return probability. -------------------------------------------------------------#
   return(prob)
   #---------------------------------------------------------------------------------------#
}#end function psn.mult
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#       Auxiliary function, to find the root of the cumulative probability function.       #
#------------------------------------------------------------------------------------------#
psn.root <<- function(x,z,n,p,shift="location",...){
   short = substring(tolower(shift),1,2)
   if ( short == "lo"){
      ans = psn.mult(z,n,location=x,...) - p
   }else if ( short == "sc" ){
      ans = psn.mult(z,n,scale=x,...) - p
   }else if ( short == "sh" ){
      ans = psn.mult(z,n,shape=x,...) - p
   }else{
      cat(" - Invalid statistics request! Shift has been set to ",shift,"\n")
      cat("   Acceptable options are location, scale, and shape","\n")
      stop("Invalid statistics!")
   }#end if

   #----- Return value. -------------------------------------------------------------------#
   rm(short)
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function finds the location parameter that makes the probability of z with n    #
# independent approaches to become p.                                                      #
#------------------------------------------------------------------------------------------#
qsn.mult <<- function( z
                     , n
                     , p
                     , location
                     , scale
                     , shape
                     , shift    = "location"
                     , lower
                     , upper
                     , nsample  = 10000
                     , ...
                     ){





   #----- Find the basic statistics in terms of the original parameters. ------------------#
   delta         = shape / sqrt(1. + shape^2)
   old.mean      = location     + scale * delta * sqrt(2/pi)
   old.sdev      = scale * sqrt( 1 - 2*delta^2/pi)
   old.skew      = 0.5 * (4 - pi) * (delta*sqrt(2/pi) / sqrt(1-2*delta^2/pi))^3
   #---------------------------------------------------------------------------------------#



   #----- Find the location parameter that solves P(Z) = p. -------------------------------#
   tol        = 0.1 / nsample
   short      = substring(tolower(shift),1,2)
   if ( short == "lo"){
      sol        = uniroot(f=psn.root,lower=lower,upper=upper,z=z,n=n,p=p,scale=scale
                          ,shape=shape,shift=shift,nsample=nsample,tol=tol,maxiter = 1000)


      #----- Find the change in relative and absolute mean and location. ------------------#
      new.location  = sol$root
      new.mean      = new.location + scale * delta * sqrt(2/pi)
      sol$dlocation = ( new.location - location)
      sol$dmean     = ( new.mean     - old.mean)
      #------------------------------------------------------------------------------------#


      #----- Free memory. -----------------------------------------------------------------#
      rm(new.location,new.mean)
      #------------------------------------------------------------------------------------#
   }else if (short == "sc"){
      sol        = uniroot(f=psn.root,lower=lower,upper=upper,z=z,n=n,p=p,location=location
                          ,shape=shape,shift=shift,nsample=nsample,tol=tol,maxiter = 1000)


      #----- Find the change in relative and absolute std. dev. and scale. ----------------#
      new.scale      = sol$root
      new.sdev       = new.scale * sqrt( 1 - 2*delta^2/pi)
      sol$dscale     = ( new.scale - scale)
      sol$dsdev      = ( new.sdev - old.sdev)
      #------------------------------------------------------------------------------------#


      #----- Free memory. -----------------------------------------------------------------#
      rm(new.scale,new.sdev)
      #------------------------------------------------------------------------------------#
   }else if (short == "sh"){
      sol        = uniroot(f=psn.root,lower=lower,upper=upper,z=z,n=n,p=p,location=location
                          ,scale=scale,shift=shift,nsample=nsample,tol=tol,maxiter = 1000)

      #----- Find the change in relative and absolute std. dev. and scale. ----------------#
      new.shape      = sol$root
      new.delta      = new.shape / sqrt(1. + new.shape^2)
      new.skew       = 0.5 * (4 - pi) * (new.delta*sqrt(2/pi) / sqrt(1-2*new.delta^2/pi))^3
      sol$dshape     = ( new.shape - shape   )
      sol$dskew      = ( new.skew  - old.skew)
      #------------------------------------------------------------------------------------#


      #----- Free memory. -----------------------------------------------------------------#
      rm(new.shape,new.delta,new.skew)
      #------------------------------------------------------------------------------------#
   }else{
      cat(" - Invalid statistics request! Shift has been set to ",shift,"\n")
      cat("   Acceptable options are location, scale, and shape","\n")
      stop("Invalid statistics!")
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Free memory. --------------------------------------------------------------------#
   rm(delta,old.mean,old.sdev,old.skew)
   #---------------------------------------------------------------------------------------#


   #----- Return solution. ----------------------------------------------------------------#
   ans        = unlist(sol)
   names(ans) = names(sol)
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end lsn
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      This function is the equivalent to maximum likelihood estimator for least squares   #
# in the normal distribution, but assuming that the residuals have skew normal             #
# distribution with mean = 0.                                                              #
#------------------------------------------------------------------------------------------#
sn.lsq <<- function(r,skew=FALSE){

   #---------------------------------------------------------------------------------------#
   #     If all data are missing, return NA.                                               #
   #---------------------------------------------------------------------------------------#
   if (all(is.na(r))){
      ans = NA
      return(ans)
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Keep only the good data. --------------------------------------------------------#
   r.fine = r[! is.na(r)]
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Decide which PDF to use (normal or skew normal.                                   #
   #---------------------------------------------------------------------------------------#
   if (skew){
      rstats    = sn.stats(r.fine)
      rscale    = rstats["scale"]
      rshape    = rstats["shape"]
      rdelta    = rshape / sqrt(1. + rshape^2)
      #----- This is the location that produces mean 0. -----------------------------------#
      rlocation = - rscale * rdelta * sqrt(2/pi)
      #------------------------------------------------------------------------------------#


      #------ Find the log-likelihood of the distribution. --------------------------------#
      if (R.major <= 2){
         ans = sum(x=dsn(x=r.fine,location=rlocation,scale=rscale,shape=rshape,log=TRUE))
      }else{
         ans = sum(x=dsn(x=r.fine,xi=rlocation,omega=rscale,alpha=rshape,log=TRUE))
      }#end if (R.major <= 2)
      #------------------------------------------------------------------------------------#

      #----- Free memory. -----------------------------------------------------------------#
      rm(rstats,rscale,rshape,rdelta,rlocation)
      #------------------------------------------------------------------------------------#
   }else{
      #------ Find the log-likelihood of the distribution. --------------------------------#
      ans = sum(x=dnorm(x=r.fine,mean=0,sd=sd(r.fine),log=TRUE))
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Free memory. --------------------------------------------------------------------#
   rm(r.fine)
   #---------------------------------------------------------------------------------------#


   #----- Free memory and return solution. ------------------------------------------------#
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end sn.lsq
#==========================================================================================#
#==========================================================================================#
