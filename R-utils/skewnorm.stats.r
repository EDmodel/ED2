n.min.skew.norm <<- 15


#==========================================================================================#
#==========================================================================================#
#    This function finds all three parameters for skew normal distribution.                #
#------------------------------------------------------------------------------------------#
sn.stats = function(x,na.rm=FALSE,maxit=9999){
   #----- Stop if package fGarch isn't loaded. --------------------------------------------#
   if (! "package:sn" %in% search()){
      stop("Function sn.stats requires package sn!")
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- If na.rm is TRUE, we select only the data that is finite. -----------------------#
   if (na.rm){
      sel = is.finite(x)
      xsel = x[sel]
   }else{
      xsel = x
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     The default is to try the maximum likelihood estimator.  It is fast, but it       #
   # doesn't always converge.  In case it fails, we fall back to EM algorithm, which is    #
   # more robust but a lot slower.  If none of them work, then falls back to Gaussian, as  #
   # this usually happens when there are very few points, or all points being the same.    #
   #---------------------------------------------------------------------------------------#
   myfit.mle = try(sn.mle(y=xsel,plot.it=FALSE,control=list(maxit=maxit)),silent=TRUE)
   if ("try-error" %in% is(myfit.mle)){
      warning(" - MLE failed, falling back to EM")
      converged = FALSE
   }else{
      dp         = cp.to.dp(myfit.mle$cp)
      ans        = c(dp,0)
      converged  = myfit.mle$optim$convergence == 0
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Run EM if it failed.                                                              #
   #---------------------------------------------------------------------------------------#
   if (! converged){
      myfit.em = try(sn.em(y=xsel),silent=TRUE)
      if ("try-error" %in% is(myfit.em)){
         #---------------------------------------------------------------------------------#
         #      Not enough point, assume Gaussian and warn the user...                     #
         #---------------------------------------------------------------------------------#
         warning("EM failed too, using mean instead.")
         ans = c(mean(xsel,na.rm=TRUE),sd(xsel,na.rm=TRUE),0.,1.)
         #---------------------------------------------------------------------------------#
      }else{
         ans = c(myfit.em$dp,0)
      }#end if
   }#end if
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
sn.location = function(x,na.rm=FALSE,maxit=9999){
   #----- Stop if package fGarch isn't loaded. --------------------------------------------#
   if (! "package:sn" %in% search()){
      stop("Function sn.location requires package sn!")
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- If na.rm is TRUE, we select only the data that is finite. -----------------------#
   if (na.rm){
      sel = is.finite(x)
      xsel = x[sel]
   }else{
      xsel = x
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     The default is to try the maximum likelihood estimator.  It is fast, but it       #
   # doesn't always converge.  In case it fails, we fall back to EM algorithm, which is    #
   # more robust but a lot slower.  If none of them work, then falls back to Gaussian, as  #
   # this usually happens when there are very few points, or all points being the same.    #
   #---------------------------------------------------------------------------------------#
   myfit.mle = try(sn.mle(y=xsel,plot.it=FALSE,control=list(maxit=maxit)),silent=TRUE)
   if ("try-error" %in% is(myfit.mle)){
      warning(" - MLE failed, falling back to EM")
      converged = FALSE
   }else{
      dp        = cp.to.dp(myfit.mle$cp)
      ans       = dp[1]
      converged = myfit.mle$optim$convergence == 0
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Run EM if it failed.                                                              #
   #---------------------------------------------------------------------------------------#
   if (! converged){
      myfit.em = try(sn.em(y=xsel),silent=TRUE)
      if ("try-error" %in% is(myfit.em)){
         #---------------------------------------------------------------------------------#
         #      Not enough point, assume Gaussian and warn the user...                     #
         #---------------------------------------------------------------------------------#
         warning("EM failed too, using mean instead.")
         ans = mean(xsel,na.rm=TRUE)
         #---------------------------------------------------------------------------------#
      }else{
         ans = myfit.em$dp[1]
      }#end if
   }#end if
   #---------------------------------------------------------------------------------------#


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
sn.scale = function(x,na.rm=FALSE,maxit=9999){
   #----- Stop if package fGarch isn't loaded. --------------------------------------------#
   if (! "package:sn" %in% search()){
      stop("Function sn.scale requires package sn!")
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- If na.rm is TRUE, we select only the data that is finite. -----------------------#
   if (na.rm){
      sel = is.finite(x)
      xsel = x[sel]
   }else{
      xsel = x
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     The default is to try the maximum likelihood estimator.  It is fast, but it       #
   # doesn't always converge.  In case it fails, we fall back to EM algorithm, which is    #
   # more robust but a lot slower.  If none of them work, then falls back to Gaussian, as  #
   # this usually happens when there are very few points, or all points being the same.    #
   #---------------------------------------------------------------------------------------#
   myfit.mle = try(sn.mle(y=xsel,plot.it=FALSE,control=list(maxit=maxit)),silent=TRUE)
   if ("try-error" %in% is(myfit.mle)){
      warning(" - MLE failed, falling back to EM")
      converged = FALSE
   }else{
      dp        = cp.to.dp(myfit.mle$cp)
      ans       = dp[2]
      converged = myfit.mle$optim$convergence == 0
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Run EM if it failed.                                                              #
   #---------------------------------------------------------------------------------------#
   if (! converged){
      myfit.em = try(sn.em(y=xsel),silent=TRUE)
      if ("try-error" %in% is(myfit.em)){
         #---------------------------------------------------------------------------------#
         #      Not enough point, assume Gaussian and warn the user...                     #
         #---------------------------------------------------------------------------------#
         warning("EM failed too, using standard deviation instead.")
         ans = sd(xsel,na.rm=TRUE)
         #---------------------------------------------------------------------------------#
      }else{
         ans = myfit.em$dp[2]
      }#end if
   }#end if
   #---------------------------------------------------------------------------------------#


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
sn.shape = function(x,na.rm=FALSE,maxit=9999){
   #----- Stop if package fGarch isn't loaded. --------------------------------------------#
   if (! "package:sn" %in% search()){
      stop("Function sn.shape requires package sn!")
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- If na.rm is TRUE, we select only the data that is finite. -----------------------#
   if (na.rm){
      sel = is.finite(x)
      xsel = x[sel]
   }else{
      xsel = x
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     The default is to try the maximum likelihood estimator.  It is fast, but it       #
   # doesn't always converge.  In case it fails, we fall back to EM algorithm, which is    #
   # more robust but a lot slower.  If none of them work, then falls back to Gaussian, as  #
   # this usually happens when there are very few points, or all points being the same.    #
   #---------------------------------------------------------------------------------------#
   myfit.mle = try(sn.mle(y=xsel,plot.it=FALSE,control=list(maxit=maxit)),silent=TRUE)
   if ("try-error" %in% is(myfit.mle)){
      warning(" - MLE failed, falling back to EM")
      converged = FALSE
   }else{
      dp        = cp.to.dp(myfit.mle$cp)
      ans       = dp[3]
      converged = myfit.mle$optim$convergence == 0
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Run EM if it failed.                                                              #
   #---------------------------------------------------------------------------------------#
   if (! converged){
      myfit.em = try(sn.em(y=xsel),silent=TRUE)
      if ("try-error" %in% is(myfit.em)){
         #---------------------------------------------------------------------------------#
         #      Not enough point, assume Gaussian and warn the user...                     #
         #---------------------------------------------------------------------------------#
         warning("EM failed too, using Gaussian (shape 0) instead.")
         ans = 0.
         #---------------------------------------------------------------------------------#
      }else{
         ans = myfit.em$dp[3]
      }#end if
   }#end if
   #---------------------------------------------------------------------------------------#


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
sn.converged = function(x,na.rm=FALSE){
   #----- Stop if package fGarch isn't loaded. --------------------------------------------#
   if (! "package:sn" %in% search()){
      stop("Function sn.converged requires package sn!")
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     The very least number of valid input data for skewed normal is 3.                 #
   #---------------------------------------------------------------------------------------#
   if (sum(is.finite(x)) < n.min.skew.norm){
      ans = FALSE
      return(ans)
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- If na.rm is TRUE, we select only the data that is finite. -----------------------#
   if (na.rm){
      sel = is.finite(x)
      xsel = x[sel]
   }else{
      xsel = x
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Fit the skewed normal distribution.                                              #
   #---------------------------------------------------------------------------------------#
   if (length(xsel) > n.min.skew.norm & sum(is.finite(xsel)) > n.min.skew.norm){
      #----- Possible to run the maximum likelihood method. -------------------------------#
      myfit = sn.em(y=xsel)
      ans   = myfit$logL
   }else{
      warning("Not enough data to determine skewed normal distribution.")
      ans = FALSE
   }#end if
   #---------------------------------------------------------------------------------------#


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
skew2normal = function(x,location,scale,shape,idx=rep(1,times=length(x))){
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
   for (n in 1:nidx){
      sel          = is.finite(x) & idx == unique.idx[n]
      stats.ok     = is.finite(location[n]) && is.finite(scale[n]) && is.finite(shape[n])
      if (any(sel) && stats.ok){
         cdf.skew[sel] = psn(x[sel],location=location[n],scale=scale[n],shape=shape[n])
      }#end if
   }#end for
   #---------------------------------------------------------------------------------------#



   #----- Find the equivalent quantile for the Gaussian distribution. ---------------------#
   xnorm    = qnorm(p=cdf.skew,mean=0,sd=1)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Return the answer.                                                                #
   #---------------------------------------------------------------------------------------#
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
normal2skew = function(xnorm,location,scale,shape,idx=rep(1,times=length(x))){
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
      if (any(sel) && stats.ok){
         x[sel] = qsn(p=cdf.skew[sel],location=location[n],scale=scale[n],shape=shape[n])
      }#end if
   }#end for
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Return the answer.                                                                #
   #---------------------------------------------------------------------------------------#
   return(x)
   #---------------------------------------------------------------------------------------#
}#end function
#==========================================================================================#
#==========================================================================================#
