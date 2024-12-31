#==========================================================================================#
#==========================================================================================#
#    Density function of the left truncated Weibull distribution.                          #
#------------------------------------------------------------------------------------------#
dlt.weibull <<- function(x,left=0,shape,scale=1,log=FALSE){
   #------ Aliases for non-dimensional variables. -----------------------------------------#
   xnorm =     x / scale
   lnorm =  left / scale
   #---------------------------------------------------------------------------------------#

   dens = ifelse( x %lt% left
                , 0.
                , shape / scale * xnorm ^ (shape-1) * exp( lnorm ^ shape - xnorm ^ shape )
                )#end ifelse
   if (log) dens = log(dens)
   return(dens)
}#end dlt.weibull
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#    Cumulative function of the left truncated Weibull distribution.                       #
#------------------------------------------------------------------------------------------#
plt.weibull <<- function(q,left=0,shape,scale=1,lower.tail=TRUE,log.p=FALSE){
   #------ Aliases for non-dimensional variables. -----------------------------------------#
   qnorm =     q / scale
   lnorm =  left / scale
   #---------------------------------------------------------------------------------------#


   prob = ifelse( q %lt% left, 0., 1. - exp( lnorm^shape - qnorm^shape ) )
   if (lower.tail) prob = 1. - prob
   if (log.p     ) prob = log(prob)
   return(prob)
}#end plt.weibull
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#    Quantile of the left truncated Weibull distribution.                                  #
#------------------------------------------------------------------------------------------#
qlt.weibull <<- function(p,left=0,shape,scale=1,lower.tail=TRUE,log.p=FALSE){
   if (log.p     ) p = exp(p)
   if (lower.tail) p = 1 - p
   quant = scale * ( (left/scale)^shape - log(1-p))^(1/shape)
   return(quant)
}#end qlt.weibull
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#    Random number generator of the left truncated Weibull distribution.                   #
#------------------------------------------------------------------------------------------#
rlt.weibull <<- function(n,left=0,shape,scale=1){
   p     = runif(n,min=0,max=1)
   quant = qlt.weibull(p=p,left=left,shape=shape,scale=scale)
   return(quant)
}#end rlt.weibull
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#    Maximum likelihood estimator for left-truncated Weibull distribution.                 #
#------------------------------------------------------------------------------------------#
mle.lt.weibull <<- function(x,left=0,start=list(shape=1,scale=1),count.min=5){
   #----- Check whether the mle will work this time. --------------------------------------#
   fine = length(unique(x)) > count.min
   if (fine){
      now  = try( fitdistr(x=x,densfun=dlt.weibull,start=start,left=left),silent=TRUE)
      fine = ! ("try-error" %in% is(now))
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- In case of success, return the answer, otherwise return NA. ---------------------#
   if (fine){
      ans = unlist(now)
   }else{
      r   = rlt.weibull(n=100,shape=2,scale=left,left=left)
      now = try( fitdistr(x=r,densfun=dlt.weibull,start=start,left=left),silent=TRUE)
      if ("try-error" %in% is(now)) browser()
      ans = unlist(now) * NA
   }#end if
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end mle.lt.weibull
#==========================================================================================#
#==========================================================================================#
