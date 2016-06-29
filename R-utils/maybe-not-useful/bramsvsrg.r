#------------------------------------------------------------------------------------------#
#     This function will build the statistics for the BRAMS dataset, based on the rain     #
# gauge network.                                                                           #
#------------------------------------------------------------------------------------------#
bramsvsrg = function(bramsgrid,bramsraw,obs){
   statout   = list()

   #----- Convert precipitation dataset to a matrix. --------------------------------------#
   prate     = matrix(bramsraw$prate,nrow=bramsraw$tmax,ncol=bramsraw$xmax*bramsraw$ymax)
   prate     = prate[,bramsgrid$nnrg]

   print ("   - Finding the precipitation rates...")
   dathour = hours(bramsraw$gtime)
   hday     = sort(unique(dathour))
   xhday    = paste("X",hday,sep="")
   nday     = length(hday)

   rainlog             = log(prate)
   itrains             = is.finite(rainlog)
   rainlog[! itrains ] = NA
   raindata            = is.finite(prate)
   #----- Compute the diurnal cycle of some statistics. -----------------------------------#
   print ("   - Computing the average diurnal cycle...")
   dmean = qapply(mat=rainlog ,bycol=TRUE,index=dathour,func=mean,na.rm=TRUE)
   print ("   - Computing the standard deviation of the diurnal cycle...")
   dsdev = qapply(mat=rainlog ,bycol=TRUE,index=dathour,func=sd  ,na.rm=TRUE)
   print ("   - Counting how many data we used for the diunal cycle...")
   dccnt = qapply(mat=raindata,bycol=TRUE,index=dathour,func=sum ,na.rm=TRUE)
   print ("   - Computing the diurnal cycle of probability of precipitation...")
   dprob = qapply(mat=itrains ,bycol=TRUE,index=dathour,func=sum ,na.rm=TRUE) / dccnt

   #----- Discard hours that we can't compare. --------------------------------------------#
   print ("   - Discarding hours that are not available in TRMM...")
   usehr  = match(obs$rg$xhday,xhday)
   xhday  = xhday[usehr]
   dmean  = dmean[usehr,]
   dsdev  = dsdev[usehr,]
   dccnt  = dccnt[usehr,]
   dprob  = dprob[usehr,]

   #----- Standardise all -Inf to NA. -----------------------------------------------------#
   dmean[! is.finite(dmean) | dccnt == 0] = NA
   dsdev[! is.finite(dsdev) | dccnt == 0] = NA
   dprob[! is.finite(dprob) | dccnt == 0] = NA
   dmean = exp(dmean)
   dsdev = exp(dsdev)

   #----- Transpose the matrices. ---------------------------------------------------------#
   statout$dmean = t(dmean)
   statout$dsdev = t(dsdev)
   statout$dprob = t(dprob)
   statout$dccnt = t(dccnt)

   dimnames(statout$dmean) = list(NULL,xhday)
   dimnames(statout$dsdev) = list(NULL,xhday)
   dimnames(statout$dprob) = list(NULL,xhday)
   dimnames(statout$dccnt) = list(NULL,xhday)

   return(statout)
}#end function
#------------------------------------------------------------------------------------------#

