#------------------------------------------------------------------------------------------#
#     This function will build the statistics for the BRAMS dataset, based on TRMM grid.   #
#------------------------------------------------------------------------------------------#
bramsvstrmm = function(bramsgrid,bramsraw,obs){
   statout   = list()

   #----- Convert precipitation dataset to a matrix. --------------------------------------#
   prate     = matrix(bramsraw$prate,nrow=bramsraw$tmax,ncol=bramsraw$xmax*bramsraw$ymax)

   print ("   - Finding the precipitation rate statistics in each TRMM grid area...")
   aux     = qapply(mat=prate,index=bramsgrid$maptrmm,func=mean,bycol=FALSE,na.rm=TRUE)
   dathour = hours(bramsraw$gtime)
   hday     = sort(unique(dathour))
   xhday    = paste("X",hday,sep="")
   nday     = length(hday)

   rainlog             = log(aux)
   itrains             = is.finite(rainlog)
   rainlog[! itrains ] = NA
   raindata            = is.finite(aux  )
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
   usehr  = match(obs$trmm$xhday,xhday)
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
   dmean = t(dmean)
   dsdev = t(dsdev)
   dprob = t(dprob)
   dccnt = t(dccnt)

   br2tr = sort(unique(bramsgrid$maptrmm[is.finite(bramsgrid$maptrmm)]))

   #---------------------------------------------------------------------------------------#
   #     Make the matrices with the same size as TRMM, initialise them with NA, and fill   #
   # only those with valid BRAMS grid points in it.                                        #
   #---------------------------------------------------------------------------------------#
   statout$dmean = matrix(NA,nrow=dim(obs$trmm$dmean)[1],ncol=dim(obs$trmm$dmean)[2])
   statout$dsdev = matrix(NA,nrow=dim(obs$trmm$dsdev)[1],ncol=dim(obs$trmm$dsdev)[2])
   statout$dprob = matrix(NA,nrow=dim(obs$trmm$dprob)[1],ncol=dim(obs$trmm$dprob)[2])
   statout$dccnt = matrix(NA,nrow=dim(obs$trmm$dccnt)[1],ncol=dim(obs$trmm$dccnt)[2])

   dimnames(statout$dmean) = list(NULL,xhday)
   dimnames(statout$dsdev) = list(NULL,xhday)
   dimnames(statout$dprob) = list(NULL,xhday)
   dimnames(statout$dccnt) = list(NULL,xhday)

   statout$dmean[br2tr,] = dmean
   statout$dsdev[br2tr,] = dsdev
   statout$dprob[br2tr,] = dprob
   statout$dccnt[br2tr,] = dccnt


   return(statout)
}#end function
#------------------------------------------------------------------------------------------#

