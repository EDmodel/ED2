#------------------------------------------------------------------------------------------#
#     This function will build the statistics for the BRAMS dataset, based on the auto-    #
# matic weather station network.                                                           #
#------------------------------------------------------------------------------------------#
bramsvsaws = function(bramsgrid,bramsraw,obs){
   statout   = list()

   for (vari in obs$aws$vars){

      print (paste("   - Finding the statistics for",vari,"..."))


      #----- Convert the current variable to a matrix. ------------------------------------#
      this    = matrix(bramsraw[[vari]],nrow=bramsraw$tmax,ncol=bramsraw$xmax*bramsraw$ymax)
      this   = this[,bramsgrid$nnaws]
      isvalid = is.finite(this)

      dathour = hours(bramsraw$gtime)
      hday    = sort(unique(dathour))
      xhday   = paste("X",hday,sep="")
      nday    = length(hday)


      #----- Compute the diurnal cycle of some statistics. --------------------------------#
      print ("     * Computing the average diurnal cycle...")
      dmean = qapply(mat=this ,bycol=TRUE,index=dathour,func=mean,na.rm=TRUE)
      print ("     * Computing the standard deviation of the diurnal cycle...")
      dsdev = qapply(mat=this ,bycol=TRUE,index=dathour,func=sd  ,na.rm=TRUE)
      print ("     * Counting how many data we used for the diunal cycle...")
      dccnt = qapply(mat=isvalid,bycol=TRUE,index=dathour,func=sum ,na.rm=TRUE)

      #----- Discard hours that we can't compare. -----------------------------------------#
      print ("   - Discarding hours that are not available in TRMM...")
      usehr  = match(obs$rg$xhday,xhday)
      xhday  = xhday[usehr]
      dmean  = dmean[usehr,]
      dsdev  = dsdev[usehr,]
      dccnt  = dccnt[usehr,]

      #----- Standardise all -Inf to NA. --------------------------------------------------#
      if (any(dccnt == 0) || any(! is.finite(dmean))){
         dmean[! is.finite(dmean) | dccnt == 0] = NA
      }
      if (any(dccnt == 0) || any(! is.finite(dsdev))){
         dsdev[! is.finite(dsdev) | dccnt == 0] = NA
      }

      #----- Transpose the matrices. ------------------------------------------------------#
      statout[[vari]]$dmean = t(dmean)
      statout[[vari]]$dsdev = t(dsdev)
      statout[[vari]]$dccnt = t(dccnt)

      dimnames(statout[[vari]]$dmean) = list(obs$aws$places,xhday)
      dimnames(statout[[vari]]$dsdev) = list(obs$aws$places,xhday)
      dimnames(statout[[vari]]$dccnt) = list(obs$aws$places,xhday)

   }#end for (vari in obs$aws$vars)
   return(statout)
}#end function
#------------------------------------------------------------------------------------------#

