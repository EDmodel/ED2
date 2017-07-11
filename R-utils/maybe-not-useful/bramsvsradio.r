#------------------------------------------------------------------------------------------#
#     This function will build the statistics for the BRAMS dataset, based on the          #
# radiosonde network.                                                                      #
#------------------------------------------------------------------------------------------#
bramsvsradio = function(bramsgrid,bramsraw,obs){
   statout   = list()

   for (vari in obs$radio$vars){

      print (paste("   - Finding the statistics for",vari,"..."))
      vardim=dim(bramsraw[[vari]])

      #----- Save dimensions in scalar variables. -----------------------------------------#
      tmax    = vardim[1]
      zmax    = vardim[2]
      xymax   = vardim[3]*vardim[4]
      ilmax   = obs$radio$nplevs
      places  = obs$radio$places
      nplaces = length(places)
      
      #----- Convert the current variable to a matrix. ------------------------------------#
      this   = array(bramsraw[[vari]],dim=c(tmax,zmax,xymax))
      this   = this[,,bramsgrid$nnradio]
      presin = array(bramsraw$pres,dim=c(tmax,zmax,xymax))
      presin = presin[,,bramsgrid$nnradio]

      #----- Switch the order of the indices so the pressure is the leftmost one. ---------#
      this    = aperm(a=this,perm=c(2,1,3))
      presin  = aperm(a=presin,perm=c(2,1,3))
      presout = obs$radio$plevs

      #----- Find the interpolated array. -------------------------------------------------#
      inter   = interpol(varin=this,presin=presin,presout=presout)

      #------------------------------------------------------------------------------------#
      #     Swap the pressure and time again, and collapse it to a matrix, so the          #
      # computation of the diurnal cycle statistics can be done using qapply.              #
      #------------------------------------------------------------------------------------#
      inter   = aperm(a=inter,perm=c(2,1,3))
      inter   = matrix(inter,nrow=tmax,ncol=ilmax*nplaces)
      isvalid = is.finite(inter)

      dathour = hours(bramsraw$gtime)
      hday    = sort(unique(dathour))
      xhday   = paste("X",hday,sep="")
      nday    = length(hday)

      #----- Compute the diurnal cycle of some statistics. --------------------------------#
      print ("     * Computing the average diurnal cycle...")
      dmean = qapply(mat=inter ,bycol=TRUE,index=dathour,func=mean,na.rm=TRUE)
      print ("     * Computing the standard deviation of the diurnal cycle...")
      dsdev = qapply(mat=inter ,bycol=TRUE,index=dathour,func=sd  ,na.rm=TRUE)
      print ("     * Counting how many data we used for the diunal cycle...")
      dccnt = qapply(mat=isvalid,bycol=TRUE,index=dathour,func=sum ,na.rm=TRUE)

      #----- Discard hours that we can't compare. -----------------------------------------#
      print ("   - Discarding hours that are not available in the radiosonde data...")
      usehr  = match(obs$radio$xhday,xhday)
      nxday  = length(usehr)



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

      #----- Expand the matrices back to arrays. ------------------------------------------#
      statout[[vari]]$dmean = array(dmean,dim=c(nxday,ilmax,nplaces))
      statout[[vari]]$dsdev = array(dsdev,dim=c(nxday,ilmax,nplaces))
      statout[[vari]]$dccnt = array(dccnt,dim=c(nxday,ilmax,nplaces))

      dimnames(statout[[vari]]$dmean) = list(xhday,presout,places)
      dimnames(statout[[vari]]$dsdev) = list(xhday,presout,places)
      dimnames(statout[[vari]]$dccnt) = list(xhday,presout,places)

   }#end for (vari in obs$radio$vars)
   return(statout)
}#end function
#------------------------------------------------------------------------------------------#

