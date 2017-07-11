#------------------------------------------------------------------------------------------#
#     This function will build the support, or log of the likelihood, based on radiosonde  #
# statistics.  For each variable:                                                          #
#                                                                                          #
#    Support is computed as:                                                               #
#                                                                                          #
#  S = k + n * { - ln(SS) - [s^2+(m-MM)^2]/[2*SS^2]}                                       #
#                                                                                          #
# where:                                                                                   #
#                                                                                          #
# k  - arbitrary constant (it will be assigned later, here it will be 0.                   #
# n  - number of observations.                                                             #
# m  - mean log of non-zero precipitation, based on observations.                          #
# s  - standard deviation of the log of non-zero precipitation, based on observations.     #
# NN - number of points used by the model statistics.                                      #
# MM - model mean of log of non-zero precipitation.                                        #
# SS - model standard deviation of log of non-zero precipitation.                          #
#------------------------------------------------------------------------------------------#
lnlikeradio = function(brams,obs){

   lnlike = list()

   for (v in 1:obs$radio$nvars){
      vari = obs$radio$vars[v]

      if (vari == "rvap"){
         mfac = 0.001
         afac = 0.000
      }else{
         mfac = 1.000
         afac = 0.000
      }#end if

      #----- For the likelihood, we are assuming normal distribution. ---------------------#
      n  = obs$radio[[vari]]$dccnt
      m  = obs$radio[[vari]]$dmean
      s  = obs$radio[[vari]]$dsdev
      NN = brams$radio[[vari]]$dccnt
      MM = mfac * brams$radio[[vari]]$dmean + afac
      SS = mfac * brams$radio[[vari]]$dsdev + afac

      #----- Make the values above 300 hPa NA, since the observations are not good. -------#
      if (vari == "rvap"){
         sel = obs$radio$plevs < 300.
         n [,sel,] = 0
         m [,sel,] = NA
         s [,sel,] = NA
         NN[,sel,] = 0
         MM[,sel,] = NA
         SS[,sel,] = NA
      }
      #------------------------------------------------------------------------------------#
      #   Now find the likelihood for each point independently.                            #
      #------------------------------------------------------------------------------------#
      lnlike[[vari]] = n * (- log(SS) - (s^2 + (m - MM)^2)/(2*SS^2))
      lnlike[[vari]] = lnlike[[vari]] - max(lnlike[[vari]],na.rm=TRUE)
   }#end if

   return(lnlike)
}#end function
#------------------------------------------------------------------------------------------#

