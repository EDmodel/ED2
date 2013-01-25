#------------------------------------------------------------------------------------------#
#     This function will build the support, or log of the likelihood, based on automatic   #
# weather station statistics.  For each variable:                                          #
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
lnlikeaws = function(brams,obs){

   lnlike = list()

   for (v in 1:obs$aws$nvars){
      vari = obs$aws$vars[v]

      if (vari == "psfc"){
         mfac = 100.0
         afac = 0.000
      }else if (vari == "rvsfc"){
         mfac = 0.001
         afac = 0.000
      }else{
         mfac = 1.000
         afac = 0.000
      }#end if

      #----- For the likelihood, we are assuming normal distribution. ---------------------#
      n  = t(obs$aws[[vari]]$dccnt)
      m  = t(obs$aws[[vari]]$dmean)
      s  = t(obs$aws[[vari]]$dsdev)
      NN = brams$aws[[vari]]$dccnt
      MM = mfac * brams$aws[[vari]]$dmean + afac
      SS = mfac * brams$aws[[vari]]$dsdev + afac


      #---------------------------------------------------------------------------------------#
      #   Now find the likelihood for each point independently.                               #
      #---------------------------------------------------------------------------------------#
      lnlike[[vari]] = n * (- log(SS) - (s^2 + (m - MM)^2)/(2*SS^2))
      lnlike[[vari]] = lnlike[[vari]] - max(lnlike[[vari]],na.rm=TRUE)
   }#end if

   return(lnlike)
}#end function
#------------------------------------------------------------------------------------------#

