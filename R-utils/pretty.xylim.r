#==========================================================================================#
#==========================================================================================#
#     This function creates a expanded range for a variable, expanding the limits by a     #
# certain factor so it fits a legend.                                                      #
#------------------------------------------------------------------------------------------#
pretty.xylim <<- function(u,fracexp=0.00,is.log=FALSE){
   #----- Make fracexp a vector o length[2]. ----------------------------------------------#
   if (length(fracexp) == 1){
      if (fracexp < 0.0){
         fracexp = c(-fracexp,0.0)
      }else{
         fracexp = c(0.0,fracexp)
      }#end if
   }else{
      fracexp = abs(fracexp)
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Set lnudge, which is used to create limits when variable "u" is constant. -------#
   if (any(fracexp != 0.0)){
      lnudge = 0.5 * max(fracexp)
   }else{
      lnudge = 0.20
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Save the original warning, so it won't issue warnings during the selection. -----#
   warn.orig = getOption("warn")
   options(warn=-1)
   #---------------------------------------------------------------------------------------#



   #----- Select the data that we consider for the range. ---------------------------------#
   vec.u  = c(u)
   if (is.log) vec.u[! (vec.u %gt% 0)] = NA
   ulimit  = range(vec.u,finite=TRUE)
   #---------------------------------------------------------------------------------------#


   #----- Revert the warning option back to the original. ---------------------------------#
   options(warn=warn.orig)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Expand the axis.                                                                  #
   #---------------------------------------------------------------------------------------#
   if (is.log) ulimit = log(ulimit)
   dulimit = diff(ulimit)
   if (any(! is.finite(ulimit)) || (ulimit[1] == ulimit[2] && ulimit[1] == 0)){
      ulimit = c(-1,1)
   }else if (ulimit[1] == ulimit[2] ){
      ulimit[1] = ulimit[1] * ( 1. - sign(ulimit[1]) * lnudge)
      ulimit[2] = ulimit[2] * ( 1. + sign(ulimit[2]) * lnudge)
   }else{
      ulimit[1] = ulimit[1] - fracexp[1] * dulimit
      ulimit[2] = ulimit[2] + fracexp[2] * dulimit
   }#end if
   if (is.log) ulimit = exp(ulimit)
   #---------------------------------------------------------------------------------------#

   return(ulimit)
}#end pretty.xylim
#==========================================================================================#
#==========================================================================================#
