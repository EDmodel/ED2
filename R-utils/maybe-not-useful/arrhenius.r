#==========================================================================================#
#==========================================================================================#
#      This is the Arrhenius function, written as a function of a reference temperature,   #
# given by tarrh8.  The output variable will have the same unit as the pre-factor          #
# coefficient.                                                                             #
#------------------------------------------------------------------------------------------#
arrhenius = function(temp,prefactor,expcoeff){

   #---------------------------------------------------------------------------------------#
   #     Find the term that goes to the exponential term, and check its size.  This is to  #
   # avoid floating point exceptions due to overflow or underflow.                         #
   #---------------------------------------------------------------------------------------#
   lnexp = expcoeff * (tarrhi - 1.0/temp)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     If the exponential factor is tiny, make it zero, otherwise compute the actual     #
   # function.                                                                             #
   #---------------------------------------------------------------------------------------#
   sel      = lnexp >= lnexp.min
   arr      = lnexp
   arr[ sel] = prefactor * exp(lnexp[sel])
   arr[!sel] = 0.0
   #---------------------------------------------------------------------------------------#

   return(arr)
}#end function arrhenius
#==========================================================================================#
#==========================================================================================#
