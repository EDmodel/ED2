#==========================================================================================#
#==========================================================================================#
#      This is the generic function used for a handful of parameters in Collatz et al.     #
# (1991), written as a function of a reference temperature, given by tcollatz.  The output #
# variable will have the same unit as the pre-factor coefficient.                          #
#------------------------------------------------------------------------------------------#
collatz <<- function(temp,refval,q10){
   #---------------------------------------------------------------------------------------#
   #     If the exponential factor is tiny, make it zero, otherwise compute the actual     #
   # function.                                                                             #
   #---------------------------------------------------------------------------------------#
   coll  = refval  * q10 ^ (fcoll * (temp - tphysref))
   #---------------------------------------------------------------------------------------#

   return(coll)
}#end function collatz
#==========================================================================================#
#==========================================================================================#
