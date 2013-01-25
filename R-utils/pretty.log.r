#==========================================================================================#
#==========================================================================================#
#     Function that creates a pretty scale in log co-ordinates.                            #
#------------------------------------------------------------------------------------------#
pretty.log = function(x,base=10,n=10){
   vlevels     = pretty(x=log(x,base=base),n=n)
   nact        = length(vlevels)
   vlevels     = base^vlevels

   power.base  = base^(floor(log(x=vlevels,base=base)))
   npow        = length(unique(power.base))
   ndigits     = ceiling(log(nact/npow,base=base))
   vlevels     = round(vlevels / power.base,digits=ndigits) * power.base
   return(vlevels)
}#end function
#==========================================================================================#
#==========================================================================================#

