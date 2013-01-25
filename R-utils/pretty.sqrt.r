#==========================================================================================#
#==========================================================================================#
#     Function that creates a pretty scale in square root.                                 #
#------------------------------------------------------------------------------------------#
pretty.sqrt    = function(x,n=10){
   vlevels     = pretty(x=sign(x)*sqrt(abs(x)),n=n)
   nact        = length(vlevels)
   vlevels     = abs(vlevels)*vlevels

   power.base  = 10^(floor(log10(x=abs(vlevels))))
   npow        = length(unique(power.base))
   ndigits     = ceiling(log10(nact/npow))
   vlevels     = ifelse( vlevels==0
                       , 0.
                       , round(vlevels / power.base,digits=ndigits) * power.base
                       )#end ifelse
   return(vlevels)
}#end function
#==========================================================================================#
#==========================================================================================#

