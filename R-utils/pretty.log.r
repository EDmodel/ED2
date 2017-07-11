#==========================================================================================#
#==========================================================================================#
#     Function that creates a pretty scale in log co-ordinates.                            #
#------------------------------------------------------------------------------------------#
pretty.log = function(x,base=10,n=5,forcelog=FALSE){
   log.neat  = pretty(x=log(x,base=base),n=n)
   dlog.neat = median(diff(log.neat))
   neat      = base^log.neat
   nact      = length(neat)

   #---------------------------------------------------------------------------------------#
   #     In case it is a base 10, make it even prettier.  Also, in case the log scale is   #
   # very close to the actual scale, forget about the log scale and use regular pretty, it #
   # gives a more legible scale.                                                           #
   #---------------------------------------------------------------------------------------#
   if (base == 10 && dlog.neat %wr% c(0.1,0.5) && (! forcelog)){
      sel  = abs(log.neat - as.integer(log.neat)) %<% (0.5 * dlog.neat)
      tens = sort(c(neat[sel],base^(c(floor(min(log.neat)),ceiling(max(log.neat))))))

      #------------------------------------------------------------------------------------#
      #      Fix scale so it looks nicer if dtens is 2 or 3.                               #
      #------------------------------------------------------------------------------------#
      if (dlog.neat %<% 0.15){
         mult    = c(1,2,3,5,7)
      }else if (dlog.neat %<% 0.35){
         mult    = c(1,2,5)
      }else{
         mult    = c(1,3)
      }#end if
      vlevels   = sort(unique(c(mult %o% tens)))
      aa        = min(which(vlevels > min(neat)))-1
      zz        = max(which(vlevels < max(neat)))+1
      vlevels   = vlevels[aa:zz]
      still.log = TRUE
      #------------------------------------------------------------------------------------#
   }else if(dlog.neat %<% 0.1 && (! forcelog)){
      #-----  The plot is hardly log, use normal units instead. ---------------------------#
      vlevels = pretty(x=x,n=n)
      still.log = FALSE
      #------------------------------------------------------------------------------------#
   }else{
      vlevels   = neat
      still.log = TRUE
   }#end if
   #---------------------------------------------------------------------------------------#



   if (still.log){
      power.base  = base^(floor(log(x=vlevels,base=base)))
      npow        = length(unique(power.base))
      ndigits     = ceiling(log(nact/npow,base=base))
      vlevels     = round(vlevels / power.base,digits=ndigits) * power.base
   }#end if

   return(vlevels)
}#end function
#==========================================================================================#
#==========================================================================================#
