#==========================================================================================#
#==========================================================================================#
#     Function DENSITY.SAFE -- this is just a wrapper for "density" function.  In case of  #
#                              a single element, we assign NAs.  We don't fix the values   #
#                              of anything other than $y, this is just a quick and cheap   #
#                              fix.                                                        #
#------------------------------------------------------------------------------------------#
density.safe <<- function(x,weights=NULL,...){
   sel     = is.finite(x)
   x.use   = x[sel]

   n.x.use = length(x.use)

   if (n.x.use == 0){
      x.use       = c(-1,1)
      weights.use = c(0.5,0.5)
      ans         = density(x=x.use,weights=weights.use,...)
   }else if (n.x.use == 1){
      x.use       = rep(x.use,times=2)
      weights.use = c(1,0)
      ans         = density(x=x.use,weights=weights.use,...)
   }else if (! is.null(weights)){
      weights.use = weights[sel] / sum(weights[sel])
      ans         = density(x=x.use,weights=weights.use,...)
   }else{
      weights.use = NULL
      ans         = density(x=x.use,...)
   }#end if

   if (n.x.use == 0) ans$y = ans$y + NA
   return(ans)
}#end function
#==========================================================================================#
#==========================================================================================#
