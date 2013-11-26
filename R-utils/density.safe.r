#==========================================================================================#
#==========================================================================================#
#     Function DENSITY.SAFE -- this is just a wrapper for "density" function.  In case of  #
#                              a single element, we assign NAs.  We don't fix the values   #
#                              of anything other than $y, this is just a quick and cheap   #
#                              fix.                                                        #
#------------------------------------------------------------------------------------------#
density.safe <<- function(x,weights=NULL,...){
   nx = length(x)

   if (! is.null(weights)) weights = weights / sum(weights)

   if (nx == 1){
     x = rep(x,times=2)
     if (! is.null(weights)) weights=rep(weights,times=2)
   }#end if

   ans = density(x=x,weights=weights,...)
   if (nx == 1) ans$y = ans$y + NA
   return(ans)
}#end function
#==========================================================================================#
#==========================================================================================#
