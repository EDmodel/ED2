#==========================================================================================#
#==========================================================================================#
#     Function that extracts a certain bit position of a numeric value.                    #
#------------------------------------------------------------------------------------------#
bit.extract <<- function(x,n){
   dummy = stopifnot(is.numeric(x))
   dummy = stopifnot(n >= 1L & n <= 32L)
   nx    = length(x)
   ans   = rep(NA,times=nx)
   for (i in sequence(nx)) ans[i] = as.integer(intToBits(x=x[i])[n])
   return(ans)
}#end function bit.extract
bit.extract <<- cmpfun(bit.extract)
#==========================================================================================#
#==========================================================================================#
