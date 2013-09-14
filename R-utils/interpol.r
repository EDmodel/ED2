#==========================================================================================#
#==========================================================================================#
#    This function will perform interpolation of multi-dimensional arrays.  The result     #
# will be an array with the same size for dimensions other than the first, which will      #
# become the same size as the number of the output y.                                      #
#------------------------------------------------------------------------------------------#
interpol <<- function(x,y,xout,along=1,is.log=FALSE,...){

   #----- A bit of sanity check. ----------------------------------------------------------#
   remember = match.call()
   error    = NULL
   if (! ( is.array(y) || is.matrix(y) || is.data.frame(y)) ){
      error = c(error," 'y' must be an array, matrix, or data.frame!")
   }else if (! along %in% seq_along(dim(y))){
      error = c(error," 'along' must be a valid dimension for 'y'!")
   }else if (dim(y)[along] < 1){
      error = c(error," size of 'y' at dimension 'along' must at least 2!")
   }#end if
   if (is.array(x) || is.matrix(x) || is.data.frame(x)){
      if (any(dim(y) != dim(x))){
         error = c(error," 'x' and 'y' dimensions must match!")
      }#end if
      #------------------------------------------------------------------------------------#
   }else{
      #----- Coerce x to a vector. --------------------------------------------------------#
      x = unlist(x)
      #------------------------------------------------------------------------------------#
   }#end if
   if (along %in% seq_along(dim(y))){
      if (length(x) != dim(y)[along]){
         error = c(error," Size of 'x' doesn't match along-th dimension of y")
      }#end if
      #------------------------------------------------------------------------------------#
   }#end if
   if (length(xout) < 1){
      error = c(error," Size of 'xout' must be greater than 2")
   }#end if
   if (length(error) > 0){
      cat("Call: ",remember,"\n")
      cat("Errors: ","\n")
      for (e in seq_along(error)) cat(error[e],"\n")
      stop("Problems in function interpol")
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Save instructions to re-organise y. ---------------------------------------------#
   orig.dims = seq_along(dim(y))
   temp.dims = c(along,orig.dims[-along])
   back.dims = match(orig.dims,temp.dims)
   #---------------------------------------------------------------------------------------#




   #----- Send interp. dimension to first dimension, adjust names accordingly. ------------#
   y.temp             = aperm(a=y,perm=temp.dims)
   dim.yout           = c(length(xout),dim(y.temp)[-1])
   dimnames.yout      = dimnames(y.temp)
   dimnames.yout[[1]] = xout
   nr                 = dim(y.temp)[1]
   nc                 = prod(dim(y.temp)[-1])
   #---------------------------------------------------------------------------------------#


   #----- If we should use log to interpolate, transform variables now. -------------------#
   if (is.log){
      x.temp    = log(x)
      xout.temp = log(xout)
   }else{
      x.temp    = x
      xout.temp = xout
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Grab the arguments from the ellipsis, and append "xout".                           #
   #---------------------------------------------------------------------------------------#
   dotdotdot = modifyList(x=list(...),val=list(xout=xout.temp))
   #---------------------------------------------------------------------------------------#



   #----- Make sure "x.temp" has the same size as "y.temp". -------------------------------#
   if (length(x.temp) != length(y.temp)){
      x.temp = array(x.temp,dim=dim(y.temp),dimnames=dimnames(y.temp))
   }#end if
   #---------------------------------------------------------------------------------------#




   #----- Transform data into a matrix, then into a list. ---------------------------------#
   y.temp = matrix(y.temp,nrow=nr,ncol=nc)
   x.temp = matrix(x.temp,nrow=nr,ncol=nc)
   f.temp = col(y.temp)
   y.temp = split(x=y.temp,f=f.temp)
   x.temp = split(x=x.temp,f=f.temp)
   #---------------------------------------------------------------------------------------#




   #----- Perform interpolation. ----------------------------------------------------------#
   yo.temp = mapply( FUN      = interpolfun
                   , x        = x.temp
                   , y        = y.temp
                   , MoreArgs = dotdotdot
                   , SIMPLIFY = TRUE
                   )#end mapply
   #---------------------------------------------------------------------------------------#



   #----- Transform data back to the array. -----------------------------------------------#
   if (dim(yo.temp)[2] == length(xout)) yo.temp = t(yo.temp)
   yo.temp = array(data=yo.temp,dim=dim.yout,dimnames=dimnames.yout)
   yo.temp = aperm(a = yo.temp, perm = back.dims)
   #---------------------------------------------------------------------------------------#



   #----- Anything that went wrong becomes . -----------------------------------------------#
   if (dim(yo.temp)[2] == length(xout)) yo.temp = t(yo.temp)
   yo.temp = array(data=yo.temp,dim=dim.yout,dimnames=dimnames.yout)
   yo.temp = aperm(a = yo.temp, perm = back.dims)
   #---------------------------------------------------------------------------------------#



   #----- Free memory. --------------------------------------------------------------------#
   rm(x.temp,f.temp.y.temp)
   #---------------------------------------------------------------------------------------#


   #----- Free memory. --------------------------------------------------------------------#
   return(yo.temp)
   #---------------------------------------------------------------------------------------#
}# end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      This is a auxiliary function that runs spline only if it is safe to do so.          #
#------------------------------------------------------------------------------------------#
interpolfun <<- function(x,y,xout,silent=TRUE,...){
   sel = is.finite(x) & is.finite(y)
   if (any(sel)){
      yout = splinefun(x=x[sel],y=y[sel],...)(xout)
   }else{
      if (! silent) warning("Zero non-NA points")
      yout = xout + NA
   }#end interpol
   return(yout)
}#end function interpolfun
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      This is a auxiliary function that finds the best value of x given y using spline    #
# to approximate the function.                                                             #
#------------------------------------------------------------------------------------------#
interpol.root <<- function( x
                          , y
                          , yout
                          , extrapolate = FALSE
                          , silent      = TRUE
                          , tol         = sqrt(.Machine$double.eps)
                          , maxit       = 20
                          ,...
                          ){

   #----- Ignore iterative attempt if the user doesn't want to extrapolate. ---------------#
   if (! extrapolate) maxit = 1
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #      Check whether a solition is possible.                                            #
   #---------------------------------------------------------------------------------------#
   sel = is.finite(x) & is.finite(y)
   if (any(sel) && is.finite(yout)){
      #----- Find the span of x variable, and iterate until we get a good interval. -------#
      xrange  = range(x[sel])
      delta   = diff (xrange)
      iterate = TRUE
      it      = 0
      while (iterate && it < maxit){
         #----- Update iteration count and the range. -------------------------------------#
         it         = it + 1
         xrange.now = xrange + (it != 1)*c(-1,1)*2^(it-1)*delta
         #---------------------------------------------------------------------------------#


         #----- Fail-safe root finding. ---------------------------------------------------#
         ans        = try( expr   = uniroot( f        = splinefun.root
                                           , interval = xrange.now
                                           , yout     = yout
                                           , xin      = x[sel]
                                           , yin      = y[sel]
                                           , tol      = tol
                                           , ...
                                           )#end uniroot
                         , silent = silent
                         )#end try
         iterate     = "try-error" %in% is(ans)
         #---------------------------------------------------------------------------------#
      }#end while
      #------------------------------------------------------------------------------------#

      #------------------------------------------------------------------------------------#
      #     If "iterate" is TRUE, then it means we failed finding a solution.  In case     #
      # extrapolate is FALSE, we coerce to the edges, but return a failure sign, otherwise #
      # we return NA.                                                                      #
      #------------------------------------------------------------------------------------#
      if (iterate && ! extrapolate){
         #----- Evaluate the root function at the edges, and pick the closest one. --------#
         off = splinefun.root(xout=xrange,yout=yout,xin=x[sel],yin=y[sel],...)
         if (! silent) warning("No good solution was found, coercing to the best edge")
         xout = xrange[which.min(abs(off))]
         #---------------------------------------------------------------------------------#
      }else if (iterate){
         #----- Failed finding answer, and extrapolation was off, return NA. --------------#
         if (! silent) warning("No good solution was found, returning NA")
         xout = yout + NA
         #---------------------------------------------------------------------------------#
      }else{
         #----- Found an answer, use it. --------------------------------------------------#
         xout = ans$root
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#
   }else{
      #----- Everything is NA, return NA as well. -----------------------------------------#
      if ( (! silent) && (! any(sel))        ) warning("Zero non-NA points")
      if ( (! silent) && (! is.finite(yout)) ) warning("Output y is not defined")
      xout = yout + NA
      #------------------------------------------------------------------------------------#
   }#end if (any(sel))
   #---------------------------------------------------------------------------------------#
   return(xout)
}#end function interpolfun
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      This is the inverse function of splinefun.                                          #
#------------------------------------------------------------------------------------------#
splinefun.root <<- function(xout,yout,xin,yin,...){
   ans = splinefun(x=xin,y=yin,...)(xout)-yout
   return(ans)
}#end function splinefun.root
#==========================================================================================#
#==========================================================================================#
