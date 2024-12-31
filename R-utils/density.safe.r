#==========================================================================================#
#==========================================================================================#
#     Function DENSITY.SAFE.                                                               #
#                                                                                          #
#     This is very similar to R's original "density" function, but it is less likely to    #
# crash.  Unlike regular density, this function allows incorporating the effect of         #
# uncertainty in x entries (variable e).  Currently only normal and log-normal             #
# distributions are possible.                                                              #
#     In case the series doesn't have any valid entries, we return NA in the y axis.       #
#------------------------------------------------------------------------------------------#
density.safe <<- function( x
                         , e       = NULL
                         , weights = NULL
                         , ne      = 1000
                         , dist    = c("norm","lnorm")
                         , xmin    = -Inf
                         , xmax    = +Inf
                         , n       = 512
                         ,...
                         ){


   #----- Standardise distribution. -------------------------------------------------------#
   if (! is.null(e)){
      dist = match.arg(dist)
      if (ne < 30){
         stop(" Variable 'ne' must be at least 30 if running density.safe with errors.")
      }#end if (ne < 30)
   }#end if (! is.null(e))
   #---------------------------------------------------------------------------------------#



   #----- Select only valid entries. ------------------------------------------------------#
   if (is.null(e)){
      #----- Select the valid x values. ---------------------------------------------------#
      sel     = is.finite(x)
      x.use   = x[sel]
      n.x.use = sum(sel)
      #------------------------------------------------------------------------------------#

      #------ Copy weights a new data set. ------------------------------------------------#
      if (is.null(weights)){
         weights.use = NULL
      }else{
         weights.use = weights[sel] / sum(weights[sel])
      }#end if (is.null(weights))
      #------------------------------------------------------------------------------------#

   }else{
      #----- Select the valid x and e values. ---------------------------------------------#
      sel     = (x %gt% 0) & (e %ge% 0)
      n.x.use = sum(sel)
      #------------------------------------------------------------------------------------#

      #------------------------------------------------------------------------------------#
      #     Decide which distribution to use.                                              #
      #------------------------------------------------------------------------------------#
      if (dist %in% "norm"){
         x.use = rnorm(n=n.x.use*ne,mean=x[sel],sd=e[sel])
      }else if (dist %in% "lnorm"){
         x.use = x[sel] * rlnorm(n=n.x.use*ne,meanlog=0,sdlog=sqrt(1 + (e[sel]/x[sel])^2))
      }#end if (dist %in% "rnorm")
      x.use = pmax(xmin,pmin(xmax,x.use))
      #------------------------------------------------------------------------------------#

      #------ Copy weights a new data set. ------------------------------------------------#
      if (is.null(weights)){
         weights.use = NULL
      }else{
         weights.use = rep(weights[sel],times=ne)
         weights.use = weights.use / sum(weights.use)
      }#end if (is.null(weights))
      #------------------------------------------------------------------------------------#
   }#end if (is.null(e))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Decide what to leave in the output based on the number of valid entries.          #
   #---------------------------------------------------------------------------------------#
   if (n.x.use == 0){
      x.use       = c(-1,1)
      weights.use = c(0.5,0.5)
      ans         = density(x=x.use,weights=weights.use,n=n,...)
   }else if (n.x.use == 1 && is.null(e)){
      x.use       = rep(x.use,times=2)
      weights.use = c(1,0)
      ans         = density(x=x.use,weights=weights.use,n=n,...)
   }else if (! is.null(weights)){
      ans         = density(x=x.use,weights=weights.use,n=n,...)
   }else{
      weights.use = NULL
      ans         = density(x=x.use,n=n,...)
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Ensure that the values outside the acceptable range are transferred to the        #
   # acceptable range.  This is to avoid density to be non-zero in non-physical ranges     #
   # (for example, negative CO2 concentration).                                            #
   #---------------------------------------------------------------------------------------#
   if (is.finite(xmin)){
      bye         = ans$x %lt% xmin
      xadd        = min(ans$x[! bye],na.rm=TRUE)
      iadd        = which(ans$x %eq% xadd)
      ans$y[iadd] = ans$y[iadd] + sum(ans$y[bye])
      ans$y[bye ] = 0.
   }#end if(is.finite(xmin))
   if (is.finite(xmax)){
      bye         = ans$x %gt% xmax
      xadd        = max(ans$x[! bye],na.rm=TRUE)
      iadd        = which(ans$x %eq% xadd)
      ans$y[iadd] = ans$y[iadd] + sum(ans$y[bye])
      ans$y[bye ] = 0.
   }#end if(is.finite(xmax))
   #---------------------------------------------------------------------------------------#


   #----- Return NA in case there weren't enough values to determine the KDE. -------------#
   if (n.x.use == 0) ans$y = ans$y + NA_real_
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function
#==========================================================================================#
#==========================================================================================#
