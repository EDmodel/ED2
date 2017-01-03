#==========================================================================================#
#==========================================================================================#
#     This function plots a scatter plot adding error bars.                                #
#------------------------------------------------------------------------------------------#
error.bar <<- function( x
                      , y
                      , xerr     = NULL
                      , yerr     = NULL
                      , xlow     = NULL
                      , xhigh    = NULL
                      , ylow     = NULL
                      , yhigh    = NULL
                      , xlim     = NULL
                      , ylim     = NULL
                      , cap      = 0.015
                      , main     = NULL
                      , sub      = NULL
                      , xlab     = as.character(substitute(x))
                      , ylab     = if (is.factor(x) || is.character(x)){
                                      ""
                                   }else{
                                      as.character(substitute(y))
                                   }#end if
                      , add      = FALSE
                      , lty      = 1
                      , type     = "p"
                      , lwd      = 1
                      , pch      = 16
                      , col      = "black"
                      , err.col  = col
                      , err.type = c("l","e")
                      , edensity = 60
                      , eangle   = runif(n=length(x),min=0,max=180)
                      ,...){


   #---------------------------------------------------------------------------------------#
   #     Define the error type.                                                            #
   #---------------------------------------------------------------------------------------#
   err.type = match.arg(err.type)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Make sure arguments for density have the same length as x, y, xerr, yerr.         #
   #---------------------------------------------------------------------------------------#
   nxy = length(x)
   if (length(err.col ) == 1) err.col  = rep(err.col ,times=nxy)
   if (length(edensity) == 1) edensity = rep(edensity,times=nxy)
   if (length(eangle  ) == 1) eangle   = rep(eangle  ,times=nxy)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Make sure that errors settings are not inconsistent.                             #
   #---------------------------------------------------------------------------------------#
   error  = FALSE
   skip.x = FALSE
   skip.y = FALSE
   #----- Check errors along 'x' axis. ----------------------------------------------------#
   if (is.null(xerr)){
      #------------------------------------------------------------------------------------#
      #    When 'xerr' is NULL xlow and xhigh must be consistent.  Either they must be     #
      # both NULL (so no error is plotted in the x axis), or they must be both defined.    #
      #------------------------------------------------------------------------------------#
      if ( is.null(xlow) && is.null(xhigh) ){
         #----- No error in the x axis. ---------------------------------------------------#
         skip.x = TRUE
         #---------------------------------------------------------------------------------#
      }else if (is.null(xlow) != is.null(xhigh)){
         #----- 'xerr' is NULL so xlow and xhigh must be both NULL or both defined. -------#
         cat(" You cannot define only one of 'xlow' and 'xhigh'...","\n")
         error = TRUE 
         #---------------------------------------------------------------------------------#
      }#end if (is.null(xlow) && is.null(xhigh))
      #------------------------------------------------------------------------------------#
   }else{
      #------------------------------------------------------------------------------------#
      #    When 'xerr' is defined, then both xlow and xhigh must be NULL.  They will be    #
      # defined here so for the remaining of this function xerr is not called.             #
      #------------------------------------------------------------------------------------#
      if (is.null(xlow) && is.null(xhigh)){
         #----- Define xlow and xhigh based on xerr. --------------------------------------#
         xlow  = x - xerr
         xhigh = x + xerr
         #---------------------------------------------------------------------------------#
      }else{
         #----- One of xlow/xhigh was provided, which shouldn't happen. -------------------#
         cat(" If 'xerr' is provided, then both 'xlow' and 'xhigh' must be NULL...","\n")
         error = TRUE
         #---------------------------------------------------------------------------------#
      }#end if (! (is.null(xlow) && is.null(xhigh)) )
      #------------------------------------------------------------------------------------#
   }#end if (is.null(xerr))
   #----- Check errors along 'y' axis. ----------------------------------------------------#
   if (is.null(yerr)){
      #------------------------------------------------------------------------------------#
      #    When 'yerr' is NULL ylow and yhigh must be consistent.  Either they must be     #
      # both NULL (so no error is plotted in the x axis), or they must be both defined.    #
      #------------------------------------------------------------------------------------#
      if ( is.null(ylow) && is.null(yhigh) ){
         #----- No error in the x axis. ---------------------------------------------------#
         skip.y = TRUE
         #---------------------------------------------------------------------------------#
      }else if (is.null(ylow) != is.null(yhigh)){
         #----- 'yerr' is NULL so ylow and yhigh must be both NULL or both defined. -------#
         cat(" You cannot define only one of 'ylow' and 'yhigh'...","\n")
         error = TRUE 
         #---------------------------------------------------------------------------------#
      }#end if (is.null(ylow) && is.null(yhigh))
      #------------------------------------------------------------------------------------#
   }else{
      #------------------------------------------------------------------------------------#
      #    When 'yerr' is defined, then both ylow and yhigh must be NULL.  They will be    #
      # defined here so for the remaining of this function yerr is not called.             #
      #------------------------------------------------------------------------------------#
      if (is.null(ylow) && is.null(yhigh)){
         #----- Define ylow and yhigh based on yerr. --------------------------------------#
         ylow  = y - yerr
         yhigh = y + yerr
         #---------------------------------------------------------------------------------#
      }else{
         #----- One of ylow/yhigh was provided, which shouldn't happen. -------------------#
         cat(" If 'yerr' is provided, then both 'ylow' and 'yhigh' must be NULL...","\n")
         error = TRUE
         #---------------------------------------------------------------------------------#
      }#end if (! (is.null(ylow) && is.null(yhigh)) )
      #------------------------------------------------------------------------------------#
   }#end if (is.null(yerr))
   if (error) stop("Inconsistent error settings...")
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Decide whether to start a new plot window or put over an existing plot.           #
   #---------------------------------------------------------------------------------------#
   if (! add){
      #----- Find range for x and y in case none has been provided. -----------------------#
      if (is.null(xlim)) xlim = range(c(x,xlow,xhigh),na.rm=TRUE)
      if (is.null(ylim)) ylim = range(c(y,ylow,yhigh),na.rm=TRUE)
      #------------------------------------------------------------------------------------#


      #----- New plot ---------------------------------------------------------------------#
      plot(x=x,y=y,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,type="n",...)
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Define how to plot the error type.                                                 #
   #---------------------------------------------------------------------------------------#
   if (err.type %in% "l"){
       #-----------------------------------------------------------------------------------#
       #     Plot the X error bars.                                                        #
       #-----------------------------------------------------------------------------------#
       if (! skip.x){
          #----- Scale the cap with the total scale. --------------------------------------#
          ycoord = par()$usr[3:4]
          smidge = 0.5 * cap * (ycoord[2] - ycoord[1])
          #--------------------------------------------------------------------------------#



          #----- Decide the scale depending on whether the y scale is linear or log. ------#
          if (par()$ylog){
             ycapa  = y * 10^(-smidge)
             ycapz  = y * 10^( smidge)
          }else{
             ycapa  = y - smidge
             ycapz  = y + smidge
          }#end if
          #--------------------------------------------------------------------------------#


          #---- Plot the main stem. -------------------------------------------------------#
          segments(xlow, y, xhigh, y, lty = lty, lwd = lwd,col=err.col)
          #--------------------------------------------------------------------------------#



          #---- Plot the caps. ------------------------------------------------------------#
          segments(xlow ,ycapa,xlow ,ycapz,lty=lty,lwd=lwd,col=err.col)
          segments(xhigh,ycapa,xhigh,ycapz,lty=lty,lwd=lwd,col=err.col)
          #--------------------------------------------------------------------------------#
       }#end if (! skip.x)
       #-----------------------------------------------------------------------------------#



       #-----------------------------------------------------------------------------------#
       #     Plot the Y error bars.                                                        #
       #-----------------------------------------------------------------------------------#
       if (! skip.y){
          #----- Scale the cap with the total scale. --------------------------------------#
          xcoord = par()$usr[1:2]
          smidge = 0.5 * cap * (xcoord[2] - xcoord[1])
          #--------------------------------------------------------------------------------#



          #----- Decide the scale depending on whether the x scale is linear or log. ------#
          if (par()$xlog) {
             xcapa = x * 10^(-smidge)
             xcapz = x * 10^( smidge)
          }else{
             xcapa = x - smidge
             xcapz = x + smidge
          }#end if
          #--------------------------------------------------------------------------------#


          #---- Plot the main stem. -------------------------------------------------------#
          segments(x,ylow,x,yhigh,lty = lty,lwd=lwd,col=err.col)
          #--------------------------------------------------------------------------------#



          #---- Plot the caps. ------------------------------------------------------------#
          segments(xcapa,ylow ,xcapz,ylow ,lwd=lwd,lty=lty,col=err.col)
          segments(xcapa,yhigh,xcapz,yhigh,lwd=lwd,lty=lty,col=err.col)
          #--------------------------------------------------------------------------------#
       }#end if (! skip.y)
       #-----------------------------------------------------------------------------------#
   }else if (err.type %in% "e"){
       #------ Make shaded ellipses. ------------------------------------------------------#
       mapply( FUN      = err.ellipse
             , x        = x
             , y        = y
             , xlow     = xlow
             , xhigh    = xhigh
             , ylow     = ylow
             , yhigh    = yhigh
             , col      = err.col
             , border   = err.col
             , density  = edensity
             , angle    = eangle
             )#end mapply
       #-----------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#


   #---- Points go after everything else. -------------------------------------------------#
   points(x, y, pch = pch, type = type, col=col,lwd=lwd,lty=lty,...)
   #---------------------------------------------------------------------------------------#

   return(invisible())
}#end function
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     This routine plots an elliptical feature to show the error area.  It interpolates    #
# the radius using spline, so it allows different errors on the x and y axis, as well as   #
# differences in the lower and upper bounds of errors.                                     #
#------------------------------------------------------------------------------------------#
err.ellipse <<- function( x
                        , y
                        , xlow
                        , xhigh
                        , ylow
                        , yhigh
                        , nv      = 360
                        , lty     = 1
                        , lwd     = 1
                        , ...
                        ){


   #------ Define radius dimensions and angles. -------------------------------------------#
   a.plus  = xhigh - x
   a.minus = x  - xlow
   b.plus  = yhigh - y
   b.minus = y  - ylow
   theta   = seq(from=0,to=2*pi,length.out=nv+1)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Interpolate radius.                                                              #
   #---------------------------------------------------------------------------------------#
   aa = ifelse( test = theta < pi/2 | theta >= 3*pi/2, yes = a.plus, no = a.minus)
   bb = ifelse( test = theta > 0.   & theta <= pi    , yes = b.plus, no = b.minus)
   xx = aa * cos(theta)
   yy = bb * sin(theta)
   ee = atan2(yy,xx)
   rr = sqrt(xx^2 + yy^2)
   xp = x + rr * cos(ee)
   yp = y + rr * sin(ee)
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #     Plot polygon.                                                                     #
   #---------------------------------------------------------------------------------------#
   polygon(x=xp,y=yp,...)
   #---------------------------------------------------------------------------------------#


   #----- Nothing is returned. ------------------------------------------------------------#
   invisible(NULL)
   #---------------------------------------------------------------------------------------#
}#end err.ellipse
#==========================================================================================#
#==========================================================================================#
