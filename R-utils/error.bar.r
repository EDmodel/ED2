#------------------------------------------------------------------------------------------#
#     This function plots a scatter plot adding error bars.                                #
#------------------------------------------------------------------------------------------#
error.bar <<- function(x, y,xerr=NULL,yerr=NULL,xlow=NULL,xhigh=NULL,ylow=NULL,yhigh=NULL
                      ,xlim=NULL,ylim=NULL,cap=0.015,main=NULL,sub=NULL
                      ,xlab = as.character(substitute(x))
                      ,ylab = if (is.factor(x) || is.character(x)){
                                 ""
                              }else{
                                 as.character(substitute(y))
                              }#end if
                      ,add = FALSE, lty = 1, type = "p", lwd = 1
                      ,pch = 16,col="black",err.col=col
                      ,...){

    #--------------------------------------------------------------------------------------#
    #     Decide whether to start a new plot window or put over an existing plot.          #
    #--------------------------------------------------------------------------------------#
    if (add){
       #---- Points. ----------------------------------------------------------------------#
       points(x, y, pch = pch, type = type, col=col,lwd=lwd,lty=lty,...)
       #-----------------------------------------------------------------------------------#
    }else{
       #----- Find range for x and y in case none has been provided. ----------------------#
       if (is.null(xlim)){
          if ((! is.null(xlow)) & (! is.null(xhigh))){
             xlim = range(c(x,xlow,xhigh),na.rm=TRUE)
           }else{
             xlim = range(c(x,x+xerr,x-xerr),na.rm=TRUE)
          }#end if
       }#end if
       if (is.null(ylim)){
          if ((! is.null(ylow)) & (! is.null(yhigh))){
             ylim = range(c(y,ylow,yhigh),na.rm=TRUE)
          }else{
             ylim = range(c(y,y+yerr,y-yerr),na.rm=TRUE)
          }#end if
       }#end if
       #-----------------------------------------------------------------------------------#


       #----- New plot --------------------------------------------------------------------#
       plot(x=x,y=y,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,type=type
           ,col=col,lwd=lwd,lty=lty,...)
       #-----------------------------------------------------------------------------------#
    }#end if
    #--------------------------------------------------------------------------------------#



    #--------------------------------------------------------------------------------------#
    #     Decide the bounds for X and Y when low and high aren't provided.                 #
    #--------------------------------------------------------------------------------------#
    if ( (! is.null(xerr)) && (is.null(xlow) | is.null(xhigh))){
       xlow  = x - xerr
       xhigh = x + xerr
    }else if (is.null(xerr) && (! is.null(xlow)) && (! is.null(xhigh))){
       xerr  = TRUE
    }#end if
    if ( (! is.null(yerr)) && (is.null(ylow) | is.null(yhigh))){
       ylow  = y - yerr
       yhigh = y + yerr
    }else if (is.null(yerr) && (! is.null(ylow)) && (! is.null(yhigh))){
       yerr  = TRUE
    }#end if
    #--------------------------------------------------------------------------------------#



    #--------------------------------------------------------------------------------------#
    #     Plot the X error bars.                                                           #
    #--------------------------------------------------------------------------------------#
    if (! is.null(xerr)){
       #----- Scale the cap with the total scale. -----------------------------------------#
       ycoord = par()$usr[3:4]
       smidge = 0.5 * cap * (ycoord[2] - ycoord[1])
       #-----------------------------------------------------------------------------------#



       #----- Decide the scale depending on whether the y scale is linear or log. ---------#
       if (par()$ylog){
          ycapa  = y * 10^(-smidge)
          ycapz  = y * 10^( smidge)
       }else{
          ycapa  = y - smidge
          ycapz  = y + smidge
       }#end if
       #-----------------------------------------------------------------------------------#


       #---- Plot the main stem. ----------------------------------------------------------#
       segments(xlow, y, xhigh, y, lty = lty, lwd = lwd,col=err.col)
       #-----------------------------------------------------------------------------------#



       #---- Plot the caps. ---------------------------------------------------------------#
       segments(xlow ,ycapa,xlow ,ycapa,lty=lty,lwd=lwd,col=err.col)
       segments(xhigh,ycapz,xhigh,ycapz,lty=lty,lwd=lwd,col=err.col)
       #-----------------------------------------------------------------------------------#
    }#end if
    #--------------------------------------------------------------------------------------#



    #--------------------------------------------------------------------------------------#
    #     Plot the Y error bars.                                                           #
    #--------------------------------------------------------------------------------------#
    if (! is.null(yerr)){
       #----- Scale the cap with the total scale. -----------------------------------------#
       xcoord = par()$usr[1:2]
       smidge = cap * (xcoord[2] - xcoord[1])/2
       #-----------------------------------------------------------------------------------#



       #----- Decide the scale depending on whether the x scale is linear or log. ---------#
       if (par()$xlog) {
          xcapa = x * 10^(-smidge)
          xcapz = x * 10^( smidge)
       }else{
          xcapa = x - smidge
          xcapz = x + smidge
       }#end if
       #-----------------------------------------------------------------------------------#


       #---- Plot the main stem. ----------------------------------------------------------#
       segments(x,ylow,x,yhigh,lty = lty,lwd=lwd,col=err.col)
       #-----------------------------------------------------------------------------------#



       #---- Plot the caps. ---------------------------------------------------------------#
       segments(xcapa,ylow ,xcapz,ylow ,lwd=lwd,lty=lty,col=err.col)
       segments(xcapa,yhigh,xcapz,yhigh,lwd=lwd,lty=lty,col=err.col)
       #-----------------------------------------------------------------------------------#
    }#end if
    #--------------------------------------------------------------------------------------#


    return(invisible())
}#end function
#------------------------------------------------------------------------------------------#
