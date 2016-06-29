#==========================================================================================#
#==========================================================================================#
# Function sombreado                                                                       #
#                                                                                          #
#    This function is exactly the same as filled.contour, but without the black lines      #
# between the colours in the scale...                                                      #
#------------------------------------------------------------------------------------------#
sombreado <<- function(x=seq(from=0,to=1,len=nrow(z)),y=seq(from=0,to=1,len=ncol(z)),z
                      ,xlim=range(x,finite=TRUE),ylim=range(y,finite=TRUE)
                      ,zlim=range(z,finite=TRUE)
                      ,levels = if (key.log){
                                   sort(unique(pretty.log(x=z,n=nlevels,forcelog=TRUE)))
                                }else{
                                   sort(unique(pretty(x=z,n=nlevels)))
                                }#end if
                      ,nlevels=100,colour.palette=color.palette,color.palette=cm.colors
                      ,col=colour.palette(length(levels)-1)
                      ,plot.title,plot.axes,xlog=FALSE,ylog=FALSE,key.title,key.axes
                      ,key.log=FALSE,asp=NA,interp=TRUE
                      ,xaxs="i",yaxs="i",las=1,axes=TRUE,frame.plot=axes,useRaster=TRUE
                      ,...){

   #----- Check which kind of input was given. --------------------------------------------#
   if (missing(z)) {
      #----- No z was given x must be a list or the user didn't provide any axis... -------# 
      if (!missing(x)) {
         if (is.list(x)) {
            #----- X is a list, copy the elements to variables. ---------------------------#
            z = x$z
            y = x$y
            x = x$x
         }else{
            #----- x is an array, make up some x axis. ------------------------------------#
            z = x
            x = seq(from = 0, to = 1, length.out = nrow(z))
            y = seq(from = 0, to = 1, length.out = ncol(z))
         }#end if
       }else{
          #----- Bad setting. -------------------------------------------------------------#
          stop("no `z' matrix specified")
       }#end if
   }else if (is.list(x)) {
       #----- Z is there, just need to check whether x and y were given as a list... ------#
       y = x$y
       x = x$x
   }#end if

   #----- Check whether the z matrix makes sense or not. ----------------------------------#
   if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1){
      stop("no proper `z' matrix specified")
   }else if(!is.double(z)){
      storage.mode(z) = "double"
   }#end if

   #----- No messed-up axes are allowed, they must increase. ------------------------------#
   if (any(diff(x) <= 0) || any(diff(y) <= 0)){
       stop("increasing x and y values expected")
   }#end if

   #----- Save the margins to avoid losing the data. --------------------------------------#
   par.orig = par(c("mar", "las", "mfrow"))
   mar.orig = par.orig$mar
   on.exit(par(par.orig))

   #----- Split the screen into 2. --------------------------------------------------------#
   w = (3 + mar.orig[2]) * par("csi") * 2.54
   layout(matrix(c(2, 1), nc = 2), widths = c(1, lcm(w)))
   par(las = las)
   mar = mar.orig
   mar[4] = mar[2]
   mar[2] = 1

   #----- First plot: the key scale. ------------------------------------------------------#
   par(mar = mar)
   plot.new()
   #----- Decide whether the scale is logarithmic or not. ---------------------------------#
   if (key.log){
      plot.window(xlim=c(0,1),ylim=range(levels),xaxs="i",yaxs="i",log="y")
   }else{
      plot.window(xlim=c(0,1),ylim=range(levels),xaxs="i",yaxs="i")
   }#end if

   #----- Draw the colour bar. ------------------------------------------------------------#
   rect(xleft=0,ybottom=levels[-length(levels)],xright=1,ytop=levels[-1],col=col
       ,border=col)

   #----- Check whether there are specific instructions for plotting the key axis. --------#
   if (missing(key.axes)) {
      if (axes) axis(side=4,...)
   }else{
      key.axes
   }#end if

   #----- Draw box. -----------------------------------------------------------------------#
   box()

   #----- Plot the title. -----------------------------------------------------------------#
   if (!missing(key.title)) key.title

   #----- Make the log variable for the main window. --------------------------------------#
   if (xlog) xlim = log(xlim)
   if (ylog) ylim = log(ylim)

   #----- Now we plot the filled contour. -------------------------------------------------#
   mar    = mar.orig
   mar[4] = 1
   par(mar = mar)
   plot.new()
   plot.window(xlim=xlim,ylim=ylim, xaxs = xaxs, yaxs = yaxs, asp = asp)



   #---------------------------------------------------------------------------------------#
   #      We use image to plot, so it looks nice in PDF.                                   #
   #---------------------------------------------------------------------------------------#
   #----- Make x and y dimensionless. -----------------------------------------------------#
   if (xlog){
      xx    = log(as.numeric(x))
   }else{
      xx    = as.numeric(x)
   }#end if
   if (ylog){
      yy    = log(as.numeric(y))
   }else{
      yy    = as.numeric(y)
   }#end if
   nx    = length(xx)
   ny    = length(yy)
   xlow  = min(xx)
   xhigh = max(xx)
   ylow  = min(yy)
   yhigh = max(yy)
   #----- Scale x and y. ------------------------------------------------------------------#
   xxx    = rep((xx-xlow)/(xhigh-xlow),times=length(yy))
   yyy    = rep((yy-ylow)/(yhigh-ylow),each =length(xx))
   sss    = is.finite(z)
   xo     = seq(from=0,to=1,length.out=10*length(xx))
   yo     = seq(from=0,to=1,length.out=10*length(yy))
   zint   = interp(x=xxx[sss],y=yyy[sss],z=z[sss],xo=xo,yo=yo)
   sint   = interp(x=xxx     ,y=yyy     ,z=sss   ,xo=xo,yo=yo)
   sint$z = sint$z > twothirds
   zint$z = ifelse(sint$z,zint$z,NA)
   zint$x = xlow + zint$x * (xhigh - xlow)
   zint$y = ylow + zint$y * (yhigh - ylow)
   image(zint,breaks=levels,col=col,add=TRUE,useRaster=useRaster)
   #---------------------------------------------------------------------------------------#



   #----- Check whether there are especial instructions for plotting the axes. ------------#
   if (missing(plot.axes)) {
       if (axes) {
           if (xlog){
              xlab = pretty.log(x)
              xat  = log(xlab)
              axis(side=1,at=xat,labels=xlab)
           }else{
              axis(side=1)
           }#end if
           if (ylog){
              ylab = pretty.log(y)
              yat  = log(ylab)
              axis(side=2,las=1,at=yat,labels=ylab)
           }else{
              axis(side=2,las=1)
           }#end if
       }#end if
   }else{
      plot.axes
   }#end if
   #---------------------------------------------------------------------------------------#


   if (frame.plot) box()
   #----- Check whether there are especial instructions for plotting the title. -----------#
   if (missing(plot.title)){
       if (axes) title(main = "", xlab = "", ylab = "",...)
   }else{
       plot.title
   }#end if

   invisible()
}#end function sombreado
#==========================================================================================#
#==========================================================================================#
