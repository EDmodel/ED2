#==========================================================================================#
#==========================================================================================#
# Function colourmap                                                                       #
#                                                                                          #
#    Given the x and y coordinates, this function will plot the vegetation values (modis)  #
# with a colour scheme but no interpolation...                                             #
#------------------------------------------------------------------------------------------#
mcd12q1.map = function( x         = seq(from=0,to=1,len=nrow(z))
                      , y         = seq(from=0,to=1,len=ncol(z))
                      , z
                      , xlim      = range(x,finite=TRUE)
                      , ylim      = range(y,finite=TRUE)
                      , zlim      = range(v,finite=TRUE)
                      , highlight = matrix(data=rep(FALSE,times=length(z))
                                          ,nrow=nrow(z),ncol=ncol(z))
                      , na.col="gray94",plot.title,plot.axes
                      , col.highlight="red3"
                      , asp=NA,xaxs="i",yaxs="i",las=1,axes=TRUE,frame.plot=axes,pch=15
                      , cex=1.0,...){


   modis.col = c("midnightblue","steelblue3","#004E00","deepskyblue","chartreuse4"
                ,"orange1","lawngreen","lightgoldenrod4","darkseagreen3","gray30"
                ,"lightskyblue3","saddlebrown",col.highlight,"gray80")
   modis.leg = c("H2O","ENT","EBT","DNT","DBT","SHR","GRS","CER","BLC","URB","ICE","DES"
                ,"USE","UND")
   n.modis   = length(modis.leg)
   modis.val = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13)
   modis.bks = c(modis.val[1]-0.5,modis.val+0.5)
   modis.at  = sequence(n.modis)
   z[! z %in% modis.val] = 13
   z[highlight]          = 12


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
            x = seq(0, 1, len = length(z))
         }#end if
       }else{
          #----- Bad setting. -------------------------------------------------------------#
          stop("no `z' variable specified")
       }#end if
   }else if (is.list(x)) {
       #----- Z is there, just need to check whether x and y were given as a list... ------#
       y = x$y
       x = x$x
   }#end if

   #----- Save the margins to avoid losing the data. --------------------------------------#
   mar.orig = (par.orig = par(c("mar", "las", "mfrow")))$mar
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
   plot.window(xlim=c(0,1),ylim=c(0.5,n.modis+0.5),xaxs="i",yaxs="i")

   #----- Draw the colour bar. ------------------------------------------------------------#
   rect(xleft=0,ybottom=seq(0.5,n.modis-0.5,1),xright=1,ytop=seq(1.5,n.modis+0.5,1)
       ,col=modis.col,border="black")

   #----- Check whether there are specific instructions for plotting the key axis. --------#
   if (axes){
      axis(side=4,at=modis.at,labels=modis.leg)
   }#end if

   #----- Draw box. -----------------------------------------------------------------------#
   box()

   #----- Now we plot the filled contour. -------------------------------------------------#
   mar    = mar.orig
   mar[4] = 1
   par(mar = mar)
   plot.new()
   plot.window(xlim=xlim,ylim=ylim,log="", xaxs = xaxs, yaxs = yaxs, asp = asp)

   #---------------------------------------------------------------------------------------#
   #    Split z into the breaks defined by the colour palette.                             #
   #---------------------------------------------------------------------------------------#
   zcol = modis.col[match(z,modis.val)]
   zcol[is.na(zcol)] = na.col

   #----- Call the function that actually plots the data. ---------------------------------#
   image(x=x,y=y,z=z,add=TRUE,col=modis.col,breaks=modis.bks,...)

   #----- Check whether there are especial instructions for plotting the axes. ------------#
   if (missing(plot.axes)) {
       if (axes) {
           axis(1)
           axis(2)
       }
   }else{
      plot.axes
   }#end if
   
   if (frame.plot) box()
   #----- Check whether there are especial instructions for plotting the title. -----------#
   if (missing(plot.title)){
       if (axes) title(main = "", xlab = "", ylab = "",...)
   }else{
       plot.title
   }#end if

   invisible()
}#end function simbrazil
#==========================================================================================#
#==========================================================================================#
