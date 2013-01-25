#==========================================================================================#
#==========================================================================================#
# Function colourmap                                                                       #
#                                                                                          #
#    Given the x and y coordinates, this function will plot the z values with a colour     #
# scheme but no interpolation...                                                           #
#------------------------------------------------------------------------------------------#
colourmap  = function( x
                     , y
                     , z
                     , xlim             = range(x,finite=TRUE)
                     , ylim             = range(y,finite=TRUE)
                     , zlim             = range(z,finite=TRUE)
                     , levels           = if (key.log){
                                             pretty.log(x=zlim,n=nlevels)
                                          }else{
                                            pretty(x=zlim,n=nlevels)
                                          }#end if
                     , nlevels          = 20
                     , colour.palette   = cm.colors
                     , col              = colour.palette(length(levels)-1)
                     , na.col           = "grey94"
                     , plot.title
                     , plot.axes
                     , key.title
                     , key.axes
                     , key.log          = FALSE
                     , axes             = TRUE
                     , frame.plot       = axes
                     , pch              = 15
                     , cex              = 1.0
                     , ...
                     ){

   #---------------------------------------------------------------------------------------#
   #      All three coordinates must be given.                                             #
   #---------------------------------------------------------------------------------------#
   if (missing(x) || missing(y) || missing(z)){
      print(paste(" X is missing: ",missing(x)))
      print(paste(" Y is missing: ",missing(y)))
      print(paste(" Z is missing: ",missing(z)))
      stop ("At least one of the data points is missing...")
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Check whether x, y, and z are the same type of data.                             #
   #---------------------------------------------------------------------------------------#
   same.kind = (is.list(x) == is.list(y) && is.list(x) == is.list(y))
   if (! same.kind){
      print(paste(" X is list: ",is.list(x)))
      print(paste(" Y is list: ",is.list(y)))
      print(paste(" Z is list: ",is.list(z)))
      stop ("X, Y, and Z must be of the same kind...")
   }else if (!is.list(x)){
      #----- Convert x, y, and z to lists. ------------------------------------------------#
      x       = list(x)
      y       = list(y)
      z       = list(z)
      pch     = list(pch)
      cex     = list(cex)
      npanels = 1
   }else{
      npanels = length(x)
      if (! is.list(pch)){
         orig.pch = pch
         pch      = list()
         for (p in 1:npanels) pch[[p]] = orig.pch
      }#end if
      if (! is.list(cex)){
         orig.cex = cex
         cex      = list()
         for (p in 1:npanels) cex[[p]] = orig.cex
      }#end if
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Save the margins to avoid losing the data. --------------------------------------#
   mar.orig = (par.orig = par(c("mar", "las", "mfrow")))$mar
   on.exit(par(par.orig))
   #---------------------------------------------------------------------------------------#



   #----- Split the screen into 3, the two panels and the scale. --------------------------#
   if (npanels == 1){
      w = (3 + mar.orig[2]) * par("csi") * 2.54
      layout(matrix(c(2, 1), nc = 2), widths = c(1, lcm(w)))
      mar = mar.orig
      mar[4] = mar[2]
      mar[2] = 1
      key.vertical = TRUE
   }else{
      h = (1 + mar.orig[3]) * par("csi") * 2.54
      layout( mat     =rbind(seq(from=2,to=npanels+1),rep(1,times=npanels))
            , heights = c(1, lcm(h))
            )#end layout
      mar = mar.orig
      mar[1] = 2.6
      mar[3] = 2.1
      key.vertical = FALSE
   }#end if
   #---------------------------------------------------------------------------------------#
   #=======================================================================================#
   #=======================================================================================#





   #=======================================================================================#
   #=======================================================================================#
   #      First plot: the key scale.                                                       #
   #---------------------------------------------------------------------------------------#
      par(mar = mar)
      plot.new()
      #------------------------------------------------------------------------------------#
      #     Plot in the horizontal or vertical depending on where the scale is going to    #
      # be plotted.                                                                        #
      #------------------------------------------------------------------------------------#
      if (key.vertical){
         #----- Decide whether the scale is logarithmic or not. ---------------------------#
         if (key.log){
            plot.window(xlim=c(0,1),ylim=range(levels),xaxs="i",yaxs="i",log="y")
         }else{
            plot.window(xlim=c(0,1),ylim=range(levels),xaxs="i",yaxs="i")
         }#end if
         #---------------------------------------------------------------------------------#

         #----- Draw the colour bar. ------------------------------------------------------#
         rect(xleft=0,ybottom=levels[-length(levels)],xright=1,ytop=levels[-1],col=col
             ,border=col)
         #---------------------------------------------------------------------------------#

         #----- Check whether there are specific instructions for plotting the key axis. --#
         if (missing(key.axes)) {
            if (axes) axis(side=4,las=1,...)
         }else{
            key.axes
         }#end if
         #---------------------------------------------------------------------------------#
      }else{
         #----- Decide whether the scale is logarithmic or not. ---------------------------#
         if (key.log){
            plot.window(xlim=range(levels),ylim=c(0,1),xaxs="i",yaxs="i",las=1,log="y")
         }else{
            plot.window(xlim=range(levels),ylim=c(0,1),xaxs="i",yaxs="i",las=1)
         }#end if
         #---------------------------------------------------------------------------------#


         #----- Draw the colour bar. ------------------------------------------------------#
         rect(xleft=levels[-length(levels)],ybottom=0,xright=levels[-1],ytop=1
             ,col=col,border=col)
         #---------------------------------------------------------------------------------#


         #----- Check whether there are specific instructions for plotting the key axis. --#
         if (missing(key.axes)) {
            if (axes) axis(side=1,...)
         }else{
            key.axes
         }#end if
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#


      #----- Draw box. --------------------------------------------------------------------#
      box()
      #------------------------------------------------------------------------------------#


      #----- Plot the title. --------------------------------------------------------------#
      if (!missing(key.title)) key.title
      #------------------------------------------------------------------------------------#
   #=======================================================================================#
   #=======================================================================================#





   #=======================================================================================#
   #=======================================================================================#
   #      Now we plot the other panels.                                                    #
   #---------------------------------------------------------------------------------------#
   for (p in 1:npanels){
      #----- Set the window. --------------------------------------------------------------#
      mar    = mar.orig
      if (! key.vertical) mar[1] = 4.1
      par(mar = mar)
      plot.new()
      plot.window(xlim=xlim,ylim=ylim,...)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    Split zleft into the breaks defined by the colour palette.                      #
      #------------------------------------------------------------------------------------#
      zcut              = cut(z[[p]],breaks=levels)
      zlev              = levels(zcut)
      zcol              = col[match(zcut,zlev)]
      zcol[is.na(zcol)] = na.col
      #------------------------------------------------------------------------------------#



      #----- Call the function that actually plots the data. ------------------------------#
      points(x=x[[p]],y=y[[p]],pch=pch[[p]],cex=cex[[p]],col=zcol,...)
      #------------------------------------------------------------------------------------#



      #----- Check whether there are especial instructions for plotting the axes. ---------#
      if (missing(plot.axes)) {
          if (axes) {
              axis(1)
              axis(2)
          }
      }else{
         if (is.list(plot.axes)){
            if (! is.null(plot.axes[[p]]$x.axis)) do.call("axis",plot.axes[[p]]$x.axis)
            if (! is.null(plot.axes[[p]]$y.axis)) do.call("axis",plot.axes[[p]]$y.axis)
            other = which( ! names(plot.axes[[p]]) %in% c("x.axis","y.axis"))
            if (length(other) > 0){
               for (o in other){
                  do.call(names(plot.axes[[p]])[o],plot.axes[[p]][[o]])
               }#end for
            }#end if
         }else{
            plot.axes
         }#end if
      }#end if
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    Plot the frame and the tiles.                                                   #
      #------------------------------------------------------------------------------------#
      if (frame.plot) box()
      #----- Check whether there are especial instructions for plotting the title. --------#
      if (missing(plot.title)){
          if (axes) title(main = "", xlab = "", ylab = "",...)
      }else if (is.list(plot.title)){
          do.call("title",plot.title[[p]])
      }else{
          plot.title
      }#end if
      #------------------------------------------------------------------------------------#
   }#end for
   #=======================================================================================#
   #=======================================================================================#

   invisible()
}#end function colourmap
#==========================================================================================#
#==========================================================================================#
