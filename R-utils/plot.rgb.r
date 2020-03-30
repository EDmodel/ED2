#==========================================================================================#
#==========================================================================================#
#     This function plots a graph as a function of 3 parameters that add to a constant,    #
# using a RGB scale.                                                                       #
#------------------------------------------------------------------------------------------#
plot.rgb <<- function( x
                     , y
                     , r
                     , g
                     , b
                     , dx               = ifelse( is.list(x)
                                                , lapply(lapply(x,diff),median,na.rm=TRUE)
                                                , median(diff(x),na.rm=TRUE)
                                                )#end ifelse
                     , dy               = ifelse( is.list(y)
                                                , lapply(lapply(y,diff),median,na.rm=TRUE)
                                                , median(diff(y),na.rm=TRUE)
                                                )#end ifelse
                     , xlim             = range(unlist(x),finite=TRUE)
                     , ylim             = range(unlist(y),finite=TRUE)
                     , zlim             = range(unlist(z),finite=TRUE)
                     , rgb.delta        = 0.10
                     , rgb.at           = c(0.1,0.9)
                     , rgb.axis.labels  = list(l="Red",b="Green",r="Blue")
                     , rgb.show.axis    = TRUE
                     , rgb.shift.labels = ! rgb.show.axis
                     , rgbsum.fixed     = TRUE
                     , main.vlab        = c(0,1)
                     , v.delta          = 0.04
                     , v.axis.options   = NULL
                     , x.axis.options   = NULL
                     , y.axis.options   = NULL
                     , sub.options      = NULL
                     , main.title       = NULL
                     , main.xlab        = NULL
                     , main.ylab        = NULL
                     , rgb.title        = NULL
                     , v.title          = NULL
                     , plot.after       = NULL
                     , key.right        = TRUE
                     , f.key            = 1/5
                     , na.col           = "grey94"
                     , ...
                     ){

   #---------------------------------------------------------------------------------------#
   #     Find out whether x, y, r, g, and b are single values or lists.                    #
   #---------------------------------------------------------------------------------------#
   if (missing(x) | missing(y) | missing(r) | missing(g) | missing (b)){
      cat(" - x is missing: ",missing(x),"\n")
      cat(" - y is missing: ",missing(y),"\n")
      cat(" - r is missing: ",missing(r),"\n")
      cat(" - g is missing: ",missing(g),"\n")
      cat(" - b is missing: ",missing(b),"\n")
      stop(" x, y, r, g, and b must be given")
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Check whether x, y, and z are the same type of data.                             #
   #---------------------------------------------------------------------------------------#
   same.kind = (is.list(x) == is.list(y) && is.list(x) == is.list(y) &&
                is.list(x) == is.list(r) && is.list(x) == is.list(g) &&
                is.list(x) == is.list(b))
   if (! same.kind){
      cat(" X is list: ",is.list(x),"\n")
      cat(" Y is list: ",is.list(y),"\n")
      cat(" R is list: ",is.list(r),"\n")
      cat(" G is list: ",is.list(g),"\n")
      cat(" B is list: ",is.list(b),"\n")
      stop ("X, Y, R, G, and B must be of the same kind...")
   }else if (!is.list(x)){
      #----- Convert x, y, and z to lists. ------------------------------------------------#
      x              = list(x )
      y              = list(y )
      r              = list(r )
      g              = list(g )
      b              = list(b )
      dx             = list(dx)
      dy             = list(dy)
      if (! missing(x.axis.options)) x.axis.options = list(x.axis.options)
      if (! missing(y.axis.options)) y.axis.options = list(y.axis.options)
      if (! missing(sub.options   )) sub.options    = list(sub.options   )
      npanels = 1
   }else{
      npanels = length(x)
   }#end if
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     x and y should be axes for r, g, and b.  Check whether they are or not.           #
   #---------------------------------------------------------------------------------------#
   nx = sapply(X = x, FUN = length)
   ny = sapply(X = y, FUN = length)
   nr = sapply(X = r, FUN = length)
   ng = sapply(X = g, FUN = length)
   nb = sapply(X = b, FUN = length)
   if ( any(nx != ny) ||  any(nx != nr) ||  any(nx != ng) ||  any(nx != nb)){
      cat(" - length(x): ",paste(nx,sep=" "),"\n")
      cat(" - length(y): ",paste(ny,sep=" "),"\n")
      cat(" - length(r): ",paste(nr,sep=" "),"\n")
      cat(" - length(g): ",paste(ng,sep=" "),"\n")
      cat(" - length(b): ",paste(nb,sep=" "),"\n")
      stop(" x, y, r, g, and b must have all the same length")
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Find the colour scaling factor. -------------------------------------------------#
   if (rgbsum.fixed){
      rgb.max = max(unlist(r)+unlist(g)+unlist(b),na.rm=TRUE)
   }else{
      rgb.max = max(c(unlist(r),unlist(g),unlist(b)),na.rm=TRUE)
   }#end if (rgbsum.fixed)
   #---------------------------------------------------------------------------------------#



   #----- Correct colours. ----------------------------------------------------------------#
   pcol    = mapply( FUN      = rgb.norm
                   , red      = r
                   , green    = g
                   , blue     = b
                   , MoreArgs = list(maxColorValue=rgb.max,na.col=na.col)
                   , SIMPLIFY = FALSE
                   )#end mapply
   #---------------------------------------------------------------------------------------#


   #----- Save the margins to avoid losing the data. --------------------------------------#
   par.orig     = par(no.readonly=TRUE)
   par.here     = par(las=1,xpd=FALSE,pty="s",mar=par.orig$mar)
   par(par.here)
   par(oma=c(0,0,0,0))
   on.exit(par(par.orig))
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #     Split the window according to whether to the sought position of the legends, and  #
   # whether to plot the vibrancy scale in addition to the triangle.                       #
   #---------------------------------------------------------------------------------------#
   lo.box = pretty.box(n=npanels)
   if (rgbsum.fixed && key.right){
      #------ Fixed RGB sum, legend to be to the right of the plots. ----------------------#
      dummy = layout( mat    = cbind(lo.box$mat.off,rep(1,times=lo.box$nrow))
                    , width  = c(rep((1.-f.key)/lo.box$ncol,times=lo.box$ncol),f.key)
                    )#end layout
      #------------------------------------------------------------------------------------#
   }else if (rgbsum.fixed){
      #------ Fixed RGB sum, legend to be underneath the plots. ---------------------------#
      dummy = layout( mat    = rbind(lo.box$mat.off,rep(1,times=lo.box$ncol))
                    , height = c(rep((1.-f.key)/lo.box$nrow,times=lo.box$nrow),f.key)
                    )#end layout
      #------------------------------------------------------------------------------------#
   }else if (key.right){
      #----- Make sure we have an even number of rows. ------------------------------------#
      if ((lo.box$nrow%%2) == 0){
         matnow = lo.box$mat.off2
      }else{
         matnow = apply(X=lo.box$mat.off2,MARGIN=2,FUN=rep,each=2)
      }#end if ((lo.box$nrow%%2) == 0)
      nrmat  = nrow(matnow)
      ncmat  = ncol(matnow)
      #------------------------------------------------------------------------------------#

      #------ Variable RGB sum, legends to be to the right of the plots. ------------------#
      dummy = layout( mat   = cbind(matnow,rep(c(1,2),each=nrmat/2))
                    , width = c(rep((1.-f.key)/ncmat,times=ncmat),f.key)
                    )#end layout
      #------------------------------------------------------------------------------------#
   }else{
      #----- Make sure we have an even number of columns. ---------------------------------#
      if ((lo.box$ncol%%2) == 0){
         matnow = lo.box$mat.off2
      }else{
         matnow = t(apply(X=lo.box$mat.off2,MARGIN=1,FUN=rep,each=2))
      }#end if ((lo.box$nrow%%2) == 0)
      nrmat  = nrow(matnow)
      ncmat  = ncol(matnow)
      #------------------------------------------------------------------------------------#

      #------ Variable RGB sum, legends to be underneath the plots. -----------------------#
      dummy = layout( mat    = rbind(matnow,rep(c(1,2),each=ncmat/2))
                    , height = c(rep((1.-f.key)/nrmat,times=nrmat),f.key)
                    )#end layout
      #------------------------------------------------------------------------------------#
   }#end if (rgbsum.fixed && leg.horizontal)
   #---------------------------------------------------------------------------------------#



   #=======================================================================================#
   #=======================================================================================#
   #    Plot the vibrancy palette (only in case RGB sum is not fixed).                     #
   #---------------------------------------------------------------------------------------#
   if (! rgbsum.fixed){
      nv.delta = round(1./v.delta)
      v.delta  = 1./nv.delta
      v.show   = seq(from=0,to=1.,by=v.delta)
      v.col    = hsv(h=0.5,s=0.0,v=mid.points(v.show))
      v.lwr    = v.show[-nv.delta]
      v.upr    = v.show[-1]
      if (key.right){
         v.side   = 4
         v.mar    = c(0.6,2.6,3.1,3.1)
      }else{
         v.side   = 1
         v.mar    = c(3.1,0.6,3.1,3.1)
      }#end if (key.right)

      #----- Plot legend (check whether to place it to the right or beneath the plot). ----#
      par(mar=v.mar)
      plot.new()
      plot.window(xlim=c(0,1),ylim=c(0,1),xaxs="i",yaxs="i")
      if (key.right){
         rect(xleft=0,ybottom=v.lwr,xright=1,ytop=v.upr,col=v.col,border=v.col)
      }else{
         rect(xleft=v.lwr,ybottom=0,xright=v.upr,ytop=1,col=v.col,border=v.col)
      }#end if (key.right)
      box()
      #------------------------------------------------------------------------------------#


      #----- Check whether there are specific instructions for plotting the legend axis. --#
      if (missing(v.axis.options)) {
         v.axis.now = list(side=v.side,las=1,...)
      }else{
         v.axis.now = modifyList(x=key.axis.options,val=list(side=v.side,las=1))
      }#end if
      do.call (what="axis",args=v.axis.now)
      #------------------------------------------------------------------------------------#


      #----- Plot the title. --------------------------------------------------------------#
      if (! is.null(v.title)) do.call(what="title",args=v.title)
      #------------------------------------------------------------------------------------#
   }#end if (! rgbsum.fixed)
   #=======================================================================================#
   #=======================================================================================#




   #=======================================================================================#
   #=======================================================================================#
   #    Plot the RGB palette.                                                              #
   #---------------------------------------------------------------------------------------#
      #------------------------------------------------------------------------------------#
      #     Make sure rgb.delta yields an integer number of bins.  In case it doesn't,     #
      # find the nearest number.                                                           #
      #------------------------------------------------------------------------------------#
      delta = 1./round(1./rgb.delta)
      #------------------------------------------------------------------------------------#


      #---- Find the triangle list. -------------------------------------------------------#
      red.span = seq(from=1,to=0,by=-delta)
      nrl=length(red.span)
      tri.val = list()
      tri.col = list()
      for (n in sequence(nrl)){
         green.l  = seq(from=0,to=1-red.span[n],by=0.5*delta)
         blue.l   = seq(from=0,to=1-red.span[n],by=0.5*delta)
         rgb.l    = expand.grid(red=red.span[n],green=green.l,blue=blue.l)
         keep     = ( rowSums(rgb.l) %>=% (1-sqrt(.Machine$double.eps))
                    & rowSums(rgb.l) %<=% (1+sqrt(.Machine$double.eps)) )
         rgb.l    = rgb.l[keep,] / rowSums(rgb.l[keep,])

         #---------------------------------------------------------------------------------#
         #     Append colours.                                                             #
         #---------------------------------------------------------------------------------#
         tri.val[[n]] = list(red=rgb.l$red,green=rgb.l$green,blue=rgb.l$blue)
         if (rgbsum.fixed){
            tri.col[[n]] = rgb( red           = rgb.l$red
                              , green         = rgb.l$green
                              , blue          = rgb.l$blue
                              , maxColorValue = 1
                              )#end rgb
         }else{
            hsv.l        = data.frame(t(rgb2hsv(r=rgb.l$red,g=rgb.l$green,b=rgb.l$blue)))
            hsv.l$v      = 1.0 + 0. * hsv.l$v
            tri.col[[n]] = hsv( h             = hsv.l$h
                              , s             = hsv.l$s
                              , v             = hsv.l$v
                              )#end hsv
         }#end if (rgbsum.fixed)
         #---------------------------------------------------------------------------------#
      }#end for
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #    Plot the colour legend.                                                         #
      #------------------------------------------------------------------------------------#
      triax.leg( at           = rgb.at
               , axis.labels  = rgb.axis.labels
               , mar          = c(0,0.5,0,0.5)
               , cex.axis     = 0.8
               , cex.ticks    = 0.5
               , show.axis    = rgb.show.axis
               , shift.labels = rgb.shift.labels
               , show.grid    = TRUE
               )#end triax.leg
      triax.fill(col=tri.col)
      #------------------------------------------------------------------------------------#


      #----- Plot the title. --------------------------------------------------------------#
      if (! is.null(rgb.title)) do.call(what="title",args=rgb.title)
      #------------------------------------------------------------------------------------#
   #=======================================================================================#
   #=======================================================================================#





   #=======================================================================================#
   #=======================================================================================#
   #      Now we plot the other panels.                                                    #
   #---------------------------------------------------------------------------------------#
   for (p in sequence(npanels)){
      #----- Set the window. --------------------------------------------------------------#
      
      plot.new()
      mar.orig    = par.here$mar
      mar.orig[4] = 1.1
      par(par.here)
      plot.window(xlim=xlim,ylim=ylim,...)
      #------------------------------------------------------------------------------------#



      #----- Find the corners for the rectangles. -----------------------------------------#
      nx      = length(x[[p]])
      ny      = length(y[[p]])
      xleft   = x[[p]] - 0.5 * dx[[p]]
      xright  = x[[p]] + 0.5 * dx[[p]]
      ybottom = y[[p]] - 0.5 * dy[[p]]
      ytop    = y[[p]] + 0.5 * dy[[p]]
      rect(xleft=xleft,ybottom=ybottom,xright=xright,ytop=ytop
          ,col=pcol[[p]],border=pcol[[p]])
      #------------------------------------------------------------------------------------#




      #---- Plot the X axis. --------------------------------------------------------------#
      if (! is.null(x.axis.options)){
         x.axis.now = modifyList(x=x.axis.options[[p]],val=list(side=1))
      }else{
         x.axis.now = list(side=1,las=1)
      }#end if
      do.call(what="axis",args=x.axis.now)
      #------------------------------------------------------------------------------------#




      #---- Plot the Y axis. --------------------------------------------------------------#
      if (! is.null(y.axis.options)){
         y.axis.now = modifyList(x=y.axis.options[[p]],val=list(side=2))
      }else{
         y.axis.now = list(side=2,las=1)
      }#end if
      do.call(what="axis",args=y.axis.now)
      #------------------------------------------------------------------------------------#



      #---- Plot the title. ---------------------------------------------------------------#
      if (! is.null(sub.options)){
         do.call(what="title",args=sub.options[[p]])
      }#end if
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Plot other options.                                                            #
      #------------------------------------------------------------------------------------#
      n.after = length(plot.after)
      for (a in sequence(n.after)){
          do.call(what=names(plot.after)[a],args=plot.after[[a]])
      }#end for
      #------------------------------------------------------------------------------------#


      #---- Lastly, the box. --------------------------------------------------------------#
      box()
      #------------------------------------------------------------------------------------#

   }#end for
   #=======================================================================================#
   #=======================================================================================#





   #=======================================================================================#
   #=======================================================================================#
   #     Plot the global title.                                                            #
   #---------------------------------------------------------------------------------------#
   par(las=0)
   if (! is.null(main.xlab)){
      #----- Make sure we get the main text. ----------------------------------------------#
      if (! is.list(main.xlab)){
         main.xlab=list(text=main.xlab)
      }else if (! "text" %in% names(main.xlab)){
         names(main.xlab)[[1]] = "text"
      }#end if
      #----- Outer must be set to TRUE, overwrite if needed be. ---------------------------#
      main.xlab$outer = TRUE
      if (! "side" %in% names(main.xlab)) main.xlab$side = 1
      if (! "padj" %in% names(main.xlab)) main.xlab$padj = -4.75
      do.call("mtext",main.xlab)
   }#end if
   #---------------------------------------------------------------------------------------#



   if (! is.null(main.ylab)){
      #----- Make sure we get the main text. ----------------------------------------------#
      if (! is.list(main.ylab)){
         main.ylab = list(text = main.ylab)
      }else if (! "text" %in% names(main.ylab)){
         names(main.ylab)[[1]] = "text"
      }#end if
      #----- Outer must be set to TRUE, overwrite if needed be. ---------------------------#
      main.ylab$outer = TRUE
      if (! "side" %in% names(main.ylab)) main.ylab$side = 2
      if (! "padj" %in% names(main.ylab)) main.ylab$padj = -0.75
      do.call("mtext",main.ylab)
   }#end if
   if (! is.null(main.title)){
      #----- Make sure we get the main text. ----------------------------------------------#
      if (! is.list(main.title)){
         main.title=list(text=main.title)
      }else if (! "text" %in% names(main.title)){
         names(main.title)[[1]] = "text"
      }#end if
      #----- Outer must be set to TRUE, overwrite if needed be. ---------------------------#
      main.title$outer = TRUE
      if (! "side" %in% names(main.title)) main.xlab$side = 3
      if (! "padj" %in% names(main.title)) main.xlab$padj = 0
      if (! "cex"  %in% names(main.title)) main.xlab$cex  = 1.1
      if (! "font" %in% names(main.title)) main.xlab$font = 2
      do.call("mtext",main.title)
   }#end if
   #---------------------------------------------------------------------------------------#

   invisible()
}#end function plot.rgb
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#    This function is borrowed from package plotrix.  It plots the axes for barycentric    #
# plots.                                                                                   #
#------------------------------------------------------------------------------------------#
triax.leg <<- function( x            = NULL
                      , main         = ""
                      , at           = seq(from=0.1,to=0.9,by=0.1)
                      , axis.labels  = NULL
                      , show.axis    = TRUE
                      , tick.labels  = NULL
                      , col.axis     = "black"
                      , cex.axis     = 1
                      , cex.ticks    = 1
                      , shift.labels = FALSE
                      , show.grid    = FALSE
                      , col.grid     = "grey"
                      , lty.grid     = par("lty")
                      , show.legend  = FALSE
                      , label.points = FALSE
                      , point.labels = NULL
                      , col.symbols  = "black"
                      , pch          = par("pch")
                      , mar          = c(5,2,4,2)
                      , no.add       = TRUE
                      , ...
                      ){

    #---- Save and restore the par settings once we leave the function. -------------------#
    oldpar = par(new=FALSE,no.readonly=TRUE)
    on.exit(par(oldpar))
    #--------------------------------------------------------------------------------------#


    #----- In case axes are not defined, use the column names. ----------------------------#
    if (is.null(axis.labels)) axis.labels <- colnames(x)[sequence(3)]
    #--------------------------------------------------------------------------------------#

    #----- Settings for the barycentric plot. ---------------------------------------------#
    par.send=par(pty="s",xpd = TRUE, mar = mar)
    #--------------------------------------------------------------------------------------#


    #----- Initialise the plot window. ----------------------------------------------------#
    plot.new()
    plot.window(xlim=c(0,1),ylim=c(0,1))
    #--------------------------------------------------------------------------------------#



    #----- Plot the barycentric frame. ----------------------------------------------------#
    triax.mesh( at           = at
              , axis.labels  = axis.labels
              , tick.labels  = tick.labels
              , col.axis     = col.axis
              , cex.axis     = cex.axis
              , cex.ticks    = cex.ticks
              , shift.labels = shift.labels
              , show.axis    = show.axis
              , show.grid    = show.grid
              , col.grid     = col.grid
              , lty.grid     = lty.grid
              )#end triax.mesh
    #--------------------------------------------------------------------------------------#



    #--------------------------------------------------------------------------------------#
    #    Add points in case they are needed.                                               #
    #--------------------------------------------------------------------------------------#
    if (is.null(x)){
        xypos = NULL
    }else{
        xypos = triax.points( x
                            , show.legend  = show.legend
                            , label.points = label.points
                            , point.labels = point.labels
                            , col.symbols  = col.symbols
                            , pch          = pch
                            , cc.axes      = FALSE
                            , par.send
                            ,...
                            )#end triax.points
    }#end if (is.null(x))
    #--------------------------------------------------------------------------------------#


    #----- Invisible output. --------------------------------------------------------------#
    invisible(list(xypos = xypos, oldpar = oldpar))
    #--------------------------------------------------------------------------------------#
}#end function triax.leg
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     This is pretty much triax.frame from plotrix, with a few tweaks to plot the colour   #
# ramp axes where they peak.                                                               #
#------------------------------------------------------------------------------------------#
triax.mesh <<- function ( at           = seq(from=0.1,to=0.9,by=0.1)
                        , tick.labels  = NULL
                        , axis.labels  = NULL
                        , col.axis     = "black"
                        , cex.axis     = 1
                        , cex.ticks    = 1
                        , show.axis    = TRUE
                        , show.grid    = FALSE
                        , shift.labels = FALSE
                        , col.grid     = "gray"
                        , lty.grid     = par("lty")
                        ){ #end triax.mesh

   #----- Define some handy constants. ----------------------------------------------------#
   sin60  = sin(pi/3)
   sin120 = sin(2*pi/3)
   #---------------------------------------------------------------------------------------#



   #----- Define tick marks for the bottom axis. ------------------------------------------#
   bx1 = at
   bx2 = bx1 + 0.01
   by1 = rep(x=0,times=9)
   by2 = rep(-0.02 * sin60, 9)
   #---------------------------------------------------------------------------------------#



   #----- Define tick marks for the left axis. --------------------------------------------#
   ly1 = at * sin60
   lx1 = bx1 * 0.5
   lx2 = lx1 - 0.02
   ly2 = ly1
   #---------------------------------------------------------------------------------------#



   #----- Define tick marks for the right axis. -------------------------------------------#
   rx1 = at * 0.5 + 0.5
   rx2 = rx1 + 0.01
   ry1 = rev(ly1)
   ry2 = rev(ly2) + 0.02 * sin60
   #---------------------------------------------------------------------------------------#



   #----- Draw grid in case the grid is sought. -------------------------------------------#
   if (show.grid){
      par(fg = col.grid)
      segments(x0=bx1     ,y0=by1     ,x1=lx1     ,y1=ly1     ,lty=lty.grid)
      segments(x0=lx1     ,y0=ly1     ,x1=rev(rx1),y1=rev(ry1),lty=lty.grid)
      segments(x0=rx1     ,y0=ry1     ,x1=bx1     ,y1=by1     ,lty=lty.grid)
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Plot tick labels. ---------------------------------------------------------------#
   par(fg = col.axis, xpd = TRUE)
   if (is.null(tick.labels)) tick.labels <- list(l = at, r = at, b = at)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Plot the axis labels.                                                              #
   #---------------------------------------------------------------------------------------#
   #----- Make sure that axes are in the right order. -------------------------------------#
   if (all(c("l","r","b") %in% names(axis.labels))){
      axis.labels = axis.labels[c("l","r","b")]
   }else{
      axis.labels = axis.labels[sequence(3)]
      names(axis.labels) = c("l","r","b")
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Plot left axis text. ------------------------------------------------------------#
   par(srt = 60)
   lnow = axis.labels["l"]
   if (shift.labels){
      xnow = 0.42
      ynow = 1.00
      anow = 1.00
   }else{
      xnow = 0.13
      ynow = 0.50
      anow = 0.50
   }#end (shift.labels)
   text(x=xnow,y=ynow,labels=lnow,adj=anow,cex=cex.axis,xpd=TRUE)
   #---------------------------------------------------------------------------------------#


   #----- Check whether to plot axis labels. ----------------------------------------------#
   if (show.axis){
      par(srt = 0)
      xoffset = 0.05
      yoffset = 0
      text(x=lx1-xoffset,y=ly1+yoffset,labels=tick.labels$l,cex=cex.ticks,xpd=TRUE)
   }#end if (show.axis)
   #---------------------------------------------------------------------------------------#



   #----- Plot right axis text. -----------------------------------------------------------#
   par(srt=300)
   lnow = axis.labels["r"]
   if (shift.labels){
      xnow = 1.13
      ynow = 0.00
      anow = 1.00
   }else{
      xnow = 0.86
      ynow = 0.52
      anow = 0.50
   }#end (shift.labels)
   text(x=xnow,y=ynow,labels=lnow,adj=anow,cex=cex.axis,xpd=TRUE)
   #---------------------------------------------------------------------------------------#


   #----- Check whether to plot axis labels. ----------------------------------------------#
   if (show.axis){
      par(srt = 60)
      xoffset = 0.015
      yoffset = 0.045
      text(rx2 + xoffset, ry1 + yoffset, tick.labels$r, cex = cex.ticks,xpd=TRUE)
   }#end if 
   #---------------------------------------------------------------------------------------#



   #----- Plot the bottom axis text. ------------------------------------------------------#
   par(srt = 0)
   lnow = axis.labels["b"]
   if (shift.labels){
      xnow =  0.00
      ynow = -0.14
      anow =  0.00
   }else{
      xnow =  0.50
      ynow = -0.14
      anow =  0.50
   }#end (shift.labels)
   text(x=xnow,y=ynow,labels=lnow,adj=anow,cex=cex.axis,xpd=TRUE)
   #---------------------------------------------------------------------------------------#


   #----- Check whether to plot axis labels. ----------------------------------------------#
   if (show.axis){
      par(srt = 300)
      xoffset = 0.03
      text(x=bx1+xoffset,y=by1-0.05,labels=rev(tick.labels$b),cex=cex.ticks,xpd=TRUE)
   }#end if 
   #---------------------------------------------------------------------------------------#


   #------ Plot the box. ------------------------------------------------------------------#
   x1 = c(0, 0, 0.5)
   x2 = c(1, 0.5, 1)
   y1 = c(0, 0, sin60)
   y2 = c(0, sin60, 0)
   par(fg = col.axis)
   segments(x0=x1 ,y0=y1 ,x1=x2 ,y1=y2 )
   segments(x0=bx1,y0=by1,x1=bx2,y1=by2)
   segments(x0=lx1,y0=ly1,x1=lx2,y1=ly2)
   segments(x0=rx1,y0=ry1,x1=rx2,y1=ry2)
   #---------------------------------------------------------------------------------------#
}#end function triax.mesh
#==========================================================================================#
#==========================================================================================#



#==========================================================================================#
#==========================================================================================#
#    Function that sets the colour scale, but replaces missing values with the NA colour.  #
#------------------------------------------------------------------------------------------#
rgb.norm <<- function(red,green,blue,maxColorValue=255,na.col="grey94"){

   #----- Normalise colours. --------------------------------------------------------------#
   red   = red   / maxColorValue
   green = green / maxColorValue
   blue  = blue  / maxColorValue
   fade  = ! (is.finite(red) & is.finite(green) & is.finite(blue) )
   #---------------------------------------------------------------------------------------#


   #----- Find the RGB for NA entries. ----------------------------------------------------#
   na.rgb = data.frame(t(col2rgb(na.col))) / 255
   #---------------------------------------------------------------------------------------#


   #----- Replace colours of missing points with NA colour. -------------------------------#
   red   = ifelse(test=fade,yes=na.rgb$red  ,no=red  )
   green = ifelse(test=fade,yes=na.rgb$green,no=green)
   blue  = ifelse(test=fade,yes=na.rgb$blue ,no=blue )
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Generate colour.  Make sure to set maxColorValue to 1 (they are already          #
   # normalised).                                                                          #
   #---------------------------------------------------------------------------------------#
   ans = rgb(red=red,green=green,blue=blue,maxColorValue=1)
   #---------------------------------------------------------------------------------------#


   return(ans)
}#end function rgb.norm
#==========================================================================================#
#==========================================================================================#
