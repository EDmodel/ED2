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
                     , rgb.at           = c(0.1,0.9)
                     , rgb.axis.labels  = c("Green","Red","Blue")
                     , x.axis.options   = NULL
                     , y.axis.options   = NULL
                     , sub.options      = NULL
                     , main.title       = NULL
                     , main.xlab        = NULL
                     , main.ylab        = NULL
                     , key.title        = NULL
                     , plot.after       = NULL
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


   #---------------------------------------------------------------------------------------#
   #     Find the colours.                                                                 #
   #---------------------------------------------------------------------------------------#
   rgb.max = max(unlist(r)+unlist(g)+unlist(b),na.rm=TRUE)
   pcol    = mapply( FUN      = rgb
                   , red      = r
                   , green    = g
                   , blue     = b
                   , MoreArgs = list(maxColorValue=rgb.max)
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





   #---- Split the window into two. -------------------------------------------------------#
   lo.box = pretty.box(n=npanels)
   layout( mat    = cbind(1+lo.box$mat,rep(1,times=lo.box$nrow))
         , width  = c(rep(7/lo.box$ncol,times=lo.box$ncol),2)
         )#end layout
   #---------------------------------------------------------------------------------------#




   #---- Find the triangle list. ----------------------------------------------------------#
   red.span = seq(from=1,to=0,by=-0.1)
   nrl=length(red.span)
   tri.val = list()
   tri.col = list()
   for (n in 1:nrl){
      green.l  = seq(from=0,to=1-red.span[n],by=0.05)
      blue.l   = seq(from=0,to=1-red.span[n],by=0.05)
      rgb.l    = expand.grid(red=red.span[n],green=green.l,blue=blue.l)
      keep     = ( rowSums(rgb.l) >= 1-sqrt(.Machine$double.eps)
                 & rowSums(rgb.l) <= 1+sqrt(.Machine$double.eps) )
      rgb.l    = rgb.l[keep,] / rowSums(rgb.l[keep,])
      tri.val[[n]] = list(red=rgb.l$red,green=rgb.l$green,blue=rgb.l$blue)
      tri.col[[n]] = rgb(red=rgb.l$red,green=rgb.l$green,blue=rgb.l$blue,maxColorValue=1)
   }#end for
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Plot the colour legend.                                                            #
   #---------------------------------------------------------------------------------------#
   triax.leg( at          = rgb.at
            , axis.labels = rgb.axis.labels
            , mar         = c(0,0.5,0,0.5)
            , cex.axis    = 0.8
            , cex.ticks   = 0.5
            , show.grid   = TRUE
            )#end triax.leg
   triax.fill(col=tri.col)
   #---------------------------------------------------------------------------------------#





   #=======================================================================================#
   #=======================================================================================#
   #      Now we plot the other panels.                                                    #
   #---------------------------------------------------------------------------------------#
   for (p in 1:npanels){
      #----- Set the window. --------------------------------------------------------------#
      
      plot.new()
      mar.orig    = par.here$mar
      mar.orig[4] = 1.1
      par(par.here)
      plot.window(xlim=xlim,ylim=ylim,...)
      box()
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
triax.leg <<- function (x = NULL, main = "", at = seq(0.1, 0.9, by = 0.1),
    axis.labels = NULL, tick.labels = NULL, col.axis = "black", 
    cex.axis = 1, cex.ticks = 1, align.labels = TRUE, show.grid = FALSE, 
    col.grid = "gray", lty.grid = par("lty"), cc.axes = FALSE, 
    show.legend = FALSE, label.points = FALSE, point.labels = NULL, 
    col.symbols = "black", pch = par("pch"), mar = c(5, 2, 4, 
        2), no.add = TRUE, ...) 
{
    oldpar <- par(new=FALSE,no.readonly=TRUE)
    on.exit(par(oldpar))
    if (is.null(axis.labels)) axis.labels <- colnames(x)[1:3]
    par.send=par(pty="s",xpd = TRUE, mar = mar)
    plot(0.5, type = "n", axes = FALSE, xlim = c(0, 1), ylim = c(0, 
        1), main = main, xlab = "", ylab = "")
    triax.frame(at = at, axis.labels = axis.labels, tick.labels = tick.labels, 
        col.axis = col.axis, cex.axis = cex.axis, cex.ticks = cex.ticks, 
        align.labels = align.labels, show.grid = show.grid, col.grid = col.grid, 
        lty.grid = lty.grid, cc.axes = cc.axes)
    if (is.null(x)) 
        xypos <- NULL
    else xypos <- triax.points(x, show.legend = show.legend, 
        label.points = label.points, point.labels = point.labels, 
        col.symbols = col.symbols, pch = pch, cc.axes = cc.axes,par.send,...)
    invisible(list(xypos = xypos, oldpar = oldpar))
}
#==========================================================================================#
#==========================================================================================#
