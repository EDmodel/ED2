#==========================================================================================#
#==========================================================================================#
#     This function plots a graph as a function of 3 parameters, with the colour scheme    #
# given.                                                                                   #
#------------------------------------------------------------------------------------------#
image.map <<- function( x
                      , y
                      , z
                      , dx               = if (is.list(x) & length(x) > 1){
                                              lapply(lapply(x,diff),median,na.rm=TRUE)
                                           }else{
                                              median( diff(sort(unique(unlist(x))))
                                                    , na.rm=TRUE )
                                           }#end if
                      , dy               = if (is.list(y) & length(y) > 1){
                                              lapply(lapply(y,diff),median,na.rm=TRUE)
                                           }else{
                                              median( diff(sort(unique(unlist(y))))
                                                    , na.rm=TRUE )
                                           }#end if
                      , xlim             = range(unlist(x),finite=TRUE)
                      , ylim             = range(unlist(y),finite=TRUE)
                      , zlim             = range(unlist(z),finite=TRUE)
                      , levels           = if (key.log){
                                              pretty.log(x=zlim,n=nlevels)
                                           }else{
                                             pretty(x=zlim,n=nlevels)
                                           }#end if
                      , nlevels          = 20
                      , colour.palette   = cm.colors
                      , col              = colour.palette(length(levels)-1)
                      , na.col           = "grey94"
                      , key.log          = FALSE
                      , key.vertical     = TRUE
                      , x.axis.options   = NULL
                      , y.axis.options   = NULL
                      , key.axis.options = NULL
                      , key.options      = NULL
                      , sub.options      = NULL
                      , main.title       = NULL
                      , main.xlab        = NULL
                      , main.ylab        = NULL
                      , key.title        = NULL
                      , plot.after       = NULL
                      , matrix.plot      = FALSE
                      , edge.axes        = FALSE
                      , oma              = if (is.list(x) && length(x) > 1){
                                              c( 2.0 - 2.0 * is.null(main.xlab)
                                               , 2.0 - 2.0 * is.null(main.ylab)
                                               , 2.0 - 2.0 * is.null(main.title)
                                               , 0.0
                                               )#end c
                                           }else{
                                              c(0.,0.,0.,0.)
                                           }#end if
                      , f.key            = 1/6
                      , xaxs             = "i"
                      , yaxs             = "i"
                      , ...
                      ){


   #---------------------------------------------------------------------------------------#
   #     Find out whether x, y, r, g, and b are single values or lists.                    #
   #---------------------------------------------------------------------------------------#
   if (missing(x) || missing(y) || missing(z)){
      cat(" - x is missing: ",missing(x),"\n")
      cat(" - y is missing: ",missing(y),"\n")
      cat(" - z is missing: ",missing(z),"\n")
      stop(" x, y, and z must be given")
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Check whether x, y, and z are the same type of data.                             #
   #---------------------------------------------------------------------------------------#
   same.kind = (is.list(x) == is.list(y) && is.list(x) == is.list(y))
   if (! same.kind){
      cat(" X is list: ",is.list(x),"\n")
      cat(" Y is list: ",is.list(y),"\n")
      cat(" Z is list: ",is.list(z),"\n")
      stop ("X, Y, and Z must be of the same kind...")
   }else if (!is.list(x)){
      #----- Convert x, y, and z to lists. ------------------------------------------------#
      x              = list(x )
      y              = list(y )
      z              = list(z )
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


   #----- Find the box structure for the panels. ------------------------------------------#
   lo.panel = pretty.box(npanels)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     x and y should be axes for r, g, and b.  Check whether they are or not.           #
   #---------------------------------------------------------------------------------------#
   nx = sapply(X = x, FUN = length)
   ny = sapply(X = y, FUN = length)
   nz = sapply(X = z, FUN = length)
   if ( any(nz != nx) || any(nz != ny)){
      cat(" - length(x): ",paste(nx,sep=" "),"\n")
      cat(" - length(y): ",paste(ny,sep=" "),"\n")
      cat(" - length(z): ",paste(nz,sep=" "),"\n")
      stop(" x, y, and z must have the same length")
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Save the margins to avoid losing the data. --------------------------------------#
   par.orig = par(no.readonly=TRUE)
   mar.orig = par.orig$mar
   on.exit(par(par.orig))
   par(par.user)
   par(oma=oma)
   #---------------------------------------------------------------------------------------#



   #----- Split the screen into 3, the two panels and the scale. --------------------------#
   f.panel = 1 - f.key
   if (npanels == 1 && key.vertical){
      layout(matrix(c(2, 1), nc = 2), widths = c(f.panel,f.key))
      mar.key = c(4.1,1.0,4.1,4.1)
   }else if (matrix.plot && key.vertical){
      layout( mat     = cbind(lo.panel$mat.off,rep(1,times=lo.panel$nrow))
            , widths  = c(rep(f.panel/lo.panel$ncol,times=lo.panel$ncol),f.key)
            )#end layout
      mar.key = c(4.1,1.0,4.1,4.1)
   }else if (matrix.plot){
      layout( mat     = rbind(lo.panel$mat.off,rep(1,times=lo.panel$ncol))
            , heights = c(rep(f.panel,times=lo.panel$nrow),f.key)
            )#end layout
      mar.key = c(2.6,4.1,2.1,2.1)
   }else{
      layout( mat     =rbind(seq(from=2,to=npanels+1),rep(1,times=npanels))
            , heights = c(f.panel,f.key)
            )#end layout
      mar.key = c(2.6,4.1,2.1,2.1)
   }#end if
   #---------------------------------------------------------------------------------------#
   #=======================================================================================#
   #=======================================================================================#





   #=======================================================================================#
   #=======================================================================================#
   #      First plot: the key scale.                                                       #
   #---------------------------------------------------------------------------------------#
      par(mar = lo.panel$mar.key)
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
             ,border="transparent")
         #---------------------------------------------------------------------------------#

         #----- Check whether there are specific instructions for plotting the key axis. --#
         if (missing(key.axis.options)) {
            key.now = list(side=4,las=1,...)
         }else{
            key.now = modifyList(x=key.axis.options,val=list(side=4,las=1))
         }#end if
         do.call (what="axis",args=key.now)
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
             ,col=col,border="transparent")
         #---------------------------------------------------------------------------------#


         #----- Check whether there are specific instructions for plotting the key axis. --#
         if (missing(key.axis.options)) {
            key.now = list(side=1,las=1,...)
         }else{
            key.now = modifyList(x=key.axis.options,val=list(side=1,las=1))
         }#end if
         do.call (what="axis",args=key.now)
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
      #----- Find out where the box goes, and set up axes and margins. --------------------#
      left    = (p %% lo.panel$ncol) == 1
      right   = (p %% lo.panel$ncol) == 0
      top     = p <= lo.panel$ncol
      bottom  = p > (lo.panel$nrow - 1) * lo.panel$ncol
      #------------------------------------------------------------------------------------#


      #----- Set the window. --------------------------------------------------------------#
      if (matrix.plot & edge.axes){
         if (left && right){
            mar.left  = 4.1
            mar.right = 2.1
         }else if (left){
            mar.left  = 3.1
            mar.right = 0.1
         }else if (right){
            mar.left  = 0.1
            mar.right = 3.1
         }else{
            mar.left  = 1.6
            mar.right = 1.6
         }#end if
         if (bottom && top){
            mar.bottom = 5.1
            mar.top    = 4.1
         }else if (bottom){
            mar.bottom = 3.1
            mar.top    = 1.1
         }else if (top){
            mar.bottom = 1.1
            mar.top    = 3.1
         }else{
            mar.bottom = 1.1
            mar.top    = 1.1
         }#end if
         mar.now = c(mar.bottom,mar.left,mar.top,mar.right)
         #---------------------------------------------------------------------------------#
      }else if (matrix.plot){
         mar.left   = 3.1 + 1.0 * (npanels == 1)
         mar.right  = 1.1 + 1.0 * (npanels == 1)
         mar.bottom = 4.1 + 1.0 * (npanels == 1)
         mar.top    = 3.1 + 1.0 * (npanels == 1)
         mar.now = c(mar.bottom,mar.left,mar.top,mar.right)
      }else{
         #----- Find out where the box goes, and set up axes and margins. -----------------#
         left    = TRUE
         right   = TRUE
         bottom  = TRUE
         top     = TRUE
         mar.now = mar.orig
         if (! key.vertical) mar[1] = 4.1
         #---------------------------------------------------------------------------------#
      }#end if
      par(mar = mar.now)
      plot.new()
      plot.window(xlim=xlim,ylim=ylim,xaxs=xaxs,yaxs=yaxs,...)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    Split zleft into the breaks defined by the colour palette.                      #
      #------------------------------------------------------------------------------------#
      zcut              = cut(z[[p]],breaks=levels)
      zlev              = levels(zcut)
      zcol              = col[match(zcut,zlev)]
      zcol[is.na(zcol)] = na.col
      #------------------------------------------------------------------------------------#



      #----- Find the corners for the rectangles. -----------------------------------------#
      nx      = length(x[[p]])
      ny      = length(y[[p]])
      xleft   = x[[p]] - 0.5 * dx[[p]]
      xright  = x[[p]] + 0.5 * dx[[p]]
      ybottom = y[[p]] - 0.5 * dy[[p]]
      ytop    = y[[p]] + 0.5 * dy[[p]]
      rect(xleft=xleft,ybottom=ybottom,xright=xright,ytop=ytop,col=zcol,border=zcol)
      #------------------------------------------------------------------------------------#




      #---- Plot the X axis. --------------------------------------------------------------#
      if (bottom){
         if (! is.null(x.axis.options)){
            x.axis.now = modifyList(x=x.axis.options[[p]],val=list(side=1))
         }else{
            x.axis.now = list(side=1,las=1)
         }#end if
         do.call(what="axis",args=x.axis.now)
      }#end if
      #------------------------------------------------------------------------------------#




      #---- Plot the Y axis. --------------------------------------------------------------#
      if (left){
         if (! is.null(y.axis.options)){
            y.axis.now = modifyList(x=y.axis.options[[p]],val=list(side=2))
         }else{
            y.axis.now = list(side=2,las=1)
         }#end if
         do.call(what="axis",args=y.axis.now)
      }#end if
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


      #----- Lastly, add the box (so it stays on top). ------------------------------------#
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
      main.xlab$outer = par("oma")[1] > 0
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
      main.ylab$outer = par("oma")[2] > 0
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
      main.title$outer = par("oma")[3] > 0
      if (! "side" %in% names(main.title)) main.xlab$side = 3
      if (! "padj" %in% names(main.title)) main.xlab$padj = 0
      if (! "cex"  %in% names(main.title)) main.xlab$cex  = 1.1
      if (! "font" %in% names(main.title)) main.xlab$font = 2
      do.call("mtext",main.title)
   }#end if
   #---------------------------------------------------------------------------------------#




   invisible()
}#end function image.map
#==========================================================================================#
#==========================================================================================#
