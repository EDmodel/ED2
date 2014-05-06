#==========================================================================================#
#==========================================================================================#
#     This function plots a graph as a function of 3 parameters, with the colour scheme    #
# given.                                                                                   #
#------------------------------------------------------------------------------------------#
image.map <<- function( x
                      , y
                      , z
                      , dx               = if (is.list(x) & length(x) > 1){
                                              if (xlog){
                                                 lapply(lapply(lapply(lapply(lapply(x,sort),unique),log),diff),median,na.rm=TRUE)
                                              }else{
                                                 lapply(lapply(lapply(lapply(x,sort),unique),diff),median,na.rm=TRUE)
                                              }#end if
                                           }else if(xlog){
                                              median( diff(log(sort(unique(unlist(x)))))
                                                    , na.rm=TRUE )
                                           }else{
                                              median( diff(sort(unique(unlist(x))))
                                                    , na.rm=TRUE )
                                           }#end if
                      , dy               = if (is.list(y) & length(y) > 1){
                                              if (ylog){
                                                 lapply(lapply(lapply(lapply(lapply(y,sort),unique),log),diff),median,na.rm=TRUE)
                                              }else{
                                                 lapply(lapply(lapply(lapply(y,sort),unique),diff),median,na.rm=TRUE)
                                              }#end if
                                           }else if (ylog){
                                              median( diff(log(sort(unique(unlist(y)))))
                                                    , na.rm=TRUE )
                                           }else{
                                              median( diff(sort(unique(unlist(y))))
                                                    , na.rm=TRUE )
                                           }#end if
                      , xlim             = range(unlist(x),finite=TRUE)
                      , ylim             = range(unlist(y),finite=TRUE)
                      , zlim             = range(unlist(z),finite=TRUE)
                      , xlog             = FALSE
                      , ylog             = FALSE
                      , levels           = if (key.log){
                                              sort(unique(pretty.log(x=zlim,n=nlevels
                                                                    ,forcelog=TRUE)))
                                           }else{
                                              sort(unique(pretty    (x=zlim,n=nlevels)))
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
                      , legend.options   = NULL
                      , edge.axes        = FALSE
                      , oma              = NULL
                      , f.key            = 1/6
                      , f.leg            = 1/6
                      , off.xlab         = NULL
                      , off.right        = NULL
                      , xaxs             = "i"
                      , yaxs             = "i"
                      , smidgen          = 0
                      , interp.xyz       = FALSE
                      , useRaster        = TRUE
                      , byrow            = TRUE
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
      if (! is.null(x.axis.options)) x.axis.options = list(x.axis.options)
      if (! is.null(y.axis.options)) y.axis.options = list(y.axis.options)
      if (! is.null(sub.options   )) sub.options    = list(sub.options   )
      npanels = 1
   }else{
      npanels = length(x)
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     If legend is to be plotted, key.vertical has to be TRUE.  In case the user said   #
   # otherwise, return a warning.  Also, define offsets for X and Y according to the       #
   # legends and keys.                                                                     #
   #---------------------------------------------------------------------------------------#
   plot.legend = ! is.null(legend.options)
   if ( plot.legend && (! key.vertical)){
      warning(" key.vertical=FALSE ignored due to the legend.")
      key.vertical = TRUE
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Find the box structure for the panels. ------------------------------------------#
   lo.panel = pretty.box(npanels,byrow=byrow)
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
   par.orig = par(no.readonly=TRUE )
   mar.orig = par.orig$mar
   on.exit(par(par.orig))
   par(par.user)





   #----- Check for outer margins. --------------------------------------------------------#
   if (is.null(oma) && npanels == 1){
      par(oma=c(0,0,0,0))
   }else if (is.null(oma)){
      omd = c(0.02,1.00,0.01,0.93)
      par(omd=omd)
   }else{
      par(oma=oma)
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Find offset for x axis label and right, based on legends and keys and outer      #
   # margin .                                                                              #
   #---------------------------------------------------------------------------------------#
   par.tout = par(no.readonly=FALSE)
   #----- Bottom margin. ------------------------------------------------------------------#
   if (is.null(off.xlab)){
      if (plot.legend && key.vertical){
         off.xlab = with(par.tout,( omi[1] + f.leg * (din[2]-omi[1]-omi[3]) ) / din[2])
      }else if (key.vertical){
         off.xlab = with(par.tout,omi[1] / din[2])
      }else{
         off.xlab = with(par.tout,( omi[1] + f.key * (din[2]-omi[1]-omi[3]) ) / din[2])
      }#end if
   }#end if
   #----- Right margin. -------------------------------------------------------------------#
   if (is.null(off.right)){
      if (key.vertical){
         off.right = with(par.tout,( omi[4] + f.key * (din[1]-omi[2]-omi[4]) ) / din[1])
      }else if (plot.legend){
         off.right = with(par.tout,( omi[4] + f.leg * (din[1]-omi[2]-omi[4]) ) / din[1])
      }else{
         off.right = with(par.tout,omi[4] / din[1])
      }#end if
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Split the screen into multiple pieces (legend, key, plots...) -------------------#
   fh.panel = 1. - f.key
   fv.panel = 1. - f.leg
   if (npanels == 1 && plot.legend){
      layout( mat     = rbind(c(3, 2),c(1,0))
            , heights = c(fv.panel,f.leg)
            , widths  = c(fh.panel,f.key)
            )#end layout
   }else if (matrix.plot && plot.legend){
      layout( mat     = rbind( cbind(lo.panel$mat.off2,rep(2,times=lo.panel$nrow))
                             , c(rep(1,times=lo.panel$ncol),0)
                             )#end rbind
            , heights = c(rep(fv.panel/lo.panel$nrow,times=lo.panel$nrow),f.leg)
            , widths  = c(rep(fh.panel/lo.panel$ncol,times=lo.panel$ncol),f.key)
            )#end layout
   }else if (npanels == 1 && key.vertical){
      layout(mat=cbind(2, 1), widths = c(fh.panel,f.key))
   }else if (matrix.plot && key.vertical){
      layout( mat     = cbind(lo.panel$mat.off,rep(1,times=lo.panel$nrow))
            , widths  = c(rep(fh.panel/lo.panel$ncol,times=lo.panel$ncol),f.key)
            )#end layout
   }else if (matrix.plot){
      layout( mat     = rbind(lo.panel$mat.off,rep(1,times=lo.panel$ncol))
            , heights = c(rep(fh.panel/lo.panel$nrow,times=lo.panel$nrow),f.key)
            )#end layout
   }else{
      layout( mat     =rbind(1+sequence(npanels),rep(1,times=npanels))
            , heights = c(fh.panel,f.key)
            )#end layout
   }#end if
   #---------------------------------------------------------------------------------------#
   #=======================================================================================#
   #=======================================================================================#







   #=======================================================================================#
   #=======================================================================================#
   #      Check whether the "plot.after" list is shared by all plots or if it is one list  #
   # for each sub-plot.                                                                    #
   #---------------------------------------------------------------------------------------#
   if (is.null(plot.after)){
      same.for.all = TRUE
   }else{
      pa.names = names(plot.after)
      named    = ! is.null(pa.names)
      if (named){
         same.for.all = all(mapply(FUN=exists,x=pa.names,MoreArgs=list(mode="function")))
      }else{
         same.for.all = FALSE
      }#end if
   }#end if
   #=======================================================================================#
   #=======================================================================================#






   #=======================================================================================#
   #=======================================================================================#
   #      First plot: the legend.                                                          #
   #---------------------------------------------------------------------------------------#
   if (plot.legend){
      par(mar = c(0.1,0.1,0.1,0.1))
      plot.new()
      plot.window(xlim=c(0,1),ylim=c(0,1))
      do.call(what="legend",args=legend.options)
   }#end if
   #=======================================================================================#
   #=======================================================================================#






   #=======================================================================================#
   #=======================================================================================#
   #      Second plot: the key scale.                                                      #
   #---------------------------------------------------------------------------------------#
      if (key.vertical){
         par(mar = lo.panel$mar.key)
      }else{
         par(mar = c(2.1,4.6,1.6,2.1))
      }#end if
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
            plot.window(xlim=range(levels),ylim=c(0,1),xaxs="i",yaxs="i",las=1,log="x")
         }else{
            plot.window(xlim=range(levels),ylim=c(0,1),xaxs="i",yaxs="i",las=1)
         }#end if
         #---------------------------------------------------------------------------------#


         #----- Draw the colour bar. ------------------------------------------------------#
         rect(xleft=levels[-length(levels)],ybottom=0,xright=levels[-1],ytop=1
             ,col=col,border=col)
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
      if (! is.null(key.title)) do.call(what="title",args=key.title)
      #------------------------------------------------------------------------------------#
   #=======================================================================================#
   #=======================================================================================#





   #=======================================================================================#
   #=======================================================================================#
   #      Now we plot the other panels.                                                    #
   #---------------------------------------------------------------------------------------#
   for (p in sequence(npanels)){
      #----- Find out where the box goes, and set up axes and margins. --------------------#
      left    = lo.panel$left  [p]
      right   = lo.panel$right [p]
      top     = lo.panel$top   [p]
      bottom  = lo.panel$bottom[p]
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
         if (! key.vertical) mar.now[1] = 4.1
         #---------------------------------------------------------------------------------#
      }#end if
      plog = ""
      if (xlog) plog = paste(plog,"x",sep="")
      if (ylog) plog = paste(plog,"y",sep="")
      par(mar = mar.now)
      plot.new()
      plot.window(xlim=xlim,ylim=ylim,log=plog,xaxs=xaxs,yaxs=yaxs,...)
      #------------------------------------------------------------------------------------#



      #----- Find the corners for the rectangles. -----------------------------------------#
      if (interp.xyz){

         #----------------------------------------------------------------------------------#
         #      We use image to plot, so it looks nice in PDF.                              #
         #----------------------------------------------------------------------------------#
         useRaster.now = useRaster && (! xlog) && (! ylog) 
         #----- Make x and y dimensionless. ------------------------------------------------#
         if (xlog){
            xx    = log(as.numeric(x[[p]]))
         }else{
            xx    = as.numeric(x[[p]])
         }#end if
         if (ylog){
            yy    = log(as.numeric(y[[p]]))
         }else{
            yy    = as.numeric(y[[p]])
         }#end if
         zz    = z[[p]]
         nx    = length(xx)
         ny    = length(yy)
         xlow  = min(xx)
         xhigh = max(xx)
         ylow  = min(yy)
         yhigh = max(yy)
         #----- Scale x and y. -------------------------------------------------------------#
         xxx    = ( xx - xlow ) / ( xhigh - xlow )
         yyy    = ( yy - ylow ) / ( yhigh - ylow )
         sss    = is.finite(zz)
         xo     = seq(from=0,to=1,length.out=10*length(unique(xx)))
         yo     = seq(from=0,to=1,length.out=10*length(unique(yy)))

         if (any(is.finite(zz[sss]))){
            zint   = interp(x=xxx[sss],y=yyy[sss],z=zz[sss],xo=xo,yo=yo)
            sint   = try(interp(x=xxx     ,y=yyy     ,z=sss    ,xo=xo,yo=yo))
         }else{
            zint   = list(x=xo,y=yo,z=matrix(nrow=length(xo),ncol=length(yo)))
            sint   = list(x=xo,y=yo,z=matrix(nrow=length(xo),ncol=length(yo)))
         }#end if
         sint$z = sint$z %>% twothirds
         zint$z = ifelse(sint$z,zint$z,NA)
         zint$x = xlow + zint$x * (xhigh - xlow)
         zint$y = ylow + zint$y * (yhigh - ylow)

         if (xlog) zint$x = exp(zint$x)
         if (ylog) zint$y = exp(zint$y)

         image(zint,breaks=levels,col=col,add=TRUE,useRaster=useRaster.now)
      }else{
         #---------------------------------------------------------------------------------#
         #    Split zleft into the breaks defined by the colour palette.                   #
         #---------------------------------------------------------------------------------#
         zcut              = try(cut(as.numeric(z[[p]]),breaks=levels))
         zlev              = levels(zcut)
         zcol              = col[match(zcut,zlev)]
         zcol[is.na(zcol)] = na.col
         #---------------------------------------------------------------------------------#

         nx      = length(x[[p]])
         ny      = length(y[[p]])
         if (xlog){
            xleft   = exp(log(x[[p]]) - 0.5 * (1. + smidgen) * dx[[p]])
            xright  = exp(log(x[[p]]) + 0.5 * (1. + smidgen) * dx[[p]])
         }else{
            xleft   = x[[p]] - 0.5 * (1. + smidgen) * dx[[p]]
            xright  = x[[p]] + 0.5 * (1. + smidgen) * dx[[p]]
         }#end if
         if (ylog){
            ybottom = exp(log(y[[p]]) - 0.5 * (1. + smidgen) * dy[[p]])
            ytop    = exp(log(y[[p]]) + 0.5 * (1. + smidgen) * dy[[p]])
         }else{
            ybottom = y[[p]] - 0.5 * (1. + smidgen) * dy[[p]]
            ytop    = y[[p]] + 0.5 * (1. + smidgen) * dy[[p]]
         }#end if
         rect( xleft   = xleft
             , ybottom = ybottom
             , xright  = xright
             , ytop    = ytop
             , col     = zcol
             , border  = zcol
             , xpd     = FALSE
             )#end rect
      }#end if
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
      if (! is.null(sub.options) && npanels != 1){
         do.call(what="title",args=sub.options[[p]])
      }#end if
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Plot other options.  Check use a shared list, or one list for each sub-plot.   #
      #------------------------------------------------------------------------------------#
      if (same.for.all){
         n.after = length(plot.after)
         for (a in sequence(n.after)){
             do.call(what=names(plot.after)[a],args=plot.after[[a]])
         }#end for
      }else{
         n.after = length(plot.after[[p]])
         for (a in sequence(n.after)){
            do.call(what=names(plot.after[[p]])[a],args=plot.after[[p]][[a]])
         }#end for
      }#end if
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
   if (! is.null(main.title)){
      #----- Make sure we get the main text. ----------------------------------------------#
      if (! is.list(main.title)){
         main.title=list(main=main.title)
      }else if (! "main" %in% names(main.title)){
         names(main.title)[[1]] = "main"
      }#end if
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #       Check whether to use title or gtitle.                                        #
      #------------------------------------------------------------------------------------#
      if (npanels == 1){
         #------ Convert subtitle options into true subtitle options. ---------------------#
         if (! is.null(sub.options)){
            names(sub.options[[1]]) = gsub( pattern     = "main"
                                          , replacement = "sub"
                                          , x           = names(sub.options[[1]])
                                          )#end gsub
         }#end if (! is.null(sub.options))
         main.title              = modifyList(x=main.title,val=list(sub.options[[1]]))
         do.call(what="title",args=main.title)
      }else{
         main.title = modifyList( x   = main.title
                                , val = list(off.xlab=off.xlab,off.right=off.right)
                                )#end modifyList
         do.call(what="gtitle",args=main.title)
      }#end if
      #------------------------------------------------------------------------------------#
   }#end if
   #=======================================================================================#
   #=======================================================================================#




   invisible()
}#end function image.map
#==========================================================================================#
#==========================================================================================#
