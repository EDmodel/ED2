#==========================================================================================#
#==========================================================================================#
#     This function plots a density point cloud that represents the point density at any   #
# given area of the graph.                                                                 #
#------------------------------------------------------------------------------------------#
xy.density <<- function( x
                       , y
                       , xlim             = if (xlog){
                                               range(pretty.log(x))
                                            }else{
                                               range(pretty(x))
                                            }#end if (xlog)
                       , ylim             = if (ylog){
                                               range(pretty.log(y))
                                            }else{
                                               range(pretty(y))
                                            }#end if (xlog)
                       , xlevels          = NULL
                       , ylevels          = NULL
                       , zlim             = NULL
                       , xlog             = FALSE
                       , ylog             = FALSE
                       , nbins            = 80
                       , colour.palette   = cm.colors
                       , nlevels          = 20
                       , plot.key         = TRUE
                       , key.log          = FALSE
                       , key.vertical     = TRUE
                       , x.axis.options   = NULL
                       , y.axis.options   = NULL
                       , key.axis.options = NULL
                       , key.options      = NULL
                       , sub.options      = NULL
                       , main.title       = NULL
                       , key.title        = NULL
                       , plot.after       = NULL
                       , legend.options   = NULL
                       , edge.axes        = FALSE
                       , oma              = NULL
                       , omd              = NULL
                       , f.key            = 1/6
                       , f.leg            = 1/6
                       , off.xlab         = NULL
                       , off.right        = NULL
                       , xaxs             = "i"
                       , yaxs             = "i"
                       , method           = c("table","density")
                       , zignore          = 0.0001
                       , mar.main         = c(4.1,4.1,4.1,1.1)
                       , mar.key          = NULL
                       , useRaster        = ! (xlog || ylog)
                       , reparse          = TRUE
                       , add              = FALSE
                       , ...
                       ){


   #---------------------------------------------------------------------------------------#
   #     Find out whether x and y are both provided.                                       #
   #---------------------------------------------------------------------------------------#
   if (missing(x) || missing(y)){
      cat(" - x is missing: ",missing(x),"\n")
      cat(" - y is missing: ",missing(y),"\n")
      stop(" Both x and y must be provided.")
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Standardise method.                                                               #
   #---------------------------------------------------------------------------------------#
   method = match.arg(method)
   off    = (method %in% "density") 
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    In case this plot is using and existing plot.window, make sure we use the window   #
   # settings.                                                                             #
   #---------------------------------------------------------------------------------------#
   if (add){
      #------ Retrieve direct settings. ---------------------------------------------------#
      xlog   = par("xlog")
      ylog   = par("ylog")
      xa     = par("usr")[1]
      xz     = par("usr")[2]
      ya     = par("usr")[3]
      yz     = par("usr")[4]
      xaxs   = par("xaxs")
      yaxs   = par("yaxs")
      xfac   = ifelse(test=xaxs %in% "r",yes=0.04,no=0.00)
      yfac   = ifelse(test=yaxs %in% "r",yes=0.04,no=0.00)
      #------------------------------------------------------------------------------------#



      #----- Find bounds in the X direction. ----------------------------------------------#
      xlim = c( (1.+xfac)*xa + xfac*xz, xfac*xa + (1+xfac)*xz ) / (1.+2.*xfac)
      ylim = c( (1.+yfac)*ya + yfac*yz, yfac*ya + (1+yfac)*yz ) / (1.+2.*yfac)
      if (xlog) xlim=10^xlim
      if (ylog) ylim=10^ylim
      #------------------------------------------------------------------------------------#
   }#end if 
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #    Split the domain into bins, and count points.                                      #
   #---------------------------------------------------------------------------------------#
   if (is.null(xlevels) && xlog){
      xlwr    = log(xlim[1])-eps()*diff(log(xlim))
      xupr    = log(xlim[2])+eps()*diff(log(xlim))
      xlevels = exp(seq(from=xlwr,to=xupr,length.out=nbins+off))
      xdens   = exp(mid.points(log(xlevels)))
   }else if (is.null(xlevels)){
      xlwr    = xlim[1]-eps()*diff(xlim)
      xupr    = xlim[2]+eps()*diff(xlim)
      xlevels = seq(from=xlwr,to=xupr,length.out=nbins+off)
      xdens   = mid.points(xlevels)
   }else if (xlog){
      xlwr    = min(log(xlevels),na.rm=TRUE)
      xupr    = max(log(xlevels),na.rm=TRUE)
      xdens   = exp(mid.points(log(xlevels)))
   }else{
      xlwr    = min(xlevels,na.rm=TRUE)
      xupr    = max(xlevels,na.rm=TRUE)
      xdens   = mid.points(xlevels)
   }#end if (is.null(xlevels))
   if (is.null(ylevels) && ylog){
      ylwr    = log(ylim[1])-eps()*diff(log(ylim))
      yupr    = log(ylim[2])+eps()*diff(log(ylim))
      ylevels = exp(seq(from=ylwr,to=yupr,length.out=nbins+off))
      ydens   = exp(mid.points(log(ylevels)))
   }else if (is.null(ylevels)){
      ylwr    = ylim[1]-eps()*diff(ylim)
      yupr    = ylim[2]+eps()*diff(ylim)
      ylevels = seq(from=ylwr,to=yupr,length.out=nbins+off)
      ydens   = mid.points(ylevels)
   }else if (ylog){
      ylwr    = min(log(ylevels),na.rm=TRUE)
      yupr    = max(log(ylevels),na.rm=TRUE)
      ydens   = exp(mid.points(log(ylevels)))
   }else{
      ylwr    = min(ylevels,na.rm=TRUE)
      yupr    = max(ylevels,na.rm=TRUE)
      ydens   = mid.points(ylevels)
   }#end if (is.null(ylevels))
   #---------------------------------------------------------------------------------------#


   #------ Cut x and y points into the bins, then use table to count occurrences. ---------#
   if (method %in% "table"){
      xcut       = as.integer(cut(x,breaks=xlevels,labels=xdens))
      ycut       = as.integer(cut(y,breaks=ylevels,labels=ydens))
      ztable     = table(xcut,ycut)
      idx        = cbind( row = as.integer(rownames(ztable)[row(ztable)])
                        , col = as.integer(colnames(ztable)[col(ztable)])
                        )#end cbind
      zdens      = matrix(data=0,nrow=length(xdens),ncol=length(ydens))
      zdens[idx] = c(as.matrix(ztable))
      zdens      = 100. * zdens / sum(zdens)
   }else{
      xkde    = if(xlog){log(x)}else{x}
      ykde    = if(ylog){log(y)}else{y}
      zdens   = kde2d(x=xkde,y=ykde,n=nbins,lims=c(xlwr,xupr,ylwr,yupr))
      zdens   = 100. * zdens[[3]] / sum(zdens[[3]])
   }#end if (method %in% "table")
   zlwr    = zignore * max(zdens,na.rm=TRUE)
   zdens   = 0. * zdens + ifelse(test=zdens %ge% zlwr,yes=zdens,no=0.)
   #---------------------------------------------------------------------------------------#



   #------ Find colour levels. ------------------------------------------------------------#
   if (key.log){
      if (is.null(zlim)){
         zlim = range(pretty.log(zdens))
      }#end if (is.null(zlim))
      clevels = sort(unique(pretty.log(x=zlim,n=nlevels,forcelog=TRUE)))
   }else{
      if (is.null(zlim)){
         zlim = range(pretty(zdens))
      }#end if (is.null(zlim))
      clevels = sort(unique(pretty(x=zlim,n=nlevels)))
   }#end if
   ccolours = colour.palette(length(clevels)-1)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Skip this part in case we are adding another plot.                                #
   #---------------------------------------------------------------------------------------#
   if (! add){

      #------------------------------------------------------------------------------------#
      #     If legend is to be plotted, key.vertical has to be TRUE.  In case the user     #
      # said otherwise, return a warning.  Also, define offsets for X and Y according to   #
      # the legends and keys.                                                              #
      #------------------------------------------------------------------------------------#
      plot.legend = ! is.null(legend.options)
      if ( plot.legend && (! key.vertical)){
         warning(" key.vertical=FALSE ignored due to the legend.")
         key.vertical = TRUE
      }#end if
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Find key margins.                                                              #
      #------------------------------------------------------------------------------------#
      if (key.vertical && is.null(mar.key)){
         mar.key = c(4.1,0.1,4.1,4.1)
      }else if (is.null(mar.key)){
         mar.key = c(4.1,4.1,0.6,1.1)
      }#end if
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Coerce x, y, and key axis options, and key and main title options into lists.  #
      #------------------------------------------------------------------------------------#
      if (is.null(x.axis.options)){
         x.axis.options = list(side=1,las=1)
      }else{
         x.axis.options = as.list(x.axis.options)
         if (! "side" %in% names(x.axis.options)){
            x.axis.options = modifyList(x=x.axis.options,val=list(side=1))
         }#end if (! "side" %in% names(x.axis.options))
      }#end if
      if (is.null(y.axis.options)){
         y.axis.options = list(side=2,las=1)
      }else{
         y.axis.options = as.list(y.axis.options)
         if (! "side" %in% names(y.axis.options)){
            y.axis.options = modifyList(x=y.axis.options,val=list(side=2))
         }#end if (! "side" %in% names(y.axis.options))
      }#end if
      if (is.null(key.axis.options)){
         if (key.log){
            zat     = pretty.log(zlim)
            zlabels = sprintf("%g",zat)
         }else{
            zat     = pretty(zlim)
            zlabels = sprintf("%g",zat)
         }#end if
         key.axis.options = list(side=ifelse(key.vertical,4,1),las=1,at=zat,labels=zlabels)
      }else{
         key.axis.options = as.list(key.axis.options)
         if (! "side" %in% names(y.axis.options)){
            key.axis.options = modifyList( x   = key.axis.options
                                         , val = list(side=ifelse(key.vertical,4,1))
                                         )#end modifyList
         }#end if (! "side" %in% names(y.axis.options))
      }#end if
      if (! is.null(key.title )) key.title  = as.list(key.title )
      if (! is.null(main.title)) main.title = as.list(main.title)
      #------------------------------------------------------------------------------------#



      #----- Save the margins to avoid losing the data. -----------------------------------#
      par.orig = par(no.readonly=TRUE )
      mar.orig = par.orig$mar
      on.exit(par(par.orig))
      par(par.user)
      #------------------------------------------------------------------------------------#




      #----- Check for outer margins. -----------------------------------------------------#
      if ( (! is.null(oma)) && (! is.null(omd))){
         stop ("You cannot provide both oma and omd!")
      }else if (is.null(oma) && is.null(omd)){
         par(oma=c(0,0,0,0))
      }else if (is.null(omd)){
         par(oma=oma)
      }else{
         par(omd=omd)
      }#end if
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Find offset for x axis label and right, based on legends and keys and outer   #
      # margin .                                                                           #
      #------------------------------------------------------------------------------------#
      par.tout = par(no.readonly=FALSE)
      #----- Bottom margin. ---------------------------------------------------------------#
      if (is.null(off.xlab)){
         if (plot.legend && key.vertical){
            off.xlab = with(par.tout,( omi[1] + f.leg * (din[2]-omi[1]-omi[3]) ) / din[2])
         }else if (key.vertical){
            off.xlab = with(par.tout,omi[1] / din[2])
         }else{
            off.xlab = with(par.tout,( omi[1] + f.key * (din[2]-omi[1]-omi[3]) ) / din[2])
         }#end if
      }#end if
      #----- Right margin. ----------------------------------------------------------------#
      if (is.null(off.right)){
         if (key.vertical){
            off.right = with(par.tout,( omi[4] + f.key * (din[1]-omi[2]-omi[4]) ) / din[1])
         }else if (plot.legend){
            off.right = with(par.tout,( omi[4] + f.leg * (din[1]-omi[2]-omi[4]) ) / din[1])
         }else{
            off.right = with(par.tout,omi[4] / din[1])
         }#end if
      }#end if
      #------------------------------------------------------------------------------------#



      #----- Split the screen into multiple pieces (legend, key, plots...) ----------------#
      fh.panel = 1. - f.key
      fv.panel = 1. - f.leg
      if (plot.legend && plot.key){
         layout( mat     = rbind(c(3, 2),c(1,0))
               , heights = c(fv.panel,f.leg)
               , widths  = c(fh.panel,f.key)
               )#end layout
      }else if (plot.key){
         if (key.vertical){
            layout(mat=cbind(2, 1), widths  = c(fh.panel,f.key))
         }else{
            layout(mat=rbind(2, 1), heights = c(fh.panel,f.key))
         }#end (if key.vertical)
      }else if (plot.legend){
         layout(mat=rbind(2,1), heights=c(fv.panel,f.leg))
      }#end if (plot.legend)
      #------------------------------------------------------------------------------------#






      #====================================================================================#
      #====================================================================================#
      #      First plot: the legend.                                                       #
      #------------------------------------------------------------------------------------#
      if (plot.legend){
         par(mar = c(0.1,0.1,0.1,0.1))
         plot.new()
         plot.window(xlim=c(0,1),ylim=c(0,1))
         do.call(what="legend",args=legend.options)
      }#end if
      #====================================================================================#
      #====================================================================================#






      #====================================================================================#
      #====================================================================================#
      #      Second plot: the key scale.                                                   #
      #------------------------------------------------------------------------------------#
      if (plot.key){
         par(mar = mar.key)
         plot.new()
         #---------------------------------------------------------------------------------#
         #     Plot in the horizontal or vertical depending on where the scale is going to #
         # be plotted.                                                                     #
         #---------------------------------------------------------------------------------#
         if (key.vertical){
            #----- Decide whether the scale is logarithmic or not. ------------------------#
            if (key.log){
               plot.window(xlim=c(0,1),ylim=range(clevels),xaxs="i",yaxs="i",log="y")
            }else{
               plot.window(xlim=c(0,1),ylim=range(clevels),xaxs="i",yaxs="i")
            }#end if
            #------------------------------------------------------------------------------#

            #----- Draw the colour bar. ---------------------------------------------------#
            rect( xleft   = 0
                , ybottom = clevels[-length(clevels)]
                , xright  = 1
                , ytop    = clevels[-1]
                , col     = ccolours
                , border  = ccolours
                )#end rect
            #------------------------------------------------------------------------------#
         }else{
            #----- Decide whether the scale is logarithmic or not. ------------------------#
            if (key.log){
               plot.window(xlim=range(clevels),ylim=c(0,1),xaxs="i",yaxs="i",las=1,log="x")
            }else{
               plot.window(xlim=range(clevels),ylim=c(0,1),xaxs="i",yaxs="i",las=1)
            }#end if
            #------------------------------------------------------------------------------#


            #----- Draw the colour bar. ---------------------------------------------------#
            rect( xleft   = clevels[-length(clevels)]
                , ybottom = 0
                , xright  = clevels[-1]
                , ytop    = 1
                , col     = ccolours
                , border  = ccolours
                )#end rect
            #------------------------------------------------------------------------------#
         }#end if
         #---------------------------------------------------------------------------------#


         #----- Plot the key axis. --------------------------------------------------------#
         do.call (what="axis",args=key.axis.options)
         #---------------------------------------------------------------------------------#


         #----- Draw box. -----------------------------------------------------------------#
         box()
         #---------------------------------------------------------------------------------#


         #----- Plot the title. -----------------------------------------------------------#
         if (! is.null(key.title)) do.call(what="title",args=key.title)
         #---------------------------------------------------------------------------------#
      }#end if (plot.key)
      #====================================================================================#
      #====================================================================================#



      #====================================================================================#
      #====================================================================================#
      #     Plot the main panel.                                                           #
      #------------------------------------------------------------------------------------#
      plog = paste0(ifelse(xlog,"x",""),ifelse(ylog,"y",""))
      par(mar = mar.main)
      plot.new()
      plot.window(xlim=xlim,ylim=ylim,log=plog,xaxs=xaxs,yaxs=yaxs,...)
      zupr  = zlim[1] + (1.-sqrt(.Machine$double.eps))*diff(zlim)
      zdens = pmin(zupr,zdens) + ifelse(zdens %gt% 0,0,NA) + 0. * zdens
      xyz   = list(x=xdens,y=ydens,z=zdens)
      image(x=xyz,zlim=zlim,col=ccolours,breaks=clevels,add=TRUE,useRaster=useRaster)
      #====================================================================================#
      #====================================================================================#
   }else{
      #====================================================================================#
      #====================================================================================#
      #     Plot the main panel on existing box.                                           #
      #------------------------------------------------------------------------------------#
      zupr  = zlim[1] + (1.-sqrt(.Machine$double.eps))*diff(zlim)
      zdens = pmin(zupr,zdens) + ifelse(zdens %gt% 0,0,NA) + 0. * zdens
      xyz   = list(x=xdens,y=ydens,z=zdens)
      image(x=xyz,zlim=zlim,col=ccolours,breaks=clevels,add=TRUE,useRaster=useRaster)
      #====================================================================================#
      #====================================================================================#
   }#end if (! add)
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Plot other options.  Check use a shared list, or one list for each sub-plot.      #
   #---------------------------------------------------------------------------------------#
   n.after = length(plot.after)
   for (a in sequence(n.after)){
      a.fun = names(plot.after)[a]
      a.args = plot.after[[a]]
      if (a.fun %in% "text" && reparse) a.args$labels = parse(text=a.args$labels)
      do.call(what=a.fun,args=a.args)
   }#end for
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Plot axes and title annotation.                                                   #
   #---------------------------------------------------------------------------------------#
   if (! add){
      do.call(what="axis" ,args=x.axis.options)
      do.call(what="axis" ,args=y.axis.options)
      do.call(what="title",args=main.title    )
      #------------------------------------------------------------------------------------#

      #----- Lastly, add the box (so it stays on top). ------------------------------------#
      box()
      #------------------------------------------------------------------------------------#
   }#end if (! add)
   #---------------------------------------------------------------------------------------#


   invisible()
}#end function xy.density
#==========================================================================================#
#==========================================================================================#
