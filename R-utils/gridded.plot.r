#==========================================================================================#
#==========================================================================================#
#     Function gridded.plot                                                                #
#                                                                                          #
#     This function plots a gridded plot, with the colour scheme given.  This is similar   #
# to sombreado, but structured more like image.map.  Currently only one panel is allowed,  #
# though this may change soon.                                                             #
#------------------------------------------------------------------------------------------#
gridded.plot <<- function( x                = seq(from=0,to=1,len=nrow(z))
                         , y                = seq(from=0,to=1,len=ncol(z))
                         , z
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
                         , main.title       = NULL
                         , key.title        = NULL
                         , plot.after       = NULL
                         , matrix.plot      = FALSE
                         , legend.options   = NULL
                         , edge.axes        = FALSE
                         , mar              = NULL
                         , oma              = NULL
                         , omd              = NULL
                         , f.key            = 1/6
                         , f.leg            = 1/6
                         , xaxs             = "i"
                         , yaxs             = "i"
                         , smidgen          = 0
                         , useRaster        = TRUE
                         , ...
                         ){

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
         #---------------------------------------------------------------------------------#
      }else{
         #----- Bad setting. -------------------------------------------------------------#
         stop("no `z' matrix specified")
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#
   }else if (is.list(x)) {
      #----- Z is there, just need to check whether x and y were given as a list... -------#
      y = x$y
      x = x$x
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#

   #----- Check whether the z matrix makes sense or not. ----------------------------------#
   if (! (is.matrix(z) || is.data.frame(z))){
      stop("no proper `z' matrix specified")
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- No messed-up axes are allowed, they must increase. ------------------------------#
   if (any(diff(x) %<=% 0) || any(diff(y) %<=% 0)){
       stop("increasing x and y values expected")
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


   #----- Save the margins to avoid losing the data. --------------------------------------#
   par.orig = par(no.readonly=TRUE)
   mar.orig = par.orig$mar
   on.exit(par(par.orig))
   par(par.user)
   #---------------------------------------------------------------------------------------#





   #----- Check for outer margins. --------------------------------------------------------#
   if ( (! is.null(oma)) && (! is.null(omd))){
      stop ("You cannot provide both oma and omd!")
   }else if (is.null(oma) && is.null(omd)){
      par(oma=c(0,0,0,0))
   }else if (is.null(omd)){
      par(oma=oma)
   }else{
      par(omd=omd)
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Split the screen into multiple pieces (legend, key, plots...) -------------------#
   fh.panel = 1. - f.key
   fv.panel = 1. - f.leg
   if (plot.legend){
      layout( mat     = rbind(c(3, 2),c(1,0))
            , heights = c(fv.panel,f.leg)
            , widths  = c(fh.panel,f.key)
            )#end layout
   }else if (key.vertical){
      layout(mat=cbind(2, 1), widths = c(fh.panel,f.key))
   }else{
      layout( mat     = rbind(2,1)
            , heights = c(fh.panel,f.key)
            )#end layout
   }#end if
   #---------------------------------------------------------------------------------------#
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
         par(mar = pretty.box(n=1)$mar.key)
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


   #----- Set the window. -----------------------------------------------------------------#
   if (! is.null(mar)){
      mar.now = mar
   }else if (key.vertical){
      mar.now = c(5.1,4.1,4.1,2.1)
   }else{
      mar.now = c(4.1,4.1,4.1,2.1)
   }#end if
   plog = ""
   if (xlog) plog = paste(plog,"x",sep="")
   if (ylog) plog = paste(plog,"y",sep="")
   par(mar = mar.now)
   plot.new()
   plot.window(xlim=xlim,ylim=ylim,log=plog,xaxs=xaxs,yaxs=yaxs,...)
   #---------------------------------------------------------------------------------------#


   #----- Plot the field. -----------------------------------------------------------------#
   xyz = list(x=x,y=y,z=z)
   image(xyz,breaks=levels,col=col,add=TRUE,useRaster=TRUE)
   #---------------------------------------------------------------------------------------#




   #---- Plot the X axis. -----------------------------------------------------------------#
   if (! is.null(x.axis.options)){
      x.axis.now = modifyList(x=x.axis.options,val=list(side=1))
   }else{
      x.axis.now = list(side=1,las=1)
   }#end if
   do.call(what="axis",args=x.axis.now)
   #---------------------------------------------------------------------------------------#




   #---- Plot the Y axis. -----------------------------------------------------------------#
   if (! is.null(y.axis.options)){
      y.axis.now = modifyList(x=y.axis.options,val=list(side=2))
   }else{
      y.axis.now = list(side=2,las=1)
   }#end if
   do.call(what="axis",args=y.axis.now)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Plot other options.                                                               #
   #---------------------------------------------------------------------------------------#
   n.after = length(plot.after)
   for (a in sequence(n.after)){
       do.call(what=names(plot.after)[a],args=plot.after[[a]])
   }#end for
   #---------------------------------------------------------------------------------------#


   #----- Lastly, add the box (so it stays on top). ---------------------------------------#
   box()
   #---------------------------------------------------------------------------------------#



   #----- Make sure we get the main text. -------------------------------------------------#
   if (! is.list(main.title)){
      main.title=list(main=main.title)
   }else if (! "main" %in% names(main.title)){
      names(main.title)[[1]] = "main"
   }#end if
   do.call(what="title",args=main.title)
   #---------------------------------------------------------------------------------------#


   #----- Don't bother the user with messages. --------------------------------------------#
   invisible()
   #---------------------------------------------------------------------------------------#
}#end function gridded.plot
#==========================================================================================#
#==========================================================================================#
