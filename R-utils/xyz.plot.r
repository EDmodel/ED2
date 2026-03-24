#==========================================================================================#
#==========================================================================================#
# Function xyz.plot                                                                        #
#                                                                                          #
#    Given the x and y coordinates, this function will plot the xy scatter plot with the   #
# colour given by z...   This will generate as many plots as the number of lists in x, y,  #
# and z, and add a legend (if legend is not NULL), and a colour palette.                   #
#------------------------------------------------------------------------------------------#
xyz.plot <<- function( x
                   , y
                   , z
                   , fixed.xlim     = FALSE
                   , fixed.ylim     = FALSE
                   , xy.log         = ""
                   , xlim           = if (is.list(x) & ! fixed.xlim){
                                         lapply( X      = x
                                               , FUN    = pretty.xylim
                                               , is.log = regexpr("x",tolower(xy.log)) > 0
                                               )#end lapply
                                      }else{
                                         pretty.xylim( u      = unlist(x)
                                               , is.log = regexpr("x",tolower(xy.log)) > 0)
                                      }#end if
                   , ylim           = if (is.list(y) & ! fixed.ylim){
                                         lapply( X      = y
                                               , FUN    = pretty.xylim
                                               , is.log = regexpr("y",tolower(xy.log)) > 0
                                               )#end lapply
                                      }else{
                                         pretty.xylim( u      = unlist(y)
                                               , is.log = regexpr("y",tolower(xy.log)) > 0)
                                      }#end if
                   , zlim           = pretty.xylim(u=unlist(z),is.log=key.log)
                   , pch            = 15
                   , cex            = 1.0
                   , levels         = if (key.log){
                                         pretty.log(x=zlim,n=nlevels,forcelog=TRUE)
                                      }else{
                                         pretty(x=zlim,n=nlevels)
                                      }#end if
                   , nlevels        = 20
                   , colour.palette = cm.colors
                   , col            = colour.palette(length(levels)-1)
                   , na.col         = "grey94"
                   , xyz.title      = NULL
                   , xyz.sub        = if (length(x) == 1) {""} else {names(x)}
                   , xyz.legend     = NULL
                   , xyz.xaxis      = NULL
                   , xyz.yaxis      = NULL
                   , xyz.more       = NULL
                   , xyz.before     = xyz.more
                   , xyz.after      = NULL
                   , edge.axes      = FALSE
                   , key.title      = NULL
                   , key.log        = FALSE
                   , key.axis       = NULL
                   , key.width      = if (is.list(x)){
                                         if (length(x) > 1){7}else{5}
                                      }else{
                                         4.5
                                      }#end if
                   , leg.height     = if (is.list(x)){
                                         if (length(x) > 2){5}else{6}
                                      }else{
                                         6
                                      }#end if
                   , shuffle        = TRUE
                   , ...
                   ){



   #---------------------------------------------------------------------------------------#
   #      All three coordinates must be given.                                             #
   #---------------------------------------------------------------------------------------#
   if (missing(x) || missing(y) || missing(z)){
      cat (" X is missing: ",missing(x),"\n")
      cat (" Y is missing: ",missing(y),"\n")
      cat (" Z is missing: ",missing(z),"\n")
      stop("At least one of the data points is missing...")
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
   }else if (! is.list(x)){
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
         for (p in sequence(npanels)) pch[[p]] = orig.pch
      }#end if
      if (! is.list(cex)){
         orig.cex = cex
         cex      = list()
         for (p in sequence(npanels)) cex[[p]] = orig.cex
      }#end if
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Save the margins to avoid losing the data. --------------------------------------#
   par.orig = par(no.readonly = TRUE)
   mar.orig = par.orig$mar
   on.exit(par(par.orig))
   #---------------------------------------------------------------------------------------#



   #----- Check whether to add outer margins (we add it only if npanels > 1). -------------#
   par(oma = c(0.2,3,4.5,0)*(npanels > 1))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Split the screen into multiple blocks, plus one extra line for the legend and    #
   # one extra row for the colour bar.                                                     #
   #---------------------------------------------------------------------------------------#
   lo.box = pretty.box(npanels)
   if (is.null(xyz.legend)){
      emat   = cbind(lo.box$mat.off,rep(1,times=lo.box$nrow))
      layout( mat     = emat
            , heights = rep(1/lo.box$nrow,times=lo.box$nrow)
            , widths  = c(rep(key.width/lo.box$ncol,times=lo.box$ncol),1)
            )#end layout
      off.xlab  = 0
      off.right = 1/(key.width+1)
   }else{
      emat   = rbind( cbind(lo.box$mat.off2,rep(2,times=lo.box$nrow))
                    , c(rep(1,times=lo.box$ncol),0)
                    )#end rbind
      layout( mat     = emat
            , heights = c(rep(leg.height/lo.box$nrow,times=lo.box$nrow),1)
            , widths  = c(rep(key.width /lo.box$ncol,times=lo.box$ncol),1)
            )#end layout
      off.xlab  = 1/(leg.height+1)
      off.right = 1/(key.width +1)
   }#end if
   #---------------------------------------------------------------------------------------#
   #=======================================================================================#
   #=======================================================================================#





   #=======================================================================================#
   #=======================================================================================#
   #     If xyz.legend is not NULL, plot the legend first.                                 #
   #---------------------------------------------------------------------------------------#
   if (! is.null(xyz.legend)){
      par(mar=c(0.1,0.1,0.1,0.1))
      plot.new()
      plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
      do.call("legend",xyz.legend)
   }#end if
   #=======================================================================================#
   #=======================================================================================#





   #=======================================================================================#
   #=======================================================================================#
   #      Next plot (or first plot): the key scale.                                        #
   #---------------------------------------------------------------------------------------#
   par(mar = lo.box$mar.key)
   plot.new()
   #---------------------------------------------------------------------------------------#
   #     Plot in the horizontal or vertical depending on where the scale is going to       #
   # be plotted.                                                                           #
   #---------------------------------------------------------------------------------------#
   #----- Decide whether the scale is logarithmic or not. ---------------------------------#
   if (key.log){
      plot.window(xlim=c(0,1),ylim=range(levels),xaxs="i",yaxs="i",log="y")
   }else{
      plot.window(xlim=c(0,1),ylim=range(levels),xaxs="i",yaxs="i")
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Draw the colour bar. ------------------------------------------------------------#
   rect(xleft=0,ybottom=levels[-length(levels)],xright=1,ytop=levels[-1],col=col,border=col)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Check whether there are specific instructions for plotting the key axis.          #
   #---------------------------------------------------------------------------------------#
   if (is.null(key.axis)) {
      axis(side=4,las=1)
   }else{
      if (! "side" %in% names(key.axis)) key.axis$side = 4
      do.call("axis",key.axis)
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Draw box. -----------------------------------------------------------------------#
   box()
   #---------------------------------------------------------------------------------------#


   #----- Plot the title. -----------------------------------------------------------------#
   if (! is.null(key.title)) do.call("title",key.title)
   #---------------------------------------------------------------------------------------#

   #=======================================================================================#
   #=======================================================================================#





   #=======================================================================================#
   #=======================================================================================#
   #      Now we plot the other panels.                                                    #
   #---------------------------------------------------------------------------------------#
   for (p in sequence(npanels)){
      #------ Decide the margins based upon the XY axes limits. ---------------------------#
      if (edge.axes || (! fixed.xlim) || (! fixed.ylim)){
         left    = TRUE
         right   = TRUE
         top     = TRUE
         bottom  = TRUE
         mar.now = lo.box$mar0
      }else{
         left    = lo.panel$left  [p]
         right   = lo.panel$right [p]
         top     = lo.panel$top   [p]
         bottom  = lo.panel$bottom[p]
         mar.now = lo.panel$mar   [p,]
      }#end if (edge.axes || (! fixed.xlim) || (! fixed.ylim))
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Find out whether xlim and ylim are lists.                                      #
      #------------------------------------------------------------------------------------#
      xlim.now = if(is.list(xlim)) { xlim[[p]] }else{ xlim }
      ylim.now = if(is.list(ylim)) { ylim[[p]] }else{ ylim }
      #------------------------------------------------------------------------------------#



      #----- Set the window. --------------------------------------------------------------#
      par(mar = mar.now)
      plot.new()
      plot.window(xlim=xlim.now,ylim=ylim.now,log=xy.log,...)
      box()
      if (npanels != 1){
         title(main=xyz.sub[p],xlab="",ylab="",line=0.5)
      }#end if
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #    Split zleft into the breaks defined by the colour palette.                      #
      #------------------------------------------------------------------------------------#
      zcut              = cut(z[[p]],breaks=levels)
      zlev              = levels(zcut)
      zcol              = col[match(zcut,zlev)]
      zcol[is.na(zcol)] = na.col
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Check whether there are especial instructions for plotting the axes.           #
      #------------------------------------------------------------------------------------#
      if (is.null(xyz.xaxis) && bottom ){
         axis(side=1)
      }else if ( bottom ){
         if (! "side" %in% names(xyz.xaxis)) xyz.xaxis$side = 1
         do.call("axis",xyz.xaxis)
      }#end if
      if (is.null(xyz.yaxis) && left ){
         axis(side=2,las=1)
      }else if ( left ){
         if (! "side" %in% names(xyz.yaxis)) xyz.yaxis$side = 2
         do.call("axis",xyz.yaxis)
      }#end if
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Check whether there are additional instructions to plot. 
      #------------------------------------------------------------------------------------#
      if (! is.null(xyz.before)) {
         for (m in seq_along(xyz.before)){
            do.call(names(xyz.before)[m],xyz.before[[m]])
         }#end for
      }#end if
      #------------------------------------------------------------------------------------#



      #----- Call the function that actually plots the data. ------------------------------#
      if (shuffle){
         shf   = sample(x=length(x  [[p]]))
      }else{
         shf   = seq_along(x[[p]])
      }#end if
      shf.x    = pmin(shf,length(x  [[p]]))
      shf.y    = pmin(shf,length(y  [[p]]))
      shf.pch  = pmin(shf,length(pch[[p]]))
      shf.cex  = pmin(shf,length(cex[[p]]))
      shf.zcol = pmin(shf,length(zcol    ))
      points(x=x[[p]][shf.x],y=y[[p]][shf.y],pch=pch[[p]][shf.pch]
            ,cex=cex[[p]][shf.cex],col=zcol[shf.zcol],...)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Check whether there are additional instructions to plot. 
      #------------------------------------------------------------------------------------#
      if (! is.null(xyz.after)) {
         for (m in seq_along(xyz.after)){
            do.call(names(xyz.after)[m],xyz.after[[m]])
         }#end for
      }#end if
      #------------------------------------------------------------------------------------#
   }#end for
   #=======================================================================================#
   #=======================================================================================#



   #---------------------------------------------------------------------------------------#
   #     Plot the global title.                                                            #
   #---------------------------------------------------------------------------------------#
   if (! is.null(xyz.title)){
      #----- Make sure we get the main text. ----------------------------------------------#
      if (! is.list(xyz.title)){
         xyz.title=list(main=xyz.title)
      }else if (! "main" %in% names(xyz.title)){
         names(xyz.title)[[1]] = "main"
      }#end if
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #       Check whether to use title or gtitle.                                        #
      #------------------------------------------------------------------------------------#
      if (npanels == 1){
         xyz.title = modifyList(x=xyz.title,val=list(sub=xyz.sub[1]))
         do.call(what="title",args=xyz.title)
      }else{
         xyz.title = modifyList( x   = xyz.title
                               , val = list(off.xlab=off.xlab,off.right=off.right)
                               )#end modifyList
         do.call(what="gtitle",args=xyz.title)
      }#end if
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#

   invisible()
}#end function xyz.plot
#==========================================================================================#
#==========================================================================================#
