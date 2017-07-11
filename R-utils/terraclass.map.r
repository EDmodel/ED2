#==========================================================================================#
#==========================================================================================#
# Function terraclass.map                                                                  #
#                                                                                          #
#    Given the x and y coordinates, this function will plot the vegetation values using    #
# the simplified TerraClass categories along with a colour scheme (no interpolation).      #
#------------------------------------------------------------------------------------------#
terraclass.map <<- function( x
                           , y
                           , z
                           , h                = if (is.list(z) & length(z) > 1){
                                                   lapply( z
                                                         , function(x) rep(FALSE,length(x))
                                                         )#end lapply
                                                }else{
                                                   rep(FALSE,length(z))
                                                }#end if
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
                           , plot.highlight   = TRUE
                           , col.highlight    = "white"
                           , dens.highlight   = 2
                           , pch.highlight    = "."
                           , cex.highlight    = 0.7
                           , na.col           = "grey94"
                           , plot.key         = TRUE
                           , key.vertical     = TRUE
                           , x.axis.options   = NULL
                           , y.axis.options   = NULL
                           , key.options      = NULL
                           , sub.options      = NULL
                           , main.title       = NULL
                           , main.xlab        = NULL
                           , main.ylab        = NULL
                           , key.title        = NULL
                           , plot.after       = NULL
                           , matrix.plot      = FALSE
                           , edge.axes        = FALSE
                           , oma              = NULL
                           , f.key            = 1/6
                           , f.leg            = 1/6
                           , off.xlab         = NULL
                           , off.right        = NULL
                           , xaxs             = "i"
                           , yaxs             = "i"
                           , smidgen          = 0
                           , main.mar         = NULL
                           , key.mar          = NULL
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
   same.kind = (   is.list(x) == is.list(y)
                && is.list(x) == is.list(z) 
                && is.list(x) == is.list(h))
   if (! same.kind){
      cat(" X is list: ",is.list(x),"\n")
      cat(" Y is list: ",is.list(y),"\n")
      cat(" Z is list: ",is.list(z),"\n")
      cat(" H is list: ",is.list(h),"\n")
      stop ("X, Y, Z, and H must be of the same kind...")
   }else if (!is.list(x)){
      #----- Convert x, y, and z to lists. ------------------------------------------------#
      x              = list(x )
      y              = list(y )
      z              = list(z )
      h              = list(h )
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
   if (is.null(off.xlab )) off.xlab  = 0.0
   if (is.null(off.right)) off.right = f.key * plot.key
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



   #----- Make sure colours are defined for missing values. -------------------------------#
   z = lapply(X=z,FUN=tif2terraclass)
   #---------------------------------------------------------------------------------------#

   #----- Select only colours that appear in the map. -------------------------------------#
   tcnow.show = sort(unique(unlist(z)))
   tcnow.leg  = tc2012.leg[match(tcnow.show,tc2012.val)]
   tcnow.clss = tc2012.col[match(tcnow.show,tc2012.val)]
   tcnow.bks  = c(0,seq_along(tcnow.show))
   tcnow.at   = mid.points(tcnow.bks)
   #---------------------------------------------------------------------------------------#




   #----- Check for outer margins. --------------------------------------------------------#
   if (is.null(oma)) oma = c(0.2,3,4.5,0)*(npanels > 1)
   #---------------------------------------------------------------------------------------#



   #----- Save the margins to avoid losing the data. --------------------------------------#
   par.orig = par(no.readonly=TRUE)
   mar.orig = par.orig$mar
   on.exit(par(par.orig))
   par(par.user)
   par(oma = oma)
   #---------------------------------------------------------------------------------------#



   #----- Split the screen into multiple pieces (legend, key, plots...) -------------------#
   if (plot.key){
      fh.panel = 1. - f.key
      if ((! plot.key) && (npanels > 1)){
      }else if (npanels == 1 && key.vertical){
         layout(mat=cbind(2, 1), widths = c(fh.panel,f.key))
      }else if (matrix.plot && key.vertical){
         layout( mat     = cbind(lo.panel$mat.off,rep(1,times=lo.panel$nrow))
               , widths  = c(rep(fh.panel/lo.panel$ncol,times=lo.panel$ncol),f.key)
               )#end layout
      }else if (matrix.plot){
         layout( mat     = rbind(lo.panel$mat.off,rep(1,times=lo.panel$ncol))
               , heights = c(rep(fh.panel,times=lo.panel$nrow),f.key)
               )#end layout
      }else{
         layout( mat     =rbind(1+sequence(npanels),rep(1,times=npanels))
               , heights = c(fh.panel,f.key)
               )#end layout
      }#end if
   }else if (npanels > 1){
      layout(mat=lo.panel$mat)
   }#end if (plot.key)
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
   #      First plot: the key scale.  Plot in the horizontal or vertical depending on      #
   # where the scale is going to be plotted.                                               #
   #---------------------------------------------------------------------------------------#
   if (plot.key){
      if (key.vertical){
         if (is.null(key.mar)){
            par(mar = lo.panel$mar.key)
         }else{
            par(mar = key.mar)
         }#end if
         plot.new()

         #----- Set plot window. ----------------------------------------------------------#
         plot.window(xlim=c(0,1),ylim=range(tcnow.bks),xaxs="i",yaxs="i")
         #---------------------------------------------------------------------------------#



         #----- Draw the colour bar. ------------------------------------------------------#
         rect( xleft   = 0
             , ybottom = tcnow.bks[-length(tcnow.bks)]
             , xright  = 1
             , ytop    = tcnow.bks[-1]
             , col     = tcnow.clss
             ,border   = foreground
             )#end rect
         #---------------------------------------------------------------------------------#



         #----- Check whether there are specific instructions for plotting the key axis. --#
         key.now = list(side=4,at=tcnow.at,labels=tcnow.leg,las=1,...)
         do.call (what="axis",args=key.now)
         #---------------------------------------------------------------------------------#
      }else{
         if (is.null(key.mar)){
            par(mar = c(4.1,4.1,0.1,2.1))
         }else{
            par(mar = key.mar)
         }#end if
         plot.new()

         #----- Set plot window. ----------------------------------------------------------#
         plot.window(xlim=range(tcnow.bks),ylim=c(0,1),xaxs="i",yaxs="i",las=1)
         #---------------------------------------------------------------------------------#


         #----- Draw the colour bar. ------------------------------------------------------#
         rect( xleft   = tcnow.bks[-length(tcnow.bks)]
             , ybottom = 0
             , xright  = tcnow.bks[-1]
             , ytop    = 1
             , col     = tcnow.clss
             , border  = foreground
             )#end rect
         #---------------------------------------------------------------------------------#


         #----- Check whether there are specific instructions for plotting the key axis. --#
         key.now = list(side=1,las=1,at=tcnow.at,labels=tcnow.leg,...)
         do.call (what="axis",args=key.now)
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#


      #----- Draw box. --------------------------------------------------------------------#
      box()
      #------------------------------------------------------------------------------------#


      #----- Plot the title. --------------------------------------------------------------#
      if (!missing(key.title)) do.call(what="title",args=key.title)
      #------------------------------------------------------------------------------------#
   }#end if (plot.key)
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
      if (! is.null(main.mar)){
         mar.now = main.mar
      }else if (matrix.plot & edge.axes){
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
      par(mar = mar.now)
      plot.new()
      plot.window(xlim=xlim,ylim=ylim,xaxs=xaxs,yaxs=yaxs,...)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    Split zleft into the breaks defined by the colour palette.                      #
      #------------------------------------------------------------------------------------#
      zcol = tc2012.col[match(z[[p]],tc2012.val)]
      zcol[is.na(zcol)] = na.col
      #------------------------------------------------------------------------------------#



      #----- Find the corners for the rectangles. -----------------------------------------#
      nx      = length(x[[p]])
      ny      = length(y[[p]])
      xleft   = x[[p]] - 0.5 * (1. + smidgen) * dx[[p]]
      xright  = x[[p]] + 0.5 * (1. + smidgen) * dx[[p]]
      ybottom = y[[p]] - 0.5 * (1. + smidgen) * dy[[p]]
      ytop    = y[[p]] + 0.5 * (1. + smidgen) * dy[[p]]
      rect( xleft   = xleft
          , ybottom = ybottom
          , xright  = xright
          , ytop    = ytop
          , col     = zcol
          , border  = zcol
          , xpd     = FALSE
          )#end rect
      #------------------------------------------------------------------------------------#



      #----- Plot stippling points on highlighted areas. ----------------------------------#
      if (plot.highlight){
         stipples( xleft   = xleft  [h[[p]]]
                 , ybottom = ybottom[h[[p]]]
                 , xright  = xright [h[[p]]]
                 , ytop    = ytop   [h[[p]]]
                 , col     = col.highlight
                 , density = dens.highlight
                 , pch     = pch.highlight
                 , cex     = cex.highlight
                 )#end rect
      }#end if (plot.highlight)
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
}#end function terraclass.map
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      This function transforms the classes read in the GeoTIFF into indices for           #
# the simplified Terra Class categories.                                                   #
#------------------------------------------------------------------------------------------#
tif2terraclass <<- function(x){
   xo10 = round(x/10)
   xo10 = ifelse(xo10 %in% seq_along(tc2012.t2c),xo10,11)
   ans  = tc2012.t2c[xo10]
   return(ans)
}#end function tif2terraclass
#==========================================================================================#
#==========================================================================================#
