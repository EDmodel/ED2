#==========================================================================================#
#==========================================================================================#
#     This function plots a point cloud in 3-D, with colours representing the Z dimension. #
#------------------------------------------------------------------------------------------#
plot.ptcloud <<- function( pt.cloud
                         , col.variable     = c("z","intensity","pt.class")
                         , cex.intensity    = TRUE
                         , xlim             = NULL
                         , ylim             = NULL
                         , zlim             = NULL
                         , clim             = NULL
                         , xlog             = FALSE
                         , ylog             = FALSE
                         , zlog             = FALSE
                         , xaxt             = c("s","t","l","n")
                         , yaxt             = c("s","t","l","n")
                         , zaxt             = c("s","t","l","n")
                         , levels           = NULL
                         , nlevels          = 100
                         , colour.palette   = cm.colors
                         , col              = NULL
                         , na.col           = "grey94"
                         , key.log          = if(col.variable %in% "z"){zlog}else{FALSE}
                         , key.vertical     = TRUE
                         , key.axis.options = NULL
                         , key.options      = NULL
                         , plot.title       = NULL
                         , xlab             = NULL
                         , ylab             = NULL
                         , zlab             = NULL
                         , key.title        = NULL
                         , plot.after       = NULL
                         , f.key            = 1/9
                         , theta            = 315.
                         , phi              = 30.
                         , expand           = 0.5
                         , ticktype         = "detailed"
                         , shade            = 0.125
                         , ltheta           = -210.
                         , pch              = 16
                         , cex              = if(cex.intensity){c(0.1,0.8)}else{0.5}
                         , floor.col        = "grey94"
                         , plot.density     = TRUE
                         , plot.metrics     = TRUE
                         , zzdens           = NULL
                         , mhdens           = NULL
                         , n.dens           = 512
                         , from.dens        = NULL
                         , to.dens          = NULL
                         , skip.dens        = NULL
                         , plot.zmah        = TRUE
                         , plot.peaks       = TRUE
                         , col.dens         = c("grey10","grey60")
                         , lwd.dens         = c(2,2)
                         , lty.dens         = c("solid","solid")
                         , grid.dens        = TRUE
                         , pch.peaks        = c(4,3)
                         , lwd.peaks        = c(2,2)
                         , col.peaks        = c("midnightblue","deepskyblue")
                         , draw.box         = FALSE
                         , quant.metrics    = c(0.10,0.25,0.50,0.75,0.90,0.95)
                         , abline.metrics   = FALSE
                         , cex.metrics      = 1.5
                         , col.metrics      = c("midnightblue","deepskyblue")
                         , lty.metrics      = c("dotted")
                         , ...
                         ){



   #----- Standardise the colour and cex variables. ---------------------------------------#
   col.variable = match.arg(col.variable)
   if (col.variable %in% "pt.class") key.log = FALSE
   #---------------------------------------------------------------------------------------#



   #----- Check which kind of input was given. --------------------------------------------#
   if (missing(pt.cloud)) {
      #----- The user didn't provide anything. --------------------------------------------#
      cat  ("---------------------------------------------------------------------","\n")
      cat  ("  Variable pt.cloud must be provided!"                                ,"\n")
      cat  ("  This variable must be a data frame, a list, or a matrix."           ,"\n")
      cat  ("---------------------------------------------------------------------","\n")
      stop ("Point cloud missing")
      #------------------------------------------------------------------------------------#
   }else if (is.list(pt.cloud) || is.data.frame(pt.cloud) || is.matrix(pt.cloud)){

      #----- Turn matrix into a data frame. -----------------------------------------------#
      if (is.matrix(pt.cloud)){
         #----- Give names to the columns in case the matrix doesn't have names. ----------#
         if (is.null(colnames(pt.cloud))){
            colnames(pt.cloud) = c("x","y","z","intensity")[sequence(ncol(pt.cloud))]
         }#end if
         #---------------------------------------------------------------------------------#
         pt.cloud = as.data.frame(pt.cloud)
      }#end if
      #------------------------------------------------------------------------------------#


      #----- Make names case insensitive. -------------------------------------------------#
      names(pt.cloud) = tolower(names(pt.cloud))
      x = pt.cloud$x
      y = pt.cloud$y
      z = pt.cloud$z
      i = pt.cloud$intensity
      p = pt.cloud$pt.class
      #------------------------------------------------------------------------------------#




      #----- Check that x, y, and z were given. -------------------------------------------#
      if (is.null(x) || is.null(y) || is.null(z)){
         #----- The pt.cloud provided doesn't have all variables. -------------------------#
         cat  ("---------------------------------------------------------------------","\n")
         cat  ("  Variable pt.cloud must be provided!"                                ,"\n")
         cat  ("  Variable pt.cloud must contain variables x, y, and z."              ,"\n")
         cat  ("  In case you want to colour by intensity, provide intensity as well.","\n")
         cat  ("---------------------------------------------------------------------","\n")
         stop ("Coordinates are missing from point cloud")
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#




      #----- Intensity must be given in case it will be used to colour or scale points. ---#
      col.intensity = col.variable %in% "intensity"
      if (is.null(i) && (cex.intensity || plot.density)){
         cat  ("------------------------------------------------------------","\n")
         cat  ("     Variable pt.cloud does not have intensity."             ,"\n")
         cat  ("     Provide intensity as well, or set col.variable to z or" ,"\n")
         cat  (" pt.class, and cex.intensity and plot.density to FALSE."     ,"\n")
         cat  ("------------------------------------------------------------","\n")
         stop ("Intensity is missing from point cloud")
      }else if (is.null(i)){
         #----- Make a dummy intensity. ---------------------------------------------------#
         i = 0*x + 1
         #---------------------------------------------------------------------------------#
      }#end if (is.null(i))
      #------------------------------------------------------------------------------------#




      #----- pt.class must be given in case it will be used to colour or scale points. ----#
      col.pt.class = col.variable %in% "pt.class"
      if (is.null(p) && col.pt.class){
         cat  ("--------------------------------------------------------","\n")
         cat  ("     Variable pt.cloud does not have pt.class. "         ,"\n")
         cat  ("     Provide pt.class as well, or set col.variable to"   ,"\n")
         cat  (" z or intensity."                                        ,"\n")
         cat  ("--------------------------------------------------------","\n")
         stop ("Point class is missing from point cloud")
      }else if (is.null(i)){
         #----- Make a dummy intensity. ---------------------------------------------------#
         p = 0*x + 1
         #---------------------------------------------------------------------------------#
      }#end if (is.null(i))
      #------------------------------------------------------------------------------------#


      #----- Make sure the length of all four variables is the same. ----------------------#
      npts = unique(c(length(x),length(y),length(z),length(i),length(p)))
      if (length(npts) != 1){
         cat  ("--------------------------------------------------------","\n")
         cat  ("  LENGTH(x) = ",length(x)                                ,"\n")
         cat  ("  LENGTH(y) = ",length(y)                                ,"\n")
         cat  ("  LENGTH(z) = ",length(z)                                ,"\n")
         cat  ("  LENGTH(i) = ",length(i)                                ,"\n")
         cat  ("  LENGTH(p) = ",length(p)                                ,"\n")
         cat  ("  All variables above must have the same length."        ,"\n")
         cat  ("--------------------------------------------------------","\n")
         stop (" Variable pt.cloud has multiple lengths for elements.")
      }#end if
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Standardise xaxt, yaxt, and zaxt. -----------------------------------------------#
   xaxt = match.arg(xaxt)
   yaxt = match.arg(yaxt)
   zaxt = match.arg(zaxt)
   #---------------------------------------------------------------------------------------#


   #----- Use height to colour the points in case no field is given. ----------------------#
   o  = order(x=z,na.last=FALSE,decreasing=FALSE)
   x  = x [o]
   y  = y [o]
   z  = z [o]
   i  = i [o]
   p  = p [o]
   p  = ifelse(p %in% asprs.val,p,1)
   if (col.intensity){cc = i}else{cc = z}
   #---------------------------------------------------------------------------------------#



   #----- cex must have length 2 or one value for each point in case cex.intensity=TRUE. --#
   if ( (! length(cex) %in% c(2,npts)) && cex.intensity){
      warning("cex.intensity will be ignored since length(x) is not 2")
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Set default limits for x, y, z.                                                   #
   #---------------------------------------------------------------------------------------#
   if (is.null(xlim)) xlim = range(x ,finite=TRUE)
   if (is.null(ylim)) ylim = range(y ,finite=TRUE)
   if (is.null(zlim)) zlim = range(z ,finite=TRUE)
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #    Split colour field into the breaks defined by the colour palette.                  #
   #---------------------------------------------------------------------------------------#
   if (col.pt.class){
      #----- List classes that were actually used. ----------------------------------------#
      class.show = sort(unique(p))
      class.leg  = asprs.leg[match(class.show,asprs.val)]
      col        = asprs.col[match(class.show,asprs.val)]
      levels     = c(0,seq_along(class.show))
      class.at   = mid.points(levels)
      #------------------------------------------------------------------------------------#

      #----- Set colour classes. ----------------------------------------------------------#
      ccol       = asprs.col[match(p,asprs.val)]
      #------------------------------------------------------------------------------------#
   }else{
      #----- Set colour palette for colour field. -----------------------------------------#
      if (is.null(clim)) clim = range(cc,finite=TRUE)
      if (is.null(levels) && key.log){
         levels = sort(unique(pretty.log(x=clim,n=nlevels,forcelog=TRUE)))
      }else if (is.null(levels)){
         levels = sort(unique(pretty    (x=clim,n=nlevels)))
      }#end if
      if (is.null(col)) col = colour.palette(length(levels)-1)
      #------------------------------------------------------------------------------------#


      #----- Split the colours according to the cc field. ---------------------------------#
      ccut              = cut(cc,breaks=levels)
      clev              = levels(ccut)
      ccol              = col[match(ccut,clev)]
      ccol[is.na(ccol)] = na.col
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #    Define the expansion factor.                                                       #
   #---------------------------------------------------------------------------------------#
   if (length(cex) == 1){
      cex = rep(x=cex,times=npts)
   }else if (length(cex) == 2){
      if (length(unique(i)) == 1 || (! cex.intensity)){
         cex = rep(x=mean(cex),times=npts)
      }else{
         #----- Scale the data using intensity. -------------------------------------------#
         imin = range(i,finite=TRUE)[1]
         imax = range(i,finite=TRUE)[2]
         cmin = min(cex)
         cmax = max(cex)
         #---------------------------------------------------------------------------------#


         #----- Find the scale, and bound the sizes in case of infinite values. -----------#
         cex             = cmin + (i - imin) * (cmax - cmin) / (imax - imin)
         cex[is.na(cex)] = 0.5 * ( cmin + cmax )
         cex             = pmax(cmin,pmin(cmax,cex))
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Save the margins to avoid losing the data. --------------------------------------#
   par.orig = par(no.readonly=TRUE )
   mar.orig = par.orig$mar
   on.exit(par(par.orig))
   par(par.user)
   #---------------------------------------------------------------------------------------#



   #----- Split window into two blocks. ---------------------------------------------------#
   if (plot.density){
      layout(cbind(3,2,1),widths=c(c(1/3,2/3)*(1-f.key),f.key))
   }else{
      layout(cbind(2,1),widths=c(1-f.key,f.key))
   }#end if
   #---------------------------------------------------------------------------------------#






   #=======================================================================================#
   #=======================================================================================#
   #      Firs plot: the key scale.                                                        #
   #---------------------------------------------------------------------------------------#
      if (key.vertical){
         par(mar = c(5.1,0.6,4.1,4.1))
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
         if (col.pt.class && is.null(key.axis.options)){
            key.now = list(side=4,las=1,at=class.at,labels=class.leg)
         }else if (is.null(key.axis.options)){
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
         if (col.pt.class && is.null(key.axis.options)){
            key.now = list(side=1,las=1,at=class.at,labels=class.leg)
         }else if (is.null(key.axis.options)){
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



      #------------------------------------------------------------------------------------#
      #     Plot the key title.                                                            #
      #------------------------------------------------------------------------------------#
      if (! is.null(key.title)){
         #----- Make sure we get the main text. -------------------------------------------#
         if (! is.list(key.title)){
            key.title=list(main=key.title)
         }else if (! "main" %in% names(key.title)){
            names(key.title)[[1]] = "main"
         }#end if
         do.call(what="title",args=key.title)
         #---------------------------------------------------------------------------------#
      }#end if (! is.null(key.title))
      #------------------------------------------------------------------------------------#
   #=======================================================================================#
   #=======================================================================================#





   #=======================================================================================#
   #=======================================================================================#
   #      Now we build the floor and grid information for the 3-D plot.                    #
   #---------------------------------------------------------------------------------------#
   if (prod(zlim) < 0.){
      floor3d = 0.
   }else if (zlim[1] >= 0.){
      floor3d = zlim[1]
   }else{
      floor3d = zlim[2]
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Define the grid information for the 3-D plot. -----------------------------------#
   xpretty = if(xlog){pretty.log(xlim,n=5)}else{pretty(xlim,n=5)}
   ypretty = if(ylog){pretty.log(ylim,n=5)}else{pretty(ylim,n=5)}
   zpretty = if(zlog){pretty.log(zlim,n=5)}else{pretty(zlim,n=5)}
   xat     = if(xlog){log(xpretty)}else{xpretty}
   yat     = if(ylog){log(ypretty)}else{ypretty}
   zat     = if(zlog){log(zpretty)}else{zpretty}
   xlabels = sprintf("%g",xpretty)
   ylabels = sprintf("%g",ypretty)
   zlabels = sprintf("%g",zpretty)
   xlimit = range(x=xat)
   ylimit = range(x=yat)
   zlimit = range(x=zat)
   xfloor  = seq(from=xlimit[1],to=xlimit[2],length.out=16)
   yfloor  = seq(from=ylimit[1],to=ylimit[2],length.out=16)
   zfloor  = matrix(floor3d,nrow=length(xfloor),ncol=length(yfloor))
   #---------------------------------------------------------------------------------------#
   #=======================================================================================#
   #=======================================================================================#




   #=======================================================================================#
   #=======================================================================================#
   #     Plot the 3-D plot.                                                                #
   #---------------------------------------------------------------------------------------#
   if (plot.density && all(c(xaxt,yaxt,zaxt) %in% "n")){
      par(mar=c(0.1,0.1,4.1,0.1))
   }else if (plot.density){
      par(mar=c(1.1,3.1,4.1,1.1))
   }else{
      par(mar=c(1.1,1.1,4.1,1.1))
   }#end if
   pout = perspx( x         = xfloor
                , y         = yfloor
                , z         = zfloor
                , xlim      = xlimit
                , ylim      = ylimit
                , zlim      = zlimit
                , theta     = theta
                , phi       = phi
                , col       = floor.col
                , expand    = expand
                , box       = draw.box
                , ticktype  = ticktype
                , border    = NA
                , shade     = shade
                , ltheta    = ltheta
                , cex.main  = 0.8*cex.ptsz
                , axes      = FALSE
                )#end perspx
   #---------------------------------------------------------------------------------------#



   #----- Add X axes (or not). ------------------------------------------------------------#
   if (! xaxt %in% "n"){
      paxis3d(edge="X--",pmat=pout,at=xat,cex=0.9*cex.ptsz,labels=xlabels)
      #----- Plot label only if xlab is not NULL and xaxt is not "n". ---------------------#
      if (! is.null(xlab)){
         mtext3d(edge="X--",pmat=pout,labels=xlab,cex=cex.ptsz,srt=theta+90)
      }#end if
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Add Y axes (or not). ------------------------------------------------------------#
   if (! yaxt %in% "n"){
      paxis3d(edge="Y--",pmat=pout,at=yat,cex=0.9*cex.ptsz,labels=ylabels)
      #----- Plot label only if xlab is not NULL and xaxt is not "n". ---------------------#
      if (! is.null(ylab)){
         mtext3d(edge="Y--",pmat=pout,labels=ylab,cex=cex.ptsz,srt=theta)
      }#end if
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Add Z axes (or not). ------------------------------------------------------------#
   if (! zaxt %in% "n"){
      paxis3d(edge="Z-+",pmat=pout,at=zat,cex=0.9*cex.ptsz,labels=zlabels)
      #----- Plot label only if xlab is not NULL and xaxt is not "n". ---------------------#
      if (! is.null(zlab)){
         mtext3d(edge="Z-+",pmat=pout,labels=zlab,cex=cex.ptsz,srt=-75)
      }#end if
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Add the point clouds. ----------------------------------------------------------------#
   points(trans3d(x=x,y=y,z=z,pmat=pout),type="p",pch=pch,cex=cex,col=ccol)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Plot other options.  This should be a list, one list element for each sub-plot.   #
   #---------------------------------------------------------------------------------------#
   n.after = length(plot.after)
   for (a in sequence(n.after)){
      plot.now  = plot.after[[a]]
      plot.name = names(plot.after[a])
      if (plot.name %in% c("points","lines","text")){
         #---- Retrieve the coordinates. --------------------------------------------------#
         if (! all(c("x","y","z") %in% names(plot.now))){
            xyz = plot.now[[1]]
            if ((is.list(xyz) || is.data.frame(xyz)) && length(xyz) >= 3){
               if (all(c("x","y","z") %in% names(xyz))){
                  x3d = xyz$x
                  y3d = xyz$y
                  z3d = xyz$z
               }else{
                  x3d = xyz[[1]]
                  y3d = xyz[[2]]
                  z3d = xyz[[3]]
               }#end if
            }else if (is.matrix(xyz) && ncol(xyz) >= 3){
               x3d = xyz[,1]
               y3d = xyz[,2]
               z3d = xyz[,3]
            }else{
               cat ("------------------------------------------------------","\n",sep="")
               cat (" In plot.after, command ",plot.name,":"                ,"\n",sep="")
               cat (" Missing arguments! x, y, and z must be provided!"     ,"\n",sep="")
               cat ("------------------------------------------------------","\n",sep="")
               stop(" Missing arguments")
            }#end if
         }else{
            x3d = plot.now$x
            y3d = plot.now$y
            z3d = plot.now$z
         }#end if
         #---------------------------------------------------------------------------------#


         #------ Transform the coordinates. -----------------------------------------------#
         xyz=trans3d(x=x3d,y=y3d,z=z3d,pmat=pout)
         plot.now = modifyList(x=plot.now,val=list(x=xyz$x,y=xyz$y,z=NULL))
         do.call(what=plot.name,args=plot.now)
         #---------------------------------------------------------------------------------#
      }else{
         do.call(what=plot.name,args=plot.now)
      }#end if
   }#end for (a in sequence(n.after))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Plot the global title.                                                            #
   #---------------------------------------------------------------------------------------#
   if (! is.null(plot.title)){
      #----- Make sure we get the main text. ----------------------------------------------#
      if (! is.list(plot.title)){
         plot.title=list(main=plot.title)
      }else if (! "main" %in% names(plot.title)){
         names(plot.title)[[1]] = "main"
      }#end if
      do.call(what="title",args=plot.title)
      #------------------------------------------------------------------------------------#
   }#end if (! is.null(plot.title))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Plot the density title in case plot.density is TRUE.                             #
   #---------------------------------------------------------------------------------------#
   if (plot.density){
      #----- Make sure plot parameters have length 2. -------------------------------------#
      if (length(col.dens   ) == 1) col.dens    = rep(x=col.dens   ,times=2)
      if (length(lwd.dens   ) == 1) lwd.dens    = rep(x=lwd.dens   ,times=2)
      if (length(pch.peaks  ) == 1) pch.peaks   = rep(x=pch.peaks  ,times=2)
      if (length(lwd.peaks  ) == 1) lwd.peaks   = rep(x=lwd.peaks  ,times=2)
      if (length(col.peaks  ) == 1) col.peaks   = rep(x=col.peaks  ,times=2)
      if (length(col.metrics) == 1) col.metrics = rep(x=col.metrics,times=2)
      #------------------------------------------------------------------------------------#


      #----- Find the density functions. --------------------------------------------------#
      if (is.null(to.dens  )) to.dens   = zlim[2]
      if (is.null(from.dens) && zlim[1] > 0){
         from.dens = zlim[1]
      }else if (is.null(from.dens)){
         from.dens = to.dens / n.dens
      }#end if
      if (is.null(skip.dens)) skip.dens = quantile(x=z,probs=0.10,na.rm=TRUE)
      if (is.null(zzdens) && zlog){
         zdens  = density(x=log(z),n=n.dens,from=log(from.dens),to=log(to.dens))
         xdens  = zdens$y
         ydens  = exp(zdens$x)
      }else if (is.null(zzdens)){
         zdens  = density(x=z,n=n.dens,from=from.dens,to=to.dens)
         xdens  = zdens$y
         ydens  = zdens$x
      }else if (zlog){
         zdens  = zzdens
         xdens  = zzdens$y
         ydens  = exp(zzdens$x)
      }else{
         zdens  = zzdens
         xdens  = zzdens$y
         ydens  = zzdens$x
      }#end if
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #    Apply the MacArthur and Horn (1969) correction to the profile.                  #
      #------------------------------------------------------------------------------------#
      if (plot.zmah && is.null(mhdens)){
         mhdens = macarthur.horn( pt.cloud = pt.cloud
                                , zl       = from.dens
                                , zh       = to.dens
                                , zo       = skip.dens
                                , nz       = n.dens
                                )#end macarthur.horn
         if (zlog){
            zi     = jitter( x      = sample( x       = mhdens$x
                                            , size    = sum(i)
                                            , replace =TRUE
                                            , prob    = mhdens$y
                                            )#end sample
                           , amount = 0.5 * mean(diff(mhdens$x))
                           )#end jitter
            zidens = density(x=log(zi),n=n.dens,from=log(from.dens),to=log(to.dens))
            xidens = zidens$y
            yidens = exp(zidens$x)
         }else{
            xidens = mhdens$y
            yidens = mhdens$x
         }#end if
      }else if (plot.zmah){
         xidens = mhdens$y
         yidens = mhdens$x
      }else{
         xidens = rep(NA,2)
         yidens = rep(NA,2)
      }#end if
      #------------------------------------------------------------------------------------#




      #----- Fix limits. ------------------------------------------------------------------#
      xdlim    = pretty.xylim(u=c(xdens,xidens),fracexp=0.0,is.log=FALSE)
      ydlim    = pretty.xylim(u=c(ydens,yidens),fracexp=0.2,is.log=zlog )
      xdat     = pretty(c(xdens,xidens),n=5)
      xdlabels = sprintf("%g",xdat)
      if (zlog){
         ydat  = pretty.log(c(ydens,yidens),n=5)
      }else{
         ydat  = pretty(c(ydens,yidens),n=5)
      }#end if
      ydlabels = sprintf("%g",ydat)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Find the peaks of the distribution.                                            #
      #------------------------------------------------------------------------------------#
      if (plot.peaks){
         pk.dens  = peaks(xdens )
         #---- Check whether we can find peaks for MacArthur-Horn correction. -------------#
         if (plot.zmah){
            pk.idens = peaks(xidens)
         }else{
            pk.idens = rep(FALSE,times=length(xidens))
         }#end if (plot.zmah)
         #---------------------------------------------------------------------------------#
      }else{
         pk.dens  = rep(FALSE,times=length(xdens ))
         pk.idens = rep(FALSE,times=length(xidens))
      }#end if (plot.peaks)
      #------------------------------------------------------------------------------------#



      #----- Open the plotting area.  Note that y is the height and x is the density. -----#
      par(mar = c(5.1,4.1,4.1,1.1))
      plot.new()
      plot.window(xlim=xdlim,ylim=ydlim,log=if(zlog){"y"}else{""})
      if (grid.dens) abline(h=ydat,v=xdat,col=grid.colour,lty="dotted",lwd=0.75)
      #------------------------------------------------------------------------------------#



      #----- Plot metrics. ----------------------------------------------------------------#
      if (plot.metrics){
        z.quant = quantile(x=z,probs=quant.metrics,na.rm=TRUE)
        z.idx   = mapply(FUN=which.closest,x=z.quant,MoreArgs=list(A=ydens))
        text( x      = xdens[z.idx]
            , y      = ydens[z.idx]
            , labels = round(100*quant.metrics,1)
            , col    = col.metrics[2]
            , font   = 2
            , cex    = 1.1
            )#end text
        
#        abline(h=z.quant,col=col.metrics[2],lty=lty.metrics,lwd=2)
      }#end if
      #------------------------------------------------------------------------------------#



      #----- Unscaled density function. ---------------------------------------------------#
      lines ( x    = xdens
            , y    = ydens
            , type = "l"
            , col  = col.dens[1]
            , lwd  = lwd.dens[1]
            , lty  = lty.dens[1]
            )#end lines
      points( x    = xdens[pk.dens]
            , y    = ydens[pk.dens]
            , type = "p"
            , col  = col.peaks[1]
            , lwd  = lwd.peaks[1]
            , pch  = pch.peaks
            )#end points
      #------------------------------------------------------------------------------------#



      #----- Scaled density function. -----------------------------------------------------#
      if (plot.zmah){
         lines ( x    = xidens
               , y    = yidens
               , type = "l"
               , col  = col.dens[2]
               , lwd  = lwd.dens[2]
               , lty  = lty.dens[2]
               )#end lines
         points( x    = xidens[pk.idens]
               , y    = yidens[pk.idens]
               , type = "p"
               , col  = col.peaks[2]
               , lwd  = lwd.peaks[2]
               , pch  = pch.peaks
               )#end points
      }#end if (plot.zmah)
      #------------------------------------------------------------------------------------#



      #----- Plot metrics. ----------------------------------------------------------------#
      if (plot.metrics){
        z.mean = mean(x=z,na.rm=TRUE)
        z.sdev = sd  (x=z,na.rm=TRUE)
        error.bar( x    = mean(xdlim)
                 , y    = z.mean
                 , yerr = z.sdev
                 , col  = col.metrics[2]
                 , pch  = 15
                 , cex  = 2
                 , add  = TRUE
                 )#end error.bar
        text     ( x      = xdlim[1] + c(0.4,0.5,0.5)*diff(xdlim)
                 , y      = c(z.mean,z.mean+c(-0.5,0.5)*z.sdev)
                 , labels = parse(text=c("mu[z]","sigma[z]"))
                 , col    = col.metrics[1]
                 )#end text
      }#end if
      #------------------------------------------------------------------------------------#



      #----- Plot legend. -----------------------------------------------------------------#
      if (plot.zmah){
         legend( x       = "topright"
               , inset   = 0.001
               , legend  = c("Raw","MacArthur-Horn")
               , col     = col.dens
               , lty     = lty.dens
               , lwd     = lwd.dens
               , bg      = background
               , cex     = 0.6
               , bty     = "n"
               )#end legend
      }#end if (plot.zmah)
      #------------------------------------------------------------------------------------#


      #----- Plot annotation. -------------------------------------------------------------#
      axis (side=1,las=1,at=xdat,labels=xdlabels)
      axis (side=2,las=1,at=ydat,labels=ydlabels)
      title(xlab=desc.unit(desc="Density function",unit=untab$empty))
      title(ylab=zlab)
      box()
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Nothing to return. --------------------------------------------------------------#
   invisible()
   #---------------------------------------------------------------------------------------#
}#end function plot.ptcloud
#==========================================================================================#
#==========================================================================================#
