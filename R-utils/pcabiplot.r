#==========================================================================================#
#==========================================================================================#
#     This function finds the principal component analysis biplot.                         #
#------------------------------------------------------------------------------------------#
pcabiplot <<- function( pca
                      , nshow          = 20
                      , vnames         = NULL
                      , f.leg          = 1/6
                      , plot.vec.axes  = TRUE
                      , colour.palette = cm.colors
                      , pt.col         = "grey89"
                      , pt.pch         = 16
                      , pt.cex         = 0.5
                      , arr.length     = 0.10
                      , arr.lwd        = 2
                      , arr.lty        = "solid"
                      , cex.leg        = 0.8
                      , main           = ""
                      , ...
                      ){
   #----- Make sure x is a prcomp object. -------------------------------------------------#
   stopifnot(inherits(pca,"prcomp"))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Fix vnames in case it is NULL.                                                   #
   #---------------------------------------------------------------------------------------#
   if (is.null(vnames)){
      vnames = rownames(pca$rotation)
   }else if (length(vnames) != nrow(pca$rotation)){
      cat0("------------------------------------------------------------------------")
      cat0("    length(vnames):     ",length(vnames)    )
      cat0("    # of PCA variables: ",nrow(pca$rotation))
      cat0("------------------------------------------------------------------------")
      stop(" 'vnames' must have the same length as number of variables for PCA!")
   }#end if
   #---------------------------------------------------------------------------------------#




   #----- Load PCA, and find the strongest components. ------------------------------------#
   summ.pca = summary(pca)
   scores   = pca$x
   rota     = pca$rotation
   nvars    = ncol(scores)
   npca     = nrow(scores)
   lam      = pca$sdev * sqrt(npca)
   pca.pts  = t(t(scores) / lam)
   pca.vec  = t(t(rota  ) * lam)
   explain  = round(100*summ.pca$importance[2,],1)
   mag.vec  = apply(X=pca.vec,MARGIN=1,FUN=sum2)
   ord.vec  = order(mag.vec,decreasing=TRUE)
   nshow    = min(nshow,nvars)
   oidx     = ord.vec[sequence(nshow)]
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Create a colour palette to pick some colours.                                      #
   #---------------------------------------------------------------------------------------#
   sample.palette   = colour.palette(nshow)
   vars.colour      = colour.palette(nshow)
   pca.frmt         = ceiling(log10(nshow))
   pca.frmt         = paste("%",pca.frmt,".",pca.frmt,"i",sep="")
   vars.idx         = sprintf(pca.frmt,oidx)
   vars.labels      = paste0(vars.idx,": ",vnames[oidx])
   #---------------------------------------------------------------------------------------#




   #----- Split the panel. ----------------------------------------------------------------#
   par(par.user)
   par(oma=c(0,0,2,0))
   layout(mat= rbind(2,1),heights=c(1-f.leg,f.leg))
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     First, the legend.                                                                #
   #---------------------------------------------------------------------------------------#
   par(mar=c(0.1,0.1,0.1,0.1))
   plot.new()
   plot.window(xlim=c(0,1),ylim=c(0,1))
   legend( x      = "center"
         , inset  = 0.0
         , legend = vars.labels
         , fill   = vars.colour
         , bty    = "n"
         , xpd    = TRUE
         , cex    = cex.leg * cex.ptsz
         , ncol   = min(4,pretty.box(nshow)$ncol)
         , title  = ifelse( nshow < nvars
                          , paste(nshow," statistics with most variance")
                          , paste("Variables sorted by decreasing variance")
                          )#end ifelse
         )#end legend
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #   Make axes.                                                                          #
   #---------------------------------------------------------------------------------------#
   pvar.explained = paste("Component",sequence(npca)
                         ," - ", sprintf("%.1f",explain),"%",sep="")
   #---------------------------------------------------------------------------------------#



   #------ Find limits. -------------------------------------------------------------------#
   pxlimit  = max(abs(pca.pts[,1])) * c(-1,1)
   pylimit  = max(abs(pca.pts[,2])) * c(-1,1)
   vxlimit  = max(abs(pca.vec[,1])) * c(-1.20,1.20)
   vylimit  = max(abs(pca.vec[,2])) * c(-1.20,1.20)
   #---------------------------------------------------------------------------------------#


   #----- Plot the PCA points. ------------------------------------------------------------#
   par(mar=c(4.1,4.1,5.1,2.1))
   plot.new()
   plot.window(xlim=pxlimit,ylim=pylimit)
   abline(h=0,v=0,col=foreground,lty="solid",lwd=2)
   axis(side=1)
   axis(side=2)
   title( main = main
        , xlab = pvar.explained[1]
        , ylab = pvar.explained[2]
        , ...
        )#end title
   points(x=pca.pts[,1],y=pca.pts[,2],col=pt.col,cex=pt.cex,pch=pt.pch)
   #---------------------------------------------------------------------------------------#


   #----- Plot the PCA vectors. -----------------------------------------------------------#
   plot.window(xlim=vxlimit,ylim=vylimit)
   if (plot.vec.axes){
      axis(side=3,col.ticks=firebrick.fg,col.axis=firebrick.mg)
      axis(side=4,col.ticks=firebrick.fg,col.axis=firebrick.mg)
   }#end if
   box()
   for (u in sequence(nshow)){
      onow = oidx[u]
      arrows( x0     = 0
            , y0     = 0
            , x1     = pca.vec[onow,1]
            , y1     = pca.vec[onow,2]
            , col    = vars.colour[u]
            , length = arr.length
            , lwd    = arr.lwd
            , lty    = arr.lty
            )#end arrows
   }#end for
   for (u in sequence(nshow)){
      onow = oidx[u]
      rot  = 180*atan2(pca.vec[onow,2],pca.vec[onow,1])/pi
      if (rot >  90 & rot <=  180) rot = rot + 180
      if (rot < -90 & rot >= -180) rot = rot + 180
      mult = 1.00 + runif(n=1,min=0.040,max=0.125)
      text( x      = mult*pca.vec[onow,1]
          , y      = mult*pca.vec[onow,2]
          , labels = vars.idx[u]
          , col    = vars.colour[u]
          , cex    = 0.9
          , srt    = rot
          )#end text
   }#end for
   #---------------------------------------------------------------------------------------#


   #----- Discreet exit. ------------------------------------------------------------------#
   invisible()
   #---------------------------------------------------------------------------------------#
}#end pcabiplot
#==========================================================================================#
#==========================================================================================#
