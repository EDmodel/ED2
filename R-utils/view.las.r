#==========================================================================================#
#==========================================================================================#
#     This function plots the point cloud in an interactive way, using the rgl package.    #
#                                                                                          #
# Function based on script provided by Yasir Kaheil.                                       #
# https://stat.ethz.ch/pipermail/r-help/2008-May/161202.html                               #
#------------------------------------------------------------------------------------------#
view.las <<- function( x
                     , off.x          = TRUE
                     , off.y          = TRUE
                     , off.z          = FALSE
                     , colour.palette = cm.colors
                     , nlevels        = 20
                     , new            = TRUE
                     , axes           = TRUE
                     , pt.cex         = 1.0
                     , nshow          = nrow(pt.cloud)
                     , fcol           = c("z","intensity","pt.class")
                     , ...
                     ){


   #---------------------------------------------------------------------------------------#
   #     Standardise field colour.                                                         #
   #---------------------------------------------------------------------------------------#
   fcol = match.arg(fcol)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Find out whether pt.cloud is a data frame or a character.  In case it is a        #
   # character, we read in the file.                                                       #
   #---------------------------------------------------------------------------------------#
   if (is.character(x)){
      ptc.file = x

      #----- Crash in case this file doesn't exist. ---------------------------------------#
      dummy = stopifnot(file.exists(ptc.file))
      #------------------------------------------------------------------------------------#


      #----- Find out whether the file is compressed or not. ------------------------------#
      pt.cloud = read.las(lasfile=ptc.file)
      #------------------------------------------------------------------------------------#
   }else if (is.data.frame(x)){
      pt.cloud = x
   }#end if (is.character(pt.cloud))
   #---------------------------------------------------------------------------------------#



   #------ Remove minimum values to change origin to minimum. -----------------------------#
   if (off.x) pt.cloud$x = pt.cloud$x - min(pt.cloud$x,na.rm=TRUE)
   if (off.y) pt.cloud$y = pt.cloud$y - min(pt.cloud$y,na.rm=TRUE)
   if (off.z) pt.cloud$z = pt.cloud$z - min(pt.cloud$z,na.rm=TRUE)
   #---------------------------------------------------------------------------------------#


   #------ Define variable scale. ---------------------------------------------------------#
   x.range = range(pt.cloud$x        ,finite=TRUE)
   y.range = range(pt.cloud$y        ,finite=TRUE)
   z.range = range(pt.cloud$z        ,finite=TRUE)
   i.range = range(pt.cloud$intensity,finite=TRUE)
   x.show  = (pt.cloud$x - x.range[1])/diff(x.range)
   y.show  = (pt.cloud$y - y.range[1])/diff(y.range)
   z.show  = (pt.cloud$z - z.range[1])/diff(z.range)
   i.show  = (pt.cloud$i - i.range[1])/diff(i.range)
   #---------------------------------------------------------------------------------------#


   #------ Find colour levels. ------------------------------------------------------------#
   if (fcol %in% "z"){
      vlevels  = pretty(x=c(0,1),n=nlevels)
      vclevels = colour.palette(length(vlevels)-1)
      vcut     = cut(z.show,breaks=nlevels)
      vcol     = vclevels[match(vcut,levels(vcut))]
   }else if (fcol %in% "intensity"){
      vlevels  = pretty(x=c(0,1),n=nlevels)
      vclevels = colour.palette(length(vlevels)-1)
      vcut     = cut(i.show,breaks=nlevels)
      vcol     = vclevels[match(vcut,levels(vcut))]
   }else if (fcol %in% "pt.class"){
      vval     = match(pt.cloud$pt.class,asprs.val)
      vcol     = asprs.col[vval]
   }#end if
   #---------------------------------------------------------------------------------------#


   #-----  Find out which points to show. -------------------------------------------------#
   idx = sample(nrow(pt.cloud),size=pmin(nshow,nrow(pt.cloud)),replace=FALSE)
   #---------------------------------------------------------------------------------------#

   #----- Start new image in case the user wants it. --------------------------------------#
   if (new) rgl.clear("all")
   #---------------------------------------------------------------------------------------#


   #------ Define axis labels. ------------------------------------------------------------#
   if (axes){
      xlabels = pretty(pt.cloud$x)
      ylabels = pretty(pt.cloud$y)
      zlabels = pretty(pt.cloud$z)
      xat     = (xlabels - x.range[1])/diff(x.range)
      yat     = (ylabels - y.range[1])/diff(y.range)
      zat     = (zlabels - z.range[1])/diff(z.range)
      rgl.lines(x=c(0,1.1),y=0       ,z=0       )
      rgl.lines(x=0       ,y=c(0,1.1),z=0       )
      rgl.lines(x=0       ,y=0       ,z=c(0,1.1)) 
      rgl.texts(x=xat  ,y=-0.05,z=-0.05,text=xlabels)
      rgl.texts(x=-0.05,y=yat  ,z=-0.05,text=ylabels)
      rgl.texts(x=-0.05,y=-0.05,z=zat  ,text=zlabels)
      rgl.texts( x    = c(  0.5,-0.15,-0.15)
               , y    = c(-0.15,  0.5,-0.15)
               , z    = c(-0.15,-0.15,  0.5)
               , text = c("x","y","z")
               )#end rgl.texts
   }#end if
   #---------------------------------------------------------------------------------------#

   rgl.points( x      = x.show[idx]
             , y      = y.show[idx]
             , z      = z.show[idx]
             , color  = vcol[idx]
             , lit    = FALSE
             , fog    = FALSE
             , size   = 3.0 *pt.cex
             ,...
              )#end rgl.spheres
   invisible()
}#end function view.las
#==========================================================================================#
#==========================================================================================#
