#==========================================================================================#
#==========================================================================================#
#      This function plots stipples within a given set of coordinates.                     #
#------------------------------------------------------------------------------------------#
stipples <<- function(xleft,ybottom,xright,ytop,density=2,pch=".",cex=0.7,col="black",...){

   #----- Check that coordinates are all provided. ----------------------------------------#
   stopifnot(! ( missing(xleft) | missing(ybottom) | missing(xright) | missing(ytop)) )
   stopifnot(density %>% 0)
   #---------------------------------------------------------------------------------------#


   #----- Get arguments, and append points, cex, and col. ---------------------------------#
   dotdotdot = list(...)
   dotdotdot = modifyList(x=dotdotdot,val=list(pch=pch,cex=cex,col=col))
   #---------------------------------------------------------------------------------------#


   #----- Transform x and y in case this is a log plot. -----------------------------------#
   xlog = par("xlog")
   ylog = par("ylog")
   if (xlog){
      xleft  = log(xleft )
      xright = log(xright)
   }#end if
   if (ylog){
      ybottom = log(ybottom)
      ytop    = log(ytop   )
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- n is approximately how many points to plot. -------------------------------------#
   n          = round(sqrt(density)) + 1
   stipple.xy = function(a,z,n) mid.points(seq(from=a,to=z,length.out=n))
   xpts       = mapply( FUN      = stipple.xy
                      , a        = xleft
                      , z        = xright
                      , MoreArgs = list(n=n)
                      , SIMPLIFY = FALSE
                      )#end mapply
   ypts       = mapply( FUN      = stipple.xy
                      , a        = ybottom
                      , z        = ytop
                      , MoreArgs = list(n=n)
                      , SIMPLIFY = FALSE
                      )#end mapply
   #---------------------------------------------------------------------------------------#


   #----- Convert x and y back to scale in case of log plot. ------------------------------#
   if (xlog) xpts = lapply(X=xpts,exp)
   if (ylog) ypts = lapply(X=ypts,exp)
   xypts          = mapply(FUN=expand.grid,x=xpts,y=ypts,SIMPLIFY=TRUE)
   xypts          = apply(X=xypts,MARGIN=1,FUN=unlist)
   xypts          = as.data.frame(xypts)
   #---------------------------------------------------------------------------------------#


   #------ Plot points. -------------------------------------------------------------------#
   dotdotdot = modifyList(x=dotdotdot,val=list(x=xypts$x,y=xypts$y))
   do.call(what="points",args=dotdotdot)
   #---------------------------------------------------------------------------------------#

   invisible()
}#end function stipples
#==========================================================================================#
#==========================================================================================#
