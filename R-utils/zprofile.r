#==========================================================================================#
#==========================================================================================#
# Function zprofile                                                                        #
#                                                                                          #
#    This function is based on filled.contour, but without the black lines between the     #
# colours in the scale... It also plots additional information at the bottom of the panel  #
#------------------------------------------------------------------------------------------#
zprofile = function(lon=seq(from=0,to=1,len=nrow(vari)),lev=seq(from=0,to=1,len=ncol(vari))
                   ,vari,lai=rep(0.,times=length(lon)),land=rep(0.5,times=length(lon))
                   ,limlon=range(lon,finite=TRUE),limlev=range(lev,finite=TRUE)
                   ,limvari=range(vari,finite=TRUE)
                   ,levels = if (key.log){
                                pretty.log(x=limvari,n=nlevels)
                             }else{
                                pretty(x=limvari,n=nlevels)
                             }#end if
                    ,nlevels=100,color.palette=cm.colors,col=color.palette(length(levels)-1)
                    ,plot.title,plot.axes,foot.axes,foot.title,key.title,key.axes
                    ,key.log=FALSE,asp=NA
                    ,xaxs="i",yaxs="i",las=0,axes=TRUE,frame.plot=axes
                    ,lty.foot="h",lwd.foot=1,col.foot="olivedrab",bg.foot="steelblue",...){

   #----- Check whether the vari matrix makes sense or not. -------------------------------#
   if (!is.matrix(vari) || nrow(vari) <= 1 || ncol(vari) <= 1){
      stop("no proper `vari' matrix specified")
   }else if(!is.double(vari)){
      storage.mode(vari) = "double"
   }#end if

   #----- No messed-up axes are allowed, they must increase. ------------------------------#
   if (any(diff(lon) <= 0) || any(diff(lev) <= 0)){
       stop("increasing lon and lev values expected")
   }#end if

   #----- Save the margins to avoid losing the data. --------------------------------------#
   mar.orig = (par.orig = par(c("mar", "las", "mfrow")))$mar
   mar.orig[2] = mar.orig[2] + 0.5
   on.exit(par(par.orig))

   #----- Split the screen into 2. --------------------------------------------------------#
   w = (3 + mar.orig[2]) * par("csi") * 2.54
   h = (5 + mar.orig[1]) * par("csi") * 2.54
   layout(mat=matrix(c(2,3,1,0), nc = 2), widths = c(1, lcm(w)),heights=c(1,lcm(h)))
   par(las = las)
   mar = mar.orig
   mar[4] = mar[2]
   mar[2] = 1
   mar[1] = 1./3.

   #----- First plot: the key scale. ------------------------------------------------------#
   par(mar = mar)
   plot.new()
   #----- Decide whether the scale is logarithmic or not. ---------------------------------#
   if (key.log){
      plot.window(xlim=c(0,1),ylim=range(levels),xaxs="i",yaxs="i",log="y")
   }else{
      plot.window(xlim=c(0,1),ylim=range(levels),xaxs="i",yaxs="i")
   }#end if

   #----- Draw the colour bar. ------------------------------------------------------------#
   rect(xleft=0,ybottom=levels[-length(levels)],xright=1,ytop=levels[-1],col=col
       ,border=col)

   #----- Check whether there are specific instructions for plotting the key axis. --------#
   if (missing(key.axes)) {
      if (axes) axis(side=4,...)
   }else{
      key.axes
   }#end if

   #----- Draw box. -----------------------------------------------------------------------#
   box()

   #----- Plot the title. -----------------------------------------------------------------#
   if (!missing(key.title)) key.title

   #----- Now we plot the filled contour. -------------------------------------------------#
   mar    = mar.orig
   mar[1] = 1./3.
   mar[4] = 1
   par(mar = mar)
   plot.new()
   plot.window(xlim=limlon,ylim=limlev,log="", xaxs = xaxs, yaxs = yaxs, asp = asp)

   #----- Call the function that actually plots the data. ---------------------------------#
   .Internal(filledcontour(as.double(lon), as.double(lev), vari, as.double(levels)
            ,col = col))

   #----- Check whether there are especial instructions for plotting the axes. ------------#
   if (missing(plot.axes)) {
       if (axes) {
           axis(2)
       }
   }else{
      plot.axes
   }#end if
   
   if (frame.plot) box()
   #----- Check whether there are especial instructions for plotting the title. -----------#
   if (missing(plot.title)){
       if (axes) title(main = "", xlab = "", ylab = "",...)
   }else{
       plot.title
   }#end if


   #----- Plot the footnote plot. ---------------------------------------------------------#
   mar    = mar.orig
   mar[4] = 1
   mar[3] = 1./3.
   par(mar = mar)
   plot.new()
   plot.window(xlim=limlon,ylim=c(0,max(lai,na.rm=TRUE)),log="",xaxs=xaxs, yaxs = yaxs
              ,asp = asp)
   if (missing(foot.axes)) {
       if (axes) {
           axis(1)
           axis(2)
       }
   }else{
      foot.axes
   }#end if
   points(x=lon,y=lai,type=lty.foot,lwd=lwd.foot,col=col.foot)
   points(x=lon,y=lai*(1-land),type="h",lwd=lwd.foot,col=bg.foot)

   #----- Check whether there are especial instructions for plotting the title. -----------#
   if (missing(foot.title)){
       if (axes) title(main = "", xlab = "", ylab = "",...)
   }else{
       foot.title
   }#end if

   if (frame.plot) box()

   invisible()
}#end function zprofile
#==========================================================================================#
#==========================================================================================#
