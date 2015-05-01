#------------------------------------------------------------------------------------------#
#   Function that creates an atlas-like colour scheme.                                     #
#------------------------------------------------------------------------------------------#
atlas <<- function(n){
   green  = c("#D4FFCC","#73E69A","#36B260","#137F38","#004D1A")
   orange = c("#FFF2B3","#FFCC66","#CB7A51","#994317","#4C1900")

   nodes     = c(rev(green),orange)
   nodes     = data.frame(t(col2rgb(nodes)))
   pivot     = round(seq(from=1,to=n,length.out=nrow(nodes)),digits=0)
   rgb.out   = data.frame(t(mapply(FUN=spline,y=nodes,MoreArgs=list(x=pivot,n=n))))$y
   rgb.out   = lapply(X=rgb.out,FUN=as.integer)
   rgb.out   = lapply(X=rgb.out,FUN=pmax,  0)
   rgb.out   = lapply(X=rgb.out,FUN=pmin,255)
   rgb.out   = rgb(r=rgb.out$red,g=rgb.out$green,b=rgb.out$blue,maxColorValue=255)
   return(rgb.out)
}#end function atlas
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#   Function that creates an atlas-like colour scheme.                                     #
#------------------------------------------------------------------------------------------#
iatlas <<- function(n){
   rgb.out   = rev(atlas(n=n))
   return(rgb.out)
}#end function iatlas
#------------------------------------------------------------------------------------------#
