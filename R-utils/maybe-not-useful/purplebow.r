#------------------------------------------------------------------------------------------#
#   Function that creates a white to purple colour scheme.                                 #
#------------------------------------------------------------------------------------------#
purplebow <<- function(n){
   nodes     = c("#FFFFFF","#DEDEFF","#B7AFE3","#AA80FF","#8C41D8","#5A009A")
   nodes     = data.frame(t(col2rgb(nodes)))
   pivot     = round(seq(from=1,to=n,length.out=nrow(nodes)),digits=0)
   rgb.out   = data.frame(t(mapply(FUN=spline,y=nodes,MoreArgs=list(x=pivot,n=n))))$y
   rgb.out   = lapply(X=rgb.out,FUN=as.integer)
   rgb.out   = lapply(X=rgb.out,FUN=pmax,  0)
   rgb.out   = lapply(X=rgb.out,FUN=pmin,255)
   rgb.out   = rgb(r=rgb.out$red,g=rgb.out$green,b=rgb.out$blue,maxColorValue=255)
   return(rgb.out)
}#end function greenbow
#------------------------------------------------------------------------------------------#






#------------------------------------------------------------------------------------------#
#   Function that creates a purple to white colour scheme.                                 #
#------------------------------------------------------------------------------------------#
ipurplebow <<- function(n){
   rgb.out   = rev(purplebow(n=n))
   return(rgb.out)
}#end function igreenbow
#------------------------------------------------------------------------------------------#
