#------------------------------------------------------------------------------------------#
#   Function that creates a white to green colour scheme.                                  #
#------------------------------------------------------------------------------------------#
greenbow <<- function(n){
   nodes     = c("#FFFFFF","#EFFFCF","#C1E573","#89CC14","#4E7F0D","#254400")
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
#   Function that creates a green to white colour scheme.                                  #
#------------------------------------------------------------------------------------------#
igreenbow <<- function(n){
   rgb.out   = rev(greenbow(n=n))
   return(rgb.out)
}#end function igreenbow
#------------------------------------------------------------------------------------------#
