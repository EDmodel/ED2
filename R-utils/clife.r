#------------------------------------------------------------------------------------------#
#   Function that creates a purple to green colour scheme.                                 #
#------------------------------------------------------------------------------------------#
clife <<- function(n){
   purple = c("#DEE1FF","#7E76FF","#581BA2","#2A0053")
#   green  = c("#BFF684","#85C040","#408010","#274E08")
   green  = c("#C1E573","#89CC14","#4E7F0D","#152600")
   nodes  = c(rev(purple),green)
#   rrr       = c(  32,  96,  96, 212, 160,  32,   0)
#   ggg       = c(   0,   0, 128, 212, 255, 192,  48)
#   bbb       = c(  64, 255, 255, 212,   0,   0,   0)

#   rrr       = c(  60,  60, 100, 140, 180, 200, 150, 100,  50)
#   ggg       = c(   0,   0,  80, 140, 200, 240, 240, 170, 100)
#   bbb       = c( 120, 180, 240, 240, 180, 120,  60,   0,   0)
#   nodes     = mapply(FUN=rgb,red=rrr,green=ggg,blue=bbb,MoreArgs=list(maxColorValue=255))

#   nodes     = c("#3F1368","purple2","slateblue","lightslateblue","#C0ACCF"
#                ,"darkolivegreen1","olivedrab3","chartreuse2","forestgreen","#004000")
   nodes     = data.frame(t(col2rgb(nodes)))
   pivot     = round(seq(from=1,to=n,length.out=nrow(nodes)),digits=0)
   rgb.out   = data.frame(t(mapply(FUN=spline,y=nodes,MoreArgs=list(x=pivot,n=n))))$y
   rgb.out   = lapply(X=rgb.out,FUN=as.integer)
   rgb.out   = lapply(X=rgb.out,FUN=pmax,  0)
   rgb.out   = lapply(X=rgb.out,FUN=pmin,255)
   rgb.out   = rgb(r=rgb.out$red,g=rgb.out$green,b=rgb.out$blue,maxColorValue=255)
   return(rgb.out)
}#end function clife
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#   Function that creates a nice colour scheme.                                            #
#------------------------------------------------------------------------------------------#
iclife <<- function(n){
   rgb.out   = rev(clife(n=n))
   return(rgb.out)
}#end function iclife
#------------------------------------------------------------------------------------------#
