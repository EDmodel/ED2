#------------------------------------------------------------------------------------------#
#   Function that re-creates QGIS "magma" colour scheme.                                   #
#------------------------------------------------------------------------------------------#
magma <<- function(n,sat=1.00){
   rrr   = c(   0,  26,  74, 121, 170, 217, 247, 254, 253, 252)*sat #- Red ----------------#
   ggg   = c(   0,  16,  16,  34,  51,  70, 114, 170, 226, 253)*sat #- Green --------------#
   bbb   = c(   4,  66, 121, 130, 125, 107,  92, 116, 163, 191)*sat #- Blue ---------------#
   pivot = round(seq(from=1,to=n,by=(n-1)/(length(rrr)-1)),digits=0)

   red   = pmax(0,pmin(255,as.integer(spline(x=pivot,y=rrr,n=n)$y)))
   green = pmax(0,pmin(255,as.integer(spline(x=pivot,y=ggg,n=n)$y)))
   blue  = pmax(0,pmin(255,as.integer(spline(x=pivot,y=bbb,n=n)$y)))

   mycolsch = rgb(r=red,g=green,b=blue,maxColorValue=255)
   return(mycolsch)
}#end function atlas
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#   Function that re-creates QGIS "magma" colour scheme, but in inverse order.             #
#------------------------------------------------------------------------------------------#
imagma <<- function(n,sat=1.00){
   mycolsch = rev(magma(n,sat=sat))
   return(mycolsch)
}#end function iatlas
#------------------------------------------------------------------------------------------#
