#------------------------------------------------------------------------------------------#
#   Function that re-creates QGIS "inferno" colour scheme.                                 #
#------------------------------------------------------------------------------------------#
inferno <<- function(n,sat=1.00){
   rrr   = c(  16,  40, 102, 161, 214, 248, 252, 252)*sat #- Red --------------------------#
   ggg   = c(   0,   6,  12,  38,  71, 125, 194, 255)*sat #- Green ------------------------#
   bbb   = c(  42,  85, 112,  99,  61,   0,   0, 159)*sat #- Blue -------------------------#
   pivot = round(seq(from=1,to=n,by=(n-1)/(length(rrr)-1)),digits=0)

   red   = pmax(0,pmin(255,as.integer(spline(x=pivot,y=rrr,n=n)$y)))
   green = pmax(0,pmin(255,as.integer(spline(x=pivot,y=ggg,n=n)$y)))
   blue  = pmax(0,pmin(255,as.integer(spline(x=pivot,y=bbb,n=n)$y)))

   mycolsch = rgb(r=red,g=green,b=blue,maxColorValue=255)
   return(mycolsch)
}#end function atlas
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#   Function that re-creates QGIS "inferno" colour scheme, but in inverse order.           #
#------------------------------------------------------------------------------------------#
iinferno <<- function(n,sat=1.00){
   mycolsch = rev(inferno(n,sat=sat))
   return(mycolsch)
}#end function iatlas
#------------------------------------------------------------------------------------------#
