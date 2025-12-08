#------------------------------------------------------------------------------------------#
#   Function that re-creates QGIS "PuBuGn" colour scheme.                                  #
#------------------------------------------------------------------------------------------#
pubugn <<- function(n,sat=1.00){
   rrr   = c( 246, 214, 177, 128,  81,  38,  16,   1)*sat #- Red --------------------------#
   ggg   = c( 239, 218, 197, 178, 162, 148, 129, 108)*sat #- Green ------------------------#
   bbb   = c( 247, 235, 223, 212, 192, 161, 126,  89)*sat #- Blue -------------------------#
   pivot = round(seq(from=1,to=n,by=(n-1)/(length(rrr)-1)),digits=0)

   red   = pmax(0,pmin(255,as.integer(spline(x=pivot,y=rrr,n=n)$y)))
   green = pmax(0,pmin(255,as.integer(spline(x=pivot,y=ggg,n=n)$y)))
   blue  = pmax(0,pmin(255,as.integer(spline(x=pivot,y=bbb,n=n)$y)))

   mycolsch = rgb(r=red,g=green,b=blue,maxColorValue=255)
   return(mycolsch)
}#end function pubugn
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#   Function that re-creates QGIS "PuBuGn" colour scheme, but in inverse order.            #
#------------------------------------------------------------------------------------------#
ipubugn <<- function(n,sat=1.00){
   mycolsch = rev(pubugn(n,sat=sat))
   return(mycolsch)
}#end function ipubugn
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#   Function that re-creates QGIS "BuPu" colour scheme.                                    #
#------------------------------------------------------------------------------------------#
bupu <<- function(n,sat=1.00){
   rrr   = c( 237, 204, 174, 151, 139, 137, 133, 129)*sat #- Red --------------------------#
   ggg   = c( 248, 224, 197, 166, 132,  95,  55,  15)*sat #- Green ------------------------#
   bbb   = c( 251, 238, 223, 207, 189, 172, 149, 124)*sat #- Blue -------------------------#
   pivot = round(seq(from=1,to=n,by=(n-1)/(length(rrr)-1)),digits=0)

   red   = pmax(0,pmin(255,as.integer(spline(x=pivot,y=rrr,n=n)$y)))
   green = pmax(0,pmin(255,as.integer(spline(x=pivot,y=ggg,n=n)$y)))
   blue  = pmax(0,pmin(255,as.integer(spline(x=pivot,y=bbb,n=n)$y)))

   mycolsch = rgb(r=red,g=green,b=blue,maxColorValue=255)
   return(mycolsch)
}#end function pubugn
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#   Function that re-creates QGIS "BuPu" colour scheme, but in inverse order.              #
#------------------------------------------------------------------------------------------#
ibupu <<- function(n,sat=1.00){
   mycolsch = rev(bupu(n,sat=sat))
   return(mycolsch)
}#end function ipubugn
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#   Function that re-creates QGIS "YlGnBu" colour scheme.                                  #
#------------------------------------------------------------------------------------------#
ylgnbu <<- function(n,sat=1.00){
   rrr   = c( 255, 202, 147,  92,  59,  47,  41,  37)*sat #- Red --------------------------#
   ggg   = c( 255, 234, 213, 193, 166, 135,  95,  52)*sat #- Green ------------------------#
   bbb   = c( 204, 191, 182, 192, 193, 186, 169, 148)*sat #- Blue -------------------------#
   pivot = round(seq(from=1,to=n,by=(n-1)/(length(rrr)-1)),digits=0)

   red   = pmax(0,pmin(255,as.integer(spline(x=pivot,y=rrr,n=n)$y)))
   green = pmax(0,pmin(255,as.integer(spline(x=pivot,y=ggg,n=n)$y)))
   blue  = pmax(0,pmin(255,as.integer(spline(x=pivot,y=bbb,n=n)$y)))

   mycolsch = rgb(r=red,g=green,b=blue,maxColorValue=255)
   return(mycolsch)
}#end function ylgnbu
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#   Function that re-creates QGIS "YlGnBu" colour scheme, but in inverse order.            #
#------------------------------------------------------------------------------------------#
iylgnbu <<- function(n,sat=1.00){
   mycolsch = rev(ylgnbu(n,sat=sat))
   return(mycolsch)
}#end function iylgnbu
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#   Function that re-creates QGIS "BuGn" colour scheme.                                    #
#------------------------------------------------------------------------------------------#
bugn <<- function(n,sat=1.00){
   rrr   = c( 194, 122,  99,  51,  17,   0)*sat #- Red ------------------------------------#
   ggg   = c( 210, 177, 166, 128,  89,  51)*sat #- Green ----------------------------------#
   bbb   = c( 242, 204, 166, 102,  41,   0)*sat #- Blue -----------------------------------#
   pivot = round(seq(from=1,to=n,by=(n-1)/(length(rrr)-1)),digits=0)

   red   = pmax(0,pmin(255,as.integer(spline(x=pivot,y=rrr,n=n)$y)))
   green = pmax(0,pmin(255,as.integer(spline(x=pivot,y=ggg,n=n)$y)))
   blue  = pmax(0,pmin(255,as.integer(spline(x=pivot,y=bbb,n=n)$y)))

   mycolsch = rgb(r=red,g=green,b=blue,maxColorValue=255)
   return(mycolsch)
}#end function bugn
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#   Function that re-creates QGIS "BuGn" colour scheme, but in inverse order.              #
#------------------------------------------------------------------------------------------#
ibugn <<- function(n,sat=1.00){
   mycolsch = rev(bugn(n,sat=sat))
   return(mycolsch)
}#end function ibugn
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#   Function that re-creates QGIS "OrRd" colour scheme.                                    #
#------------------------------------------------------------------------------------------#
orrd <<- function(n,sat=1.00){
   rrr   = c( 254, 254, 242, 204, 166, 128)*sat #- Red ------------------------------------#
   ggg   = c( 250, 231, 187, 117,  52,   0)*sat #- Green ----------------------------------#
   bbb   = c( 228, 183, 131,  73,  30,   0)*sat #- Blue -----------------------------------#
   pivot = round(seq(from=1,to=n,by=(n-1)/(length(rrr)-1)),digits=0)

   red   = pmax(0,pmin(255,as.integer(spline(x=pivot,y=rrr,n=n)$y)))
   green = pmax(0,pmin(255,as.integer(spline(x=pivot,y=ggg,n=n)$y)))
   blue  = pmax(0,pmin(255,as.integer(spline(x=pivot,y=bbb,n=n)$y)))

   mycolsch = rgb(r=red,g=green,b=blue,maxColorValue=255)
   return(mycolsch)
}#end function orrd
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#   Function that re-creates QGIS "OrRd" colour scheme, but in inverse order.              #
#------------------------------------------------------------------------------------------#
iorrd <<- function(n,sat=1.00){
   mycolsch = rev(bugn(n,sat=sat))
   return(mycolsch)
}#end function iorrd
#------------------------------------------------------------------------------------------#
