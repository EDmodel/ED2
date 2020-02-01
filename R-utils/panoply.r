#------------------------------------------------------------------------------------------#
#   Function that creates a nice colour scheme.                                            #
#------------------------------------------------------------------------------------------#
panoply <<- function(n,sat=1.00){
   rrr   = c(   0,   0,  70, 180, 230, 255, 200, 150)*sat #- Red --------------------------#
   ggg   = c(  50, 130, 180, 230, 230, 180,  90,   0)*sat #- Green ------------------------#
   bbb   = c( 100, 200, 255, 255, 180,  60,  10,   0)*sat #- Blue -------------------------#
   pivot = round(seq(from=1,to=n,by=(n-1)/(length(rrr)-1)),digits=0)

   red   = as.integer(spline(x=pivot,y=rrr,n=n)$y)
   green = as.integer(spline(x=pivot,y=ggg,n=n)$y)
   blue  = as.integer(spline(x=pivot,y=bbb,n=n)$y)

   red  [red   > 255] = 255; red  [red   < 0] = 0
   green[green > 255] = 255; green[green < 0] = 0
   blue [blue  > 255] = 255; blue [blue  < 0] = 0
   mycolsch = rgb(r=red,g=green,b=blue,maxColorValue=255)
   return(mycolsch)
}#end function panoply
#------------------------------------------------------------------------------------------#






#------------------------------------------------------------------------------------------#
#   Function that creates a nice colour scheme.                                            #
#------------------------------------------------------------------------------------------#
ipanoply <<- function(n,sat=1.0){
  mycolsch = rev(panoply(n,sat=sat))
  return(mycolsch)
}#end ivisible
#------------------------------------------------------------------------------------------#
