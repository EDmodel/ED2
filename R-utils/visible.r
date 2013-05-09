#------------------------------------------------------------------------------------------#
#   Function that creates a nice colour scheme.                                            #
#------------------------------------------------------------------------------------------#
visible <<- function(n,sat=0.92){
  rrr   = c(  64, 128,  64,   0,   0,   0,  64, 128, 255, 255, 255, 128)*sat #- Red -------#
  ggg   = c(   0,   0,   0, 128, 255, 255, 255, 255, 255, 128,   0,   0)*sat #- Green -----#
  bbb   = c( 128, 255, 255, 255, 255, 128,   0,   0,   0,   0,   0,   0)*sat #- Blue ------#
  pivot = round(seq(from=1,to=n,by=(n-1)/(length(rrr)-1)),digits=0)

  red   = as.integer(spline(x=pivot,y=rrr,n=n)$y)
  green = as.integer(spline(x=pivot,y=ggg,n=n)$y)
  blue  = as.integer(spline(x=pivot,y=bbb,n=n)$y)

  red  [red   > 255] = 255; red  [red   < 0] = 0
  green[green > 255] = 255; green[green < 0] = 0
  blue [blue  > 255] = 255; blue [blue  < 0] = 0
  mycolsch = rgb(r=red,g=green,b=blue,maxColorValue=255)
  return(mycolsch)
}
#------------------------------------------------------------------------------------------#






#------------------------------------------------------------------------------------------#
#   Function that creates a nice colour scheme.                                            #
#------------------------------------------------------------------------------------------#
ivisible <<- function(n,sat=1.0){
  mycolsch = rev(visible(n,sat=sat))
  return(mycolsch)
}#end ivisible
#------------------------------------------------------------------------------------------#
