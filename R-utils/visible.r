#------------------------------------------------------------------------------------------#
#   Function that creates a nice colour scheme.                                            #
#------------------------------------------------------------------------------------------#
visible <<- function(n,sat=1.00){
  rrr   = c(  80,  40,  70, 230, 215, 178)*sat #- Red -------#
  ggg   = c(  30, 100, 195, 235, 125,   0)*sat #- Green -----#
  bbb   = c( 160, 250,  40,   0,   0,   0)*sat #- Blue ------#
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
