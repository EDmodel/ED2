#------------------------------------------------------------------------------------------#
#   Function that creates a nice colour scheme.                                            #
#------------------------------------------------------------------------------------------#
colsuccess = function(n){
  rrr   = c( 212, 177, 141, 106,  70,  35,   0) #---- Red pivots. -------------------------#
  ggg   = c( 255, 221, 186, 152, 117,  83,  48) #---- Green pivots. -----------------------#
  bbb   = c(   0,  32,  64,  48,  32,  16,   0) #---- Blue pivots. ------------------------#
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
icolsuccess = function(n){
  rrr   = c(   0,  35,  70, 106, 141, 177, 212) #---- Red pivots. -------------------------#
  ggg   = c(  48,  83, 117, 152, 186, 221, 255) #---- Green pivots. -----------------------#
  bbb   = c(   0,  16,  32,  48,  64,  32,   0) #---- Blue pivots. ------------------------#
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
