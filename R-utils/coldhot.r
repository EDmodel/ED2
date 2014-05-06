#------------------------------------------------------------------------------------------#
#   Function that creates a nice cold colour scheme.                                       #
#------------------------------------------------------------------------------------------#
hue.cold <<- function(n,sat=1.00){
  rrr   = c( 192, 132,  98,  46,   0,  74,  41)*sat #- Red --------------------------------#
  ggg   = c( 237, 221, 147, 111,  68,   0,   0)*sat #- Green ------------------------------#
  bbb   = c( 254, 253, 254, 253, 254, 230, 130)*sat #- Blue -------------------------------#
  pivot = round(seq(from=1,to=n,by=(n-1)/(length(rrr)-1)),digits=0)

  red   = as.integer(spline(x=pivot,y=rrr,n=n)$y)
  green = as.integer(spline(x=pivot,y=ggg,n=n)$y)
  blue  = as.integer(spline(x=pivot,y=bbb,n=n)$y)

  red  [red   > 255] = 255; red  [red   < 0] = 0
  green[green > 255] = 255; green[green < 0] = 0
  blue [blue  > 255] = 255; blue [blue  < 0] = 0
  mycolsch = rgb(r=red,g=green,b=blue,maxColorValue=255)
  return(mycolsch)
}#end function
#------------------------------------------------------------------------------------------#





#------------------------------------------------------------------------------------------#
#   Function that creates a nice hot colour scheme.                                        #
#------------------------------------------------------------------------------------------#
hue.hot <<- function(n,sat=1.00){
  rrr   = c( 255, 254, 253, 255, 217, 110,  70)*sat #- Red --------------------------------#
  ggg   = c( 252, 209, 165,  83,  11,   5,   0)*sat #- Green ------------------------------#
  bbb   = c( 171, 100,  49,   8,   0,   0,   0)*sat #- Blue -------------------------------#
  pivot = round(seq(from=1,to=n,by=(n-1)/(length(rrr)-1)),digits=0)

  red   = as.integer(spline(x=pivot,y=rrr,n=n)$y)
  green = as.integer(spline(x=pivot,y=ggg,n=n)$y)
  blue  = as.integer(spline(x=pivot,y=bbb,n=n)$y)

  red  [red   > 255] = 255; red  [red   < 0] = 0
  green[green > 255] = 255; green[green < 0] = 0
  blue [blue  > 255] = 255; blue [blue  < 0] = 0
  mycolsch = rgb(r=red,g=green,b=blue,maxColorValue=255)
  return(mycolsch)
}#end function
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#   Combined palettes.                                                                     #
#------------------------------------------------------------------------------------------#
coldhot   <<- function(n,sat=1.0){
    #--------------------------------------------------------------------------------------#
    n.each     = floor(n/2)
    cold       = ihue.cold(n=n.each)
    tepid      = rgb(r=227,g=227,b=227,maxColorValue=255)
    hot        = hue.hot (n=n.each)
    mycolsch   = c(cold,rep(tepid,times=(n%%2)),hot)
    return(mycolsch)
    #--------------------------------------------------------------------------------------#
}#end coldhot
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#   Inverse colour palettes.                                                               #
#------------------------------------------------------------------------------------------#
ihue.cold <<- function(n,sat=1.0) rev(hue.cold(n,sat=sat))
ihue.hot  <<- function(n,sat=1.0) rev(hue.hot (n,sat=sat))
hotcold   <<- function(n,sat=1.0) rev(coldhot (n,sat=sat))
#------------------------------------------------------------------------------------------#
