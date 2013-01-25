#------------------------------------------------------------------------------------------#
#   Function that creates a nice colour scheme.                                            #
#------------------------------------------------------------------------------------------#
muitas <- function(n){
  rrr   <- c(   0,   0,   0, 212, 255, 192,  64) #---- Red pivots. ------------------------#
  ggg   <- c(   0,  64, 192, 212, 128,  32,   0) #---- Green pivots. ----------------------#
  bbb   <- c(  64, 255,  64,   0,   0,   0,   0) #---- Blue pivots. -----------------------#
  pivot <- round(seq(from=1,to=n,by=(n-1)/(length(rrr)-1)),digits=0)

  red   <- as.integer(spline(x=pivot,y=rrr,n=n)$y)
  green <- as.integer(spline(x=pivot,y=ggg,n=n)$y)
  blue  <- as.integer(spline(x=pivot,y=bbb,n=n)$y)

  red  [red   > 255] <- 255; red  [red   < 0] <- 0
  green[green > 255] <- 255; green[green < 0] <- 0
  blue [blue  > 255] <- 255; blue [blue  < 0] <- 0
  mycolsch <- rgb(r=red,g=green,b=blue,maxColorValue=255)
  return(mycolsch)
}
#------------------------------------------------------------------------------------------#






#------------------------------------------------------------------------------------------#
#   Function that creates a nice colour scheme.                                            #
#------------------------------------------------------------------------------------------#
imuitas <- function(n){
  rrr   <- c(  64, 192, 255, 212,   0,   0,   0) #---- Red pivots. ------------------------#
  ggg   <- c(   0,  32, 128, 212, 192,  64,   0) #---- Green pivots. ----------------------#
  bbb   <- c(   0,   0,   0,   0,  64, 255,  64) #---- Blue pivots. -----------------------#
  pivot <- round(seq(from=1,to=n,by=(n-1)/(length(rrr)-1)),digits=0)

  red   <- as.integer(spline(x=pivot,y=rrr,n=n)$y)
  green <- as.integer(spline(x=pivot,y=ggg,n=n)$y)
  blue  <- as.integer(spline(x=pivot,y=bbb,n=n)$y)

  red  [red   > 255] <- 255; red  [red   < 0] <- 0
  green[green > 255] <- 255; green[green < 0] <- 0
  blue [blue  > 255] <- 255; blue [blue  < 0] <- 0
  mycolsch <- rgb(r=red,g=green,b=blue,maxColorValue=255)
  return(mycolsch)
}
#------------------------------------------------------------------------------------------#
