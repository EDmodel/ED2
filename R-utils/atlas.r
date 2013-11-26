#------------------------------------------------------------------------------------------#
#   Function that creates a nice colour scheme.                                            #
#------------------------------------------------------------------------------------------#
atlas <<- function(n){
  rrr   <- c(  18,  42, 149, 237, 160, 102,  55) #---- Red pivots. ------------------------#
  ggg   <- c(  63, 158, 202, 186, 121,  80,  40) #---- Green pivots. ----------------------#
  bbb   <- c(  20,  33,  56,  33,   5,  17,  28) #---- Blue pivots. -----------------------#
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
iatlas <<- function(n){
  rrr   <- c(  55, 102, 160, 237, 149,  42,  18) #---- Red pivots. ------------------------#
  ggg   <- c(  40,  80, 121, 186, 202, 158,  63) #---- Green pivots. ----------------------#
  bbb   <- c(  28,  17,   5,  33,  56,  33,  20) #---- Blue pivots. -----------------------#
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
