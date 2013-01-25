#------------------------------------------------------------------------------------------#
#   Function that creates a nice grayscale colour scheme.                                  #
#------------------------------------------------------------------------------------------#
peb <- function(n){
  rrr   <- c(  32,  64,  96, 128, 159, 191, 223) #---- Red pivots. ------------------------#
  ggg   <- c(  32,  64,  96, 128, 159, 191, 223) #---- Green pivots. ----------------------#
  bbb   <- c(  32,  64,  96, 128, 159, 191, 223) #---- Blue pivots. -----------------------#
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
#   Function that creates a nice grayscale colour scheme.                                  #
#------------------------------------------------------------------------------------------#
bep <- function(n){
  rrr   <- rev(c(  32,  64,  96, 128, 159, 191, 223)) #---- Red pivots. -------------------#
  ggg   <- rev(c(  32,  64,  96, 128, 159, 191, 223)) #---- Green pivots. -----------------#
  bbb   <- rev(c(  32,  64,  96, 128, 159, 191, 223)) #---- Blue pivots. ------------------#
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
