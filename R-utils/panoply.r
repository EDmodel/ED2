#------------------------------------------------------------------------------------------#
#   Function that creates a nice colour scheme.                                            #
#------------------------------------------------------------------------------------------#
panoply <<- function(n,sat=1.00){
   rrr      = c(   5, 146, 230, 244, 202)*sat # Red
   ggg      = c( 113, 197, 230, 165,   0)*sat # Green
   bbb      = c( 176, 222, 230, 130,  32)*sat # Blue

   pivot    = round(seq(from=1,to=n,by=(n-1)/(length(rrr)-1)),digits=0)

   red      = as.integer(spline(x=pivot,y=rrr,n=n)$y)
   green    = as.integer(spline(x=pivot,y=ggg,n=n)$y)
   blue     = as.integer(spline(x=pivot,y=bbb,n=n)$y)

   red      = pmax(0L,pmin(255L,red  ))
   green    = pmax(0L,pmin(255L,green))
   blue     = pmax(0L,pmin(255L,blue ))
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
