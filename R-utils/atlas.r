#------------------------------------------------------------------------------------------#
#   Function that creates an atlas-like colour scheme.                                     #
#------------------------------------------------------------------------------------------#
atlas <<- function(n,sat=1.00){
   rrr   = c( 140, 191, 223, 246, 235, 199, 128,  53,   1)*sat #- Red ---------------------#
   ggg   = c(  81, 129, 194, 232, 245, 234, 205, 151, 102)*sat #- Green -------------------#
   bbb   = c(  10,  45, 125, 195, 245, 229, 193, 143,  94)*sat #- Blue --------------------#
   pivot = round(seq(from=1,to=n,by=(n-1)/(length(rrr)-1)),digits=0)

   red   = as.integer(spline(x=pivot,y=rrr,n=n)$y)
   green = as.integer(spline(x=pivot,y=ggg,n=n)$y)
   blue  = as.integer(spline(x=pivot,y=bbb,n=n)$y)

   red  [red   > 255] = 255; red  [red   < 0] = 0
   green[green > 255] = 255; green[green < 0] = 0
   blue [blue  > 255] = 255; blue [blue  < 0] = 0
   mycolsch = rgb(r=red,g=green,b=blue,maxColorValue=255)
   return(mycolsch)
}#end function atlas
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#   Function that creates an atlas-like colour scheme.                                     #
#------------------------------------------------------------------------------------------#
iatlas <<- function(n,sat=1.00){
   mycolsch = rev(atlas(n,sat=sat))
   return(mycolsch)
}#end function iatlas
#------------------------------------------------------------------------------------------#
