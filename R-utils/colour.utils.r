#==========================================================================================#
#==========================================================================================#
#      This function is just a convenience function to convert RGB into colour names, with #
# default value being 255.                                                                 #
#------------------------------------------------------------------------------------------#
RGB <<- function(R,G,B) rgb(red=R,green=G,blue=B,maxColorValue=255)
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#      This function is just a convenience function to convert HSV into colour names, with #
# default value being 0-360 for hue, and 0-100 for saturation and value.                   #
#------------------------------------------------------------------------------------------#
HSV <<- function(H,S,V) hsv(h=(H/360)%%1,s=0.01*S,v=0.01*V)
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#      This function is just a convenience function to convert colour names into the hsv   #
# space.                                                                                   #
#------------------------------------------------------------------------------------------#
col2hsv <<- function(col,alpha=FALSE) rgb2hsv(col2rgb(col=col,alpha=alpha))
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#      This function creates a washed out version of the input colours (lower saturation,  #
# higher vibrancy, same hue).                                                              #
#------------------------------------------------------------------------------------------#
washout <<- function(colour,f=c(s=0.50,v=0.50),...){
   #----- Make sure the washout factor makes sense. ---------------------------------------#
   if (length(f) == 1){
      f = c(s=f,v=f)
   }else if ((length(f) == 2) && is.null(names(f))){
      names(f) = c("s","v")
   }else if (! all(c("s","v") %in% names(f))){
      cat0(" Invalid input variable \"f\"!")
      cat0(" - Length of \"f\": ",length(f),".")
      cat0(" - Element 's' is found: ","s" %in% names(f),".")
      cat0(" - Element 'v' is found: ","v" %in% names(f),".")
      cat0(" Variable \"f\" must be a scalar, unnamed vector with 2 elements,")
      cat0("   or a vector with named elements 's' and 'v'.")
      stop(" Fix variable \"f\".")
   }#end if
   #---------------------------------------------------------------------------------------#

   c.hsv   = data.frame(t(col2hsv(colour,...)))
   c.hsv$s = pmax(0,pmin(1,c.hsv$s * (1. - f["s"])         ))
   c.hsv$v = pmax(0,pmin(1,c.hsv$v * (1. - f["v"]) + f["v"]))
   ans     = hsv(h=c.hsv$h,s=c.hsv$s,v=c.hsv$v,...)
   return(ans)
}#end function washout
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#      This function creates a brighter version of the input colours (higher saturation,   #
# higher vibrancy, same hue).                                                              #
#------------------------------------------------------------------------------------------#
brighten <<- function(colour,f=c(s=0.50,v=0.50),...){
   #----- Make sure the brighten factor makes sense. --------------------------------------#
   if (length(f) == 1){
      f = c(s=f,v=f)
   }else if ((length(f) == 2) && is.null(names(f))){
      names(f) = c("s","v")
   }else if (! all(c("s","v") %in% names(f))){
      cat0(" Invalid input variable \"f\"!")
      cat0(" - Length of \"f\": ",length(f),".")
      cat0(" - Element 's' is found: ","s" %in% names(f),".")
      cat0(" - Element 'v' is found: ","v" %in% names(f),".")
      cat0(" Variable \"f\" must be a scalar, unnamed vector with 2 elements,")
      cat0("   or a vector with named elements 's' and 'v'.")
      stop(" Fix variable \"f\".")
   }#end if
   #---------------------------------------------------------------------------------------#

   c.hsv   = data.frame(t(col2hsv(colour,...)))
   c.hsv$s = pmax(0,pmin(1,c.hsv$s * (1. - f["s"]) + f["s"]))
   c.hsv$v = pmax(0,pmin(1,c.hsv$v * (1. - f["v"]) + f["v"]))
   ans     = hsv(h=c.hsv$h,s=c.hsv$s,v=c.hsv$v,...)
   return(ans)
}#end function brighten
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#      This function creates a darker version of the input colours (higher saturation,     #
# lower vibrancy, same hue).                                                               #
#------------------------------------------------------------------------------------------#
darken <<- function(colour,f=c(s=0.50,v=0.50),...){
   #----- Make sure the brighten factor makes sense. --------------------------------------#
   if (length(f) == 1){
      f = c(s=f,v=f)
   }else if ((length(f) == 2) && is.null(names(f))){
      names(f) = c("s","v")
   }else if (! all(c("s","v") %in% names(f))){
      cat0(" Invalid input variable \"f\"!")
      cat0(" - Length of \"f\": ",length(f),".")
      cat0(" - Element 's' is found: ","s" %in% names(f),".")
      cat0(" - Element 'v' is found: ","v" %in% names(f),".")
      cat0(" Variable \"f\" must be a scalar, unnamed vector with 2 elements,")
      cat0("   or a vector with named elements 's' and 'v'.")
      stop(" Fix variable \"f\".")
   }#end if
   #---------------------------------------------------------------------------------------#

   c.hsv   = data.frame(t(col2hsv(colour,...)))
   c.hsv$s = pmax(0,pmin(1,c.hsv$s * (1. - f["s"]) + f["s"]))
   c.hsv$v = pmax(0,pmin(1,c.hsv$v * (1. - f["v"])))
   ans     = hsv(h=c.hsv$h,s=c.hsv$s,v=c.hsv$v,...)
   return(ans)
}#end function darken
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#      This function creates a blackened version of the input colours (lower saturation,   #
# lower vibrancy, same hue).                                                               #
#------------------------------------------------------------------------------------------#
blacken <<- function(colour,f=c(s=0.50,v=0.50),...){
   #----- Make sure the brighten factor makes sense. --------------------------------------#
   if (length(f) == 1){
      f = c(s=f,v=f)
   }else if ((length(f) == 2) && is.null(names(f))){
      names(f) = c("s","v")
   }else if (! all(c("s","v") %in% names(f))){
      cat0(" Invalid input variable \"f\"!")
      cat0(" - Length of \"f\": ",length(f),".")
      cat0(" - Element 's' is found: ","s" %in% names(f),".")
      cat0(" - Element 'v' is found: ","v" %in% names(f),".")
      cat0(" Variable \"f\" must be a scalar, unnamed vector with 2 elements,")
      cat0("   or a vector with named elements 's' and 'v'.")
      stop(" Fix variable \"f\".")
   }#end if
   #---------------------------------------------------------------------------------------#

   c.hsv   = data.frame(t(col2hsv(colour,...)))
   c.hsv$s = pmax(0,pmin(1,c.hsv$s * (1. - f["s"])))
   c.hsv$v = pmax(0,pmin(1,c.hsv$v * (1. - f["v"])))
   ans     = hsv(h=c.hsv$h,s=c.hsv$s,v=c.hsv$v,...)
   return(ans)
}#end function blacken
#==========================================================================================#
#==========================================================================================#
