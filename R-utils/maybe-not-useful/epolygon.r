#==========================================================================================#
#==========================================================================================#
#      This function is the same as the standard polygon function in R, with the exception #
# that it also allows hatch on log scale.                                                  #
#------------------------------------------------------------------------------------------#
epolygon <<- function(x, y = NULL, density = NULL, angle = 45, border = NULL
                     ,col = NA, lty = par("lty"), ..., fillOddEven = FALSE){

   #----- Check whether x or y are in log scale. -----------------------------------------#
   plog = ""
   xlog = par("xlog")
   ylog = par("ylog")
   if (xlog) plog = paste(plog,"x",sep="")
   if (ylog) plog = paste(plog,"y",sep="")
   #--------------------------------------------------------------------------------------#
   xy = xy.coords(x, y,log=plog)

   #----- Delete density if the input value doesn't make sense... ------------------------#
   if (is.numeric(density) && all(is.na(density) | density < 0)) density = NULL
   #--------------------------------------------------------------------------------------#



   #--------------------------------------------------------------------------------------#
   #     Plot the hatch in case both angle and density are given.                         #
   #--------------------------------------------------------------------------------------#
   if (!is.null(angle) && !is.null(density)) {
      if (missing(col) || is.null(col) || is.na(col)) col = par("fg")
      if (is.null(border)) border = col
      if (is.logical(border)){
         if (! is.na(border) && border){ 
            border = col
         }else{
            border = NA
         }#end if
      }#end if
      start        = 1
      ends         = c(seq_along(xy$x)[is.na(xy$x) | is.na(xy$y)],length(xy$x) + 1)
      num.polygons = length(ends)
      col          = rep(col, length.out = num.polygons)
      border       = rep(border, length.out = num.polygons)
      lty          = rep(lty, length.out = num.polygons)
      density      = rep(density, length.out = num.polygons)
      angle        = rep(angle, length.out = num.polygons)
      i = 1
      for (end in ends) {
         if (end > start) {
            den = density[i]
            if (is.na(den) || den < 0){
               if (R.Version()$major == "3"){
                  polygon( x       = xy$x[start:(end - 1)]
                         , y       = xy$y[start:(end - 1)]
                         , col     = col[i]
                         , border  = NA
                         , density = den
                         , lty     = lty[i]
                         , ...
                         )#end polygon
               }else{
                  .Internal(polygon( xy$x[start:(end - 1)]
                                   , xy$y[start:(end - 1)], col[i], NA, lty[i], ...))
               }#end if
            }else if (den > 0) {
               epolygon.fullhatch( x           = xy$x[start:(end - 1)]
                                 , y           = xy$y[start:(end - 1)]
                                 , xlog        = xlog
                                 , ylog        = ylog
                                 , col         = col[i]
                                 , lty         = lty[i]
                                 , density     = density[i]
                                 , angle       = angle[i]
                                 , fillOddEven = fillOddEven
                                 , ...
                                 )#end polygon.fullhatch
            }#end if
            i = i + 1
         }#end if
         start = end + 1
      }#end for
      if (R.Version()$major == "3"){
         polygon( x       = xy$x
                , y       = xy$y
                , density = 0
                , col     = NA
                , border  = border
                , lty     = lty
                , ...
                )#end polygon
      }else{
         .Internal(polygon(xy$x, xy$y, NA, border, lty, ...))
      }#end if
   }else{
      if (is.logical(border)) {
         if (!is.na(border) && border){
            border = par("fg")
         }else{
            border = NA
         }#end if
      }#end if
      
      if (R.Version()$major == "3"){
         polygon( x       = xy$x
                , y       = xy$y
                , density = NULL
                , col     = col
                , border  = border
                , lty     = lty
                , ...
                )#end polygon
      }else{
         .Internal(polygon(xy$x, xy$y, col, border, lty, ...))
      }#end if
   }#end if
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      Function the draws on line.                                                         #
#------------------------------------------------------------------------------------------#
epolygon.onehatch <<- function(x, y, x0, y0, xd, yd, xlog,ylog, fillOddEven = FALSE, ...){

   halfplane  = as.integer(xd * (y - y0) - yd * (x - x0) <= 0)
   cross      = halfplane[-1L] - halfplane[-length(halfplane)]
   does.cross = cross != 0
   if (! any(does.cross)) return()


   x1 = x[-length(x)][does.cross]
   y1 = y[-length(y)][does.cross]
   x2 = x[-1L][does.cross]
   y2 = y[-1L][does.cross]

   tt = ( ( (x1 - x0) * (y2 - y1) - (y1 - y0) * (x2 - x1) )
        / (xd * (y2 - y1) - yd * (x2 - x1)) )
   oo = order(tt)
   tsort = tt[oo]
   crossings = cumsum(cross[does.cross][oo])

   if (fillOddEven) crossings = crossings%%2

   drawline = crossings != 0
   lx = x0 + xd * tsort
   ly = y0 + yd * tsort

   lx1 = lx[-length(lx)][drawline]
   ly1 = ly[-length(ly)][drawline]
   lx2 = lx[-1L][drawline]
   ly2 = ly[-1L][drawline]


   if (xlog){
      lx1 = 10^(lx1)
      lx2 = 10^(lx2)
   }#end if

   if (ylog){
      ly1 = 10^(ly1)
      ly2 = 10^(ly2)
   }#end if
   segments(lx1, ly1, lx2, ly2, ...)
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
epolygon.fullhatch <<- function(x, y, xlog,ylog, density, angle, fillOddEven = FALSE,...) {
   
   #----- Convert x and y to log scale if needed be. --------------------------------------#
   x = c(x, x[1L])
   if (xlog) x = log10(x)
   y = c(y, y[1L])
   if (ylog) y = log10(y)
   #---------------------------------------------------------------------------------------#
   
   angle = angle%%180

   usr = par("usr")
   pin = par("pin")

   upix = ( usr[2L] - usr[1L] ) / pin[1L]
   upiy = ( usr[4L] - usr[3L] ) / pin[2L]
   upi  = c(upix,upiy)
   
   if (upi[1L] < 0) angle = 180 - angle
   if (upi[2L] < 0) angle = 180 - angle

   upi = abs(upi)
   xd  = cos(angle/180 * pi) * upi[1L]
   yd  = sin(angle/180 * pi) * upi[2L]


   if (angle < 45 || angle > 135) {
       if (angle < 45) {
         first.x = max(x)
         last.x  = min(x)
       }else{
         first.x = min(x)
         last.x  = max(x)
       }#end if
       y.shift = upi[2L]/density/abs(cos(angle/180 * pi))
       x0      = 0
       y0      = floor((min(y) - first.x * yd/xd)/y.shift) * y.shift
       y.end   = max(y) - last.x * yd/xd
       while (y0 < y.end) {
         epolygon.onehatch(x, y, x0, y0, xd, yd , xlog, ylog,fillOddEven,...)
         y0 = y0 + y.shift
       }#end while
   }else{
       if (angle < 90){
         first.y = max(y)
         last.y = min(y)
       }else{
         first.y = min(y)
         last.y = max(y)
       }#end if
       x.shift = upi[1L]/density/abs(sin(angle/180 * pi))
       x0 = floor((min(x) - first.y * xd/yd)/x.shift) * x.shift
       y0 = 0
       x.end = max(x) - last.y * xd/yd
       while (x0 < x.end) {
         epolygon.onehatch(x, y, x0, y0, xd, yd,xlog,ylog,fillOddEven,...)
         x0 = x0 + x.shift
       }#end while
   }#end if
}#end function
#==========================================================================================#
#==========================================================================================#
