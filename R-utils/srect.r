#==========================================================================================#
#==========================================================================================#
#      This function is the same as the standard rect function in R, with the exception    #
# that it also allows hatch on log scale.                                                  #
#------------------------------------------------------------------------------------------#
srect <<- function( xleft
                  , ybottom
                  , xright
                  , ytop
                  , density = NULL
                  , angle   = 45
                  , col     = NA
                  , border  = NULL
                  , lty     = par("lty")
                  , lwd     = par("lwd")
                  , ...
                  ){


    #----- Discard density if it is NA or negative. ---------------------------------------#
    if (is.numeric(density) && all(is.na(density) | density < 0)) density = NULL
    #--------------------------------------------------------------------------------------#



    #--------------------------------------------------------------------------------------#
    #      Check whether to use hatched filling or full colour filling.                    #
    #--------------------------------------------------------------------------------------#
    if (! is.null(density) && ! is.null(angle)){
       #-----------------------------------------------------------------------------------#
       #      Hatched filling, prepare to call epolygon.                                   #
       #-----------------------------------------------------------------------------------#


       #----- Check border colour. --------------------------------------------------------#
       if (is.logical(border) && ! is.na(border)){
          if (border){
             border = col
          }else{
             border = NA
          }#end if (border)
          #--------------------------------------------------------------------------------#
        }#end if (is.logical(border) && ! is.na(border))
       #-----------------------------------------------------------------------------------#



       #----- Ensure that all corners are given. ------------------------------------------#
       n = range(length(xleft), length(xright), length(ybottom),length(ytop))
       if (n[1L] == 0) stop("invalid rectangle specification")
       n = n[2L]
       #-----------------------------------------------------------------------------------#


       #----- Create multiple polygons separating them with NA. ---------------------------#
       x = rbind(rep.int(NA, n), xleft, xright, xright, xleft)[-1L]
       y = rbind(rep.int(NA, n), ybottom, ybottom, ytop, ytop)[-1L]
       #-----------------------------------------------------------------------------------#


       #----- Call "epolygon", which allows hatching in the log scale. --------------------#
       epolygon( x       = c(x)
               , y       = c(y)
               , col     = col
               , border  = border
               , lty     = lty
               , lwd     = lwd
               , density = density
               , angle   = angle
               , ...
               )#end epolygon
       #-----------------------------------------------------------------------------------#
    }else{
       #----- No hatches, call the standard rect function. --------------------------------#
       rect( xleft   = xleft
           , ybottom = ybottom
           , xright  = xright
           , ytop    = ytop
           , density = NULL
           , col     = col
           , border  = border
           , lty     = lty
           , lwd     = lwd
           , ...
           )#end rect
       #-----------------------------------------------------------------------------------#
    }#end if
    #--------------------------------------------------------------------------------------#
}#end function srect
#==========================================================================================#
#==========================================================================================#
