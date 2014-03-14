#==========================================================================================#
#==========================================================================================#
#     Function radial.flex                                                                 #
#                                                                                          #
#   This function is almost the same as radial.plot (package plotrix), but it also allows  #
# the axes to be expressions.  It also has the option to set the background colour for     #
# boxed labels.                                                                            #
#------------------------------------------------------------------------------------------#
radial.flex <<- function ( lengths
                         , radial.pos       = NULL
                         , labels           = NA
                         , label.pos        = NULL
                         , radlab           = FALSE
                         , lab.col          = par("fg")
                         , lab.bg           = par("bg")
                         , lab.cex          = par("cex.axis")
                         , start            = 0
                         , clockwise        = FALSE
                         , rp.type          = "r"
                         , pt.type          = "p"
                         , label.prop       = 1.15
                         , main             = ""
                         , xlab             = ""
                         , ylab             = ""
                         , line.col         = par("fg")
                         , lty              = par("lty")
                         , lwd              = par("lwd")
                         , mar              = c(2, 2, 3, 2)
                         , show.grid        = TRUE
                         , show.grid.labels = 4
                         , show.radial.grid = TRUE
                         , show.radial.edge = FALSE
                         , grid.col         = "grey"
                         , grid.bg          = "transparent"
                         , grid.left        = FALSE
                         , grid.unit        = NULL
                         , point.symbols    = NULL
                         , point.col        = NULL
                         , show.centroid    = FALSE
                         , radial.lim       = NULL
                         , radial.labels    = NULL
                         , radial.col       = lab.col
                         , radial.bg        = lab.bg
                         , radial.cex       = par("cex.lab")
                         , boxed.radial     = TRUE
                         , poly.col         = NULL
                         , add              = FALSE
                         , ...
                         ){



   #---------------------------------------------------------------------------------------#
   #     Save the current par, we will revert back upon exit.                              #
   #---------------------------------------------------------------------------------------#
   oldpar = par("xpd", "mar", "pty")
   on.exit(par(oldpar))
   #---------------------------------------------------------------------------------------#



   #----- Define default radial limit in case none has been given. ------------------------#
   if (is.null(radial.lim)) radial.lim = range(lengths,na.rm=TRUE)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Find out whether the data are a vector or a matrix.  We will plot each column as  #
   # a separate vector in case it is a matrix.                                             #
   #---------------------------------------------------------------------------------------#
   length.dim = dim(lengths)
   if (is.null(length.dim)) {
       npoints = length(lengths)
       nsets   = 1
       lengths = matrix(lengths, nrow = 1)
   }else{
       npoints = length.dim[2]
       nsets   = length.dim[1]
       lengths = as.matrix(lengths)
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Make sure that none of the lengths are negative.                                   #
   #---------------------------------------------------------------------------------------#
   lengths              = lengths - radial.lim[1]
   lengths[lengths < 0] = NA
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #      If the azimuths have not been given, assume equal spread around the circle.      #
   #---------------------------------------------------------------------------------------#
   if (is.null(radial.pos[1])){
       radial.pos = seq(from = 0, to = pi * (2 - 2/npoints), length.out = npoints)
   }#end if
   #---- Make sure that the radial position is also a matrix. -----------------------------#
   radial.pos.dim = dim(radial.pos)
   if (is.null(radial.pos.dim)){
      radial.pos = matrix(rep(radial.pos, nsets), nrow = nsets, byrow = TRUE)
   }else{
      radial.pos = as.matrix(radial.pos)
   }#end if
   #---- Revert direction in case the radial plot is to go clockwise. ---------------------#
   if (clockwise) radial.pos = -radial.pos
   #---- Offset the first position in case it is not supposed to start at angle 0. --------#
   if (start) radial.pos = radial.pos + start
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Set up grid properties.                                                           #
   #---------------------------------------------------------------------------------------#
   if (show.grid){
       #---- If radial.lim is just the limit, find the nice places to add grid. -----------#
       if (length(radial.lim) < 3){
          grid.pos = pretty(radial.lim)
       }else{
          grid.pos = radial.lim
       }#end if
       if (grid.pos[1] < radial.lim[1]) grid.pos = grid.pos[-1]
       maxlength = max(grid.pos - radial.lim[1])
       angles    = seq(from=0,to=360,by=0.5) * pio180
   }else{
       grid.pos  = NA
       maxlength = diff(radial.lim)
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Open the plotting window in case we aren't simply adding more points.            #
   #---------------------------------------------------------------------------------------#
   if (! add) {
       par(mar = mar, pty = "s")
       plot.new()
       plot.window( xlim = c(-maxlength, maxlength), ylim = c(-maxlength, maxlength),...)
       title(main = main, xlab = xlab, ylab = ylab)
       #---- Add the grid in case it is sought. -------------------------------------------#
       if (show.grid){
          for (i in rev(seq_along(grid.pos))){
             xpos = cos(angles) * (grid.pos[i] - radial.lim[1])
             ypos = sin(angles) * (grid.pos[i] - radial.lim[1])
             polygon(xpos, ypos, border = grid.col, col = grid.bg)
          }#end for
       }#end if
   }#end if (! add)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Make sure the spider settings all have the same length.                          #
   #---------------------------------------------------------------------------------------#
   par(xpd = TRUE)
   if (length(line.col     ) < nsets) line.col      = sequence(nsets)
   if (length(rp.type      ) < nsets) rp.type       = rep(rp.type      ,length.out=nsets)
   if (length(point.symbols) < nsets) point.symbols = rep(point.symbols,length.out=nsets)
   if (length(point.col    ) < nsets) point.col     = rep(point.col    ,length.out=nsets)
   if (length(poly.col     ) < nsets) poly.col      = rep(poly.col     ,length.out=nsets)
   if (length(lty          ) < nsets) lty           = rep(lty          ,length.out=nsets)
   if (length(lwd          ) < nsets) lwd           = rep(lwd          ,length.out=nsets)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Loop over all elements.                                                           #
   #---------------------------------------------------------------------------------------#
   for (i in 1:nsets){
      if (nsets > 1){
         linecol      = line.col     [i]
         polycol      = poly.col     [i]
         pointcol     = point.col    [i]
         pointsymbols = point.symbols[i]
         ltype        = lty          [i]
         lwidth       = lwd          [i]
      }else{
         linecol      = line.col
         polycol      = poly.col
         pointcol     = point.col
         pointsymbols = point.symbols
         ltype        = lty
         lwidth       = lwd
      }#end if

      #----- Decide which type of radial plot, and set up the point symbols and colour. ---#
      rptype = unlist(strsplit(rp.type[i], ""))
      if (match("s", rptype, 0)) {
          if (is.null(pointsymbols)) pointsymbols = i
          if (is.null(pointcol))     pointcol     = i
      }#end if
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Translate radius and angle into X and Y.                                       #
      #------------------------------------------------------------------------------------#
      xpos = cos(radial.pos[i, ]) * lengths[i, ]
      ypos = sin(radial.pos[i, ]) * lengths[i, ]
      if (match("r", rptype, 0)){
         segments(0,0,xpos,ypos,col=linecol,lty=ltype,lwd=lwidth,...)
      }#end if
      if (match("p", rptype, 0)){
         polygon(xpos,ypos,border=linecol,col=polycol,lty=ltype,lwd=lwidth,...)
      }#end if
      if (match("s", rptype, 0)){
         points(xpos,ypos,pch=pointsymbols,col=pointcol,type=pt.type,...)
      }#end if
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Plot the centre of mass (roughly).                                             #
      #------------------------------------------------------------------------------------#
      if (show.centroid){
         #---------------------------------------------------------------------------------#
         #      Check whether this is a point centroid or not.                             #
         #---------------------------------------------------------------------------------#
         if (match("p", rptype, 0)) {
            nvertices   = length(xpos)
            polygonarea = xpos[nvertices] * ypos[1] - xpos[1] * ypos[nvertices]
            for (vertex in 1:(nvertices - 1)){
               polygonarea = ( polygonarea 
                             + xpos[vertex  ] * ypos[vertex+1]
                             - xpos[vertex+1] * ypos[vertex  ]
                             )
            }#end for

            polygonarea = polygonarea/2
            centroidx = ( (xpos[nvertices] + xpos[1])
                        * (xpos[nvertices] * ypos[1] - xpos[1] * ypos[nvertices]) )
            centroidy = ( (ypos[nvertices] + ypos[1]) 
                        * (xpos[nvertices] * ypos[1] - xpos[1] * ypos[nvertices]) )
            for (vertex in 1:(nvertices - 1)){
               centroidx = ( centroidx 
                           + ( xpos[vertex  ] + xpos[vertex+1] )
                           * ( xpos[vertex  ] * ypos[vertex+1]
                             - xpos[vertex+1] * ypos[vertex  ] ) )
               centroidy = ( centroidy
                           + ( ypos[vertex  ] + ypos[vertex+1] )
                           * ( xpos[vertex  ] * ypos[vertex+1]
                             - xpos[vertex+1] * ypos[vertex  ] ) )
            }#end for

            points( x   = centroidx/(6 * polygonarea)
                  , y   = centroidy/(6 * polygonarea)
                  , col = point.col[i]
                  , pch = point.symbols[i]
                  , cex = 2
                  , ...
                  )#end points
         }else{
            points(x=mean(xpos),y=mean(ypos),col=pointcol,pch=pointsymbols,cex=2,...)
         }#end if
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Check whether this is a new plot, and add labels.                                 #
   #---------------------------------------------------------------------------------------#
   if (! add){
      #---- Find out whether labels are expressions or regular text. ----------------------#
      if (is.expression(labels[1])){
         is.exp = TRUE
      }else{
         is.exp = FALSE
         #----- Create default angles (in degrees) for labels. ----------------------------#
         if (is.na(labels[1])){
            label.deg = seq(from = 0, to = 330, by = 30)
            label.pos = label.deg * pio180
            labels    = as.character(round(label.deg, 0))
         }#end if
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Find where to place the labels. 
      #------------------------------------------------------------------------------------#
      if (is.null(label.pos[1])) {
          lablen = length(labels)
          label.pos = seq(0, pi * (2 - 2/lablen), length.out = lablen)
      }#end if
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Adjust the labels in case the start is not 0 or the plot is clockwise, then   #
      # find the (x;y) coordinates associated with the labels and plot the lines.          #
      #------------------------------------------------------------------------------------#
      if (clockwise) label.pos = - label.pos
      if (start    ) label.pos =   label.pos + start
      xpos = cos(label.pos) * maxlength
      ypos = sin(label.pos) * maxlength
      if (show.radial.grid) segments(0, 0, xpos, ypos, col = grid.col)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Plot the angle labels.                                                         #
      #------------------------------------------------------------------------------------#
      xpos = cos(label.pos) * maxlength * label.prop
      ypos = sin(label.pos) * maxlength * label.prop
      if (radlab) {
         for (label in sequence(length(labels))){
            lab.srt = ( (180 * label.pos[label]/pi)
                      +  180 * (label.pos[label] > pi/2 && label.pos[label] < 3 * pi/2) )
            text( x      = xpos[label]
                , y      = ypos[label]
                , labels = if (is.exp){as.expression(labels[label])}else{labels[label]}
                , cex    = lab.cex
                , srt    = lab.srt
                , col    = lab.col
                , bg     = lab.bg
                )#end text
          }#end for (label in sequence(length(labels)))
      }else{
         boxed.labels( x      = xpos
                     , y      = ypos
                     , labels = if (is.exp){as.expression(labels)}else{labels}
                     , ypad   = 0.7
                     , border = FALSE
                     , col    = lab.col
                     , bg     = lab.bg
                     , cex    = lab.cex
                     )#end boxed.labels
      }#end if (radlab)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Plot the radius labels.                                                        #
      #------------------------------------------------------------------------------------#
      if (show.grid.labels) {
         #---- Get the X and Y coordinates to add the labels. -----------------------------#
         if ( show.grid.labels %% 2){
            ypos = grid.pos - radial.lim[1]
            xpos = rep(0, length(grid.pos))
            if (show.grid.labels == 1) ypos = -ypos
         }else{
            xpos = grid.pos - radial.lim[1]
            ypos = rep(0, length(grid.pos))
            if (show.grid.labels == 2) xpos = -xpos
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Make the labels, and append the units to the last element if the user       #
         # provided it.                                                                    #
         #---------------------------------------------------------------------------------#
         if (is.null(radial.labels)) radial.labels = as.character(grid.pos)
         n.pos            = length(grid.pos)
         if (! is.null(grid.unit)){
            radial.labels[n.pos] = paste(radial.labels[n.pos],grid.unit)
         }else if (! show.radial.edge){
            xpos          = xpos         [-n.pos]
            ypos          = ypos         [-n.pos]
            radial.labels = radial.labels[-n.pos]
            n.pos         = n.pos-1
         }#end if
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Plot the labels using either boxed or regular labels.                       #
         #---------------------------------------------------------------------------------#
         if (boxed.radial){
            boxed.labels( x      = xpos
                        , y      = ypos
                        , labels = radial.labels
                        , border = FALSE
                        , cex    = radial.cex
                        , col    = radial.col
                        , bg     = radial.bg
                        )#end boxed.labels
         }else{
            text( x      = xpos
                , y      = ypos
                , labels = radial.labels
                , col    = radial.col
                , bg     = radial.bg
                , cex    = radial.cex
                )#end text
         }#end if (boxed.radial)
         #---------------------------------------------------------------------------------#
      }#end if (show.grid.labels)
      #------------------------------------------------------------------------------------#
   }#end if (! add)
   #---------------------------------------------------------------------------------------#


   invisible()
}#end function radial.flex
#==========================================================================================#
#==========================================================================================#
