#==========================================================================================#
#==========================================================================================#
#     Function taylor.plot.  This is similar to the taylor.diagram (package plotrix), with #
# some adaptations to make it more flexible.                                               #
#                                                                                          #
# 1.  The size of the labels is controlled by an input variable: cex.label (default is the #
#     same size as the original).                                                          #
# 2.  The code now checks whether 'mod' is a list, vector, or array, and standardises      #
#     as a list.  This allows multiple models to be dumped to the plot at once.            #
# 3.  Observations are allowed to contain NA, in which case values will be ignored.        #
#                                                                                          #
#------------------------------------------------------------------------------------------#
taylor.plot <<- function ( obs
                         , mod
                         , maxsd           = NULL
                         , add             = FALSE
                         , pt.col          = "red"
                         , pt.bg           = pt.col
                         , pt.pch          = 19
                         , pt.cex          = 1.0
                         , pt.lwd          = 2.0
                         , pos.corr        = NA
                         , plot.obs        = TRUE
                         , obs.col         = "black"
                         , obs.cex         = 2.0
                         , xlab            = "Standard deviation (Residuals)"
                         , ylab            = "Standard deviation (Observations and Model)"
                         , zlab            = "Correlation"
                         , main            = "Taylor Diagram"
                         , show.gamma      = TRUE
                         , ngamma          = 3
                         , gamma.col       = "purple3"
                         , gamma.bg        = "white"
                         , gamma.lty       = "solid"
                         , gamma.lwd       = 1.2
                         , sd.arcs         = TRUE
                         , sd.col          = "grey40"
                         , sd.lty          = "dotted"
                         , sd.lwd          = 1.2
                         , plot.sd.obs     = TRUE
                         , sd.obs.col      = "orchid"
                         , sd.obs.lty      = "solid"
                         , sd.obs.lwd      = 2.0
                         , corr.grid       = TRUE
                         , corr.label      = c(0, .2, .4, .6, .8, .9, .95, .99)
                         , corr.col        = "black"
                         , corr.lty        = "dotdash"
                         , corr.lwd        = 1.0
                         , cex.axis        = 1.0
                         , cex.label       = 1.0
                         , normalize
                         , normalise       = if (! missing(normalize)){
                                                normalize
                                             }else{
                                                FALSE
                                             }#end if
                         , mar             = c(5,5,4,4) + 0.1
                         , ...
                         ){
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Make sure that 'obs' is a vector, and that 'mod' is a list.                        #
   #---------------------------------------------------------------------------------------#
   #----- Reference. ----------------------------------------------------------------------#
   obs   = unlist   (obs)
   keep  = is.finite(obs)
   n.obs = length   (obs)
   #----- Check which type of variable "mod" is. ------------------------------------------#
   if (is.list(mod)){
      #----- Make sure that dimensions match, otherwise leave "mod" as a list. ------------#
      n.mod = sapply(X=mod,FUN=length)
      if (any(n.mod) != n.obs && sum(n.mod) != n.obs){
         cat(" - Length(mod): ",paste(n.mod,collapse=" "),"\n")
         cat(" - Length(obs):",n.obs,"\n")
         stop(" Dimensions of 'mod' and 'obs' must match!")
      }else if(sum(n.mod) == n.obs){
         warning(" Converting 'mod' from multiple lists to a single list...")
         mod = list(unlist(mod))
      }#end if
      #------------------------------------------------------------------------------------#
   }else if(is.array(mod) || is.data.frame(mod)){
      #------------------------------------------------------------------------------------#
      #     "mod" is a matrix or array, find the matching dimension convert it to list.    #
      #------------------------------------------------------------------------------------#
      dim.mod = dim(mod)
      n.mod   = length(mod)
      n.dim   = length(dim.mod)
      o.dim   = which(dim.mod == n.obs)
      if (length(o.dim) == 0 && n.obs != n.obs){
         cat(" - Dim(mod)   : ",paste(dim.mod,collapse=" "),"\n")
         cat(" - Length(mod): ",n.mod                      ,"\n")
         cat(" - Length(obs): ",n.obs                      ,"\n")
         stop(" Either the length or one dimension of 'mod' must match 'obs' size!")
      }else if(length(o.dim) == 0 && n.mod == n.obs){
         warning(" Converting array 'mod' to a vector...")
         mod = list(unlist(mod[keep]))
      }else if(length(o.dim) > 1){
         cat (" - Dim(mod)   : ",paste(dim.mod,collapse=" "),"\n")
         cat (" - Length(mod): ",n.mod                      ,"\n")
         cat (" - Length(obs): ",n.obs                      ,"\n")
         cat (" Ambiguous: 2 or more dimensions of 'mod' match the 'obs' length...","\n")
         stop(" Hint: split the array 'mod' into a list and try again...")
      }else{
         #----- Success! Split the array into lists. --------------------------------------#
         mod = aperm(a = mod, perm = c(sequence(n.dim)[-o.dim],o.dim))
         mod = t(matrix(mod,nrow=prod(dim.mod[-o.dim]),ncol=dim.mod[o.dim]))
         mod = split(mod, col(mod))
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#
   }else{
      #------------------------------------------------------------------------------------#
      #      "mod" is something else (probably a  vector).  Convert it to a list and hope  #
      # for the best...                                                                    #
      #------------------------------------------------------------------------------------#
      mod   = unlist(mod)
      n.mod = length(mod)
      if (n.mod != n.obs){
         cat (" - Length(mod): ",n.mod,"\n")
         cat (" - Length(obs): ",n.obs,"\n")
         stop(" Dimensions of 'mod' and 'obs' must match...")
      }else{
         mod = list(unlist(mod))
      }#end if
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Keep only the valid observations. -----------------------------------------------#
   # obs  = obs[keep]
   # mod  = lapply(X=mod,FUN="[",keep)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     In case the user asked positive correlation only but there are negative           #
   # correlation points, warn the user.                                                    #
   #---------------------------------------------------------------------------------------#
   R      = sapply(X = mod, FUN = cor, x = obs, use = "pairwise")
   if (is.na(pos.corr)){
      pos.corr = all(R %>=% 0. | is.na(R))
   }else if (pos.corr && any(R %<% 0)){
      warning(" There are negative correlations, but you've asked positive side only!")
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Find the coefficient of correlation and standard deviations.                       #
   #---------------------------------------------------------------------------------------#
   sd.obs.orig = sd    (x=obs,na.rm=TRUE)
   sd.mod.orig = sapply(X=mod,FUN=sd,na.rm=TRUE)
   if (normalise) {
       sd.mod = sd.mod.orig / sd.obs.orig
       sd.obs = 1.0
   }else{
       sd.mod = sd.mod.orig
       sd.obs = sd.obs.orig
   }#end if
   if (is.null(maxsd)){
      if (all(is.na(c(sd.mod,sd.obs)))){
         maxsd = c(0,1)
      }else{
         maxsd = 1.2 * max(unlist(c(sd.mod, sd.obs)),na.rm=TRUE)
      }#end if
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Check whether this is a new plot or not.  If it is a new plot, create the         #
   # diagram background before adding the new information.                                 #
   #---------------------------------------------------------------------------------------#
   if (! add){
       #----- Save previous PAR settings. -------------------------------------------------#
       par.all  = par(no.readonly=FALSE)
       par.orig = par(no.readonly=TRUE )
       on.exit(par.orig)
       #-----------------------------------------------------------------------------------#



       #-----------------------------------------------------------------------------------#
       #      Define some settings differently depending on whether this is a positive-    #
       # -only or full correlation Taylor plot.                                            #
       #-----------------------------------------------------------------------------------#
       quarter    = seq(from=0,to= 90,by=1) * pio180
       half       = seq(from=0,to=180,by=1) * pio180
       if (pos.corr){
          xlim           = c(0,maxsd)
          ylim           = c(0,maxsd)
          x.xlab         = sd.obs
          y.xlab         = -0.15 * maxsd
          srt.xlab       =  0
          x.ylab         = -0.15 * maxsd
          y.ylab         =   0.5 * maxsd
          srt.ylab       = 90
          x0.axis        = c(     0,     0)
          y0.axis        = c(     0,     0)
          x1.axis        = c( maxsd,     0)
          y1.axis        = c(     0, maxsd)
          r.gamma        = half
          r.sdarc        = quarter
          r.corr         = quarter
          corr.label     = corr.label[corr.label != 0]
          corr.at        = acos(corr.label)
          corr.angle     = 45. * pio180
          corr.angle.lab = acos(cos(corr.angle) * par.all$din[2] / par.all$din[1])
       }else{
          xlim           = c(-maxsd,maxsd)
          ylim           = c(     0,maxsd)
          x.xlab         = sd.obs
          y.xlab         = -0.20 * maxsd
          srt.xlab       =  0
          x.ylab         =  0
          y.ylab         = -0.20 * maxsd
          srt.ylab       =  0
          x0.axis        = c(-maxsd,     0)
          y0.axis        = c(     0,     0)
          x1.axis        = c( maxsd,     0)
          y1.axis        = c(     0, maxsd)
          r.gamma        = half
          r.sdarc        = half
          r.corr         = half
          corr.label     = sort(unique(c(-corr.label,0,corr.label)))
          corr.at        = acos(corr.label)
          corr.angle     = 90. * pio180
          corr.angle.lab = corr.angle
       }#end if
       #-----------------------------------------------------------------------------------#




       #----- Start up a new plotting window. ---------------------------------------------#
       par(mar=mar,xpd=TRUE)
       plot.new()
       plot.window(xlim=xlim, ylim = ylim,xaxs="i",yaxs="i",...)
       title(main=main,cex=cex.axis)
       par(xpd = TRUE)
       text(x=x.xlab,y=y.xlab,labels=xlab,srt=srt.xlab,cex=cex.axis,col=gamma.col,font=2)
       text(x=x.ylab,y=y.ylab,labels=ylab,srt=srt.ylab,cex=cex.axis,col=sd.col,font=2)
       #-----------------------------------------------------------------------------------#




       #----- Plot the correlation grid. --------------------------------------------------#
       if (corr.grid) {
          for (coat in corr.at){
             lines( x   = c(0, maxsd * cos(coat))
                  , y   = c(0, maxsd * sin(coat))
                  , col = corr.col
                  , lty = corr.lty
                  , lwd = corr.lwd
                  )#end lines
          }#end for
          #--------------------------------------------------------------------------------#
       }#end if
       #-----------------------------------------------------------------------------------#



       #----- Plot the axis lines and labels. ---------------------------------------------#
       segments( x0  = x0.axis, y0 = y0.axis, x1 = x1.axis, y1 = y1.axis
               , col = par.orig$fg, lty = par.orig$lty, lwd = par.orig$lwd)
       axis.ticks = pretty(xlim)
       axis.ticks = axis.ticks[abs(axis.ticks) <= maxsd]
       axis( side = 1 + pos.corr, at = axis.ticks, las = 1, labels = abs(axis.ticks)
           , col.axis  = sd.col, cex.axis = cex.axis)
       #-----------------------------------------------------------------------------------#




       #-----------------------------------------------------------------------------------#
       #      Plot the standard deviation arcs if the user wants it.                       #
       #-----------------------------------------------------------------------------------#
       plot.sd.arcs = (  ( is.logical(sd.arcs) && any(sd.arcs)      )
                      || ( is.numeric(sd.arcs) && any(sd.arcs != 0) ) )
       if (plot.sd.arcs) {
           #---- Grab the axis ticks to plot the curves. ----------------------------------#
           if (is.logical(sd.arcs)) sd.arcs = axis.ticks[axis.ticks >= 0]
           #-------------------------------------------------------------------------------#



           #-------------------------------------------------------------------------------#
           #   Plot all arcs.                                                              #
           #-------------------------------------------------------------------------------#
           for (sdarc in sd.arcs) {
              x.curve = cos(r.sdarc) * sdarc
              y.curve = sin(r.sdarc) * sdarc
              lines(x = x.curve, y = y.curve, col = sd.col, lty = sd.lty, lwd = sd.lwd)
           }#end for (sdarc in sd.arcs)
           #-------------------------------------------------------------------------------#
       }# if (plot.sd.arcs)
       #-----------------------------------------------------------------------------------#




       #-----------------------------------------------------------------------------------#
       #     Check whether to show the concentric circles around the target.               #
       #-----------------------------------------------------------------------------------#
       plot.gamma = (  ( is.logical(show.gamma) && any(show.gamma     ) )
                    || ( is.numeric(show.gamma) && any(show.gamma != 0) ) )

       if (plot.gamma) {
           #----- Use the default circles. ------------------------------------------------#
           if (is.logical(show.gamma)){
              show.gamma = pretty(c(0, maxsd), n = ngamma)[-1]
           }#end if
           ngamma     = length(show.gamma)
           #-------------------------------------------------------------------------------#


           #-------------------------------------------------------------------------------#
           #     Plot the curves.                                                          #
           #-------------------------------------------------------------------------------#
           for (gg in 1:ngamma){
              x.curve  = cos(r.gamma) * show.gamma[gg] + sd.obs
              y.curve  = sin(r.gamma) * show.gamma[gg]
              r2.curve = x.curve * x.curve + y.curve * y.curve
              if (pos.corr){
                 bye      = r2.curve > maxsd * maxsd | x.curve < 0
              }else{
                 bye      = r2.curve > maxsd * maxsd
              }#end if
              x.curve[bye] = NA
              y.curve[bye] = NA
              idx          = floor(median(which(! bye)))
              lines( x   = x.curve
                   , y   = y.curve
                   , col = gamma.col
                   , lwd = gamma.lwd
                   , lty = gamma.lty
                   )#end lines
              #----- Labels. --------------------------------------------------------------#
              boxed.labels( x      = x.curve[idx]
                          , y      = y.curve[idx]
                          , labels = show.gamma[gg]
                          , col    = gamma.col
                          , bg     = gamma.bg
                          , cex    = cex.label
                          , border = FALSE
                          )#end border.labels
              #----------------------------------------------------------------------------#
           }#end for (gg in 1:ngamma)
       }#end if (plot.gamma)
       #-----------------------------------------------------------------------------------#



       #-----------------------------------------------------------------------------------#
       #     Plot the outer curves and tick marks ("correlation axis").                    #
       #-----------------------------------------------------------------------------------#
       #------ Outer curve. ---------------------------------------------------------------#
       x.curve = cos(r.corr) * maxsd
       y.curve = sin(r.corr) * maxsd
       lines( x   = x.curve
            , y   = y.curve
            , col = par.orig$fg
            , lty = par.orig$lty
            , lwd = par.orig$lwd
            )#end lines
       x0.at   = 1.00  * cos(corr.at) * maxsd
       x1.at   = 1.025 * cos(corr.at) * maxsd
       y0.at   = 1.00  * sin(corr.at) * maxsd
       y1.at   = 1.025 * sin(corr.at) * maxsd
       #------ Tick marks. ----------------------------------------------------------------#
       segments( x0  = x0.at
               , y0  = y0.at
               , x1  = x1.at
               , y1  = y1.at
               , col = par.orig$fg
               , lty = par.orig$lty
               , lwd = par.orig$lwd
               )#end segments
       #------ Axis tick labels. ----------------------------------------------------------#
       for (cc in 1:length(corr.at)){
          text ( x      = 1.05 * cos(corr.at[cc]) * maxsd
               , y      = 1.05 * sin(corr.at[cc]) * maxsd
               , labels = corr.label[cc]
               , col    = corr.col
               , cex    = 1.1 * cex.label
               , srt    = 180 * corr.at[cc] / pi - 90
               )#end text
       }#end for (cc in 1:length(corr.at))
       #------ Axis label. ----------------------------------------------------------------#
       text( x      = cos(corr.angle) * maxsd * 1.15
           , y      = sin(corr.angle) * maxsd * 1.15
           , labels = zlab
           , col    = corr.col
           , cex    = cex.axis
           , srt    = 180. * corr.angle.lab / pi - 90.
           , font   = 2
           )#end text
       #-----------------------------------------------------------------------------------#


       #-----------------------------------------------------------------------------------#
       #    Check whether to plot the reference standard deviation arc.                    #
       #-----------------------------------------------------------------------------------#
       if (plot.sd.obs) {
           x.curve = cos(r.sdarc) * sd.obs
           y.curve = sin(r.sdarc) * sd.obs
           lines(x.curve, y.curve,col=sd.obs.col,lty=sd.obs.lty,lwd=sd.obs.lwd)
       }#end if
       #-----------------------------------------------------------------------------------#


       #-----------------------------------------------------------------------------------#
       #    Check whether to plot a point showing the sweetest spot.                       #
       #-----------------------------------------------------------------------------------#
       if (plot.obs) {
          points( x = sd.obs, y = 0, pch=16, cex = 2/3 * obs.cex, col = obs.col)
          points( x = sd.obs, y = 0, pch=21, cex =       obs.cex, col = obs.col)
       }#end if
       #-----------------------------------------------------------------------------------#
   }#end if (! add)
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Plot the points.                                                                  #
   #---------------------------------------------------------------------------------------#
   points( x   = sd.mod * R
         , y   = sd.mod * sin(acos(R))
         , pch = pt.pch
         , col = pt.col
         , bg  = pt.bg
         , cex = pt.cex
         , lwd = pt.lwd
         )#end points
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Return the basic statistics.                                                       #
   #---------------------------------------------------------------------------------------#
   ans = list( corr = R, sd.obs = sd.obs.orig, sd.mod = sd.mod.orig)
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function
#==========================================================================================#
#==========================================================================================#
