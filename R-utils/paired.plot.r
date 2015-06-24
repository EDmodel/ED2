#==========================================================================================#
#==========================================================================================#
#     This function is very similar to pairs, from which most of the code was taken. The   #
# main difference is that paired.plot allows for different x and y limits for upper and    #
# lower panels, which may be handy.  Another minor differences:                            #
# 1.  x and y labels are drawn on every corner plot.                                       #
# 2.  Unless diag.box = TRUE or diag.panel is a function, it does not plot the box in the  #
#     diagonal.                                                                            #
# 3.  It allows for xlab, ylab, and more specific labels for lower and upper labels.       #
#------------------------------------------------------------------------------------------#
paired.plot <<- function ( x
                         , diag.labels
                         , panel       = points
                         , lower.panel = panel
                         , upper.panel = panel
                         , diag.panel  = NULL
                         , text.panel  = textPanel
                         , xlim        = NULL
                         , ylim        = NULL
                         , lower.xlim  = NULL
                         , lower.ylim  = NULL
                         , upper.xlim  = NULL
                         , upper.ylim  = NULL
                         , label.pos   = 0.5 + has.diag/3
                         , line.main   = 3
                         , line.xlab   = 3
                         , line.ylab   = 3
                         , cex.diag    = NULL
                         , font.diag   = 1
                         , row1attop   = TRUE
                         , gap         = 1
                         , xlog        = FALSE
                         , lower.xlog  = xlog
                         , upper.xlog  = xlog
                         , ylog        = FALSE
                         , lower.ylog  = ylog
                         , upper.ylog  = ylog
                         , lower.xlab  = NULL
                         , upper.xlab  = NULL
                         , lower.ylab  = NULL
                         , upper.ylab  = NULL
                         , lower.xaxis = NULL
                         , upper.xaxis = NULL
                         , lower.yaxis = NULL
                         , upper.yaxis = NULL
                         , diag.box    = ! is.null(diag.panel)
                         , ...
                         ){

   #---------------------------------------------------------------------------------------#
   #     Save some Boolean flags.                                                          #
   #---------------------------------------------------------------------------------------#
   doText    = missing(text.panel)
   has.lower = ! is.null(lower.panel)
   has.upper = ! is.null(upper.panel)
   has.diag  = ! is.null(diag.panel )
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Define generic functions for plotting boxes.                                      #
   #---------------------------------------------------------------------------------------#
   if (doText || is.function(text.panel)){ 
      textPanel = function(x = 0.5, y = 0.5, txt, cex, font){
         text(x,y, txt, cex = cex, font = font)
      }#end textPanel
   }#end if
   localAxis       = function(side,bg,col=NULL,main,oma,...)       axis(side=side,...)
   localLowerPanel = function(..., main, oma, font.main, cex.main) lower.panel(...)
   localUpperPanel = function(..., main, oma, font.main, cex.main) upper.panel(...)
   localDiagPanel  = function(..., main, oma, font.main, cex.main) diag.panel (...)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Copy dots to a list.                                                               #
   #---------------------------------------------------------------------------------------#
   dots   = list(...)
   nmdots = names(dots)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Make sure x is a nice numeric matrix.                                              #
   #---------------------------------------------------------------------------------------#
   if (! is.matrix(x)) {
      x   = as.data.frame(x)
      nnx = length(names(x))
      #----- Loop over columns. -----------------------------------------------------------#
      for (i in sequence(nnx)){
         #----- Don't make a fuss if the variable is logical of factorial. ----------------#
         if (is.factor(x[[i]]) || is.logical(x[[i]])){
            x[[i]] <- as.numeric(x[[i]])
         }#end if (is.factor(x[[i]]) || is.logical(x[[i]]))
         if (! is.numeric(unclass(x[[i]]))){ 
            stop("non-numeric argument to 'paired.plot'")
         }#end if (! is.numeric(unclass(x[[i]])))
         #---------------------------------------------------------------------------------#
      }#end for (i in seq_along(names(x)))
      #------------------------------------------------------------------------------------#
   }else if (! is.numeric(x)){
      stop("non-numeric argument to 'paired.plot'")
   }#end if (! is.numeric(x))
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Get functions.                                                                    #
   #---------------------------------------------------------------------------------------#
   panel <- match.fun(panel)
   if (has.lower && (! missing(lower.panel))) lower.panel = match.fun(lower.panel)
   if (has.upper && (! missing(upper.panel))) upper.panel = match.fun(upper.panel)
   if (has.diag  && (! missing(diag.panel ))) diag.panel  = match.fun(diag.panel )
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Set x and y limits.                                                                #
   #---------------------------------------------------------------------------------------#
   if (is.null(xlim      )) xlim       = range(c(x),finite=TRUE)
   if (is.null(lower.xlim)) lower.xlim = xlim
   if (is.null(lower.xlog)) lower.xlog = xlog
   if (is.null(upper.xlim)) upper.xlim = xlim
   if (is.null(upper.xlog)) upper.xlog = xlog
   if (is.null(ylim      )) ylim       = range(c(x),finite=TRUE)
   if (is.null(lower.ylim)) lower.ylim = ylim
   if (is.null(lower.ylog)) lower.ylog = ylog
   if (is.null(upper.ylim)) upper.ylim = ylim
   if (is.null(upper.ylog)) upper.ylog = ylog
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Swap lower and upper functions in case row #1 is at the top.                       #
   #---------------------------------------------------------------------------------------#
   if (row1attop) {
      placeholder = lower.panel
      lower.panel = upper.panel
      upper.panel = placeholder

      placeholder = has.lower
      has.lower   = has.upper
      has.upper   = placeholder

      placeholder = lower.xlim
      lower.xlim  = upper.xlim
      upper.xlim  = placeholder

      placeholder = lower.xlog
      lower.xlog  = upper.xlog
      upper.xlog  = placeholder

      placeholder = lower.ylim
      lower.ylim  = upper.ylim
      upper.ylim  = placeholder

      placeholder = lower.ylog
      lower.ylog  = upper.ylog
      upper.ylog  = placeholder
   
      right.ylab  = upper.ylab
      left.ylab   = lower.ylab
   }else{
      right.ylab  = lower.ylab
      left.ylab   = upper.ylab
   }#end if (row1attop)
   #---------------------------------------------------------------------------------------#



   #----- Number of columns. --------------------------------------------------------------#
   ncolx <- ncol(x)
   if (ncolx < 2) stop("only one column in the argument to 'paired.plot'")
   #---------------------------------------------------------------------------------------#



   #----- Check labels for diagonal elements. ---------------------------------------------#
   if (doText){
      if (missing(diag.labels)) {
         diag.labels = colnames(x)
         if (is.null(diag.labels)) diag.labels = paste("var", sequence(ncolx))
      }else if (is.null(diag.labels)){ 
         doText = FALSE
      }#end if (missing(labels))
   }#end if (doText)
   #---------------------------------------------------------------------------------------#



   #----- Grab the main if it has been provided. ------------------------------------------#
   if ("main" %in% nmdots){
      main = dots$main
   }else{
      main = NULL
   }#end if ("main" %in% nmdots)
   #---------------------------------------------------------------------------------------#



   #----- Number of columns. --------------------------------------------------------------#
   if ("oma" %in% nmdots){
      oma = dots$oma
   }else{
      oma = c( 4 + (1 - is.null(lower.xlab))
             , 4 + (1 - is.null(left.ylab ))
             , 4 + (1 - is.null(upper.ylab)) + 2 * (1 - is.null(main))
             , 4 + (1 - is.null(right.ylab))
             )#end c
   }#end if ("oma" %in% nmdots)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Save some configuration for later.                                                #
   #---------------------------------------------------------------------------------------#
   opar <- par(mfrow = c(ncolx, ncolx), mar = rep.int(gap/2, 4), oma = oma)
   on.exit(par(opar))
   #---------------------------------------------------------------------------------------#


   #------ Make sure things pop up in the end only. ---------------------------------------#
   dev.hold()
   on.exit(dev.flush(), add = TRUE)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #       Vectors with loop elements (rows and columns).                                  #
   #---------------------------------------------------------------------------------------#
   if (row1attop){
      row.loop = sequence(ncolx)
   }else{
      row.loop = rev(sequence(ncolx))
   }#end if (row1attop)
   col.loop = sequence(ncolx)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Outer loop: rows.                                                                 #
   #---------------------------------------------------------------------------------------#
   for (rr in row.loop){
      #------------------------------------------------------------------------------------#
      #     Inner loop: columns.                                                           #
      #------------------------------------------------------------------------------------#
      for (cc in col.loop) {
         #---------------------------------------------------------------------------------#
         #     Copy some settings depending on whether this is upper or lower plot.        #
         #---------------------------------------------------------------------------------#
         if (rr <= cc){
            xlog = lower.xlog
            ylog = lower.ylog
            xlim = lower.xlim
            ylim = lower.ylim
         }else{
            xlog = upper.xlog
            ylog = upper.ylog
            xlim = upper.xlim
            ylim = upper.ylim
         }#end if (rr <= cc)
         boxed = (rr != cc) || diag.box
         plog  = paste0(ifelse(xlog,"x",""),ifelse(ylog,"y",""))
         #---------------------------------------------------------------------------------#
       

         #---------------------------------------------------------------------------------#
         #     Open the plotting window.                                                   #
         #---------------------------------------------------------------------------------#
         plot.new()
         plot.window(xlim=xlim,ylim=ylim,log=plog,...)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Plot this or skip it.                                                       #
         #---------------------------------------------------------------------------------#
         if (rr == cc || ((rr < cc) && has.lower) || ((rr > cc) && has.upper)){
         
            #----- Check whether to plot box. ---------------------------------------------#
            if (rr != cc || diag.box){
               box()
            }#end if
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #    Decide whether to plot lower of upper axis.                               #
            #------------------------------------------------------------------------------#
            if (row1attop){
               draw.x.upper = (rr == 1)     && (rr != cc) && has.upper
               draw.x.lower = (rr == ncolx) && (rr != cc) && has.lower
               draw.y.upper = (cc == ncolx) && (rr != cc) && has.upper
               draw.y.lower = (cc == 1)     && (rr != cc) && has.lower
            }else{
               draw.x.upper = (rr == ncolx) && (rr != cc) && has.upper
               draw.x.lower = (rr == 1    ) && (rr != cc) && has.lower
               draw.y.upper = (cc == 1    ) && (rr != cc) && has.upper
               draw.y.lower = (cc == ncolx) && (rr != cc) && has.lower
            }#end if (row1attop)
            if (draw.x.upper){
               #----- Upper x axis. -------------------------------------------------------#
               if (is.null(upper.xaxis)){
                  localAxis(side=3,...)
               }else{
                  do.call( what = "localAxis"
                         , args = modifyList(x=upper.xaxis,val=list(side=3))
                         )#end do.call
               }#end if (is.null(upper.xaxis))
               #---------------------------------------------------------------------------#
            }#end if
            if (draw.x.lower){
               #----- Lower x axis. -------------------------------------------------------#
               if (is.null(lower.xaxis)){
                  localAxis(side=1,...)
               }else{
                  do.call( what = "localAxis"
                         , args = modifyList(x=lower.xaxis,val=list(side=1))
                         )#end do.call
               }#end if (is.null(upper.xaxis))
               #---------------------------------------------------------------------------#
            }#end if
            if (draw.y.upper){
               #----- Upper x axis. -------------------------------------------------------#
               sss = ifelse(row1attop,4,2)
               if (is.null(upper.yaxis)){
                  localAxis(side=sss,...)
               }else{
                  do.call( what = "localAxis"
                         , args = modifyList(x=upper.yaxis,val=list(side=sss))
                         )#end do.call
               }#end if (is.null(upper.xaxis))
               #---------------------------------------------------------------------------#
            }#end if
            if (draw.y.lower){
               #----- Upper x axis. -------------------------------------------------------#
               sss = ifelse(row1attop,2,4)
               if (is.null(upper.yaxis)){
                  localAxis(side=sss,...)
               }else{
                  do.call( what = "localAxis"
                         , args = modifyList(x=lower.yaxis,val=list(side=sss))
                         )#end do.call
               }#end if (is.null(upper.xaxis))
               #---------------------------------------------------------------------------#
            }#end if
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Plot the stuff.                                                         #
            #------------------------------------------------------------------------------#
            mfg = par("mfg")
            if (rr == cc) {
               #----- Plot the diagonal plot. ---------------------------------------------#
               if (has.diag) localDiagPanel(as.vector(x[, rr]), ...)
               #---------------------------------------------------------------------------#

               #----- Plot labels. --------------------------------------------------------#
               if (doText){
                  par(usr=c(0,1,0,1))
                  #----- Fix labels so they don't overflow. -------------------------------#
                  if (is.null(cex.diag)){
                    l.wid      = strwidth(diag.labels, "user")
                    cex.diag = max(0.8, min(2, 0.9/max(l.wid)))
                  }#end if(is.null(cex.diag))
                  #------------------------------------------------------------------------#


                  #------------------------------------------------------------------------#
                  #     Find the position for text.                                        #
                  #------------------------------------------------------------------------#
                  xlp = if (xlog){10^0.5}else{0.5}
                  ylp = if (ylog){10^label.pos}else{label.pos}
                  #------------------------------------------------------------------------#


                  #------------------------------------------------------------------------#
                  #     Plot text.                                                         #
                  #------------------------------------------------------------------------#
                  text.panel( x      = xlp
                            , y      = ylp
                            , txt    = diag.labels[rr]
                            , cex    = cex.diag
                            , font   = font.diag
                            )#end text.panel
                  #------------------------------------------------------------------------#
               }#end if (doText)
               #---------------------------------------------------------------------------#
            }else if (rr < cc){
               #----- Plot lower panel. ---------------------------------------------------#
               localLowerPanel(as.vector(x[, cc]),as.vector(x[,rr]),...)
               #---------------------------------------------------------------------------#
            }else{
               localUpperPanel(as.vector(x[, cc]),as.vector(x[,rr]),...)
            }#end if (rr == cc)
            #------------------------------------------------------------------------------#

            #----- Only low-level functions are allowed. ----------------------------------#
            if (any(par("mfg") != mfg))  stop("the 'panel' function made a new plot")
            #------------------------------------------------------------------------------#
         }#end if
         #---------------------------------------------------------------------------------#
      }#end for (cc in col.loop)
      #------------------------------------------------------------------------------------#
   }#end for (rr in row.loop)
   #---------------------------------------------------------------------------------------#
   


  #----------------------------------------------------------------------------------------#
  #     Plot axis labels.                                                                  #
  #----------------------------------------------------------------------------------------#
  if (! is.null(main)) {
     if ("font.main" %in% nmdots){ 
        font.main = dots$font.main
     }else{ 
        font.main = par("font.main")
     }#end if ("font.main" %in% nmdots)
     if ("cex.main" %in% nmdots){
        cex.main = dots$cex.main
     }else{
        cex.main = par("cex.main")
     }#end if ("cex.main" %in% nmdots)
     mtext( text  = main
          , side  = 3
          , line  = line.main + (1 - is.null(upper.xlab)) * line.xlab
          , outer = TRUE
          , at    = 0.5
          , cex   = cex.main
          , font  = font.main
          )#end mtext
  }#end if (! is.null(main))
  #----------------------------------------------------------------------------------------#
   


  #----------------------------------------------------------------------------------------#
  #     Plot Lower X-axis labels.                                                          #
  #----------------------------------------------------------------------------------------#
  if (! is.null(lower.xlab)){
     if ("font.labels" %in% nmdots){ 
        font.labels = dots$font.labels
     }else{ 
        font.labels = par("font.labels")
     }#end if ("font.labels" %in% nmdots)
     if ("cex.labels" %in% nmdots){
        cex.labels = dots$cex.labels
     }else{
        cex.labels = par("cex.labels")
     }#end if ("cex.labels" %in% nmdots)
     mtext( text  = lower.xlab
          , side  = 1
          , line  = line.xlab
          , outer = TRUE
          , at    = 0.5
          , cex   = cex.labels
          , font  = font.labels
          )#end mtext
  }#end if (! is.null(lower.xlab))
  #----------------------------------------------------------------------------------------#
   


  #----------------------------------------------------------------------------------------#
  #     Plot Upper X-axis labels.                                                          #
  #----------------------------------------------------------------------------------------#
  if (! is.null(upper.xlab)){
     if ("font.labels" %in% nmdots){ 
        font.labels = dots$font.labels
     }else{ 
        font.labels = par("font.labels")
     }#end if ("font.labels" %in% nmdots)
     if ("cex.labels" %in% nmdots){
        cex.labels = dots$cex.labels
     }else{
        cex.labels = par("cex.labels")
     }#end if ("cex.labels" %in% nmdots)
     mtext( text  = upper.xlab
          , side  = 3
          , line  = line.xlab
          , outer = TRUE
          , at    = 0.5
          , cex   = cex.labels
          , font  = font.labels
          )#end mtext
  }#end if (! is.null(labels))
  #----------------------------------------------------------------------------------------#
   


  #----------------------------------------------------------------------------------------#
  #     Plot axis labels.                                                                  #
  #----------------------------------------------------------------------------------------#
  if (! is.null(left.ylab)){
     if ("font.labels" %in% nmdots){ 
        font.labels = dots$font.labels
     }else{ 
        font.labels = par("font.labels")
     }#end if ("font.labels" %in% nmdots)
     if ("cex.labels" %in% nmdots){
        cex.labels = dots$cex.labels
     }else{
        cex.labels = par("cex.labels")
     }#end if ("cex.labels" %in% nmdots)
     mtext( text  = left.ylab
          , side  = 2
          , line  = line.ylab
          , outer = TRUE
          , at    = ifelse(row1attop,(ncolx-1)/(2*ncolx),(ncolx+1)/(2*ncolx))
          , cex   = cex.labels
          , font  = font.labels
          )#end mtext
  }#end if (! is.null(labels))
  #----------------------------------------------------------------------------------------#
   


  #----------------------------------------------------------------------------------------#
  #     Plot axis labels.                                                                  #
  #----------------------------------------------------------------------------------------#
  if (! is.null(right.ylab)){
     if ("font.labels" %in% nmdots){ 
        font.labels = dots$font.labels
     }else{ 
        font.labels = par("font.labels")
     }#end if ("font.labels" %in% nmdots)
     if ("cex.labels" %in% nmdots){
        cex.labels = dots$cex.labels
     }else{
        cex.labels = par("cex.labels")
     }#end if ("cex.labels" %in% nmdots)
     mtext( text  = right.ylab
          , side  = 4
          , line  = line.ylab
          , outer = TRUE
          , at    = ifelse(row1attop,(ncolx+1)/(2*ncolx),(ncolx-1)/(2*ncolx))
          , cex   = cex.labels
          , font  = font.labels
          )#end mtext
  }#end if (! is.null(labels))
  #----------------------------------------------------------------------------------------#



  invisible(NULL)
}#end function(paired.plot)
#==========================================================================================#
#==========================================================================================#
