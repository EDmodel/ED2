#==========================================================================================#
#==========================================================================================#
#     This function creates a Gantt chart.                                                 #
#------------------------------------------------------------------------------------------#
gantt <<- function( task
                  , start         = NULL
                  , end           = NULL
                  , links         = NULL
                  , main          = ""
                  , xlab          = "Time frame"
                  , ylab          = ""
                  , arr.off       = 0.05
                  , y.off         = 0.2
                  , mar           = c(4.1,12.1,3.1,0.1)
                  , cex.task      = NULL
                  , xgrid         = NULL
                  , grid.options  = NULL
                  , title.options = list()
                  , axis.options  = list()
                  , arrow.options = list()
                  , ...
                  ){


   #----- Grab info from data frame if only one object has been given. --------------------#
   if (is.null(start) && is.null(end)){
      start  = task$start
      end    = task$end
      links  = task$links
      task   = task$task
   }#end if
   #---------------------------------------------------------------------------------------#



   #------ Find a margin that fits all information. ---------------------------------------#
   if (is.null(cex.task)){
      width    = max(nchar(task))
      mleft    = mar[2] - 0.1
      cex.task = mleft * 18 / ( width * 12 )
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Number of tasks. ----------------------------------------------------------------#
   ntasks = length(task)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Check what type of object links is.  Standardise it to list as tasks may link to #
   # multiple tasks.                                                                       #
   #---------------------------------------------------------------------------------------#
   if (is.null(links)){
      links = sapply(X=rep(NA,ntasks),FUN=list)
   }else if(is.vector(links)){
      links = sapply(X=links,FUN=list)
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Find the ranges. ----------------------------------------------------------------#
   xlimit = range(c(start,end))
   ylimit = c(-ntasks,0)
   deltax = arr.off * diff(ylimit)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Make pretty x axes but first check whether the objects are chron or not.         #
   #---------------------------------------------------------------------------------------#
   if (is.null(xgrid)){
      if (is.time(start) || is.time(end)){
         xat     = pretty.time(chron(xlimit))
         xlabels = chron(xat)
      }else{
         xat     = pretty(xlimit)
         xlabels = sprintf("%g",xat)
      }#end if
   }else{
      xat     = xgrid
      if (is.time(xgrid)){
         xlabels = paste(xgrid)
      }else{
         xlabels = sprintf("%g",xat)
      }#end if
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Axes for y.                                                                       #
   #---------------------------------------------------------------------------------------#
   yat     = - sequence(ntasks) + 0.5
   ytop    = - sequence(ntasks) + (1.0-y.off)
   ybottom = - sequence(ntasks) + y.off
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #       Open plotting device.                                                           #
   #---------------------------------------------------------------------------------------#
   par(par.user)
   par(mar=mar)
   plot.new()
   plot.window(xlim=xlimit,ylim=ylimit)
   do.call    ( what = "title"
              , args = modifyList(x=title.options,val=list(main=main,xlab=xlab,ylab=ylab))
              )#end do.call
   do.call    ( what = "axis"
              , args = modifyList(x=axis.options,val=list(side=1,at=xat,labels=xlabels))
              )#end do.call
   do.call    ( what = "axis"
              , args = modifyList( x   = axis.options
                                 , val = list( side     = 2
                                             , las      = 1
                                             , at       = yat
                                             , labels   = task
                                             , cex.axis = cex.task
                                             )
                                 )#end modifyList
              )#end do.call
   if (! is.null(xgrid)){
      do.call( what = "abline", args = modifyList(x=grid.options,val=list(v=xgrid)))
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Plot the boxes.                                                                  #
   #---------------------------------------------------------------------------------------#
   rect.options = modifyList( x   = list(...)
                            , val = list(xleft=start,ybottom=ybottom,xright=end,ytop=ytop)
                            )#end modifyList
   do.call( what = "rect", args = rect.options)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Plot links.                                                                      #
   #---------------------------------------------------------------------------------------#
   neach  = sapply(X=links,FUN=length)
   from   = unlist(mapply(FUN=rep,x=sequence(ntasks),each=neach))
   to     = unlist(links)
   sel    = is.finite(to)
   from   = from[sel]
   to     = to  [sel]
   if (any(sel)){
      #------ Find the beginning and end of arrows. ---------------------------------------#
      wstart = as.numeric(start)
      wend   = as.numeric(end  )
      
      wa.from = wstart[from]
      wz.from = wend  [from]
      wa.to   = wstart[to  ]
      wz.to   = wend  [to  ]
      
      yfrom  = - from + ifelse(wz.to > wz.from, 0.5, y.off)
      yto    = - to + (1.0-y.off)
      xfrom  = ifelse( wz.to > wz.from, wz.from, 0.5*(pmax(wa.from,wa.to)+wz.to) + deltax)
      xto    = ifelse( wz.to > wz.from, pmax(wz.from,wa.to) + deltax, xfrom)
      hey    = cbind( from  = from
                    , to    = to
                    , wa.from = wstart[from]
                    , wz.from = wend  [from]
                    , wa.to   = wstart[to  ]
                    , wz.to   = wend  [to  ]
                    , xfrom   = round(xfrom,1)
                    , xto     = round(xto  ,1)
                    , yfrom   = round(yfrom,1)
                    , yto     = round(yto  ,1)
                    )
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Arguments for segments.                                                        #
      #------------------------------------------------------------------------------------#
      argnow = modifyList( x   = arrow.options
                         , val = list( x0     = xfrom
                                     , y0     = yfrom
                                     , x1     = xto
                                     , y1     = yfrom
                                     , length = NULL
                                     , angle  = NULL
                                     , code   = NULL
                                     , lty    = "solid"
                                     )#end list
                         )#end modifyList
      do.call( what = "segments", args = argnow )
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Arguments for arrows.                                                          #
      #------------------------------------------------------------------------------------#
      argnow = modifyList( x   = arrow.options
                         , val = list( x0     = xto
                                     , y0     = yfrom
                                     , x1     = xto
                                     , y1     = yto
                                     , code   = 2
                                     , lty    = "solid"
                                     )#end list
                         )#end modifyList
      do.call( what = "arrows", args = argnow )
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Plot the boxes.                                                                  #
   #---------------------------------------------------------------------------------------#
   rect.options = modifyList( x   = list(...)
                            , val = list(xleft=start,ybottom=ybottom,xright=end,ytop=ytop)
                            )#end modifyList
   do.call( what = "rect", args = rect.options)
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Plot dashed lines.                                                                #
   #---------------------------------------------------------------------------------------#
   if (any(sel)){
      #------------------------------------------------------------------------------------#
      #     Arguments for segments.                                                        #
      #------------------------------------------------------------------------------------#
      argnow = modifyList( x   = arrow.options
                         , val = list( x0     = xfrom
                                     , y0     = yfrom
                                     , x1     = xto
                                     , y1     = yfrom
                                     , length = NULL
                                     , angle  = NULL
                                     , code   = NULL
                                     , lty    = "dotted"
                                     )#end list
                         )#end modifyList
      do.call( what = "segments", args = argnow )
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Arguments for arrows.                                                          #
      #------------------------------------------------------------------------------------#
      argnow = modifyList( x   = arrow.options
                         , val = list( x0     = xto
                                     , y0     = yfrom
                                     , x1     = xto
                                     , y1     = yto
                                     , length = NULL
                                     , angle  = NULL
                                     , code   = NULL
                                     , lty    = "dotted"
                                     )#end list
                         )#end modifyList
      do.call( what = "segments", args = argnow )
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Draw a box. ---------------------------------------------------------------------#
   box()
   #---------------------------------------------------------------------------------------#
}#end function gantt
#==========================================================================================#
#==========================================================================================#
