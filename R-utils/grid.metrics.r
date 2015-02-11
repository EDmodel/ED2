#==========================================================================================#
#==========================================================================================#
#     This function returns the cloud metrics from a point cloud.                          #
#------------------------------------------------------------------------------------------#
grid.metrics <<- function( x
                         , pixres   = 1
                         , min.pts  = 4 * pixres^2
                         , maxblock = 50
                         , verbose  = FALSE
                         , ...
                         ){

   #---------------------------------------------------------------------------------------#
   #     Make sure some cloud.metrics options are properly set for grid metrics.           #
   #---------------------------------------------------------------------------------------#
   dotdotdot = list(...)
   warn = 0
   #----- mat.out must be FALSE. ----------------------------------------------------------#
   if (! is.null(dotdotdot$mat.out)){
      #----- Warn only in case mat.out is TRUE. -------------------------------------------#
      if (dotdotdot$mat.out){
         warning(" mat.out cannot be true in grid.metrics, setting it to false")
         warn = warn+1
      }#end if
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#

   #----- mhdens must be NULL. ------------------------------------------------------------#
   if (! is.null(dotdotdot$mhdens)){
      warning(" Ignoring mhdens as it varies for each profile")
      warn = warn+1
   }#end if
   #---------------------------------------------------------------------------------------#

   #----- zzdens must be NULL. ------------------------------------------------------------#
   if (! is.null(dotdotdot$zzdens)){
      warning(" Ignoring zzdens as it varies for each profile")
      warn = warn+1
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    If warn is not zero, then call grid metrics with the correct settings.             #
   #---------------------------------------------------------------------------------------#
   if (warn > 0){
      dotdotdot = modifyList( x   = dotdotdot
                            , val = list( x        = x
                                        , pixres   = pixres
                                        , min.pts  = min.pts
                                        , maxblock = maxblock
                                        , verbose  = verbose
                                        , mat.out  = FALSE
                                        , zzdens   = NULL
                                        , mhdens   = NULL
                                        )
                            )#end modifyList
      ans       = do.call(what="grid.metrics",args=dotdotdot)
      return(ans)
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- "Degrade" point cloud resolution to the desired resolution. ---------------------#
   when.round = proc.time()
   x$x = round(x$x / pixres) * pixres
   x$y = round(x$y / pixres) * pixres
   when.round = proc.time() - when.round
   if (verbose) cat("    + Rounding took ",when.round[3],"s... \n",sep="")
   #---------------------------------------------------------------------------------------#


   #----- "Degrade" point cloud resolution to the desired resolution. ---------------------#
   when.index = proc.time()
   xuniq      = sort(unique(x$x))
   yuniq      = sort(unique(x$y))
   nx         = length(xuniq)
   ny         = length(yuniq)
   x$idx      = match(x$x,xuniq) + nx * (match(x$y,yuniq))
   when.index = proc.time() - when.index
   if (verbose) cat("    + Finding indices took ",when.index[3],"s... \n",sep="")
   #---------------------------------------------------------------------------------------#



   #----- Split data frame using 1-D index. -----------------------------------------------#
   when.split1d = proc.time()
   pc.list      = split(x=x,f=x$idx)
   when.split1d = proc.time() - when.split1d
   if (verbose) cat("    + Splitting point cloud ",when.split1d[3],"s... \n",sep="")
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Split list in smaller blocks.                                                     #
   #---------------------------------------------------------------------------------------#
   nlist = length(pc.list)
   nloop = ceiling(nlist/maxblock)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Run the first element to initialise the output matrix.                           #
   #---------------------------------------------------------------------------------------#
   if (verbose) cat("    + Initialising matrix...","\n",sep="")
   if (nlist >= 1){
      dotdotdot    = list(...)
      dotdotdot    = modifyList( x   = dotdotdot
                               , val = list(x=pc.list[[1]],min.pts=min.pts,mat.out=FALSE)
                               )#end modifyList
      template     = do.call(what="cloud.metrics",args=dotdotdot)
      cm.table     = matrix( ncol     = nlist
                           , nrow     = length(template)
                           , dimnames = list(names(template),names(pc.list))
                           )#end matrix
      cm.table[,1] = template
      rm(dotdotdot,template)
   }#end if (nlist >= 1)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Loop over remaining elements using for loop.                                     #
   #---------------------------------------------------------------------------------------#
   if (verbose) cat("    + Aggregating data: ","\n",sep="")
   nshow = sort(c(1,nlist,pretty(c(1,nlist),n=10)))
   nshow = nshow[nshow %wr% c(1,nlist)]
   when.block = proc.time()
   for (n in sequence(nlist)){
      if (verbose && (n %in% nshow)){
         cat("      - Processing pixel ",n,"/",nlist,"... \n",sep=" ")
      }#end if (verbose && (n %in% nshow))
      dotdotdot    = list(...)
      dotdotdot    = modifyList( x   = dotdotdot
                               , val = list(x=pc.list[[n]],min.pts=min.pts,mat.out=FALSE)
                               )#end modifyList
      cm.table[,n] = do.call(what="cloud.metrics",args=dotdotdot)
      rm(dotdotdot)
   }#end for (n in sequence(nloop))
   when.block = proc.time() - when.block
   if (verbose) cat("    + Running cloud metrics took ",when.block[3],"s... \n",sep="")
   #---------------------------------------------------------------------------------------#



   #----- Return the mean, maximum, and the standard deviation of all metrics. ------------#
   when.aggr              = proc.time()
   mean.cm.table          = apply(X=cm.table,MARGIN=1,FUN=mean  ,na.rm=TRUE)
   sdev.cm.table          = apply(X=cm.table,MARGIN=1,FUN=sd    ,na.rm=TRUE)
   median.cm.table        = apply(X=cm.table,MARGIN=1,FUN=median,na.rm=TRUE)
   max.cm.table           = apply(X=cm.table,MARGIN=1,FUN=max   ,na.rm=TRUE)
   names(mean.cm.table  ) = paste("mean"  ,rownames(cm.table),sep=".")
   names(sdev.cm.table  ) = paste("sdev"  ,rownames(cm.table),sep=".")
   names(median.cm.table) = paste("median",rownames(cm.table),sep=".")
   names(max.cm.table   ) = paste("max"   ,rownames(cm.table),sep=".")
   ans                    = c(mean.cm.table,sdev.cm.table,median.cm.table,max.cm.table)
   when.aggr              = proc.time() - when.aggr
   if (verbose) cat("    + Aggregating data took ",when.aggr[3],"s... \n",sep="")
   #---------------------------------------------------------------------------------------#



   #----- Free memory and return table. ---------------------------------------------------#
   rm(cm.table,mean.cm.table,sdev.cm.table,median.cm.table,max.cm.table)
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end grid.metrics
#==========================================================================================#
#==========================================================================================#
