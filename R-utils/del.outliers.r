#==========================================================================================#
#==========================================================================================#
#      This function deletes outliers from a time series.  It uses the time to find        #
# statistics as a function of the time of the day.                                         #
#------------------------------------------------------------------------------------------#
del.outliers <<- function( x                  # Vector to be evaluated
                         , when               # Time vector
                         , out.hour   = TRUE  # Check mean diel
                         , out.all    = TRUE  # Check time series
                         , ncheck.min = 100   # Minimum vector size to check
                         , bw         = 2     # One-sided bandwidth for spike check
                         , spike.min  = 0.25  # Minimum difference to be considered a spike
                         , f.unreal   = 0.25  # Factor above outlier limit to be considered
                                              #    unrealistic (to be removed regardless of
                                              #    whether it is a spike or not).
                         ){

   #----- Copy x to a local variable. -----------------------------------------------------#
   thisvar = x
   nx      = length(x)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Check whether there is any valid dataset in the input data.  If not, there is    #
   # nothing to be done.                                                                   #
   #---------------------------------------------------------------------------------------#
   if (! any(is.finite(thisvar))) return(thisvar)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Decide whether to discard outliers based on the hour of the day.                   #
   #---------------------------------------------------------------------------------------#
   if (out.hour){
      cat0("     * Diel-based:")
      #------------------------------------------------------------------------------------#
      #    First check: we discard the instantaneous data that are considered weird (i.e.  #
      # normalised variable that is unacceptably far from 0 or a spike that is far from 0  #
      # and surrounded by reasonable values.  We keep iterating it until we have no more   #
      # points removed.                                                                    #
      #------------------------------------------------------------------------------------#
      hh      = hours  (when)
      mm      = minutes(when)
      hhmm    = paste0(sprintf("%2.2i",hh),sprintf("%2.2i",mm))
      n       = 0
      nremain = sum(is.finite(thisvar))
      iterate = nremain >= ncheck.min
      while (iterate){
         n = n + 1

         #---------------------------------------------------------------------------------#
         #     Find the normalised values that correspond to the extreme values (minimum   #
         # and maximum) that would be normal if it happenned only once and the variable    #
         # had a normal distribution.                                                      #
         #---------------------------------------------------------------------------------#
         valid       = length(is.finite(thisvar))
         pfine       = (valid - 2) / valid
         max.fine    = max(3.0,qnorm(pfine, mean = 0., sd = 1.0))
         max.real    = (1.+f.unreal) * max.fine
         #---------------------------------------------------------------------------------#


         #----- Find the mean diurnal cycle and the mean variability of the diel. ---------#
         snstats.dcycle  = tapply(X=thisvar,INDEX=hhmm,FUN=sn.stats    ,na.rm=TRUE)
         snstats.dcycle  = as.data.frame(do.call(what=rbind,args=snstats.dcycle))
         location.dcycle = snstats.dcycle$location
         scale.dcycle    = snstats.dcycle$scale
         shape.dcycle    = snstats.dcycle$shape
         #---------------------------------------------------------------------------------#



         #----- Match the full time series with the average hour. -------------------------#
         unique.hhmm = rownames(snstats.dcycle)
         idx         = match(hhmm,unique.hhmm)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Find the normalised variable.                                               #
         #---------------------------------------------------------------------------------#
         thisnorm    = skew2normal( x        = thisvar 
                                  , location = location.dcycle
                                  , scale    = scale.dcycle
                                  , shape    = shape.dcycle
                                  , idx      = idx
                                  )#end skew2normal
         NEIGH.NORM  = neighbour.mat( x      = thisnorm
                                    , bw     = bw
                                    , cyclic = FALSE
                                    )#end neighbour.mat
         neigh.norm = apply(X=NEIGH.NORM,MARGIN=1,FUN=mean,na.rm=TRUE)
         neigh.norm = ifelse(test=is.finite(neigh.norm),neigh.norm,NA)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Discard suspicious data.                                                    #
         #---------------------------------------------------------------------------------#
         is.infty            = is.infinite(thisnorm)
         is.unreal           = abs(thisnorm  )          %gt% max.real
         is.outlier          = abs(thisnorm  )          %gt% max.fine
         is.spike            = abs(thisnorm-neigh.norm) %ge% spike.min
         weird               = is.infty | is.unreal | (is.outlier & is.spike)
         thisvar[weird]      = NA_real_
         nweird              = sum(weird)
         nremain             = sum(is.finite(thisvar))
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Decide whether we should continue filtering.                               #
         #---------------------------------------------------------------------------------#
         iterate             = (nweird > 0) & (nremain >= ncheck.min)
         cat0("       > Iteration : ",n
             ,";   # of weird hours: ",nweird
             ,";   # of valid hours: ",nremain
             ,";   max.fine = ",sprintf("%.2f",max.fine),"."
             )#end cat0
         #---------------------------------------------------------------------------------#

      }#end while
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Decide whether to discard outliers based on the full time series (to detect weird  #
   # measurements that weren't spikes).                                                    #
   #---------------------------------------------------------------------------------------#
   if (out.all){
      #------------------------------------------------------------------------------------#
      #    Second check: we discard the instantaneous data that are considered weird (i.e. #
      # normalised variable that is unacceptably far from 0 or a spike that is far from 0  #
      # and surrounded by reasonable values.  We keep iterating it until we have no points #
      # considered outliers.                                                               #
      #------------------------------------------------------------------------------------#
      cat0("     * Time-series based:")
      today   = dates(when)
      n       = 0
      nremain = sum(is.finite(thisvar))
      iterate = nremain >= ncheck.min
      while (iterate){
         n = n + 1

         #---------------------------------------------------------------------------------#
         #     Find the normalised values that correspond to the extreme values            #
         # (minimum and maximum) that would be normal if it happenned only once and the    #
         # variable had a normal distribution.                                             #
         #---------------------------------------------------------------------------------#
         valid       = length(is.finite(thisvar))
         pfine       = (valid - 2) / valid
         max.fine    = max(3.0,qnorm(pfine, mean = 0., sd = 1.0))
         max.real    = (1.+f.unreal) * max.fine
         #---------------------------------------------------------------------------------#


         #----- Find the mean diurnal cycle and the mean variability of the diel. ---------#
         snstats.all  = sn.stats(x=thisvar,na.rm=TRUE)
         location.all = snstats.all["location"]
         scale.all    = snstats.all["scale"   ]
         shape.all    = snstats.all["shape"   ]
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Find the normalised variable.                                               #
         #---------------------------------------------------------------------------------#
         thisnorm    = skew2normal( x        = thisvar 
                                  , location = location.all
                                  , scale    = scale.all
                                  , shape    = shape.all
                                  )#end skew2normal
         NEIGH.NORM  = neighbour.mat( x      = thisnorm
                                    , bw     = bw
                                    , cyclic = FALSE
                                    )#end neighbour.mat
         neigh.norm = apply(X=NEIGH.NORM,MARGIN=1,FUN=mean,na.rm=TRUE)
         neigh.norm = ifelse(test=is.finite(neigh.norm),neigh.norm,NA)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Discard suspicious data.                                                    #
         #---------------------------------------------------------------------------------#
         is.infty            = is.infinite(thisnorm)
         is.unreal           = abs(thisnorm  )          %gt% max.real
         is.outlier          = abs(thisnorm  )          %gt% max.fine
         is.spike            = abs(thisnorm-neigh.norm) %ge% spike.min
         weird               = is.infty | is.unreal | (is.outlier & is.spike)
         thisvar[weird]      = NA_real_
         nweird              = sum(weird)
         nremain             = sum(is.finite(thisvar))
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Decide whether we should continue filtering.                               #
         #---------------------------------------------------------------------------------#
         iterate             = (nweird > 0) & (nremain >= ncheck.min)
         cat0("       > Iteration : ",n
             ,";   # of weird hours: ",nweird
             ,";   # of valid hours: ",nremain
             ,";   max.fine = ",sprintf("%.2f",max.fine),"."
             )#end cat0
         #---------------------------------------------------------------------------------#
      }#end while (iterate)
      #------------------------------------------------------------------------------------#
   }#end if (out.all)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Return the clean time series.                                                     #
   #---------------------------------------------------------------------------------------#
   return(thisvar)
   #---------------------------------------------------------------------------------------#
}#end function del.outliers
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#      This function creates a matrix in which each row 'r' has values of vector x near    #
# the element 'r' of vector x, but excluding element 'r' itself.  The number of neighbours #
# is defined by bandwidth variable bw.  If cyclic is TRUE, it assumes vector x to be       #
# cyclic, otherwise it puts NA in the elements near the edge.                              #
#------------------------------------------------------------------------------------------#
neighbour.mat <<- function(x,bw=1,cyclic=FALSE){

   #----- Find the vector size. -----------------------------------------------------------#
   nx  = length(x)
   i   = sequence(nx)
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #      Find the vector with datum.                                                      #
   #---------------------------------------------------------------------------------------#
   idat = rep(i,times=2*bw+2)
   iuse = sequence((nx+1)*(2*bw+1))
   I    = matrix(idat[iuse],nrow=nx+1,ncol=2*bw+1)[-(nx+1),,drop=FALSE]
   I    = 1 + (I-bw-1)%%nx
   X    = matrix(data=x[I],nrow=nx,ncol=2*bw+1)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     In case result should not be cyclic, reset elements to NA.                        #
   #---------------------------------------------------------------------------------------#
   if (! cyclic){
      bye1   = matrix((col(I)+row(I)-1) <= bw,nrow=nx,ncol=2*bw+1)
      bye2   = bye1[rev(sequence(nx)),rev(sequence(2*bw+1))]
      X      = ifelse(test=bye1|bye2,yes=NA,no=X)
   }#end if (! cyclic)
   #---------------------------------------------------------------------------------------#

   #----- Remove middle column: it is the index itself. -----------------------------------#
   X    = X[,-(bw+1),drop=FALSE]
   #---------------------------------------------------------------------------------------#
   return(X)
}#end function neighbour.mat
#==========================================================================================#
#==========================================================================================#
