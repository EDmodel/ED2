#==========================================================================================#
#==========================================================================================#
#      This function deletes outliers from a time series.  It uses the time to find        #
# statistics as a function of the time of the day.                                         #
#------------------------------------------------------------------------------------------#
del.outliers <<- function(x,when,out.hour=TRUE,out.all=TRUE){

   thisvar = x
   nx      = length(x)

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
      cat("     * Hourly data...","\n")
      #------------------------------------------------------------------------------------#
      #    First check: we discard the instantaneous data that are considered weird (i.e.  #
      # normalised variable that is unacceptably far from 0 or a spike that is far from 0  #
      # and surrounded by reasonable values.  We keep iterating it until we have no more   #
      # points removed.                                                                    #
      #------------------------------------------------------------------------------------#
      hh    = hours  (when)
      mm    = minutes(when)
      hhmm  = paste(sprintf("%2.2i",hh),sprintf("%2.2i",mm),sep="")
      n       = 0
      iterate = TRUE
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
         #---------------------------------------------------------------------------------#


         #----- Find the mean diurnal cycle and the mean variability of the diel. ---------#
         location.dcycle = tapply(X=thisvar,INDEX=hhmm,FUN=sn.location ,na.rm=TRUE)
         scale.dcycle    = tapply(X=thisvar,INDEX=hhmm,FUN=sn.scale    ,na.rm=TRUE)
         shape.dcycle    = tapply(X=thisvar,INDEX=hhmm,FUN=sn.shape    ,na.rm=TRUE)
         #---------------------------------------------------------------------------------#



         #----- Match the full time series with the average hour. -------------------------#
         unique.hhmm = names(location.dcycle)
         idx         = match(hhmm,unique.hhmm)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Find the normalised variable.                                               #
         #---------------------------------------------------------------------------------#
         thisnorm      = skew2normal( x        = thisvar 
                                    , location = location.dcycle
                                    , scale    = scale.dcycle
                                    , shape    = shape.dcycle
                                    , idx      = idx
                                    )
         thisnormprev  = c(thisnorm[nx],thisnorm[seq(from=1,to=nx-1,by=1)])
         thisnormnext  = c(thisnorm[seq(from=2,to=nx,by=1)],thisnorm[1])
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Discard suspicious data.                                                    #
         #---------------------------------------------------------------------------------#
         unacceptable        = abs(thisnorm) > max.fine
         weird               = unacceptable
         weird[is.na(weird)] = FALSE
         
         thisvar[weird]      = NA
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Decide whether we should continue filtering.                               #
         #---------------------------------------------------------------------------------#
         iterate             = sum(weird) > 0
         cat("       > Iteration :",n,"# of weird hours:",sum(weird)
                                   ,"; max.fine=",sprintf("%.2f",max.fine),"...","\n")
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
      cat("     * Daily data...","\n")
      today = dates  (when)
      n       = 0
      iterate = TRUE
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
         #---------------------------------------------------------------------------------#


         #----- Find the mean diurnal cycle and the mean variability of the diel. ---------#
         location.all = sn.location (x=thisvar,na.rm=TRUE)
         scale.all    = sn.scale    (x=thisvar,na.rm=TRUE)
         shape.all    = sn.shape    (x=thisvar,na.rm=TRUE)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Find the normalised variable.                                               #
         #---------------------------------------------------------------------------------#
         thisnorm      = skew2normal( x        = thisvar 
                                    , location = location.all
                                    , scale    = scale.all
                                    , shape    = shape.all
                                    )
         thisnormprev  = c(thisnorm[nx],thisnorm[seq(from=1,to=nx-1,by=1)])
         thisnormnext  = c(thisnorm[seq(from=2,to=nx,by=1)],thisnorm[1])
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Discard suspicious data.                                                    #
         #---------------------------------------------------------------------------------#
         unacceptable        = abs(thisnorm) > max.fine
         weird               = unacceptable
         weird[is.na(weird)] = FALSE
         
         thisvar[weird]      = NA
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Decide whether we should continue filtering.                               #
         #---------------------------------------------------------------------------------#
         iterate             = sum(weird) > 0
         #---------------------------------------------------------------------------------#

         cat("       > Iteration :",n,"# of weird days:",sum(weird)
                                   ,"; max.fine=",sprintf("%.2f",max.fine),"...","\n")
      }#end while
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Return the clean time series.                                                     #
   #---------------------------------------------------------------------------------------#
   return(thisvar)
   #---------------------------------------------------------------------------------------#
}#end function
#==========================================================================================#
#==========================================================================================#
