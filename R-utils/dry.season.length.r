#==========================================================================================#
#==========================================================================================#
#     dry.season.length -- This function finds the length of the dry season, using a very  #
#                          simple threshold on a previously smoothed time series.          #
#                                                                                          #
#     Input variables:                                                                     #
#     zmean.rain     -- Filtered and gap-filled rainfall rate.  Daily values, although the #
#                       units must be mm/month                                             #
#     today          -- The dates.  It can go for over a year, but only the first year     #
#                       will be analysed. The first year must be complete                  #
#     dry.thresh     -- Minimum rainfall below which dry season starts.  In mm/month       #
#     wet.thresh     -- Minimum rainfall above which wet season starts.  In mm/month       #
#     min.dsl        -- Minimum dry season length to be considered an additional dry       #
#                       season.  In days.                                                  #
#------------------------------------------------------------------------------------------#
dry.season.length <<- function(zmean.rain,today,dry.thresh=90,wet.thresh=110,min.dsl=15){

   #---------------------------------------------------------------------------------------#
   #      Find a suitable zero.                                                            #
   #---------------------------------------------------------------------------------------#
   yeara               = min(numyears(today))
   zero                = chron(paste(12,31,yeara-1,sep="/"))
   doy                 = as.numeric(today-zero)
   maxbeg              = as.numeric(chron(paste(12,31,yeara,sep="/"))-zero)
   first               = numyears(today) == yeara
   rain.first          = zmean.rain
   rain.first[! first] = NA
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      The dry season must start in the current year.  Otherwise the year has no dry    #
   # season.  Also, there must be at least one "wet" day, otherwise, the year has a        #
   # 12-month long dry season.                                                             #
   #---------------------------------------------------------------------------------------#
   wet     = which(zmean.rain > wet.thresh)
   dry     = which(zmean.rain < dry.thresh)
   memory  = 0
   if (length(wet) > 0){
      if (any(dry <= wet[1]) && wet[1] <= maxbeg) memory = wet[1] - 1
      dry = dry[dry > wet[1]]
   }#end if
   sel.1st = dry <= maxbeg
   if (length(wet) == 0 || wet[1] > maxbeg){
      #------------------------------------------------------------------------------------#
      #      Not a single day qualify as wet season.  The entire year is dry season.       #
      #------------------------------------------------------------------------------------#
      doy.beg = min(doy[first])
      doy.end = max(doy[first])
      #------------------------------------------------------------------------------------#
   }else if (length(dry[sel.1st]) == 0){
      #------------------------------------------------------------------------------------#
      #      Not a single day qualify as dry season.  Set the beginning of the dry season  #
      # as the driest day of the year, but make the length zero.                           #
      #------------------------------------------------------------------------------------#
      doy.beg = doy[which.min(rain.first)]
      doy.end = doy[which.min(rain.first)]
   }else if (all(dry[sel.1st] < min(wet[1])) ){
      #------------------------------------------------------------------------------------#
      #      There is a trailing dry season from previous year, but no dry season starts   #
      # in this year.                                                                      #
      #------------------------------------------------------------------------------------#
      sel           = doy < min(wet[1])
      rain.use      = rain.first
      rain.use[sel] = NA
      doy.beg       = doy[which.min(rain.use)]
      doy.end       = doy[which.min(rain.use)]
   }else{
      #------------------------------------------------------------------------------------#
      #      There is a dry season that is not 12-month long that starts in the first      #
      # year.  Now we check whether there is any trailing dry season from previous year.   #
      # If there is, then we remove the trailing dry season and start looking for days     #
      # past the wet season.  In case more than one dry season occurs in the same year, we #
      # pick all of them.                                                                  #
      #------------------------------------------------------------------------------------#
      doy.beg    = NULL
      doy.end    = NULL
      tt         = 0
      iterate    = TRUE
      while ((dry[1] <= maxbeg) && (length(dry) > 0)){
         tt     = tt + 1
         ndry   = length(dry)

         #----- Crop the trailing wet season. ---------------------------------------------#
         wet         = wet[wet > dry[1]]
         #---------------------------------------------------------------------------------#


         if (length(wet) == 0){
            #----- There isn't any wet season in site, make everything dry season. --------#
            doy.beg.now = dry[1]
            doy.end.now = max(dry)
            #------------------------------------------------------------------------------#
         }else{
            #----- Dry season is the period before the next wet season begins. ------------#
            doy.beg.now = dry[1]
            doy.end.now = wet[1]-1
            #------------------------------------------------------------------------------#
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Append the dry season to the list of dry seasons.                           #
         #---------------------------------------------------------------------------------#
         doy.beg = c(doy.beg,doy.beg.now)
         doy.end = c(doy.end,doy.end.now)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Delete the dry season points that were already used.                        #
         #---------------------------------------------------------------------------------#
         dry = dry[dry > doy.end.now]
         #---------------------------------------------------------------------------------#
      }#end while
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     If the last day of the dry season is the last day of all, then the next year is a #
   # completely dry year.  We trim the day so it doesn't go all the way to the other year, #
   # because we don't want to double count.                                                #
   #---------------------------------------------------------------------------------------#
   if (doy.end[length(doy.end)] == max(doy)) doy.end[length(doy.end)] = maxbeg
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the dry season length in days, and remove those that are too short.   If all #
   # of them are short, then we keep only one, the longest one.                            #
   #---------------------------------------------------------------------------------------#
   dsl  = pmin(doy.end,maxbeg+1) - doy.beg
   keep = which(dsl > min.dsl)
   if (any(keep)){
      doy.beg      = doy.beg[keep]
      doy.end      = doy.end[keep]
      dsl          = dsl    [keep]
      dry.beg      = chron(zero + doy.beg)
      dry.end      = chron(zero + doy.end)
      imx          = which.max(dsl)
      long.doy.beg = doy.beg[imx]
      long.doy.end = doy.end[imx]
      long.dsl     = dsl    [imx]
      long.dry.beg = dry.beg[imx]
      long.dry.end = dry.end[imx]
      dsl.total    = sum(dsl) + memory
   }else{
      doy.beg      = NA
      doy.end      = NA
      dsl          = 0
      dry.beg      = NA
      dry.end      = NA
      long.doy.beg = NA
      long.doy.end = NA
      long.dsl     = 0
      long.dry.beg = NA
      long.dry.end = NA
      dsl.total    = memory
   }#end if
   #---------------------------------------------------------------------------------------#


   ans = list( doy.beg      = doy.beg
             , doy.end      = doy.end
             , dry.beg      = dry.beg
             , dry.end      = dry.end
             , dsl          = dsl
             , long.doy.beg = long.doy.beg
             , long.doy.end = long.doy.end
             , long.dry.beg = long.dry.beg
             , long.dry.end = long.dry.end
             , long.dsl     = long.dsl
             , dsl.total    = dsl.total
             )#end list
   return(ans)
}#end function dry.season.length
#==========================================================================================#
#==========================================================================================#
