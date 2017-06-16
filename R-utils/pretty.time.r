#==========================================================================================#
#==========================================================================================#
#     This function creates a pretty time scale.  It is loosely based on pretty, but here  #
# we make extensive use of the chron functions, and define the suitable scales in a        #
# different way as time has a non-decimal scale.                                           #
#  The result is a list containing the levels, and nice labels for the plots.              #
#------------------------------------------------------------------------------------------#
pretty.time <<- function(when,n=10,...){

   #----- Find the 1st and last time. -----------------------------------------------------#
   whena = min(when,na.rm=TRUE)
   whenz = max(when,na.rm=TRUE)

   #----- Table for accepted bases for months, hours, minutes, and seconds. ---------------#
   base.months = c(1,2,3,4,6)
   base.days   = c(1,2,4,7,15)
   base.hours  = c(1,2,3,4,6,12)
   base.minsec = c(1,2,5,10,15,20,30)

   #----- Convert time to seconds. --------------------------------------------------------#
   when.sec = as.numeric(when) * day.sec

   #---------------------------------------------------------------------------------------#
   #    Find a first guess of the step size, so we decide whether to use years, months,    #
   # days, hours, minutes, or seconds.                                                     #
   #---------------------------------------------------------------------------------------#
   wstep.1st  = mean(diff(pretty(when.sec,n)))
   if (wstep.1st == 0){
      myunit=NA
      #---- Whatever, use just what comes out from the regular pretty. --------------------#
      vlevels = chron(pretty(when,n))
      vlabels = as.character(vlevels)
      padj    = rep(0,times=length(vlabels))
   }else if(wstep.1st / yr.sec > 0.8){
      myunit="years"
      #------------------------------------------------------------------------------------#
      #     Years are the best scale for the plots.                                        #
      #------------------------------------------------------------------------------------#
      yrrange = numyears(when)
      vlevels = pretty(yrrange,n)
      vlevels = dates(x=paste(1,1,vlevels,sep="/"))

      vlabels = paste(months(vlevels),years(vlevels),sep="-")
      padj    = rep(0,times=length(vlabels))
   }else if(wstep.1st / (30. * day.sec) > 0.8){
      myunit="months"
      #------------------------------------------------------------------------------------#
      #     Months are the best scale for the plots.                                       #
      #------------------------------------------------------------------------------------#
      #----- Find the time step that is the closest to the base. --------------------------#
      wstep     = wstep.1st / (30. * day.sec)
      whichbase = base.months[which.min(abs(wstep-base.months))]

      #----- Find the list of years to plot. ----------------------------------------------#
      allyears = numyears(when)
      yeara    = min(allyears,na.rm=TRUE)
      yearz    = max(allyears,na.rm=TRUE)+1
      vlevels  = seq.dates(from = paste(1,1,yeara,sep="/")
                          ,to   = paste(1,1,yearz,sep="/")
                          ,by   = "months")
      mon1st   = nummonths(vlevels)
      monlevs  = seq(from=1,to=12,by=whichbase)
      #------------------------------------------------------------------------------------#

      #----- Find the limits that will keep the labels not too far from the data. ---------#
      vlevels  = vlevels[mon1st %in% monlevs]
      vlev.a   = which.max(vlevels[vlevels <= whena])
      vlev.z   = length(vlevels) + 1 - which.min(rev(vlevels[vlevels >= whenz]))
      keep     = seq(from=vlev.a,to=vlev.z,by=1)
      #------------------------------------------------------------------------------------#

      #----- Format output. ---------------------------------------------------------------#
      vlevels = dates(vlevels[keep],out.format="m/d/year")
      vlabels = paste(months(vlevels),years(vlevels),sep="-")
      padj    = rep(0,times=length(vlabels))
      #------------------------------------------------------------------------------------#

   }else if(wstep.1st / day.sec > 0.8){
      myunit="days"

      #------------------------------------------------------------------------------------#
      #     Days are the best scale for the plots, but we keep them tethered to months,    #
      # even if the grid becomes slightly irregular.                                       #
      #------------------------------------------------------------------------------------#
      #----- Find the time step that is the closest to the base. --------------------------#
      wstep     = wstep.1st / day.sec
      whichbase = base.days[which.min(abs(wstep-base.days))]
      #----- Find the list of years to plot. ----------------------------------------------#
      allyears = numyears(when)
      yeara    = min(allyears,na.rm=TRUE)
      yearz    = max(allyears,na.rm=TRUE)+1
      #------------------------------------------------------------------------------------#
      #     Impose the list of months to be from January to December, we will trim the     #
      # numbers later.                                                                     #
      #------------------------------------------------------------------------------------#
      montha   = 1
      monthz   = 12
      #------------------------------------------------------------------------------------#
      #     Impose the list of months to be from January to December, we will trim the     #
      # numbers later.                                                                     #
      #------------------------------------------------------------------------------------#
      daylevs=seq(from=1,to=31-whichbase+1,by=whichbase)
      #----- First guess for the levels. --------------------------------------------------#
      vlevels  = chron( seq.dates( from = paste(1,1,yeara,sep="/")
                                 , to   = paste(1,1,yearz,sep="/")
                                 , by   = "days")
                      )#end chron 
      day1st   = numdays(vlevels)
      #------------------------------------------------------------------------------------#

      #----- Find the limits that will keep the labels not too far from the data. ---------#
      vlevels  = vlevels[day1st %in% daylevs]
      vlev.a   = which.max(vlevels[vlevels <= whena])
      vlev.z   = length(vlevels) + 1 - which.min(rev(vlevels[vlevels >= whenz]))
      keep     = seq(from=vlev.a,to=vlev.z,by=1)
      #------------------------------------------------------------------------------------#

      #----- Format output. ---------------------------------------------------------------#
      vlevels  = dates(vlevels[keep],out.format="m/d/y")
      vlabels  = paste(months(vlevels),days(vlevels),sep="/")
      padj     = rep(0,times=length(vlabels))

      #----- First day of the year. -------------------------------------------------------#
      sel          = numdays(vlevels) == 1 & nummonths(vlevels) == 1
      vlabels[sel] = paste(months(vlevels[sel]),"/",days(vlevels[sel]),"\n"
                          ,years(vlevels[sel]),sep="")
      padj   [sel] = 0.5
      #------------------------------------------------------------------------------------#


      #----- First label.  Always include year in it. -------------------------------------#
      vlabels[1] = paste(months(vlevels[1]),"/",days(vlevels[1]),"\n",years(vlevels[1])
                        ,sep="")
      padj   [1] = 0.5
      #------------------------------------------------------------------------------------#

   }else if(wstep.1st / hr.sec > 0.8){
      myunit="hours"
      #------------------------------------------------------------------------------------#
      #     Hours are the best scale for the plots.                                        #
      #------------------------------------------------------------------------------------#
      #----- Find the time step that is the closest to the base. --------------------------#
      wstep     = wstep.1st / hr.sec
      whichbase = base.hours[which.min(abs(wstep-base.hours))]
      #----- Find the list of days to plot. -----------------------------------------------#
      when1st  = dates(min(when  ,na.rm=TRUE))
      whenlast = dates(max(when+1,na.rm=TRUE))
      mydates  = seq.dates(from=when1st,to=whenlast,by="days")
      mytimes  = times(seq(from=0,to=day.sec-1,by=whichbase*hr.sec)) / day.sec
      ndates   = length(mydates)
      ntimes   = length(mytimes)
      #----- First guess for the levels. --------------------------------------------------#
      vlevels  = chron(dates=rep(x=mydates,each=ntimes),times=rep(x=mytimes,times=ndates))
      #------------------------------------------------------------------------------------#

      #----- Find the limits that will keep the labels not too far from the data. ---------#
      vlev.a   = which.max(vlevels[vlevels <= whena])
      vlev.z   = length(vlevels) + 1 - which.min(rev(vlevels[vlevels >= whenz]))
      keep     = seq(from=vlev.a,to=vlev.z,by=1)
      #------------------------------------------------------------------------------------#

      #----- Format output. ---------------------------------------------------------------#
      vlevels      = vlevels[keep]
      vlabels      = paste(sprintf("%2.2i",hours(vlevels))
                          ,sprintf("%2.2i",minutes(vlevels)),sep=":")
      padj         = rep(0,times=length(vlabels))
      #----- First time of the day. -------------------------------------------------------#
      sel          = hours(vlevels) == 0
      vlabels[sel] = paste(sprintf("%2.2i",hours(vlevels[sel])),":"
                          ,sprintf("%2.2i",minutes(vlevels[sel])),"\n"
                          ,months(vlevels[sel]),"-",days(vlevels[sel]),sep="")
      padj[sel]    = 0.5
      #----- First time or first time of the year. ----------------------------------------#
      sel          = ( vlevels == vlevels[1] 
                     |  ( nummonths(vlevels) == 1 & numdays(vlevels) == 1 
                        & hours(vlevels) == 0 ))
      vlabels[sel] = paste(sprintf("%2.2i",hours(vlevels[sel])),":"
                          ,sprintf("%2.2i",minutes(vlevels[sel])),"\n"
                          ,months(vlevels[sel]),"-",days(vlevels[sel]),"\n"
                          ,years(vlevels[sel]),sep="")
      padj[sel]    = 0.5
      #------------------------------------------------------------------------------------#


   }else if(wstep.1st / min.sec > 0.8){
      myunit="minutes"
      #------------------------------------------------------------------------------------#
      #     Minutes are the best scale for the plots.                                      #
      #------------------------------------------------------------------------------------#
      #----- Find the time step that is the closest to the base. --------------------------#
      wstep     = wstep.1st / min.sec
      whichbase = base.minsec[which.min(abs(wstep-base.minsec))]
      #----- Find the list of days to plot. -----------------------------------------------#
      when1st  = dates(min(when  ,na.rm=TRUE))
      whenlast = dates(max(when+1,na.rm=TRUE))
      mydates  = seq.dates(from=when1st,to=whenlast,by="days")
      mytimes  = times(seq(from=0,to=day.sec-1,by=whichbase*min.sec)) / day.sec
      ndates   = length(mydates)
      ntimes   = length(mytimes)
      #----- First guess for the levels. --------------------------------------------------#
      vlevels  = chron(dates=rep(x=mydates,each=ntimes),times=rep(x=mytimes,times=ndates))
      #------------------------------------------------------------------------------------#

      #----- Find the limits that will keep the labels not too far from the data. ---------#
      vlev.a   = which.max(vlevels[vlevels <= whena])
      vlev.z   = length(vlevels) + 1 - which.min(rev(vlevels[vlevels >= whenz]))
      keep     = seq(from=vlev.a,to=vlev.z,by=1)
      #------------------------------------------------------------------------------------#

      #----- Format output. ---------------------------------------------------------------#
      vlevels      = vlevels[keep]
      vlabels      = paste(sprintf("%2.2i",hours(vlevels))
                          ,sprintf("%2.2i",minutes(vlevels)),sep=":")
      padj         = rep(0,times=length(vlabels))
      #----- First time of the day. -------------------------------------------------------#
      sel          = hours(vlevels) == 0 & minutes(vlevels) == 0
      vlabels[sel] = paste(sprintf("%2.2i",hours(vlevels[sel])),":"
                          ,sprintf("%2.2i",minutes(vlevels[sel])),"\n"
                          ,months(vlevels[sel]),"-",days(vlevels[sel]),sep="")
      padj[sel]    = 0.5
      #----- First time of the year. ------------------------------------------------------#
      sel          = ( vlevels == vlevels[1] 
                     |  ( nummonths(vlevels) == 1 & numdays(vlevels) == 1 
                        & hours(vlevels) == 0 & minutes(vlevels) == 0))
      vlabels[sel] = paste(sprintf("%2.2i",hours(vlevels[sel])),":"
                          ,sprintf("%2.2i",minutes(vlevels[sel])),"\n"
                          ,months(vlevels[sel]),"-",days(vlevels[sel]),"\n"
                          ,years(vlevels[sel]),sep="")
      padj[sel]    = 0.5
      #------------------------------------------------------------------------------------#



   }else{
      myunit="seconds"
      #------------------------------------------------------------------------------------#
      #     Seconds are the best scale for the plots.                                      #
      #------------------------------------------------------------------------------------#
      #----- Find the time step that is the closest to the base. --------------------------#
      wstep     = wstep.1st
      whichbase = base.minsec[which.min(abs(wstep-base.minsec))]
      #----- Find the list of days to plot. -----------------------------------------------#
      when1st  = dates(min(when  ,na.rm=TRUE))
      whenlast = dates(max(when+1,na.rm=TRUE))
      mydates  = seq.dates(from=when1st,to=whenlast,by="days")
      mytimes  = times(seq(from=0,to=day.sec-1,by=whichbase)) / day.sec
      ndates   = length(mydates)
      ntimes   = length(mytimes)
      #----- First guess for the levels. --------------------------------------------------#
      vlevels  = chron(dates=rep(x=mydates,each=ntimes),times=rep(x=mytimes,times=ndates))
      #------------------------------------------------------------------------------------#

      #----- Find the limits that will keep the labels not too far from the data. ---------#
      vlev.a   = which.max(vlevels[vlevels <= whena])
      vlev.z   = length(vlevels) + 1 - which.min(rev(vlevels[vlevels >= whenz]))
      keep     = seq(from=vlev.a,to=vlev.z,by=1)
      #------------------------------------------------------------------------------------#

      #----- Format output. ---------------------------------------------------------------#
      vlevels      = vlevels[keep]
      vlabels      = paste(sprintf("%2.2i",hours(vlevels))
                          ,sprintf("%2.2i",minutes(vlevels))
                          ,sprintf("%2.2i",seconds(vlevels)),sep=":")
      padj         = rep(0,times=length(vlabels))
      #----- First time of the day. -------------------------------------------------------#
      sel          = hours(vlevels) == 0 & minutes(vlevels) == 0 & seconds(vlevels) == 0
      vlabels[sel] = paste(sprintf("%2.2i",hours(vlevels[sel])),":"
                          ,sprintf("%2.2i",minutes(vlevels[sel])),":"
                          ,sprintf("%2.2i",seconds(vlevels[sel])),"\n"
                          ,months(vlevels[sel]),"-",days(vlevels[sel]),sep="")
      padj[sel]    = 0.5
      #----- First time of the year. ------------------------------------------------------#
      sel          = ( vlevels == vlevels[1] 
                     |  ( nummonths(vlevels) == 1 & numdays(vlevels) == 1 
                        & hours(vlevels) == 0  & minutes(vlevels) == 0 
                        & seconds(vlevels) == 0))
      vlabels[sel] = paste(sprintf("%2.2i",hours(vlevels[sel])),":"
                          ,sprintf("%2.2i",minutes(vlevels[sel])),":"
                          ,sprintf("%2.2i",seconds(vlevels[sel])),"\n"
                          ,months(vlevels[sel]),"-",days(vlevels[sel]),"\n"
                          ,years(vlevels[sel]),sep="")
      padj[sel]    = 0.5
      #------------------------------------------------------------------------------------#
   }#end if

   vresult=list(levels=vlevels,labels=vlabels,n=length(vlevels),scale=myunit,padj=padj)
   return(vresult)
}#end function
#==========================================================================================#
#==========================================================================================#







#==========================================================================================#
#==========================================================================================#
#      This function finds pretty axis scales for elapsed time, using bases typical of     #
# time (12, 24, 60) instead of decimal.                                                    #
#------------------------------------------------------------------------------------------#
pretty.elapsed <<- function(x,base,n=5,...){

   #------ Make sure the base is known. ---------------------------------------------------#
   if (! base %in% c(12,24,60)){
      stop( paste0(" Invalid base (",base,")!"
                  ,"  It must be 12 (months), 24 (hours), or 60 (minutes/seconds)."
                  )#end paste0
          )#end stop
   }#end if  (! base %in% c(12,24,60))
   #---------------------------------------------------------------------------------------#


   #------ Run instrinsic function pretty. ------------------------------------------------#
   x.base    = x / base
   neat.base = pretty(x.base,n=n,...)
   diff.base = median(diff(neat.base))
   xlwr      = min(x,na.rm=TRUE)
   xupr      = max(x,na.rm=TRUE)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Re-scale data based on the pretty limits.                                        #
   #---------------------------------------------------------------------------------------#
   if (diff.base %<=% 0.075 || diff.base %>=% 0.75){
      dbase = diff.base * base
   }else if (diff.base %<% 0.15){
      dbase = base / 12
   }else if (diff.base %<% 0.35){
      dbase = base / 4
   }else{ 
      dbase = base / 2
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Build sequence. -----------------------------------------------------------------#      
   nlwr = dbase * floor  (xlwr/dbase)
   nupr = dbase * ceiling(xupr/dbase)
   neat = seq(from=nlwr,to=nupr,by=dbase)
   #---------------------------------------------------------------------------------------#

   return(neat)
}#end pretty.elapsed
#==========================================================================================#
#==========================================================================================#
