#==========================================================================================#
#==========================================================================================#
#     This function deletes some bad radiation data.  The input must be a data.frame or    #
# list with all the radiation components already with standard name.  The output will be   #
# the same structure, but with a nicer and cleaner radiation.                              #
#------------------------------------------------------------------------------------------#
del.bad.rshort <<- function(dat,rshort.day.min=0.,par.frac.min=0.80,alb.max=0.30){
   #----- Remove suspicious incoming shortwave radiation. ---------------------------------#
   cat0("   - Remove negative daytime shortwave radiation.")
   suspect             = as.integer(dat$daytime & dat$rshort.in  < rshort.day.min)
   testdays            = dates(sort(unique(dat$today)))
   weird               = testdays[tapply(X=suspect,INDEX=dat$today,FUN=sum,na.rm=TRUE) > 0]
   del                 = dat$today %in% weird
   dat$rshort.in [del] = NA
   #---------------------------------------------------------------------------------------#



   #----- Remove suspicious outgoing shortwave radiation. ---------------------------------#
   cat0("   - Remove suspicious outgoing shortwave radiation.")
   suspect             = as.integer(dat$daytime & dat$rshort.out  < rshort.day.min)
   testdays            = dates(sort(unique(dat$today)))
   weird               = testdays[tapply(X=suspect,INDEX=dat$today,FUN=sum,na.rm=TRUE) > 0]
   del                 = dat$today %in% weird
   dat$rshort.out[del] = NA
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Remove nocturnal data so we can check the diurnal cycle statistics.              #
   #---------------------------------------------------------------------------------------#
   dat$rshort.in [dat$nighttime] = NA
   dat$rshort.out[dat$nighttime] = NA
   dat$par.in    [dat$nighttime] = NA
   dat$par.out   [dat$nighttime] = NA
   #---------------------------------------------------------------------------------------#



   #----- Find the mean diurnal cycle and the mean variability of the diel. ---------------#
   cat0("   - Find probability of each measurement.")
   for (this in c("rshort.in","rshort.out","par.in","par.out")){
      cat0("    > ",this,".")
      hhmm = paste0(sprintf("%2.2i",dat$hour),sprintf("%2.2i",dat$minu))
      #----- Find the mean diurnal cycle. -------------------------------------------------#
      location.this = tapply(X=dat[[this]] ,INDEX=hhmm,FUN=sn.location ,na.rm=TRUE)
      scale.this    = tapply(X=dat[[this]] ,INDEX=hhmm,FUN=sn.scale    ,na.rm=TRUE)
      shape.this    = tapply(X=dat[[this]] ,INDEX=hhmm,FUN=sn.shape    ,na.rm=TRUE)
      #------------------------------------------------------------------------------------#


      #----- Match the full time series with the average hour. ----------------------------#
      unique.hhmm   = names(location.this)
      idx           = match(hhmm,unique.hhmm)
      #------------------------------------------------------------------------------------#


      #----- Normalise the data. ----------------------------------------------------------#
      norm.this     = skew2normal( x        = dat[[this]]
                                 , location = location.this
                                 , scale    = scale.this
                                 , shape    = shape.this
                                 , idx      = idx
                                 )#skew2normal
      prob.this     = dnorm(norm.this)
      assign(paste("prob",this,sep="."),prob.this)
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      PAR should be always close to 45-55% of the total radiation.  Here check whether #
   # the fraction goes outside the range by a large amount (33.3-66.7%).  For these times, #
   # we discard the value with the least probability that we found using a skewed normal   #
   # distribution for the hour of the day.                                                 #
   #---------------------------------------------------------------------------------------#
   cat0("   - Discard data that has weird PAR/SW ratio.")
   weird.in            = ( dat$par.in  < onethird  * dat$rshort.in
                         | dat$par.in  > twothirds * dat$rshort.in  )
   weird.out           = ( dat$par.out < onethird  * dat$rshort.out
                         | dat$par.out > twothirds * dat$rshort.out )
   #----- Delete suspicious incoming shortwave radiation. ---------------------------------#
   del                 = weird.in & prob.rshort.in < prob.par.in
   del[is.na(del)]     = FALSE
   dat$rshort.in[del]  = NA
   #----- Delete suspicious incoming photosynthetically active radiation. -----------------#
   del                 = weird.in & prob.par.in    < prob.rshort.in
   del[is.na(del)]     = FALSE
   dat$par.in[del]     = NA
   #----- Delete suspicious outgoing shortwave radiation. ---------------------------------#
   del                 = weird.out & prob.rshort.out < prob.par.out
   del[is.na(del)]     = FALSE
   dat$rshort.out[del] = NA
   #----- Delete suspicious outgoing photosynthetically active radiation. -----------------#
   del                 = weird.out & prob.par.out    < prob.rshort.out
   del[is.na(del)]     = FALSE
   dat$par.out[del]    = NA
   #---------------------------------------------------------------------------------------#



   #----- Nighttime radiation should be zero. ---------------------------------------------#
   cat0("   - Ensure night time irradiance components are all zero.")
   dat$rshort.in [dat$nighttime] = 0.
   dat$rshort.out[dat$nighttime] = 0.
   dat$par.in    [dat$nighttime] = 0.
   dat$par.out   [dat$nighttime] = 0.
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     We don't let radiation exceed maximum.                                            #
   #---------------------------------------------------------------------------------------#
   cat ("   - Ensure daytime radiation does not exceed maximum.")
   dat$rshort.in[dat$daytime] = pmin(dat$rshort.in[dat$daytime],dat$rshort.pot[dat$daytime])
   dat$par.in   [dat$daytime] = pmin(dat$par.in   [dat$daytime],dat$par.pot   [dat$daytime])
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   return(dat)
   #---------------------------------------------------------------------------------------#
}#end function del.bad.rshort
#==========================================================================================#
#==========================================================================================#
