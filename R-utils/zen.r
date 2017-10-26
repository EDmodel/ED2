#==========================================================================================#
#==========================================================================================#
#    This function finds the cosine of the zenith angle either for the right instant, or   #
# to the interval between two consecutive times.                                           #
#                                                                                          #
# Input variables:                                                                         #
#    - lon       - longitude of the point.  Mandatory, one value only.                     #
#    - lat       - latitude of the point.  Mandatory, one value only.                      #
#    - when      - time.  Mandatory, one point only or a vector.                           #
#    - ed21      - I shall use ED-2.1 method (TRUE/FALSE).  Default is TRUE                #
#    - zeronight - The cosine of zenith angle shall be set to zero at night.               #
#                  Default is FALSE                                                        #
#    - meanval   - I shall find the mean cosine of the integration time.  The beginning    #
#                  and the end are given by variable imetavg.  Default is FALSE.  In case  #
#                  it is TRUE but "when" has just one point, this flag will be ignored and #
#                  it will be solved as instantaneous.                                     #
#    - imetavg   - Which kind of time average was used?                                    #
#                  1 - averages ending at the reference time;                              #
#                  2 - averages beginning at the reference time;                           #
#                  3 - averages centred at the reference time.                             #
#    - nmean     - Number of intermediate points for the average                           #
#                                                                                          #
# The output is going to be a list with the following values:                              #
#    - cosz      - Cosine of zenith angle                                                  #
#    - zen       - The zenith angle in degrees                                             #
#    - height    - The sun height in degrees                                               #
#    - declin    - Declination in degrees                                                  #
#    - day       - Daytime (Sun above horizon)                                             #
#    - night     - Night time (Sun below 6 degrees below the horizon                       #
#    (N.B. When both day and night are false, we consider it twilight.                     #
#------------------------------------------------------------------------------------------#
ed.zen <<- function (lon,lat,when,ed21=TRUE,zeronight=FALSE,meanval=FALSE,imetavg=1
                    ,nmean=120,...){
   #------ Constants. ---------------------------------------------------------------------#
   dcoeff   = c( 0.006918, -0.399912,  0.070257, -0.006758,  0.000907, -0.002697,  0.001480)
   #---------------------------------------------------------------------------------------#

   #------ Find the number of elements. ---------------------------------------------------#
   ntimes  = length(when)
   if ((! meanval) | ntimes == 1) nmean = 1
   #---------------------------------------------------------------------------------------#


   #------ Make matrix of times to make the results time averages if needed be. -----------#
   if (nmean > 1){
      #------------------------------------------------------------------------------------#
      #     The minimum difference is safer than the mean in case the time series has      #
      # gaps.                                                                              #
      #------------------------------------------------------------------------------------#
      dwhen = diff(as.numeric(when))
      sel   = is.finite(dwhen)
      dwhen = dwhen[sel]
      dwhen = min(dwhen[dwhen > 0])
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    Decide the beginning and ending times depending on imetavg.                     #
      #------------------------------------------------------------------------------------#
      if (imetavg == 1){
         #----- Averages ending at the reference time. ------------------------------------#
         na = 1 - nmean
         nz = 0
      }else if (imetavg == 2){
         #----- Averages starting at the reference time. ----------------------------------#
         na = 0
         nz = nmean - 1
      }else if (imetavg == 3){
         #---------------------------------------------------------------------------------#
         #     Averages centered at the reference time. The initial and ending times na    #
         # and nz will be slightly different depending on whether the number of mean       #
         # points is odd or even.                                                          #
         #---------------------------------------------------------------------------------#
         nz = floor(nmean/2) + 0.5 * ((nmean %% 2) - 1.0)
         na = - nz
      }else{
         cat(" ---> In function ed.zen: imetavg =",imetavg,".","\n",sep="")
         stop ("Invalid imetavg, it must be 1, 2, or 3!")
      }#end if
      #------------------------------------------------------------------------------------#



      #----- Averages ending at the reference time. ---------------------------------------#
      dtidx = seq(from=na,to=nz,by=1) / (nz - na + 1)
      WHEN  = chron( matrix(as.numeric(when),ncol=nmean,nrow=ntimes)
                   + matrix(dtidx,ncol=nmean,nrow=ntimes,byrow=TRUE) * dwhen)
      #------------------------------------------------------------------------------------#
   }else{
      #----- Single time, use only the instantaneous value. -------------------------------#
      WHEN  = matrix(as.numeric(when),ncol=nmean,nrow=ntimes)
      #------------------------------------------------------------------------------------#
   }#end if
   empty = as.numeric(WHEN) * NA
   #---------------------------------------------------------------------------------------#



   #------ Find the day of year, list of leap year times, and sun hour. -------------------#
   doy     = matrix(dayofyear(when)           ,ncol=nmean,nrow=ntimes)
   leap    = matrix(is.leap  (when)           ,ncol=nmean,nrow=ntimes)
   fracday = matrix(hms2frac (as.vector(WHEN)),ncol=nmean,nrow=ntimes)
   sunhr   = (fracday * day.hr + lon / 15. + day.hr) %% day.hr
   #---------------------------------------------------------------------------------------#



   #------ Find the hour angle and its cosine. --------------------------------------------#
   hrangle = 15 * (sunhr - 12) * pio180
   chra    = cos(hrangle)
   #---------------------------------------------------------------------------------------#



   #------ Find the declination
   if (ed21){
      doyfun = empty
      doyfun[!leap] = 2 * pi * (doy[!leap] - shsummer) / 365.
      doyfun[ leap] = 2 * pi * (doy[ leap] - shsummer) / 366.

      declin = capri * cos(doyfun)
   }else{
      doyfun = empty
      doyfun[!leap] = 2 * pi * (doy[!leap] - 1) / 365.
      doyfun[ leap] = 2 * pi * (doy[ leap] - 1) / 366.

      declin = ( dcoeff[1]
               + dcoeff[2] * cos(1.*doyfun) + dcoeff[3] * sin(1.*doyfun)
               + dcoeff[4] * cos(2.*doyfun) + dcoeff[5] * sin(2.*doyfun)
               + dcoeff[6] * cos(3.*doyfun) + dcoeff[7] * sin(3.*doyfun) )
   }#end if
   #---------------------------------------------------------------------------------------#

   #------ Find the cosine and sine of latitude and declination. --------------------------#
   clat = cos(pio180*lat)
   slat = sin(pio180*lat)
   cdec = matrix(cos(declin),ncol=nmean,nrow=ntimes)
   sdec = matrix(sin(declin),ncol=nmean,nrow=ntimes)

   #------ Find the cosine of the zenith angle, the zenith angle, and day/night flag. -----#
   cosz   = rowMeans(slat * sdec + clat * cdec * chra,...)
   zen    = acos(cosz) / pio180
   hgt    = 90. - zen
   declin = declin / pio180
   night  = cosz <  cosz.twilight
   day    = cosz >= cosz.min

   if (zeronight){
      cosz[night] =  0.
      hgt [night] =  0.
      zen [night] = 90.
   }#end if
   ans = data.frame(cosz=cosz,zen=zen,hgt=hgt,declin=declin,day=day,night=night)
   return(ans)
}#end function ed.zen
#==========================================================================================#
#==========================================================================================#
