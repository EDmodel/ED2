#==========================================================================================#
#==========================================================================================#
#     Function that checks whether an object is chron.                                     #
#------------------------------------------------------------------------------------------#
is.dates <<- function(x) inherits(x,"dates")
is.times <<- function(x) inherits(x,"times")
is.time  <<- function(x){ is.chron(x) || is.dates(x) || is.times(x)}
#------------------------------------------------------------------------------------------#


#==========================================================================================#
#==========================================================================================#
#      Function that determines whether the year is leap or not.                           #
#------------------------------------------------------------------------------------------#
is.leap <<- function(when){

   wit = is(when)
   
   if ("dates" %in% wit || "chron" %in% wit){
      year = numyears(when)
   }else{
      year = when
   }#end if
   leaptf = year %% 400 == 0 | (year %% 4 == 0 & year %% 100 != 0)
   return(leaptf)
} #end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      Function that determines the number of days in a given month.                       #
#------------------------------------------------------------------------------------------#
daymax <<- function(month,year){
  mmm  = c(31,28,31,30,31,30,31,31,30,31,30,31)

  if (missing(year)){
     wit = is(month)
     if ("dates" %in% wit || "chron" %in% wit){
        when  = month
        year  = numyears (when)
        month = nummonths(when)
     }else{
        stop("  No year given and month is not time")
     }#end if
  }#end if

  mday         = mmm[month]
  addone       = month == 2 & is.leap(year)
  mday[addone] = mday[addone] + 1

  return(mday)
} #end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      Function that determines the number of the month given the character name.          #
#------------------------------------------------------------------------------------------#
mmm2mon <<- function(mmm,lang="English"){
  lang = substring(tolower(lang),1,2)
  if (lang %in% c("en")){
     m3l  = tolower(month.abb)
  }else if(lang %in% c("po","pt")){
     m3l  = c("jan","fev","mar","abr","mai","jun","jul","ago","set","out","nov","dez")
  }else if(lang %in% c("fr")){
     m3l  = c("jan","fev","mar","avr","mai","jun","jul","aou","sep","oct","nov","dec")
  }#end if

  mmmloc = tolower(substring(as.character(mmm),1,3))
  monout = match(mmmloc,m3l)
  return(monout)
} #end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      Function that determines 3-letter name of the month given their number.             #
#------------------------------------------------------------------------------------------#
mon2mmm <<- function(mon,lang="English",cap1=FALSE){
  lang = substring(tolower(lang),1,2)
  if (lang %in% c("en")){
     m3l  = month.abb
  }else if(lang %in% c("po","pt")){
     m3l  = c("Jan","Fev","Mar","Abr","Mai","Jun","Jul","Ago","Set","Out","Nov","Dez")
  }else if(lang %in% c("fr")){
     m3l  = c("Jan","Fev","Mar","Avr","Mai","Jun","Jul","Aou","Sep","Oct","Nov","Dec")
  }#end if

  if (! cap1) m3l = tolower(m3l)
  monout = m3l[mon]
  return(monout)
} #end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      Function that converts a chron object to numeric years.                             #
#------------------------------------------------------------------------------------------#
numyears <<- function(when){
   yrs    = years(when)
   lyrs   = levels(yrs)
   yrout  = as.numeric(lyrs[match(yrs,lyrs)])
   return(yrout)
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      Function that converts a chron object to numeric months.                            #
#------------------------------------------------------------------------------------------#
nummonths <<- function(when){
   mos    = months(when)
   lmos   = levels(mos)
   moout  = match(mos,lmos)
   return(moout)
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      Function that converts a chron object to numeric fortnight count.                   #
# Fortnight index goes from 1 to 26, where 1 includes January 1-14, 2 is January 15-28,    #
# and so on.  Some periods get 15 days just to make sure that nothing is lost.             #
#------------------------------------------------------------------------------------------#
numfortnight <<- function(when){
   doy   = floor(dayofyear(when))
   dmax  = 365 + is.leap(when)
   fnfac = dmax / 26
   fnout = ceiling(doy/fnfac)
   return(fnout)
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      Function that converts a chron object to numeric days.                              #
#------------------------------------------------------------------------------------------#
numdays <<- function(when){
   dys    = days(when)
   ldys   = levels(dys)
   dyout  = match(dys,ldys)
   return(dyout)
} #end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      Function that returns the dates as characters.                                      #
#------------------------------------------------------------------------------------------#
chardates <<- function(when){
   mymonth = substring(100   + nummonths(when),2,3)
   myday   = substring(100   + numdays  (when),2,3)
   myyear  = substring(10000 + numyears (when),2,5)
   mydate  = paste(mymonth,myday,myyear,sep="/")
  return(mydate)
} #end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      Function that returns the dates as characters.                                      #
#------------------------------------------------------------------------------------------#
label.dates <<- function(when,add.hours=TRUE){
   mymonth = substring(100   + nummonths(when),2,3)
   myday   = substring(100   + numdays  (when),2,3)
   myyear  = substring(10000 + numyears (when),2,5)
   mydate  = paste(myyear,mymonth,myday,sep="-")

   if (add.hours){
      mytime  = paste(substring(100 + hours  (when),2,3)
                     ,substring(100 + minutes(when),2,3)
                     ,substring(100 + seconds(when),2,3)
                     ,sep="")
      mylabel = paste(mydate,mytime,sep="-")
   }else{
      mylabel = mydate
   }#end if

  return(mylabel)
} #end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      Function that returns the times as characters.                                      #
#------------------------------------------------------------------------------------------#
chartimes <<- function(when){
   myhour = substring(100 + hours  (when),2,3)
   myminu = substring(100 + minutes(when),2,3)
   myseco = substring(100 + seconds(when),2,3)
   mytime = paste(myhour,myminu,myseco,sep=":")
  return(mytime)
} #end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      Function that finds the fraction of the day.                                        #
#------------------------------------------------------------------------------------------#
hms2frac <<- function(when){
   thishour  = hours    (when)
   thismin   = minutes  (when)
   thissec   = seconds  (when)

   elapsed = thishour / day.hr + thismin / day.min + thissec / day.sec
   return(elapsed)
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      Function that finds the numeric version of the days.                                #
#------------------------------------------------------------------------------------------#
dayofyear <<- function(when){
   offdays   = c(0, 31,59,90,120,151,181,212,243,273,304,334,365)

   thisday   = numdays  (when)
   thismonth = nummonths(when)
   thisyear  = numyears (when)
   thisfrac  = hms2frac (when)
   
   addone    = as.integer(thismonth > 2 & is.leap(when))

   doy =  thisday + offdays[thismonth] + addone + thisfrac
   return(doy)
} #end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function appends several time-related variables for a given data frame.         #
#------------------------------------------------------------------------------------------#
alltimes <<- function(datin,lon,lat,ed21=TRUE,zeronight=FALSE,meanval=FALSE,imetavg=1
                   ,nmean=120,...){
   #------ Copy the input data frame, and call the other functions. -----------------------#
   datout = datin
   datout$year       = numyears (datout$when)
   datout$month      = nummonths(datout$when)
   datout$day        = numdays  (datout$when)
   datout$hour       = hours    (datout$when)
   datout$minu       = minutes  (datout$when)
   datout$today      = dates    (datout$when)
   datout$tomonth    = chron(paste(datout$month,1,datout$year,sep="/"))
   datout$doy        = dayofyear(datout$when)
   zenith            = ed.zen   (when=datout$when,lon=lon,lat=lat,ed21=ed21
                                ,zeronight=zeronight,meanval=meanval,imetavg=imetavg
                                ,nmean=nmean,...)
   datout$cosz       =   zenith$cosz
   datout$zen        =   zenith$zen
   datout$sunhgt     =   zenith$hgt
   datout$nighttime  =   zenith$night
   datout$daytime    =   zenith$day
   datout$twilight   = (! zenith$night) & (! zenith$day)
   datout$diel       = as.integer(! datout$nighttime) + as.integer(datout$daytime)
   datout$highsun    = zenith$cosz >= cosz.highsun
   datout$riseset    = zenith$cosz >= cosz.twilight & zenith$cosz < cosz.highsun
   datout$notdaytime = ! zenith$day
   datout$season     = sign(-lat) * cos( 2.0 * pi * (dayofyear(datout$when) - 1)
                                       / (365 + is.leap(datout$year)) )

   return(datout)
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      List of trimestral seasons.                                                         #
#------------------------------------------------------------------------------------------#
season <<- function(when,add.year=FALSE,dec.next=TRUE){


   #----- Get the year and month. ---------------------------------------------------------#
   year = numyears (when)
   mon  = nummonths(when)
   #---------------------------------------------------------------------------------------#



   #----- We don't give summer/winter, instead we make generic season names. --------------#
   sidx = season.index
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Assign the season names depending on the month and year.                          #
   #---------------------------------------------------------------------------------------#
   if (add.year){
      #----- Add year before the season. --------------------------------------------------#
      seasout      = paste(year,sprintf("%2.2i",sidx[mon]),sep="")
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     December, January, and February have two years, decide whether to bind         #
      # December with next year, of January and February with previous year.               #
      #------------------------------------------------------------------------------------#
      if (dec.next){
         mp1          = mon == 12
         seasout[mp1] = paste( sprintf("%4.4i",year[mp1]+1)
                             , sprintf("%2.2i",sidx[mon[mp1]])
                             , sep = "")
      }else{
         mm1          = mon %in%  c(1,2)
         seasout[mp1] = paste( sprintf("%4.4i",year[mm1]-1)
                             , sprintf("%2.2i",sidx[mon[mm1]])
                             , sep = "")
      }#end if
      #------------------------------------------------------------------------------------#
   }else{
      #----- No year to be tagged. --------------------------------------------------------#
      seasout = sidx[mon]
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Return variable. ----------------------------------------------------------------#
   return(seasout)
   #---------------------------------------------------------------------------------------#
}#end for
#==========================================================================================#
#==========================================================================================#







#==========================================================================================#
#==========================================================================================#
#      List of periods of the day.  This is going to split the hours into blocks, with     #
# 00 UTC being always at the last block.                                                   #
#------------------------------------------------------------------------------------------#
period.day <<- function(when,dtblock=3){


   #----- Get the year and month. ---------------------------------------------------------#
   this.hour = hours  (when)
   this.minu = minutes(when)
   #----- Eliminate minutes, and always make it go to the next hour. ----------------------#
   this.hour = (this.hour + ceiling(this.minu/60)) %% 24
   #---------------------------------------------------------------------------------------#



   #----- Split the day into blocks. ------------------------------------------------------#
   shift.hour = ( this.hour - 1 ) %% 24
   break.hour = seq(from=-0.5,to=23.5,by=dtblock)
   cut.hour   = cut(x=shift.hour,breaks=break.hour)
   idx.period = match(cut.hour,levels(cut.hour))
   #---------------------------------------------------------------------------------------#



   #----- Return variable. ----------------------------------------------------------------#
   return(idx.period)
   #---------------------------------------------------------------------------------------#
}#end for
#==========================================================================================#
#==========================================================================================#







#==========================================================================================#
#==========================================================================================#
#      Create a time stamp that is a chron object from fortnight index and year.  The      #
# default is the middle point of the fortnight period.                                     #
#------------------------------------------------------------------------------------------#
fnyear.2.chron <<- function(fortnight,year,loc="centre"){

   loc    = tolower(substring(loc,1,1))


   #---------------------------------------------------------------------------------------#
   #     If fortnight is actually a chron object, find the fortnight index.                #
   #---------------------------------------------------------------------------------------#
   if (missing(year)){
      wit = is(fortnight)
      if ("dates" %in% wit || "chron" %in% wit){
         when      = fortnight
         fortnight = numfortnight(when)
         year      = numyears (when)
      }else{
         stop("  No year given and fortnight is not time")
      }#end if
   }#end if
   #---------------------------------------------------------------------------------------#


   dmax  = ifelse(is.leap(year),366,365)
   fnfac = dmax / 26
   zero  = chron(paste(12,31,year-1,sep="/"))


   if ( loc %in% c("l","b","s")){
      when.out = chron(floor(zero+(fortnight-1)*fnfac + 1))
   }else if ( loc %in% c("c","m")){
      when.out = chron(round(zero+(fortnight-0.5)*fnfac))
   }else{
      when.out = chron(floor(zero+fortnight*fnfac))
   }#end if

   return(when.out)
}#end fnyear.2.chron
#==========================================================================================#
#==========================================================================================#




#==========================================================================================#
#==========================================================================================#
#      List with the names of months and seasons.                                          #
#------------------------------------------------------------------------------------------#
mlist          <<- month.name
mlist3         <<- month.abb

if (! "mseason.1st" %in% ls()) mseason.1st = 12


if (mseason.1st %in% c(1,4,7,10)){
   season.list  <<- c("JFM","AMJ","JAS","OND","ALL")
   season.full  <<- c("January-March","April-June","July-September","October-December"
                     ,"All months")
   season.index <<- c(1,1,1,2,2,2,3,3,3,4,4,4)
   mon.add1     <<- 13
}else if (mseason.1st %in% c(2,5,8,11)){
   season.list  <<- c("NDJ","FMA","MJJ","ASO","ALL")
   season.full  <<- c("November-January","February-April","May-July","August-October"
                     ,"All months")
   season.index <<- c(1,2,2,2,3,3,3,4,4,4,1,1)
   mon.add1     <<- 11
}else if (mseason.1st %in% c(3,6,9,12)){
   season.list  <<- c("DJF","MAM","JJA","SON","ALL")
   season.full  <<- c("December-February","March-May","June-August","September-November"
                       ,"All months")
   mon.add1     <<- 12
   season.index <<- c(1,1,2,2,2,3,3,3,4,4,4,1)
}else{
   stop (paste("Invalid mseason.1st (",mseason.1st,")",sep=""))
}#end if

nseasons       <<- length(season.list)
season.cols    <<- c("#3B24B3","#2996CC","#D9BE36","#990F0F",all.colour)
season.label   <<- paste(sprintf("%2.2i",sequence(nseasons)),season.list,sep="-")
#==========================================================================================#
#==========================================================================================#
