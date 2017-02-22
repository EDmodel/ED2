#------------------------------------------------------------------------------------------#
# Function bintemp                                                                         #
# Developed by Marcos Longo - EPS/Harvard University                                       #
#                                                                                          #
#   This function writes a vector containing the binary files to be read in case of temp-  #
# late ctl file.                                                                           #
#------------------------------------------------------------------------------------------#
bintemp <<- function(binary,gtime){
   mnames = tolower(month.abb)
   seps   = 2^-24

   #----- Add a tiny perturbation to avoid rounding to the wrong time ---------------------#
   if (length(gtime) > 1){
     dtime = rep(mean(diff(gtime))*seps,times=length(gtime))
   }else{
     dtime = 1/86400 # 1 second. Since GrADS finest resolution is 1 minute...
   } #end if (length(gtime) > 1)
   #---------------------------------------------------------------------------------------#


   #----- Extract all the information that may be used. -----------------------------------#
   gdtime   = gtime+dtime
   levyears = as.integer(levels(years(gdtime)))
   y4       = levyears[as.integer(years(gdtime))]
   y2       = sprintf("%2.2i",as.numeric(y4) %% 100)
   mc       = tolower(as.character(months(gtime+dtime)))
   m1       = match(mc,table=mnames)
   m2       = sprintf("%2.2i",as.numeric(m1))
   d1       = days(gtime+dtime)
   d2       = sprintf("%2.2i",as.numeric(d1))
   h1       = hours(gtime+dtime)
   h2       = sprintf("%2.2i",as.numeric(h1))
   n1       = minutes(gtime+dtime)
   n2       = sprintf("%2.2i",as.numeric(n1))
   bintemp  = rep(binary,times=length(y4))
   #---------------------------------------------------------------------------------------#



   #----- Loop until all "%" are removed from the binary name. ----------------------------#
   timestr = gregexpr("%",bintemp[1])
   pos     = timestr[[1]][1]
   count   = 0
   while (pos > 0){
      count = count+1
      now   = tolower(substr(bintemp[1],pos,pos+2))
      bef   = pos-1
      aft   = pos+3
      #---- Check which part of the template I'm dealing with -----------------------------#
      if(now=="%y4") bintemp = paste(substr(bintemp,1,bef),y4,substring(bintemp,aft),sep="")
      if(now=="%y2") bintemp = paste(substr(bintemp,1,bef),y2,substring(bintemp,aft),sep="")
      if(now=="%mc") bintemp = paste(substr(bintemp,1,bef),mc,substring(bintemp,aft),sep="")
      if(now=="%m1") bintemp = paste(substr(bintemp,1,bef),m1,substring(bintemp,aft),sep="")
      if(now=="%m2") bintemp = paste(substr(bintemp,1,bef),m2,substring(bintemp,aft),sep="")
      if(now=="%d1") bintemp = paste(substr(bintemp,1,bef),d1,substring(bintemp,aft),sep="")
      if(now=="%d2") bintemp = paste(substr(bintemp,1,bef),d2,substring(bintemp,aft),sep="")
      if(now=="%h1") bintemp = paste(substr(bintemp,1,bef),h1,substring(bintemp,aft),sep="")
      if(now=="%h2") bintemp = paste(substr(bintemp,1,bef),h2,substring(bintemp,aft),sep="")
      if(now=="%n2") bintemp = paste(substr(bintemp,1,bef),n2,substring(bintemp,aft),sep="")
      #------------------------------------------------------------------------------------#

      #----- Look for more "%"... ---------------------------------------------------------#
      timestr = gregexpr("%",bintemp[1])
      pos     = timestr[[1]][1]   
      #------------------------------------------------------------------------------------#
   }#end while (pos > 0)
   return(bintemp)
}#end function bintemp
#------------------------------------------------------------------------------------------#
