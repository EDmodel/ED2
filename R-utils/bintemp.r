#------------------------------------------------------------------------------------------#
# Function bintemp                                                                         #
# Developed by Marcos Longo - EPS/Harvard University                                       #
#                                                                                          #
#   This function writes a vector containing the binary files to be read in case of temp-  #
# late ctl file.                                                                           #
#------------------------------------------------------------------------------------------#

bintemp = function(binary,gtime){
  mnames = c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
#----- Adding a tiny perturbation to avoid rounding to the wrong time ---------------------#
  if (length(gtime) > 1){
    dtime = rep(mean(diff(gtime))*eps(),times=length(gtime))
  }else{
    dtime = 1/86400 # 1 second. Since GrADS finest resolution is 1 minute...
  } #end if (length(gtime) > 1)

#----- Extracting all the information that I might need -----------------------------------#
  gdtime   = gtime+dtime
  levyears = as.integer(levels(years(gdtime)))
  y4       = levyears[as.integer(years(gdtime))]
  y2       = substring(100 + (as.numeric(y4) %% 100),2,3)
  mc       = tolower(as.character(months(gtime+dtime)))
  m1       = match(mc,table=mnames)
  m2       = substring(100+as.numeric(m1),2,3)
  d1       = days(gtime+dtime)
  d2       = substring(100+as.numeric(d1),2,3)
  h1       = hours(gtime+dtime)
  h2       = substring(100+as.numeric(h1),2,3)
  n1       = minutes(gtime+dtime)
  n2       = substring(100+as.numeric(n1),2,3)
  bintemp  = rep(binary,times=length(y4))

#----- Looping until I remove all "%" from the binary name --------------------------------#
  timestr = gregexpr("%",bintemp[1])
  pos = timestr[[1]][1]
  count = 0
  while (pos > 0){
    count = count+1
    now = tolower(substr(bintemp[1],pos,pos+2))
    bef = pos-1
    aft = pos+3
#---- Checking which part of the template I'm dealing with --------------------------------#
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
#----- Looking for more "%"... ------------------------------------------------------------#
    timestr = gregexpr("%",bintemp[1])
    pos     = timestr[[1]][1]   
  } #end while (pos > 0)
  return(bintemp)
} #end function bintemp
#------------------------------------------------------------------------------------------#
