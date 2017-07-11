#==========================================================================================#
#==========================================================================================#
#     This function creates a pretty time scale for annual cycle using fortnightly means.  #
#------------------------------------------------------------------------------------------#
pretty.ftnight <<- function(when,add.year=FALSE,...){
   #----- Find months and years. ----------------------------------------------------------#
   nwhen    = length(when)
   mons     = nummonths(when)
   yrs      = numyears (when)

   mon.app  = mons[nwhen]%%12 + 1
   yr.app   = yrs [nwhen] + (mon.app == 1)
   mons     = c(mons,mon.app)
   yrs      = c(yrs, yr.app)
   vlevels  = chron(paste(mons,1,yrs,sep="/"))
   vlevels  = unique(vlevels)
   mons     = nummonths(vlevels)
   yrs      = numyears(vlevels)


   add.year = add.year & ( mons == 1 | vlevels == vlevels[1] )
   if (any(add.year)){
      vlabels = ifelse(add.year,paste(month.abb[mons],yrs,sep="\n"),month.abb[mons])
      padj    = ifelse(add.year,0.5,0.0)
   }else{
      vlabels = substring(month.abb[mons],1,1)
      padj    = rep(0,times=length(mons))
   }#end if
   vresult=list(levels=vlevels,labels=vlabels,n=length(vlevels),scale="months",padj=padj)
   return(vresult)
}#end function
#==========================================================================================#
#==========================================================================================#

