#----- This should be called before anything else, don't define stuff before this line. ---#
rm(list=ls())

#----- workdir is the user-defined variable section. -----------------------------------------#
main             = "pathhere"
histomain        = "paththere"

polyg  = "thispoly"
queue  = "thisqueue"
output = paste(main,polyg,sep="/") # Current directory.

library(hdf5)

#----- Retrieve the last history file. ----------------------------------------------------#
histodir    = paste(histomain,polyg,"histo",NULL,sep="/")
histolist   = dir(histodir)
bye         =-grep("-Z-",histolist)
if (length(bye) > 0) histolist   = histolist[bye]
nhisto      = length(histolist)

#----- Check whether there is a last history file or not... -------------------------------#
if (nhisto > 0){
   latesthisto = histolist[nhisto]
   pathhisto   = paste(histodir,histolist[nhisto],sep="/")
   print (paste("Polygon",polyg," - Last:",latesthisto,"..."))

   #----- Open the last history and check how many cohorts exist. -------------------------#
   mydata = hdf5load(file=pathhisto,load=FALSE,verbosity=0,tidy=TRUE)

   fl=nchar(latesthisto)

   yyyy = substring(latesthisto,fl-23,fl-20)
   mm   = substring(latesthisto,fl-18,fl-17)
   dd   = substring(latesthisto,fl-15,fl-14)
   hhhh = substring(latesthisto,fl-12,fl-9)

   #----- Check whether the job is still running. -----------------------------------------#
   endrun  = paste(main,polyg,"serial_lsf.out",sep="/")
   hasrun  = paste(main,polyg,"serial_out.out",sep="/")

   running = file.exists(hasrun) && ! file.exists(endrun)

   #---------------------------------------------------------------------------------------#
   #     Here we check whether the simulation is actually running, or if it has finished   #
   # or crashed.  the existence of serial_lsf.out does not tell the entire story, because  #
   # the polygon may be running under the "unrestricted_parallel" queue.                   #
   #---------------------------------------------------------------------------------------#
   if (file.exists(hasrun)){
      simout   = scan(hasrun,what="raw")
      simout   = paste(simout,collapse=" ")
      crashed  = length(grep("FATAL ERROR",simout))       > 0
      finished = length(grep("ED execution ends",simout)) > 0
   }else{
      crashed  = FALSE
      finished = FALSE
   }#end if


   #---------------------------------------------------------------------------------------#
   #    In order to qualify to be killed, a polygon must be bare and it can't be run-      #
   # ning on unrestricted_parallel, because the jobs have other jobs in it.                #
   #---------------------------------------------------------------------------------------#
   desert   = mydata$NCOHORTS.GLOBAL == 0
   agb      = sum(mydata$AGB)
   lai      = mydata$LAI

   killable = queue != "unrestricted_parallel" 

   if (desert){
      status    = paste(polyg,yyyy,mm,dd,hhhh,"EXTINCT",agb,lai,sep=" ")
   }else if(crashed){
      status    = paste(polyg,yyyy,mm,dd,hhhh,"CRASHED",agb,lai,sep=" ")
   }else if(finished){
      status    = paste(polyg,yyyy,mm,dd,hhhh,"THE_END",agb,lai,sep=" ")
   }else{
      status    = paste(polyg,yyyy,mm,dd,hhhh,"HISTORY",agb,lai,sep=" ")
   }#end if

}else{
   status    = paste(polyg,"1500","01","01","0000","INITIAL","NA","NA",sep=" ")
}#end if (nhisto > 0)

statusout    = paste(output,"statusrun.txt",sep="/")
dum          = write(x=status,file=statusout,append=FALSE)

q("no")
