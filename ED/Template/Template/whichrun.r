#----- This should be called before anything else, don't define stuff before this line. ---#
rm(list=ls())
graphics.off()
#------------------------------------------------------------------------------------------#


#----- The user-defined variable section. -------------------------------------------------#
main             = "pathhere"
histomain        = "pathhere"
srcdir           = "thisrscpath"
polyg            = "thispoly"
queue            = "thisqueue"
yeara            = thisyeara
montha           = thismontha
datea            = thisdatea
timea            = thistimea
output = paste(main,polyg,sep="/") # Current directory.
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#   No need to change anything beyond this point unless you are developing the script.     #
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#



#----- Load some useful scripts and packages. ---------------------------------------------#
isok = require(hdf5)
#------------------------------------------------------------------------------------------#




#----- Retrieve the last history file. ----------------------------------------------------#
histodir    = paste(histomain,polyg,"histo",NULL,sep="/")
histolist   = dir(histodir)
bye         =-grep("-Z-",histolist)
if (length(bye) > 0) histolist   = histolist[bye]
nhisto      = length(histolist)
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#      Check whether there is a last history file or not...                                #
#------------------------------------------------------------------------------------------#
if (nhisto > 0){
   latesthisto = histolist[nhisto]
   pathhisto   = paste(histodir,histolist[nhisto],sep="/")
   print (paste("Polygon",polyg," - Last:",latesthisto,"..."))

   #----- Open the last history and check how many cohorts exist. -------------------------#
   mydata = hdf5load(file=pathhisto,load=FALSE,verbosity=0,tidy=TRUE)

   fl=nchar(latesthisto)

   #------ Determine the time. ------------------------------------------------------------#
   yyyy = substring(latesthisto,fl-23,fl-20)
   mm   = substring(latesthisto,fl-18,fl-17)
   dd   = substring(latesthisto,fl-15,fl-14)
   hhhh = substring(latesthisto,fl-12,fl- 9)
   #---------------------------------------------------------------------------------------#



   #----- Find the standard output files. -------------------------------------------------#
   endrun  = paste(main,polyg,"serial_lsf.out",sep="/")
   hasrun  = paste(main,polyg,"serial_out.out",sep="/")
   haserr  = paste(main,polyg,"serial_out.err",sep="/")
   #---------------------------------------------------------------------------------------#


   #----- Check whether it is still running. ----------------------------------------------#
   running = file.exists(hasrun) && ! file.exists(endrun)
   #---------------------------------------------------------------------------------------#


   #----- Check whether it doesn't have segmentation violation. ---------------------------#
   if ( file.exists(haserr)){
      errout  = scan(haserr,what="raw")
      errout  = paste(errout,collapse=" ")
      sigsegv = ( length(grep("sigsegv"           ,errout,ignore.case=TRUE)) > 0
                | length(grep("segmentation fault",errout,ignore.case=TRUE)) > 0
                )#end sigsegv
   }else{
      sigsegv = FALSE
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Here we check whether the simulation is actually running, or if it has finished   #
   # or crashed.  the existence of serial_lsf.out does not tell the entire story, because  #
   # the polygon may be running under the "unrestricted_parallel" queue.                   #
   #---------------------------------------------------------------------------------------#
   running = file.exists(hasrun) && ! file.exists(endrun)
   if (file.exists(hasrun) && ! sigsegv){
      simout   = scan(hasrun,what="raw")
      simout   = paste(simout,collapse=" ")
      metmiss  = ( length(grep("Cannot open met driver input file",simout)) > 0
                 | length(grep("Specify ED_MET_DRIVER_DB properly",simout)) > 0 )
      crashed  = length(grep("IFLAG1 problem."                    ,simout)) > 0
      stopped  = length(grep("FATAL ERROR"                        ,simout)) > 0
      finished = length(grep("ED execution ends"                  ,simout)) > 0
   }else{
      metmiss  = FALSE
      crashed  = FALSE
      stopped  = FALSE
      finished = FALSE
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Sum the AGB and basal area for each PFT.                                           #
   #---------------------------------------------------------------------------------------#
   agb.pft  = colSums(mydata$AGB.PY       [1,,],na.rm=TRUE)
   bsa.pft  = colSums(mydata$BASAL.AREA.PY[1,,],na.rm=TRUE)
   lai.pft  = colSums(mydata$LAI.PY       [1,,],na.rm=TRUE)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    In order to qualify to be killed, a polygon must be bare and it can't be run-      #
   # ning on unrestricted_parallel, because the jobs have other jobs in it.                #
   #---------------------------------------------------------------------------------------#
   desert   = mydata$NCOHORTS.GLOBAL == 0
   agb      = sprintf("%7.3f",sum(agb.pft))
   bsa      = sprintf("%7.3f",sum(bsa.pft))
   lai      = sprintf("%7.3f",sum(lai.pft))
   #---------------------------------------------------------------------------------------#


   #----- Find out whether this job can be killed. ----------------------------------------#
   killable = queue != "unrestricted_parallel" 
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Determine the status (be pessimistic and check for all possible errors before     #
   # convincing that the run is working).                                                  #
   #---------------------------------------------------------------------------------------#
   if (sigsegv){
      status    = paste(polyg,yyyy,mm,dd,hhhh,"SIGSEGV",agb,bsa,lai,sep=" ")
   }else if (metmiss){
      status    = paste(polyg,yyyy,mm,dd,hhhh,"METMISS",agb,bsa,lai,sep=" ")
   }else if(crashed){
      status    = paste(polyg,yyyy,mm,dd,hhhh,"CRASHED",agb,bsa,lai,sep=" ")
   }else if(stopped){
      status    = paste(polyg,yyyy,mm,dd,hhhh,"STOPPED",agb,bsa,lai,sep=" ")
   }else if(finished){
      status    = paste(polyg,yyyy,mm,dd,hhhh,"THE_END",agb,bsa,lai,sep=" ")
   }else{
      status    = paste(polyg,yyyy,mm,dd,hhhh,"HISTORY",agb,bsa,lai,sep=" ")
   }#end if
   #---------------------------------------------------------------------------------------#

}else{
   #---------------------------------------------------------------------------------------#
   #     No history file, the run hasn't started yet...                                    #
   #---------------------------------------------------------------------------------------#
   yyyy      = sprintf("%4.4i",yeara )
   mm        = sprintf("%2.2i",montha)
   dd        = sprintf("%2.2i",datea )
   hhhh      = sprintf("%4.4i",timea )
   agb       = sprintf("%7.3f",NA    )
   bsa       = sprintf("%7.3f",NA    )
   lai       = sprintf("%7.3f",NA    )
   status    = paste(polyg,yyyy,mm,dd,hhhh,"INITIAL",agb,bsa,lai,sep=" ")
   #---------------------------------------------------------------------------------------#
}#end if (nhisto > 0)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Dump the output to a file.                                                           #
#------------------------------------------------------------------------------------------#
statusout    = paste(output,"statusrun.txt",sep="/")
dum          = write(x=status,file=statusout,append=FALSE)
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#      Quit R.                                                                             #
#------------------------------------------------------------------------------------------#
q("no")
#------------------------------------------------------------------------------------------#
