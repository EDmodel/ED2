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
check.steady     = thischecksteady
metcyca          = thismetcyc1
metcycz          = thismetcycf
nyear.min        = thisnyearmin
ststcrit         = thisststcrit
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


#----- Output path. -----------------------------------------------------------------------#
output           = file.path(main,polyg)
#------------------------------------------------------------------------------------------#


#----- Load some useful scripts and packages. ---------------------------------------------#
isok = require(hdf5)
#------------------------------------------------------------------------------------------#



#----- Find the standard output files. ----------------------------------------------------#
endrun     = file.path(main,polyg,"serial_lsf.out" )
hasrun     = file.path(main,polyg,"serial_out.out" )
haserr     = file.path(main,polyg,"serial_out.err" )
hascrashed = file.path(main,polyg,"crashed_out.out")
hassigsegv = file.path(main,polyg,"sigsegv_out.out")
hasbad.met = file.path(main,polyg,"bad_met_out.out")
hasmetmiss = file.path(main,polyg,"metmiss_out.out")
hasstopped = file.path(main,polyg,"stopped_out.out")
#------------------------------------------------------------------------------------------#


#----- Check whether it doesn't have segmentation violation. ------------------------------#
if ( file.exists(haserr)){
   errout  = scan(haserr,what="raw")
   errout  = paste(errout,collapse=" ")
   sigsegv = ( length(grep("sigsegv"           ,errout,ignore.case=TRUE)) > 0
             | length(grep("segmentation fault",errout,ignore.case=TRUE)) > 0
             )#end sigsegv
}else{
   sigsegv = FALSE
}#end if
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Here we check whether the simulation is actually running, or if it has finished      #
# or crashed.  the existence of serial_lsf.out does not tell the entire story, because     #
# the polygon may be running under the "unrestricted_parallel" queue.                      #
#------------------------------------------------------------------------------------------#
if (file.exists(hasrun) && ! sigsegv){
   simout   = scan(hasrun,what="raw")
   simout   = paste(simout,collapse=" ")
   metmiss  = ( length(grep("Cannot open met driver input file",simout)) > 0
              | length(grep("Specify ED_MET_DRIVER_DB properly",simout)) > 0 )
   crashed  = length(grep("IFLAG1 problem."                    ,simout)) > 0
   bad.met  = length(grep("Meteorological forcing has issues"  ,simout)) > 0
   stopped  = length(grep("FATAL ERROR"                        ,simout)) > 0
   finished = length(grep("ED execution ends"                  ,simout)) > 0
   running  = ! (metmiss || crashed || stopped || finished)
}else if (file.exists(hasrun) && sigsegv){
   metmiss  = FALSE
   crashed  = FALSE
   bad.met  = FALSE
   stopped  = FALSE
   finished = FALSE
   running  = FALSE
}else{
   metmiss  = file.exists(hasmetmiss)
   bad.met  = file.exists(hasbad.met)
   crashed  = file.exists(hascrashed)
   stopped  = file.exists(hasstopped)
   sigsegv  = file.exists(hassigsegv)
   finished = FALSE
   running  = FALSE
}#end if
#------------------------------------------------------------------------------------------#




#----- Retrieve the last history file. ----------------------------------------------------#
histodir    = paste(histomain,polyg,"histo",NULL,sep="/")
histolist   = dir(histodir)
bye         =-grep("-Z-",histolist)
if (length(bye) > 0) histolist   = histolist[bye]
nhisto      = length(histolist)
hasoutput   = nhisto > 0
#------------------------------------------------------------------------------------------#





#------------------------------------------------------------------------------------------#
#      Check whether there is a last history file or not...                                #
#------------------------------------------------------------------------------------------#
if (hasoutput){
   latesthisto = histolist[nhisto]
   pathhisto   = file.path(histodir,histolist[nhisto])
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



   #---------------------------------------------------------------------------------------#
   #    In order to qualify to be killed, a polygon must be bare and it can't be run-      #
   # ning on unrestricted_parallel, because the jobs have other jobs in it.                #
   #---------------------------------------------------------------------------------------#
   desert   = mydata$NCOHORTS.GLOBAL == 0
   if (desert){
      agb   = 0.
      bsa   = 0.
      lai   = 0.
   }else{
      areapa = mydata$AREA*rep(mydata$AREA.SI,times=mydata$SIPA.N)
      areapa = areapa / sum(areapa)
      area   = rep(areapa,times=mydata$PACO.N)
      agb    = sum(mydata$NPLANT * mydata$AGB.CO * area)
      bsa    = sum(mydata$NPLANT * mydata$BA.CO  * area)
      lai    = sum(                mydata$LAI.CO * area)
   }#end if
   scb = mydata$FAST.SOIL.C.PY+mydata$STRUCT.SOIL.C.PY+mydata$SLOW.SOIL.C.PY
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Check whether we have reached steady state (if the user wants it, of course...    #
   #---------------------------------------------------------------------------------------#
   if (check.steady && yyyy > yeara + nyear.min){
      #------------------------------------------------------------------------------------#
      #    If all plants went extinct, this has become a desert.  No need to check for     #
      # steady state.                                                                      #
      #------------------------------------------------------------------------------------#
      if (desert){
         extinct = TRUE
         ststate = FALSE
      }else{
         #----- Not a desert, check for steady state. -------------------------------------#
         extinct = FALSE

         #---- Get the meteorological cycle. ----------------------------------------------#
         metcycle   = metcycz - metcyca
         years.test = seq(from=yeara,to=yeara+nhisto-1,by=metcycle)
         ntest      = length(years.test)
         npft       = dim(mydata$AGB.PY)[3]
         #---------------------------------------------------------------------------------#


         #----- Get AGB and soil carbon time series. --------------------------------------#
         agb.pft    = matrix(NA,nrow=ntest,ncol=npft)
         scb.scp    = matrix(NA,nrow=ntest,ncol=   3)
         for (y in 1:ntest){
            thishisto    = grep(paste("-",years.test[y],"-",sep=""),histolist)
            pathhisto    = paste(histodir,histolist[thishisto],sep="/")
            olddata      = hdf5load(file=pathhisto,load=FALSE,verbosity=0,tidy=TRUE)
            agb.pft[y,]  = colSums(olddata$AGB.PY[1,,])
            scb.scp[y,1] = olddata$FAST.SOIL.C.PY
            scb.scp[y,2] = olddata$STRUCT.SOIL.C.PY
            scb.scp[y,3] = olddata$SLOW.SOIL.C.PY
         }#end for
         #---------------------------------------------------------------------------------#

         #---------------------------------------------------------------------------------#
         #     Extinct PFTs by definition have reached the equilibrium.  Likewise, if any  #
         # soil carbon is zero, we cannot check so assume equilibrium.  From this point,   #
         # we assume we are at steady state.  If any PFT or soil carbon pool fail, we stop #
         # checking.                                                                       #
         #---------------------------------------------------------------------------------#
         pft.check = which(agb.pft[ntest,] > 0)
         scp.check = which(scp.sft[ntest,] > 0)
         steady    = TRUE
         #----- Check AGB for each PFT. ---------------------------------------------------#
         for (p in pft.check){
            #----- Skip check of no longer steady state. ----------------------------------#
            if (steady){
               #----- Find the AGB rate of change. ----------------------------------------#
               rate    = c(NA,diff(agb.pft[,p])) / max(agb.pft[,p])
               data.in = data.frame(rate,years.test)
               #---------------------------------------------------------------------------#


               #----- Smooth the rate of change, so we remove the noise. ------------------#
               fit.in  = loess(formula = rate~years.test,data=data.in)
               pred.in = predict(fit.in,data=data.in)
               #---------------------------------------------------------------------------#


               #----- Check whether the last point has low enough rate. -------------------#
               npred   = length(pred.in)
               steady  = ( abs(pred.in[npred])   <= ststcrit   &&
                           abs(pred.in[npred-1]) <= ststcrit     )
               #---------------------------------------------------------------------------#
            }#end if (steady)
            #------------------------------------------------------------------------------#
         }#end for (p in pft.check)
         #----- Check Soil carbon for each soil carbon pool. ------------------------------#
         for (p in scp.check){
            #----- Skip check of no longer steady state. ----------------------------------#
            if (steady){
               #----- Find the AGB rate of change. ----------------------------------------#
               rate    = c(NA,diff(scb.scp[,p])) / max(scb.scp[,p])
               data.in = data.frame(rate,years.test)
               #---------------------------------------------------------------------------#


               #----- Smooth the rate of change, so we remove the noise. ------------------#
               fit.in  = loess(formula = rate~years.test,data=data.in)
               pred.in = predict(fit.in,data=data.in)
               #---------------------------------------------------------------------------#


               #----- Check whether the last point has low enough rate. -------------------#
               npred   = length(pred.in)
               steady  = ( abs(pred.in[npred])   <= ststcrit   &&
                           abs(pred.in[npred-1]) <= ststcrit     )
               #---------------------------------------------------------------------------#
            }#end if (steady)
         }#end for (p in scp.check)
         #---------------------------------------------------------------------------------#
      }#end if (desert)
      #------------------------------------------------------------------------------------#
   }else{
      extinct = FALSE
      ststate = FALSE
   }#end if (check.steady && yyyy > yeara + nyear.min)
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
   scb       = sprintf("%7.3f",NA    )
   extinct   = FALSE
   ststate   = FALSE
   #---------------------------------------------------------------------------------------#
}#end if (nhisto > 0)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Determine the status (be pessimistic and check for all possible errors before        #
# convincing that the run is working).                                                     #
#------------------------------------------------------------------------------------------#
if (running && hasoutput){
   if (extinct){
      status    = paste(polyg,yyyy,mm,dd,hhhh,"EXTINCT",agb,bsa,lai,scb,sep=" ")
   }else if(ststate){
      status    = paste(polyg,yyyy,mm,dd,hhhh,"STSTATE",agb,bsa,lai,scb,sep=" ")
   }else{
      status    = paste(polyg,yyyy,mm,dd,hhhh,"HISTORY",agb,bsa,lai,scb,sep=" ")
   }#end if
}else if(running){
   status    = paste(polyg,yyyy,mm,dd,hhhh,"INITIAL",agb,bsa,lai,scb,sep=" ")
}else if (sigsegv){
   status    = paste(polyg,yyyy,mm,dd,hhhh,"SIGSEGV",agb,bsa,lai,scb,sep=" ")
}else if(crashed){
   status    = paste(polyg,yyyy,mm,dd,hhhh,"CRASHED",agb,bsa,lai,scb,sep=" ")
}else if(bad.met){
   status    = paste(polyg,yyyy,mm,dd,hhhh,"BAD_MET",agb,bsa,lai,scb,sep=" ")
}else if (metmiss){
   status    = paste(polyg,yyyy,mm,dd,hhhh,"METMISS",agb,bsa,lai,scb,sep=" ")
}else if(stopped){
   status    = paste(polyg,yyyy,mm,dd,hhhh,"STOPPED",agb,bsa,lai,scb,sep=" ")
}else if(finished){
   status    = paste(polyg,yyyy,mm,dd,hhhh,"THE_END",agb,bsa,lai,scb,sep=" ")
}else if(hasoutput){
   if (extinct){
      status    = paste(polyg,yyyy,mm,dd,hhhh,"EXTINCT",agb,bsa,lai,scb,sep=" ")
   }else if(ststate){
      status    = paste(polyg,yyyy,mm,dd,hhhh,"STSTATE",agb,bsa,lai,scb,sep=" ")
   }else{
      status    = paste(polyg,yyyy,mm,dd,hhhh,"HISTORY",agb,bsa,lai,scb,sep=" ")
   }#end if
}else{
   status    = paste(polyg,yyyy,mm,dd,hhhh,"INITIAL",agb,bsa,lai,scb,sep=" ")
}#end if
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
