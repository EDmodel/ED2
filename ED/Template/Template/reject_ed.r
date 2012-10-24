#----- Here is the user-defined variable section. -----------------------------------------#
here           = "thispath" # Current directory.
srcdir         = "/n/moorcroft_data/mlongo/util/Rsc"      # Source  directory.
outroot        = "thisoutroot"
myplaces       = c("thispoly")
outform        = "thisoutform"          # Formats for output file.  Supported formats are:
                                #   - "X11" - for printing on screen
                                #   - "eps" - for postscript printing
                                #   - "png" - for PNG printing
                                #   - "pdf" - for PDF printing

ntopoffe       =  6             # Number of top "offenders" that go to the plot.
coloffe        = c("firebrick","goldenrod","chartreuse"
                  ,"forestgreen","steelblue","midnightblue")

depth          = 96             # PNG resolution, in pixels per inch
paper          = "letter"       # Paper size, to define the plot shape
ptsz           = 14             # Font size.
lwidth         = 2.5            # Line width
plotgrid       = TRUE           # Should I plot the grid in the background? 
scalleg        = 0.25           # Increase in y scale to fit the legend.
ncolshov       = 200            # Target number of colours for Hovmoller diagrams.
hovgrid        = TRUE           # Should I include a grid on the Hovmoller plots?
offegrid       = TRUE           # Should I include a grid on the top offender plots?

#----- Comparative plot of monthly means of errors. ---------------------------------------#
nbox = 14
bplot01 = list(vnam="can.theta"    ,desc="Canopy air pot. temperature"   ,unit="",plt=T)
bplot02 = list(vnam="can.shv"      ,desc="Canopy air specific humidity"  ,unit="",plt=T)
bplot03 = list(vnam="can.co2"      ,desc="Canopy air CO2 mixing ratio"   ,unit="",plt=T)
bplot04 = list(vnam="can.prss"     ,desc="Canopy air pressure"           ,unit="",plt=T)
bplot05 = list(vnam="veg.energy"   ,desc="Vegetation energy/temperature" ,unit="",plt=T)
bplot06 = list(vnam="veg.water"    ,desc="Vegetation temperature"        ,unit="",plt=T)
bplot07 = list(vnam="soil.ener.16" ,desc="Top soil temperature"          ,unit="",plt=T)
bplot08 = list(vnam="soil.water.16",desc="Top soil moisture"             ,unit="",plt=T)
bplot09 = list(vnam="sfcw.ener.01" ,desc="Surface water energy"          ,unit="",plt=T)
bplot10 = list(vnam="sfcw.mass.01" ,desc="Surface water mass"           ,unit="",plt=T)
bplot11 = list(vnam="co2b.storage" ,desc="Carbon storage"               ,unit="",plt=T)
bplot12 = list(vnam="co2b.loss2atm",desc="Carbon loss to atmosphere"    ,unit="",plt=T)
bplot13 = list(vnam="eb.loss2atm"  ,desc="Energy loss to atmosphere"    ,unit="",plt=T)
bplot14 = list(vnam="watb.loss2atm",desc="Water loss to atmosphere"     ,unit="",plt=T)

#----- Similar to Hovmoller diagrams. -----------------------------------------------------#
nhov = 14
hovdi01 = list(vnam="can.theta"     ,desc="Canopy air pot. temperature"     ,unit=""  
                                    ,csch="muitas"                          ,plt=T)
hovdi02 = list(vnam="can.shv"       ,desc="Canopy air specific humidity"    ,unit=""  
                                    ,csch="muitas"                          ,plt=T)
hovdi03 = list(vnam="can.co2"       ,desc="Canopy air CO2 mixing ratio"     ,unit=""  
                                    ,csch="muitas"                          ,plt=T)
hovdi04 = list(vnam="can.prss"      ,desc="Canopy air pressure"             ,unit="" 
                                    ,csch="muitas"                          ,plt=T)
hovdi05 = list(vnam="veg.energy"    ,desc="Vegetation energy/temperature"   ,unit=""  
                                    ,csch="muitas"                          ,plt=T)
hovdi06 = list(vnam="veg.water"     ,desc="Vegetation temperature"          ,unit=""    
                                    ,csch="muitas"                          ,plt=T)
hovdi07 = list(vnam="soil.ener.16"  ,desc="Top soil temperature"            ,unit=""    
                                    ,csch="muitas"                          ,plt=T)
hovdi08 = list(vnam="soil.water.16" ,desc="Top soil moisture"               ,unit=""   
                                    ,csch="muitas"                          ,plt=T)
hovdi09 = list(vnam="sfcw.ener.01"  ,desc="Surface water energy"            ,unit=""
                                    ,csch="muitas"                          ,plt=T)
hovdi10 = list(vnam="sfcw.mass.01"  ,desc="Surface water mass"              ,unit=""    
                                    ,csch="muitas"                          ,plt=T)
hovdi11 = list(vnam="co2b.storage"  ,desc="'Carbon storage"                 ,unit=""    
                                    ,csch="muitas"                          ,plt=T)
hovdi12 = list(vnam="co2b.loss2atm" ,desc="'Carbon loss to atmosphere"      ,unit=""     
                                    ,csch="muitas"                          ,plt=T)
hovdi13 = list(vnam="eb.loss2atm"   ,desc="'Energy loss to atmosphere"      ,unit=""     
                                    ,csch="muitas"                          ,plt=T)
hovdi14 = list(vnam="watb.loss2atm" ,desc="'Water loss to atmosphere"       ,unit=""     
                                    ,csch="muitas"                          ,plt=T)
#------------------------------------------------------------------------------------------#


#------- Name of all variables for the top offenders. -------------------------------------#
nvad=61
vades01 = list(vnam="year"          ,desc="Year"                        )
vades02 = list(vnam="mon"           ,desc="Month"                       )
vades03 = list(vnam="day"           ,desc="Day"                         )
vades04 = list(vnam="can.theiv"     ,desc="Canopy Theta_Eiv"            )
vades05 = list(vnam="can.theta"     ,desc="Canopy Pot. Temp."           )
vades06 = list(vnam="can.shv"       ,desc="Canopy Specific humidity"    )
vades07 = list(vnam="can.temp"      ,desc="Canopy temperature"          )
vades08 = list(vnam="can.prss"      ,desc="Canopy pressure"             )
vades09 = list(vnam="can.co2"       ,desc="Canopy CO2"                  )
vades10 = list(vnam="veg.water"     ,desc="Leaf water"                  )
vades11 = list(vnam="veg.energy"    ,desc="Leaf energy/temperature"     )
vades12 = list(vnam="virt.heat"     ,desc="Virtual energy/temperature"  )
vades13 = list(vnam="virt.water"    ,desc="Virtual water"               )
vades14 = list(vnam="co2b.storage"  ,desc="CO2 storage"                 )
vades15 = list(vnam="co2b.loss2atm" ,desc="CO2 loss to atmosphere"      )
vades16 = list(vnam="eb.loss2atm"   ,desc="Energy loss to atmosphere"   )
vades17 = list(vnam="watb.loss2atm" ,desc="Water loss to atmosphere"    )
vades18 = list(vnam="enb.loss2dra"  ,desc="Energy drainage"             )
vades19 = list(vnam="watb.loss2dra" ,desc="Water drainage"              )
vades20 = list(vnam="enb.storage"   ,desc="Energy storage"              )
vades21 = list(vnam="watb.storage"  ,desc="Water storage"               )
vades22 = list(vnam="soil.water.01" ,desc="Soil water  (z= 1)"          )
vades23 = list(vnam="soil.water.02" ,desc="Soil water  (z= 2)"          )
vades24 = list(vnam="soil.water.03" ,desc="Soil water  (z= 3)"          )
vades25 = list(vnam="soil.water.04" ,desc="Soil water  (z= 4)"          )
vades26 = list(vnam="soil.water.05" ,desc="Soil water  (z= 5)"          )
vades27 = list(vnam="soil.water.06" ,desc="Soil water  (z= 6)"          )
vades28 = list(vnam="soil.water.07" ,desc="Soil water  (z= 7)"          )
vades29 = list(vnam="soil.water.08" ,desc="Soil water  (z= 8)"          )
vades30 = list(vnam="soil.water.09" ,desc="Soil water  (z= 9)"          )
vades31 = list(vnam="soil.water.10" ,desc="Soil water  (z=10)"          )
vades32 = list(vnam="soil.water.11" ,desc="Soil water  (z=11)"          )
vades33 = list(vnam="soil.water.12" ,desc="Soil water  (z=12)"          )
vades34 = list(vnam="soil.water.13" ,desc="Soil water  (z=13)"          )
vades35 = list(vnam="soil.water.14" ,desc="Soil water  (z=14)"          )
vades36 = list(vnam="soil.water.15" ,desc="Soil water  (z=15)"          )
vades37 = list(vnam="soil.water.16" ,desc="Soil water  (z=16)"          )
vades38 = list(vnam="soil.ener.01"  ,desc="Soil energy (z= 1)"          )
vades39 = list(vnam="soil.ener.02"  ,desc="Soil energy (z= 2)"          )
vades40 = list(vnam="soil.ener.03"  ,desc="Soil energy (z= 3)"          )
vades41 = list(vnam="soil.ener.04"  ,desc="Soil energy (z= 4)"          )
vades42 = list(vnam="soil.ener.05"  ,desc="Soil energy (z= 5)"          )
vades43 = list(vnam="soil.ener.06"  ,desc="Soil energy (z= 6)"          )
vades44 = list(vnam="soil.ener.07"  ,desc="Soil energy (z= 7)"          )
vades45 = list(vnam="soil.ener.08"  ,desc="Soil energy (z= 8)"          )
vades46 = list(vnam="soil.ener.09"  ,desc="Soil energy (z= 9)"          )
vades47 = list(vnam="soil.ener.10"  ,desc="Soil energy (z=10)"          )
vades48 = list(vnam="soil.ener.11"  ,desc="Soil energy (z=11)"          )
vades49 = list(vnam="soil.ener.12"  ,desc="Soil energy (z=12)"          )
vades50 = list(vnam="soil.ener.13"  ,desc="Soil energy (z=13)"          )
vades51 = list(vnam="soil.ener.14"  ,desc="Soil energy (z=14)"          )
vades52 = list(vnam="soil.ener.15"  ,desc="Soil energy (z=15)"          )
vades53 = list(vnam="soil.ener.16"  ,desc="Soil energy (z=16)"          )
vades54 = list(vnam="sfcw.ener.01"  ,desc="Sfc. water energy (z= 1)"    )
vades55 = list(vnam="sfcw.ener.02"  ,desc="Sfc. water energy (z= 2)"    )
vades56 = list(vnam="sfcw.ener.03"  ,desc="Sfc. water energy (z= 3)"    )
vades57 = list(vnam="sfcw.ener.04"  ,desc="Sfc. water energy (z= 4)"    )
vades58 = list(vnam="sfcw.mass.01"  ,desc="Sfc. water mass   (z= 1)"    )
vades59 = list(vnam="sfcw.mass.02"  ,desc="Sfc. water mass   (z= 2)"    )
vades60 = list(vnam="sfcw.mass.03"  ,desc="Sfc. water mass   (z= 3)"    )
vades61 = list(vnam="sfcw.mass.04"  ,desc="Sfc. water mass   (z= 4)"    )
#------------------------------------------------------------------------------------------#



#----- Loading some packages. -------------------------------------------------------------#
library(hdf5)
library(chron)
library(scatterplot3d)
library(lattice)
library(maps)
library(mapdata)
library(akima)
#------------------------------------------------------------------------------------------#



#----- In case there is some graphic still opened. ----------------------------------------#
graphics.off()
#------------------------------------------------------------------------------------------#



#----- Setting how many formats we must output. -------------------------------------------#
outform = tolower(outform)
nout = length(outform)
#------------------------------------------------------------------------------------------#



#----- Names of all kinds of errors to be analysed. ---------------------------------------#
errvar  = c("errmax","sanchk")
errdesc = c("Maximum Error","Sanity Check")
nerror  = length(errvar)
#------------------------------------------------------------------------------------------#



#----- Avoiding unecessary and extremely annoying beeps. ----------------------------------#
options(locatorBell=FALSE)
#------------------------------------------------------------------------------------------#



#----- Loading some files with functions. -------------------------------------------------#
source(paste(srcdir,"atlas.r",sep="/"))
source(paste(srcdir,"globdims.r",sep="/"))
source(paste(srcdir,"locations.r",sep="/"))
source(paste(srcdir,"muitas.r",sep="/"))
source(paste(srcdir,"pretty.log.r",sep="/"))
source(paste(srcdir,"pretty.time.r",sep="/"))
source(paste(srcdir,"plotsize.r",sep="/"))
source(paste(srcdir,"qapply.r",sep="/"))
source(paste(srcdir,"rconstants.r",sep="/"))
source(paste(srcdir,"sombreado.r",sep="/"))
source(paste(srcdir,"southammap.r",sep="/"))
source(paste(srcdir,"timeutils.r",sep="/"))
#------------------------------------------------------------------------------------------#



#----- Defining plot window size ----------------------------------------------------------#
size = plotsize(proje=FALSE,paper=paper)
#------------------------------------------------------------------------------------------#



#----- Box plots --------------------------------------------------------------------------#
bplot = list()
for (s in 1:nbox){
  sss  = substring(100+s,2,3)
  bpbp = paste("bplot",sss,sep="")
  if (s == 1){
     bplot   = get(bpbp)
  }else{
     bpmerge = get(bpbp)
     for (n in union(names(bplot),names(bpmerge))){
        bplot[[n]] = c(bplot[[n]],bpmerge[[n]])
     } #end for
  } #end if
} #end for
#------------------------------------------------------------------------------------------#



#----- Hovmoller diagram ------------------------------------------------------------------#
hovdi = list()
for (s in 1:nhov){
  sss  = substring(100+s,2,3)
  hdhd = paste("hovdi",sss,sep="")
  if (s == 1){
     hovdi   = get(hdhd)
  }else{
     hdmerge = get(hdhd)
     for (n in union(names(hovdi),names(hdmerge))){
        hovdi[[n]] = c(hovdi[[n]],hdmerge[[n]])
     } #end for
  } #end if
} #end for
#------------------------------------------------------------------------------------------#



#----- Variable description  --------------------------------------------------------------#
vades = list()
for (v in 1:nvad){
  vvv  = substring(100+v,2,3)
  vdvd = paste("vades",vvv,sep="")
  if (v == 1){
     vades   = get(vdvd)
  }else{
     vdmerge = get(vdvd)
     for (n in union(names(vades),names(vdmerge))){
        vades[[n]] = c(vades[[n]],vdmerge[[n]])
     } #end for
  } #end if
} #end for
#------------------------------------------------------------------------------------------#



#---- Create the main output directory in case there is none. -----------------------------#
if (! file.exists(outroot)) dir.create(outroot)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Big place loop starts here...                                                        #
#------------------------------------------------------------------------------------------#
for (place in myplaces){

   #----- Retrieve default information about this place and set up some variables. --------#
   thispoi    = locations(where=place,here=here)
   inpref     = paste(here,place,sep="/")
   outpref    = paste(outroot,place,sep="/")
   lieu       = thispoi$lieu
   suffix     = thispoi$iata
   errmaxfile = paste(inpref,"error_max_count.txt",sep="/")
   sanchkfile = paste(inpref,"sanity_check_count.txt",sep="/")
   errmaxout  = paste(outpref,"errmax",sep="/")
   sanchkout  = paste(outpref,"sanchk",sep="/")
   #---------------------------------------------------------------------------------------#


   #----- Print the banner to entretain the user. -----------------------------------------#
   print (paste("  > ",thispoi$lieu,"...",sep=""))


   #----- Check whether the output directories exist. -------------------------------------#
   if (! file.exists(outpref)) dir.create(outpref)
   if (! file.exists(errmaxout)) dir.create(errmaxout)
   if (! file.exists(sanchkout)) dir.create(sanchkout)
   #---------------------------------------------------------------------------------------#



   #----- Read the first line of the errmax file, just to grab the variable names. --------#
   vnames   = scan(file=errmaxfile,what="raw",nlines=1,quiet=TRUE)
   nvars    = length(vnames)
   for (v in 1:nvars){
      aux          = tolower(vnames[v])
      saux         = strsplit(aux,split="")[[1]]
      uscore       = which(saux == "_")
      saux[uscore] = "."
      vnames[v]    = paste(saux,collapse="")
   }#end for
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Read the error max file, this time reading all the data and skipping the first    #
   # line.                                                                                 #
   #---------------------------------------------------------------------------------------#
   aux                 = as.numeric(scan(file=errmaxfile,what="numeric",skip=1,quiet=TRUE))
   errmaxraw           = matrix(aux,ncol=nvars,byrow=TRUE)
   dimnames(errmaxraw) = list(NULL,vnames)
   yyyymm              = paste(errmaxraw[,"year"]
                              ,substring(100+errmaxraw[,"mon"],2,3),sep="-")
   #----- Monthly totals. -----------------------------------------------------------------#
   tserrmax            = qapply(X=errmaxraw,INDEX=yyyymm,DIM=1,FUN=sum)
   tserrmax            = data.frame(tserrmax)
   yyyymm              = unique(yyyymm)
   tserrmax$year       = as.numeric(substring(yyyymm,1,4))
   tserrmax$mon        = as.numeric(substring(yyyymm,6,7))
   tserrmax$day        = 15
   #----- Monthly means. ------------------------------------------------------------------#
   mon                 = tserrmax$mon
   cyerrmax            = qapply(X=tserrmax,INDEX=mon,DIM=1,FUN=mean)
   cyerrmax            = data.frame(cyerrmax)
   mon                 = unique(mon)
   cyerrmax$mon        = mon
   #----- Total, and rank the variables by the number of rejections. ----------------------#
   tterrmax            = colSums(errmaxraw)
   tterrmax[1:3]       = 0
   errmaxoffe          = order(tterrmax,decreasing=TRUE)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Read the sanity check file, this time reading all the data and skipping the first #
   # line.                                                                                 #
   #---------------------------------------------------------------------------------------#
   aux                 = as.numeric(scan(file=sanchkfile,what="numeric",skip=1,quiet=TRUE))
   sanchkraw           = matrix(aux,ncol=nvars,byrow=TRUE)
   dimnames(sanchkraw) = list(NULL,vnames)
   yyyymm              = paste(sanchkraw[,"year"]
                              ,substring(100+sanchkraw[,"mon"],2,3),sep="-")
   #----- Monthly totals. -----------------------------------------------------------------#
   tssanchk            = qapply(X=sanchkraw,INDEX=yyyymm,DIM=1,FUN=sum)
   tssanchk            = data.frame(tssanchk)
   yyyymm              = unique(yyyymm)
   tssanchk$year       = as.numeric(substring(yyyymm,1,4))
   tssanchk$mon        = as.numeric(substring(yyyymm,6,7))
   tssanchk$day        = 15
   #----- Monthly means. ------------------------------------------------------------------#
   mon                 = tssanchk$mon
   cysanchk            = qapply(X=tssanchk,INDEX=mon,DIM=1,FUN=mean)
   cysanchk            = data.frame(cysanchk)
   mon                 = unique(mon)
   cysanchk$mon        = mon
   #----- Total, and rank the variables by the number of rejections. ----------------------#
   ttsanchk            = colSums(sanchkraw)
   ttsanchk[1:3]       = 0
   sanchkoffe          = order(ttsanchk,decreasing=TRUE)
   #---------------------------------------------------------------------------------------#



   #----- Define months with three letters. -----------------------------------------------#
   errmaxmon3             = as.factor(mon2mmm(tserrmax$mon))
   sanchkmon3             = as.factor(mon2mmm(tssanchk$mon))
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #   Plot the monthly boxplots.                                                          #
   #---------------------------------------------------------------------------------------#
   for (e in 1:nerror){
      #----- Retrieve the error we are about to plot and associated variables. ------------#
      thiserr = get(paste("ts",errvar[e],sep=""))
      thiseds = errdesc[e]
      thisout = get(paste(errvar[e],"out",sep=""))
      thism3  = get(paste(errvar[e],"mon3",sep=""))

      print (paste("    #",errvar[e]," seasonal cycle..."))

      for (v in 1:nbox){

         #----- Retrieve variable information from the list. ------------------------------#
         vnam        = bplot$vnam[v]
         description = bplot$desc[v]
         unit        = bplot$unit[v]
         plotit      = bplot$plt[v]

         if (plotit){
            #------------------------------------------------------------------------------#
            #     Check if the directory exists.  If not, create it.                       #
            #------------------------------------------------------------------------------#
            outdir  =  paste(thisout,"boxplot",sep="/")
            if (! file.exists(outdir)) dir.create(outdir)
            print (paste("      +",description,"box plot..."))

            #----- Load this variable into "thisvar". -------------------------------------#
            thisvar = thiserr[[vnam]]

            for (o in 1:nout){
               fichier = paste(outdir,"/",vnam,"-",suffix,".",outform[o],sep="")
               if (outform[o] == "x11"){
                  X11(width=size$width,height=size$height,pointsize=ptsz)
               }else if(outform[o] == "png"){
                  png(filename=fichier,width=size$width*depth,height=size$height*depth
                     ,pointsize=ptsz,res=depth)
               }else if(outform[o] == "eps"){
                  postscript(file=fichier,width=size$width,height=size$height
                            ,pointsize=ptsz,paper=size$paper)
               }else if(outform[o] == "pdf"){
                  pdf(file=fichier,onefile=FALSE
                     ,width=size$width,height=size$height,pointsize=ptsz,paper=size$paper)
               }#end if
               ylimit  = range(thisvar, na.rm=TRUE)
               
               letitre = paste(thiseds,description,lieu,sep=" - ")
               plot(thism3,thisvar,main=letitre,ylim=ylimit
                   ,xlab="Time",ylab=paste("Error counts [Error/month]",sep=""))

               if (outform[o] == "x11"){
                  locator(n=1)
                  dev.off()
               }else{
                  dev.off()
               }#end if
            } #end for outform
         }#end if
      }#end for nbox
   }#end for error
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #   Plot the Hovmoller diagrams.                                                        #
   #---------------------------------------------------------------------------------------#
   for (e in 1:nerror){
      #----- Retrieve the error we are about to plot and associated variables. ------------#
      thiserr = get(paste("ts",errvar[e],sep=""))
      thiseds = errdesc[e]
      thisout = get(paste(errvar[e],"out",sep=""))
      thism3  = get(paste(errvar[e],"mon3",sep=""))

      print (paste("    #",errvar[e]," seasonal cycle..."))

      for (v in 1:nhov){

         #----- Retrieve variable information from the list. ------------------------------#
         vnam        = hovdi$vnam[v]
         description = hovdi$desc[v]
         unit        = hovdi$unit[v]
         vcscheme    = hovdi$csch[v]
         plotit      = hovdi$plt[v]

         #---------------------------------------------------------------------------------#
         #     Find the last full year.  This will be the actual last year only if the     #
         # year is complete, otherwise it will stop at the previous year.  If this         #
         # simulation has not completed two years, then this plot will be ignored.         #
         #---------------------------------------------------------------------------------#
         ntimes   = length(thiserr[,1])
         yearzz   = max(thiserr$year)
         twoyears = ntimes >= 24
         if (ntimes %% 12 == 0){
            sel      = thiserr$year  <= yearzz
         }else{
            sel      = thiserr$year  <  yearzz
         }#end if

         if (plotit && twoyears){
            #------------------------------------------------------------------------------#
            #     Check if the directory exists.  If not, create it.                       #
            #------------------------------------------------------------------------------#
            outdir  =  paste(thisout,"hovmoller",sep="/")
            if (! file.exists(outdir)) dir.create(outdir)
            print (paste("      +",description,"Hovmoller time series..."))

            #----- Load this variable into "thisvar". -------------------------------------#
            thisvar = thiserr[[vnam]]

            #----- Find the number of rows and columns, and the axes. ---------------------#
            monaxis = sort(unique(thiserr$mon[sel]))
            yraxis  = sort(unique(thiserr$year[sel]))
            nmon    = length(monaxis)
            nyear   = length(yraxis)

            #----- Save the meaningful months and years. ----------------------------------#
            monat   = 1:12
            monlab  = c("J","F","M","A","M","J","J","A","S","O","N","D")
            yrat    = pretty(yraxis)

            #----- Convert the vector data into an array. ---------------------------------#
            vararr  = array(thisvar[sel],c(nmon,nyear))

            #----- Convert the vector data into an array. ---------------------------------#
            vararr  = array(thisvar[sel],c(nmon,nyear))

            #----- Copy Decembers ans Januaries to make the edges buffered. ---------------#
            january  = vararr[1,]
            january  = c(january,january[nyear],january[nyear])

            december = vararr[12,]
            december = c(december[1],december[1],december)

            #----- Bind first and last year to the array, to make the edges buffered. -----#
            varbuff  = cbind(vararr[,1],vararr,vararr[,nyear])
            varbuff  = rbind(december,varbuff,january)

            #----- Expand the month and year axes. ----------------------------------------#
            monaxis = c(min(monaxis)-1,monaxis,max(monaxis)+1)
            yraxis  = c(min(yraxis)-1,yraxis,max(yraxis)+1)

            vrange  = range(varbuff,na.rm=TRUE)
            vlevels = pretty(x=vrange,n=ncolshov)
            vnlev   = length(vlevels)

            for (o in 1:nout){
               fichier = paste(outdir,"/",vnam,"-",suffix,".",outform[o],sep="")
               if (outform[o] == "x11"){
                  X11(width=size$width,height=size$height,pointsize=ptsz)
               }else if(outform[o] == "png"){
                  png(filename=fichier,width=size$width*depth,height=size$height*depth
                     ,pointsize=ptsz,res=depth)
               }else if(outform[o] == "eps"){
                  postscript(file=fichier,width=size$width,height=size$height
                            ,pointsize=ptsz,paper=size$paper)
               }else if(outform[o] == "pdf"){
                  pdf(file=fichier,onefile=FALSE
                     ,width=size$width,height=size$height,pointsize=ptsz,paper=size$paper)
               }#end if

               letitre = paste(description," - ",lieu,sep="")
               sombreado(x=monaxis,y=yraxis,z=varbuff,levels=vlevels,nlevels=vnlev
                        ,color.palette=get(vcscheme)
                        ,plot.title=title(main=letitre,xlab="Month",ylab="Year")
                        ,key.title=title(main="[Error/Month]",cex.main=0.8)
                        ,plot.axes={axis(side=1,at=monat,labels=monlab)
                                    axis(side=2,at=yrat)
                                    if (hovgrid){
                                       for (yl in yrat){
                                          abline(h=yl,col="lightgray",lty="dotted")
                                       } #end for yl
                                       for (ml in monat){
                                          abline(v=ml,col="lightgray",lty="dotted")
                                       } #end for ml
                                    }#end if hovgrid
                                   }#end plot.axes
                        )

               if (outform[o] == "x11"){
                  locator(n=1)
                  dev.off()
               }else{
                  dev.off()
               }#end if
            } #end for outform
         }#end if
      }#end for nbox
   }#end for error
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #   Plot the 6 most recurrent offenders.                                                #
   #---------------------------------------------------------------------------------------#
   for (e in 1:nerror){
      #----- Retrieve the error we are about to plot and associated variables. ------------#
      thiserr  = get(paste("cy",errvar[e],sep=""))
      thisrank = get(paste(errvar[e],"offe",sep=""))
      thiseds = errdesc[e]
      thisout = get(paste(errvar[e],"out",sep=""))

      print (paste("    #",errvar[e]," top offender seasonal cycle..."))

      #----- First loop, just to get the range. -------------------------------------------#
      ylimit = NULL
      lege   = NULL
      colleg = NULL
      for (v in 1:ntopoffe){

         offe = thisrank[v]

         #----- Retrieve variable information from the list. ------------------------------#
         vnam        = vades$vnam[offe]
         description = vades$desc[offe]

         #----- Retrieve the variable and count the number of offenses. -------------------#
         thisvar     = thiserr[vnam]

         #----- Find the range. -----------------------------------------------------------#
         if (any(thisvar > 0)){
            ylimit = range(c(ylimit,range(thisvar[thisvar >0])))
         }#end if

         #----- Append the description and colour for plots. ------------------------------#
         lege   = c(lege,description)
         colleg = c(colleg,coloffe[v])
      }#end for

      #----- Make the range larger to make sure that the legend will fit. -----------------#
      ylimit[2] = ylimit[2] * (ylimit[2]/ylimit[1])^scalleg

      #----- Save the meaningful months and years. ----------------------------------------#
      monat   = 1:12
      monlab  = c("J","F","M","A","M","J","J","A","S","O","N","D")
      errat   = pretty.log(ylimit)

      outdir  =  paste(thisout,"topoffe",sep="/")
      if (! file.exists(outdir)) dir.create(outdir)

      #----- Loop over formats. -----------------------------------------------------------#
      for (o in 1:nout){
         fichier = paste(outdir,"/topoffe","-",suffix,".",outform[o],sep="")
         if (outform[o] == "x11"){
            X11(width=size$width,height=size$height,pointsize=ptsz)
         }else if(outform[o] == "png"){
            png(filename=fichier,width=size$width*depth,height=size$height*depth
               ,pointsize=ptsz,res=depth)
         }else if(outform[o] == "eps"){
            postscript(file=fichier,width=size$width,height=size$height
                      ,pointsize=ptsz,paper=size$paper)
         }else if(outform[o] == "pdf"){
            pdf(file=fichier,onefile=FALSE
               ,width=size$width,height=size$height,pointsize=ptsz,paper=size$paper)
         }#end if

         letitre = paste("Mean annual cycle - ",lieu,sep="")
         plot(x=thiserr$mon,y=thiserr$mon,type="n",log="y",ylim=ylimit,xlab="Month"
             ,ylab="Error count [error/month]",main=letitre,xaxt="n",yaxt="n")
         axis(side=1,at=monat,labels=monlab)
         axis(side=2,at=errat,labels=TRUE)
         if (offegrid) abline(v=monat,h=errat,lty="dashed",col="gray66")
         for (v in 1:ntopoffe){


            offe = thisrank[v]
            #----- Retrieve variable information from the list. ---------------------------#
            vnam        = vades$vnam[offe]

            #----- Retrieve the variable and count the number of offenses. ----------------#
            thisvar     = thiserr[vnam]
            points(x=thiserr$mon,y=t(thisvar),type="o",pch=16,col=coloffe[v],lwd=2,cex=1.2)
         }#end for
         legend(x="top",inset=0.02,legend=lege,col=colleg,pch=16,lty="solid",lwd=2
               ,pt.cex=1.2,cex=0.8,bg="white",ncol=3
               ,title=paste(ntopoffe,"commonest causes for rejection"))

         if (outform[o] == "x11"){
            locator(n=1)
            dev.off()
         }else{
            dev.off()
         }#end if
      } #end for outform
   }#end for error
   #---------------------------------------------------------------------------------------#
}#end for
#------------------------------------------------------------------------------------------#
