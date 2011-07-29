#----- Here is the user-defined variable section. -----------------------------------------#
here           = "thispath"                               # Current directory.
srcdir         = "/n/Moorcroft_Lab/Users/mlongo/util/Rsc" # Property directory.
when           = c("thismontha/thisdatea/thisyeara","thishoura:thisminua:00")               # Time to grab the history
outroot        = "thisoutroot"
myplaces       = c("thispoly")                            # Places to find patch properties


outform        = "png"           # Formats for output file.  Supported formats are:
                                 #   - "X11" - for printing on screen
                                 #   - "eps" - for postscript printing
                                 #   - "png" - for PNG printing
depth          = 96              # PNG resolution, in pixels per inch
paper          = "letter"        # Paper size, to define the plot shape
ptsz           = 14              # Font size.
plotgrid       = TRUE            # Should I plot the grid in the background? 
exfac          = 0.30            # Expansion factor to plot legend.
legwhere       = "topleft"       # Where to put the PFT legend
inset          = 0.01            # Inset
legbg          = "white"         # Colour for legend background
#------------------------------------------------------------------------------------------#
#     List of possible plots. In case you don't want some of them, simply switch plt to F. #
#------------------------------------------------------------------------------------------#
#----- Patch property plots. --------------------------------------------------------------#
npatprop = 11
patprop01 = list(vnam="can.depth"     ,desc="Canopy depth"       ,unit="m"        ,lwd=2
                ,colour="midnightblue",plt=TRUE)
patprop02 = list(vnam="can.age"       ,desc="Patch age"          ,unit="yr"       ,lwd=2
                ,colour="goldenrod"   ,plt=TRUE)
patprop03 = list(vnam="can.area"      ,desc="Patch area"         ,unit="%"        ,lwd=2
                ,colour="lawngreen"   ,plt=TRUE)
patprop04 = list(vnam="patch.lai"     ,desc="Patch LAI"          ,unit="m2/m2"    ,lwd=2
                ,colour="forestgreen" ,plt=TRUE)
patprop05 = list(vnam="patch.wai"     ,desc="Patch WAI"          ,unit="m2/m2"    ,lwd=2
                ,colour="sienna"      ,plt=TRUE)
patprop06 = list(vnam="patch.tai"     ,desc="Patch TAI"          ,unit="m2/m2"    ,lwd=2
                ,colour="darkorange3" ,plt=TRUE)
patprop07 = list(vnam="patch.agb"     ,desc="Patch AGB"          ,unit="kgC/m2"   ,lwd=2
                ,colour="olivedrab"   ,plt=TRUE)
patprop08 = list(vnam="patch.ba"      ,desc="Patch BA"           ,unit="m2"    ,lwd=2
                ,colour="olivedrab"   ,plt=TRUE)
patprop09 = list(vnam="veg.height"    ,desc="Vegetation height"  ,unit="m"     ,lwd=2
                ,colour="midnightblue",plt=TRUE)
patprop10 = list(vnam="veg.displace"  ,desc="Displacement height",unit="m"     ,lwd=2
                ,colour="steelblue"   ,plt=TRUE)
patprop11 = list(vnam="veg.rough"     ,desc="Roughness"          ,unit="m"     ,lwd=2
                ,colour="deepskyblue" ,plt=TRUE)
#----- Cohort property plots. -------------------------------------------------------------#
ncohprop = 8
cohprop01 = list(vnam="lai"   ,desc="LAI"          ,unit="m2/m2"    ,plog=F,lwd=2,plt=T)
cohprop02 = list(vnam="wai"   ,desc="WAI"          ,unit="m2/m2"    ,plog=F,lwd=2,plt=T)
cohprop03 = list(vnam="tai"   ,desc="TAI"          ,unit="m2/m2"    ,plog=F,lwd=2,plt=T)
cohprop04 = list(vnam="cai"   ,desc="CAI"          ,unit="m2/m2"    ,plog=F,lwd=2,plt=T)
cohprop05 = list(vnam="agb"   ,desc="AGB"          ,unit="kgC/plant",plog=T,lwd=2,plt=T)
cohprop06 = list(vnam="height",desc="Height"       ,unit="m"        ,plog=F,lwd=2,plt=T)
cohprop07 = list(vnam="nplant",desc="Plant density",unit="plant/m2" ,plog=T,lwd=2,plt=T)
cohprop08 = list(vnam="ba"    ,desc="Basal area"   ,unit="m2"       ,plog=T,lwd=2,plt=T)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#       NO NEED TO CHANGE ANYTHING BEYOND THIS POINT UNLESS YOU ARE DEVELOPING THE CODE.   #
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#


#----- In case there is some graphic still opened. ----------------------------------------#
graphics.off()
#------------------------------------------------------------------------------------------#



#----- Setting how many formats and places we must output. --------------------------------#
outform = tolower(outform)
nout = length(outform)
nplaces = length(myplaces)
#------------------------------------------------------------------------------------------#

#----- Define some default legend colours and names. --------------------------------------#
pftnames = c("C4 Grass","Early Tropical","Mid Tropical","Late Tropical","Temp. C3 Grass"
             ,"North Pine","South Pine","Late Conifer","Early Temperate","Mid Temperate"
             ,"Late Temperate","C3 Pasture","C3 Crop","C4 Pasture","C4 Crop"
             ,"C3 Grass","Araucaria","Total")
pftcols  = c("gold","chartreuse","limegreen","darkgreen","purple3"
            ,"deepskyblue","aquamarine","midnightblue","darkorange3","sienna"
            ,"firebrick","orchid","coral","gray45","olivedrab"
            ,"goldenrod","steelblue","gray22")
#------------------------------------------------------------------------------------------#


#----- Loading some packages. -------------------------------------------------------------#
library(hdf5)
library(chron)
library(scatterplot3d)
library(lattice)
library(maps)
library(mapdata)
library(akima)
library(fields)
#------------------------------------------------------------------------------------------#



#----- Avoiding unecessary and extremely annoying beeps. ----------------------------------#
options(locatorBell=FALSE)
#------------------------------------------------------------------------------------------#



#----- Loading some files with functions. -------------------------------------------------#
source(paste(srcdir,"atlas.r"     ,sep="/"))
source(paste(srcdir,"globdims.r"  ,sep="/"))
source(paste(srcdir,"locations.r" ,sep="/"))
source(paste(srcdir,"muitas.r"    ,sep="/"))
source(paste(srcdir,"plotsize.r"  ,sep="/"))
source(paste(srcdir,"rconstants.r",sep="/"))
source(paste(srcdir,"sombreado.r" ,sep="/"))
source(paste(srcdir,"southammap.r",sep="/"))
source(paste(srcdir,"timeutils.r" ,sep="/"))
#------------------------------------------------------------------------------------------#



#----- Define plot window size ------------------------------------------------------------#
size = plotsize(proje=FALSE,paper=paper)
#------------------------------------------------------------------------------------------#


#----- Properties for the patch -----------------------------------------------------------#
patprop      = list()
namespprop = NULL
for (s in 1:npatprop){
  ppp          = substring(100+s,2,3)
  prpr         = paste("patprop",ppp,sep="")
  namespprop   = c(namespprop,prpr)
  patprop[[s]] = get(prpr)
} #end for
names(patprop) = namespprop
#------------------------------------------------------------------------------------------#


#----- Properties for the patch -----------------------------------------------------------#
cohprop    = list()
namescprop = NULL
for (s in 1:ncohprop){
  ccc          = substring(100+s,2,3)
  prpr         = paste("cohprop",ccc,sep="")
  namescprop   = c(namescprop,prpr)
  cohprop[[s]] = get(prpr)
} #end for
names(cohprop) = namescprop
#------------------------------------------------------------------------------------------#



#----- Defining the history time to use. --------------------------------------------------#
when = chron(dates=when[1],times=when[2])
#------------------------------------------------------------------------------------------#


for (ipy in 1:nplaces){

   place = myplaces[ipy]

   #----- Retrieve default information about this place and set up some variables. --------#
   thispoi  = locations(where=place,here=here)
   inpref   = paste(here,place,"histo",place,sep="/")
   outpref  = thispoi$pathout
   lieu     = thispoi$lieu
   suffix   = thispoi$iata
   poilon   = thispoi$lon
   poilat   = thispoi$lat

   print(paste(" + Place: ",lieu,"...",sep=""))
   #---------------------------------------------------------------------------------------#



   #---- Extract time information for file name... ----------------------------------------#
   whenoff = when + 1.e-7
   year    = numyears (whenoff)
   month   = nummonths(whenoff)
   day     = numdays  (whenoff)
   hour    = hours    (whenoff)
   minu    = minutes  (whenoff)
   seco    = seconds  (whenoff)
   #---------------------------------------------------------------------------------------#



   #----- Make the time information as fixed-length characters. ---------------------------#
   cyear  = substring(10000+year,2,5)
   cmonth = substring(100+month,2,3)
   cday   = substring(100+day,2,3)
   chour  = substring(100+hour,2,3)
   cminu  = substring(100+minu,2,3)
   cseco  = substring(100+seco,2,3)
   #---------------------------------------------------------------------------------------#



   #----- Build the history file name. ----------------------------------------------------#
   myfile = paste(inpref,"-S-",cyear,"-",cmonth,"-",cday,"-",chour,cminu,cseco,"-g01.h5"
                 ,sep="")
   #---------------------------------------------------------------------------------------#


   #----- Read in the file. ---------------------------------------------------------------#
   print(paste("   - Reading in file: ",myfile,"...",sep=""))
   myhist = hdf5load(myfile,load=FALSE,verbosity=0,tidy=TRUE)
   #---------------------------------------------------------------------------------------#



   #---- Grab properties. -----------------------------------------------------------------#
   print(paste("   - Grabbing properties...",sep=""))
   npatches     = myhist$SIPA.N
   can.depth    = myhist$CAN.DEPTH
   can.area     = myhist$AREA * 100.
   can.age      = myhist$AGE
   veg.height   = myhist$VEG.HEIGHT
   veg.displace = myhist$VEG.DISPLACE
   veg.rough    = myhist$VEG.ROUGH
   patch.lai    = rep(NA,times=npatches)
   patch.wai    = rep(NA,times=npatches)
   patch.tai    = rep(NA,times=npatches)
   patch.agb    = rep(NA,times=npatches)
   patch.ba     = rep(NA,times=npatches)
   patchidx     = seq(from=1,to=npatches,by=1)

   #---- Get cohort properties. -----------------------------------------------------------#
   ncohorts  = myhist$PACO.N
   firstcoh  = myhist$PACO.ID
   cpatch    = list()
   for (ipa in patchidx){
      cpatch[[ipa]] = list()
      if (ncohorts[ipa] > 0){
         icoa = firstcoh[ipa]
         icoz = icoa + ncohorts[ipa] - 1 
         cpatch[[ipa]]$cohidx = 1:ncohorts[ipa]
         cpatch[[ipa]]$pft    = myhist$PFT[icoa:icoz]
         cpatch[[ipa]]$nplant = myhist$NPLANT[icoa:icoz]
         cpatch[[ipa]]$lai    = myhist$LAI.CO[icoa:icoz]
         cpatch[[ipa]]$wai    = myhist$WAI.CO[icoa:icoz]
         cpatch[[ipa]]$tai    = cpatch[[ipa]]$lai + cpatch[[ipa]]$wai
         cpatch[[ipa]]$cai    = myhist$CROWN.AREA[icoa:icoz]
         cpatch[[ipa]]$agb    = myhist$AGB.CO[icoa:icoz]
         cpatch[[ipa]]$ba     = myhist$BA.CO[icoa:icoz]
         cpatch[[ipa]]$height = myhist$HITE[icoa:icoz]
         patch.lai[ipa] = sum(cpatch[[ipa]]$lai)
         patch.wai[ipa] = sum(cpatch[[ipa]]$wai)
         patch.tai[ipa] = sum(cpatch[[ipa]]$tai)
         patch.agb[ipa] = sum(cpatch[[ipa]]$nplant * cpatch[[ipa]]$agb)
         patch.ba [ipa] = sum(cpatch[[ipa]]$nplant * cpatch[[ipa]]$ba)
      }else{
         patch.lai[ipa] = 0.
         patch.wai[ipa] = 0.
         patch.tai[ipa] = 0.
         patch.agb[ipa] = 0.
      }#end if
   }#end for
   #---------------------------------------------------------------------------------------#



   #----- Create the output directory if there is none. -----------------------------------#
   if (! file.exists(outroot)) dir.create(outroot)
   outmain = paste(outroot,place,sep="/")
   if (! file.exists(outmain)) dir.create(outmain)
   outdir = paste(outmain,"patchprop",sep="/")
   if (!file.exists(outdir)) dir.create(outdir)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Loop over the properties, and plot them.                                          #
   #---------------------------------------------------------------------------------------#
   print(paste("   - Plotting patch properties...",sep=""))
   for (pp in 1: npatprop){
      #----- Retrieve variable information from the list. ---------------------------------#
      ppropnow  = patprop[[pp]]
      pvnam     = ppropnow$vnam
      pdesc     = ppropnow$desc
      punit     = ppropnow$unit
      plwd      = ppropnow$lwd
      pcolour   = ppropnow$colour
      plotit    = ppropnow$plt

      #----- Check whether to plot it or skip it. -----------------------------------------#
      if (plotit){

         print(paste("     * ",pdesc,"...",sep=""))

         #----- Get the variable. ---------------------------------------------------------#
         thisvar = get(pvnam)

         #----- Loop over formats. --------------------------------------------------------#
         for (o in 1:nout){
            #----- Make the file name. ----------------------------------------------------#
            fichier = paste(outdir,paste(pvnam,"-",suffix,".",outform[o],sep=""),sep="/")
            if(outform[o] == "x11"){
               X11(width=size$width,height=size$height,pointsize=ptsz)
            }else if(outform[o] == "png"){
               png(filename=fichier,width=size$width*depth,height=size$height*depth
                  ,pointsize=ptsz,res=depth)
            }else if(outform[o] == "eps"){
               postscript(file=fichier,width=size$width,height=size$height
                         ,pointsize=ptsz,paper=paper)
            }#end if
            #------------------------------------------------------------------------------#



            #----- Make the title and the labels. -----------------------------------------#
            letitre = paste(pdesc,lieu,sep=" - ")
            leylab  = paste(pdesc," [",punit,"]",sep="")
            lexlab  = "Patch number"
            #------------------------------------------------------------------------------#



            #----- Plot the graph. --------------------------------------------------------#
            plot(x=patchidx,y=thisvar,type="n",main=letitre,xlab=lexlab,ylab=leylab)
            grid(nx=0,ny=NULL,col="gray59",lty="dotted")
            abline(v=patchidx,col="gray59",lty="dotted")
            points(x=patchidx,y=thisvar,type="h",lwd=plwd,col=pcolour)
            points(x=patchidx,y=thisvar,type="p",pch=16,cex=1.5,col=pcolour)
            #------------------------------------------------------------------------------#



            #----- Close the device. ------------------------------------------------------#
            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
            #------------------------------------------------------------------------------#

         }#end for
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Loop over the cohort properties, and plot them.                                          #
   #---------------------------------------------------------------------------------------#
   print(paste("   - Plotting cohort properties...",sep=""))
   for (cp in 1: ncohprop){
      #----- Retrieve variable information from the list. ---------------------------------#
      cpropnow  = cohprop[[cp]]
      pvnam     = cpropnow$vnam
      pdesc     = cpropnow$desc
      punit     = cpropnow$unit
      plwd      = cpropnow$lwd
      if (cpropnow$plog){
         plog      = "y"
      }else{
         plog      = ""
      }#end if
      plotit    = cpropnow$plt

      #----- Check whether to plot it or skip it. -----------------------------------------#
      if (plotit){

         print(paste("     * ",pdesc,"...",sep=""))

         for (ipa in 1:npatches){
            if (ncohorts[ipa] > 0){

               cipa   = paste("patch_",substring(10000+ipa,2,5),sep="")
               outpat = paste(outdir,cipa,sep="/")
               if (! file.exists(outpat)) dir.create(outpat)

               #----- Get the variable. ---------------------------------------------------#
               thisvar = cpatch[[ipa]][[pvnam]]
               thisidx = cpatch[[ipa]]$cohidx
               pcolour = pftcols[cpatch[[ipa]]$pft]

               #----- Find which PFTs go to the legend. -----------------------------------#
               selpft = seq(from=1,to=npft,by=1) %in% cpatch[[ipa]]$pft
               legpft = pftnames[selpft]
               colpft = pftcols[selpft]

               if (plog == "y"){
                  ylimit  = range(thisvar,na.rm=TRUE)
                  ylimit  = c(ylimit[1]
                             ,exp(log(ylimit[2]) + exfac * log(ylimit[2]/ylimit[1])))
               }else{
                  ylimit  = range(thisvar,na.rm=TRUE)
                  ylimit  = c(ylimit[1],ylimit[2] + exfac * (ylimit[2]-ylimit[1]))
               }#end if


               #----- Loop over formats. --------------------------------------------------#
               for (o in 1:nout){
                  #----- Make the file name. ----------------------------------------------#
                  fichier = paste(outpat,paste(pvnam,"-",suffix,".",outform[o],sep="")
                                        ,sep="/")
                  if(outform[o] == "x11"){
                     X11(width=size$width,height=size$height,pointsize=ptsz)
                  }else if(outform[o] == "png"){
                     png(filename=fichier,width=size$width*depth,height=size$height*depth
                        ,pointsize=ptsz,res=depth)
                  }else if(outform[o] == "eps"){
                     postscript(file=fichier,width=size$width,height=size$height
                               ,pointsize=ptsz,paper=paper)
                  }#end if
                  #------------------------------------------------------------------------#



                  #----- Make the title and the labels. -----------------------------------#
                  letitre = paste(pdesc,paste("Patch",ipa,sep=""),lieu,sep=" - ")
                  leylab  = paste(pdesc," [",punit,"]",sep="")
                  lexlab  = "Cohort number"
                  #------------------------------------------------------------------------#



                  #----- Plot the graph. --------------------------------------------------#
                  plot(x=thisidx,y=thisvar,type="n",main=letitre,log=plog
                      ,xlab=lexlab,ylab=leylab,ylim=ylimit)
                  abline(v=thisidx,h=axTicks(side=2),col="gray59",lty="dotted")
                  points(x=thisidx,y=thisvar,type="h",lwd=plwd,col=pcolour)
                  points(x=thisidx,y=thisvar,type="p",pch=16,cex=1.5,col=pcolour)
                  legend(x=legwhere,inset=inset,bg=legbg,pch=16,pt.cex=1.5,legend=legpft
                        ,cex=0.7,col=colpft,lwd=plwd)
                  #------------------------------------------------------------------------#



                  #----- Close the device. ------------------------------------------------#
                  if (outform[o] == "x11"){
                     locator(n=1)
                     dev.off()
                  }else{
                     dev.off()
                  }#end if
                  #------------------------------------------------------------------------#

               }#end for
               #---------------------------------------------------------------------------#
            }#end if
            #------------------------------------------------------------------------------#
         }#end for
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#
}#end for
#------------------------------------------------------------------------------------------#
