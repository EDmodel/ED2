#==========================================================================================#
#==========================================================================================#
#     Leave these commands at the beginning.  They will refresh the session.               #
#------------------------------------------------------------------------------------------#
rm(list=ls())
graphics.off()
#==========================================================================================#
#==========================================================================================#



#==========================================================================================#
#==========================================================================================#
#      Here is the user defined variable section.                                          #
#------------------------------------------------------------------------------------------#

#----- Paths. -----------------------------------------------------------------------------#
here        = getwd()                     # Current directory.
there       = getwd()                     # Directory where analyses/history are 
srcdir      = "/n/home00/mlongo/util/Rsc" # Source  directory.
#------------------------------------------------------------------------------------------#



#----- Output path (background-dependent). ------------------------------------------------#
ibackground = 0                           # Target background colour:
                                          #   0 -- white
                                          #   1 -- black
                                          #   2 -- dark grey
outroot     = file.path(here,paste("agbdlen_ibg",sprintf("%2.2i",ibackground),sep=""))
#------------------------------------------------------------------------------------------#




#----- Time options. ----------------------------------------------------------------------#
monthbeg       = 01           # First month to use
yearbeg        = 1967         # First year to consider
yearend        = 2011         # Maximum year to consider
#------------------------------------------------------------------------------------------#



#----- Plot options. ----------------------------------------------------------------------#
outform        = c("pdf")               # Formats for output file.  Supported formats are:
                                        #   - "X11" - for printing on screen
                                        #   - "eps" - for postscript printing
                                        #   - "png" - for PNG printing
                                        #   - "pdf" - for PDF printing
depth          = 96                     # PNG resolution, in pixels per inch
paper          = "letter"               # Paper size, to define the plot shape
ptsz           = 16                     # Font size.
#------------------------------------------------------------------------------------------#


#----- Name of the simulations. -----------------------------------------------------------#
site             = c("gyf","s67")
site.key         = paste("t",site,sep="")
site.desc        = c("Paracou (GYF)","Santarem (S67)")
drain            = seq(from=0.0,to=-1.6,by=-0.2)
drain.key        = paste("r",sprintf("%+3.3i",100*drain),sep="")
drain.desc       = paste("dR = ",drain,"S",sep="")
dtemp            = seq(from=0.0,to=+3.0,by=+1.0)
dtemp.key        = paste("t",sprintf("%+3.3i",100*dtemp),sep="")
dtemp.desc       = paste("dT = ",dtemp,"K",sep="")
iphen            = c(-1,2)
iphen.key        = paste("iphen",sprintf("%+2.2i",iphen),sep="")
iphen.desc       = paste("Phenology:",c("Evergreen","Drought deciduous"))
realisation      = seq(from=0,to=15,by=1)
n.realisation    = length(realisation)
realisation.key  = paste("real",sprintf("%2.2i",realisation),sep="-")
realisation.desc = paste("Realisation:",realisation)
stext            = c(2,6,8,16,11)
stext.key        = paste("stext",sprintf("%2.2i",stext),sep="")
stext.desc       = c("Loamy Sand","Sandy Clay Loam","Clay Loam","Clayey Sand","Clay")

isite            = 2
idrain           = 6
idtemp           = 1
iiphen           = 2
istext           = 4

myplaces         = paste(site.key[isite],drain.key[idrain],dtemp.key[idtemp]
                        ,realisation.key,iphen.key[iiphen],stext.key[istext]
                        ,sep="_")
simpath          = paste(site.key[isite]  ,drain.key[idrain],dtemp.key[idtemp]
                        ,iphen.key[iiphen],stext.key[istext]
                        ,sep="_")

line.1st         = paste(site.desc[isite],iphen.desc[iiphen],stext.desc[istext],sep="  -  ")
line.2nd         = paste(drain.desc[idrain],realisation.desc,sep="  -  ")
letitre          = paste(line.1st,line.2nd,sep="\n")
#------------------------------------------------------------------------------------------#



#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#      NO NEED TO CHANGE ANYTHING BEYOND THIS POINT UNLESS YOU ARE DEVELOPING THE CODE...  #
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#



#----- Load some packages and scripts. ----------------------------------------------------#
source(file.path(srcdir,"load.everything.r"))
#------------------------------------------------------------------------------------------#


#----- Set how many formats we must output. -----------------------------------------------#
outform = tolower(outform)
nout    = length (outform)
#------------------------------------------------------------------------------------------#

#----- Define plot window size ------------------------------------------------------------#
size = plotsize(proje=FALSE,paper=paper)
#------------------------------------------------------------------------------------------#



#---- Create the main output directory in case there is none. -----------------------------#
outpath = file.path(outroot,simpath)
if (! file.exists(outroot)) dir.create(outroot)
if (! file.exists(outpath)) dir.create(outpath)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Big place loop starts here...                                                        #
#------------------------------------------------------------------------------------------#
for (r in sequence(n.realisation)){

   place = myplaces[r]
   cat(" + Simulation: ",place,"...","\n")


   #----- Retrieve default information about this place and set up some variables. --------#
   thispoi = locations(where=place,here=there,yearbeg=yearbeg,yearend=yearend
                      ,monthbeg=monthbeg)
   inpref  = thispoi$pathin
   lieu    = thispoi$lieu
   iata    = thispoi$iata
   suffix  = thispoi$iata
   #---------------------------------------------------------------------------------------#




   #----- Load the modelled dataset. ------------------------------------------------------#
   cat("   - Loading previous session...","\n")
   ed22.rdata = file.path(here,place,"rdata_month",paste(place,"RData",sep="."))
   load(ed22.rdata)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Load variables of interest.                                                      #
   #---------------------------------------------------------------------------------------#
   when        = datum$when
   agb         = datum$emean$agb
   dlen        = datum$emean$nmon.wdef / 12
   xlimit      = range(when)
   ylimit.agb  = range(agb)
   ylimit.dlen = range(dlen)

   xlabels     = pretty(numyears(when))
   xat         = chron(paste(1,1,xlabels,sep="/"))

   #----- Find the tick marks for both axes. ----------------------------------------------#
   yat.agb     = pretty( agb,n=5)
   yat.dlen    = pretty(dlen,n=5)
   n.yat.agb   = length(yat.agb)
   n.yat.dlen  = length(yat.dlen)
   n.yat       = max(n.yat.agb,n.yat.dlen)
   #----- Make sure that both axes have the same number of tick marks. --------------------#
   if (n.yat.agb > n.yat.dlen){
      dyat.dlen = mean(diff(yat.dlen))
      n.add     = n.yat.agb - n.yat.dlen
      yedge     = rep(range(yat.dlen),times=n.add)[1:n.add]
      ynudge    = (c(-dyat.dlen,dyat.dlen) * rep(sequence(n.add),each=2))[1:n.add]
      yat.dlen  = sort(c(yat.dlen,yedge+ynudge))
   }else if (n.yat.dlen > n.yat.agb){
      dyat.agb  = mean(diff(yat.agb))
      n.add     = n.yat.dlen - n.yat.agb
      yedge     = rep(range(yat.agb),times=n.add)[1:n.add]
      ynudge    = (c(-dyat.agb,dyat.agb) * rep(sequence(n.add),each=2))[1:n.add]
      yat.agb   = sort(c(yat.agb,yedge+ynudge))
   }#end if
   #------ Correct the limits to make both match perfectly. -------------------------------#
   ylimit.agb   = range(yat.agb    )
   ylimit.dlen  = range(yat.dlen   )
   dylimit.agb  = diff (ylimit.agb )
   dylimit.dlen = diff (ylimit.dlen)
   #---------------------------------------------------------------------------------------#



   #----- Convert dry season length to AGB scale. -----------------------------------------#
   dlen.plot = ylimit.agb[1] + (dlen - ylimit.dlen[1] ) * dylimit.agb / dylimit.dlen
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Make plot annotation.                                                            #
   #---------------------------------------------------------------------------------------#
   lex = desc.unit(desc="Year",unit=untab$empty)
   ley.agb  = desc.unit(desc="Above-Ground Biomass",unit=untab$kgcom2)
   ley.dlen = desc.unit(desc="Drought length"      ,unit=untab$yr    )
   #---------------------------------------------------------------------------------------#




   #----- Loop over output formats. -------------------------------------------------------#
   cat("   - Plotting AGB and drought length time series...","\n")
   for (o in sequence(nout)){
      #----- Open file. -------------------------------------------------------------------#
      fichier = file.path(outpath,paste("ts_agb-dlen_",realisation.key[r],"_"
                                       ,site.key[isite],".",outform[o],sep=""))
      if(outform[o] == "x11"){
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
      #------------------------------------------------------------------------------------#



      #----- Plot window based on AGB. ----------------------------------------------------#
      par(par.user)
      par(mar=c(4.1,4.6,4.1,4.6))
      plot.new()
      plot.window(xlim=xlimit,ylim=ylimit.agb)
      abline(h=yat.agb,v=xat,col=grid.colour,lty="dotted")
      axis(side=1,at=xat,labels=xlabels)
      axis(side=2,at=yat.agb,labels=yat.agb ,col.axis=chartreuse.fg)
      axis(side=4,at=yat.agb,labels=yat.dlen,col.axis=indigo.fg    )
      title(main=letitre[r],xlab=lex)
      mtext(text=ley.agb      ,side=2,line=3,col=chartreuse.fg)
      mtext(text=ley.dlen     ,side=4,line=3,col=indigo.fg    )
      lines(x=when,y=dlen.plot,type="l",lwd=2.0,col=indigo.fg    )
      lines(x=when,y=agb      ,type="l",lwd=2.0,col=chartreuse.fg)
      box()
      #------------------------------------------------------------------------------------#



      #----- Close the device. ------------------------------------------------------------#
      if (outform[o] == "x11"){
         locator(n=1)
         dev.off()
      }else{
         dev.off()
      }#end if
      clean.tmp()
      #------------------------------------------------------------------------------------#

   }#end for (o in sequence(nout))
   #---------------------------------------------------------------------------------------#
}#end for (r in sequence(n.realisation))
#------------------------------------------------------------------------------------------#
