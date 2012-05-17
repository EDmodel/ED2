#----- Here is the user-defined variable section. -----------------------------------------#
here           = "thispath"    # Current directory.
there          = "thatpath"    # Directory where analyses/history are 
srcdir         = "/n/moorcroft_data/mlongo/util/Rsc"      # Source  directory.
outroot        = "thisoutroot"
monthbeg       = thismontha
yearbeg        = thisyeara         # First year to consider
yearend        = thisyearz         # Maximum year to consider
season.mona    = thisseasonmona
myplaces       = c("thispoly")
outform        = "png"           # Formats for output file.  Supported formats are:
                                 #   - "X11" - for printing on screen
                                 #   - "eps" - for postscript printing
                                 #   - "png" - for PNG printing

byeold         = TRUE           # Remove old files of the given format?

depth          = 96             # PNG resolution, in pixels per inch
paper          = "letter"       # Paper size, to define the plot shape
ptsz           = 14             # Font size.
lwidth         = 2.5            # Line width
plotgrid       = TRUE           # Should I plot the grid in the background? 

sasfixlimits   = FALSE          # Should I use a fixed scale for size and age-structure
                                # plots? (FALSE will set a suitable scale for each year)

ncolsfc        = 200            # Target number of colours for filled contour plots.
fcgrid         = TRUE           # Should I include a grid on the filled contour plots?

ncolshov       = 200            # Target number of colours for Hovmoller diagrams.
hovgrid        = TRUE           # Should I include a grid on the Hovmoller plots?

legwhere       = "topleft"      # Where should I place the legend?
inset          = 0.01           # inset distance between legend and edge of plot region.
legbg          = "white"        # Legend background colour.
scalleg        = 0.20
cex.main       = 0.8             # Scale coefficient for the title

theta           = 315.                    # Azimuth for perspective projection
phi             = 30.                     # Vertical angle for perspective projection
ltheta          = -210.                   # Azimuth angle for light
shade           = 0.125                   # Shade intensity
expz            = 0.5                     # Expansion factor for Z axis
gcol            = c("lightblue","white")  # Colours for the 50's style floor
cexmin          = 0.5                     # Minimum "head" size of the lollipop
cexmax          = 3.0                     # Maximum "head" size of the lollipop

fracexp         = 0.40                    # Nudging factor for ylimit
ptype           = "l"                     # Type of plot
ptyped          = "p"                     # Type of plot
ptypeb          = "o"                     # Type of plot



#----- Loading some packages. -------------------------------------------------------------#
library(hdf5)
library(chron)
library(scatterplot3d)
library(lattice)
library(maps)
library(mapdata)
library(akima)
library(Hmisc)

#----- In case there is some graphic still opened. ----------------------------------------#
graphics.off()

#----- Set how many formats we must output. -----------------------------------------------#
outform = tolower(outform)
nout = length(outform)

#----- Avoid unecessary and extremely annoying beeps. -------------------------------------#
options(locatorBell=FALSE)


#----- Load some files with functions. ----------------------------------------------------#
source(paste(srcdir,"atlas.r"           ,sep="/"))
source(paste(srcdir,"charutils.r"       ,sep="/"))
source(paste(srcdir,"census.r"          ,sep="/"))
source(paste(srcdir,"cloudy.r"          ,sep="/"))
source(paste(srcdir,"colourmap.r"       ,sep="/"))
source(paste(srcdir,"error.bar.r"       ,sep="/"))
source(paste(srcdir,"globdims.r"        ,sep="/"))
source(paste(srcdir,"locations.r"       ,sep="/"))
source(paste(srcdir,"muitas.r"          ,sep="/"))
source(paste(srcdir,"plotsize.r"        ,sep="/"))
source(paste(srcdir,"pretty.log.r"      ,sep="/"))
source(paste(srcdir,"pretty.time.r"     ,sep="/"))
source(paste(srcdir,"qapply.r"          ,sep="/"))
source(paste(srcdir,"rconstants.r"      ,sep="/"))
source(paste(srcdir,"soilutils.r"       ,sep="/"))
source(paste(srcdir,"sombreado.r"       ,sep="/"))
source(paste(srcdir,"southammap.r"      ,sep="/"))
source(paste(srcdir,"thermlib.r"        ,sep="/"))
source(paste(srcdir,"timeutils.r"       ,sep="/"))
#----- These should be called after the others. --------------------------------------------#
source(paste(srcdir,"pft.coms.r"        ,sep="/"))
source(paste(srcdir,"pyearly_varlist.r" ,sep="/"))


#----- Load observations. -----------------------------------------------------------------#
obsrfile = paste(srcdir,"LBA_MIP.v6.RData",sep="/")
load(file=obsrfile)

#----- Define plot window size ------------------------------------------------------------#
size = plotsize(proje=FALSE,paper=paper)
#------------------------------------------------------------------------------------------#



#---- Create the main output directory in case there is none. -----------------------------#
if (! file.exists(outroot)) dir.create(outroot)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Big place loop starts here...                                                        #
#------------------------------------------------------------------------------------------#
for (place in myplaces){

   #----- Retrieve default information about this place and set up some variables. --------#
   thispoi = locations(where=place,here=there,yearbeg=yearbeg,yearend=yearend
                      ,monthbeg=monthbeg)
   inpref  = thispoi$pathin
   outmain = paste(outroot,place,sep="/")
   outpref = paste(outmain,"yearly",sep="/")
   lieu    = thispoi$lieu
   iata    = thispoi$iata
   suffix  = thispoi$iata
   yeara   = thispoi$yeara
   yearz   = thispoi$yearz
   meszz   = thispoi$monz



   #----- Create the directories in case they don't exist. --------------------------------#
   if (! file.exists(outmain)) dir.create(outmain)
   if (! file.exists(outpref)) dir.create(outpref)

   #----- Print a banner to entretain the user. -------------------------------------------#
   print(paste(" + Post-processing output from ",lieu,"...",sep=""))


   #---------------------------------------------------------------------------------------#
   #     Flush all variables that will hold the data.                                      #
   #---------------------------------------------------------------------------------------#
   totmon      = (yearz-yeara-1)*12+meszz+(12-monthbeg+1)
   #----- Polygon level vectors. ----------------------------------------------------------#
   fast.soil.c     = rep(NA,times=totmon)
   slow.soil.c     = rep(NA,times=totmon)
   struct.soil.c   = rep(NA,times=totmon)
   nep             = rep(NA,times=totmon)
   nee             = rep(NA,times=totmon)
   gpp             = rep(NA,times=totmon)
   npp             = rep(NA,times=totmon)
   plresp          = rep(NA,times=totmon)
   leaf.resp       = rep(NA,times=totmon)
   root.resp       = rep(NA,times=totmon)
   growth.resp     = rep(NA,times=totmon)
   hetresp         = rep(NA,times=totmon)
   reco            = rep(NA,times=totmon)
   mco             = rep(NA,times=totmon)
   npp             = rep(NA,times=totmon)
   cba             = rep(NA,times=totmon)
   ldrop           = rep(NA,times=totmon)
   cflxca          = rep(NA,times=totmon)
   cflxst          = rep(NA,times=totmon)
   evap            = rep(NA,times=totmon)
   transp          = rep(NA,times=totmon)
   ustar           = rep(NA,times=totmon)
   atm.vels        = rep(NA,times=totmon)
   atm.prss        = rep(NA,times=totmon)
   atm.temp        = rep(NA,times=totmon)
   can.prss        = rep(NA,times=totmon)
   can.temp        = rep(NA,times=totmon)
   atm.co2         = rep(NA,times=totmon)
   can.co2         = rep(NA,times=totmon)
   leaf.temp       = rep(NA,times=totmon)
   wood.temp       = rep(NA,times=totmon)
   atm.shv         = rep(NA,times=totmon)
   can.shv         = rep(NA,times=totmon)
   can.co2         = rep(NA,times=totmon)
   hflxca          = rep(NA,times=totmon)
   qwflxca         = rep(NA,times=totmon)
   wflxca          = rep(NA,times=totmon)
   agb             = rep(NA,times=totmon)
   basarea         = rep(NA,times=totmon)
   nplant          = rep(NA,times=totmon)
   lai             = rep(NA,times=totmon)
   wai             = rep(NA,times=totmon)
   tai             = rep(NA,times=totmon)
   rain            = rep(NA,times=totmon)
   gnd.temp        = rep(NA,times=totmon)
   gnd.shv         = rep(NA,times=totmon)
   fs.open         = rep(NA,times=totmon)
   hflxgc          = rep(NA,times=totmon)
   hflxlc          = rep(NA,times=totmon)
   hflxwc          = rep(NA,times=totmon)
   wflxgc          = rep(NA,times=totmon)
   wflxlc          = rep(NA,times=totmon)
   wflxwc          = rep(NA,times=totmon)
   rshort          = rep(NA,times=totmon)
   rshort.beam     = rep(NA,times=totmon)
   rshort.diff     = rep(NA,times=totmon)
   rlong           = rep(NA,times=totmon)
   rshort.gnd      = rep(NA,times=totmon)
   rlong.gnd       = rep(NA,times=totmon)
   rlongup         = rep(NA,times=totmon)
   rlong.albedo    = rep(NA,times=totmon)
   demand          = rep(NA,times=totmon)

   n               = 0
   tomonth         = rep(NA,times=totmon)

   now.month = monthbeg - 1 + 12 * as.integer(monthbeg == 1)
   now.year  = yeara        -      as.integer(monthbeg == 1)

   #----- Loop over years. ----------------------------------------------------------------#
   for (m in 1:totmon){
      now.month = (now.month %% 12) + 1
      now.year  = now.year + as.integer(now.month == 1)

      #----- Build the file name. ---------------------------------------------------------#
      cmonth = sprintf("%2.2i",now.month)
      cyear  = sprintf("%4.4i",now.year )
      ddd    = daymax(now.month,now.year)
      myfile = paste(inpref,"-Q-",cyear,"-",cmonth,"-00-000000-g01.h5",sep="")
      print (paste("    - Reading data from file ",basename(myfile),"...",sep=""))
      #------------------------------------------------------------------------------------#



      #----- Read data and close connection immediately after. ----------------------------#
      mymont = hdf5load(file=myfile,load=FALSE,verbosity=0,tidy=TRUE)
      #------------------------------------------------------------------------------------#


      #----- Build the time. --------------------------------------------------------------#
      tomonth [m] = as.numeric(chron(dates=paste(now.month,1,now.year,sep="/")
                                    ,times="0:0:0") )
      #------------------------------------------------------------------------------------#



      #----- Load the simple variables. ---------------------------------------------------#
      gpp           [m] =   mymont$MMEAN.GPP
      plresp        [m] =   mymont$MMEAN.PLRESP
      leaf.resp     [m] =   mymont$MMEAN.LEAF.RESP
      root.resp     [m] =   mymont$MMEAN.ROOT.RESP
      growth.resp   [m] =   mymont$MMEAN.GROWTH.RESP
      hetresp       [m] =   mymont$MMEAN.RH
      reco          [m] =   plresp[m] + hetresp[m]
      nep           [m] =   mymont$MMEAN.NEP
      nee           [m] = ( - mymont$MMEAN.CARBON.AC + mymont$MMEAN.CARBON.ST ) * yr.sec
      cflxca        [m] = - mymont$MMEAN.CARBON.AC
      cflxst        [m] =   mymont$MMEAN.CARBON.ST
      hflxca        [m] = - mymont$MMEAN.SENSIBLE.AC
      hflxlc        [m] =   mymont$MMEAN.SENSIBLE.LC
      hflxwc        [m] =   mymont$MMEAN.SENSIBLE.WC
      hflxgc        [m] =   mymont$MMEAN.SENSIBLE.GC
      qwflxca       [m] = - mymont$MMEAN.VAPOR.AC        * alvli(mymont$MMEAN.CAN.TEMP)
      wflxca        [m] = - mymont$MMEAN.VAPOR.AC        * day.sec
      wflxlc        [m] =   mymont$MMEAN.VAPOR.LC        * day.sec
      wflxwc        [m] =   mymont$MMEAN.VAPOR.WC        * day.sec
      wflxgc        [m] =   mymont$MMEAN.VAPOR.GC        * day.sec
      evap          [m] =   mymont$MMEAN.EVAP            * day.sec
      transp        [m] =   mymont$MMEAN.TRANSP          * day.sec
      ustar         [m] =   mymont$MMEAN.USTAR
      atm.vels      [m] =   mymont$MMEAN.ATM.VELS
      atm.prss      [m] =   mymont$MMEAN.ATM.PRSS        * 0.01
      atm.temp      [m] =   mymont$MMEAN.ATM.TEMP        - t00
      atm.shv       [m] =   mymont$MMEAN.ATM.SHV         * kg2g
      atm.co2       [m] =   mymont$MMEAN.ATM.CO2           
      can.prss      [m] =   mymont$MMEAN.CAN.PRSS        * 0.01
      can.temp      [m] =   mymont$MMEAN.CAN.TEMP        - t00
      can.shv       [m] =   mymont$MMEAN.CAN.SHV         * kg2g
      can.co2       [m] =   mymont$MMEAN.CAN.CO2           
      gnd.temp      [m] =   mymont$MMEAN.GND.TEMP        - t00
      gnd.shv       [m] =   mymont$MMEAN.GND.SHV         * kg2g
      leaf.temp     [m] =   mymont$MMEAN.LEAF.TEMP        -t00
      wood.temp     [m] =   mymont$MMEAN.WOOD.TEMP        -t00
      rain          [m] =   mymont$MMEAN.PCPG*ddd        * day.sec
      fs.open       [m] =   mymont$MMEAN.FS.OPEN
      rshort        [m] =   mymont$MMEAN.RSHORT
      rshort.beam   [m] =   mymont$MMEAN.RSHORT          - mymont$MMEAN.RSHORT.DIFF
      rshort.diff   [m] =   mymont$MMEAN.RSHORT.DIFF
      rlong         [m] =   mymont$MMEAN.RLONG
      rshort.gnd    [m] =   mymont$MMEAN.RSHORT.GND
      rlong.gnd     [m] =   mymont$MMEAN.RLONG.GND
      rlongup       [m] =   mymont$MMEAN.RLONGUP
      rlong.albedo  [m] =   mymont$MMEAN.RLONG.ALBEDO
      #------------------------------------------------------------------------------------#


      #---- Read in the site-level area. --------------------------------------------------#
      areasi     = mymont$AREA.SI
      npatches   = mymont$SIPA.N
      #------------------------------------------------------------------------------------#


      #----- Read a few patch-level variables. --------------------------------------------#
      areapa     = mymont$AREA * rep(areasi,times=npatches)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Get the soil carbon.                                                           #
      #------------------------------------------------------------------------------------#
      fast.soil.c   [m] = sum(mymont$FAST.SOIL.C       * areapa)
      slow.soil.c   [m] = sum(mymont$SLOW.SOIL.C       * areapa)
      struct.soil.c [m] = sum(mymont$STRUCTURAL.SOIL.C * areapa)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Read the cohort-level variables.  Because empty patchs do exist (deserts),     #
      # we must check whether there is any cohort to be read.  If not, assign NA to        #
      # all variables.                                                                     #
      #------------------------------------------------------------------------------------#
      ncohorts   = mymont$PACO.N
      if (any (ncohorts > 0)){
         areaconow  = rep(areapa,times=ncohorts)
         nplantconow     = mymont$NPLANT
         baconow         = mymont$BA.CO
         agbconow        = mymont$AGB.CO
         laiconow        = mymont$LAI.CO
         waiconow        = mymont$WAI.CO
         taiconow        = laiconow + waiconow
         gppconow        = mymont$MMEAN.GPP.CO
         leafrespconow   = mymont$MMEAN.LEAF.RESP.CO
         rootrespconow   = mymont$MMEAN.ROOT.RESP.CO
         growthrespconow = mymont$MMEAN.GROWTH.RESP.CO
         respconow       = ( mymont$MMEAN.LEAF.RESP.CO   + mymont$MMEAN.ROOT.RESP.CO
                           + mymont$MMEAN.GROWTH.RESP.CO + mymont$MMEAN.STORAGE.RESP.CO
                           + mymont$MMEAN.VLEAF.RESP.CO  )
         nppconow        = gppconow-respconow
         cbaconow        = mymont$MMEAN.CB
         mcostconow      = ( mymont$MMEAN.LEAF.MAINTENANCE
                           + mymont$MMEAN.ROOT.MAINTENANCE )
         ldropconow      = mymont$MMEAN.LEAF.DROP.CO
         demandconow     = mymont$MMEAN.PSI.OPEN * day.sec

         #---------------------------------------------------------------------------------#
      }else{
         #----- Make everything NA. -------------------------------------------------------#
         areaconow       = NA
         nplantconow     = NA
         baconow         = NA
         agbconow        = NA
         laiconow        = NA
         waiconow        = NA
         taiconow        = NA
         gppconow        = NA
         leafrespconow   = NA
         rootrespconow   = NA
         growthrespconow = NA
         respconow       = NA
         nppconow        = NA
         cbaconow        = NA
         mcostconow      = NA
         ldropconow      = NA
         demandconow     = NA
      }#end if
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     The following two variables are used to scale "intensive" properties           #
      # (whatever/plant) to "extensive" (whatever/m2).  Sometimes they may be used to      #
      # build weighted averages.                                                           #
      #------------------------------------------------------------------------------------#
      w.nplant = nplantconow * areaconow
      w.lai    = laiconow    * areaconow
      w.wai    = laiconow    * areaconow
      w.tai    = laiconow    * areaconow
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #       Build the derived variables.                                                 #
      #------------------------------------------------------------------------------------#
      npp     [m] = sum( w.nplant  * nppconow   )
      mco     [m] = sum( w.nplant  * mcostconow )
      cba     [m] = sum( w.nplant  * cbaconow   )
      nplant  [m] = sum( w.nplant               )
      lai     [m] = sum( w.lai                  )
      wai     [m] = sum( w.wai                  )
      tai     [m] = sum( w.tai                  )
      agb     [m] = sum( w.nplant  * agbconow   )
      basarea [m] = sum( w.nplant  * baconow    )
      demand  [m] = sum( w.lai     * demandconow)
      #------------------------------------------------------------------------------------#
   }#end for, year
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Convert tomonth back to chron.                                                     #
   #---------------------------------------------------------------------------------------#
   print (paste("    - Finding times...",sep=""))
   tomonth = chron(tomonth)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Create a list with unique years.                                                   #
   #---------------------------------------------------------------------------------------#
   print (paste("    - Finding the years and seasons...",sep=""))
   year.all = unique(numyears(tomonth))
   if (season.mona == 1){
      nyears    = length(year.all) + 1
      year4     = c(year.all,2*year.all[nyears-1] - year.all[nyears-2])
   }else{
      nyears    = length(year.all)
      year4     = year.all
   }#end if
   year.desc = paste(year4-c(diff(year4)[1],diff(year4)),year4,sep="-")
   year.col  = eft.col[match(year4,eft.year)]
   #----- Year for seasonal means. --------------------------------------------------------#
   yr3mon  = year.all
   nyr3mon = length(yr3mon)
   yr3mon.desc = paste("Dec",sprintf("%2.2i"
                                    ,(yr3mon-c(diff(yr3mon)[1],diff(yr3mon))) %% 100 )
                      ,"-Nov",sprintf("%2.2i",yr3mon %% 100),sep="")
   yr3mon.col  = eft.col[match(yr3mon,eft.year)]
   yr3mon.pch  = eft.pch[match(yr3mon,eft.year)]
   #---------------------------------------------------------------------------------------#


   dwhena     = as.numeric(chron(paste(season.mona,1,year4[2],sep="/")))
   dwhenz     = as.numeric(chron(paste(season.mona,1,year4[3],sep="/")))
   sel        = tomonth >= dwhena & tomonth <= dwhenz
   month.when = tomonth[sel]
   
   for (tt in 1:ntvar){
      thisvar = tvar[[tt]]
      vname   = thisvar$vnam
      desc    = thisvar$desc
      if (vname == "rain"){
         unit    = "mm/month"
      }else{
         unit    = thisvar$unit
      }#end if
      add0             = thisvar$add
      mult0            = thisvar$mult
      colmean          = thisvar$colmean
      colerr           = thisvar$colerr
      coledge          = thisvar$coledge
      plotit           = thisvar$plt
      cumul            = thisvar$cumul
      if (cumul){
         plttype = "l"
      }else{
         plttype = "o"
      }#end if




      if (plotit){
         print(paste("    * ",desc,"..."))


         #----- Create the directories. ---------------------------------------------------#
         outtser  = paste(outpref,"tseries",sep="/")
         outbplot = paste(outpref,"bplot"  ,sep="/")
         if (! file.exists(outtser )) dir.create(outtser )
         if (! file.exists(outbplot)) dir.create(outbplot)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Create the list of data.                                                    #
         #---------------------------------------------------------------------------------#
         outplot = list()
         ylimit  = NULL
         for (y in 2:nyears){
            whena.now      = chron(paste(season.mona,1,year4[y-1],sep="/"))
            whenz.now      = chron(paste(season.mona,1,year4[y]  ,sep="/"))
            sel            = tomonth >= whena.now & tomonth <= whenz.now
            nsel           = sum(sel)
            outplot[[y]]   = list()
            outplot[[y]]$x = month.when
            if (cumul){
               datum               = get(vname)[sel]
               datum[is.na(datum)] = 0
               outplot[[y]]$y = cumsum(datum)
            }else{
               outplot[[y]]$y = get(vname)[sel]
            }#end if

            #----- Update range. ----------------------------------------------------------#
            if (any(is.finite(outplot[[y]]$y)) ){
               y.min  = min(outplot[[y]]$y,na.rm=TRUE)
               y.max  = max(outplot[[y]]$y,na.rm=TRUE)
               ylimit = range(c(ylimit,y.min,y.max))
            }#end if
            #------------------------------------------------------------------------------#
         }#end for
         #----- Make a dummy limit in case everything is empty. ---------------------------#
         if (is.null(ylimit)){
            ylimit = c(-1,1)
         }else{
            ylimit[2] = ylimit[2] + fracexp * (ylimit[2] - ylimit[1])
         }#end if
         #---------------------------------------------------------------------------------#

         #---------------------------------------------------------------------------------#
         #      Make the time axis.                                                        #
         #---------------------------------------------------------------------------------#
         whenplot        = pretty.time(chron(c(dwhena,dwhenz)),n=8)
         whenplot$labels = substring(mlist[nummonths(whenplot$levels)],1,3)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Plot the data.                                                             #
         #---------------------------------------------------------------------------------#
         for (o in 1:nout){
            fichier = paste(outtser,"/ycomp-",vname,"-",suffix,".",outform[o],sep="")
            if(outform[o] == "x11"){
               X11(width=size$width,height=size$height,pointsize=ptsz)
            }else if(outform[o] == "png"){
               png(filename=fichier,width=size$width*depth,height=size$height*depth
                  ,pointsize=ptsz,res=depth)
            }else if(outform[o] == "eps"){
               postscript(file=fichier,width=size$width,height=size$height,pointsize=ptsz
                         ,paper=paper)
            }#end if
            
            letitre = paste(lieu," \n",desc,sep="")
            lex     = "Months"
            ley     = paste(desc," [",unit,"]")


            plot(x=outplot[[2]]$x,y=outplot[[2]]$y,type="n",main=letitre,xlab=lex,ylab=ley
                ,ylim=ylimit,cex.main=cex.main,xaxt="n")
            axis(side=1,at=whenplot$levels,labels=whenplot$labels,padj=whenplot$padj)
            if (plotgrid){ 
               abline(v=whenplot$levels,h=axTicks(side=2),col="gray62",lty="solid")
            }#end if
            for (y in 2:nyears){
               lines(x=outplot[[y]]$x,y=outplot[[y]]$y,type=plttype,pch=16,cex=1.0
                    ,lwd=2.5,col=year.col[y])
            }#end for
            legend(x="topleft",inset=0.01,legend=year.desc[2:nyears]
                  ,lwd=2.5,col=year.col[2:nyears],title="Period",cex=0.9,ncol=3)


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



         #---------------------------------------------------------------------------------#
         #    Split the data into seasonal means.                                          #
         #---------------------------------------------------------------------------------#
         print(paste("    * Seasonal bar plot: ",desc,"..."))
         this.season = season(tomonth,add.year=TRUE)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #    Make the last December become part of the first year.                        #
         #---------------------------------------------------------------------------------#
         sel              = this.season == paste(max(yr3mon)+1,"01",sep="")
         this.season[sel] = paste(min(yr3mon),"01",sep="")
         yr.season        = as.numeric(substring(this.season,1,4))
         #---------------------------------------------------------------------------------#

         #---------------------------------------------------------------------------------#
         #     Find the seasonality matrix.                                                #
         #---------------------------------------------------------------------------------#
         this.var = get(vname)
         if (cumul){
            season.vec = tapply(X=this.var,INDEX=this.season,FUN=sum)
         }else{
            season.vec = tapply(X=this.var,INDEX=this.season,FUN=mean)
         }#end if
         season.mat = matrix( data     = season.vec
                            , ncol     = 4
                            , nrow     = nyr3mon
                            , dimnames = list(yr3mon.desc,season.list)
                            , byrow    = TRUE
                            )#end matrix
         if (vname == "rain"){
            ylimit = c(0,1.25 * max(season.vec,na.rm=TRUE))
         }else{
            yl.try = range(season.vec,na.rm=TRUE)
            if (any(! is.finite(yl.try)) || ( yl.try[1] == yl.try[2] && yl.try[1] == 0)){
               ylimit    = c(-1,1)
               leg.pos   = "topright"
            }else if (yl.try[1] == yl.try[2]){
               ylimit    = yl.try
               ylimit[1] = yl.try[1] - 0.30 * abs(yl.try[1])
               ylimit[2] = yl.try[2] + 0.30 * abs(yl.try[2])
               leg.pos   = "topright"
            }else if(yl.try[2] > 0){
               ylimit    = yl.try
               ylimit[1] = yl.try[1] - 0.05 * (yl.try[2] - yl.try[1])
               ylimit[2] = yl.try[2] + 0.40 * (yl.try[2] - yl.try[1])
               leg.pos   = "topright"
            }else{
               ylimit    = yl.try
               ylimit[1] = yl.try[1] - 0.40 * (yl.try[2] - yl.try[1])
               ylimit[2] = yl.try[2] + 0.05 * (yl.try[2] - yl.try[1])
               leg.pos   = "bottomright"
            }#end if
         }#end if
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Plot the bar plot.                                                          #
         #---------------------------------------------------------------------------------#
         for (o in 1:nout){
            fichier = paste(outbplot,"/bp_season_",vname,"-",suffix,".",outform[o],sep="")
            if(outform[o] == "x11"){
               X11(width=size$width,height=size$height,pointsize=ptsz)
            }else if(outform[o] == "png"){
               png(filename=fichier,width=size$width*depth,height=size$height*depth
                  ,pointsize=ptsz,res=depth)
            }else if(outform[o] == "eps"){
               postscript(file=fichier,width=size$width,height=size$height,pointsize=ptsz
                         ,paper=paper)
            }#end if
            
            letitre = paste(lieu," \n Year comparison: ",desc,sep="")
            lex     = "Season"
            ley     = paste(desc," [",unit,"]")


            barplot(season.mat,col=yr3mon.col,main=letitre,xlab=lex,ylab=ley
                   ,cex.main=cex.main,ylim=ylimit,legend.text=FALSE,beside=TRUE
                   ,border="gray23",xpd=FALSE)
            box()
            legend(x=leg.pos,inset=0.01,legend=yr3mon.desc,fill=yr3mon.col
                  ,title="Period",cex=0.9,ncol=2)


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
   #     Plot the monthly mean variables as functions of other 2 environment variables.    #
   #---------------------------------------------------------------------------------------#
   print(paste(" + Plotting parameter space: ",lieu,"...",sep=""))
   for (xyz in 1:nxyzvar){
      this.xyz = xyzvar[[xyz]]
      zvname   = this.xyz$zvname
      zdesc    = this.xyz$zdesc
      zkey     = this.xyz$zkey
      zunit    = this.xyz$zunit
      xvname   = this.xyz$xvname
      xdesc    = this.xyz$xdesc
      xunit    = this.xyz$xunit
      xleg     = this.xyz$xleg
      yvname   = this.xyz$yvname
      ydesc    = this.xyz$ydesc
      yunit    = this.xyz$yunit
      yleg     = this.xyz$yleg

     #----- Create the directories. -------------------------------------------------------#
     outxyzp = paste(outpref,"xyzplot",sep="/")
     outzvar = paste(outxyzp,zvname   ,sep="/")
     if (! file.exists(outxyzp )) dir.create(outxyzp )
     if (! file.exists(outzvar )) dir.create(outzvar )
     #-------------------------------------------------------------------------------------#


      #----- Number of x and y axes. ------------------------------------------------------#
      nxvars   = length(xvname)
      nyvars   = length(yvname)
      #------------------------------------------------------------------------------------#

      #------------------------------------------------------------------------------------#
      #     Loop over all x and y axes.                                                    #
      #------------------------------------------------------------------------------------#
      print(paste("    * Z: ",zdesc,"..."))
      for (y in 1:nyvars){
         for (x in 1:nxvars){
            print(paste("       ~ X: ",xdesc[x]," Y: ",ydesc[y],"..."))

            #----- Title. -----------------------------------------------------------------#
            letitre    = paste(lieu,zdesc,sep="\n")
            #----- Attribute symbols according to the year. -------------------------------#
            this.pch  = eft.pch[match(yr.season,eft.year)]
            #----- Expand the edges of the x axis. ----------------------------------------#
            xvar      = get(xvname[x])
            lex       = paste(xdesc[x]," [",xunit[x],"]",sep="")
            xrange    = range(xvar,na.rm=TRUE)
            xlimit    = xrange
            xlimit[1] = xrange[1] - 0.05 * (xrange[2] - xrange[1])
            xlimit[2] = xrange[2] + 0.05 * (xrange[2] - xrange[1])
            #----- Expand the edges of the y axis. ----------------------------------------#
            yvar      = get(yvname[y])
            ley       = paste(ydesc[y]," [",yunit[y],"]",sep="")
            yrange    = range(yvar,na.rm=TRUE)
            ylimit    = yrange
            ylimit[1] = yrange[1] - 0.05 * (yrange[2] - yrange[1])
            ylimit[2] = yrange[2] + 0.40 * (yrange[2] - yrange[1])
            #----- Annotation for the colour map ("Z" axis). ------------------------------#
            zvar      = get(zvname)
            lez       = paste(zkey,"\n [",zunit,"]",sep="")
            #----- Find the position to plot the legend. ----------------------------------#
            leg.pos   = paste(yleg[y],xleg[x],sep="")
            #------------------------------------------------------------------------------#
            #     Plot the bar plot.                                                       #
            #------------------------------------------------------------------------------#
            for (o in 1:nout){
               fichier = paste(outzvar,"/cmap_x_",xvname[x],"_y_",yvname[y]
                                       ,"_z_",zvname,"-",suffix,".",outform[o],sep="")
               if(outform[o] == "x11"){
                  X11(width=size$width,height=size$height,pointsize=ptsz)
               }else if(outform[o] == "png"){
                  png(filename=fichier,width=size$width*depth,height=size$height*depth
                     ,pointsize=ptsz,res=depth)
               }else if(outform[o] == "eps"){
                  postscript(file=fichier,width=size$width,height=size$height
                            ,pointsize=ptsz,paper=paper)
               }#end if



               #----- Plot the parameter space. -------------------------------------------#
               colourmap(x=xvar,y=yvar,z=zvar,xlim=xlimit,ylim=ylimit
                        ,colour.palette=muitas,cex=1.6,pch=this.pch,lwd=2
                        ,plot.title=title(main=letitre,xlab=lex,ylab=ley,cex.main=cex.main)
                        ,key.title=title(main=lez,cex.main=0.8)
                        ,plot.axes={ axis(side=1)
                                     axis(side=2)
                                     grid(col="gray62",lty="solid")
                                     legend( x      = leg.pos
                                           , inset  = 0.01
                                           , legend = yr3mon.desc
                                           , col    = "black"
                                           , bg     = "white"
                                           , pch    = yr3mon.pch
                                           , title  = "Period"
                                           , ncol   = 2
                                           , pt.cex = 1./0.9
                                           , cex    = 0.9
                                           )#end legend
                                   }#end plot.axes
                        )#end colourmap
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #     Close the device.                                                     #
               #---------------------------------------------------------------------------#
               if (outform[o] == "x11"){
                  locator(n=1)
                  dev.off()
               }else{
                  dev.off()
               }#end if
               #---------------------------------------------------------------------------#
            }#end for outform
            #------------------------------------------------------------------------------#
         }#end for
         #---------------------------------------------------------------------------------#
      }#end for
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#
}#end for places
#------------------------------------------------------------------------------------------#

#q("no")
