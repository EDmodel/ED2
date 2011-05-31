#----- Here is the user-defined variable section. -----------------------------------------#
here           = "thispath" # Current directory.
srcdir         = "/n/Moorcroft_Lab/Users/mlongo/util/Rsc"
outroot        = "thisoutroot"
frqsum         = 3600.

yearbeg        = thisyeara         # First year to consider
monthbeg       = thismontha
daybeg         = thisdatea
hourbeg        = thishoura
minubeg        = thisminua
secobeg        = 0

yearend        = thisyearz          # Maximum year to consider
monthend       = thismonthz
dayend         = thisdatez
hourend        = thishourz
minuend        = thisminuz
secoend        = 0

region         = "thispoly"   # Region name.
myplaces       = c("thispoly")
outform        = "png"          # Formats for output file.  Supported formats are:
                                 #   - "X11" - for printing on screen
                                 #   - "eps" - for postscript printing
                                 #   - "png" - for PNG printing
ptype          = "l"                  # Type of plot
ptyped         = "p"                  # Type of plot

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
plotgrid        = TRUE          # Should I include a grid on the plots?

legwhere       = "topleft"      # Where should I place the legend?
inset          = 0.05           # inset distance between legend and edge of plot region.
legbg          = "white"        # Legend background colour.
scalleg        = 0.32           # Increase in y scale to fit the legend.

theta           = 315.
phi             = 30.
ltheta          = -210.
shade           = 0.125
expz            = 0.5
gcol            = c("lightblue","white")
cexmin          = 0.5
cexmax          = 3.0
fullonly        = FALSE

#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#     List of possible plots. In case you don't want some of them, simply switch plt to F. #
#------------------------------------------------------------------------------------------#
#----- Time series plots. -----------------------------------------------------------------#
nhov = 9
hovdi01 = list(vnam   = c("gpp","plresp","hetresp")
              ,desc   = c("GPP","Plant resp.","Het. resp.")
              ,colour = c("forestgreen","chartreuse","sienna")
              ,lwd    = c(1.5,1.5,1.5)
              ,type   = ptype
              ,prefix = "carbflux"
              ,theme  = "Ecosystem carbon fluxes"
              ,unit   = "umol/m2/s"
              ,legpos = "topleft"
              ,plt    = TRUE)
hovdi02 = list(vnam   = c("rsnet","rlong","rlongup","qwflxca","hflxca")
              ,desc   = c("Net SW","LW down","LW up","Latent","Sensible")
              ,colour = c("deepskyblue","gray33","goldenrod","midnightblue","firebrick")
              ,lwd    = c(1.5,1.5,1.5,1.5,1.5)
              ,type   = ptype
              ,prefix = "eneflux"
              ,theme  = "Energy fluxes"
              ,unit   = "W/m2"
              ,legpos = "topleft"
              ,plt    = TRUE)
hovdi03 = list(vnam   = c("qwflxgc","qwflxac","qwflxvc","qtransp","qdewgnd")
              ,desc   = c("Ground->Canopy","Air->Canopy","Leaf->Canopy","Transpiration","Dew")
              ,colour = c("firebrick","midnightblue","chartreuse"
                         ,"darkolivegreen","deepskyblue")
              ,lwd    = c(1.5,1.5,1.5,1.5,1.5)
              ,type   = ptype
              ,prefix = "h2oflux"
              ,theme  = "Water fluxes"
              ,unit   = "W/m2"
              ,legpos = "topleft"
              ,plt    = TRUE)
hovdi04 = list(vnam   = c("hflxgc","hflxac","hflxvc")
              ,desc   = c("Ground->Canopy","Air->Canopy","Leaf->Canopy")
              ,colour = c("firebrick","midnightblue","chartreuse")
              ,lwd    = c(1.5,1.5,1.5)
              ,type   = ptype
              ,prefix = "sensflux"
              ,theme  = "Sensible heat fluxes"
              ,unit   = "W/m2"
              ,legpos = "topleft"
              ,plt    = TRUE)
hovdi05 = list(vnam   = c("atm.temp","can.temp","veg.temp","soil.temp")
              ,desc   = c("Atmosphere","Canopy air","Leaf","Soil (Top)")
              ,colour = c("deepskyblue","gray21","chartreuse","sienna")
              ,lwd    = c(1.5,1.5,1.5,1.5)
              ,type   = ptype
              ,prefix = "temperature"
              ,theme  = "Temperature"
              ,unit   = "degC"
              ,legpos = "topleft"
              ,plt    = TRUE)
hovdi06 = list(vnam   = c("atm.shv","can.shv")
              ,desc   = c("Atmosphere","Canopy air")
              ,colour = c("deepskyblue","gray21")
              ,lwd    = c(1.5,1.5)
              ,type   = ptype
              ,prefix = "h2ovapour"
              ,theme  = "Water vapour mixing ratio"
              ,unit   = "g/kg"
              ,legpos = "topleft"
              ,plt    = TRUE)
hovdi07 = list(vnam   = c("atm.co2","can.co2")
              ,desc   = c("Atmosphere","Canopy air")
              ,colour = c("chartreuse","darkolivegreen")
              ,lwd    = c(1.5,1.5)
              ,type   = ptype
              ,prefix = "co2con"
              ,theme  = "CO2 mixing ratio"
              ,unit   = "umol/mol"
              ,legpos = "topleft"
              ,plt    = TRUE)
hovdi08 = list(vnam   = c("soilcp","soilwp","sfldcap","slmsts","soil.water")
              ,desc   = c("Dry soil","Wilting point","Field capacity",
                          "Saturation","Soil Moisture")
              ,colour = c("firebrick","goldenrod","steelblue","midnightblue","olivedrab")
              ,lwd    = c(1.5,1.5,1.5,1.5,1.5)
              ,type   = ptype
              ,prefix = "smoist"
              ,theme  = "Soil moisture at top layer"
              ,unit   = "m3/m3"
              ,legpos = "topleft"
              ,plt    = TRUE)
hovdi09 = list(vnam   = c("prec","intercept","throughfall","wshed")
              ,desc   = c("Precipitation","Interception","Throughfall","Shedding")
              ,colour = c("midnightblue","forestgreen","dodgerblue","aquamarine")
              ,lwd    = c(2.5,1.5,1.5,1.5)
              ,type   = ptype
              ,prefix = "prec"
              ,theme  = "Precipitation rate"
              ,unit   = "mm/hr"
              ,legpos = "topleft"
              ,plt    = TRUE)


#----- Loading some packages. -------------------------------------------------------------#
library(hdf5)
library(chron)
library(scatterplot3d)
library(lattice)
library(maps)
library(mapdata)
library(akima)
library(fields)

#----- In case there is more than one point to be plotted. --------------------------------#

#----- Setting how many formats we must output. -------------------------------------------#
outform = tolower(outform)
nout    = length(outform)
nplaces = length(myplaces)

#----- Avoiding unecessary and extremely annoying beeps. ----------------------------------#
options(locatorBell=FALSE)

#----- Loading some files with functions. -------------------------------------------------#
source(paste(srcdir,"atlas.r"      ,sep="/"))
source(paste(srcdir,"globdims.r"   ,sep="/"))
source(paste(srcdir,"locations.r"  ,sep="/"))
source(paste(srcdir,"muitas.r"     ,sep="/"))
source(paste(srcdir,"plotsize.r"   ,sep="/"))
source(paste(srcdir,"pretty.log.r" ,sep="/"))
source(paste(srcdir,"pretty.time.r",sep="/"))
source(paste(srcdir,"rconstants.r" ,sep="/"))
source(paste(srcdir,"soil_coms.r"  ,sep="/"))
source(paste(srcdir,"sombreado.r"  ,sep="/"))
source(paste(srcdir,"southammap.r" ,sep="/"))
source(paste(srcdir,"timeutils.r"  ,sep="/"))


#----- Convert frqsum to days, which is the standard unit for time in chron. --------------#
frqsum = frqsum / day.sec


#----- Define some default legend colours and names. --------------------------------------#
pftnames = c("C4 Grass","Early Tropical","Mid Tropical","Late Tropical","Temp. C3 Grass"
             ,"North Pine","South Pine","Late Conifer","Early Temperate","Mid Temperate"
             ,"Late Temperate","C3 Pasture","C3 Crop","C4 Pasture","C4 Crop"
             ,"C3 Grass","Araucaria","Total")
pftcols  = c("gold","chartreuse","limegreen","darkgreen","purple3"
            ,"deepskyblue","aquamarine","midnightblue","darkorange3","sienna"
            ,"firebrick","orchid","coral","gray45","olivedrab"
            ,"goldenrod","steelblue","gray22")

lunames = c("Agricultural","Secondary","Primary","Total")
lucols  = c("goldenrod","chartreuse","darkgreen","gray22")

distnames = c("Agr->Agr" ,"2nd->Agr" ,"Prim->Agr"
              ,"Agr->2nd" ,"2nd->2nd" ,"Prim->2nd"
              ,"Agr->Prim","2nd->Prim","Prim->Prim")
distcols  = c("gold","darkorange2","firebrick"
              ,"lightskyblue","turquoise","steelblue"
              ,"palegreen","chartreuse","forestgreen")

#----- Define plot window size ------------------------------------------------------------#
size = plotsize(proje=FALSE,paper=paper)

#------------------------------------------------------------------------------------------#
#     Use brute force to concatenate the plotting lists.                                   #
#------------------------------------------------------------------------------------------#
#----- Time series plots ------------------------------------------------------------------#
hovdi      = list()
nameshovdi = NULL
for (s in 1:nhov){
  sss        = substring(100+s,2,3)
  hdhd       = paste("hovdi",sss,sep="")
  nameshovdi = c(nameshovdi,hdhd)
  hovdi[[s]] = get(hdhd)
} #end for
names(hovdi) = nameshovdi
#------------------------------------------------------------------------------------------#



#----- Retrieve the information about the regional run. -----------------------------------#
regio    = locations(where=region,here=here,yearbeg=yearbeg,yearend=yearend,
                     monthbeg=monthbeg,daybeg=daybeg,monthly=FALSE,fullonly=fullonly)
inpref   = regio$pathin
rstpref  = regio$pathrst
yeara    = regio$yeara
yearz    = regio$yearz
monaa    = regio$mona
monzz    = regio$monz
dayaa    = regio$daya
dayzz    = regio$dayz


#------------------------------------------------------------------------------------------#
#     Flush all variables that will hold the data.                                         #
#------------------------------------------------------------------------------------------#
whena       = chron(dates=paste(monthbeg,daybeg,yearbeg,sep="/")
                   ,times=paste(hourbeg,minubeg,secobeg,sep=":"))
whenz       = chron(dates=paste(monthend,dayend,yearend,sep="/")
                   ,times=paste(hourend,minuend,secoend,sep=":"))
ntimes      = as.integer(round((whenz-whena)/frqsum)+1)





#---- Create the main output directory in case there is none. -----------------------------#
if (! file.exists(outroot)) dir.create(outroot)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Initialise list...                                                                   #
#------------------------------------------------------------------------------------------#
poi = list()

#----- Printing a banner to entretain the user. -------------------------------------------#
print(paste(" + Initialising the POI structures...",sep=""))

for (ipy in 1:nplaces){
   place = myplaces[ipy]
   p     = list()

   #----- Retrieve default information about this place and set up some variables. --------#
   thispoi  = locations(where=place,here=here,yearbeg=yearbeg,yearend=yearend,
                        monthbeg=monthbeg,daybeg=daybeg,monthly=FALSE,fullonly=TRUE)
   p$pathmain = paste(outroot,place,sep="/")
   p$outpref  = paste(p$pathmain,"fast",sep="/")
   p$lieu     = thispoi$lieu
   p$suffix   = thispoi$iata
   p$poilon   = thispoi$lon
   p$poilat   = thispoi$lat

   print(paste("   - Place: ",p$lieu,"...",sep=""))

   if (! file.exists(p$pathmain)) dir.create(p$pathmain)
   if (! file.exists(p$outpref))  dir.create(p$outpref)

   #----- Polygon level vectors. ----------------------------------------------------------#
   p$gpp          = NULL
   p$plresp       = NULL
   p$hetresp      = NULL
   p$rsnet        = NULL
   p$rlong        = NULL
   p$rlongup      = NULL
   p$qwflxca      = NULL
   p$hflxca       = NULL
   p$qwflxgc      = NULL
   p$qwflxac      = NULL
   p$qwflxvc      = NULL
   p$qtransp      = NULL
   p$qdewgnd      = NULL
   p$hflxgc       = NULL
   p$hflxac       = NULL
   p$hflxvc       = NULL
   p$atm.temp     = NULL
   p$can.temp     = NULL
   p$veg.temp     = NULL
   p$soil.temp    = NULL
   p$atm.shv      = NULL
   p$can.shv      = NULL
   p$atm.co2      = NULL
   p$can.co2      = NULL
   p$prec         = NULL
   p$intercept    = NULL
   p$throughfall  = NULL
   p$wshed        = NULL
   p$soilcp       = NULL
   p$soilwp       = NULL
   p$sfldcap      = NULL
   p$slmsts       = NULL
   p$soil.water   = NULL

   #----- Copy back the list to the poi structure. ----------------------------------------#
   poi[[place]] = p
}#end for

print(" + Looping over time...")
when       = whena - frqsum
n          = 0
d          = 0
thistime   = NULL
firstcall  = TRUE
#----- Looping over time. -----------------------------------------------------------------#
for (n in 1:ntimes){
   #---- Add time, and make a time with a tiny offset to avoid the 24:00:00 bug... --------#
   when    = when + frqsum
   whenoff = when + 1.e-7


   year  = numyears (whenoff)
   month = nummonths(whenoff)
   day   = numdays  (whenoff)
   hour  = hours    (whenoff)
   minu  = minutes  (whenoff)
   seco  = seconds  (whenoff)

   #----- Building the file name. ---------------------------------------------------------#
   cyear  = substring(10000+year,2,5)
   cmonth = substring(100+month,2,3)
   cday   = substring(100+day,2,3)
   chour  = substring(100+hour,2,3)
   cminu  = substring(100+minu,2,3)
   cseco  = substring(100+seco,2,3)

   myfile = paste(inpref,"-I-",cyear,"-",cmonth,"-",cday,"-",chour,cminu,cseco
                 ,"-g01.h5",sep="")

   print(paste ("   - Reading the fast file: ",basename(myfile),"...",sep=""))

   #----- Read data and close connection immediately after. -------------------------------#
   myfast = hdf5load(file=myfile,load=FALSE,verbosity=0,tidy=TRUE)
   

   #----- If this is the first call for this place, find the closest location. ------------#
   if (firstcall){
      for (l in 1:nplaces){
         place = myplaces[l]
         reglon   = as.vector(myfast$LONGITUDE)
         reglat   = as.vector(myfast$LATITUDE)
         distance = rdist.earth(x1=cbind(poi[[place]]$poilon,poi[[place]]$poilat)
                               ,x2=cbind(reglon,reglat),miles=FALSE)
         poi[[place]]$ipy  = which.min(distance)
         poi[[place]]$lon  = myfast$LONGITUDE[poi[[place]]$ipy]
         poi[[place]]$lat  = myfast$LATITUDE [poi[[place]]$ipy]
         firstcall = FALSE
      }#end for
   }#end if


   #----- Building the time. --------------------------------------------------------------#
   thistime = c(thistime,chron(dates=paste(month,day,year,sep="/")
                              ,times=paste(hour,minu,seco,sep=":")))

   #----- Define the number of soil layers. -----------------------------------------------#
   nzg       = length(myfast$AVG.SOIL.TEMP[ipy,])


   #---------------------------------------------------------------------------------------#
   #   Loop over places.                                                                   #
   #---------------------------------------------------------------------------------------#
   for (l in 1:nplaces){

      place = myplaces[l]
      p     = poi[[place]]
      ipy   = p$ipy

      print (paste("      # Retrieving data from ",p$lieu,"...",sep=""))

      p$gpp         = c(p$gpp        , myfast$AVG.GPP         [ipy    ]           )
      p$plresp      = c(p$plresp     , myfast$AVG.PLANT.RESP  [ipy    ]           )
      p$hetresp     = c(p$hetresp    , myfast$AVG.HTROPH.RESP [ipy    ]           )
      p$rsnet       = c(p$rsnet      , myfast$AVG.RSHORT      [ipy    ]
                                     * (1. - myfast$AVG.ALBEDT[ipy    ])          )
      p$rlong       = c(p$rlong      , myfast$AVG.RLONG       [ipy    ]           )
      p$rlongup     = c(p$rlongup    , myfast$AVG.RLONGUP     [ipy    ]           )
      p$qwflxca     = c(p$qwflxca    ,-myfast$AVG.VAPOR.AC    [ipy    ] * alvl    )
      p$hflxca      = c(p$hflxca     ,-myfast$AVG.SENSIBLE.AC [ipy    ]           )
      p$qwflxgc     = c(p$qwflxgc    , myfast$AVG.VAPOR.GC    [ipy    ] * alvl    )
      p$qwflxac     = c(p$qwflxac    , myfast$AVG.VAPOR.AC    [ipy    ] * alvl    )
      p$qwflxvc     = c(p$qwflxvc    , myfast$AVG.VAPOR.VC    [ipy    ] * alvl    )
      p$qtransp     = c(p$qtransp    , myfast$AVG.TRANSP      [ipy    ] * alvl    )
      p$qdewgnd     = c(p$qdewgnd    ,-myfast$AVG.DEW.CG      [ipy    ] * alvl    )
      p$hflxgc      = c(p$hflxgc     , myfast$AVG.SENSIBLE.GC [ipy    ]           )
      p$hflxac      = c(p$hflxac     , myfast$AVG.SENSIBLE.AC [ipy    ]           )
      p$hflxvc      = c(p$hflxvc     , myfast$AVG.SENSIBLE.VC [ipy    ]           )
      p$atm.temp    = c(p$atm.temp   , myfast$AVG.ATM.TMP     [ipy    ] - t00     )
      p$can.temp    = c(p$can.temp   , myfast$AVG.CAN.TEMP    [ipy    ] - t00     )
      p$veg.temp    = c(p$veg.temp   , myfast$AVG.VEG.TEMP    [ipy    ] - t00     )
      p$soil.temp   = c(p$soil.temp  , myfast$AVG.SOIL.TEMP   [ipy,nzg] - t00     )
      p$atm.shv     = c(p$atm.shv    , myfast$AVG.ATM.SHV     [ipy    ] * 1000.   )
      p$can.shv     = c(p$can.shv    , myfast$AVG.CAN.SHV     [ipy    ] * 1000.   )
      p$soil.water  = c(p$soil.water , myfast$AVG.SOIL.WATER  [ipy,nzg]           )
      p$atm.co2     = c(p$atm.co2    , myfast$AVG.ATM.CO2     [ipy    ]           )
      p$can.co2     = c(p$can.co2    , myfast$AVG.CAN.CO2     [ipy    ]           )
      p$prec        = c(p$prec       , myfast$AVG.PCPG        [ipy    ] * 3600.   )
      p$intercept   = c(p$intercept  , myfast$AVG.INTERCEPTED [ipy    ] * 3600.   )
      p$throughfall = c(p$throughfall, myfast$AVG.INTERCEPTED [ipy    ] * 3600.   )
      p$wshed       = c(p$wshed      , myfast$AVG.WSHED.VG    [ipy    ] * 3600.   )

      #------ Collecting the properties of this soil type. --------------------------------#
      nsoil = myfast$NTEXT.SOIL[ipy]
      p$soilcp     = c(p$soilcp   , soil[[nsoil]]$soilcp  )
      p$soilwp     = c(p$soilwp   , soil[[nsoil]]$soilwp  )
      p$sfldcap    = c(p$sfldcap  , soil[[nsoil]]$sfldcap )
      p$slmsts     = c(p$slmsts   , soil[[nsoil]]$slmsts  )

      poi[[place]] = p
   }#end for, place
}#end for, year


#----- Finding the date variables. --------------------------------------------------------#
thistime = chron(thistime,out.format=c(dates="day-mon-yr",times="h:m:s"))
tlast = length(thistime)

#------------------------------------------------------------------------------------------#
#      Define a suitable scale for those time series that uses thismonth...                #
#------------------------------------------------------------------------------------------#
whenplot = pretty.time(thistime,n=8)

#------------------------------------------------------------------------------------------#
#      Plotting section begins here...                                                     #
#------------------------------------------------------------------------------------------#
print ("  + Plotting figures...")
#------------------------------------------------------------------------------------------#

#----- Find which PFTs, land uses and transitions we need to consider ---------------------#
for (l in 1:nplaces){
   place = myplaces[l]
   p     = poi[[place]]

   print (paste("    - Place: ",p$lieu,"...",sep=""))



   #---------------------------------------------------------------------------------------#
   #   Plot the time series diagrams showing months and years.                             #
   #---------------------------------------------------------------------------------------#
   for (hh in 1:nhov){

      #----- Retrieve variable information from the list. ---------------------------------#
      hovdinow    = hovdi[[hh]]
      vnames      = hovdinow$vnam  
      description = hovdinow$desc  
      lcolours    = hovdinow$colour
      llwd        = hovdinow$lwd
      ltype       = hovdinow$type
      prefix      = hovdinow$prefix
      theme       = hovdinow$theme 
      unit        = hovdinow$unit  
      legpos      = hovdinow$legpos
      plotit      = hovdinow$plt   
 
      if (plotit){


         #----- Define the number of layers. ----------------------------------------------#
         nlayers = length(vnames)
         ylimit  = range(p[vnames],na.rm=TRUE)
         if (ylimit[1] == ylimit[2]  & ylimit[1] == 0){
            ylimit[1] = -1
            ylimit[2] =  1
         }else if (ylimit[1] == ylimit[2] & ylimit[1] > 0){
            ylimit[2] = (1.0+scalleg) * ylimit[1]
         }else if (ylimit[1] == ylimit[2] & ylimit[1] < 0){
            ylimit[2] = (1.0-scalleg) * ylimit[1]
         }else{
            ylimit[2] = ylimit[2] + scalleg * (ylimit[2] - ylimit[1])
         }#end if


         #---------------------------------------------------------------------------------#
         #     Check if the directory exists.  If not, create it.                          #
         #---------------------------------------------------------------------------------#
         outdir = paste(p$outpref,"tsvar",sep="/")
         if (! file.exists(outdir)) dir.create(outdir)
         print (paste("      +",theme,"Time series ..."))

         #----- Loop over formats. --------------------------------------------------------#
         for (o in 1:nout){
            fichier = paste(outdir,"/",prefix,"-",p$suffix,".",outform[o],sep="")
            if(outform[o] == "x11"){
               X11(width=size$width,height=size$height,pointsize=ptsz)
            }else if(outform[o] == "png"){
               png(filename=fichier,width=size$width*depth,height=size$height*depth
                  ,pointsize=ptsz,res=depth)
            }else if(outform[o] == "eps"){
               postscript(file=fichier,width=size$width,height=size$height
                         ,pointsize=ptsz,paper=paper)
            }#end if

            letitre = paste(theme," - ",p$lieu,
                            " \n"," Time series of hourly averages: ",theme,sep="")

            plot(x=thistime,y=p[[vnames[1]]],type="n",main=letitre,xlab="Day",
                 ylim=ylimit,ylab=paste("[",unit,"]",sep=""),xaxt="n")
            axis(side=1,at=whenplot$levels,labels=whenplot$labels,padj=whenplot$padj)
            if (plotgrid){
               abline(v=whenplot$levels,h=axTicks(side=2),col="lightgray",lty="solid")
            }#end if
            for (l in 1:nlayers){
               points(x=thistime,y=p[[vnames[l]]],pch=16,col=lcolours[l],cex=1.2
                     ,type=ptype,lwd=llwd[l])
            }#end for
            legend(x=legpos,inset=0.05,legend=description,col=lcolours,pch=16,lwd=2)
            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
         } #end for outform
      }#end if plotit
   }#end for nhov
   #---------------------------------------------------------------------------------------#
}#end for places

#q("no")
