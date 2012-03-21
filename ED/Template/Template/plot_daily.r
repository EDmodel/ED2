#----- Here is the user-defined variable section. -----------------------------------------#
here           = "thispath"      # Current directory.
srcdir         = "/n/moorcroft_data/mlongo/util/Rsc"     
outroot        = "thisoutroot"
daybeg         = thisdatea
monthbeg       = thismontha
yearbeg        = thisyeara         # First year to consider
yearend        = thisyearz         # Maximum year to consider
region         = "thispoly"   # Region name.
myplaces       = c("thispoly")
sasday         = 15
outform        = "png"          # Formats for output file.  Supported formats are:
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
inset          = 0.05           # inset distance between legend and edge of plot region.
legbg          = "white"        # Legend background colour.

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
#----- Time series per PFT. ---------------------------------------------------------------#
ntspft   = 4
tsplot01 = list(vnam="laipft"   ,desc="Leaf area index"           ,unit="m2/m2"     ,plt=T)
tsplot02 = list(vnam="gpppft"   ,desc="Gross primary productivity",unit="kgC/m2/yr" ,plt=T)
tsplot03 = list(vnam="npppft"   ,desc="Net primary productivity"  ,unit="kgC/m2/yr" ,plt=T)
tsplot04 = list(vnam="mcopft"   ,desc="Maintenance costs"         ,unit="kgC/m2/yr" ,plt=T)

#----- Size (DBH) and age structure of cohort level variables. ----------------------------#
npsas  = 10
psas01 = list(vnam="lightbeamco",desc="Downward direct light"   ,unit="--"          , plt=T)
psas02 = list(vnam="gppco"      ,desc="Gross primary product."  ,unit="kgC/plant/yr", plt=T)
psas03 = list(vnam="respco"     ,desc="Total plant respiration" ,unit="kgC/plant/yr", plt=T)
psas04 = list(vnam="nppco"      ,desc="Net primary productivity",unit="kgC/plant/yr", plt=T)
psas05 = list(vnam="mcostco"    ,desc="Maintenance costs"       ,unit="kgC/plant/yr", plt=T)
psas06 = list(vnam="ncbmortco"  ,desc="Mortality due to Neg. CB",unit="1/yr"        , plt=T)
psas07 = list(vnam="agbco"      ,desc="Above-ground biomass"    ,unit="kgC/plant"   , plt=T)
psas08 = list(vnam="fsoco"      ,desc="Fraction of open stomata",unit="--"          , plt=T)
psas09 = list(vnam="nplantco"   ,desc="Plant density"           ,unit="plant/m2"    , plt=T)
psas10 = list(vnam="laico"      ,desc="Leaf area index"         ,unit="m2/m2"       , plt=T)
#----- Similar to Hovmoller diagrams. -----------------------------------------------------#
nhov = 4
hovdi01 = list(vnam   = c("gpp","plresp","hetresp","nep")
              ,desc   = c("GPP","Plant resp.","Het. resp.","NEP")
              ,colour = c("forestgreen","chartreuse","sienna","darkolivegreen")
              ,prefix = "carbflux"
              ,theme  = "Ecosystem carbon fluxes"
              ,unit   = "kgC/m2/yr"
              ,legpos = "topleft"
              ,plt    = TRUE)
hovdi02 = list(vnam   = c("sens","evap","transp")
              ,desc   = c("Surface sensible heat","Surface evaporation","Transpiration")
              ,colour = c("firebrick","midnightblue","chartreuse")
              ,prefix = "surfflux"
              ,theme  = "Surface fluxes"
              ,unit   = "W/m2"
              ,legpos = "topleft"
              ,plt    = TRUE)
hovdi03 = list(vnam   = c("atm.temp","can.temp","leaf.temp","wood.temp","soil.temp")
              ,desc   = c("Atmosphere","Canopy air","Leaf","Wood","Soil (Top)")
              ,colour = c("deepskyblue","gray21","chartreuse","goldenrod","sienna")
              ,prefix = "temperature"
              ,theme  = "Temperature"
              ,unit   = "degC"
              ,legpos = "topleft"
              ,plt    = TRUE)
hovdi04 = list(vnam   = c("atm.shv","can.shv")
              ,desc   = c("Atmosphere","Canopy air")
              ,colour = c("deepskyblue","gray21")
              ,prefix = "h2ovapour"
              ,theme  = "Water vapour mixing ratio"
              ,unit   = "g/kg"
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

#----- Close all graphics devices. --------------------------------------------------------#
graphics.off()

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
source(paste(srcdir,"sombreado.r"  ,sep="/"))
source(paste(srcdir,"southammap.r" ,sep="/"))
source(paste(srcdir,"thermlib.r"   ,sep="/"))
source(paste(srcdir,"timeutils.r"  ,sep="/"))

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
#----- Time series plot. ------------------------------------------------------------------#
tspft = list()
for (s in 1:ntspft){
  sss  = substring(100+s,2,3)
  tsts = paste("tsplot",sss,sep="")
  if (s == 1){
     tspft = get(tsts)
  }else{
     tsmerge = get(tsts)
     for (n in union(names(tspft),names(tsmerge))){
        tspft[[n]] = c(tspft[[n]],tsmerge[[n]])
     } #end for
  }# end if
} #end for
#----- Size (DBH) and age structure cohort level variables. -------------------------------#
sasplot = list()
for (s in 1:npsas){
  sss  = substring(100+s,2,3)
  psps = paste("psas",sss,sep="")
  if (s == 1){
     sasplot = get(psps)
  }else{
     sasmerge = get(psps)
     for (n in union(names(sasplot),names(sasmerge))){
        sasplot[[n]] = c(sasplot[[n]],sasmerge[[n]])
     } #end for
  } #end if
} #end for
#----- Hovmoller diagram ------------------------------------------------------------------#
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
                     monthbeg=monthbeg,daybeg=daybeg,filetype="D",fullonly=fullonly)
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
whena       = chron(dates=paste(monaa,dayaa,yeara,sep="/"))
whenz       = chron(dates=paste(monzz,dayzz,yearz,sep="/"))
totday      = as.integer(round(whenz-whena+1))



#---- Create the main output directory in case there is none. -----------------------------#
if (! file.exists(outroot) dir.create(outroot)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Initialise list...                                                                   #
#------------------------------------------------------------------------------------------#
poi = list()
#------------------------------------------------------------------------------------------#



#----- Printing a banner to entretain the user. ----------------------------------------#
print(paste(" + Initialising the POI structures...",sep=""))
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Big place loop starts here...                                                        #
#------------------------------------------------------------------------------------------#
for (ipy in 1:nplaces){
   place = myplaces[ipy]
   p     = list()

   #----- Retrieve default information about this place and set up some variables. -----#
   thispoi  = locations(where=place,here=here,yearbeg=yearbeg,yearend=yearend,
                        monthbeg=monthbeg,daybeg=daybeg,monthly=FALSE,fullonly=TRUE)
   p$pathmain = paste(outroot,place,sep="/")
   p$outpref  = paste(p$pathmain,"daily",sep="/")
   p$lieu     = thispoi$lieu
   p$suffix   = thispoi$iata
   p$poilon   = thispoi$lon
   p$poilat   = thispoi$lat

   print(paste("   - Place: ",p$lieu,"...",sep=""))

   if (! file.exists(p$pathmain)) dir.create(p$pathmain)
   if (! file.exists(p$outpref))  dir.create(p$outpref)

   #----- PFT arrays. ---------------------------------------------------------------------#
   p$agbpft       = matrix(data=0,nrow=totday,ncol=npft)
   p$laipft       = matrix(data=0,nrow=totday,ncol=npft)
   p$gpppft       = matrix(data=0,nrow=totday,ncol=npft)
   p$npppft       = matrix(data=0,nrow=totday,ncol=npft)
   p$mcopft       = matrix(data=0,nrow=totday,ncol=npft)
   #----- Polygon level vectors. ----------------------------------------------------------#
   p$gpp          = NULL
   p$plresp       = NULL
   p$hetresp      = NULL
   p$npp          = NULL
   p$nep          = NULL
   p$evap         = NULL
   p$transp       = NULL
   p$atm.temp     = NULL
   p$can.temp     = NULL
   p$leaf.temp    = NULL
   p$wood.temp    = NULL
   p$atm.shv      = NULL
   p$can.shv      = NULL
   p$sens         = NULL
   p$agb          = NULL
   p$lai          = NULL
   p$area         = NULL
   p$rain         = NULL
   p$soil.temp    = NULL
   p$soil.moist   = NULL
   #----- Cohort level lists. -------------------------------------------------------------#
   p$lightbeamco  = list()
   p$gppco        = list()
   p$respco       = list()
   p$nppco        = list()
   p$mcostco      = list()
   p$ncbmortco    = list()
   p$agbco        = list()
   p$fsoco        = list()
   p$nplantco     = list()
   p$pftco        = list()
   p$dbhco        = list()
   p$laico        = list()
   p$ageco        = list()
   p$areaco       = list()

   #----- Copy back the list to the poi structure. ----------------------------------------#
   poi[[place]] = p
}#end for

print(" + Looping over years, months, and days...")
n          = 0
d          = 0
thisday    = NULL
firstcall  = TRUE
#----- Looping over years. ----------------------------------------------------------------#
for (year in yeara:yearz){
    if (year == yeara){
       firstmonth = monaa
    }else{
       firstmonth = 1
    }#end if
    if (year == yearz){
       lastmonth = monzz
    }else{
       lastmonth = 12
    }#end if

    #----- Looping over months. -----------------------------------------------------------#
    for (month in firstmonth:lastmonth){

       #----- Define the first day with data. ---------------------------------------------#
       if (month == monaa && year == yeara){
          firstday = dayaa
       }else{
          firstday = 1
       }#end if

       #----- Define the last day with data. ----------------------------------------------#
       if (month == monzz && year == yearz){
          lastday = dayzz
       }else{
          lastday = daymax(month,year)
       }#end if

       refday = min(max(firstday,sasday),lastday)


       #----- Define the inverse of the number of days. -----------------------------------#
       ndaysi = 1./ (lastday - firstday + 1)


       #----- Looping over days. ----------------------------------------------------------#
       print("   - Looping over days...")
       for (day in firstday:lastday){

          d = d + 1

          #----- Building the file name. --------------------------------------------------#
          cmonth = substring(100+month,2,3)
          cday   = substring(100+day,2,3)
          myfile = paste(inpref,"-D-",year,"-",cmonth,"-",cday,"-000000-g01.h5",sep="")
          print (paste("    - Reading data from ",month,"/",day,"/",year,"...",sep=""))

          #----- Read data and close connection immediately after. ------------------------#
          myday = hdf5load(file=myfile,load=FALSE,verbosity=0,tidy=TRUE)

          #----- If this is the first call for this place, find the closest location. -----#
          if (firstcall){
             for (l in 1:nplaces){
                place = myplaces[l]
                reglon   = as.vector(myday$LONGITUDE)
                reglat   = as.vector(myday$LATITUDE)
                distance = rdist.earth(x1=cbind(poi[[place]]$poilon,poi[[place]]$poilat)
                                      ,x2=cbind(reglon,reglat),miles=FALSE)
                poi[[place]]$ipy  = which.min(distance)
                firstcall = FALSE
             }#end for
          }#end if


          #----- Building the time. -------------------------------------------------------#
          thisday = c(thisday,chron(paste(month,day,year,sep="/"),times="0:0:0"))

          #----- Define the number of soil layers. ----------------------------------------#
          nzg       = length(myday$AVG.SOIL.TEMP[ipy,])


          #--------------------------------------------------------------------------------#
          #   Loop over places.                                                            #
          #--------------------------------------------------------------------------------#
          for (l in 1:nplaces){

             place = myplaces[l]
             p     = poi[[place]]
             ipy   = p$ipy

             print (paste("      # Retrieving data from ",p$lieu,"...",sep=""))

             #-----------------------------------------------------------------------------#
             #    Define the position of the current polygon for cohort-level vari-        #
             # ables.                                                                      #
             #-----------------------------------------------------------------------------#
             apa        = myday$SIPA.ID[ipy]
             npatches   = myday$SIPA.N [ipy]
             zpa        = apa + npatches - 1
             aco        = myday$PACO.ID[apa]
             ncohorts   = myday$PACO.N[apa:zpa]
             cohcount   = sum(ncohorts)
             zco        = aco + cohcount -1 

             #----- Loading the simple variables. -----------------------------------------#
             #p$gpp        = c (p$gpp    , myday$DMEAN.GPP         [ipy]            )
             #p$plresp     = c (p$plresp , myday$DMEAN.PLRESP      [ipy]            )
             #p$hetresp    = c (p$hetresp, myday$DMEAN.RH          [ipy]            )
             #p$nep        = c (p$nep    , myday$DMEAN.NEP         [ipy]            )
             #p$npp        = c (p$npp    , myday$DMEAN.GPP         [ipy] - 
             #                             myday$DMEAN.LEAF.RESP   [ipy] - 
             #                             myday$DMEAN.ROOT.RESP   [ipy] -
             #                             myday$DMEAN.GROWTH.RESP [ipy] - 
             #                             myday$DMEAN.STORAGE.RESP[ipy] -
             #                             myday$DMEAN.VLEAF.RESP  [ipy]            )

             p$sens       = c (p$sens  , myday$DMEAN.SENSIBLE.GC    [ipy] 
                                       + myday$DMEAN.SENSIBLE.LC    [ipy]
                                       + myday$DMEAN.SENSIBLE.WC    [ipy])
             p$evap       = c (p$evap  , myday$DMEAN.EVAP           [ipy]
                                       * alvli(myday$DMEAN.CAN.TEMP [ipy])         )
             p$transp     = c (p$transp, myday$DMEAN.TRANSP         [ipy]
                                       * alvli(myday$DMEAN.CAN.TEMP [ipy])         )

             p$atm.temp   = c (p$atm.temp,myday$DMEAN.ATM.TEMP [ipy] - t00         )
             p$atm.shv    = c (p$atm.shv ,myday$DMEAN.ATM.SHV  [ipy] * 1000.       )

             p$can.temp   = c (p$can.temp,myday$DMEAN.CAN.TEMP [ipy] - t00         )
             p$can.shv    = c (p$can.shv ,myday$DMEAN.CAN.SHV  [ipy] * 1000.       )

             p$leaf.temp  = c (p$leaf.temp,myday$DMEAN.LEAF.TEMP [ipy] - t00       )
             p$wood.temp  = c (p$wood.temp,myday$DMEAN.WOOD.TEMP [ipy] - t00       )
             p$rain       = c (p$rain,myday$DMEAN.PCPG         [ipy] * day.sec     )

             #----- Read soil data, and only retrieve the top layer.  ---------------------#
             p$soil.temp  = c (p$soil.temp , myday$DMEAN.SOIL.TEMP [ipy,nzg] - t00)
             p$soil.moist = c (p$soil.moist, myday$DMEAN.SOIL.WATER[ipy,nzg]      )

             #----- Creating a vector of areas the same size as the cohort level. ---------#
             areapa     = myday$AREA[apa:zpa]
             areaconow  = rep(areapa,times=ncohorts)

             #----- Defining the DBH classes. ---------------------------------------------#
             dbhconow   = myday$DBH[aco:zco]

             #----- Defining the age classes. ---------------------------------------------#
             agepa      = myday$AGE[apa:zpa]
             ageconow   = rep(x=agepa,times=ncohorts)

             #----- Load the heterotrophic respiration at the patch level. ----------------#
             hetresppanow  = myday$DMEAN.RH.PA[apa:zpa]

             #----- Load the other cohort-level variables. --------------------------------#
             pftconow       = myday$PFT[aco:zco]
             nplantconow    = myday$NPLANT[aco:zco]
             agbconow       = myday$AGB.CO[aco:zco]
             laiconow       = myday$LAI.CO[aco:zco]
             gppconow       = myday$DMEAN.GPP.CO[aco:zco]
             respconow      = ( myday$DMEAN.LEAF.RESP.CO   [aco:zco]
                              + myday$DMEAN.ROOT.RESP.CO   [aco:zco]
                              + growth.resp.fac[myday$PFT[aco:zco]]
                              * ( myday$DMEAN.GPP.CO      [aco:zco]
                                - myday$DMEAN.LEAF.RESP.CO[aco:zco]
                                - myday$DMEAN.ROOT.RESP.CO[aco:zco]  )
                              )
             nppconow       = gppconow-respconow
             mcostconow     = myday$LEAF.MAINTENANCE [aco:zco] + 
                              myday$ROOT.MAINTENANCE [aco:zco] 
             ncbmortconow   = myday$MORT.RATE.CO[aco:zco,2]
             fsoconow       = myday$DMEAN.FS.OPEN.CO[aco:zco]
             lightbeamconow = myday$DMEAN.LIGHTBEAM.LEVEL[aco:zco]
             #-----------------------------------------------------------------------------#



             #----- Compute the carbon fluxes from the cohort/patch levels. ---------------#
             gpppynow     = sum(gppconow*nplantconow*areaconow ,na.rm=TRUE)
             plresppynow  = sum(respconow*nplantconow*areaconow,na.rm=TRUE)
             hetresppynow = sum(hetresppanow*areapa            ,na.rm=TRUE)
             npppynow     = gpppynow - plresppynow
             neppynow     = npppynow - hetresppynow

             p$gpp        = c (p$gpp    , gpppynow    )
             p$plresp     = c (p$plresp , plresppynow )
             p$hetresp    = c (p$hetresp, hetresppynow)
             p$nep        = c (p$nep    , neppynow    )
             p$npp        = c (p$npp    , npppynow    )
             #-----------------------------------------------------------------------------#


             #-----------------------------------------------------------------------------#
             #     Building the PFT arrays.                                                #
             #-----------------------------------------------------------------------------#
             for (n in 1:npft){
                 sel      = pftconow == n
                 if (any(sel)){
                    p$laipft [d,n]    = p$laipft [d,n] + 
                                        sum(laiconow [sel]    * areaconow[sel])
                    p$agbpft [d,n]    = p$agbpft [d,n] + 
                                        sum(nplantconow[sel] * agbconow [sel]   * 
                                            areaconow[sel])
                    p$gpppft [d,n]    = p$gpppft [d,n] + 
                                        sum(nplantconow[sel] * gppconow [sel]   * 
                                            areaconow[sel])
                    p$npppft [d,n]    = p$npppft [d,n] +
                                        sum(nplantconow[sel] * nppconow [sel]   * 
                                            areaconow[sel])
                    p$mcopft [d,n]    = p$mcopft [d,n] +
                                        sum(nplantconow[sel] * mcostconow [sel] * 
                                            areaconow[sel])
                 }#end if
             }#end if
             #-----------------------------------------------------------------------------#

             #-----------------------------------------------------------------------------#
             #       Building the derived variables.                                       #
             #-----------------------------------------------------------------------------#
             p$lai   = c(p$lai ,sum(p$laipft[d,]))
             p$agb   = c(p$agb ,sum(p$agbpft[d,]))
             #-----------------------------------------------------------------------------#


             #-----------------------------------------------------------------------------#
             #      Building the cohort-level average lists if this is the right month.    #
             #-----------------------------------------------------------------------------#
             monyear  = paste("m",cmonth,"y",year,sep="")
             if (day == firstday){
                #----- Binding the current cohorts. ---------------------------------------#
                p$lightbeamco[[monyear]] = lightbeamconow * ndaysi
                p$gppco      [[monyear]] = gppconow       * ndaysi
                p$respco     [[monyear]] = respconow      * ndaysi
                p$nppco      [[monyear]] = nppconow       * ndaysi
                p$mcostco    [[monyear]] = mcostconow     * ndaysi
                p$ncbmortco  [[monyear]] = ncbmortconow   * ndaysi
                p$agbco      [[monyear]] = agbconow       * ndaysi
                p$fsoco      [[monyear]] = fsoconow       * ndaysi
                p$nplantco   [[monyear]] = nplantconow    * ndaysi
                p$dbhco      [[monyear]] = dbhconow       * ndaysi
                p$laico      [[monyear]] = laiconow       * ndaysi
                #----- The following variables are not averaged. --------------------------#
                p$pftco      [[monyear]] = pftconow
                p$ageco      [[monyear]] = ageconow
                p$areaco     [[monyear]] = areaconow
             }else{
                p$lightbeamco[[monyear]] = p$lightbeamco[[monyear]] + lightbeamconow*ndaysi
                p$gppco      [[monyear]] = p$gppco      [[monyear]] + gppconow      *ndaysi
                p$respco     [[monyear]] = p$respco     [[monyear]] + respconow     *ndaysi
                p$nppco      [[monyear]] = p$nppco      [[monyear]] + nppconow      *ndaysi
                p$mcostco    [[monyear]] = p$mcostco    [[monyear]] + mcostconow    *ndaysi
                p$ncbmortco  [[monyear]] = p$ncbmortco  [[monyear]] + ncbmortconow  *ndaysi
                p$agbco      [[monyear]] = p$agbco      [[monyear]] + agbconow      *ndaysi
                p$fsoco      [[monyear]] = p$fsoco      [[monyear]] + fsoconow      *ndaysi
                p$nplantco   [[monyear]] = p$nplantco   [[monyear]] + nplantconow   *ndaysi
                p$dbhco      [[monyear]] = p$dbhco      [[monyear]] + dbhconow      *ndaysi
                p$laico      [[monyear]] = p$laico      [[monyear]] + laiconow      *ndaysi
             } #end if month=sasmonth
             #-----------------------------------------------------------------------------#

             poi[[place]] = p
          }#end for, place
      }#end for, day
   }#end for, month
}#end for, year


#----- Finding the date variables. --------------------------------------------------------#
thisday = chron(thisday,out.format=c(dates="day-mon-yr",times=NULL))
tlast = length(thisday)

#------------------------------------------------------------------------------------------#
#      Define a suitable scale for those time series that uses thismonth...                #
#------------------------------------------------------------------------------------------#
whenplot = pretty.time(thisday,n=8)

#------------------------------------------------------------------------------------------#
#      Plotting section begins here...                                                     #
#------------------------------------------------------------------------------------------#
print ("  + Plotting figures...")
#------------------------------------------------------------------------------------------#

#----- Finding which PFTs, land uses and transitions we need to consider ---------------#
for (l in 1:nplaces){
   place = myplaces[l]
   p     = poi[[place]]

   print (paste("    - Place: ",p$lieu,"...",sep=""))
   pftave  = colMeans(p$agbpft,na.rm=TRUE)
   selpft  = pftave  > 0.



   #---------------------------------------------------------------------------------------#
   #      Time series by PFT.                                                              #
   #---------------------------------------------------------------------------------------#
   for (v in 1:ntspft){
      vnam        = tspft$vnam[v]
      description = tspft$desc[v]
      unit        = tspft$unit[v]
      plotit      = tspft$plt[v]

      #----- Check whether the user wants to have this variable plotted. ------------------#
      if (plotit){

         #---------------------------------------------------------------------------------#
         #    Checking whether the time series directory exists.  If not, create it.       #
         #---------------------------------------------------------------------------------#
         p$outdir = paste(p$outpref,"tspft",sep="/")
         if (! file.exists(p$outdir)) dir.create(p$outdir)
         print (paste("      +",description,"time series for all PFTs..."))

         #----- Load variable -------------------------------------------------------------#
         thisvar = p[[vnam]]

         #----- Loop over output formats. -------------------------------------------------#
         for (o in 1:nout){
            fichier = paste(p$outdir,"/",vnam,"-",p$suffix,".",outform[o],sep="")
            if(outform[o] == "x11"){
               X11(width=size$width,height=size$height,pointsize=ptsz)
            }else if(outform[o] == "png"){
               png(filename=fichier,width=size$width*depth,height=size$height*depth
                  ,pointsize=ptsz,res=depth)
            }else if(outform[o] == "eps"){
               postscript(file=fichier,width=size$width,height=size$height
                         ,pointsize=ptsz,paper=paper)
            }#end if
            ylimit  = range(thisvar[,selpft],na.rm=TRUE)
            letitre = paste(description,p$lieu,sep=" - ")
            cols    = pftcols[selpft]
            legs    = pftnames[selpft]
            plot(thisday,thisvar[,1],type="n",main=letitre,ylim=ylimit
                ,xlab="Time",xaxt="n",ylab=unit)
            axis(side=1,at=whenplot$levels,labels=whenplot$labels,padj=whenplot$padj)
            if (plotgrid){
               abline(v=whenplot$levels,h=axTicks(side=2),col="lightgray",lty="solid")
            }#end if
            for (n in 1:npft){
               if (selpft[n]){
                  lines(thisday,thisvar[,n],type="l",col=pftcols[n],lwd=lwidth)
               }#end if
            }#end for
            legend(x=legwhere,inset=inset,bg=legbg,legend=legs,col=cols,lwd=lwidth)

            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
         } #end for outform
      }#end if (tseragbpft)
   } #end for tseries
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #    Plotting the 3-D size and age structure of light level.                            #
   #---------------------------------------------------------------------------------------#
   for (v in 1:npsas){
      #----- Retrieve variable information from the list. ---------------------------------#
      vnam        = sasplot$vnam[v]
      description = sasplot$desc[v]
      unit        = sasplot$unit[v]
      plotit      = sasplot$plt[v]

      #----- If this variable is to be plotted, then go through this if block. ------------#
      if (plotit){

         print (paste("      +",description,"size and age structure plot..."))

         #---------------------------------------------------------------------------------#
         #     Check if the directory exists.  If not, create it.                          #
         #---------------------------------------------------------------------------------#
         sasdir = paste(p$outpref,"sas",sep="/")
         if (! file.exists(sasdir)) dir.create(sasdir)
         outdir = paste(sasdir,vnam,sep="/")
         if (! file.exists(outdir)) dir.create(outdir)

         #----- Load this list into "thislist". -------------------------------------------#
         varco  = p[[vnam]]
         nameco = names(varco)

         for (m2y4 in nameco){

            #----- Find which year we are plotting. ---------------------------------------#
            cmonth   = substring(m2y4,2,3)
            mm       = as.numeric(cmonth)
            thisyear = substring(m2y4,5,8)

            #----- Retrieve variable list, age, DBH, and PFT for this year. ---------------#
            agem2y4   =    p$ageco[[m2y4]]
            dbhm2y4   =    p$dbhco[[m2y4]]
            pftm2y4   =    p$pftco[[m2y4]]
            varm2y4   =      varco[[m2y4]]
            popm2y4   = p$nplantco[[m2y4]] * p$areaco[[m2y4]]

            #------------------------------------------------------------------------------#
            #      Find the range.  If the user wants the range to be fixed, then use the  #
            # global range, otherwise, simply use the range for this year.                 #
            #------------------------------------------------------------------------------#
            if (sasfixlimits){
               xlimit  = range(p$ageco            ,na.rm=TRUE)
               ylimit  = range(p$dbhco            ,na.rm=TRUE)
               zlimit  = range(p$varco            ,na.rm=TRUE)
               popmin  = min  (p$nplant * p$areaco,na.rm=TRUE)
               popmax  = max  (p$nplant * p$areaco,na.rm=TRUE)
            }else{
               xlimit  = range(agem2y4  ,na.rm=TRUE)
               ylimit  = range(dbhm2y4  ,na.rm=TRUE)
               zlimit  = range(varm2y4  ,na.rm=TRUE)
               popmin  = min  (popm2y4  ,na.rm=TRUE)
               popmax  = max  (popm2y4  ,na.rm=TRUE)
            }#end if

            #----- Define the scale-dependent population size. ----------------------------#
            cexm2y4 = cexmin + (cexmax - cexmin) * log(popm2y4/popmin) / log(popmax/popmin)

            if (zlimit[1] == zlimit[2]) zlimit = sort(c(0.9,1.1)*zlimit[1])
            if ((zlimit[1] > 0) != (zlimit[2] > 0)){
               floor3d = 0.
            }else if (zlimit[1] > 0){
               floor3d = zlimit[1]
            }else{
               floor3d = zlimit[2]
            }#end if

            #----- Define the grid information for the 3-D plot. --------------------------#
            ageaxis   = pretty(xlimit,n=20)
            dbhaxis   = pretty(ylimit,n=20)
            xlimit    = range(ageaxis)
            ylimit    = range(dbhaxis)
            flooraxis = matrix(floor3d,nrow=length(ageaxis),ncol=length(dbhaxis))

            #----- Expand the lines to make the lollipops. --------------------------------#
            ncohnow  = length(varm2y4)
            agem2y4  = rep(agem2y4,each=3)
            dbhm2y4  = rep(dbhm2y4,each=3)
            pftm2y4  = rep(pftm2y4,each=3)
            varm2y4  = as.vector(rbind(rep(floor3d,times=ncohnow)
                                      ,varco[[m2y4]]
                                      ,rep(NA,times=ncohnow)))
            pchm2y4  = rep(c(NA,16,NA),times=ncohnow)
            cexm2y4  = rep(cexm2y4,each=3)
            colm2y4  = pftcols[pftm2y4]

            pftin   = sort(unique(p$pftco[[m2y4]]))
            colleg  = pftcols[pftin]
            pftleg  = pftnames[pftin]

            #----- Loop over output formats. ----------------------------------------------#
            for (o in 1:nout){
               fichier = paste(outdir,"/",vnam,"-",thisyear,"-",cmonth,"-",
                               p$suffix,".",outform[o],sep="")
               if (outform[o] == "x11"){
                  X11(width=size$width,height=size$height,pointsize=ptsz)
               }else if(outform[o] == "png"){
                  png(filename=fichier,width=size$width*depth,height=size$height*depth
                     ,pointsize=ptsz,res=depth)
               }else if(outform[o] == "eps"){
                  postscript(file=fichier,width=size$width,height=size$height
                            ,pointsize=ptsz,paper=paper)
               }#end if

               #----- Define the labels and colours. --------------------------------------#
               stcol   = pftcols[pftm2y4]
               letitre = paste(description," - ",p$lieu,
                              "\n Time :",mlist[mm],"-",thisyear,sep=" ")
               lezlab  = paste(description," [",unit,"]",sep="")

               #----- First plot: the box. ------------------------------------------------#
               pout = persp(x=ageaxis,y=dbhaxis,z=flooraxis,xlim=xlimit,ylim=ylimit
                           ,zlim=zlimit,theta=theta,phi=phi,col=gcol,expand=expz
                           ,ticktype="detailed",border=NA,xlab="Gap age[yr]"
                           ,ylab="DBH [cm]",zlab=lezlab,shade=shade,ltheta=ltheta
                           ,main=letitre)
               #----- Second plot, the actual data (aka my lollipop trees). ---------------#
               lines (trans3d(x=agem2y4,y=dbhm2y4,z=varm2y4,pmat=pout),type="l"
                     ,col="gray29",lwd=2)
               points(trans3d(x=agem2y4,y=dbhm2y4,z=varm2y4,pmat=pout),type="p"
                     ,pch=pchm2y4,col=colm2y4,cex=cexm2y4)
               legend(x="bottomright",inset=0.01,legend=pftleg,fill=colleg
                     ,ncol=1,bg="white",cex=0.9)
               if (outform[o] == "x11"){
                  locator(n=1)
                  dev.off()
               }else{
                  dev.off()
               }#end if
            } #end for outform
         }#end for nameco
      } #end if
   }#end for npsas
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #   Plot the time series diagrams showing months and years.                             #
   #---------------------------------------------------------------------------------------#
   for (hh in 1:nhov){

      #----- Retrieve variable information from the list. ---------------------------------#
      hovdinow    = hovdi[[hh]]
      vnames      = hovdinow$vnam  
      description = hovdinow$desc  
      lcolours    = hovdinow$colour
      prefix      = hovdinow$prefix
      theme       = hovdinow$theme 
      unit        = hovdinow$unit  
      legpos      = hovdinow$legpos
      plotit      = hovdinow$plt   
 
      if (plotit){


         #----- Define the number of layers. ----------------------------------------------#
         nlayers = length(vnames)
         ylimit  = range(p[vnames],na.rm=TRUE)


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
                            " \n"," Time series of daily averages: ",theme,sep="")

            plot(x=thisday,y=p[[vnames[1]]],type="n",main=letitre,xlab="Day",
                 ylim=ylimit,ylab=paste("[",unit,"]",sep=""))
            if (hovgrid) grid(col="lightgray",lty="dotted")
            for (l in 1:nlayers){
               points(x=thisday,y=p[[vnames[l]]],pch=16,col=lcolours[l],cex=1.2,type="o"
                     ,lwd=2)
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
