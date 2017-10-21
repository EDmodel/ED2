#==========================================================================================#
#==========================================================================================#
#     Leave these commands at the beginning.  They will refresh the session.               #
#------------------------------------------------------------------------------------------#
rm(list=ls())
graphics.off()
#==========================================================================================#
#==========================================================================================#




#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#      Here is the user defined variable section.                                          #
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
here    = getwd()                                #   Current directory
srcdir  = "/prj/prjidfca/marcosl/Util/Rsc"       #   Script directory
outroot = file.path(here,"radprof_comp")         #   Output directory
ascroot = file.path(here,"ascii_radprof")        #   ASCII output directory

#------ Plot options. ---------------------------------------------------------------------#
outform        = c("pdf")              # Formats for output file.  Supported formats are:
                                       #   - "X11" - for printing on screen
                                       #   - "eps" - for postscript printing
                                       #   - "png" - for PNG printing
                                       #   - "pdf" - for PDF printing
depth          = 96                    # PNG resolution, in pixels per inch
paper          = "square"              # Paper size, to define the plot shape
wpaper         = "letter"              # Wide paper size, to define the plot shape
ptsz           = 18                    # Font size.
ibackground    = 0                     # Make figures compatible to which background?
                                       # 0 -- white
                                       # 1 -- black
                                       # 2 -- dark grey
#------------------------------------------------------------------------------------------#




#----- Scale shortwave variables? ---------------------------------------------------------#
scale.sw = FALSE
#------------------------------------------------------------------------------------------#


#----- Profile options. -------------------------------------------------------------------#
my.hours   = sequence(24)-1
prof.hours = c(0,3,6,9,12,15,18,21)                  # Hours to try
prof.zen   = c(0,5,20,30,40,50,60,70,75,80,84,88,92) # Zenith angle breaks
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     List with all sites.  List elements are:                                             #
#     - iata: site code                                                                    #
#     - desc: full name                                                                    #
#     - pch : symbol for this site (currently ignored)                                     #
#     - col : default colour for this site (currently ignored)                             #
#     - year: year to process                                                              #
#                                                                                          #
# The use.sites key allows to subset.  It can be either the number or the iata code.  For  #
#     convenience, you can also set it to TRUE, which will load all sites (FALSE will skip #
#     all of them, kind of useless...)                                                     #
#------------------------------------------------------------------------------------------#
sites      = list()
sites[[1]] = list(iata="gyf",desc="Paracou"       ,yeara=2004,yearz=2004)
sites[[2]] = list(iata="s67",desc="Santarem km 67",yeara=2004,yearz=2004)
sites[[3]] = list(iata="s83",desc="Santarem km 83",yeara=2001,yearz=2001)
sites[[4]] = list(iata="pdg",desc="Pe-de-Gigante" ,yeara=2003,yearz=2003)
sites[[5]] = list(iata="rja",desc="Rebio Jaru"    ,yeara=2001,yearz=2001)
sites[[6]] = list(iata="m34",desc="Manaus K34"    ,yeara=2004,yearz=2004)
sites[[7]] = list(iata="pnz",desc="Petrolina"     ,yeara=2004,yearz=2004)
sites[[8]] = list(iata="ban",desc="Bananal"       ,yeara=2005,yearz=2005)
sites[[9]] = list(iata="dcm",desc="Santarem km 67",yeara=2009,yearz=2009)
use.sites  = "s67" # "rja"
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Height characteristics for output.                                                   #
#------------------------------------------------------------------------------------------#
zdart.bot     = 1
zdart.top     = 53
dzdart        = 1
dart.crit.par = 225 # Critical PAR above which the light utilisation decreases
                    # Value in umol/m2leaf/s, the code will convert to W/m2 if needed
#------------------------------------------------------------------------------------------#





#----- Info on hourly data. ---------------------------------------------------------------#
reload.profile  = TRUE
rdata.path      = file.path(here,"RData_profile")
rdata.suffix    = paste0("radprof_ed22.RData")
finished.suffix = paste0("radprof_ed22.txt")
#------------------------------------------------------------------------------------------#



#----- Info on hourly data. ---------------------------------------------------------------#
reload.summary  = TRUE
rdata.summary   = paste0("summary_radprof_ed22.RData")
#------------------------------------------------------------------------------------------#



#----- Flags to control whether to plot site-specific and multi-site comparison. ----------#
plot.site    = c(FALSE,TRUE)[2]
plot.month   = c(FALSE,TRUE)[2]
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     List all simulations.                                                                #
#------------------------------------------------------------------------------------------#
n           = 0
simul       = list()
n           = n + 1
simul[[ n]] = list( suff           = "ivdyn00_ihrzrad00"
                  , desc           = "HRZ-OFF"
                  , colour         = "#3B24B3"
                  , bgcol          = "#3B24B3"
                  , mult.down      = 1.
                  , mult.beam.down = 1.
                  , mult.diff.down = 1. # 0.97
                  , mult.up        = 1.
                  , angle          = -45
                  , density        = 25
                  )#end list

n           = n + 1
simul[[ n]] = list( suff           = "ivdyn00_ihrzrad01"
                  , desc           = "HRZ-ON"
                  , colour         = "#E65C17"
                  , bgcol          = "#E65C17"
                  , mult.down      = 1. # 0.97
                  , mult.beam.down = 1. # 0.97
                  , mult.diff.down = 1.
                  , mult.up        = 1.015
                  , angle          = +45
                  , density        = 25
                  )#end list
#------------------------------------------------------------------------------------------#


#----- List of bands. ---------------------------------------------------------------------#
band = list()
band[[1]] =  list( key  = "sol"
                 , ed2  = "RSHORT"
                 , desc = "Total Solar"
                 , unit = "wom2"
                 , col  = "#0742C3"
                 )#end list
band[[2]] =  list( key  = "par"
                 , ed2  = "PAR"
                 , desc = "Photosynthetically Active"
                 , unit = "umolom2os"
                 , col  = "#46FF32"
                 )#end list
band[[3]] =  list( key  = "nir"
                 , ed2  = "NIR"
                 , desc = "Near Infrared"
                 , unit = "wom2"
                 , col  = "#FF5700"
                 )#end list
band[[4]] =  list( key  = "tir"
                 , ed2  = "RLONG"
                 , desc = "Thermal Infrared"
                 , unit = "wom2"
                 , col  = "#70000E"
                 )#end list
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#      NO NEED TO CHANGE ANYTHING BEYOND THIS POINT UNLESS YOU ARE DEVELOPING THE CODE...  #
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#



#----- Load some packages. ----------------------------------------------------------------#
source(file.path(srcdir,"load.everything.r"))
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#      Convert lists to a data frames.  For sites only, we may skip several of them using  #
# use.sites.                                                                               #
#------------------------------------------------------------------------------------------#
sites = list.2.data.frame(sites)
if (is.logical(use.sites)){
   use.sites = rep(use.sites,times=nrow(sites))
}else if (is.character(use.sites)){
   use.sites = match(tolower(use.sites),sites$iata)
   use.sites = use.sites[is.finite(use.sites)]
}#end if
sites = sites[use.sites,]
simul = list.2.data.frame(simul)
band  = list.2.data.frame(band )
#------------------------------------------------------------------------------------------#



#----- Number of sites and simulations. ---------------------------------------------------#
nsites  = nrow(sites)
nsimul  = nrow(simul)
nband   = nrow(band )
nprof   = length(prof.hours)
nhours  = length(my.hours)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Height characteristics for output.                                                   #
#------------------------------------------------------------------------------------------#
zdart.above = zdart.top + dzdart
zdart       = seq(from=zdart.above,to=zdart.bot,by=-dzdart)
nzdart      = length(zdart)
#------------------------------------------------------------------------------------------#



#----- Number of output formats. ----------------------------------------------------------#
outform = tolower(outform)
nout    = length(outform)
#------------------------------------------------------------------------------------------#



#----- Define plot window size ------------------------------------------------------------#
size  = plotsize(proje=FALSE,paper=paper )
wsize = plotsize(proje=FALSE,paper=wpaper)
tsize = plotsize(proje=FALSE,stdheight=8.5*6/5,stdwidth=11)
#------------------------------------------------------------------------------------------#




#----- Define sub-plot instructions for simulation comparison. ----------------------------#
lo.simul = pretty.box(n=nsimul)
#------------------------------------------------------------------------------------------#


#----- Critical PAR in W/m2. --------------------------------------------------------------#
dart.crit.par = dart.crit.par * Ein.2.Watts * 1.0e-6
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#      Create all output directories, separated by format.                                 #
#------------------------------------------------------------------------------------------#
if (! file.exists(rdata.path)) dir.create(rdata.path)
if (! file.exists(outroot   )) dir.create(outroot   )
if (! file.exists(ascroot   )) dir.create(ascroot   )
#------------------------------------------------------------------------------------------#


#----- Summary file. ----------------------------------------------------------------------#
summary.fullname = file.path(rdata.path,rdata.summary)
#------------------------------------------------------------------------------------------#


#----- These files are quite heavy, so we load them only if they are necessary. -----------#
if ( (! file.exists(summary.fullname)) || plot.site){
   #---------------------------------------------------------------------------------------#
   #      Create all output directories, separated by format.                              #
   #---------------------------------------------------------------------------------------#
   eft        = list()
   loop.sites = integer(length=0)
   for (p in sequence(nsites)){
      #----- Retrieve site sinformation. --------------------------------------------------#
      iata     = sites$iata[p]
      ip       = match(iata,poilist$iata)
      yeara    = sites$yeara[p]
      yearz    = sites$yearz[p]
      lon      = poilist$lon     [ip]
      lat      = poilist$lat     [ip]
      desc     = poilist$longname[ip]
      short    = poilist$short   [ip]
      #------------------------------------------------------------------------------------#



      #----- Find file name. --------------------------------------------------------------#
      rdata.iata = file.path(rdata.path,paste(iata,rdata.suffix,sep="_"))
      if (file.exists(rdata.iata)){
         #----- Reload data and copy to the general list. ---------------------------------#
         cat(" + Loading data from ",basename(rdata.iata),"...","\n",sep="")
         dummy       = load(file=rdata.iata)
         eft[[iata]] = model
         rm(model)
         #---------------------------------------------------------------------------------#
      }else{
         #---------------------------------------------------------------------------------#
         #      Find all days that we will grab data.                                      #
         #---------------------------------------------------------------------------------#
         cat(" - Generating profiles for ",desc,"...","\n")
         alldays   = chron(seq( from = as.numeric(chron(paste( 1, 1,yeara,sep="/")) )
                              , to   = as.numeric(chron(paste(12,31,yearz,sep="/")) )
                              , by   = 1
                              )#end seq
                          )#end chron
         ndates    = length(alldays)
         when      = chron( dates = rep(paste(alldays),each=nhours)
                          , times = paste(rep(my.hours,times=ndates),0,0,sep=":")
                          )#end chron
         w.year    = numyears (when)
         w.month   = nummonths(when)
         w.day     = numdays  (when)
         w.hour    = hours    (when)
         w.tomonth = chron(paste(w.month,1,w.year,sep="/"))
         ntimes    = sort(unique(w.tomonth))


         zenith   = mapply( FUN      = function(bef,now,...){
                                          zz  = ed.zen(when=c(bef,now),...)
                                          ans = c( zen   = zz$zen  [2]
                                                 , day   = zz$day  [2]
                                                 , night = zz$night[2]
                                                 )#end c
                                          return(ans)
                                       }#end function
                          , bef      = when - 1/24
                          , now      = when
                          , MoreArgs = list( lon       = lon
                                           , lat       = lat
                                           , ed21      = TRUE
                                           , zeronight = FALSE
                                           , meanval   = TRUE
                                           , imetavg   = 1
                                           , nmean     = 15
                                           )#end list
                          , SIMPLIFY = TRUE
                          )#end mapply

         zen   = tapply( X     = zenith["zen",]
                       , INDEX = list(w.tomonth,w.hour)
                       , FUN   = mean
                       , na.rm = TRUE
                       )#end tapply
         day   = tapply( X     = zenith["day",]
                       , INDEX = list(w.tomonth,w.hour)
                       , FUN   = commonest
                       , na.rm = TRUE
                       )#end tapply
         night = tapply( X     = zenith["night",]
                       , INDEX = list(w.tomonth,w.hour)
                       , FUN   = commonest
                       , na.rm = TRUE
                       )#end tapply

         model    = list( iata         = iata
                        , short        = short
                        , longname     = desc
                        , lon          = lon
                        , lat          = lat
                        , yeara        = yeara
                        , yearz        = yearz
                        , sresume      = 1
                        , tomonth      = sort(unique(w.tomonth))
                        , month        = nummonths(sort(unique(w.tomonth)))
                        , year         = numyears (sort(unique(w.tomonth)))
                        , hour         = sort(unique(w.hour))
                        , ntimes       = length(unique(w.tomonth))
                        , ndcycle      = length(unique(w.hour))
                        , nzdart       = nzdart
                        , zdart        = zdart
                        , zen          = zen
                        , diel         = round( 1 + day + night) + 0.*night
                        , gloom        = matrix( data = FALSE
                                               , ncol = ncol(night)
                                               , nrow = nrow(night)
                                               )#end matrix
                        )#end list
         #---------------------------------------------------------------------------------#

         eft[[iata]] = model

      }#end if
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#





   #=======================================================================================#
   #=======================================================================================#
   #      Loop over sites.                                                                 #
   #---------------------------------------------------------------------------------------#
   for (p in sequence(nsites)){

      #------ Grab data for this site. ----------------------------------------------------#
      iata      = sites$iata[p]
      model     = eft[[iata]]
      yeara     = model$yeara
      yearz     = model$yearz
      lon       = model$lon
      lat       = model$lat
      longname  = model$longname
      short     = model$short
      sresume   = model$sresume
      ntimes    = model$ntimes
      ndcycle   = model$ndcycle
      tomonth   = model$tomonth
      #------------------------------------------------------------------------------------#


      #----- File name for this site. -----------------------------------------------------#
      rdata.iata = file.path(rdata.path,paste(iata,rdata.suffix,sep="_"))
      #------------------------------------------------------------------------------------#

      #---- Get the times that must be read. ----------------------------------------------#
      if (sresume > nsimul){
         loop.simul = numeric(0)
      }else{
         loop.simul = seq(from=sresume,to=nsimul,by=1)
      }#end if
      #------------------------------------------------------------------------------------#







      #------------------------------------------------------------------------------------#
      #    Loop over simulations.                                                          #
      #------------------------------------------------------------------------------------#
      for (s in loop.simul){

         cat("   * Simulation ",simul$desc[s],":","\n",sep="")

         #----- Simulation prefix. --------------------------------------------------------#
         simplace = paste("t",iata,"_",simul$suff[s],sep="")
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Loop over dates.                                                            #
         #---------------------------------------------------------------------------------#
         for (w in sequence(ntimes)){


            #----- Grab date/time information on the data set. ----------------------------#
            now   = tomonth[w]
            yyyy = sprintf("%4.4i",numyears (now))
            mm   = sprintf("%2.2i",nummonths(now))
            dd   = sprintf("%2.2i",numdays  (now))
            hh   = sprintf("%2.2i",hours    (now))
            cat("     > ",paste(now),"\n",sep="")
            #------------------------------------------------------------------------------#



            #----- Build file suffix. -----------------------------------------------------#
            tsuff = paste(yyyy,"-",mm,"-00-000000-g01.h5",sep="")
            #------------------------------------------------------------------------------#

            #----- Build the file name. ---------------------------------------------------#
            filenow  = paste(simplace,"Q",tsuff,sep="-")
            #------------------------------------------------------------------------------#


            #----- Make file name as well as the compressed names. ------------------------#
            fullfile     = file.path(here,simplace,"analy",filenow)
            fullfile.gz  = paste(fullfile,"gz" ,sep=".")
            fullfile.bz2 = paste(fullfile,"bz2",sep=".")
            #------------------------------------------------------------------------------#



            #----- Open the file. ---------------------------------------------------------#
            if (file.exists(fullfile)){
               mymont    = hdf5load(file=fullfile,load=FALSE,tidy=TRUE,verbosity=0)
            }else if (file.exists(fullfile.gz )){
               temp.file = file.path(tempdir(),basename(fullfile))
               dummy     = gunzip(filename=fullfile.gz,destname=temp.file,remove=FALSE)
               mymont    = hdf5load(file=temp.file,load=FALSE,verbosity=0,tidy=TRUE)
               dummy     = file.remove(temp.file)
            }else if (file.exists(fullfile.bz2)){
               temp.file = file.path(tempdir(),basename(fullfile))
               dummy     = bunzip2(filename=fullfile.bz2,destname=temp.file,remove=FALSE)
               mymont    = hdf5load(file=temp.file,load=FALSE,verbosity=0,tidy=TRUE)
               dummy     = file.remove(temp.file)
            }else{
               stop ("File not found")
            }#end if
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Append some auxiliary variables to mymont.                              #
            #------------------------------------------------------------------------------#
            mymont$SIPA.N                   = length(mymont$PACO.ID)
            mymont$QMEAN.ATM.NIR.PY         = ( mymont$QMEAN.ATM.RSHORT.PY
                                              - mymont$QMEAN.ATM.PAR.PY         )
            mymont$QMEAN.ATM.NIR.DIFF.PY    = ( mymont$QMEAN.ATM.RSHORT.DIFF.PY
                                              - mymont$QMEAN.ATM.PAR.DIFF.PY    )
            mymont$QMEAN.ATM.RSHORT.BEAM.PY = ( mymont$QMEAN.ATM.RSHORT.PY
                                              - mymont$QMEAN.ATM.RSHORT.DIFF.PY )
            mymont$QMEAN.ATM.PAR.BEAM.PY    = ( mymont$QMEAN.ATM.PAR.PY
                                              - mymont$QMEAN.ATM.PAR.DIFF.PY    )
            mymont$QMEAN.ATM.NIR.BEAM.PY    = ( mymont$QMEAN.ATM.NIR.PY
                                              - mymont$QMEAN.ATM.NIR.DIFF.PY    )
            mymont$QMEAN.ATM.RLONG.DIFF.PY  = mymont$QMEAN.ATM.RLONG.PY
            mymont$QMEAN.ATM.RLONG.BEAM.PY  = 0. * mymont$QMEAN.ATM.RLONG.PY
            mymont$QMEAN.ATM.NIR.SI         = ( mymont$QMEAN.ATM.RSHORT.SI
                                              - mymont$QMEAN.ATM.PAR.SI         )
            mymont$QMEAN.ATM.NIR.DIFF.SI    = ( mymont$QMEAN.ATM.RSHORT.DIFF.SI
                                              - mymont$QMEAN.ATM.PAR.DIFF.SI    )
            mymont$QMEAN.ATM.RSHORT.BEAM.SI = ( mymont$QMEAN.ATM.RSHORT.SI
                                              - mymont$QMEAN.ATM.RSHORT.DIFF.SI )
            mymont$QMEAN.ATM.PAR.BEAM.SI    = ( mymont$QMEAN.ATM.PAR.SI
                                              - mymont$QMEAN.ATM.PAR.DIFF.SI    )
            mymont$QMEAN.ATM.NIR.BEAM.SI    = ( mymont$QMEAN.ATM.NIR.SI
                                              - mymont$QMEAN.ATM.NIR.DIFF.SI    )
            mymont$QMEAN.ATM.RLONG.DIFF.SI  = mymont$QMEAN.ATM.RLONG.SI
            mymont$QMEAN.ATM.RLONG.BEAM.SI  = 0. * mymont$QMEAN.ATM.RLONG.SI
            mymont$QMEAN.RSHORTUP.PA        = ( mymont$QMEAN.PARUP.PA 
                                              + mymont$QMEAN.NIRUP.PA )
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Make the radiation at the top of the atmosphere a patch variables.       #
            #------------------------------------------------------------------------------#
            mymont$QMEAN.PAR.BEAM.PA    = with( mymont
                                              , matrix( rep( x    = QMEAN.ATM.PAR.BEAM.SI
                                                           , each = SIPA.N
                                                           )#end rep
                                                      , nrow = SIPA.N, ncol = ndcycle
                                                      )#end matrix
                                              )#end with
            mymont$QMEAN.PAR.DIFF.PA    = with( mymont
                                              , matrix( rep( x    = QMEAN.ATM.PAR.DIFF.SI
                                                           , each = SIPA.N
                                                           )#end rep
                                                      , nrow = SIPA.N, ncol = ndcycle
                                                      )#end matrix
                                              )#end with
            mymont$QMEAN.NIR.BEAM.PA    = with( mymont
                                              , matrix( rep( x    = QMEAN.ATM.NIR.BEAM.SI
                                                           , each = SIPA.N
                                                           )#end rep
                                                      , nrow = SIPA.N, ncol = ndcycle
                                                      )#end matrix
                                              )#end with
            mymont$QMEAN.NIR.DIFF.PA    = with( mymont
                                              , matrix( rep( x    = QMEAN.ATM.NIR.DIFF.SI
                                                           , each = SIPA.N
                                                           )#end rep
                                                      , nrow = SIPA.N, ncol = ndcycle
                                                      )#end matrix
                                              )#end with
            mymont$QMEAN.RSHORT.BEAM.PA = with( mymont
                                              , matrix( rep( x    = QMEAN.ATM.RSHORT.BEAM.SI
                                                        , each = SIPA.N
                                                        )#end rep
                                                      , nrow = SIPA.N, ncol = ndcycle
                                                      )#end matrix
                                              )#end with
            mymont$QMEAN.RSHORT.DIFF.PA = with( mymont
                                              , matrix( rep( x    = QMEAN.ATM.RSHORT.DIFF.SI
                                                        , each = SIPA.N
                                                        )#end rep
                                                      , nrow = SIPA.N, ncol = ndcycle
                                                      )#end matrix
                                              )#end with
            mymont$QMEAN.RLONG.BEAM.PA  = with( mymont
                                              , matrix( rep( x    = QMEAN.ATM.RLONG.BEAM.SI
                                                         , each = SIPA.N
                                                         )#end rep
                                                      , nrow = SIPA.N, ncol = ndcycle
                                                      )#end matrix
                                              )#end with
            mymont$QMEAN.RLONG.DIFF.PA  = with( mymont
                                              , matrix( rep( x    = QMEAN.ATM.RLONG.DIFF.SI
                                                         , each = SIPA.N
                                                         )#end rep
                                                      , nrow = SIPA.N, ncol = ndcycle
                                                      )#end matrix
                                              )#end with
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Correct profiles.                                                       #
            #------------------------------------------------------------------------------#
            mymont$QMEAN.PAR.BEAM.CO    =   mymont$QMEAN.RAD.PROFILE[,, 1]
            mymont$QMEAN.PAR.DIFF.CO    =   mymont$QMEAN.RAD.PROFILE[,, 3]
            mymont$QMEAN.PARUP.CO       = ( mymont$QMEAN.RAD.PROFILE[,, 2]
                                          + mymont$QMEAN.RAD.PROFILE[,, 4] )
            mymont$QMEAN.NIR.BEAM.CO    =   mymont$QMEAN.RAD.PROFILE[,, 5]
            mymont$QMEAN.NIR.DIFF.CO    =   mymont$QMEAN.RAD.PROFILE[,, 7]
            mymont$QMEAN.NIRUP.CO       = ( mymont$QMEAN.RAD.PROFILE[,, 6]
                                          + mymont$QMEAN.RAD.PROFILE[,, 8] )
            mymont$QMEAN.RLONG.BEAM.CO  =   mymont$QMEAN.RAD.PROFILE[,, 9] * 0.
            mymont$QMEAN.RLONG.DIFF.CO  =   mymont$QMEAN.RAD.PROFILE[,, 9]
            mymont$QMEAN.RLONGUP.CO     =   mymont$QMEAN.RAD.PROFILE[,,10]
            mymont$QMEAN.RSHORT.BEAM.CO = ( mymont$QMEAN.PAR.BEAM.CO
                                          + mymont$QMEAN.NIR.BEAM.CO
                                          )#end
            mymont$QMEAN.RSHORT.DIFF.CO = ( mymont$QMEAN.PAR.DIFF.CO
                                          + mymont$QMEAN.NIR.DIFF.CO
                                          )#end
            mymont$QMEAN.RSHORTUP.CO    = ( mymont$QMEAN.PARUP.CO
                                          + mymont$QMEAN.NIRUP.CO
                                          )#end
            mymont$QMEAN.NIR.L.CO       = ( mymont$QMEAN.RSHORT.L.CO
                                          - mymont$QMEAN.PAR.L.CO
                                          )#end
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Initialise the list with simulation data.                                #
            #------------------------------------------------------------------------------#
            if (w == 1 && s == 1){
               #----- List with model simulations. ----------------------------------------#
               model$simul = replicate(nsimul,list())
               names(model$simul) = simul$suff
               #---------------------------------------------------------------------------#


               #----- Interpolated data to DART levels. -----------------------------------#
               dart.ai     = array (data= NA,dim=c(nsimul,ntimes        ,nzdart))
               dart.atm    = array (data= NA,dim=c(nsimul,ntimes,ndcycle       ))
               dart.prof   = array (data= NA,dim=c(nsimul,ntimes,ndcycle,nzdart))
               dart.list   = list  ( down.beam = dart.prof
                                   , down.diff = dart.prof
                                   , down.tot  = dart.prof
                                   , up.tot    = dart.prof
                                   , abs.lyr   = dart.prof
                                   , abs.wood  = dart.prof
                                   , abs.leaf  = dart.prof
                                   , eff.leaf  = dart.prof
                                   , cum.lyr   = dart.prof
                                   , cum.wood  = dart.prof
                                   , cum.leaf  = dart.prof
                                   , cum.eff   = dart.prof
                                   , atm.tot   = dart.atm
                                   , atm.beam  = dart.atm
                                   , atm.diff  = dart.atm
                                   )#end list
               model$dart  = list  ( par      = dart.list
                                   , nir      = dart.list
                                   , sol      = dart.list
                                   , tir      = dart.list
                                   , lai      = dart.ai
                                   , cumlai   = dart.ai
                                   )#end list
               #---------------------------------------------------------------------------#
            }#end if
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Load the patch and cohort information.                                  #
            #------------------------------------------------------------------------------#
            if (w == 1){
               nsites    = mymont$PYSI.N
               npatches  = mymont$SIPA.N
               ncohorts  = mymont$PACO.N
               a.ico     = mymont$PACO.ID
               z.ico     = a.ico + ncohorts - 1

               areapa    = mymont$AREA * rep(mymont$AREA.SI,times=npatches)
               areaco    = rep(areapa,times=ncohorts  )
               pareaco   = rep(areapa,times=ncohorts+1)

               isi       = rep(sequence(nsites),times=npatches)
               ipa       = sequence(npatches)
               ipaco     = rep(ipa,times=ncohorts  )
               ppaco     = rep(ipa,times=ncohorts+1)
               ico       = sequence(ncohorts)
               pco       = sequence(ncohorts+1)-1
               nlayer    = sum(ncohorts+1)


               #----- Initialise data regarding this polygon. -----------------------------#
               empty.ai    = array (data= NA,dim=c(ntimes        ,nlayer))
               empty.atm   = array (data= NA,dim=c(ntimes,ndcycle       ))
               empty.prof  = array (data= NA,dim=c(ntimes,ndcycle,nlayer))
               empty.list  = list  ( down.beam = empty.prof
                                   , down.diff = empty.prof
                                   , down.tot  = empty.prof
                                   , up.tot    = empty.prof
                                   , abs.lyr   = empty.prof
                                   , abs.tot   = empty.prof
                                   , abs.leaf  = empty.prof
                                   , abs.wood  = empty.prof
                                   , atm.tot   = empty.atm
                                   , atm.beam  = empty.atm
                                   , atm.diff  = empty.atm
                                   )#end list
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #     Simul is a template for each simulation.                              #
               #---------------------------------------------------------------------------#
               model$simul[[s]] = list( npatches = npatches
                                      , ncohorts = ncohorts
                                      , desert   = sum(ncohorts) == 0
                                      , nlayer   = nlayer
                                      , ipa      = ipa
                                      , ipaco    = ipaco
                                      , ppaco    = ppaco
                                      , ico      = ico
                                      , pco      = pco
                                      , areapa   = areapa
                                      , areaco   = areaco
                                      , pareaco  = pareaco
                                      , height   = empty.ai
                                      , lai      = empty.ai
                                      , wai      = empty.ai
                                      , tai      = empty.ai
                                      , par      = empty.list
                                      , nir      = empty.list
                                      , sol      = empty.list
                                      , tir      = empty.list
                                      )#end list
               #---------------------------------------------------------------------------#
            }#end if (w == 1)
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Retrieve some useful data from this run.                                 #
            #------------------------------------------------------------------------------#
            simnow   = model$simul[[s]]
            npatches = simnow$npatches
            ipa      = simnow$ipa
            ipaco    = simnow$ipaco
            ico      = simnow$ico
            ppaco    = simnow$ppa
            pco      = simnow$pco
            desert   = simnow$desert
            #------------------------------------------------------------------------------#



            #------ Variables to fill in constant boundary conditions. --------------------#
            zero        = rep(0               ,times=npatches)
            veritas     = rep(TRUE            ,times=npatches)
            zabove.pa   = rep(zdart.above     ,times=npatches)
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Make the cumulative LAI profile.                                         #
            #------------------------------------------------------------------------------#
            if (desert){
               lai.prof = zero
               wai.prof = zero
               tai.prof = zero
               hgt.prof = zabove.pa
               use.lyr  = veritas
            }else{
               lai.prof = append.patch( ipa   = ipa
                                      , ipaco = ipaco
                                      , xpa   = zero
                                      , xco   = mymont$MMEAN.LAI.CO
                                      )#end append.path
               wai.prof = append.patch( ipa   = ipa
                                      , ipaco = ipaco
                                      , xpa   = zero
                                      , xco   = mymont$WAI.CO
                                      )#end append.path
               tai.prof = lai.prof + wai.prof



               #---------------------------------------------------------------------------#
               #     Map the height to the closest DART level.                             #
               #---------------------------------------------------------------------------#
               ipft     = mymont$PFT
               ed.hgt   = mymont$HITE
               dart.hgt = ( zdart.above * (ed.hgt            - pft$hgt.min[ipft])
                                        / (pft$hgt.max[ipft] - pft$hgt.min[ipft]) )
               dart.hgt = mapply(FUN=closest,x=dart.hgt,MoreArgs=list(A=zdart[-1]))
               #---------------------------------------------------------------------------#


               #------ Height profile and flag for resolved cohort. -----------------------#
               hgt.prof = append.patch(ipa=ipa,ipaco=ipaco,xpa=zabove.pa,xco=dart.hgt)
               use.lyr  = append.patch( ipa   = ipa
                                      , ipaco = ipaco
                                      , xpa   = veritas
                                      , xco   = mymont$MMEAN.RLONG.L.CO != 0.
                                      )#end append.patch
               #---------------------------------------------------------------------------#

            }#end if
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Append layer characteristics for this time.                              #
            #------------------------------------------------------------------------------#
            simnow$lai[w,] = lai.prof
            simnow$wai[w,] = wai.prof
            simnow$tai[w,] = tai.prof
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Loop over radiation bands.                                               #
            #------------------------------------------------------------------------------#
            for (b in sequence(nband)){
               bnd   = band$key [b]
               ed2   = band$ed2 [b]
               bdesc = band$desc[b]


               #----- Get the components. -------------------------------------------------#
               down.beam.top = mymont[[paste("QMEAN.",ed2,".BEAM.PA",sep="")]]
               down.diff.top = mymont[[paste("QMEAN.",ed2,".DIFF.PA",sep="")]]
               up.tot.top    = mymont[[paste("QMEAN.",ed2,"UP.PA"   ,sep="")]]
               down.beam.coh = mymont[[paste("QMEAN.",ed2,".BEAM.CO",sep="")]]
               down.diff.coh = mymont[[paste("QMEAN.",ed2,".DIFF.CO",sep="")]]
               up.tot.coh    = mymont[[paste("QMEAN.",ed2,"UP.CO"   ,sep="")]]
               abs.leaf.coh  = mymont[[paste("QMEAN.",ed2,".L.CO"   ,sep="")]]

               atm.beam      = mymont[[paste("QMEAN.ATM.",ed2,".BEAM.PY",sep="")]]
               atm.diff      = mymont[[paste("QMEAN.ATM.",ed2,".DIFF.PY",sep="")]]
               atm.tot       = atm.beam + atm.diff
               #---------------------------------------------------------------------------#




               #---------------------------------------------------------------------------#
               #    Append data by time of day.                                            #
               #---------------------------------------------------------------------------#
               for (h in sequence(ndcycle)){

                  cat("       - ",bdesc,": ",sprintf("%2.2i",h)," UTC","\n",sep="")
                  
                  #------------------------------------------------------------------------#
                  #     Add top boundary condition to each cohort.                         #
                  #------------------------------------------------------------------------#
                  if (desert){
                     #----- Desert run, no cohorts. ---------------------------------------#
                     down.beam = down.beam.top[,h]
                     down.diff = down.diff.top[,h]
                     up.tot    = up.top       [,h]
                     abs.leaf  = zero
                     
                     #---------------------------------------------------------------------#
                  }else{
                     #----- Append the patch-level top boundary condition. ----------------#
                     down.beam = append.patch( ipa   = ipa
                                             , ipaco = ipaco
                                             , xpa   = down.beam.top[,h]
                                             , xco   = down.beam.coh[,h]
                                             )#end append.patch
                     down.diff = append.patch( ipa   = ipa
                                             , ipaco = ipaco
                                             , xpa   = down.diff.top[,h]
                                             , xco   = down.diff.coh[,h]
                                             )#end append.patch
                     up.tot    = append.patch( ipa   = ipa
                                             , ipaco = ipaco
                                             , xpa   = up.tot.top   [,h]
                                             , xco   = up.tot.coh   [,h]
                                             )#end append.patch
                     abs.leaf  = append.patch( ipa   = ipa
                                             , ipaco = ipaco
                                             , xpa   = zero
                                             , xco   = abs.leaf.coh [,h]
                                             )#end append.patch
                     #---------------------------------------------------------------------#



                     #----- Fill in cohorts that haven't been resolved. -------------------#
                     down.beam = fill.unresolved(ipaco=ppaco,x=down.beam,use=use.lyr)
                     down.diff = fill.unresolved(ipaco=ppaco,x=down.diff,use=use.lyr)
                     up.tot    = fill.unresolved(ipaco=ppaco,x=up.tot   ,use=use.lyr)
                     #---------------------------------------------------------------------#
                  }#end if
                  #------------------------------------------------------------------------#



                  #----- Find the derived quantities. -------------------------------------#
                  down.tot  = down.beam + down.diff
                  abs.lyr   = layer.absorption(ipaco=ppaco,down=down.tot,up=up.tot)
                  abs.wood  = abs.lyr - abs.leaf
                  cum.lyr   = unlist(tapply(X=abs.lyr ,INDEX=ppaco,FUN=cumsum))
                  cum.leaf  = unlist(tapply(X=abs.leaf,INDEX=ppaco,FUN=cumsum))
                  cum.wood  = unlist(tapply(X=abs.wood,INDEX=ppaco,FUN=cumsum))
                  #------------------------------------------------------------------------#



                  #----- Update profile. --------------------------------------------------#
                  simnow[[bnd]]$down.beam [w,h,] = down.beam
                  simnow[[bnd]]$down.diff [w,h,] = down.diff
                  simnow[[bnd]]$down.tot  [w,h,] = down.tot
                  simnow[[bnd]]$up.tot    [w,h,] = up.tot
                  simnow[[bnd]]$abs.lyr   [w,h,] = abs.lyr
                  simnow[[bnd]]$abs.tot   [w,h,] = cum.lyr
                  simnow[[bnd]]$abs.leaf  [w,h,] = cum.leaf
                  simnow[[bnd]]$abs.wood  [w,h,] = cum.wood
                  #------------------------------------------------------------------------#


                  #------------------------------------------------------------------------#
                  #    Eliminate duplicates within each patch.  We keep only the last      #
                  # light level for each group of duplicates, the sum of absorption and    #
                  # LAI, TRUE if any of the use.lyr is TRUE, and the commonest height and  #
                  # ppaco.                                                                 #
                  #------------------------------------------------------------------------#
                  #------ Split data. -----------------------------------------------------#
                  unq.index     = 10000*ppaco - hgt.prof
                  unq.ppaco     = tapply( X     = ppaco
                                        , INDEX = unq.index
                                        , FUN   = commonest
                                        , na.rm = TRUE
                                        )#end tapply
                  unq.hgt.prof  = tapply( X     = hgt.prof
                                        , INDEX = unq.index
                                        , FUN   = commonest
                                        , na.rm = TRUE
                                        )#end tapply
                  unq.down.beam = tapply( X     = down.beam
                                        , INDEX = unq.index
                                        , FUN   = get.last
                                        , na.rm = TRUE
                                        )#end tapply
                  unq.down.diff = tapply( X     = down.diff
                                        , INDEX = unq.index
                                        , FUN   = get.last
                                        , na.rm = TRUE
                                        )#end tapply
                  unq.up.tot    = tapply( X     = up.tot
                                        , INDEX = unq.index
                                        , FUN   = get.last
                                        , na.rm = TRUE
                                        )#end tapply
                  unq.use.lyr   = tapply( X     = use.lyr
                                        , INDEX = unq.index
                                        , FUN   = any
                                        , na.rm = TRUE
                                        )#end tapply
                  unq.abs.leaf  = tapply( X     = abs.leaf
                                        , INDEX = unq.index
                                        , FUN   = sum
                                        , na.rm = TRUE
                                        )#end tapply
                  unq.lai.prof  = tapply( X     = lai.prof
                                        , INDEX = unq.index
                                        , FUN   = sum
                                        , na.rm = TRUE
                                        )#end tapply
                  #------ Re-combine data. ------------------------------------------------#
                  names(unq.ppaco    ) = NULL
                  names(unq.hgt.prof ) = NULL
                  names(unq.down.beam) = NULL
                  names(unq.down.diff) = NULL
                  names(unq.up.tot   ) = NULL
                  names(unq.use.lyr  ) = NULL
                  names(unq.abs.leaf ) = NULL
                  names(unq.lai.prof ) = NULL
                  #------------------------------------------------------------------------#



                  #----- Initialise matrix profiles. --------------------------------------#
                  dart.down.beam = matrix(data=     0,nrow=nzdart,ncol=npatches)
                  dart.down.diff = matrix(data=     0,nrow=nzdart,ncol=npatches)
                  dart.up.tot    = matrix(data=     0,nrow=nzdart,ncol=npatches)
                  dart.use.lyr   = matrix(data= FALSE,nrow=nzdart,ncol=npatches)
                  dart.abs.leaf  = matrix(data=     0,nrow=nzdart,ncol=npatches)
                  dart.lai       = matrix(data=     0,nrow=nzdart,ncol=npatches)
                  dart.area      = matrix(data=areapa,nrow=nzdart,ncol=npatches,byrow=TRUE)
                  #------------------------------------------------------------------------#


                  #----- Create height/patch index. ---------------------------------------#
                  idx                  = cbind(match(unq.hgt.prof,zdart),unq.ppaco)
                  dart.down.beam [idx] = unq.down.beam
                  dart.down.diff [idx] = unq.down.diff
                  dart.up.tot    [idx] = unq.up.tot
                  dart.use.lyr   [idx] = unq.use.lyr
                  dart.abs.leaf  [idx] = unq.abs.leaf
                  dart.lai       [idx] = unq.lai.prof
                  #------------------------------------------------------------------------#



                  #----- Convert matrices to data frame (so mapply works). ----------------#
                  dart.down.beam = as.data.frame(dart.down.beam)
                  dart.down.diff = as.data.frame(dart.down.diff)
                  dart.up.tot    = as.data.frame(dart.up.tot   )
                  dart.use.lyr   = as.data.frame(dart.use.lyr  )
                  dart.abs.leaf  = as.data.frame(dart.abs.leaf )
                  dart.lai       = as.data.frame(dart.lai      )
                  dart.area      = as.data.frame(dart.area     )
                  #------------------------------------------------------------------------#




                  #----- Find effective light absorption. ---------------------------------#
                  if (bnd %in% "par"){
                     dart.eff.leaf  = ( 0 * dart.abs.leaf
                                      + pmin(dart.crit.par*dart.lai,dart.abs.leaf)
                                      )#end dart.eff.lyr
                  }else{
                     dart.eff.leaf  = dart.abs.leaf
                  }#end if (bnd %in% "par")
                  #------------------------------------------------------------------------#




                  #----- Convert matrices to data frame (so mapply works). ----------------#
                  dart.down.beam = mapply( FUN      = fill.unresolved
                                         , x        = dart.down.beam
                                         , use      = dart.use.lyr
                                         , MoreArgs = list(ipaco=rep(1,nzdart))
                                         )#end mapply
                  dart.down.diff = mapply( FUN      = fill.unresolved
                                         , x        = dart.down.diff
                                         , use      = dart.use.lyr
                                         , MoreArgs = list(ipaco=rep(1,nzdart))
                                         )#end mapply
                  dart.up.tot    = mapply( FUN      = fill.unresolved
                                         , x        = dart.up.tot
                                         , use      = dart.use.lyr
                                         , MoreArgs = list(ipaco=rep(1,nzdart))
                                         )#end mapply
                  dart.down.beam = as.data.frame(dart.down.beam)
                  dart.down.diff = as.data.frame(dart.down.diff)
                  dart.up.tot    = as.data.frame(dart.up.tot   )
                  #------------------------------------------------------------------------#




                  #----- Find the derived quantities. -------------------------------------#
                  dart.down.tot  = dart.down.beam + dart.down.diff
                  dart.abs.lyr   = mapply( FUN      = layer.absorption
                                         , down     = dart.down.tot
                                         , up       = dart.up.tot
                                         , MoreArgs = list(ipaco=rep(1,nzdart))
                                         )#end mapply
                  dart.abs.lyr   = as.data.frame(dart.abs.lyr)
                  dart.abs.wood  = dart.abs.lyr - dart.abs.leaf
                  #------------------------------------------------------------------------#


                  #------------------------------------------------------------------------#
                  #     Find the average light level and absorption.                       #
                  #------------------------------------------------------------------------#
                  dart.lai       = rowSums(dart.lai      *dart.area)
                  dart.down.beam = rowSums(dart.down.beam*dart.area)
                  dart.down.diff = rowSums(dart.down.diff*dart.area)
                  dart.down.tot  = rowSums(dart.down.tot *dart.area)
                  dart.up.tot    = rowSums(dart.up.tot   *dart.area)
                  dart.abs.lyr   = rowSums(dart.abs.lyr  *dart.area)
                  dart.abs.wood  = rowSums(dart.abs.wood *dart.area)
                  dart.abs.leaf  = rowSums(dart.abs.leaf *dart.area)
                  dart.eff.leaf  = rowSums(dart.eff.leaf *dart.area)
                  dart.cumlai    = cumsum(dart.lai     )
                  dart.cum.lyr   = cumsum(dart.abs.lyr )
                  dart.cum.wood  = cumsum(dart.abs.wood)
                  dart.cum.leaf  = cumsum(dart.abs.leaf)
                  dart.cum.eff   = cumsum(dart.eff.leaf)
                  #------------------------------------------------------------------------#


                  #----- Update profile. --------------------------------------------------#
                  model$dart[[bnd]]$down.beam [s,w,h,] = dart.down.beam
                  model$dart[[bnd]]$down.diff [s,w,h,] = dart.down.diff
                  model$dart[[bnd]]$down.tot  [s,w,h,] = dart.down.tot
                  model$dart[[bnd]]$up.tot    [s,w,h,] = dart.up.tot
                  model$dart[[bnd]]$abs.lyr   [s,w,h,] = dart.abs.lyr
                  model$dart[[bnd]]$abs.wood  [s,w,h,] = dart.abs.wood
                  model$dart[[bnd]]$abs.leaf  [s,w,h,] = dart.abs.leaf
                  model$dart[[bnd]]$eff.leaf  [s,w,h,] = dart.eff.leaf
                  model$dart[[bnd]]$cum.lyr   [s,w,h,] = dart.cum.lyr
                  model$dart[[bnd]]$cum.wood  [s,w,h,] = dart.cum.wood
                  model$dart[[bnd]]$cum.leaf  [s,w,h,] = dart.cum.leaf
                  model$dart[[bnd]]$cum.eff   [s,w,h,] = dart.cum.eff
                  if (h == 1 && b == 1){
                     model$dart$lai   [s,w,] = dart.lai
                     model$dart$cumlai[s,w,] = dart.cumlai
                  }#end if
                  #------------------------------------------------------------------------#



                  #----- Free memory. -----------------------------------------------------#
                  rm(unq.index,unq.ppaco,unq.hgt.prof,unq.down.beam,uq.down.diff)
                  rm(unq.up.tot,unq.use.lyr,unq.abs.leaf,unq.lai.prof)
                  rm(dart.lai,dart.down.beam,dart.down.diff,dart.down.tot,dart.up.tot)
                  rm(dart.abs.lyr,dart.abs.wood,dart.abs.leaf,dart.eff.leaf)
                  rm(dart.cumlai,dart.cum.lyr,dart.cum.wood,dart.cum.leaf,dart.cum.eff)
                  rm(dart.area,idx,down.tot,abs.lyr,abs.wood,cum.lyr,cum.leaf,cum.wood)
                  rm(down.beam,down.diff,up.tot,abs.leaf)
                  #------------------------------------------------------------------------#
               }#end for (h in sequence(ndcycle))
               #---------------------------------------------------------------------------#



               #----- Update ToC radiation. -----------------------------------------------#
               simnow[[bnd]]$atm.beam    [  w,] = atm.beam
               simnow[[bnd]]$atm.diff    [  w,] = atm.diff
               simnow[[bnd]]$atm.tot     [  w,] = atm.tot
               model$dart[[bnd]]$atm.beam[s,w,] = atm.beam
               model$dart[[bnd]]$atm.diff[s,w,] = atm.diff
               model$dart[[bnd]]$atm.tot [s,w,] = atm.tot
               rm(atm.beam,atm.diff,atm.tot,down.beam.top,down.diff.top,up.tot.top)
               rm(down.beam.coh,down.diff.coh,up.tot.coh,abs.leaf.coh)
               #---------------------------------------------------------------------------#
            }#end for (b in sequence(nband))
            #------------------------------------------------------------------------------#
         }#end for (w in sequence(nwhen))
         #---------------------------------------------------------------------------------#

         #----- Reload data and copy to the general list. ---------------------------------#
         rdata.iata     = file.path(rdata.path,paste(iata,rdata.suffix,sep="_"))
         rdata.finished = file.path(rdata.path
                                   ,paste("loaded",iata,finished.suffix,sep="_"))
         cat(" + Saving data to ",basename(rdata.iata),"...","\n",sep="")
         model$sresume  = s+1
         eft[[iata]]    = model
         dummy          = save(model,file=rdata.iata)
         dummy          = write(x=paste(Sys.time()),file=rdata.finished)
         rm(mymont,simnow,lai.prof,wai.prof,tai.prof,ed.hgt,dart.hgt,hgt.prof,use.lyr)
         #---------------------------------------------------------------------------------#
      }#end for (s in sequence(nsimul))
      #------------------------------------------------------------------------------------#
   }#end for (p in sequence(nsites))
   #=======================================================================================#
   #=======================================================================================#
}#end if ( ( (! file.exists(summary.fullname)) && plot.patch ) || plot.site )
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     Site-specific plots.                                                                 #
#------------------------------------------------------------------------------------------#
cat(" + Entering plot.site...","\n")
if (plot.site){
   #---------------------------------------------------------------------------------------#
   #    Site loop.                                                                         #
   #---------------------------------------------------------------------------------------#
   for (p in sequence(nsites)){
      #----- Retrieve site information. ---------------------------------------------------#
      iata     = sites$iata[p]
      model    = eft[[iata]]
      yeara    = model$yeara
      yearz    = model$yearz
      lon      = model$lon
      lat      = model$lat
      desc     = model$longname
      dart     = model$dart
      cat(" - Generating profiles for ",desc,"...","\n")
      #------------------------------------------------------------------------------------#


      #----- Make output path for this place. ---------------------------------------------#
      sitepath = file.path(outroot,model$short)
      if (! file.exists(sitepath)) dir.create(sitepath)
      #------------------------------------------------------------------------------------#



      #----- Choose the patches to plot based on the LAI. ---------------------------------#
      cat("   > LAI by month...","\n")
      if (plot.month){
         loop.mon = sequence(12)
      }else{
         loop.mon = numeric(0)
      }#end if (plot.month)
      #------------------------------------------------------------------------------------#



      #----- Make output path for this place. ---------------------------------------------#
      profilepath  = file.path(sitepath,"profile" )
      dielpath     = file.path(sitepath,"diel"    )
      reldielpath  = file.path(sitepath,"rel_diel")
      laipath      = file.path(sitepath,"lai"     )
      if (! file.exists(profilepath)) dir.create(profilepath)
      if (! file.exists(dielpath   )) dir.create(dielpath   )
      if (! file.exists(reldielpath)) dir.create(reldielpath)
      if (! file.exists(laipath    )) dir.create(laipath    )
      #------------------------------------------------------------------------------------#




      #----- Loop over months and hours. --------------------------------------------------#
      for (m in loop.mon){
         mm   = sprintf("%2.2i",m)
         ms   = model$month %in% m
         cat("     > ",month.name[m],"\n")


         #---------------------------------------------------------------------------------#
         #     Load LAI/WAI/TAI profiles.                                                  #
         #---------------------------------------------------------------------------------#
         lai     = dart$lai
         cum.lai = dart$cumlai
         #---------------------------------------------------------------------------------#


         #----- Fix the limits for the Y axis. --------------------------------------------#
         lai.lim = range(-cum.lai,finite=TRUE)
         lai.at   = pretty(lai.lim)
         lai.lab  = -lai.at
         #---------------------------------------------------------------------------------#


         zen.mon = colMeans(model$zen[ms,,drop=FALSE])


         #---------------------------------------------------------------------------------#
         #     Generate table with indices and hours.                                      #
         #---------------------------------------------------------------------------------#
         for (b in sequence(nband)){
            bnd      = band$key [b]
            bnd.desc = band$desc[b]
            cat("       * ",bnd.desc,"\n")

            mult     = if (bnd %in% c("par")){Watts.2.Ein * 1.e6}else{1.0}
            prof     = dart[[bnd]]



            #------------------------------------------------------------------------------#
            #     Loop over each simulation.                                               #
            #------------------------------------------------------------------------------#
            for (ss in sequence(nsimul)){
               s.suff = simul$suff[ss]
               s.desc = simul$desc[ss]
               cat("         - ",s.desc,"\n")


               #----- Find the averages for the hour and month. ---------------------------#
               down.beam = apply( X      = prof$down.beam[ss,ms,,,drop=FALSE] * mult
                                , MARGIN = c(1,3,4)
                                , FUN    = mean
                                )#end apply
               down.diff = apply( X      = prof$down.diff[ss,ms,,,drop=FALSE] * mult
                                , MARGIN = c(1,3,4)
                                , FUN    = mean
                                )#end apply
               down.tot  = apply( X      = prof$down.tot [ss,ms,,,drop=FALSE] * mult
                                , MARGIN = c(1,3,4)
                                , FUN    = mean
                                )#end apply
               up.tot    = apply( X      = prof$up.tot   [ss,ms,,,drop=FALSE] * mult
                                , MARGIN = c(1,3,4)
                                , FUN    = mean
                                )#end apply
               abs.lyr   = apply( X      = prof$abs.lyr  [ss,ms,,,drop=FALSE] * mult
                                , MARGIN = c(1,3,4)
                                , FUN    = mean
                                )#end apply
               abs.wood  = apply( X      = prof$abs.wood [ss,ms,,,drop=FALSE] * mult
                                , MARGIN = c(1,3,4)
                                , FUN    = mean
                                )#end apply
               abs.leaf  = apply( X      = prof$abs.leaf [ss,ms,,,drop=FALSE] * mult
                                , MARGIN = c(1,3,4)
                                , FUN    = mean
                                )#end apply
               eff.leaf  = apply( X      = prof$eff.leaf [ss,ms,,,drop=FALSE] * mult
                                , MARGIN = c(1,3,4)
                                , FUN    = mean
                                )#end apply
               cum.lyr   = apply( X      = prof$cum.lyr  [ss,ms,,,drop=FALSE] * mult
                                , MARGIN = c(1,3,4)
                                , FUN    = mean
                                )#end apply
               cum.wood  = apply( X      = prof$cum.wood [ss,ms,,,drop=FALSE] * mult
                                , MARGIN = c(1,3,4)
                                , FUN    = mean
                                )#end apply
               cum.leaf  = apply( X      = prof$cum.leaf [ss,ms,,,drop=FALSE] * mult
                                , MARGIN = c(1,3,4)
                                , FUN    = mean
                                )#end apply
               cum.eff   = apply( X      = prof$cum.eff  [ss,ms,,,drop=FALSE] * mult
                                , MARGIN = c(1,3,4)
                                , FUN    = mean
                                )#end apply
               down.beam = aperm(a=down.beam,perm=c(3,2,1))[,,1]
               down.diff = aperm(a=down.diff,perm=c(3,2,1))[,,1]
               down.tot  = aperm(a=down.tot ,perm=c(3,2,1))[,,1]
               up.tot    = aperm(a=up.tot   ,perm=c(3,2,1))[,,1]
               abs.lyr   = aperm(a=abs.lyr  ,perm=c(3,2,1))[,,1]
               abs.wood  = aperm(a=abs.wood ,perm=c(3,2,1))[,,1]
               abs.leaf  = aperm(a=abs.leaf ,perm=c(3,2,1))[,,1]
               eff.leaf  = aperm(a=eff.leaf ,perm=c(3,2,1))[,,1]
               cum.lyr   = aperm(a=cum.lyr  ,perm=c(3,2,1))[,,1]
               cum.wood  = aperm(a=cum.wood ,perm=c(3,2,1))[,,1]
               cum.leaf  = aperm(a=cum.leaf ,perm=c(3,2,1))[,,1]
               cum.eff   = aperm(a=cum.eff  ,perm=c(3,2,1))[,,1]
               #---------------------------------------------------------------------------#



               #----- Create data frame that will be written to a csv file. ---------------#
               out.df = data.frame( hour      = model$hour[col(abs.leaf)]
                                  , zen.sun   = zen.mon   [col(abs.leaf)]
                                  , height    = rep(model$zdart    ,times=nhours)
                                  , lai.lyr   = rep(lai    [ss,ms,],times=nhours)
                                  , cum.lai   = rep(cum.lai[ss,ms,],times=nhours)
                                  , down.beam = c(down.beam)
                                  , down.diff = c(down.diff)
                                  , down.tot  = c(down.tot)
                                  , up.tot    = c(up.tot)
                                  , abs.lyr   = c(abs.lyr)
                                  , abs.wood  = c(abs.wood)
                                  , abs.leaf  = c(abs.leaf)
                                  , eff.leaf  = c(eff.leaf)
                                  , cum.lyr   = c(cum.lyr)
                                  , cum.wood  = c(cum.wood)
                                  , cum.leaf  = c(cum.leaf)
                                  , cum.eff   = c(cum.eff )
                                  )#end data.frame
               #---------------------------------------------------------------------------#

               #---------------------------------------------------------------------------#
               #     Save data frame.                                                      #
               #---------------------------------------------------------------------------#
               basefile  = paste0(iata,"_",bnd,"_",mm,"-",tolower(month.abb[m])
                                 ,"_",s.suff,".csv")
               asciifile = file.path(ascroot,basefile)
               dummy     = write.csv(x=format(out.df),file=asciifile,quote=FALSE
                                    ,na="NaN",row.names=FALSE)
               #---------------------------------------------------------------------------#
            }#end for (ss in sequence(nsimul))
            #------------------------------------------------------------------------------#

         }#end for (b in sequence(nband))
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Generate profiles for each hour of the day.                                 #
         #---------------------------------------------------------------------------------#
         cat("       ~ Average profile by hour of the day...","\n")
         for (h in prof.hours){
            hh      = sprintf("%2.2i",h)
            hs      = model$hour %in% h
            zen.bar = mean(model$zen[ms,hs])



            #------ Create path for this hour. --------------------------------------------#
            profhrpath  = file.path(profilepath,paste(hh,"utc",sep=""))
            if (! file.exists(profhrpath)) dir.create(profhrpath)
            #------------------------------------------------------------------------------#



            #----- Band loop. -------------------------------------------------------------#
            for (b in sequence(nband)){
               bnd      = band$key [b]
               bnd.desc = band$desc[b]

               mult     = if (bnd %in% c("par")){Watts.2.Ein * 1.e6}else{1.0}
               prof     = dart[[bnd]]

               #----- Find the averages for the hour and month. ---------------------------#
               down.beam = apply( X      = prof$down.beam[,ms,hs,,drop=FALSE] * mult
                                , MARGIN = c(1,4)
                                , FUN    = mean
                                )#end apply
               down.diff = apply( X      = prof$down.diff[,ms,hs,,drop=FALSE] * mult
                                , MARGIN = c(1,4)
                                , FUN    = mean
                                )#end apply
               down.tot  = apply( X      = prof$down.tot [,ms,hs,,drop=FALSE] * mult
                                , MARGIN = c(1,4)
                                , FUN    = mean
                                )#end apply
               up.tot    = apply( X      = prof$up.tot   [,ms,hs,,drop=FALSE] * mult
                                , MARGIN = c(1,4)
                                , FUN    = mean
                                )#end apply
               abs.lyr   = apply( X      = prof$abs.lyr  [,ms,hs,,drop=FALSE] * mult
                                , MARGIN = c(1,4)
                                , FUN    = mean
                                )#end apply
               abs.wood  = apply( X      = prof$abs.wood [,ms,hs,,drop=FALSE] * mult
                                , MARGIN = c(1,4)
                                , FUN    = mean
                                )#end apply
               abs.leaf  = apply( X      = prof$abs.leaf [,ms,hs,,drop=FALSE] * mult
                                , MARGIN = c(1,4)
                                , FUN    = mean
                                )#end apply
               eff.leaf  = apply( X      = prof$eff.leaf [,ms,hs,,drop=FALSE] * mult
                                , MARGIN = c(1,4)
                                , FUN    = mean
                                )#end apply
               cum.lyr   = apply( X      = prof$cum.lyr  [,ms,hs,,drop=FALSE] * mult
                                , MARGIN = c(1,4)
                                , FUN    = mean
                                )#end apply
               cum.wood  = apply( X      = prof$cum.wood [,ms,hs,,drop=FALSE] * mult
                                , MARGIN = c(1,4)
                                , FUN    = mean
                                )#end apply
               cum.leaf  = apply( X      = prof$cum.leaf [,ms,hs,,drop=FALSE] * mult
                                , MARGIN = c(1,4)
                                , FUN    = mean
                                )#end apply
               cum.eff   = apply( X      = prof$cum.eff  [,ms,hs,,drop=FALSE] * mult
                                , MARGIN = c(1,4)
                                , FUN    = mean
                                )#end apply
               #---------------------------------------------------------------------------#




               #---------------------------------------------------------------------------#
               #     Shift the results a bit so we can see when two curves are on top      #
               # of each other.                                                            #
               #---------------------------------------------------------------------------#
               if ( ! bnd %in% "tir"){
                  down.beam = apply(X = down.beam,MARGIN=2,FUN="*",simul$mult.beam.down)
                  down.diff = apply(X = down.diff,MARGIN=2,FUN="*",simul$mult.diff.down)
                  up.tot    = apply(X = up.tot   ,MARGIN=2,FUN="*",simul$mult.up       )
               }#end if
               #---------------------------------------------------------------------------#





               #---------------------------------------------------------------------------#
               #     Skip plot if it is solar radiation and it is nighttime.               #
               #---------------------------------------------------------------------------#
               if (max(down.tot,na.rm=TRUE) > 1.){


                  #----- Retrieve the data. -----------------------------------------------#
                  if (! bnd %in% "tir"){
                     down.lim = pretty.xylim(u=c(down.beam,down.diff))
                  }else{
                     down.lim = pretty.xylim(u=down.diff)
                  }#end if
                  up.lim   = pretty.xylim(u=up.tot )
                  abs.lim  = pretty.xylim(u=c(0,cum.eff,cum.leaf))
                  #------------------------------------------------------------------------#


                  for (o in sequence(nout)){
                     #----- Open file. ----------------------------------------------------#
                     fichier = file.path(profhrpath
                                        ,paste0("profile-",iata,"-",bnd,"-"
                                               ,mm,"-",hh,"utc",".",outform[o]))
                     if (bnd %in% "tir"){nsize=wsize}else{nsize=size}
                     if       (outform[o] == "x11"){
                       X11(width=nsize$width,height=nsize$height,pointsize=ptsz)
                     }else if (outform[o] == "png"){
                       png(filename=fichier,width=nsize$width*depth
                          ,height=nsize$height*depth,pointsize=ptsz,res=depth)
                     }else if (outform[o] == "eps"){
                       postscript(file=fichier,width=nsize$width,height=nsize$height
                                 ,pointsize=ptsz,paper=nsize$paper)
                     }else if (outform[o] == "pdf"){
                       pdf(file=fichier,onefile=FALSE,width=nsize$width
                          ,height=nsize$height,pointsize=ptsz,paper=nsize$paper)
                     }#end if
                     #---------------------------------------------------------------------#


                     #----- Set up the plotting window. -----------------------------------#
                     par(par.user)
                     par(oma=c(0,1,4,0))
                     if (bnd %in% "tir"){
                        layout(mat=rbind(c(2,3,4),c(1,1,1)),heights=c(6,1))
                     }else{
                        layout(mat=rbind(c(2,4),c(3,5),c(1,1)),heights=c(3,3,1))
                     }#end if
                     #---------------------------------------------------------------------#




                     #----- First plot: the legend. ---------------------------------------#
                     par(mar=c(0.1,4.1,0.1,2.1))
                     plot.new()
                     plot.window(xlim=c(0,1),ylim=c(0,1))
                     legend( x      = "bottom"
                           , inset  = 0.0
                           , legend = simul$desc
                           , col    = simul$colour
                           , lty    = "solid"
                           , lwd    = 2.5
                           , cex    = 0.9 * cex.ptsz
                           , ncol   = nsimul
                           , xpd    = TRUE
                           )#end legend
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #      Direct radiation is not plotted for TIR.                       #
                     #---------------------------------------------------------------------#
                     if (! bnd %in% "tir"){
                        #------ Second plot: Downward direct. -----------------------------#
                        par(mar=c(3.1,3.1,3.1,1.1))
                        plot.new()
                        plot.window(xlim=down.lim,ylim=lai.lim)
                        axis(side=1,cex.axis=0.8,padj=-0.8)
                        axis(side=2,at=lai.at,labels=lai.lab,las=1)
                        grid(col=grid.colour,lwd=1.0,lty="solid")
                        box()
                        title(main="Downward - Direct",font.main=1,line=1)
                        for (s in sequence(nsimul)){
                           lines(x=down.beam[s,],y=-cum.lai[s,ms,]
                                ,col=simul$colour[s],lwd=2.5,lty="solid")
                        }#end for
                        #------------------------------------------------------------------#
                     }#end if
                     #---------------------------------------------------------------------#




                     #------ Third plot: Downward diffuse. --------------------------------#
                     par(mar=c(4.1,3.1,2.1,1.1))
                     plot.new()
                     plot.window(xlim=down.lim,ylim=lai.lim)
                     axis(side=1,cex.axis=0.8,padj=-0.8)
                     axis(side=2,at=lai.at,labels=lai.lab,las=1)
                     grid(col=grid.colour,lwd=1.0,lty="solid")
                     box()
                     title(main="Downward - Diffuse",font.main=1,line=1)
                     for (s in sequence(nsimul)){
                        lines(x=down.diff[s,],y=-cum.lai[s,ms,]
                             ,col=simul$colour[s],lwd=2.5,lty="solid")
                     }#end for
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Fourth plot: Upward radiation.  Margins depend on the band.     #
                     #---------------------------------------------------------------------#
                     if (bnd %in% "tir"){
                        par(mar=c(4.1,3.1,2.1,1.1))
                     }else{
                        par(mar=c(3.1,2.1,3.1,2.1))
                     }#end if
                     plot.new()
                     plot.window(xlim=up.lim,ylim=lai.lim)
                     axis(side=1,cex.axis=0.8,padj=-0.8)
                     axis(side=2,at=lai.at,labels=lai.lab,las=1)
                     grid(col=grid.colour,lwd=1.0,lty="solid")
                     box()
                     title(main="Upward",font.main=1,line=1)
                     for (s in sequence(nsimul)){
                        lines(x=up.tot[s,],y=-cum.lai[s,ms,]
                             ,col=simul$colour[s],lwd=2.5,lty="solid")
                     }#end for
                     #---------------------------------------------------------------------#



                     #------ Fifth plot: Absorbed radiation. ------------------------------#
                     par(mar=c(4.1,2.1,2.1,2.1))
                     plot.new()
                     plot.window(xlim=abs.lim,ylim=lai.lim)
                     axis(side=1,cex.axis=0.8,padj=-0.8)
                     axis(side=2,at=lai.at,labels=lai.lab,las=1)
                     grid(col=grid.colour,lwd=1.0,lty="solid")
                     box()
                     title(main="Cumulative absorption (dashed = effective)"
                          ,font.main=1,line=1)
                     for (s in sequence(nsimul)){
                        lines(x=cum.eff [s,],y=-cum.lai[s,ms,]
                             ,col=simul$colour[s],lwd=2.0,lty="dotdash")
                        lines(x=cum.leaf[s,],y=-cum.lai[s,ms,]
                             ,col=simul$colour[s],lwd=2.5,lty="solid")
                     }#end for
                     #---------------------------------------------------------------------#



                     #---------------------------------------------------------------------#
                     #     Plot the general annotation.                                    #
                     #---------------------------------------------------------------------#
                     lex     = desc.unit( desc = paste(band$desc[b],"Irradiance")
                                        , unit = untab[[band$unit[b]]]
                                        )#end desc.unit
                     ley     = desc.unit(desc="Cumulative LAI",unit=untab$m2lom2)
                     whenlab = paste(month.name[m],"-",hh,"UTC")
                     letitre = paste(desc," - Canopy profiles","  -  ",whenlab,sep="")
                     laimon  = max(cum.lai[,ms,],na.rm=TRUE)
                     lelai   = desc.unit( desc = paste("Average LAI"
                                                      ,sprintf("%6.2f",laimon)
                                                      )#end paste
                                        , unit = untab$m2lom2
                                        , bracket = FALSE
                                        )#end desc.unit
                     lezen   = desc.unit( desc    = paste("Sun zenith angle"
                                                         ,sprintf("%5.2f",zen.bar)
                                                         )#end paste
                                        , unit    = untab$deg
                                        , bracket = FALSE
                                        )#end desc.unit
                     lesub   = parse(text = paste(lelai,lezen,sep=" - "))
                     gtitle( main      = letitre
                           , sub       = lesub
                           , xlab      = lex
                           , ylab      = ley
                           , line.xlab = 3.2
                           , line.ylab = 2.6
                           , line.sub  = 4.2
                           , off.sub   = 2/21
                           , cex.main  = 1.0*cex.ptsz
                           , cex.xlab  = 0.9*cex.ptsz
                           , cex.ylab  = 0.9*cex.ptsz
                           , cex.sub   = 0.7*cex.ptsz
                           )#end gtitle
                     #---------------------------------------------------------------------#



                     #----- Close the device. ---------------------------------------------#
                     if (outform[o] == "x11"){
                        locator(n=1)
                        dev.off()
                     }else{
                        dev.off()
                     }#end if
                     clean.tmp()
                     #---------------------------------------------------------------------#
                  }#end for (o in sequence(nout))
                  #------------------------------------------------------------------------#
               }#end if (max(tot.down) > 1.)
               #---------------------------------------------------------------------------#
            }#end for (b in sequence(nband))
            #------------------------------------------------------------------------------#
         }#end for (h in prof.hours)
         #---------------------------------------------------------------------------------#
      }#end for (m in sequence(12))
      #------------------------------------------------------------------------------------#
   }#end for (p in sequence(nsites))
#------------------------------------------------------------------------------------------#
}#end if (plot.site)
#==========================================================================================#
#==========================================================================================#
