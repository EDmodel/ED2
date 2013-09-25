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
srcdir  = "/n/home00/mlongo/util/Rsc"            #   Script directory
outroot = file.path(here,"radprof_comp")         #   Output directory

#------ Plot options. ---------------------------------------------------------------------#
outform        = c("pdf")              # Formats for output file.  Supported formats are:
                                       #   - "X11" - for printing on screen
                                       #   - "eps" - for postscript printing
                                       #   - "png" - for PNG printing
                                       #   - "pdf" - for PDF printing
depth          = 96                    # PNG resolution, in pixels per inch
paper          = "letter"              # Paper size, to define the plot shape
wpaper         = "legal"               # Wide paper size, to define the plot shape
ptsz           = 18                    # Font size.
ibackground    = 0                     # Make figures compatible to which background?
                                       # 0 -- white
                                       # 1 -- black
                                       # 2 -- dark grey
#------------------------------------------------------------------------------------------#


#----- Profile options. -------------------------------------------------------------------#
my.hours   = sequence(24)-1
prof.hours = c(0,3,6,9,12,15,18,21)                  # Hours to try
prof.zen   = c(0,5,20,30,40,50,60,70,75,80,84,88,92) # Zenith angle breaks
#------------------------------------------------------------------------------------------#


#----- Default break for LAI classes. -----------------------------------------------------#
lai.defbrk  = c(  1.5,  2.0,  2.5,  3.0,  3.5,  4.5,  5.5,  7.0)
add.others  = FALSE
report      = list()
#---------------------------------------   1.5   2.0   2.5   3.0   3.5   4.5   5.5   7.0  -#
report[[1]] = list( solar   = rbind( min = c(  NA ,  NA , 22.0,  8.9,  5.4,  NA ,  NA )
                                   , max = c(  NA ,  NA , 45.4, 11.3,  6.7,  NA ,  NA )
                                   )#end rbind
                  , par     = rbind( min = c(  NA ,  NA ,  NA ,  NA ,  NA ,  NA ,  NA )
                                   , max = c(  NA ,  NA ,  NA ,  NA ,  NA ,  NA ,  NA )
                                   )#end rbind
                  , nir     = rbind( min = c(  NA ,  NA ,  NA ,  NA ,  NA ,  NA ,  NA )
                                   , max = c(  NA ,  NA ,  NA ,  NA ,  NA ,  NA ,  NA )
                                   )#end rbind
                  , colour  = "steelblue1"
                  , angle   = 20
                  , density = 40
                  , lwd     = 2
                  , desc    = "Torquebiau and Walter (1987)"
                  )#end report
#---------------------------------------   1.5   2.0   2.5   3.0   3.5   4.5   5.5   7.0  -#
report[[2]] = list( solar   = rbind( min = c(  NA ,  NA ,  NA ,  NA ,  NA ,  NA ,  NA )
                                   , max = c(  NA ,  NA ,  NA ,  NA ,  NA ,  NA ,  NA )
                                   )#end rbind
                  , par     = rbind( min = c(  NA ,  NA , 10.8,  2.6, 0.35,  NA ,  NA )
                                   , max = c(  NA ,  NA , 22.2,  3.3, 1.04,  NA ,  NA )
                                   )#end rbind
                  , nir     = rbind( min = c(  NA ,  NA ,  NA ,  NA ,  NA ,  NA ,  NA )
                                   , max = c(  NA ,  NA ,  NA ,  NA ,  NA ,  NA ,  NA )
                                   )#end rbind
                  , colour  = "mediumpurple1"
                  , angle   = 70
                  , density = 20
                  , lwd     = 2
                  , desc    = "Torquebiau (1988)"
                  )#end report
#---------------------------------------   1.5   2.0   2.5   3.0   3.5   4.5   5.5   7.0  -#
report[[3]] = list( solar   = rbind( min = c(  NA ,  NA ,  NA ,  NA ,  NA ,  NA , 0.89)
                                   , max = c(  NA ,  NA ,  NA ,  NA ,  NA ,  NA , 1.61)
                                   )#end rbind
                  , par     = rbind( min = c(  NA ,  NA ,  NA ,  NA ,  NA ,  NA ,  NA )
                                   , max = c(  NA ,  NA ,  NA ,  NA ,  NA ,  NA ,  NA )
                                   )#end rbind
                  , nir     = rbind( min = c(  NA ,  NA ,  NA ,  NA ,  NA ,  NA ,  NA )
                                   , max = c(  NA ,  NA ,  NA ,  NA ,  NA ,  NA ,  NA )
                                   )#end rbind
                  , colour  = "sandybrown"
                  , angle   = -20
                  , density =  40
                  , lwd     =   2
                  , desc    = "Shuttleworth (1984)"
                  )#end report
#------------------------------------------------------------------------------------------#





#----- Info on hourly data. ---------------------------------------------------------------#
reload.profile  = TRUE
rdata.path      = file.path(here,"RData_profile")
rdata.suffix    = "radprof_ed22.RData"
finished.suffix = "radprof_ed22.txt"
#------------------------------------------------------------------------------------------#



#----- Info on hourly data. ---------------------------------------------------------------#
reload.summary  = TRUE
rdata.summary   = "summary_radprof_ed22.RData"
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
sites[[1]] = list(iata="gyf",desc="Paracou"       ,lai=6.0,pch= 2,col="#520485",year=2004)
sites[[2]] = list(iata="s67",desc="Santarem km 67",lai=5.0,pch= 5,col="#46FF32",year=2003)
sites[[3]] = list(iata="s83",desc="Santarem km 83",lai=5.0,pch= 9,col="#FF5700",year=2001)
sites[[4]] = list(iata="pdg",desc="Pe-de-Gigante" ,lai=3.5,pch=13,col="#A00014",year=2003)
sites[[5]] = list(iata="rja",desc="Rebio Jaru"    ,lai=4.0,pch= 1,col="#006715",year=2001)
sites[[6]] = list(iata="m34",desc="Manaus K34"    ,lai=5.8,pch= 6,col="#0742C3",year=2004)
sites[[7]] = list(iata="pnz",desc="Petrolina"     ,lai=2.0,pch= 4,col="#B49ED2",year=2004)
sites[[8]] = list(iata="ban",desc="Bananal"       ,lai=3.0,pch= 8,col="#F5C858",year=2005)
use.sites  = "m34" # "rja"
#------------------------------------------------------------------------------------------#



#----- Flags to control whether to plot site-specific and multi-site comparison. ----------#
plot.site   = c(FALSE,TRUE)[2]
plot.patch  = c(FALSE,TRUE)[1]
plot.month  = c(FALSE,TRUE)[2]
plot.zenith = c(FALSE,TRUE)[2]
plot.gloom  = c(FALSE,TRUE)[2]
dlai        = 0.25
ci.estimate = "six.summary"    # How to define the statistics?
                               # "boot.summary" -- Use bootstrap to get CI
                               # "six.summary"  -- Use t distribution to get CI
n.boot      = 1000
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



#----- Load profiles from ABRACOS. --------------------------------------------------------#
load(file.path(srcdir,"ABRACOS.light.RData"))
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     List all simulations.                                                                #
#------------------------------------------------------------------------------------------#
simul       = list()
simul[[ 1]] = list( suff           = "icanrad00_crown01"
                  , desc           = "DM06"
                  , colour         = "#0059B3"
                  , bgcol          = "#537CA6"
                  , mult.down      = 1.
                  , mult.beam.down = 1.
                  , mult.diff.down = 0.97
                  , mult.up        = 1.
                  , angle          = -30
                  , density        = 25
                  )#end list
simul[[ 2]] = list( suff           = "icanrad00_crown00"
                  , desc           = "DC08"
                  , colour         = "#4D0099"
                  , bgcol          = "#8659B3"
                  , mult.down      = 1.
                  , mult.beam.down = 1.
                  , mult.diff.down = 1.03
                  , mult.up        = 0.985
                  , angle          = 30
                  , density        = 25
                  )#end list
simul[[ 3]] = list( suff           = "icanrad01_crown00"
                  , desc           = "ZQ05"
                  , colour         = "#B3B300"
                  , bgcol          = "#A6A653"
                  , mult.down      = 1.03
                  , mult.beam.down = 1.03
                  , mult.diff.down = 1.
                  , mult.up        = 1.
                  , angle          = -60
                  , density        = 25
                  )#end list
simul[[ 4]] = list( suff           = "icanrad02_crown00"
                  , desc           = "ED22"
                  , colour         = "#446600"
                  , bgcol          = "#AACC66"
                  , mult.down      = 0.97
                  , mult.beam.down = 0.97
                  , mult.diff.down = 1.
                  , mult.up        = 1.015
                  , angle          = 60
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
nreport = length(report)
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




#------------------------------------------------------------------------------------------#
#      Create all output directories, separated by format.                                 #
#------------------------------------------------------------------------------------------#
if (! file.exists(rdata.path)) dir.create(rdata.path)
if (! file.exists(outroot   )) dir.create(outroot   )
#------------------------------------------------------------------------------------------#


#----- Summary file. ----------------------------------------------------------------------#
summary.fullname = file.path(rdata.path,rdata.summary)
#------------------------------------------------------------------------------------------#



#----- Zenith angle bins. -----------------------------------------------------------------#
nzb       = length(prof.zen)
zen.key   = paste("zen_",sprintf("%2.2i",prof.zen[-nzb]),"-"
                        ,sprintf("%2.2i",prof.zen[-1]),sep="")
zen.desc  = paste("Mean zenith angle ",sprintf("%.1f",mid.points(prof.zen)),sep="")
nzen.prof = length(zen.key)
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
      yeara    = sites$year[p]
      yearz    = sites$year[p]
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
         alldays  = chron(seq( from = as.numeric(chron(paste( 1, 1,yeara,sep="/")) )
                             , to   = as.numeric(chron(paste(12,31,yearz,sep="/")) )
                             , by   = 1
                             )#end seq
                         )#end chron
         ndates   = length(alldays)
         when     = chron( dates = rep(paste(alldays),each=nhours)
                         , times = paste(rep(my.hours,times=ndates),0,0,sep=":")
                         )#end chron
         nwhen    = length(when)
         zenith   = mapply( FUN      = function(bef,now,...){
                                          zz = ed.zen(when=c(bef,now),...)
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
         model    = list( iata         = iata
                        , short        = short
                        , longname     = desc
                        , lon          = lon
                        , lat          = lat
                        , yeara        = yeara
                        , yearz        = yearz
                        , tresume      = 1
                        , nwhen        = nwhen
                        , when         = when
                        , year         = numyears (when)
                        , month        = nummonths(when)
                        , day          = numdays  (when)
                        , hour         = hours    (when)
                        , zen          = zenith["zen",]
                        , diel         = round( 1 + zenith["day",] - zenith["night",])
                        , gloom        = rep(FALSE,times=nwhen)
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
      tresume   = model$tresume
      nwhen     = model$nwhen
      when      = model$when
      #------------------------------------------------------------------------------------#


      #----- File name for this site. -----------------------------------------------------#
      rdata.iata = file.path(rdata.path,paste(iata,rdata.suffix,sep="_"))
      #------------------------------------------------------------------------------------#

      #---- Get the times that must be read. ----------------------------------------------#
      if (tresume > nwhen){
         loop.when = numeric(0)
      }else{
         loop.when = seq(from=tresume,to=nwhen,by=1)
      }#end if
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Loop over dates.                                                               #
      #------------------------------------------------------------------------------------#
      for (w in loop.when){


         #----- Grab date/time information on the data set. -------------------------------#
         now   = when[w]
         y     = numyears (now); yyyy = sprintf("%4.4i",y)
         m     = nummonths(now); mm   = sprintf("%2.2i",m)
         d     = numdays  (now); dd   = sprintf("%2.2i",d)
         h     = hours    (now); hh   = sprintf("%2.2i",h)


         #---------------------------------------------------------------------------------#
         #      Save data every 4 months, to avoid R collapse.                             #
         #---------------------------------------------------------------------------------#
         time.to.save = ( ! (w %in% c(tresume,nwhen)) && m %in% c(1,5,9)
                        &&   d %in% c(1) && h %in% c(0) )
         if (time.to.save){
            #----- Reload data and copy to the general list. ------------------------------#
            rdata.iata = file.path(rdata.path,paste(iata,rdata.suffix,sep="_"))
            cat(" + Save partial data to ",basename(rdata.iata),"...","\n",sep="")
            model$tresume = w
            dummy         = save(model,file=rdata.iata)
            rm(model)
            q("no")
            #------------------------------------------------------------------------------#
         }else{
            cat("   * Time ",paste(dates(when[w]))," ",hh,"UTC","\n",sep="")
         }#end if
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#



         #----- Build file suffix. --------------------------------------------------------#
         tsuff = paste(yyyy,"-",mm,"-",dd,"-",hh,"0000-g01.h5",sep="")
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #    Loop over simulations.                                                       #
         #---------------------------------------------------------------------------------#
         for (s in sequence(nsimul)){

            #----- Build the file name. ---------------------------------------------------#
            simplace = paste("t",iata,"_",simul$suff[s],sep="")
            filenow  = paste(simplace,"I",tsuff,sep="-")
            #------------------------------------------------------------------------------#


            #----- Make file name as well as the compressed names. ------------------------#
            fullfile     = file.path(here,simplace,"analy",filenow)
            fullfile.gz  = paste(fullfile,"gz" ,sep=".")
            fullfile.bz2 = paste(fullfile,"bz2",sep=".")
            #------------------------------------------------------------------------------#



            #----- Open the file. ---------------------------------------------------------#
            if (file.exists(fullfile)){
               myinst    = hdf5load(file=fullfile,load=FALSE,tidy=TRUE,verbosity=0)
            }else if (file.exists(fullfile.gz )){
               temp.file = file.path(tempdir(),basename(fullfile))
               dummy     = gunzip(filename=fullfile.gz,destname=temp.file,remove=FALSE)
               myinst    = hdf5load(file=temp.file,load=FALSE,verbosity=0,tidy=TRUE)
               dummy     = file.remove(temp.file)
            }else if (file.exists(fullfile.bz2)){
               temp.file = file.path(tempdir(),basename(fullfile))
               dummy     = bunzip2(filename=fullfile.bz2,destname=temp.file,remove=FALSE)
               myinst    = hdf5load(file=temp.file,load=FALSE,verbosity=0,tidy=TRUE)
               dummy     = file.remove(temp.file)
            }else{
               stop ("File not found")
            }#end if
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Append some auxiliary variables to myinst.                              #
            #------------------------------------------------------------------------------#
            myinst$SIPA.N                   = length(myinst$PACO.ID)
            myinst$FMEAN.ATM.NIR.PY         = ( myinst$FMEAN.ATM.RSHORT.PY
                                              - myinst$FMEAN.ATM.PAR.PY         )
            myinst$FMEAN.ATM.NIR.DIFF.PY    = ( myinst$FMEAN.ATM.RSHORT.DIFF.PY
                                              - myinst$FMEAN.ATM.PAR.DIFF.PY    )
            myinst$FMEAN.ATM.RSHORT.BEAM.PY = ( myinst$FMEAN.ATM.RSHORT.PY
                                              - myinst$FMEAN.ATM.RSHORT.DIFF.PY )
            myinst$FMEAN.ATM.PAR.BEAM.PY    = ( myinst$FMEAN.ATM.PAR.PY
                                              - myinst$FMEAN.ATM.PAR.DIFF.PY    )
            myinst$FMEAN.ATM.NIR.BEAM.PY    = ( myinst$FMEAN.ATM.NIR.PY
                                              - myinst$FMEAN.ATM.NIR.DIFF.PY    )
            myinst$FMEAN.ATM.RLONG.DIFF.PY  = myinst$FMEAN.ATM.RLONG.PY
            myinst$FMEAN.ATM.RLONG.BEAM.PY  = 0. * myinst$FMEAN.ATM.RLONG.PY
            myinst$FMEAN.ATM.NIR.SI         = ( myinst$FMEAN.ATM.RSHORT.SI
                                              - myinst$FMEAN.ATM.PAR.SI         )
            myinst$FMEAN.ATM.NIR.DIFF.SI    = ( myinst$FMEAN.ATM.RSHORT.DIFF.SI
                                              - myinst$FMEAN.ATM.PAR.DIFF.SI    )
            myinst$FMEAN.ATM.RSHORT.BEAM.SI = ( myinst$FMEAN.ATM.RSHORT.SI
                                              - myinst$FMEAN.ATM.RSHORT.DIFF.SI )
            myinst$FMEAN.ATM.PAR.BEAM.SI    = ( myinst$FMEAN.ATM.PAR.SI
                                              - myinst$FMEAN.ATM.PAR.DIFF.SI    )
            myinst$FMEAN.ATM.NIR.BEAM.SI    = ( myinst$FMEAN.ATM.NIR.SI
                                              - myinst$FMEAN.ATM.NIR.DIFF.SI    )
            myinst$FMEAN.ATM.RLONG.DIFF.SI  = myinst$FMEAN.ATM.RLONG.SI
            myinst$FMEAN.ATM.RLONG.BEAM.SI  = 0. * myinst$FMEAN.ATM.RLONG.SI
            myinst$FMEAN.RSHORTUP.PA        = ( myinst$FMEAN.PARUP.PA 
                                              + myinst$FMEAN.NIRUP.PA )
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Make the radiation at the top of the atmosphere a patch variables.       #
            #------------------------------------------------------------------------------#
            myinst$FMEAN.PAR.BEAM.PA    = rep( x     = myinst$FMEAN.ATM.PAR.BEAM.SI   
                                             , times = myinst$SIPA.N
                                             )#end rep
            myinst$FMEAN.PAR.DIFF.PA    = rep( x     = myinst$FMEAN.ATM.PAR.DIFF.SI   
                                             , times = myinst$SIPA.N
                                             )#end rep
            myinst$FMEAN.NIR.BEAM.PA    = rep( x     = myinst$FMEAN.ATM.NIR.BEAM.SI   
                                             , times = myinst$SIPA.N
                                             )#end rep
            myinst$FMEAN.NIR.DIFF.PA    = rep( x     = myinst$FMEAN.ATM.NIR.DIFF.SI   
                                             , times = myinst$SIPA.N
                                             )#end rep
            myinst$FMEAN.RSHORT.DIFF.PA = rep( x     = myinst$FMEAN.ATM.RSHORT.DIFF.SI
                                             , times = myinst$SIPA.N
                                             )#end rep
            myinst$FMEAN.RSHORT.BEAM.PA = rep( x     = myinst$FMEAN.ATM.RSHORT.BEAM.SI
                                             , times = myinst$SIPA.N
                                             )#end rep
            myinst$FMEAN.RLONG.DIFF.PA  = rep( x     = myinst$FMEAN.ATM.RLONG.DIFF.SI 
                                             , times = myinst$SIPA.N
                                             )#end rep
            myinst$FMEAN.RLONG.BEAM.PA  = rep( x     = myinst$FMEAN.ATM.RLONG.BEAM.SI 
                                             , times = myinst$SIPA.N
                                             )#end rep
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Correct profiles.                                                       #
            #------------------------------------------------------------------------------#
            if (sum(myinst$PACO.N) > 1){
               myinst$FMEAN.PAR.BEAM.CO    = myinst$FMEAN.RAD.PROFILE[,1]
               myinst$FMEAN.PAR.DIFF.CO    = myinst$FMEAN.RAD.PROFILE[,3]
               myinst$FMEAN.PARUP.CO       = rowSums(myinst$FMEAN.RAD.PROFILE[,c(2,4)])
               myinst$FMEAN.NIR.BEAM.CO    = myinst$FMEAN.RAD.PROFILE[,5]
               myinst$FMEAN.NIR.DIFF.CO    = myinst$FMEAN.RAD.PROFILE[,7]
               myinst$FMEAN.NIRUP.CO       = rowSums(myinst$FMEAN.RAD.PROFILE[,c(6,8)])
               myinst$FMEAN.RLONG.BEAM.CO  = 0. * myinst$FMEAN.RAD.PROFILE[, 9]
               myinst$FMEAN.RLONG.DIFF.CO  = myinst$FMEAN.RAD.PROFILE[, 9]
               myinst$FMEAN.RLONGUP.CO     = myinst$FMEAN.RAD.PROFILE[,10]
               myinst$FMEAN.RSHORT.BEAM.CO = ( myinst$FMEAN.PAR.BEAM.CO
                                             + myinst$FMEAN.NIR.BEAM.CO
                                             )#end
               myinst$FMEAN.RSHORT.DIFF.CO = ( myinst$FMEAN.PAR.DIFF.CO
                                             + myinst$FMEAN.NIR.DIFF.CO
                                             )#end
               myinst$FMEAN.RSHORTUP.CO    = ( myinst$FMEAN.PARUP.CO
                                             + myinst$FMEAN.NIRUP.CO
                                             )#end
               myinst$FMEAN.NIR.L.CO       = ( myinst$FMEAN.RSHORT.L.CO
                                             - myinst$FMEAN.PAR.L.CO
                                             )#end
            }else if (sum(myinst$PACO.N)== 1){
               myinst$FMEAN.PAR.BEAM.CO    = myinst$FMEAN.RAD.PROFILE[1]
               myinst$FMEAN.PAR.DIFF.CO    = myinst$FMEAN.RAD.PROFILE[3]
               myinst$FMEAN.PARUP.CO       = sum(myinst$FMEAN.RAD.PROFILE[c(2,4)])
               myinst$FMEAN.NIR.BEAM.CO    = myinst$FMEAN.RAD.PROFILE[5]
               myinst$FMEAN.NIR.DIFF.CO    = myinst$FMEAN.RAD.PROFILE[7]
               myinst$FMEAN.NIRUP.CO       = sum(myinst$FMEAN.RAD.PROFILE[c(6,8)])
               myinst$FMEAN.RLONG.BEAM.CO  = 0. * myinst$FMEAN.RAD.PROFILE[9]
               myinst$FMEAN.RLONG.DIFF.CO  = myinst$FMEAN.RAD.PROFILE[9]
               myinst$FMEAN.RLONGUP.CO     = myinst$FMEAN.RAD.PROFILE[10]
               myinst$FMEAN.RSHORT.BEAM.CO = ( myinst$FMEAN.PAR.BEAM.CO
                                             + myinst$FMEAN.NIR.BEAM.CO
                                             )#end
               myinst$FMEAN.RSHORT.DIFF.CO = ( myinst$FMEAN.PAR.DIFF.CO
                                             + myinst$FMEAN.NIR.DIFF.CO
                                             )#end
               myinst$FMEAN.RSHORTUP.CO    = ( myinst$FMEAN.PARUP.CO
                                             + myinst$FMEAN.NIRUP.CO
                                             )#end
               myinst$FMEAN.NIR.L.CO       = ( myinst$FMEAN.RSHORT.L.CO
                                             - myinst$FMEAN.PAR.L.CO
                                             )#end
            }#end if
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #      Load the patch and cohort information.                                  #
            #------------------------------------------------------------------------------#
            if (w == 1 && s == 1){
               nsites    = myinst$PYSI.N
               npatches  = myinst$SIPA.N
               ncohorts  = myinst$PACO.N
               a.ico     = myinst$PACO.ID
               z.ico     = a.ico + ncohorts - 1

               isi       = rep(sequence(nsites),times=npatches)
               ipa       = sequence(npatches)
               ipaco     = rep(ipa,times=ncohorts  )
               ppaco     = rep(ipa,times=ncohorts+1)
               ico       = sequence(ncohorts)
               pco       = sequence(ncohorts+1)-1

               nlayer    = sum(ncohorts+1)

               #----- Initialise data regarding this polygon. -----------------------------#
               empty.atm   = array (data= NA,dim=c(nsimul,nwhen))
               empty.prof  = array (data= NA,dim=c(nsimul,nwhen,nlayer))
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
               #      Append initial conditions to this polygon.                           #
               #---------------------------------------------------------------------------#
               model = modifyList( x   = model
                                 , val = list( npatches = npatches
                                             , ncohorts = ncohorts
                                             , desert   = sum(ncohorts) == 0
                                             , nlayer   = nlayer
                                             , ipa      = ipa
                                             , ipaco    = ipaco
                                             , ppaco    = ppaco
                                             , ico      = ico
                                             , pco      = pco
                                             , lai      = empty.prof
                                             , wai      = empty.prof
                                             , tai      = empty.prof
                                             , par      = empty.list
                                             , nir      = empty.list
                                             , sol      = empty.list
                                             , tir      = empty.list
                                             )#end list
                                 )#end modifyList
               #---------------------------------------------------------------------------#
            }#end if (w == 1 && s == 1)
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Retrieve some useful data from this run.                                 #
            #------------------------------------------------------------------------------#
            npatches = model$npatches
            ipa      = model$ipa
            ipaco    = model$ipaco
            ico      = model$ico
            ppaco    = model$ppa
            pco      = model$pco
            desert   = model$desert
            #------------------------------------------------------------------------------#



            #------ Variables to fill in constant boundary conditions. --------------------#
            zero     = rep(0   ,times=npatches)
            veritas  = rep(TRUE,times=npatches)
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Make the cumulative LAI profile.                                         #
            #------------------------------------------------------------------------------#
            if (desert){
               lai.prof = zero
               wai.prof = zero
               tai.prof = zero
               use.lyr  = veritas
            }else{
               lai.prof = append.patch(ipaco=ipaco,xpa=zero,xco=myinst$LAI.CO)
               wai.prof = append.patch(ipaco=ipaco,xpa=zero,xco=myinst$WAI.CO)
               tai.prof = lai.prof + wai.prof
               use.lyr  = append.patch( ipaco = ipaco
                                      , xpa   = veritas
                                      , xco   = myinst$FMEAN.RLONG.DIFF.CO != 0.
                                      )#end append.patch
            }#end if
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #     Append layer characteristics for this time.                              #
            #------------------------------------------------------------------------------#
            model$lai      [s,w,] = lai.prof
            model$wai      [s,w,] = wai.prof
            model$tai      [s,w,] = tai.prof
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Loop over radiation bands.                                               #
            #------------------------------------------------------------------------------#
            for (b in sequence(nband)){
               bnd  = band$key [b]
               ed2  = band$ed2 [b]



               #----- Get the components. -------------------------------------------------#
               down.beam.top = myinst[[paste("FMEAN.",ed2,".BEAM.PA",sep="")]]
               down.diff.top = myinst[[paste("FMEAN.",ed2,".DIFF.PA",sep="")]]
               up.tot.top    = myinst[[paste("FMEAN.",ed2,"UP.PA"   ,sep="")]]
               down.beam.coh = myinst[[paste("FMEAN.",ed2,".BEAM.CO",sep="")]]
               down.diff.coh = myinst[[paste("FMEAN.",ed2,".DIFF.CO",sep="")]]
               up.tot.coh    = myinst[[paste("FMEAN.",ed2,"UP.CO"   ,sep="")]]
               abs.leaf.coh  = myinst[[paste("FMEAN.",ed2,".L.CO"   ,sep="")]]

               atm.beam      = myinst[[paste("FMEAN.ATM.",ed2,".BEAM.PY",sep="")]]
               atm.diff      = myinst[[paste("FMEAN.ATM.",ed2,".DIFF.PY",sep="")]]
               atm.tot       = atm.beam + atm.diff
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Add top boundary condition to each cohort.                            #
               #---------------------------------------------------------------------------#
               if (desert){
                  #----- Desert run, no cohorts. ------------------------------------------#
                  down.beam = down.beam.top
                  down.diff = down.diff.top
                  up.tot    = up.top
                  abs.leaf  = zero
                  #------------------------------------------------------------------------#
               }else{
                  #----- Append the patch-level top boundary condition. -------------------#
                  down.beam = append.patch(ipaco=ipaco,xpa=down.beam.top,xco=down.beam.coh)
                  down.diff = append.patch(ipaco=ipaco,xpa=down.diff.top,xco=down.diff.coh)
                  up.tot    = append.patch(ipaco=ipaco,xpa=up.tot.top   ,xco=up.tot.coh   )
                  abs.leaf  = append.patch(ipaco=ipaco,xpa=zero         ,xco=abs.leaf.coh )
                  #------------------------------------------------------------------------#



                  #----- Fill in cohorts that haven't been resolved. ----------------------#
                  down.beam = fill.unresolved(ipaco=ppaco,x=down.beam,use=use.lyr)
                  down.diff = fill.unresolved(ipaco=ppaco,x=down.diff,use=use.lyr)
                  up.tot    = fill.unresolved(ipaco=ppaco,x=up.tot   ,use=use.lyr)
                  #------------------------------------------------------------------------#
               }#end if
               #---------------------------------------------------------------------------#



               #----- Find the derived quantities. ----------------------------------------#
               down.tot  = down.beam + down.diff
               abs.lyr   = layer.absorption(ipaco=ppaco,down=down.tot,up=up.tot)
               abs.wood  = abs.lyr - abs.leaf
               abs.tot   = unlist(tapply(X=abs.lyr ,INDEX=ppaco,FUN=cumsum))
               abs.leaf  = unlist(tapply(X=abs.leaf,INDEX=ppaco,FUN=cumsum))
               abs.wood  = unlist(tapply(X=abs.wood,INDEX=ppaco,FUN=cumsum))
               #---------------------------------------------------------------------------#



               #----- Update profile. -----------------------------------------------------#
               model[[bnd]]$down.beam [s,w,] = down.beam
               model[[bnd]]$down.diff [s,w,] = down.diff
               model[[bnd]]$down.tot  [s,w,] = down.tot
               model[[bnd]]$up.tot    [s,w,] = up.tot
               model[[bnd]]$abs.lyr   [s,w,] = abs.lyr
               model[[bnd]]$abs.tot   [s,w,] = abs.tot
               model[[bnd]]$abs.leaf  [s,w,] = abs.leaf
               model[[bnd]]$abs.wood  [s,w,] = abs.wood
               #---------------------------------------------------------------------------#



               #----- Update ToC radiation. -----------------------------------------------#
               model[[bnd]]$atm.beam  [s,w] =  atm.beam
               model[[bnd]]$atm.diff  [s,w] =  atm.diff
               model[[bnd]]$atm.tot   [s,w] =  atm.tot
               #---------------------------------------------------------------------------#
            }#end for (b in sequence(nband))
            #------------------------------------------------------------------------------#



            #------ Save gloomy times. ----------------------------------------------------#
            model$gloom[w] = (model$diel[w] == 2) & (model$sol$atm.beam[w] == 0.)
            #------------------------------------------------------------------------------#
         }#end for (s in sequence(nsimul))
         #---------------------------------------------------------------------------------#
      }#end for (w in sequence(nwhen))
      #------------------------------------------------------------------------------------#

      if (tresume <= nwhen){
         #---------------------------------------------------------------------------------#
         #     Now we find which patch has the most, the least and the median LAI.         #
         #---------------------------------------------------------------------------------#
         model$lai.pa   = qapply(X=model$lai,INDEX =model$ppaco,DIM=3,FUN=sum ,na.rm=TRUE)
         lai.bar        = apply (X=model$lai.pa,MARGIN=3,FUN=mean,na.rm=TRUE)
         model$lai.plot = c( max = which.max(lai.bar)
                           , mid = which.closest(x=median(lai.bar),A=lai.bar)
                           , min = which.min(lai.bar)
                           )#end c
         #---------------------------------------------------------------------------------#



         #----- Reload data and copy to the general list. ---------------------------------#
         rdata.iata     = file.path(rdata.path,paste(iata,rdata.suffix,sep="_"))
         rdata.finished = file.path(rdata.path
                                   ,paste("loaded",iata,finished.suffix,sep="_"))
         cat(" + Saving data to ",basename(rdata.iata),"...","\n",sep="")
         model$tresume = nwhen+1
         dummy         = save(model,file=rdata.iata)
         dummy         = write(x=paste(Sys.time()),file=rdata.finished)
         eft[[iata]]   = model
         rm(model)
         #---------------------------------------------------------------------------------#
      }#end if (tresume <= nwhen)
   }#end for (p in sequence(nsites))
   #=======================================================================================#
   #=======================================================================================#
}#end if ( ( (! file.exists(summary.fullname)) && plot.patch ) || plot.site )
#==========================================================================================#
#==========================================================================================#








#==========================================================================================#
#==========================================================================================#
#     Patch-level summary for all sites.                                                   #
#------------------------------------------------------------------------------------------#
if (! file.exists(summary.fullname) && plot.patch){
   summ = list()


   for (p in sequence(nsites)){
      #----- Retrieve site information. ---------------------------------------------------#
      iata     = sites$iata[p]
      model    = eft[[iata]]
      yeara    = model$yeara
      yearz    = model$yearz
      lon      = model$lon
      lat      = model$lat
      desc     = model$longname
      cat(" - Generating patch-level summary for ",desc,"...","\n")
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #     Load LAI/WAI/TAI profiles.                                                     #
      #------------------------------------------------------------------------------------#
      lai     = qapply(X=model$lai,INDEX=model$month,DIM=2,FUN=mean)
      wai     = qapply(X=model$wai,INDEX=model$month,DIM=2,FUN=mean)
      tai     = qapply(X=model$tai,INDEX=model$month,DIM=2,FUN=mean)
      tot.lai = qapply(X=lai      ,INDEX=model$ppaco,DIM=3,FUN=sum )
      tot.wai = qapply(X=wai      ,INDEX=model$ppaco,DIM=3,FUN=sum )
      tot.tai = qapply(X=tai      ,INDEX=model$ppaco,DIM=3,FUN=sum )
      #------------------------------------------------------------------------------------#



      #----- Initialise arrays that will hold the mean diel. ------------------------------#
      npatches = model$npatches
      empty    = array( data     = NA
                      , dim      = c(nband,nsimul,12,npatches)
                      , dimnames = list(band$key,simul$suff,month.abb,sequence(npatches))
                      )#end array
      flux     = list( down.beam = empty
                     , down.diff = empty
                     , down.tot  = empty
                     , up.tot    = empty
                     )#end list
      #------------------------------------------------------------------------------------#




      #----- Band loop. -------------------------------------------------------------------#
      for (b in sequence(nband)){
         bnd      = band$key [b]
         bnd.desc = band$desc[b]
         prof     = model[[bnd]]


         #----- Find the annual averages for each cohort. ---------------------------------#
         down.beam = qapply(X=prof$down.beam,INDEX=model$month,DIM=2,FUN=mean)
         down.diff = qapply(X=prof$down.diff,INDEX=model$month,DIM=2,FUN=mean)
         down.tot  = qapply(X=prof$down.tot ,INDEX=model$month,DIM=2,FUN=mean)
         up.tot    = qapply(X=prof$up.tot   ,INDEX=model$month,DIM=2,FUN=mean)
         #---------------------------------------------------------------------------------#

         #----- Get only the bottom of the cohort. ----------------------------------------#
         flux$down.beam[b,,,] = qapply(X=down.beam,INDEX=model$ppaco,DIM=3,FUN=rel.last)
         flux$down.diff[b,,,] = qapply(X=down.diff,INDEX=model$ppaco,DIM=3,FUN=rel.last)
         flux$down.tot [b,,,] = qapply(X=down.tot ,INDEX=model$ppaco,DIM=3,FUN=rel.last)
         flux$up.tot   [b,,,] = qapply(X=up.tot   ,INDEX=model$ppaco,DIM=3,FUN=rel.last)
         #---------------------------------------------------------------------------------#
      }#end for (b in sequence(nband))
      #------------------------------------------------------------------------------------#



      #----- Save the summary for this site. ----------------------------------------------#
      summ[[iata]] = list( lai = tot.lai, wai = tot.wai, tai = tot.tai, flux = flux )
      #------------------------------------------------------------------------------------#

   }#end for (p in sequence(nsites))
   #---------------------------------------------------------------------------------------#

   cat(" - Saving patch-level summary...","\n")
   dummy = save(summ,file=summary.fullname)

}else{
   cat(" - Retrieving patch-level summary...","\n")
   dummy = load(summary.fullname)
}#end if (! file.exists(summary.fullname))
#==========================================================================================#
#==========================================================================================#


















#==========================================================================================#
#==========================================================================================#
#     Site-specific plots.                                                                 #
#------------------------------------------------------------------------------------------#
if (plot.patch){

   #----- Create path if it doesn't exist. ------------------------------------------------#
   outpatch = file.path(outroot,"all_sites")
   if (! file.exists(outpatch)) dir.create(outpatch)
   #---------------------------------------------------------------------------------------#



   #----- Initialise arrays that will light extinction data. ------------------------------#
   tot.lai  = array( data     = NA
                   , dim      = c(nsimul,12,0)
                   , dimnames = list(simul$suff,month.abb,NULL)
                   )#end array
   empty    = array( data     = NA
                   , dim      = c(nband,nsimul,12,0)
                   , dimnames = list(band$key,simul$suff,month.abb,NULL)
                   )#end array
   flux     = list( down.beam = empty
                  , down.diff = empty
                  , down.tot  = empty
                  , up.tot    = empty
                  )#end list
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Site loop.                                                                         #
   #---------------------------------------------------------------------------------------#
   for (p in sequence(nsites)){
      #----- Retrieve site information. ---------------------------------------------------#
      iata     = sites$iata[p]
      plai     = summ[[iata]]$lai
      this     = summ[[iata]]$flux
      #------------------------------------------------------------------------------------#



      #----- Append site. -----------------------------------------------------------------#
      tot.lai        = abind(tot.lai       ,plai          ,along=3)
      flux$down.beam = abind(flux$down.beam,this$down.beam,along=4)
      flux$down.diff = abind(flux$down.diff,this$down.diff,along=4)
      flux$down.tot  = abind(flux$down.tot ,this$down.tot ,along=4)
      flux$up.tot    = abind(flux$up.tot   ,this$up.tot   ,along=4)
      #------------------------------------------------------------------------------------#
   }#end for (p in sequence(nsites)
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #     Plot the total extinction of light as a function of cumulative LAI and month.     #
   #---------------------------------------------------------------------------------------#
   cat(" + Multi-site comparison of light transmissivity as a function of patch LAI","\n")


   #---------------------------------------------------------------------------------------#
   #     We display only solar radiation, since thermal radiation scale is completely      #
   # different and has even less information than solar radiation.                         #
   #---------------------------------------------------------------------------------------#
   bwhich     = which(! band$key %in% "tir")
   nbwhich    = length(bwhich)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Define LAI classes.  Make sure to remove empty classes but to include all         #
   # patches.                                                                              #
   #---------------------------------------------------------------------------------------#
   if (add.others){
      lai.extremes = range(pretty.log(range(c(tot.lai)),n=5))
      lai.breaks   = sort(c(lai.defbrk,lai.extremes[! lai.extremes %wr% lai.defbrk]))
   }else{
      lai.breaks   = lai.defbrk
   }#end if
   discard    = ! (tot.lai %wr% lai.breaks)
   tot.lai    = ifelse(discard,mean(lai.breaks),tot.lai)
   lai.labels = pair.paste(x=format(lai.breaks),collapse="-")
   lai.cut    = as.integer(cut(x=c(tot.lai),breaks=lai.breaks))
   idx.cut    = as.integer(names(table(lai.cut)))
   lai.cut    = array(match(lai.cut,idx.cut),dim=dim(tot.lai),dimnames=dimnames(tot.lai))
   n.cut      = length(idx.cut)
   lai.grid   = 0.5 + nbwhich*(sequence(n.cut+1)-1)
   lai.at     = mid.points(lai.grid)
   lai.limit  = range(lai.grid)
   lai.labels = lai.labels[idx.cut]
   mycols     = rep(band$col[bwhich],times=n.cut)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Irradiance component.                                                             #
   #---------------------------------------------------------------------------------------#
   nflux = length(flux)
   for (f in sequence(nflux)){
      #----- Get this flux component. -----------------------------------------------------#
      flux.key   = names(flux)[f]
      flux.desc  = c("Downward (direct)","Downward (diffuse)"
                    ,"Downward (total)" ,"Upward")[f]
      this       = ifelse(rep(discard,each=nband),NA,flux[[f]]) + 0 * flux[[f]]
      this.limit = range(c(this[bwhich,,,]),finite=TRUE)
      y.at       = pretty.log(this.limit,n=4)
      y.labels   = sprintf("%g",100.*y.at)
      #------------------------------------------------------------------------------------#


      #----- Check whether to plot the reports or not. ------------------------------------#
      plot.report = flux.key %in% c("down.tot") && ! add.others
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Adjust legend stuff based on the whether we are plotting reported values.     #
      #------------------------------------------------------------------------------------#
      if (plot.report){
         nrleg                    = max(nbwhich,nreport)
         ncleg                    = 2
         leg.name                 = character(length=nrleg*ncleg)
         leg.fill                 = rep("transparent",times=nrleg*ncleg)
         leg.border               = rep("transparent",times=nrleg*ncleg)
         leg.density              = rep(-1,times=nrleg*ncleg)
         leg.angle                = numeric(length=nrleg*ncleg)
         leg.name   [1:nbwhich]   = c(band$desc[bwhich])
         leg.fill   [1:nbwhich]   = c(band$col [bwhich])
         leg.border [1:nbwhich]   = foreground
         for (r in sequence(nreport)){
             leg.name   [nrleg+r] = report[[r]]$desc
             leg.fill   [nrleg+r] = report[[r]]$colour
             leg.border [nrleg+r] = report[[r]]$colour
             leg.density[nrleg+r] = report[[r]]$density
             leg.angle  [nrleg+r] = report[[r]]$angle
         }#end for
      }else{
         ncleg                    = 1
         leg.density              = rep(-1,times=nbwhich)
         leg.angle                = numeric(length=nbwhich)
         leg.name                 = c(band$desc[bwhich])
         leg.fill                 = c(band$col [bwhich])
         leg.border               = rep(foreground,times=nbwhich)
      }#end if
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Loop over formats.                                                            #
      #------------------------------------------------------------------------------------#
      for (o in sequence(nout)){
         #----- Open file. ----------------------------------------------------------------#
         fichier = file.path(outpatch
                            ,paste("beneath_relative-",flux.key,".",outform[o],sep=""))
         if       (outform[o] == "x11"){
           X11(width=tsize$width,height=tsize$height,pointsize=ptsz)
         }else if (outform[o] == "png"){
           png(filename=fichier,width=tsize$width*depth
              ,height=tsize$height*depth,pointsize=ptsz,res=depth)
         }else if (outform[o] == "eps"){
           postscript(file=fichier,width=tsize$width,height=tsize$height
                     ,pointsize=ptsz,paper=tsize$paper)
         }else if (outform[o] == "pdf"){
           pdf(file=fichier,onefile=FALSE,width=tsize$width
              ,height=tsize$height,pointsize=ptsz,paper=tsize$paper)
         }#end if
         #---------------------------------------------------------------------------------#



         #----- Set up the plotting window. -----------------------------------------------#
         par(par.user)
         par(oma=c(0.2,1,4,0.2))
         layout( mat     = rbind(lo.simul$mat.off,rep(1,times=lo.simul$ncol))
               , heights = c(rep(5/lo.simul$nrow,lo.simul$nrow),1)
               )#end layout
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #      First plot: the legend.  Skip TIR.                                         #
         #---------------------------------------------------------------------------------#
         par(mar=c(0.1,4.1,0.1,2.1))
         plot.new()
         plot.window(xlim=c(0,1),ylim=c(0,1))
         legend( x       = "bottom"
               , inset   = 0.0
               , legend  = leg.name
               , fill    = leg.fill
               , border  = leg.border
               , angle   = leg.angle
               , density = leg.density
               , cex     = cex.ptsz
               , ncol    = ncleg
               , xpd     = TRUE
               )#end legend
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #      Loop over the simulations, then split data by band and LAI class.          #
         #---------------------------------------------------------------------------------#
         for (s in sequence(nsimul)){
            #------ Open plotting device. -------------------------------------------------#
            par(mar=lo.simul$mar0)
            plot.new()
            plot.window(xlim = lai.limit,ylim = this.limit,log="y")
            abline(h=y.at    ,col=grid.colour,lty="solid")
            abline(v=lai.grid,col=grey.fg    ,lty="solid")
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Add some of the reported values in the background.                       #
            #------------------------------------------------------------------------------#
            if (plot.report){
               xleft  = sequence(nbwhich*n.cut) - 0.45
               xright = sequence(nbwhich*n.cut) + 0.45

               #----- Organise the reported values and plot rectangles. -------------------#
               for (r in sequence(nreport)){
                   rnow    = report[[r]]
                   ybottom = c(rbind(rnow$solar[1,],rnow$par[1,],rnow$nir[1,]))*0.01
                   ytop    = c(rbind(rnow$solar[2,],rnow$par[2,],rnow$nir[2,]))*0.01
                   srect( xleft   = xleft
                        , ybottom = ybottom
                        , xright  = xright
                        , ytop    = ytop
                        , density = rnow$density
                        , col     = rnow$colour
                        , border  = rnow$colour
                        , angle   = rnow$angle
                        , lwd     = rnow$lwd
                        )#end srect
               }#end for
               #---------------------------------------------------------------------------------------#
            }#end if
            #------------------------------------------------------------------------------#


            #----- Create the lists to generate the box-and-whisker plots. ----------------#
            this.now  = this[bwhich,s,,]
            f.band    = arrayInd(seq_along(this.now),.dim=dim(this.now))[,1]
            f.lai     = c(rep(lai.cut[s,,],each=nbwhich))
            this.list = split(x=c(this.now),f=list(f.band,f.lai))
            boxplot(x=this.list,col=mycols,axes=FALSE,add=TRUE)
            #------------------------------------------------------------------------------#


            #----- Plot the final information. --------------------------------------------#
            axis   (side=2,at=y.at,labels=y.labels,las=1)
            axis.rt(side=1,las=5,at=lai.at,labels=lai.labels,off=0.075)
            box()
            title(main=simul$desc[s],font.main=1,line=1)
            #------------------------------------------------------------------------------#
         }#end for (s in sequence(nsimul))
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Plot the general annotation.                                                #
         #---------------------------------------------------------------------------------#
         lex     = desc.unit(desc="Patch LAI class",unit=untab$m2lom2)
         ley     = desc.unit(desc="Relative Irradiance",unit=untab$pc)
         letitre = paste("Relative irradiance beneath the canopy"
                        ,"\n",flux.desc,sep="")
         gtitle( main      = letitre
               , xlab      = lex
               , ylab      = ley
               , line.xlab = 4.1
               , line.ylab = 2.6
               , cex.main  = 0.8*cex.ptsz
               , cex.axis  = 0.8*cex.ptsz
               , off.xlab  = 1/9
               )#end gtitle
         #---------------------------------------------------------------------------------#



         #----- Close the device. ---------------------------------------------------------#
         if (outform[o] == "x11"){
            locator(n=1)
            dev.off()
         }else{
            dev.off()
         }#end if
         clean.tmp()
         #---------------------------------------------------------------------------------#
      }#end for (o in sequence(nout))
      #------------------------------------------------------------------------------------#
   }#end for (f in sequence(nflux))
   #---------------------------------------------------------------------------------------#

}#end if (plot.patch)
#==========================================================================================#
#==========================================================================================#







#==========================================================================================#
#==========================================================================================#
#     Site-specific plots.                                                                 #
#------------------------------------------------------------------------------------------#
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
      cat(" - Generating profiles for ",desc,"...","\n")
      #------------------------------------------------------------------------------------#


      #----- Make output path for this place. ---------------------------------------------#
      sitepath = file.path(outroot,model$short)
      if (! file.exists(sitepath)) dir.create(sitepath)
      #------------------------------------------------------------------------------------#


      #----- Find the zenith classes.  Make "NA" a class after valid classes. -------------#
      model$zen.class = as.integer(cut(model$zen,breaks=prof.zen))
      #------------------------------------------------------------------------------------#


      #----- Choose the patches to plot based on the LAI. ---------------------------------#
      cat("   > LAI by month...","\n")
      if (plot.month){
         lai.mon  = qapply(X = model$lai.pa, INDEX  = model$month, DIM = 2, FUN = mean)
         lai.mon  = apply (X = lai.mon     , MARGIN = c(2,3), FUN = mean)
         loop.mon = sequence(12)
      }else{
         loop.mon = numeric(0)
      }#end if (plot.month)
      #------------------------------------------------------------------------------------#



      #----- Choose the patches to plot based on the LAI. ---------------------------------#
      cat("   > LAI by zenith angle class...","\n")
      if (plot.zenith){
         lai.1st  = qapply(X = model$lai.pa, INDEX  = model$zen.class, DIM = 2, FUN = mean)
         lai.1st  = apply (X = lai.1st     , MARGIN = c(2,3), FUN = mean)

         #----- Copy lai class to a full matrix. ------------------------------------------#
         idx          = cbind( row = as.numeric(rownames(lai.1st)[row(lai.1st)])
                             , col = as.numeric(colnames(lai.1st)[col(lai.1st)])
                             )#end cbind
         lai.zen      = matrix(data=NA,nrow=nzen.prof,ncol=model$npatches)
         lai.zen[idx] = lai.1st
         rm(lai.1st)
         #---------------------------------------------------------------------------------#

         loop.zen = sequence(nzen.prof)
      }else{
         loop.zen = numeric(0)
      }#end if
      #------------------------------------------------------------------------------------#



      #----- Select patches to be plotted. ------------------------------------------------#
      lai.bar      = apply (X=model$lai.pa,MARGIN=3,FUN=mean,na.rm=TRUE)
      use.lai.plot = c( ref = which.closest( x    = sites$lai[p]
                                           , A    = lai.bar
                                           , mask = lai.bar-sites$lai[p]>0
                                           )#end which.closest
                      , max = which.max(lai.bar)
                      , mid = which.closest(x=median(lai.bar),A=lai.bar)
                      , min = which.min(lai.bar)
                      )#end c
      use.lai.key  = c("lairef","laimax","laimed","laimin")
      use.lai.desc = c("Default","Greenest","Medium","Brownest")
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




      #------------------------------------------------------------------------------------#
      #     Loop over the LAI profiles we will plot.                                       #
      #------------------------------------------------------------------------------------#
      for (l in seq_along(use.lai.plot)){
         pnow      = use.lai.plot[l]
         lai.key   = use.lai.key [l]
         lai.desc  = use.lai.desc[l]
         pp        = model$ppaco == pnow
         ptb       = c(which.min(which(pp)),which.max(which(pp)))

         cat("   * ",lai.desc," patch","\n")


         #----- Loop over months and hours. -----------------------------------------------#
         for (m in loop.mon){
            mm = sprintf("%2.2i",m)
            cat("     > ",month.name[m],"\n")
            cat("       ~ Average profile by hour of the day...","\n")


            for (h in prof.hours){
               hh      = sprintf("%2.2i",h)
               mh      = model$month == m & model$hour == h
               zen.bar = mean(model$zen[mh])



               #------ Create path for this hour. -----------------------------------------#
               profhrpath  = file.path(profilepath,paste(hh,"utc",sep=""))
               if (! file.exists(profhrpath)) dir.create(profhrpath)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Load LAI/WAI/TAI profiles.                                            #
               #---------------------------------------------------------------------------#
               lai     = apply(X=model$lai[,mh,pp,drop=FALSE],MARGIN=c(1,3),FUN=mean)
               wai     = apply(X=model$wai[,mh,pp,drop=FALSE],MARGIN=c(1,3),FUN=mean)
               tai     = apply(X=model$tai[,mh,pp,drop=FALSE],MARGIN=c(1,3),FUN=mean)
               cum.lai = t(apply(X=lai,MARGIN=1,FUN=cumsum))
               cum.wai = t(apply(X=wai,MARGIN=1,FUN=cumsum))
               cum.tai = t(apply(X=tai,MARGIN=1,FUN=cumsum))
               #---------------------------------------------------------------------------#


               #----- Fix the limits for the Y axis. --------------------------------------#
               lai.lim = range(-cum.lai)
               lai.at   = pretty(lai.lim)
               lai.lab  = -lai.at
               #---------------------------------------------------------------------------#



               #----- Band loop. ----------------------------------------------------------#
               for (b in sequence(nband)){
                  bnd      = band$key [b]
                  bnd.desc = band$desc[b]

                  mult     = if (bnd %in% c("par")){Watts.2.Ein * 1.e6}else{1.0}
                  prof     = model[[bnd]]

                  #----- Find the averages for the hour and month. ------------------------#
                  down.beam = apply( X      = prof$down.beam[,mh,pp,drop=FALSE] * mult
                                   , MARGIN = c(1,3)
                                   , FUN    = mean
                                   )#end apply
                  down.diff = apply( X      = prof$down.diff[,mh,pp,drop=FALSE] * mult
                                   , MARGIN = c(1,3)
                                   , FUN    = mean
                                   )#end apply
                  down.tot  = apply( X      = prof$down.tot [,mh,pp,drop=FALSE] * mult
                                   , MARGIN = c(1,3)
                                   , FUN    = mean
                                   )#end apply
                  up.tot    = apply( X      = prof$up.tot   [,mh,pp,drop=FALSE] * mult
                                   , MARGIN = c(1,3)
                                   , FUN    = mean
                                   )#end apply
                  abs.lyr   = apply( X      = prof$abs.lyr  [,mh,pp,drop=FALSE] * mult
                                   , MARGIN = c(1,3)
                                   , FUN    = mean
                                   )#end apply
                  abs.tot   = apply( X      = prof$abs.tot  [,mh,pp,drop=FALSE] * mult
                                   , MARGIN = c(1,3)
                                   , FUN    = mean
                                   )#end apply
                  abs.leaf  = apply( X      = prof$abs.leaf [,mh,pp,drop=FALSE] * mult
                                   , MARGIN = c(1,3)
                                   , FUN    = mean
                                   )#end apply
                  abs.wood  = apply( X      = prof$abs.wood [,mh,pp,drop=FALSE] * mult
                                   , MARGIN = c(1,3)
                                   , FUN    = mean
                                   )#end apply
                  #------------------------------------------------------------------------#




                  #------------------------------------------------------------------------#
                  #     Shift the results a bit so we can see when two curves are on top   #
                  # of each other.                                                         #
                  #------------------------------------------------------------------------#
                  if ( ! bnd %in% "tir"){
                     down.beam = apply(X = down.beam,MARGIN=2,FUN="*",simul$mult.beam.down)
                     down.diff = apply(X = down.diff,MARGIN=2,FUN="*",simul$mult.diff.down)
                     up.tot    = apply(X = up.tot   ,MARGIN=2,FUN="*",simul$mult.up       )
                  }#end if
                  #------------------------------------------------------------------------#





                  #------------------------------------------------------------------------#
                  #     Skip plot if it is solar radiation and it is nighttime.            #
                  #------------------------------------------------------------------------#
                  if (max(down.tot) > 1.){


                     #----- Retrieve the data. --------------------------------------------#
                     if (! bnd %in% "tir"){
                        down.lim = pretty.xylim(u=c(down.beam,down.diff))
                     }else{
                        down.lim = pretty.xylim(u=down.diff)
                     }#end if
                     up.lim   = pretty.xylim(u=up.tot )
                     abs.lim  = pretty.xylim(u=c(0,abs.tot,abs.leaf))
                     #---------------------------------------------------------------------#


                     for (o in sequence(nout)){
                        #----- Open file. -------------------------------------------------#
                        fichier = file.path(profhrpath
                                           ,paste("profile-",iata,"-",bnd,"-",mm,"-"
                                                 ,hh,"utc-",lai.key,".",outform[o],sep=""))
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
                        #------------------------------------------------------------------#


                        #----- Set up the plotting window. --------------------------------#
                        par(par.user)
                        par(oma=c(0,1,4,0))
                        if (bnd %in% "tir"){
                           layout(mat=rbind(c(2,3,4),c(1,1,1)),heights=c(6,1))
                        }else{
                           layout(mat=rbind(c(2,4),c(3,5),c(1,1)),heights=c(3,3,1))
                        }#end if
                        #------------------------------------------------------------------#




                        #----- First plot: the legend. ------------------------------------#
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
                              , ncol   = 2
                              , xpd    = TRUE
                              )#end legend
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #      Direct radiation is not plotted for TIR.                    #
                        #------------------------------------------------------------------#
                        if (! bnd %in% "tir"){
                           #------ Second plot: Downward direct. --------------------------#
                           par(mar=c(3.1,3.1,3.1,1.1))
                           plot.new()
                           plot.window(xlim=down.lim,ylim=lai.lim)
                           axis(side=1,cex.axis=0.8,padj=-0.8)
                           axis(side=2,at=lai.at,labels=lai.lab,las=1)
                           grid(col=grid.colour,lwd=1.0,lty="solid")
                           box()
                           title(main="Downward - Direct",font.main=1,line=1)
                           for (s in sequence(nsimul)){
                              lines(x=down.beam[s,],y=-cum.lai[s,]
                                   ,col=simul$colour[s],lwd=2.5,lty="solid")
                           }#end for
                           #---------------------------------------------------------------#
                        }#end if
                        #------------------------------------------------------------------#




                        #------ Third plot: Downward diffuse. -----------------------------#
                        par(mar=c(4.1,3.1,2.1,1.1))
                        plot.new()
                        plot.window(xlim=down.lim,ylim=lai.lim)
                        axis(side=1,cex.axis=0.8,padj=-0.8)
                        axis(side=2,at=lai.at,labels=lai.lab,las=1)
                        grid(col=grid.colour,lwd=1.0,lty="solid")
                        box()
                        title(main="Downward - Diffuse",font.main=1,line=1)
                        for (s in sequence(nsimul)){
                           lines(x=down.diff[s,],y=-cum.lai[s,]
                                ,col=simul$colour[s],lwd=2.5,lty="solid")
                        }#end for
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #     Fourth plot: Upward radiation.  Margins depend on the band.  #
                        #------------------------------------------------------------------#
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
                           lines(x=up.tot[s,],y=-cum.lai[s,]
                                ,col=simul$colour[s],lwd=2.5,lty="solid")
                        }#end for
                        #------------------------------------------------------------------#



                        #------ Fifth plot: Absorbed radiation. ---------------------------#
                        par(mar=c(4.1,2.1,2.1,2.1))
                        plot.new()
                        plot.window(xlim=abs.lim,ylim=lai.lim)
                        axis(side=1,cex.axis=0.8,padj=-0.8)
                        axis(side=2,at=lai.at,labels=lai.lab,las=1)
                        grid(col=grid.colour,lwd=1.0,lty="solid")
                        box()
                        title(main="Cumulative absorption (dashed = leaves)"
                             ,font.main=1,line=1)
                        for (s in sequence(nsimul)){
                           lines(x=abs.leaf[s,],y=-cum.lai[s,]
                                ,col=simul$colour[s],lwd=2.0,lty="dotdash")
                           lines(x=abs.tot [s,],y=-cum.lai[s,]
                                ,col=simul$colour[s],lwd=2.5,lty="solid")
                        }#end for
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #     Plot the general annotation.                                 #
                        #------------------------------------------------------------------#
                        lex     = desc.unit( desc = paste(band$desc[b],"Irradiance")
                                           , unit = untab[[band$unit[b]]]
                                           )#end desc.unit
                        ley     = desc.unit(desc="Cumulative LAI",unit=untab$m2lom2)
                        whenlab = paste(month.name[m],"-",hh,"UTC")
                        letitre = paste(desc," - Canopy profiles","  -  ",whenlab,sep="")
                        lelai   = desc.unit( desc = paste("Average LAI"
                                                         ,sprintf("%6.2f",lai.mon[m,pnow])
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
                              , cex.main  = 1.1*cex.ptsz
                              , cex.xlab  = 0.9*cex.ptsz
                              , cex.ylab  = 0.9*cex.ptsz
                              , cex.sub   = 0.7*cex.ptsz
                              )#end gtitle
                        #------------------------------------------------------------------#



                        #----- Close the device. ------------------------------------------#
                        if (outform[o] == "x11"){
                           locator(n=1)
                           dev.off()
                        }else{
                           dev.off()
                        }#end if
                        clean.tmp()
                        #------------------------------------------------------------------#
                     }#end for (o in sequence(nout))
                     #---------------------------------------------------------------------#
                  }#end if (max(tot.down) > 1.)
                  #------------------------------------------------------------------------#
               }#end for (b in sequence(nband))
               #---------------------------------------------------------------------------#
            }#end for (h in prof.hours)
            #------------------------------------------------------------------------------#





            #------------------------------------------------------------------------------#
            #      Create diel plots for radiation profile at the bottom and top of the    #
            # canopy.                                                                      #
            #------------------------------------------------------------------------------#
            cat("       ~ Diel plots above and beneath canopy.","\n")

            mn = model$month == m

            #----- Initialise arrays that will hold the mean diel. ------------------------#
            empty     = array( data     = NA
                             , dim      = c(nband,nsimul,nhours,2)
                             , dimnames = list(band$key,simul$suff
                                              ,my.hours,c("top","bottom"))
                             )#end array
            flux      = list( down.beam = empty
                            , down.diff = empty
                            , down.tot  = empty
                            , up.tot    = empty
                            )#end list
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Find the mean diurnal flag for each hour.  We will skip the nighttime   #
            # for relative components, to avoid weird results.                             #
            #------------------------------------------------------------------------------#
            sel.day = tapply(X=model$diel[mn],INDEX=model$hour[mn],FUN=commonest)
            sel.day = sel.day %==% 2
            #------------------------------------------------------------------------------#



            #----- Band loop. -------------------------------------------------------------#
            for (b in sequence(nband)){
               bnd      = band$key [b]
               bnd.desc = band$desc[b]
               prof     = model[[bnd]]


               #----- Find the averages for the hour and month. ---------------------------#
               flux$down.beam[b,,,] = qapply( X      = prof$down.beam[,mn,ptb,drop=FALSE]
                                            , INDEX  = model$hour[mn]
                                            , DIM    = 2
                                            , FUN    = mean
                                            )#end apply
               flux$down.diff[b,,,] = qapply( X      = prof$down.diff[,mn,ptb,drop=FALSE]
                                            , INDEX  = model$hour[mn]
                                            , DIM    = 2
                                            , FUN    = mean
                                            )#end apply
               flux$down.tot [b,,,] = qapply( X      = prof$down.tot [,mn,ptb,drop=FALSE]
                                            , INDEX  = model$hour[mn]
                                            , DIM    = 2
                                            , FUN    = mean
                                            )#end apply
               flux$up.tot   [b,,,] = qapply( X      = prof$up.tot   [,mn,ptb,drop=FALSE]
                                            , INDEX  = model$hour[mn]
                                            , DIM    = 2
                                            , FUN    = mean
                                            )#end apply
               #---------------------------------------------------------------------------#
            }#end for (b in sequence(nband))
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #     Define range for x and y.                                                #
            #------------------------------------------------------------------------------#
            hour.limit = range(my.hours)
            hour.at    = prof.hours
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #      Loop over the components.  Since shortwave and longwave are very        #
            # different, we only plot shortwave components, and we only use the hours that #
            # are "sufficiently diurnal".                                                  #
            #------------------------------------------------------------------------------#
            nflux = length(flux)
            for (f in sequence(nflux)){
               #----- Get this flux component. --------------------------------------------#
               flux.key   = names(flux)[f]
               flux.desc  = c("Downward (direct)","Downward (diffuse)"
                             ,"Downward (total)" ,"Upward")[f]
               if (flux.key %in% c("down.beam","down.diff")){
                  bwhich     = which(! band$key %in% "tir")
                  nbwhich    = length(bwhich)
               }else{
                  bwhich     = sequence(nband)
                  nbwhich    = nband
               }#end if
               this       = flux[[f]]
               this.limit = range(c(this[bwhich,,,]))
               #---------------------------------------------------------------------------#




               #---------------------------------------------------------------------------#
               #      Loop over formats.                                                   #
               #---------------------------------------------------------------------------#
               for (o in sequence(nout)){
                  #----- Open file. -------------------------------------------------------#
                  fichier = file.path(dielpath
                                     ,paste("diel-",iata,"-",flux.key,"-",mm,"-",lai.key
                                           ,".",outform[o],sep="") )
                  if       (outform[o] == "x11"){
                    X11(width=size$width,height=size$height,pointsize=ptsz)
                  }else if (outform[o] == "png"){
                    png(filename=fichier,width=size$width*depth
                       ,height=size$height*depth,pointsize=ptsz,res=depth)
                  }else if (outform[o] == "eps"){
                    postscript(file=fichier,width=size$width,height=size$height
                              ,pointsize=ptsz,paper=size$paper)
                  }else if (outform[o] == "pdf"){
                    pdf(file=fichier,onefile=FALSE,width=size$width
                       ,height=size$height,pointsize=ptsz,paper=size$paper)
                  }#end if
                  #------------------------------------------------------------------------#



                  #----- Set up the plotting window. --------------------------------------#
                  par(par.user)
                  par(oma=c(0,1,4,0))
                  layout( mat     = rbind(lo.simul$mat.off,rep(1,times=lo.simul$ncol))
                        , heights = c(rep(6/lo.simul$nrow,lo.simul$nrow),1)
                        )#end layout
                  #------------------------------------------------------------------------#




                  #------------------------------------------------------------------------#
                  #      First plot: the legend.  No reason to plot thermal IR radiation   #
                  # for direct component.                                                  #
                  #------------------------------------------------------------------------#
                  if (flux.key %in% c("down.beam","down.diff")){
                     #------ Select the bands that should be used. ------------------------#
                     par(mar=c(0.1,4.1,0.1,2.1))
                     plot.new()
                     plot.window(xlim=c(0,1),ylim=c(0,1))
                     legend( x      = "bottom"
                           , inset  = 0.0
                           , legend = c(band$desc[bwhich],"","Top","Bottom")
                           , fill   = c(band$col [bwhich],rep("transparent",3))
                           , border = c(band$col [bwhich],rep("transparent",3))
                           , col    = c(rep("transparent",nband),rep(foreground,2))
                           , lty    = c(rep("solid",nband),"solid","dashed")
                           , lwd    = 2.5
                           , cex    = 0.9 * cex.ptsz
                           , ncol   = 3
                           , xpd    = TRUE
                           )#end legend
                     #---------------------------------------------------------------------#
                  }else{
                     #------ All bands shall be used. -------------------------------------#
                     par(mar=c(0.1,4.1,0.1,2.1))
                     plot.new()
                     plot.window(xlim=c(0,1),ylim=c(0,1))
                     legend( x      = "bottom"
                           , inset  = 0.0
                           , legend = c(band$desc,"Top","Bottom")
                           , fill   = c(band$col ,rep("transparent",2))
                           , border = c(band$col ,rep("transparent",2))
                           , col    = c(rep("transparent",nband),rep(foreground,2))
                           , lty    = c(rep("solid",nband),"solid","dashed")
                           , lwd    = 2.5
                           , cex    = 0.9 * cex.ptsz
                           , ncol   = 3
                           , xpd    = TRUE
                           )#end legend
                     #---------------------------------------------------------------------#
                  }#end if
                  #------------------------------------------------------------------------#




                  #------------------------------------------------------------------------#
                  #      Loop over the simulations, and the bands.                         #
                  #------------------------------------------------------------------------#
                  for (s in sequence(nsimul)){
                     #------ Open plotting device. ----------------------------------------#
                     par(mar=lo.simul$mar0)
                     plot.new()
                     plot.window(xlim = hour.limit,ylim = this.limit)
                     axis(side=1,at=hour.at)
                     axis(side=2,las=1)
                     abline(v=hour.at,h=axTicks(2),col=grid.colour,lty="dotted")
                     #---------------------------------------------------------------------#



                     #----- Loop over bands to be plotted. --------------------------------#
                     for (b in bwhich){
                         lines(x=my.hours,y=this[b,s,,1],col=band$col[b],type="o",pch=16
                              ,lty="solid" ,lwd=2.0)
                         lines(x=my.hours,y=this[b,s,,2],col=band$col[b],type="o",pch=16
                              ,lty="dashed",lwd=2.0)
                     }#end for (b in bwhich)
                     #---------------------------------------------------------------------#


                     #----- Plot the final information. -----------------------------------#
                     box()
                     title(main=simul$desc[s],font.main=1,line=1)
                     #---------------------------------------------------------------------#
                  }#end for (s in sequence(nsimul))
                  #------------------------------------------------------------------------#


                  #------------------------------------------------------------------------#
                  #     Plot the general annotation.                                       #
                  #------------------------------------------------------------------------#
                  lex     = desc.unit(desc="Hour"      ,unit=untab$utc )
                  ley     = desc.unit(desc="Irradiance",unit=untab$wom2)
                  lesub   = desc.unit( desc    = paste("Average LAI"
                                                      ,sprintf("%6.2f",lai.mon[m,pnow])
                                                      )#end paste
                                     , unit    = untab$m2lom2
                                     , bracket = FALSE
                                     )#end desc.unit
                  letitre = paste(desc," - Mean diel of ",flux.desc,"  -  ",month.name[m]
                                 ,sep="")
                  gtitle( main      = letitre
                        , sub       = lesub
                        , xlab      = lex
                        , ylab      = ley
                        , line.xlab = 3.2
                        , line.ylab = 2.6
                        , line.sub  = 4.2
                        , off.sub   = 2/21
                        , cex.main  = 1.1*cex.ptsz
                        , cex.xlab  = 0.9*cex.ptsz
                        , cex.ylab  = 0.9*cex.ptsz
                        , cex.sub   = 0.7*cex.ptsz
                        )#end gtitle
                  #------------------------------------------------------------------------#



                  #----- Close the device. ------------------------------------------------#
                  if (outform[o] == "x11"){
                     locator(n=1)
                     dev.off()
                  }else{
                     dev.off()
                  }#end if
                  clean.tmp()
                  #------------------------------------------------------------------------#
               }#end for (o in sequence(nout))
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     For relative plots, we use only daytime and only solar radiation,     #
               # since thermal radiation scale is completely different.                    #
               #---------------------------------------------------------------------------#
               bwhich     = which(! band$key %in% "tir")
               nbwhich    = length(bwhich)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      We only plot the daytime hours.  In case 00UTC is daytime and this   #
               # is not a polar day, we put the last bits of the day in the end.           #
               #---------------------------------------------------------------------------#
               dayhours = my.hours[sel.day]
               if ( 0 %in% dayhours & any(! sel.day)){
                  #----- Find the local midnight. -----------------------------------------#
                  midnight   = mean(my.hours[!sel.day])
                  dayhours   = c(      dayhours[dayhours > midnight]
                                , 24 + dayhours[dayhours < midnight] )
                  hn         = match(dayhours %% 24 , my.hours)
                  #------------------------------------------------------------------------#
               }#end if
               #---------------------------------------------------------------------------#


               #----- Find the hours to plot, labels and axis. ----------------------------#
               hn         = match(dayhours %% 24 , my.hours)
               day.at     = dayhours
               day.labels = day.at %% 24
               day.limit  = pretty.xylim(u=dayhours,fracexp=0.0,is.log=FALSE)
               #---------------------------------------------------------------------------#



               #----- Discard the nighttime data, and find the range. ---------------------#
               this       = this[,,hn,,drop=FALSE]
               this       = apply(X=this,MARGIN=c(1,2,3),FUN=ediff)
               this.limit = range(c(this[bwhich,,]),finite=TRUE)
               #---------------------------------------------------------------------------#




               #---------------------------------------------------------------------------#
               #      Loop over formats.                                                   #
               #---------------------------------------------------------------------------#
               for (o in sequence(nout)){
                  #----- Open file. -------------------------------------------------------#
                  fichier = file.path(reldielpath
                                     ,paste("rel_diel-",iata,"-",flux.key,"-",mm,"-"
                                           ,lai.key,".",outform[o],sep=""))
                  if       (outform[o] == "x11"){
                    X11(width=size$width,height=size$height,pointsize=ptsz)
                  }else if (outform[o] == "png"){
                    png(filename=fichier,width=size$width*depth
                       ,height=size$height*depth,pointsize=ptsz,res=depth)
                  }else if (outform[o] == "eps"){
                    postscript(file=fichier,width=size$width,height=size$height
                              ,pointsize=ptsz,paper=size$paper)
                  }else if (outform[o] == "pdf"){
                    pdf(file=fichier,onefile=FALSE,width=size$width
                       ,height=size$height,pointsize=ptsz,paper=size$paper)
                  }#end if
                  #------------------------------------------------------------------------#



                  #----- Set up the plotting window. --------------------------------------#
                  par(par.user)
                  par(oma=c(0,1,4,0))
                  layout( mat     = rbind(lo.simul$mat.off,rep(1,times=lo.simul$ncol))
                        , heights = c(rep(6/lo.simul$nrow,lo.simul$nrow),1)
                        )#end layout
                  #------------------------------------------------------------------------#




                  #------------------------------------------------------------------------#
                  #      First plot: the legend.  Skip TIR.                                #
                  #------------------------------------------------------------------------#
                  par(mar=c(0.1,4.1,0.1,2.1))
                  plot.new()
                  plot.window(xlim=c(0,1),ylim=c(0,1))
                  legend( x      = "bottom"
                        , inset  = 0.0
                        , legend = c(band$desc[bwhich])
                        , fill   = c(band$col [bwhich])
                        , border = c(band$col [bwhich])
                        , cex    = 0.9 * cex.ptsz
                        , ncol   = 3
                        , xpd    = TRUE
                        )#end legend
                  #------------------------------------------------------------------------#




                  #------------------------------------------------------------------------#
                  #      Loop over the simulations, and the bands.                         #
                  #------------------------------------------------------------------------#
                  for (s in sequence(nsimul)){
                     #------ Open plotting device. ----------------------------------------#
                     par(mar=lo.simul$mar0)
                     plot.new()
                     plot.window(xlim = day.limit,ylim = this.limit)
                     axis(side=1,at=day.at,labels=day.labels)
                     axis(side=2,las=1)
                     abline(v=hour.at,h=axTicks(2),col=grid.colour,lty="dotted")
                     #---------------------------------------------------------------------#



                     #----- Loop over bands to be plotted. --------------------------------#
                     for (b in bwhich){
                         lines(x=dayhours,y=this[b,s,],col=band$col[b],type="o",pch=16
                              ,lty="solid" ,lwd=2.0)
                     }#end for (b in bwhich)
                     #---------------------------------------------------------------------#


                     #----- Plot the final information. -----------------------------------#
                     box()
                     title(main=simul$desc[s],font.main=1,line=1)
                     #---------------------------------------------------------------------#
                  }#end for (s in sequence(nsimul))
                  #------------------------------------------------------------------------#


                  #------------------------------------------------------------------------#
                  #     Plot the general annotation.                                       #
                  #------------------------------------------------------------------------#
                  lex     = desc.unit(desc="Hour"      ,unit=untab$utc )
                  ley     = desc.unit(desc="Relative Irradiance",unit=untab$empty)
                  lesub   = desc.unit( desc    = paste("Average LAI"
                                                      ,sprintf("%6.2f",lai.mon[m,pnow])
                                                      )#end paste
                                     , unit    = untab$m2lom2
                                     , bracket = FALSE
                                     )#end desc.unit
                  letitre = paste(desc," - Relative irradiance beneath the canopy"
                                 ,"\n",flux.desc,"  -  ",month.name[m],sep="")
                  gtitle( main      = letitre
                        , sub       = lesub
                        , xlab      = lex
                        , ylab      = ley
                        , line.xlab = 3.2
                        , line.ylab = 2.6
                        , line.sub  = 4.2
                        , off.sub   = 2/21
                        , cex.main  = 1.1*cex.ptsz
                        , cex.xlab  = 0.9*cex.ptsz
                        , cex.ylab  = 0.9*cex.ptsz
                        , cex.sub   = 0.7*cex.ptsz
                        )#end gtitle
                  #------------------------------------------------------------------------#



                  #----- Close the device. ------------------------------------------------#
                  if (outform[o] == "x11"){
                     locator(n=1)
                     dev.off()
                  }else{
                     dev.off()
                  }#end if
                  clean.tmp()
                  #------------------------------------------------------------------------#
               }#end for (o in sequence(nout))
               #---------------------------------------------------------------------------#

            }#end for (f in sequence(nflux))
            #------------------------------------------------------------------------------#
         }#end for (m in sequence(12))
         #---------------------------------------------------------------------------------#









         #---------------------------------------------------------------------------------#
         #     Profile by zenith angle class.                                              #
         #---------------------------------------------------------------------------------#
         for (z in loop.zen){
            cat("     > ",zen.desc[z],"\n")
            if (iata %in% names(obs.light)){
               sel      = ( obs.light[[iata]]$highsun 
                          & is.finite(obs.light[[iata]]$par.tot.top)
                          )#end sel
               obs.mon  = unique(floor(obs.light[[iata]]$doy[sel]))
               
               use      = floor(dayofyear(when=model$when)) %in% obs.mon
               zc       = use & model$zen.class %in% z
            }else{
               zc       = model$zen.class %in% z
            }#end if
            if (any(zc)){


               #------ Create path for this hour. -----------------------------------------#
               profzenpath  = file.path(profilepath,"zenith")
               if (! file.exists(profzenpath)) dir.create(profzenpath)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Load LAI/WAI/TAI profiles.                                            #
               #---------------------------------------------------------------------------#
               lai     = apply(X=model$lai[,zc,pp,drop=FALSE],MARGIN=c(1,3),FUN=mean)
               wai     = apply(X=model$wai[,zc,pp,drop=FALSE],MARGIN=c(1,3),FUN=mean)
               tai     = apply(X=model$tai[,zc,pp,drop=FALSE],MARGIN=c(1,3),FUN=mean)
               cum.lai = t(apply(X=lai,MARGIN=1,FUN=cumsum))
               cum.wai = t(apply(X=wai,MARGIN=1,FUN=cumsum))
               cum.tai = t(apply(X=tai,MARGIN=1,FUN=cumsum))
               #---------------------------------------------------------------------------#


               #----- Fix the limits for the Y axis. --------------------------------------#
               tai.lim = range(-cum.tai)
               tai.at   = pretty(tai.lim)
               tai.lab  = -tai.at
               #---------------------------------------------------------------------------#



               #----- Band loop. ----------------------------------------------------------#
               for (b in sequence(nband)){
                  bnd      = band$key [b]
                  bnd.desc = band$desc[b]
                  prof     = model[[bnd]]


                  #------------------------------------------------------------------------#
                  #     Shift the results a bit so we can see when two curves are on top   #
                  # of each other.                                                         #
                  #------------------------------------------------------------------------#
                  if ( ! bnd %in% "tir"){
                     #----- Scale all variables. ------------------------------------------#
                     down.tot  = 100. * apply ( X        = prof$down.tot [,zc,pp,drop=F]
                                              , MARGIN   = c(1,2)
                                              , FUN      = scal.first
                                              )#end apply
                     down.beam = 100. * apply ( X        = prof$down.beam[,zc,pp,drop=F]
                                              , MARGIN   = c(1,2)
                                              , FUN      = scal.first
                                              )#end apply
                     down.diff = 100. * apply ( X        = prof$down.diff[,zc,pp,drop=F]
                                              , MARGIN   = c(1,2)
                                              , FUN      = scal.first
                                              )#end apply
                     up.tot    = 100. * apply ( X        = prof$up.tot   [,zc,pp,drop=F]
                                              , MARGIN   = c(1,2)
                                              , FUN      = scal.first
                                              )#end apply

                     abs.leaf  = prof$abs.leaf[,zc,pp,drop=FALSE]
                     aidx      = seq_along(abs.leaf)
                     abs.leaf  = split( x = abs.leaf
                                      , f = list( arrayInd(aidx,dim(abs.leaf))[,1]
                                                , arrayInd(aidx,dim(abs.leaf))[,2]
                                                )#end list
                                      )#end split
                     abs.tot  = prof$abs.tot[,zc,pp,drop=FALSE]
                     aidx     = seq_along(abs.tot)
                     abs.tot  = split( x = abs.tot
                                      , f = list( arrayInd(aidx,dim(abs.tot))[,1]
                                                , arrayInd(aidx,dim(abs.tot))[,2]
                                                )#end list
                                      )#end split


                     abs.leaf  = 100. * mapply( FUN      = scal.last
                                              , x        = abs.leaf
                                              , y        = abs.tot
                                              , SIMPLIFY = TRUE
                                              )#end mapply
                     abs.tot   = 100. * apply ( X      = prof$abs.tot  [,zc,pp,drop=F]
                                              , MARGIN = c(1,2)
                                              , FUN    = scal.last
                                              )#end apply
                     abs.leaf  = c(abs.leaf) + 0 * abs.tot
                     #---------------------------------------------------------------------#


                     #----- Organise the mess. --------------------------------------------#
                     down.tot  = aperm(a=down.tot ,perm=c(2,3,1))
                     down.beam = aperm(a=down.beam,perm=c(2,3,1))
                     down.diff = aperm(a=down.diff,perm=c(2,3,1))
                     up.tot    = aperm(a=up.tot   ,perm=c(2,3,1))
                     abs.tot   = aperm(a=abs.tot  ,perm=c(2,3,1))
                     abs.leaf  = aperm(a=abs.leaf ,perm=c(2,3,1))
                     #---------------------------------------------------------------------#
                  }else{
                     down.tot  = prof$down.tot [,zc,pp,drop=F]
                     down.beam = prof$down.beam[,zc,pp,drop=F]
                     down.diff = prof$down.diff[,zc,pp,drop=F]
                     up.tot    = prof$up.tot   [,zc,pp,drop=F]
                     abs.tot   = prof$abs.tot  [,zc,pp,drop=F]
                     abs.leaf  = prof$abs.leaf [,zc,pp,drop=F]
                  }#end if
                  #------------------------------------------------------------------------#



                  #----- Find the averages for the hour and month. ------------------------#
                  if (tolower(ci.estimate) %in% "six.summary"){
                     down.beam = apply(X=down.beam,MARGIN=c(1,3),FUN=six.summary)
                     down.diff = apply(X=down.diff,MARGIN=c(1,3),FUN=six.summary)
                     down.tot  = apply(X=down.tot ,MARGIN=c(1,3),FUN=six.summary)
                     up.tot    = apply(X=up.tot   ,MARGIN=c(1,3),FUN=six.summary)
                     abs.tot   = apply(X=abs.tot  ,MARGIN=c(1,3),FUN=six.summary)
                     abs.leaf  = apply(X=abs.leaf ,MARGIN=c(1,3),FUN=six.summary)
                  }else{
                     down.beam = apply(X=down.beam,MARGIN=c(1,3),FUN=boot.six.summary,R=n.boot)
                     down.diff = apply(X=down.diff,MARGIN=c(1,3),FUN=boot.six.summary,R=n.boot)
                     down.tot  = apply(X=down.tot ,MARGIN=c(1,3),FUN=boot.six.summary,R=n.boot)
                     up.tot    = apply(X=up.tot   ,MARGIN=c(1,3),FUN=boot.six.summary,R=n.boot)
                     abs.tot   = apply(X=abs.tot  ,MARGIN=c(1,3),FUN=boot.six.summary,R=n.boot)
                     abs.leaf  = apply(X=abs.leaf ,MARGIN=c(1,3),FUN=boot.six.summary,R=n.boot)
                  }#end if
                  #------------------------------------------------------------------------#




                  #------------------------------------------------------------------------#
                  #     Keep only mean, q025, and q975.                                    #
                  #------------------------------------------------------------------------#
                  down.tot  = down.tot [c(1,5,6),,]
                  down.beam = down.beam[c(1,5,6),,]
                  down.diff = down.diff[c(1,5,6),,]
                  up.tot    = up.tot   [c(1,5,6),,]
                  abs.tot   = abs.tot  [c(1,5,6),,]
                  abs.leaf  = abs.leaf [c(1,5,6),,]
                  #------------------------------------------------------------------------#



                  #------------------------------------------------------------------------#
                  #     Shift the results a bit so we can see when two curves are on top   #
                  # of each other.                                                         #
                  #------------------------------------------------------------------------#
                  if ( ! bnd %in% "tir"){
                     down.tot  = apply(X=down.tot ,c(1,3),"^",simul$mult.down     )
                     down.beam = apply(X=down.beam,c(1,3),"^",simul$mult.beam.down)
                     down.diff = apply(X=down.diff,c(1,3),"^",simul$mult.diff.down)
                     up.tot    = apply(X=up.tot   ,c(1,3),"^",simul$mult.up       )
                     down.tot  = aperm(a=down.tot ,perm=c(2,1,3))
                     down.beam = aperm(a=down.beam,perm=c(2,1,3))
                     down.diff = aperm(a=down.diff,perm=c(2,1,3))
                     up.tot    = aperm(a=up.tot   ,perm=c(2,1,3))
                  }#end if
                  #------------------------------------------------------------------------#




                  #------------------------------------------------------------------------#
                  #     Get data for this site if the site has data.                       #
                  #------------------------------------------------------------------------#
                  is.iata = iata %in% names(obs.light)
                  obs = list()
                  if (is.iata){
                     is.tot  = bnd %in% names(obs.light[[iata]]$zmean)
                     is.diff = bnd %in% names(obs.light[[iata]]$mmean)
                     obs$lai = obs.light[[iata]]$cum.lai
                     zmean   = obs.light[[iata]]$zmean[[bnd]]
                     mmean   = obs.light[[iata]]$mmean[[bnd]]
                  }else{
                     is.tot  = FALSE
                     is.diff = FALSE
                     obs$lai = NA
                     zmean   = list()
                     mmean   = list()
                  }#end if
                  #------ Build an standardised summary. ----------------------------------#
                  stats  = c("expected","ci.lower","ci.upper")
                  k      = match(stats,six.summary.names)
                  nstats = length(stats)
                  nlai   = length(obs$lai)
                  empty  = matrix(data=NA,nrow=nstats,ncol=nlai,dimnames=list(stats,NULL))
                  #------------------------------------------------------------------------#



                  #------ Initialise lists without data. ----------------------------------#
                  obs$down.tot  = empty
                  obs$down.beam = empty
                  obs$down.diff = empty
                  obs$up.tot    = empty
                  obs$abs.tot   = empty
                  obs$abs.leaf  = empty
                  if ("down.tot"  %in% names(zmean)) obs$down.tot  = zmean$down.tot [k,z,]
                  if ("up.tot"    %in% names(zmean)) obs$up.tot    = zmean$up.tot   [k,z,]
                  if ("down.diff" %in% names(mmean)) obs$down.diff = mmean$down.diff[k,  ]
                  if ("abs.tot"   %in% names(zmean)) obs$abs.tot   = zmean$abs.tot  [k,z,]
                  #------------------------------------------------------------------------#




                  #------------------------------------------------------------------------#
                  #     Skip plot if it is solar radiation and it is nighttime.            #
                  #------------------------------------------------------------------------#
                  if (max(down.tot) > 1.){


                     #----- Near infrarred cannot be in log scale. ------------------------#
                     is.log = ! bnd %in% "tir"
                     if (is.log){xylog = "x"}else{xylog=""}
                     #---------------------------------------------------------------------#


                     #----- Find limits. --------------------------------------------------#
                     down.tot.lim  = pretty.xylim(u=down.tot           ,is.log=is.log)
                     down.diff.lim = pretty.xylim(u=down.diff          ,is.log=is.log)
                     up.tot.lim    = pretty.xylim(u=up.tot             ,is.log=is.log)
                     abs.tot.lim   = pretty.xylim(u=c(abs.leaf,abs.tot),is.log=is.log)
                     #---------------------------------------------------------------------#


                     #----- Find axis ticks. ----------------------------------------------#
                     if (is.log){
                        down.tot.at   = pretty.log(down.tot.lim )
                        down.diff.at  = pretty.log(down.diff.lim)
                        up.tot.at     = pretty.log(up.tot.lim   )
                        abs.tot.at    = pretty.log(abs.tot.lim  )
                     }else{
                        down.tot.at   = pretty(down.tot.lim )
                        down.diff.at  = pretty(down.diff.lim)
                        up.tot.at     = pretty(up.tot.lim   )
                        abs.tot.at    = pretty(abs.tot.lim  )
                     }#end if
                     #---------------------------------------------------------------------#


                     for (o in sequence(nout)){
                        #----- Open file. -------------------------------------------------#
                        fichier = file.path(profzenpath
                                           ,paste("profile-",iata,"-",bnd,"-",zen.key[z]
                                                 ,"-",lai.key,".",outform[o],sep=""))
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
                        #------------------------------------------------------------------#


                        #----- Set up the plotting window. --------------------------------#
                        par(par.user)
                        par(oma=c(0,1,4,0))
                        if (bnd %in% "tir"){
                           layout(mat=rbind(c(2,3,4),c(1,1,1)),heights=c(6,1))
                        }else{
                           layout(mat=rbind(c(2,4),c(3,5),c(1,1)),heights=c(3,3,1))
                        }#end if
                        #------------------------------------------------------------------#




                        #----- First plot: the legend. ------------------------------------#
                        par(mar=c(0.1,4.1,0.1,2.1))
                        plot.new()
                        plot.window(xlim=c(0,1),ylim=c(0,1))
                        legend( x       = "bottom"
                              , inset   = 0.0
                              , legend  = c(simul$desc,"ABRACOS","")
                              , fill    = c(simul$bgcol ,grey.bg,"transparent")
                              , border  = c(simul$bgcol ,grey.bg,"transparent")
                              , col     = c(simul$colour,grey.fg,"transparent")
                              , density = c(simul$density,25,25)
                              , angle   = c(simul$angle  , 0, 0)
                              , lty     = "solid"
                              , lwd     = 2.5
                              , cex     = 0.9 * cex.ptsz
                              , ncol    = 3
                              , xpd     = TRUE
                              )#end legend
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #      Direct radiation is not plotted for TIR.                    #
                        #------------------------------------------------------------------#
                        if (! bnd %in% "tir"){
                           #------ Second plot: Downward direct. --------------------------#
                           par(mar=c(3.1,3.1,3.1,1.1))
                           plot.new()
                           plot.window(xlim=down.tot.lim,ylim=tai.lim,log=xylog)
                           axis(side=1,at=down.tot.at,cex.axis=0.8,padj=-0.8)
                           axis(side=2,at=tai.at,labels=tai.lab,las=1)
                           abline(v=down.tot.at,h=tai.at,col=grid.colour,lty="solid")
                           box()
                           title(main="Downward - Total",font.main=1,line=1)
                           #----- Confidence intervals. -----------------------------------#
                           for (s in sequence(nsimul)){
                              epolygon( x = c(down.tot[2,s,],rev(down.tot[3,s,]))
                                      , y = c(-cum.tai[  s,],rev(-cum.tai[  s,]))
                                      , density = simul$density[s]
                                      , angle   = simul$angle  [s]
                                      , col     = simul$bgcol  [s]
                                      )#end epolygon
                           }#end for
                           epolygon( x = c( obs$down.tot[2,],rev( obs$down.tot[3,]))
                                   , y = c(-obs$lai         ,rev(-obs$lai         ))
                                   , density = 25
                                   , angle   = 0
                                   , col     = grey.bg
                                   )#end epolygon
                           #----- Expected values. ----------------------------------------#
                           for (s in sequence(nsimul)){
                              lines(x=down.tot[1,s,],y=-cum.tai[s,]
                                   ,col=simul$colour[s],lwd=2.5,lty="solid")
                           }#end for
                           lines(x=obs$down.tot[1,],y=-obs$lai,type="o",pch=16,col=grey.fg)
                           #---------------------------------------------------------------#
                        }#end if
                        #------------------------------------------------------------------#




                        #------ Third plot: Downward diffuse. -----------------------------#
                        par(mar=c(4.1,3.1,2.1,1.1))
                        plot.new()
                        plot.window(xlim=down.diff.lim,ylim=tai.lim,log=xylog)
                        axis(side=1,at=down.diff.at,cex.axis=0.8,padj=-0.8)
                        axis(side=2,at=tai.at,labels=tai.lab,las=1)
                        abline(v=down.diff.at,h=tai.at,col=grid.colour,lty="solid")

                        box()
                        title(main="Downward - Diffuse",font.main=1,line=1)
                        #----- Confidence intervals. --------------------------------------#
                        for (s in sequence(nsimul)){
                           epolygon( x = c(down.diff[2,s,],rev(down.diff[3,s,]))
                                   , y = c(-cum.tai [  s,],rev(-cum.tai [  s,]))
                                   , density = simul$density[s]
                                   , angle   = simul$angle  [s]
                                   , col     = simul$bgcol  [s]
                                   )#end epolygon
                        }#end for
                        epolygon( x = c( obs$down.diff[2,],rev( obs$down.diff[3,]))
                                , y = c(-obs$lai          ,rev(-obs$lai          ))
                                , density = 25
                                , angle   = 0
                                , col     = grey.bg
                                )#end epolygon
                        #----- Expected values. -------------------------------------------#
                        for (s in sequence(nsimul)){
                           lines(x=down.diff[1,s,],y=-cum.tai[s,]
                                ,col=simul$colour[s],lwd=2.5,lty="solid")
                        }#end for
                        lines(x=obs$down.diff[1,],y=-obs$lai,type="o",pch=16,col=grey.fg)
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #     Fourth plot: Upward radiation.  Margins depend on the band.  #
                        #------------------------------------------------------------------#
                        if (bnd %in% "tir"){
                           par(mar=c(4.1,3.1,2.1,1.1))
                        }else{
                           par(mar=c(3.1,2.1,3.1,2.1))
                        }#end if
                        plot.new()
                        plot.window(xlim=up.tot.lim,ylim=tai.lim,log=xylog)
                        axis(side=1,at=up.tot.at,cex.axis=0.8,padj=-0.8)
                        axis(side=2,at=tai.at,labels=tai.lab,las=1)
                        abline(v=up.tot.at,h=tai.at,col=grid.colour,lty="solid")
                        box()
                        title(main="Upward",font.main=1,line=1)
                        #----- Confidence intervals. --------------------------------------#
                        for (s in sequence(nsimul)){
                           epolygon( x = c(up.tot   [2,s,],rev(up.tot   [3,s,]))
                                   , y = c(-cum.tai [  s,],rev(-cum.tai [  s,]))
                                   , density = simul$density[s]
                                   , angle   = simul$angle  [s]
                                   , col     = simul$bgcol  [s]
                                   )#end epolygon
                        }#end for
                        epolygon( x = c( obs$up.tot   [2,],rev( obs$up.tot   [3,]))
                                , y = c(-obs$lai          ,rev(-obs$lai          ))
                                , density = 25
                                , angle   = 0
                                , col     = grey.bg
                                )#end epolygon
                        #----- Expected values. -------------------------------------------#
                        for (s in sequence(nsimul)){
                           lines(x=up.tot[1,s,],y=-cum.tai[s,]
                                ,col=simul$colour[s],lwd=2.5,lty="solid")
                        }#end for
                        lines(x=obs$up.tot[1,],y=-obs$lai,type="o",pch=16,col=grey.fg)
                        #------------------------------------------------------------------#



                        #------ Fifth plot: Absorbed radiation. ---------------------------#
                        par(mar=c(4.1,2.1,2.1,2.1))
                        plot.new()
                        plot.window(xlim=abs.tot.lim,ylim=tai.lim,log=xylog)
                        axis(side=1,at=abs.tot.at,cex.axis=0.8,padj=-0.8)
                        axis(side=2,at=tai.at,labels=tai.lab,las=1)
                        abline(v=abs.tot.at,h=tai.at,col=grid.colour,lty="solid")
                        box()
                        title(main="Cumulative absorption (dashed = leaves)"
                             ,font.main=1,line=1)
                        #----- Confidence intervals. --------------------------------------#
                        for (s in sequence(nsimul)){
                           epolygon( x       = c(abs.tot  [2,s,],rev(abs.tot  [3,s,]))
                                   , y       = c(-cum.tai [  s,],rev(-cum.tai [  s,]))
                                   , density = simul$density[s]
                                   , angle   = simul$angle  [s]
                                   , col     = simul$bgcol  [s]
                                   )#end epolygon
                        }#end for
                        epolygon( x       = c( obs$abs.tot  [2,],rev( obs$abs.tot  [3,]))
                                , y       = c(-obs$lai          ,rev(-obs$lai          ))
                                , density = 25
                                , angle   = 0
                                , col     = grey.bg
                                )#end epolygon
                        #----- Expected values. -------------------------------------------#
                        for (s in sequence(nsimul)){
                           lines(x=abs.tot [1,s,],y=-cum.tai[s,]
                                ,col=simul$colour[s],lwd=2.5,lty="solid")
                        }#end for
                        for (s in sequence(nsimul)){
                           lines(x=abs.leaf[1,s,],y=-cum.tai[s,]
                                ,col=simul$colour[s],lwd=2.0,lty="dotdash")
                        }#end for
                        lines(x=obs$abs.tot[1,],y=-obs$lai,type="o",pch=16,col=grey.fg)
                        #------------------------------------------------------------------#



                        #------------------------------------------------------------------#
                        #     Plot the general annotation.                                 #
                        #------------------------------------------------------------------#
                        if (bnd %in% "tir"){
                           lex     = desc.unit( desc = paste(band$desc[b],"Irradiance")
                                              , unit = untab[[band$unit[b]]]
                                              )#end desc.unit
                        }else{
                           lex     = desc.unit( desc = paste(band$desc[b]
                                                            ,"Relative Irradiance")
                                              , unit = untab$empty
                                              )#end desc.unit
                        }#end if
                        ley     = desc.unit(desc="Cumulative TAI",unit=untab$m2lom2)
                        letitre = paste(desc," - Canopy profiles","\n",zen.desc[z],sep="")
                        lesub   = desc.unit( desc = paste("Average LAI"
                                                         ,sprintf("%6.2f",lai.zen[z,pnow])
                                                         )#end paste
                                           , unit = untab$m2lom2
                                           , bracket = FALSE
                                           )#end desc.unit
                        gtitle( main      = letitre
                              , sub       = lesub
                              , xlab      = lex
                              , ylab      = ley
                              , line.xlab = 3.2
                              , line.ylab = 2.6
                              , line.sub  = 4.2
                              , off.sub   = 2/21
                              , cex.main  = 1.1*cex.ptsz
                              , cex.xlab  = 0.9*cex.ptsz
                              , cex.ylab  = 0.9*cex.ptsz
                              , cex.sub   = 0.7*cex.ptsz
                              )#end gtitle
                        #------------------------------------------------------------------#



                        #----- Close the device. ------------------------------------------#
                        if (outform[o] == "x11"){
                           locator(n=1)
                           dev.off()
                        }else{
                           dev.off()
                        }#end if
                        clean.tmp()
                        #------------------------------------------------------------------#
                     }#end for (o in sequence(nout))
                     #---------------------------------------------------------------------#
                  }#end if (max(tot.down) > 1.)
                  #------------------------------------------------------------------------#


                  rm(down.tot,down.beam,down.diff,up.tot,abs.tot,abs.leaf,obs)

               }#end for (b in sequence(nband)){
               #---------------------------------------------------------------------------#
            }#end if (any(zc))
         }#end for (z in sequence(nzen.prof))
         #---------------------------------------------------------------------------------#









         #---------------------------------------------------------------------------------#
         #     Profile of diffuse radiation for gloomy periods.                            #
         #---------------------------------------------------------------------------------#
         cat("     > Gloomy times...","\n")
         if (iata %in% names(obs.light)){
            sel      = ( obs.light[[iata]]$highsun 
                       & is.finite(obs.light[[iata]]$par.tot.top)
                       )#end sel
            obs.mon  = unique(floor(obs.light[[iata]]$doy[sel]))
            
            use      = floor(dayofyear(when=model$when)) %in% obs.mon
            gm       = use & model$gloom
         }else{
            gm       = model$gloom
         }#end if
         gm = model$gloom
         if (any(gm) & plot.gloom){



            #------ Create path for gloomy periods. ---------------------------------------#
            profgloompath  = file.path(profilepath,"gloom")
            if (! file.exists(profgloompath)) dir.create(profgloompath)
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Load LAI/WAI/TAI profiles.                                               #
            #------------------------------------------------------------------------------#
            lai     = apply(X=model$lai[,gm,pp,drop=FALSE],MARGIN=c(1,3),FUN=mean)
            wai     = apply(X=model$wai[,gm,pp,drop=FALSE],MARGIN=c(1,3),FUN=mean)
            tai     = apply(X=model$tai[,gm,pp,drop=FALSE],MARGIN=c(1,3),FUN=mean)
            cum.lai = t(apply(X=lai,MARGIN=1,FUN=cumsum))
            cum.wai = t(apply(X=wai,MARGIN=1,FUN=cumsum))
            cum.tai = t(apply(X=tai,MARGIN=1,FUN=cumsum))
            #------------------------------------------------------------------------------#



            #----- Fix the limits for the Y axis. -----------------------------------------#
            tai.lim = range(-cum.tai)
            tai.at   = pretty(tai.lim)
            tai.lab  = -tai.at
            #------------------------------------------------------------------------------#


            #----- Select the statistics to keep. -----------------------------------------#
            stats          = c("expected","ci.lower","ci.upper")
            k              = match(stats,six.summary.names)
            nstats         = length(stats)
            #------------------------------------------------------------------------------#




            #----- Build the place holder array. ------------------------------------------#
            ncohort       = length(pp)
            mod.down.diff = array( data = NA
                                , dim  = c(nband,nstats,nsimul,ncohort)
                                , dimnames = list( band.key,stats,simul$key,NULL)
                                )#end array
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #     Get data for this site if the site has data.                             #
            #------------------------------------------------------------------------------#
            is.iata = iata %in% names(obs.light)
            obs = list()
            if (is.iata){
               obs$lai = obs.light[[iata]]$cum.lai
               mmean   = obs.light[[iata]]$mmean[[bnd]]
            }else{
               obs$lai = NA
               mmean   = list()
            }#end if
            #------ Build an standardised summary. ----------------------------------------#
            nlai           = length(obs$lai)
            obs$down.diff  = array( data     = NA
                                  , dim      = c(nband,nstats,nlai)
                                  , dimnames = list(band$key,stats,NULL)
                                  )#end array
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Loop over all bands and grab data.                                      #
            #------------------------------------------------------------------------------#
            for (b in sequence(nband)){
               #----- Grab information from this band. ------------------------------------#
               bnd      = band$key [b]
               bnd.desc = band$desc[b]
               prof     = model[[bnd]]
               #---------------------------------------------------------------------------#

               #---------------------------------------------------------------------------#
               #     Scale the radiation by the top layer.
               #---------------------------------------------------------------------------#
               if ( ! bnd %in% "tir"){
                  down.diff = 100. * apply ( X        = prof$down.diff[,gm,pp,drop=F]
                                           , MARGIN   = c(1,2)
                                           , FUN      = scal.first
                                           )#end apply


                  #----- Organise the mess. -----------------------------------------------#
                  down.diff = aperm(a=down.diff,perm=c(2,3,1))
                  #------------------------------------------------------------------------#
               }else{
                  down.diff = prof$down.diff[,zc,pp,drop=F]
               }#end if
               #---------------------------------------------------------------------------#



               #----- Find the averages for the hour and month. ---------------------------#
               if (tolower(ci.estimate) %in% "six.summary"){
                  down.diff = apply(X=down.diff,MARGIN=c(1,3),FUN=six.summary)
               }else{
                  down.diff = apply(X=down.diff,MARGIN=c(1,3),FUN=boot.six.summary,R=n.boot)
               }#end if
               #---------------------------------------------------------------------------#




               #---------------------------------------------------------------------------#
               #     Keep just mean, and confidence intervals.                             #
               #---------------------------------------------------------------------------#
               k         = match(c("expected","ci.lower","ci.upper"),six.summary.names)
               down.diff = down.diff[k,,]
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Shift the results a bit so we can see when two curves are on top      #
               # of each other.                                                            #
               #---------------------------------------------------------------------------#
               if ( ! bnd %in% "tir"){
                  down.diff = apply(X=down.diff,c(1,3),"^",simul$mult.diff.down)
                  down.diff = aperm(a=down.diff,perm=c(2,1,3))
               }#end if
               #---------------------------------------------------------------------------#



               #------ Check whether there are observed values. ---------------------------#
               if ("down.diff" %in% names(mmean)) obs$down.diff[b,,] = mmean$down.diff[k,]
               #---------------------------------------------------------------------------#




               #------ Save the data for this band. ---------------------------------------#
               mod.down.diff[b,,,] = down.diff
               #---------------------------------------------------------------------------#
            }#end for (b in sequence(nband))
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #     Plot the profiles.                                                       #
            #------------------------------------------------------------------------------#
            lo.band = pretty.box(n=nband)
            for (o in sequence(nout)){
               #----- Open file. ----------------------------------------------------------#
               fichier = file.path(profgloompath
                                  ,paste("profile-",iata,"-diffuse-gloomy"
                                        ,"-",lai.key,".",outform[o],sep=""))
               if       (outform[o] == "x11"){
                 X11(width=size$width,height=size$height,pointsize=ptsz)
               }else if (outform[o] == "png"){
                 png(filename=fichier,width=size$width*depth
                    ,height=size$height*depth,pointsize=ptsz,res=depth)
               }else if (outform[o] == "eps"){
                 postscript(file=fichier,width=size$width,height=size$height
                           ,pointsize=ptsz,paper=size$paper)
               }else if (outform[o] == "pdf"){
                 pdf(file=fichier,onefile=FALSE,width=size$width
                    ,height=size$height,pointsize=ptsz,paper=size$paper)
               }#end if
               #---------------------------------------------------------------------------#


               #----- Set up the plotting window. -----------------------------------------#
               par(par.user)
               par(oma=c(0,1,4,0))
               layout( mat     = rbind(lo.band$mat.off,rep(1,lo.band$ncol))
                     , heights = c(rep(6/lo.band$nrow,lo.band$nrow),1)
                     )#end layout
               #---------------------------------------------------------------------------#

               #------ Force logarithmic scale. -------------------------------------------#
               xylog  = "x"
               is.log = TRUE
               #---------------------------------------------------------------------------#




               #----- Find limits. --------------------------------------------------------#
               down.diff.lim = pretty.xylim(u=c(mod.down.diff,obs.down.diff),is.log=is.log)
               down.diff.at  = pretty.log(down.diff.lim)
               #---------------------------------------------------------------------------#




               #----- First plot: the legend. ---------------------------------------------#
               par(mar=c(0.1,4.1,0.1,2.1))
               plot.new()
               plot.window(xlim=c(0,1),ylim=c(0,1))
               legend( x       = "bottom"
                     , inset   = 0.0
                     , legend  = c(simul$desc,"ABRACOS","")
                     , fill    = c(simul$bgcol ,grey.bg,"transparent")
                     , border  = c(simul$bgcol ,grey.bg,"transparent")
                     , col     = c(simul$colour,grey.fg,"transparent")
                     , density = c(simul$density,25,25)
                     , angle   = c(simul$angle  , 0, 0)
                     , lty     = "solid"
                     , lwd     = 2.5
                     , cex     = 0.9 * cex.ptsz
                     , ncol    = 3
                     , xpd     = TRUE
                     )#end legend
               #---------------------------------------------------------------------------#






               #---------------------------------------------------------------------------#
               #      Band loop.                                                           #
               #---------------------------------------------------------------------------#
               for (b in sequence(nband)){
                  #------ Third plot: Downward diffuse. -----------------------------------#
                  par(mar=c(4.1,3.1,2.1,1.1))
                  plot.new()
                  plot.window(xlim=down.diff.lim,ylim=tai.lim,log=xylog)
                  axis(side=1,at=down.diff.at,cex.axis=0.8,padj=-0.8)
                  axis(side=2,at=tai.at,labels=tai.lab,las=1)
                  abline(v=down.diff.at,h=tai.at,col=grid.colour,lty="solid")

                  box()
                  title(main=band$desc[b],font.main=1,line=1)
                  #----- Confidence intervals. --------------------------------------------#
                  for (s in sequence(nsimul)){
                     epolygon( x = c(mod.down.diff[b,2,s,],rev(mod.down.diff[b,3,s,]))
                             , y = c(-cum.tai [s,],rev(-cum.tai [s,]))
                             , density = simul$density[s]
                             , angle   = simul$angle  [s]
                             , col     = simul$bgcol  [s]
                             )#end epolygon
                  }#end for
                  epolygon( x = c( obs$down.diff[b,2,],rev( obs$down.diff[b,3,]))
                          , y = c(-obs$lai,rev(-obs$lai))
                          , density = 25
                          , angle   = 0
                          , col     = grey.bg
                          )#end epolygon
                  #----- Expected values. -------------------------------------------------#
                  for (s in sequence(nsimul)){
                     lines(x=mod.down.diff[b,1,s,],y=-cum.tai[s,]
                          ,col=simul$colour[s],lwd=2.5,lty="solid")
                  }#end for
                  lines(x=obs$down.diff[b,1,],y=-obs$lai,type="o",pch=16,col=grey.fg)
                  #------------------------------------------------------------------------#
               }#end for (b in sequence(nband))
               #---------------------------------------------------------------------------#

               #---------------------------------------------------------------------------#
               #     Plot the general annotation.                                          #
               #---------------------------------------------------------------------------#
               lex     = desc.unit(desc="Relative Irradiance",unit=untab$empty)
               ley     = desc.unit(desc="Cumulative TAI",unit=untab$m2lom2)
               letitre = paste(desc," - Canopy profiles (diffuse)",sep="")
               lesub   = desc.unit( desc = paste("Average LAI"
                                                ,sprintf("%6.2f",lai.zen[z,pnow])
                                                )#end paste
                                  , unit = untab$m2lom2
                                  , bracket = FALSE
                                  )#end desc.unit
               gtitle( main      = letitre
                     , sub       = lesub
                     , xlab      = lex
                     , ylab      = ley
                     , line.xlab = 3.2
                     , line.ylab = 2.6
                     , line.sub  = 4.2
                     , off.sub   = 2/21
                     , cex.main  = 1.1*cex.ptsz
                     , cex.xlab  = 0.9*cex.ptsz
                     , cex.ylab  = 0.9*cex.ptsz
                     , cex.sub   = 0.7*cex.ptsz
                     )#end gtitle
               #---------------------------------------------------------------------------#



               #----- Close the device. ---------------------------------------------------#
               if (outform[o] == "x11"){
                  locator(n=1)
                  dev.off()
               }else{
                  dev.off()
               }#end if
               clean.tmp()
               #---------------------------------------------------------------------------#
            }#end for (o in sequence(nout))
            #------------------------------------------------------------------------------#
            rm(mod.down.diff,obs)
         }#end if (any(gm))
         #---------------------------------------------------------------------------------#



      }#end for (l in seq_along(model$lai.plot)
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #     Plot the total extinction of light as a function of cumulative LAI and month.  #
      #------------------------------------------------------------------------------------#
      cat("     > Relative light beneath canopy as a function of patch LAI","\n")





      #------------------------------------------------------------------------------------#
      #     Load LAI and flux for this site.                                               #
      #------------------------------------------------------------------------------------#
      tot.lai = summ[[iata]]$lai
      flux    = summ[[iata]]$flux
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Define range for x and y.                                                      #
      #------------------------------------------------------------------------------------#
      lai.limit = range(tot.lai)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     We display only solar radiation, since thermal radiation scale is completely   #
      # different and has even less information than solar radiation.                      #
      #------------------------------------------------------------------------------------#
      bwhich     = which(! band$key %in% "tir")
      nbwhich    = length(bwhich)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Irradiance component.                                                          #
      #------------------------------------------------------------------------------------#
      nflux = length(flux)
      for (f in sequence(nflux)){
         #----- Get this flux component. --------------------------------------------------#
         flux.key   = names(flux)[f]
         flux.desc  = c("Downward (direct)","Downward (diffuse)"
                       ,"Downward (total)" ,"Upward")[f]
         this       = flux[[f]]
         this.limit = range(c(this[bwhich,,,]),finite=TRUE)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Loop over formats.                                                         #
         #---------------------------------------------------------------------------------#
         for (o in sequence(nout)){
            #----- Open file. -------------------------------------------------------------#
            fichier = file.path(laipath
                               ,paste("beneath_relative-",iata,"-",flux.key
                                     ,".",outform[o],sep=""))
            if       (outform[o] == "x11"){
              X11(width=size$width,height=size$height,pointsize=ptsz)
            }else if (outform[o] == "png"){
              png(filename=fichier,width=size$width*depth
                 ,height=size$height*depth,pointsize=ptsz,res=depth)
            }else if (outform[o] == "eps"){
              postscript(file=fichier,width=size$width,height=size$height
                        ,pointsize=ptsz,paper=size$paper)
            }else if (outform[o] == "pdf"){
              pdf(file=fichier,onefile=FALSE,width=size$width
                 ,height=size$height,pointsize=ptsz,paper=size$paper)
            }#end if
            #------------------------------------------------------------------------------#



            #----- Set up the plotting window. --------------------------------------------#
            par(par.user)
            par(oma=c(0,1,4,0))
            layout( mat     = rbind(lo.simul$mat.off,rep(1,times=lo.simul$ncol))
                  , heights = c(rep(6/lo.simul$nrow,lo.simul$nrow),1)
                  )#end layout
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #      First plot: the legend.  Skip TIR.                                      #
            #------------------------------------------------------------------------------#
            par(mar=c(0.1,4.1,0.1,2.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1))
            legend( x      = "bottom"
                  , inset  = 0.0
                  , legend = c(band$desc[bwhich])
                  , col    = c(band$col [bwhich])
                  , pch    = 16
                  , pt.cex = 1.2
                  , cex    = 0.9 * cex.ptsz
                  , ncol   = 3
                  , xpd    = TRUE
                  )#end legend
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #      Loop over the simulations, and the bands.                               #
            #------------------------------------------------------------------------------#
            for (s in sequence(nsimul)){
               #------ Open plotting device. ----------------------------------------------#
               par(mar=lo.simul$mar0)
               plot.new()
               plot.window(xlim = lai.limit,ylim = this.limit,log="y")
               axis(side=1)
               axis(side=2,las=1)
               grid(col=grid.colour,lty="dotted")
               #---------------------------------------------------------------------------#



               #----- Loop over bands to be plotted. --------------------------------------#
               for (b in bwhich){
                  this.bar = tapply( X     = c(this[b,s,,])
                                   , INDEX = c(round(tot.lai[s,,],1))
                                   , FUN   = mean
                                   , na.rm = TRUE
                                   )#end if
                  lai.bar  = as.numeric(names(this.bar))
                  points(x=lai.bar,y=this.bar,col=band$col[b],type="p",pch=16,cex=0.8)
                  
               }#end for (b in bwhich)
               #---------------------------------------------------------------------------#


               #----- Plot the final information. -----------------------------------------#
               box()
               title(main=simul$desc[s],font.main=1,line=1)
               #---------------------------------------------------------------------------#
            }#end for (s in sequence(nsimul))
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Plot the general annotation.                                             #
            #------------------------------------------------------------------------------#
            lex     = desc.unit(desc="Patch LAI",unit=untab$m2lom2)
            ley     = desc.unit(desc="Relative Irradiance",unit=untab$empty)
            letitre = paste(desc," - ","Mean relative irradiance beneath the canopy"
                           ,"\n",flux.desc,sep="")
            gtitle( main      = letitre
                  , xlab      = lex
                  , ylab      = ley
                  , line.xlab = 2.6
                  , line.ylab = 2.6
                  , cex.main  = 0.8*cex.ptsz
                  , cex.xlab  = 0.8*cex.ptsz
                  , cex.ylab  = 0.8*cex.ptsz
                  )#end gtitle
            #------------------------------------------------------------------------------#



            #----- Close the device. ------------------------------------------------------#
            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
            clean.tmp()
            #------------------------------------------------------------------------------#
         }#end for (o in sequence(nout))
         #---------------------------------------------------------------------------------#
      }#end for (f in sequence(nflux))
      #------------------------------------------------------------------------------------#
   }#end for (p in sequence(nsites))
#------------------------------------------------------------------------------------------#
}#end if (plot.site)
#==========================================================================================#
#==========================================================================================#
