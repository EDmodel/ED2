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
here           = "thispath"    # Current directory.
there          = "thatpath"    # Directory where analyses/history are 
srcdir         = "thisrscpath" # Source  directory.
outroot        = "thisoutroot" # Directory for figures
pov.incs       = "thispovincs" # Path with POV-Ray include files
#------------------------------------------------------------------------------------------#


#----- Time options. ----------------------------------------------------------------------#
monthbeg    = thismontha   # First month to use
yearbeg     = thisyeara    # First year to consider
yearend     = thisyearz    # Maximum year to consider
reload.data = TRUE         # Should I reload partially loaded data?
pov.month   = 5            # Months for POV-Ray plots
pop.scale   = 1.0          # Scaling factor to REDUCE displayed population.
sasmonth    = sequence(12)
#------------------------------------------------------------------------------------------#



#----- Name of the simulations. -----------------------------------------------------------#
myplaces       = c("thispoly")
#------------------------------------------------------------------------------------------#



#----- Plot options. ----------------------------------------------------------------------#
depth          = 1200                   # PNG resolution, in pixels per inch
paper          = "letter"               # Paper size, to define the plot shape
ptsz           = 14                     # Font size.
ibackground    = mybackground           # Background settings (check load_everything.r)
#------------------------------------------------------------------------------------------#


#------ Miscellaneous settings. -----------------------------------------------------------#
slz.min        = -5.0         # The deepest depth that trees access water.
idbh.type      = myidbhtype   # Type of DBH class
                              # 1 -- Every 10 cm until 100cm; > 100cm
                              # 2 -- 0-10; 10-20; 20-35; 35-50; 50-70; > 70 (cm)
                              # 3 -- 0-10; 10-35; 35-70; > 70 (cm)
klight         = myklight     # Weighting factor for maximum carbon balance
iallom         = myallom      # Allometry
isoil.hydro    = myslhydro    # Soil hydrology method
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



#----- Loading some packages and scripts. -------------------------------------------------#
source(file.path(srcdir,"load.everything.r"))
#------------------------------------------------------------------------------------------#


#----- Avoid unecessary and extremely annoying beeps. -------------------------------------#
options(locatorBell=FALSE)
#------------------------------------------------------------------------------------------#



#----- Load observations. -----------------------------------------------------------------#
obsrfile = paste(srcdir,"LBA_MIP.v9.RData",sep="/")
load(file=obsrfile)
#------------------------------------------------------------------------------------------#



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
   bnpref  = basename(inpref)
   outmain = paste(outroot,place,sep="/")
   outpref = paste(outmain,"povray",sep="/")
   lieu    = thispoi$lieu
   iata    = thispoi$iata
   suffix  = thispoi$iata
   yeara   = thispoi$yeara
   yearz   = thispoi$yearz
   meszz   = thispoi$monz
   #---------------------------------------------------------------------------------------#



   #----- Create the directories in case they don't exist. --------------------------------#
   if (! file.exists(outmain)) dir.create(outmain)
   if (! file.exists(outpref)) dir.create(outpref)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the total number of months that can be loaded this time.                     #
   #---------------------------------------------------------------------------------------#
   ntimes     = (yearz-yeara-1)*12+meszz+(12-monthbeg+1)
   #---------------------------------------------------------------------------------------#



   #----- Print a banner to entretain the user. -------------------------------------------#
   cat(" + Post-processing output from ",lieu,"...","\n")
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Make the RData file name, then we check whether we must read the files again     #
   # or use the stored RData.                                                              #
   #---------------------------------------------------------------------------------------#
   path.data  = paste(here,place,"rdata_month",sep="/")
   if (! file.exists(path.data)) dir.create(path.data)
   ed22.rdata = file.path(path.data,paste0(place,".RData"))
   if (reload.data && file.exists(ed22.rdata)){
      #----- Load the modelled dataset. ---------------------------------------------------#
      cat("   - Loading previous session...","\n")
      load(ed22.rdata)
      tresume = datum$ntimes + 1
      datum   = update.monthly( new.ntimes = ntimes 
                              , old.datum  = datum
                              , montha     = monthbeg
                              , yeara      = yeara
                              , inpref     = inpref
                              , slz.min    = slz.min
                              )#end update.monthly
   }else{
      cat("   - Starting new session...","\n")
      tresume    = 1
      datum      = create.monthly( ntimes  = ntimes
                                 , montha  = monthbeg
                                 , yeara   = yeara
                                 , inpref  = inpref
                                 , slz.min = slz.min
                                 )#end create.monthly
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Check whether we have anything to update.                                         #
   #---------------------------------------------------------------------------------------#
   complete = tresume > ntimes
   #---------------------------------------------------------------------------------------#



   #----- Copy some dimensions to scalars. ------------------------------------------------#
   nzg        = datum$nzg
   nzs        = datum$nzs
   ndcycle    = datum$ndcycle
   isoilflg   = datum$isoilflg
   slz        = datum$slz
   slxsand    = datum$slxsand
   slxclay    = datum$slxclay
   ntext      = datum$ntext
   soil.prop  = datum$soil.prop
   dslz       = datum$dslz
   soil.depth = datum$soil.depth
   soil.dry   = datum$soil.dry
   soil.poro  = datum$soil.poro
   ka         = datum$ka
   kz         = datum$kz
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Loop over all times in case there is anything new to be read.                     #
   #---------------------------------------------------------------------------------------#
   if (! complete){

      #------------------------------------------------------------------------------------#
      #     This function will read the files.                                             #
      #------------------------------------------------------------------------------------#
      datum = read.q.files(datum=datum,ntimes=ntimes,tresume=tresume,sasmonth=sasmonth)
      #------------------------------------------------------------------------------------#

      #------ Save the data to the R object. ----------------------------------------------#
      cat(" + Saving data to ",basename(ed22.rdata),"...","\n")
      save(datum,file=ed22.rdata)
      #------------------------------------------------------------------------------------#
   }#end if (! complete)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Copy some data to local variables.                                                #
   #---------------------------------------------------------------------------------------#
   patch  = datum$patch
   cohort = datum$cohort
   #---------------------------------------------------------------------------------------#


   years   = sort(unique(datum$year))
   pclabs  = paste("y",sprintf("%4.4i",years),"m",sprintf("%2.2i",pov.month),sep="")
   pcwhens = paste(mon2mmm(pov.month,cap1=TRUE),years,sep="-")
   pcouts  = paste(sprintf("%4.4i",years),sprintf("%2.2i",pov.month),sep="-")
   pcloop  = sequence(length(pclabs))

   for (w in pcloop){
      #----- Grab this time. --------------------------------------------------------------#
      pclab  = pclabs [w]
      pcwhen = pcwhens[w]
      pcout  = pcouts [w]
      #------------------------------------------------------------------------------------#
      cat (" + Creating POV-Ray image for ",pcwhen,"...","\n")



      #----- Copy some patch variables to local variables. --------------------------------#
      ipa      = patch$ipa [[pclab]]
      areapa   = patch$area[[pclab]]
      #------------------------------------------------------------------------------------#



      #----- Copy some cohort variables to local variables. -------------------------------#
      ipaco    = cohort$ipa[[pclab]]
      icoco    = cohort$ico[[pclab]]
      nplantco = round(cohort$nplant[[pclab]] * pov.total.area * pop.scale)
      dbhco    = cohort$dbh   [[pclab]]
      pftco    = cohort$pft   [[pclab]]
      #------------------------------------------------------------------------------------#


      #----- Remove small plants to reduce the clutter. -----------------------------------#
      keep     = is.finite(dbhco) & dbhco >= pov.dbh.min
      ipaco    = ipaco   [keep]
      icoco    = icoco   [keep]
      nplantco = nplantco[keep]
      dbhco    = dbhco   [keep]
      pftco    = pftco   [keep]
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Find the tree coordinates.                                                     #
      #------------------------------------------------------------------------------------#
      if (sum(keep) > 0){
         #----- Determine the quadricules where trees of each patch can go. ---------------#
         npatches    = length(areapa)
         nquadpa.1st = round( pov.nxy.patch * areapa )
         totquad     = sum(nquadpa.1st)
         off         = pov.nxy.patch - totquad
         nfix        = sample( c( rep( x = 0        , times = npatches-abs(off))
                                , rep( x = sign(off), times = abs(off)         )
                                )#end c
                             )#end sample
         nquadpa     = nquadpa.1st + nfix
         if (sum(nquadpa) != pov.nxy.patch){
            cat (" - NQUADPA.1ST: ",sum(nquadpa.1st),"\n")
            cat (" - NQUADPA    : ",sum(nquadpa)    ,"\n")
            cat (" - NXY.POV    : ",pov.nxy.patch   ,"\n")
            stop(" Not the correct number of patches")
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Create a list with the quadricules for each patch.  We double each count   #
         # to make sure each patch gets at least 2 numbers (so we can safely use 'sample'. #
         #---------------------------------------------------------------------------------#
         ipa.quad    = unlist(mapply(FUN=rep,x=ipa,each=2*nquadpa))
         quad        = split(x=rep(sample(pov.nxy.patch),each=2),f=ipa.quad)
         names(quad) = paste("patch",sprintf("%3.3i",sort(unique(ipa))),sep="_")
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Expand population variables.                                                #
         #---------------------------------------------------------------------------------#
         ipaco  = unlist(mapply(FUN=rep,x=ipaco,each=nplantco))
         dbhco  = unlist(mapply(FUN=rep,x=dbhco,each=nplantco))
         pftco  = unlist(mapply(FUN=rep,x=pftco,each=nplantco))
         nco    = sum(nplantco)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Randomly assign the cohorts to the patches where they belong.               #
         #---------------------------------------------------------------------------------#
         nco.pa.pop  = table(ipaco)
         lab         = paste("patch",sprintf("%3.3i",as.integer(names(nco.pa.pop))),sep="_")
         idx         = match(lab,names(quad))
         nco.pa      = rep(0,times=length(quad))
         nco.pa[idx] = nco.pa.pop
         quadco      = unlist( mapply( FUN      = sample
                                     , x        = quad
                                     , size     = nco.pa
                                     , MoreArgs = list(replace=TRUE)
                                     )#end mapply
                             )#end unlist
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Assign coordinates for all cohorts.                                         #
         #---------------------------------------------------------------------------------#
         xco    = ( sample(x=0.1*(sequence(10*pov.patch.xmax)-0.5),size=nco,replace=TRUE)
                  + pov.x0[quadco] )
         yco    = ( sample(x=0.1*(sequence(10*pov.patch.ymax)-0.5),size=nco,replace=TRUE)
                  + pov.y0[quadco] )
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Make the labels.                                                            #
         #---------------------------------------------------------------------------------#
         povplant = rbind(       "//----- The plants. ---------------------------------//"
                         , rbind( paste("plant(", unlist( mapply( FUN = paste
                                                                , sprintf("%2i"  ,iallom)
                                                                , sprintf("%7.2f",dbhco )
                                                                , sprintf("%2i"  ,pftco )
                                                                , sprintf("%7.2f",xco   )
                                                                , sprintf("%7.2f",yco   )
                                                                , MoreArgs = list(sep=",")
                                                                )#end mapply
                                                        )#end unlist
                                       ,")"
                                       ,sep=""
                                       )#end paste
                                )#rbind
                         ,       "//---------------------------------------------------//"
                         )#end rbind
         #---------------------------------------------------------------------------------#
      }else{
         povplant = rbind(       "//----- No plants. ----------------------------------//"
                        ,        "//---------------------------------------------------//"
                         )#end rbind
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Copy the POV-ray template file to the working directory and append the set-   #
      # tings for this time.                                                               #
      #------------------------------------------------------------------------------------#
      povscript = file.path(here,place,"polygon.pov")
      dummy     = file.copy(from=file.path(srcdir,"povray_template.pov"),to=povscript)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Decide the colour of the text depending on the background.                    #
      #------------------------------------------------------------------------------------#
      if (ibackground == 0){
         pigment = "       pigment  { color rgb <0. ,0. ,0. >}   "
      }else if (ibackground == 1){
         pigment = "       pigment  { color rgb <1. ,1. ,1. >}   "
      }else if (ibackground == 2){
         pigment = "       pigment  { color rgb <1. ,1. ,1. >}   "
      }#end if
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #      Append the title and time stamp.                                              #
      #------------------------------------------------------------------------------------#
      lesim    = unlist(strsplit(lieu,split="\n"))
      povtitle = NULL
      for (n in sequence(length(lesim))){
         yt       = sprintf("%7.1f",200 - 20 * (n-1))
         povtitle = rbind( povtitle
                         ,       "//----- The header. --------------------------------//"
                         , paste("text { ttf \"cyrvetic.ttf\" \"",lesim[n],"\" 5,0",sep="")
                         ,               pigment
                         ,       "       finish{ ambient 1.0 diffuse 0.0}"
                         ,       "       scale     <   12.0,   12.0,    0.1>"
                         ,       "       rotate    <   28.8,    0.0,    0.0>"
                         ,       "       rotate    <    0.0,   45.0,    0.0>"
                         , paste("       translate < -150.0,", yt,",  200.0>",sep="")
                         ,       "     }// end text"
                         ,       "//--------------------------------------------------//"
                         ,       " "
                         ,       " "
                         )#end rbind
      }#end for (n in 1:lesim)
      povstamp = rbind(       "//----- The time stamp. ----------------------------//"
                      , paste("text { ttf \"cyrvetic.ttf\" \"",pcwhen,"\" 5,0",sep="")
                      ,        pigment
                      ,       "       finish{ ambient 1.0 diffuse 0.0}"
                      ,       "       scale     <   12.0,   12.0,    0.1>"
                      ,       "       rotate    <   28.8,    0.0,    0.0>"
                      ,       "       rotate    <    0.0,   45.0,    0.0>"
                      ,       "       translate <  100.0,  180.0, -150.0>"
                      ,       "     }// end text"
                      ,       "//--------------------------------------------------//"
                      ,       " "
                      ,       " "
                      )#end rbind
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Append the title, stamp, and plants...                                        #
      #------------------------------------------------------------------------------------#
      write(x = povtitle, file = povscript, ncolumns = 1, append = TRUE)
      write(x = povstamp, file = povscript, ncolumns = 1, append = TRUE)
      write(x = povplant, file = povscript, ncolumns = 1, append = TRUE)
      #------------------------------------------------------------------------------------#

      #------------------------------------------------------------------------------------#
      #     Call POV-Ray (must use system, though...)                                      #
      #------------------------------------------------------------------------------------#
      povray  = Sys.which("povray")
      outtemp = file.path(tempdir(),paste(place,".png",sep=""))
      outfile = file.path(outpref,paste(bnpref,"-",pcout,".png",sep=""))

      povopts = paste("-D"
                     ,"-V"
                     ,"+UA"
                     ,paste0("+L",pov.incs)
                     ,paste0("+W",round(size$width*depth ))
                     ,paste0("+H",round(size$height*depth))
                     ,paste0("+O",outtemp)
                     ,sep = " "
                     )#end paste
      dummy   = system( command       = paste(povray,povopts,povscript,sep=" ")
                      , intern        = TRUE
                      , ignore.stdout = TRUE
                      , ignore.stderr = TRUE
                      )#end system
      dummy   = file.copy(from=outtemp,to=outfile,overwrite=TRUE)
      dummy   = file.remove(outtemp,povscript)
      #------------------------------------------------------------------------------------#
   }#end for (w in 1:pclabs)
   #---------------------------------------------------------------------------------------#
}#end for (place in myplaces)
#==========================================================================================#
#==========================================================================================#
