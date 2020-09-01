#==========================================================================================#
#==========================================================================================#
#     This function retrieves the information of a given site using its name to determine  #
# the properties.  To add more sites, go to the end of this file and adjust the lists      #
# accordingly.                                                                             #
#------------------------------------------------------------------------------------------#
locations <<- function(where,here=getwd(),yearbeg=1500,yearend=2008,monthbeg=1,daybeg=1
                      ,filetype="Q",fullonly=FALSE){
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     This is the name of the place, made case-insensitive.  Output directories will be #
   # always in lower case.                                                                 #
   #---------------------------------------------------------------------------------------#
   ici = tolower(where)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Now we check which place or simulation we are running.                           #
   #---------------------------------------------------------------------------------------#
   if (ici %in% poilist$short){
      #----- Regular POI. -----------------------------------------------------------------#
      pp       = match(ici,poilist$short)
      lieu     = poilist$longname[pp]
      iata     = poilist$iata[pp]
      wmo      = poilist$wmo [pp]
      pathroot = file.path(here,ici)
      pathin   = file.path(here,ici,"analy",ici)
      pathrst  = file.path(here,ici,"histo",ici)
      pathout  = file.path(here,ici,"epost")
      lon      = poilist$lon[pp]
      lat      = poilist$lat[pp]

   }else if( ( substring(ici,1,1) %in% c("t","s","e","z") )
           & ( ( substring(ici,5,5) == "_" ) | ( nchar(ici) == 4 ) ) ){
      #---- IATA-based name with model configuration. -------------------------------------#
      metflag = substring(ici,1,1)
      if (tolower(metflag) == "s"){
         thismet = " (Sheffield)"
      }else if (tolower(metflag) == "t"){
         thismet = " (Tower)"
      }else if (tolower(metflag) %in% c("e","z")){
         thismet = ""
      }#end if

      #----- Grab the IATA code. ----------------------------------------------------------#
      iata    = substring(ici,2,4)
      pp      = which(poilist$iata == iata)

      #----- Put the name of the place and the meteorological forcing. --------------------#
      testpoi = paste0(poilist$longname[pp],thismet)
      lon     = poilist$lon[pp]
      lat     = poilist$lat[pp]
      wmo     = poilist$wmo[pp]

      lieu     = simul.description(ici,testpoi,iata=TRUE)
      pathroot = file.path(here,ici)
      pathin   = file.path(pathroot,"analy",ici)
      pathrst  = file.path(here,ici,"histo",ici)
      pathout  = file.path(pathroot,"epost")

   }else if( substring(ici,4,4) == "_" & substring(ici,9,9) == "_"){
      #---- Convert back to upper case. ---------------------------------------------------#
      ici      = paste0(toupper(substring(ici,1,3)),substring(ici,4))

      #---- Regional polygon. -------------------------------------------------------------#
      pnumber  = as.numeric(substring(ici,5,8))
      testpoi  = paste("Polygon",pnumber)

      #----- Grab the "IATA" code. --------------------------------------------------------#
      iata    = substring(ici,1,8)

      #----- Put the name of the place and the meteorological forcing. --------------------#
      lon     = as.numeric(substring(ici,13,18))
      lat     = as.numeric(substring(ici,23,28))
      wmo     = NA

      lieu     = simul.description(ici,testpoi,iata=FALSE)
      pathroot = file.path(here,ici)
      pathin   = file.path(pathroot,"analy",ici)
      pathrst  = file.path(here,ici,"histo",ici)
      pathout  = file.path(pathroot,"epost")

   }else if(ici == "rondonia")         {
       lieu     = "Rondonia"
       iata     = "RO"
       wmo      = NA
       pathroot = file.path(here,ici)
       pathin   = file.path(here,"analy",ici)
       pathrst  = file.path(here,"histo",ici)
       pathout  = file.path(here,"epost")
       lon      = NA
       lat      = NA

   }else if(ici == "potveg")         {
       lieu     = "South America"
       iata     = "potveg"
       wmo      = NA
       pathroot = file.path(here,ici)
       pathin   = file.path(here,"analy",ici)
       pathrst  = file.path(here,"histo",ici)
       pathout  = file.path(here,"epost")
       lon      = NA
       lat      = NA

   }else if(ici == "coupled")         {
       lieu     = "Coupled model"
       iata     = "coupled"
       wmo      = NA
       pathroot = here
       pathin   = file.path(here,"ecoss","coupled")
       pathrst  = file.path(here,"ecort","coupled")
       pathout  = file.path(here,"epost")
       lon      = NA
       lat      = NA

   }else{
       lieu     = NA
       iata     = NA
       wmo      = NA
       pathroot = NA
       pathin   = NA
       pathrst  = NA
       pathout  = NA
       lon      = NA
       lat      = NA
   } #end if
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #     Obtain information on how far this simulation has run.                            #
   #---------------------------------------------------------------------------------------#
   if (filetype == "Q"){
      #------- Now we loop over years and months to determine what is the last file... ----#
      yeara = yearbeg
      yearz = yearbeg - 1
      yr    = yearend + 1
      while (yr > yeara){
         yr = yr - 1
         if (yr == yeara){
            mona=monthbeg
         }else{
            mona=1
         }#end if
         cmon=sprintf("%2.2i",mona)
         filename     = paste0(pathin,"-Q-",yr,"-",cmon,"-00-000000-g01.h5"    )
         filename.bz2 = paste0(pathin,"-Q-",yr,"-",cmon,"-00-000000-g01.h5.bz2")
         filename.gz  = paste0(pathin,"-Q-",yr,"-",cmon,"-00-000000-g01.h5.gz" )
         if (any(file.exists(c(filename,filename.bz2,filename.gz)))){
            yearz = yr
            yr    = yeara #----- This will make it leave the loop. ------------------------#
         } #end if
      } #end while
      monz   = mona
      mo     = 13
      while (mo > mona){
         mo           = mo -1
         cmon         = sprintf("%2.2i",mo)
         filename     = paste0(pathin,"-Q-",yearz,"-",cmon,"-00-000000-g01.h5"    )
         filename.bz2 = paste0(pathin,"-Q-",yearz,"-",cmon,"-00-000000-g01.h5.bz2")
         filename.gz  = paste0(pathin,"-Q-",yearz,"-",cmon,"-00-000000-g01.h5.gz" )
         if (any(file.exists(c(filename,filename.bz2,filename.gz)))){
            monz = mo
            mo   = 1
         }#end if
      }#end while
      dayz = daymax(monz,yearz)
   }else if (filetype == "E"){
      #------- Now we loop over years and months to determine what is the last file... ----#
      yeara = yearbeg
      yearz = yearbeg - 1
      yr    = yearend + 1
      while (yr > yeara){
         yr = yr - 1
         if (yr == yeara){
            mona=monthbeg
         }else{
            mona=1
         }#end if
         cmon         = sprintf("%2.2i",mona)
         filename     = paste0(pathin,"-E-",yr,"-",cmon,"-00-000000-g01.h5"    )
         filename.bz2 = paste0(pathin,"-E-",yr,"-",cmon,"-00-000000-g01.h5.bz2")
         filename.gz  = paste0(pathin,"-E-",yr,"-",cmon,"-00-000000-g01.h5.gz" )
         if (any(file.exists(c(filename,filename.bz2,filename.gz)))){
            yearz = yr
            yr    = yeara #----- This will make it leave the loop. ------------------------#
         } #end if
      } #end while
      monz   = mona
      mo     = 13
      while (mo > mona){
         mo           = mo -1
         cmon         = sprintf("%2.2i",mo)
         filename     = paste0(pathin,"-E-",yearz,"-",cmon,"-00-000000-g01.h5"    )
         filename.bz2 = paste0(pathin,"-E-",yearz,"-",cmon,"-00-000000-g01.h5.bz2")
         filename.gz  = paste0(pathin,"-E-",yearz,"-",cmon,"-00-000000-g01.h5.gz" )
         if (any(file.exists(c(filename,filename.bz2,filename.gz)))){
            monz = mo
            mo   = 1
         }#end if
      }#end while
      dayz = daymax(monz,yearz)
   }else if (filetype == "S"){
      #------- Now we loop over years and months to determine what is the last file... ----#
      yeara = yearbeg
      yearz = yearbeg - 1
      yr    = yearend + 1
      while (yr > yeara){
         yr = yr - 1
         if (yr == yeara){
            mona=monthbeg
         }else{
            mona=1
         }#end if
         cmon=sprintf("%2.2i",mona)
         filename     = paste0(pathrst,"-S-",yr,"-",cmon,"-01-000000-g01.h5"    )
         filename.bz2 = paste0(pathrst,"-S-",yr,"-",cmon,"-01-000000-g01.h5.bz2")
         filename.gz  = paste0(pathrst,"-S-",yr,"-",cmon,"-01-000000-g01.h5.gz" )
         if (any(file.exists(c(filename,filename.bz2,filename.gz)))){
            yearz = yr
            yr    = yeara #----- This will make it leave the loop. ------------------------#
         } #end if
      } #end while
      monz   = mona
      mo     = 13
      while (mo > mona){
         mo           = mo -1
         cmon         = sprintf("%2.2i",mo)
         filename     = paste0(pathrst,"-S-",yearz,"-",cmon,"-01-000000-g01.h5"    )
         filename.bz2 = paste0(pathrst,"-S-",yearz,"-",cmon,"-01-000000-g01.h5.bz2")
         filename.gz  = paste0(pathrst,"-S-",yearz,"-",cmon,"-01-000000-g01.h5.gz" )
         if (any(file.exists(c(filename,filename.bz2,filename.gz)))){
            monz = mo
            mo   = 1
         }#end if
      }#end while
      dayz = daymax(monz,yearz)
   }else if (filetype == "D"){
      #------- Now we loop over years and months to determine what is the last file... ----#
      yeara = yearbeg
      yearz = yearbeg - 1
      yr    = yearend + 1
      while (yr > yeara){
         yr = yr - 1
         if (yr == yeara){
            mona=monthbeg
            daya=daybeg
         }else{
            mona=1
            daya=1
         }#end if
         cmon         = sprintf("%2.2i",mona)
         cday         = sprintf("%2.2i",daya)
         filename     = paste0(pathin,"-D-",yr,"-",cmon,"-",cday,"-000000-g01.h5"    )
         filename.bz2 = paste0(pathin,"-D-",yr,"-",cmon,"-",cday,"-000000-g01.h5.bz2")
         filename.gz  = paste0(pathin,"-D-",yr,"-",cmon,"-",cday,"-000000-g01.h5.gz" )
         if (any(file.exists(c(filename,filename.bz2,filename.gz)))){
            yearz = yr
            yr    = yeara #----- This will make it leave the loop. ------------------------#
         } #end if
      } #end while

      monz   = mona
      mo     = 13
      while (mo > mona){
         mo       = mo -1
         if (yearz == yeara && mo == mona){
            daya = daybeg
         }else{
            daya = 1
         }#end if
         cmon         = sprintf("%2.2i",mo  )
         cday         = sprintf("%2.2i",daya)
         filename     = paste0(pathin,"-D-",yearz,"-",cmon,"-",cday,"-000000-g01.h5"    )
         filename.bz2 = paste0(pathin,"-D-",yearz,"-",cmon,"-",cday,"-000000-g01.h5.bz2")
         filename.gz  = paste0(pathin,"-D-",yearz,"-",cmon,"-",cday,"-000000-g01.h5.gz" )
         if (any(file.exists(c(filename,filename.bz2,filename.gz)))){
            monz = mo
            mo   = 1
         }#end if
      }#end while

      dayz = daya
      dy   = daymax(monz,yearz)
      cmon = sprintf("%2.2i",monz)
      while (dy > daya){
         dy           = dy - 1
         cday         = sprintf("%2.2i",dy)
         filename     = paste0(pathin,"-D-",yearz,"-",cmon,"-",cday,"-000000-g01.h5"    )
         filename.bz2 = paste0(pathin,"-D-",yearz,"-",cmon,"-",cday,"-000000-g01.h5.bz2")
         filename.gz  = paste0(pathin,"-D-",yearz,"-",cmon,"-",cday,"-000000-g01.h5.gz" )
         if (any(file.exists(c(filename,filename.bz2,filename.gz)))){
            dayz = dy
            dy   = daya #----- This will make it leave the loop. -------------------------#
         } #end if
      }#end while

      if (fullonly && dayz < daymax(monz,yearz)){
         monz = monz -1
         if (monz == 0){
            monz  = 12
            yearz = yearz - 1
         }#end if
         dayz = daymax(monz,yearz)
      }#end if
   }else if (filetype == "I"){
      #------- Now we loop over years and months to determine what is the last file... ----#
      yeara = yearbeg
      yearz = yearbeg - 1
      yr    = yearend + 1
      while (yr > yeara){
         yr = yr - 1
         if (yr == yeara){
            mona=monthbeg
            daya=daybeg
         }else{
            mona=1
            daya=1
         }#end if
         cmon         = sprintf("%2.2i",mona)
         cday         = sprintf("%2.2i",daya)
         filename     = paste0(pathin,"-I-",yr,"-",cmon,"-",cday,"-000000-g01.h5"    )
         filename.bz2 = paste0(pathin,"-I-",yr,"-",cmon,"-",cday,"-000000-g01.h5.bz2")
         filename.gz  = paste0(pathin,"-I-",yr,"-",cmon,"-",cday,"-000000-g01.h5.gz" )
         if (any(file.exists(c(filename,filename.bz2,filename.gz)))){
            yearz = yr
            yr    = yeara #----- This will make it leave the loop. ------------------------#
         } #end if
      } #end while

      monz   = mona
      mo     = 13
      while (mo > mona){
         mo       = mo -1
         if (yearz == yeara && mo == mona){
            daya = daybeg
         }else{
            daya = 1
         }#end if
         cmon         = sprintf("%2.2i",mo  )
         cday         = sprintf("%2.2i",daya)
         filename     = paste0(pathin,"-I-",yearz,"-",cmon,"-",cday,"-000000-g01.h5"    )
         filename.bz2 = paste0(pathin,"-I-",yearz,"-",cmon,"-",cday,"-000000-g01.h5.bz2")
         filename.gz  = paste0(pathin,"-I-",yearz,"-",cmon,"-",cday,"-000000-g01.h5.gz" )
         if (any(file.exists(c(filename,filename.bz2,filename.gz)))){
            monz = mo
            mo   = 1
         }#end if
      }#end while

      dayz = daya
      dy   = daymax(monz,yearz)
      cmon = sprintf("%2.2i",monz)
      while (dy > daya){
         dy           = dy - 1
         cday         = sprintf("%2.2i",dy)
         filename     = paste0(pathin,"-I-",yearz,"-",cmon,"-",cday,"-000000-g01.h5"    )
         filename.bz2 = paste0(pathin,"-I-",yearz,"-",cmon,"-",cday,"-000000-g01.h5.bz2")
         filename.gz  = paste0(pathin,"-I-",yearz,"-",cmon,"-",cday,"-000000-g01.h5.gz" )
         if (any(file.exists(c(filename,filename.bz2,filename.gz)))){
            dayz = dy
            dy   = daya #----- This will make it leave the loop. -------------------------#
         } #end if
      }#end while

      if (fullonly && dayz < daymax(monz,yearz)){
         monz = monz -1
         if (monz == 0){
            monz  = 12
            yearz = yearz - 1
         }#end if
         dayz = daymax(monz,yearz)
      }#end if

   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Return the information about this place in a neatly organised list.                #
   #---------------------------------------------------------------------------------------#
   thisplace = list(lieu=lieu,iata=iata,lon=lon,lat=lat,
                    pathroot=pathroot,pathin=pathin,pathrst=pathrst,pathout=pathout,
                    yeara=yeara,yearz=yearz,mona=monthbeg,monz=monz,daya=daybeg,dayz=dayz)
   return(thisplace)
   #---------------------------------------------------------------------------------------#
} #end function
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     This function finds a good name to describe the simulation.                          #
#------------------------------------------------------------------------------------------#
simul.description <<- function(ici,testpoi,iata=TRUE,max.char=66){

   #---------------------------------------------------------------------------------------#
   #     This block contains names for the commonest settings for variables that are       #
   # controlled by flags.                                                                  #
   #---------------------------------------------------------------------------------------#
   flagvar = list()
   #----- imetrad is the radiation splitting method. --------------------------------------# 
   flagvar[["imetrad"]]         = list( descr   = "SW Radiation splitting"
                                      , numeric = TRUE
                                      , values  = seq(from=0,to=5,by=1)
                                      , names   = c( "Input files"
                                                   , "SiB"
                                                   , "Weissman-Norman"
                                                   , "100% diffuse"
                                                   , "100% direct"
                                                   , "Clearness Index"
                                                   )#end c
                                      )#end list
   #----- icanrad is the canopy radiation model. ------------------------------------------#
   flagvar[["icanrad"]]         = list( descr   = "Canopy radiation method"
                                      , numeric = TRUE
                                      , values  = seq(from=0,to=2,by=1)
                                      , names   = c("Medvigy","Zhao and Qualls","Liou")
                                      )#end list
   #----- ihrzrad is the horizontal shading model. ----------------------------------------#
   flagvar[["ihrzrad"]]         = list( descr   = "Horizontal shading"
                                      , numeric = TRUE
                                      , values  = seq(from=0,to=4,by=1)
                                      , names   = c("OFF","GAP","PIXEL","DUMMY","OFF")
                                      )#end list
   #----- imetavg is the averaging method for met driver. ---------------------------------#
   flagvar[["imetavg"]]         = list( descr   = "Met driver average"
                                      , numeric = TRUE
                                      , values  = seq(from=-1,to=3,by=1)
                                      , names   = c("Linear"
                                                  ,"Cos(Zen) - Instantaneous"
                                                  ,"Cos(Zen) - Forward"
                                                  ,"Cos(Zen) - Backward"
                                                  ,"Cos(Zen) - Centred")
                                      )#end list
   #----- iallom is the allometry method. -------------------------------------------------#
   flagvar[["iallom"]]          = list( descr   = "Allometry"
                                      , numeric = TRUE
                                      , values  = seq(from=0,to=3,by=1)
                                      , names   = c("ED-1.0","ED-2.0","ED-2.1","ED-2.2")
                                      )#end list
   #----- ieconomics is the trait trade-off method. ---------------------------------------#
   flagvar[["ieconomics"]]      = list( descr   = "Economics"
                                      , numeric = TRUE
                                      , values  = seq(from=0,to=1,by=1)
                                      , names   = c("ED-2.0","TNRG data sets")
                                      )#end list
   #----- h2o.plant.limit is the water limitation method. ---------------------------------#
   flagvar[["h2o.plant.limit"]] = list( descr   = "Water limitation"
                                      , numeric = TRUE
                                      , values  = seq(from=0,to=3,by=1)
                                      , names   = c("No limitation"
                                                  ,"Demand/Supply"
                                                  ,"Matric potl. on fsw"
                                                  ,"Matric potl. on gsw")
                                      )#end list
   #----- h2o.plant.limit is the water limitation method. ---------------------------------#
   flagvar[["dd.mort.control"]] = list( descr   = "Dens.-dep. mortality"
                                      , numeric = TRUE
                                      , values  = seq(from=2,to=4,by=1)
                                      , names   = c("Light"
                                                   ,"Light+soil water"
                                                   ,"Light+water(air/soil)")
                                      )#end list
   #----- ibranch.thermo is the branch thermodynamics. ------------------------------------#
   flagvar[["ibranch.thermo"]]  = list( descr   = "Branch thermodynamics"
                                      , numeric = TRUE
                                      , values  = seq(from=0,to=2,by=1)
                                      , names   = c("No branches"
                                                  ,"Wood+Leaf together"
                                                  ,"Wood and leaf (separate)")
                                      )#end list

   #----- isfclyrm is the Surface layer model. --------------------------------------------#
   flagvar[["isfclyrm"]]        = list( descr   = "Sfc. Lyr. Model"
                                      , numeric = TRUE
                                      , values  = seq(from=1,to=4,by=1)
                                      , names   = c("Louis (1979)"
                                                  ,"Oncley and Dudhia (1995)"
                                                  ,"Beljaars and Holtslag (1991)"
                                                  ,"CLM (2004)")
                                      )#end list
   #----- icanturb is the roughness model. ------------------------------------------------#
   flagvar[["icanturb"]]        = list( descr   = "Roughness Model"
                                      , numeric = TRUE
                                      , values  = seq(from=0,to=4,by=1)
                                      , names   = c("Leuning (1995)"
                                                  ,"ED-2.1 default"
                                                  ,"Massman (1997)"
                                                  ,"Massman and Weil (1999)"
                                                  ,"Leuning (1995) + CLM")
                                      )#end list
   #----- igrndvap is the ground vapour model. --------------------------------------------#
   flagvar[["igrndvap"]]        = list( descr   = "Gnd. Vapour Model"
                                      , numeric = TRUE
                                      , values  = seq(from=0,to=4,by=1)
                                      , names   = c("Lee and Pielke (1991)"
                                                  ,"Mahfouf and Noilhan (1991) - T1"
                                                  ,"Mahfouf and Noilhan (1991) - T2"
                                                  ,"Mahfouf and Noilhan (1991) - T3"
                                                  ,"Mahfouf and Noilhan (1991) - T4")
                                      )#end list
   #----- isoil.text is the soil texture class. -------------------------------------------#
   flagvar[["isoil.text"]]      = list( descr   = "Soil texture type"
                                      , numeric = TRUE
                                      , values  = seq(from=0,to=17,by=1)
                                      , names   = c("Site default"
                                                  ,"Sand"
                                                  ,"Loamy sand"
                                                  ,"Sandy loam"
                                                  ,"Silty loam"
                                                  ,"Loam"
                                                  ,"Sandy clay loam"
                                                  ,"Silty clay loam"
                                                  ,"Clayey loam"
                                                  ,"Sandy clay"
                                                  ,"Silty clay"
                                                  ,"Clay"
                                                  ,"Peat"
                                                  ,"Bedrock"
                                                  ,"Silt"
                                                  ,"Heavy clay"
                                                  ,"Clayey sand"
                                                  ,"Clayey silt")
                                      )#end list
   #----- crown.mod is the crown model for radiation. -------------------------------------#
   flagvar[["crown.mod"]]       = list( descr   = "Crown model"
                                      , numeric = TRUE
                                      , values  = seq(from=0,to=1,by=1)
                                      , names   = c("ED-2.1 original","Finite crown")
                                      )#end list
   #----- ivegt.dynamics is the flag for vegetation dynamics. -----------------------------#
   flagvar[["ivegt.dynamics"]]  = list( descr   = "Vegetation dynamics"
                                      , numeric = TRUE
                                      , values  = seq(from=0,to=2,by=1)
                                      , names   = c("OFF","ON","Multi")
                                      )#end list
   #----- iphen.scheme is the phenology scheme for tropical broadleaf trees. --------------#
   flagvar[["iphen.scheme"]]    = list( descr   = "Phenology scheme"
                                      , numeric = TRUE
                                      , values  = seq(from=-1,to=3,by=1)
                                      , names   = c("Evergreen"
                                                  ,"Drought - ED-2.0"
                                                  ,"Prescribed"
                                                  ,"Drought - ED-2.1"
                                                  ,"Light + drought")
                                      )#end list
   #----- idiversity is an outside flag to tell how many PFTs were allowed. ---------------#
   flagvar[["idiversity"]]      = list( descr   = "PFTs used"
                                      , numeric = TRUE
                                      , values  = c(2,5)
                                      , names   = c("1 Grass + 1 Tree"
                                                   ,"2 Grasses + 3 Trees"
                                                   )#end names
                                      )#end list
   #----- ipatch is an outside flag to tell how many patches existed. ---------------------#
   flagvar[["ipatch"]]      = list( descr   = "Patch structure"
                                  , numeric = TRUE
                                  , values  = c(0,1,25,40)
                                  , names   = c("Multiple patches"
                                               ,"Single patch"
                                               ,"Multiple patches"
                                               ,"Multiple patches"
                                               )#end names
                                      )#end list
   #----- idrought is the pseudo-drought flag. --------------------------------------------#
   flagvar[["idrought"]]        = list( descr   = "Met forcing"
                                      , numeric = TRUE
                                      , values  = seq(from=0,to=1,by=1)
                                      , names   = c("Tower (2004)","Pseudo Drought")
                                      )#end list
   #----- itemp.drought is the temperature response. --------------------------------------#
   flagvar[["imet.drought"]]    = list( descr   = "Tower"
                                      , numeric = TRUE
                                      , values  = c(66,67,77,83)
                                      , names   = c("Km 66","Km 67","Km 77","Km 83")
                                      )#end list
   #----- itemp.drought is the temperature response. --------------------------------------#
   flagvar[["itemp.drought"]]   = list( descr   = "Temperature"
                                      , numeric = TRUE
                                      , values  = seq(from=0,to=1,by=1)
                                      , names   = c("Normal","2C warmer")
                                      )#end list
   #----- irad.drought is the radiation response. -----------------------------------------#
   flagvar[["irad.drought"]]    = list( descr   = "Radiation"
                                      , numeric = TRUE
                                      , values  = seq(from=0,to=1,by=1)
                                      , names   = c("Normal","15% Brighter")
                                      )#end list
   #----- ihum.drought is the humidity response. ------------------------------------------#
   flagvar[["ihum.drought"]]    = list( descr   = "Humidity"
                                      , numeric = TRUE
                                      , values  = seq(from=0,to=1,by=1)
                                      , names   = c("Same Specific humidity"
                                                  ,"Same Relative humidity")
                                      )#end list
   #----- irad.drought is the rainfall response. ------------------------------------------#
   flagvar[["irain.drought"]]    = list( descr  = "Rainfall"
                                      , numeric = TRUE
                                       , values = seq(from=0,to=1,by=1)
                                       , names  = c("Normal","2010 drought")
                                       )#end list
   #----- Photosynthesis parameters. ------------------------------------------------------#
   flagvar[["iphysiol"]]         = list( descr  = "Photosynthesis"
                                      , numeric = TRUE
                                       , values = seq(from=0,to=3,by=1)
                                       , names  = c("Arrhenius (no Jmax/TPmax)"
                                                   ,"Arrenhius (with Jmax/TPmax)"
                                                   ,"Q10 function (no Jmax/TPmax)"
                                                   ,"Q10 function (with Jmax/TPmax)"
                                                   )#end c
                                       )#end list
   #----- Grass type. ---------------------------------------------------------------------#
   flagvar[["igrass"]]           = list( descr  = "Grass scheme"
                                      , numeric = TRUE
                                       , values = c(0,1)
                                       , names  = c("ED-1","Swann")
                                       )#end list
   #----- Grass type. ---------------------------------------------------------------------#
   flagvar[["isoilbc"]]          = list( descr  = "Soil Bnd. Cond."
                                      , numeric = TRUE
                                       , values = c(0,1,2,3,4)
                                       , names  = c( "Bedrock"
                                                   , "Free drainage"
                                                   , "Sinkhole"
                                                   , "Field capacity"
                                                   , "Aquifer"
                                                   , "Lateral drainage")
                                       )#end list
   #----- Grass type. ---------------------------------------------------------------------#
   flagvar[["ipercol"]]          = list( descr  = "Percolation"
                                      , numeric = TRUE
                                       , values = c(0,1,2)
                                       , names  = c( "Walko et al. (2000)"
                                                   , "Anderson (1976)"
                                                   , "Anderson + CLM")
                                       )#end list
   #----- Grass type. ---------------------------------------------------------------------#
   flagvar[["cgrass"]]           = list( descr  = "Grass type"
                                      , numeric = FALSE
                                       , values = c("c3","c4","cb")
                                       , names  = c("C3","C4","C3 and C4")
                                       )#end list
   #----- Fire model. ---------------------------------------------------------------------#
   flagvar[["include.fire"]]     = list( descr   = "Fire model"
                                       , numeric = TRUE
                                       , values  = c(0,1,2,3)
                                       , names   = c( "Off"
                                                    , "ED-1.0 type"
                                                    , "Default ED-2.1"
                                                    , "Water-deficit based"
                                                    )#end c
                                       )#end list
   #----- Energy budget type. -------------------------------------------------------------#
   flagvar[["iloss"]]            = list( descr  = "Soil internal energy lost for transp."
                                      , numeric = TRUE
                                       , values = c(0,1)
                                       , names  = c("Disappears","Goes to leaves")
                                       )#end list
   #----- Revision. -----------------------------------------------------------------------#
   flagvar[["revision"]]         = list( descr  = "Revision"
                                      , numeric = TRUE
                                       , values = c(78,280)
                                       , names  = c("Mainline - 78"
                                                   ,"Marcos   - 280")
                                       )#end list
   #----- isoil.text is the soil texture class. -------------------------------------------#
   flagvar[["pft"]]             = list( descr   = "PFT used"
                                      , numeric = TRUE
                                      , values  = seq(from=1,to=17,by=1)
                                      , names   = c("C4 grass"
                                                   ,"Early tropical"
                                                   ,"Mid tropical"
                                                   ,"Late tropical"
                                                   ,"Temp. C3 grass"
                                                   ,"North pine"
                                                   ,"South pine"
                                                   ,"Late conifer"
                                                   ,"Early hardwood"
                                                   ,"Mid hardwood"
                                                   ,"Late hardwood"
                                                   ,"C3 pasture"
                                                   ,"C3 crop"
                                                   ,"C4 pasture"
                                                   ,"C4 crop"
                                                   ,"C3 grass"
                                                   ,"Araucaria")
                                      )#end list
   #----- Vegetation structure. -----------------------------------------------------------#
   flagvar[["isas"]]             = list( descr  = "Structure"
                                       , numeric = FALSE
                                       , values = c("sas","ble")
                                       , names  = c("Size structure"
                                                   ,"Single size")
                                       )#end list
   #----- Density-dependent mortality. ----------------------------------------------------#
   flagvar[["iddmort.scheme"]]  = list( descr  = "Carbon balance method"
                                       , numeric = TRUE
                                       , values = c(0,1)
                                       , names  = c("Rates"
                                                   ,"Rates and storage"
                                                   )#end names
                                       )#end list
   #----- Respiration scheme mortality. ---------------------------------------------------#
   flagvar[["resp.scheme"]]     = list( descr  = "Respiration method"
                                       , numeric = TRUE
                                       , values = c(0,1)
                                       , names  = c("ED-2.1"
                                                   ,"Based on ED-1.0"
                                                   )#end names
                                       )#end list
   #----- Integration scheme. -------------------------------------------------------------#
   flagvar[["integ.scheme"]]    = list( descr  = "Integration scheme"
                                       , numeric = TRUE
                                       , values = c(0,1,2,3)
                                       , names  = c("Forward Euler"
                                                   ,"4th order Runge-Kutta"
                                                   ,"Heun"
                                                   ,"Hybrid"
                                                   )#end names
                                       )#end list
   #----- Drought scenario. ---------------------------------------------------------------#
   flagvar[["drought.ext"]]      = list( descr  = "Drought extent"
                                       , numeric = TRUE
                                       , values = c(0,1,2,3)
                                       , names  = c("Same as default"
                                                   ,"+1 month"
                                                   ,"+2 months"
                                                   ,"+3 months"
                                                   )#end names
                                       )#end list
   flagvar[["drought.seq"]]      = list( descr  = "Drought sequence"
                                       , numeric = TRUE
                                       , values = c(0,1,2)
                                       , names  = c("Same as default"
                                                   ,"Consecutive"
                                                   ,"Random"
                                                   )#end names
                                       )#end list
   flagvar[["isoilcol"]]         = list( descr  = "Soil colour"
                                       , numeric = TRUE
                                       , values = seq(1,21,1)
                                       , names  = c("01"
                                                   ,"02 - Bright"
                                                   ,"03"
                                                   ,"04"
                                                   ,"05"
                                                   ,"06 - Medium bright"
                                                   ,"07"
                                                   ,"08"
                                                   ,"09"
                                                   ,"10 - Medium"
                                                   ,"11"
                                                   ,"12"
                                                   ,"13"
                                                   ,"14 - Medium dark"
                                                   ,"15"
                                                   ,"16"
                                                   ,"17"
                                                   ,"18 - Dark"
                                                   ,"19"
                                                   ,"20 - Very dark"
                                                   ,"21 - ED-2.1"
                                                   )#end names
                                       )#end list
   flagvar[["iscenario"]]        = list( descr  = "Scenario"
                                       , numeric = TRUE
                                       , values = c(-150,-125,-100, -75, -50, -25
                                                   , -10,   0,  10
                                                   ,  25,  50,  75, 100, 125, 150
                                                   )#end c
                                       , names  = c("Drier:  A = 3/2"
                                                   ,"Drier:  A = 5/4"
                                                   ,"Drier:  A =   1"
                                                   ,"Drier:  A = 3/4"
                                                   ,"Drier:  A = 1/2"
                                                   ,"Drier:  A = 1/4"
                                                   ,"INMET period"
                                                   ,"EFT period only"
                                                   ,"Uniform sampling"
                                                   ,"Wetter: A = 1/4"
                                                   ,"Wetter: A = 1/2"
                                                   ,"Wetter: A = 3/4"
                                                   ,"Wetter: A =   1"
                                                   ,"Wetter: A = 5/4"
                                                   ,"Wetter: A = 3/2"
                                                   )#end names
                                       )#end list
   flagvar[["met.forcing"]]      = list( descr  = "Meteorological forcing"
                                       , numeric = FALSE
                                       , values = c("eft","shr","wmo","rag")
                                       , names  = c("Eddy flux tower"
                                                   ,"Sub-hourly"
                                                   ,"WMO-based"
                                                   ,"Average rainfall"
                                                   )#end names
                                       )#end list
   flagvar[["init.mode"]]        = list( descr  = "Initial Cond.:"
                                       , numeric = TRUE
                                       , values = c(-1,0,1,2,3,4,5,6)
                                       , names  = c("Bare Ground"
                                                   ,"Near bare ground"
                                                   ,"ED-1.0"
                                                   ,"ED-2.0"
                                                   ,"ED-2.0 + Hydro"
                                                   ,"ED-2.1"
                                                   ,"ED-2.2"
                                                   ,"Biometry"
                                                   )#end names
                                       )#end list
   flagvar[["struct"]]           = list( descr  = "Structure"
                                       , numeric = FALSE
                                       , values = c("ble_iage30_pft02"
                                                   ,"ble_iage30_pft05"
                                                   ,"sas_iage01_pft02"
                                                   ,"sas_iage01_pft05"
                                                   ,"sas_iage30_pft02"
                                                   ,"sas_iage30_pft05"
                                                   ,"ble_age30_pft02"
                                                   ,"ble_age30_pft05"
                                                   ,"sas_age01_pft02"
                                                   ,"sas_age01_pft05"
                                                   ,"sas_age30_pft02"
                                                   ,"sas_age30_pft05"
                                                   )#end values
                                       , names  = c("2 PFTs"
                                                   ,"5 PFTs"
                                                   ,"Size + 2 PFTs"
                                                   ,"Size + 5 PFTs"
                                                   ,"Size + Age + 2 PFTs"
                                                   ,"Size + Age + 5 PFTs"
                                                   ,"2 PFTs"
                                                   ,"5 PFTs"
                                                   ,"Size + 2 PFTs"
                                                   ,"Size + 5 PFTs"
                                                   ,"Size + Age + 2 PFTs"
                                                   ,"Size + Age + 5 PFTs"
                                                   )#end names
                                       )#end list
   flagvar[["iustar"]]           = list( descr  = "Friction Vel."
                                       , numeric = TRUE
                                       , values = c(0,1)
                                       , names  = c("Prescribed","Predicted")
                                       )#end list
   flagvar[["idimort"]]          = list( descr  = "DI Mortality"
                                       , numeric = TRUE
                                       , values = c(0,1)
                                       , names  = c("ED-1.0","ED-2.2")
                                       )#end list
   flagvar[["iplastic"]]         = list( descr  = "Plasticity"
                                       , numeric = TRUE
                                       , values = c(0,1,2)
                                       , names  = c("Off","height-based","LAI-based")
                                       )#end list
   flagvar[["repro.scheme" ]]    = list( descr  = "Seed allocation"
                                       , numeric = TRUE
                                       , values = c(0,1,2,3)
                                       , names  = c("No allocation"
                                                   ,"Big-bang within site"
                                                   ,"Big-bang across site"
                                                   ,"Asymptotic"
                                                   )#end c
                                       )#end list
   flagvar[["ianth.disturb"]]    = list( descr  = "Anthropogenic disturbance"
                                       , numeric = TRUE
                                       , values = c(0,1,2)
                                       , names  = c("OFF","ON","Logging")
                                       )#end list
   flagvar[["decomp.scheme"]]    = list( descr  = "Decomposition"
                                       , numeric = TRUE
                                       , values = c(0,1,2,3,4,5)
                                       , names  = c( "Moorcroft et al. (2001)"
                                                   , "Lloyd and Taylor (1994)"
                                                   , "Longo et al. (2019)"
                                                   , "Moorcroft+Moyano"
                                                   , "Lloyd_Taylor+Moyano"
                                                   , "Bolker et al. (1998)"
                                                   )#end c
                                       )#end list
   flagvar[["ihydro"       ]]    = list( descr  = "Plant hydraulics"
                                       , numeric = TRUE
                                       , values = c(0,1,2)
                                       , names  = c( "Static"
                                                   , "Xu/Christoffersen"
                                                   , "Xu et al (2016)"
                                                   )#end c
                                       )#end list
   flagvar[["hydrodyn.set" ]]    = list( descr  = "Model settings"
                                       , numeric = TRUE
                                       , values = c(0,1,2)
                                       , names  = c( "ED-2.2"
                                                   , "Hybrid"
                                                   , "Xu"
                                                   )#end c
                                       )#end list
   flagvar[["ianth.dataset"]]    = list( descr  = "LULCC dataset"
                                       , numeric = FALSE
                                       , values = c("lum-off"
                                                   ,"glu-331"
                                                   ,"glu-sa1"
                                                   ,"glu-sag"
                                                   ,"glu-sa2"
                                                   ,"lurcp26"
                                                   ,"lurcp45"
                                                   ,"lurcp60"
                                                   ,"lurcp85"
                                                   )#end c
                                       , names  = c("No LU model"
                                                   ,"GLM-3.3.1"
                                                   ,"GLM-3.3.1 + SimAmazonia1 (BAU)"
                                                   ,"GLM-3.3.1 + SimAmazonia1 (GOV)"
                                                   ,"GLM-3.3.1 + SimAmazonia2 (BAU)"
                                                   ,"LUH-1.1 - RCP2.6 (IMAGE)"
                                                   ,"LUH-1.1 - RCP4.5 (MiniCAM)"
                                                   ,"LUH-1.1 - RCP6.0 (AIM)"
                                                   ,"LUH-1.1 - RCP8.5 (MESSAGE)"
                                                   )#end c
                                       )#end list
   flagvar[["user"]]             = list( descr  = "Settings:"
                                       , numeric = FALSE
                                       , values = c("mlo","rgk","nml","kzh","als","mcd"
                                                   ,"dmm","prm")
                                       , names  = c("Marcos Longo"
                                                   ,"Ryan Knox"
                                                   ,"Naomi Levine"
                                                   ,"Ke Zhang"
                                                   ,"Abby Swann"
                                                   ,"Mike Dietze"
                                                   ,"David Medvigy"
                                                   ,"Paul Moorcroft"
                                                   )
                                       )#end list
   flagvar[["npatch"]]           = list( descr  = "Initial patch count"
                                       , numeric = TRUE
                                       , values = c(1,4,25,100,2500)
                                       , names  = c("1","4","25","100","2500")
                                       )#end list
   flagvar[["tf.slope"]]         = list( descr  = "Mortality slope"
                                       , numeric = FALSE
                                       , values = c("-","+")
                                       , names  = c("0.05","0.15")
                                       )#end list
   flagvar[["iszpft"  ]]         = list( descr  = "Initial condition"
                                       , numeric = TRUE
                                       , values = c(0,1,2)
                                       , names  = c("Filled","Inventory","Airborne Lidar")
                                       )#end list
   flagvar[["teff"  ]]           = list( descr  = "Temperature increase"
                                       , numeric = TRUE
                                       , values = c(0,1)
                                       , names  = c("OFF","ON")
                                       )#end list
   flagvar[["soil.hydro"]]       = list( descr   = "Soil hydraulics"
                                       , numeric = TRUE
                                       , values  = c(0,1,2)
                                       , names   = c( "Brooks-Corey-Cosby"
                                                    , "Brooks-Cory-Tomasella"
                                                    , "van Genuchten-Hodnett"
                                                    )#end c
                                       )#end list
   flagvar[["dhist" ]]           = list( descr  = "LCLU history"
                                       , numeric = FALSE
                                       , values = c("int","ril","cl1","cl2","lt1"
                                                   ,"bn1","bn2","bn3","bn6"
                                                   ,"lb1","lb2","blb"
                                                   ,"sec","sb2","lwr","upr")
                                       , names  = c("Intact"
                                                   ,"Reduced-Impact Logging"
                                                   ,"Conventional Logging"
                                                   ,"Logged twice"
                                                   ,"Logged and thinned"
                                                   ,"Burned once"
                                                   ,"Burned twice"
                                                   ,"Burned three times"
                                                   ,"Burned six times"
                                                   ,"Logged and burned once"
                                                   ,"Logged and burned twice"
                                                   ,"Burned, logged, burned"
                                                   ,"Secondary growth"
                                                   ,"Secondary and burned twice"
                                                   ,"Bottomland"
                                                   ,"Plateau"
                                                   )#end c
                                       )#end list
   flagvar[["sl.scale"     ]]    = list( descr  = "Logging scale"
                                       , numeric = TRUE
                                       , values = c(0,1)
                                       , names  = c("Working unit","Landscape")
                                       )#end list
   flagvar[["logging.type" ]]    = list( descr  = "Logging type"
                                       , numeric = FALSE
                                       , values = c("ril","cvl","scl")
                                       , names  = c("Reduced-Impact Logging"
                                                   ,"Conventional Logging"
                                                   ,"Stand clearing"
                                                   )#end c
                                       )#end list
   flagvar[["icompile" ]]        = list( descr  = "Compilation type"
                                       , numeric = FALSE
                                       , values = c("openmp","serial")
                                       , names  = c("OpenMP","Serial")
                                       )#end list
   flagvar[["version"  ]]        = list( descr  = "ED version"
                                       , numeric = FALSE
                                       , values = c("ed20","ed21","ed22")
                                       , names  = c("ED-2.0.12","ED-2.1","ED-2.2")
                                       )#end list
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     The following variables are numeric variables that may need conversion from       #
   # integer to numeric.                                                                   #
   #---------------------------------------------------------------------------------------#
   numvar = list()
   numvar[["kw.grass"       ]] = list( descr = "Kw[grass]"          
                                     , unit  = "m2/yr/kgC_rt"          
                                     , fmt   = "%.0f"          
                                     , off   =    0.0
                                     , mult  =    1.0  )
   numvar[["kw.tree"        ]] = list( descr = "Kw[tree]"                   
                                     , unit  = "m2/yr/kgC_rt"          
                                     , fmt   = "%.0f"          
                                     , off   =    0.0
                                     , mult  =    1.0  )
   numvar[["kw"             ]] = list( descr = "Kw"                     
                                     , unit  = "m2/yr/kgC_rt"          
                                     , fmt   = "%.0f"          
                                     , off   =    0.0
                                     , mult  =    1.0  )
   numvar[["kwfact"         ]] = list( descr = "Kw scale factor"
                                     , unit  = ""
                                     , fmt   = "%.3f"
                                     , off   =    0.0
                                     , mult  =    0.001  )
   numvar[["d0fact"         ]] = list( descr = "D0 scale factor"
                                     , unit  = ""
                                     , fmt   = "%.3f"
                                     , off   =    0.0
                                     , mult  =    0.001  )
   numvar[["fclump"         ]] = list( descr = "Clumping factor"
                                     , unit  = ""
                                     , fmt   = "%.3f"
                                     , off   =    0.0
                                     , mult  =    0.001  )
   numvar[["mphoto.trc3"    ]] = list( descr = "Stomatal slope (C3)"
                                     , unit  = ""
                                     , fmt   = "%.1f"
                                     , off   =    0.0
                                     , mult  =    0.1  )
   numvar[["soil.depth"     ]] = list( descr = "Soil depth"                 
                                     , unit  = "m"                     
                                     , fmt   = "%.1f"          
                                     , off   =    0.0
                                     , mult  =    0.01 )
   numvar[["atm.co2"        ]] = list( descr = "Air CO2"                    
                                     , unit  = "umol/mol"              
                                     , fmt   = "%.0f"          
                                     , off   =    0.0
                                     , mult  =    1.0  )
   numvar[["gamma.c3"       ]] = list( descr = "gamma[C3]"                  
                                     , unit  = ""                      
                                     , fmt   = "%.3f"          
                                     , off   =    0.0
                                     , mult  =    0.0001)
   numvar[["gamma.c4"       ]] = list( descr = "gamma[C4]"                  
                                     , unit  = ""                      
                                     , fmt   = "%.3f"          
                                     , off   =    0.0
                                     , mult  =    0.001)
   numvar[["vm0.c3"         ]] = list( descr = "Vm0fac[C3]"                 
                                     , unit  = ""                      
                                     , fmt   = "%.3f"          
                                     , off   =    0.0
                                     , mult  =    0.001 )
   numvar[["vm0.c4"         ]] = list( descr = "Vm0fac[C4]"                 
                                     , unit  = ""                      
                                     , fmt   = "%.3f"          
                                     , off   =    0.0
                                     , mult  =    0.001 )
   numvar[["vm0"            ]] = list( descr = "Vm0fac"                   
                                     , unit  = ""                      
                                     , fmt   = "%.1f"          
                                     , off   =    0.0
                                     , mult  =    0.01 )
   numvar[["leaf.resp"      ]] = list( descr = "Leaf Resp. Factor"                   
                                     , unit  = ""                      
                                     , fmt   = "%.4f"          
                                     , off   =    0.0
                                     , mult  =    0.0001 )
   numvar[["grow.resp"      ]] = list( descr = "Grow. Resp. Factor"                   
                                     , unit  = ""                      
                                     , fmt   = "%.3f"          
                                     , off   =    0.0
                                     , mult  =    0.001 )
   numvar[["root.resp"      ]] = list( descr = "Root Resp. Factor"                   
                                     , unit  = ""                      
                                     , fmt   = "%.3f"          
                                     , off   =    0.0
                                     , mult  =    0.001 )
   numvar[["lw.grass"       ]] = list( descr = "Leaf width [grass]"         
                                     , unit  = "cm"                    
                                     , fmt   = "%.1f"          
                                     , off   =    0.0
                                     , mult  =    1.0  )
   numvar[["lw.tree"        ]] = list( descr = "Leaf width [tree]"          
                                     , unit  = "cm"                    
                                     , fmt   = "%.1f"          
                                     , off   =    0.0
                                     , mult  =    1.0  )
   numvar[["orient.tree"    ]] = list( descr = "Orient fact.[tree]"         
                                     , unit  = ""                      
                                     , fmt   = "%.2f"          
                                     , off   =    0.0
                                     , mult  =    0.01 )
   numvar[["orient.grass"   ]] = list( descr = "Orient fact.[grass]"        
                                     , unit  = ""                      
                                     , fmt   = "%.2f"          
                                     , off   =    0.0
                                     , mult  =    0.01 )
   numvar[["klowco2"        ]] = list( descr = "Klow"                       
                                     , unit  = "mol/umol"                    
                                     , fmt   = "%.0f"          
                                     , off   =    0.0
                                     , mult  = 1000.   )
   numvar[["d0.c3"          ]] = list( descr = "D0[c3]"                     
                                     , unit  = ""                      
                                     , fmt   = "%.3f"          
                                     , off   =    0.0
                                     , mult  =    0.001)
   numvar[["d0"             ]] = list( descr = "D0"
                                     , unit  = "mol/mol"
                                     , fmt   = "%.3f"
                                     , off   =    0.0
                                     , mult  =    0.001)
   numvar[["ltover"         ]] = list( descr = "Leaf turnover"
                                     , unit  = "1/yr"
                                     , fmt   = "%.1f"
                                     , off   =    0.0
                                     , mult  =    0.1)
   numvar[["runoff"         ]] = list( descr = "Runoff time"
                                     , unit  = "h"
                                     , fmt   = "%.0f"
                                     , off   =    0.0
                                     , mult  = 3600.0)
   numvar[["leaf.drop"      ]] = list( descr = "Leaf drop threshold"
                                     , unit  = "m"
                                     , fmt   = "%.0f"
                                     , off   =    0.0
                                     , mult  =    1.0)
   numvar[["sm.fire"        ]] = list( descr = "Fire threshold"
                                     , unit  = "m"
                                     , fmt   = "%.2f"
                                     , off   =    0.0
                                     , mult  =    0.01)
   numvar[["fire.parameter" ]] = list( descr = "Fire parameter"
                                     , unit  = "--"
                                     , fmt   = "%.2f"
                                     , off   =    0.0
                                     , mult  =    0.1)
   numvar[["mslope"         ]] = list( descr = "Stomatal slope"
                                     , unit  = "--"
                                     , fmt   = "%.1f"
                                     , off   =    0.0
                                     , mult  =    0.1)
   numvar[["ubmin"          ]] = list( descr = "Min. Wind"
                                     , unit  = "m/s"
                                     , fmt   = "%.2f"
                                     , off   =    0.0
                                     , mult  =    0.01)
   numvar[["ustmin"         ]] = list( descr = "Min. u*"
                                     , unit  = "m/s"
                                     , fmt   = "%.2f"
                                     , off   =    0.0
                                     , mult  =    0.01)
   numvar[["ddmort.fac"     ]] = list( descr = "DD mortality factor"
                                     , unit  = "--"
                                     , fmt   = "%.2f"
                                     , off   =    0.0
                                     , mult  =    1.00)
   numvar[["sand"           ]] = list( descr = "Sand"
                                     , unit  = "%"
                                     , fmt   = "%.1f"
                                     , off   =    0.0
                                     , mult  =    0.1)
   numvar[["clay"           ]] = list( descr = "Clay"
                                     , unit  = "%"
                                     , fmt   = "%.1f"
                                     , off   =    0.0
                                     , mult  =    0.1)
   numvar[["treefall"       ]] = list( descr = "TF Disturbance"
                                     , unit  = "%/yr"
                                     , fmt   = "%.3f"
                                     , off   =    0.0
                                     , mult  =   0.01)
   numvar[["leaf.absorb.nir"]] = list( descr = "Leaf NIR absorptance"
                                     , unit  = ""
                                     , fmt   = "%.3f"
                                     , off   =    0.0
                                     , mult  =   0.001)
   numvar[["leaf.absorb.vis"]] = list( descr = "Leaf PAR absorptance"
                                     , unit  = ""
                                     , fmt   = "%.3f"
                                     , off   =    0.0
                                     , mult  =   0.001)
   numvar[["idrain.scen"]]     = list( descr = "Change in rainfall"
                                     , unit  = "s.d."
                                     , fmt   = "%.2f"
                                     , off   =    0.0
                                     , mult  =   0.01)
   numvar[["idtemp.scen"]]     = list( descr = "Change in temperature"
                                     , unit  = "K"
                                     , fmt   = "%.2f"
                                     , off   =    0.0
                                     , mult  =   0.01)
   numvar[["realisation"]]     = list( descr = "Realisation"
                                     , unit  = ""
                                     , fmt   = "%2.2i"
                                     , off   =    0.0
                                     , mult  =    1.0)
   numvar[["iage"]]            = list( descr = "Target # of patches"
                                     , unit  = ""
                                     , fmt   = "%3i"
                                     , off   =   0.0
                                     , mult  =   1.0)
   numvar[["yeara"]]           = list( descr = "Initial year"
                                     , unit  = ""
                                     , fmt   = "%4.4i"
                                     , off   =   0.0
                                     , mult  =   1.0)
   numvar[["yearz"]]           = list( descr = "Final year"
                                     , unit  = ""
                                     , fmt   = "%4.4i"
                                     , off   =   0.0
                                     , mult  =   1.0)
   numvar[["lon"  ]]           = list( descr = "Longitude"
                                     , unit  = ""
                                     , fmt   = "%.2f"
                                     , off   =   0.0
                                     , mult  =   1.0)
   numvar[["lat"  ]]           = list( descr = "Latitude"
                                     , unit  = ""
                                     , fmt   = "%.2f"
                                     , off   =   0.0
                                     , mult  =   1.0)
   numvar[["nzs"  ]]           = list( descr = "Snow layers"
                                     , unit  = ""
                                     , fmt   = "%3i"
                                     , off   =   0.0
                                     , mult  =   1.0)
   numvar[["rcp"  ]]           = list( descr = "RCP"
                                     , unit  = "W/m2"
                                     , fmt   = "%.1f"
                                     , off   =   0.0
                                     , mult  =   0.1)
   numvar[["drain"]]           = list( descr = "Rainfall change"
                                     , unit  = "%"
                                     , fmt   = "%.1f"
                                     , off   =   0.0
                                     , mult  =   0.1)
   numvar[["bharvest"]]        = list( descr = "Target harvest"
                                     , unit  = "kgC/m2"
                                     , fmt   = "%.3f"
                                     , off   = 0.0
                                     , mult  = 0.001)
   numvar[["logging.cycle"]]   = list( descr = "Logging cycle"
                                     , unit  = "yrs"
                                     , fmt   = "%4i"
                                     , off   = 0.0
                                     , mult  = 1.0)
   numvar[["topsoil" ]]        = list( descr = "Top soil thickness"
                                     , unit  = "cm"
                                     , fmt   = "%.1f"
                                     , off   = 0.0
                                     , mult  = -0.1)
   numvar[["rrffact" ]]        = list( descr = "Root respiration factor"
                                     , unit  = ""
                                     , fmt   = "%.3f"
                                     , off   = 0.0
                                     , mult  = 0.001)
   numvar[["theta.crit"]]      = list( descr = "Leaf drop"
                                     , unit  = "MPa"
                                     , fmt   = "%.2f"
                                     , off   =    0.0
                                     , mult  =    0.01)
   numvar[["simul"   ]]        = list( descr = "Simulation"
                                     , unit  = ""
                                     , fmt   = "%4.4i"
                                     , off   = 0.0
                                     , mult  = 1.0)
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #         The following list is currently hardcoded, it could be better designed in the #
   # future so the script would recognise the terms automatically.                         #
   #---------------------------------------------------------------------------------------#
   if (iata){
      lenici = nchar(ici)
      if (lenici == 4){
         nparms = 0
         params = character(0)
         na     = integer(0)
         nz     = integer(0)
      }else if (lenici == 8){
         nparms = 1
         param  = c("dhist")
         na     = c(      6)
         nz     = c(      8)
      }else if (lenici == 9 && grepl(pattern="_ed2",x=ici)){
         nparms = 1
         param  = c("version")
         na     = c(      6)
         nz     = c(      9)
      }else if (lenici == 9){
         nparms = 1
         param  = c("simul")
         na     = c(      7)
         nz     = c(      9)
      }else if (lenici == 10){
         nparms = 1
         param  = c("nzs")
         na     = c(    9)
         nz     = c(   10)
      }else if (lenici == 11){
         nparms = 1
         param  = c("icompile")
         na     = c(         6)
         nz     = c(        11)
      }else if (lenici == 12){
         nparms = 2
         param  = c("met.forcing","isas")
         na     = c(            6,    10)
         nz     = c(            8,    12)
      }else if (lenici == 13 && grepl(pattern="irepro",x=ici)){
         nparms = 1
         param  = c("repro.scheme")
         na     = c(            12)
         nz     = c(            13)
      }else if (lenici == 13 && grepl(pattern="iphen",x=ici)){
         nparms = 1
         param  = c("iphen.scheme")
         na     = c(            11)
         nz     = c(            13)
      }else if (lenici == 13 && grepl(pattern="istext",x=ici)){
         nparms = 1
         param  = c("isoil.text")
         na     = c(          12)
         nz     = c(          13)
      }else if (lenici == 13 && grepl(pattern="tfall",x=ici)){
         nparms = 1
         param  = c("treefall")
         na     = c(          11)
         nz     = c(          13)
      }else if (lenici == 14 && grepl(pattern="cotwo",x=ici)){
         nparms = 1
         param  = c("atm.co2")
         na     = c(       11)
         nz     = c(       14)
      }else if (lenici == 14 && grepl(pattern="ibranch",x=ici)){
         nparms = 1
         param  = c("ibranch.thermo")
         na     = c(              13)
         nz     = c(              14)
      }else if (lenici == 14 && grepl(pattern="mslope",x=ici)){
         nparms = 1
         param  = c("mslope")
         na     = c(        12)
         nz     = c(        14)
      }else if (lenici == 15 && grepl(pattern="icanrad",x=ici)){
         nparms = 1
         param  = c("icanrad")
         na     = c(        13)
         nz     = c(        15)
      }else if (lenici == 15 && grepl(pattern="iplastic",x=ici)){
         nparms = 1
         param  = c("iplastic")
         na     = c(        14)
         nz     = c(        15)
      }else if (lenici == 15 && grepl(pattern="ivegtdyn",x=ici)){
         nparms = 1
         param  = c("ivegt.dynamics")
         na     = c(        14)
         nz     = c(        15)
      }else if (lenici == 15 && grepl(pattern="iphysiol",x=ici)){
         nparms = 1
         param  = c("iphysiol")
         na     = c(        14)
         nz     = c(        15)
      }else if (lenici == 15 && grepl(pattern="icanturb",x=ici)){
         nparms = 1
         param  = c("icanturb")
         na     = c(        14)
         nz     = c(        15)
      }else if (lenici == 15 && grepl(pattern="h2olimit",x=ici)){
         nparms = 1
         param  = c("h2o.plant.limit")
         na     = c(               14)
         nz     = c(               15)
      }else if (lenici == 15 && grepl(pattern="kwfact",x=ici)){
         nparms = 1
         param  = c("kwfact")
         na     = c(      12)
         nz     = c(      15)
      }else if (lenici == 15 && grepl(pattern="fclump",x=ici)){
         nparms = 1
         param  = c("fclump")
         na     = c(      12)
         nz     = c(      15)
      }else if (lenici == 15 && grepl(pattern="ustmin",x=ici)){
         nparms = 1
         param  = c("ustmin")
         na     = c(      12)
         nz     = c(      15)
      }else if (lenici == 16 && grepl(pattern="rrffact",x=ici)){
         nparms = 1
         param  = c("rrffact")
         na     = c(       13)
         nz     = c(       16)
      }else if (lenici == 16 && grepl(pattern="ifire",x=ici)){
         nparms = 2
         param  = c("dhist","include.fire")
         na     = c(      6,            15)
         nz     = c(      8,            16)
      }else if (lenici == 17 && grepl(pattern="iphen",x=ici)){
         nparms = 2
         param  = c("dhist","iphen.scheme")
         na     = c(      6,            15)
         nz     = c(      8,            17)
      }else if (lenici == 17){
         nparms = 2
         param  = c("met.forcing","isoil.text")
         na     = c(            6,          16)
         nz     = c(            8,          17)
      }else if (lenici == 18){
         nparms = 2
         param  = c("tf.slope","treefall")
         na     = c(        14,        15)
         nz     = c(        14,        18)
      }else if (lenici == 19){
         nparms = 3
         param  = c("met.forcing","isas","iage")
         na     = c(            6,    10,    18)
         nz     = c(            8,    12,    19)
      }else if (lenici == 20){
         nparms = 2
         param  = c("include.fire","isoil.text")
         na     = c(            11,          19)
         nz     = c(            12,          20)
      }else if (lenici == 21){
         nparms = 3
         param  = c("isas","iage","idiversity")
         na     = c(     6,    14,          20)
         nz     = c(     8,    15,          21)
      }else if (lenici == 22 && grepl(pattern="iszpft",x=ici)){
         nparms = 2
         param  = c("iszpft","iphen.scheme")
         na     = c(      12,            20)
         nz     = c(      13,            22)
      }else if (lenici == 22 && grepl(pattern="istext",x=ici)){
         nparms = 2
         param  = c("isoil.text","iphen.scheme")
         na     = c(          12,            20)
         nz     = c(          13,            22)
      }else if (lenici == 22 && grepl(pattern="iallom",x=ici)){
         nparms = 2
         param  = c("iallom","igrass")
         na     = c(      12,      21)
         nz     = c(      13,      22)
      }else if (lenici == 22 && grepl(pattern="iustar",x=ici)){
         nparms = 2
         param  = c("iustar","icanturb")
         na     = c(      14,      21)
         nz     = c(      15,      22)
      }else if (lenici == 22 && grepl(pattern="ihrzrad",x=ici)){
         nparms = 2
         param  = c("ihrzrad","include.fire")
         na     = c(       13,            21)
         nz     = c(       14,            22)
      }else if (lenici == 23 && grepl(pattern="iustar",x=ici)){
         nparms = 2
         param  = c("iustar","icanturb")
         na     = c(      12,        22)
         nz     = c(      13,        23)
      }else if (lenici == 23 && grepl(pattern="idecomp",x=ici)){
         nparms = 2
         param  = c("decomp.scheme","iphen.scheme")
         na     = c(             13,            21)
         nz     = c(             14,            23)
      }else if (lenici == 24 && grepl(pattern="ihydro",x=ici)){
         nparms = 2
         param  = c("ihydro","ivegt.dynamics")
         na     = c(      12,              23)
         nz     = c(      13,              24)
      }else if (lenici == 24 && grepl(pattern="icanturb",x=ici)){
         nparms = 2
         param  = c("iphen.scheme","icanturb")
         na     = c(            11,        23)
         nz     = c(            13,        24)
      }else if (lenici == 24 && grepl(pattern="thcrit",x=ici)){
         nparms = 2
         param  = c("iphen.scheme","theta.crit")
         na     = c(            11,          21)
         nz     = c(            13,          24)
      }else if (lenici == 24 && grepl(pattern="ihrzrad",x=ici)){
         nparms = 2
         param  = c("ihrzrad","imetrad")
         na     = c(       13,       23)
         nz     = c(       14,       24)
      }else if (lenici == 24 && grepl(pattern="iallom",x=ici)){
         nparms = 2
         param  = c("iallom","iphysiol")
         na     = c(      12,      23)
         nz     = c(      13,      24)
      }else if (lenici == 25 && grepl(pattern="cotwo",x=ici)){
         nparms = 2
         param  = c("soil.depth","atm.co2")
         na     = c(          12,       22)
         nz     = c(          15,       25)
      }else if (lenici == 25){
         nparms = 3
         param  = c("iphen.scheme","d0","include.fire")
         na     = c(            10,  16,            24)
         nz     = c(            12,  18,            25)
      }else if (lenici == 26 & grep(pattern="islhydro",x=ici)){
         nparms = 2
         param  = c("ivegt.dynamics","soil.hydro")
         na     = c(              14,          25)
         nz     = c(              15,          26)
      }else if (lenici == 26){
         nparms = 3
         param  = c("ianth.disturb","logging.type", "bharvest")
         na     = c(             11,            14,         23)
         nz     = c(             12,            16,         26)
      }else if (lenici == 27 && grepl(pattern="ihydrodyn",x=ici)){
         nparms = 2
         param  = c("ivegt.dynamics","hydrodyn.set")
         na     = c(              14,            26)
         nz     = c(              15,            27)
      }else if (lenici == 27 && grepl(pattern="idecomp",x=ici)){
         nparms = 3
         param  = c("met.forcing","decomp.scheme","iphen.scheme")
         na     = c(            6,             17,            25)
         nz     = c(            8,             18,            27)
      }else if (lenici == 27 && grepl(pattern="ihrzrad",x=ici)){
         nparms = 2
         param  = c("ivegt.dynamics","ihrzrad")
         na     = c(             16,        26)
         nz     = c(             17,        27)
      }else if (lenici == 28){
         nparms = 3
         param  = c("iphen.scheme", "isoil.text","treefall")
         na     = c(            11,           20,        25)
         nz     = c(            13,           21,        28)
      }else if (lenici == 29 && grepl(pattern="d0x",x=ici)){
         nparms = 3
         param  = c("iphen.scheme","kwfact","d0fact")
         na     = c(            11,      18,      26)
         nz     = c(            13,      21,      29)
      }else if (lenici == 29){
         nparms = 3
         param  = c("ivegt.dynamics","ihrzrad","leaf.absorb.vis")
         na     = c(             11,        18,               26)
         nz     = c(             12,        19,               29)
      }else if (lenici == 30 && grepl(pattern="kwx",x=ici)){
         nparms = 3
         param  = c("iphen.scheme","kwfact","mphoto.trc3")
         na     = c(            11,      18,           29)
         nz     = c(            13,      21,           30)
      }else if (lenici == 30){
         nparms = 3
         param  = c("include.fire","isoil.text","treefall")
         na     = c(            11,          19,        27)
         nz     = c(            12,          20,        30)
      }else if (lenici == 31){
         nparms = 3
         param  = c("iallom","idimort","include.fire")
         na     = c(      12,       22,            30)
         nz     = c(      13,       23,            31)
      }else if (lenici == 32 && grepl(pattern="iecon",x=ici)){
         nparms = 3
         param  = c("iallom","ieconomics","iplastic")
         na     = c(      12,          20,        31)
         nz     = c(      13,          21,        32)
      }else if (lenici == 32 && grepl(pattern="imetrad",x=ici)){
         nparms = 3
         param  = c("iage","ihrzrad","imetrad")
         na     = c(    10,       21,       31)
         nz     = c(    12,       22,       32)
      }else if (lenici == 34 && grepl(pattern="imetrad",x=ici)){
         nparms = 3
         param  = c("ihrzrad","imetrad","vm0")
         na     = c(       13,       23,   32)
         nz     = c(       14,       24,   34)
      }else if (lenici == 34 && grepl(pattern="bharv",x=ici)){
         nparms = 4
         param  = c("ianth.disturb","logging.type","logging.cycle","bharvest")
         na     = c(             11,            14,             22,        31)
         nz     = c(             12,            16,             24,        34)
      }else if (lenici == 35){
         nparms = 4
         param  = c("include.fire","ianth.disturb","ianth.dataset","isoil.text")
         na     = c(            11,             19,             22,          34)
         nz     = c(            12,             20,             28,          35)
      }else if (lenici == 36 && grepl(pattern="idecomp",x=ici)){
         nparms = 3
         param  = c("ivegt.dynamics","decomp.scheme","hydrodyn.set")
         na     = c(              14,             24,            35)
         nz     = c(              15,             25,            36)
      }else if (lenici == 36){
         nparms = 3
         param  = c("iustar","icanturb","topsoil")
         na     = c(      12,        22,       32)
         nz     = c(      13,        23,       36)
      }else if (lenici == 37){
         nparms = 4
         param  = c("yeara","iphen.scheme","isoil.text","treefall")
         na     = c(      9,            18,          27,        35)
         nz     = c(     11,            20,          28,        37)
      }else if (lenici == 39){
         nparms = 3
         param  = c("struct","iphen.scheme","include.fire")
         na     = c(      10,            30,            38)
         nz     = c(      24,            32,            39)
      }else if (lenici == 40){
         nparms = 4
         param  = c("grow.resp","vm0","leaf.resp","root.resp")
         na     = c(         10,   18,         28,         37)
         nz     = c(         13,   21,         30,         40)
      }else if (lenici == 41 && grepl(pattern="stext",x=ici)){
         nparms = 5
         param  = c("idrain.scen","idtemp.scen","realisation","iphen.scheme","isoil.text")
         na     = c(            7,           13,           23,            31,          40)
         nz     = c(           10,           16,           24,            33,          41)
      }else if (lenici == 41 && grepl(pattern="iecon",x=ici)){
         nparms = 4
         param  = c("iallom","ieconomics","repro.scheme","iplastic")
         na     = c(      12,          20,           29,         40)
         nz     = c(      13,          21,           30,         41)
      }else if (lenici == 50 && grepl(pattern="bharv",x=ici)){
         nparms = 6
         param  = c("ianth.disturb","logging.type",    "sl.scale"
                   ,"logging.cycle",    "bharvest","include.fire")
         na     = c(             11,            14,            23
                   ,             30,            39,            49)
         nz     = c(             12,            16,            24
                   ,             32,            42,            50)
      }else if (lenici == 50){
         nparms = 5
         param  = c("idrain.scen","idtemp.scen","realisation","iphen.scheme","struct")
         na     = c(            7,           13,           23,            31,      35)
         nz     = c(           10,           16,           24,            33,      50)
      }#end if
   }else{
      lenici = nchar(ici)
      if (lenici == 28){
         nparms = 2
         param  = c("lon","lat")
         na     = c(   13,   23)
         nz     = c(   18,   28)
      }#end if
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Loop over all parameters and concatenate the description.                         #
   #---------------------------------------------------------------------------------------#
   description = testpoi
   nchar.line  = nchar(description)
   for (p in sequence(nparms)){
      #----- Retrieve the parameter value or flag. ----------------------------------------#
      mycharval  = substring(ici,na[p],nz[p])
      myval      = as.numeric(mycharval)
      #----- Check whether the parameter is a value or a flag. ----------------------------#
      mytest = param[p]
      if (mytest %in% names(flagvar)){
         #----- Flag variable, find the number that corresponds to the setting. -----------#
         if (flagvar[[mytest]]$numeric){
            n = match(myval,flagvar[[mytest]]$values)
         }else{
            n = match(mycharval,flagvar[[mytest]]$values)
         }#end if
         if (is.na(n)){
            browser()
            stop(paste0(" Option ",myval," doesn't exist for",flagvar[[mytest]]$descr,"!"))
         }else{
            thisdesc = paste0(flagvar[[mytest]]$descr,": ",flagvar[[mytest]]$names[n])
         }
      }else if (mytest %in% names(numvar)){
         #----- Numeric variable, put the value in a nicer format. ------------------------#
         val.pretty = sprintf(numvar[[mytest]]$fmt
                             ,numvar[[mytest]]$off + myval*numvar[[mytest]]$mult)

         thisdesc = paste0(numvar[[mytest]]$descr," = ",val.pretty,numvar[[mytest]]$unit)
      }else{
         stop(paste0(" Invalid parameter test: ",mytest,"!"))
      }#end if
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Check whether to break the line or continue in the same one.                   #
      #------------------------------------------------------------------------------------#
      if (nchar.line + nchar(thisdesc) > max.char){
         description = paste(description,thisdesc,sep="\n")
         nchar.line  = nchar(thisdesc)
      }else{
         description = paste(description,thisdesc,sep="    -    ")
         nchar.line  = nchar.line + nchar(thisdesc)
      }#end if
   }#end for
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Return the description of this simulation.                                        #
   #---------------------------------------------------------------------------------------#
   return(description)
   #---------------------------------------------------------------------------------------#
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     File poilist.csv replaces the original list here, to make it more manageable.  Feel  #
# free to add more sites in that file.                                                     #
#------------------------------------------------------------------------------------------#
poilist <<- read.csv( file             = file.path(srcdir,"poilist.csv")
                    , header           = TRUE
                    , stringsAsFactors = FALSE
                    )#end read.csv
npoi    <<- nrow(poilist)
#==========================================================================================#
#==========================================================================================#
