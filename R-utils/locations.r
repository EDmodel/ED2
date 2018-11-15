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
   }else if(nchar(ici) == 5 & substring(ici,1,2) == "kz")     {
      #----- Ke's list of polygons. -------------------------------------------------------#
       testpoi  = as.numeric(substring(ici,3,5))

       lieu     = paste0("KZ test - polygon ",testpoi)
       iata     = ici
       wmo      = NA
       pathroot = file.path(here,ici)
       pathin   = file.path(pathroot,"analy",ici)
       pathrst  = file.path(here,ici,"histo",ici)
       pathout  = file.path(pathroot,"epost")
       lon      = kzlist[[testpoi]]$lon
       lat      = kzlist[[testpoi]]$lat

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
                                      , values  = seq(from=0,to=1,by=1)
                                      , names   = c("OFF","ON")
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
   #----- Density-dependent mortality. ----------------------------------------------------#
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
                                       , values = c("eft","shr","wmo")
                                       , names  = c("Eddy flux tower"
                                                   ,"Sub-hourly"
                                                   ,"WMO-based"
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
   flagvar[["iunder"]]           = list( descr  = "Initial condition"
                                       , numeric = TRUE
                                       , values = c(0,1)
                                       , names  = c("Inventory","Equilibrium")
                                       )#end list
   flagvar[["teff"  ]]           = list( descr  = "Temperature increase"
                                       , numeric = TRUE
                                       , values = c(0,1)
                                       , names  = c("OFF","ON")
                                       )#end list
   flagvar[["dhist" ]]           = list( descr  = "LCLU history"
                                       , numeric = FALSE
                                       , values = c("int","ril","cl1","cl2","bn1","bn2"
                                                   ,"bn3","bn6","lb1","lb2","blb"
                                                   ,"sec","sb2","lwr","upr")
                                       , names  = c("Intact"
                                                   ,"Reduced-Impact Logging"
                                                   ,"Conventional Logging"
                                                   ,"Logged twice"
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
   numvar[["fclump"         ]] = list( descr = "Clumping factor"
                                     , unit  = ""
                                     , fmt   = "%.3f"
                                     , off   =    0.0
                                     , mult  =    0.001  )
   numvar[["soil.depth"     ]] = list( descr = "Soil depth"                 
                                     , unit  = "m"                     
                                     , fmt   = "%.1f"          
                                     , off   =    0.0
                                     , mult  =    0.1  )
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
         nparms = 1
         param  = c("iphen")
         na     = c(     15)
         nz     = c(     17)
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
      }else if (lenici == 23){
         nparms = 2
         param  = c("iustar","icanturb")
         na     = c(      12,        22)
         nz     = c(      13,        23)
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
      }else if (lenici == 25){
         nparms = 3
         param  = c("iphen.scheme","d0","include.fire")
         na     = c(            10,  16,            24)
         nz     = c(            12,  18,            25)
      }else if (lenici == 26){
         nparms = 3
         param  = c("ianth.disturb","logging.type", "bharvest")
         na     = c(             11,            14,         23)
         nz     = c(             12,            16,         26)
      }else if (lenici == 27){
         nparms = 2
         param  = c("ivegt.dynamics","ihrzrad")
         na     = c(             16,        26)
         nz     = c(             17,        27)
      }else if (lenici == 28){
         nparms = 3
         param  = c("iphen.scheme", "isoil.text","treefall")
         na     = c(            11,           20,        25)
         nz     = c(            13,           21,        28)
      }else if (lenici == 29){
         nparms = 2
         param  = c("ivegt.dynamics","ihrzrad","leaf.absorb.vis")
         na     = c(             11,        18,               26)
         nz     = c(             12,        19,               29)
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
#     This list has the commonest locations, feel free to add new places.                  #
#------------------------------------------------------------------------------------------#
u           = 0
poitmp      = list()
u           = u + 1
poitmp[[u]] = list( short           = "allpahuayo"        
                  , longname        = "Allpahuayo, PER"             
                  , iata            = "alp"
                  , lon             = -73.437
                  , lat             =  -3.953
                  , alt             = 143.
                  , wmo             = NA
                  , isoilflg        = 2
                  , ntext           = 1
                  , sand            = 0.937
                  , clay            = 0.026
                  , depth           = "D"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "alta_floresta"     
                  , longname        = "Alta Floresta, MT"           
                  , iata            = "afl"
                  , lon             = -56.100
                  , lat             =  -9.867
                  , alt             = 269.
                  , wmo             = 82965
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "andiroba"    
                  , longname        = "Fazenda Andiroba, PA"
                  , iata            = "and"
                  , lon             = -47.517   # -46.832
                  , lat             =  -3.315   # -2.551
                  , alt             = 122       # 105.
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "H"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = +2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "angra_dos_reis"    
                  , longname        = "Angra dos Reis, RJ"          
                  , iata            = "aei"
                  , lon             = -44.300
                  , lat             = -22.970
                  , alt             = 4.
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "araguaiana"        
                  , longname        = "Araguaiana, MT"              
                  , iata            = "ayx"
                  , lon             = -51.810
                  , lat             = -15.710
                  , alt             = 284.
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "araracuara"        
                  , longname        = "Araracuara, COL"             
                  , iata            = "arc"
                  , lon             = -72.398
                  , lat             =  -0.601
                  , alt             = 243
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "asuncion"          
                  , longname        = "Asuncion, PRY"               
                  , iata            = "asu"
                  , lon             = -57.560
                  , lat             = -25.300
                  , alt             = 127
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "bananal"           
                  , longname        = "Bananal Island, TO"          
                  , iata            = "ban"
                  , lon             = -50.159
                  , lat             =  -9.824
                  , alt             = 173
                  , wmo             = NA
                  , isoilflg        = 2
                  , ntext           = 8
                  , sand            = 0.240
                  , clay            = 0.370
                  , depth           = "C"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "Bananal"
                  , yeara           = 1999
                  , yearz           = 2007
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "barro_colorado"    
                  , longname        = "Barro Colorado, PAN"         
                  , iata            = "bci"
                  , lon             = -79.850
                  , lat             =   9.160
                  , alt             = 145
                  , wmo             = NA
                  , isoilflg        = 2
                  , ntext           = 17
                  , sand            = 0.200
                  , clay            = 0.420
                  , depth           = "D"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "belem"             
                  , longname        = "Belem, PA"                   
                  , iata            = "bel"
                  , lon             = -48.480
                  , lat             =  -1.380
                  , alt             = 15
                  , wmo             = 82193
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "belo_horizonte"    
                  , longname        = "Belo Horizonte, MG"          
                  , iata            = "cnf"
                  , lon             = -43.950
                  , lat             = -19.850
                  , alt             = NA
                  , wmo             = 83566
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "belterra"             
                  , longname        = "Belterra, PA"                   
                  , iata            = "bte"
                  , lon             = -54.944
                  , lat             =  -2.646
                  , alt             = 171
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "benjamin_constant"             
                  , longname        = "Benjamin Constant, AM"                   
                  , iata            = "qav"
                  , lon             = -70.030
                  , lat             =  -4.380
                  , alt             = 65
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "blumenau"          
                  , longname        = "Blumenau, SC"                
                  , iata            = "bnu"
                  , lon             = -49.060
                  , lat             = -26.920
                  , alt             = 15
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "boa_vista"         
                  , longname        = "Boa Vista, RR"               
                  , iata            = "bvb"
                  , lon             = -60.610
                  , lat             =   2.920
                  , alt             = 82
                  , wmo             = 82022
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "bogota"            
                  , longname        = "Bogota, COL"                 
                  , iata            = "bog"
                  , lon             = -74.100
                  , lat             =   4.650
                  , alt             = 2556
                  , wmo             = 80222
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "bom_jesus"          
                  , longname        = "Bom Jesus da Lapa, BA"                
                  , iata            = "laz"
                  , lon             = -43.417
                  , lat             = -13.267
                  , alt             = 440
                  , wmo             = 83288
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "bonal"    
                  , longname        = "Fazenda Bonal, AC"
                  , iata            = "bon"
                  , lon             = -67.286
                  , lat             = -9.874
                  , alt             = 165.
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "H"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = +2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "brasilia"          
                  , longname        = "Brasilia, DF"                
                  , iata            = "bsb"
                  , lon             = -47.713
                  , lat             = -15.601
                  , alt             = 1023
                  , wmo             = 83378
                  , isoilflg        = 2
                  , ntext           = 17
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "H"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "Brasilia"
                  , yeara           = 2006
                  , yearz           = 2012
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "bridgetown"          
                  , longname        = "Bridgetown, BRB"                
                  , iata            = "bgi"
                  , lon             = -59.500
                  , lat             =  13.067
                  , alt             = 56
                  , wmo             = 78954
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "Sheffield"
                  , yeara           = 1969
                  , yearz           = 2008
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "cabo_frio"         
                  , longname        = "Cabo Frio, RJ"               
                  , iata            = "cfb"
                  , lon             = -42.070
                  , lat             = -22.920
                  , alt             = 6
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "cacoal_grande"
                  , longname        = "Cacoal Grande (Embrapa), PA"
                  , iata            = "ecc"
                  , lon             = -54.329
                  , lat             =  -2.389
                  , alt             = 11
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "caico"        
                  , longname        = "Caico, RN"              
                  , iata            = "cqo"
                  , lon             = -37.085
                  , lat             =  -6.467
                  , alt             = 176
                  , wmo             = 82690
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "Caico"
                  , yeara           = 1998
                  , yearz           = 2008
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "cajazeiras"        
                  , longname        = "Cajazeiras, PB"              
                  , iata            = "cjz"
                  , lon             = -38.570
                  , lat             =  -6.900
                  , alt             = 306
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "calabozo"          
                  , longname        = "Calabozo, VEN"               
                  , iata            = "clz"
                  , lon             = -67.420
                  , lat             =   8.920
                  , alt             = 109
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "campo_grande"          
                  , longname        = "Campo Grande, MS"                
                  , iata            = "cgr"
                  , lon             = -54.538
                  , lat             = -20.438
                  , alt             = 677
                  , wmo             = 83612
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "Campo_Grande"
                  , yeara           = 2002
                  , yearz           = 2013
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "canarana"          
                  , longname        = "Canarana, MT"                
                  , iata            = "qnr"
                  , lon             = -52.250
                  , lat             = -13.560
                  , alt             = 395
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "carajas"           
                  , longname        = "Parauapebas (Carajas), PA"   
                  , iata            = "cks"
                  , lon             = -50.722
                  , lat             =  -5.786
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "cardoso"           
                  , longname        = "Cardoso Island, SP"          
                  , iata            = "czi"
                  , lon             = -48.010
                  , lat             = -25.096
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 2
                  , ntext           = 1
                  , sand            = 0.950
                  , clay            = 0.010
                  , depth           = "C"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "carolina"           
                  , longname        = "Carolina, MA"          
                  , iata            = "cln"
                  , lon             = -47.467
                  , lat             =  -7.333
                  , alt             = 193.
                  , wmo             = 82765
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "cauaxi"
                  , longname        = "Fazenda Cauaxi, PA"
                  , iata            = "cau"
                  , lon             = -47.517 # -48.483
                  , lat             =  -3.315 #  -3.748
                  , alt             = 122     # 140
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.00
                  , clay            = -1.00
                  , depth           = "H"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = +2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "cauaxi_und"
                  , longname        = "Faz. Cauaxi (Undist.), PA"
                  , iata            = "ca3"
                  , lon             = -48.483
                  , lat             =  -3.748
                  , alt             = 140
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.00
                  , clay            = -1.00
                  , depth           = "H"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = +2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "cauaxi_ril"
                  , longname        = "Faz. Cauaxi (Logged), PA"
                  , iata            = "ca6"
                  , lon             = -48.483
                  , lat             =  -3.748
                  , alt             = 140
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.00
                  , clay            = -1.00
                  , depth           = "H"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = +2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "caxiuana"          
                  , longname        = "Caxiuana, PA"                
                  , iata            = "cax"
                  , lon             = -51.458
                  , lat             =  -1.720
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 2
                  , ntext           = 16
                  , sand            = 0.380
                  , clay            = 0.440
                  , depth           = "D"
                  , isoilbc         = 2
                  , sldrain         = 6.
                  , scolour         = 14
                  , met.driver      = "Caxiuana"
                  , yeara           = 1994
                  , yearz           = 2004
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "cayenne"          
                  , longname        = "Cayenne, GUF"                
                  , iata            = "cay"
                  , lon             = -52.367
                  , lat             =   4.833
                  , alt             = 105.
                  , wmo             = 81405
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "chaiten"           
                  , longname        = "Chaiten, CHL"                
                  , iata            = "wch"
                  , lon             = -72.500
                  , lat             = -42.500
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "Sheffield"
                  , yeara           = 1969
                  , yearz           = 2008
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "ciudad_bolivar"           
                  , longname        = "Ciudad Bolivar, VEN"                
                  , iata            = "cbl"
                  , lon             = -63.550
                  , lat             =   8.150
                  , alt             = 48
                  , wmo             = 80444
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "ciudad_guayana"           
                  , longname        = "Ciudad Guayana, VEN"                
                  , iata            = "cgu"
                  , lon             = -62.762
                  , lat             =   8.289
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "cochabamba"        
                  , longname        = "Cochabamba, BOL"             
                  , iata            = "cbb"
                  , lon             = -66.170
                  , lat             = -17.420
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "cuiaba"            
                  , longname        = "Cuiaba, MT"                  
                  , iata            = "cgb"
                  , lon             = -56.070
                  , lat             = -15.555
                  , alt             = 185
                  , wmo             = 83362
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "Cuiaba"
                  , yeara           = 2002
                  , yearz           = 2013
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "curacao"            
                  , longname        = "Curacao, ANT"                  
                  , iata            = "cur"
                  , lon             = -68.967
                  , lat             =  12.200
                  , alt             = 9
                  , wmo             = 78988
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "curitiba"          
                  , longname        = "Curitiba, PR"                
                  , iata            = "cwb"
                  , lon             = -49.230
                  , lat             = -25.410
                  , alt             = 908
                  , wmo             = 83840
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "curua_una"
                  , longname        = "Curua-Una (Sudam), PA"
                  , iata            = "scr"
                  , lon             = -54.091
                  , lat             =  -2.544
                  , alt             = 19
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "cruzeiro_do_sul"
                  , longname        = "Cruzeiro do Sul, AC"
                  , iata            = "czs"
                  , lon             = -72.667
                  , lat             =  -7.663
                  , alt             = 170
                  , wmo             = 82705
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "diamantino"        
                  , longname        = "Diamantino, MT"              
                  , iata            = "dmt"
                  , lon             = -56.620
                  , lat             = -14.370
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "douglas_tnf"
                  , longname        = "Tapajos Natl Forest, PA"
                  , iata            = "dcm"
                  , lon             = -54.959
                  , lat             =  -2.857
                  , alt             = 198.
                  , wmo             = NA
                  , isoilflg        = 2
                  , ntext           = 16
                  , sand            = 0.390
                  , clay            = 0.590
                  , depth           = "H"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "Santarem_Km67"
                  , yeara           = 1996
                  , yearz           = 2012
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "dourados"          
                  , longname        = "Dourados, MS"                
                  , iata            = "dou"
                  , lon             = -54.810
                  , lat             = -22.220
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "ducke"            
                  , longname        = "Ducke , AM"                  
                  , iata            = "duc"
                  , lon             = -59.940
                  , lat             =  -2.951
                  , alt             = 117
                  , wmo             = NA
                  , isoilflg        = 2
                  , ntext           = 8
                  , sand            = 0.433
                  , clay            = 0.366
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "east_belterra"            
                  , longname        = "Belterra (East), PA"                  
                  , iata            = "ebt"
                  , lon             = -54.888
                  , lat             =  -3.156
                  , alt             = 147
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "east_feliz_natal"            
                  , longname        = "Feliz Natal (East), MT"                  
                  , iata            = "fn2"
                  , lon             = -54.197
                  , lat             = -11.862
                  , alt             = 336
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "east_sao_felix"   
                  , longname        = "Sao Felix (East), PA"      
                  , iata            = "sx2"
                  , lon             = -51.793
                  , lat             =  -6.601
                  , alt             = 224
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "el_triunfo"        
                  , longname        = "El Triunfo, BOL"             
                  , iata            = "etf"
                  , lon             = -67.000
                  , lat             = -13.500
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "el_zafire"         
                  , longname        = "El Zafire, COL"              
                  , iata            = "zar"
                  , lon             = -69.902
                  , lat             =  -4.007
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 2
                  , ntext           = 2
                  , sand            = 0.748
                  , clay            = 0.006
                  , depth           = "B"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "emas"         
                  , longname        = "Parque das Emas, GO"              
                  , iata            = "jti"
                  , lon             = -52.19
                  , lat             = -18.12
                  , alt             = 848
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "erechim"           
                  , longname        = "Erechim, RS"                 
                  , iata            = "erm"
                  , lon             = -52.240
                  , lat             = -27.610
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "eunapolis"         
                  , longname        = "Eunapolis, BA"               
                  , iata            = "enp"
                  , lon             = -39.580
                  , lat             = -16.330
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "fazendans"         
                  , longname        = "Fazenda Nsa. Sra., RO"       
                  , iata            = "fns"
                  , lon             = -62.357
                  , lat             = -10.762
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 2
                  , ntext           = 2
                  , sand            = 0.800
                  , clay            = 0.100
                  , depth           = "G"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "Fazenda_Nossa_Senhora"
                  , yeara           = 1994
                  , yearz           = 2003
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "feliz_natal"
                  , longname        = "Feliz Natal, MT"
                  , iata            = "fnz"
                  , lon             = -54.692
                  , lat             = -12.146
                  , alt             = 350.
                  , wmo             = NA
                  , isoilflg        = 2
                  , ntext           = 6
                  , sand            = 0.641
                  , clay            = 0.247
                  , depth           = "I"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "floriano"         
                  , longname        = "Floriano, PI"               
                  , iata            = "flb"
                  , lon             = -43.020
                  , lat             =  -6.770
                  , alt             = 123
                  , wmo             = 82678
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "fortaleza"         
                  , longname        = "Fortaleza, CE"               
                  , iata            = "for"
                  , lon             = -38.530
                  , lat             =  -3.780
                  , alt             = 26
                  , wmo             = 82397
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "guarana"
                  , longname        = "Guarana, PA"
                  , iata            = "gun"
                  , lon             = -54.325
                  , lat             =  -2.677
                  , alt             = 138
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "harvard"           
                  , longname        = "Harvard Forest, MA"          
                  , iata            = "hvd"
                  , lon             = -72.170
                  , lat             =  42.540
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 2
                  , ntext           = 2
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "E"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "Harvard"
                  , yeara           = 1987
                  , yearz           = 2003
                  , iphen           = 1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "humaita"
                  , longname        = "Fazenda Humaita, AC"
                  , iata            = "hum"
                  , lon             = -67.645
                  , lat             =  -9.756
                  , alt             = 172.
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "iguape"            
                  , longname        = "Iguape, SP"                  
                  , iata            = "igp"
                  , lon             = -47.590
                  , lat             = -24.630
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "imperatriz"        
                  , longname        = "Imperatriz, MA"              
                  , iata            = "imp"
                  , lon             = -47.460
                  , lat             =  -5.530
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "iquique"           
                  , longname        = "Iquique, CHL"                
                  , iata            = "iqq"
                  , lon             = -69.970
                  , lat             = -20.240
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "iquitos"           
                  , longname        = "Iquitos, PER"                
                  , iata            = "iqt"
                  , lon             = -73.250
                  , lat             =  -3.750
                  , alt             = 125
                  , wmo             = 84378
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "itabaiana"         
                  , longname        = "Itabaiana, SE"               
                  , iata            = "ibn"
                  , lon             = -37.420
                  , lat             = -10.680
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "itapeva"           
                  , longname        = "Itapeva, SP"                 
                  , iata            = "ipv"
                  , lon             = -48.880
                  , lat             = -23.980
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "itirapina"           
                  , longname        = "Itirapina, SP"                 
                  , iata            = "ity"
                  , lon             = -47.85
                  , lat             = -22.22
                  , alt             = 810
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "jacareacanga"      
                  , longname        = "Jacareacanga, PA"            
                  , iata            = "jcr"
                  , lon             = -57.777
                  , lat             =  -6.233
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "jamaraqua"
                  , longname        = "Jamaraqua, PA"
                  , iata            = "jmq"
                  , lon             = -55.036
                  , lat             =  -2.806
                  , alt             = 15
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "jamari"
                  , longname        = "Jamari Natl. Forest, RO"
                  , iata            = "jam"
                  , lon             = -63.007
                  , lat             =  -9.116
                  , alt             = 112.
                  , wmo             = NA
                  , isoilflg        = 2
                  , ntext           = 6
                  , sand            = 0.458
                  , clay            = 0.328
                  , depth           = "I"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "jiparana"          
                  , longname        = "Ji-Parana, RO"               
                  , iata            = "jpr"
                  , lon             = -61.980
                  , lat             = -10.860
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "joao_pessoa"       
                  , longname        = "Joao Pessoa, PB"             
                  , iata            = "jpa"
                  , lon             = -34.910
                  , lat             =  -7.100
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "joinville"       
                  , longname        = "Joinville, SC"             
                  , iata            = "joi"
                  , lon             = -48.858
                  , lat             = -26.252
                  , alt             = 48
                  , wmo             = 83905
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "Joinville"
                  , yeara           = 2004
                  , yearz           = 2013
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "kenia"             
                  , longname        = "Kenia, BOL"                  
                  , iata            = "qea"
                  , lon             = -62.730
                  , lat             = -16.010
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 2
                  , ntext           = 3
                  , sand            = 0.760
                  , clay            = 0.160
                  , depth           = "D"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "las_gaviotas"      
                  , longname        = "Las Gaviotas, COL"           
                  , iata            = "lgt"
                  , lon             = -70.970
                  , lat             =   4.550
                  , alt             = 167
                  , wmo             = 80241
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "la_esmeralda"      
                  , longname        = "La Esmeralda, VEN"           
                  , iata            = "lfe"
                  , lon             = -65.540
                  , lat             =   3.170
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "la_lorena"         
                  , longname        = "La Lorena, COL"              
                  , iata            = "lor"
                  , lon             = -69.991
                  , lat             =  -3.056
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 2
                  , ntext           = 9
                  , sand            = 0.380
                  , clay            = 0.310
                  , depth           = "A"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "la_planada"        
                  , longname        = "La Planada, COL"             
                  , iata            = "lpn"
                  , lon             = -77.994
                  , lat             =   1.116
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "H"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "la_selva"          
                  , longname        = "La Selva, CRI"               
                  , iata            = "lse"
                  , lon             = -84.010
                  , lat             =  10.430
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 2
                  , ntext           = 8
                  , sand            = 0.388
                  , clay            = 0.329
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "labrea"            
                  , longname        = "Labrea, AM"                  
                  , iata            = "lbr"
                  , lon             = -64.770
                  , lat             =  -7.280
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "lencois"           
                  , longname        = "Lencois, BA"                 
                  , iata            = "lec"
                  , lon             = -41.350
                  , lat             = -12.480
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "leticia"           
                  , longname        = "Leticia, COL"                 
                  , iata            = "let"
                  , lon             = -69.950
                  , lat             = -4.167
                  , alt             = 84
                  , wmo             = 80398
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "lima"           
                  , longname        = "Lima, PER"                 
                  , iata            = "lim"
                  , lon             = -77.117
                  , lat             = -12.000
                  , alt             = 12
                  , wmo             = 84629
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "linden"           
                  , longname        = "Linden, GUY"                 
                  , iata            = "lyd"
                  , lon             = -58.302
                  , lat             =   6.015
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "llochegua"           
                  , longname        = "Llochegua, PER"                 
                  , iata            = "llo"
                  , lon             = -73.908
                  , lat             = -12.410
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "macapa"            
                  , longname        = "Macapa, AP"                  
                  , iata            = "mcp"
                  , lon             = -51.090
                  , lat             =   0.330
                  , alt             = 16
                  , wmo             = 82099
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "malalcahuello"     
                  , longname        = "Malalcahuello, CHL"          
                  , iata            = "zmh"
                  , lon             = -71.580
                  , lat             = -38.470
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "Sheffield"
                  , yeara           = 1969
                  , yearz           = 2008
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "manaus"            
                  , longname        = "Manaus, AM"                  
                  , iata            = "mao"
                  , lon             = -60.020
                  , lat             =  -3.110
                  , alt             = 84
                  , wmo             = 82332
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "manaus_km34"       
                  , longname        = "Manaus - Km34, AM"           
                  , iata            = "m34"
                  , lon             = -60.209
                  , lat             =  -2.609
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 2
                  , ntext           = 11
                  , sand            = 0.200
                  , clay            = 0.680
                  , depth           = "H"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "Manaus_Km34"
                  , yeara           = 1994
                  , yearz           = 2007
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "manaus_bdffp"      
                  , longname        = "Manaus - BDFFP, AM"          
                  , iata            = "bdf"
                  , lon             = -60.209
                  , lat             =  -2.609
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "manicore"          
                  , longname        = "Manicore, AM"                
                  , iata            = "mnx"
                  , lon             = -61.280
                  , lat             =  -5.820
                  , alt             = 49
                  , wmo             = 82532
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "maracay"           
                  , longname        = "Maracay, VEN"                 
                  , iata            = "myc"
                  , lon             = -67.650
                  , lat             =  10.250
                  , alt             = 437
                  , wmo             = 80413
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "maringa"           
                  , longname        = "Maringa, PR"                 
                  , iata            = "mgf"
                  , lon             = -52.010
                  , lat             = -23.470
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "mariscal_estigarribia"           
                  , longname        = "Mariscal Estigarribia, PRY"
                  , iata            = "esg"
                  , lon             = -60.624
                  , lat             = -22.043
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "maracarume"           
                  , longname        = "Maracarume, MA"                 
                  , iata            = "zme"
                  , lon             = -45.954
                  , lat             =  -2.041
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "montes_claros"     
                  , longname        = "Montes Claros, MG"           
                  , iata            = "moc"
                  , lon             = -43.820
                  , lat             = -16.710
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "mojui_dos_campos"
                  , longname        = "Mojui dos Campos, PA"
                  , iata            = "mjs"
                  , lon             = -54.579
                  , lat             =  -2.767
                  , alt             = 128
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "natal"           
                  , longname        = "Natal, RN"                
                  , iata            = "nat"
                  , lon             = -35.206
                  , lat             = -5.837
                  , alt             = 58
                  , wmo             = 82599
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "Natal"
                  , yeara           = 2004
                  , yearz           = 2013
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "neuquen"           
                  , longname        = "Neuquen, ARG"                
                  , iata            = "nqn"
                  , lon             = -68.133
                  , lat             = -38.950
                  , alt             = 270
                  , wmo             = 87715
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "Sheffield"
                  , yeara           = 1969
                  , yearz           = 2008
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "nouragues"
                  , longname        = "Nouragues, GUF"
                  , iata            = "nrf"
                  , lon             = -52.683
                  , lat             =   4.033
                  , alt             = 77
                  , wmo             = NA
                  , isoilflg        = 2
                  , ntext           = 6
                  , sand            = 0.562
                  , clay            = 0.345
                  , depth           = "H"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "Nouragues"
                  , yeara           = 2012
                  , yearz           = 2018
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "oeiras"            
                  , longname        = "Oeiras, PI"                  
                  , iata            = "oei"
                  , lon             = -42.160
                  , lat             =  -7.020
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "palmas"            
                  , longname        = "Palmas, TO"                  
                  , iata            = "pmw"
                  , lon             = -48.360
                  , lat             = -10.290
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "panama"            
                  , longname        = "Panama, PAN"                  
                  , iata            = "how"
                  , lon             = -79.600
                  , lat             =   8.917
                  , alt             = 16
                  , wmo             = 78807
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "paramaribo"        
                  , longname        = "Paramaribo, SUR"             
                  , iata            = "pbm"
                  , lon             = -55.150
                  , lat             =   5.830
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "paracou"           
                  , longname        = "Paracou, GUF"                
                  , iata            = "gyf"
                  , lon             = -52.912
                  , lat             =   5.282
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 2
                  , ntext           = 6
                  , sand            = 0.562
                  , clay            = 0.345
                  , depth           = "H"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "Paracou"
                  , yeara           = 1999
                  , yearz           = 2015
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "neonita"
                  , longname        = "Neonita, PA"
                  , iata            = "par"
                  , lon             = -47.517
                  , lat             =  -3.315
                  , alt             = 122
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.00
                  , clay            = -1.00
                  , depth           = "H"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = +2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "paragominas"
                  , longname        = "Paragominas, PA"
                  , iata            = "prg"
                  , lon             = -47.792
                  , lat             =  -3.020
                  , alt             = 81
                  , wmo             = NA
                  , isoilflg        = 2
                  , ntext           = 6
                  , sand            = 0.471
                  , clay            = 0.328
                  , depth           = "I"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "Paragominas"
                  , yeara           = 2004
                  , yearz           = 2015
                  , iphen           = +2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "pedegigante"       
                  , longname        = "Pe-de-Gigante, SP"           
                  , iata            = "pdg"
                  , lon             = -47.650
                  , lat             = -21.619
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 2
                  , ntext           = 2
                  , sand            = 0.850
                  , clay            = 0.030
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "Pe-de-Gigante"
                  , yeara           = 1996
                  , yearz           = 2004
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "petrolina"
                  , longname        = "Petrolina, PE"
                  , iata            = "pnz"
                  , lon             = -40.370
                  , lat             =  -9.165
                  , alt             = 383
                  , wmo             = 82983
                  , isoilflg        = 2
                  , ntext           = 2
                  , sand            = 0.821
                  , clay            = 0.052
                  , depth           = "C"
                  , isoilbc         = 2
                  , sldrain         = 15.
                  , scolour         = 14
                  , met.driver      = "Petrolina"
                  , yeara           = 1999
                  , yearz           = 2012
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "piura"         
                  , longname        = "Piura, PER"               
                  , iata            = "piu"
                  , lon             = -80.617
                  , lat             =  -5.207
                  , alt             = 49
                  , wmo             = 84416
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "port_of_spain"      
                  , longname        = "Port-of-Spain, TTO"            
                  , iata            = "pos"
                  , lon             = -61.350
                  , lat             =  10.583
                  , alt             = 15
                  , wmo             = 78970
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "porto_de_moz"      
                  , longname        = "Porto de Moz, PA"            
                  , iata            = "ptq"
                  , lon             = -52.236
                  , lat             =  -1.741
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "porto_velho"      
                  , longname        = "Porto Velho, RO"            
                  , iata            = "pvh"
                  , lon             = -63.917
                  , lat             =  -8.767
                  , alt             = 102
                  , wmo             = 82824
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "pucallpa"          
                  , longname        = "Pucallpa, PER"               
                  , iata            = "pcl"
                  , lon             = -74.570
                  , lat             =  -8.380
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "puerto_suarez"     
                  , longname        = "Puerto Suarez, BOL"          
                  , iata            = "psz"
                  , lon             = -58.090
                  , lat             = -18.580
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "quibdo"            
                  , longname        = "Quibdo, COL"                 
                  , iata            = "uib"
                  , lon             = -76.640
                  , lat             =   5.690
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "rebio_jaru"        
                  , longname        = "Rebio Jaru, RO"              
                  , iata            = "rja"
                  , lon             = -61.931
                  , lat             = -10.083
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 2
                  , ntext           = 2
                  , sand            = 0.800
                  , clay            = 0.100
                  , depth           = "C"
                  , isoilbc         = 2
                  , sldrain         = 6
                  , scolour         = 14
                  , met.driver      = "Rebio_Jaru"
                  , yeara           = 1994
                  , yearz           = 2003
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "recife"            
                  , longname        = "Recife, PE"                  
                  , iata            = "rec"
                  , lon             = -34.910
                  , lat             =  -8.070
                  , alt             = 7
                  , wmo             = 82900
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "redencao"          
                  , longname        = "Redencao, PA"                
                  , iata            = "rdc"
                  , lon             = -49.980
                  , lat             =  -8.030
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "ribeirao_preto"    
                  , longname        = "Ribeirao Preto, SP"          
                  , iata            = "rao"
                  , lon             = -47.780
                  , lat             = -21.140
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "rio_branco"        
                  , longname        = "Rio Branco, AC"              
                  , iata            = "rbr"
                  , lon             = -67.752
                  , lat             =  -9.824
                  , alt             = 199.
                  , wmo             = 82917
                  , isoilflg        = 2
                  , ntext           = 8
                  , sand            = 0.382
                  , clay            = 0.308
                  , depth           = "I"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "riohacha"        
                  , longname        = "Riohacha, COL"              
                  , iata            = "rch"
                  , lon             = -72.930
                  , lat             =  11.530
                  , alt             = 4
                  , wmo             = 80035
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "salta"             
                  , longname        = "Salta, ARG"                  
                  , iata            = "sla"
                  , lon             = -65.483
                  , lat             = -24.850
                  , alt             = 1238
                  , wmo             = 87047
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "salvador"             
                  , longname        = "Salvador, BA"                  
                  , iata            = "ssa"
                  , lon             = -38.517
                  , lat             = -13.017
                  , alt             = 91
                  , wmo             = 83229
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "san_andres"             
                  , longname        = "San Andres, COL"                  
                  , iata            = "adz"
                  , lon             = -81.717
                  , lat             =  12.583
                  , alt             = 1
                  , wmo             = 80001
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "san_antonio_del_tachira"             
                  , longname        = "San Antonio del Tachira, VEN"                  
                  , iata            = "svz"
                  , lon             = -72.450
                  , lat             =  7.850
                  , alt             = 378
                  , wmo             = 80447
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "san_jose"             
                  , longname        = "San Jose, CRI"                  
                  , iata            = "sjo"
                  , lon             = -84.217
                  , lat             =  9.983
                  , alt             = 939
                  , wmo             = 78762
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "san_fernando_de_apure"             
                  , longname        = "San Fernando de Apure, VEN"                  
                  , iata            = "sfd"
                  , lon             = -67.417
                  , lat             =  7.867
                  , alt             = 48
                  , wmo             = 80450
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "san_pedro"         
                  , longname        = "San Pedro, ARG"              
                  , iata            = "zpe"
                  , lon             = -54.110
                  , lat             = -26.630
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "sao_carlos"             
                  , longname        = "Sao Carlos, SP"                  
                  , iata            = "qsc"
                  , lon             = -47.83
                  , lat             = -21.97
                  , alt             = 904
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "sao_felix_araguaia"
                  , longname        = "Sao Felix do Araguaia, MT"   
                  , iata            = "sxo"
                  , lon             = -50.690
                  , lat             = -11.630
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "sao_felix_xingu"
                  , longname        = "Sao Felix do Xingu, PA"
                  , iata            = "sfx"
                  , lon             = -52,349
                  , lat             =  -6.507
                  , alt             = 197.
                  , wmo             = NA
                  , isoilflg        = 2
                  , ntext           = 6
                  , sand            = 0.451
                  , clay            = 0.309
                  , depth           = "I"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "sao_gabriel"       
                  , longname        = "Sao Gabriel da Cachoeira, AM"
                  , iata            = "sjl"
                  , lon             = -66.980
                  , lat             =  -0.140
                  , alt             = 90
                  , wmo             = 82107
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "sao_jorge"       
                  , longname        = "Sao Jorge, PA"
                  , iata            = "tsj"
                  , lon             = -54.971
                  , lat             =  -3.135
                  , alt             = 111
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "sao_luis"          
                  , longname        = "Sao Luis, MA"                
                  , iata            = "slz"
                  , lon             = -44.236
                  , lat             =  -2.586
                  , alt             = 53
                  , wmo             = 82281
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "sao_martinho"          
                  , longname        = "Sao Martinho da Serra, RS"                
                  , iata            = "mhs"
                  , lon             = -53.823
                  , lat             = -29.443
                  , alt             = 489
                  , wmo             = 83937
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "Sao_Martinho"
                  , yeara           = 2000
                  , yearz           = 2013
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "santa_fe"          
                  , longname        = "Santa Fe, ARG"                
                  , iata            = "sfn"
                  , lon             = -60.809
                  , lat             = -31.712
                  , alt             = 198.
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "santarem"          
                  , longname        = "Santarem, PA"                
                  , iata            = "stm"
                  , lon             = -54.740
                  , lat             =  -2.419
                  , alt             = 198.
                  , wmo             = 82244
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "santarem_km67"     
                  , longname        = "Santarem - Km67, PA"         
                  , iata            = "s67"
                  , lon             = -54.959
                  , lat             =  -2.857
                  , alt             = 198.
                  , wmo             = NA
                  , isoilflg        = 2
                  , ntext           = 16
                  , sand            = 0.390
                  , clay            = 0.590
                  , depth           = "I"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "Santarem_Km67"
                  , yeara           = 1996
                  , yearz           = 2012
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "santarem_km77"     
                  , longname        = "Santarem - Km 77, PA"        
                  , iata            = "s77"
                  , lon             = -54.537
                  , lat             =  -3.012
                  , alt             = 33.
                  , wmo             = NA
                  , isoilflg        = 2
                  , ntext           = 16
                  , sand            = 0.390
                  , clay            = 0.590
                  , depth           = "I"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "Santarem_Km77"
                  , yeara           = 1996
                  , yearz           = 2006
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "santarem_km83"     
                  , longname        = "Santarem - Km 83, PA"        
                  , iata            = "s83"
                  , lon             = -54.971
                  , lat             =  -3.018
                  , alt             = 191
                  , wmo             = NA
                  , isoilflg        = 2
                  , ntext           = 16
                  , sand            = 0.390
                  , clay            = 0.590
                  , depth           = "I"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "Santarem_Km83"
                  , yeara           = 1995
                  , yearz           = 2004
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "santarem_km117"
                  , longname        = "Santarem - Km 117, PA"
                  , iata            = "117"
                  , lon             = -54.936
                  , lat             =  -3.350
                  , alt             = 186
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "saraca-taquera"     
                  , longname        = "Saraca-Taquera Natl. Forest, PA"
                  , iata            = "fst"
                  , lon             = -56.223
                  , lat             = -1.622
                  , alt             = 104.
                  , wmo             = NA
                  , isoilflg        = 2
                  , ntext           = 6
                  , sand            = 0.588
                  , clay            = 0.272
                  , depth           = "I"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "serra_do_navio"             
                  , longname        = "Serra do Navio, AP"               
                  , iata            = "svi"
                  , lon             = -52.000
                  , lat             =   0.906
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "sinop"             
                  , longname        = "Sinop, MT"                   
                  , iata            = "ops"
                  , lon             = -55.322
                  , lat             = -11.421
                  , alt             = 423
                  , wmo             = NA
                  , isoilflg        = 2
                  , ntext           = 2
                  , sand            = 0.840
                  , clay            = 0.120
                  , depth           = "E"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "sobral"             
                  , longname        = "Sobral, CE"               
                  , iata            = "qbx"
                  , lon             = -40.337
                  , lat             =  -3.679
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "sorriso"             
                  , longname        = "Sorriso, MT"                   
                  , iata            = "smt"
                  , lon             = -55.723
                  , lat             = -12.555
                  , alt             = 379
                  , wmo             = NA
                  , isoilflg        = 2
                  , ntext           = 2
                  , sand            = 0.840
                  , clay            = 0.120
                  , depth           = "E"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "Sorriso"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "tabatinga"         
                  , longname        = "Tabatinga, AM"              
                  , iata            = "tbt"
                  , lon             = -69.667
                  , lat             = -3.667
                  , alt             = 85
                  , wmo             = 82411
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "talisma"   
                  , longname        = "Fazenda Talisma, PA"      
                  , iata            = "tal"
                  , lon             = -67.977
                  , lat             = -10.255
                  , alt             = 189
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "tambopata"         
                  , longname        = "Tambopata, PER"              
                  , iata            = "tam"
                  , lon             = -69.271
                  , lat             = -12.830
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 2
                  , ntext           = 11
                  , sand            = 0.400
                  , clay            = 0.430
                  , depth           = "B"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "tanguro"     
                  , longname        = "Fazenda Tanguro, MT"
                  , iata            = "tan"
                  , lon             = -52.409
                  , lat             = -13.081
                  , alt             = 349.
                  , wmo             = NA
                  , isoilflg        = 2
                  , ntext           = 6
                  , sand            = 0.664
                  , clay            = 0.246
                  , depth           = "I"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "Tanguro"
                  , yeara           = 2008
                  , yearz           = 2018
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "tanguro"     
                  , longname        = "Fazenda Tanguro, MT"
                  , iata            = "tb0"
                  , lon             = -52.409
                  , lat             = -13.081
                  , alt             = 349.
                  , wmo             = NA
                  , isoilflg        = 2
                  , ntext           = 6
                  , sand            = 0.664
                  , clay            = 0.246
                  , depth           = "I"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "Tanguro"
                  , yeara           = 2008
                  , yearz           = 2018
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "tapajos"     
                  , longname        = "Tapajos Natl. Forest, PA"
                  , iata            = "tni"
                  , lon             = -54.955
                  , lat             =  -2.857
                  , alt             = 199.
                  , wmo             = NA
                  , isoilflg        = 2
                  , ntext           = 16
                  , sand            = 0.390
                  , clay            = 0.590
                  , depth           = "H"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "Santarem_Km67"
                  , yeara           = 1996
                  , yearz           = 2012
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "tapajos"     
                  , longname        = "Tapajos Natl. Forest, PA"
                  , iata            = "tnf"
                  , lon             = -54.941
                  , lat             = -3.129
                  , alt             = 110.
                  , wmo             = NA
                  , isoilflg        = 2
                  , ntext           = 16
                  , sand            = 0.433
                  , clay            = 0.430
                  , depth           = "I"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "Santarem_Km67"
                  , yeara           = 1996
                  , yearz           = 2012
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "tarauaca"          
                  , longname        = "Tarauaca, AC"                
                  , iata            = "trq"
                  , lon             = -70.781
                  , lat             =  -8.157
                  , alt             = 197
                  , wmo             = 82807
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "tefe"              
                  , longname        = "Tefe, AM"                    
                  , iata            = "tff"
                  , lon             = -64.720
                  , lat             =  -3.380
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "teresina"          
                  , longname        = "Teresina, PI"                
                  , iata            = "the"
                  , lon             = -42.800
                  , lat             =  -5.090
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "tirios"            
                  , longname        = "Tirios, PA"                  
                  , iata            = "obi"
                  , lon             = -55.940
                  , lat             =   2.220
                  , alt             = 325
                  , wmo             = 82026
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "tolhuin"
                  , longname        = "Tolhuin, ARG"
                  , iata            = "tqh"
                  , lon             = -67.222
                  , lat             = -54.502
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "Sheffield"
                  , yeara           = 1969
                  , yearz           = 2008
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "tome-acu"   
                  , longname        = "Tome-Acu, PA"      
                  , iata            = "tac"
                  , lon             = -48.314
                  , lat             =  -2.461
                  , alt             = 50
                  , wmo             = NA
                  , isoilflg        = 2
                  , ntext           = 6
                  , sand            = 0.504
                  , clay            = 0.294
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "tonzi"
                  , longname        = "Tonzi, CA"
                  , iata            = "tzi"
                  , lon             = -120.894
                  , lat             =   38.427
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 2
                  , ntext           = 5
                  , sand            = 0.447
                  , clay            = 0.191
                  , depth           = "D"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "Tonzi"
                  , yeara           = 1964
                  , yearz           = 2009
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "vina_del_mar"      
                  , longname        = "Vina del Mar, CHL"           
                  , iata            = "kna"
                  , lon             = -71.480
                  , lat             = -32.950
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "vila_franca"
                  , longname        = "Vila Franca, PA"
                  , iata            = "vfr"
                  , lon             = -55.029
                  , lat             =  -2.349
                  , alt             = 13
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "vilhena"           
                  , longname        = "Vilhena, RO"             
                  , iata            = "bvh"
                  , lon             = -60.100
                  , lat             = -12.730
                  , alt             = 612
                  , wmo             = 83208
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "vitoria"           
                  , longname        = "Vitoria, ES"                 
                  , iata            = "vix"
                  , lon             = -40.390
                  , lat             = -20.310
                  , alt             = 3
                  , wmo             = 83649
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "west_feliz_natal"
                  , longname        = "Feliz Natal (Vitoria), MT"
                  , iata            = "fna"
                  , lon             = -55.019
                  , lat             = -12.501
                  , alt             = 392.
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "west_sao_felix"   
                  , longname        = "Sao Felix (West), PA"      
                  , iata            = "sx1"
                  , lon             = -52.904
                  , lat             =  -6.410
                  , alt             = 213
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "xingu"             
                  , longname        = "Xingu National Park, MT"     
                  , iata            = "xgu"
                  , lon             = -52.636
                  , lat             =  -9.692
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 1
                  , ntext           = 1
                  , sand            = -1.000
                  , clay            = -1.000
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = 2
                  )#end list
u           = u + 1
poitmp[[u]] = list( short           = "yasuni"            
                  , longname        = "Yasuni, ECU"                 
                  , iata            = "ysn"
                  , lon             = -76.396
                  , lat             =  -0.686
                  , alt             = NA
                  , wmo             = NA
                  , isoilflg        = 2
                  , ntext           = 4
                  , sand            = 0.259
                  , clay            = 0.255
                  , depth           = "F"
                  , isoilbc         = 1
                  , sldrain         = 90.
                  , scolour         = 14
                  , met.driver      = "MERRA2_MSWEP2"
                  , yeara           = 1974
                  , yearz           = 2017
                  , iphen           = -1
                  )#end list
#----- Make the table global. -------------------------------------------------------------#
npoi <<- length(poitmp)
poilist = list()
for (pp in sequence(npoi)){
   for (n in union(names(poilist),names(poitmp[[pp]]))){
      poilist[[n]] = c(poilist[[n]],poitmp[[pp]][[n]])
   }#end for
}#end for
poilist <<- data.frame(poilist,stringsAsFactors=FALSE)
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This list has the longitude and latitude of Ke's test polygons.                      #
#------------------------------------------------------------------------------------------#
u     = 0
kztmp = list()
u     = u + 1; kztmp[[u]] = list(lon= -60.5,lat=   5.5)
u     = u + 1; kztmp[[u]] = list(lon= -61.5,lat=   4.5)
u     = u + 1; kztmp[[u]] = list(lon= -59.5,lat=   4.5)
u     = u + 1; kztmp[[u]] = list(lon= -60.5,lat=   3.5)
u     = u + 1; kztmp[[u]] = list(lon= -58.5,lat=   3.5)
u     = u + 1; kztmp[[u]] = list(lon= -56.5,lat=   3.5)
u     = u + 1; kztmp[[u]] = list(lon= -59.5,lat=   2.5)
u     = u + 1; kztmp[[u]] = list(lon= -57.5,lat=   2.5)
u     = u + 1; kztmp[[u]] = list(lon= -55.5,lat=   2.5)
u     = u + 1; kztmp[[u]] = list(lon= -58.5,lat=   1.5)
u     = u + 1; kztmp[[u]] = list(lon= -56.5,lat=   1.5)
u     = u + 1; kztmp[[u]] = list(lon= -54.5,lat=   1.5)
u     = u + 1; kztmp[[u]] = list(lon= -57.5,lat=   0.5)
u     = u + 1; kztmp[[u]] = list(lon= -55.5,lat=   0.5)
u     = u + 1; kztmp[[u]] = list(lon= -53.5,lat=   0.5)
u     = u + 1; kztmp[[u]] = list(lon= -54.5,lat=  -0.5)
u     = u + 1; kztmp[[u]] = list(lon= -52.5,lat=  -0.5)
u     = u + 1; kztmp[[u]] = list(lon= -53.5,lat=  -1.5)
u     = u + 1; kztmp[[u]] = list(lon= -51.5,lat=  -1.5)
u     = u + 1; kztmp[[u]] = list(lon= -52.5,lat=  -2.5)
u     = u + 1; kztmp[[u]] = list(lon= -50.5,lat=  -2.5)
u     = u + 1; kztmp[[u]] = list(lon= -72.5,lat= -10.5)
u     = u + 1; kztmp[[u]] = list(lon= -71.5,lat= -10.5)
u     = u + 1; kztmp[[u]] = list(lon= -70.5,lat= -10.5)
u     = u + 1; kztmp[[u]] = list(lon= -69.5,lat= -10.5)
u     = u + 1; kztmp[[u]] = list(lon= -68.5,lat= -10.5)
u     = u + 1; kztmp[[u]] = list(lon= -67.5,lat= -10.5)
u     = u + 1; kztmp[[u]] = list(lon= -66.5,lat= -10.5)
u     = u + 1; kztmp[[u]] = list(lon= -65.5,lat= -10.5)
u     = u + 1; kztmp[[u]] = list(lon= -64.5,lat= -10.5)
u     = u + 1; kztmp[[u]] = list(lon= -63.5,lat= -10.5)
u     = u + 1; kztmp[[u]] = list(lon= -62.5,lat= -10.5)
u     = u + 1; kztmp[[u]] = list(lon= -61.5,lat= -10.5)
u     = u + 1; kztmp[[u]] = list(lon= -60.5,lat= -10.5)
u     = u + 1; kztmp[[u]] = list(lon= -59.5,lat= -10.5)
u     = u + 1; kztmp[[u]] = list(lon= -58.5,lat= -10.5)
u     = u + 1; kztmp[[u]] = list(lon= -57.5,lat= -10.5)
u     = u + 1; kztmp[[u]] = list(lon= -56.5,lat= -10.5)
u     = u + 1; kztmp[[u]] = list(lon= -55.5,lat= -10.5)
u     = u + 1; kztmp[[u]] = list(lon= -54.5,lat= -10.5)
u     = u + 1; kztmp[[u]] = list(lon= -53.5,lat= -10.5)
u     = u + 1; kztmp[[u]] = list(lon= -52.5,lat= -10.5)
u     = u + 1; kztmp[[u]] = list(lon= -51.5,lat= -10.5)
u     = u + 1; kztmp[[u]] = list(lon= -50.5,lat= -10.5)
u     = u + 1; kztmp[[u]] = list(lon= -49.5,lat= -10.5)
u     = u + 1; kztmp[[u]] = list(lon= -72.5,lat= -11.5)
u     = u + 1; kztmp[[u]] = list(lon= -71.5,lat= -11.5)
u     = u + 1; kztmp[[u]] = list(lon= -70.5,lat= -11.5)
u     = u + 1; kztmp[[u]] = list(lon= -69.5,lat= -11.5)
u     = u + 1; kztmp[[u]] = list(lon= -68.5,lat= -11.5)
u     = u + 1; kztmp[[u]] = list(lon= -67.5,lat= -11.5)
u     = u + 1; kztmp[[u]] = list(lon= -66.5,lat= -11.5)
u     = u + 1; kztmp[[u]] = list(lon= -65.5,lat= -11.5)
u     = u + 1; kztmp[[u]] = list(lon= -64.5,lat= -11.5)
u     = u + 1; kztmp[[u]] = list(lon= -63.5,lat= -11.5)
u     = u + 1; kztmp[[u]] = list(lon= -62.5,lat= -11.5)
u     = u + 1; kztmp[[u]] = list(lon= -61.5,lat= -11.5)
u     = u + 1; kztmp[[u]] = list(lon= -60.5,lat= -11.5)
u     = u + 1; kztmp[[u]] = list(lon= -59.5,lat= -11.5)
u     = u + 1; kztmp[[u]] = list(lon= -58.5,lat= -11.5)
u     = u + 1; kztmp[[u]] = list(lon= -57.5,lat= -11.5)
u     = u + 1; kztmp[[u]] = list(lon= -56.5,lat= -11.5)
u     = u + 1; kztmp[[u]] = list(lon= -55.5,lat= -11.5)
u     = u + 1; kztmp[[u]] = list(lon= -54.5,lat= -11.5)
u     = u + 1; kztmp[[u]] = list(lon= -53.5,lat= -11.5)
u     = u + 1; kztmp[[u]] = list(lon= -52.5,lat= -11.5)
u     = u + 1; kztmp[[u]] = list(lon= -51.5,lat= -11.5)
u     = u + 1; kztmp[[u]] = list(lon= -50.5,lat= -11.5)
u     = u + 1; kztmp[[u]] = list(lon= -49.5,lat= -11.5)
u     = u + 1; kztmp[[u]] = list(lon= -72.5,lat= -12.5)
u     = u + 1; kztmp[[u]] = list(lon= -71.5,lat= -12.5)
u     = u + 1; kztmp[[u]] = list(lon= -70.5,lat= -12.5)
u     = u + 1; kztmp[[u]] = list(lon= -69.5,lat= -12.5)
u     = u + 1; kztmp[[u]] = list(lon= -68.5,lat= -12.5)
u     = u + 1; kztmp[[u]] = list(lon= -67.5,lat= -12.5)
u     = u + 1; kztmp[[u]] = list(lon= -66.5,lat= -12.5)
u     = u + 1; kztmp[[u]] = list(lon= -65.5,lat= -12.5)
u     = u + 1; kztmp[[u]] = list(lon= -64.5,lat= -12.5)
u     = u + 1; kztmp[[u]] = list(lon= -63.5,lat= -12.5)
u     = u + 1; kztmp[[u]] = list(lon= -62.5,lat= -12.5)
u     = u + 1; kztmp[[u]] = list(lon= -61.5,lat= -12.5)
u     = u + 1; kztmp[[u]] = list(lon= -60.5,lat= -12.5)
u     = u + 1; kztmp[[u]] = list(lon= -59.5,lat= -12.5)
u     = u + 1; kztmp[[u]] = list(lon= -58.5,lat= -12.5)
u     = u + 1; kztmp[[u]] = list(lon= -57.5,lat= -12.5)
u     = u + 1; kztmp[[u]] = list(lon= -56.5,lat= -12.5)
u     = u + 1; kztmp[[u]] = list(lon= -55.5,lat= -12.5)
u     = u + 1; kztmp[[u]] = list(lon= -54.5,lat= -12.5)
u     = u + 1; kztmp[[u]] = list(lon= -53.5,lat= -12.5)
u     = u + 1; kztmp[[u]] = list(lon= -52.5,lat= -12.5)
u     = u + 1; kztmp[[u]] = list(lon= -51.5,lat= -12.5)
u     = u + 1; kztmp[[u]] = list(lon= -50.5,lat= -12.5)
u     = u + 1; kztmp[[u]] = list(lon= -49.5,lat= -12.5)
u     = u + 1; kztmp[[u]] = list(lon= -72.5,lat= -13.5)
u     = u + 1; kztmp[[u]] = list(lon= -71.5,lat= -13.5)
u     = u + 1; kztmp[[u]] = list(lon= -70.5,lat= -13.5)
u     = u + 1; kztmp[[u]] = list(lon= -69.5,lat= -13.5)
u     = u + 1; kztmp[[u]] = list(lon= -68.5,lat= -13.5)
u     = u + 1; kztmp[[u]] = list(lon= -67.5,lat= -13.5)
u     = u + 1; kztmp[[u]] = list(lon= -66.5,lat= -13.5)
u     = u + 1; kztmp[[u]] = list(lon= -65.5,lat= -13.5)
u     = u + 1; kztmp[[u]] = list(lon= -64.5,lat= -13.5)
u     = u + 1; kztmp[[u]] = list(lon= -63.5,lat= -13.5)
u     = u + 1; kztmp[[u]] = list(lon= -62.5,lat= -13.5)
u     = u + 1; kztmp[[u]] = list(lon= -61.5,lat= -13.5)
u     = u + 1; kztmp[[u]] = list(lon= -60.5,lat= -13.5)
u     = u + 1; kztmp[[u]] = list(lon= -59.5,lat= -13.5)
u     = u + 1; kztmp[[u]] = list(lon= -58.5,lat= -13.5)
u     = u + 1; kztmp[[u]] = list(lon= -57.5,lat= -13.5)
u     = u + 1; kztmp[[u]] = list(lon= -56.5,lat= -13.5)
u     = u + 1; kztmp[[u]] = list(lon= -55.5,lat= -13.5)
u     = u + 1; kztmp[[u]] = list(lon= -54.5,lat= -13.5)
u     = u + 1; kztmp[[u]] = list(lon= -53.5,lat= -13.5)
u     = u + 1; kztmp[[u]] = list(lon= -52.5,lat= -13.5)
u     = u + 1; kztmp[[u]] = list(lon= -51.5,lat= -13.5)
u     = u + 1; kztmp[[u]] = list(lon= -50.5,lat= -13.5)
u     = u + 1; kztmp[[u]] = list(lon= -49.5,lat= -13.5)
u     = u + 1; kztmp[[u]] = list(lon= -62.5,lat= -16.5)
u     = u + 1; kztmp[[u]] = list(lon= -60.5,lat= -2.5 )
u     = u + 1; kztmp[[u]] = list(lon= -54.5,lat= -2.5 )
u     = u + 1; kztmp[[u]] = list(lon= -54.5,lat= -3.5 )
u     = u + 1; kztmp[[u]] = list(lon= -50.5,lat= -9.5 )
u     = u + 1; kztmp[[u]] = list(lon= -69.5,lat= -3.5 )
u     = u + 1; kztmp[[u]] = list(lon= -69.5,lat= -4.5 )
u     = u + 1; kztmp[[u]] = list(lon= -52.5,lat= -5.5 )
u     = u + 1; kztmp[[u]] = list(lon= -73.5,lat= -5.5 )
#----- Make the table global. -------------------------------------------------------------#
nkzpoly <<- length(kztmp)
kzlist = list()
for (kz in sequence(nkzpoly)){
   for (n in union(names(kzlist),names(kztmp[[kz]]))){
      kzlist[[n]] = c(kzlist[[n]],kztmp[[kz]][[n]])
   }#end for
}#end for
kzlist <<- kzlist  
#------------------------------------------------------------------------------------------#
