#==========================================================================================#
#==========================================================================================#
#     Reset session.                                                                       #
#------------------------------------------------------------------------------------------#
rm(list=ls())
graphics.off()
#------------------------------------------------------------------------------------------#





#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#      Here is the user defined variable section.                                          #
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
here       = getwd()                                  # Current directory
srcdir     = c("/prj/prjidfca/marcosl/Util/Rsc"       # Possible paths with libraries
              ,"/prj/bramsolam/marcos.longo/Util/Rsc" #    R will select the first one that
              ,"/n/home00/mlongo/util/Rsc"            #    is found, or issue an error
              ,"/Users/mlongo/Util/Rsc"               #    message in case none of them
              ,"/home/b14275/Util/Rsc"                #    exist.
              )#end c                                 #
outfile    = file.path(here,"joborder.txt")           # Job order
defjob     = FALSE                                    # Generate the default job order?
append.job = FALSE                                    # Append job? (FALSE means new file)
lonlat     = NULL                                     # NULL - define runs locally 
                                                      #        (varrun/varlabel)
#lonlat  = file.path(here,"lonlat_input.txt")         # Not NULL - read lon/lat from file,
                                                      #    and finish up settings below
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Variables that will vary (check the list below for all possibilities).              #
#------------------------------------------------------------------------------------------#
if (! defjob && is.null(lonlat)){
   #----- Site-run simulation, define the sites and settings. -----------------------------#
   varrun   = list( iata      = c("gyf","m34","s67","s83","ban","rja","pdg","cax","pnz")
                  , iscenario = c("wmo")
                  , ibigleaf  = c(0,1)
                  , iage      = c(1,25)
                  , isizepft  = c(2 ,5)
                  , ivegtdyn  = 0
                  , queue     = "general"
                  )#end list
   varlabel = list( iata      = varrun$iata
                  , iscenario = varrun$iscenario
                  , ibigleaf  = c("sas","ble")
                  , iage      = paste("iage" , sprintf( "%2.2i",varrun$iage     ), sep="" )
                  , isizepft  = paste("pft"  , sprintf( "%2.2i",varrun$isizepft ), sep="" )
                  , ivegtdyn  = 0
                  , queue     = varrun$queue
                  )#end list
}else if(! defjob){
   #----- Regional run.  Read in the list of polygons, and complete the settings. ---------#
   myruns  = read.table( file             = lonlat
                       , header           = TRUE
                       , stringsAsFactors = FALSE
                       , comment.char     = ""
                       , colClasses       = c("character",rep("numeric",times=2))
                       )#end read.table
   nruns             = nrow(myruns)
   myruns$iscenario  = rep("default"  ,times=nruns)
   myruns$met.driver = rep("Sheffield",times=nruns)
   myruns$init.mode  = rep(          0,times=nruns)
   myruns$iphen      = rep(          2,times=nruns)
   myruns$fire       = rep(          3,times=nruns)
   myruns$iage       = rep(         25,times=nruns)
   myruns$ibigleaf   = rep(          0,times=nruns)
   myruns$yeara      = rep(       1500,times=nruns)
   myruns$yearz      = rep(       2011,times=nruns)

   nvars             = ncol(myruns)
}#end if
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




#----- Load some packages and scripts. ----------------------------------------------------#
if (any(file.exists(srcdir))){
   srcdir = (srcdir[file.exists(srcdir)])[1]
   source(file.path(srcdir,"load.everything.r"))
}else{
   stop("None of the paths provided in variable \"srcdir\" exist.")
}#end if (any(file.exists(srcdir)))
#------------------------------------------------------------------------------------------#


#----- Check whether there is a joborder.  In case there is one, back it up. --------------#
backupfile = file.path(dirname(outfile),paste("old",basename(outfile),sep="_"))
if (file.exists(outfile)) dummy = file.copy(from=outfile,to=backupfile)
#------------------------------------------------------------------------------------------#



#----- Check that varrun and varlabel have the same length and dimensions. ----------------#
if (! defjob && is.null(lonlat)){
   if (length(varrun) != length(varlabel)){
      stop(" Varrun and varlabel must have the same number of variables!")
   }else if (any(names(varrun) != names(varlabel))){
      stop(" Variable names of varrun and varlabel must match an be in the same order!")
   }else if (any(sapply(X=varrun,FUN=length) != sapply(X=varlabel,FUN=length))){
      length.matrix = cbind( run   = sapply(X=varrun  ,FUN=length)
                           , label = sapply(X=varlabel,FUN=length)
                           )#end cbind
      print(length.matrix)
      stop(" Length of all variables in varrun and varlabel must match!")
   }#end if
}#end if
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Default properties.                                                                  #
#------------------------------------------------------------------------------------------#
default = list( run           = "unnamed"
              , iata          = "xxx"
              , lon           = 0.00
              , lat           = 0.00
              , yeara         = "1967"
              , montha        = "01"
              , daya          = "01"
              , timea         = "0000"
              , yearz         = "2013"
              , monthz        = "01"
              , dayz          = "01"
              , timez         = "0000"
              , init.mode     = 6
              , iscenario     = "default"
              , isizepft      = 0
              , iage          = 30
              , imaxcohort    = 50
              , isoilflg      = 1
              , istext        = 1
              , sand          = -1.0
              , clay          = -1.0
              , depth         = "F"
              , isoilbc       = 1
              , sldrain       = 90.
              , scolour       = 16
              , slzres        = 0
              , queue         = "linux.q"
              , met.driver    = "tower"
              , dtlsm         = 600.
              , vmfact.c3     = 1.00
              , vmfact.c4     = 1.00
              , mphoto.trc3   = 9.0
              , mphoto.tec3   = 7.2
              , mphoto.c4     = 5.2
              , bphoto.blc3   = 10000.
              , bphoto.nlc3   = 1000.
              , bphoto.c4     = 10000.
              , kw.grass      = 900.
              , kw.tree       = 600.
              , gamma.c3      = 0.020
              , gamma.c4      = 0.040
              , d0.grass      = 0.016
              , d0.tree       = 0.016
              , alpha.c3      = 0.080
              , alpha.c4      = 0.055
              , klowco2       = 4000.
              , decomp.scheme = 2
              , rrffact       = 1.000
              , growthresp    = 0.400
              , lwidth.grass  = 0.05
              , lwidth.bltree = 0.10
              , lwidth.nltree = 0.05
              , q10.c3        = 2.4
              , q10.c4        = 2.4
              , h2o.limit     = 2
              , imort.scheme  = 1
              , ddmort.const  = 0.8
              , cbr.scheme    = 0
              , isfclyrm      = 4
              , icanturb      = 2
              , ubmin         = 1.00
              , ugbmin        = 0.40
              , ustmin        = 0.10
              , gamm          = 13.0
              , gamh          = 13.0
              , tprandtl      = 0.74
              , ribmax        = 0.50
              , atmco2        = 400.
              , thcrit        = -1.20
              , sm.fire       = -1.45
              , ifire         = 0
              , fire.parm     = 0.5
              , ipercol       = 0
              , runoff.time   = 3600.
              , imetrad       = 2
              , ibranch       = 1
              , icanrad       = 2
              , ihrzrad       = 0
              , crown.mod     = 0
              , ltrans.vis    = 0.050
              , lreflect.vis  = 0.100
              , ltrans.nir    = 0.200
              , lreflect.nir  = 0.400
              , orient.tree   = 0.100
              , orient.grass  = 0.000
              , clump.tree    = 0.800
              , clump.grass   = 1.000
              , igoutput      = 0
              , ivegtdyn      = 1
              , igndvap       = 0
              , iphen         = -1
              , iallom        = 4
              , ibigleaf      = 0
              , irepro        = 2
              , treefall      = -0.0140
              , ianth.disturb = 0
              , ianth.dataset = "glu-331"
              ) #end list
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Create the full combination of runs (if reading lonlat, this has been already        #
# generated).                                                                              #
#------------------------------------------------------------------------------------------#
if (defjob){
   myruns  = data.frame(iata=poilist$iata,stringsAsFactors=FALSE)
   defname = poilist$short
   nvars     = ncol(myruns)
   nruns     = nrow(myruns)
}else if (is.null(lonlat)){
   myruns  = expand.grid(varrun,stringsAsFactors=FALSE)
   nvars     = ncol(myruns)
   nruns     = nrow(myruns)
}#end if
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Create a table with room for all simulations.  Check whether this is a special case. #
#------------------------------------------------------------------------------------------#
joborder = data.frame(sapply(X=default,FUN=rep,times=nruns),stringsAsFactors=FALSE)
for (n in sequence(nvars)){
   name.now = names(myruns)[n]
   if (name.now %in% "labsorb.vis"){
      joborder$ltrans.vis   = 1. * (1. - myruns$labsorb.vis) / 3.
      joborder$lreflect.vis = 2. * (1. - myruns$labsorb.vis) / 3.
   }else if (name.now %in% "labsorb.nir"){
      joborder$ltrans.nir   = 1. * (1. - myruns$labsorb.nir) / 3.
      joborder$lreflect.nir = 2. * (1. - myruns$labsorb.nir) / 3.
   }else if (name.now %in% names(joborder)){
      joborder[[name.now]] = myruns[[name.now]]
   }else{
      stop(paste(" Variable ",name.now," is not a valid joborder variable!",sep=""))
   }#end if
   #---------------------------------------------------------------------------------------#
}#end for (n in sequence(nvars))
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Remove forbidden combinations.   Big leaf cannot be run with a single patch, because #
# each PFT is allocated to one patch.  Also, if ianthropogenic disturbance is turned off,  #
# we only need to run with one land use data set, since they would be always the same.     #
#------------------------------------------------------------------------------------------#
forbidden = joborder$iage == 1 & joborder$ibigleaf == 1
if ( "ianth.dataset" %in% names(varrun)){
   lu.1st    = varrun$ianth.dataset[1]
   redundant = joborder$ianth.disturb == 0 & joborder$ianth.dataset != lu.1st
   forbidden = forbiden | redundant
}#end if ( "ianth.dataset" %in% names(varrun))
joborder     = joborder[! forbidden,]
myruns       = myruns[! forbidden,]
nruns        = nrow(myruns)
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Check whether to replace information by site-specific information if no default has  #
# been provided (only if lonlat is NULL).                                                  #
#------------------------------------------------------------------------------------------#
if (is.null(lonlat)){
   #---------------------------------------------------------------------------------------#
   #     Create a table with room for all simulations.                                     #
   #---------------------------------------------------------------------------------------#
   poitout = poilist[match(joborder$iata,poilist$iata),]
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Replace met driver when it is supposed to be tower data.                          #
   #---------------------------------------------------------------------------------------#
   is.tower                      = joborder$met.driver == "tower"
   joborder$met.driver[is.tower] = poitout$met.driver[is.tower]
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Fix the initial year in case it was a negative number.                             #
   #---------------------------------------------------------------------------------------#
   if ("yeara" %in% names(joborder)){
      joborder$yeara = ifelse( test = joborder$yeara > 0
                             , yes  = joborder$yeara
                             , no   = poitout$yeara + 5 + joborder$yeara
                             )#end ifelse
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Replace other POI-specific variables as long as they are not to be specified by   #
   # the user settings.                                                                    #
   #---------------------------------------------------------------------------------------#
   keep    = ( names(poitout) %in% names(joborder) 
             & ( ! names(poitout) %in% names(myruns) )
             & ( ! names(poitout) %in% c("iata","met.driver") )
             )#end keep
   poidata = poitout[,keep]
   npois   = ncol(poidata)
   for (p in sequence(npois)){
      joborder[[names(poidata)[p]]] = poidata[[names(poidata)[p]]]
   }#end for (p in sequence(npois))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Make sure that soil texture is standard when the user provides soil texture and   #
   # the value is not zero.                                                                #
   #---------------------------------------------------------------------------------------#
   if ("istext" %in% names(varrun)){
      sel                      = joborder$istext != 0
      joborder$sand    [  sel] = -1.0
      joborder$clay    [  sel] = -1.0
      joborder$istext  [! sel] = poidata$istext[! sel]
      joborder$isoilflg[  sel] = 2
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Build the name of the simulations.  The polygon name always stays, even if the run #
   # is for one site only.                                                                 #
   #---------------------------------------------------------------------------------------#
   if (defjob){
      runname = defname[! forbidden]
      metname = ""
   }else{
      stay     = which(names(varlabel) %in% c("iata"))
      bye      = which(sapply(X=varlabel,FUN=length) == 1)
      bye      = bye[! bye %in% stay]
      if (length(bye) > 0) for (b in sort(bye,decreasing=TRUE)) varlabel[[b]] = NULL
      runname  = apply( X        = expand.grid(varlabel, stringsAsFactors = FALSE)
                      , MARGIN   = 1
                      , FUN      = paste
                      , collapse = "_"
                      )#end apply
      runname  = runname[! forbidden]
      metname  = ifelse(is.tower,"t","s")
   }#end if
   joborder$run      = paste(metname,runname,sep="")
   #---------------------------------------------------------------------------------------#
}else{
   #---------------------------------------------------------------------------------------#
   #     Job name by appending longitude and latitude.                                     #
   #---------------------------------------------------------------------------------------#
   joborder$run = paste( joborder$iata
                       , "_lon",sprintf("%+06.2f",joborder$lon)
                       , "_lat",sprintf("%+06.2f",joborder$lat)
                       , sep = ""
                       )#end paste
   #---------------------------------------------------------------------------------------#
}#end if
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     In case we must append the file, we read the previous file and merge them before     #
# writing, so they follow a single format.                                                 #
#------------------------------------------------------------------------------------------#
if (append.job){
   jobprev    = read.table( file             = outfile
                          , skip             = 3
                          , stringsAsFactors = FALSE
                          )#end read.table
   names(jobprev) = names(joborder)
   joborder       = rbind(jobprev,joborder)

}else{

}#end append.job
#------------------------------------------------------------------------------------------#


#----- Make sure days, months, and hours are properly formatted. --------------------------#
joborder$montha = sprintf("%2.2i",as.numeric(joborder$montha))
joborder$daya   = sprintf("%2.2i",as.numeric(joborder$daya  ))
joborder$timea  = sprintf("%4.4i",as.numeric(joborder$timea ))
joborder$monthz = sprintf("%2.2i",as.numeric(joborder$monthz))
joborder$dayz   = sprintf("%2.2i",as.numeric(joborder$dayz  ))
joborder$timez  = sprintf("%4.4i",as.numeric(joborder$timez ))
#------------------------------------------------------------------------------------------#



#----- Create new joborder file. ----------------------------------------------------------#
dash       = sapply( FUN      = paste
                   , X        = mapply( FUN      = rep
                                      , times    = nchar(names(joborder))
                                      , MoreArgs = list(x="-")
                                      )#end mapply
                   , collapse = "")
jobneat    = format(rbind(dash,toupper(names(joborder)),dash,joborder),justify="right"
                   ,width=nchar(names(joborder)))
jobneat    = apply(X=jobneat,MARGIN=1,FUN=paste,collapse="  ")
jobneat[1] = gsub(pattern=" ",replacement="-",x=jobneat[1])
jobneat[3] = gsub(pattern=" ",replacement="-",x=jobneat[3])
dum        = write.table( x         = jobneat
                        , file      = outfile
                        , quote     = FALSE
                        , row.names = FALSE
                        , col.names = FALSE
                        )#end write.table
#------------------------------------------------------------------------------------------#
