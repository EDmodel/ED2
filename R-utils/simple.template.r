#==========================================================================================#
#==========================================================================================#
#     Create the simplified monthly mean structure to be filled by plot_monthly.r          #
#------------------------------------------------------------------------------------------#
create.simple <<- function(ntimes,montha,yeara,inpref,slz.min){

   #----- Read the first HDF5 to grab some simulation-dependent dimensions. ---------------#
   cyear        = sprintf("%4.4i",yeara )
   cmonth       = sprintf("%2.2i",montha)
   h5first      = paste(inpref,"-Q-",cyear,"-",cmonth,"-00-000000-g01.h5"    ,sep="")
   h5first.bz2  = paste(inpref,"-Q-",cyear,"-",cmonth,"-00-000000-g01.h5.bz2",sep="")
   h5first.gz   = paste(inpref,"-Q-",cyear,"-",cmonth,"-00-000000-g01.h5.gz" ,sep="")
   if ( file.exists(h5first) ){
      mymont    = hdf5load(file=h5first,load=FALSE,verbosity=0,tidy=TRUE)

   }else if ( file.exists(h5first.bz2) ){
      temp.file = file.path(tempdir(),basename(h5first))
      dummy     = bunzip2(filename=h5first.bz2,destname=temp.file,remove=FALSE)
      mymont    = hdf5load(file=temp.file,load=FALSE,verbosity=0,tidy=TRUE)
      dummy     = file.remove(temp.file)

   }else if ( file.exists(h5first.gz) ){
      temp.file = file.path(tempdir(),basename(h5first))
      dummy     = gunzip(filename=h5first.gz,destname=temp.file,remove=FALSE)
      mymont    = hdf5load(file=temp.file,load=FALSE,verbosity=0,tidy=TRUE)
      dummy     = file.remove(temp.file)

   }else{
      cat (" Path: ",dirname (h5first),"\n")
      cat (" File: ",basename(h5first),"\n")
      stop(" File not found...")

   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Start up the list.                                                                #
   #---------------------------------------------------------------------------------------#
   ed      = list()
   #---------------------------------------------------------------------------------------#


   #----- Define the number of soil layers. -----------------------------------------------#
   ed$nzg        = mymont$NZG
   ed$nzs        = mymont$NZS
   ed$ndcycle    = mymont$NDCYCLE
   ed$ntimes     = ntimes
   #---------------------------------------------------------------------------------------#



   #----- Find which soil are we solving, and save properties into soil.prop. -------------#
   ed$isoilflg   = mymont$ISOILFLG
   ed$slz        = mymont$SLZ
   ed$slxsand    = mymont$SLXSAND
   ed$slxclay    = mymont$SLXCLAY
   ed$ntext      = mymont$NTEXT.SOIL[ed$nzg]
   #---------------------------------------------------------------------------------------#



   #----- Derive the soil properties. -----------------------------------------------------#
   ed$soil.prop  = soil.params(ed$ntext,ed$isoilflg,ed$slxsand,ed$slxclay)
   ed$dslz       = diff(c(ed$slz,0))
   ed$soil.depth = rev(cumsum(rev(ed$dslz)))
   ed$soil.dry   = rev(cumsum(rev(ed$soil.prop$soilcp * wdns * ed$dslz)))
   ed$soil.poro  = rev(cumsum(rev(ed$soil.prop$slmsts * wdns * ed$dslz)))
   #---------------------------------------------------------------------------------------#


   #----- Find the layers we care about. --------------------------------------------------#
   sel        = ed$slz < slz.min
   if (any(sel)){
      ed$ka      = which.max(ed$slz[sel])
   }else{
      ed$ka      = 1
   }#end if
   ed$kz         = ed$nzg
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Find all time information.                                                        #
   #---------------------------------------------------------------------------------------#
   runmonths         = montha + sequence(ntimes) - 1
   ed$month         = 1 + (runmonths-1) %% 12
   ed$year          = yeara - 1 + ceiling(runmonths/12)
   ed$when          = chron(paste(ed$month,1,ed$year,sep="/"))
   ed$tomonth       = chron(ed$when,out.format=c(dates="day-mon-yr",times=NULL))
   ed$toyear        = sort(unique(ed$year))
   #----- Count the number of months for each month. --------------------------------------#
   ed$montable      = rep(0,times=12)
   montable            = table(ed$month)
   idx                 = as.numeric(names(montable))
   ed$montable[idx] = montable
   ed$moncnt        = matrix( data = rep(ed$montable,times=ed$ndcycle)
                               , ncol = ed$ndcycle
                               , nrow = 12
                               )#end matrix
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find all the file names.                                                          #
   #---------------------------------------------------------------------------------------#
   cmonth   = sprintf("%2.2i",ed$month)
   cyear    = sprintf("%4.4i",ed$year )
   ed$input = paste(inpref,"-Q-",cyear,"-",cmonth,"-00-000000-g01.h5",sep="")
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Make a copy of the dimensions to avoid clutter.                                    #
   #---------------------------------------------------------------------------------------#
   ndcycle  = ed$ndcycle
   nzg      = ed$nzg
   nzs      = ed$nzs
   #---------------------------------------------------------------------------------------#



   #=======================================================================================#
   #=======================================================================================#
   #     Flush all variables that will hold the data.  For convenience they are split into #
   # multiple lists.                                                                       #
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   # emean -- variables that we can either compare directly with observations, or are      #
   #          or that may be used to draw time series.   They don't need to be really      #
   #          monthly means, but you should put only the variables that make sense to be   #
   #          plotted in simple time series (with no PFT or DBH information).              #
   #---------------------------------------------------------------------------------------#
   emean = list()
   emean$fast.soil.c             = rep(NA,times=ntimes)
   emean$slow.soil.c             = rep(NA,times=ntimes)
   emean$struct.soil.c           = rep(NA,times=ntimes)
   emean$het.resp                = rep(NA,times=ntimes)
   emean$cwd.resp                = rep(NA,times=ntimes)
   emean$gpp                     = rep(NA,times=ntimes)
   emean$npp                     = rep(NA,times=ntimes)
   emean$plant.resp              = rep(NA,times=ntimes)
   emean$leaf.resp               = rep(NA,times=ntimes)
   emean$root.resp               = rep(NA,times=ntimes)
   emean$growth.resp             = rep(NA,times=ntimes)
   emean$reco                    = rep(NA,times=ntimes)
   emean$nep                     = rep(NA,times=ntimes)
   emean$nee                     = rep(NA,times=ntimes)
   emean$ustar                   = rep(NA,times=ntimes)
   emean$cflxca                  = rep(NA,times=ntimes)
   emean$cflxst                  = rep(NA,times=ntimes)
   emean$ustar                   = rep(NA,times=ntimes)
   emean$atm.vels                = rep(NA,times=ntimes)
   emean$atm.prss                = rep(NA,times=ntimes)
   emean$atm.temp                = rep(NA,times=ntimes)
   emean$atm.shv                 = rep(NA,times=ntimes)
   emean$atm.co2                 = rep(NA,times=ntimes)
   emean$atm.vpd                 = rep(NA,times=ntimes)
   emean$can.prss                = rep(NA,times=ntimes)
   emean$can.temp                = rep(NA,times=ntimes)
   emean$can.shv                 = rep(NA,times=ntimes)
   emean$can.co2                 = rep(NA,times=ntimes)
   emean$can.vpd                 = rep(NA,times=ntimes)
   emean$gnd.temp                = rep(NA,times=ntimes)
   emean$gnd.shv                 = rep(NA,times=ntimes)
   emean$leaf.temp               = rep(NA,times=ntimes)
   emean$leaf.water              = rep(NA,times=ntimes)
   emean$leaf.vpd                = rep(NA,times=ntimes)
   emean$wood.temp               = rep(NA,times=ntimes)
   emean$hflxca                  = rep(NA,times=ntimes)
   emean$qwflxca                 = rep(NA,times=ntimes)
   emean$hflxgc                  = rep(NA,times=ntimes)
   emean$hflxlc                  = rep(NA,times=ntimes)
   emean$hflxwc                  = rep(NA,times=ntimes)
   emean$wflxca                  = rep(NA,times=ntimes)
   emean$wflxgc                  = rep(NA,times=ntimes)
   emean$wflxlc                  = rep(NA,times=ntimes)
   emean$wflxwc                  = rep(NA,times=ntimes)
   emean$runoff                  = rep(NA,times=ntimes)
   emean$intercepted             = rep(NA,times=ntimes)
   emean$wshed                   = rep(NA,times=ntimes)
   emean$evap                    = rep(NA,times=ntimes)
   emean$transp                  = rep(NA,times=ntimes)
   emean$et                      = rep(NA,times=ntimes)
   emean$rain                    = rep(NA,times=ntimes)
   emean$sm.stress               = rep(NA,times=ntimes)
   emean$rshort                  = rep(NA,times=ntimes)
   emean$rshort.beam             = rep(NA,times=ntimes)
   emean$rshort.diff             = rep(NA,times=ntimes)
   emean$rshortup                = rep(NA,times=ntimes)
   emean$rlong                   = rep(NA,times=ntimes)
   emean$rshort.gnd              = rep(NA,times=ntimes)
   emean$rlong.gnd               = rep(NA,times=ntimes)
   emean$rlongup                 = rep(NA,times=ntimes)
   emean$par.tot                 = rep(NA,times=ntimes)
   emean$par.beam                = rep(NA,times=ntimes)
   emean$par.diff                = rep(NA,times=ntimes)
   emean$par.gnd                 = rep(NA,times=ntimes)
   emean$parup                   = rep(NA,times=ntimes)
   emean$rnet                    = rep(NA,times=ntimes)
   emean$albedo                  = rep(NA,times=ntimes)
   emean$albedo.par              = rep(NA,times=ntimes)
   emean$albedo.nir              = rep(NA,times=ntimes)
   emean$rlong.albedo            = rep(NA,times=ntimes)
   emean$leaf.gbw                = rep(NA,times=ntimes)
   emean$leaf.gsw                = rep(NA,times=ntimes)
   emean$wood.gbw                = rep(NA,times=ntimes)
   emean$mco                     = rep(NA,times=ntimes)
   emean$ldrop                   = rep(NA,times=ntimes)
   emean$paw                     = rep(NA,times=ntimes)
   emean$smpot                   = rep(NA,times=ntimes)
   emean$workload                = rep(NA,times=ntimes)
   emean$specwork                = rep(NA,times=ntimes)
   emean$agb                     = rep(NA,times=ntimes)
   emean$lai                     = rep(NA,times=ntimes)
   emean$wai                     = rep(NA,times=ntimes)
   emean$basarea                 = rep(NA,times=ntimes)
   emean$nmon.lt.090             = rep(NA,times=ntimes)
   emean$nmon.lt.100             = rep(NA,times=ntimes)
   emean$nmon.lt.110             = rep(NA,times=ntimes)
   emean$nmon.lt.120             = rep(NA,times=ntimes)
   emean$nmon.wdef               = rep(NA,times=ntimes)
   emean$nmon.mdef               = rep(NA,times=ntimes)
   #----- Soil variables. -----------------------------------------------------------------#
   emean$soil.water              = matrix(data=0,nrow=ntimes,ncol=nzg)
   emean$soil.temp               = matrix(data=0,nrow=ntimes,ncol=nzg)
   emean$soil.mstpot             = matrix(data=0,nrow=ntimes,ncol=nzg)
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   # emsqu -- mean sum of squares of polygon-level variable.                               #
   #---------------------------------------------------------------------------------------#
   emsqu                   = list()
   emsqu$gpp               = rep(NA,times=ntimes)
   emsqu$plant.resp        = rep(NA,times=ntimes)
   emsqu$leaf.resp         = rep(NA,times=ntimes)
   emsqu$root.resp         = rep(NA,times=ntimes)
   emsqu$het.resp          = rep(NA,times=ntimes)
   emsqu$cwd.resp          = rep(NA,times=ntimes)
   emsqu$reco              = rep(NA,times=ntimes)
   emsqu$cflxca            = rep(NA,times=ntimes)
   emsqu$cflxst            = rep(NA,times=ntimes)
   emsqu$hflxca            = rep(NA,times=ntimes)
   emsqu$hflxlc            = rep(NA,times=ntimes)
   emsqu$hflxwc            = rep(NA,times=ntimes)
   emsqu$hflxgc            = rep(NA,times=ntimes)
   emsqu$wflxca            = rep(NA,times=ntimes)
   emsqu$qwflxca           = rep(NA,times=ntimes)
   emsqu$wflxlc            = rep(NA,times=ntimes)
   emsqu$wflxwc            = rep(NA,times=ntimes)
   emsqu$wflxgc            = rep(NA,times=ntimes)
   emsqu$evap              = rep(NA,times=ntimes)
   emsqu$transp            = rep(NA,times=ntimes)
   emsqu$ustar             = rep(NA,times=ntimes)
   emsqu$albedo            = rep(NA,times=ntimes)
   emsqu$rshortup          = rep(NA,times=ntimes)
   emsqu$rlongup           = rep(NA,times=ntimes)
   emsqu$parup             = rep(NA,times=ntimes)
   emsqu$rnet              = rep(NA,times=ntimes)
   #---------------------------------------------------------------------------------------#









   #---------------------------------------------------------------------------------------#
   # QMEAN -- Polygon-level variables, containing the mean diel (diurnal cycle).           #
   #---------------------------------------------------------------------------------------#
   qmean                = list()
   qmean$gpp            = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$npp            = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$plant.resp     = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$leaf.resp      = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$root.resp      = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$growth.resp    = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$het.resp       = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$cwd.resp       = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$nep            = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$nee            = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$reco           = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$cflxca         = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$cflxst         = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$hflxca         = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$hflxlc         = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$hflxwc         = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$hflxgc         = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$qwflxca        = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$wflxca         = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$wflxlc         = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$wflxwc         = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$wflxgc         = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$runoff         = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$intercepted    = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$wshed          = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$evap           = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$transp         = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$atm.temp       = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$can.temp       = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$leaf.temp      = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$leaf.water     = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$wood.temp      = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$gnd.temp       = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$atm.shv        = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$can.shv        = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$gnd.shv        = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$atm.vpd        = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$can.vpd        = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$leaf.vpd       = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$atm.co2        = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$can.co2        = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$atm.prss       = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$can.prss       = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$atm.vels       = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$ustar          = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$sm.stress      = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$rain           = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$rshort         = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$rshort.beam    = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$rshort.diff    = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$rshort.gnd     = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$rshortup       = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$rlong          = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$rlong.gnd      = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$rlongup        = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$par.tot        = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$par.beam       = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$par.diff       = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$par.gnd        = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$parup          = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$rnet           = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$albedo         = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$albedo.par     = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$albedo.nir     = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$rlong.albedo   = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$leaf.gsw       = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$leaf.gbw       = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$wood.gbw       = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   #---------------------------------------------------------------------------------------#









   #---------------------------------------------------------------------------------------#
   # QMSQU -- Polygon-level variables, containing the mean sum of squares for the diel     #
   #          (diurnal cycle).                                                             #
   #---------------------------------------------------------------------------------------#
   qmsqu                = list()
   qmsqu$gpp            = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmsqu$npp            = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmsqu$plant.resp     = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmsqu$leaf.resp      = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmsqu$root.resp      = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmsqu$het.resp       = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmsqu$cwd.resp       = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmsqu$nep            = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmsqu$reco           = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmsqu$cflxca         = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmsqu$cflxst         = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmsqu$hflxca         = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmsqu$hflxlc         = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmsqu$hflxwc         = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmsqu$hflxgc         = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmsqu$qwflxca        = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmsqu$wflxca         = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmsqu$wflxlc         = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmsqu$wflxwc         = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmsqu$wflxgc         = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmsqu$transp         = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmsqu$ustar          = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmsqu$albedo         = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmsqu$rshortup       = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmsqu$rlongup        = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmsqu$parup          = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmsqu$rnet           = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   #---------------------------------------------------------------------------------------#





   #----- Copy the polygon-level variable to the main structure. --------------------------#
   ed$emean  = emean
   ed$emsqu  = emsqu
   ed$qmean  = qmean
   ed$qmsqu  = qmsqu
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   return(ed)
   #---------------------------------------------------------------------------------------#

}#end create.simple
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     Expand the monthly array so it fits the new times.                                   #
#------------------------------------------------------------------------------------------#
update.simple <<- function(new.ntimes,old.datum,montha,yeara,inpref,slz.min){

   #----- Create the new data set. --------------------------------------------------------#
   new.datum = create.monthly(new.ntimes,montha,yeara,inpref,slz.min)
   #---------------------------------------------------------------------------------------#


   #----- Find out which times to copy. ---------------------------------------------------#
   idx = match(old.datum$when,new.datum$when)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   # emean -- variables that we can either compare directly with observations, or are      #
   #          or that may be used to draw time series.   They don't need to be really      #
   #          monthly means, but you should put only the variables that make sense to be   #
   #          plotted in simple time series (with no PFT or DBH information).              #
   #---------------------------------------------------------------------------------------#
   new.datum$emean$fast.soil.c    [idx ] = old.datum$emean$fast.soil.c
   new.datum$emean$slow.soil.c    [idx ] = old.datum$emean$slow.soil.c
   new.datum$emean$struct.soil.c  [idx ] = old.datum$emean$struct.soil.c
   new.datum$emean$het.resp       [idx ] = old.datum$emean$het.resp
   new.datum$emean$cwd.resp       [idx ] = old.datum$emean$cwd.resp
   new.datum$emean$gpp            [idx ] = old.datum$emean$gpp
   new.datum$emean$npp            [idx ] = old.datum$emean$npp
   new.datum$emean$plant.resp     [idx ] = old.datum$emean$plant.resp
   new.datum$emean$leaf.resp      [idx ] = old.datum$emean$leaf.resp
   new.datum$emean$root.resp      [idx ] = old.datum$emean$root.resp
   new.datum$emean$growth.resp    [idx ] = old.datum$emean$growth.resp
   new.datum$emean$reco           [idx ] = old.datum$emean$reco
   new.datum$emean$nep            [idx ] = old.datum$emean$nep
   new.datum$emean$nee            [idx ] = old.datum$emean$nee
   new.datum$emean$ustar          [idx ] = old.datum$emean$ustar
   new.datum$emean$cflxca         [idx ] = old.datum$emean$cflxca
   new.datum$emean$cflxst         [idx ] = old.datum$emean$cflxst
   new.datum$emean$ustar          [idx ] = old.datum$emean$ustar
   new.datum$emean$atm.vels       [idx ] = old.datum$emean$atm.vels
   new.datum$emean$atm.prss       [idx ] = old.datum$emean$atm.prss
   new.datum$emean$atm.temp       [idx ] = old.datum$emean$atm.temp
   new.datum$emean$atm.shv        [idx ] = old.datum$emean$atm.shv
   new.datum$emean$atm.co2        [idx ] = old.datum$emean$atm.co2
   new.datum$emean$atm.vpd        [idx ] = old.datum$emean$atm.vpd
   new.datum$emean$can.prss       [idx ] = old.datum$emean$can.prss
   new.datum$emean$can.temp       [idx ] = old.datum$emean$can.temp
   new.datum$emean$can.shv        [idx ] = old.datum$emean$can.shv
   new.datum$emean$can.co2        [idx ] = old.datum$emean$can.co2
   new.datum$emean$can.vpd        [idx ] = old.datum$emean$can.vpd
   new.datum$emean$gnd.temp       [idx ] = old.datum$emean$gnd.temp
   new.datum$emean$gnd.shv        [idx ] = old.datum$emean$gnd.shv
   new.datum$emean$leaf.temp      [idx ] = old.datum$emean$leaf.temp
   new.datum$emean$leaf.water     [idx ] = old.datum$emean$leaf.water
   new.datum$emean$leaf.vpd       [idx ] = old.datum$emean$leaf.vpd
   new.datum$emean$wood.temp      [idx ] = old.datum$emean$wood.temp
   new.datum$emean$hflxca         [idx ] = old.datum$emean$hflxca
   new.datum$emean$qwflxca        [idx ] = old.datum$emean$qwflxca
   new.datum$emean$hflxgc         [idx ] = old.datum$emean$hflxgc
   new.datum$emean$hflxlc         [idx ] = old.datum$emean$hflxlc
   new.datum$emean$hflxwc         [idx ] = old.datum$emean$hflxwc
   new.datum$emean$wflxca         [idx ] = old.datum$emean$wflxca
   new.datum$emean$wflxgc         [idx ] = old.datum$emean$wflxgc
   new.datum$emean$wflxlc         [idx ] = old.datum$emean$wflxlc
   new.datum$emean$wflxwc         [idx ] = old.datum$emean$wflxwc
   new.datum$emean$runoff         [idx ] = old.datum$emean$runoff
   new.datum$emean$intercepted    [idx ] = old.datum$emean$intercepted
   new.datum$emean$wshed          [idx ] = old.datum$emean$wshed
   new.datum$emean$evap           [idx ] = old.datum$emean$evap
   new.datum$emean$transp         [idx ] = old.datum$emean$transp
   new.datum$emean$et             [idx ] = old.datum$emean$et
   new.datum$emean$rain           [idx ] = old.datum$emean$rain
   new.datum$emean$sm.stress      [idx ] = old.datum$emean$sm.stress
   new.datum$emean$rshort         [idx ] = old.datum$emean$rshort
   new.datum$emean$rshort.beam    [idx ] = old.datum$emean$rshort.beam
   new.datum$emean$rshort.diff    [idx ] = old.datum$emean$rshort.diff
   new.datum$emean$rshortup       [idx ] = old.datum$emean$rshortup
   new.datum$emean$rlong          [idx ] = old.datum$emean$rlong
   new.datum$emean$rshort.gnd     [idx ] = old.datum$emean$rshort.gnd
   new.datum$emean$rlong.gnd      [idx ] = old.datum$emean$rlong.gnd
   new.datum$emean$rlongup        [idx ] = old.datum$emean$rlongup
   new.datum$emean$par.tot        [idx ] = old.datum$emean$par.tot
   new.datum$emean$par.beam       [idx ] = old.datum$emean$par.beam
   new.datum$emean$par.diff       [idx ] = old.datum$emean$par.diff
   new.datum$emean$par.gnd        [idx ] = old.datum$emean$par.gnd
   new.datum$emean$parup          [idx ] = old.datum$emean$parup
   new.datum$emean$rnet           [idx ] = old.datum$emean$rnet
   new.datum$emean$albedo         [idx ] = old.datum$emean$albedo
   new.datum$emean$albedo.par     [idx ] = old.datum$emean$albedo.par
   new.datum$emean$albedo.nir     [idx ] = old.datum$emean$albedo.nir
   new.datum$emean$rlong.albedo   [idx ] = old.datum$emean$rlong.albedo
   new.datum$emean$leaf.gbw       [idx ] = old.datum$emean$leaf.gbw
   new.datum$emean$leaf.gsw       [idx ] = old.datum$emean$leaf.gsw
   new.datum$emean$wood.gbw       [idx ] = old.datum$emean$wood.gbw
   new.datum$emean$mco            [idx ] = old.datum$emean$mco
   new.datum$emean$ldrop          [idx ] = old.datum$emean$ldrop
   new.datum$emean$paw            [idx ] = old.datum$emean$paw
   new.datum$emean$smpot          [idx ] = old.datum$emean$smpot
   new.datum$emean$workload       [idx ] = old.datum$emean$workload
   new.datum$emean$specwork       [idx ] = old.datum$emean$specwork
   new.datum$emean$agb            [idx ] = old.datum$emean$agb
   new.datum$emean$lai            [idx ] = old.datum$emean$lai
   new.datum$emean$wai            [idx ] = old.datum$emean$wai
   new.datum$emean$basarea        [idx ] = old.datum$emean$basarea
   new.datum$emean$nmon.lt.090    [idx ] = old.datum$emean$nmon.lt.090
   new.datum$emean$nmon.lt.100    [idx ] = old.datum$emean$nmon.lt.100
   new.datum$emean$nmon.lt.110    [idx ] = old.datum$emean$nmon.lt.110
   new.datum$emean$nmon.lt.120    [idx ] = old.datum$emean$nmon.lt.120
   new.datum$emean$nmon.wdef      [idx ] = old.datum$emean$nmon.wdef
   new.datum$emean$nmon.mdef      [idx ] = old.datum$emean$nmon.mdef
   #----- Soil variables. -----------------------------------------------------------------#
   new.datum$emean$soil.water     [idx,] = old.datum$emean$soil.water
   new.datum$emean$soil.temp      [idx,] = old.datum$emean$soil.temp
   new.datum$emean$soil.mstpot    [idx,] = old.datum$emean$soil.mstpot
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   # emsqu -- mean sum of squares of polygon-level variable.                               #
   #---------------------------------------------------------------------------------------#
   new.datum$emsqu$gpp            [idx] = old.datum$emsqu$gpp
   new.datum$emsqu$plant.resp     [idx] = old.datum$emsqu$plant.resp
   new.datum$emsqu$het.resp       [idx] = old.datum$emsqu$het.resp
   new.datum$emsqu$cwd.resp       [idx] = old.datum$emsqu$cwd.resp
   new.datum$emsqu$cflxca         [idx] = old.datum$emsqu$cflxca
   new.datum$emsqu$cflxst         [idx] = old.datum$emsqu$cflxst
   new.datum$emsqu$hflxca         [idx] = old.datum$emsqu$hflxca
   new.datum$emsqu$hflxlc         [idx] = old.datum$emsqu$hflxlc
   new.datum$emsqu$hflxwc         [idx] = old.datum$emsqu$hflxwc
   new.datum$emsqu$hflxgc         [idx] = old.datum$emsqu$hflxgc
   new.datum$emsqu$wflxca         [idx] = old.datum$emsqu$wflxca
   new.datum$emsqu$qwflxca        [idx] = old.datum$emsqu$qwflxca
   new.datum$emsqu$wflxlc         [idx] = old.datum$emsqu$wflxlc
   new.datum$emsqu$wflxwc         [idx] = old.datum$emsqu$wflxwc
   new.datum$emsqu$wflxgc         [idx] = old.datum$emsqu$wflxgc
   new.datum$emsqu$transp         [idx] = old.datum$emsqu$transp
   new.datum$emsqu$ustar          [idx] = old.datum$emsqu$ustar
   new.datum$emsqu$albedo         [idx] = old.datum$emsqu$albedo
   new.datum$emsqu$rshortup       [idx] = old.datum$emsqu$rshortup
   new.datum$emsqu$rlongup        [idx] = old.datum$emsqu$rlongup
   new.datum$emsqu$parup          [idx] = old.datum$emsqu$parup
   new.datum$emsqu$rnet           [idx] = old.datum$emsqu$rnet
   #---------------------------------------------------------------------------------------#









   #---------------------------------------------------------------------------------------#
   # QMEAN -- Polygon-level variables, containing the mean diel (diurnal cycle).           #
   #---------------------------------------------------------------------------------------#
   new.datum$qmean$gpp           [idx,] = old.datum$qmean$gpp
   new.datum$qmean$npp           [idx,] = old.datum$qmean$npp
   new.datum$qmean$plant.resp    [idx,] = old.datum$qmean$plant.resp
   new.datum$qmean$leaf.resp     [idx,] = old.datum$qmean$leaf.resp
   new.datum$qmean$root.resp     [idx,] = old.datum$qmean$root.resp
   new.datum$qmean$growth.resp   [idx,] = old.datum$qmean$growth.resp
   new.datum$qmean$het.resp      [idx,] = old.datum$qmean$het.resp
   new.datum$qmean$cwd.resp      [idx,] = old.datum$qmean$cwd.resp
   new.datum$qmean$nep           [idx,] = old.datum$qmean$nep
   new.datum$qmean$nee           [idx,] = old.datum$qmean$nee
   new.datum$qmean$reco          [idx,] = old.datum$qmean$reco
   new.datum$qmean$cflxca        [idx,] = old.datum$qmean$cflxca
   new.datum$qmean$cflxst        [idx,] = old.datum$qmean$cflxst
   new.datum$qmean$hflxca        [idx,] = old.datum$qmean$hflxca
   new.datum$qmean$hflxlc        [idx,] = old.datum$qmean$hflxlc
   new.datum$qmean$hflxwc        [idx,] = old.datum$qmean$hflxwc
   new.datum$qmean$hflxgc        [idx,] = old.datum$qmean$hflxgc
   new.datum$qmean$qwflxca       [idx,] = old.datum$qmean$qwflxca
   new.datum$qmean$wflxca        [idx,] = old.datum$qmean$wflxca
   new.datum$qmean$wflxlc        [idx,] = old.datum$qmean$wflxlc
   new.datum$qmean$wflxwc        [idx,] = old.datum$qmean$wflxwc
   new.datum$qmean$wflxgc        [idx,] = old.datum$qmean$wflxgc
   new.datum$qmean$runoff        [idx,] = old.datum$qmean$runoff
   new.datum$qmean$intercepted   [idx,] = old.datum$qmean$intercepted
   new.datum$qmean$wshed         [idx,] = old.datum$qmean$wshed
   new.datum$qmean$evap          [idx,] = old.datum$qmean$evap
   new.datum$qmean$transp        [idx,] = old.datum$qmean$transp
   new.datum$qmean$atm.temp      [idx,] = old.datum$qmean$atm.temp
   new.datum$qmean$can.temp      [idx,] = old.datum$qmean$can.temp
   new.datum$qmean$leaf.temp     [idx,] = old.datum$qmean$leaf.temp
   new.datum$qmean$leaf.water    [idx,] = old.datum$qmean$leaf.water
   new.datum$qmean$wood.temp     [idx,] = old.datum$qmean$wood.temp
   new.datum$qmean$gnd.temp      [idx,] = old.datum$qmean$gnd.temp
   new.datum$qmean$atm.shv       [idx,] = old.datum$qmean$atm.shv
   new.datum$qmean$can.shv       [idx,] = old.datum$qmean$can.shv
   new.datum$qmean$gnd.shv       [idx,] = old.datum$qmean$gnd.shv
   new.datum$qmean$atm.vpd       [idx,] = old.datum$qmean$atm.vpd
   new.datum$qmean$can.vpd       [idx,] = old.datum$qmean$can.vpd
   new.datum$qmean$leaf.vpd      [idx,] = old.datum$qmean$leaf.vpd
   new.datum$qmean$atm.co2       [idx,] = old.datum$qmean$atm.co2
   new.datum$qmean$can.co2       [idx,] = old.datum$qmean$can.co2
   new.datum$qmean$atm.prss      [idx,] = old.datum$qmean$atm.prss
   new.datum$qmean$can.prss      [idx,] = old.datum$qmean$can.prss
   new.datum$qmean$atm.vels      [idx,] = old.datum$qmean$atm.vels
   new.datum$qmean$ustar         [idx,] = old.datum$qmean$ustar
   new.datum$qmean$sm.stress     [idx,] = old.datum$qmean$sm.stress
   new.datum$qmean$rain          [idx,] = old.datum$qmean$rain
   new.datum$qmean$rshort        [idx,] = old.datum$qmean$rshort
   new.datum$qmean$rshort.beam   [idx,] = old.datum$qmean$rshort.beam
   new.datum$qmean$rshort.diff   [idx,] = old.datum$qmean$rshort.diff
   new.datum$qmean$rshort.gnd    [idx,] = old.datum$qmean$rshort.gnd
   new.datum$qmean$rshortup      [idx,] = old.datum$qmean$rshortup
   new.datum$qmean$rlong         [idx,] = old.datum$qmean$rlong
   new.datum$qmean$rlong.gnd     [idx,] = old.datum$qmean$rlong.gnd
   new.datum$qmean$rlongup       [idx,] = old.datum$qmean$rlongup
   new.datum$qmean$par.tot       [idx,] = old.datum$qmean$par.tot
   new.datum$qmean$par.beam      [idx,] = old.datum$qmean$par.beam
   new.datum$qmean$par.diff      [idx,] = old.datum$qmean$par.diff
   new.datum$qmean$par.gnd       [idx,] = old.datum$qmean$par.gnd
   new.datum$qmean$parup         [idx,] = old.datum$qmean$parup
   new.datum$qmean$rnet          [idx,] = old.datum$qmean$rnet
   new.datum$qmean$albedo        [idx,] = old.datum$qmean$albedo
   new.datum$qmean$albedo.par    [idx,] = old.datum$qmean$albedo.par
   new.datum$qmean$albedo.nir    [idx,] = old.datum$qmean$albedo.nir
   new.datum$qmean$rlong.albedo  [idx,] = old.datum$qmean$rlong.albedo
   new.datum$qmean$leaf.gsw      [idx,] = old.datum$qmean$leaf.gsw
   new.datum$qmean$leaf.gbw      [idx,] = old.datum$qmean$leaf.gbw
   new.datum$qmean$wood.gbw      [idx,] = old.datum$qmean$wood.gbw
   #---------------------------------------------------------------------------------------#









   #---------------------------------------------------------------------------------------#
   # QMSQU -- Polygon-level variables, containing the mean sum of squares for the diel     #
   #          (diurnal cycle).                                                             #
   #---------------------------------------------------------------------------------------#
   new.datum$qmsqu$gpp           [idx,] = old.datum$qmsqu$gpp
   new.datum$qmsqu$npp           [idx,] = old.datum$qmsqu$npp
   new.datum$qmsqu$plant.resp    [idx,] = old.datum$qmsqu$plant.resp
   new.datum$qmsqu$het.resp      [idx,] = old.datum$qmsqu$het.resp
   new.datum$qmsqu$cwd.resp      [idx,] = old.datum$qmsqu$cwd.resp
   new.datum$qmsqu$nep           [idx,] = old.datum$qmsqu$nep
   new.datum$qmsqu$cflxca        [idx,] = old.datum$qmsqu$cflxca
   new.datum$qmsqu$cflxst        [idx,] = old.datum$qmsqu$cflxst
   new.datum$qmsqu$hflxca        [idx,] = old.datum$qmsqu$hflxca
   new.datum$qmsqu$hflxlc        [idx,] = old.datum$qmsqu$hflxlc
   new.datum$qmsqu$hflxwc        [idx,] = old.datum$qmsqu$hflxwc
   new.datum$qmsqu$hflxgc        [idx,] = old.datum$qmsqu$hflxgc
   new.datum$qmsqu$qwflxca       [idx,] = old.datum$qmsqu$qwflxca
   new.datum$qmsqu$wflxca        [idx,] = old.datum$qmsqu$wflxca
   new.datum$qmsqu$wflxlc        [idx,] = old.datum$qmsqu$wflxlc
   new.datum$qmsqu$wflxwc        [idx,] = old.datum$qmsqu$wflxwc
   new.datum$qmsqu$wflxgc        [idx,] = old.datum$qmsqu$wflxgc
   new.datum$qmsqu$transp        [idx,] = old.datum$qmsqu$transp
   new.datum$qmsqu$ustar         [idx,] = old.datum$qmsqu$ustar
   new.datum$qmsqu$albedo        [idx,] = old.datum$qmsqu$albedo
   new.datum$qmsqu$rshortup      [idx,] = old.datum$qmsqu$rshortup
   new.datum$qmsqu$rlongup       [idx,] = old.datum$qmsqu$rlongup
   new.datum$qmsqu$parup         [idx,] = old.datum$qmsqu$parup
   new.datum$qmsqu$rnet          [idx,] = old.datum$qmsqu$rnet
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Send the data back.                                                               #
   #---------------------------------------------------------------------------------------#
   return(new.datum)
   #---------------------------------------------------------------------------------------#
}#end function update.simple
#==========================================================================================#
#==========================================================================================#
