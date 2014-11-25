#==========================================================================================#
#==========================================================================================#
#     Create the monthly mean structure to be filled by plot_monthly.r                     #
#------------------------------------------------------------------------------------------#
create.monthly <<- function(ntimes,montha,yeara,inpref,slz.min){

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
   emean$wood.dens               = rep(NA,times=ntimes)
   emean$fast.soil.c             = rep(NA,times=ntimes)
   emean$slow.soil.c             = rep(NA,times=ntimes)
   emean$struct.soil.c           = rep(NA,times=ntimes)
   emean$het.resp                = rep(NA,times=ntimes)
   emean$cwd.resp                = rep(NA,times=ntimes)
   emean$soil.resp               = rep(NA,times=ntimes)
   emean$gpp                     = rep(NA,times=ntimes)
   emean$last.1yr.gpp            = rep(NA,times=ntimes)
   emean$last.2yr.gpp            = rep(NA,times=ntimes)
   emean$last.3yr.gpp            = rep(NA,times=ntimes)
   emean$plant.resp              = rep(NA,times=ntimes)
   emean$last.1yr.plresp         = rep(NA,times=ntimes)
   emean$last.2yr.plresp         = rep(NA,times=ntimes)
   emean$last.3yr.plresp         = rep(NA,times=ntimes)
   emean$cue                     = rep(NA,times=ntimes)
   emean$last.1yr.cue            = rep(NA,times=ntimes)
   emean$last.2yr.cue            = rep(NA,times=ntimes)
   emean$last.3yr.cue            = rep(NA,times=ntimes)
   emean$ecue                    = rep(NA,times=ntimes)
   emean$last.1yr.ecue           = rep(NA,times=ntimes)
   emean$last.2yr.ecue           = rep(NA,times=ntimes)
   emean$last.3yr.ecue           = rep(NA,times=ntimes)
   emean$leaf.resp               = rep(NA,times=ntimes)
   emean$root.resp               = rep(NA,times=ntimes)
   emean$froot.resp              = rep(NA,times=ntimes)
   emean$croot.resp              = rep(NA,times=ntimes)
   emean$stem.resp               = rep(NA,times=ntimes)
   emean$growth.resp             = rep(NA,times=ntimes)
   emean$storage.resp            = rep(NA,times=ntimes)
   emean$reco                    = rep(NA,times=ntimes)
   emean$assim.light             = rep(NA,times=ntimes)
   emean$assim.rubp              = rep(NA,times=ntimes)
   emean$assim.co2               = rep(NA,times=ntimes)
   emean$assim.ratio             = rep(NA,times=ntimes)
   emean$mco                     = rep(NA,times=ntimes)
   emean$cba                     = rep(NA,times=ntimes)
   emean$last.1yr.cba            = rep(NA,times=ntimes)
   emean$last.2yr.cba            = rep(NA,times=ntimes)
   emean$last.3yr.cba            = rep(NA,times=ntimes)
   emean$cbamax                  = rep(NA,times=ntimes)
   emean$cbalight                = rep(NA,times=ntimes)
   emean$cbamoist                = rep(NA,times=ntimes)
   emean$cbarel                  = rep(NA,times=ntimes)
   emean$ldrop                   = rep(NA,times=ntimes)
   emean$nep                     = rep(NA,times=ntimes)
   emean$nee                     = rep(NA,times=ntimes)
   emean$cflxca                  = rep(NA,times=ntimes)
   emean$cflxst                  = rep(NA,times=ntimes)
   emean$ustar                   = rep(NA,times=ntimes)
   emean$atm.vels                = rep(NA,times=ntimes)
   emean$atm.prss                = rep(NA,times=ntimes)
   emean$atm.temp                = rep(NA,times=ntimes)
   emean$atm.shv                 = rep(NA,times=ntimes)
   emean$atm.vpd                 = rep(NA,times=ntimes)
   emean$atm.co2                 = rep(NA,times=ntimes)
   emean$can.prss                = rep(NA,times=ntimes)
   emean$can.temp                = rep(NA,times=ntimes)
   emean$can.co2                 = rep(NA,times=ntimes)
   emean$can.shv                 = rep(NA,times=ntimes)
   emean$can.vpd                 = rep(NA,times=ntimes)
   emean$can.depth               = rep(NA,times=ntimes) 
   emean$can.area                = rep(NA,times=ntimes)
   emean$gnd.temp                = rep(NA,times=ntimes)
   emean$gnd.shv                 = rep(NA,times=ntimes)
   emean$leaf.temp               = rep(NA,times=ntimes)
   emean$phap.ltemp              = rep(NA,times=ntimes)
   emean$last.1yr.ltemp          = rep(NA,times=ntimes)
   emean$last.2yr.ltemp          = rep(NA,times=ntimes)
   emean$last.3yr.ltemp          = rep(NA,times=ntimes)
   emean$leaf.vpd                = rep(NA,times=ntimes)
   emean$phap.lvpd               = rep(NA,times=ntimes)
   emean$last.1yr.lvpd           = rep(NA,times=ntimes)
   emean$last.2yr.lvpd           = rep(NA,times=ntimes)
   emean$last.3yr.lvpd           = rep(NA,times=ntimes)
   emean$leaf.water              = rep(NA,times=ntimes)
   emean$phap.lwater             = rep(NA,times=ntimes)
   emean$last.1yr.lwater         = rep(NA,times=ntimes)
   emean$last.2yr.lwater         = rep(NA,times=ntimes)
   emean$last.3yr.lwater         = rep(NA,times=ntimes)
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
   emean$last.1yr.evap           = rep(NA,times=ntimes)
   emean$last.2yr.evap           = rep(NA,times=ntimes)
   emean$last.3yr.evap           = rep(NA,times=ntimes)
   emean$npp                     = rep(NA,times=ntimes)
   emean$last.1yr.npp            = rep(NA,times=ntimes)
   emean$last.2yr.npp            = rep(NA,times=ntimes)
   emean$last.3yr.npp            = rep(NA,times=ntimes)
   emean$dcbadt                  = rep(NA,times=ntimes)
   emean$last.1yr.dcbadt         = rep(NA,times=ntimes)
   emean$last.2yr.dcbadt         = rep(NA,times=ntimes)
   emean$last.3yr.dcbadt         = rep(NA,times=ntimes)
   emean$et                      = rep(NA,times=ntimes)
   emean$last.1yr.et             = rep(NA,times=ntimes)
   emean$last.2yr.et             = rep(NA,times=ntimes)
   emean$last.3yr.et             = rep(NA,times=ntimes)
   emean$transp                  = rep(NA,times=ntimes)
   emean$last.1yr.transp         = rep(NA,times=ntimes)
   emean$last.2yr.transp         = rep(NA,times=ntimes)
   emean$last.3yr.transp         = rep(NA,times=ntimes)
   emean$etue                    = rep(NA,times=ntimes)
   emean$last.1yr.etue           = rep(NA,times=ntimes)
   emean$last.2yr.etue           = rep(NA,times=ntimes)
   emean$last.3yr.etue           = rep(NA,times=ntimes)
   emean$wue                     = rep(NA,times=ntimes)
   emean$last.1yr.wue            = rep(NA,times=ntimes)
   emean$last.2yr.wue            = rep(NA,times=ntimes)
   emean$last.3yr.wue            = rep(NA,times=ntimes)
   emean$rain                    = rep(NA,times=ntimes)
   emean$last.1yr.rain           = rep(NA,times=ntimes)
   emean$last.2yr.rain           = rep(NA,times=ntimes)
   emean$last.3yr.rain           = rep(NA,times=ntimes)
   emean$nmon.lt.090             = rep(NA,times=ntimes)
   emean$nmon.lt.100             = rep(NA,times=ntimes)
   emean$nmon.lt.110             = rep(NA,times=ntimes)
   emean$nmon.lt.120             = rep(NA,times=ntimes)
   emean$nmon.wdef               = rep(NA,times=ntimes)
   emean$nmon.mdef               = rep(NA,times=ntimes)
   emean$rue                     = rep(NA,times=ntimes)
   emean$last.1yr.rue            = rep(NA,times=ntimes)
   emean$last.2yr.rue            = rep(NA,times=ntimes)
   emean$last.3yr.rue            = rep(NA,times=ntimes)
   emean$sm.stress               = rep(NA,times=ntimes)
   emean$phap.sms                = rep(NA,times=ntimes)
   emean$last.1yr.sms            = rep(NA,times=ntimes)
   emean$last.2yr.sms            = rep(NA,times=ntimes)
   emean$last.3yr.sms            = rep(NA,times=ntimes)
   emean$rshort                  = rep(NA,times=ntimes)
   emean$last.1yr.rshort         = rep(NA,times=ntimes)
   emean$last.2yr.rshort         = rep(NA,times=ntimes)
   emean$last.3yr.rshort         = rep(NA,times=ntimes)
   emean$rshort.beam             = rep(NA,times=ntimes)
   emean$rshort.diff             = rep(NA,times=ntimes)
   emean$rshortup                = rep(NA,times=ntimes)
   emean$rshort.gnd              = rep(NA,times=ntimes)
   emean$rlong                   = rep(NA,times=ntimes)
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
   emean$nplant                  = rep(NA,times=ntimes)
   emean$agb                     = rep(NA,times=ntimes)
   emean$biomass                 = rep(NA,times=ntimes)
   emean$lai                     = rep(NA,times=ntimes)
   emean$wai                     = rep(NA,times=ntimes)
   emean$tai                     = rep(NA,times=ntimes)
   emean$area                    = rep(NA,times=ntimes)
   emean$workload                = rep(NA,times=ntimes)
   emean$specwork                = rep(NA,times=ntimes)
   emean$demand                  = rep(NA,times=ntimes)
   emean$supply                  = rep(NA,times=ntimes)
   emean$paw                     = rep(NA,times=ntimes)
   emean$smpot                   = rep(NA,times=ntimes)
   emean$last.1yr.smpot          = rep(NA,times=ntimes)
   emean$last.2yr.smpot          = rep(NA,times=ntimes)
   emean$last.3yr.smpot          = rep(NA,times=ntimes)
   emean$npat.global             = rep(NA,times=ntimes)
   emean$ncoh.global             = rep(NA,times=ntimes)
   emean$water.deficit           = rep(NA,times=ntimes)
   emean$last.1yr.mwd            = rep(NA,times=ntimes)
   emean$last.2yr.mwd            = rep(NA,times=ntimes)
   emean$last.3yr.mwd            = rep(NA,times=ntimes)
   emean$malhi.deficit           = rep(NA,times=ntimes)
   emean$leaf.gsw                = rep(NA,times=ntimes)
   emean$phap.lgsw               = rep(NA,times=ntimes)
   emean$last.1yr.lgsw           = rep(NA,times=ntimes)
   emean$last.2yr.lgsw           = rep(NA,times=ntimes)
   emean$last.3yr.lgsw           = rep(NA,times=ntimes)
   emean$leaf.gbw                = rep(NA,times=ntimes)
   emean$phap.lgbw               = rep(NA,times=ntimes)
   emean$wood.gbw                = rep(NA,times=ntimes)
   emean$i.gpp                   = rep(NA,times=ntimes)
   emean$i.npp                   = rep(NA,times=ntimes)
   emean$i.plant.resp            = rep(NA,times=ntimes)
   emean$i.mco                   = rep(NA,times=ntimes)
   emean$i.cba                   = rep(NA,times=ntimes)
   emean$i.cbamax                = rep(NA,times=ntimes)
   emean$i.cbalight              = rep(NA,times=ntimes)
   emean$i.cbamoist              = rep(NA,times=ntimes)
   emean$i.transp                = rep(NA,times=ntimes)
   emean$i.wflxlc                = rep(NA,times=ntimes)
   emean$i.hflxlc                = rep(NA,times=ntimes)
   emean$f.gpp                   = rep(NA,times=ntimes)
   emean$f.plant.resp            = rep(NA,times=ntimes)
   emean$f.npp                   = rep(NA,times=ntimes)
   emean$f.mco                   = rep(NA,times=ntimes)
   emean$f.cba                   = rep(NA,times=ntimes)
   emean$f.bstorage              = rep(NA,times=ntimes)
   emean$f.bleaf                 = rep(NA,times=ntimes)
   emean$f.bstem                 = rep(NA,times=ntimes)
   emean$f.broot                 = rep(NA,times=ntimes)
   emean$f.bseeds                = rep(NA,times=ntimes)
   emean$f.dcbadt                = rep(NA,times=ntimes)
   emean$leaf.par                = rep(NA,times=ntimes)
   emean$phap.lpar               = rep(NA,times=ntimes)
   emean$last.1yr.lpar           = rep(NA,times=ntimes)
   emean$last.2yr.lpar           = rep(NA,times=ntimes)
   emean$last.3yr.lpar           = rep(NA,times=ntimes)
   emean$leaf.rshort             = rep(NA,times=ntimes)
   emean$leaf.rlong              = rep(NA,times=ntimes)
   emean$last.1yr.growth         = rep(NA,times=ntimes)
   emean$last.2yr.growth         = rep(NA,times=ntimes)
   emean$last.3yr.growth         = rep(NA,times=ntimes)
   emean$last.1yr.recr           = rep(NA,times=ntimes)
   emean$last.2yr.recr           = rep(NA,times=ntimes)
   emean$last.3yr.recr           = rep(NA,times=ntimes)
   emean$last.1yr.mort           = rep(NA,times=ntimes)
   emean$last.2yr.mort           = rep(NA,times=ntimes)
   emean$last.3yr.mort           = rep(NA,times=ntimes)
   emean$last.1yr.dimort         = rep(NA,times=ntimes)
   emean$last.2yr.dimort         = rep(NA,times=ntimes)
   emean$last.3yr.dimort         = rep(NA,times=ntimes)
   emean$last.1yr.ncbmort        = rep(NA,times=ntimes)
   emean$last.2yr.ncbmort        = rep(NA,times=ntimes)
   emean$last.3yr.ncbmort        = rep(NA,times=ntimes)
   emean$agb.change              = rep(NA,times=ntimes)
   emean$last.1yr.change         = rep(NA,times=ntimes)
   emean$last.2yr.change         = rep(NA,times=ntimes)
   emean$last.3yr.change         = rep(NA,times=ntimes)
   #----- Soil variables. -----------------------------------------------------------------#
   emean$soil.water              = matrix(data=0,nrow=ntimes,ncol=nzg)
   emean$soil.temp               = matrix(data=0,nrow=ntimes,ncol=nzg)
   emean$soil.mstpot             = matrix(data=0,nrow=ntimes,ncol=nzg)
   emean$soil.extracted          = matrix(data=0,nrow=ntimes,ncol=nzg)
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   # emsqu -- mean sum of squares of polygon-level variable.                               #
   #---------------------------------------------------------------------------------------#
   emsqu                   = list()
   emsqu$gpp               = rep(NA,times=ntimes)
   emsqu$plant.resp        = rep(NA,times=ntimes)
   emsqu$leaf.resp         = rep(NA,times=ntimes)
   emsqu$froot.resp        = rep(NA,times=ntimes)
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
   # SZPFT -- Size (DBH) and plant functional type (PFT) array.  An extra level is         #
   #          appended to the end, which will hold the sum of all categories.              #
   #                                                                                       #
   #          The initial value depends on the type of variable:                           #
   #          - If absence is equal to zero, then the initial value must be zero. This is  #
   #            normally the case for state variables (e.g. biomass, demographic density,  #
   #            LAI, etc.)                                                                 #
   #          - If absence makes the variable meaningless, then the initial value must     #
   #            be NA.  This is true for plant-derived properties (temperature, mortality, #
   #            gsw, etc.)                                                                 #
   #---------------------------------------------------------------------------------------#
   szpft = list()
   #----- Initial value should be zero. ---------------------------------------------------#
   szpft$agb               = array(data=0 ,dim=c(ntimes,ndbh+1,npft+1))
   szpft$biomass           = array(data=0 ,dim=c(ntimes,ndbh+1,npft+1))
   szpft$lai               = array(data=0 ,dim=c(ntimes,ndbh+1,npft+1))
   szpft$wai               = array(data=0 ,dim=c(ntimes,ndbh+1,npft+1))
   szpft$tai               = array(data=0 ,dim=c(ntimes,ndbh+1,npft+1))
   szpft$ba                = array(data=0 ,dim=c(ntimes,ndbh+1,npft+1))
   szpft$nplant            = array(data=0 ,dim=c(ntimes,ndbh+1,npft+1))
   szpft$bdead             = array(data=0 ,dim=c(ntimes,ndbh+1,npft+1))
   szpft$balive            = array(data=0 ,dim=c(ntimes,ndbh+1,npft+1))
   szpft$bleaf             = array(data=0 ,dim=c(ntimes,ndbh+1,npft+1))
   szpft$bstem             = array(data=0 ,dim=c(ntimes,ndbh+1,npft+1))
   szpft$broot             = array(data=0 ,dim=c(ntimes,ndbh+1,npft+1))
   szpft$bfroot            = array(data=0 ,dim=c(ntimes,ndbh+1,npft+1))
   szpft$bcroot            = array(data=0 ,dim=c(ntimes,ndbh+1,npft+1))
   szpft$bsapwood          = array(data=0 ,dim=c(ntimes,ndbh+1,npft+1))
   szpft$bstorage          = array(data=0 ,dim=c(ntimes,ndbh+1,npft+1))
   szpft$bseeds            = array(data=0 ,dim=c(ntimes,ndbh+1,npft+1))
   szpft$census.lai        = array(data=0 ,dim=c(ntimes,ndbh+1,npft+1))
   szpft$census.wai        = array(data=0 ,dim=c(ntimes,ndbh+1,npft+1))
   szpft$census.tai        = array(data=0 ,dim=c(ntimes,ndbh+1,npft+1))
   szpft$census.agb        = array(data=0 ,dim=c(ntimes,ndbh+1,npft+1))
   szpft$census.ba         = array(data=0 ,dim=c(ntimes,ndbh+1,npft+1))
   #----- Initial value should be NA. -----------------------------------------------------#
   szpft$wood.dens         = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$gpp               = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$npp               = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$dcbadt            = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$leaf.resp         = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$root.resp         = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$froot.resp        = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$croot.resp        = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$stem.resp         = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$growth.resp       = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$storage.resp      = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$plant.resp        = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$assim.light       = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$assim.rubp        = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$assim.co2         = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$assim.ratio       = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$mco               = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$cba               = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$cbamax            = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$cbalight          = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$cbamoist          = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$cbarel            = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$ldrop             = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$sm.stress         = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$phap.sms          = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$leaf.gbw          = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$phap.lgbw         = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$leaf.gsw          = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$phap.lgsw         = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$wood.gbw          = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$leaf.temp         = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$phap.ltemp        = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$leaf.water        = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$phap.lwater       = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$wood.temp         = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$leaf.vpd          = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$phap.lvpd         = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$demand            = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$supply            = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$mort              = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$dimort            = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$ncbmort           = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$growth            = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$recr              = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$change            = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$agb.mort          = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$agb.dimort        = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$agb.ncbmort       = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$agb.growth        = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$agb.recr          = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$agb.change        = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$bsa.mort          = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$bsa.dimort        = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$bsa.ncbmort       = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$bsa.growth        = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$bsa.recr          = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$bsa.change        = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$hflxlc            = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$wflxlc            = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$i.gpp             = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$i.npp             = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$i.plant.resp      = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$i.mco             = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$i.cba             = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$i.cbamax          = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$i.cbalight        = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$i.cbamoist        = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$i.transp          = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$i.wflxlc          = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$i.hflxlc          = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$f.gpp             = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$f.plant.resp      = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$f.npp             = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$f.mco             = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$f.cba             = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$f.bstorage        = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$f.dcbadt          = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$f.bleaf           = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$f.bstem           = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$f.broot           = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$f.bseeds          = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$leaf.par          = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$phap.lpar         = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$leaf.rshort       = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$leaf.rlong        = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$transp            = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$wue               = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$cue               = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$ecue              = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$etue              = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   szpft$rue               = array(data=NA,dim=c(ntimes,ndbh+1,npft+1))
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   # LU -- Polygon-level variables, split by land use type.  One extra dimension is        #
   #       appended to the end, which will hold the sum of all land use types.             #
   #                                                                                       #
   #       The initial value depends on the type of variable:                              #
   #       - If absence is equal to zero, then the initial value must be zero. This is     #
   #         normally the case for state variables (e.g. biomass, demographic density,     #
   #         LAI, etc.)                                                                    #
   #       - If absence makes the variable meaningless, then the initial value must        #
   #         be NA.  This is true for plant-derived properties (temperature, mortality,    #
   #         gsw, etc.)                                                                    #
   #---------------------------------------------------------------------------------------#
   lu           = list()
   #----- Initial value should be zero. ---------------------------------------------------#
   lu$agb       = matrix(data=0 ,nrow=ntimes,ncol=nlu+1)
   lu$biomass   = matrix(data=0 ,nrow=ntimes,ncol=nlu+1)
   lu$lai       = matrix(data=0 ,nrow=ntimes,ncol=nlu+1)
   lu$area      = matrix(data=0 ,nrow=ntimes,ncol=nlu+1)
   lu$ba        = matrix(data=0 ,nrow=ntimes,ncol=nlu+1)
   #----- Initial value should be NA. -----------------------------------------------------#
   lu$gpp       = matrix(data=NA,nrow=ntimes,ncol=nlu+1)
   lu$npp       = matrix(data=NA,nrow=ntimes,ncol=nlu+1)
   lu$f.agb     = matrix(data=NA,nrow=ntimes,ncol=nlu+1)
   lu$f.biomass = matrix(data=NA,nrow=ntimes,ncol=nlu+1)
   lu$f.lai     = matrix(data=NA,nrow=ntimes,ncol=nlu+1)
   lu$f.gpp     = matrix(data=NA,nrow=ntimes,ncol=nlu+1)
   lu$f.npp     = matrix(data=NA,nrow=ntimes,ncol=nlu+1)
   lu$f.ba      = matrix(data=NA,nrow=ntimes,ncol=nlu+1)
   lu$dist      = array (data=NA,dim=c(ntimes,nlu,nlu))
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
   qmean$froot.resp     = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$het.resp       = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$cwd.resp       = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$soil.resp      = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$assim.light    = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$assim.rubp     = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$assim.co2      = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
   qmean$assim.ratio    = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
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
   qmsqu$froot.resp     = matrix(data=NA,nrow=ntimes,ncol=ndcycle)
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




   #---------------------------------------------------------------------------------------#
   #  PATCH -- patch level variables, we save as lists because the dimensions vary.    #
   #---------------------------------------------------------------------------------------#
   patch               = list()
   patch$ipa           = list()
   patch$age           = list()
   patch$area          = list()
   patch$lu            = list()
   patch$nep           = list()
   patch$het.resp      = list()
   patch$soil.resp     = list()
   patch$can.temp      = list()
   patch$gnd.temp      = list()
   patch$can.shv       = list()
   patch$gnd.shv       = list()
   patch$can.vpd       = list()
   patch$can.co2       = list()
   patch$can.prss      = list()
   patch$cflxca        = list()
   patch$cflxst        = list()
   patch$nee           = list()
   patch$hflxca        = list()
   patch$hflxgc        = list()
   patch$qwflxca       = list()
   patch$wflxca        = list()
   patch$wflxgc        = list()
   patch$ustar         = list()
   patch$albedo        = list()
   patch$rshortup      = list()
   patch$rlongup       = list()
   patch$parup         = list()
   patch$rnet          = list()
   patch$lai           = list()
   patch$wai           = list()
   patch$tai           = list()
   patch$agb           = list()
   patch$ba            = list()
   patch$wood.dens     = list()
   patch$can.depth     = list()
   patch$can.area      = list()
   patch$leaf.temp     = list()
   patch$leaf.water    = list()
   patch$leaf.vpd      = list()
   patch$wood.temp     = list()
   patch$gpp           = list()
   patch$npp           = list()
   patch$plant.resp    = list()
   patch$reco          = list()
   patch$hflxlc        = list()
   patch$hflxwc        = list()
   patch$wflxlc        = list()
   patch$wflxwc        = list()
   patch$transp        = list()
   patch$fast.soil.c   = list()
   patch$slow.soil.c   = list()
   patch$struct.soil.c = list()
   #---------------------------------------------------------------------------------------#




   #----- Cohort level, we save as lists because the dimensions vary. ---------------------#
   cohort                = list()
   cohort$ipa            = list()
   cohort$ico            = list()
   cohort$area           = list()
   cohort$lu             = list()
   cohort$dbh            = list()
   cohort$age            = list()
   cohort$pft            = list()
   cohort$nplant         = list()
   cohort$height         = list()
   cohort$ba             = list()
   cohort$agb            = list()
   cohort$biomass        = list()
   cohort$lai            = list()
   cohort$wai            = list()
   cohort$tai            = list()
   cohort$gpp            = list()
   cohort$leaf.resp      = list()
   cohort$root.resp      = list()
   cohort$froot.resp     = list()
   cohort$croot.resp     = list()
   cohort$stem.resp      = list()
   cohort$growth.resp    = list()
   cohort$storage.resp   = list()
   cohort$plant.resp     = list()
   cohort$assim.light    = list()
   cohort$assim.rubp     = list()
   cohort$assim.co2      = list()
   cohort$assim.ratio    = list()
   cohort$npp            = list()
   cohort$cba            = list()
   cohort$cbamax         = list()
   cohort$cbalight       = list()
   cohort$cbamoist       = list()
   cohort$cbarel         = list()
   cohort$dcbadt         = list()
   cohort$mcost          = list()
   cohort$ldrop          = list()
   cohort$sm.stress      = list()
   cohort$phap.sms       = list()
   cohort$light          = list()
   cohort$light.beam     = list()
   cohort$light.diff     = list()
   cohort$balive         = list()
   cohort$bdead          = list()
   cohort$bleaf          = list()
   cohort$bstem          = list()
   cohort$broot          = list()
   cohort$bfroot         = list()
   cohort$bcroot         = list()
   cohort$bsapwood       = list()
   cohort$bstorage       = list()
   cohort$bseeds         = list()
   cohort$hflxlc         = list()
   cohort$wflxlc         = list()
   cohort$transp         = list()
   cohort$wue            = list()
   cohort$cue            = list()
   cohort$ecue           = list()
   cohort$etue           = list()
   cohort$demand         = list()
   cohort$supply         = list()
   cohort$mort           = list()
   cohort$dimort         = list()
   cohort$ncbmort        = list()
   cohort$recruit        = list()
   cohort$growth         = list()
   cohort$agb.growth     = list()
   cohort$bsa.growth     = list()
   cohort$f.gpp          = list()
   cohort$f.plant.resp   = list()
   cohort$f.npp          = list()
   cohort$f.mco          = list()
   cohort$f.cba          = list()
   cohort$f.bstorage     = list()
   cohort$f.dcbadt       = list()
   cohort$f.bleaf        = list()
   cohort$f.bstem        = list()
   cohort$f.broot        = list()
   cohort$f.bseeds       = list()
   cohort$leaf.par       = list()
   cohort$phap.lpar      = list()
   cohort$leaf.rshort    = list()
   cohort$leaf.rlong     = list()
   cohort$rue            = list()
   #---------------------------------------------------------------------------------------#





   #----- Copy the polygon-level variable to the main structure. --------------------------#
   ed$emean  = emean
   ed$emsqu  = emsqu
   ed$qmean  = qmean
   ed$qmsqu  = qmsqu
   ed$lu     = lu
   ed$szpft  = szpft
   ed$patch  = patch
   ed$cohort = cohort
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   return(ed)
   #---------------------------------------------------------------------------------------#

}#end create.monthly
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     Expand the monthly array so it fits the new times.                                   #
#------------------------------------------------------------------------------------------#
update.monthly <<- function(new.ntimes,old.datum,montha,yeara,inpref,slz.min){

   #----- Create the new data set. --------------------------------------------------------#
   new.datum = create.monthly(new.ntimes,montha,yeara,inpref,slz.min)
   #---------------------------------------------------------------------------------------#


   #----- Find out which times to copy. ---------------------------------------------------#
   sel = old.datum$when %in% new.datum$when
   idx = match(old.datum$when[sel],new.datum$when)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   # emean -- variables that we can either compare directly with observations, or are      #
   #          or that may be used to draw time series.   They don't need to be really      #
   #          monthly means, but you should put only the variables that make sense to be   #
   #          plotted in simple time series (with no PFT or DBH information).              #
   #---------------------------------------------------------------------------------------#
   new.datum$emean$wood.dens         [idx ] = old.datum$emean$wood.dens           [sel ]
   new.datum$emean$fast.soil.c       [idx ] = old.datum$emean$fast.soil.c         [sel ]
   new.datum$emean$slow.soil.c       [idx ] = old.datum$emean$slow.soil.c         [sel ]
   new.datum$emean$struct.soil.c     [idx ] = old.datum$emean$struct.soil.c       [sel ]
   new.datum$emean$het.resp          [idx ] = old.datum$emean$het.resp            [sel ]
   new.datum$emean$cwd.resp          [idx ] = old.datum$emean$cwd.resp            [sel ]
   new.datum$emean$soil.resp         [idx ] = old.datum$emean$soil.resp           [sel ]
   new.datum$emean$gpp               [idx ] = old.datum$emean$gpp                 [sel ]
   new.datum$emean$last.1yr.gpp      [idx ] = old.datum$emean$last.1yr.gpp        [sel ]
   new.datum$emean$last.2yr.gpp      [idx ] = old.datum$emean$last.2yr.gpp        [sel ]
   new.datum$emean$last.3yr.gpp      [idx ] = old.datum$emean$last.3yr.gpp        [sel ]
   new.datum$emean$cue               [idx ] = old.datum$emean$cue                 [sel ]
   new.datum$emean$last.1yr.cue      [idx ] = old.datum$emean$last.1yr.cue        [sel ]
   new.datum$emean$last.2yr.cue      [idx ] = old.datum$emean$last.2yr.cue        [sel ]
   new.datum$emean$last.3yr.cue      [idx ] = old.datum$emean$last.3yr.cue        [sel ]
   new.datum$emean$ecue              [idx ] = old.datum$emean$ecue                [sel ]
   new.datum$emean$last.1yr.ecue     [idx ] = old.datum$emean$last.1yr.ecue       [sel ]
   new.datum$emean$last.2yr.ecue     [idx ] = old.datum$emean$last.2yr.ecue       [sel ]
   new.datum$emean$last.3yr.ecue     [idx ] = old.datum$emean$last.3yr.ecue       [sel ]
   new.datum$emean$npp               [idx ] = old.datum$emean$npp                 [sel ]
   new.datum$emean$last.1yr.npp      [idx ] = old.datum$emean$last.1yr.npp        [sel ]
   new.datum$emean$last.2yr.npp      [idx ] = old.datum$emean$last.2yr.npp        [sel ]
   new.datum$emean$last.3yr.npp      [idx ] = old.datum$emean$last.3yr.npp        [sel ]
   new.datum$emean$dcbadt            [idx ] = old.datum$emean$dcbadt              [sel ]
   new.datum$emean$last.1yr.dcbadt   [idx ] = old.datum$emean$last.1yr.dcbadt     [sel ]
   new.datum$emean$last.2yr.dcbadt   [idx ] = old.datum$emean$last.2yr.dcbadt     [sel ]
   new.datum$emean$last.3yr.dcbadt   [idx ] = old.datum$emean$last.3yr.dcbadt     [sel ]
   new.datum$emean$plant.resp        [idx ] = old.datum$emean$plant.resp          [sel ]
   new.datum$emean$last.1yr.plresp   [idx ] = old.datum$emean$last.1yr.plresp     [sel ]
   new.datum$emean$last.2yr.plresp   [idx ] = old.datum$emean$last.2yr.plresp     [sel ]
   new.datum$emean$last.3yr.plresp   [idx ] = old.datum$emean$last.3yr.plresp     [sel ]
   new.datum$emean$leaf.resp         [idx ] = old.datum$emean$leaf.resp           [sel ]
   new.datum$emean$root.resp         [idx ] = old.datum$emean$root.resp           [sel ]
   new.datum$emean$froot.resp        [idx ] = old.datum$emean$froot.resp          [sel ]
   new.datum$emean$croot.resp        [idx ] = old.datum$emean$croot.resp          [sel ]
   new.datum$emean$stem.resp         [idx ] = old.datum$emean$stem.resp           [sel ]
   new.datum$emean$growth.resp       [idx ] = old.datum$emean$growth.resp         [sel ]
   new.datum$emean$storage.resp      [idx ] = old.datum$emean$storage.resp        [sel ]
   new.datum$emean$reco              [idx ] = old.datum$emean$reco                [sel ]
   new.datum$emean$assim.light       [idx ] = old.datum$emean$assim.light         [sel ]
   new.datum$emean$assim.rubp        [idx ] = old.datum$emean$assim.rubp          [sel ]
   new.datum$emean$assim.co2         [idx ] = old.datum$emean$assim.co2           [sel ]
   new.datum$emean$assim.ratio       [idx ] = old.datum$emean$assim.ratio         [sel ]
   new.datum$emean$mco               [idx ] = old.datum$emean$mco                 [sel ]
   new.datum$emean$cba               [idx ] = old.datum$emean$cba                 [sel ]
   new.datum$emean$last.1yr.cba      [idx ] = old.datum$emean$last.1yr.cba        [sel ]
   new.datum$emean$last.2yr.cba      [idx ] = old.datum$emean$last.2yr.cba        [sel ]
   new.datum$emean$last.3yr.cba      [idx ] = old.datum$emean$last.3yr.cba        [sel ]
   new.datum$emean$cbamax            [idx ] = old.datum$emean$cbamax              [sel ]
   new.datum$emean$cbalight          [idx ] = old.datum$emean$cbalight            [sel ]
   new.datum$emean$cbamoist          [idx ] = old.datum$emean$cbamoist            [sel ]
   new.datum$emean$cbarel            [idx ] = old.datum$emean$cbarel              [sel ]
   new.datum$emean$ldrop             [idx ] = old.datum$emean$ldrop               [sel ]
   new.datum$emean$nep               [idx ] = old.datum$emean$nep                 [sel ]
   new.datum$emean$nee               [idx ] = old.datum$emean$nee                 [sel ]
   new.datum$emean$cflxca            [idx ] = old.datum$emean$cflxca              [sel ]
   new.datum$emean$cflxst            [idx ] = old.datum$emean$cflxst              [sel ]
   new.datum$emean$runoff            [idx ] = old.datum$emean$runoff              [sel ]
   new.datum$emean$intercepted       [idx ] = old.datum$emean$intercepted         [sel ]
   new.datum$emean$wshed             [idx ] = old.datum$emean$wshed               [sel ]
   new.datum$emean$evap              [idx ] = old.datum$emean$evap                [sel ]
   new.datum$emean$last.1yr.evap     [idx ] = old.datum$emean$last.1yr.evap       [sel ]
   new.datum$emean$last.2yr.evap     [idx ] = old.datum$emean$last.2yr.evap       [sel ]
   new.datum$emean$last.3yr.evap     [idx ] = old.datum$emean$last.3yr.evap       [sel ]
   new.datum$emean$ustar             [idx ] = old.datum$emean$ustar               [sel ]
   new.datum$emean$atm.vels          [idx ] = old.datum$emean$atm.vels            [sel ]
   new.datum$emean$atm.prss          [idx ] = old.datum$emean$atm.prss            [sel ]
   new.datum$emean$atm.temp          [idx ] = old.datum$emean$atm.temp            [sel ]
   new.datum$emean$can.prss          [idx ] = old.datum$emean$can.prss            [sel ]
   new.datum$emean$can.temp          [idx ] = old.datum$emean$can.temp            [sel ]
   new.datum$emean$atm.co2           [idx ] = old.datum$emean$atm.co2             [sel ]
   new.datum$emean$can.co2           [idx ] = old.datum$emean$can.co2             [sel ]
   new.datum$emean$can.depth         [idx ] = old.datum$emean$can.depth           [sel ]
   new.datum$emean$can.area          [idx ] = old.datum$emean$can.area            [sel ]
   new.datum$emean$leaf.temp         [idx ] = old.datum$emean$leaf.temp           [sel ]
   new.datum$emean$phap.ltemp        [idx ] = old.datum$emean$phap.ltemp          [sel ]
   new.datum$emean$last.1yr.ltemp    [idx ] = old.datum$emean$last.1yr.ltemp      [sel ]
   new.datum$emean$last.2yr.ltemp    [idx ] = old.datum$emean$last.2yr.ltemp      [sel ]
   new.datum$emean$last.3yr.ltemp    [idx ] = old.datum$emean$last.3yr.ltemp      [sel ]
   new.datum$emean$leaf.water        [idx ] = old.datum$emean$leaf.water          [sel ]
   new.datum$emean$phap.lwater       [idx ] = old.datum$emean$phap.lwater         [sel ]
   new.datum$emean$last.1yr.lwater   [idx ] = old.datum$emean$last.1yr.lwater     [sel ]
   new.datum$emean$last.2yr.lwater   [idx ] = old.datum$emean$last.2yr.lwater     [sel ]
   new.datum$emean$last.3yr.lwater   [idx ] = old.datum$emean$last.3yr.lwater     [sel ]
   new.datum$emean$wood.temp         [idx ] = old.datum$emean$wood.temp           [sel ]
   new.datum$emean$atm.shv           [idx ] = old.datum$emean$atm.shv             [sel ]
   new.datum$emean$can.shv           [idx ] = old.datum$emean$can.shv             [sel ]
   new.datum$emean$atm.vpd           [idx ] = old.datum$emean$atm.vpd             [sel ]
   new.datum$emean$can.vpd           [idx ] = old.datum$emean$can.vpd             [sel ]
   new.datum$emean$leaf.vpd          [idx ] = old.datum$emean$leaf.vpd            [sel ]
   new.datum$emean$phap.lvpd         [idx ] = old.datum$emean$phap.lvpd           [sel ]
   new.datum$emean$last.1yr.lvpd     [idx ] = old.datum$emean$last.1yr.lvpd       [sel ]
   new.datum$emean$last.2yr.lvpd     [idx ] = old.datum$emean$last.2yr.lvpd       [sel ]
   new.datum$emean$last.3yr.lvpd     [idx ] = old.datum$emean$last.3yr.lvpd       [sel ]
   new.datum$emean$can.co2           [idx ] = old.datum$emean$can.co2             [sel ]
   new.datum$emean$hflxca            [idx ] = old.datum$emean$hflxca              [sel ]
   new.datum$emean$qwflxca           [idx ] = old.datum$emean$qwflxca             [sel ]
   new.datum$emean$wflxca            [idx ] = old.datum$emean$wflxca              [sel ]
   new.datum$emean$agb               [idx ] = old.datum$emean$agb                 [sel ]
   new.datum$emean$biomass           [idx ] = old.datum$emean$biomass             [sel ]
   new.datum$emean$nplant            [idx ] = old.datum$emean$nplant              [sel ]
   new.datum$emean$lai               [idx ] = old.datum$emean$lai                 [sel ]
   new.datum$emean$wai               [idx ] = old.datum$emean$wai                 [sel ]
   new.datum$emean$tai               [idx ] = old.datum$emean$tai                 [sel ]
   new.datum$emean$area              [idx ] = old.datum$emean$area                [sel ]
   new.datum$emean$et                [idx ] = old.datum$emean$et                  [sel ]
   new.datum$emean$last.1yr.et       [idx ] = old.datum$emean$last.1yr.et         [sel ]
   new.datum$emean$last.2yr.et       [idx ] = old.datum$emean$last.2yr.et         [sel ]
   new.datum$emean$last.3yr.et       [idx ] = old.datum$emean$last.3yr.et         [sel ]
   new.datum$emean$transp            [idx ] = old.datum$emean$transp              [sel ]
   new.datum$emean$last.1yr.transp   [idx ] = old.datum$emean$last.1yr.transp     [sel ]
   new.datum$emean$last.2yr.transp   [idx ] = old.datum$emean$last.2yr.transp     [sel ]
   new.datum$emean$last.3yr.transp   [idx ] = old.datum$emean$last.3yr.transp     [sel ]
   new.datum$emean$wue               [idx ] = old.datum$emean$wue                 [sel ]
   new.datum$emean$last.1yr.wue      [idx ] = old.datum$emean$last.1yr.wue        [sel ]
   new.datum$emean$last.2yr.wue      [idx ] = old.datum$emean$last.2yr.wue        [sel ]
   new.datum$emean$last.3yr.wue      [idx ] = old.datum$emean$last.3yr.wue        [sel ]
   new.datum$emean$etue              [idx ] = old.datum$emean$etue                [sel ]
   new.datum$emean$last.1yr.etue     [idx ] = old.datum$emean$last.1yr.etue       [sel ]
   new.datum$emean$last.2yr.etue     [idx ] = old.datum$emean$last.2yr.etue       [sel ]
   new.datum$emean$last.3yr.etue     [idx ] = old.datum$emean$last.3yr.etue       [sel ]
   new.datum$emean$rain              [idx ] = old.datum$emean$rain                [sel ]
   new.datum$emean$last.1yr.rain     [idx ] = old.datum$emean$last.1yr.rain       [sel ]
   new.datum$emean$last.2yr.rain     [idx ] = old.datum$emean$last.2yr.rain       [sel ]
   new.datum$emean$last.3yr.rain     [idx ] = old.datum$emean$last.3yr.rain       [sel ]
   new.datum$emean$nmon.lt.090       [idx ] = old.datum$emean$nmon.lt.090         [sel ]
   new.datum$emean$nmon.lt.100       [idx ] = old.datum$emean$nmon.lt.100         [sel ]
   new.datum$emean$nmon.lt.110       [idx ] = old.datum$emean$nmon.lt.110         [sel ]
   new.datum$emean$nmon.lt.120       [idx ] = old.datum$emean$nmon.lt.120         [sel ]
   new.datum$emean$nmon.wdef         [idx ] = old.datum$emean$nmon.wdef           [sel ]
   new.datum$emean$nmon.mdef         [idx ] = old.datum$emean$nmon.mdef           [sel ]
   new.datum$emean$rue               [idx ] = old.datum$emean$rue                 [sel ]
   new.datum$emean$last.1yr.rue      [idx ] = old.datum$emean$last.1yr.rue        [sel ]
   new.datum$emean$last.2yr.rue      [idx ] = old.datum$emean$last.2yr.rue        [sel ]
   new.datum$emean$last.3yr.rue      [idx ] = old.datum$emean$last.3yr.rue        [sel ]
   new.datum$emean$gnd.temp          [idx ] = old.datum$emean$gnd.temp            [sel ]
   new.datum$emean$gnd.shv           [idx ] = old.datum$emean$gnd.shv             [sel ]
   new.datum$emean$workload          [idx ] = old.datum$emean$workload            [sel ]
   new.datum$emean$specwork          [idx ] = old.datum$emean$specwork            [sel ]
   new.datum$emean$sm.stress         [idx ] = old.datum$emean$sm.stress           [sel ]
   new.datum$emean$phap.sms          [idx ] = old.datum$emean$phap.sms            [sel ]
   new.datum$emean$last.1yr.sms      [idx ] = old.datum$emean$last.1yr.sms        [sel ]
   new.datum$emean$last.2yr.sms      [idx ] = old.datum$emean$last.2yr.sms        [sel ]
   new.datum$emean$last.3yr.sms      [idx ] = old.datum$emean$last.3yr.sms        [sel ]
   new.datum$emean$demand            [idx ] = old.datum$emean$demand              [sel ]
   new.datum$emean$supply            [idx ] = old.datum$emean$supply              [sel ]
   new.datum$emean$hflxgc            [idx ] = old.datum$emean$hflxgc              [sel ]
   new.datum$emean$hflxlc            [idx ] = old.datum$emean$hflxlc              [sel ]
   new.datum$emean$hflxwc            [idx ] = old.datum$emean$hflxwc              [sel ]
   new.datum$emean$wflxgc            [idx ] = old.datum$emean$wflxgc              [sel ]
   new.datum$emean$wflxlc            [idx ] = old.datum$emean$wflxlc              [sel ]
   new.datum$emean$wflxwc            [idx ] = old.datum$emean$wflxwc              [sel ]
   new.datum$emean$rshort            [idx ] = old.datum$emean$rshort              [sel ]
   new.datum$emean$last.1yr.rshort   [idx ] = old.datum$emean$last.1yr.rshort     [sel ]
   new.datum$emean$last.2yr.rshort   [idx ] = old.datum$emean$last.2yr.rshort     [sel ]
   new.datum$emean$last.3yr.rshort   [idx ] = old.datum$emean$last.3yr.rshort     [sel ]
   new.datum$emean$rshort.beam       [idx ] = old.datum$emean$rshort.beam         [sel ]
   new.datum$emean$rshort.diff       [idx ] = old.datum$emean$rshort.diff         [sel ]
   new.datum$emean$rshortup          [idx ] = old.datum$emean$rshortup            [sel ]
   new.datum$emean$rshort.gnd        [idx ] = old.datum$emean$rshort.gnd          [sel ]
   new.datum$emean$rlong             [idx ] = old.datum$emean$rlong               [sel ]
   new.datum$emean$rlong.gnd         [idx ] = old.datum$emean$rlong.gnd           [sel ]
   new.datum$emean$rlongup           [idx ] = old.datum$emean$rlongup             [sel ]
   new.datum$emean$par.tot           [idx ] = old.datum$emean$par.tot             [sel ]
   new.datum$emean$par.beam          [idx ] = old.datum$emean$par.beam            [sel ]
   new.datum$emean$par.diff          [idx ] = old.datum$emean$par.diff            [sel ]
   new.datum$emean$par.gnd           [idx ] = old.datum$emean$par.gnd             [sel ]
   new.datum$emean$parup             [idx ] = old.datum$emean$parup               [sel ]
   new.datum$emean$rnet              [idx ] = old.datum$emean$rnet                [sel ]
   new.datum$emean$albedo            [idx ] = old.datum$emean$albedo              [sel ]
   new.datum$emean$albedo.par        [idx ] = old.datum$emean$albedo.par          [sel ]
   new.datum$emean$albedo.nir        [idx ] = old.datum$emean$albedo.nir          [sel ]
   new.datum$emean$rlong.albedo      [idx ] = old.datum$emean$rlong.albedo        [sel ]
   new.datum$emean$paw               [idx ] = old.datum$emean$paw                 [sel ]
   new.datum$emean$smpot             [idx ] = old.datum$emean$smpot               [sel ]
   new.datum$emean$last.1yr.smpot    [idx ] = old.datum$emean$last.1yr.smpot      [sel ]
   new.datum$emean$last.2yr.smpot    [idx ] = old.datum$emean$last.2yr.smpot      [sel ]
   new.datum$emean$last.3yr.smpot    [idx ] = old.datum$emean$last.3yr.smpot      [sel ]
   new.datum$emean$npat.global       [idx ] = old.datum$emean$npat.global         [sel ]
   new.datum$emean$ncoh.global       [idx ] = old.datum$emean$ncoh.global         [sel ]
   new.datum$emean$water.deficit     [idx ] = old.datum$emean$water.deficit       [sel ]
   new.datum$emean$last.1yr.mwd      [idx ] = old.datum$emean$last.1yr.mwd        [sel ]
   new.datum$emean$last.2yr.mwd      [idx ] = old.datum$emean$last.2yr.mwd        [sel ]
   new.datum$emean$last.3yr.mwd      [idx ] = old.datum$emean$last.3yr.mwd        [sel ]
   new.datum$emean$malhi.deficit     [idx ] = old.datum$emean$malhi.deficit       [sel ]
   new.datum$emean$i.gpp             [idx ] = old.datum$emean$i.gpp               [sel ]
   new.datum$emean$i.npp             [idx ] = old.datum$emean$i.npp               [sel ]
   new.datum$emean$i.plant.resp      [idx ] = old.datum$emean$i.plant.resp        [sel ]
   new.datum$emean$i.mco             [idx ] = old.datum$emean$i.mco               [sel ]
   new.datum$emean$i.cba             [idx ] = old.datum$emean$i.cba               [sel ]
   new.datum$emean$i.cbamax          [idx ] = old.datum$emean$i.cbamax            [sel ]
   new.datum$emean$i.cbalight        [idx ] = old.datum$emean$i.cbalight          [sel ]
   new.datum$emean$i.cbamoist        [idx ] = old.datum$emean$i.cbamoist          [sel ]
   new.datum$emean$i.transp          [idx ] = old.datum$emean$i.transp            [sel ]
   new.datum$emean$i.wflxlc          [idx ] = old.datum$emean$i.wflxlc            [sel ]
   new.datum$emean$i.hflxlc          [idx ] = old.datum$emean$i.hflxlc            [sel ]
   new.datum$emean$f.gpp             [idx ] = old.datum$emean$f.gpp               [sel ]
   new.datum$emean$f.plant.resp      [idx ] = old.datum$emean$f.plant.resp        [sel ]
   new.datum$emean$f.npp             [idx ] = old.datum$emean$f.npp               [sel ]
   new.datum$emean$f.mco             [idx ] = old.datum$emean$f.mco               [sel ]
   new.datum$emean$f.cba             [idx ] = old.datum$emean$f.cba               [sel ]
   new.datum$emean$f.bstorage        [idx ] = old.datum$emean$f.bstorage          [sel ]
   new.datum$emean$f.bleaf           [idx ] = old.datum$emean$f.bleaf             [sel ]
   new.datum$emean$f.bstem           [idx ] = old.datum$emean$f.bstem             [sel ]
   new.datum$emean$f.broot           [idx ] = old.datum$emean$f.broot             [sel ]
   new.datum$emean$f.bseeds          [idx ] = old.datum$emean$f.bseeds            [sel ]
   new.datum$emean$f.dcbadt          [idx ] = old.datum$emean$f.dcbadt            [sel ]
   new.datum$emean$leaf.gsw          [idx ] = old.datum$emean$leaf.gsw            [sel ]
   new.datum$emean$phap.lgsw         [idx ] = old.datum$emean$phap.lgsw           [sel ]
   new.datum$emean$last.1yr.lgsw     [idx ] = old.datum$emean$last.1yr.lgsw       [sel ]
   new.datum$emean$last.2yr.lgsw     [idx ] = old.datum$emean$last.2yr.lgsw       [sel ]
   new.datum$emean$last.3yr.lgsw     [idx ] = old.datum$emean$last.3yr.lgsw       [sel ]
   new.datum$emean$leaf.gbw          [idx ] = old.datum$emean$leaf.gbw            [sel ]
   new.datum$emean$phap.lgbw         [idx ] = old.datum$emean$phap.lgbw           [sel ]
   new.datum$emean$wood.gbw          [idx ] = old.datum$emean$wood.gbw            [sel ]
   new.datum$emean$leaf.par          [idx ] = old.datum$emean$leaf.par            [sel ]
   new.datum$emean$phap.lpar         [idx ] = old.datum$emean$phap.lpar           [sel ]
   new.datum$emean$last.1yr.lpar     [idx ] = old.datum$emean$last.1yr.lpar       [sel ]
   new.datum$emean$last.2yr.lpar     [idx ] = old.datum$emean$last.2yr.lpar       [sel ]
   new.datum$emean$last.3yr.lpar     [idx ] = old.datum$emean$last.3yr.lpar       [sel ]
   new.datum$emean$leaf.rshort       [idx ] = old.datum$emean$leaf.rshort         [sel ]
   new.datum$emean$leaf.rlong        [idx ] = old.datum$emean$leaf.rlong          [sel ]
   new.datum$emean$soil.water        [idx,] = old.datum$emean$soil.water          [sel,]
   new.datum$emean$soil.temp         [idx,] = old.datum$emean$soil.temp           [sel,]
   new.datum$emean$soil.mstpot       [idx,] = old.datum$emean$soil.mstpot         [sel,]
   new.datum$emean$soil.extracted    [idx,] = old.datum$emean$soil.extracted      [sel,]
   new.datum$emean$last.1yr.growth   [idx ] = old.datum$emean$last.1yr.growth     [sel ]
   new.datum$emean$last.2yr.growth   [idx ] = old.datum$emean$last.2yr.growth     [sel ]
   new.datum$emean$last.3yr.growth   [idx ] = old.datum$emean$last.3yr.growth     [sel ]
   new.datum$emean$last.1yr.mort     [idx ] = old.datum$emean$last.1yr.mort       [sel ]
   new.datum$emean$last.2yr.mort     [idx ] = old.datum$emean$last.2yr.mort       [sel ]
   new.datum$emean$last.3yr.mort     [idx ] = old.datum$emean$last.3yr.mort       [sel ]
   new.datum$emean$last.1yr.dimort   [idx ] = old.datum$emean$last.1yr.dimort     [sel ]
   new.datum$emean$last.2yr.dimort   [idx ] = old.datum$emean$last.2yr.dimort     [sel ]
   new.datum$emean$last.3yr.dimort   [idx ] = old.datum$emean$last.3yr.dimort     [sel ]
   new.datum$emean$last.1yr.ncbmort  [idx ] = old.datum$emean$last.1yr.ncbmort    [sel ]
   new.datum$emean$last.2yr.ncbmort  [idx ] = old.datum$emean$last.2yr.ncbmort    [sel ]
   new.datum$emean$last.3yr.ncbmort  [idx ] = old.datum$emean$last.3yr.ncbmort    [sel ]
   new.datum$emean$agb.change        [idx ] = old.datum$emean$agb.change          [sel ]
   new.datum$emean$last.1yr.change   [idx ] = old.datum$emean$last.1yr.change     [sel ]
   new.datum$emean$last.2yr.change   [idx ] = old.datum$emean$last.2yr.change     [sel ]
   new.datum$emean$last.3yr.change   [idx ] = old.datum$emean$last.3yr.change     [sel ]
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   # emsqu -- mean sum of squares of polygon-level variable.                               #
   #---------------------------------------------------------------------------------------#
   new.datum$emsqu$gpp            [idx] = old.datum$emsqu$gpp            [sel]
   new.datum$emsqu$plant.resp     [idx] = old.datum$emsqu$plant.resp     [sel]
   new.datum$emsqu$het.resp       [idx] = old.datum$emsqu$het.resp       [sel]
   new.datum$emsqu$cwd.resp       [idx] = old.datum$emsqu$cwd.resp       [sel]
   new.datum$emsqu$soil.resp      [idx] = old.datum$emsqu$soil.resp      [sel]
   new.datum$emsqu$cflxca         [idx] = old.datum$emsqu$cflxca         [sel]
   new.datum$emsqu$cflxst         [idx] = old.datum$emsqu$cflxst         [sel]
   new.datum$emsqu$hflxca         [idx] = old.datum$emsqu$hflxca         [sel]
   new.datum$emsqu$hflxlc         [idx] = old.datum$emsqu$hflxlc         [sel]
   new.datum$emsqu$hflxwc         [idx] = old.datum$emsqu$hflxwc         [sel]
   new.datum$emsqu$hflxgc         [idx] = old.datum$emsqu$hflxgc         [sel]
   new.datum$emsqu$wflxca         [idx] = old.datum$emsqu$wflxca         [sel]
   new.datum$emsqu$qwflxca        [idx] = old.datum$emsqu$qwflxca        [sel]
   new.datum$emsqu$wflxlc         [idx] = old.datum$emsqu$wflxlc         [sel]
   new.datum$emsqu$wflxwc         [idx] = old.datum$emsqu$wflxwc         [sel]
   new.datum$emsqu$wflxgc         [idx] = old.datum$emsqu$wflxgc         [sel]
   new.datum$emsqu$transp         [idx] = old.datum$emsqu$transp         [sel]
   new.datum$emsqu$ustar          [idx] = old.datum$emsqu$ustar          [sel]
   new.datum$emsqu$albedo         [idx] = old.datum$emsqu$albedo         [sel]
   new.datum$emsqu$rshortup       [idx] = old.datum$emsqu$rshortup       [sel]
   new.datum$emsqu$rlongup        [idx] = old.datum$emsqu$rlongup        [sel]
   new.datum$emsqu$parup          [idx] = old.datum$emsqu$parup          [sel]
   new.datum$emsqu$rnet           [idx] = old.datum$emsqu$rnet           [sel]
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   # SZPFT -- Size (DBH) and plant functional type (PFT) array.  An extra level is         #
   #          appended to the end, which will hold the sum of all categories.              #
   #---------------------------------------------------------------------------------------#
   new.datum$szpft$agb            [idx,,] = old.datum$szpft$agb             [sel,,]
   new.datum$szpft$biomass        [idx,,] = old.datum$szpft$biomass         [sel,,]
   new.datum$szpft$lai            [idx,,] = old.datum$szpft$lai             [sel,,]
   new.datum$szpft$wai            [idx,,] = old.datum$szpft$wai             [sel,,]
   new.datum$szpft$tai            [idx,,] = old.datum$szpft$tai             [sel,,]
   new.datum$szpft$ba             [idx,,] = old.datum$szpft$ba              [sel,,]
   new.datum$szpft$gpp            [idx,,] = old.datum$szpft$gpp             [sel,,]
   new.datum$szpft$npp            [idx,,] = old.datum$szpft$npp             [sel,,]
   new.datum$szpft$dcbadt         [idx,,] = old.datum$szpft$dcbadt          [sel,,]
   new.datum$szpft$wood.dens      [idx,,] = old.datum$szpft$wood.dens       [sel,,]
   new.datum$szpft$leaf.resp      [idx,,] = old.datum$szpft$leaf.resp       [sel,,]
   new.datum$szpft$root.resp      [idx,,] = old.datum$szpft$root.resp       [sel,,]
   new.datum$szpft$froot.resp     [idx,,] = old.datum$szpft$froot.resp      [sel,,]
   new.datum$szpft$croot.resp     [idx,,] = old.datum$szpft$croot.resp      [sel,,]
   new.datum$szpft$stem.resp      [idx,,] = old.datum$szpft$stem.resp       [sel,,]
   new.datum$szpft$growth.resp    [idx,,] = old.datum$szpft$growth.resp     [sel,,]
   new.datum$szpft$storage.resp   [idx,,] = old.datum$szpft$storage.resp    [sel,,]
   new.datum$szpft$plant.resp     [idx,,] = old.datum$szpft$plant.resp      [sel,,]
   new.datum$szpft$assim.light    [idx,,] = old.datum$szpft$assim.light     [sel,,]
   new.datum$szpft$assim.rubp     [idx,,] = old.datum$szpft$assim.rubp      [sel,,]
   new.datum$szpft$assim.co2      [idx,,] = old.datum$szpft$assim.co2       [sel,,]
   new.datum$szpft$mco            [idx,,] = old.datum$szpft$mco             [sel,,]
   new.datum$szpft$cba            [idx,,] = old.datum$szpft$cba             [sel,,]
   new.datum$szpft$cbamax         [idx,,] = old.datum$szpft$cbamax          [sel,,]
   new.datum$szpft$cbalight       [idx,,] = old.datum$szpft$cbalight        [sel,,]
   new.datum$szpft$cbamoist       [idx,,] = old.datum$szpft$cbamoist        [sel,,]
   new.datum$szpft$cbarel         [idx,,] = old.datum$szpft$cbarel          [sel,,]
   new.datum$szpft$ldrop          [idx,,] = old.datum$szpft$ldrop           [sel,,]
   new.datum$szpft$sm.stress      [idx,,] = old.datum$szpft$sm.stress       [sel,,]
   new.datum$szpft$phap.sms       [idx,,] = old.datum$szpft$phap.sms        [sel,,]
   new.datum$szpft$leaf.gbw       [idx,,] = old.datum$szpft$leaf.gbw        [sel,,]
   new.datum$szpft$phap.lgbw      [idx,,] = old.datum$szpft$phap.lgbw       [sel,,]
   new.datum$szpft$leaf.gsw       [idx,,] = old.datum$szpft$leaf.gsw        [sel,,]
   new.datum$szpft$phap.lgsw      [idx,,] = old.datum$szpft$phap.lgsw       [sel,,]
   new.datum$szpft$wood.gbw       [idx,,] = old.datum$szpft$wood.gbw        [sel,,]
   new.datum$szpft$leaf.temp      [idx,,] = old.datum$szpft$leaf.temp       [sel,,]
   new.datum$szpft$phap.ltemp     [idx,,] = old.datum$szpft$phap.ltemp      [sel,,]
   new.datum$szpft$leaf.water     [idx,,] = old.datum$szpft$leaf.water      [sel,,]
   new.datum$szpft$phap.lwater    [idx,,] = old.datum$szpft$phap.lwater     [sel,,]
   new.datum$szpft$wood.temp      [idx,,] = old.datum$szpft$wood.temp       [sel,,]
   new.datum$szpft$leaf.vpd       [idx,,] = old.datum$szpft$leaf.vpd        [sel,,]
   new.datum$szpft$phap.lvpd      [idx,,] = old.datum$szpft$phap.lvpd       [sel,,]
   new.datum$szpft$demand         [idx,,] = old.datum$szpft$demand          [sel,,]
   new.datum$szpft$supply         [idx,,] = old.datum$szpft$supply          [sel,,]
   new.datum$szpft$nplant         [idx,,] = old.datum$szpft$nplant          [sel,,]
   new.datum$szpft$mort           [idx,,] = old.datum$szpft$mort            [sel,,]
   new.datum$szpft$dimort         [idx,,] = old.datum$szpft$dimort          [sel,,]
   new.datum$szpft$ncbmort        [idx,,] = old.datum$szpft$ncbmort         [sel,,]
   new.datum$szpft$growth         [idx,,] = old.datum$szpft$growth          [sel,,]
   new.datum$szpft$recr           [idx,,] = old.datum$szpft$recr            [sel,,]
   new.datum$szpft$change         [idx,,] = old.datum$szpft$change          [sel,,]
   new.datum$szpft$agb.mort       [idx,,] = old.datum$szpft$agb.mort        [sel,,]
   new.datum$szpft$agb.dimort     [idx,,] = old.datum$szpft$agb.dimort      [sel,,]
   new.datum$szpft$agb.ncbmort    [idx,,] = old.datum$szpft$agb.ncbmort     [sel,,]
   new.datum$szpft$agb.growth     [idx,,] = old.datum$szpft$agb.growth      [sel,,]
   new.datum$szpft$agb.recr       [idx,,] = old.datum$szpft$agb.recr        [sel,,]
   new.datum$szpft$agb.change     [idx,,] = old.datum$szpft$agb.change      [sel,,]
   new.datum$szpft$bsa.mort       [idx,,] = old.datum$szpft$bsa.mort        [sel,,]
   new.datum$szpft$bsa.dimort     [idx,,] = old.datum$szpft$bsa.dimort      [sel,,]
   new.datum$szpft$bsa.ncbmort    [idx,,] = old.datum$szpft$bsa.ncbmort     [sel,,]
   new.datum$szpft$bsa.growth     [idx,,] = old.datum$szpft$bsa.growth      [sel,,]
   new.datum$szpft$bsa.recr       [idx,,] = old.datum$szpft$bsa.recr        [sel,,]
   new.datum$szpft$bsa.change     [idx,,] = old.datum$szpft$bsa.change      [sel,,]
   new.datum$szpft$bdead          [idx,,] = old.datum$szpft$bdead           [sel,,]
   new.datum$szpft$balive         [idx,,] = old.datum$szpft$balive          [sel,,]
   new.datum$szpft$bleaf          [idx,,] = old.datum$szpft$bleaf           [sel,,]
   new.datum$szpft$bstem          [idx,,] = old.datum$szpft$bstem           [sel,,]
   new.datum$szpft$broot          [idx,,] = old.datum$szpft$broot           [sel,,]
   new.datum$szpft$bfroot         [idx,,] = old.datum$szpft$bfroot          [sel,,]
   new.datum$szpft$bcroot         [idx,,] = old.datum$szpft$bcroot          [sel,,]
   new.datum$szpft$bsapwood       [idx,,] = old.datum$szpft$bsapwood        [sel,,]
   new.datum$szpft$bstorage       [idx,,] = old.datum$szpft$bstorage        [sel,,]
   new.datum$szpft$bseeds         [idx,,] = old.datum$szpft$bseeds          [sel,,]
   new.datum$szpft$hflxlc         [idx,,] = old.datum$szpft$hflxlc          [sel,,]
   new.datum$szpft$wflxlc         [idx,,] = old.datum$szpft$wflxlc          [sel,,]
   new.datum$szpft$census.lai     [idx,,] = old.datum$szpft$census.lai      [sel,,]
   new.datum$szpft$census.wai     [idx,,] = old.datum$szpft$census.wai      [sel,,]
   new.datum$szpft$census.tai     [idx,,] = old.datum$szpft$census.tai      [sel,,]
   new.datum$szpft$census.agb     [idx,,] = old.datum$szpft$census.agb      [sel,,]
   new.datum$szpft$census.ba      [idx,,] = old.datum$szpft$census.ba       [sel,,]
   new.datum$szpft$i.gpp          [idx,,] = old.datum$szpft$i.gpp           [sel,,]
   new.datum$szpft$i.npp          [idx,,] = old.datum$szpft$i.npp           [sel,,]
   new.datum$szpft$i.plant.resp   [idx,,] = old.datum$szpft$i.plant.resp    [sel,,]
   new.datum$szpft$i.mco          [idx,,] = old.datum$szpft$i.mco           [sel,,]
   new.datum$szpft$i.cba          [idx,,] = old.datum$szpft$i.cba           [sel,,]
   new.datum$szpft$i.cbamax       [idx,,] = old.datum$szpft$i.cbamax        [sel,,]
   new.datum$szpft$i.cbalight     [idx,,] = old.datum$szpft$i.cbalight      [sel,,]
   new.datum$szpft$i.cbamoist     [idx,,] = old.datum$szpft$i.cbamoist      [sel,,]
   new.datum$szpft$i.transp       [idx,,] = old.datum$szpft$i.transp        [sel,,]
   new.datum$szpft$i.wflxlc       [idx,,] = old.datum$szpft$i.wflxlc        [sel,,]
   new.datum$szpft$i.hflxlc       [idx,,] = old.datum$szpft$i.hflxlc        [sel,,]
   new.datum$szpft$f.gpp          [idx,,] = old.datum$szpft$f.gpp           [sel,,]
   new.datum$szpft$f.plant.resp   [idx,,] = old.datum$szpft$f.plant.resp    [sel,,]
   new.datum$szpft$f.npp          [idx,,] = old.datum$szpft$f.npp           [sel,,]
   new.datum$szpft$f.mco          [idx,,] = old.datum$szpft$f.mco           [sel,,]
   new.datum$szpft$f.cba          [idx,,] = old.datum$szpft$f.cba           [sel,,]
   new.datum$szpft$f.bstorage     [idx,,] = old.datum$szpft$f.bstorage      [sel,,]
   new.datum$szpft$f.bleaf        [idx,,] = old.datum$szpft$f.bleaf         [sel,,]
   new.datum$szpft$f.bstem        [idx,,] = old.datum$szpft$f.bstem         [sel,,]
   new.datum$szpft$f.broot        [idx,,] = old.datum$szpft$f.broot         [sel,,]
   new.datum$szpft$f.bseeds       [idx,,] = old.datum$szpft$f.bseeds        [sel,,]
   new.datum$szpft$f.dcbadt       [idx,,] = old.datum$szpft$f.dcbadt        [sel,,]
   new.datum$szpft$leaf.par       [idx,,] = old.datum$szpft$leaf.par        [sel,,]
   new.datum$szpft$phap.lpar      [idx,,] = old.datum$szpft$phap.lpar       [sel,,]
   new.datum$szpft$leaf.rshort    [idx,,] = old.datum$szpft$leaf.rshort     [sel,,]
   new.datum$szpft$leaf.rlong     [idx,,] = old.datum$szpft$leaf.rlong      [sel,,]
   new.datum$szpft$transp         [idx,,] = old.datum$szpft$transp          [sel,,]
   new.datum$szpft$wue            [idx,,] = old.datum$szpft$wue             [sel,,]
   new.datum$szpft$cue            [idx,,] = old.datum$szpft$cue             [sel,,]
   new.datum$szpft$ecue           [idx,,] = old.datum$szpft$ecue            [sel,,]
   new.datum$szpft$etue           [idx,,] = old.datum$szpft$etue            [sel,,]
   new.datum$szpft$rue            [idx,,] = old.datum$szpft$rue             [sel,,]
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   # LU -- Polygon-level variables, split by land use type.  One extra dimension is        #
   #       appended to the end, which will hold the sum of all land use types.             #
   #---------------------------------------------------------------------------------------#
   new.datum$lu$agb               [idx, ] = old.datum$lu$agb               [sel, ]
   new.datum$lu$biomass           [idx, ] = old.datum$lu$biomass           [sel, ]
   new.datum$lu$lai               [idx, ] = old.datum$lu$lai               [sel, ]
   new.datum$lu$gpp               [idx, ] = old.datum$lu$gpp               [sel, ]
   new.datum$lu$npp               [idx, ] = old.datum$lu$npp               [sel, ]
   new.datum$lu$area              [idx, ] = old.datum$lu$area              [sel, ]
   new.datum$lu$ba                [idx, ] = old.datum$lu$ba                [sel, ]
   new.datum$lu$dist              [idx,,] = old.datum$lu$dist              [sel,,]
   new.datum$lu$f.agb             [idx, ] = old.datum$lu$f.agb             [sel, ]
   new.datum$lu$f.biomass         [idx, ] = old.datum$lu$f.biomass         [sel, ]
   new.datum$lu$f.lai             [idx, ] = old.datum$lu$f.lai             [sel, ]
   new.datum$lu$f.gpp             [idx, ] = old.datum$lu$f.gpp             [sel, ]
   new.datum$lu$f.npp             [idx, ] = old.datum$lu$f.npp             [sel, ]
   new.datum$lu$f.ba              [idx, ] = old.datum$lu$f.ba              [sel, ]
   #---------------------------------------------------------------------------------------#









   #---------------------------------------------------------------------------------------#
   # QMEAN -- Polygon-level variables, containing the mean diel (diurnal cycle).           #
   #---------------------------------------------------------------------------------------#
   new.datum$qmean$gpp           [idx,] = old.datum$qmean$gpp            [sel,]
   new.datum$qmean$npp           [idx,] = old.datum$qmean$npp            [sel,]
   new.datum$qmean$plant.resp    [idx,] = old.datum$qmean$plant.resp     [sel,]
   new.datum$qmean$leaf.resp     [idx,] = old.datum$qmean$leaf.resp      [sel,]
   new.datum$qmean$root.resp     [idx,] = old.datum$qmean$root.resp      [sel,]
   new.datum$qmean$froot.resp    [idx,] = old.datum$qmean$froot.resp     [sel,]
   new.datum$qmean$het.resp      [idx,] = old.datum$qmean$het.resp       [sel,]
   new.datum$qmean$cwd.resp      [idx,] = old.datum$qmean$cwd.resp       [sel,]
   new.datum$qmean$soil.resp     [idx,] = old.datum$qmean$soil.resp      [sel,]
   new.datum$qmean$assim.light   [idx,] = old.datum$qmean$assim.light    [sel,]
   new.datum$qmean$assim.rubp    [idx,] = old.datum$qmean$assim.rubp     [sel,]
   new.datum$qmean$assim.co2     [idx,] = old.datum$qmean$assim.co2      [sel,]
   new.datum$qmean$assim.ratio   [idx,] = old.datum$qmean$assim.ratio    [sel,]
   new.datum$qmean$nep           [idx,] = old.datum$qmean$nep            [sel,]
   new.datum$qmean$nee           [idx,] = old.datum$qmean$nee            [sel,]
   new.datum$qmean$reco          [idx,] = old.datum$qmean$reco           [sel,]
   new.datum$qmean$cflxca        [idx,] = old.datum$qmean$cflxca         [sel,]
   new.datum$qmean$cflxst        [idx,] = old.datum$qmean$cflxst         [sel,]
   new.datum$qmean$hflxca        [idx,] = old.datum$qmean$hflxca         [sel,]
   new.datum$qmean$hflxlc        [idx,] = old.datum$qmean$hflxlc         [sel,]
   new.datum$qmean$hflxwc        [idx,] = old.datum$qmean$hflxwc         [sel,]
   new.datum$qmean$hflxgc        [idx,] = old.datum$qmean$hflxgc         [sel,]
   new.datum$qmean$qwflxca       [idx,] = old.datum$qmean$qwflxca        [sel,]
   new.datum$qmean$wflxca        [idx,] = old.datum$qmean$wflxca         [sel,]
   new.datum$qmean$wflxlc        [idx,] = old.datum$qmean$wflxlc         [sel,]
   new.datum$qmean$wflxwc        [idx,] = old.datum$qmean$wflxwc         [sel,]
   new.datum$qmean$wflxgc        [idx,] = old.datum$qmean$wflxgc         [sel,]
   new.datum$qmean$runoff        [idx,] = old.datum$qmean$runoff         [sel,]
   new.datum$qmean$intercepted   [idx,] = old.datum$qmean$intercepted    [sel,]
   new.datum$qmean$wshed         [idx,] = old.datum$qmean$wshed          [sel,]
   new.datum$qmean$evap          [idx,] = old.datum$qmean$evap           [sel,]
   new.datum$qmean$transp        [idx,] = old.datum$qmean$transp         [sel,]
   new.datum$qmean$atm.temp      [idx,] = old.datum$qmean$atm.temp       [sel,]
   new.datum$qmean$can.temp      [idx,] = old.datum$qmean$can.temp       [sel,]
   new.datum$qmean$leaf.temp     [idx,] = old.datum$qmean$leaf.temp      [sel,]
   new.datum$qmean$leaf.water    [idx,] = old.datum$qmean$leaf.water     [sel,]
   new.datum$qmean$wood.temp     [idx,] = old.datum$qmean$wood.temp      [sel,]
   new.datum$qmean$gnd.temp      [idx,] = old.datum$qmean$gnd.temp       [sel,]
   new.datum$qmean$atm.shv       [idx,] = old.datum$qmean$atm.shv        [sel,]
   new.datum$qmean$can.shv       [idx,] = old.datum$qmean$can.shv        [sel,]
   new.datum$qmean$gnd.shv       [idx,] = old.datum$qmean$gnd.shv        [sel,]
   new.datum$qmean$atm.vpd       [idx,] = old.datum$qmean$atm.vpd        [sel,]
   new.datum$qmean$can.vpd       [idx,] = old.datum$qmean$can.vpd        [sel,]
   new.datum$qmean$leaf.vpd      [idx,] = old.datum$qmean$leaf.vpd       [sel,]
   new.datum$qmean$atm.co2       [idx,] = old.datum$qmean$atm.co2        [sel,]
   new.datum$qmean$can.co2       [idx,] = old.datum$qmean$can.co2        [sel,]
   new.datum$qmean$atm.prss      [idx,] = old.datum$qmean$atm.prss       [sel,]
   new.datum$qmean$can.prss      [idx,] = old.datum$qmean$can.prss       [sel,]
   new.datum$qmean$atm.vels      [idx,] = old.datum$qmean$atm.vels       [sel,]
   new.datum$qmean$ustar         [idx,] = old.datum$qmean$ustar          [sel,]
   new.datum$qmean$sm.stress     [idx,] = old.datum$qmean$sm.stress      [sel,]
   new.datum$qmean$rain          [idx,] = old.datum$qmean$rain           [sel,]
   new.datum$qmean$rshort        [idx,] = old.datum$qmean$rshort         [sel,]
   new.datum$qmean$rshort.beam   [idx,] = old.datum$qmean$rshort.beam    [sel,]
   new.datum$qmean$rshort.diff   [idx,] = old.datum$qmean$rshort.diff    [sel,]
   new.datum$qmean$rshort.gnd    [idx,] = old.datum$qmean$rshort.gnd     [sel,]
   new.datum$qmean$rshortup      [idx,] = old.datum$qmean$rshortup       [sel,]
   new.datum$qmean$rlong         [idx,] = old.datum$qmean$rlong          [sel,]
   new.datum$qmean$rlong.gnd     [idx,] = old.datum$qmean$rlong.gnd      [sel,]
   new.datum$qmean$rlongup       [idx,] = old.datum$qmean$rlongup        [sel,]
   new.datum$qmean$par.tot       [idx,] = old.datum$qmean$par.tot        [sel,]
   new.datum$qmean$par.beam      [idx,] = old.datum$qmean$par.beam       [sel,]
   new.datum$qmean$par.diff      [idx,] = old.datum$qmean$par.diff       [sel,]
   new.datum$qmean$par.gnd       [idx,] = old.datum$qmean$par.gnd        [sel,]
   new.datum$qmean$parup         [idx,] = old.datum$qmean$parup          [sel,]
   new.datum$qmean$rnet          [idx,] = old.datum$qmean$rnet           [sel,]
   new.datum$qmean$albedo        [idx,] = old.datum$qmean$albedo         [sel,]
   new.datum$qmean$albedo.par    [idx,] = old.datum$qmean$albedo.par     [sel,]
   new.datum$qmean$albedo.nir    [idx,] = old.datum$qmean$albedo.nir     [sel,]
   new.datum$qmean$rlong.albedo  [idx,] = old.datum$qmean$rlong.albedo   [sel,]
   new.datum$qmean$leaf.gsw      [idx,] = old.datum$qmean$leaf.gsw       [sel,]
   new.datum$qmean$leaf.gbw      [idx,] = old.datum$qmean$leaf.gbw       [sel,]
   new.datum$qmean$wood.gbw      [idx,] = old.datum$qmean$wood.gbw       [sel,]
   #---------------------------------------------------------------------------------------#









   #---------------------------------------------------------------------------------------#
   # QMSQU -- Polygon-level variables, containing the mean sum of squares for the diel     #
   #          (diurnal cycle).                                                             #
   #---------------------------------------------------------------------------------------#
   new.datum$qmsqu$gpp           [idx,] = old.datum$qmsqu$gpp           [sel,] 
   new.datum$qmsqu$npp           [idx,] = old.datum$qmsqu$npp           [sel,] 
   new.datum$qmsqu$plant.resp    [idx,] = old.datum$qmsqu$plant.resp    [sel,] 
   new.datum$qmsqu$leaf.resp     [idx,] = old.datum$qmsqu$leaf.resp     [sel,] 
   new.datum$qmsqu$root.resp     [idx,] = old.datum$qmsqu$root.resp     [sel,] 
   new.datum$qmsqu$froot.resp    [idx,] = old.datum$qmsqu$froot.resp    [sel,] 
   new.datum$qmsqu$het.resp      [idx,] = old.datum$qmsqu$het.resp      [sel,] 
   new.datum$qmsqu$cwd.resp      [idx,] = old.datum$qmsqu$cwd.resp      [sel,] 
   new.datum$qmsqu$soil.resp     [idx,] = old.datum$qmsqu$soil.resp     [sel,] 
   new.datum$qmsqu$nep           [idx,] = old.datum$qmsqu$nep           [sel,] 
   new.datum$qmsqu$cflxca        [idx,] = old.datum$qmsqu$cflxca        [sel,] 
   new.datum$qmsqu$cflxst        [idx,] = old.datum$qmsqu$cflxst        [sel,] 
   new.datum$qmsqu$hflxca        [idx,] = old.datum$qmsqu$hflxca        [sel,] 
   new.datum$qmsqu$hflxlc        [idx,] = old.datum$qmsqu$hflxlc        [sel,] 
   new.datum$qmsqu$hflxwc        [idx,] = old.datum$qmsqu$hflxwc        [sel,] 
   new.datum$qmsqu$hflxgc        [idx,] = old.datum$qmsqu$hflxgc        [sel,] 
   new.datum$qmsqu$qwflxca       [idx,] = old.datum$qmsqu$qwflxca       [sel,] 
   new.datum$qmsqu$wflxca        [idx,] = old.datum$qmsqu$wflxca        [sel,] 
   new.datum$qmsqu$wflxlc        [idx,] = old.datum$qmsqu$wflxlc        [sel,] 
   new.datum$qmsqu$wflxwc        [idx,] = old.datum$qmsqu$wflxwc        [sel,] 
   new.datum$qmsqu$wflxgc        [idx,] = old.datum$qmsqu$wflxgc        [sel,] 
   new.datum$qmsqu$transp        [idx,] = old.datum$qmsqu$transp        [sel,] 
   new.datum$qmsqu$ustar         [idx,] = old.datum$qmsqu$ustar         [sel,] 
   new.datum$qmsqu$albedo        [idx,] = old.datum$qmsqu$albedo        [sel,] 
   new.datum$qmsqu$rshortup      [idx,] = old.datum$qmsqu$rshortup      [sel,] 
   new.datum$qmsqu$rlongup       [idx,] = old.datum$qmsqu$rlongup       [sel,] 
   new.datum$qmsqu$parup         [idx,] = old.datum$qmsqu$parup         [sel,] 
   new.datum$qmsqu$rnet          [idx,] = old.datum$qmsqu$rnet          [sel,] 
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #  PATCH -- patch level variables, we save as lists because the dimensions vary.    #
   #---------------------------------------------------------------------------------------#
   new.datum$patch$ipa           = old.datum$patch$ipa
   new.datum$patch$age           = old.datum$patch$age
   new.datum$patch$area          = old.datum$patch$area
   new.datum$patch$lu            = old.datum$patch$lu
   new.datum$patch$nep           = old.datum$patch$nep
   new.datum$patch$het.resp      = old.datum$patch$het.resp
   new.datum$patch$soil.resp     = old.datum$patch$soil.resp
   new.datum$patch$can.temp      = old.datum$patch$can.temp
   new.datum$patch$gnd.temp      = old.datum$patch$gnd.temp
   new.datum$patch$can.shv       = old.datum$patch$can.shv
   new.datum$patch$gnd.shv       = old.datum$patch$gnd.shv
   new.datum$patch$can.vpd       = old.datum$patch$can.vpd
   new.datum$patch$can.co2       = old.datum$patch$can.co2
   new.datum$patch$can.prss      = old.datum$patch$can.prss
   new.datum$patch$cflxca        = old.datum$patch$cflxca
   new.datum$patch$cflxst        = old.datum$patch$cflxst
   new.datum$patch$nee           = old.datum$patch$nee
   new.datum$patch$hflxca        = old.datum$patch$hflxca
   new.datum$patch$hflxgc        = old.datum$patch$hflxgc
   new.datum$patch$qwflxca       = old.datum$patch$qwflxca
   new.datum$patch$wflxca        = old.datum$patch$wflxca
   new.datum$patch$wflxgc        = old.datum$patch$wflxgc
   new.datum$patch$ustar         = old.datum$patch$ustar
   new.datum$patch$albedo        = old.datum$patch$albedo
   new.datum$patch$rshortup      = old.datum$patch$rshortup
   new.datum$patch$rlongup       = old.datum$patch$rlongup
   new.datum$patch$parup         = old.datum$patch$parup
   new.datum$patch$rnet          = old.datum$patch$rnet
   new.datum$patch$lai           = old.datum$patch$lai
   new.datum$patch$wai           = old.datum$patch$wai
   new.datum$patch$tai           = old.datum$patch$tai
   new.datum$patch$agb           = old.datum$patch$agb
   new.datum$patch$ba            = old.datum$patch$ba
   new.datum$patch$wood.dens     = old.datum$patch$wood.dens
   new.datum$patch$can.depth     = old.datum$patch$can.depth
   new.datum$patch$can.area      = old.datum$patch$can.area
   new.datum$patch$leaf.temp     = old.datum$patch$leaf.temp
   new.datum$patch$leaf.water    = old.datum$patch$leaf.water
   new.datum$patch$leaf.vpd      = old.datum$patch$leaf.vpd
   new.datum$patch$wood.temp     = old.datum$patch$wood.temp
   new.datum$patch$gpp           = old.datum$patch$gpp
   new.datum$patch$npp           = old.datum$patch$npp
   new.datum$patch$plant.resp    = old.datum$patch$plant.resp
   new.datum$patch$reco          = old.datum$patch$reco
   new.datum$patch$hflxlc        = old.datum$patch$hflxlc
   new.datum$patch$hflxwc        = old.datum$patch$hflxwc
   new.datum$patch$wflxlc        = old.datum$patch$wflxlc
   new.datum$patch$wflxwc        = old.datum$patch$wflxwc
   new.datum$patch$transp        = old.datum$patch$transp
   new.datum$patch$fast.soil.c   = old.datum$patch$fast.soil.c
   new.datum$patch$slow.soil.c   = old.datum$patch$slow.soil.c
   new.datum$patch$struct.soil.c = old.datum$patch$struct.soil.c
   #---------------------------------------------------------------------------------------#




   #----- Cohort level, we save as lists because the dimensions vary. ---------------------#
   new.datum$cohort$ipa              = old.datum$cohort$ipa
   new.datum$cohort$ico              = old.datum$cohort$ico
   new.datum$cohort$area             = old.datum$cohort$area
   new.datum$cohort$lu               = old.datum$cohort$lu
   new.datum$cohort$dbh              = old.datum$cohort$dbh
   new.datum$cohort$age              = old.datum$cohort$age
   new.datum$cohort$pft              = old.datum$cohort$pft
   new.datum$cohort$nplant           = old.datum$cohort$nplant
   new.datum$cohort$height           = old.datum$cohort$height
   new.datum$cohort$ba               = old.datum$cohort$ba
   new.datum$cohort$agb              = old.datum$cohort$agb
   new.datum$cohort$biomass          = old.datum$cohort$biomass
   new.datum$cohort$lai              = old.datum$cohort$lai
   new.datum$cohort$wai              = old.datum$cohort$wai
   new.datum$cohort$tai              = old.datum$cohort$tai
   new.datum$cohort$gpp              = old.datum$cohort$gpp
   new.datum$cohort$leaf.resp        = old.datum$cohort$leaf.resp
   new.datum$cohort$root.resp        = old.datum$cohort$root.resp
   new.datum$cohort$froot.resp       = old.datum$cohort$froot.resp
   new.datum$cohort$croot.resp       = old.datum$cohort$croot.resp
   new.datum$cohort$stem.resp        = old.datum$cohort$stem.resp
   new.datum$cohort$storage.resp     = old.datum$cohort$storage.resp
   new.datum$cohort$plant.resp       = old.datum$cohort$plant.resp
   new.datum$cohort$assim.light      = old.datum$cohort$assim.light
   new.datum$cohort$assim.rubp       = old.datum$cohort$assim.rubp
   new.datum$cohort$assim.co2        = old.datum$cohort$assim.co2
   new.datum$cohort$assim.ratio      = old.datum$cohort$assim.ratio
   new.datum$cohort$npp              = old.datum$cohort$npp
   new.datum$cohort$cba              = old.datum$cohort$cba
   new.datum$cohort$cbamax           = old.datum$cohort$cbamax
   new.datum$cohort$cbalight         = old.datum$cohort$cbalight
   new.datum$cohort$cbamoist         = old.datum$cohort$cbamoist
   new.datum$cohort$cbarel           = old.datum$cohort$cbarel
   new.datum$cohort$dcbadt           = old.datum$cohort$dcbadt 
   new.datum$cohort$mcost            = old.datum$cohort$mcost
   new.datum$cohort$ldrop            = old.datum$cohort$ldrop
   new.datum$cohort$sm.stress        = old.datum$cohort$sm.stress
   new.datum$cohort$phap.sms         = old.datum$cohort$phap.sms
   new.datum$cohort$light            = old.datum$cohort$light
   new.datum$cohort$lightbeam        = old.datum$cohort$lightbeam
   new.datum$cohort$lightdiff        = old.datum$cohort$lightdiff
   new.datum$cohort$balive           = old.datum$cohort$balive
   new.datum$cohort$bdead            = old.datum$cohort$bdead
   new.datum$cohort$bleaf            = old.datum$cohort$bleaf
   new.datum$cohort$bstem            = old.datum$cohort$bstem
   new.datum$cohort$broot            = old.datum$cohort$broot
   new.datum$cohort$bfroot           = old.datum$cohort$bfroot
   new.datum$cohort$bcroot           = old.datum$cohort$bcroot
   new.datum$cohort$bsapwood         = old.datum$cohort$bsapwood
   new.datum$cohort$bstorage         = old.datum$cohort$bstorage
   new.datum$cohort$bseeds           = old.datum$cohort$bseeds
   new.datum$cohort$hflxlc           = old.datum$cohort$hflxlc
   new.datum$cohort$wflxlc           = old.datum$cohort$wflxlc
   new.datum$cohort$transp           = old.datum$cohort$transp
   new.datum$cohort$wue              = old.datum$cohort$wue
   new.datum$cohort$cue              = old.datum$cohort$cue
   new.datum$cohort$ecue             = old.datum$cohort$ecue
   new.datum$cohort$etue             = old.datum$cohort$etue
   new.datum$cohort$demand           = old.datum$cohort$demand
   new.datum$cohort$supply           = old.datum$cohort$supply
   new.datum$cohort$mort             = old.datum$cohort$mort
   new.datum$cohort$dimort           = old.datum$cohort$dimort
   new.datum$cohort$ncbmort          = old.datum$cohort$ncbmort
   new.datum$cohort$recruit          = old.datum$cohort$recruit
   new.datum$cohort$growth           = old.datum$cohort$growth
   new.datum$cohort$agb.growth       = old.datum$cohort$agb.growth
   new.datum$cohort$bsa.growth       = old.datum$cohort$bsa.growth
   new.datum$cohort$f.gpp            = old.datum$cohort$f.gpp
   new.datum$cohort$f.plant.resp     = old.datum$cohort$f.plant.resp
   new.datum$cohort$f.npp            = old.datum$cohort$f.npp
   new.datum$cohort$f.mco            = old.datum$cohort$f.mco
   new.datum$cohort$f.cba            = old.datum$cohort$f.cba
   new.datum$cohort$f.bstorage       = old.datum$cohort$f.bstorage
   new.datum$cohort$f.bleaf          = old.datum$cohort$f.bleaf
   new.datum$cohort$f.bstem          = old.datum$cohort$f.bstem
   new.datum$cohort$f.broot          = old.datum$cohort$f.broot
   new.datum$cohort$f.bseeds         = old.datum$cohort$f.bseeds
   new.datum$cohort$f.dcbadt         = old.datum$cohort$f.dcbadt 
   new.datum$cohort$leaf.par         = old.datum$cohort$leaf.par
   new.datum$cohort$phap.lpar        = old.datum$cohort$phap.lpar
   new.datum$cohort$leaf.rshort      = old.datum$cohort$leaf.rshort
   new.datum$cohort$leaf.rlong       = old.datum$cohort$leaf.rlong
   new.datum$cohort$rue              = old.datum$cohort$rue
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Send the data back.                                                               #
   #---------------------------------------------------------------------------------------#
   return(new.datum)
   #---------------------------------------------------------------------------------------#
}#end update
#==========================================================================================#
#==========================================================================================#
