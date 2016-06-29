#==========================================================================================#
#==========================================================================================#
#     This function returns the cloud metrics from a point cloud.                          #
#------------------------------------------------------------------------------------------#
cloud.metrics <<- function( x
                          , zmah.stats = TRUE
                          , z1st.stats = TRUE
                          , zveg.stats = TRUE
                          , zvfr.stats = TRUE
                          , probs      = c(0.01,0.05,0.10,0.25,0.50,0.75,0.90,0.95,0.99)
                          , zbreaks    = c(1.0,2.5,5.0,7.5,10.0,15.0,20.0,25.0,30.0)
                          , n.dens     = 256
                          , zl.dens    = 1
                          , zh.dens    = 50
                          , zo.dens    = 5
                          , zzdens     = NULL
                          , mhdens     = NULL
                          , zmin       = 0.
                          , zmax       = Inf
                          , intmin     = 0
                          , intmax     = Inf
                          , mat.out    = FALSE
                          , min.pts    = 500
                          , summ.only  = FALSE
                          ){

   #----- Remove missing values. ----------------------------------------------------------#
   z     = x$z
   i     = x$intensity
   ptc   = x$pt.class
   keep  = is.finite(z) & z >= zmin & z <= zmax & is.finite(i) & i >= intmin & i <= intmax
   z     = z[keep]
   i     = i[keep]
   ptc   = ptc[keep]
   x     = x[keep,]
   nz    = length(z)
   nzmah = sum(i)
   #---------------------------------------------------------------------------------------#


   #----- Return a dummy vector in case there aren't enough points. -----------------------#
   if (nz < min.pts){
      #----- Generate a dummy point cloud. ------------------------------------------------#
      dummy = data.frame( x                = runif(n=min.pts*2)
                        , y                = runif(n=min.pts*2)
                        , z                = runif(n=min.pts*2,min=zl.dens,max=zh.dens)
                        , intensity        = rep  (x=1     ,times=min.pts*2)
                        , retn.number      = rep  (x=c(1,2),each=min.pts)
                        , number.retn.gp   = rep  (x=1     ,times=min.pts*2)
                        , scan.dir.flag    = rep  (x=FALSE ,times=min.pts*2)
                        , edge.flight.line = rep  (x=FALSE ,times=min.pts*2)
                        , pt.class         = rep  (x=4     ,times=min.pts*2)
                        , synthetic        = rep  (x=FALSE ,times=min.pts*2)
                        , key.point        = rep  (x=FALSE ,times=min.pts*2)
                        , withheld         = rep  (x=FALSE ,times=min.pts*2)
                        , scan.anlge.rank  = rep  (x=1     ,times=min.pts*2)
                        , user.data        = rep  (x=0     ,times=min.pts*2)
                        , pt.source.ID     = rep  (x=1     ,times=min.pts*2)
                        , gpstime          = rep  (x=1     ,times=min.pts*2)
                        )#end data.frame
      #------------------------------------------------------------------------------------#



      #----- Call the cloud metrics for the dummy point cloud. ----------------------------#
      ans        = cloud.metrics( x          = dummy
                                , zmah.stats = zmah.stats
                                , z1st.stats = z1st.stats
                                , zveg.stats = zveg.stats
                                , zvfr.stats = zvfr.stats
                                , probs      = probs
                                , zbreaks    = zbreaks
                                , n.dens     = n.dens
                                , zl.dens    = zl.dens
                                , zh.dens    = zh.dens
                                , zo.dens    = zo.dens
                                , zzdens     = zzdens
                                , mhdens     = mhdens
                                , zmin       = zmin
                                , zmax       = zmax
                                , intmin     = intmin
                                , intmax     = intmax
                                , mat.out    = mat.out
                                , min.pts    = min.pts
                                , summ.only  = summ.only
                                )#end cloud.metrics
      #------------------------------------------------------------------------------------#



      #----- Throw away all information, except for the number of points. -----------------#
      ans        = ans * NA
      if (mat.out){
         ans[1,]           = nz
      }else{
         ans["elev.count"] = nz
         if (zmah.stats) ans["zmah.count"] = nz
         if (z1st.stats) ans["z1st.count"] = nz
         if (zveg.stats) ans["zveg.count"] = nz
         if (zvfr.stats) ans["zvfr.count"] = nz
      }#end if
      return(ans)
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Find the general metrics. -------------------------------------------------------#
   if (summ.only){
      every = .Int.summary.metrics( z     = z
                                  , ptc   = ptc
                                  , pref  = "elev"
                                  , probs = probs
                                  )#end .Int.cloud.metrics
   }else{
      every = .Int.cloud.metrics( z       = z
                                , ptc     = ptc
                                , pref    = "elev"
                                , probs   = probs
                                , zbreaks = zbreaks
                                , n.dens  = n.dens
                                , zl.dens = zl.dens
                                , zh.dens = zh.dens
                                , zzdens  = zzdens
                                , zmin    = zmin
                                , zmax    = zmax
                                )#end .Int.cloud.metrics
   }#end summ.only
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find height metrics weighted by intensity.                                        #
   #---------------------------------------------------------------------------------------#
   if (zmah.stats){
      #------------------------------------------------------------------------------------#
      #     Create a height dataset that follows the PDF given by the mhdens distribution  #
      # function.                                                                          #
      #------------------------------------------------------------------------------------#
      if (is.null(mhdens)){
         mhdens = macarthur.horn(pt.cloud=x,zl= zl.dens,zh=zh.dens,zo=zo.dens,nz=n.dens)
      }#end if (is.null(mhdens))
      zmah   = jitter( x      = sample(x=mhdens$x,size=nzmah,replace=TRUE,prob=mhdens$y)
                     , amount = 0.5*mean(diff(mhdens$x)))
      ptcmah = rep(5,times=nzmah)
      #------------------------------------------------------------------------------------#


      #----- Find the Mac-Arthurn Horn metrics. -------------------------------------------#
      if (summ.only){
         mah = .Int.summary.metrics( z     = zmah
                                   , ptc   = ptcmah
                                   , pref  = "zmah"
                                   , probs = probs
                                   )#end .Int.cloud.metrics
      }else{
         mah = .Int.cloud.metrics( z       = zmah
                                 , ptc     = ptcmah
                                 , pref    = "zmah"
                                 , probs   = probs
                                 , zbreaks = zbreaks
                                 , n.dens  = n.dens
                                 , zl.dens = zl.dens
                                 , zh.dens = zh.dens
                                 , zzdens  = NULL
                                 , zmin    = zmin
                                 , zmax    = zmax
                                 )#end .Int.cloud.metrics
      }#end summ.only
      #------------------------------------------------------------------------------------#

   }else{
      mah = NULL
   }#end if (zmah.stats)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find height metrics for first return only.                                        #
   #---------------------------------------------------------------------------------------#
   if (z1st.stats){
      #------------------------------------------------------------------------------------#
      #     Keep only the first returns.                                                   #
      #------------------------------------------------------------------------------------#
      sel1st = x$retn.number %in% min(x$retn.number)
      #----- Check whether to use the subset or the entire thing (to generate NA). --------#
      if (sum(sel1st) >= min.pts){
         z1st   = x$z[sel1st]
         ptc1st = x$pt.class[sel1st]
      }else{
         z1st   = x$z
         ptc1st = x$pt.class
      }#end if
      #------------------------------------------------------------------------------------#



      #----- Find the Mac-Arthurn Horn metrics. -------------------------------------------#
      if (summ.only){
         first = .Int.summary.metrics( z       = z1st
                                     , ptc     = ptc1st
                                     , pref    = "z1st"
                                     , probs   = probs
                                     )#end .Int.cloud.metrics
      }else{
         first = .Int.cloud.metrics( z       = z1st
                                   , ptc     = ptc1st
                                   , pref    = "z1st"
                                   , probs   = probs
                                   , zbreaks = zbreaks
                                   , n.dens  = n.dens
                                   , zl.dens = zl.dens
                                   , zh.dens = zh.dens
                                   , zzdens  = NULL
                                   , zmin    = zmin
                                   , zmax    = zmax
                                   )#end .Int.cloud.metrics
      }#end summ.only
      #------------------------------------------------------------------------------------#

      #------ Discard results in case the point cloud is too thin. ------------------------#
      if (sum(sel1st) < min.pts){
         first = first + NA
      }#end if (length(z1st) > min.pts)
      #------------------------------------------------------------------------------------#

   }else{
      first = NULL
   }#end if (z1st.stats)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find metrics only for those points flagged as vegetation.                         #
   #---------------------------------------------------------------------------------------#
   if (zveg.stats){
      #------------------------------------------------------------------------------------#
      #     Keep only the vegetation returns.                                              #
      #------------------------------------------------------------------------------------#
      selveg = x$pt.class %in% c(0,1,3,4,5)
      #----- Check whether to use the subset or the entire thing (to generate NA). --------#
      if (sum(selveg) >= min.pts){
         zveg   = x$z[selveg]
         ptcveg = x$pt.class[selveg]
      }else{
         zveg   = x$z
         ptcveg = x$pt.class
      }#end if
      #------------------------------------------------------------------------------------#



      #----- Find the vegetation metrics. -------------------------------------------------#
      if (summ.only){
         plant = .Int.summary.metrics( z       = zveg
                                     , ptc     = ptcveg
                                     , pref    = "zveg"
                                     , probs   = probs
                                     )#end .Int.cloud.metrics
      }else{
         plant = .Int.cloud.metrics( z       = zveg
                                   , ptc     = ptcveg
                                   , pref    = "zveg"
                                   , probs   = probs
                                   , zbreaks = zbreaks
                                   , n.dens  = n.dens
                                   , zl.dens = zl.dens
                                   , zh.dens = zh.dens
                                   , zzdens  = NULL
                                   , zmin    = zmin
                                   , zmax    = zmax
                                   )#end .Int.cloud.metrics
      }#end summ.only
      #------------------------------------------------------------------------------------#

      #------ Discard results in case the point cloud is too thin. ------------------------#
      if (sum(selveg) < min.pts){
         plant = plant + NA
      }#end if (length(zveg) > min.pts)
      #------------------------------------------------------------------------------------#

   }else{
      plant = NULL
   }#end if (zveg.stats)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find metrics only for those points flagged as vegetation AND first return.        #
   #---------------------------------------------------------------------------------------#
   if (zvfr.stats){
      #------------------------------------------------------------------------------------#
      #     Keep only the first returns that are vegetation.                               #
      #------------------------------------------------------------------------------------#
      selvfr = (x$pt.class %in% c(0,1,3,4,5)) & (x$retn.number %in% min(x$retn.number))
      #----- Check whether to use the subset or the entire thing (to generate NA). --------#
      if (sum(selvfr) >= min.pts){
         zvfr   = x$z[selvfr]
         ptcvfr = x$pt.class[selvfr]
      }else{
         zvfr   = x$z
         ptcvfr = x$pt.class
      }#end if
      #------------------------------------------------------------------------------------#



      #----- Find the Mac-Arthurn Horn metrics. -------------------------------------------#
      if (summ.only){
         pl1st = .Int.summary.metrics( z       = zvfr
                                     , ptc     = ptcvfr
                                     , pref    = "zvfr"
                                     , probs   = probs
                                     )#end .Int.cloud.metrics
      }else{
         pl1st = .Int.cloud.metrics( z       = zvfr
                                   , ptc     = ptcvfr
                                   , pref    = "zvfr"
                                   , probs   = probs
                                   , zbreaks = zbreaks
                                   , n.dens  = n.dens
                                   , zl.dens = zl.dens
                                   , zh.dens = zh.dens
                                   , zzdens  = NULL
                                   , zmin    = zmin
                                   , zmax    = zmax
                                   )#end .Int.cloud.metrics
      }#end summ.only
      #------------------------------------------------------------------------------------#



      #------ Discard results in case the point cloud is too thin. ------------------------#
      if (sum(selvfr) < min.pts){
         pl1st = pl1st + NA
      }#end if (length(zvfr) > min.pts)
      #------------------------------------------------------------------------------------#

   }else{
      pl1st = NULL
   }#end if (zvfr.stats)
   #---------------------------------------------------------------------------------------#



   #----- Return the answer as a vector or a matrix, depending on the user's choice. ------#
   if (mat.out){
      ans = cbind(elev = every, zmah = mah, z1st = first, zveg = plant, zvfr = pl1st)
   }else{
      ans = c(every,mah,first,plant,pl1st)
   }#end if
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end cloud.metrics
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      This is the internal routine that computes the metrics for a given dataset.  This   #
# routine cannot be called directly, it can only be called by cloud.metrics.               #
#------------------------------------------------------------------------------------------#
.Int.cloud.metrics <<- function(z,ptc,pref,probs,zbreaks,n.dens,zl.dens,zh.dens,zzdens
                               ,zmin,zmax){


   #----- Make sure this function has been called by cloud.metrics or grid.metrics. -------#
   patt  = "^([A-Za-z0-9]+)(\\({1})(.*)(\\){1})$"
   repl  = "\\1"
   n     = 0
   mess  = TRUE
   top   = FALSE
   wcm   = list()
   while (! top){
      n = n + 1
      wcm[[n]] = try( gsub( pattern     = patt
                          , replacement = repl
                          , x           = deparse(sys.call(-n))
                          )#end gsub
                    , silent = TRUE
                    )#end try
      if ("try-error" %in% is(wcm[[n]])){
         wcm[[n]] = NA
         top      = TRUE
      }else{
         #----- Not an error.  Check whether this has been called by a friend function. ---#
         wcm[[n]] = paste(wcm[[n]],collapse="")
         top      = substring(wcm[[n]],1,4) %==% "NULL"
         mess     = mess && ! ( grepl("cloud.metrics",wcm[[n]]) ||
                                grepl("grid.metrics" ,wcm[[n]])  )
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#
   }#end while
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Do not allow this function to continue in case it was an externall call.          #
   #---------------------------------------------------------------------------------------#
   if (mess){
      wcm =  sapply(X=wcm,FUN=rbind)
      print(wcm)
      bye =  paste0( " Function .Int.cloud.metrics is internal,"
                   , " and can only be called by cloud.metrics or grid.metrics"
                   )#end paste
      stop(bye)
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Lend of the point cloud. --------------------------------------------------------#
   nz = length(z)
   #---------------------------------------------------------------------------------------#


   #----- Labels for mode probability. ----------------------------------------------------#
   plist = c("prob","pmah","p1st","pveg","pvfr")
   zlist = c("elev","zmah","z1st","zveg","zvfr")
   pmode = plist[match(pref,zlist)]
   #---------------------------------------------------------------------------------------#


   #----- Labels for probabilities. -------------------------------------------------------#
   prob.names = paste("p",sprintf(ifelse(probs==1,"%3.3i","%2.2i"),round(100*probs)),sep="")
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Use heights to break the return signal in classes.  We only keep the breaks that #
   # are between bounds, and force the bounds to be break points.                          #
   #---------------------------------------------------------------------------------------#
   zbreaks = sort(unique(c(zmin,zbreaks,zmax)))
   zbreaks = zbreaks[zbreaks >= zmin & zbreaks <= zmax]
   nbreaks = length(zbreaks)
   llabels = ifelse(is.finite(zbreaks),sprintf("%.1f",zbreaks),zbreaks)
   zlabels = paste0(llabels[-nbreaks],".to.",llabels[-1],".m")
   alabels = paste0("above.",llabels[-nbreaks],".m")
   #---------------------------------------------------------------------------------------#



   #----- Initialise the vector with some basic statistics. -------------------------------#
   ans            = list()
   ans$elev.count = nz
   ans$elev.mean  = mean(z)
   ans$elev.sdev  = sd  (z)
   ans$elev.skew  = skew(z)
   ans$elev.kurt  = kurt(z)
   ans$elev.max   = max (z)
   #---------------------------------------------------------------------------------------#


   #----- Find quantiles. -----------------------------------------------------------------#
   quant        = as.list(quantile(x=z,probs=probs,names=FALSE))
   names(quant) = paste("elev",prob.names,sep=".")
   ans          = modifyList(x=ans,val=quant)
   #---------------------------------------------------------------------------------------#



   #----- Find the inter-quantile range. --------------------------------------------------#
   ans$elev.iqr = diff(quantile(x=z,probs=c(0.25,0.75),names=FALSE))
   #---------------------------------------------------------------------------------------#



   #----- Find the density function, and use it to retrieve the mode. ---------------------#
   if (is.null(zzdens)) zzdens = density.safe(x=z,n=n.dens,from=zl.dens,to=zh.dens)
   zdens             = data.frame(x=zzdens$x,y=zzdens$y)
   if (any(is.finite(zzdens$y))){
      dz                = mean(diff(zdens$x))
      spk               = peaks(zdens$y) & zdens$y %>% 1.e-10
      zpeaks            = zdens[spk,]
      o                 = order(zpeaks$y,decreasing=TRUE)
      zpeaks            = zpeaks[o,]
      top               = which.max(zpeaks$x)
      bot               = which.min(zpeaks$x)
      npeaks            = nrow(zpeaks)
      ans$elev.1st.mode = if (npeaks >= 1){zpeaks$x[  1]     }else{NA}
      ans$prob.1st.mode = if (npeaks >= 1){zpeaks$y[  1] * dz}else{NA}
      ans$elev.2nd.mode = if (npeaks >= 2){zpeaks$x[  2]     }else{NA}
      ans$prob.2nd.mode = if (npeaks >= 2){zpeaks$y[  2] * dz}else{NA}
      ans$elev.3rd.mode = if (npeaks >= 3){zpeaks$x[  3]     }else{NA}
      ans$prob.3rd.mode = if (npeaks >= 3){zpeaks$y[  3] * dz}else{NA}
      ans$elev.top.mode = if (npeaks >= 1){zpeaks$x[top]     }else{NA}
      ans$prob.top.mode = if (npeaks >= 1){zpeaks$y[top] * dz}else{NA}
      ans$elev.bot.mode = if (npeaks >= 1){zpeaks$x[bot]     }else{NA}
      ans$prob.bot.mode = if (npeaks >= 1){zpeaks$y[bot] * dz}else{NA}
   }else{
      ans$elev.1st.mode = NA
      ans$prob.1st.mode = NA
      ans$elev.2nd.mode = NA
      ans$prob.2nd.mode = NA
      ans$elev.3rd.mode = NA
      ans$prob.3rd.mode = NA
      ans$elev.top.mode = NA
      ans$prob.top.mode = NA
      ans$elev.bot.mode = NA
      ans$prob.bot.mode = NA
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Check whether we can compute tree cover.                                         #
   #---------------------------------------------------------------------------------------#
   if (is.null(ptc)){
      #----- Skip canopy height calculation. ----------------------------------------------#
      find.can = FALSE
      #------------------------------------------------------------------------------------#
   }else{
      #----- Discard data that are not classified as vegetation. --------------------------#
      zveg     = ifelse( ! ptc %in% 2, z,NA)
      find.can = any(is.finite(zveg))
      #------------------------------------------------------------------------------------#
   }#end if (is.null(ptc))
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #      Find tree cover if possible.                                                     #
   #---------------------------------------------------------------------------------------#
   if (find.can){
      #----- Use heights to break the return signal in classes. ---------------------------#
      zcut         = cut(x=zveg,breaks=zbreaks,right=FALSE)
      zprop        = lapply( X   = tapply(X=zveg,INDEX=zcut,FUN=length,simplify=FALSE)
                           , FUN = '/'
                           , e2  = nz
                           )#end lapply
      zprop        = lapply( X   = zprop
                           , FUN = function(x) if(length(x) == 0){x = 0}else{x=x}
                           )#end lapply
      names(zprop) = paste("fcan.elev.",zlabels,sep="")
      ans          = modifyList(x=ans,val=zprop)
      #------------------------------------------------------------------------------------#




      #----- Use proportion to find the canopy fraction. ----------------------------------#
      zfcan         = as.list(rev(cumsum(rev(unlist(zprop)))))
      if (length(zfcan) != length(alabels)) browser()
      names(zfcan)  = paste("fcan.elev.",alabels,sep="")
      ans           = modifyList(x=ans,val=zfcan)
      #------------------------------------------------------------------------------------#
   }else{
      #----- Point cloud classes have not been provided, make variables undefined. --------#
      zprop              = replicate(n=length(zlabels),list(NA))
      names(zprop)       = paste("fcan.elev.",zlabels,sep="")
      ans                = modifyList(x=ans,val=zprop)
      zfcan              = replicate(n=length(alabels),list(NA))
      names(zfcan)       = paste("fcan.elev.",alabels,sep="")
      ans                = modifyList(x=ans,val=zfcan)
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Coerce ans to a vector. ---------------------------------------------------------#
   ans          = sapply(X=ans,FUN=c)
   #---------------------------------------------------------------------------------------#


   #----- Fix the names. ------------------------------------------------------------------#
   names(ans) = gsub(pattern="elev",replacement=pref ,x=names(ans))
   names(ans) = gsub(pattern="prob",replacement=pmode,x=names(ans))
   #---------------------------------------------------------------------------------------#


   #----- Send metrics back to cloud.metrics. ---------------------------------------------#
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function .Int.cloud.metrics
#==========================================================================================#
#==========================================================================================#







#==========================================================================================#
#==========================================================================================#
#      This is the internal routine that computes the simple set of metrics.               #
#------------------------------------------------------------------------------------------#
.Int.summary.metrics <<- function(z,ptc,pref,probs){


   #----- Make sure this function has been called by cloud.metrics or grid.metrics. -------#
   patt  = "^([A-Za-z0-9]+)(\\({1})(.*)(\\){1})$"
   repl  = "\\1"
   n     = 0
   mess  = TRUE
   top   = FALSE
   wcm   = list()
   while (! top){
      n = n + 1
      wcm[[n]] = try( gsub( pattern     = patt
                          , replacement = repl
                          , x           = deparse(sys.call(-n))
                          )#end gsub
                    , silent = TRUE
                    )#end try
      if ("try-error" %in% is(wcm[[n]])){
         wcm[[n]] = NA
         top      = TRUE
      }else{
         #----- Not an error.  Check whether this has been called by a friend function. ---#
         wcm[[n]] = paste(wcm[[n]],collapse="")
         top      = substring(wcm[[n]],1,4) %==% "NULL"
         mess     = mess && ! ( grepl("cloud.metrics",wcm[[n]]) ||
                                grepl("grid.metrics" ,wcm[[n]])  )
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#
   }#end while
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Do not allow this function to continue in case it was an externall call.          #
   #---------------------------------------------------------------------------------------#
   if (mess){
      wcm =  sapply(X=wcm,FUN=rbind)
      print(wcm)
      bye =  paste0( " Function .Int.summary.metrics is internal,"
                   , " and can only be called by cloud.metrics or grid.metrics"
                   )#end paste
      stop(bye)
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Lend of the point cloud. --------------------------------------------------------#
   nz = length(z)
   #---------------------------------------------------------------------------------------#


   #----- Labels for probabilities. -------------------------------------------------------#
   prob.names = paste("p",sprintf(ifelse(probs==1,"%3.3i","%2.2i"),round(100*probs)),sep="")
   #---------------------------------------------------------------------------------------#



   #----- Initialise the vector with some basic statistics. -------------------------------#
   ans            = list( elev.count = nz
                        , elev.mean  = mean(z)
                        , elev.max   = max (z)
                        )#end list
   #---------------------------------------------------------------------------------------#


   #----- Find quantiles. -----------------------------------------------------------------#
   quant        = as.list(quantile(x=z,probs=probs,names=FALSE))
   names(quant) = paste("elev",prob.names,sep=".")
   ans          = modifyList(x=ans,val=quant)
   #---------------------------------------------------------------------------------------#



   #----- Coerce ans to a vector. ---------------------------------------------------------#
   ans          = sapply(X=ans,FUN=c)
   #---------------------------------------------------------------------------------------#


   #----- Fix the names. ------------------------------------------------------------------#
   names(ans) = gsub(pattern="elev",replacement=pref ,x=names(ans))
   #---------------------------------------------------------------------------------------#


   #----- Send metrics back to cloud.metrics. ---------------------------------------------#
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function .Int.summary.metrics
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#      This function finds the canopy opening as a function of height breaks.              #
#------------------------------------------------------------------------------------------#
open.fcan <<- function( pt.cloud
                      , zabove    = seq(from=0,to=35,by=5)
                      ){

   #----- Breaks.  Ensure the last break is Infinity. -------------------------------------#
   if (any(! is.finite(zabove))){
      stop("zabove contains non-finite entries, which is forbidden!")
   }#end if
   zbreaks = sort(c(zabove,Inf))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #       Check whether the point cloud has any data in it.                               #
   #---------------------------------------------------------------------------------------#
   if (nrow(pt.cloud) == 0){
      ans = rep(NA,times=length(zabove))
   }else{
      #----- Discard data that are not classified as vegetation. --------------------------#
      zveg    = ifelse( pt.cloud$pt.class %in% c(0,1,3,4,5) & pt.cloud$z %>=% zabove[1]
                      , pt.cloud$z
                      , NA
                      )#end ifelse
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     If everything became NA, then all returns were ground returns, assume no       #
      # canopy cover.                                                                      #
      #------------------------------------------------------------------------------------#
      if (all(is.na(zveg))){
         #----- Set all classes to zero. --------------------------------------------------#
         ans = rep(0,times=length(zabove))
         #---------------------------------------------------------------------------------#
      }else{

         #----- Split elevation in bins. --------------------------------------------------#
         zcut    = as.integer(cut(x=zveg,breaks=zbreaks,labels=zabove,right=FALSE))
         #---------------------------------------------------------------------------------#


         #----- Count layers, and return the fraction above each height layer. ------------#
         tabnow   = table(zcut)
         ans      = rep(0,times=length(zabove))
         idx      = as.numeric(names(tabnow))
         ans[idx] = tabnow
         ans      = rev(cumsum(rev(ans))) / length(zcut)
         #---------------------------------------------------------------------------------#
      }#end if (all(is.na(zveg)))
      #------------------------------------------------------------------------------------#
   }#end if (nrow(pt.cloud) == 0)
   #---------------------------------------------------------------------------------------#


   #----- Make names based on zabove, then return the answer. -----------------------------#
   fmt        = ceiling(log10(max(abs(zabove))*(1.+10*.Machine$double.eps)))
   if (all(zabove == as.integer(zabove))){
      fmt        = paste("%0",fmt,".",fmt,"i",sep="")
   }else{
      fmt        = paste("%0",fmt+3,".",2,"f",sep="")
   }#end if
   names(ans) = paste("fcan",sprintf(fmt,zabove),sep=".")
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end frac.can.lidar
#==========================================================================================#
#==========================================================================================#
