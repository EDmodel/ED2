#==========================================================================================#
#==========================================================================================#
#     This function returns the cloud metrics from a point cloud.                          #
#------------------------------------------------------------------------------------------#
cloud.metrics <<- function( x
                          , zmah.stats = TRUE
                          , z1st.stats = TRUE
                          , tree.stats = TRUE
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
                          ){

   #----- Remove missing values. ----------------------------------------------------------#
   z     = x$z
   i     = x$intensity
   keep  = is.finite(z) & z >= zmin & z <= zmax & is.finite(i) & i >= intmin & i <= intmax
   z     = z[keep]
   i     = i[keep]
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
                                , tree.stats = tree.stats
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
         if (tree.stats) ans["tree.count"] = nz
      }#end if
      return(ans)
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Find the general metrics. -------------------------------------------------------#
   every = .Int.cloud.metrics( z     = z
                             , pref  = "elev"
                             , probs = probs
                             , zbreaks = zbreaks
                             , n.dens  = n.dens
                             , zl.dens = zl.dens
                             , zh.dens = zh.dens
                             , zzdens  = zzdens
                             , zmin    = zmin
                             , zmax    = zmax
                             )#end .Int.cloud.metrics
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
      zmah  = jitter( x      = sample(x=mhdens$x,size=nzmah,replace=TRUE,prob=mhdens$y)
                    , amount = 0.5*mean(diff(mhdens$x)))
      #------------------------------------------------------------------------------------#


      #----- Find the general metrics. ----------------------------------------------------#
      mah = .Int.cloud.metrics( z       = zmah
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
      z1st = x$z[x$retn.number %in% min(x$retn.number)]
      #------------------------------------------------------------------------------------#


      #----- Find the general metrics. ----------------------------------------------------#
      if (length(z1st) > min.pts){
         first = .Int.cloud.metrics( z       = z1st
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
      }else{
         first = .Int.cloud.metrics( z       = x$z
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
         first = first + NA
      }#end if (length(z1st) > min.pts)
      #------------------------------------------------------------------------------------#

   }else{
      first = NULL
   }#end if (z1st.stats)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find metrics for points flagged as medium and high vegetation only.               #
   #---------------------------------------------------------------------------------------#
   if (tree.stats){
      #------------------------------------------------------------------------------------#
      #     Keep only the first returns.                                                   #
      #------------------------------------------------------------------------------------#
      zveg = x$z[x$pt.class %in% c(4,5)]
      #------------------------------------------------------------------------------------#


      #----- Find the general metrics. ----------------------------------------------------#
      if (length(zveg) > min.pts){
         tree = .Int.cloud.metrics( z       = zveg
                                  , pref    = "tree"
                                  , probs   = probs
                                  , zbreaks = zbreaks
                                  , n.dens  = n.dens
                                  , zl.dens = zl.dens
                                  , zh.dens = zh.dens
                                  , zzdens  = NULL
                                  , zmin    = zmin
                                  , zmax    = zmax
                                  )#end .Int.cloud.metrics
      }else{
         tree = .Int.cloud.metrics( z       = x$z
                                  , pref    = "tree"
                                  , probs   = probs
                                  , zbreaks = zbreaks
                                  , n.dens  = n.dens
                                  , zl.dens = zl.dens
                                  , zh.dens = zh.dens
                                  , zzdens  = NULL
                                  , zmin    = zmin
                                  , zmax    = zmax
                                  )#end .Int.cloud.metrics
         tree = tree + NA
      }#end if
      #------------------------------------------------------------------------------------#

   }else{
      tree = NULL
   }#end if (tree.stats)
   #---------------------------------------------------------------------------------------#



   #----- Return the answer as a vector or a matrix, depending on the user's choice. ------#
   if (mat.out){
      ans = cbind(elev = every, zmah = mah, z1st = first, tree = tree)
   }else{
      ans = c(every,mah,first,tree)
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
.Int.cloud.metrics <<- function(z,pref,probs,zbreaks,n.dens,zl.dens,zh.dens,zzdens
                               ,zmin,zmax){


   #----- Make sure this function has been called by summnum.  Otherwise, stop. -----------#
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
         mess     = mess && ! ( substring(wcm[[n]],1,13) %==% "cloud.metrics" ||
                                substring(wcm[[n]],1,12) %==% "grid.metrics"  )
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
      bye =  paste( " Function .Int.cloud.metrics is internal,"
                  , " and can only be called by cloud.metrics or grid.metrics"
                  , sep = ""
                  )#end paste
      stop(bye)
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Lend of the point cloud. --------------------------------------------------------#
   nz = length(z)
   #---------------------------------------------------------------------------------------#


   #----- Labels for mode probability. ----------------------------------------------------#
   pmode = c("prob","pmah","p1st","pveg")[match(pref,c("elev","zmah","z1st","tree"))]
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
   zlabels = ifelse(is.finite(zbreaks),sprintf("%.1f",zbreaks),zbreaks)
   zlabels = paste(zlabels[-nbreaks],".to.",zlabels[-1],".m",sep="")
   #---------------------------------------------------------------------------------------#



   #----- Initialise the vector with some basic statistics. -------------------------------#
   ans            = list()
   ans$elev.count = nz
   ans$elev.mean  = mean(z)
   ans$elev.sdev  = sd  (z)
   ans$elev.skew  = skew(z)
   ans$elev.kurt  = kurt(z)
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
   if (is.null(zzdens)) zzdens = density(x=z,n=n.dens,from=zl.dens,to=zh.dens)
   zdens             = data.frame(x=zzdens$x,y=zzdens$y)
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
   #---------------------------------------------------------------------------------------#



   #----- Use heights to break the return signal in classes. ------------------------------#
   zcut         = cut(x=z,breaks=zbreaks,right=FALSE)
   zprop        = lapply( X   = tapply(X=z,INDEX=zcut,FUN=length,simplify=FALSE)
                        , FUN = '/'
                        , e2  = nz
                        )#end lapply
   zprop        = lapply( X   = zprop
                        , FUN = function(x) if(length(x) == 0){x = 0}else{x=x}
                        )#end lapply
   names(zprop) = paste("fret.elev.",zlabels,sep="")
   ans          = modifyList(x=ans,val=zprop)
   #---------------------------------------------------------------------------------------#




   #----- Use proportion to find the canopy fraction. -------------------------------------#
   zfcan         = as.list(rev(cumsum(rev(unlist(zprop)))))
   if (length(zfcan) != length(zlabels)) browser()
   names(zfcan)  = paste("fcan.elev.",zlabels,sep="")
   ans           = modifyList(x=ans,val=zfcan)
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
      zveg    = ifelse( pt.cloud$pt.class %in% c(3,4,5) & pt.cloud$pt.class %>=% zabove[1]
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
