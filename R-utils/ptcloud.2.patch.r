#==========================================================================================#
#==========================================================================================#
#     This function determines cohort distribution based on a point cloud, and target      #
# values for total aboveground biomass, basal area, stem (demographic) density and a look- #
# -up table of relative PFT composition as a function of height.                           #
#------------------------------------------------------------------------------------------#
ptcloud.2.patch <<- function( pt.cloud
                            , pname
                            , zo
                            , zh
                            , nz
                            , lai.goal
                            , bsa.goal
                            , agb.goal
                            , acd.goal = agb.goal
                            , npl.goal
                            , w.lai
                            , w.bsa
                            , w.agb
                            , w.npl
                            , lookup
                            , necro
                            , dist.type
                            , dist.age
                            , syear
                            , tall.at.zh     = TRUE
                            , use.intensity  = FALSE
                            , lai.pst        = 2.0
                            , pft.pst        = 1
                            , use.dbh0.std   = TRUE
                            , dbh0.std       = 10.0
                            , dbh0.min       = 5.0
                            , pft.def        = 3
                            , use.lookup     = FALSE
                            , use.net.method = c("linear","ratio","fixed")
                            , f.max.oth      = 10^0.50
                            , f.max.npl      = 10^0.75
                            , b.min.oth      = 1. / f.max.oth
                            , b.min.npl      = 1. / f.max.npl
                            , f.net.def      = 1.0
                            , wfsim          = TRUE
                            ){


   #---------------------------------------------------------------------------------------#
   #    Standardise the net method. "Linear" is default as it seems less biased for        #
   # biomass and basal area (but more biased for stem density).                            #
   #---------------------------------------------------------------------------------------#
   use.net.method = match.arg(use.net.method)
   #---------------------------------------------------------------------------------------#




   #------ Run MacArthur-Horn (1969) correction to profiles. ------------------------------#
   mh.now = macarthur.horn( pt.cloud      = pt.cloud
                          , zo            = zo
                          , zh            = zh
                          , nz            = nz
                          , tall.at.zh    = tall.at.zh
                          , use.intensity = use.intensity
                          , wfsim         = wfsim
                          )#end macarthur.horn 
   #---------------------------------------------------------------------------------------#




   #----- Swap order of profiles. ---------------------------------------------------------#
   nzprof     = nrow(mh.now)
   t2b        = rev(sequence(nzprof))
   mh.now     = mh.now[t2b,]
   zmcprof    = mh.now$z
   dzprof     = 0*zmcprof - median(diff(zmcprof))
   hgtprof    = zmcprof / pft$b1Mh[pft.def]
   dbhprof    = h2dbh(h=hgtprof,ipft=pft.def)
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #   Find the most similar profile based on the lookup table.                            #
   #---------------------------------------------------------------------------------------#
   ks.stat    = apply(X=lookup$cdf,MARGIN=1,FUN=max.abs.diff,y=mh.now$cdf)
   o          = which.min(ks.stat)
   if (use.dbh0.std){
      dbh.cutoff = dbh0.std
   }else{
      dbh.cutoff = lookup$dbh.cutoff[o]
   }#end if (use.dbh0.std)
   pftprof    = lookup$pft[o,,]
   ladprof    = lookup$lad[o,]
   indiv      = lookup$idv
   #---------------------------------------------------------------------------------------#


   #----- Find necromass associated with standing dead trees. -----------------------------#
   acd.necro = max(0.,acd.goal-agb.goal) / pft$agf.bs[pft.def]
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Check whether to use the standard look-up table or to optimise scales.            #
   #---------------------------------------------------------------------------------------#
   if (use.lookup){
      #------------------------------------------------------------------------------------#
      #      Find cohort properties.                                                       #
      #------------------------------------------------------------------------------------#
      dsel   = indiv$identifier %in% rownames(lookup$lad)[o]
      i.now  = indiv[dsel,,drop=FALSE]
      ni.now = nrow(i.now)
      
      dtype  = dist.type
      dage   = dist.age
      if (ni.now > 0){
         #----- Sort cohorts. -------------------------------------------------------------#
         oi     = order(-i.now$DBH,i.now$pft)
         i.now  = i.now[oi,,drop=FALSE]
         #---------------------------------------------------------------------------------#



         #----- Create PSS and CSS structures. --------------------------------------------#
         cssnow = data.frame( time             = rep(syear,times=ni.now)
                            , patch            = rep(pname,times=ni.now)
                            , cohort           = sequence(ni.now)
                            , dbh              = i.now$DBH
                            , hite             = i.now$Htot
                            , pft              = i.now$pft
                            , n                = i.now$nplant
                            , bdead            = i.now$bdead
                            , balive           = i.now$balive
                            , lai              = i.now$lai
                            , stringsAsFactors = FALSE
                            )#end data.frame
         pssnow = data.frame( time             = syear
                            , patch            = pname
                            , trk              = dist.type - 1
                            , age              = dist.age
                            , area             = 1.0
                            , water            = 0.
                            , fsc              = necro$fsc
                            , stsc             = necro$stsc + acd.necro
                            , stsl             = necro$stsl + acd.necro
                            , ssc              = necro$ssc
                            , lai              = sum(cssnow$lai)
                            , isn              = necro$isn
                            , fsn              = necro$fsn
                            , msc              = necro$msc
                            , gpp              = 0.
                            , rh               = 0.
                            , stringsAsFactors = FALSE
                            )#end data.frame
         #---------------------------------------------------------------------------------#
      }else{
         #----- Create empty CSS and a PSS data frame. ------------------------------------#
         cssnow = data.frame( time             = numeric(0)
                            , patch            = character(0)
                            , cohort           = numeric(0)
                            , dbh              = numeric(0)
                            , hite             = numeric(0)
                            , pft              = numeric(0)
                            , n                = numeric(0)
                            , bdead            = numeric(0)
                            , balive           = numeric(0)
                            , lai              = numeric(0)
                            , stringsAsFactors = FALSE
                            )#end data.frame
         pssnow = data.frame( time             = syear
                            , patch            = pname
                            , trk              = dist.type - 1
                            , age              = dist.age
                            , area             = 1.0
                            , water            = 0.
                            , fsc              = necro$fsc
                            , stsc             = necro$stsc + acd.necro
                            , stsl             = necro$stsl + acd.necro
                            , ssc              = necro$ssc
                            , lai              = 0.
                            , isn              = necro$isn
                            , fsn              = necro$fsn
                            , msc              = necro$msc
                            , gpp              = 0.
                            , rh               = 0.
                            , stringsAsFactors = FALSE
                            )#end data.frame
      }#end if (nidx > 0)
      #------------------------------------------------------------------------------------#
   }else{

      #------------------------------------------------------------------------------------#
      #     Retrieve PFTs to consider.                                                     #
      #------------------------------------------------------------------------------------#
      pft.keys = toupper(dimnames(pftprof)[[2]])
      #------------------------------------------------------------------------------------#





      #------------------------------------------------------------------------------------#
      #     Grab current pdf profile, and find maxima, minima, and inflection points.      #
      #------------------------------------------------------------------------------------#
      lsmnow = if(wfsim){mh.now$lad}else{mh.now$pdf}
      cfz = curve.features(x=lsmnow,span=5)
      if (sum(c(cfz$max,cfz$iph)) == 0){
      
         #----- Reduce window to make peaks more likely. ----------------------------------#
         cfz = curve.features(x=lsmnow,span=3)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     In case a narrow failed to capture peaks, use function peaks.               #
         #---------------------------------------------------------------------------------#
         if (sum(c(cfz$max,cfz$iph)) == 0){
            iterate = TRUE
            nspan   = 17
            while (iterate){
               nspan = nspan - 2
               pktry = peaks(series=lsmnow,span=nspan)
               iterate = sum(pktry) == 0 & nspan > 3
            }#end while
         
            if (sum(pktry) == 0){
               stop("Failed finding any peaks in this distribution!")
            }else{
               #----- Make dummy data frame. ----------------------------------------------#
               cfz = data.frame( max = pktry
                               , min = rep(FALSE,times=nzprof)
                               , iph = rep(FALSE,times=nzprof)
                               , ipv = rep(FALSE,times=nzprof)
                               )#end data.frame
               #---------------------------------------------------------------------------#
            }#end if (sum(pktry) == 0)
            #------------------------------------------------------------------------------#
         }#end if (sum(c(cfz$max,cfz$iph)) == 0)
         #---------------------------------------------------------------------------------#
      }#end if (sum(c(cfz$max,cfz$iph)) == 0)
      #------------------------------------------------------------------------------------#



      #----- Indices of interesting points. -----------------------------------------------#
      imx = which(cfz$max)
      imn = which(cfz$min)
      iih = which(cfz$iph)
      iiv = which(cfz$ipv)
      #------------------------------------------------------------------------------------#



      #----- Find the peaks of the PDF distribution. --------------------------------------#
      ipk  = sort(unique(c(imx,iih)))
      ncoh = length(ipk)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #   Make sure that at least one peak is defined.                                     #
      #------------------------------------------------------------------------------------#
      ibnd = c( mapply( FUN      = left.bnd
                      , k        = sequence(ncoh)
                      , MoreArgs = list(ipk=ipk,imn=imn,iiv=iiv)
                      )#end mapply
              , nzprof
              )#end c
      #------------------------------------------------------------------------------------#
      dz   = diff(c(0,ibnd[-1]))
      pcoh = rep(ipk,times=dz)
      icoh = rep(sequence(ncoh),times=dz)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Run the procedure twice, one including only levels with DBH > dbh.cutoff, and  #
      # thesecond for the entire profile.                                                  #
      #------------------------------------------------------------------------------------#
      for (cstep in c("calibration","prediction")){

         #----- Select layers to use. -----------------------------------------------------#
         if (cstep %in% "calibration"){
            hgt.cutoff = dbh2h(dbh=dbh.cutoff,ipft=pft.def)
            hgt.bottom = dbh2h(dbh=dbh0.min  ,ipft=pft.def)
         }else{
            hgt.cutoff = dbh2h(dbh=dbh0.min  ,ipft=pft.def)
            hgt.bottom = hgt.cutoff
         }#end if (cstep %in% "calibration")
         k.sel = which(hgtprof[ipk] %ge% hgt.cutoff)
         k.bot = which(hgtprof[ipk] %ge% hgt.bottom)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #    In case all peaks are lower than the minimum standard height, assume that we #
         # cannot correct biases for this plot.                                            #
         #---------------------------------------------------------------------------------#
         if ((length(k.sel) == 0) && (length(k.bot) == 0)){
            #------------------------------------------------------------------------------#
            #      All the signal is underneath the minimum "visible" height.  Turn cohort #
            # into a grassland.                                                            #
            #------------------------------------------------------------------------------#
            sel  = rep(FALSE,times=length(hgtprof))
            bel  = rep(FALSE,times=length(hgtprof))
            #------------------------------------------------------------------------------#
         }else if (length(k.sel) == 0){
            #----- Tallest cohort is outside range of DBH.  Don't constrain. --------------#
            k.idx = ibnd[max(k.bot) + 1]
            sel   = hgtprof %ge% hgtprof[k.idx]
            bel   = sel
            #------------------------------------------------------------------------------#
         }else{
            #----- Keep only cohorts whose peak is above minimum DBH. ---------------------#
            k.idx = ibnd[max(k.sel) + 1]
            b.idx = ibnd[max(k.bot) + 1]
            sel   = hgtprof %ge% hgtprof[k.idx]
            bel   = hgtprof %ge% hgtprof[b.idx]
            #------------------------------------------------------------------------------#
         }#end if ((length(k.sel) == 0) && (length(k.bot) == 0))
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Split layers by cohort, in case any cohort is qualified.                   #
         #---------------------------------------------------------------------------------#
         if (any(sel)){
            #------------------------------------------------------------------------------#
            #      Split layers by cohort.                                                 #
            #------------------------------------------------------------------------------#
            c.ulai = split(x=mh.now$lad[sel]*dzprof[sel] ,f=icoh[sel])
            c.uipk = split(x=pcoh[sel]                   ,f=icoh[sel])
            b.ulai = split(x=mh.now$lad[bel]*dzprof[bel] ,f=icoh[bel])
            b.uipk = split(x=pcoh[bel]                   ,f=icoh[bel])
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #    Sometimes the density function creates blocks with LAI=0.  Drop them as   #
            # they are meaningless.                                                        #
            #------------------------------------------------------------------------------#
            #------ Selected data. --------------------------------------------------------#
            ulai = sapply(X=c.ulai,FUN=sum      )
            uipk = sapply(X=c.uipk,FUN=commonest)
            keep = ulai %gt% 0.
            if (any(keep)){
               ulai = c(ulai [keep],0)
               uipk = uipk   [keep]
               uhgt = c(hgtprof[uipk],pft$hgt.max[pft.pst])
               fpft = pftprof[uipk,,drop=FALSE]
               fpft = cbind(0,fpft)
               fpft = rbind(fpft,0)
            }else{
               ulai = 0.
               uhgt = pft$hgt.max[pft.pst]
               fpft = cbind(0,0 * pftprof[1,,drop=FALSE])
            }#end if (length(fpft) > 0)
            #------------------------------------------------------------------------------#


            #------ Total data. -----------------------------------------------------------#
            blai = sapply(X=b.ulai,FUN=sum      )
            bipk = sapply(X=b.uipk,FUN=commonest)
            beep = blai %gt% 0.
            if (any(beep)){
               blai = c(blai [beep],0)
               bipk = bipk   [beep]
               bhgt = c(hgtprof[bipk],pft$hgt.max[pft.pst])
               bpft = pftprof[bipk,,drop=FALSE]
               bpft = cbind(0,bpft)
               bpft = rbind(bpft,0)
            }else{
               blai = 0.
               bhgt = pft$hgt.max[pft.pst]
               bpft = cbind(0,0 * pftprof[1,,drop=FALSE])
            }#end if (length(fpft) > 0)
            #------------------------------------------------------------------------------#
         }else{
            #------------------------------------------------------------------------------#
            #    Nothing to include, turn this into a grassland.                           #
            #------------------------------------------------------------------------------#
            keep = FALSE
            ulai = 0.
            uhgt = pft$hgt.max[pft.pst]
            fpft = cbind(0,0 * pftprof[1,,drop=FALSE])
            beep = FALSE
            blai = 0.
            bhgt = pft$hgt.max[pft.pst]
            bpft = cbind(0,0 * pftprof[1,,drop=FALSE])
            #------------------------------------------------------------------------------#
         }#end if (any(sel))
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Append layer for tall grasses (set to zero by default)                      #
         #---------------------------------------------------------------------------------#
         ncoh = sum(keep) + 1
         bcoh = sum(beep) + 1
         #---------------------------------------------------------------------------------#


         #----- Append grasses to the pft list. -------------------------------------------#
         mypfts   = c(pft.pst,match(pft.keys,pft$key))
         npfts    = length(mypfts)
         #---------------------------------------------------------------------------------#


         #----- Rename fpft. --------------------------------------------------------------#
         colnames(fpft) = c(pft$key[mypfts])
         colnames(bpft) = c(pft$key[mypfts])
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Find uncalibrated properties.  Because ED2 has a cap in height and leaf    #
         # biomass, we find the equivalent DBH that would produce the same above-ground    #
         # biomass when the lidar height exceeds ED2 maximum height.                       #
         #---------------------------------------------------------------------------------#
         udbx     = h2dbh(h=uhgt,ipft=pft.def)              + 0. * uhgt
         udbh     = size2de(dbh=udbx,hgt=uhgt,ipft=pft.def) + 0. * uhgt
         uhgt.bnd = dbh2h(dbh=udbh,ipft=pft.def)            + 0. * uhgt
         unpl     = ulai / la.SL(dbh=udbh,height=uhgt.bnd)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Find the uncalibrated biomass and basal area.                              #
         #---------------------------------------------------------------------------------#
         npl.pft    = rep(unpl           , times= npfts)      * fpft
         dbh.pft    = rep(udbh           , times= npfts) + 0. * fpft
         hgt.pft    = rep(uhgt.bnd       , times= npfts) + 0. * fpft
         ipft.pft   = rep(mypfts         , each = ncoh ) + 0L * fpft
         wdns.pft   = rep(pft$rho[mypfts], each = ncoh ) + 0. * fpft
         sla.pft    = rep(pft$SLA[mypfts], each = ncoh ) + 0. * fpft
         bleaf.pft  = ( size2bl(dbh=dbh.pft,hgt=hgt.pft,sla=sla.pft,ipft=ipft.pft)
                      + 0. * fpft
                      )#end bleaf.pft
         bdead.pft  = size2bd(dbh=dbh.pft,hgt=hgt.pft,ipft=ipft.pft)  + 0. * fpft
         broot.pft  = rep(pft$qroot[mypfts],each = ncoh) * bleaf.pft
         bsw.pft    = rep(pft$qsw  [mypfts],each = ncoh) * hgt.pft   * bleaf.pft
         bbark.pft  = rep(pft$qbark[mypfts],each = ncoh) * hgt.pft   * bleaf.pft
         balive.pft = bleaf.pft + broot.pft + bsw.pft + bbark.pft
         agb.pft    = npl.pft * agb.SL(dbh=dbh.pft,height=hgt.pft,wdens=wdns.pft)
         bsa.pft    = npl.pft * pi/4 * dbh.pft^2
         lai.pft    = npl.pft * sla.pft * bleaf.pft
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Find uncalibrated properties.  Because ED2 has a cap in height and leaf    #
         # biomass, we find the equivalent DBH that would produce the same above-ground    #
         # biomass when the lidar height exceeds ED2 maximum height.                       #
         #---------------------------------------------------------------------------------#
         bdbx     = h2dbh(h=bhgt,ipft=pft.def)              + 0. * bhgt
         bdbh     = size2de(dbh=bdbx,hgt=bhgt,ipft=pft.def) + 0. * bhgt
         bhgt.bnd = dbh2h(dbh=bdbh,ipft=pft.def)            + 0. * bhgt
         bnpl     = blai / la.SL(dbh=bdbh,height=bhgt.bnd)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Find the uncalibrated biomass and basal area.                              #
         #---------------------------------------------------------------------------------#
         npl.bft    = rep(bnpl           , times= npfts)      * bpft
         dbh.bft    = rep(bdbh           , times= npfts) + 0. * bpft
         hgt.bft    = rep(bhgt.bnd       , times= npfts) + 0. * bpft
         ipft.bft   = rep(mypfts         , each = bcoh ) + 0L * bpft
         wdns.bft   = rep(pft$rho[mypfts], each = bcoh ) + 0. * bpft
         sla.bft    = rep(pft$SLA[mypfts], each = bcoh ) + 0. * bpft
         bleaf.bft  = ( size2bl(dbh=dbh.bft,hgt=hgt.bft,sla=sla.bft,ipft=ipft.bft)
                      + 0. * bpft
                      )#end bleaf.bft
         bdead.bft  = size2bd(dbh=dbh.bft,hgt=hgt.bft,ipft=ipft.bft)  + 0. * bpft
         broot.bft  = rep(pft$qroot[mypfts],each = bcoh) * bleaf.bft
         bsw.bft    = rep(pft$qsw  [mypfts],each = bcoh) * hgt.bft   * bleaf.bft
         bbark.bft  = rep(pft$qbark[mypfts],each = bcoh) * hgt.bft   * bleaf.bft
         balive.bft = bleaf.bft + broot.bft + bsw.bft + bbark.bft
         agb.bft    = npl.bft * agb.SL(dbh=dbh.bft,height=hgt.bft,wdens=wdns.bft)
         bsa.bft    = npl.bft * pi/4 * dbh.bft^2
         lai.bft    = npl.bft * sla.bft * bleaf.bft
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #    Different steps depending on whether this is the first or second run.        #
         #---------------------------------------------------------------------------------#
         suspicious = FALSE
         if (cstep %in% "calibration"){
            #------------------------------------------------------------------------------#
            #     Find calibration factors for pft and bsa, based on cross-validation      #
            # algorithms.                                                                  #
            #------------------------------------------------------------------------------#
            #----- Stem density.  Check for singularities before finding the global f. ----#
            if (sum(npl.pft) == 0){
               f.npl = f.net.def
               w.npl = 0.
            }else{
               f.npl = npl.goal / sum(npl.pft)
            }#end if (sum(npl.pft) == 0)
            #----- Basal area.  Check for singularities before finding the global f. ------#
            if (sum(bsa.pft) == 0){
               f.bsa = f.net.def
               w.bsa = 0.
            }else{
               f.bsa = bsa.goal / sum(bsa.pft)
            }#end if (sum(bsa.pft) == 0)
            #----- Leaf area index.  Check for singularities before finding the global f. -#
            if (sum(lai.pft) == 0){
               f.lai = f.net.def
               w.lai = 0.
            }else{
               f.lai = lai.goal / sum(lai.pft)
            }#end if (sum(bsa.pft) == 0)
            #----- Biomass.  Check for singularities before finding the global f. ---------#
            if (sum(agb.pft) == 0){
               f.agb = f.net.def
               w.agb = 0.
            }else{
               f.agb = agb.goal / sum(agb.pft)
            }#end if (sum(agb.pft) == 0)
            #----- Net, check for singularities. ------------------------------------------#
            w.all = w.npl + w.lai + w.bsa + w.agb
            if ( (use.net.method %in% "fixed") ||  (w.all == 0)){
               f.net = f.net.def
            }else if (use.net.method %in% "ratio"){
               s.bsa = f.bsa^2
               s.lai = f.lai^2
               s.agb = f.agb^2
               s.npl = f.npl^2
               f.net = ( f.npl * f.lai * f.bsa * f.agb
                       * ( w.npl*f.lai*f.bsa*f.agb + f.npl*w.lai*f.bsa*f.agb
                         + f.npl*f.lai*w.bsa*f.agb + f.npl*f.lai*f.bsa*w.agb )
                       / ( w.npl*s.lai*s.bsa*s.agb + s.npl*w.lai*s.bsa*s.agb
                         + s.npl*s.lai*w.bsa*s.agb + s.npl*s.lai*s.bsa*w.agb )
                       )#end f.net
            }else{
               f.net = (w.npl*f.npl+w.lai*f.lai+w.bsa*f.bsa+w.agb*f.agb) / w.all
            }#end if ( (use.net.method %in% "fixed") ||  (w.all == 0))
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Find calibration factors if we include everything.                       #
            #------------------------------------------------------------------------------#
            #----- Stem density.  Check for singularities before finding the global f. ----#
            if (sum(npl.bft) == 0){
               b.npl = f.net.def
               m.npl = 0.
            }else{
               b.npl = npl.goal / sum(npl.bft)
               m.npl = w.npl
            }#end if (sum(npl.pft) == 0)
            #----- Basal area.  Check for singularities before finding the global f. ------#
            if (sum(bsa.bft) == 0){
               b.bsa = f.net.def
               m.bsa = 0.
            }else{
               b.bsa = bsa.goal / sum(bsa.bft)
               m.bsa = w.bsa
            }#end if (sum(bsa.bft) == 0)
            #----- Leaf area index.  Check for singularities before finding the global f. -#
            if (sum(lai.bft) == 0){
               b.lai = f.net.def
               m.lai = 0.
            }else{
               b.lai = lai.goal / sum(lai.bft)
               m.lai = w.lai
            }#end if (sum(bsa.bft) == 0)
            #----- Biomass.  Check for singularities before finding the global f. ---------#
            if (sum(agb.bft) == 0){
               b.agb = f.net.def
               m.agb = 0.
            }else{
               b.agb = agb.goal / sum(agb.bft)
               m.agb = w.agb
            }#end if (sum(agb.pft) == 0)
            #----- Net, check for singularities. ------------------------------------------#
            m.all = m.npl + m.lai + m.bsa + m.agb
            if ( (use.net.method %in% "fixed") ||  (m.all == 0)){
               b.net = f.net.def
            }else if (use.net.method %in% "ratio"){
               z.bsa = b.bsa^2
               z.lai = b.lai^2
               z.agb = b.agb^2
               z.npl = b.npl^2
               b.net = ( b.npl * b.lai * b.bsa * b.agb
                       * ( m.npl*b.lai*b.bsa*b.agb + b.npl*m.lai*b.bsa*b.agb
                         + b.npl*b.lai*m.bsa*b.agb + b.npl*b.lai*b.bsa*m.agb )
                       / ( m.npl*z.lai*z.bsa*z.agb + z.npl*m.lai*z.bsa*z.agb
                         + z.npl*z.lai*m.bsa*z.agb + z.npl*z.lai*z.bsa*m.agb )
                       )#end b.net
            }else{
               b.net = (m.npl*b.npl+m.lai*b.lai+m.bsa*b.bsa+m.agb*b.agb) / m.all
            }#end if ( (use.net.method %in% "fixed") ||  (m.all == 0))
            #------------------------------------------------------------------------------#


            #----- Check for degenerate solution. -----------------------------------------#
            suspicious = (sum(lai.pft*f.net) %gt% 12)
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Ignore correction in case the scaling factor above threshold has very    #
            # different magnitude than scaling factor if we considered the entire profile. #
            # This indicates that most of the returns are coming from below the            #
            # calibration layer (DBH > 10cm) so it is really unconstrained.                #
            #------------------------------------------------------------------------------#
            x.bsa = (f.bsa %ge% f.max.oth) && (b.bsa %le% b.min.oth)
            x.lai = (f.lai %ge% f.max.oth) && (b.lai %le% b.min.oth)
            x.agb = (f.agb %ge% f.max.oth) && (b.agb %le% b.min.oth)
            x.npl = (f.npl %ge% f.max.npl) && (b.npl %le% b.min.npl)
            x.net = (f.net %ge% f.max.oth) && (b.net %le% b.min.oth)
            x.any = x.npl || x.bsa || x.lai || x.agb || x.net
            f.bsa = ifelse(test=x.any,yes=f.net.def,no=f.bsa)
            f.lai = ifelse(test=x.any,yes=f.net.def,no=f.lai)
            f.agb = ifelse(test=x.any,yes=f.net.def,no=f.agb)
            f.npl = ifelse(test=x.any,yes=f.net.def,no=f.npl)
            f.net = ifelse(test=x.any,yes=f.net.def,no=f.net)
            #------------------------------------------------------------------------------#

         }else{

            #------------------------------------------------------------------------------#
            #      Scale plant count and LAI.                                              #
            #------------------------------------------------------------------------------#
            npl.pft = npl.pft * f.net
            lai.pft = lai.pft * f.net
            #------------------------------------------------------------------------------#

            #------------------------------------------------------------------------------#
            #      In case this is a pasture and has trees, make it abandoned.             #
            #------------------------------------------------------------------------------#
            if (any(is.na(lai.pft)) || is.na(dist.type)) browser()
            if (sum(lai.pft) > 0 && dist.type == 1) dist.type = 5
            grass.lai          = max(0.,lai.pst - sum(lai.pft))
            igrs               = which(colnames(lai.pft) %in% pft$key[pft.pst])
            lai.pft[ncoh,igrs] = grass.lai
            npl.pft[ncoh,igrs] = grass.lai / (pft$SLA[pft.pst] * bleaf.pft[ncoh,igrs])
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Find cohort properties.                                                 #
            #------------------------------------------------------------------------------#
            use    = which(npl.pft > 0)
            idx    = use[order(dbh.pft[use],decreasing=TRUE)]
            nidx   = length(idx)
            dtype  = dist.type
            dage   = dist.age
            if (nidx > 0){
               cssnow = data.frame( time             = rep(syear,times=nidx)
                                  , patch            = rep(pname,times=nidx)
                                  , cohort           = sequence(nidx)
                                  , dbh              = dbh.pft   [idx]
                                  , hite             = hgt.pft   [idx]
                                  , pft              = ipft.pft  [idx]
                                  , n                = npl.pft   [idx]
                                  , bdead            = bdead.pft [idx]
                                  , balive           = balive.pft[idx]
                                  , lai              = lai.pft   [idx]
                                  , stringsAsFactors = FALSE
                                  )#end data.frame
               pssnow = data.frame( time             = syear
                                  , patch            = pname
                                  , trk              = dist.type - 1
                                  , age              = dist.age
                                  , area             = 1.0
                                  , water            = 0.
                                  , fsc              = necro$fsc
                                  , stsc             = necro$stsc + acd.necro
                                  , stsl             = necro$stsl + acd.necro
                                  , ssc              = necro$ssc
                                  , lai              = sum(cssnow$lai)
                                  , msn              = necro$isn
                                  , fsn              = necro$fsn
                                  , msc              = necro$msc
                                  , gpp              = 0.
                                  , rh               = 0.
                                  , stringsAsFactors = FALSE
                                  )#end data.frame
            }else{
               cssnow = data.frame( time             = numeric(0)
                                  , patch            = character(0)
                                  , cohort           = numeric(0)
                                  , dbh              = numeric(0)
                                  , hite             = numeric(0)
                                  , pft              = numeric(0)
                                  , n                = numeric(0)
                                  , bdead            = numeric(0)
                                  , balive           = numeric(0)
                                  , lai              = numeric(0)
                                  , stringsAsFactors = FALSE
                                  )#end data.frame
               pssnow = data.frame( time             = syear
                                  , patch            = pname
                                  , trk              = dist.type - 1
                                  , age              = dist.age
                                  , area             = 1.0
                                  , water            = 0.
                                  , fsc              = necro$fsc
                                  , stsc             = necro$stsc + acd.necro
                                  , stsl             = necro$stsl + acd.necro
                                  , ssc              = necro$ssc
                                  , lai              = 0.
                                  , isn              = necro$isn
                                  , fsn              = necro$fsn
                                  , msc              = necro$msc
                                  , gpp              = 0.
                                  , rh               = 0.
                                  , stringsAsFactors = FALSE
                                  )#end data.frame
            }#end if (nidx > 0)
            #------------------------------------------------------------------------------#



            #----- Check for degenerate solution. -----------------------------------------#
            suspicious = suspicious || (pssnow$lai %gt% 12)
            if (suspicious){
               if (nidx > 0){
                  #----- Aggregate information of the suspicious patch. -------------------#
                  css.agf    = pft$agf.bs[cssnow$pft]
                  css.sla    = pft$SLA[cssnow$pft]
                  css.bleaf  = with(cssnow,size2bl(dbh=dbh,hgt=hite,sla=css.sla,ipft=pft))
                  css.bsap   = pft$qsw   [cssnow$pft] * cssnow$hite * css.bleaf
                  css.bbark  = pft$qbark [cssnow$pft] * cssnow$hite * css.bleaf
                  css.agb    = ( css.bleaf
                               + css.agf * ( css.bsap + css.bbark + cssnow$bdead )
                               )#end css.agb
                  css.bsa    = 0.25 * pi * cssnow$dbh * cssnow$dbh
                  css.use    = cssnow$dbh %ge% 10.0
                  npl.show   = 10000. * sum(cssnow$n[css.use])
                  lai.show   = sum(cssnow$lai)
                  agb.show   = sum(cssnow$n[css.use]*css.agb[css.use])
                  bsa.show   = sum(cssnow$n[css.use]*css.bsa[css.use])
                  #------------------------------------------------------------------------#
               }else{
                  #----- Empty patch. -----------------------------------------------------#
                  npl.show   = 0.
                  lai.show   = 0.
                  agb.show   = 0.
                  bsa.show   = 0.
                  #------------------------------------------------------------------------#
               }#end if (nidx > 0)
               cat0(" ")
               cat0(" ")
               cat0(" ")
               cat0(" ")
               cat0(" ")
               cat0(" ")
               cat0("-------------------------------------------------------------------")
               cat0("       SUSPICIOUS PATCH! SUSPICIOUS PATCH! SUSPICIOUS PATCH!       ")
               cat0("-------------------------------------------------------------------")
               cat0("  Patch information: "                                              )
               cat0("  + Name                      -- ",pssnow$patch                     )
               cat0("  + Disturbance type          -- ",dist.type                        )
               cat0("  + Multiplication factor     -- ",f.net                            )
               cat0("  + Age                       -- ",dist.age                         )
               cat0("  + Stem density [DBH > 10cm] -- ",npl.show                         )
               cat0("  + Leaf area index           -- ",lai.show                         )
               cat0("  + AG Biomass   [DBH > 10cm] -- ",agb.show                         )
               cat0("  + Basal Area   [DBH > 10cm] -- ",bsa.show                         )
               cat0("-------------------------------------------------------------------")
               cat0(" ")
               cat0(" ")
               cat0(" ")
               cat0(" ")
               cat0(" ")
               cat0(" ")
            }#end if
            #------------------------------------------------------------------------------#
         }#end if (cstep %in% "calibration")
         #---------------------------------------------------------------------------------#
      }#end for (cstep in c("calibration","prediction"))
      #------------------------------------------------------------------------------------#
   }#end if (use.lookup)
   #---------------------------------------------------------------------------------------#



   #----- Return PSS and CSS components. --------------------------------------------------#
   ans = list(pss = pssnow, css = cssnow)
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function ptcloud.2.patch
#==========================================================================================#
#==========================================================================================#
