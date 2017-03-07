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
                            , bsa.goal
                            , agb.goal
                            , acd.goal
                            , npl.goal
                            , w.bsa
                            , w.agb
                            , w.npl
                            , lookup
                            , necro
                            , dist.type
                            , dist.age
                            , syear
                            , tall.at.zh    = TRUE
                            , use.intensity = FALSE
                            , lai.pst       = 2.0
                            , pft.pst       = 1
                            , dbh0.min      = 5.0
                            , pft.def       = 3
                            , use.lookup    = FALSE
                            ){

   #------ Run MacArthur-Horn (1969) correction to profiles. ------------------------------#
   mh.now = macarthur.horn( pt.cloud = pt.cloud
                          , zo       = zo
                          , zh       = zh
                          , nz       = nz
                          , tall.at.zh = tall.at.zh
                          , use.intensity = use.intensity
                          )#end macarthur.horn 
   #---------------------------------------------------------------------------------------#

   #----- Swap order of profiles. ---------------------------------------------------------#
   nzprof     = nrow(mh.now)
   t2b        = rev(sequence(nzprof))
   mh.now     = mh.now[t2b,]
   zmcprof    = mh.now$z
   dzprof     = 0*zmcprof - median(diff(zmcprof))
   hgtprof    = zmcprof / pft$b1Mh[pft.def]
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #   Find the most similar profile based on the lookup table.                            #
   #---------------------------------------------------------------------------------------#
   ks.stat    = apply(X=lookup$cdf,MARGIN=1,FUN=max.abs.diff,y=mh.now$cdf)
   o          = which.min(ks.stat)
   dbh.cutoff = lookup$dbh.cutoff[o]
   pftprof    = lookup$pft[o,,]
   ladprof    = lookup$lad[o,]
   indiv      = lookup$idv
   #---------------------------------------------------------------------------------------#


   #----- Find necromass associated with standing dead trees. -----------------------------#
   acd.necro = max(0.,acd.goal-agb.goal) / pft$agf.bs[pft.def]
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   # 
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
                            , msn              = necro$msn
                            , fsn              = necro$fsn
                            , nep              = 0.
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
                            , msn              = necro$msn
                            , fsn              = necro$fsn
                            , nep              = 0.
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
      cfz = curve.features(x=mh.now$pdf,span=5)
      if (sum(c(cfz$max,cfz$iph)) == 0){
      
         #----- Reduce window to make peaks more likely. ----------------------------------#
         cfz = curve.features(x=mh.now$pdf,span=3)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     In case a narrow failed to capture peaks, use function peaks.               #
         #---------------------------------------------------------------------------------#
         if (sum(c(cfz$max,cfz$iph)) == 0){
            iterate = TRUE
            nspan   = 17
            while (iterate){
               nspan = nspan - 2
               pktry = peaks(series=mh.now$pdf,span=nspan)
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
            sel = hgtprof %>=% dbh2h(dbh=dbh.cutoff,ipft=pft.def)
         }else{
            sel = hgtprof %>=% dbh2h(dbh=dbh0.min  ,ipft=pft.def)
         }#end if (cstep %in% "calibration") 
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Split layers by cohort.                                                    #
         #---------------------------------------------------------------------------------#
         c.ulai = split(x=mh.now$lad[sel]*dzprof[sel] ,f=icoh[sel])
         c.uipk = split(x=pcoh[sel]                   ,f=icoh[sel])
         #---------------------------------------------------------------------------------#
      
      
         #---------------------------------------------------------------------------------#
         #    Sometimes the density function creates blocks with LAI=0.  Drop them as      #
         # they are meaningless.                                                           #
         #---------------------------------------------------------------------------------#
         ulai = sapply(X=c.ulai,FUN=sum      )
         uipk = sapply(X=c.uipk,FUN=commonest)
         keep = ulai %>% 0.
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
            fpft = matrix(data=0,nrow=1,ncol=length(pft.keys)+1)
         }#end if (length(fpft) > 0)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Append layer for tall grasses (set to zero by default)                      #
         #---------------------------------------------------------------------------------#
         ncoh = sum(keep) + 1
         #---------------------------------------------------------------------------------#


         #----- Append grasses to the pft list. -------------------------------------------#
         mypfts   = c(pft.pst,match(pft.keys,pft$key))
         npfts    = length(mypfts)
         #---------------------------------------------------------------------------------#


         #----- Rename fpft. --------------------------------------------------------------#
         colnames(fpft) = c(pft$key[mypfts])
         #---------------------------------------------------------------------------------#



         #----- Find uncalibrated properties. ---------------------------------------------#
         uhgt.bnd = pmin(uhgt,pft$hgt.max[pft.def])
         udbh     = ( h2dbh(h=uhgt.bnd,ipft=pft.def)
                    * sqrt(pmax(uhgt,pft$hgt.max[pft.def])/pft$hgt.max[pft.def])
                    )#end udbh
         unpl     = ulai / (pft$SLA[pft.def] * dbh2bl(dbh=udbh,ipft=pft.def))
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Find the uncalibrated biomass and basal area.                              #
         #---------------------------------------------------------------------------------#
         npl.pft    = rep(unpl           , times= npfts)      * fpft
         dbh.pft    = rep(udbh           , times= npfts) + 0. * fpft
         hgt.pft    = rep(uhgt.bnd       , times= npfts) + 0. * fpft  
         ipft.pft   = rep(mypfts         , each = ncoh ) + 0L * fpft
         wdns.pft   = rep(pft$rho[mypfts], each = ncoh ) + 0. * fpft
         bleaf.pft  = dbh2bl(dbh=dbh.pft,ipft=ipft.pft)
         broot.pft  = rep(pft$qroot[mypfts],each = ncoh) * bleaf.pft
         bsw.pft    = rep(pft$qsw  [mypfts],each = ncoh) * hgt.pft   * bleaf.pft
         balive.pft = bleaf.pft + broot.pft + bsw.pft
         bdead.pft  = dbh2bd(dbh=dbh.pft,ipft=ipft.pft)
         agb.pft    = npl.pft * agb.SL(dbh=dbh.pft,height=hgt.pft,wdens=wdns.pft)
         bsa.pft    = npl.pft * pi/4 * dbh.pft^2
         lai.pft    = npl.pft * rep(pft$SLA[mypfts], each = ncoh ) * bleaf.pft
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #    Different steps depending on whether this is the first or second run.        #
         #---------------------------------------------------------------------------------#
         if (cstep %in% "calibration"){
            #------------------------------------------------------------------------------#
            #     Find calibration factors for pft and bsa, based on cross-validation      #
            # algorithms.                                                                  #
            #------------------------------------------------------------------------------#
            #----- Stem density.  Check for singularities before finding the global f. ----#
            if (sum(npl.pft) == 0){
               f.npl = 1.
               w.npl = 0.
            }else{
               f.npl = npl.goal / sum(npl.pft)
            }#end if (sum(npl.pft) == 0)
            #----- Basal area.  Check for singularities before finding the global f. ------#
            if (sum(bsa.pft) == 0){
               f.bsa = 1.
               w.bsa = 0.
            }else{
               f.bsa = bsa.goal / sum(bsa.pft)
            }#end if (sum(bsa.pft) == 0)
            #----- Biomass.  Check for singularities before finding the global f. ---------#
            if (sum(agb.pft) == 0){
               f.agb = 1.
               w.agb = 0.
            }else{
               f.agb = agb.goal / sum(agb.pft)
            }#end if (sum(agb.pft) == 0)
            #----- Net, check for singularities. ------------------------------------------#
            if ( (w.npl+w.bsa+w.agb) == 0){
               f.net = 1.
            }else{
               s.bsa = f.bsa^2
               s.agb = f.agb^2
               s.npl = f.npl^2
               f.net = ( f.npl * f.bsa * f.agb
                       * (w.npl*f.bsa*f.agb+f.npl*w.bsa*f.agb+f.npl*f.bsa*w.agb)
                       / (w.npl*s.bsa*s.agb+s.npl*w.bsa*s.agb+s.npl*s.bsa*w.agb)
                       )#end f.net
            }#end if( (w.npl+w.bsa+w.agb) == 0)
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
                                  , msn              = necro$msn
                                  , fsn              = necro$fsn
                                  , nep              = 0.
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
                                  , msn              = necro$msn
                                  , fsn              = necro$fsn
                                  , nep              = 0.
                                  , gpp              = 0.
                                  , rh               = 0.
                                  , stringsAsFactors = FALSE
                                  )#end data.frame
            }#end if (nidx > 0)
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
