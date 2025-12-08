#==========================================================================================#
#==========================================================================================#
#      This function applies the bias correction based on Cannon et al. (2015), assuming   #
# skew-normal distribution.  The default method is the 'QDM', which preserves the ratio    #
# and does not exaggerate much the occurrence of extremes.                                 #
#                                                                                          #
# References:                                                                              #
#                                                                                          #
# A. J. Cannon, S. R. Sobie, and T. Q. Murdock. Bias correction of GCM precipitation by    #
#     quantile mapping: How well do methods preserve changes in quantiles and extremes?    #
#     J. Climate, 28(17):6938-6959, Sep 2015. doi:10.1175/JCLI-D-14-00754.1.               #
#                                                                                          #
# S. Hempel, K. Frieler, L. Warszawski, J. Schewe, and F. Piontek. A trend-preserving bias #
#     correction --- the ISI-MIP approach. Earth Syst. Dynam., 4(2):219-236, Jul 2013.     #
#     doi:10.5194/esd-4-219-2013.                                                          #
#                                                                                          #
# ---------------------------------------------------------------------------------------- #
# Inputs                                                                                   #
# ---------------------------------------------------------------------------------------- #
# x.mp  - Data from the model projection (i.e. data to be corrected)                       #
# x.mh  - Modelled data from the historical period                                         #
# x.oh  - Observed data from the historical period                                         #
# distr - Which distribution to fit. Possible values are:                                  #
#         empirical: use ECDF and spline to interpolate between values                     #
#         sn (or skewnormal): use skew-normal distribution                                 #
#         normal (or Gaussian): use normal distribution                                    #
# method - Which method to apply. Possible values are:                                     #
#          QM  - simple quantile mapping, which relies on distributions from historical    #
#                period only.                                                              #
#          DQM - detrended quantile mapping, which takes into account the change in the    #
#                mean values between modelled historical and modelled projected data.      #
#          QDM - quantile delta mapping.  The method developed by Cannon et al. (2015),    #
#                which preserves relative changes in quantiles in the modelled             #
#                projections.                                                              #
#          H13 - Hempel et al. (2013) bias-correction approach.  Simple scaling factor     #
#                based on the ratio of the means.                                          #
# xlwr   - Lowest value for x (at which CDF should be zero). Used only for empirical       #
#          distributions.                                                                  #
# xupr   - Uppermost value of x (at which CDF should be one). Used only for empirical      #
#          distributions.                                                                  #
# na.rm  - Should NA's be ignored when deriving statistics?  Used when                     #
#------------------------------------------------------------------------------------------#
cannon.biascorr <<- function( x.mp
                            , x.mh
                            , x.oh
                            , distr  = c("empirical","sn","skewnormal","normal","Gaussian")
                            , method = c("QDM","DQM","QM","H13")
                            , xlwr   = 0
                            , xupr   = 1.1*max(c(x.mp,x.mh,x.oh),na.rm=TRUE)
                            , cmax   = 10.
                            , na.rm  = TRUE
                            ){

   #----- Find method to use. -------------------------------------------------------------#
   distr  = match.arg(distr )
   method = match.arg(method)
   #---------------------------------------------------------------------------------------#



   #----- Crash in case this is an old sn version. ----------------------------------------#
   if ((distr %in% c("sn","skewnormal")) && (sn.version == 0)){
      stop( paste0( " Function cannon.biascorr was called with distr set to \"sn\"\n"
                  , " (\"skewnormal\").  Your \"sn\" package is too old, update it\n"
                  , " to a version > 1.0."
                  )#end paste0
          )#end stop
   }#end if (sn.version == 0)
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #      Decide which distribution to use.                                                #
   #---------------------------------------------------------------------------------------#
   if (method %in% "H13"){
      #----- Hempel et al. (2013) method.  No distribution is needed. ---------------------#
      x.oh.bar = mean(x.oh,na.rm = na.rm)
      x.mh.bar = mean(x.mh,na.rm = na.rm)
      if (x.oh.bar == 0.){
         corr = 0.
      }else{
         corr = x.oh.bar / pmax(x.oh.bar / cmax,x.mh.bar)
      }#end if 
      ans = corr * x.mp
      #------------------------------------------------------------------------------------#
   }else if (distr %in% "empirical"){
      #------------------------------------------------------------------------------------#
      #     Make sure that xlwr and xupr are outside the bounds of x.oh, x.mh, and x.mp.   #
      #------------------------------------------------------------------------------------#
      fine.oh = is.na(x.oh) | ( x.oh %wr% c(xlwr,xupr) )
      fine.mh = is.na(x.mh) | ( x.mh %wr% c(xlwr,xupr) )
      fine.mp = is.na(x.mp) | ( x.mp %wr% c(xlwr,xupr) )
      if (! all(c(fine.oh,fine.mh,fine.mp)) ){
         #------ Inform the error. --------------------------------------------------------#
         cat0("--------------------------------------------------------------------------")
         cat0(" FATAL ERROR!"                                                             )
         cat0("--------------------------------------------------------------------------")
         cat0(" - XLWR          = ",xlwr                                                  )
         cat0(" - XUPR          = ",xupr                                                  )
         cat0(" - Range of x.oh = ",min(x.oh,na.rm=TRUE)," ",max(x.oh,na.rm=TRUE)         )
         cat0(" - Range of x.mh = ",min(x.mh,na.rm=TRUE)," ",max(x.mh,na.rm=TRUE)         )
         cat0(" - Range of x.mp = ",min(x.mp,na.rm=TRUE)," ",max(x.mp,na.rm=TRUE)         )
         cat0(" >>> All 'x' variables must be within range of 'xlwr' and 'xupr'!"         )
         cat0("--------------------------------------------------------------------------")
         stop(" Invalid range using emprirical distribution.")
         #---------------------------------------------------------------------------------#
      }#end if (! all(c(fine.oh,fine.mh,fine.mp)) )
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Define the empirical functions for observed and modelled (historical).        #
      #------------------------------------------------------------------------------------#
      #----- Observed (historical) data. --------------------------------------------------#
      ecdf.oh = ecdf(sort(c(x.oh,xupr)))
      q.oh    = c(xlwr,x.oh,xupr)
      tau.oh  = ecdf.oh(q.oh)
      pfun.oh = splinefun(x=q.oh  ,y=tau.oh,method="monoH.FC")
      qfun.oh = splinefun(x=tau.oh,y=q.oh  ,method="monoH.FC")
      #----- Modelled (historical) data. --------------------------------------------------#
      ecdf.mh = ecdf(sort(c(x.mh,xupr)))
      q.mh    = c(xlwr,x.mh,xupr)
      tau.mh  = ecdf.mh(q.mh)
      pfun.mh = splinefun(x=q.mh  ,y=tau.mh,method="monoH.FC")
      qfun.mh = splinefun(x=tau.mh,y=q.mh  ,method="monoH.FC")
      #------------------------------------------------------------------------------------#

      #------------------------------------------------------------------------------------#
      #      Empirical distribution.  Pick the method.                                     #
      #------------------------------------------------------------------------------------#
      if (method %in% "QM"){
         #----- Quantile mapping. ---------------------------------------------------------#
         tau = pfun.mh(x.mp)
         ans = qfun.oh(tau )
         #---------------------------------------------------------------------------------#
      }else if (method %in% "DQM"){
         #----- Detrended quantile mapping. -----------------------------------------------#
         mhomp = mean(x.mh,na.rm=na.rm) / mean(x.mp,na.rm=na.rm)
         tau   = pfun.mh(mhomp*x.mp)
         ans   = qfun.oh(tau       ) / mhomp
         #---------------------------------------------------------------------------------#
      }else if (method %in% "QDM"){
         #---------------------------------------------------------------------------------#
         #     Quantile delta mapping.                                                     #
         #---------------------------------------------------------------------------------#


         #----- Find functions for modelled (projection). ---------------------------------#
         ecdf.mp = ecdf(sort(c(x.mp,xupr)))
         q.mp    = c(xlwr,x.mp,xupr)
         tau.mp  = ecdf.mp(q.mp)
         pfun.mp = splinefun(x=q.mp  ,y=tau.mp,method="monoH.FC")
         qfun.mp = splinefun(x=tau.mp,y=q.mp  ,method="monoH.FC")
         #---------------------------------------------------------------------------------#


         #----- Find the projection CDF. --------------------------------------------------#
         tau.mp  = pfun.mp(x.mp)
         delta.m = x.mp    / qfun.mh(tau.mp)
         ans     = delta.m * qfun.oh(tau.mp)
         #---------------------------------------------------------------------------------#
      }#end if (method %in% "QM")
      #------------------------------------------------------------------------------------#
   }else if (distr %in% c("sn","skewnormal")){

      #----- Find skew-normal statistics for both observed and modelled (historical). -----#
      stats.oh = sn.stats(x.oh,na.rm=na.rm)[sequence(3)]
      stats.mh = sn.stats(x.mh,na.rm=na.rm)[sequence(3)]
      #------------------------------------------------------------------------------------#

      #------------------------------------------------------------------------------------#
      #      Skew-normal distribution.  Decide which method to use.                        #
      #------------------------------------------------------------------------------------#
      if (method %in% "QM"){
         #----- Quantile mapping. ---------------------------------------------------------#
         tau = psn(x=x.mp,dp=stats.mh,engine="T.Owen")
         ans = qsn(p=tau ,dp=stats.oh,engine="T.Owen",solver="RFB")
         #---------------------------------------------------------------------------------#
      }else if (method %in% "DQM"){
         #----- Detrended quantile mapping. -----------------------------------------------#
         mhomp = mean(x.mh,na.rm=na.rm) / mean(x.mp,na.rm=na.rm)
         tau   = psn(x=mhomp*x.mp,dp=stats.mh,engine="T.Owen")
         ans   = qsn(p=tau       ,dp=stats.oh,engine="T.Owen",solver="RFB") / mhomp
         #---------------------------------------------------------------------------------#

      }else if (method %in% "QDM"){
         #---------------------------------------------------------------------------------#
         #     Quantile delta mapping.                                                     #
         #---------------------------------------------------------------------------------#


         #----- Find statistics for projection. -------------------------------------------#
         stats.mp = sn.stats(x.mp,na.rm=na.rm)[sequence(3)]
         #---------------------------------------------------------------------------------#


         #----- Find the projection CDF. --------------------------------------------------#
         tau.mp  = psn(x=x.mp,dp=stats.mp,engine="T.Owen")
         delta.m = x.mp    / qsn(p=tau.mp,dp=stats.mh,engine="T.Owen",solver="RFB")
         ans     = delta.m * qsn(p=tau.mp,dp=stats.oh,engine="T.Owen",solver="RFB")
         #---------------------------------------------------------------------------------#
      }#end if (method %in% "QM")
      #------------------------------------------------------------------------------------#
   }else if (distr %in% c("normal","Gaussian")){

      #----- Find Gaussian statistics for both observed and modelled (historical). --------#
      stats.oh = c(mu=mean(x.oh,na.rm=na.rm),sigma=sd(x.oh,na.rm=na.rm))
      stats.mh = c(mu=mean(x.mh,na.rm=na.rm),sigma=sd(x.mh,na.rm=na.rm))
      #------------------------------------------------------------------------------------#

      #------------------------------------------------------------------------------------#
      #      Skew-normal distribution.  Decide which method to use.                        #
      #------------------------------------------------------------------------------------#
      if (method %in% "QM"){
         #----- Quantile mapping. ---------------------------------------------------------#
         tau = pnorm(x=x.mp,mean=stats.mh[1],sd=stats.mh[2])
         ans = qnorm(p=tau ,mean=stats.oh[1],sd=stats.mh[2])
         #---------------------------------------------------------------------------------#
      }else if (method %in% "DQM"){
         #----- Detrended quantile mapping. -----------------------------------------------#
         mhomp = mean(x.mh,na.rm=na.rm) / mean(x.mp,na.rm=na.rm)
         tau   = pnorm(x=mhomp*x.mp,mean=stats.mh[1],sd=stats.mh[2])
         ans   = qnorm(p=tau       ,mean=stats.oh[1],sd=stats.oh[2]) / mhomp
         #---------------------------------------------------------------------------------#

      }else if (method %in% "QDM"){
         #---------------------------------------------------------------------------------#
         #     Quantile delta mapping.                                                     #
         #---------------------------------------------------------------------------------#


         #----- Find statistics for projection. -------------------------------------------#
         stats.mp = c(mu=mean(x.mp,na.rm=na.rm),sigma=sd(x.mp,na.rm=na.rm))
         #---------------------------------------------------------------------------------#


         #----- Find the projection CDF. --------------------------------------------------#
         tau.mp  = pnorm(x=x.mp,mean=stats.mp[1],sd=stats.mp[2])
         delta.m = x.mp    / qnorm(p=tau.mp,mean=stats.mh[1],sd=stats.mh[2])
         ans     = delta.m * qnorm(p=tau.mp,mean=stats.oh[1],sd=stats.oh[2])
         #---------------------------------------------------------------------------------#
      }#end if (method %in% "QM")
      #------------------------------------------------------------------------------------#
   }#end if (distr %in% "empirical")
   #---------------------------------------------------------------------------------------#



   #----- Answer is the bias-corrected model projection. ----------------------------------#
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end cannon.sn.biascorr
#==========================================================================================#
#==========================================================================================#
