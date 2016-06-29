#==========================================================================================#
#==========================================================================================#
#      This function predicts the longwave radiation using the "CDO" scheme as presented   #
# at:                                                                                      #
#                                                                                          #
# Marthews, T. R., Y. Malhi, H. Iwata, 2012: Calculating downward longwave radiation under #
#     clear and cloudy conditions over a tropical lowland forest site: an evaluation of    #
#     model schemes for hourly data.  Theor. Appl. Climatol., 107, 461-477.                #
#                                                                                          #
# x is the set of parameters that are optimised for the dataset.                           #
#------------------------------------------------------------------------------------------#
rlong.in.alogit <<- function(x,datum){
   a = x
   #---------------------------------------------------------------------------------------#
   #     Estimate clear-sky emissivity.                                                    #
   #---------------------------------------------------------------------------------------#
   emiss.csky = a[1] + (1. - a[1]) * inv.logit(a[2] * (datum$atm.ph2o+a[3]))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Estimate the cloud cover.                                                         #
   #---------------------------------------------------------------------------------------#
   f.cloud  = 1. - inv.logit(a[4] * (datum$tau+a[5]))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Estimate the effective emissivity.                                                #
   #---------------------------------------------------------------------------------------#
   emiss.eff = emiss.csky * f.cloud
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Estimate cloud temperature.  This is really a really simple way of doing it, we   #
   # simply penalise the drier days by making the TLCL cooler.                             #
   #---------------------------------------------------------------------------------------#
   cld.tmp = (a[6] + (1. - a[6]) * datum$atm.rhv) * datum$atm.tlcl
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Estimate longwave radiation.                                                      #
   #---------------------------------------------------------------------------------------#
   rlong  = (                 emiss.csky  * stefan * datum$atm.tmp  ^ 4 
            + f.cloud * (1. - emiss.csky) * stefan * cld.tmp        ^ 4 )
   #---------------------------------------------------------------------------------------#

   ans = list(rlong.in=rlong,emiss.csky=emiss.csky,f.cloud=f.cloud,cld.tmp=cld.tmp)
   return(ans)
}#end predict
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#       This function finds the sum of the squares, which is the log likelihood if we      #
# asume the errors to be independent and normally distributed (a big assumption).          #
#==========================================================================================#
#==========================================================================================#
rlong.in.alogit.support <<- function(x,sigma,datum){

   rlong.in.try = rlong.in.alogit(x,datum)$rlong.in
   residual     = rlong.in.try - datum$rlong.in
   chi.square   = sum((residual/sigma)^2,na.rm=TRUE)

   support    = - chi.square
   return(support)
}#end function rlong.in.mmi.lnlike
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#       This function optimises the longwave radiation.  To make it simple, we use only    #
# the parametrisations that work day and night as shown in Table 2 of:                     #
#                                                                                          #
# Marthews, T. R., Y. Malhi, H. Iwata, 2012: Calculating downward longwave radiation under #
#     clear and cloudy conditions over a tropical lowland forest site: an evaluation of    #
#     model schemes for hourly data.  Theor. Appl. Climatol., 107, 461-477.                #
#                                                                                          #
#==========================================================================================#
#==========================================================================================#
rlong.in.alogit.optim <<- function(datum,keep.day=TRUE,keep.ngt=TRUE,verbose=FALSE){


   #---------------------------------------------------------------------------------------#
   #     Estimate some auxiliary thermodynamic and radiation properties.                   #
   #---------------------------------------------------------------------------------------#
   #----- Precipitable water. -------------------------------------------------------------#
   cat ("   - Precipitable water...","\n")
   datum$atm.ph2o = 4.65 * datum$atm.pvap / datum$atm.tmp
   #----- LCL (assumed to be the cloud base). ---------------------------------------------#
   cat ("   - LCL temperature (cloud base)...","\n")
   atm.theta      = datum$atm.tmp * (p00 / datum$atm.prss)^rocp
   lcl            = lcl.il(thil=atm.theta,pres=datum$atm.prss,temp=datum$atm.tmp
                          ,hum=datum$atm.pvap,type.hum="pvap")
   datum$atm.tlcl = lcl$temp
   #----- Tau, defined in terms of potential radiation at canopy. -------------------------#
   cat ("   - Fraction of shortwave (TOA) that reaches the ground ...","\n")
   datum$tau  = datum$rshort.in / datum$rshort.pot
   cat ("     * Gap-fill night time tau using linear interpolation...","\n")
   datum$tau[ ! datum$highsun ] = NA
   datum$tau  = na.fill(na.approx(datum$tau,na.rm=FALSE),"extend")
   #----- Split the time series into diel and season indices. -----------------------------#
   cat ("   - Splitting time series into hours and seasons...","\n")
   diel       = datum$hour
   un.diel    = unique(diel)
   idx.diel   = match(diel,un.diel)
   mon2season = c( 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 1)
   season     = mon2season[datum$month]
   un.season  = sort(unique(season))
   idx.season = match(season,un.season)
   #----- Find the standard deviation for weighted least squares. -------------------------#
   cat ("   - S.D. of incoming longwave, by hour and season ...","\n")
   sigma.hrse = tapply(X=datum$rlong.in,INDEX=list(diel,season),FUN=sd,na.rm=TRUE)
   idx        = cbind(idx.diel,idx.season)
   sigma      = sigma.hrse[idx]
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #     Data selection and total number of parameters.  Here we must decide whether we    #
   # are optimising daytime or nighttime.                                                  #
   #---------------------------------------------------------------------------------------#
   use   = ( is.finite(datum$rlong.in)  & is.finite(datum$atm.ph2o )
           & is.finite(datum$atm.tlcl)  & is.finite(datum$rshort.in)
           & is.finite(datum$atm.rhv)   )
   if (keep.day & keep.ngt){
      pre   = rep(TRUE,times=length(datum$rlong.in))
   }else if (keep.day){
      use   = use  & datum$highsun
      pre   = datum$highsun
   }else if (keep.ngt){
      use   = use & (! datum$highsun)
      pre   = ! datum$highsun
   }else{
      stop(" At least one between keep.day and keep.ngt must be TRUE...")
   }#end if
   n.use = sum(use)
   x.1st = c(0.7,0.1,-40.0,6.0,-0.6,0.6)
   n.par = 6
   #---------------------------------------------------------------------------------------#




   #----- 1st. guess, use Marthews et al. (2012) numbers for CDO. -------------------------#
   cat ("   - Run the optimiser...","\n")
   opt  = optim( par         = x.1st
               , fn          = rlong.in.alogit.support
               , sigma       = sigma[use]
               , datum       = datum[use,]
               , control     = list( trace   = verbose
                                   , fnscale = -1
                                   , maxit   = 20000
                                   , REPORT  = 1
                                   )#end list
               , hessian     = TRUE
               )#end optim
   if (opt$convergence != 0)  stop (" Solution of rlong.in.mmi.optim didn't converge...")
   #---------------------------------------------------------------------------------------#

   ans               = list()
   ans$hessian       = opt$hessian
   ans$df            = n.use - n.par
   ans$coefficients  = opt$par
   ans$std.err       = sqrt(diag(solve(-ans$hessian)))
   ans$t.value       = ans$coefficients / ans$std.err
   ans$p.value       = 2.0 * pt(-abs(ans$t.value),df=ans$df)
   ans$first.guess   = x.1st
   ans$support       = opt$value
   #----- Save the fitted values and the residuals. ---------------------------------------#
   myfit             = rlong.in.alogit(x=ans$coefficients,datum=datum)
   ans$fitted.values = myfit$rlong.in[pre]
   ans$emiss.csky    = myfit$emiss.csky ; ans$emiss.csky    [! use] = NA
   ans$f.cloud       = myfit$f.cloud    ; ans$f.cloud       [! use] = NA
   ans$cld.tmp       = myfit$cld.tmp    ; ans$cld.tmp       [! use] = NA
   ans$atm.ph2o      = datum$atm.ph2o   ; ans$atm.ph2o      [! use] = NA
   ans$atm.tlcl      = datum$atm.tlcl   ; ans$atm.tlcl      [! use] = NA
   ans$atm.tquant    = datum$atm.tquant ; ans$atm.tquant    [! use] = NA
   ans$tau           = datum$tau        ; ans$tau           [! use] = NA
   ans$lsq.weights   = 1. / sigma^2     ; ans$lsq.weights   [! use] = 0
   ans$residuals     = ans$fitted.values[use] - datum$rlong.in[use]

   #---------------------------------------------------------------------------------------#
   #     Estimate the global standard error and R2.                                        #
   #---------------------------------------------------------------------------------------#
   ss.err       = sum(ans$residuals^2)
   df.err       = ans$df
   mean.y       = mean(datum$rlong.in[use])
   ss.tot       = sum((datum$rlong.in[use]-mean.y)^2)
   df.tot       = n.use - 1
   ans$r.square = 1.0 - ss.err * df.tot / ( ss.tot * df.err )
   ans$sigma    = sqrt(ss.err / ans$df)
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end function rlong.in.mmi.lnlike
#==========================================================================================#
#==========================================================================================#
