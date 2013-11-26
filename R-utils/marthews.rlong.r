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
rlong.in.mmi.predict <<- function(x,scheme,datum){

   #---------------------------------------------------------------------------------------#
   #     Estimate emissivity, based on the chosen scheme.                                  #
   #---------------------------------------------------------------------------------------#
   if (scheme == "cij"){

      #----- Idso and Jackson (1969). -----------------------------------------------------#
      emiss.csky = NA * datum$atm.ph2o
      f.cloud    = NA * datum$tau
      emiss.eff  = 1.0 + x[1]  * exp(x[2] * (273 - datum$atm.tmp)^2)
      cld.tmp    = datum$atm.tlcl
      #------------------------------------------------------------------------------------#

   }else if (scheme == "cid"){

      #----- Idso (1981). -----------------------------------------------------------------#
      emiss.csky = NA * datum$atm.ph2o
      f.cloud    = NA * datum$tau
      emiss.eff  = x[1] + x[2] * datum$atm.pvap * exp (x[3] / datum$atm.tmp)
      cld.tmp    = datum$atm.tlcl
      #------------------------------------------------------------------------------------#

   }else if (scheme == "cdo"){

      #----- Dilley and O'Brien (1998), model B. ------------------------------------------#
      emiss.csky = NA * datum$atm.ph2o
      f.cloud    = NA * datum$tau
      emiss.eff  = ( (x[1] + x[2] * (datum$atm.tmp/t3ple)^6
                           + x[3] * sqrt(datum$atm.ph2o/25) )
                   / (stefan * datum$atm.tlcl^4 ) )
      cld.tmp    = datum$atm.tlcl
      #------------------------------------------------------------------------------------#

   }else if (scheme == "cmu"){

      #----- Monteith and Unsworth (2008). ------------------------------------------------#
      emiss.csky = NA * datum$atm.ph2o
      f.cloud    = NA * datum$tau
      emiss.eff  = ( x[1] + x[2] * stefan * datum$atm.tmp^4) / (stefan * datum$atm.tlcl^4)
      cld.tmp    = datum$atm.tlcl
      #------------------------------------------------------------------------------------#


   }else if (scheme == "amu"){

      #----- Clear sky, from Monteith and Unsworth (2008). --------------------------------#
      emiss.csky = ( x[1] + x[2] * stefan * datum$atm.tmp^4) / (stefan * datum$atm.tlcl^4)
      #------------------------------------------------------------------------------------#

      #----- Cloud cover from Black (1956). -----------------------------------------------#
      k.crit  = x[5] - x[4]^2 / (4. * x[3])
      use.kappa = pmin(k.crit,datum$kappa)
      f.cloud = ( x[4] - sqrt(x[4]^2 - 4. * x[3] * ( x[5] - use.kappa)) ) / (2. * x[3])
      f.cloud = pmax(0.,pmin(1,f.cloud))
      #------------------------------------------------------------------------------------#


      #----- Correction for overcast emissivity (Monteith and Unsworth 2008). -------------#
      oc.fac    = 1.0 - 4.0 * (datum$atm.tmp - datum$atm.tlcl) / datum$atm.tmp
      emiss.eff = (1.0 - oc.fac * f.cloud) * emiss.csky + oc.fac * f.cloud
      #------------------------------------------------------------------------------------#


      #----- Cloud temperature, assumed to be the LCL. ------------------------------------#
      cld.tmp    = datum$atm.tlcl
      #------------------------------------------------------------------------------------#
   }else if (scheme == "alm"){

      #----- Idso (1981). -----------------------------------------------------------------#
      emiss.csky = x[1] + x[2] * datum$atm.pvap * exp (x[3] / datum$atm.tmp)
      #------------------------------------------------------------------------------------#


      #----- Stockli et al. (2007), the LBA-MIP. ------------------------------------------#
      f.cloud    = (1.0 - datum$tau)^2
      emiss.eff  = emiss.csky * (1.0 + x[4] * f.cloud)
      #------------------------------------------------------------------------------------#


      #----- Cloud temperature, assumed to be the LCL. ------------------------------------#
      cld.tmp    = datum$atm.tlcl
      #------------------------------------------------------------------------------------#

   }else if (scheme == "adk"){

      #----- Dilley and O'Brien (1998), model B. ------------------------------------------#
      emiss.csky  = ( (x[1] + x[2] * (datum$atm.tmp/t3ple)^6
                    + x[3] * sqrt(datum$atm.ph2o/25) ) / (stefan * datum$atm.tlcl^4 ) )
      #------------------------------------------------------------------------------------#

      #----- Cloud cover by Kasten and Czeplak (1980). ------------------------------------#
      f.cloud  = (x[4] * (1.0 - datum$tau))^x[5]
      f.cloud  = pmax(0.0,pmin(1.0,f.cloud))
      #------------------------------------------------------------------------------------#

      #----- Dilley-Kimball scheme from Flerchinger et al. (2009). ------------------------#
      emiss.8z  = x[6] + x[7] * datum$atm.pvap^2 * exp (x[8]/datum$atm.tmp)
      f.8i      = x[9] + x[10] * datum$atm.tlcl + x[11] * datum$atm.tlcl^2
      emiss.eff = ( emiss.csky
                  + (1.0 - emiss.8z * (1 + x[12]*(1 - emiss.8z))) * f.cloud * f.8i )
      #------------------------------------------------------------------------------------#


      #----- Cloud temperature, assumed to be the LCL. ------------------------------------#
      cld.tmp    = datum$atm.tlcl
      #------------------------------------------------------------------------------------#
   }else if (scheme == "aml"){
      #------------------------------------------------------------------------------------#
      #    Not yet published method by Longo et al., based on Monteith and Unsworth (2008) #
      # and Marthews et al. (2012), but applying a relative humidity correction to the     #
      # cloud temperature, and using the inverse of logistic function to predict the clear #
      # sky emissivity and cloud cover (so the functions are always bounded).              #
      #------------------------------------------------------------------------------------#


      #---- Estimate clear-sky emissivity. ------------------------------------------------#
      emiss.csky = x[1] + (1. - x[1]) * inv.logit(x[2] * (datum$atm.ph2o+x[3]))
      #------------------------------------------------------------------------------------#



      #---- Estimate the cloud cover. -----------------------------------------------------#
      f.cloud  = 1. - inv.logit(x[4] * (datum$kappa+x[5]))
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Estimate cloud temperature.  This is really a really simple way of doing it,   #
      # we simply penalise the drier days by making the TLCL cooler.                       #
      #------------------------------------------------------------------------------------#
      cld.tmp = (x[6] + (1. - x[6]) * datum$atm.rhv) * datum$atm.tlcl
      #------------------------------------------------------------------------------------#



      #---- Estimate the effective emissivity. --------------------------------------------#
      emiss.eff = emiss.csky * (datum$atm.tmp/cld.tmp)^4 + f.cloud * (1.0 - emiss.csky)
      #------------------------------------------------------------------------------------#


   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Estimate longwave radiation.                                                      #
   #---------------------------------------------------------------------------------------#
   rlong.in  = emiss.eff * stefan * cld.tmp  ^ 4
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Update effective emissivity so it is based on air temperature.                    #
   #---------------------------------------------------------------------------------------#
   emiss.eff = emiss.eff * ( cld.tmp / datum$atm.tmp ) ^ 4
   #---------------------------------------------------------------------------------------#

   ans = list( rlong.in   = rlong.in
             , emiss.csky = emiss.csky
             , emiss.eff  = emiss.eff
             , f.cloud    = f.cloud
             , cld.tmp    = cld.tmp
             )#end list
   return(ans)
}#end predict
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#       This function finds the support assuming that the residuals are independent and    #
# normally distributed, which is quite a big assumption...                                 #
#------------------------------------------------------------------------------------------#
rlong.in.mmi.support <<- function(x,scheme,sigma,datum){

   #----- Solve long wave for this set of parameters. -------------------------------------#
   rlong.in.try  = rlong.in.mmi.predict(x,scheme,datum)$rlong.in
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Determine the residuals and find the standard deviation.  We do not find the mean #
   # residual (bias), because we want to find the parameters that are the closest to       #
   # unbiased distributions.                                                               #
   #---------------------------------------------------------------------------------------#
   residual      = rlong.in.try - datum$rlong.in
   sdev.residual = sd(residual,na.rm=TRUE)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Since we don't know the distribution of the residual a priori, we must include   #
   # the variance in the support, but force the mean to be 0.                              #
   #---------------------------------------------------------------------------------------#
   support = sum(dnorm(residual,mean=0,sd=sdev.residual,log=TRUE),na.rm=TRUE)
   #---------------------------------------------------------------------------------------#

   return(support)
}#end function rlong.in.mmi.support
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#       This function models the incoming longwave radiation using some of the approaches  #
# by:                                                                                      #
#                                                                                          #
# Marthews, T. R., Y. Malhi, H. Iwata, 2012: Calculating downward longwave radiation under #
#     clear and cloudy conditions over a tropical lowland forest site: an evaluation of    #
#     model schemes for hourly data.  Theor. Appl. Climatol., 107, 461-477.                #
#                                                                                          #
#------------------------------------------------------------------------------------------#
rlong.in.mmi.optim <<- function(datum,run.optim=FALSE,keep.day=TRUE,keep.ngt=TRUE
                               ,scheme=NULL,verbose=FALSE){




   #----- List with all schemes available. ------------------------------------------------#
   avail.schemes       = list()
   avail.schemes[[ 1]] = list( name  = "cij"
                             , desc  = "Idso and Jackson (1969)"
                             , n.par = 2
                             , x.def = c(-0.261, -7.77e-4)
                             )#end list
   avail.schemes[[ 2]] = list( name  = "cid"
                             , desc  = "Idso (1981)"
                             , n.par = 3
                             , x.def = c(0.7, 5.95e-7,1500.)
                             )#end list
   avail.schemes[[ 3]] = list( name  = "cdo"
                             , desc  = "Dilley and O'Brien (1998), model B"
                             , n.par = 3
                             , x.def = c(59.38, 113.7,96.96)
                             )#end.list
   avail.schemes[[ 4]] = list( name  = "cmu"
                             , desc  = "Clear sky - Monteith and Unsworth (2008)"
                             , n.par = 2
                             , x.def = c(-119.0,1.06)
                             )#end.list
   avail.schemes[[ 5]] = list( name  = "amu"
                             , desc  = "All skies - Monteith and Unsworth (2008)"
                             , n.par = 5
                             , x.def = c(-119.0,1.06,-0.458,0.34,0.803)
                             )#end.list
   avail.schemes[[ 6]] = list( name  = "alm"
                             , desc  = "All skies - Stockli et al. (2007), LBA-MIP"
                             , n.par = 4
                             , x.def = c(0.7,5.95e-7,1500.,0.3)
                             )#end.list
   avail.schemes[[ 7]] = list( name  = "adk"
                             , desc  = "All skies - Flerchinger et al. (2009), DK scheme"
                             , n.par = 12
                             , x.def = c(59.38,113.7,96.96,4./3.,1/3.4,0.24,2.98e-12,3000.
                                        ,-0.6732,6.24e-3,-9.14e-6,0.4)
                             )#end.list
   avail.schemes[[ 8]] = list( name  = "aml"
                             , desc  = "All skies - Longo et al. (2012)"
                             , n.par = 6
                             , x.def = c(0.7161,0.02572,-46.63,9.205,-0.500,0.324)
                             )#end.list
   n.avail.schemes = length(avail.schemes)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Find out which schemes to use.                                                    #
   #---------------------------------------------------------------------------------------#
   name.schemes = unlist(sapply(X=avail.schemes,FUN=sapply,c,simplify=TRUE)[1,])
   if (is.null(scheme)){
      sel.schemes = sequence(n.avail.schemes)
   }else{
      sel.schemes  = match(scheme,name.schemes)
   }#end if
   #---------------------------------------------------------------------------------------#



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
   #----- Tau, defined here as the ratio between actual and clear sky. --------------------#
   cat ("   - Tau, rshort/rshort(clear)...","\n")
   datum$tau                  = datum$rshort.in / pmax(datum$rshort.in,datum$rshort.pot)
   datum$tau[! datum$highsun] = NA
   datum$tau                  = na.fill(na.approx(datum$tau,na.rm=FALSE),"extend")
   datum$tau                  = pmin(1.,pmax(0.,datum$tau))
   #----- Kappa, defined here as the ratio between actual radiation and TOA shortwave. ----#
   cat ("   - Kappa, rshort/rshort(TOA)...","\n")
   datum$kappa                  = ( datum$rshort.in 
                                  / pmax(datum$rshort.in,solar * datum$cosz) )
   datum$kappa[! datum$highsun] = NA
   datum$kappa                  = na.fill(na.approx(datum$kappa,na.rm=FALSE),"extend")
   datum$kappa                  = pmin(1.,pmax(0.,datum$kappa))
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
   #      Test only with the actual data.  We only try to estimate the goodness of the fit #
   # in case longwave radiation has some actual data.                                      #
   #---------------------------------------------------------------------------------------#
   if (any(datum$gfflg.rlong.in == 0 & is.finite(datum$rlong.in),na.rm=TRUE)){
      fit.use   = ( datum$gfflg.rlong.in == 0  & is.finite(datum$rlong.in)
                  & is.finite(datum$atm.prss)  & is.finite(datum$atm.tmp)
                  & is.finite(datum$atm.pvap)  & is.finite(datum$rshort.in) )
   }else{
      if (run.optim){
         warning(" No longwave radiation, skip fitting and use defaults...")
      }#end if
      fit.use   = ( is.finite(datum$atm.prss)  & is.finite(datum$atm.tmp)
                  & is.finite(datum$atm.pvap)  & is.finite(datum$rshort.in) )
      run.optim = FALSE
   }#end if
   if (keep.day & keep.ngt){
      pred.use = rep(TRUE,times=length(datum$rlong.in))
   }else if (keep.day){
      fit.use  = fit.use  & datum$highsun
      pred.use = datum$highsun
   }else if (keep.ngt){
      fit.use  = fit.use & (! datum$highsun)
      pred.use = ! datum$highsun
   }else{
      stop(" At least one between keep.day and keep.ngt must be TRUE...")
   }#end if
   n.fit.use = sum(fit.use)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Loop over all schemes that we should solve.                                       #
   #---------------------------------------------------------------------------------------#
   ans             = list()
   attributes(ans) = list( run.optim = run.optim
                         , keep.day  = keep.day
                         , keep.ngt  = keep.ngt
                         , scheme    = name.schemes[sel.schemes]
                         )#end list
   for (s in sel.schemes){
      this         = avail.schemes[[s]]

      #------------------------------------------------------------------------------------#
      #     Check whether to optimise the scheme or not.                                   #
      #------------------------------------------------------------------------------------#
      if (run.optim){
         cat ("   - Run the optimiser for ",this$desc,"...","\n")
         datum.fit = datum[fit.use,]
         opt   = optim( par         = this$x.def
                      , fn          = rlong.in.mmi.support
                      , datum       = datum[fit.use,]
                      , sigma       = sigma
                      , scheme      = this$name
                      , control     = list( trace   = verbose
                                          , fnscale = -1
                                          , maxit   = 30000
                                          , rel.tol = 1.e-9
                                          , REPORT  = 1
                                          )#end list
                      , hessian     = TRUE
                      )#end optim
         fine  = opt$convergence == 0
         if (! fine) cat ("     > Optimiser didn't converge, use default instead...","\n")
      }else{
         cat ("   - Check goodness of fit for ",this$desc,"...","\n")
         fine = FALSE
      }#end if
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Check whether the optimiser worked or not.  If it didn't, or if the user       #
      # preferred to use the default values, then use the default coefficients.            #
      #------------------------------------------------------------------------------------#
      if (! fine){
         this$convergence  = FALSE
         this$coefficients = this$x.def
         this$support      = rlong.in.mmi.support( x      = this$coefficients
                                                 , scheme = this$name
                                                 , sigma  = sigma
                                                 , datum  = datum[fit.use,]
                                                 )#end rlong.in.mmi.support
      }else{
         this$convergence  = TRUE
         this$coefficients = opt$par
         this$support      = opt$value
      }#end if
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      The remaining values are common for the convergence and non convergence.      #
      #------------------------------------------------------------------------------------#
      this$hessian      = hessian( func   = rlong.in.mmi.support
                                 , x      = this$coefficients
                                 , scheme = this$name
                                 , sigma  = sigma
                                 , datum  = datum[fit.use,]
                                 , method.args = list( eps = 1.e-4
                                                     , d   = 0.0001
                                                     , r   = 14
                                                     , v   = 2
                                                     , zero.tol=sqrt(.Machine$double.eps)
                                                     )
                                 )#end hessian
      this$df           = n.fit.use - this$n.par
      attempt           = try(solve(- this$hessian),silent=TRUE)
      if ("try-error" %in% is(attempt)){
         warning("     > Cannot solve the Information matrix...")
         this$std.err   = NA * this$coefficients
      }else{
         this$std.err   = sqrt(this$coefficients^2 * diag(solve(- this$hessian)))
      }#end if
      #------------------------------------------------------------------------------------#

      this$t.value       = this$coefficients / this$std.err
      this$p.value       = 2.0 * pt(-abs(this$t.value),df=this$df)
      myfit              = rlong.in.mmi.predict( x      = this$coefficients
                                               , scheme = this$name
                                               , datum  = datum
                                               )#end rlong.in.mmi.predict
      this$fitted.values = myfit$rlong.in
      this$emiss.csky    = myfit$emiss.csky
      this$emiss.eff     = myfit$emiss.eff
      this$f.cloud       = myfit$f.cloud
      this$cld.tmp       = myfit$cld.tmp
      this$atm.ph2o      = datum$atm.ph2o
      this$atm.tlcl      = datum$atm.tlcl
      this$tau           = datum$tau
      this$kappa         = datum$kappa
      this$lsq.weights   = 1. / sigma^2
      this$residuals     = this$fitted.values - datum$rlong.in
      this$fit.sel       = fit.use
      this$predict.sel   = pred.use
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    Eliminate residual values that were not used for fitting and fitted values for  #
      # points in which this fitting should not be applied.                                #
      #------------------------------------------------------------------------------------#
      this$fitted.values [ ! pred.use] = NA
      this$residuals     [ ! fit.use ] = NA
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Estimate the global standard error and R2.                                     #
      #------------------------------------------------------------------------------------#
      ss.err        = sum(this$residuals[fit.use]^2)
      df.err        = this$df
      mean.y        = mean(datum$rlong.in[fit.use])
      ss.tot        = sum((datum$rlong.in[fit.use]-mean.y)^2)
      df.tot        = n.fit.use - 1
      this$r.square = 1.0 - ss.err * df.tot / ( ss.tot * df.err )
      this$sigma    = sqrt(ss.err / this$df)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Save this result to the global list.                                           #
      #------------------------------------------------------------------------------------#
      ans[[this$name]] = this
      #------------------------------------------------------------------------------------#
   }#end for (s in 1:n.schemes)
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end function rlong.in.mmi.lnlike
#==========================================================================================#
#==========================================================================================#
