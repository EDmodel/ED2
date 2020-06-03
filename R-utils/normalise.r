





#==========================================================================================#
#==========================================================================================#
#     This function normalises a vector.  Results are just the normalised values, mean and #
# standard deviation are in attributes.                                                    #
#------------------------------------------------------------------------------------------#
normalise <<- function(x,mu,sigma,distr=c("best","normal","lognormal","skewnormal")){

   #----- Standardise distribution. -------------------------------------------------------#
   distr       = match.arg(distr)
   #---------------------------------------------------------------------------------------#


   #----- Exclude missing points. ---------------------------------------------------------#
   xfit        = x[is.finite(x)]
   nxfit       = length(xfit)
   lntry       = all(xfit %>% 0)
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #     Issue an error message in case distr should be lognormal but there are negative   #
   # numbers.                                                                              #
   #---------------------------------------------------------------------------------------#
   if ( (distr %in% "lognormal") && (! lntry) ){
      stop("Log-normal distribution was requested, but data have non-positive numbers.")
   }#end if ( (distr %in% "lognormal") && (! lntry) )
   #---------------------------------------------------------------------------------------#


   #----- Use Gaussian distribution in case mu and sigma are provided. --------------------#
   force.gauss = ( distr %in% "normal"    ) || ( (! missing(mu)) && (! missing(sigma)) )
   force.skewn = ( distr %in% "skewnormal")
   force.lnorm = ( distr %in% "lognormal" )
   #---------------------------------------------------------------------------------------#


   if (nxfit == 0){
      #------ Do nothing. It should return an empty vector. -------------------------------#
      use.sn   = FALSE
      use.norm = FALSE
      #------------------------------------------------------------------------------------#
   }else if (force.gauss){
      #------ Make sure that both mu and sigma are defined. -------------------------------#
      if (missing(mu   )) mu    = mean(xfit)
      if (missing(sigma)) sigma = sd  (xfit)
      #------------------------------------------------------------------------------------#


      #------ Fit normal distribution. ----------------------------------------------------#
      lnlike    = -0.5*nxfit*(log(2*pi) + log(sigma^2)) - 0.5*sum((xfit-mu)^2)/sigma^2
      normfit   = list( estimate = c(mean=mu,sd=sigma)
                      , loglik   = lnlike
                      , bic      = log(nxfit)*2 - 2. * lnlike
                      )#end list
      use.sn    = FALSE
      use.norm  = TRUE
      #------------------------------------------------------------------------------------#
   }else{
      #------ Fit normal distribution. ----------------------------------------------------#
      normfit     = fitdistr(x=xfit,densfun="normal")
      normfit$bic = log(nxfit)*2 - 2. * normfit$loglik
      #------------------------------------------------------------------------------------#



      #------ Find parameters for skew-normal distribution (first guess). -----------------#
      sn3       = sn.stats(x)
      xi        = sn3["location"]
      omega     = sn3["scale"   ]
      alpha     = sn3["shape"   ]
      #------------------------------------------------------------------------------------#



      #------ Try to fit skew-normal distribution, then check whether it worked. ----------#
      snfit     = try( expr   = fitdistr( x       = xfit
                                        , densfun = dsn
                                        , start   = list(xi=xi,omega=omega,alpha=alpha)
                                        )#end fitdistr
                     , silent = TRUE
                     )#end try
      #------------------------------------------------------------------------------------#



      #------ Verify that snfit actually worked. ------------------------------------------#
      if ("try-error" %in% is(snfit)){
         #----- Fitting failed, copy normfit, but set likelihood to infinity to skip it. --#
         snfit     = normfit
         if (! force.skewn) snfit$bic = +Inf
         #---------------------------------------------------------------------------------#
      }else{
         #---------------------------------------------------------------------------------#
         #     Fit worked, find BIC.  Impose a penalty of 6 to select it only when the     #
         # improvement is really worth.                                                    #
         #---------------------------------------------------------------------------------#
         snfit$bic = log(nxfit)*3 - 2. * snfit$loglik + 6
         #---------------------------------------------------------------------------------#
      }#end if ("try-error" %in% is(snfit))
      #------------------------------------------------------------------------------------#



      #------ Fit log-normal distribution if it is possible. ------------------------------#
      if (lntry){
         lnormfit     = fitdistr(x=xfit,densfun="log-normal")
         lnormfit$bic = log(nxfit)*2 - 2. * lnormfit$loglik + 6
      }else{
         #----- Copy normfit, but set BIC to +infinity so it won't be chosen. -------------#
         lnormfit     = normfit
         lnormfit$bic = +Inf
         #---------------------------------------------------------------------------------#
      }#end if (lntry)
      #------------------------------------------------------------------------------------#



      #------ Choose the best distribution to normalise the data. -------------------------#
      use.sn   = force.skewn || ( (snfit$bic < normfit$bic) && (snfit$bic < lnormfit$bic) )
      use.norm = (! force.lnorm) && (! use.sn) && (normfit$bic < lnormfit$bic)
      #------------------------------------------------------------------------------------#



   }#end if (force.gauss)
   #---------------------------------------------------------------------------------------#



   #----- Normalise. ----------------------------------------------------------------------#
   if (nxfit == 0){
      #----- Make a dummy output, with nothing. -------------------------------------------#
      normal             = rep(NA_real_,times=length(x))
      attributes(normal) = list( n     = nxfit
                               , distr = NA_character_
                               )#end list
      #------------------------------------------------------------------------------------#
   }else if (use.sn){
      #------ Skew-normal distribution. ---------------------------------------------------#
      if ("xi" %in% names(snfit$estimate)){
         xi              = snfit$estimate["xi"   ]
         omega           = snfit$estimate["omega"]
         alpha           = snfit$estimate["alpha"]
      }else{
         xi              = snfit$estimate["mean" ]
         omega           = snfit$estimate["sd"   ]
         alpha           = 0.
      }#end if ("xi" %in names(snfit$estimate))
      normal             = skew2normal(x=x,location=xi,scale=omega,shape=alpha)
      normal             = ifelse(test=is.finite(normal),yes=normal,no=NA_real_)
      attributes(normal) = list( location = xi
                               , scale    = omega
                               , shape    = alpha
                               , n        = nxfit
                               , distr    = "sn"
                               )#end list
      #------------------------------------------------------------------------------------#
   }else if (use.norm){
      #------ Normal distribution. --------------------------------------------------------#
      mu                 = normfit$estimate["mean"]
      sigma              = normfit$estimate["sd"  ]
      normal             = (x - mu) / sigma
      normal             = ifelse(test=is.finite(normal),yes=normal,no=NA_real_)
      attributes(normal) = list( mean  = mu
                               , sdev  = sigma
                               , n     = nxfit
                               , distr = "normal"
                               )#end list
      #------------------------------------------------------------------------------------#
   }else{
      #------ Log-normal distribution. ----------------------------------------------------#
      lnmu               = lnormfit$estimate["meanlog"]
      lnsigma            = lnormfit$estimate["sdlog"  ]
      normal             = (log(x) - lnmu) / lnsigma
      normal             = ifelse(test=is.finite(normal),yes=normal,no=NA_real_)
      attributes(normal) = list( meanlog = lnmu
                               , sdevlog = lnsigma
                               , n       = nxfit
                               , distr   = "log-normal"
                               )#end list
      #------------------------------------------------------------------------------------#
   }#end if (use.sn)
   #---------------------------------------------------------------------------------------#

   return(normal)
}#end function normalise
#==========================================================================================#
#==========================================================================================#
