#==========================================================================================#
#==========================================================================================#
#     This function computes several metrics to test the goodness of fit.  The list is     #
# not complete, contributions are welcome.                                                 #
#                                                                                          #
#  INPUT:                                                                                  #
#  ~ x.mod          -- the modelled values of x                                            #
#  ~ x.obs          -- the observed values of x                                            #
#  ~ n.parameters   -- number of parameters (if the number of parameters is unknown, we    #
#                      don't assume any parameters).                                       #
#------------------------------------------------------------------------------------------#
test.goodness <<- function(x.mod,x.obs,n.parameters=NULL){


   #---- Crash if the x.mod and x.obs don't have the same size and class. -----------------#
   dlength = (length(x.mod) - length(x.obs)) != 0
   dclass  = any(sort(is(x.mod)) != sort(is(x.obs)))
   if (dlength || dclass){
      stop (" x.mod and x.obs must have the same size and class","\n")
   }#end if
   #---------------------------------------------------------------------------------------#



   #---- Find some general variables that will be used during this function. --------------#
   sel      = is.finite(x.mod) & is.finite(x.obs)
   n.ok     = sum(sel)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the number of degrees of freedom for both total and explained variance.      #
   #---------------------------------------------------------------------------------------#
   df.tot = n.ok - 1
   if (is.null(n.parameters)){
      df.err = n.ok - 1
   }else{
      df.err = n.ok - n.parameters - 1
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Check whether we have enough valid observations before proceeding.                #
   #---------------------------------------------------------------------------------------#
   if (df.err > 1){
      #------------------------------------------------------------------------------------#
      #      Find the total and model residuals.                                           #
      #------------------------------------------------------------------------------------#
      x.mod.ok = x.mod[sel]
      x.obs.ok = x.obs[sel]
      x.tot.ok = x.obs.ok - mean(x.obs.ok)
      x.res.ok = x.obs.ok - x.mod.ok
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Find the mean bias, standard deviation of the residuals, and the support for   #
      # the errors being normally distributed around the mean.                             #
      #------------------------------------------------------------------------------------#
      bias        = mean(-x.res.ok)
      sigma       = sd  (-x.res.ok)
      lsq.lnlike  = sum(dnorm(x.res.ok,mean=0,sd=sigma,log=TRUE))
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Find the mean square error of the estimator, and the root mean square error.   #
      #------------------------------------------------------------------------------------#
      mse         = bias^2 + sigma^2
      rmse        = sqrt(mse)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Find the coefficient of determination (R2).  If the number of parameters is    #
      # given, then we find the adjusted R^2 (the one that penalises due to the number of  #
      # parameters).  Otherwise we can't correct, so we just compare the sum of the        #
      # variances.                                                                         #
      #------------------------------------------------------------------------------------#
      ss.tot    = sum(x.tot.ok^2)
      ss.err    = sum(x.res.ok^2)
      r.squared = 1. - df.tot * ss.err / ( df.err * ss.tot )
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Find the estimator's fraction of variance unexplained (FVU).  Because we use   #
      # MSE instead of ss.err, this is not 1 - R2.                                         #
      #------------------------------------------------------------------------------------#
      fvue = mse / var(x.obs.ok)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Run the Shapiro-Wilk test to determine whether the residuals are normally      #
      # distributed or not.                                                                #
      #------------------------------------------------------------------------------------#
      sw.test      = sw.test(x.res.ok)
      sw.statistic = sw.test$statistic
      sw.p.value   = sw.test$p.value
      #------------------------------------------------------------------------------------#
   }else{
      #------------------------------------------------------------------------------------#
      #    Not enough data points.                                                         #
      #------------------------------------------------------------------------------------#
      bias         = NA
      sigma        = NA
      lsq.lnlike   = NA
      mse          = NA
      rmse         = NA
      ss.tot       = NA
      ss.err       = NA
      r.squared    = NA
      fvue         = NA
      sw.statistic = NA
      sw.p.value   = NA
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Return everything to the user.                                                    #
   #---------------------------------------------------------------------------------------#
   ans = list ( n            = n.ok 
              , df.tot       = df.tot
              , df.err       = df.err
              , bias         = bias
              , sigma        = sigma
              , lsq.lnlike   = lsq.lnlike
              , mse          = mse
              , rmse         = rmse
              , ss.tot       = ss.tot
              , ss.err       = ss.err
              , r.squared    = r.squared
              , fvue         = fvue
              , sw.statistic = sw.statistic
              , sw.p.value   = sw.p.value
              )#end list
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function
#==========================================================================================#
#==========================================================================================#
