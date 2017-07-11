#==========================================================================================#
#==========================================================================================#
#      This function uses an iterative technique to determine the u* filter for a given    #
# CO2 flux dataset.  The method is based on:                                               #
#                                                                                          #
#    Gu, L.; E. M. Falge, T. Boden, D. D. Baldocchi, T. A. Black, S. R. Saleska, T. Suni,  #
#        S. B. Verma, T. Vesala, S. C. Wofsy, L. Xu, 2005: Objective threshold determina-  #
#        tion for nighttime eddy flux filtering.  Ag. Forest Meteorol., 128, 179-197.      #
#------------------------------------------------------------------------------------------#
ustar.filter = function(dat                # The dataset, it must be a data.frame
                       ,frac.mw    = 0.05  # Relative window size for t-test
                       ,min.nmw    = 20    # Minimum window size
                       ,max.nmw    = Inf   # Minimum window size
                       ,pval.ustar = 0.05  # p-value to reject null hypothesis
                       ,quant.unst = 0.20  # Minimum fraction of valid nighttime nee
                       ,retain.min = 1/3   # Minimum amount of data to retain (warning)
                       ,toler      = 0.01  # Tolerance for convergence test
                       ,up.limit   = TRUE  # Filter the upper limit
                       ,verbose    = FALSE # Print additional information
                       ,test.stor  = TRUE  # Should we include storage in the test?
                       ){


   #---------------------------------------------------------------------------------------#
   #     Number of points.                                                                 #
   #---------------------------------------------------------------------------------------#
   ndat   = length(dat$when)
   ndatm1 = ndat - 1
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Save the flux into a temporary array.                                             #
   #---------------------------------------------------------------------------------------#
   if (verbose) print("   + Copying the flux...")
   if (test.stor){
      dat.flux = dat$fco2 + dat$storco2
   }else{
      dat.flux = dat$fco2
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Make the night time filter, which remains constant throughout the process.        #
   #---------------------------------------------------------------------------------------#
   if (verbose) print("   + Flagging night time data...")
   nighttime = dat$cosz < cosz.twilight
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Tonight is the night associated with 12:01 am of a particular night.  So the      #
   # evening of Jan 26 will considered the same night as the dawn of Jan 27.               #
   #---------------------------------------------------------------------------------------#
   if (verbose) print("   + Counting number of nights...")
   dat.tonight             = cumsum(nighttime & nighttime != c(FALSE,nighttime[1:ndatm1]))
   dat.tonight[!nighttime] = NA
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Make the night time filter, which remains constant throughout the process.        #
   #---------------------------------------------------------------------------------------#
   if (verbose) print("   + Flagging valid data...")
   use = is.finite(dat$ustar)   & is.finite(dat.flux)    & is.finite(dat$atm.tmp)
   #---------------------------------------------------------------------------------------#


   #----- Initialise the low and high u* filter (outer loop). -----------------------------#
   ustlo.out = 0.
   usthi.out = Inf
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     This is the outer iteration loop.                                                 #
   #---------------------------------------------------------------------------------------#
   if (verbose) print("   + Outer loop...")
   out.iter = TRUE
   oo       = 0
   while (out.iter){

      #------------------------------------------------------------------------------------#
      #    Counter.                                                                        #
      #------------------------------------------------------------------------------------#
      oo = oo + 1
      if (verbose) print(paste("     - Iteration: ",oo,":",sep=""))
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Selections.  The "sel" flag keeps only the night time data that has finite     #
      # values, whereas the "sel.fit" flag includes only those "sel" values whose ustar    #
      # is between the outer loop ustar limits.                                            #
      #------------------------------------------------------------------------------------#
      if (verbose) print(paste("       * Selecting data...",sep=""))
      sel      = use & nighttime
      sel.fit  = sel & dat$ustar >= ustlo.out & dat$ustar <= usthi.out
      nsel     = sum(sel)
      nsel.fit = sum(sel.fit)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    Retain the data that we will consider in this part.                             #
      #------------------------------------------------------------------------------------#
      reco     = dat.flux   [sel    ]
      temp     = dat$atm.tmp[sel    ]
      ustar    = dat$ustar  [sel    ]
      night    = dat.tonight[sel    ]
      reco.fit = dat.flux   [sel.fit]
      temp.fit = dat$atm.tmp[sel.fit]
      #------------------------------------------------------------------------------------#



      #----- Find a linear fit between temperature and respiration. -----------------------#
      if (verbose) print(paste("       * Fitting a temperature line...",sep=""))
      fit      = lm(reco.fit ~ temp.fit)
      summ.fit = summary(fit)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Now we check whether there is an obvious relationship between the two.         #
      # Normally there will be a relationship in temperate and subpolar sites, but not in  #
      # equatorial sites.                                                                  #
      #------------------------------------------------------------------------------------#
      a0         = summ.fit$coefficients[1,1]
      a1         = summ.fit$coefficients[2,1]
      r2         = summ.fit$adj.r.squared
      pval.inter = summ.fit$coefficients[1,4]
      pval.slope = summ.fit$coefficients[2,4]
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Determine whether a linear fit is a good fit.  To be considered one, the slope #
      # has to be statistically significant, and the fit should explain something...       #
      #------------------------------------------------------------------------------------#
      lindep = r2 >= r2.min && pval.slope <= pval.max
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Find the expected values (either the mean or the linear fit), then normalise   #
      # the respiration flux by dividing the ecosystem respiration by the expected values. #
      # At this point we use the full night time data set.                                 #
      #------------------------------------------------------------------------------------#
      if (lindep){
         if (verbose) print(paste("       * Linear fit significant! Normalising data..."
                                 ,sep=""))
         reco.exp = a0 + a1 * temp
         reco.sde = summ.fit$sigma
      }else{
         if (verbose) print(paste("       * Linear fit non-significant! Normalising data..."
                                 ,sep=""))
         reco.exp = mean(reco.fit) + 0. * temp
         reco.sde = sd(reco.fit)   + 0. * temp
      }#end if
      recon = (reco - reco.exp) / (reco.sde)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Conduct an outlier detection test.  Now that the data are normalised, this    #
      # means that we only retain values between -3. and 3.                                #
      #------------------------------------------------------------------------------------#
      if (verbose) print(paste("       * Discarding outliers...",sep=""))
      ok    = abs(recon) <= 3.
      reco  = reco  [ok]
      recon = recon [ok]
      temp  = temp  [ok]
      ustar = ustar [ok]
      night = night [ok]
      #------------------------------------------------------------------------------------#


      #----- Initialise the low and high u* filter (inner loop). --------------------------#
      ustlo.in = 0.
      usthi.in = Inf
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     This is the inner iteration loop.                                              #
      #------------------------------------------------------------------------------------#
      in.iter = TRUE
      ii = 0
      while (in.iter){

         #---------------------------------------------------------------------------------#
         #    Counter.                                                                     #
         #---------------------------------------------------------------------------------#
         ii = ii + 1
         if (verbose) print(paste("       * Inner loop, iteration ",ii,":",sep=""))
         #---------------------------------------------------------------------------------#



         #----- Find the median ustar for each night. -------------------------------------#
         if (verbose) print(paste("         > Discarding too stable nights...",sep=""))
         ust.qunst   = tapply(X=ustar,INDEX=night,FUN=quantile,prob=(1.0-quant.unst)
                             ,na.rm=TRUE)
         night.qunst = as.integer(names(ust.qunst))
         #---------------------------------------------------------------------------------#



         #----- Find the nights whose median u* is above ustlo.in -------------------------#
         turb    = night %in% night.qunst[ust.qunst >= ustlo.in]
         sel.in  = ustar >= ustlo.in & ustar <= usthi.in
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #    Keep data from turbulent nights only.                                        #
         #---------------------------------------------------------------------------------#
         turb.reco  = reco  [turb & sel.in]
         turb.recon = recon [turb & sel.in]
         turb.temp  = temp  [turb & sel.in]
         turb.ustar = ustar [turb & sel.in]
         turb.night = night [turb & sel.in] 
         #---------------------------------------------------------------------------------#



         #----- Find the sample size and the size of the moving window. -------------------#
         nturb      = length(turb.reco)
         nmw        = min(max.nmw,max(min.nmw,floor(frac.mw * nturb)))
         if (verbose) print(paste("         > Sample size = ",nturb
                                 ," Window size = ",nmw,"...",sep=""))
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #    Rank dataset from lowest to highest u*, and from highest to lowest u*.       #
         #---------------------------------------------------------------------------------#
         if (verbose) print(paste("         > Ranking data by u*...",sep=""))
         lohi       = order(turb.ustar,na.last=TRUE,decreasing=FALSE)
         hilo       = order(turb.ustar,na.last=TRUE,decreasing=TRUE )
         reco.lohi  = turb.reco [lohi]
         recon.lohi = turb.recon[lohi]
         temp.lohi  = turb.temp [lohi]
         ustar.lohi = turb.ustar[lohi]
         night.lohi = turb.night[lohi]
         reco.hilo  = turb.reco [hilo]
         recon.hilo = turb.recon[hilo]
         temp.hilo  = turb.temp [hilo]
         ustar.hilo = turb.ustar[hilo]
         night.hilo = turb.night[hilo]
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #   Initialise the indices for the moving window.                                 #
         #---------------------------------------------------------------------------------#
         nloa = 1
         nloz = nmw
         nhia = 1
         nhiz = nmw
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #  Here we will check whether the following null hypothesis can be rejected.      #
         #                                                                                 #
         #  H0:  The mean flux of the sub-sample of ranked fluxes, beginning at nloa and   #
         #       ending at nloz is greater than or equal to the mean flux of the entire    #
         #       series.                                                                   #
         #                                                                                 #
         #  If the hypothesis is rejected, we shift nloa and nloz until we find a case in  #
         #  which we can't reject.                                                         #
         #---------------------------------------------------------------------------------#
         if (verbose) print(paste("         > Finding the filter for low end...",sep=""))
         reject = TRUE
         while (reject){
            mytest = t.test(x=recon.lohi[nloa:nloz],y=recon.lohi
                           ,alternative="less")
            reject = mytest$p.value <= pval.ustar
            if (reject){
               nloa = nloa + 1
               nloz = nloz + 1
            }#end if
         }#end while
         if (nloa != 1){
            new.ustlo.in = ustar.lohi[nloa]
         }else{
            new.ustlo.in = ustlo.in
         }#end if
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #  Here we will check whether the following null hypothesis can be rejected.      #
         #  (only if up.limit is true)                                                     #
         #                                                                                 #
         #  H0:  The mean flux of the sub-sample of ranked fluxes, beginning at nhia and   #
         #       ending at nhiz is less than or equal to the mean flux of the entire       #
         #       series.                                                                   #
         #                                                                                 #
         #  If the hypothesis is rejected, we shift nhia and nhiz until we find a case in  #
         #  which we can't reject.                                                         #
         #---------------------------------------------------------------------------------#
         if (up.limit){
            if (verbose) print(paste("         > Finding the filter for high end..."
                                    ,sep=""))
            reject = TRUE
            while (reject){
               mytest = t.test(x=recon.hilo[nhia:nhiz],y=recon.hilo
                              ,alternative="greater")
               reject = mytest$p.value <= pval.ustar
               if (reject){
                  nhia = nhia + 1
                  nhiz = nhiz + 1
               }#end if
            }#end while
            if (nhia != 1){
               new.usthi.in = ustar.hilo[nhia]
            }else{
               new.usthi.in = usthi.in
            }#end if
         }else{
            new.usthi.in = usthi.in
         }#end if
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Check whether the inner loop converged or not.                              #
         #---------------------------------------------------------------------------------#
         if (verbose) print(paste("         > Checking convergence...",sep=""))
         #----- Low end. ------------------------------------------------------------------#
         if (new.ustlo.in == 0. && ustlo.in == 0.){
            error.lo     = 0.
         }else{
            error.lo  = 2.0 * ( abs(new.ustlo.in - ustlo.in) 
                              / (abs(new.ustlo.in) + abs(ustlo.in)) )
         }#end if
         #----- High end. -----------------------------------------------------------------#
         if (! up.limit){
            error.hi = 0.
         }else if (new.usthi.in == Inf && usthi.in == Inf){
            error.hi = 0.
         }else if(usthi.in == Inf){
            error.hi = 2.0
         }else{
            error.hi = 2.0 * ( abs(new.usthi.in - usthi.in) 
                             / (abs(new.usthi.in) + abs(usthi.in)) )
         }#end if
         #---------------------------------------------------------------------------------#



         #----- Update the guesses in any case. -------------------------------------------#
         ustlo.in = new.ustlo.in
         usthi.in = new.usthi.in
         #---------------------------------------------------------------------------------#



         #----- A new iteration is needed unless both ends have converged. ----------------#
         in.iter = (error.lo > toler) | (error.hi > toler)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Print errors to entretain the users if they want...                         #
         #---------------------------------------------------------------------------------#
         if (verbose){
            print(paste("           . Error (low)  = ",sprintf("%.2f",error.lo  )
                                                      ,"...",sep=""))
            print(paste("           . Error (high) = ",sprintf("%.2f",error.hi  )
                                                      ,"...",sep=""))
            print(paste("           . New low  u*  = ",sprintf("%.2f",ustlo.in  )
                                                      ,"...",sep=""))
            print(paste("           . New high u*  = ",sprintf("%.2f",usthi.in  )
                                                      ,"...",sep=""))
            print(paste("           . Converged    = ",! in.iter ,"...",sep=""))
         }#end if
         #---------------------------------------------------------------------------------#
      }#end while
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Update the guesses at the outer loop by using the updated value from the inner #
      # loop.                                                                              #
      #------------------------------------------------------------------------------------#
      new.ustlo.out = ustlo.in
      new.usthi.out = usthi.in
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #     Check whether the outer loop converged or not.                                 #
      #------------------------------------------------------------------------------------#
      if (verbose) print(paste("       * Check convergence of the outer loop...",sep=""))
      #----- Low end. ---------------------------------------------------------------------#
      if (new.ustlo.out == 0. && ustlo.out == 0.){
         error.lo     = 0.
      }else{
         error.lo  = 2.0 * ( abs(new.ustlo.out - ustlo.out) 
                           / (abs(new.ustlo.out) + abs(ustlo.out)) )
      }#end if
      #----- High end. --------------------------------------------------------------------#
      if (! up.limit){
         error.hi  = 0.
      }else if (new.usthi.out == Inf && usthi.out == Inf){
         error.hi  = 0.
      }else if(usthi.out == Inf){
         error.hi  = 2.0
      }else{
         error.hi  = 2.0 * ( abs(new.usthi.out - usthi.out) 
                           / (abs(new.usthi.out) + abs(usthi.out)) )
      }#end if
      #------------------------------------------------------------------------------------#



      #----- Update the guesses in any case. ----------------------------------------------#
      ustlo.out = new.ustlo.out
      usthi.out = new.usthi.out
      #------------------------------------------------------------------------------------#



      #----- A new iteration is needed unless both ends have converged. -------------------#
      out.iter = (error.lo > toler) | (error.hi > toler)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Print errors to entretain the users if they want...                            #
      #------------------------------------------------------------------------------------#
      if (verbose){
         print(paste("         > Error (low)  = ",sprintf("%.2f",error.lo   ),"...",sep=""))
         print(paste("         > Error (high) = ",sprintf("%.2f",error.hi   ),"...",sep=""))
         print(paste("         > New low  u*  = ",sprintf("%.2f",ustlo.out  ),"...",sep=""))
         print(paste("         > New high u*  = ",sprintf("%.2f",usthi.out  ),"...",sep=""))
         print(paste("         > Converged    = ",! out.iter ,"...",sep=""))
      }#end if
      #------------------------------------------------------------------------------------#
   }#end while
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Save the filter as a two-dimension vector, with the u* range that the NEE is       #
   # valid.                                                                                #
   #---------------------------------------------------------------------------------------#
   ustar.thre = c(ustlo.out,usthi.out)
   retain     = dat$ustar >= ustar.thre[1] & dat$ustar <= ustar.thre[2]
   nretain    = sum(retain,na.rm=TRUE)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Warn the user in case the filter cause too much data to be discarded.             #
   #---------------------------------------------------------------------------------------#
   if (nretain / ndat < retain.min){
      warning(paste(" ==> U* filter will retain only "
                   ,sprintf("%.1f",100*nretain/ndat)," of data...",sep=""))
   }else if (verbose){
      print(paste(" + U* filter will retain "
                 ,sprintf("%.1f",100*nretain/ndat)," of data...",sep=""))

   }#end if
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #    Quit, and return the vector with ustar filters.                                    #
   #---------------------------------------------------------------------------------------#
   return(ustar.thre)
   #---------------------------------------------------------------------------------------#
}#end function
#==========================================================================================#
#==========================================================================================#
