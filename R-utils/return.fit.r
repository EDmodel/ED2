#==========================================================================================#
#==========================================================================================#
#      This function predicts the change in biomass as a function of the return period     #
# from ED simulations using a logistic function.                                           #
# x is the set of parameters that are optimised for the dataset.                           #
#------------------------------------------------------------------------------------------#
predict.change <<- function(x,datum){


   #----- Copy parameters. ----------------------------------------------------------------#
   y0 = x[1]
   a  = x[2]
   b  = x[3]
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #     Estimate biomass change.                                                          #
   #---------------------------------------------------------------------------------------#
   change = y0 + a / datum$pret^b
   #---------------------------------------------------------------------------------------#
   rm(list=c("y0","a","b"))
   return(change)
}#end predict.change
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#      This function predicts the return period associated with a change in biomass using  #
# the inverse function.                                                                    #
# x is the set of parameters that are optimised for the dataset.                           #
#------------------------------------------------------------------------------------------#
predict.pret <<- function(x,y){


   #----- Copy parameters. ----------------------------------------------------------------#
   y0 = x[1]
   a  = x[2]
   b  = x[3]
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #     Estimate biomass change.                                                          #
   #---------------------------------------------------------------------------------------#
   pret = ( a / (y - y0) )^(1/b)
   #---------------------------------------------------------------------------------------#

   rm(list=c("y0","a","b"))
   return(pret)
}#end predict.change
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#       This function finds the sum of the squares, which is the log likelihood if we      #
# asume the errors to be independent and normally distributed (a big assumption).          #
#==========================================================================================#
#==========================================================================================#
support.change <<- function(x,datum,skew=FALSE){

   change.guess = predict.change(x,datum)
   residual     = change.guess - datum$change
   residual     = residual[is.finite(residual)]
   support      = sn.lsq(r=residual,skew=skew)

   rm(change.guess,residual)
   return(support)
}#end function support.change
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#       This function optimises the change as a function of the drought return period.     #
#==========================================================================================#
#==========================================================================================#
change.return.optim <<- function( datum
                                , y.crit      = 0
                                , first       = NULL
                                , skew        = FALSE
                                , tol.gain    = 0.001
                                , tol.optim   = sqrt(.Machine$double.eps)
                                , is.debug    = FALSE
                                , maxit       = 100
                                , maxit.optim = 20000
                                , verbose     = FALSE
                                ){

   #---------------------------------------------------------------------------------------#
   #     Data selection and total number of parameters.                                    #
   #---------------------------------------------------------------------------------------#
   use      = is.finite(datum$change)  & is.finite(datum$pret)
   n.use    = sum(use)
   n.par    = 3
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Check how verbose to be.                                                          #
   #---------------------------------------------------------------------------------------#
   if (is.logical(verbose)){
      optim.verbose = FALSE
   }else{
      optim.verbose = verbose >= 2
      verbose       = verbose >= 1
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Check that there are at least four valid data points, otherwise, crash!           #
   #---------------------------------------------------------------------------------------#
   if (n.use <= n.par+1){
      cat0(" - Number of valid points: ",n.use)
      cat0(" - Minimum number of valid points: ",n.par+1)
      stop(" Too few valid data points!")
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    If first is null, guess it!                                                        #
   #---------------------------------------------------------------------------------------#
   if (is.null(first)){
      #----- First guess.  Anything on the right quadrant should do. ----------------------#
      x.1st    = c( max(datum$change),-1,1)
      #------------------------------------------------------------------------------------#
   }else{
      #----- The user has the last word (well, almost, see below). ------------------------#
      x.1st    = first
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Set options to the optimiser.                                                     #
   #---------------------------------------------------------------------------------------#
   ctrl.optim = list( trace = optim.verbose
                    , REPORT  = 1
                    , fnscale = -1
                    , maxit   = maxit.optim
                    , reltol  = tol.optim
                    , ndeps   = rep(sqrt(tol.optim),times=n.par)
                    )#end list
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #    If we are going to use skew normal, we update the first guess with Gaussian, so    #
   # it is not too far from the answer.  If the residuals are too large, the skew          #
   # normal may not converge.                                                              #
   #---------------------------------------------------------------------------------------#
   opt         = optim( par     = x.1st
                      , fn      = support.change
                      , skew    = FALSE
                      , datum   = datum[use,]
                      , control = ctrl.optim
                      , hessian = TRUE
                      )#end optim
   success     = is.finite(opt$convergence) && opt$convergence == 0
   #----- Update both the first guess and the scale in case it converged. -----------------#
   if (success){
      x.1st      = opt$par
      ctrl.optim = modifyList(x=ctrl.optim,val=list(parscale=abs(x.1st)))
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Show first guess. ---------------------------------------------------------------#
   if (verbose){
      cat0("             > First guess "
          ,";   y0: ",sprintf("%.3f",x.1st[1])
          ,";   a:  ",sprintf("%.3f",x.1st[2])
          ,";   b:  ",sprintf("%.3f",x.1st[3])
          )#end cat0
   }#end if
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Optimise the parameters.                                                          #
   #---------------------------------------------------------------------------------------#
   bef.support = -Inf
   it          = 0
   iterate     = TRUE
   skew.optim  = skew
   nsteps      = 0
   while (iterate){
      it      = it + 1
      opt     = optim( par     = x.1st
                     , fn      = support.change
                     , skew    = skew.optim
                     , datum   = datum[use,]
                     , control = list( trace   = optim.verbose
                                     , fnscale = -1
                                     , maxit   = maxit.optim
                                     , REPORT  = 1
                                     )#end list
                     , hessian = TRUE
                     )#end optim
      success     = is.finite(opt$convergence) && opt$convergence == 0
      now.support = opt$value
      nsteps.now  = opt$counts["function"]
      nsteps      = nsteps + nsteps.now

      #----- Check whether we are improving things. ---------------------------------------#
      if (is.finite(bef.support)){
         gain = 200. * (now.support - bef.support) / (abs(now.support) + abs(bef.support))
      }else{
         gain = 200.
      }#end if
      #------------------------------------------------------------------------------------#


      #----- Check whether we should iterate. ---------------------------------------------#
      iterate = success && ( gain > (100. * tol.gain) ) && it < maxit
      if (verbose){
         cat0("             > Iteration ",it
             ,";   Converged: ",! iterate
             ,";   Success: "  ,success
             ,";   Support: "  ,sprintf("%.3f",now.support)
             ,";   Gain: "     ,sprintf("%.3f",gain)   
             ,";   Steps: "    ,sprintf("%6i",nsteps.now)
             )#end cat0
      }#end if
      #------------------------------------------------------------------------------------#



      #----- Check whether the optimiser worked at least once. ----------------------------#
      if (! success && it == 1){
         #---------------------------------------------------------------------------------#
         #     Skew normal usually fits better but sometimes it doesn't converge...        #
         #---------------------------------------------------------------------------------#
         if (skew.optim){
            skew.optim  = FALSE
            bef.support = -Inf
            it          = 0
            nsteps      = 0
            iterate     = TRUE
            if (verbose){
               cat0("               @ Skew normal failed, falling back to normal ")
            }else{
               warning(" Skew normal optimiser failed, falling back to normal...")
            }#end if
         }else if (is.debug){
            browser()
         }else{
            stop (" Solution of change didn't converge with normal distribution...")
         }#end if
         #---------------------------------------------------------------------------------#
      }else if (success){
         bef.support = now.support
         x.1st       = opt$par
         ctrl.optim  = modifyList(x=ctrl.optim,val=list(parscale=abs(x.1st)))
      }#end if
      if (iterate) rm(opt)
      #------------------------------------------------------------------------------------#
   }#end while
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Copy the results to ans.                                                          #
   #---------------------------------------------------------------------------------------#
   ans                 = list()
   ans$hessian         = opt$hessian
   ans$df              = n.use - n.par
   ans$coefficients    = opt$par
   ans$std.err         = sqrt(diag(solve(-ans$hessian)))
   ans$t.value         = ans$coefficients / ans$std.err
   ans$p.value         = 2.0 * pt(-abs(ans$t.value),df=ans$df)
   ans$first.guess     = x.1st
   ans$support         = opt$value
   ans$distribution    = ifelse(skew.optim,"Skew Normal","Normal")
   ans$nsteps          = nsteps
   #----- Save the fitted values and the residuals. ---------------------------------------#
   ans$fitted.values   = predict.change(x=ans$coefficients,datum=datum)
   ans$residuals       = ans$fitted.values[use] - datum$change[use]
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Estimate the global standard error and R2.                                        #
   #---------------------------------------------------------------------------------------#
   ss.err           = sum(ans$residuals^2)
   df.err           = ans$df
   mean.y           = mean(datum$change[use])
   ss.tot           = sum((datum$change[use]-mean.y)^2)
   df.tot           = n.use - 1
   ans$r.square     = 1.0 - ss.err * df.tot / ( ss.tot * df.err )
   ans$sigma        = sqrt(ss.err / ans$df)
   res.mean         = mean(ans$residuals)
   res.sdev         = sd  (ans$residuals)
   res.rmse         = sqrt(res.mean^2+res.sdev^2)
   res.skew         = skew(ans$residuals)
   res.stats        = sn.stats(ans$residuals)
   ans$res.summary  = list( mean     = res.mean
                          , sdev     = res.sdev
                          , skew     = res.skew
                          , rmse     = res.rmse
                          , location = res.stats[1]
                          , scale    = res.stats[2]
                          , shape    = res.stats[3]
                          )#end list
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Estimate the critical value given the estimates.                                  #
   #---------------------------------------------------------------------------------------#
   x.crit        = predict.pret(x=ans$coefficients,y=y.crit)
   y0            = ans$coefficients[1]
   a             = ans$coefficients[2]
   b             = ans$coefficients[3]
   z             = 1. / datum$pret[use]^b
   z.crit        = 1. / x.crit^b
   z.mean        = mean(z)
   z.sdev        = sd  (z)
   se.y          = ans$sigma * sqrt(1/n.use + (z.crit - z.mean)^2/(ans$df * z.sdev^2))
   se.x          = se.y * x.crit^(b+1) / abs(a*b)
   ans$x.crit    = x.crit
   ans$x.crit.se = se.x
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Print the result on standard output.                                              #
   #---------------------------------------------------------------------------------------#
   if (verbose){
      cat0("             > Results "
          ,";   y0: "    ,sprintf("%.3f",y0          )
          ,";   a: "     ,sprintf("%.3f",a           )
          ,";   b: "     ,sprintf("%.3f",b           )
          ,";   x.crit: ",sprintf("%.3f",x.crit      )
          )#end cat0
      cat0("             > Fit "
          ,";   R2: "  ,sprintf("%.3f",ans$r.square)
          ,";   Bias: ",sprintf("%.3f",res.mean    )
          ,";   RMSE: ",sprintf("%.3f",res.rmse    )
          ,";   Skew: ",sprintf("%.3f",res.skew    )
          )#end cat0
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Return answer -------------------------------------------------------------------#
   bye = c("use","n.use","n.par","optim.verbose","x.1st","opt","success","bef.support"
          ,"it","iterate","skew.optim","nsteps","success","now.support","nsteps.now"
          ,"gain","ss.err","df.err","mean.y","ss.tot","df.tot","res.mean","res.sdev"
          ,"res.rmse","res.skew","res.stats","x.crit","y0","a","b","z","z.crit"
          ,"z.mean","z.sdev","se.y","se.x")
   rm(list=bye)
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function change.return.optim
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#       This function optimises the change as a function of the drought return period.     #
#==========================================================================================#
#==========================================================================================#
change.return.nlrob <<- function( datum
                                , y.crit      = 0
                                , first       = NULL
                                , tol.optim   = sqrt(.Machine$double.eps)
                                , is.debug    = FALSE
                                , maxit       = 100
                                , n.boot      = 1000
                                , verbose     = FALSE
                                ){

   #---------------------------------------------------------------------------------------#
   #     Data selection and total number of parameters.                                    #
   #---------------------------------------------------------------------------------------#
   use      = which(is.finite(datum$change)  & is.finite(datum$pret))
   n.use    = length(use)
   n.par    = 3
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Check how verbose to be.                                                          #
   #---------------------------------------------------------------------------------------#
   if (is.logical(verbose)){
      optim.verbose = FALSE
   }else{
      optim.verbose = verbose >= 2
      verbose       = verbose >= 1
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Check that there are at least four valid data points, otherwise, crash!           #
   #---------------------------------------------------------------------------------------#
   if (n.use <= n.par+1){
      cat0(" - Number of valid points: ",n.use)
      cat0(" - Minimum number of valid points: ",n.par+1)
      stop(" Too few valid data points!")
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    If first is null, guess it!                                                        #
   #---------------------------------------------------------------------------------------#
   if (is.null(first)){
      #----- First guess.  Anything on the right quadrant should do. ----------------------#
      x.1st = c( max(datum$change),-1,1)
      x.1st = change.return.optim(datum,y.crit=y.crit)
      x.1st = x.1st$coefficients
      #------------------------------------------------------------------------------------#
   }else{
      #----- The user has the last word (well, almost, see below). ------------------------#
      x.1st    = first
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Show first guess. ---------------------------------------------------------------#
   if (verbose){
      cat0("             > First guess "
          ,";   y0: ",sprintf("%.3f",x.1st[1])
          ,";   a:  ",sprintf("%.3f",x.1st[2])
          ,";   b:  ",sprintf("%.3f",x.1st[3])
          )#end cat0
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      We now cheat because we don't want to take risks and send the parameters         #
   # straight to the wrong zone.  a must be negative, whilst b must be positive, so we     #
   # transform the parameters internally.                                                  #
   #---------------------------------------------------------------------------------------#
   y0.1st      = x.1st[1]
   a.prime.1st = log(-x.1st[2])
   b.prime.1st = log( x.1st[3])
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Save the data to a shorter matrix.                                               #
   #---------------------------------------------------------------------------------------#
   datum.use         = datum[use,]
   datum.pred        = datum
   datum.pred$change = datum.pred$change + NA
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Optimise the parameters.                                                          #
   #---------------------------------------------------------------------------------------#
   opt          = try( nlrob( formula   = change ~ I(y0 - exp(a.prime) / pret^exp(b.prime))
                            , data      = datum.use
                            , start     = list( y0      = y0.1st
                                              , a.prime = a.prime.1st
                                              , b.prime = b.prime.1st
                                              )#end list
                            , maxit     = maxit
                            , na.action = na.exclude
                            , tol       = tol.optim
                            , control   = list( maxiter   = maxit
                                              , tol       = tol.optim
                                              , minFactor = 1/32768
                                              )#end list
                            )#end nlfun
                     )#end try
   if ("try-error" %in% is(opt)){
      ans = change.return.htscd( datum
                               , y.crit      = y.crit
                               , first       = first
                               , tol.optim   = tol.optim
                               , is.debug    = is.debug
                               , maxit       = maxit
                               , verbose     = verbose
                               )#end change.return.optim
      return(ans)
   }else{
      summ.opt     = summary(opt)
      f.mu         = summ.opt$coeff[,1]
      f.sigma      = summ.opt$coeff[,2]
      f.expmu      = exp(f.mu+0.5*f.sigma^2)
      f.expsigma   = sqrt((exp(f.sigma^2)-1)*exp(2*f.mu+f.sigma^2))
      coeff        = c(f.mu[1],-f.expmu[2],f.expmu[3])
      se.coeff     = c(f.sigma[1],f.expsigma[2],f.expsigma[3])
      names(coeff) = c("y0","a","b")
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Copy the results to ans.                                                          #
   #---------------------------------------------------------------------------------------#
   ans                 = list()
   ans$object          = opt
   ans$df              = n.use - n.par
   ans$coefficients    = coeff
   ans$std.err         = se.coeff
   ans$t.value         = f.mu / sqrt(f.sigma)
   ans$p.value         = 2.0 * pt(-abs(ans$t.value),df=ans$df)
   ans$first.guess     = x.1st
   #----- Save the fitted values and the residuals. ---------------------------------------#
   ans$fitted.values   = predict.change(x=ans$coefficients,datum=datum)
   ans$residuals       = ans$fitted.values - datum$change
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Estimate the global standard error and R2.                                        #
   #---------------------------------------------------------------------------------------#
   ss.err           = sum(ans$residuals[use]^2)
   df.err           = ans$df
   mean.y           = mean(datum$change[use])
   ss.tot           = sum((datum$change[use]-mean.y)^2)
   df.tot           = n.use - 1
   ans$r.square     = 1.0 - ss.err * df.tot / ( ss.tot * df.err )
   ans$sigma        = sqrt(ss.err / ans$df)
   res.mean         = mean(ans$residuals)
   res.sdev         = sd  (ans$residuals)
   res.rmse         = sqrt(res.mean^2+res.sdev^2)
   res.skew         = skew(ans$residuals)
   res.stats        = sn.stats(ans$residuals)
   ans$res.summary  = list( mean     = res.mean
                          , sdev     = res.sdev
                          , skew     = res.skew
                          , rmse     = res.rmse
                          , location = res.stats[1]
                          , scale    = res.stats[2]
                          , shape    = res.stats[3]
                          )#end list
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Estimate the critical value given the estimates.                                  #
   #---------------------------------------------------------------------------------------#
   x.crit        = predict.pret(x=ans$coefficients,y=y.crit)
   y0            = ans$coefficients[1]
   a             = ans$coefficients[2]
   b             = ans$coefficients[3]
   z             = 1. / datum$pret[use]^b
   z.crit        = 1. / x.crit^b
   z.mean        = mean(z)
   z.sdev        = sd  (z)
   se.y          = ans$sigma * sqrt(1/n.use + (z.crit - z.mean)^2/(ans$df * z.sdev^2))
   se.x          = se.y * x.crit^(b+1) / abs(a*b)
   ans$x.crit    = x.crit
   ans$x.crit.se = se.x
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Print the result on standard output.                                              #
   #---------------------------------------------------------------------------------------#
   if (verbose){
      cat0("             > Results "
          ,";   y0: "    ,sprintf("%.3f",y0          )
          ,";   a: "     ,sprintf("%.3f",a           )
          ,";   b: "     ,sprintf("%.3f",b           )
          ,";   x.crit: ",sprintf("%.3f",x.crit      )
          )#end cat0
      cat0("             > Fit "
          ,";   R2: "  ,sprintf("%.3f",ans$r.square)
          ,";   Bias: ",sprintf("%.3f",res.mean    )
          ,";   RMSE: ",sprintf("%.3f",res.rmse    )
          ,";   Skew: ",sprintf("%.3f",res.skew    )
          )#end cat0
   }#end if
   #---------------------------------------------------------------------------------------#

   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function change.return.nls
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#       This function optimises the change as a function of the drought return period.     #
#==========================================================================================#
#==========================================================================================#
change.return.htscd <<- function( datum
                                , y.crit      = 0
                                , first       = NULL
                                , tol.optim   = sqrt(.Machine$double.eps)
                                , is.debug    = FALSE
                                , maxit       = 100
                                , n.boot      = 1000
                                , verbose     = FALSE
                                , err.method  = "hess"
                                ){

   #---------------------------------------------------------------------------------------#
   #     Data selection and total number of parameters.                                    #
   #---------------------------------------------------------------------------------------#
   use      = which(is.finite(datum$change)  & is.finite(datum$pret))
   n.use    = length(use)
   n.par    = 3
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Check how verbose to be.                                                          #
   #---------------------------------------------------------------------------------------#
   if (is.logical(verbose)){
      optim.verbose = FALSE
   }else{
      optim.verbose = verbose >= 2
      verbose       = verbose >= 1
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Check that there are at least four valid data points, otherwise, crash!           #
   #---------------------------------------------------------------------------------------#
   if (n.use <= n.par+1){
      cat0(" - Number of valid points: ",n.use)
      cat0(" - Minimum number of valid points: ",n.par+1)
      stop(" Too few valid data points!")
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    If first is null, guess it!                                                        #
   #---------------------------------------------------------------------------------------#
   if (is.null(first)){
      #----- First guess.  Anything on the right quadrant should do. ----------------------#
      x.1st    = c( max(datum$change),-1,1,0)
      #------------------------------------------------------------------------------------#
   }else{
      #----- The user has the last word (well, almost, see below). ------------------------#
      x.1st    = c(first,0)
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Show first guess. ---------------------------------------------------------------#
   if (verbose){
      cat0("             > First guess "
          ,";   y0: ",sprintf("%.3f",x.1st[1])
          ,";   a:  ",sprintf("%.3f",x.1st[2])
          ,";   b:  ",sprintf("%.3f",x.1st[3])
          ,";   s:  ",sprintf("%.3f",x.1st[4])
          )#end cat0
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      We now cheat because we don't want to take risks and send the parameters         #
   # straight to the wrong zone.  a must be negative, whilst b must be positive, so we     #
   # transform the parameters internally.                                                  #
   #---------------------------------------------------------------------------------------#
   y0.1st      = x.1st[1]
   a.prime.1st = log(-x.1st[2])
   b.prime.1st = log( x.1st[3])
   s0.1st      = x.1st[4]
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Save the data to a shorter matrix.                                               #
   #---------------------------------------------------------------------------------------#
   datum.use         = datum[use,]
   datum.pred        = datum
   datum.pred$change = datum.pred$change + NA
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Optimise the parameters.                                                          #
   #---------------------------------------------------------------------------------------#
   lsq.formula = "change ~ I(y0 - exp(a.prime) / pret^exp(b.prime))"
   sig.formula = "~ I(pret^s0)"
   opt         = try( optim.lsq.htscd( lsq.formula = lsq.formula
                                     , sig.formula = sig.formula
                                     , data        = datum.use
                                     , lsq.first   = list( y0      = y0.1st
                                                         , a.prime = a.prime.1st
                                                         , b.prime = b.prime.1st
                                                         )#end list
                                     , sig.first   = list( s0      = s0.1st
                                                         )#end list
                                     , err.method  = err.method
                                     , is.debug    = TRUE
                                     , n.boot      = n.boot
                                     , verbose     = optim.verbose
                                     )#end optim.lsq.htscd
                    )#end try
   if ("try-error" %in% is(opt)){
      ans = change.return.optim( datum
                               , y.crit      = y.crit
                               , first       = first
                               , tol.optim   = tol.optim
                               , is.debug    = is.debug
                               , maxit       = maxit
                               , verbose     = verbose
                               )#end change.return.optim
      return(ans)
   }else{
      summ.opt     = summary(opt)
      f.mu         = summ.opt$coefficients[,1]
      f.sigma      = summ.opt$coefficients[,2]
      f.expmu      = exp(f.mu)
      f.expsigma   = sqrt((exp(f.sigma^2)-1)*exp(2*f.mu+f.sigma^2))
      coeff        = c(f.mu[1],-f.expmu[2],f.expmu[3],f.mu[4])
      se.coeff     = c(f.sigma[1],f.expsigma[2],f.expsigma[3],f.sigma[4])
      names(coeff) = c("y0","a","b","s0")
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Copy the results to ans.                                                          #
   #---------------------------------------------------------------------------------------#
   ans                 = list()
   ans$object          = opt
   ans$df              = n.use - n.par
   ans$coefficients    = coeff
   ans$std.err         = se.coeff
   ans$t.value         = f.mu / sqrt(f.sigma)
   ans$p.value         = 2.0 * pt(-abs(ans$t.value),df=ans$df)
   ans$first.guess     = x.1st
   #----- Save the fitted values and the residuals. ---------------------------------------#
   if (err.method %in% "boot"){
      ans$fitted.values = predict(object=opt,newdata=datum,pred.boot=TRUE )$fit
   }else{
      ans$fitted.values = predict(object=opt,newdata=datum,pred.boot=FALSE)
   }#end if (err.method %in% "boot")
   ans$residuals       = ans$fitted.values - datum$change
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Estimate the global standard error and R2.                                        #
   #---------------------------------------------------------------------------------------#
   ss.err           = sum(ans$residuals[use]^2)
   df.err           = ans$df
   mean.y           = mean(datum$change[use])
   ss.tot           = sum((datum$change[use]-mean.y)^2)
   df.tot           = n.use - 1
   ans$r.square     = 1.0 - ss.err * df.tot / ( ss.tot * df.err )
   ans$sigma        = sqrt(ss.err / ans$df)
   res.mean         = mean(ans$residuals)
   res.sdev         = sd  (ans$residuals)
   res.rmse         = sqrt(res.mean^2+res.sdev^2)
   res.skew         = skew(ans$residuals)
   res.stats        = sn.stats(ans$residuals)
   ans$res.summary  = list( mean     = res.mean
                          , sdev     = res.sdev
                          , skew     = res.skew
                          , rmse     = res.rmse
                          , location = res.stats[1]
                          , scale    = res.stats[2]
                          , shape    = res.stats[3]
                          )#end list
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Estimate the critical value given the estimates.                                  #
   #---------------------------------------------------------------------------------------#
   x.crit        = predict.pret(x=ans$coefficients,y=y.crit)
   y0            = ans$coefficients[1]
   a             = ans$coefficients[2]
   b             = ans$coefficients[3]
   s0            = ans$coefficients[4]
   z             = 1. / datum$pret[use]^b
   z.crit        = 1. / x.crit^b
   z.mean        = mean(z)
   z.sdev        = sd  (z)
   se.y          = ans$sigma * sqrt(1/n.use + (z.crit - z.mean)^2/(ans$df * z.sdev^2))
   se.x          = se.y * x.crit^(b+1) / abs(a*b)
   ans$x.crit    = x.crit
   ans$x.crit.se = se.x
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Print the result on standard output.                                              #
   #---------------------------------------------------------------------------------------#
   if (verbose){
      cat0("             > Results "
          ,";   y0: "    ,sprintf("%.3f",y0          )
          ,";   a: "     ,sprintf("%.3f",a           )
          ,";   b: "     ,sprintf("%.3f",b           )
          ,";   s0: "    ,sprintf("%.3f",s0          )
          ,";   x.crit: ",sprintf("%.3f",x.crit      )
          )#end cat0
      cat0("             > Fit "
          ,";   R2: "  ,sprintf("%.3f",ans$r.square)
          ,";   Bias: ",sprintf("%.3f",res.mean    )
          ,";   RMSE: ",sprintf("%.3f",res.rmse    )
          ,";   Skew: ",sprintf("%.3f",res.skew    )
          )#end cat0
   }#end if
   #---------------------------------------------------------------------------------------#

   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function change.return.htscd
#==========================================================================================#
#==========================================================================================#
