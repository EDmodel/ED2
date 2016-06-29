#==========================================================================================#
#==========================================================================================#
#     This function fits a light response curve to given data.                             #
#------------------------------------------------------------------------------------------#
nls.light.response <<- function( par.in
                               , gpp
                               , pred.par = seq( from = min(par.in,na.rm=TRUE)
                                               , to   = max(par.in,na.rm=TRUE)
                                               , length.out = 100
                                               )#end seq
                               , first  = c(1,40,500)
                               , n.boot = 1000
                               , ...
                               ){#end function

   #----- Transform parameters to avoid singularities. ------------------------------------#
   x.1st      = list(a1 = first[1], a2 = log(first[2]), a3 = log(first[3]))
   #---------------------------------------------------------------------------------------#


   #----- Set the input dataset. ----------------------------------------------------------#
   use        = is.finite(par.in) & is.finite(gpp)
   data.in    = data.frame(par=par.in[use],gpp=gpp[use])
   data.pred  = data.frame(par=pred.par)
   n.use      = sum(use)
   n.coeff    = 3
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Main fit.                                                                         #
   #---------------------------------------------------------------------------------------#
   fit.gpp = try( nls( gpp ~ ( a1 + (exp(a2) * par) / (exp(a3) + par) )
                     , data  = data.in
                     , start = x.1st
                     , ...
                     )#end if
                , silent = TRUE
                )#end try
   if ("try-error" %in% is(fit.gpp)){
      #------------------------------------------------------------------------------------#
      #     Make all data sets NA.                                                         #
      #------------------------------------------------------------------------------------#
      coeff.gpp    = rep(NA,times=3)
      std.err.gpp  = rep(NA,times=3)
      tvalue.gpp   = rep(NA,times=3)
      pvalue.gpp   = rep(NA,times=3)
      df.gpp       = NA
      df.tot       = sum(is.finite(gpp)) - 1
      r2.gpp       = NA
      pred.gpp     = NA + pred.par
      expected.gpp = NA + pred.par
      q025.gpp     = NA + pred.par
      q975.gpp     = NA + pred.par
      #------------------------------------------------------------------------------------#
   }else{
      #------------------------------------------------------------------------------------#
      #     Save summary and predictions.                                                  #
      #------------------------------------------------------------------------------------#
      summ.gpp = summary(fit.gpp)
      pred.gpp = predict(fit.gpp,newdata=data.pred)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Copy data.                                                                     #
      #------------------------------------------------------------------------------------#
      coeff.gpp    = c(summ.gpp$coefficients[1,1],exp(summ.gpp$coefficients[2:3,1]))
      std.err.gpp  = c(summ.gpp$coefficients[1,2],exp(summ.gpp$coefficients[2:3,2]))
      tvalue.gpp   = summ.gpp$coefficients[,3]
      pvalue.gpp   = summ.gpp$coefficients[,4]
      df.gpp       = summ.gpp$df[2]
      df.tot       = sum(is.finite(gpp)) - 1
      #------------------------------------------------------------------------------------#



      #----- Find the adjusted R2. --------------------------------------------------------#
      residual.gpp   = ( predict.light.response(coeff.gpp,data.in=data.in,transf=FALSE)
                       - data.in$gpp )
      df.gpp         = n.use - n.coeff
      ss.gpp         = sum(residual.gpp^2,na.rm=TRUE) / df.gpp
      sigma.gpp      = sqrt(ss.gpp)
      ss.tot         = var(data.in$gpp,na.rm=TRUE)
      r2.gpp         = 1.0 - ss.gpp / ss.tot
      #------------------------------------------------------------------------------------#


      #----- Find predicted values and their confidence band. -----------------------------#
      expected.gpp   = pred.gpp
      q025.gpp       = NA + expected.gpp
      q975.gpp       = NA + expected.gpp
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#

   ans = list( coefficients = coeff.gpp
             , std.err      = std.err.gpp
             , tvalue       = tvalue.gpp
             , pvalue       = pvalue.gpp
             , r2           = r2.gpp
             , par          = pred.par
             , gpp          = expected.gpp
             , q025         = q025.gpp
             , q975         = q975.gpp
             )#end list
   return(ans)
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#       Bootstrap solution.                                                                #
#------------------------------------------------------------------------------------------#
boot.nls.light <<- function(data.in,index,x.1st,pred.par,...){
   fit.gpp = try( nls( gpp ~ ( a1 + (exp(a2) * par) / (exp(a3) + par) )
                     , data  = data.in[index,]
                     , start = x.1st
                     , ...
                     )#end if
                , silent = TRUE
                )#end try
   if ("try-error" %in% is(fit.gpp)){
      summ.gpp  = NA
      coeff.gpp = rep(NA,times=3)
      pred.gpp  = NA + pred.par
   }else{


      #----- Revert coefficients to intuitive units. --------------------------------------#
      summ.gpp  = summary(fit.gpp)
      coeff.gpp = summ.gpp$coefficients[,1]
      coeff.gpp[c(2,3)] = exp(coeff.gpp[c(2,3)])
      #------------------------------------------------------------------------------------#

      #----- Predict GPP. -----------------------------------------------------------------#
      pred.gpp  = predict(object = fit.gpp, newdata = data.frame(par=pred.par))
      #------------------------------------------------------------------------------------#
   }#end if


   #---------------------------------------------------------------------------------------#
   #     Bind coefficients and predicted values together.                                  #
   #---------------------------------------------------------------------------------------#
   ans = c(coeff.gpp,pred.gpp)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Discard place holders and return.                                                 #
   #---------------------------------------------------------------------------------------#
   rm(summ.gpp,coeff.gpp,pred.gpp)
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function fits a light response curve to given data.                             #
#------------------------------------------------------------------------------------------#
optim.light.response <<- function( par.in
                                 , gpp
                                 , pred.par    = seq( from = min(par.in,na.rm=TRUE)
                                                    , to   = max(par.in,na.rm=TRUE)
                                                    , length.out = 100
                                                    )#end seq
                                 , first       = NULL
                                 , skew        = FALSE
                                 , tol.optim   = sqrt(.Machine$double.eps)
                                 , n.boot      = 1000
                                 , maxit.optim = 20000
                                 , verbose     = FALSE
                                 ){#end function


   #----- Set the input dataset. ----------------------------------------------------------#
   use        = is.finite(par.in) & is.finite(gpp)
   data.in    = data.frame(par=par.in[use],gpp=gpp[use])
   n.use      = sum(use)
   n.coeff    = 3
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
   if (n.use <= n.coeff+1){
      cat (" - Number of valid points: ",n.use,"\n")
      cat (" - Minimum number of valid points: ",n.coeff+1,"\n")
      stop(" Too few valid data points!")
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    If first is null, guess it!                                                        #
   #---------------------------------------------------------------------------------------#
   if (is.null(first)){
      #----- First guess.  Anything on the right quadrant should do. ----------------------#
      x.1st    = c(0,log(40),log(500))
      #------------------------------------------------------------------------------------#
   }else{
      #----- The user has the last word (well, almost, see below). ------------------------#
      x.1st    = c(first[1],log(first[2]),log(first[3]))
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Set options to the optimiser.                                                     #
   #---------------------------------------------------------------------------------------#
   ctrl.optim = list( trace   = optim.verbose
                    , REPORT  = 1
                    , fnscale = -1
                    , maxit   = maxit.optim
                    , reltol  = tol.optim
                    , ndeps   = rep(sqrt(tol.optim),times=n.coeff)
                    )#end list
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #    If we are going to use skew normal, we update the first guess with Gaussian, so    #
   # it is not too far from the answer.  If the residuals are too large, the skew          #
   # normal may not converge.                                                              #
   #---------------------------------------------------------------------------------------#
   opt         = optim( par     = x.1st
                      , fn      = support.light.response
                      , transf  = TRUE
                      , skew    = FALSE
                      , data.in = data.in
                      , control = ctrl.optim
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
      cat("             > First guess "
         ,";   a1:  ",sprintf("%.2f",x.1st[1])
         ,";   a2:  ",sprintf("%.2f",exp(x.1st[2]))
         ,";   a3:  ",sprintf("%.2f",exp(x.1st[3]))
         ,"\n")
   }#end if
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Main fit.                                                                         #
   #---------------------------------------------------------------------------------------#
   boot.gpp = boot( data      = data.in
                  , statistic = boot.optim.light
                  , R         = n.boot
                  , x.1st     = x.1st
                  , pred.par  = pred.par
                  , skew      = skew
                  , control   = ctrl.optim
                  )#end boot
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #      Find some of the statistics.                                                     #
   #---------------------------------------------------------------------------------------#
   #----- Split bootstrap statistics into coefficients and predictions. -------------------#
   boot.coeff     = boot.gpp$t[ , sequence(n.coeff)]
   boot.pred      = boot.gpp$t[ ,-sequence(n.coeff)]

   #----- Find the coefficient estimate, t-value, p-value, and standard error. ------------#
   coeff.gpp      = boot.gpp$t0[ sequence(n.coeff)]
   #---------------------------------------------------------------------------------------#


   #----- Run t-test for parameters. ------------------------------------------------------#
   t.test.gpp           = apply (X=boot.coeff,MARGIN=2,FUN=t.test)
   t.summ.gpp           = sapply(X=t.test.gpp,FUN=rbind)
   rownames(t.summ.gpp) = names(t.test.gpp[[1]])
   tvalue.gpp           = unlist(t.summ.gpp["statistic",])
   std.err.gpp          = coeff.gpp / tvalue.gpp
   pvalue.gpp           = unlist(t.summ.gpp["p.value"  ,])
   #---------------------------------------------------------------------------------------#



   #----- Find the adjusted R2. -----------------------------------------------------------#
   residual.gpp   = ( predict.light.response(coeff.gpp,data.in=data.in,transf=FALSE)
                    - data.in$gpp )
   df.gpp         = n.use - n.coeff
   ss.gpp         = sum(residual.gpp^2,na.rm=TRUE) / df.gpp
   sigma.gpp      = sqrt(ss.gpp)
   ss.tot         = var(data.in$gpp,na.rm=TRUE)
   r2.gpp         = 1.0 - ss.gpp / ss.tot
   #---------------------------------------------------------------------------------------#


   #----- Find predicted values and their confidence band. --------------------------------#
   expected.gpp   = boot.gpp$t0[-sequence(n.coeff)]
   q025.gpp       = apply(X=boot.pred,MARGIN=2,quantile,0.025,na.rm=TRUE)
   q975.gpp       = apply(X=boot.pred,MARGIN=2,quantile,0.975,na.rm=TRUE)
   #---------------------------------------------------------------------------------------#




   #----- List with results. --------------------------------------------------------------#
   ans = list( coefficients = coeff.gpp
             , std.err      = std.err.gpp
             , tvalue       = tvalue.gpp
             , pvalue       = pvalue.gpp
             , r2           = r2.gpp
             , par          = pred.par
             , gpp          = expected.gpp
             , q025         = q025.gpp
             , q975         = q975.gpp
             )#end list
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#       Bootstrap solution.                                                                #
#------------------------------------------------------------------------------------------#
boot.optim.light <<- function(data.in,index,x.1st,pred.par,skew=FALSE,...){

   #----- Adjust curve. -------------------------------------------------------------------#
   opt         = optim( par     = x.1st
                      , fn      = support.light.response
                      , transf  = TRUE
                      , skew    = FALSE
                      , data.in = data.in[index,]
                      , ...
                      )#end optim
   #---------------------------------------------------------------------------------------#


   #----- Revert coefficients to intuitive units. -----------------------------------------#
   coeff.gpp         = opt$par
   coeff.gpp[c(2,3)] = exp(coeff.gpp[c(2,3)])
   #---------------------------------------------------------------------------------------#



   #----- Find the predictions for the current fit. ---------------------------------------#
   pred.gpp      = predict.light.response(x=coeff.gpp,data.in=data.frame(par=pred.par))
   #---------------------------------------------------------------------------------------#

   #----- Bind coefficients and predicted values to a single vector. ----------------------#
   ans = c(coeff.gpp,pred.gpp)
   #---------------------------------------------------------------------------------------#

   rm(opt,coeff.gpp,pred.gpp)
   return(ans)
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#       Predict function.                                                                  #
#------------------------------------------------------------------------------------------#
predict.light.response <<- function(x,data.in,transf=FALSE){


   #---------------------------------------------------------------------------------------#
   #     For optimisation, it is better to transform the coefficients to avoid them going  #
   # to dangerous negative values.  Otherwise, non-transformed coefficients are more       #
   # intuitive.                                                                            #
   #---------------------------------------------------------------------------------------#
   if (transf){
      a1 = x[1]
      a2 = exp(x[2])
      a3 = exp(x[3])
   }else{
      a1 = x[1]
      a2 = x[2]
      a3 = x[3]
   }#end if
   #---------------------------------------------------------------------------------------#

   pred.gpp = a1 + a2 * data.in$par / (a3 + data.in$par)

   return(pred.gpp)
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#       Support function.                                                                  #
#------------------------------------------------------------------------------------------#
support.light.response <<- function(x,data.in,transf=FALSE,skew=FALSE){

   pred.gpp = predict.light.response(x,data.in,transf=transf)
   residual = pred.gpp - data.in$gpp
   residual = residual[is.finite(residual)]
   lnlike   = sn.lsq(r=residual,skew=skew)

   rm(pred.gpp,residual)

   return(lnlike)
}#end function
#==========================================================================================#
#==========================================================================================#
