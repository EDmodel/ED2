#==========================================================================================#
#==========================================================================================#
#    Function optim.pls.  This function is a wrapper for plsr that also creates bootstrap  #
# realisations for cross validation.  It also has an option for log-transformation that    #
# back-transforms the results.                                                             #
#------------------------------------------------------------------------------------------#
optim.pls <<- function( formula
                          , data
                          , ylog          = FALSE
                          , xlog.list     = character(0)
                          , fve.tolerance = 0.99
                          , sy.data       = NULL
                          , n.boot        = 1000
                          , ci.level      = 0.95
                          , verbose       = FALSE
                          , n.syobs       = 10000
                          , boot.class    = NULL
                          , yrdm.min      = -Inf
                          , yrdm.max      = +Inf
                          ,...
                          ){
   #----- Save number of data points. -----------------------------------------------------#
   n.data    = nrow(data)
   #---------------------------------------------------------------------------------------#



   #----- Save the name of the response variable. -----------------------------------------#
   yname     = all.vars(formula)[ 1]
   xname     = all.vars(formula)[-1]
   if (all(xname %in% ".")){
      xname = names(data)
      xname = xname[! (xname %in% yname)]
   }#end if (all(xname) %in% ".")
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     If the LHS should be log-transformed for fitting, replace formula and variables.  #
   # We do this to ensure the results are back-transformed.                                #
   #---------------------------------------------------------------------------------------#
   orig.data = data
   if (ylog){
      lhs.use          = paste0("log(",yname,")")
   }else{
      lhs.use          = yname
   }#end if (ylog)
   #---------------------------------------------------------------------------------------#




   #----- Also, transform predictors that should be transformed. --------------------------#
   rhs.use   = ifelse( test = xname %in% xlog.list
                     , yes  = paste0("log(",xname,")")
                     , no   = xname
                     )#end ifelse
   rhs.use   = paste0(rhs.use,collapse=" + ")
   form.now  = as.formula(paste0(lhs.use," ~ ",rhs.use))
   #---------------------------------------------------------------------------------------#






   #----- Run PLS. ------------------------------------------------------------------------#
   if (verbose) cat0( "             > PLS (full model)")
   ans              = plsr(formula=form.now,data=data,verbose=verbose,...)
   ans$orig.formula = formula
   ans$formula      = form.now
   ans$yname        = yname
   ans$xname        = xname
   ans$xlog         = xname %in% xlog.list
   ans$ylog         = ylog
   #---------------------------------------------------------------------------------------#



   #----- Retrieve MSE in case back-transformation is sought. -----------------------------#
   resnow     = ans$residuals[,1,]
   rsumsq     = apply(X=resnow,MARGIN=2,FUN=sum2,na.rm=TRUE)
   nnow       = apply(X=resnow,MARGIN=2,FUN=function(x) sum(is.finite(x)))
   msenow     = MSEP(object=ans,"adjCV")
   ans$mse    = c(msenow$val[1,1,-1])
   ans$vary   = if(ylog){var(log(data[yname]),na.rm=TRUE)}else{var(data[yname],na.rm=TRUE)}
   ans$fve    = 1 - ans$mse / ans$vary
   max.fve    = max(ans$fve,na.rm=TRUE)
   rel.fve    = ans$fve / max.fve
   ans$ncomp  = min(which(rel.fve >= fve.tolerance))
   #---------------------------------------------------------------------------------------#





   #----- Save fitted values. -------------------------------------------------------------#
   iy         = which(names(data) %in% yname)
   ans$data   = data[,-iy,drop=FALSE]
   ans$y      = data[, iy,drop=TRUE ]
   names(ans$y) = rownames(data)
   ypred      = predict(object=ans,ncomp=ans$ncomp)
   ypred      = ypred[,1,1]
   if (ylog){
      ans$ln.pred   = ypred
      ans$ln.sigma  = sqrt(ans$mse[ans$ncomp])
      ans$predicted = exp(ypred + 0.5 * ans$mse[ans$ncomp])
      ans$sigma     = with(ans,sqrt(exp(mse[ncomp])-1)*exp(2.*ans$ln.pred + mse[ncomp]))
      plwr          = 0.5 * (1.0 - ci.level)
      pupr          = 0.5 * (1.0 + ci.level)
      ans$qlow      = qlnorm(p=plwr,meanlog=ans$ln.pred,sdlog=ans$ln.sigma)
      ans$qhigh     = qlnorm(p=pupr,meanlog=ans$ln.pred,sdlog=ans$ln.sigma)
   }else{
      ans$predicted = ypred
      ans$sigma     = sqrt(ans$mse)
      ans$qlow      = qnorm(p=plwr,mean=ans$predicted,sd=ans$sigma[ans$ncomp])
      ans$qhigh     = qnorm(p=pupr,mean=ans$predicted,sd=ans$sigma[ans$ncomp])
   }#end if (ylog)
   #---------------------------------------------------------------------------------------#



   #----- Initialise the cross validation object and matrix. ------------------------------#
   xval.boot = matrix( data     = NA_real_
                     , nrow     = n.data
                     , ncol     = n.boot
                     , dimnames = list(rownames(data),NULL)
                     )#end matrix
   pls.boot  = replicate(n=n.boot,expr=list())
   mse.boot  = rep(x=NA_real_,times=n.boot)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Get the ellipsis elements, and force the bootstrap iterations to be simpler.      #
   #---------------------------------------------------------------------------------------#
   dotdotdot = list(formula = form.now,...)
   dotdotdot = modifyList( x   = dotdotdot
                         , val = list( verbose     = FALSE
                                     , model       = FALSE
                                     , ncomp       = ans$ncomp
                                     )#end list
                         )#end modifyList
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Set data for evaluating the impact of the uncertainty in the reference           #
   # response variable on model predictions.                                               #
   #---------------------------------------------------------------------------------------#
   if (is.null(sy.data)){
      #----- Skip evaluation of the impact of errors in observations on predictions. ------#
      skip.syobs = TRUE
      n.syobs    = 0
      #------------------------------------------------------------------------------------#
   }else{
      #----- Find the impact of errors in observations on predictions. --------------------#
      skip.syobs = FALSE
      #------------------------------------------------------------------------------------#
   
      #----- Coerce sy.data to matrix. ----------------------------------------------------#
      if (is.list(sy.data)){
         #----- Transform list in data frame before converting it to matrix. --------------#
         sy.data = try(as.matrix(data.frame(sy.data)),silent=TRUE)
         #---------------------------------------------------------------------------------#
      }else if (! is.matrix(sy.data)){
         #----- Transform list in data frame before converting it to matrix. --------------#
         sy.data = try(as.matrix(sy.data),silent=TRUE)
         #---------------------------------------------------------------------------------#
      }#end if (is.list(sy.data))
      #------------------------------------------------------------------------------------#


      #----- Make sure the transformation was successful. ---------------------------------#
      if ("try-error" %in% is(sy.data)){
         stop("Argument sy.data cannot be coerced into a matrix.")
      }else if (nrow(sy.data) != nrow(data)){
         stop("Number of entries of 'sy.data' must match 'data'.")
      }#end if
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Check what to do depending on the size of sy.data.                             #
      #------------------------------------------------------------------------------------#
      #----- Check the number of rows, make sure they match data. -------------------------#
      if (ncol(sy.data) == 1){
         y.eval  = rnorm( n    = nrow(data)*n.syobs
                        , mean = rep(x=data[[yname]],times=n.syobs)
                        , sd   = rep(x=sy.data[,1]  ,times=n.syobs)
                        )#
         sy.data = matrix( data     = pmin(yrdm.max,pmax(yrdm.min,y.eval))
                         , nrow     = nrow(data)
                         , ncol     = n.syobs
                         , dimnames = list(rownames(data),NULL)
                         )#end matrix
         rm(y.eval)
      }else{
         sy.data = 0.*sy.data + pmin(yrdm.max,pmax(yrdm.min,sy.data))
         n.syobs = ncol(sy.data)
      }#end if (ncol(sy.data) == 1)
      #------------------------------------------------------------------------------------#
   }#end if (is.null(sy.data))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Decide the classes for bootstrap.                                                #
   #---------------------------------------------------------------------------------------#
   if (! is.null(boot.class)){
      uniq.class   = sort(unique(boot.class))
      n.uniq.class = length(uniq.class)
   }#end if (! is.null (boot.class))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Loop until we reach the sought number of bootstrap realisations.  In case some    #
   # realisation doesn't work, skip it.                                                    #
   #---------------------------------------------------------------------------------------#
   ib = 0
   while (ib < n.boot){
      #----- Select samples for this realisation. -----------------------------------------#
      if (is.null (boot.class)){
         idx         = sample.int(n=n.data,replace=TRUE)
         ixval       = which(! (sequence(n.data) %in% idx))
      }else{
         use.class   = lit.sample(x=uniq.class,size=n.uniq.class,replace=TRUE)
         use.sample  = mapply( FUN      = function(x,y) which(y %in% x)
                             , x        = use.class
                             , MoreArgs = list(y=boot.class)
                             )#end mapply
         use.sample  = c(unlist(use.sample))
         idx         = lit.sample(x=use.sample,size=n.data,replace=TRUE)
         ixval       = which(! boot.class %in% use.class)
      }#end if
      #------------------------------------------------------------------------------------#


      #----- Set training and testing data sets. ------------------------------------------#
      boot.data   = data[  idx,,drop=FALSE]
      xval.data   = data[ixval,,drop=FALSE]
      #------------------------------------------------------------------------------------#


      #----- Call PLS. --------------------------------------------------------------------#
      dotnow   = modifyList(x=dotdotdot,val=list(data=boot.data))
      pls.now = try(do.call(what="plsr",args=dotnow),silent=TRUE)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Check whether to append to the data set.                                      #
      #------------------------------------------------------------------------------------#
      if (! ("try-error" %in% is(pls.now))){
         ib       = ib + 1

         #----- Retrieve MSE in case back-transformation is sought. -----------------------#
         mse.now = MSEP(object=pls.now,"adjCV")
         mse.now = c(mse.now$val[1,1,ans$ncomp+1])
         #---------------------------------------------------------------------------------#



         #----- Run cross validation. -----------------------------------------------------#
         if (length(ixval) > 0){
            ypred = predict(object=pls.now,newdata=xval.data,ncomp=ans$ncomp)
            ypred = ypred[,1,1]
            if (ylog){
               xval.boot[ixval,ib] = exp(ypred + 0.5 * mse.now)
            }else{
               xval.boot[ixval,ib] = ypred
            }#end if (ylog)
         }#end if
         mse.boot[[ib]] = mse.now
         #---------------------------------------------------------------------------------#

         if (verbose){
            cat0( "             > Bootstrap  (PLS cross validation)"
                ,               "; iteration: ",ib
                )#end cat0
         }#end if (verbose)
         #---------------------------------------------------------------------------------#
      }else if (verbose){
         cat0("             > Bootstrap  realisation failed, skip it.")
      }#end if (! ("try-error" %in% is(pls.now)))
      #------------------------------------------------------------------------------------#
   }#end while (ib < n.boot)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find cross validation statistics.                                                 #
   #---------------------------------------------------------------------------------------#
   qlow        = 0.5 - 0.5 * ci.level
   qhigh       = 0.5 + 0.5 * ci.level
   xval.boot   = ifelse(is.finite(xval.boot),xval.boot,NA)
   xval.resid  = - apply(X=xval.boot ,MARGIN=2,FUN="-", data[[yname]])
   n.xval      =   apply(X=xval.boot ,MARGIN=1,FUN=function(x) sum(! is.na(x)))
   expect.xval =   apply(X=xval.boot ,MARGIN=1,FUN=mean                ,na.rm=TRUE)
   qlow.xval   =   apply(X=xval.boot ,MARGIN=1,FUN=quantile,probs=qlow ,na.rm=TRUE)
   qhigh.xval  =   apply(X=xval.boot ,MARGIN=1,FUN=quantile,probs=qhigh,na.rm=TRUE)
   bias.xval   = - apply(X=xval.resid,MARGIN=1,FUN=mean                ,na.rm=TRUE)
   sigma.xval  =   apply(X=xval.resid,MARGIN=1,FUN=sd                  ,na.rm=TRUE)
   rmse.xval   =   sqrt(bias.xval^2+sigma.xval^2)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Attach to answer.                                                                 #
   #---------------------------------------------------------------------------------------#
   ans$conf.int  = ci.level
   ans$cross.val = data.frame( n        = n.xval
                             , expected = expect.xval
                             , qlow     = qlow.xval
                             , qhigh    = qhigh.xval
                             , bias     = bias.xval
                             , sigma    = sigma.xval
                             , rmse     = rmse.xval
                             )#end data.frame
   ans$xval.mat  = xval.boot
   ans$mse.boot  = mse.boot
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Loop through all simulations to account for the observation errors.               #
   #---------------------------------------------------------------------------------------#
   if (! is.null(skip.syobs)){



      #----- Initialise the cross validation object and matrix. ---------------------------#
      keep.syobs = rep(TRUE,times=n.syobs)
      mse.syobs  = rep(x=NA_real_,times=n.syobs)
      eval.syobs = matrix( data     = NA
                         , nrow     = n.data
                         , ncol     = n.syobs
                         , dimnames = list(rownames(data),NULL)
                         )#end matrix
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Loop through simulations.                                                      #
      #------------------------------------------------------------------------------------#
      for (isy in sequence(n.syobs)){
         #----- Select samples for this realisation. --------------------------------------#
         sim.data          = data
         sim.data[[yname]] = sy.data[,isy]
         #---------------------------------------------------------------------------------#


         #----- Call PLS. -----------------------------------------------------------------#
         dotnow  = modifyList(x=dotdotdot,val=list(data=sim.data))
         pls.now = try(do.call(what="plsr",args=dotnow),silent=TRUE)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Check whether to append to the data set.                                   #
         #---------------------------------------------------------------------------------#
         if (! ("try-error" %in% is(pls.now))){

            #----- Retrieve MSE in case back-transformation is sought. --------------------#
            mse.now = MSEP(object=pls.now,"adjCV")
            mse.now = c(mse.now$val[1,1,ans$ncomp+1])
            #------------------------------------------------------------------------------#



            #----- Run cross validation. --------------------------------------------------#
            ypred = predict(object=pls.now,newdata=sim.data,ncomp=ans$icomp)
            ypred = ypred[,1,1]
            if (ylog){
               eval.syobs[,isy] = exp(ypred + 0.5 * mse.now)
            }else{
               eval.syobs[,isy] = ypred
            }#end if (ylog)
            #------------------------------------------------------------------------------#


            #----- Save model and MSE. ----------------------------------------------------#
            mse.syobs[[isy]] = mse.now
            #------------------------------------------------------------------------------#


            #----- Show banner to entertain the bored user. -------------------------------#
            if (verbose) cat0( "             > Sigma-y  (PLS); iteration: ",isy)
            #------------------------------------------------------------------------------#
         }else{
            #----- Slate this iteration to be deleted. ------------------------------------#
            keep.syobs[isy] = FALSE
            #------------------------------------------------------------------------------#

            #----- Show banner to entertain the bored user. -------------------------------#
            if (verbose) cat0("             > Sigma-y realisation failed, skip it.")
            #------------------------------------------------------------------------------#
         }#end if (! ("try-error" %in% is(pls.now)))
         #---------------------------------------------------------------------------------#


         #----- Keep only successful steps. -----------------------------------------------#
         ans$sy.data    = sy.data   [,keep.syobs,drop=FALSE]
         ans$mse.syobs  = mse.syobs [ keep.syobs]
         ans$eval.syobs = eval.syobs[,keep.syobs,drop=FALSE]
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Find cross validation statistics.                                           #
         #---------------------------------------------------------------------------------#
         qlow        = 0.5 - 0.5 * ci.level
         qhigh       = 0.5 + 0.5 * ci.level
         sy.mat      = 0 * eval.syobs + ifelse(is.finite(eval.syobs),eval.syobs,NA)
         sy.resid    = - apply(X=sy.mat ,MARGIN=2,FUN="-", data[[yname]])
         n.sy        =   apply(X=sy.mat ,MARGIN=1,FUN=function(x) sum(! is.na(x)))
         expect.sy   =   apply(X=sy.mat ,MARGIN=1,FUN=mean                ,na.rm=TRUE)
         qlow.sy     =   apply(X=sy.mat ,MARGIN=1,FUN=quantile,probs=qlow ,na.rm=TRUE)
         qhigh.sy    =   apply(X=sy.mat ,MARGIN=1,FUN=quantile,probs=qhigh,na.rm=TRUE)
         bias.sy     = - apply(X=sy.resid,MARGIN=1,FUN=mean                ,na.rm=TRUE)
         sigma.sy    =   apply(X=sy.resid,MARGIN=1,FUN=sd                  ,na.rm=TRUE)
         rmse.sy     =   sqrt(bias.sy^2+sigma.sy^2)
         ans$sy.summ = data.frame( n        = n.sy
                                 , expected = expect.sy
                                 , qlow     = qlow.sy
                                 , qhigh    = qhigh.sy
                                 , bias     = bias.sy
                                 , sigma    = sigma.sy
                                 , rmse     = rmse.sy
                                 )#end data.frame
         #---------------------------------------------------------------------------------#
      }#end for (isy in sequence(n.syobs))
      #------------------------------------------------------------------------------------#
   }#end if (! is.null(skip.sydata))
   #---------------------------------------------------------------------------------------#


   #----- Return PLS object with the cross-validation attached. ---------------------------#
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end optim.pls
#==========================================================================================#
#==========================================================================================#
