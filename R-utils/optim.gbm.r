#==========================================================================================#
#==========================================================================================#
#    Function optim.randomforest.  This function is a wrapper for RandomForest that also   #
# creates bootstrap realisations for cross validation.                                     #
#------------------------------------------------------------------------------------------#
optim.gbm <<- function( formula
                      , data
                      , sy.data  = NULL
                      , n.boot   = 1000
                      , ci.level = 0.95
                      , verbose  = FALSE
                      , n.syobs  = 10000
                      , yrdm.min = -Inf
                      , yrdm.max = +Inf
                      ,...
                      ){
   #----- Save number of data points. -----------------------------------------------------#
   n.data    = nrow(data)
   #---------------------------------------------------------------------------------------#


   #----- Save the name of the response variable. -----------------------------------------#
   yname     = all.vars(formula)[1]
   #---------------------------------------------------------------------------------------#


   #----- Run RandomForest. ---------------------------------------------------------------#
   if (verbose) cat0( "             > Generalised Boosted Model (full model)")
   dotdotdot = list(formula=formula,data=data,verbose=FALSE,...)
   ans       = do.call(what="gbm",args=dotdotdot)
   ans       = gbm(formula=formula,data=data,verbose=FALSE,...)
   #---------------------------------------------------------------------------------------#


   #----- Initialise the cross validation matrix. -----------------------------------------#
   xval.boot = matrix( data     = NA
                     , nrow     = n.data
                     , ncol     = n.boot
                     , dimnames = list(rownames(data),NULL)
                     )#end matrix
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Get the ellipsis elements, and force the bootstrap iterations to be simpler.      #
   #---------------------------------------------------------------------------------------#
   dotdotdot = list(formula=formula,...)
   dotdotdot = modifyList( x   = dotdotdot
                         , val = list( cv.folds  = 0
                                     , verbose   = FALSE
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
   #     Loop until we reach the sought number of bootstrap realisations.  In case some    #
   # realisation doesn't work, skip it.                                                    #
   #---------------------------------------------------------------------------------------#
   ib = 0
   while (ib < n.boot){
      #----- Select samples for this realisation. -----------------------------------------#
      idx         = sample.int(n=n.data,replace=TRUE)
      ixval       = which(! (sequence(n.data) %in% idx))
      boot.data   = data[  idx,,drop=FALSE]
      xval.data   = data[ixval,,drop=FALSE]
      #------------------------------------------------------------------------------------#


      #----- Call GBM. --------------------------------------------------------------------#
      dotnow   = modifyList(x=dotdotdot,val=list(data=boot.data))
      gbm.now  = try(do.call(what="gbm",args=dotnow),silent=TRUE)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Check whether to append to the data set.                                      #
      #------------------------------------------------------------------------------------#
      if (! ("try-error" %in% is(gbm.now))){
         ib = ib + 1


         #----- Run cross validation. -----------------------------------------------------#
         if (length(ixval) > 0){
            #----- Number of trees used by predict.gbm. -----------------------------------#
            if (gbm.now$train.fraction < 1) {
               n.trees = gbm.perf(gbm.now,method="test",plot.it=FALSE)
            }else if (! is.null(gbm.now$cv.error)) {
               n.trees = gbm.perf(gbm.now,method="cv"  ,plot.it=FALSE)
            }else{
               n.trees = length(gbm.now$train.error)
            }#end if (gbm.now$train.fraction < 1)
            #------------------------------------------------------------------------------#

            #----- Predict values that were left out. -------------------------------------#
            ypred               = predict( object  = gbm.now
                                         , newdata = xval.data
                                         , n.trees = n.trees
                                         , verbose = FALSE
                                         )#end ypred
            xval.boot[ixval,ib] = ypred
            #------------------------------------------------------------------------------#
         }#end if
         #---------------------------------------------------------------------------------#

         if (verbose){
            cat0( "             > Bootstrap  (Generalised Boosted Model cross validation)"
                ,               "; iteration: ",ib
                ,               "; # of trees: ",n.trees
                )#end cat0
         }#end if (verbose)
         #---------------------------------------------------------------------------------#
      }else if (verbose){
         cat0("             > Bootstrap realisation failed, skip it.")
      }#end if (! ("try-error" %in% is(gbm.now)))
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
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Loop through all simulations to account for the observation errors.               #
   #---------------------------------------------------------------------------------------#
   if (! is.null(skip.syobs)){



      #----- Initialise the cross validation object and matrix. ---------------------------#
      keep.syobs = rep(TRUE,times=n.syobs)
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


         #----- Call GBM. -----------------------------------------------------------------#
         dotnow  = modifyList(x=dotdotdot,val=list(data=sim.data))
         gbm.now = try(do.call(what="gbm",args=dotnow),silent=TRUE)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Check whether to append to the data set.                                   #
         #---------------------------------------------------------------------------------#
         if (! ("try-error" %in% is(gbm.now))){


            #----- Number of trees used by predict.gbm. -----------------------------------#
            if (gbm.now$train.fraction < 1) {
               n.trees = gbm.perf(gbm.now,method="test",plot.it=FALSE)
            }else if (! is.null(gbm.now$cv.error)) {
               n.trees = gbm.perf(gbm.now,method="cv"  ,plot.it=FALSE)
            }else{
               n.trees = length(gbm.now$train.error)
            }#end if (gbm.now$train.fraction < 1)
            #------------------------------------------------------------------------------#

            #----- Run cross validation. --------------------------------------------------#
            ypred            = predict( object  = gbm.now
                                      , newdata = sim.data
                                      , n.trees = n.trees
                                      , verbose = FALSE
                                      )#end predict
            eval.syobs[,isy] = ypred
            #------------------------------------------------------------------------------#


            #----- Show banner to entertain the bored user. -------------------------------#
            if (verbose) cat0( "             > Sigma-y  (GBM); iteration: ",isy)
            #------------------------------------------------------------------------------#
         }else{
            #----- Slate this iteration to be deleted. ------------------------------------#
            keep.syobs[isy] = FALSE
            #------------------------------------------------------------------------------#

            #----- Show banner to entertain the bored user. -------------------------------#
            if (verbose) cat0("             > Sigma-y realisation failed, skip it.")
            #------------------------------------------------------------------------------#
         }#end if (! ("try-error" %in% is(gbm.now)))
         #---------------------------------------------------------------------------------#


         #----- Keep only successful steps. -----------------------------------------------#
         ans$sy.data    = sy.data   [,keep.syobs,drop=FALSE]
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


   #----- Return generalised boosted model object with the cross-validation attached. -----#
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end optim.gbm
#==========================================================================================#
#==========================================================================================#
