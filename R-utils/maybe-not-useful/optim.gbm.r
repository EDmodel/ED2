#==========================================================================================#
#==========================================================================================#
#    Function optim.randomforest.  This function is a wrapper for RandomForest that also   #
# creates bootstrap realisations for cross validation.                                     #
#------------------------------------------------------------------------------------------#
optim.gbm <<- function(formula,data,n.boot=1000,ci.level=0.95,verbose=FALSE,...){
   #----- Save number of data points. -----------------------------------------------------#
   n.data    = nrow(data)
   #---------------------------------------------------------------------------------------#


   #----- Save the name of the response variable. -----------------------------------------#
   yname     = all.vars(formula)[1]
   #---------------------------------------------------------------------------------------#


   #----- Run RandomForest. ---------------------------------------------------------------#
   if (verbose) cat0( "             > Generalised Boosted Model (full model)")
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


      #----- Call randomForest. -----------------------------------------------------------#
      dotnow   = modifyList(x=dotdotdot,val=list(data=boot.data))
      gbm.boot = try(do.call(what="gbm",args=dotnow),silent=TRUE)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Check whether to append to the data set.                                      #
      #------------------------------------------------------------------------------------#
      if (! ("try-error" %in% is(gbm.boot))){
         ib = ib + 1

         #----- Run cross validation. -----------------------------------------------------#
         if (length(ixval) > 0){
            #----- Number of trees used by predict.gbm. -----------------------------------#
            if (gbm.boot$train.fraction < 1) {
               n.trees = gbm.perf(gbm.boot,method="test",plot.it=FALSE)
            }else if (! is.null(gbm.boot$cv.error)) {
               n.trees = gbm.perf(gbm.boot,method="cv"  ,plot.it=FALSE)
            }else{
               n.trees = length(gbm.boot$train.error)
            }#end if (gbm.boot$train.fraction < 1)
            #------------------------------------------------------------------------------#

            #----- Predict values that were left out. -------------------------------------#
            ypred               = predict( object  = gbm.boot
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
         cat0("             > Bootstrap  realisation failed, skip it.")
      }#end if (! ("try-error" %in% is(gbm.boot)))
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


   #----- Return generalised boosted model object with the cross-validation attached. -----#
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end optim.gbm
#==========================================================================================#
#==========================================================================================#
