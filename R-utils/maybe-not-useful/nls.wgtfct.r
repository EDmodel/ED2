#==========================================================================================#
#==========================================================================================#
#  Function nls.wfactor                                                                    #
#                                                                                          #
#  This function is mostly the same as the function wfct available from                    #
#     http://www.r-bloggers.com/a-weighting-function-for-nls-nlslm/                        #
#  and is intended to provide weighting factors to NLS.                                    #
#                                                                                          #
#                                                                                          #
#     The weighting function can take different variable definitions and combinations      #
# thereof:                                                                                 #
#  1.  The name of the predictor (independent) variable                                    #
#  2.  The name of the response (dependent) variable                                       #
#  3.  yhat:    the fitted values \hat{y}_i of the model                                   #
#  4.  yresid:  the residuals y_i - \hat{y}_i of the model                                 #
#  5.  xysigma: the standard deviation of the response variable as a function of the       #
#               dependent variables.                                                       #
#  6.  yysigma: similar to 5, but as a function of the response variable.                  #
#  7.  xrsigma: the standard deviation of the residuals as a function of the dependent     #
#               variables.                                                                 #
#  8.  yrsigma: similar to 7, but as a function of \hat{y}_i.                              #
#                                                                                          #
#     For fitted, resid, xrsigma and yrsigma, we fit the homoscedastic the model to        #
# estimate the fitted values and residuals, then the weights are used to fit the model     #
# again.                                                                                   #
#------------------------------------------------------------------------------------------#
nls.wgtfct <<- function(expr,nbrks=10){
   expr = deparse(expr=substitute(expr))

   #----- Create new environment. ---------------------------------------------------------#
   newEnv = new.env()
   #---------------------------------------------------------------------------------------#


   #----- Get call. -----------------------------------------------------------------------#
   sc   = sys.calls()
   fine = FALSE
   n    = 0
   while ( (! fine) && n < length(sc)){
      n      = n + 1
      mc.try = sc[[n]]
      browser()
      #----- Grab call information if this is nls. ----------------------------------------#
      if (as.character(mc.try[[1]]) %in% "nls"){
         fine = TRUE
         mc   = mc.try
      }#end if (as.character(mc.try[[1]]) %in% "nls")
      #------------------------------------------------------------------------------------#
   }#end while ( (! fine) && n < length(sc))
   #---------------------------------------------------------------------------------------#

   #----- Stop if this hasn't been called from nls. ---------------------------------------#
   if (! fine) stop(" nls.wgtfct must be called from nls!")
   #---------------------------------------------------------------------------------------#


   #----- Convert mc to a list. -----------------------------------------------------------#
   mc.list = as.list(mc)
   #---------------------------------------------------------------------------------------#


   #----- Get data and write to newEnv. ---------------------------------------------------#
   data.yo    = mc.list[["data"]]
   data.yo    = eval(data.yo)
   data.list  = as.list(data.yo)
   names.data = names(data.list)
   for (i in seq_along(data.list)){
      assign(x=names.data[i],value=data.list[[i]], envir = newEnv)
   }#end for (i in seq.along(data.list))
   #---------------------------------------------------------------------------------------#


   #----- Get parameter, response and predictor names. ------------------------------------#
   mc.formula = try(as.formula(mc.list$formula),silent=TRUE)
   if ("try-error" %in% is(mc.formula)){
      mc.formula = as.formula(eval(mc.list$formula))
   }#end if
   mc.vars = all.vars(mc.formula)
   resp    = mc.vars[ 1]
   rhs     = mc.vars[-1]
   pred    = match(rhs, names(data.list))
   pred    = names(data.list)[na.omit(pred)]
   npred   = length(pred)
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #    Make sure nbrks doesn't exceed a reasonable number.                                #
   #---------------------------------------------------------------------------------------#
   x.nbrks = min(nbrks,max(3,floor((nrow(data.now)/20)^(1/npred))))
   y.nbrks = min(nbrks,max(3,floor(nrow(data.now)/20)))
   #---------------------------------------------------------------------------------------#



   #----- Retrieve x and y. ---------------------------------------------------------------#
   x = data.list [pred]
   y = data.list[[resp]]
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Calculate the local standard variation of the response variable as a function of #
   # the predictors.                                                                       #
   #---------------------------------------------------------------------------------------#
   if (length(grep("xysigma", expr)) > 0) {
      #----- Split x into blocks, then estimate y by x classes. ---------------------------#
      xbrks       = lapply( X     = x
                          , FUN   = quantile
                          , probs = seq(from=0,to=1,length.out=x.nbrks)
                          , na.rm = TRUE
                          )#end lapply
      xcut        = mapply( FUN      = cutidx
                          , x        = x
                          , breaks   = xbrks
                          , MoreArgs = list(include.lowest=TRUE)
                          , SIMPLIFY = FALSE
                          )#end mapply
      xysigma.tab = tapply(X=y,INDEX=xcut,FUN=sd,na.rm=TRUE)
      idx         = mapply(FUN=c,xcut)
      xysigma     = c(xysigma.tab[idx])
      #------------------------------------------------------------------------------------#


      #----- Append variable to the temporary environment. --------------------------------#
      xysigma = ifelse(xysigma %>% 0 ,xysigma,NA)
      assign(x="xysigma",value=xysigma,envir=newEnv)
      #------------------------------------------------------------------------------------#
   }#end if (length(grep("xysigma", expr)) > 0)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Calculate the local standard variation of the response variable as a function of #
   # the response variable itself.                                                         #
   #---------------------------------------------------------------------------------------#
   if (length(grep("yysigma", expr)) > 0) {
      #----- Split y into blocks, then estimate sigma y for each class. -------------------#
      ybrks       = quantile(x=y,probs=seq(from=0,to=1,length.out=y.nbrks),na.rm=TRUE)
      ycut        = cutidx(x=y,breaks=ybrks,include.lowest=TRUE)
      yysigma.tab = tapply(X=y,INDEX=ycut,FUN=sd,na.rm=TRUE)
      yysigma     = ifelse(is.na(ycut),NA,yysigma.tab[ycut])
      #------------------------------------------------------------------------------------#


      #----- Append variable to the temporary environment. --------------------------------#
      yysigma = ifelse(yysigma %>% 0 ,yysigma,NA)
      assign(x="yysigma",value=yysigma,envir=newEnv)
      #------------------------------------------------------------------------------------#
   }#end if (length(grep("yysigma", expr)) > 0)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Calculate fitted, residual, and residual standard deviation as a function of      #
   # predicting values if at least one of "fitted","resid" or "rsigma" is in expression.   #
   # Since they all come from model fit, write all of them to newEnv.                      #
   #---------------------------------------------------------------------------------------#
   find.model = mapply( FUN      = grepl
                      , pattern  = c("fitted","resid","xrsigma","yrsigma")
                      , MoreArgs = list(x=expr)
                      )#end mapply
   if (any(find.model)){
      mc.w0         = mc
      mc.w0$weights = NULL
      model.w0      = eval(mc.w0)
      yhat          = predict  (model.w0,newdata=data.yo)
      yresid        = y - yhat
      #------------------------------------------------------------------------------------#



      #----- Append variables to the temporary environment. -------------------------------#
      assign(x="yhat"  , value=yhat  , envir=newEnv)
      assign(x="yresid", value=yresid, envir=newEnv)
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #      Calculate the local standard variation of the response variable as a function #
      # of the predictors.                                                                 #
      #------------------------------------------------------------------------------------#
      if (length(grep("xrsigma", expr)) > 0) {
         #----- Split x into blocks, then estimate y by x classes. ------------------------#
         xbrks       = lapply( X     = x
                             , FUN   = quantile
                             , probs = seq(from=0,to=1,length.out=x.nbrks)
                             , na.rm = TRUE
                             )#end lapply
         xcut        = mapply( FUN      = cutidx
                             , x        = x
                             , breaks   = xbrks
                             , MoreArgs = list(include.lowest=TRUE)
                             , SIMPLIFY = FALSE
                             )#end mapply
         xrsigma.tab = tapply(X=yresid,INDEX=xcut,FUN=sd,na.rm=TRUE)
         idx         = mapply(FUN=c,xcut)
         xrsigma     = c(xrsigma.tab[idx])
         #---------------------------------------------------------------------------------#


         #----- Append variable to the temporary environment. -----------------------------#
         xrsigma = ifelse(xrsigma %>% 0 ,xrsigma,NA)
         assign(x="xrsigma",value=xrsigma,envir=newEnv)
         #---------------------------------------------------------------------------------#
      }#end if (length(grep("xrsigma", expr)) > 0)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Calculate the local standard variation of the response variable as a function #
      # of the predicted variable.                                                         #
      #------------------------------------------------------------------------------------#
      if (length(grep("yrsigma", expr)) > 0) {
         #----- Split y into blocks, then estimate sigma y for each class. ----------------#
         ybrks       = quantile(x=yhat,probs=seq(from=0,to=1,length.out=y.nbrks),na.rm=TRUE)
         ycut        = cutidx(x=yhat,breaks=ybrks,include.lowest=TRUE)
         yrsigma.tab = tapply(X=yresid,INDEX=ycut,FUN=sd,na.rm=TRUE)
         yrsigma     = ifelse(is.na(ycut),NA,yrsigma.tab[ycut])
         #---------------------------------------------------------------------------------#


         #----- Append variable to the temporary environment. -----------------------------#
         yrsigma = ifelse(yrsigma %>% 0 ,yrsigma,NA)
         assign(x="yrsigma",value=yrsigma,envir=newEnv)
         #---------------------------------------------------------------------------------#
      }#end if (length(grep("yrsigma", expr)) > 0)
      #------------------------------------------------------------------------------------#
   }#end if (any(find.model))
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Evaluate the weights using the provided expression.  There shall be no undefined  #
   # values or negative numbers, if they somehow exist, replace them by zero.              #
   #---------------------------------------------------------------------------------------#
   ans = eval(expr=parse(text = expr), envir = newEnv)
   ans = ifelse(ans %>% 0, ans, 0)
   ans = c(ans)
   names(ans) = NULL
   #---------------------------------------------------------------------------------------#



   #----- Return the weight vector. -------------------------------------------------------#
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function(nls.wgtfct)
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#  cutidx -- This function decides whether "cut" is necessary or not.  It will             #
#            apply cut only if the variable is _not_ categorical, otherwise it             #
#            returns the numeric version of the categorical variable.                      #
#------------------------------------------------------------------------------------------#
cutidx <<- function(x,...){
   if (is.character(x) || is.logical(x) || is.integer(x)){
      xcut = as.numeric(as.factor(x))
   }else if (is.numeric(x)){
      xcut = as.numeric(cut(x,...))
   }else{
      stop(paste0(" x type (",typeof(x),") is unexpected!"))
   }#end if
   return(xcut)
}#end cutidx
#------------------------------------------------------------------------------------------#
