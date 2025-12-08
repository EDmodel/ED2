#==========================================================================================#
#==========================================================================================#
#     This function provides a summary for all numeric elements of a data frame or list    #
# (or matrix columns).                                                                     #
#                                                                                          #
#  + Input:                                                                                #
#    - x           -- Input data set.  Lists, data frames or matrices.                     #
#    - keep.nn     -- Keep the non-numeric elements? (TRUE/FALSE).                         #
#    - byrow       -- In case x is a matrix, should the summary be applied on rows as      #
#                     opposed to columns.                                                  #
#    - finite.only -- Use finite values only?                                              #
#    - neverlog    -- Variable names that should never be log (partial match is fine too). #
#------------------------------------------------------------------------------------------#
summnum <<- function(x,byrow=FALSE,finite.only=TRUE,neverlog=NULL,is.debug=FALSE){
   #----- In case x is matrix, turn it into a data frame. ---------------------------------#
   if (is.matrix(x)){
      if (byrow) x = t(x)
      x = as.data.frame(x,stringsAsFactors=FALSE)
   }else if (is.vector(x)){
      x = data.frame(x,stringsAsFactors=FALSE)
   }else if (! (is.list(x) || is.data.frame(x))){
      stop (" 'x' must be a data frame, list, matrix, or vector!")
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Apply the summary function. -----------------------------------------------------#
   if (is.debug){
      ans = NULL
      for (a in seq_along(x)){
         nmnow   = names(x)[[a]]
         ans.now = try(summnum.int(x=x[[a]],finite.only=finite.only),silent=TRUE)
         if ("try-error" %in% is(ans.now)){
            cat0("Problems with summnum.int.")
            cat0("Variable: ",nmnow,".")
            browser()
         }else if (is.null(ans)){
            ans          = data.frame(x=ans.now)
            names(ans)   = nmnow
         }else{
            ans[[nmnow]] = ans.now
         }#end if ("try-error" %in% is(ans.now))
      }#end for (a in seq_along(x))
   }else{
      ans = sapply(X=x,FUN=summnum.int,finite.only=finite.only)
   }#end if
   ans = data.frame(t(ans),stringsAsFactors=FALSE)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Determine whether log scale is safe.                                               #
   #---------------------------------------------------------------------------------------#
   opt.orig = options()
   options(warn=-1)
      #----- First, check whether the variable is allowed to be in logarithmic scale. -----#
      log.allowed = rep(TRUE,times=ncol(x))
      for (n in seq_along(neverlog)){
         pattnow     = neverlog[n]
         forbidden   = grepl(pattern=pattnow,x=names(x),ignore.case=TRUE)
         log.allowed = log.allowed & ! forbidden
      }#end for
      #------------------------------------------------------------------------------------#



      #----- Second, ensure that zeroes are relatively rare. ------------------------------#
      few.zeroes = ans$min %ge% 0. & ans$q100 %gt% 0
      #------------------------------------------------------------------------------------#



      #------ Third, check that the log-normal distribution is better. --------------------#
      if (is.debug){
         lnlike = NULL
         for (a in seq_along(x)){
            nmnow   = names(x)[[a]]
            ans.now = try(c(unlist(lnlike.comp(x=x[[a]]))),silent=TRUE)
            if ("try-error" %in% is(ans.now)){
               cat0("Problems with lnlike.comp.")
               cat0("Variable: ",nmnow,".")
               browser()
            }else{
               lnlike          = rbind(lnlike,ans.now)
            }#end if ("try-error" %in% is(ans.now))
         }#end for (a in seq_along(x))
         rownames(lnlike) = names(x)
      }else{
         lnlike = sapply(X=x, FUN=lnlike.comp)
      }#end if (is.debug)
      lnlike            = data.frame(t(lnlike),stringsAsFactors=FALSE)
      ln.lnorm          = with( lnlike
                              , ifelse( test = is.finite(lnorm)
                                      , yes  = signif(lnorm,3)
                                      , no   = NA_real_
                                      )#end ifelse
                              )#end with
      ln.norm           = with( lnlike
                              , ifelse( test = is.finite(norm)
                                      , yes  = signif(norm,3)
                                      , no   = NA_real_
                                      )#end ifelse
                              )#end with
      lognorm.is.better = ln.lnorm %ge% ln.norm
      #------------------------------------------------------------------------------------#


      #------ Log is safe to use only if all three conditions are met. --------------------#
      ans$lnlike.norm   = lnlike$norm
      ans$lnlike.lnorm  = lnlike$lnorm
      ans$safelog   = log.allowed & few.zeroes & lognorm.is.better
      #------------------------------------------------------------------------------------#

   options(opt.orig)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #   Return the data frame.                                                              #
   #---------------------------------------------------------------------------------------#
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     This function applies the summary to each element, returning NA in case the variable #
# is not numeric.  This shouldn't be called directly.                                      #
#------------------------------------------------------------------------------------------#
summnum.int <<- function(x,finite.only){
   #----- Make sure this function has been called by summnum.  Otherwise, stop. -----------#
   patt  = "^([A-Za-z0-9]+)(\\({1})(.*)(\\){1})$"
   repl  = "\\1"
   wcm.3 = try(gsub(pattern=patt,replacement=repl,x=deparse(sys.call(-3))),silent=TRUE)
   wcm.2 = try(gsub(pattern=patt,replacement=repl,x=deparse(sys.call(-2))),silent=TRUE)
   wcm.1 = try(gsub(pattern=patt,replacement=repl,x=deparse(sys.call(-1))),silent=TRUE)
   if ("try-error" %in% is(wcm.3)){wcm.3 = NA}else if (is.null(wcm.3)){wcm.3 = NA}
   if ("try-error" %in% is(wcm.2)){wcm.2 = NA}else if (is.null(wcm.2)){wcm.2 = NA}
   if ("try-error" %in% is(wcm.1)){wcm.1 = NA}else if (is.null(wcm.1)){wcm.1 = NA}

   #if (! all(c(wcm.1,wcm.2,wcm.3) %eq% c("lapply","sapply","summnum"))){
   #   stop(" Function summnum.int is internal, and must be called through summnum","\n")
   #}#end if
   #---------------------------------------------------------------------------------------#


   #---- Disable warning. -----------------------------------------------------------------#
   opt.orig = options()
   options(warn=-1)
   #---------------------------------------------------------------------------------------#


   #----- Make sure x is a simple vector. -------------------------------------------------#
   x  = unlist(c(x))
   if (finite.only) x[! is.finite(x)] = NA_real_
   nx = length(x)
   #---------------------------------------------------------------------------------------#


   #----- Create a summary.  If this is not numeric, create a similar vector with NA. -----#
   if (is.numeric(x)){
      ans = c( min    = min(x,na.rm=TRUE)
             , q025   = quantile(x,probs=0.025,na.rm=TRUE)
             , q100   = quantile(x,probs=0.100,na.rm=TRUE)
             , q250   = quantile(x,probs=0.250,na.rm=TRUE)
             , median = quantile(x,probs=0.500,na.rm=TRUE)
             , q750   = quantile(x,probs=0.750,na.rm=TRUE)
             , q900   = quantile(x,probs=0.900,na.rm=TRUE)
             , q975   = quantile(x,probs=0.975,na.rm=TRUE)
             , max    = max (x,na.rm=TRUE)
             , mean   = mean(x,na.rm=TRUE)
             , sdev   = sd  (x,na.rm=TRUE)
             , skew   = skew(x,na.rm=TRUE)
             , kurt   = kurt(x,na.rm=TRUE)
             , navl   = sum (! is.na(x))
             , npos   = sum (x %lt% 0)
             , ntot   = nx
             )#end c
   }else{
      ans = c( min    = NA_real_
             , q025   = NA_real_
             , q100   = NA_real_
             , q250   = NA_real_
             , median = NA_real_
             , q750   = NA_real_
             , q900   = NA_real_
             , q975   = NA_real_
             , max    = NA_real_
             , mean   = NA_real_
             , sdev   = NA_real_
             , skew   = NA_real_
             , kurt   = NA_real_
             , navl   = nx
             , npos   = nx
             , ntot   = nx
             )#end c
   }#end if
   names(ans) = c("min","q025","q100","q250","median","q750","q900","q975","max"
                 ,"mean","sdev","skew","kurt","miss","npos","ntot")
   #---------------------------------------------------------------------------------------#



   #----- Make sure the moments are finite. -----------------------------------------------#
   mfin = c("mean","sdev","skew","kurt")
   ans[mfin] = ifelse(test=is.finite(ans[mfin]),yes=ans[mfin],no=NA_real_)
   #---------------------------------------------------------------------------------------#



   #----- Name all columns, re-organise the vector, and return. ---------------------------#
   options(opt.orig)
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     This function returns TRUE if the log-normal distribution has better support than    #
# the normal distribution.                                                                 #
#------------------------------------------------------------------------------------------#
lnlike.comp <<- function(x,xlab="nothing"){
   #----- Make sure this function has been called by summnum.  Otherwise, stop. -----------#
   patt  = "^([A-Za-z0-9]+)(\\({1})(.*)(\\){1})$"
   repl  = "\\1"
   wcm.3 = try(gsub(pattern=patt,replacement=repl,x=deparse(sys.call(-3))),silent=TRUE)
   wcm.2 = try(gsub(pattern=patt,replacement=repl,x=deparse(sys.call(-2))),silent=TRUE)
   wcm.1 = try(gsub(pattern=patt,replacement=repl,x=deparse(sys.call(-1))),silent=TRUE)
   if ("try-error" %in% is(wcm.3)){wcm.3 = NA}else if (is.null(wcm.3)){wcm.3 = NA}
   if ("try-error" %in% is(wcm.2)){wcm.2 = NA}else if (is.null(wcm.2)){wcm.2 = NA}
   if ("try-error" %in% is(wcm.1)){wcm.1 = NA}else if (is.null(wcm.1)){wcm.1 = NA}

   #if (! all(c(wcm.1,wcm.2,wcm.3) %eq% c("lapply","sapply","summnum"))){
   #   stop(" Function sw.pvalue is internal and must be called through summnum","\n")
   #}#end if
   #---------------------------------------------------------------------------------------#


   #---- Disable warning. -----------------------------------------------------------------#
   opt.orig = options()
   options(warn=-1)
   #---------------------------------------------------------------------------------------#


   #----- Make sure x is a simple vector.  Keep only the positive terms. ------------------#
   x   = as.numeric(unlist(c(x)))
   lnx = ifelse(test=x %gt% 0,yes=log(x),no=NA_real_)
   sel = is.finite(lnx)
   x   = x  [sel]
   lnx = lnx[sel]
   #---------------------------------------------------------------------------------------#


   #----- Obtain the p-values. In case there aren't enough points, return NA. -------------#
   if (length(x) > 3){
      lnlike.norm   = fitdistr(x=x,densfun="normal")$loglik
      lnlike.lnorm  = fitdistr(x=x,densfun="lognormal")$loglik
   }else{
      #----- Not enough points, don't use log. --------------------------------------------#
      lnlike.norm   = NA_real_
      lnlike.lnorm  = NA_real_
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Return answer. ------------------------------------------------------------------#
   ans = c(lnlike.norm,lnlike.lnorm)
   names(ans) = c("norm","lnorm")
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function
#==========================================================================================#
#==========================================================================================#
