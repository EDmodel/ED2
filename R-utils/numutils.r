#==========================================================================================#
#==========================================================================================#
#     Function that finds the cube root, for both positive and negative numbers.           #
#------------------------------------------------------------------------------------------#
cbrt <<- function(x){
   x333      = sign(x) * abs(x)^onethird
   bad       = is.nan(x333)
   x333[bad] = NA_real_
   return(x333)
}#end function cbrt
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     Functions that finds round of the log of the number.                                 #
#------------------------------------------------------------------------------------------#
round.log   <<- function(x,base=exp(1),...) base^(round(log(x,base=base),...))
round.log10 <<- function(x,...) 10^(round(log10(x),...))
round.log2  <<- function(x,...) 2^(round(log2(x),...))
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     Base 10 exponential.                                                                 #
#------------------------------------------------------------------------------------------#
exp10 <<- function(x,...) 10^x

#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     Bound data between lower and upper limit.                                            #
#------------------------------------------------------------------------------------------#
bound <<- function(x,lwr=min(x,na.rm=TRUE),upr=max(x,na.rm=TRUE),buff=2^-23){
   #---- Return the original data in case nothing is valid. -------------------------------#
   if (! any(is.finite(x))) return(x)
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #   Don't let the bounds to be insane.                                                  #
   #---------------------------------------------------------------------------------------#
   if (lwr %>% upr){
      cat0(" - Lower limit: ",lwr)
      cat0(" - Upper limit: ",upr)
      stop(" Lower and upper limit must be finite and lower cannot be greater than upper.")
   }#end if (lwr %>% upr)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #   Don't let the bounds to be insane.                                                  #
   #---------------------------------------------------------------------------------------#
   if (! (buff %wr% c(0,1-2^-23))){
      cat0(" - buff: ",buff)
      stop(" Buffer must be between 0 (including) and 1 (excluding).")
   }#end if (lwr %>% upr)
   #---------------------------------------------------------------------------------------#


   xlwr = lwr + buff * (upr - lwr)
   xupr = upr - buff * (upr - lwr)
   ans  = pmin(xupr,pmax(xlwr,x)) + 0 * x
   return(ans)
}#end bound
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     This function finds the binary notation of an integer.                               #
#------------------------------------------------------------------------------------------#
tobin <<- function(x,reverse=TRUE){

   #----- X must be an integer, check it. -------------------------------------------------#
   if (! is.integer(x)){
      warning ("Function tobin: coercing x to an integer")
      x = as.integer(x)
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Check whether this is a single element.  If not, use recursive call. ------------#
   if (is.vector(x) && length(x) == 1){
      #----- Paste bits. ------------------------------------------------------------------#
      if (reverse){
        ans = paste(sapply(strsplit(paste(rev(intToBits(x))),""),`[[`,2),collapse="")
      }else{
        ans = paste(sapply(strsplit(paste(intToBits(x)),""),`[[`,2),collapse="")
      }#end if
      #------------------------------------------------------------------------------------#

   }else if (is.matrix(x) || is.array(x)){
      #----- Array or matrix. -------------------------------------------------------------#
      margin = length(dim(x))
      ans    = apply(X=x,MARGIN=sequence(margin),FUN=tobin,reverse=reverse)
      #------------------------------------------------------------------------------------#
   }else if (is.list(x)){
      #----- List or data frame. ----------------------------------------------------------#
      ans = sapply(X=x,FUN=tobin,reverse=reverse,simplify=is.data.frame(x))
      #------------------------------------------------------------------------------------#
   }else if (is.null(dim(x))){
      #----- Vector. ----------------------------------------------------------------------#
      ans = sapply(X=x,FUN=tobin,reverse=reverse,simplify=TRUE)
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#
   return(ans)
}#end function tobin
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     This function finds the binary notation of an integer.                               #
#------------------------------------------------------------------------------------------#
rawtoint <<- function(x){

   #----- Crash in case x is not raw. -----------------------------------------------------#
   dummy = stopifnot (is.raw(x))
   #---------------------------------------------------------------------------------------#



   #----- Pack all elements. --------------------------------------------------------------#
   ans   = as.integer(packBits(x))
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function rawtoint
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     This function computes the mean of the elements below (above) a given quantile.      #
#------------------------------------------------------------------------------------------#
qu.mean <<- function(x,p,na.rm=FALSE,lower=TRUE){
   if (na.rm) x = x[! is.na(x)]

   #----- Do the calculation only if there is anything left. ------------------------------#
   if (any(is.finite(x))){
      qu  = quantile(x,probs=p,na.rm=na.rm)
      if (lower){
         ans = mean(x[x <= qu],na.rm=na.rm)
      }else{
         ans = mean(x[x >= qu],na.rm=na.rm)
      }#end if
   }else{
      ans = NA
   }#end if

   return(ans)
}#end function qu.mean
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     This function is a quick integrator.  For more elegant ways to integrate a function, #
# check function quadrature.                                                               #
#------------------------------------------------------------------------------------------#
weighted.sum <<- function(x,w,na.rm=FALSE) sum(x*w,na.rm=na.rm)
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     These functions find weighted averages by applying the same weighting factor for     #
# each column (weighted.rowMeans) or rows (weighted.colMeans).                             #
#------------------------------------------------------------------------------------------#
weighted.rowMeans <<- function(X,wc,na.rm=FALSE){
   ans = apply( X      = X
              , MARGIN = 1
              , FUN    = weighted.mean
              , w      = wc
              , na.rm  = na.rm
              )#end apply
   return(ans)
}#end weighted.rowMeans
weighted.colMeans <<- function(X,wr,na.rm=FALSE){
   ans = apply( X      = X
              , MARGIN = 2
              , FUN    = weighted.mean
              , w      = wr
              , na.rm  = na.rm
              )#end apply
   return(ans)
}#end weighted.colMeans
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     This function finds the weighted fraction for each element of data.frame x, using    #
# weight w.  In case x is not a data frame, we will try to coerce it to a data frame,      #
# if it doesn't work then we crash it.  It will correct the final answer to provide        #
# weights that add up to one, and it will create a vector of equal chances in case weights #
# are all zeroes.                                                                          #
#------------------------------------------------------------------------------------------#
weighted.frac <<- function(x,w,na.rm=TRUE){
   #----- Make sure "x" is a data frame. --------------------------------------------------#
   if (! is.data.frame(x)){
      x = try(as.data.frame(x),silent=TRUE)
      #----- Give the bad news in case it doesn't coerce to a data frame. -----------------#
      if ("try-error" %in% is(x)){
         stop(" 'x' must be an object that can be coerced into a data frame!")
      }#end if ("try-error" %in% is(x))
      #------------------------------------------------------------------------------------#
   }#end if (! is.data.frame(x))
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Stop in case the dimensions of x and w don't match.                               #
   #---------------------------------------------------------------------------------------#
   w = c(unlist(w))
   if (length(w) != nrow(x)){
      stop(" 'x' and 'w' must have compatible dimensions (length(x) = nrow(x))!")
   }#end if (length(w) != nrow(x))
   #---------------------------------------------------------------------------------------#



   #------ In case all weights are zero, make them equal. ---------------------------------#
   if (na.rm){
      keep = rowSums(! is.finite(as.matrix(x))) == 0 & is.finite(w)
      x    = x[keep,,drop=FALSE]
      w    = w[keep]
   }else if (any(! is.finite(w))){
      #----- Return NA in case w has non-finite elements and na.rm = FALSE. ---------------#
      ans        = rep(NA,times=length(x))
      names(ans) = names(x)
      return(ans)
      #------------------------------------------------------------------------------------#
   }else if (all(w %==% 0.)){
      #----- Give equal chances in case all weights were zero. ----------------------------#
      w = rep(x=1.,times=nrow(x))
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Find the weighted mean of each element of x.                                     #
   #---------------------------------------------------------------------------------------#
   ans = sapply(X=x,FUN=weighted.mean,w=w)
   ans = ans / sum(ans)
   names(ans) = names(x)
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end weighted.frac
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     This function estimates the quantile for a table of observations x, each of which    #
# having a weight w.  This is done by finding the median of a pseudo dataset built using   #
# sample.  If the size of the resampling is not provided, then the number of samples is    #
# dependent on the range of probabilities.  By default we find the 0.50 quantile (median). #
#------------------------------------------------------------------------------------------#
weighted.quantile <<- function(x,w,qu=0.50,size.minp=10,na.rm=FALSE,out.case=FALSE){

   #----- Delete the missing values if the user asked to do it. ---------------------------#
   if (any(w <= 0, na.rm = TRUE) || any(is.infinite(w)) || any(is.na(w))){
      stop(" Weights (w) must be positive and finite, and entirely defined!")
   }else if(qu < 0. || qu > 1.){
      stop(" Quantile (qu) must be between 0 and 1. ")
   }else if(na.rm){
      keep = ! is.na(x)
      x    = x[keep]
      w    = w[keep]
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Define the probabilities by normalising the weights.                             #
   #---------------------------------------------------------------------------------------#
   prob = w / sum(w)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Sort the values by the probability.                                              #
   #---------------------------------------------------------------------------------------#
   o     = order(x,decreasing=FALSE)
   x     = x[o]
   w     = w[o]
   prob  = prob[o]
   cum   = cumsum(prob)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Sort the values by the probability.                                              #
   #---------------------------------------------------------------------------------------#
   if (qu <= cum[1]){
      qout = x[1]
      case = "minimum"
   }else if (qu >= cum[length(cum)]){
      qout = x[length(cum)]
      case = "maximum"
   }else if (any(cum == qu)){
      qout = x[which(cum == qu)]
      case = "exact"
   }else{
      below   = qu - cum ; below[below < 0] = Inf 
      above   = cum - qu ; above[above < 0] = Inf
      i.below = which.min(below)
      i.above = which.min(above)
      w.below = 1. / (below[i.below]^2)
      w.above = 1. / (above[i.above]^2)
      qout    = ( x[i.below] * w.below + x[i.above] * w.above ) / (w.below + w.above)
      case    = "interpolated"
   }#end if
   #---------------------------------------------------------------------------------------#


   #---- Decide what to return. -----------------------------------------------------------#
   if (out.case){
      ans = list(q = qout, case = case)
   }else{
      ans = qout
   }#end if (out.case)
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function weighted.quantile
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function estimates the weighted standard deviation.                             #
#------------------------------------------------------------------------------------------#
weighted.sd <<- function(x,w,M=NULL,na.rm=FALSE){

   #----- Delete the missing values if the user asked to do it. ---------------------------#
   if (any(w < 0, na.rm = TRUE) || any(is.infinite(w)) || any(is.na(w))){
      stop(" Weights (w) must be non-negative and finite, and entirely defined!")
   }else if(na.rm){
      keep = ! is.na(x)
      x    = x[keep]
      w    = w[keep]
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Assume M to be inversely proportional to the smallest weight. -------------------#
   if (is.null(M)) M = min(w)
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #      Check whether at least one weight is non-zero.                                   #
   #---------------------------------------------------------------------------------------#
   if (all(w %==% 0)){
      ans = NA
   }else{
      xwm    = weighted.mean(x=x,w=w)
      w.sum  = sum(w)
      r2.sum = sum(w*(x-xwm)^2)
      ans    = sqrt(M * r2.sum / (M * w.sum - 1))
   }#end if
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end function weighted.sd
#==========================================================================================#
#==========================================================================================#




#==========================================================================================#
#==========================================================================================#
#     This function computes the standard error.                                           #
#------------------------------------------------------------------------------------------#
se <<- function (x, na.rm = FALSE){
   #---- We follow the same convention as the standard deviation. -------------------------#
   if (is.matrix(x)) {
      msg = "se(<matrix>) is deprecated.\n Use apply(*, 2, se) instead."
      warning(paste(msg, collapse = ""), call. = FALSE, domain = NA)
      ans = apply(X = x, MARGIN = 2, FUN = se, na.rm = na.rm)
   }else if (is.list(x)){
      msg = "se(<list>) is deprecated.\n Use lapply(*, se) instead."
      warning(paste(msg, collapse = ""), call. = FALSE, domain = NA)
      ans = lapply(X = x, FUN = se, na.rm = na.rm)
   }else if (is.data.frame(x)){
      msg = "se(<data.frame>) is deprecated.\n Use sapply(*, se) instead."
      warning(paste(msg, collapse = ""), call. = FALSE, domain = NA)
      ans = sapply(X = x, se, na.rm = na.rm)
   }else{
      #----- Coerce x to a vector. --------------------------------------------------------#
      if (is.vector(x)){
         xx = x
      }else{
         xx = as.vector(x)
      }#end if
      #------------------------------------------------------------------------------------#



      #----- Decide whether to delete NA or not. ------------------------------------------#
      if (na.rm) xx = xx[! is.na(xx)]
      #------------------------------------------------------------------------------------#



      #----- Find the standard error. -----------------------------------------------------#
      nx  = length(xx)
      ans = sd(xx)/sqrt(nx)
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end function se
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     This function computes the skewness, or the third moment of a distribution.          #
#------------------------------------------------------------------------------------------#
skew <<- function (x, na.rm = FALSE){

   #---- We follow the same convention as the standard deviation. -------------------------#
   if (is.matrix(x)) {
      msg = "skew(<matrix>) is deprecated.\n Use apply(*, 2, skew) instead."
      warning(paste(msg, collapse = ""), call. = FALSE, domain = NA)
      ans = apply(X = x, MARGIN = 2, FUN = skew, na.rm = na.rm)
   }else if (is.list(x)){
      msg = "skew(<list>) is deprecated.\n Use lapply(*, skew) instead."
      warning(paste(msg, collapse = ""), call. = FALSE, domain = NA)
      ans = lapply(X = x, FUN = skew, na.rm = na.rm)
   }else if (is.data.frame(x)){
      msg = "skew(<data.frame>) is deprecated.\n Use sapply(*, skew) instead."
      warning(paste(msg, collapse = ""), call. = FALSE, domain = NA)
      ans = sapply(X = x, skew, na.rm = na.rm)
   }else{
      #----- Coerce x to a vector. --------------------------------------------------------#
      if (is.vector(x)){
         xx = x
      }else{
         xx = as.vector(x)
      }#end if
      #------------------------------------------------------------------------------------#



      #----- Decide whether to delete NA or not. ------------------------------------------#
      if (na.rm) xx = xx[! is.na(xx)]
      #------------------------------------------------------------------------------------#



      #----- Find the skewness. -----------------------------------------------------------#
      nx = length(xx)
      xx.mean = mean(xx)
      xx.sdev = sd(xx)
      ans     = sum((xx-xx.mean)^3) / (nx * xx.sdev^3)
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end function skew
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     This function computes the kurtosis, or the fourth moment of a distribution.         #
#     NOTE: This calculates the absolute kurtosis, to get the excess kurtosis you must     #
#           subtract 3.                                                                    #
#------------------------------------------------------------------------------------------#
kurt <<- function (x, na.rm = FALSE){

   #---- We follow the same convention as the standard deviation. -------------------------#
   if (is.matrix(x)) {
      msg = "kurt(<matrix>) is deprecated.\n Use apply(*, 2, kurt) instead."
      warning(paste(msg, collapse = ""), call. = FALSE, domain = NA)
      ans = apply(X = x, MARGIN = 2, FUN = kurt, na.rm = na.rm)
   }else if (is.list(x)){
      msg = "kurt(<list>) is deprecated.\n Use lapply(*, kurt) instead."
      warning(paste(msg, collapse = ""), call. = FALSE, domain = NA)
      ans = lapply(X = x, FUN = kurt, na.rm = na.rm)
   }else if (is.data.frame(x)){
      msg = "kurt(<data.frame>) is deprecated.\n Use sapply(*, kurt) instead."
      warning(paste(msg, collapse = ""), call. = FALSE, domain = NA)
      ans = sapply(X = x, kurt, na.rm = na.rm)
   }else{
      #----- Coerce x to a vector. --------------------------------------------------------#
      if (is.vector(x)){
         xx = x
      }else{
         xx = as.vector(x)
      }#end if
      #------------------------------------------------------------------------------------#



      #----- Decide whether to delete NA or not. ------------------------------------------#
      if (na.rm) xx = xx[! is.na(xx)]
      #------------------------------------------------------------------------------------#



      #----- Find the kurtosis. -----------------------------------------------------------#
      nx = length(xx)
      xx.mean = mean(xx)
      xx.sdev = sd(xx)
      ans     = sum((xx-xx.mean)^4) / (nx * xx.sdev^4)
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end function kurt
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     This function computes the inter-quartile range.                                     #
#------------------------------------------------------------------------------------------#
iqr <<- function (x, na.rm = FALSE){

   #---- We follow the same convention as the standard deviation. -------------------------#
   if (is.matrix(x)) {
      msg = "iqr(<matrix>) is deprecated.\n Use apply(*, 2, kurt) instead."
      warning(paste(msg, collapse = ""), call. = FALSE, domain = NA)
      ans = apply(X = x, MARGIN = 2, FUN = iqr, na.rm = na.rm)
   }else if (is.list(x)){
      msg = "iqr(<list>) is deprecated.\n Use lapply(*, kurt) instead."
      warning(paste(msg, collapse = ""), call. = FALSE, domain = NA)
      ans = lapply(X = x, FUN = iqr, na.rm = na.rm)
   }else if (is.data.frame(x)){
      msg = "iqr(<data.frame>) is deprecated.\n Use sapply(*, kurt) instead."
      warning(paste(msg, collapse = ""), call. = FALSE, domain = NA)
      ans = sapply(X = x, iqr, na.rm = na.rm)
   }else{
      #----- Coerce x to a vector. --------------------------------------------------------#
      if (is.vector(x)){
         xx = x
      }else{
         xx = as.vector(x)
      }#end if
      #------------------------------------------------------------------------------------#



      #----- Decide whether to delete NA or not. ------------------------------------------#
      if (na.rm) xx = xx[! is.na(xx)]
      #------------------------------------------------------------------------------------#



      #----- Find the inter-quartile range. -----------------------------------------------#
      ans        = diff(quantile(xx,probs=c(0.25,0.75)))
      names(ans) = NULL
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end function iqr
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     This function computes the first four moments of the distribution.                   #
#------------------------------------------------------------------------------------------#
four.moments <<- function (x, na.rm = FALSE){

   #---- We follow the same convention as the standard deviation. -------------------------#
   if (is.matrix(x)) {
      msg = "four.moments(<matrix>) is deprecated.\n Use apply(*, 2, four.moments) instead."
      warning(paste(msg, collapse = ""), call. = FALSE, domain = NA)
      ans = apply(X = x, MARGIN = 2, FUN = four.moments, na.rm = na.rm)
   }else if (is.list(x)){
      msg = "four.moments(<list>) is deprecated.\n Use lapply(*, four.moments) instead."
      warning(paste(msg, collapse = ""), call. = FALSE, domain = NA)
      ans = lapply(X = x, FUN = four.moments, na.rm = na.rm)
   }else if (is.data.frame(x)){
      msg = paste0("four.moments(<data.frame>) is deprecated."
                  ,"\n"
                  ,"Use sapply(*, four.moments) instead."
                  )#end paste0
      warning(paste(msg, collapse = ""), call. = FALSE, domain = NA)
      ans = sapply(X = x, four.moments, na.rm = na.rm)
   }else{
      #----- Coerce x to a vector. --------------------------------------------------------#
      if (is.vector(x)){
         xx = x
      }else{
         xx = as.vector(x)
      }#end if
      #------------------------------------------------------------------------------------#



      #----- Decide whether to delete NA or not. ------------------------------------------#
      if (na.rm) xx = xx[! is.na(xx)]
      #------------------------------------------------------------------------------------#



      #----- Find the four moments. -------------------------------------------------------#
      nx  = length(xx)
      ans = c( mean     = mean(xx)
             , variance = var (xx)
             , skewness = skew(xx)
             , kurtosis = kurt(xx)
             )#end c
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end function four.moments
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     This function computes the mean, standard deviation, and coefficient of variation    #
# of a vector.                                                                             #
#------------------------------------------------------------------------------------------#
meansdcv <<- function (x, na.rm = FALSE){

   #---- We follow the same convention as the standard deviation. -------------------------#
   if (is.matrix(x)) {
      msg = "meansdcv(<matrix>) is deprecated.\n Use apply(*, 2, meansdcv) instead."
      warning(paste(msg, collapse = ""), call. = FALSE, domain = NA)
      ans = apply(X = x, MARGIN = 2, FUN = meansdcv, na.rm = na.rm)
   }else if (is.list(x)){
      msg = "meansdcv(<list>) is deprecated.\n Use lapply(*, meansdcv) instead."
      warning(paste(msg, collapse = ""), call. = FALSE, domain = NA)
      ans = lapply(X = x, FUN = meansdcv, na.rm = na.rm)
   }else if (is.data.frame(x)){
      msg = "meansdcv(<data.frame>) is deprecated.\n Use sapply(*, meansdcv) instead."
      warning(paste(msg, collapse = ""), call. = FALSE, domain = NA)
      ans = sapply(X = x, meansdcv, na.rm = na.rm)
   }else{
      #----- Coerce x to a vector. --------------------------------------------------------#
      if (is.vector(x)){
         xx = x
      }else{
         xx = as.vector(x)
      }#end if
      #------------------------------------------------------------------------------------#



      #----- Decide whether to delete NA or not. ------------------------------------------#
      if (na.rm) xx = xx[! is.na(xx)]
      #------------------------------------------------------------------------------------#



      #----- Find the mean, standard deviation, and coefficient of variation. -------------#
      nx    = length(xx)
      mu    = mean(xx)
      sigma = sd(xx)
      cvar  = ifelse(sigma %>% 0,mu/sigma,NA)
      
      ans = c( mean = mu
             , sd   = sigma
             , cv   = cvar
             )#end c
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end function meansdcv
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#      This function finds the element index given the array index and the dimensions of   #
# the array.  This is the inverse of function arrayInd in R.                               #
#                                                                                          #
# Original author: Feng Li, Department of Statistics, Stockholm University, Sweden.        #
#                  Based on last version as of Wed Mar 14 18:03:09 CET 2012.               #
#                                                                                          #
# Modified by:     Marcos Longo.  Department of Earth and Planetary Sciences,              #
#                                 Harvard University, Cambridge, MA, USA                   #
#                  Last modified on 25 Oct 2012 - 10:56 EST                                #
#                                                                                          #
#                  The script now recognises whether the arr.ind is a vector, matrix, or   #
#                  list and call it self recursively to return the full list.              #
#------------------------------------------------------------------------------------------#
whichInd <<- function(ai, dims){

   #----- Save the number of dimensions. --------------------------------------------------#
   n.dims  = length(dims)
   #---------------------------------------------------------------------------------------#



   #---- Check the variable type. ---------------------------------------------------------#
   if (is.matrix(ai) || is.data.frame(ai)) {
      #----- Check that the dimensions match. ---------------------------------------------#
      if ( ncol(ai) != n.dims){
         cat(" - Number of columns of ai: ",ncol(ai),"\n")
         cat(" - Length of dims:          ",n.dims  ,"\n")
         stop(" Dimensions between ai and dims don't match!")
      }#end if
      #------------------------------------------------------------------------------------#

      ans = apply (X = ai, FUN = whichInd, MARGIN = 1, dims=dims)

   }else if (is.list(ai)){
      #----- Check that the dimensions match. ---------------------------------------------#
      if ( any(sapply(X=ai,FUN=length) != n.dims)){
         fail = sum(sapply(X=ai,FUN=length) != n.dims)
         cat (" - ",fail," elements of ai don't have correct dimensions","\n")
         cat (" - Length of dims: ",n.dims,"\n")
         stop(" Dimensions between some ai elements and dims don't match!")
      }#end if
      #------------------------------------------------------------------------------------#

      ans = lapply(X = ai, FUN = whichInd, dims=dims)
   }else{

      #----- Coerce the data to be a vector. ----------------------------------------------#
      if (length(ai) != length(dims)){
         cat (" - Length of ai  : ",length(ai),"\n")
         cat (" - Length of dims: ",n.dims    ,"\n")
         stop(" ai must have the same length as dims")
      }#end if
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #  Variable cumdim is the number of elements we jump every time the index increases  #
      # by 1.  Variable shif is just the array     #
      # 
      #------------------------------------------------------------------------------------#
      cumdim = c(1,cumprod(dims[-n.dims]))
      ans    = 1 + sum(cumdim*(ai-1))
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#

   return(ans)   
}#end whichInd
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function is a safe wrapper for function sample, to ensure that we always sample #
# the object given, even when it it an object of size 1.                                   #
#------------------------------------------------------------------------------------------#
lit.sample <<- function(x,size,replace=FALSE,prob=NULL){
   if (length(x) == 1){
      ans = rep(x=x,times=size)
   }else{
      ans = sample(x,size,replace,prob)
   }#end if
   return(ans)
}#end function lit.sample
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function creates a sequence that covers the entire range of data: length.out is #
# the sought length (n is a shorter option).                                               #
#------------------------------------------------------------------------------------------#
seq.range <<- function(x,n=length.out,length.out){
   xrange = range(x,finite=TRUE)
   if (any(is.finite(xrange))){
      ans = seq(from=xrange[1],to=xrange[2],length.out=n)
   }else{
      ans = rep(NA,times=n)
   }#end if
   return(ans)
}#end function seq.range
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function counts how many unique values exist in a vector.                       #
#------------------------------------------------------------------------------------------#
length.unique <<- function(x){
   ans = length(unique(x))
   return(ans)
}#end function lit.sample
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function finds the ratio of consecutive elements.  All arguments used by diff   #
# can be used here too.                                                                    #
#------------------------------------------------------------------------------------------#
ediff <<- function(x,lag=1,differences=1){
   if (differences == 1){
      n   = length(x)
      ans = x[seq(from=1+lag,to=n,by=lag)]/x[seq(from=1,to=n-lag,by=lag)]
      ans[! is.finite(ans)] = NA
   }else{
      ans = ediff(x,lag=lag,differences=differences-1)
   }#end if
   return(ans)
}#end ediff
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function finds the weighted.mean for each col, similar to colMeans.             #
#------------------------------------------------------------------------------------------#
colWgtMeans <<- function(x,w,...,na.rm=FALSE){
   #---------------------------------------------------------------------------------------#
   #     First sanity check.                                                               #
   #---------------------------------------------------------------------------------------#
   mess  = c("x must be a matrix, data.frame, or an array of dimension 2!"
            ,"w must be a matrix, data.frame, or an array of dimension 2!"
            ,"x and w must have identical sizes!"
            )#end mess
   error = c(! (is.matrix(x) || is.data.frame(x) || (is.array(x) && length(dim(x)) == 2))
            ,! (is.matrix(w) || is.data.frame(w) || (is.array(w) && length(dim(w)) == 2))
            ,any(dim(x) != dim(w))
            )#end c
   if (any(error)) stop(paste(mess,collapse="\n"))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Convert x and w into matrices                                                     #
   #---------------------------------------------------------------------------------------#
   x            = as.matrix(x)
   w            = as.matrix(w)
   x[is.nan(x)] = NA
   w[is.nan(w)] = NA
   w[is.na (w)] = 0
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Second sanity check.                                                              #
   #---------------------------------------------------------------------------------------#
   mess   = c("NA/NaN are not allowed in w (only x may contain NA)!"
             ,"All w must be finite and positive"
             ,"Some w columns have only zeroes"
             )#end c
   error  = c( any(is.na(c(w)))
             , any(! is.finite(c(w)) | (is.finite(c(w)) & c(w) < 0))
             , any(apply(X=w,MARGIN=2,FUN=sum) == 0,na.rm=TRUE)
             )#end error
   if (any(error)) stop(paste(mess,collapse="\n"))
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Find weighted means.                                                              #
   #---------------------------------------------------------------------------------------#
   ans = apply(X=x*w,MARGIN=2,FUN=sum,na.rm=na.rm) / apply(X=w,MARGIN=2,FUN=sum)
   #---------------------------------------------------------------------------------------#
   return(ans)
}#end colWgtMeans
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function finds the weighted.mean for each row, similar to rowMeans.             #
#------------------------------------------------------------------------------------------#
rowWgtMeans <<- function(x,w,...,na.rm=FALSE){
   #---------------------------------------------------------------------------------------#
   #     First sanity check.                                                               #
   #---------------------------------------------------------------------------------------#
   mess  = c("x must be a matrix, data.frame, or an array of dimension 2!"
            ,"w must be a matrix, data.frame, or an array of dimension 2!"
            ,"x and w must have identical sizes!"
            )#end mess
   error = c(! (is.matrix(x) || is.data.frame(x) || (is.array(x) && length(dim(x)) == 2))
            ,! (is.matrix(w) || is.data.frame(w) || (is.array(w) && length(dim(w)) == 2))
            ,any(dim(x) != dim(w))
            )#end c
   if (any(error)) stop(paste(mess,collapse="\n"))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Convert x and w into matrices                                                     #
   #---------------------------------------------------------------------------------------#
   x            = as.matrix(x)
   w            = as.matrix(w)
   x[is.nan(x)] = NA
   w[is.nan(w)] = NA
   w[is.na (w)] = 0
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Second sanity check.                                                              #
   #---------------------------------------------------------------------------------------#
   mess   = c("NA/NaN are not allowed in w (only x may contain NA)!"
             ,"All w must be finite and positive"
             ,"Some w columns have only zeroes"
             )#end c
   error  = c( any(is.na(c(w)))
             , any(! is.finite(c(w)) | (is.finite(c(w)) & c(w) < 0))
             , any(apply(X=w,MARGIN=1,FUN=sum) == 0,na.rm=TRUE)
             )#end error
   if (any(error)) stop(paste(mess,collapse="\n"))
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Find weighted means.                                                              #
   #---------------------------------------------------------------------------------------#
   ans = apply(X=x*w,MARGIN=1,FUN=sum,na.rm=na.rm) / apply(X=w,MARGIN=1,FUN=sum)
   #---------------------------------------------------------------------------------------#
   return(ans)
}#end colWgtMeans
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function normalises a vector.  Results are just the normalised values, mean and #
# standard deviation are in attributes.                                                    #
#------------------------------------------------------------------------------------------#
normalise <<- function(x,mu,sigma){
   if (missing(mu   )) mu    = mean(x,na.rm = TRUE)
   if (missing(sigma)) sigma = sd  (x,na.rm = TRUE)
   nx      = sum(is.finite(x))
   normal  = (x - mu) / sigma
   normal  = ifelse(is.finite(normal),normal,NA)
   attributes(normal) = list(mean = mu, sdev = sigma, n = nx)
   return(normal)
}#end normalise
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function converts data into percentil (0-100).                                  #
#------------------------------------------------------------------------------------------#
percentil <<- function(x,trim=0.0){
   qlow    = 0.0+0.5*trim
   qhigh   = 1.0-0.5*trim

   xlow    = quantile(x,probs=qlow ,na.rm=TRUE)
   xhigh   = quantile(x,probs=qhigh,na.rm=TRUE)
   xperc   = 100. * pmax(0.,pmin(1.,(x-xlow)/(xhigh-xlow)))
   xperc   = ifelse(is.finite(xperc),xperc,NA)
   return(xperc)
}#end normalise
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     A lazy option of square root of sum of squares and root mean sum of squares.         #
#------------------------------------------------------------------------------------------#
sum2           <<- function(x,...)   sqrt(x = sum          (x=x^2               ,...))
mean2          <<- function(x,...)   sqrt(x = mean         (x=x^2               ,...))
meanlog        <<- function(x,...)   exp (x = mean         (x=log(x)            ,...))
weighted.mean2 <<- function(x,w,...) sqrt(x = weighted.mean(x=x^2,w=(w/sum(w))^2,...))
mean.se        <<- function(x,...)   sqrt(x = mean(x=x^2,...) / length(x[is.finite(x)]))
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     Mean of elements that are finite and above minimum.                                  #
#------------------------------------------------------------------------------------------#
mean.above <<- function(x,xlwr,xnot=xlwr) if(any(x%>=%xlwr)){mean(x[x%>=%xlwr])}else{xnot}
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     Safe versions of common functions, which only uses finite numbers, and return either #
# finite answers of NA.                                                                    #
#------------------------------------------------------------------------------------------#
#------ Mean. -----------------------------------------------------------------------------#
mean.safe <<- function(x){
   ans = mean(x,na.rm=TRUE)
   if (! is.finite(ans)) ans = NA
   return(ans)
}#end if
#------ Median. ---------------------------------------------------------------------------#
median.safe <<- function(x){
   ans = median(x,na.rm=TRUE)
   if (! is.finite(ans)) ans = NA
   return(ans)
}#end if
#------ Standard deviation. ---------------------------------------------------------------#
sd.safe <<- function(x){
   ans = sd(x,na.rm=TRUE)
   if (! is.finite(ans)) ans = NA
   return(ans)
}#end if
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     Function to return the mid points between two consecutive points.                    #
#------------------------------------------------------------------------------------------#
mid.points <<- function(x,islog=FALSE,finite=FALSE){
   cx   = if (islog){log(c(x))}else{c(x)}
   cx   = ifelse(is.finite(cx),cx,NA)
   ix   = seq_along(cx)
   last = length(cx)
   if (finite){
      cx    = na.approx(cx,na.rm=FALSE)
      fa    = min(which(is.finite(cx)))
      fz    = max(which(is.finite(cx)))
      delta = max(c(0,mean(diff(cx),na.rm=TRUE)),na.rm=TRUE)
      #----- Append data to the left. -----------------------------------------------------#
      if (fa > 1){
         fill     = seq(from=1,to=fa-1,by=1)
         cx[fill] = cx[fa] + delta * (ix[fill]-ix[fa])
      }#end if
      if (fz < last){
         fill     = seq(from=fz+1,to=last,by=1)
         cx[fill] = cx[fz] + delta * (ix[fill]-ix[fz])
      }#end if
   }#end if

   xmid = rowMeans(cbind(cx[-1],cx[-last]))
   if(islog) xmid = exp(xmid)
   return(xmid)
}#end mid.points
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     Functions to return the fraction of points that are xxx (xxx = finite, NA, NaN).     #
#------------------------------------------------------------------------------------------#
frac.finite <<- function(x) sum(is.finite(x))/max(1,length(x))
frac.na     <<- function(x) sum(is.na    (x))/max(1,length(x))
frac.nan    <<- function(x) sum(is.nan   (x))/max(1,length(x))
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function evaluates the function only if there are enough finite points, other-  #
# wise it returns NA.                                                                      #
#------------------------------------------------------------------------------------------#
ifenough <<- function(x,f,ef.min = 1/3, ...){
   if (frac.finite(x) >= ef.min){
      ans = f(x,...)
   }else{
      ans = NA + f(runif(n=length(x),min=1,max=2),...)
   }#end if
   return(ans)
}#end ifenough
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#        For each element of vector x, this function finds the value in the "onto" vector  #
# that corresponds to the same quantile of the distribution defined by x.                  #
#------------------------------------------------------------------------------------------#
qqproject <<- function(x,onto){
   #----- Find out which elements are finite and save the position. -----------------------#
   ok      = which(! is.na(x))
   project = NA * x
   #---------------------------------------------------------------------------------------#

   #----- Keep only the non-NA elements. --------------------------------------------------#
   xf      = x[!is.na(x)]
   #---------------------------------------------------------------------------------------#


   #----- Find the corresponding cdf to which of the elements. ----------------------------#
   of      = order(xf)
   qxf     = seq(from=0,to=1,length.out=length(xf))
   qxf[of] = qxf
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   project[ok] = quantile(x=onto,probs=qxf,na.rm=TRUE)
   return(project)
   #---------------------------------------------------------------------------------------#
}#end
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function evaluates the confidence interval of a mean without assuming that the  #
# sample size is large enough.                                                             #
#------------------------------------------------------------------------------------------#
ci.mean <<- function(x,conf=0.95,finite=FALSE,...){
   if (finite) x = x[is.finite(x)]
   nx = length(x)
   if (nx <= 1){
      ci = c(NA,NA)
   }else{
      ci = mean(x) + qt(p=0.5*(1.+c(-1.,1.)*conf),df=nx-1) * se(x)
   }#end if
   return(ci)
}#end ifenough
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     These variables may be used everywhere, they define the order of the six-summary     #
# vector.                                                                                  #
#------------------------------------------------------------------------------------------#
six.summary.names <<- c("expected","variance","skewness","kurtosis","ci.lower","ci.upper")
n.six.summary     <<- length(six.summary.names)
#==========================================================================================#
#==========================================================================================#








#==========================================================================================#
#==========================================================================================#
#     This function calculates the four moments of the distribution plus the 95% C.I. of   #
# the mean using t distributioon.                                                          #
#------------------------------------------------------------------------------------------#
six.summary <<- function(x,conf=0.95,finite=TRUE,...){
   #------ Remove non-finite data in case finite=TRUE. ------------------------------------#
   if (finite) x = x[is.finite(x)]
   nx = length(x)
   #---------------------------------------------------------------------------------------#



   #------ Initialise data with NA in case the function fails. ----------------------------#
   ans = rep(NA,times=n.six.summary)
   names(ans) = six.summary.names
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Check that there are sufficient points for at least the mean.                    #
   #---------------------------------------------------------------------------------------#
   if (nx >= 1){
      #----- Find the moments and confidence interval. ------------------------------------#
      expected = mean(x)
      variance = var (x)
      skewness = skew(x)
      kurtosis = kurt(x)
      ci.lower = expected + qt(p=0.5*(1.0-conf),df=nx-1) * se(x)
      ci.upper = expected + qt(p=0.5*(1.0+conf),df=nx-1) * se(x)
      #------------------------------------------------------------------------------------#




      #----- Standardise non-finite values to NA. -----------------------------------------#
      expected = ifelse(is.finite(expected),expected,NA)
      variance = ifelse(is.finite(variance),variance,NA)
      skewness = ifelse(is.finite(skewness),skewness,NA)
      kurtosis = ifelse(is.finite(kurtosis),kurtosis,NA)
      ci.lower = ifelse(is.finite(ci.lower),ci.lower,NA)
      ci.upper = ifelse(is.finite(ci.upper),ci.upper,NA)
      #------------------------------------------------------------------------------------#



      #----- Make sure statistics go to the right place. ----------------------------------#
      ans[match("expected",six.summary.names)] = expected
      ans[match("variance",six.summary.names)] = variance
      ans[match("skewness",six.summary.names)] = skewness
      ans[match("kurtosis",six.summary.names)] = kurtosis
      ans[match("ci.lower",six.summary.names)] = ci.lower
      ans[match("ci.upper",six.summary.names)] = ci.upper
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end function six.summary
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      Thirteen-number summary:                                                            #
#                                                                                          #
# qmin - minimum                                                                           #
# q025 - 2.5 percentile                                                                    #
# q250 - first quartile                                                                    #
# q500 - median                                                                            #
# q750 - third quartile                                                                    #
# q975 - 97.5 percentile                                                                   #
# qmax - maximum                                                                           #
# mean - mean                                                                              #
# sdev - standard deviation                                                                #
# skew - skewness                                                                          #
# kurt - kurtosis                                                                          #
# ntot - total number                                                                      #
# nval - total number of valide (i.e. finite) entries)                                     #
#------------------------------------------------------------------------------------------#
thirteen.num <<- function(x){
   #----- Count total number and number of valid entries. ---------------------------------#
   ntot = length(x)
   fine = is.finite(x)
   nval = sum(fine)
   #---------------------------------------------------------------------------------------#


   #----- Discard invalid numbers. --------------------------------------------------------#
   x    = x[fine]
   #---------------------------------------------------------------------------------------#


   #----- Find quantiles. -----------------------------------------------------------------#
   quant = quantile(x=x,probs=c(0,0.025,0.25,0.50,0.75,0.975,1.000))
   names(quant) = c("qmin","q025","q250","q500","q750","q975","q100")
   #---------------------------------------------------------------------------------------#



   #----- Find mean, sd, skewness and kurtosis. -------------------------------------------#
   four  = c(mean=mean(x),sdev=sd(x),skew=skew(x),kurt=kurt(x))
   #---------------------------------------------------------------------------------------#


   #----- Append all results, and standardise weird values to NA. -------------------------#
   ans   = c(quant,four,ntot=ntot,nval=nval)
   ans[! is.finite(ans)] = NA
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end thirteen.num
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      Thirteen-number error assessment:                                                   #
#                                                                                          #
#  Input:                                                                                  #
#    x       -- Input variable                                                             #
#    e       -- Error (same units as x)                                                    #
#    min.ok  -- Minimum acceptable value                                                   #
#    max.ok  -- Maximum acceptable value                                                   #
#    nr      -- Number of replicates                                                       #
#    szmax   -- Maximum memory size to be allocated (so large vectors don't exhaust the    #
#               computer).                                                                 #
#    fun     -- Function to apply to the data set.                                         #
#    verbose -- Should the run be verbose?                                                 #
#                                                                                          #
#  Output:                                                                                 #
#    emean  -- expected mean (no replications)                                             #
#    rmean  -- average of all replications.                                                #
#    rsdev  -- standard deviation of the replicates (standard error of the mean)           #
#    rqmin  -- minimum average from replications                                           #
#    rq025  -- 2.5 percentile                                                              #
#    rq250  -- first quartile                                                              #
#    rq500  -- median                                                                      #
#    rq750  -- third quartile                                                              #
#    rq975  -- 97.5 percentile                                                             #
#    rqmax  -- maximum                                                                     #
#    ntot   -- total number                                                                #
#    nval   -- total number of valid (i.e. finite) entries)                                #
#    nrep   -- number of replicates (same as input)                                        #
#------------------------------------------------------------------------------------------#
thirteen.err <<- function( x
                         , e
                         , min.ok  = -Inf
                         , max.ok  = +Inf
                         , nr      = 10000
                         , szmax   = 100000
                         , fun     = "mean"
                         , verbose = FALSE
                         , ...
                         ){
   #----- Find the function. --------------------------------------------------------------#
   fun = match.fun(fun)
   #---------------------------------------------------------------------------------------#

   #----- Count total number and number of valid entries. ---------------------------------#
   ntot  = length(x)
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #    Select only those entries with actual error evaluation.                            #
   #---------------------------------------------------------------------------------------#
   fine  = is.finite(x) & is.finite(e)
   nval  = sum(fine)
   xfine = x[fine]
   efine = e[fine]
   #---------------------------------------------------------------------------------------#

   
   #---------------------------------------------------------------------------------------#
   #     Create vector with outcomes of each realisation.                                  #
   #---------------------------------------------------------------------------------------#
   xr = rep(NA,times=nr)
   #---------------------------------------------------------------------------------------#

   
   #---------------------------------------------------------------------------------------#
   #     Find the maximum block size.                                                      #
   #---------------------------------------------------------------------------------------#
   block.size = max(1,floor(szmax/nval))
   ia.full    = seq(from=1,to=nr,by=block.size)
   iz.full    = c(ia.full[-1]-1,nr)
   nblocks    = length(ia.full)
   if (verbose){
      cat("    -> 13.err.  Vector size: ",nval,".","\n",sep="")
      cat("    -> 13.err.  Total number of blocks: ",nblocks,".","\n",sep="")
   }#end if (verbose)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Populate outcome vector by blocks, so we run it efficiently but don't run out of  #
   # memory.                                                                               #
   #---------------------------------------------------------------------------------------#
   ishow = round(sequence(10)*nblocks/10)
   lshow = 10*sequence(10)
   for (b in sequence(nblocks)){
      if (verbose && b %in% ishow){
         lshow.now = lshow[match(b,ishow)]
         cat("       .. 13.err.  ",lshow[match(b,ishow)],"% completed...","\n",sep="")
      }#end if (verbose && b %in% ishow)
      #----- Prepare selection. -----------------------------------------------------------#
      idx  = seq(from=ia.full[b],to=iz.full[b],by=1)
      nidx = length(idx)
      #------------------------------------------------------------------------------------#

      #----- Discard bad entries. ---------------------------------------------------------#
      X = matrix(data=rep(x=xfine,times=nidx),nrow=nval,ncol=nidx)
      E = matrix(data=rep(x=efine,times=nidx),nrow=nval,ncol=nidx)
      #------------------------------------------------------------------------------------#


      #----- Create vector with random error. ---------------------------------------------#
      XR      = 0.*X + pmax(min.ok,pmin(max.ok,rnorm(n=nval*nidx,mean=X,sd=E)))
      xr[idx] = apply(X=XR,MARGIN=2,FUN=fun,...)
      rm(X,E,XR)
      #------------------------------------------------------------------------------------#
   }#end for (b in sequence(nblocks))
   #---------------------------------------------------------------------------------------#



   #----- Find mean, sd, skewness and kurtosis. -------------------------------------------#
   three = c(expct = fun(xfine,...), rmean = mean(xr), rsdev = sd  (xr))
   #---------------------------------------------------------------------------------------#


   #----- Find quantiles. -----------------------------------------------------------------#
   quant = quantile(x=xr,probs=c(0,0.025,0.25,0.50,0.75,0.975,1.000))
   names(quant) = c("rqmin","rq025","rq250","rq500","rq750","rq975","rq100")
   #---------------------------------------------------------------------------------------#


   #----- Append all results, and standardise weird values to NA. -------------------------#
   ans   = c(three,quant,ntot=ntot,nval=nval,nrep=nr)
   ans[! is.finite(ans)] = NA
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end thirteen.err
#==========================================================================================#
#==========================================================================================#


#==========================================================================================#
#==========================================================================================#
#     Function to calculate the eddy covariance.                                           # 
#------------------------------------------------------------------------------------------#
eddy.cov <<- function(x,y,...){
   xbar   = mean(x,...)
   ybar   = mean(y,...)
   xprime = x - xbar
   yprime = y - ybar
   eddy   = mean(xprime*yprime,...)
   return(eddy)
}#end function eddy.cov
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     Function to find neat breaks given the resolution.                                   # 
#------------------------------------------------------------------------------------------#
neat.breaks <<- function(x,res,...){
   bad = 0
   #----- Make sure that x is a vector with contents. -------------------------------------#
   if ((! is.vector(x)) || (! is.numeric(x)) || length(x) < 2){
      cat("x must be a numeric vector of length 2 (or greater)!","\n",sep="")
      bad = bad + 1
   }else{
      #----- Remove infinite/NA entries. --------------------------------------------------#
      x = x[is.finite(x)]
      if (length(x) < 2){
         cat("x must contain at least 2 finite values!","\n",sep="")
         bad = bad + 1
      }#end if
      #------------------------------------------------------------------------------------#
   }#end if
   #----- Make sure that res is a simple scalar. ------------------------------------------#
   if ((! is.vector(res)) || (! is.numeric(res)) || length(res) != 1){
      cat("res must be a numeric scalar!","\n",sep="")
      bad = bad + 1
      if (res <= 0){
         cat("res must be positive!","\n",sep="")
         bad = bad + 1
      }#end if
   }#end if
   if (bad > 0) stop(" x and/or res are invalid!")
   #---------------------------------------------------------------------------------------#

   xrange   = range(x,finite=TRUE)
   xextreme = c(floor(xrange[1]/res),ceiling(xrange[2]/res)) * res
   breaks   = seq(from=xextreme[1],to=xextreme[2],by=res)
   return(breaks)
}#end function neat.breaks
#==========================================================================================#
#==========================================================================================#




#==========================================================================================#
#==========================================================================================#
#     This function is just like zoo's na.approx, except that it doesn't crash if the      #
# all values are NA.                                                                       #
#------------------------------------------------------------------------------------------#
na.approx.safe <<- function(object,...){
   ans = try(na.approx(object,...),silent=TRUE)
   if ("try-error" %in% is(ans)) ans = object
   return(ans)
}#end na.approx.safe
#==========================================================================================#
#==========================================================================================#




#==========================================================================================#
#==========================================================================================#
#     These functions are the same as cumsum and cumprod, except that it shifts the        #
# cumulative values to the value until right before the point.                             #
#------------------------------------------------------------------------------------------#
left.cumsum  <<- function(x) c(0,cumsum (x)[-length(x)])
left.cumprod <<- function(x) c(0,cumprod(x)[-length(x)])
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     This function aggregates x using some function "fun" as long as the number of valid  #
# entries is greater than or equal to fmin * total number of entries.  This is normally    #
# used by rasterize, hence the dots, but na.rm will not be used.                           #
#                                                                                          #
# INPUT:                                                                                   #
# x    - the vector with points to be aggregated                                           #
# fun  - the function.  It may be the actual function or a character indicating the        #
#        function                                                                          #
# fmin - minimum number of finite entries relative to the original size.                   #
#        This must be between 0 and 1.  fmin = 1 is equivalent to na.rm =FALSE, and        #
#        fmin = 0 is equivalent to na.rm = TRUE.                                           #
# ...  - additional arguments to be passed to function fun.  Note that na.rm will not      #
#        matter because non-finite numbers will be excluded prior to calling the function. #
#------------------------------------------------------------------------------------------#
aggr.fmin <<- function(x,fun=mean,fmin=0.5,...){
   #----- Check that fmin makes sense. ----------------------------------------------------#
   if (! (fmin %>=% 0.0 & fmin %<=% 1.0)){
      stop (paste0("fmin must be between 0 and 1!  Yours is set to ",fmin,"..."))
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Make sure fun is a function. ----------------------------------------------------#
   fun = match.fun(fun)
   #---------------------------------------------------------------------------------------#


   #----- Find the minimum number of entries to actually calculate the answer. ------------#
   keep  = is.finite(x)
   nkeep = sum(keep)
   ntot  = length(x)
   nmin  = max(1,min(ntot,ceiling(fmin * ntot)))
   #---------------------------------------------------------------------------------------#

   #----- Check whether to calculate function or return NA. -------------------------------#
   if (nkeep >= nmin){
      #----- Use only valid points. -------------------------------------------------------#
      xuse         = x[is.finite(x)]
      ans          = fun(xuse,...)
      discard      = ! is.finite(ans)
      ans[discard] = NA
      #------------------------------------------------------------------------------------#
   }else if (nkeep >= 1){
      #----- Use only valid points. -------------------------------------------------------#
      xuse = lit.sample(x=x[is.finite(x)],size=nmin,replace=TRUE)
      ans  = fun(xuse,...)
      ans  = rep(NA,times=length(ans))
      #------------------------------------------------------------------------------------#
   }else{
      #----- Make dummy vector, apply the function then discard data. ---------------------#
      xuse = rnorm(n=nmin)
      ans  = fun(xuse,...)
      ans  = rep(NA,times=length(ans))
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end aggr.fmin
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function aggregates standard error assuming Gaussian distribution and           #
# independent errors.  It will return a finite result provided that the number of valid    #
# entries is greater than or equal to fmin * total number of entries.  This is normally    #
# used by rasterize, hence the dots, but na.rm will not be used.                           #
#                                                                                          #
# INPUT:                                                                                   #
# x    - the vector with standard error to be aggregated                                   #
# fmin - minimum number of finite entries relative to the original size.                   #
#        This must be between 0 and 1.  fmin = 1 is equivalent to na.rm =FALSE, and        #
#        fmin = 0 is equivalent to na.rm = TRUE.                                           #
# ...  - additional arguments to be passed to function fun.  Note that na.rm will not      #
#        matter because non-finite numbers will be excluded prior to calling the function. #
#------------------------------------------------------------------------------------------#
aggr.se <<- function(x,fmin=0.5,...){
   #----- Check that fmin makes sense. ----------------------------------------------------#
   if (! (fmin %>=% 0.0 & fmin %<=% 1.0)){
      stop (paste0("fmin must be between 0 and 1!  Yours is set to ",fmin,"..."))
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Find the minimum number of entries to actually calculate the answer. ------------#
   keep  = is.finite(x)
   nkeep = sum(keep)
   ntot  = length(x)
   nmin  = max(1,min(ntot,ceiling(fmin * ntot)))
   #---------------------------------------------------------------------------------------#

   #----- Check whether to calculate function or return NA. -------------------------------#
   if (nkeep >= nmin){
      #----- Use only valid points. -------------------------------------------------------#
      xuse         = x[is.finite(x)]
      nuse         = length(xuse)
      ans          = sqrt(mean(xuse^2)/nuse)
      discard      = ! is.finite(ans)
      ans[discard] = NA
      #------------------------------------------------------------------------------------#
   }else{
      #----- Return NA. -------------------------------------------------------------------#
      ans  = NA
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end aggr.se
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     This function finds the maximum absolute elementwise difference of two vectors.      #
#------------------------------------------------------------------------------------------#
max.abs.diff <<- function(x,y,na.rm=TRUE) max(abs(x-y),na.rm=na.rm) 
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     This function counts the number of valid entries.                                    #
#------------------------------------------------------------------------------------------#
count.valid <<- function(x,qq.rm=FALSE){
   type.x = typeof(x)
   if (type.x %in% "logical"){
      ans = sum(! is.na(x))
   }else if (type.x %in% "character"){
      ans = sum(! ( is.na(x) | ((x %in% "") & qq.rm)))
   }else{
      ans = sum(is.finite(x))
   }#end if (type.x %in% c("logical","character"))
   return(ans)
}#end count.finite
#==========================================================================================#
#==========================================================================================#
