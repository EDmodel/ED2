#==========================================================================================#
#==========================================================================================#
#     Function that finds the cubic root, for both positive and negative numbers.          #
#------------------------------------------------------------------------------------------#
cbrt <<- function(x){
   x333      = sign(x) * abs(x)^onethird
   bad       = is.nan(x333)
   x333[bad] = NA
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
#     This function estimates the quantile for a table of observations x, each of which    #
# having a weight w.  This is done by finding the median of a pseudo dataset built using   #
# sample.  If the size of the resampling is not provided, then the number of samples is    #
# dependent on the range of probabilities.  By default we find the 0.50 quantile (median). #
#------------------------------------------------------------------------------------------#
weighted.quantile <<- function(x,w,qu=0.50,size.minp=10,na.rm=FALSE){

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


   ans = list(q = qout, case = case)
   return(ans)
}#end function weighted.quantile
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function estimates the weighted standard deviation.                             #
#------------------------------------------------------------------------------------------#
weighted.sd <<- function(x,w,na.rm=FALSE){

   #----- Delete the missing values if the user asked to do it. ---------------------------#
   if (any(w < 0, na.rm = TRUE) || any(is.infinite(w)) || any(is.na(w))){
      stop(" Weights (w) must be non-negative and finite, and entirely defined!")
   }else if(na.rm){
      keep = ! is.na(x)
      x    = x[keep]
      w    = w[keep]
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Check whether at least one weight is non-zero.                                   #
   #---------------------------------------------------------------------------------------#
   if (all(w %==% 0)){
      ans = NA
   }else{
      xwm    = weighted.mean(x=x,w=w)
      M      = sum(w %>% 0)
      w.sum  = sum(w)
      r2.sum = sum(w*(x-xwm)^2)
      ans    = sqrt(M * r2.sum / ( (M-1) * w.sum))
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
      ans     = sum((xx-xx.mean)^4) / (nx * xx.sdev^4) - 3
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
      msg = "four.moment(<matrix>) is deprecated.\n Use apply(*, 2, kurt) instead."
      warning(paste(msg, collapse = ""), call. = FALSE, domain = NA)
      ans = apply(X = x, MARGIN = 2, FUN = kurt, na.rm = na.rm)
   }else if (is.list(x)){
      msg = "four.moment(<list>) is deprecated.\n Use lapply(*, kurt) instead."
      warning(paste(msg, collapse = ""), call. = FALSE, domain = NA)
      ans = lapply(X = x, FUN = kurt, na.rm = na.rm)
   }else if (is.data.frame(x)){
      msg = "four.moment(<data.frame>) is deprecated.\n Use sapply(*, kurt) instead."
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



      #----- Find the skewness. -----------------------------------------------------------#
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
