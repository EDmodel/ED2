#==========================================================================================#
#==========================================================================================#
#     Function that finds the cubic root, for both positive and negative numbers.          #
#------------------------------------------------------------------------------------------#
cbrt <<- function(x){
   x333 = NA * x

   pos  = is.finite(x) & x >= 0.
   neg  = is.finite(x) & x <  0.
   
   x333[pos] =  ( x[pos])^(1./3.)
   x333[neg] = -(-x[neg])^(1./3.)

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
#     This function finds the binary notation of an integer.                               #
#------------------------------------------------------------------------------------------#
tobin <<- function(x,reverse=TRUE){

   if (! is.integer(x)){
      stop ("x must be integer")
   }#end if

   int2bin = function(x,reverse){
      if (reverse){
        ans = paste(sapply(strsplit(paste(rev(intToBits(x))),""),`[[`,2),collapse="")
      }else{
        ans = paste(sapply(strsplit(paste(intToBits(x)),""),`[[`,2),collapse="")
      }#end if
      return(ans)
   }#end if

   if (is.list(x)){
      ans = sapply(X=x,FUN=int2bin,reverse=reverse)
   }else if (is.null(dim(x))){
      ans = sapply(X=x,FUN=int2bin,reverse=reverse,simplify=TRUE)
   }else{
      margin = length(dim(x))
      ans    = apply(X=x,MARGIN=sequence(margin),FUN=int2bin,reverse=reverse)
   }#end if
   return(ans)
}#end function tobin
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     This function computes the mean of the elements below (above) a given quantile.      #
#------------------------------------------------------------------------------------------#
qu.mean = function(x,p,na.rm=FALSE,lower=TRUE){
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
weighted.quantile = function(x,w,qu=0.50,size.minp=10,na.rm=FALSE){

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
#        Pick the element in array A that has the closest value to x.                      #
#------------------------------------------------------------------------------------------#
which.closest <<- function(x,A,mask=rep(TRUE,length(A)) ){

   #----- Create indices before applying the mask, so we retrieve the actual index. -------#
   A.idx = sequence(length(A))
   A     = A    [mask]
   A.idx = A.idx[mask]
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Check whether there is any data that 
   #---------------------------------------------------------------------------------------#
   if (! is.finite(x) || ! any(mask,na.rm=TRUE)){
      #----- x is invalid, return NA. -----------------------------------------------------#
      idx = NA
   }else if(sum(x == A,na.rm=TRUE) > 1){
      #------------------------------------------------------------------------------------#
      #      If there are multiple values of x that are the same as x, we randomly sample  #
      # one value.                                                                         #
      #------------------------------------------------------------------------------------#
      idx = sample(A.idx[x == A],size=1)
      #------------------------------------------------------------------------------------#
   }else{
      #------------------------------------------------------------------------------------#
      #      Either there is a single value that matches, or none of them match.  Either   #
      # way, we select the closest match.                                                  #
      #------------------------------------------------------------------------------------#
      idx = A.idx[which.min((A-x)^2)]
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#

   return(idx)
}#end which.closest
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



      #----- Find the skewness. -----------------------------------------------------------#
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
#     This function creates a sequence of integers that has the same length as the length  #
# of the argument.                                                                         #
#------------------------------------------------------------------------------------------#
seq.len <<- function(x){
   ans = sequence(length(x))
   return(ans)
}#end function lit.sample
#==========================================================================================#
#==========================================================================================#
