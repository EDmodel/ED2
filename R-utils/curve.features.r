#==========================================================================================#
#==========================================================================================#
#     This function determines the minima, maxima, and the "horizontal" and "vertical"     #
# inflection points for a given series.  It returns the indices of such features.          #
#                                                                                          #
# INPUTS:                                                                                  #
# x      -- vector to be analysed.                                                         #
# span   -- the window to be used                                                          #
# xscale -- scaling factor for x, so the sequence is normalised. If not provided, the      #
#           code will use by default the maximum absolute value in the series.             #
# toler  -- Tolerance for derivatives.  Values smaller than toler will be assumed zero.    #
#                                                                                          #
#------------------------------------------------------------------------------------------#
curve.features <<- function( x
                           , span   = 3L
                           , do.pad = TRUE
                           , xscale = max(abs(x),na.rm=TRUE)
                           , toler  = 1.e-4
                           ){
   #----- Make sure span is odd. ----------------------------------------------------------#
   span = as.integer(span)
   if ( ! ((span %% 2) %eq% 1 && span %ge% 3L)){
      stop(paste0(" Invalid span (",span,")! It must be an odd number (3 or greater)!"))
   }#end if ( ! ((span %% 2) %eq% 1))
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Make sure x is sufficiently long and has no NA.                                  #
   #---------------------------------------------------------------------------------------#
   nx = length(x) 
   if (nx < 4){
      stop(" Vector 'x' must have at least 4 elements!")
   }else if (! all(is.finite(x))){
      stop(" Vector 'x' cannot have non-finite elements!")
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Make sure that x can be normalised by the given scale.                            #
   #---------------------------------------------------------------------------------------#
   if (xscale %gt% 0.){
      x = x / xscale
   }else{
      stop(paste0("Invalid xscale (",xscale,").  It must be positive."))
   }#end if 
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Find the derivatives.                                                            #
   #---------------------------------------------------------------------------------------#
   i       = sequence(nx)
   ip1     = pmin(i+1,nx)
   im1     = pmax(i-1, 1)
   #----- First derivative. ---------------------------------------------------------------#
   xp      = x[ip1] - x[im1]
   xpscale = max(abs(xp),na.rm=TRUE)
   if (xpscale %gt% 0.){
      xp = ifelse(test= abs( xp / xpscale ) %ge% toler, yes= xp / xpscale, no = 0.)
   }else{
      xp = 0. * xp
   }#end if (xpscale %gt% 0.)
   #----- Second derivative. --------------------------------------------------------------#
   xpp      = x[ip1] - 2.*x[i] + x[im1]
   xpp[1]   = xpp[2]
   xpp[nx]  = xpp[nx-1] 
   xppscale = max(abs(xpp),na.rm=TRUE)
   if (xppscale %gt% 0.){
      xpp = ifelse(test= abs( xpp / xppscale ) %ge% toler, yes= xpp / xppscale, no = 0.)
   }else{
      xpp = 0. * xpp
   }#end if (xppscale %gt% 0.)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Find the padding length.                                                         #
   #---------------------------------------------------------------------------------------#
   smid = ceiling(span/2)
   soff = smid - 1
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Create an zero matrix to pad the answers.                                        #
   #---------------------------------------------------------------------------------------#
   zero    = matrix(data=0.,nrow=soff,ncol=span)
   x.mat   = t(apply(X=rbind(zero,embed(x=x  ,dimension=span),zero),MARGIN=1,FUN=rev))
   xp.mat  = t(apply(X=rbind(zero,embed(x=xp ,dimension=span),zero),MARGIN=1,FUN=rev))
   xpp.mat = t(apply(X=rbind(zero,embed(x=xpp,dimension=span),zero),MARGIN=1,FUN=rev))
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Fix edges so they can become maxima or minima in case the curve is not zero.     #
   #---------------------------------------------------------------------------------------#
   xp.mat  [ 1,] = xp.mat [   2,]
   xp.mat  [nx,] = xp.mat [nx-1,]
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Find the signs to the left, to the right and for the entire neighbourhood of     #
   # each vector element.                                                                  #
   #---------------------------------------------------------------------------------------#
   xp.left     = apply( X = xp.mat , MARGIN = 1, FUN = signblock, block = "left"  )
   xp.right    = apply( X = xp.mat , MARGIN = 1, FUN = signblock, block = "right" )
   xp.both     = apply( X = xp.mat , MARGIN = 1, FUN = signblock, block = "all"   )
   xpp.left    = apply( X = xpp.mat, MARGIN = 1, FUN = signblock, block = "left"  )
   xpp.right   = apply( X = xpp.mat, MARGIN = 1, FUN = signblock, block = "right" )
   xpp.both    = apply( X = xpp.mat, MARGIN = 1, FUN = signblock, block = "all"   )
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Find features.                                                                   #
   #---------------------------------------------------------------------------------------#
   ans = data.frame( max = xp.left*xp.right   == -1 & xpp.both         == -1
                   , min = xp.left*xp.right   == -1 & xpp.both         == +1
                   , iph = xpp.left*xpp.right == -1 & xpp.left*xp.both == -1
                   , ipv = xpp.left*xpp.right == -1 & xpp.left*xp.both == +1
                   )#end data.frame
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    It is often the case that neighbouring cells are selected. In this case we         #
   # only keep one in the neighbourhood.                                                   #
   #---------------------------------------------------------------------------------------#
   #---- Maximum. Pick the maximum in the immediate vicinity. -----------------------------#
   imax = which(ans$max)
   if (length(imax) > 1){
      xtest = NULL
      for (o in seq(from=-soff,to=soff,by=1)){
         xoff  = x[imax+o]
         xoff  = ifelse(test=(imax+o) %in% imax,yes=xoff,no=-Inf)
         xtest = cbind(xtest,xoff)
      }#end for (o in seq(from=-soff,to=soff,by=1))
      xrefr                = apply(X=xtest[,-smid],MARGIN=1,FUN=max)
      iskip                = xtest[,smid] < xrefr
      ans$max[imax[iskip]] = FALSE
   }#end if (length(imax) > 1)
   #---- Minimum. Pick the minimum in the immediate vicinity. -----------------------------#
   imin = which(ans$min)
   if (length(imin) > 1){
      xtest = NULL
      for (o in seq(from=-soff,to=soff,by=1)){
         xoff  = x[imin+o]
         xoff  = ifelse(test=(imin+o) %in% imin,yes=xoff,no=+Inf)
         xtest = cbind(xtest,xoff)
      }#end for (o in seq(from=-soff,to=soff,by=1))
      xrefr                = apply(X=xtest[,-smid],MARGIN=1,FUN=min)
      iskip                = xtest[,smid] > xrefr
      ans$min[imin[iskip]] = FALSE
   }#end if (length(imin) > 1)
   #---- Horizontal inflection point.  Pick the flattest point in the immediate vicinity. -#
   iiph = which(ans$iph)
   if (length(iiph) > 1){
      xptest = NULL
      for (o in seq(from=-soff,to=soff,by=1)){
         xpoff  = xp[iiph+o]
         xpoff  = ifelse(test=(iiph+o) %in% iiph,yes=xpoff,no=+Inf)
         xptest = cbind(xptest,xpoff)
      }#end for (o in seq(from=-soff,to=soff,by=1))
      xprefr               = apply(X=abs(xptest[,-smid]),MARGIN=1,FUN=min)
      iskip                = abs(xptest[,smid]) > xprefr
      ans$iph[iiph[iskip]] = FALSE
   }#end if (length(iiph) > 1)
   #---- Vertical inflection point.  Pick the steepest point in the immediate vicinity. ---#
   iipv = which(ans$ipv)
   if (length(iipv) > 1){
      xptest = NULL
      for (o in seq(from=-soff,to=soff,by=1)){
         xpoff  = xp[iipv+o]
         xpoff  = ifelse(test=(iipv+o) %in% iipv,yes=xpoff,no=0.)
         xptest = cbind(xptest,xpoff)
      }#end for (o in seq(from=-soff,to=soff,by=1))
      xprefr               = apply(X=abs(xptest[,-smid]),MARGIN=1,FUN=max)
      iskip                = abs(xptest[,smid]) < xprefr
      ans$ipv[iipv[iskip]] = FALSE
   }#end if (length(iipv) > 1)
   #---------------------------------------------------------------------------------------#




   #----- Flag edges as inflection points in case their first derivative is not zero. -----#
   ans$iph[ 1] = ans$max[ 1] || (xp.right[ 1] == -1)
   ans$iph[nx] = ans$max[nx] || (xp.left [nx] == +1)
   ans$ipv[ 1] = ans$max[ 1] || (xp.right[ 1] == +1)
   ans$ipv[nx] = ans$max[nx] || (xp.left [nx] == -1)
   #---------------------------------------------------------------------------------------#


   #------ Return answer. -----------------------------------------------------------------#
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end curve.features
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function finds the sign of a block, Either the left or the right side of the    #
# branch, or the entire vector.  Possible values are:                                      #
#  -1: all elements are negative                                                           #
#  +1: all elements are positive                                                           #
#   0: positive and negative elements were found.                                          #
#  NA: at least one element was not finite.                                                #
#------------------------------------------------------------------------------------------#
signblock <<- function(x,block=c("all","left","right")){
   #----- Make sure block is a valid entry. -----------------------------------------------#
   block=match.arg(block)
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #     Define the indices to look at.                                                    #
   #---------------------------------------------------------------------------------------#
   if (length(x) == 0){
      xuse = NA
   }else if (block %in% "all" || length(x) == 1){
      xuse = x[-ceiling(length(x)/2)]
   }else{
      #----- Grab indices that will be used. ----------------------------------------------#
      idx = sequence(floor(length(x)/2))
      #------------------------------------------------------------------------------------#


      #----- Find the slice to analyse. ---------------------------------------------------#
      if (block %in% "left"){
         xuse = x[idx]
      }else{
         xuse = rev(rev(x)[idx])
      }#end if (block %in% "left")
      #------------------------------------------------------------------------------------#
   }#end if (length(x) == 0)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Add the signs and divide by the length.  In case there is any sign change, the    #
   # result is always going to be zero. Also, if there is any NA, the answer will be NA.   #
   #---------------------------------------------------------------------------------------#
   ans = trunc(sum(sign(xuse))/length(xuse))
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end signblock
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     Function that defines the left boundaries between cohorts.                           #
#------------------------------------------------------------------------------------------#
left.bnd <<- function(k,ipk,imn,iiv){

   #---------------------------------------------------------------------------------------#
   #     Decide how to define the boundary.                                                #
   #---------------------------------------------------------------------------------------#
   if (k == 1){
      #----- If k is 1, then the left boundary must be 1 (first cohort). ------------------#
      ans = 1
      #------------------------------------------------------------------------------------#
   }else{
      #------------------------------------------------------------------------------------#
      #    We give preference to minima, if none exists, then use the vertical             #
      # inflection point.                                                                  #
      #------------------------------------------------------------------------------------#
      o = which(imn > ipk[k-1] & imn < ipk[k])
      if (length(o) == 0){
         o = which (iiv > ipk[k-1] & iiv < ipk[k])
         if (length(o) == 0){
            ans = round(mean(ipk[k-1],ipk[k]))
         }else{
            ans = closest(x=mean(ipk[k-1],ipk[k]),A=iiv[o])
         }#end if (length(ans) == 0)
         #---------------------------------------------------------------------------------#
      }else{
         ans = closest(x=mean(ipk[k-1],ipk[k]),A=imn[o])
      }#end if (length(ans) == 0)
      #------------------------------------------------------------------------------------#
   }#end if (k == 1){
   #---------------------------------------------------------------------------------------#

   #----- Return the answer. --------------------------------------------------------------#
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function left.bnd
#==========================================================================================#
#==========================================================================================#
