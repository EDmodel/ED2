#==========================================================================================#
#==========================================================================================#
#      Polar.periodogram.  This function applies the spatial Fourier analysis, then        #
# computes the periodogram function for different radii and theta.  The result is a data   #
# frame with three variables:                                                              #
#                                                                                          #
# radius - The wave number magnitude (in pixels)                                           #
# rsize  - Size to bin wave numbers (in pixels)                                            #
# thsize - Size to bin wave directions (in degrees)                                        #
#------------------------------------------------------------------------------------------#
polar.periodogram <<- function( X
                              , which.bin  = c("radius","length","theta"
                                              ,"angle","both","none")
                              , hemisphere = TRUE
                              , thsize     = 10
                              , thmin      = ifelse(hemisphere,0,-180)
                              , thmax      = 180-thsize
                              , rmin       = 3
                              , rmax       = 0.5
                              ){

   #---------------------------------------------------------------------------------------#
   #     Make sure X is a matrix.                                                          #
   #---------------------------------------------------------------------------------------#
   err = FALSE
   if ( is.array(X)){
      if (length(dim(X)) != 2){
         err = TRUE
      }#end if (length(dim(X)) != 2)
   }else if(! (is.matrix(X) || is.data.frame(X))){
      err = TRUE
   }else if (is.data.frame(X)){
      #----- Coerce X to a matrix. --------------------------------------------------------#
      X  = as.matrix(X)
      #------------------------------------------------------------------------------------#
   }#end if
   if (err) stop("X must be a matrix, a data frame, or an array with dimension 2")
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Find the appropriate binning term.                                                #
   #---------------------------------------------------------------------------------------#
   which.bin = match.arg(which.bin)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Apply the Fourier transform.  In R one must divide by the matrix size to obtain   #
   # the normalised components.                                                            #
   #---------------------------------------------------------------------------------------#
   Xb   = mean(c(unlist(X)),na.rm=TRUE)
   Xp   = X - Xb
   mm   = nrow (X)
   nn   = ncol (X)
   FF   = fft(Xp) / mm / nn
   II   = abs(FF) * abs(FF) * mm * nn
   sig2 = var(c(Xp),na.rm=TRUE)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Find the directional wave numbers (p and q).                                      #
   #---------------------------------------------------------------------------------------#
   np      = floor(mm/2)
   nq      = floor(nn/2)
   pp      = sequence(mm)-1
   qq      = sequence(nn)-1
   pp      = ifelse(pp > np, pp-mm, pp)
   qq      = ifelse(qq > nq, qq-nn, qq)
   op      = order(pp)
   oq      = order(qq)
   pp      = pp[op]
   qq      = qq[oq]
   FF      = FF[op,oq]
   II      = II[op,oq]
   PP      = matrix(data = rep(pp,times=nn),nrow=mm,ncol=nn)
   QQ      = matrix(data = rep(qq,each =mm),nrow=mm,ncol=nn)
   RR      = sqrt(mm*mm*PP*PP+nn*nn*QQ*QQ) / mm / nn
   LL      = 1 / RR
   THTH    = 180 * atan2(y=QQ/nn,x=PP/mm) / pi
   THTH    = ifelse(THTH > (180-0.5*thsize),THTH-360,THTH)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     We only keep the dimensions that matter.  It doesn't make sense to analyse theta  #
   # outside the 0-180 range, because the wave forms are symmetric.  Also, radii larger    #
   # than the window width are not properly sampled so they should be discarded.  In the   #
   # case of rectangular matrices, we take a conservative approach and pick the smallest   #
   # of them.                                                                              #
   #---------------------------------------------------------------------------------------#
   rsize       = 3 / max(mm,nn)
   keep.theta  = THTH > (thmin - 0.5 * thsize) & THTH <= (thmax + 0.5 * thsize)
   keep.radius = RR   > rsize & RR <= 0.5
   keep        = keep.theta & keep.radius
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #   Split radius and angles into bins.                                                  #
   #---------------------------------------------------------------------------------------#
   RR.keep   = RR  [keep]
   LL.keep   = LL  [keep]
   THTH.keep = THTH[keep]
   II.keep   = II  [keep]
   th.brks   = seq(from=thmin-0.5*thsize,to=thmax+0.5*thsize,by=thsize)
   r.brks    = seq(from=rsize,to=rsize*ceiling(0.5/rsize),by=rsize)
   R.bin     = mid.points(r.brks)
   R.cut     = cut(RR.keep,breaks=r.brks,include.lowest=FALSE,labels=R.bin)
   TH.bin    = mid.points(th.brks)
   TH.cut    = cut(THTH.keep,breaks=th.brks,include.lowest=FALSE,labels=TH.bin)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Find the variance of the retained spectrum.                                       #
   #---------------------------------------------------------------------------------------#
   sig2p = sum(II.keep)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Decide how to aggregate data.                                                    #
   #---------------------------------------------------------------------------------------#
   if (which.bin %in% "radius"){
      m.II = tapply(X=II.keep,INDEX=R.cut,FUN=mean,na.rm=TRUE)
      ans  = m.II / sig2p
   }else if (which.bin %in% c("theta","angle")){
      m.II = tapply(X=II.keep,INDEX=TH.cut,FUN=mean,na.rm=TRUE)
      ans  = m.II / sig2p
   }else if (which.bin %in% "both"){
      m.II = tapply(X=II.keep,INDEX=list(R.cut,TH.cut),FUN=mean,na.rm=TRUE)
      ans  = m.II / sig2p
   }else if (which.bin %in% "length"){
      m.II       = tapply(X=II.keep,INDEX=R.cut,FUN=mean,na.rm=TRUE)
      ans        = m.II / sig2p
      len        = tapply(X=LL.keep,INDEX=R.cut,FUN=mean,na.rm=TRUE)
      names(ans) = len
   }else if (which.bin %in% "none"){
      ans = list( x = pp, y = qq, z = ifelse(PP == 0 & QQ == 0,NA,II))
   }else{
      stop("Not sure how to aggregate stuff")
   }#end if
   #---------------------------------------------------------------------------------------#


   #------ Return matrix or list. ---------------------------------------------------------#
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end polar.periodogram
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      This function is a wrapper for the FOTO package, when using matrices as opposed to  #
# rasters. It internally creates a raster that is then sent to the package function.  This #
# copies part of the steps in the function "foto" because the package by default applies   #
# the principal component analysis and this is done separately.                            #
#------------------------------------------------------------------------------------------#
foto.periodogram <<- function(X){
   #------ Create a temporary raster to send to the FOTO function. ------------------------#
   myenv  = environment()
   nx     = nrow(X)
   ny     = ncol(X)
   yswap  = rev(sequence(ny))
   Xt     = t(X[,yswap])
   img    = raster(Xt)
   wsize  = min(nx,ny)/1
   N      = ceiling(img@nrows/wsize)
   M      = ceiling(img@ncols/wsize)
   cells  = N * M
   #---------------------------------------------------------------------------------------#



   #------ Run the FOTO function for this raster. Single window using the entire image. ---#
   Imat         = matrix(img,wsize,wsize)
   fftim        = Mod(stats::fft(Imat))^2
   offfft       = ceiling(dim(Imat)[1]/2)
   r            = sqrt((col(Imat) - offfft)^2 + (row(Imat) - offfft)^2)
   rzonal       = raster::zonal( x     = raster::raster(fftim)
                               , z     = raster::raster(r)
                               , fun   = "mean"
                               , na.rm = TRUE
                               )#end raster::zonal
   rspec        = rev(rzonal[, 2])
   names(rspec) = rev(rzonal[, 1])
   sigma        = sd(Imat,na.rm=TRUE)
   rscal        = rspec / sigma

   iuse         = rev(sequence(wsize/2)[-c(1,2)])
   ans          = rscal[iuse]
   #---------------------------------------------------------------------------------------#


   #----- Return answer. ------------------------------------------------------------------#
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end foto.periodogram
#==========================================================================================#
#==========================================================================================#
