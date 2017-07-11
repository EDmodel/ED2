#==========================================================================================#
#==========================================================================================#
#      This function fills all gaps in a time series, by using an iterative Fourier        #
# analysis of the time series.  Keep in mind that geophysical time series have random      #
# noise that are captured by the Fourier transform as weak waves, and those waves can make #
# very bad things in the middle of the gap if the gap is long.                             #
#                                                                                          #
# Input Variables                                                                          #
# ---------------                                                                          #
#                                                                                          #
# x:                                                                                       #
#    The time series to be filled.                                                         #
#                                                                                          #
# detrend.method:                                                                          #
#    Which technique to use to detrend the time series before applying the Fourier trans-  #
#    form.  Valid values are (case insensitive and the first few letters will do):         #
#    - "mean"   : don't detrend the time series, just subtract the mean                    #
#    - "linear" : (default) use a simple linear detrending                                 #
#    - "loess"  : use a polynomial surface using local fitting                             #
#                                                                                          #
# trend.back (optional):                                                                   #
#    Should the routine add back the trend? (TRUE or FALSE).  Default is TRUE.             #
#                                                                                          #
# min.signal (optional):                                                                   #
#    The minimum accumulated spectrum to retain (the count goes from most powerful to      #
#    least powerful).  It must be between 0 and signal.retain (see below).                 #
#                                                                                          #
# signal.retain (optional):                                                                #
#    The total accumulated spectrum to retain (the count goes from most powerful to least  #
#    powerful).  It must be between 0 and 1, and at least one mode will always be used.    #
#                                                                                          #
# conv.threshold:                                                                          #
#    Tolerance for change in the response for each sub-step, beyond which we move to the   #
#    next mode.                                                                            #
#                                                                                          #
# minmod:                                                                                  #
#    The minimum number of modes to use.  This avoids using too few modes when the strong- #
#    est modes are too powerful.                                                           #
#                                                                                          #
# maxmod:                                                                                  #
#    The maximum number of modes to use.  This allows using all sought modes when the      #
#    spectrum is well-defined, and lower power in case the time series is too noisy.       #
#                                                                                          #
# maxin:                                                                                   #
#    Maximum number of inner iterations before giving up convergence and moving on.        #
#                                                                                          #
# verbose:                                                                                 #
#    Prints more information, which may be useful if you want to debug or just curious.    #
#                                                                                          #
# rmse:                                                                                    #
#    Flag that tells whether to estimate the root mean square error (TRUE | FALSE).  If    #
#    true, a jackknife method will be run, otherwise a NA will be returned for error.      #
#                                                                                          #
# del.frac:                                                                                #
#    In case jackknife is to be run, this tells the fraction of valid data to be removed   #
#    each realisation.                                                                     #
#                                                                                          #
# n.jack:                                                                                  #
#    Maximum number of iterations for the jackknife method.                                #
#                                                                                          #
# Output Variables                                                                         #
# -----------------                                                                        #
#                                                                                          #
#     The output is a list containing the following variables:                             #
#                                                                                          #
# xfill:                                                                                   #
#    vector with the same length as y with the original time series where data were        #
#    available, and the gap filled value for the time series.                              #
#                                                                                          #
# error:                                                                                   #
#    Estimate of the root mean square error.  If input rmse is FALSE, this will be a NA.   #
#                                                                                          #
#------------------------------------------------------------------------------------------#
harmfill <<- function(x,detrend.method="linear",trend.back=TRUE,min.signal=0.00
                     ,signal.retain=0.80,conv.threshold=0.0001,minmod=1,maxmod=Inf
                     ,verbose=0,maxin=50,rmse=FALSE,del.frac=1/3,n.jack=100
                     ,jack.toler=0.01){

   #---------------------------------------------------------------------------------------#
   #      Harmfill requires two packages now: RSEIS and zoo.  Make sure that you have both #
   # of them installed and loaded.                                                         #
   #---------------------------------------------------------------------------------------#
   zoo.check   = "package:zoo"   %in% search()
   RSEIS.check = "package:RSEIS" %in% search() | detrend.method != "linear"
   if ( (! zoo.check) | (! RSEIS.check) ){
      cat ("    ---> ZOO:   ",zoo.check  ,"\n")
      cat ("    ---> RSEIS: ",RSEIS.check,"\n")
      stop("    Harmfill requires ZOO and RSEIS!!!")
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Copy the time series to a safe place, and never point to x again unless there is  #
   # no need for gap filling.                                                              #
   #---------------------------------------------------------------------------------------#
   x.act    = x
   nx.act   = length(x.act)
   miss.act = is.na(x.act)
   nx.inf   = sum(is.infinite(x.act),na.rm=TRUE)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Stop if there is any infinity...                                                   #
   #---------------------------------------------------------------------------------------#
   if (nx.inf != 0){
      cat ("   ---> Matt, there are infinite numbers in your time series, fix it!!!","\n")
      cat ("   ---> Number of points:                   ",nx.act,"\n")
      cat ("   ---> Number of points that are infinity: ",nx.inf,"\n")
      stop("   Time series must contain only finite numbers and NAs")
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Double the size of the time series by adding the second half before and the first #
   # half after the time series.  Harmonic filling doesn't work very well at the beginning #
   # and end of the time series, so we buffer the period.                                  #
   #---------------------------------------------------------------------------------------#
   n.brk    = floor(nx.act / 2)
   before   = sequence(n.brk)
   after    = seq(from=n.brk+1,to=nx.act,by=1)
   rxa      = length(after) + 1
   rxz      = rxa + nx.act - 1
   x.ext    = c(x.act[after],x.act,x.act[before])
   miss.ext = is.na(x.ext)
   nx.ext   = length(x.ext)
   #---------------------------------------------------------------------------------------#


   #----- Find the time series size and the Nyquist frequency. ----------------------------#
   nnyq = 1 + ceiling((nx.ext-1)/2)
   #---------------------------------------------------------------------------------------#


   #----- Find the indices to use in the FFT. ---------------------------------------------#
   fuse = sequence(nnyq)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Check and count the number of missing data. If the time series is full, we don't  #
   # need to do anything.                                                                  #
   #---------------------------------------------------------------------------------------#
   okdata =   is.finite(x.ext)
   nodata = ! is.finite(x.ext)
   nok    = sum(okdata)
   nmiss  = sum(nodata)
   if (nmiss == 0){
      if (verbose > 0) cat ("    * Time series is complete, no need to gap fill... \n")
      ans  = list(xfit=x.act,error=rep(0,times=nx.act))
      return(ans)
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Fill in the gaps with something so it doesn't go completely off during the        #
   # de-trending.  To avoid spurious trends, we keep only the true period.                 #
   #---------------------------------------------------------------------------------------#
   x.trend = na.fill(na.approx(x.ext,na.rm=FALSE),fill="extend")
   x.trend = x.trend[rxa:rxz]
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     De-trend the time series.  Notice that even though we                                                        #
   #---------------------------------------------------------------------------------------#
   if (verbose > 0) cat ("    * Detrending time series... \n")

   detmet  = substring(tolower(detrend.method),1,2)
   if (detmet == "li"){
      #----- Linear trend. ----------------------------------------------------------------#
      x0         = detrend2(x.trend)$y0
      #------------------------------------------------------------------------------------#

      #----- Make x0 the same size as x.ext -----------------------------------------------#
      x0   = c(x0[after],x0,x0[before])
      #------------------------------------------------------------------------------------#
   }else if(detmet == "lo"){
      #----- When is just a dimensionless time. -------------------------------------------#
      when  = sequence(nx.act)
      #------------------------------------------------------------------------------------#

      #----- Local fitting. ---------------------------------------------------------------#
      guess = loess(formula= x ~ when, data=data.frame(when=when,x=x.trend)
                   ,na.action="na.omit")
      x0    = predict(object=guess,when)
      #------------------------------------------------------------------------------------#

      #----- Make x0 the same size as x.ext -----------------------------------------------#
      x0   = c(x0[after],x0,x0[before])
      #------------------------------------------------------------------------------------#
   }else if (detmet == "me"){
      #----- No trend, subtract the mean and that's enough. -------------------------------#
      x0 = rep(mean(x.act,na.rm=TRUE),times=nx.ext)
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Find the deviation from the detrended time series. ------------------------------#
   xprime = x.ext - x0
   xpbar  = mean(xprime,na.rm=TRUE)
   #---------------------------------------------------------------------------------------#


   #----- Initial conditions.  Assume no wave present in missing data. --------------------#
   xnext         = xprime
   xnext[nodata] = 0.
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Find the first Fourier transform, and determine which modes to use.              #
   #---------------------------------------------------------------------------------------#
   xfftall           = fft(xnext)
   xfft              = xfftall[fuse]
   pow               = abs(xfft)^2
   npow              = order(pow,decreasing=TRUE)
   cumpow            = cumsum(pow[npow])/sum(pow)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Noutit is the number of "outer" iteractions needed to capture the meaningful      #
   # power.                                                                                #
   #---------------------------------------------------------------------------------------#
   noutit  = max(minmod,max(sum(cumpow <= min.signal)
                           ,min(nok-minmod,sum(cumpow <= signal.retain),maxmod)))
   powuse  = npow[1:noutit]
   if (verbose > 0){
      cat("    * Using ",noutit," modes out of ",nnyq,"... \n")
      cat("      ( Retained signal = ",sprintf("%.2f",cumpow[noutit]*100),"%...) \n")
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the list of iterations to print.                                             #
   #---------------------------------------------------------------------------------------#
   print.iter = unique(sort(c(1,pretty(sequence(noutit),n=10),noutit)))
   print.iter = print.iter[print.iter %in% sequence(noutit)]
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Loop over the outer iterations.                                                   #
   #---------------------------------------------------------------------------------------#
   if (verbose > 0) cat ("    * Starting the outer loop... \n")
   for (outit in 1:noutit){

      if (verbose > 0 & outit %in% print.iter){
         cat("      # Outer iteration: ",outit
            ,".  Signal = ",sprintf("%.3f",cumpow[outit]*100),"... \n")
      }#end if
      #----- Reset the values for the inner iterations. -----------------------------------#
      initerate   = TRUE
      #errest1st  = NULL
      #errestbest = NULL
      r2best      = 0
      r2prev      = 0
      init        = 0


      #------------------------------------------------------------------------------------#
      #     Conditional loop over the outer iterations.                                    #
      #------------------------------------------------------------------------------------#
      while(initerate){
         init = init + 1

         #----- Update guess. -------------------------------------------------------------#
         xnow = xnext
         #---------------------------------------------------------------------------------#



         #----- Find the Fast Fourier analysis for this guess. ----------------------------#
         xfftall = fft(xnow)
         xfft    = xfftall[fuse]
         #---------------------------------------------------------------------------------#


         #----- Keep only the powers that we should use. ----------------------------------#
         del = - powuse[1:outit]
         xfft[del] = 0+0i
         #---------------------------------------------------------------------------------#


         #----- Reconstruct the Fourier transform without the weaker powers. --------------#
         if (nx.ext %% 2 == 0){
            dseq              = seq(from=2,to=nnyq-1,by=1)
            xfft              = c(xfft[1]
                                 ,xfft[dseq]
                                 ,xfft[nnyq]
                                 ,rev(Re(xfft[dseq])) + (0-1i)*rev(Im(xfft[dseq])) )
         }else{
            dseq              = seq(from=2,to=nnyq,by=1)
            xfft              = c(xfft[1]
                                 ,xfft[dseq]
                                 ,rev(Re(xfft[dseq])) + (0-1i)*rev(Im(xfft[dseq])) )
         }#end if
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #    Find the guess by finding the inverse FFT without the trailing components.   #
         #---------------------------------------------------------------------------------#
         xfill             = Re(fft(xfft,inverse=TRUE)/nx.ext)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Find the potential next guess, which is the actual time series, with the    #
         # missing values replaced by the guess.                                           #
         #---------------------------------------------------------------------------------#
         xtry         = xprime
         xtry[nodata] = xfill[nodata]
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Simple error guess, this will tell whether the model is approaching a      #
         # place.  Although this may mean slow convergence, we assume that this is because #
         # the result is close to the best guess.                                          #
         #---------------------------------------------------------------------------------#
         r2 = ( 1. - ((nok - 1        ) * sum((xprime[okdata] - xfill[okdata])^2))
                   / ((nok - outit - 1) * sum((xprime[okdata] - xpbar        )^2)) )
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Decide whether to accept or reject the step.  We iterate in the inner loop  #
         # only if guesses are getting better by a significant amount.                     #
         #---------------------------------------------------------------------------------#
         gain       = 2.0 * (r2 - r2prev) / (r2 + r2prev)
         r2prev     = r2
         initerate  = init < maxin &&  gain > conv.threshold
         xnext      = xtry
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #    Print information if needed.                                                 #
         #---------------------------------------------------------------------------------#
         if (verbose > 1){
            cat("        > Inner iteration: ",init
               ,"; R2 = ",signif(r2,5),"; Gain = ",signif(gain,5),"...","\n")
         }#end if
         #---------------------------------------------------------------------------------#
      }#end while (initerate)
      #------------------------------------------------------------------------------------#
   }#end for outit in 1:noutit
   #---------------------------------------------------------------------------------------#


   #----- Add back the trend, and chop the time series back to the original size. ---------#
   if (trend.back){
      xfit = x0[rxa:rxz] + xnext[rxa:rxz]
   }else{
      xfit = xnext[rxa:rxz]
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     This part is the error estimate.                                                  #
   #---------------------------------------------------------------------------------------#
   if (rmse){
      if (verbose > 0) cat("    * Estimating RMSE...","\n")
      #------------------------------------------------------------------------------------# 
      #     Because this doesn't use any other data (either from other stations or other   #
      # variables), we must estimate the error using a Jackknife approach.  Each iteration #
      # will delete a fraction of valid data and run the harmonic analysis, and the root   #
      # mean square error between the predicted value and the actual value will be the     #
      # error of this realisation.  We run this n.jack times so the error estimate is more #
      # robust.                                                                            #
      #------------------------------------------------------------------------------------# 
      nj       = 0
      err.jack = 0
      #----- List indices with available data. --------------------------------------------#
      idx.avail      = which(! miss.act)
      navail         = length(idx.avail)
      iterate        = TRUE
      while (iterate){
         nj = nj + 1

         #----- Copy the time series to a scratch vector. ---------------------------------#
         ndel           = floor(del.frac * navail) 
         del            = sample(x=idx.avail,size=ndel,replace=FALSE)
         jackknife      = x.act
         jackknife[del] = NA
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Find the gap-filled time series.  This call is slightly modified so we      #
         # preserve some of the properties of the original fitting.  We impose the number  #
         # of outer iterations to be exactly the same as the result with the full time     #
         # series.  Also, we must call the function forcing the error estimate to be       #
         # FALSE, otherwise we will enter in an infinite loop.  For the error we must also #
         # add the trend back because we compare the results with the original dataset.    #
         #---------------------------------------------------------------------------------#
         realisation    = harmfill( x              = jackknife
                                  , detrend.method = detrend.method
                                  , trend.back     = TRUE
                                  , signal.retain  = 1.00
                                  , minmod         = noutit
                                  , maxmod         = noutit
                                  , verbose        = 0
                                  , rmse           = FALSE
                                  , del.frac       = del.frac
                                  , n.jack         = n.jack   )
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Find the mean square error of this realisation.                             #
         #---------------------------------------------------------------------------------#
         err.jack.now  = sqrt( sum((x.act[del] - realisation$xfit[del])^2)
                             / (ndel - 2*noutit) )
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Find the new estimate of the error.                                        #
         #---------------------------------------------------------------------------------#
         err.jack.next = (err.jack * (nj-1) + err.jack.now) / nj
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Find the new estimate of the error.                                        #
         #---------------------------------------------------------------------------------#
         gain     = 2.0 * abs(err.jack.next - err.jack) / abs(err.jack.next +  err.jack)
         iterate  = gain > jack.toler
         err.jack = err.jack.next
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Print a  banner with the information on how the method is going.           #
         #---------------------------------------------------------------------------------#
         if (verbose > 0){ 
            cat ("      # Iteration: ",nj,"; RMSE = ",signif(err.jack,4)
                ,"; GAIN = ",sprintf("%4.4f",gain),"; CONVERGE = ",! iterate,"\n")
         }#end if
         #---------------------------------------------------------------------------------#
      }#end for
      #------------------------------------------------------------------------------------#
   }else{
      #------------------------------------------------------------------------------------#
      #    Skip the error estimate and leave it missing.                                   #
      #------------------------------------------------------------------------------------#
      err.jack = NA
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Make the vector with the error estimate.  The error is 0 for full points, and the  #
   # estimate for the gap-filled points.                                                   #
   #---------------------------------------------------------------------------------------#
   error = rep(0,times=nx.act)
   error[miss.act] = err.jack
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Make the list with the output then quit.                                          #
   #---------------------------------------------------------------------------------------#
   ans = list(xfit=xfit,error=error)
   return(ans)
   #---------------------------------------------------------------------------------------#

}#end function
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     This function detrends a time series by applying a linear regression.  In case there #
# are missing data, it fills in with a linear interpolation before detrending, to make     #
# sure the data is detrended properly (i.e., the gaps contains some linear interpolation). #
#------------------------------------------------------------------------------------------#
detrend2  <<- function(y){
   #----- Fill gap with linear interpolation. ---------------------------------------------#
   yfill       = gaplin(y)
   yprime      = detrend(yfill)
   y0          = yfill - yprime
   sel         = !is.finite(y)
   yprime[sel] = NA
   ans         = list(y0=y0,yprime=yprime,yfill=yfill)
   return(ans)
}#end function
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#      This function completes the gaps with a linear interpolation between the previous   #
# and the next available points.                                                           #
#------------------------------------------------------------------------------------------#
gaplin <<- function(x){
   gap       = ! is.finite(x)
   nx        = length(x)
   if (sum(gap) > 0){
      ind       = which(gap)
      ans       = sapply(X=ind,FUN=lin.filler,dat=x)
      xout      = x
      xout[gap] = ans
   }else{
      xout      = x
   }#end if
   return(xout)
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      Auxiliary function that will fill in the gap by looking at the previous and future  #
# available data by applying a linear interpolation between the previous and next avail-   #
# able data.                                                                               #
#------------------------------------------------------------------------------------------#
lin.filler <<- function(gapind,dat){
   if (is.finite(dat[gapind])){
      #----- No need to fill gaps, data is available. -------------------------------------#
      interp = dat[gapind]
   }else{
      #------------------------------------------------------------------------------------#
      #    Interpolate.                                                                    #
      #------------------------------------------------------------------------------------#
      ndat    = length(dat)

      #----- Split the time series into before and after gapind. --------------------------#
      prevdat = dat[1:gapind]
      prevok  = which(is.finite(prevdat))
      nextdat = dat[gapind:ndat]
      nextok  = (gapind - 1) + which(is.finite(nextdat))

      #----- Size of the split time series. -----------------------------------------------#
      nprev   = length(prevok)
      nnext   = length(nextok)

      #------------------------------------------------------------------------------------#
      #    Check whether we have points on both sides (interpolation) or only on one side  #
      # (extrapolation).                                                                   #
      #------------------------------------------------------------------------------------#
      if (nprev + nnext <= 1){
         stop("Time series has either 0 or 1 valid point! No interpolation possible!")
      }else if(nprev == 0){
         inda = nextok[1]
         indz = nextok[2]
      }else if(nnext == 0){
         inda = prevok[nprev-1]
         indz = prevok[nprev]
      }else{
         inda = prevok[nprev]
         indz = nextok[1]
      }#end if

      interp = dat[inda] + (gapind - inda) * (dat[indz]-dat[inda]) / (indz - inda)
   }#end if

   return(interp)
}#end function
#==========================================================================================#
#==========================================================================================#
