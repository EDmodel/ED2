#==========================================================================================#
#==========================================================================================#
#     This function defines the density profile corrected by the return distribution,      #
# following the MH69 correction and following NM01 implementation.  To convert the         #
# discrete returns into power return, we apply a waveform simulator following BH99.  In    #
# case we should use only the regular return count, set wfsim to FALSE.                    #
#                                                                                          #
# References:                                                                              #
#                                                                                          #
# Blair JB, Hofton, MA. Modeling laser altimeter return waveforms over complex vegetation  #
#    using high- resolution elevation data. Geophys. Res. Lett. 26(16):2509--2512,         #
#    doi:10.1029/1999GL010484 (BH99).                                                      #
#                                                                                          #
# MacArthur RA, Horn HS, 1969: Foliage profile by vertical measurements, Ecology 50(5),    #
#    802--804, doi:10.2307/1933693 (MH69).                                                 #
#                                                                                          #
# Ni-Meister W, Jupp D, and Dubayah R, 2001. Modeling lidar waveforms in heterogeneous and #
#    discrete canopies. IEEE T. Geosci. Remote Sens. 39(9):1943--1958,                     #
#    doi:10.1109/36.951085 (NM01).                                                         #
#                                                                                          #
# Popescu SC, Zhao K, Neuenschwander A, Lin C. 2011. Satellite lidar vs. small footprint   #
#    airborne lidar: Comparing the accuracy of aboveground biomass estimates and forest    #
#    structure metrics at footprint level. Remote Sens. Environ. 115(11):2786--2797,       #
#    doi:10.1016/j.rse.2011.01.026 (P11).   doi:10.1016/j.rse.2011.01.026 (P11).           #
#                                                                                          #
#------------------------------------------------------------------------------------------#
macarthur.horn <<- function( pt.cloud
                           , zh            = max(pt.cloud$z,na.rm=TRUE)
                           , zo            = 0.
                           , nz            = 512
                           , rvorg         = NA_real_
                           , sigma.z       = 5. * zh / (nz-1)
                           , Gmu           = 0.5
                           , tall.at.zh    = FALSE
                           , use.intensity = FALSE
                           , wfsim         = TRUE
                           , trunc0        = sqrt(.Machine$double.eps)
                           , zpad          = 2*ceiling(sigma.z*sqrt(log(1./trunc0)))
                           , show.profiles = FALSE
                           , ...
                           ){



   #---------------------------------------------------------------------------------------#
   #      Make sure the settings make sense.                                               #
   #---------------------------------------------------------------------------------------#
   #----- Check additional values. --------------------------------------------------------#
   if (! ( zh  %>%  0                                                              &&
           nz  %>%  0                                                              &&
           zh  %>%  zo                                                             &&
           zo  %>=% 0                                                              &&
           ( is.na(rvorg) || rvorg %>% 0 )                                         &&
           ( is.logical(tall.at.zh)      )                                         &&
           ( (sigma.z %>% 0) || (! wfsim) )                                        &&
           Gmu %>%0
         ) ){
      cat0("------------------------------------------------------------------")
      cat0(" MacArthur and Horn won't run due to problems with your settings:" )
      cat0(" "                                                                 )
      cat0(" ZO         = ",zo                                                 )
      cat0(" ZH         = ",zh                                                 )
      cat0(" NZ         = ",nz                                                 )
      cat0(" RVORG      = ",rvorg                                              )
      cat0(" TALL.AT.ZH = ",tall.at.zh                                         )
      cat0(" SIGMA.Z    = ",sigma.z                                            )
      cat0(" Gmu        = ",Gmu                                                )
      cat0(" "                                                                 )
      cat0(" Please check the following:"                                      )
      cat0(" ZO must be greater than or equal to 0 and less than ZH."          )
      cat0(" NZ must be positive."                                             )
      cat0(" RVORG must be positive or NA."                                    )
      cat0(" TALL.AT.ZH must be logical."                                      )
      cat0(" SIGMA.Z must be positive when simulating waveform."               )
      cat0(" Gmu must be positive."                                            )
      cat0(" "                                                                 )
      cat0("------------------------------------------------------------------")
      stop("Invalid height settings.")
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Make sure that the point cloud has the minimum data needed.                      #
   #---------------------------------------------------------------------------------------#
   if (is.data.frame(pt.cloud) || is.list(pt.cloud) || is.matrix(pt.cloud)){

      #----- Coerce matrix to a data frame. -----------------------------------------------#
      if (is.matrix(pt.cloud)) pt.cloud = as.data.frame(pt.cloud)
      #------------------------------------------------------------------------------------#


      #----- Check that names are fine or that the data frame has three variables. --------#
      names.ok  = all(c("z","intensity","pt.class") %in% names(pt.cloud))
      three.ok  = length(names(pt.cloud)) >= 3
      #------------------------------------------------------------------------------------#



      #------ Make sure that we can derive properties. ------------------------------------#
      if (names.ok){
         #----- Names were given.  Check that the length is valid for both variables. -----#
         len.z   = length(pt.cloud$z)
         len.int = length(pt.cloud$intensity)
         len.cls = length(pt.cloud$pt.class)
         if (len.z != len.int || len.z != len.cls){
            stop(" 'z', 'intensity', and 'pt.class' elements must have the same length!")
         }#end if
         #---------------------------------------------------------------------------------#
      }else if (three.ok){
         #---------------------------------------------------------------------------------#
         #     Names were not given but there are two elements. Assume that they are       #
         # height and intensity, respectively.                                             #
         #---------------------------------------------------------------------------------#
         names(pt.cloud)[c(1,2,3)] = c("z","intensity","pt.class")
         warning(" Assumed first 3 columns are z, intensity and pt.class, respectively...")
         len.z   = length(pt.cloud$z)
         len.int = length(pt.cloud$intensity)
         len.cls = length(pt.cloud$pt.class)
         if (len.z != len.int || len.z != len.cls){
            stop(" 'z', 'intensity', and 'pt.class' elements must have the same length!")
         }#end if
         #---------------------------------------------------------------------------------#
      }else{
         #----- Cannot use this point cloud. ----------------------------------------------#
         stop("Missing elements 'z', 'intensity', and 'pt.class' from pt.cloud.")
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#
   }else{
      #----- Invalid object for a point cloud. --------------------------------------------#
      stop("Object pt.cloud must be a data frame, a list, or a matrix.")
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Find the height bin levels and width. -------------------------------------------#
   zmid     = seq(from=0,to=zh,length.out=nz)
   deltaz   = mean(diff(zmid))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     For the waveform simulation, we must some pad to ensure that the beam source is   #
   # entirely within the height domain, and that we have layers that go below ground to    #
   # fully characterize the energy.  These layers will be removed for the final output.    #
   #---------------------------------------------------------------------------------------#
   zpad       = deltaz + seq(from=0,to=zpad,by=deltaz)
   npad       = length(zpad)
   zlwr       = 0. - rev(zpad)
   zupr       = zh + zpad
   zconv      = c(zlwr,zmid,zupr)
   nzconv     = length(zconv)
   zbreaks    = c(zconv-0.5*deltaz,max(zconv)+0.5*deltaz)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     In case tall.at.zh is set to TRUE, change height.                                 #
   #---------------------------------------------------------------------------------------#
   if (tall.at.zh){
      zprune     = zo + (1.-sqrt(.Machine$double.eps))*(zh-zo)
      pt.cloud$z = pmin(pt.cloud$z,zprune)
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Assume that everything with height zero is ground. ------------------------------#
   pt.cloud$pt.class[pt.cloud$z %<=% 0] = 2
   #---------------------------------------------------------------------------------------#


   #----- Keep only the points that are within bounds. ------------------------------------#
   keep     = pt.cloud$z %>=% 0 & pt.cloud$z %<% zh & pt.cloud$pt.class %in% c(0,1,2,3,4,5)
   pt.cloud = pt.cloud[keep,,drop=FALSE]
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      In case intensity is not to be used, set intensities to 1.                       #
   #---------------------------------------------------------------------------------------#
   if (! use.intensity){
      pt.cloud$intensity = 0. * pt.cloud$intensity + 1
   }#end if (use.intensity)
   #---------------------------------------------------------------------------------------#



   #----- Generate flags for ground and vegetation returns. -------------------------------#
   veg.cloud = pt.cloud[pt.cloud$pt.class %in% c(0,1,3,4,5),,drop=FALSE]
   gnd.cloud = pt.cloud[pt.cloud$pt.class %in% c(2)        ,,drop=FALSE]
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #       In case no ground return exists, add one point with minimal intensity and       #
   # near the ground.                                                                      #
   #---------------------------------------------------------------------------------------#
   if (nrow(gnd.cloud) == 0){
      idx                      = which.min(veg.cloud$z)
      gnd.cloud                = veg.cloud[idx,]
      gnd.cloud$x              = mean(veg.cloud$x)
      gnd.cloud$y              = mean(veg.cloud$y)
      gnd.cloud$z              = 0.01
      gnd.cloud$intensity      = min(veg.cloud$intensity)
      gnd.cloud$pt.class       = 2
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Decide how to determine the vertical profile.                                      #
   #---------------------------------------------------------------------------------------#
   if (wfsim){
      #------------------------------------------------------------------------------------#
      #    Find the horizontal contribution from each layer, assuming that all returns     #
      # contribute proportionally to their intensity (i.e. assume footprint to be square   #
      # as opposed to Gaussian).                                                           #
      #------------------------------------------------------------------------------------#
      #----- Vegetation returns. ----------------------------------------------------------#
      veg.zcut    = as.integer(cut(x=veg.cloud$z,breaks=zbreaks,right=FALSE))
      veg.wh      = rep(0,times=nzconv)
      aux         = tapply(X=veg.cloud$intensity,INDEX=veg.zcut,FUN=sum)
      idx         = as.numeric(names(aux))
      veg.wh[idx] = aux
      #----- Ground returns. --------------------------------------------------------------#
      gnd.zcut    = as.integer(cut(x=gnd.cloud$z,breaks=zbreaks,right=FALSE))
      gnd.wh      = rep(0,times=nzconv)
      aux         = tapply(X=gnd.cloud$intensity,INDEX=gnd.zcut,FUN=sum)
      idx         = as.numeric(names(aux))
      gnd.wh[idx] = aux
      #----- Total returns (for debugging only). ------------------------------------------#
      tot.wh      = veg.wh + gnd.wh
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Find the vertical distribution of emitted pulses (same for vegetation and      #
      # ground returns.  We use for loop to reduce memory burden for really large chunks.  #
      #------------------------------------------------------------------------------------#
      Rvlyr = mapply( FUN = find.waveform
                    , z   = zconv
                    , MoreArgs = list( zi      = zconv
                                     , wh      = veg.wh
                                     , sigma.z = sigma.z
                                     )#end list
                    )#end mapply
      Rglyr = mapply( FUN = find.waveform
                    , z   = zconv
                    , MoreArgs = list( zi      = zconv
                                     , wh      = gnd.wh
                                     , sigma.z = sigma.z
                                     )#end list
                    )#end mapply


      #----- Obtain the cumulative energy for the MacArthur-Horn correction. --------------#
      Rv    = rev(cumsum(rev(Rvlyr)))
      Rv0   = sum(Rvlyr)
      Rg    = sum(Rglyr)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Plot curves to debug.                                                          #
      #------------------------------------------------------------------------------------#
      if (show.profiles){
         graphics.off()
         plot.new()
         plot.window(xlim=c(0,1),ylim=pretty.xylim(zmid,fracexp=c(-0.1,0.1)))
         axis(side=1,las=1)
         axis(side=2,las=2)
         abline(h=0,col="grey50",lty="dotted",lwd=1.0)
         lines(x=tot.wh/max(tot.wh),y=zconv,type="l",col="grey30" ,lwd=0.5)
         lines(x=Rglyr/max(Rglyr)  ,y=zconv,type="l",col="#CB003D",lwd=2.0)
         lines(x=Rvlyr/max(Rvlyr)  ,y=zconv,type="l",col="#107C92",lwd=2.0)
         legend( x      = "topright"
               , inset  = c(0.01,0.20)
               , legend = c( ifelse(use.intensity,"Intensity sum","Return count")
                           , "Waveform (ground)"
                           , "Waveform (vegetation)"
                           )#end legend
               , col    = c("grey30"
                           ,"#CB003D"
                           ,"#107C92")
               , lwd    = c(0.5,1.0,2.0,2.0)
               , cex    = 0.7
               , bty    = "n"
               )#end legend
         box()
         cat0(" -> Click on the plot for the next plot.")
         locator(n=1)
      }#end (show.profiles)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Remove the pads, and keep only the profile in the layers of interest.          #
      #------------------------------------------------------------------------------------#
      izero  = which.closest(0.,zconv)
      ikeep  = seq(from=izero,by=1,length.out=nz)
      veg.wh = veg.wh[ikeep]
      gnd.wh = gnd.wh[ikeep]
      Rvlyr  = Rvlyr [ikeep]
      Rglyr  = Rglyr [ikeep]
      Rv     = Rv    [ikeep]
      #------------------------------------------------------------------------------------#
   }else{



      #----- Find the total energy returned from each layer. ------------------------------#
      zcut           = cut(x=veg.cloud$z,breaks=zbreaks,right=FALSE)
      zcut           = match(zcut,levels(zcut))
      Rvlyr          = rep(0,times=nconv)
      aux            = tapply(X=veg.cloud$intensity,INDEX=zcut,FUN=sum)
      idx            = as.numeric(names(aux))
      Rvlyr[idx]     = aux
      Rvlyr[izero-1] = sum(gnd.cloud$intensity)
      #------------------------------------------------------------------------------------#


      #----- Find the total intensity returned and flagged as vegetation and as ground. ---#
      Rv  = rev(cumsum(rev(Rvlyr)))
      Rv0 = sum(veg.cloud$intensity)
      Rg  = sum(gnd.cloud$intensity)
      #------------------------------------------------------------------------------------#



      #----- Remove the pads, and keep only the profile in the layers of interest. --------#
      izero = which.closest(0.,zconv)
      ikeep = seq(from=izero,by=1,length.out=nz)
      Rvlyr = Rvlyr[ikeep]
      Rv    = Rv   [ikeep]
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Remove the pads, and keep only the profile in the layers of interest.          #
      #------------------------------------------------------------------------------------#
      if (show.profiles){
         graphics.off()
         plot.new()
         plot.window(xlim=c(0,1),ylim=range(zmid))
         axis(side=1,las=1)
         axis(side=2,las=2)
         abline(h=0,col="grey50",lty="dotted",lwd=1.0)
         lines(x=Rvlyr/max(Rvlyr),y=zmid,type="l",col="#8C510A",lwd=2)
         lines(x=Rv/max(Rv),y=zmid,type="l",col="#7EC4BC",lwd=2)
         legend( x      = "topright"
               , inset  = 0.01
               , legend = c( "Waveform (layer)"
                           , "Cumulative energy"
                           )#end legend
               , col    = c("#8C510A","#7EC4BC")
               , lwd    = c(2.0,2.0)
               , cex    = 0.7
               , bty    = "n"
               )#end legend
         box()
         cat0(" -> Click on the plot for the next plot")
         locator(n=1)
      }#end (show.profiles)
      #------------------------------------------------------------------------------------#
   }#end if(wfsim)
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Assume fraction between backscatterings in case none is provided.                 #
   #---------------------------------------------------------------------------------------#
   if (is.na(rvorg)){
     #kuse = 0.825 # Mean value by Antonarakis et al. (2014).
     kuse = 1.03   # Mean value by Tang and Dubayah (2017).
   }else{
     kuse = rvorg
   }#end if
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Find the gap fraction.                                                            #
   #---------------------------------------------------------------------------------------#
   gap.bot     = 1. - Rv / ( Rv0 + kuse * Rg)
   gap.top     = c(gap.bot[-1],1)
   gap.mid     = sqrt(gap.bot*gap.top)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Find the unscaled values for LAD.  The log ratio is from bottom to top because   #
   # we use top-down integration.                                                          #
   #---------------------------------------------------------------------------------------#
   lad    = log(gap.top / gap.bot) / (Gmu * deltaz)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     In case rvorg hasn't been provided, scale lad to unity.                           #
   #---------------------------------------------------------------------------------------#
   if (! is.na(rvorg)){
      LAI    = sum(lad * deltaz)
      lad    = lad / LAI
   }#end if
   if (! all(is.finite(lad))) browser()
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Define a waveform-related PDF.                                                    #
   #---------------------------------------------------------------------------------------#
   if (wfsim){
      #----- Use the simulated waveform to prescribe the pdf. -----------------------------#
      Rv.underground = max(0.,Rv0-sum(Rvlyr))
      x.Rvlyr        = Rvlyr
      x.Rvlyr[1]     = Rvlyr + Rv.underground
      pdfsum         = sum(x.Rvlyr*deltaz)
      lpdf           = data.frame(x=zmid,y=Rvlyr/pdfsum)
      pdfsum         = sum(lpdf$y * deltaz)
      lcdf           = rev(cumsum(rev(lpdf$y*deltaz))) / pdfsum
      #------------------------------------------------------------------------------------#

   }else{
      #----- Create a pseudo-point cloud with the MH-corrected distribution. --------------#
      nzmah = ceiling(3.*Rv0/min(veg.cloud$intensity[veg.cloud$intensity %>% 0]))
      zmah  = jitter(x= sample(x=zmid,size=nzmah,replace=TRUE,prob=lad),amount=0.5*deltaz)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Find the waveform Smooth the curve by generating a density function.  This     #
      # density function becomes the output.                                               #
      #------------------------------------------------------------------------------------#
      lpdf   = density.safe(zmah,from=0,to=zh,n=nz,...)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Find the cumulative distribution function.  Make sure the CDF goes from 0      #
      # (top) to 1 (bottom).                                                               #
      #------------------------------------------------------------------------------------#
      pdfsum = sum(lpdf$y * deltaz)
      lcdf   = rev(cumsum(rev(lpdf$y*deltaz))) / pdfsum
      #------------------------------------------------------------------------------------#
   }#end if (wfsim)
   #---------------------------------------------------------------------------------------#


   #----- Select the raw profile. ---------------------------------------------------------#
   if (wfsim){
      lraw = veg.wh + gnd.wh
   }else{
      lraw = Rvlyr
   }#end if 
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Remove points below zo.                                                           #
   #---------------------------------------------------------------------------------------#
   keep    = zmid >= zo
   zmid    = zmid   [keep]
   lad     = lad    [keep]
   deltaz  = deltaz + 0*zmid
   Rvlyr   = Rvlyr  [keep]
   Rv      = Rv     [keep]
   gap.top = gap.top[keep]
   gap.bot = gap.bot[keep]
   gap.mid = gap.mid[keep]
   lpdf    = data.frame(x=lpdf$x[keep],y=lpdf$y[keep])
   lcdf    = lcdf   [keep] / max(lcdf[keep])
   lraw    = lraw   [keep]
   #---------------------------------------------------------------------------------------#



   #----- Build data frame with full structure. -------------------------------------------#
   ans    = try( data.frame( x     = lpdf$x
                           , y     = lpdf$y
                           , z     = zmid
                           , dz    = deltaz
                           , pdf   = lpdf$y
                           , cdf   = lcdf
                           , lad   = lad
                           , raw   = lraw
                           , gap   = gap.mid
                           )#end data.frame
               , silent = TRUE
               )#end try
   if ("try-error" %in% is(ans)) browser()
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Plot curves to debug.                                                             #
   #---------------------------------------------------------------------------------------#
   if (show.profiles){
      lai = rev(cumsum(rev(lad*deltaz)))
      graphics.off()
      plot.new()
      plot.window(xlim=c(0,1),ylim=range(zmid))
      axis(side=1,las=1)
      axis(side=2,las=2)
      abline(h=0,col="grey50",lty="dotted",lwd=1.0)
      lines(x=Rvlyr/max(Rvlyr),y=zmid,type="l",col="#811F9E",lwd=2)
      lines(x=lad/max(lad)    ,y=zmid,type="l",col="#1BA2F7",lwd=2)
      lines(x=lai/max(lai)    ,y=zmid,type="l",col="#107C92",lwd=2)
      legend( x      = "topright"
            , inset  = 0.01
            , legend = c( "Waveform (vegetation)"
                        , "Leaf area density"
                        , "Leaf area index"
                        )#end legend
            , col    = c("#811F9E","#1BA2F7","#107C92")
            , lwd    = c(2.0,2.0,2.0)
            , cex    = 0.7
            , bty    = "n"
            )#end legend
      box()
      cat0(" -> Click on the plot for the next plot.")
      locator(n=1)
   }#end (show.profiles)
   #---------------------------------------------------------------------------------------#


   #----- Return the answer. --------------------------------------------------------------#
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function macarthur.horn
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#    This function computes the waveform function for any given height.                    #
#------------------------------------------------------------------------------------------#
find.waveform <<- function(wh,z,zi,sigma.z){
   #----- Find the pulse shape along the beam path. ---------------------------------------#
   wv  = exp(-(z-zi)^2/(2.*sigma.z^2)) / (sigma.z*sqrt(2*pi))
   #----- Run the convolution and obtain the waveform for each layer. ---------------------#
   ans = convolve(x=wh,y=wv,type="filter") 
   #----- Make sure the waveform is not negative. -----------------------------------------#
   ans = pmax(0,ans) + 0. * ans
   #---------------------------------------------------------------------------------------#



   #----- Return the answer. --------------------------------------------------------------#
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function find.waveform
#==========================================================================================#
#==========================================================================================#
