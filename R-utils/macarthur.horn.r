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
                           , sigma.t       = 1e-8
                           , Gmu           = 0.5
                           , tall.at.zh    = FALSE
                           , use.intensity = FALSE
                           , wfsim         = FALSE
                           , zair          = 850.
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
           ( (sigma.t %>% 0) || (! wfsim) )                                        &&
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
      cat0(" SIGMA.T    = ",sigma.t                                            )
      cat0(" Gmu        = ",Gmu                                                )
      cat0(" "                                                                 )
      cat0(" Please check the following:"                                      )
      cat0(" ZO must be greater than or equal to 0 and less than ZH."          )
      cat0(" NZ must be positive."                                             )
      cat0(" RVORG must be positive or NA."                                    )
      cat0(" TALL.AT.ZH must be logical."                                      )
      cat0(" SIGMA.T must be positive when simulating waveform."               )
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
      stop("Object pt.cloud must be a data frame, a list, or a matrix...")
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Find the height breaks in case none has been given. -----------------------------#
   zmid    = seq(from=0,to=zh,length.out=nz)
   deltaz  = mean(diff(zmid))
   zbreaks = seq(from=0-0.5*deltaz,to=zh+0.5*deltaz,length.out=nz+1)
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
   #    Decide how to determine the vertical profile.                                      #
   #---------------------------------------------------------------------------------------#
   if (wfsim){
      #------------------------------------------------------------------------------------#
      #    Find the horizontal contribution from each layer, assuming that all returns     #
      # contribute proportionally to their intensity (i.e. assume footprint to be square   #
      # as opposed to Gaussian).                                                           #
      #------------------------------------------------------------------------------------#
      zcut    = as.integer(cut(x=veg.cloud$z,breaks=zbreaks,right=FALSE))
      wh      = rep(0,times=nz)
      aux     = tapply(X=veg.cloud$intensity,INDEX=zcut,FUN=sum)
      idx     = as.numeric(names(aux))
      wh[idx] = aux
      wh[1]   = sum(gnd.cloud$intensity)
      #------------------------------------------------------------------------------------#



      #----- Find the vertical distribution of emitted pulses. ----------------------------#
      tt      = 2 * (zh-zmid) / clight
      wv      = exp(-2 * tt^2 / sigma.t^2)
      #------------------------------------------------------------------------------------#


      #----- Run the convolution. ---------------------------------------------------------#
      Rvlyr = rev(convolve(rev(wh),rev(wv)))
      Rvlyr = pmax(0,Rvlyr) + 0. * Rvlyr
      Rv    = rev(cumsum(rev(Rvlyr)))
      Rv0   = sum(Rvlyr)
      Rg    = Rvlyr[1]
      #------------------------------------------------------------------------------------#
   }else{


      #------------------------------------------------------------------------------------#
      #       In case no ground return exists, add one point with minimal intensity and    #
      # near the ground.                                                                   #
      #------------------------------------------------------------------------------------#
      if (nrow(gnd.cloud) == 0){
         idx                      = which.min(veg.cloud$z)
         gnd.cloud                = veg.cloud[idx,]
         gnd.cloud$x              = mean(veg.cloud$x)
         gnd.cloud$y              = mean(veg.cloud$y)
         gnd.cloud$z              = 0.01
         gnd.cloud$intensity      = min(veg.cloud$intensity)
         gnd.cloud$retn.number    = min(veg.cloud$retn.number)
         gnd.cloud$number.retn.gp = commonest(veg.cloud$number.retn.gp)
         gnd.cloud$pt.class       = 2
         gnd.cloud$gpstime        = max(veg.cloud$gpstime)
      }#end if
      #------------------------------------------------------------------------------------#



      #----- Find the total energy returned from each layer. ------------------------------#
      zcut       = cut(x=veg.cloud$z,breaks=zbreaks,right=FALSE)
      zcut       = match(zcut,levels(zcut))
      Rvlyr      = rep(0,times=nz)
      aux        = tapply(X=veg.cloud$intensity,INDEX=zcut,FUN=sum)
      idx        = as.numeric(names(aux))
      Rvlyr[idx] = aux
      Rvlyr[1]   = sum(gnd.cloud$intensity)
      #------------------------------------------------------------------------------------#


      #----- Find the total intensity returned and flagged as vegetation and as ground. ---#
      Rv  = rev(cumsum(rev(Rvlyr)))
      Rv0 = sum(veg.cloud$intensity)
      Rg  = sum(gnd.cloud$intensity)
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
   #     Create a pseudo point cloud using the correction by MacArthur and Horn.           #
   #---------------------------------------------------------------------------------------#
   if (! all(is.finite(lad))) browser()
   nzmah = ceiling(3.*Rv0/min(veg.cloud$intensity[veg.cloud$intensity %>% 0]))
   zmah  = jitter(x= sample(x=zmid,size=nzmah,replace=TRUE,prob=lad),amount=0.5*deltaz)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Smooth the curve by generating a density function.  This density function         #
   # becomes the output.                                                                   #
   #---------------------------------------------------------------------------------------#
   lpdf   = density.safe(zmah,from=0,to=zh,n=nz)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     In case rvorg hasn't been provided, scale lad to unity.                           #
   #---------------------------------------------------------------------------------------#
   if (! is.na(rvorg)){
      LAI    = sum(lad * deltaz)
      lad    = lad / LAI
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the total LAI given rvorg.  If rvorg hasn't been provided, then the          #
   # total LAI will be used as a normalisation factor.                                     #
   #---------------------------------------------------------------------------------------#
   pdfsum = sum(lpdf$y * deltaz)
   lcdf   = rev(cumsum(rev(lpdf$y*deltaz))) / pdfsum
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
   if (wfsim){
      wh = wh[keep]
      wv = wv[keep]
   }#end if (wfsim)
   #---------------------------------------------------------------------------------------#



   #----- Build data frame with full structure. -------------------------------------------#
   ans    = try( data.frame( x     = lpdf$x
                           , y     = lpdf$y
                           , z     = zmid
                           , dz    = deltaz
                           , pdf   = lpdf$y
                           , cdf   = lcdf
                           , lad   = lad
                           , Rvlyr = Rvlyr
                           , Rv    = Rv
                           , gap   = gap.mid
                           )#end data.frame
               , silent = TRUE
               )#end try
   if ("try-error" %in% is(ans)) browser()
   #---------------------------------------------------------------------------------------#


   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function macarthur.horn
#==========================================================================================#
#==========================================================================================#
