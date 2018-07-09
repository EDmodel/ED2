#==========================================================================================#
#==========================================================================================#
#     This function defines the density profile corrected by the pulse intensity, follow-  #
# ing the MacArthur and Horn (1969) correction.                                            #
#                                                                                          #
# Reference: MacArthur, R. A. and H. S. Horn, 1969: Foliage profile by vertical            #
#               measurements, Ecology, 50(5), 802--804.                                    #
#                                                                                          #
#------------------------------------------------------------------------------------------#
macarthur.horn <<- function( pt.cloud
                           , zh            = max(pt.cloud$z,na.rm=TRUE)
                           , zo            = 0.
                           , nz            = 512
                           , rvorg         = NA
                           , Gmu           = 0.5
                           , tall.at.zh    = FALSE
                           , use.intensity = FALSE
                           ){

   #---------------------------------------------------------------------------------------#
   #      Make sure the settings make sense.                                               #
   #---------------------------------------------------------------------------------------#
   #----- Check additional values. --------------------------------------------------------#
   if (! ( zh  %>%  0                      &&
           nz  %>%  0                      &&
           zh  %>%  zo                     &&
           zo  %>=% 0                      &&
           ( is.na(rvorg) || rvorg %>% 0 ) &&
           ( is.logical(tall.at.zh)      ) && 
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
      cat0(" Gmu        = ",Gmu                                                )
      cat0(" "                                                                 )
      cat0(" Please check the following:"                                      )
      cat0(" ZO must be greater than or equal to 0 and less than ZH."          )
      cat0(" NZ must be positive."                                             )
      cat0(" RVORG must be positive or NA."                                    )
      cat0(" TALL.AT.ZH must be logical."                                      )
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
   dzbar   = mean(diff(zmid))
   zbreaks = seq(from=0-0.5*dzbar,to=zh+0.5*dzbar,length.out=nz+1)
   deltaz  = diff(zbreaks)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     In case tall.at.zh is set to TRUE, change height.                                 #
   #---------------------------------------------------------------------------------------#
   if (tall.at.zh){
      zprune     = zo + (1.-sqrt(.Machine$double.eps))*(zh-zo)
      pt.cloud$z = pmin(pt.cloud$z,zprune)
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Keep only the points that are within bounds. ------------------------------------#
   keep     = pt.cloud$z %>=% 0 & pt.cloud$z %<% zh & pt.cloud$pt.class %in% c(0,1,2,3,4,5)
   pt.cloud = pt.cloud[keep,,drop=FALSE]
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      In case intensity is not to be used, keep only first returns, and set            #
   # intensities to 1.                                                                     #
   #---------------------------------------------------------------------------------------#
   if (! use.intensity){
      first              = unlist( tapply( X     = pt.cloud$retn.number
                                         , INDEX = pt.cloud$pulse.number
                                         , FUN   = function(x) x %==% min(x)
                                         )#end tapply
                                 )#end unlist
      pt.cloud           = pt.cloud[first,,drop=FALSE]
      pt.cloud$intensity = 0. * pt.cloud$intensity + 1
   }#end if (use.intensity)
   #---------------------------------------------------------------------------------------#


   #----- Generate flags for ground and vegetation returns. -------------------------------#
   veg.cloud = pt.cloud[pt.cloud$pt.class %in% c(0,1,3,4,5),]
   gnd.cloud = pt.cloud[pt.cloud$pt.class %in% c(2)        ,]
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #       In case no ground return exists, add one point with minimal intensity and near  #
   # the ground.                                                                           #
   #---------------------------------------------------------------------------------------#
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
   #---------------------------------------------------------------------------------------#



   #----- Find the total energy returned from each layer. ---------------------------------#
   zcut       = cut(x=veg.cloud$z,breaks=zbreaks,right=FALSE)
   zcut       = match(zcut,levels(zcut))
   Rvlyr      = rep(0,times=nz)
   aux        = tapply(X=veg.cloud$intensity,INDEX=zcut,FUN=sum)
   idx        = as.numeric(names(aux))
   Rvlyr[idx] = aux
   #---------------------------------------------------------------------------------------#


   #----- Find the total intensity returned and flagged as vegetation and as ground. ------#
   Rv  = rev(cumsum(rev(Rvlyr)))
   Rv0 = sum(veg.cloud$intensity)
   Rg  = sum(gnd.cloud$intensity)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Assume fraction between backscatterings in case none is provided.                 #
   #---------------------------------------------------------------------------------------#
   if (is.na(rvorg)){
     #kuse = 0.825 # Mean value by Antonarakis et al. (2014).
     kuse = 1.03  # Mean value by Tang and Dubayah (2017).
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
   dzbar = mean(diff(zmid))
   if (! all(is.finite(lad))) browser()
   nzmah = ceiling(3.*Rv0/min(veg.cloud$intensity[veg.cloud$intensity %>% 0]))
   zmah  = jitter(x= sample(x=zmid,size=nzmah,replace=TRUE,prob=lad),amount=0.5*dzbar)
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
   deltaz  = deltaz [keep]
   Rvlyr   = Rvlyr  [keep]
   Rv      = Rv     [keep]
   gap.top = gap.top[keep]
   gap.bot = gap.bot[keep]
   gap.mid = gap.mid[keep]
   lpdf    = data.frame(x=lpdf$x[keep],y=lpdf$y[keep])
   lcdf    = lcdf   [keep] / max(lcdf[keep])
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
