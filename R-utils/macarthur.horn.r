#==========================================================================================#
#==========================================================================================#
#     This function defines the density profile corrected by the pulse intensity, follow-  #
# ing the MacArthur and Horn (1969) correction.                                            #
#                                                                                          #
# Reference: MacArthur, R. A. and H. S. Horn, 1969: Foliage profile by vertical            #
#               measurements, Ecology, 50(5), 802--804.                                    #
#------------------------------------------------------------------------------------------#
macarthur.horn <<- function(pt.cloud,zl,zh,zo=zl,nz,rvorg=NA,Gmu=0.5){

   #----- Make sure the settings make sense. ----------------------------------------------#
   if (! ( zl  %>=% 0  &&
           zh  %>%  zl &&
           nz  %>%  0  &&
           zh  %>%  zo &&
           zo  %>=% zl &&
           ( is.na(rvorg) || rvorg %>% 0 ) && 
           Gmu %>%0
         ) ){
      cat("------------------------------------------------------------------","\n",sep="")
      cat(" MacArthur and Horn won't run due to problems with your settings:" ,"\n",sep="")
      cat(" "                                                                 ,"\n",sep="")
      cat(" ZL    = ",zl                                                      ,"\n",sep="")
      cat(" ZO    = ",zo                                                      ,"\n",sep="")
      cat(" ZH    = ",zh                                                      ,"\n",sep="")
      cat(" NZ    = ",nz                                                      ,"\n",sep="")
      cat(" RVORG = ",rvorg                                                   ,"\n",sep="")
      cat(" Gmu   = ",Gmu                                                     ,"\n",sep="")
      cat(" "                                                                 ,"\n",sep="")
      cat(" Please check the following:"                                      ,"\n",sep="")
      cat(" ZL must be greater than or equal to zero and less than ZH."       ,"\n",sep="")
      cat(" ZO must be greater than or equal to ZL and less than ZH."         ,"\n",sep="")
      cat(" NZ must be positive."                                             ,"\n",sep="")
      cat(" RVORG must be positive or NA."                                    ,"\n",sep="")
      cat(" Gmu must be positive."                                            ,"\n",sep="")
      cat(" "                                                                 ,"\n",sep="")
      cat("------------------------------------------------------------------","\n",sep="")
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
   zbreaks = sort(unique(c(0,seq(from=zl,to=zh,length.out=nz))))
   zmid    = mid.points(zbreaks)
   deltaz  = diff(zbreaks)
   #---------------------------------------------------------------------------------------#


   #----- Keep only the points that are within bounds. ------------------------------------#
   keep     = pt.cloud$z %>=% 0 & pt.cloud$z %<% zh & pt.cloud$pt.class %in% c(0,1,2,3,4,5)
   pt.cloud = pt.cloud[keep,]
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
      gnd.cloud$intensity      = ceiling(min(veg.cloud$intensity)/2)
      gnd.cloud$retn.number    = max(veg.cloud$retn.number)
      gnd.cloud$number.retn.gp = commonest(veg.cloud$number.retn.gp)
      gnd.cloud$pt.class       = 2
      gnd.cloud$gpstime        = max(veg.cloud$gpstime)
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Find the total energy returned from each layer. ---------------------------------#
   zcut    = cut(x=veg.cloud$z,breaks=zbreaks,right=FALSE)
   zcut    = match(zcut,levels(zcut))
   Rv      = rep(0,times=nz)
   aux     = tapply(X=veg.cloud$intensity,INDEX=zcut,FUN=sum)
   idx     = as.numeric(names(aux))
   Rv[idx] = aux
   #---------------------------------------------------------------------------------------#


   #----- Find the total intensity returned and flagged as vegetation and as ground. ------#
   Rv  = rev(cumsum(rev(Rv)))
   Rv0 = sum(veg.cloud$intensity)
   Rg  = sum(gnd.cloud$intensity)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Assume fraction between backscatterings in case none is provided.                 #
   #---------------------------------------------------------------------------------------#
   if (is.na(rvorg)){
     kuse = 1
   }else{
     kuse = rvorg
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the gap fraction.                                                            #
   #---------------------------------------------------------------------------------------#
   extinct.bot = 1. - Rv / ( Rv0 * (1 + kuse * Rg / Rv0 ) )
   extinct.top = c(extinct.bot[-1],1)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Find the unscaled values for LAD.  The log ratio is from bottom to top because   #
   # we use top-down integration.                                                          #
   #---------------------------------------------------------------------------------------#
   lad    = - 2 * log(extinct.bot / extinct.top) / deltaz
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Remove points below zo.                                                           #
   #---------------------------------------------------------------------------------------#
   keep     = zmid >= zo
   zmid     = zmid  [keep]
   lad      = lad   [keep]
   deltaz   = deltaz[keep]
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Smooth the curve by generating a density function.  This density function becomes #
   # the output.                                                                           #
   #---------------------------------------------------------------------------------------#
   dzbar = mean(diff(zmid))
   if (! all(is.finite(lad))) browser()
   zmah  = jitter( x      = sample(x=zmid,size=Rv0,replace=TRUE,prob=lad)
                 , amount = 0.5*dzbar
                 )#end jitter
   ans   = density.safe(zmah,from=zl,to=zh,n=nz)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     In case rvorg hasn't been provided, scale lad to unity.                           #
   #---------------------------------------------------------------------------------------#
   if (! is.na(rvorg)){
      LAI   = sum(lad * dzbar)
      ans$y = ans$y * LAI
   }#end if
   #---------------------------------------------------------------------------------------#

   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function macarthur.horn
#==========================================================================================#
#==========================================================================================#
