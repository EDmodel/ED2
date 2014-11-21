#==========================================================================================#
#==========================================================================================#
#     This function defines the density profile corrected by the pulse intensity, follow-  #
# ing the MacArthur and Horn (1969) correction.                                            #
#                                                                                          #
# Reference: MacArthur, R. A. and H. S. Horn, 1969: Foliage profile by vertical            #
#               measurements, Ecology, 50(5), 802--804.                                    #
#------------------------------------------------------------------------------------------#
macarthur.horn <<- function(pt.cloud,zl,zh,zo=zl,nz){

   #----- Make sure the settings make sense. ----------------------------------------------#
   if (! (zl %>% 0 && zh %>% zl && zh %>% zo && zo %>=% zl && nz %>% 0) ){
      cat("------------------------------------------------------------------","\n",sep="")
      cat(" MacArthur and Horn won't run due to problems with your settings:" ,"\n",sep="")
      cat(" "                                                                 ,"\n",sep="")
      cat(" ZL = ",zl                                                         ,"\n",sep="")
      cat(" ZO = ",zl                                                         ,"\n",sep="")
      cat(" ZH = ",zh                                                         ,"\n",sep="")
      cat(" NZ = ",nz                                                         ,"\n",sep="")
      cat(" "                                                                 ,"\n",sep="")
      cat(" Please check the following:"                                      ,"\n",sep="")
      cat(" ZL must be greater than zero and less than ZH."                   ,"\n",sep="")
      cat(" ZO must be greater than or equal to ZL and less than ZH."         ,"\n",sep="")
      cat(" NZ must be positive."                                             ,"\n",sep="")
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


      #----- Check that names are fine or that the data frame has two variables. ----------#
      names.ok  = all(c("z","intensity") %in% names(pt.cloud))
      two.ok    = length(names(pt.cloud)) >= 2
      #------------------------------------------------------------------------------------#



      #------ Make sure that we can derive properties. ------------------------------------#
      if (names.ok){
         #----- Names were given.  Check that the length is valid for both variables. -----#
         if (length(pt.cloud$z) != length(pt.cloud$intensity)){
            stop(" 'z' and 'intensity' elements must have the same length!")
         }#end if
         #---------------------------------------------------------------------------------#
      }else if (two.ok){
         #---------------------------------------------------------------------------------#
         #     Names were not given but there are two elements. Assume that they are       #
         # height and intensity, respectively.                                             #
         #---------------------------------------------------------------------------------#
         names(pt.cloud)[c(1,2)] = c("z","intensity")
         warning(" Assuming that 1st element is height and the 2nd is intensity...")
         if (length(pt.cloud$z) != length(pt.cloud$intensity)){
            stop(" 'z' and 'intensity' elements must have the same length!")
         }#end if
         #---------------------------------------------------------------------------------#
      }else{
         #----- Cannot use this point cloud. ----------------------------------------------#
         stop("Missing elements 'z' and 'intensity' from pt.cloud.")
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
   zbreaks = c(0,seq(from=zl,to=zh,length.out=nz))
   zmid    = mid.points(zbreaks)
   deltaz  = diff(zbreaks)
   #---------------------------------------------------------------------------------------#


   #----- Keep only the points that are within bounds. ------------------------------------#
   keep     = pt.cloud$z %>=% 0 & pt.cloud$z %<% zh
   pt.cloud = pt.cloud[keep,]
   #---------------------------------------------------------------------------------------#



   #----- Find the total energy returned from each layer. ---------------------------------#
   zcut         = cut(x=pt.cloud$z,breaks=zbreaks,right=FALSE)
   zcut         = match(zcut,levels(zcut))
   int.lev      = rep(0,times=nz)
   aux          = tapply(X=pt.cloud$intensity,INDEX=zcut,FUN=sum)
   idx          = as.numeric(names(aux))
   int.lev[idx] = aux
   #---------------------------------------------------------------------------------------#


   #----- Find the total intensity at the top of each layer. ------------------------------#
   int.tot  = sum(pt.cloud$intensity)
   int.sum  = rev(cumsum(rev(int.lev)))
   int.bot  = int.tot - int.sum
   int.top  = c(int.bot[-1],int.tot)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Find the unscaled values for both LAD and the cumulative.  Skip the first layer  #
   # to avoid infinite values. In case k is zero, then we normalise the profile.           #
   #---------------------------------------------------------------------------------------#
   lad    = log(int.top/int.bot) / deltaz
   lad    = lad   [-1]
   zmid   = zmid  [-1]
   deltaz = deltaz[-1]
   #---------------------------------------------------------------------------------------#



   #----- Scale the leaf area density to the layers we will look at. ----------------------#
   keep   = zmid >= zo
   zmid   = zmid  [keep]
   lad    = lad   [keep]
   deltaz = deltaz[keep]
   lad    = lad / sum(lad * deltaz)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Smooth the curve by generating a density function.  This density function becomes #
   # the output.                                                                           #
   #---------------------------------------------------------------------------------------#
   dzbar = mean(diff(zmid))
   zmah  = jitter( x      = sample(x=zmid,size=int.tot,replace=TRUE,prob=lad)
                 , amount = 0.5*dzbar
                 )#end jitter
   ans   = density.safe(zmah,from=zl,to=zh,n=nz)
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function macarthur.horn
#==========================================================================================#
#==========================================================================================#
