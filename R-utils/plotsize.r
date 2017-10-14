#==========================================================================================#
#==========================================================================================#
#    This function that defines the size of the figure to be plotted in case of maps.      #
# In case the plot is a map, it correct sizes so the map doesn't look distorted.           #
#------------------------------------------------------------------------------------------#
plotsize <<- function( proje                  #  Map projection? [T|F]
                     , limlon     = NULL      #  Longitude range, if proje = TRUE
                     , limlat     = NULL      #  Latitude range, if proje = TRUE
                     , deg        = TRUE      #  Are longitude and latitude in degrees?
                     , stdheight  = NULL      #  Standard height
                     , stdwidth   = NULL      #  Standard
                     , extendfc   = FALSE     #  Extend width for filled.contour [T|F]
                                              #  TRUE/FALSE  -- True = yes for longitude
                                              #  "lon","lat" -- will extend the specific
                                              #     dimension
                                              #  "both" -- will extend both dimensions
                                              #     extfactor may be a vector of two, in
                                              #     which case the first is applied to "lon"
                                              #     and the second is applied to "lat")
                     , extfactor  = 1/6       #  Factor to extend width or height
                     , paper      = "letter"  #  Paper size (ignored if stdXXX aren't NULL)
                     , landscape  = TRUE      #  Landscape? (if not swap width and height)
                     , scale.fac  = 0.8       #  Scaling factor to adjust sizes
                     ){


   null.std = is.null(stdheight) | is.null(stdwidth)

   #---------------------------------------------------------------------------------------#
   #     Check whether projection is TRUE or false.  In case it is TRUE, limlon and limlat #
   # must be given and must be a vector with dimension 2 (longitude and latitude ranges).  #
   #---------------------------------------------------------------------------------------#
   if (proje){
      if (is.null(limlon) && is.null(limlat)){
         stop("Variables limlon and limlat must be defined if proje is TRUE!")
      }else if(length(limlon) != 2 && length(limlat) !=2){
         stop("Variables limlon and limlat must be vectors of length 2 if proje is TRUE!")
      }#end if (is.null(limlon) && is.null(limlat))
      #------------------------------------------------------------------------------------#
   }#end if (proje)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Check whether stdheight and stdwidth are both available or both missing.  It is  #
   # not allowed to provide only one of them.  In case they are fine and they are          #
   # provided, override paper and make it "special".  In case they are both NULL, make the #
   # paper name lower case (so it becomes case insentitive).                               #
   #---------------------------------------------------------------------------------------#
   if (is.null(stdheight) != is.null(stdwidth)){
      cat(" Stdheight is NULL: ",is.null(stdheight),"\n")
      cat(" Stdwidth  is NULL: ",is.null(stdwidth ),"\n")
      stop(" Either set both stdheight and stdwidth to NULL or provide both...")
   }else if (! is.null(stdheight)){
      paper = "special"
   }else if (!is.null(paper)){
      paper = tolower(paper)
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #       Find the standard width and height depending on the paper.                      #
   #---------------------------------------------------------------------------------------#
   if (paper %in% "special") { 
      stdratio = max(c(stdwidth,stdheight))/min(c(stdwidth,stdheight))
   }else if (paper %in% "letter"){
      stdwidth  =  scale.fac * 11.0
      stdheight =  scale.fac *  8.5
      stdratio  = 11.0 /  8.5
   }else if (paper %in% "a4"){
      stdwidth  =  scale.fac * 29.7 / 2.54
      stdheight =  scale.fac * 21.0 / 2.54
      stdratio  = 29.7 / 21.0
   }else if (paper %in% "legal"){
      stdwidth  =  scale.fac * 14.0
      stdheight =  scale.fac *  8.5
      stdratio  = 14.0 /  8.5
   }else if (paper %in% "long"){
      stdwidth  =  scale.fac * 16.0
      stdheight =  scale.fac *  9.5
      stdratio  = 16.0 /  9.5
   }else if (paper %in% "executive"){
      stdwidth  =  scale.fac * 10.25
      stdheight =  scale.fac *  7.25
      stdratio  = 10.25 /  7.25
   }else if (paper %in% "double"){
      stdwidth  =  scale.fac * 14.0
      stdheight =  scale.fac *  7.0
      stdratio  = 14.0 / 7.0
   }else if (paper %in% "square"){
      stdwidth  =  scale.fac * 10.0
      stdheight =  scale.fac * 10.0
      stdratio  = 1.0
   }else{
      warning(paste("Unknown paper size (",paper,").  Using letter instead.",sep=""))
      stdwidth  =  scale.fac * 11.0
      stdheight =  scale.fac *  8.5
      stdratio  = 11.0 /  8.5
   }#end if 
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Correct the width and height in case this is a map.                                #
   #---------------------------------------------------------------------------------------#
   if (proje){
      #----- Find the actual ratio using the longitude and latitude. ----------------------#
      interx = max(limlon) - min(limlon)
      intery = max(limlat) - min(limlat)
      if (deg){
         ratio = interx * cos(mean(limlat) * pi / 180.) / intery
      }else{
         ratio = interx * cos(mean(limlat)) / intery
      }#end if (deg)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Fix the width or height to account for the sought ratio.                       #
      #------------------------------------------------------------------------------------#
      if (ratio >= stdratio){ 
         height = stdwidth  / ratio
         width  = stdwidth
      }else{
         height = stdheight
         width  = stdheight * ratio
      }#end if(actualratio >= stdratio)
      #------------------------------------------------------------------------------------#

   }else{
      #------------------------------------------------------------------------------------#
      #     Standard height/width, out-of-the-box.                                         #
      #------------------------------------------------------------------------------------#
      width  = stdwidth
      height = stdheight
      #------------------------------------------------------------------------------------#
   }#end if (proje)
   #---------------------------------------------------------------------------------------#



   #----- Extend the width in case this will be used for filled.contour. ------------------#
   if (is.logical(extendfc)){
      width.fac  = 1.0 + extfactor * as.numeric(extendfc)
      height.fac = 1.0
   }else if (tolower(substring(extendfc,1,2)) %in% "lo"){
      width.fac  = 1.0 + extfactor
      height.fac = 1.0
   }else if (tolower(substring(extendfc,1,2)) %in% "la"){
      width.fac  = 1.0
      height.fac = 1.0 + extfactor
   }else if (tolower(substring(extendfc,1,2)) %in% "bo"){
      width.fac  = 1.0 + extfactor[1]
      height.fac = 1.0 + extfactor[min(length(extfactor),2)]
   }else{
      width.fac  = 1.0
      height.fac = 1.0
   }#end if extendfc
   #---------------------------------------------------------------------------------------#




   #----- Not a map projection.  Use the standard size. -----------------------------------#
   height = height * height.fac
   width  = width  * width.fac
   ratio  = width  / height
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #       Swap height and width if portrait.                                              #
   #---------------------------------------------------------------------------------------#
   if (! landscape & null.std ){
      phold = height
      height = width
      width  = phold
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Append everything to a list.                                                      #
   #---------------------------------------------------------------------------------------#
   ans = list(height=height,width=width,ratio=ratio,paper="special")
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function
#==========================================================================================#
#==========================================================================================#

