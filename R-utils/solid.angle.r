#------------------------------------------------------------------------------------------#
#     This function computes the solid angle or the area associated, of a "rectangle"      #
# defined by two pairs of coordinates, representing the Southwestern and Northeaster       #
# corners.                                                                                 #
#                                                                                          #
#    INPUTS.                                                                               #
# 1. sw      = Southwestern corner (1st. column = longitude, 2nd column = latitude).       #
# 2. ne      = Northeastern corner (1st. column = longitude, 2nd column = latitude).       #
# 3. degrees = Are the longitude and latitude given in degrees? (Default: TRUE)            #
# 4. radius  = Radius of the sphere. The default, 1, gives the solid angle.                #
#------------------------------------------------------------------------------------------#
solid.angle = function(sw,ne,degrees=TRUE,radius=1.){

   #----- Assigning the multiplication factor ---------------------------------------------#
   if (degrees){
      fact = pi/180.
   }else{
      fact = 1.0
   }#end if
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #     Break the coordinates into vectors, and convert the coordinates to radians.       #
   #---------------------------------------------------------------------------------------#
   #----- 1. Southwestern corner. ---------------------------------------------------------#
   if (is.data.frame(sw) || is.matrix(sw)){
      wlon = sw[,1] * fact
      slat = sw[,2] * fact
   }else if(is.vector(sw) && length(sw) == 2){
      wlon = sw[1]  * fact
      slat = sw[2]  * fact
   }else{
      stop("sw must be a data frame, a matrix, or a vector of length 2...")
   }#end if
   #----- 2. Northeastern corner. ---------------------------------------------------------#
   if (is.data.frame(ne) || is.matrix(ne)){
      elon = ne[,1] * fact
      nlat = ne[,2] * fact
   }else if(is.vector(ne) && length(ne) == 2){
      elon = ne[1]  * fact
      nlat = ne[2]  * fact
   }else{
      stop("sw must be a data frame, a matrix, or a vector of length 2...")
   }#end if
   #---------------------------------------------------------------------------------------#






   #---------------------------------------------------------------------------------------#
   #     If the rectangle crosses the date line or Greenwich, it may happen that elon is   #
   # less than elon. Make sure that this is corrected.                                     #
   #---------------------------------------------------------------------------------------#
   sel       = elon < wlon
   elon[sel] = elon[sel] + 2 * pi

   #----- Compute the solid angle. --------------------------------------------------------#
   omega     = ( sin(nlat) - sin(slat) ) * ( elon - wlon ) * radius * radius


   return(omega)
}#end function
#------------------------------------------------------------------------------------------#
