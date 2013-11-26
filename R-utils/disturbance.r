#==========================================================================================#
#==========================================================================================#
#     This function finds the normalised area given the disturbance rates.                 #
#------------------------------------------------------------------------------------------#
find.lu.areas <<- function(dlist,minarea){

   #----- Initialise table. ---------------------------------------------------------------#
   lu.names       = c("Pasture","Secondary","Primary","Cropland")
   dist.area = matrix(NA,ncol=4,nrow=dlist$nyears+1
                          ,dimnames=list(c(dlist$years,dlist$years[dlist$nyears]+dtyear)
                                        ,lu.names))
   #---------------------------------------------------------------------------------------#


   #----- First year, everything is primary vegetation. -----------------------------------#
   dist.area[1,] = c(0,0,1,0)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Integrate area over years.                                                        #
   #---------------------------------------------------------------------------------------#
   for (y in 1:dlist$nyears){
      #----- dA contains the areas that went from type r to type c. -----------------------#
      dA      = matrix(0,nrow=4,ncol=4,dimnames=list(lu.names,lu.names))
      #------------------------------------------------------------------------------------#


      #----- Go through all types. --------------------------------------------------------#
      dA[4,1] = dist.area[y,4] * (1. - exp(-dlist$disturb[y, 1]))
      dA[1,4] = dist.area[y,1] * (1. - exp(-dlist$disturb[y, 2]))
      dA[1,3] = dist.area[y,1] * (1. - exp(-dlist$disturb[y, 3]))
      dA[3,1] = dist.area[y,3] * (1. - exp(-dlist$disturb[y, 4]))
      dA[3,4] = dist.area[y,3] * (1. - exp(-dlist$disturb[y, 5]))
      dA[4,3] = dist.area[y,4] * (1. - exp(-dlist$disturb[y, 6]))
      dA[2,4] = dist.area[y,2] * (1. - exp(-dlist$disturb[y, 7]))
      dA[4,2] = dist.area[y,4] * (1. - exp(-dlist$disturb[y, 8]))
      dA[2,1] = dist.area[y,2] * (1. - exp(-dlist$disturb[y, 9]))
      dA[1,2] = dist.area[y,1] * (1. - exp(-dlist$disturb[y,10]))
      dA[3,2] = dist.area[y,3] * (1. - exp(-dlist$disturb[y,11]))
      #------------------------------------------------------------------------------------#



      #----- Update area. -----------------------------------------------------------------#
      dist.area[y+1,1] = max(minarea,dist.area[y,1]-sum(dA[1,])+sum(dA[,1]))
      dist.area[y+1,2] = max(minarea,dist.area[y,2]-sum(dA[2,])+sum(dA[,2]))
      dist.area[y+1,3] = max(minarea,dist.area[y,3]-sum(dA[3,])+sum(dA[,3]))
      dist.area[y+1,4] = max(minarea,dist.area[y,4]-sum(dA[4,])+sum(dA[,4]))
      dist.area[y+1, ] = dist.area[y+1,] / sum(dist.area[y+1,])
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#

   return(dist.area)
}#end function find.lu.areas
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function assess the error introduced by the lambda correction.                  #
#------------------------------------------------------------------------------------------#
lambda.err <<- function(x,past,futu,int.yeara,int.yearz,minarea,err.log=FALSE){

   histo = past


   #---- Select the years to apply the additional disturbance rate. -----------------------#
   sel = histo$years >= int.yeara & histo$years <= int.yearz
   #---------------------------------------------------------------------------------------#


   #----- Select the year to test. --------------------------------------------------------#
   yrtest = which(histo$years == int.yearz + 1)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Apply lambda to the years.                                                         #
   #---------------------------------------------------------------------------------------#
   histo$disturb[sel, 1] = histo$disturb[sel, 1] + inv.logit(x[ 1])
   histo$disturb[sel, 2] = histo$disturb[sel, 2] + inv.logit(x[ 2])
   histo$disturb[sel, 4] = histo$disturb[sel, 4] + inv.logit(x[ 4])
   histo$disturb[sel, 5] = histo$disturb[sel, 5] + inv.logit(x[ 5])
   histo$disturb[sel, 7] = histo$disturb[sel, 7] + inv.logit(x[ 7])
   histo$disturb[sel, 8] = histo$disturb[sel, 8] + inv.logit(x[ 8])
   histo$disturb[sel, 9] = histo$disturb[sel, 9] + inv.logit(x[ 9])
   histo$disturb[sel,10] = histo$disturb[sel,10] + inv.logit(x[10])
   histo$disturb[sel,11] = histo$disturb[sel,11] + inv.logit(x[11])
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Find the final area.                                                               #
   #---------------------------------------------------------------------------------------#
   lu.names       = c("Pasture","Secondary","Primary","Cropland")
   for (yr in int.yeara:int.yearz){
      y   = which(histo$years == yr)
      yp1 = y+1
      #----- dA contains the areas that went from type r to type c. -----------------------#
      dA      = matrix(0,nrow=4,ncol=4,dimnames=list(lu.names,lu.names))
      #------------------------------------------------------------------------------------#


      #----- Go through all types. --------------------------------------------------------#
      dA[4,1] = histo$norm.area[y,4] * (1. - exp(-histo$disturb[y, 1]))
      dA[1,4] = histo$norm.area[y,1] * (1. - exp(-histo$disturb[y, 2]))
      dA[1,3] = histo$norm.area[y,1] * (1. - exp(-histo$disturb[y, 3]))
      dA[3,1] = histo$norm.area[y,3] * (1. - exp(-histo$disturb[y, 4]))
      dA[3,4] = histo$norm.area[y,3] * (1. - exp(-histo$disturb[y, 5]))
      dA[4,3] = histo$norm.area[y,4] * (1. - exp(-histo$disturb[y, 6]))
      dA[2,4] = histo$norm.area[y,2] * (1. - exp(-histo$disturb[y, 7]))
      dA[4,2] = histo$norm.area[y,4] * (1. - exp(-histo$disturb[y, 8]))
      dA[2,1] = histo$norm.area[y,2] * (1. - exp(-histo$disturb[y, 9]))
      dA[1,2] = histo$norm.area[y,1] * (1. - exp(-histo$disturb[y,10]))
      dA[3,2] = histo$norm.area[y,3] * (1. - exp(-histo$disturb[y,11]))
      #------------------------------------------------------------------------------------#



      #----- Update area. -----------------------------------------------------------------#
      histo$norm.area[y+1,1] = max(minarea,histo$norm.area[y,1]-sum(dA[1,])+sum(dA[,1]))
      histo$norm.area[y+1,2] = max(minarea,histo$norm.area[y,2]-sum(dA[2,])+sum(dA[,2]))
      histo$norm.area[y+1,3] = max(minarea,histo$norm.area[y,3]-sum(dA[3,])+sum(dA[,3]))
      histo$norm.area[y+1,4] = max(minarea,histo$norm.area[y,4]-sum(dA[4,])+sum(dA[,4]))
      histo$norm.area[y+1, ] = histo$norm.area[y+1,] / sum(histo$norm.area[y+1,])
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #     Determine the relative error.                                                     #
   #---------------------------------------------------------------------------------------#
   if (err.log){
      expected  = log(futu$norm.area)
      estimated = log(histo$norm.area[yrtest,])
   }else{
      expected  = futu$norm.area
      estimated = histo$norm.area[yrtest,]
   }#end if
   bias         = (estimated-expected)
   mse          = bias^2
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Return the root mean squared error.                                               #
   #---------------------------------------------------------------------------------------#
   answer = sqrt(sum(mse))
   return(answer)
   #---------------------------------------------------------------------------------------#

}#end function lambda.err
#==========================================================================================#
#==========================================================================================#
