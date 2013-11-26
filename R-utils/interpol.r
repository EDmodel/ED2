#------------------------------------------------------------------------------------------#
#    This function will perform the vertical interpolation of multi-dimensional arrays     #
# that contain the pressure as their "leftmost" dimension.  It will use the logarithm of   #
# pressure as the linear variable.  The result will be an array with the same size for     #
# dimensions other than the first, which will become the same size as the number of        #
# pressure levels.                                                                         #
#------------------------------------------------------------------------------------------#
interpol = function(varin,presin,presout,na.out=TRUE,frqval=0.5){

   if (! (is.array(varin) && is.array(presin))){
      stop("Both varin and presin must be arrays in order to use 'interpol'.")
   }#end if

   #----- Find the input variable dimensions and save them. -------------------------------#
   dimvarin  = dim(varin)              # The original variable array/matrix dimensions
   dimpresin = dim(presin)             # The original pressure array/matrix dimensions
   ndims     = length(dimvarin)        # The original number of dimensions

   npressin  = dimvarin[1]             # The input number of pressure levels
   nother    = prod(dimvarin[2:ndims]) # The product of the other dimensions of the output

   npressout = length(presout)         # The number of pressure levels of the output


   #------ Check whether the dimensions match. --------------------------------------------#
   if (any (dimvarin != dimpresin)){
      stop ("The dimensions of varin must be the same as the input pressure presin!!!")
   }#end if npressin

   #---------------------------------------------------------------------------------------#
   # "Collapse" all dimensions into two: pressure and everything else.                     #
   #---------------------------------------------------------------------------------------#
   varmat  = matrix(varin ,nrow=npressin ,ncol=nother)
   presmat = matrix(presin,nrow=npressin ,ncol=nother)
   matout  = matrix(NA    ,nrow=npressout,ncol=nother)

   #---------------------------------------------------------------------------------------#
   #      Not my favourite method, but I couldn't figure out a way without using a for     #
   # loop...                                                                               #
   #---------------------------------------------------------------------------------------#
   for (o in 1:nother){
   
      #----- Find the number of valid levels. ---------------------------------------------#
      use  = is.finite(presmat[,o])

      
      #------------------------------------------------------------------------------------#
      #     Here we determine if we have enough levels to interpolate, or if we should     #
      # simply forget about the interpolation this time and assign NA to the entire        #
      # column.                                                                            #
      #------------------------------------------------------------------------------------#
      if (sum(use) > frqval * npressout){
         out    = presout > max(presmat[use,o],na.rm=TRUE) |
                  presout < min(presmat[use,o],na.rm=TRUE)
         answer = aspline(x=log(presmat[,o]),y=varmat[,o],xout=log(presout))
         answer$y[out] = NA
         matout[,o] = answer$y
      }else{
         matout[,o] = rep(x=NA,times=npressout)
      }# end if
   }#end for (o in 1:nother)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Now we expand it back to the original size, except for the first dimension, and   #
   # that's it...                                                                          #
   #---------------------------------------------------------------------------------------#
   varout = array(matout,dim=c(npressout,dimvarin[2:ndims]))


   return(varout)
}# end function
#------------------------------------------------------------------------------------------#
