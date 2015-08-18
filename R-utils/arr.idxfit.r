#==========================================================================================#
#==========================================================================================#
#     Function arr.idxfit.r                                                                #
#                                                                                          #
#   This function fits a xy plane to the values of the matrix.  Despite the name, it also  #
# works for arrays.                                                                        #
#------------------------------------------------------------------------------------------#
arr.idxfit <<- function(A){
   #----- Make sure the matrix is an array. -----------------------------------------------#
   if (is.data.frame(A)){
      A = as.matrix(A)
   }else if (! ( is.matrix(A) | is.array(A) )){
      stop("A must be a matrix, an array, or a data frame...")
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Create the response variables.                                                    #
   #---------------------------------------------------------------------------------------#
   dfarr        = data.frame(arrayInd(seq_along(A),.dim=dim(A)))
   idnum        = ceiling(log10(length(dim(A))))
   fmt          = paste0("%",idnum,".",idnum,"i")
   names(dfarr) = paste0("x",sprintf(fmt=fmt,seq_along(dim(A))))
   dfarr$y      = c(A)
   #---------------------------------------------------------------------------------------#

   #----- Fit the plane equation. ---------------------------------------------------------#
   afit  = try(lm(formula = y ~ .,data=dfarr),silent=TRUE)
   if ("try-error" %in% is(afit)){
      apred = A * NA
   }else{
      apred = predict(object=afit,newdata=dfarr) + 0 * A
   }#end if ("try-error" %in% is(afit))
   return(apred)
   #---------------------------------------------------------------------------------------#

}#end arr.idxfit
#==========================================================================================#
#==========================================================================================#
