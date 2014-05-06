#==========================================================================================#
#==========================================================================================#
#     Function borrowed from nls2.                                                         #
#                                                                                          #
#    This function simply convert the NLS object to LM to get the names right.             #
# Effectively this function just does this:                                                #
#                                                                                          #
# lm( lhs ~ gradient - 1                                                                   #
#   , offset = fitted(object)                                                              #
#   , list(gradient = object$m$gradient()                                                  #
#   , lhs = object$m$lhs()))                                                               #
#                                                                                          #
# so most of the code is just to get the names right.                                      #
#------------------------------------------------------------------------------------------#
nls.2.lm <<- function(object, ...) {
   #----- Make sure that the object is NLS. -----------------------------------------------#
   if (! inherits(object, "nls")){
      cat(" Classes of input object: ",paste(class(object),collapse=" "),"\n",sep="")
      stop(" nls.2.lm requires NLS object, hence the name!")
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Check several names.                                                              #
   #---------------------------------------------------------------------------------------#
   gradient = object$m$gradient()
   if (is.null(colnames(gradient))) colnames(gradient) = names(object$m$getPars())

   if (length(formula(object)) == 2){
      response.name = "0"
   }else{
      response.name = as.character(formula(object)[[2]])
   }#end if
   #---------------------------------------------------------------------------------------#



   lhs = object$m$lhs()
   ell = data.frame(lhs, gradient)
   names(ell)[1] = response.name

   fo = sprintf("%s ~ %s - 1", response.name,paste(colnames(gradient), collapse = "+"))
   fo = as.formula(fo, env = as.proto.list(ell))

   #---------------------------------------------------------------------------------------#
   #    Return conversion to lm.                                                           #
   #---------------------------------------------------------------------------------------#
   fit.obj = fitted(object,...)
   ans     = do.call("lm", list(fo, offset = substitute(fit.obj)))
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function nls.2.lm
#==========================================================================================#
#==========================================================================================#
