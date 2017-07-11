#==========================================================================================#
#==========================================================================================#
#     This function converts lists into data frames, preserving original types (as long as #
# the types are defined here, you may need to add other types.                             #
#------------------------------------------------------------------------------------------#
list.2.data.frame <<- function(x){
   if (! is.list(x)){
      stop("x must be a list!")
   }else if(length(x) == 0){
      out = data.frame()
   }else{
      xlen          = sapply(X = x, FUN = length)

      #------ Check that the lists all match. ---------------------------------------------#
      if (length(unique(xlen)) != 1){
         stop("All list elements must have the same length")
      }else{
         xlen = unique(xlen)
      }#end if
      #------------------------------------------------------------------------------------#


      #------ If lists are empty, bind but make it empty. ---------------------------------#
      if (xlen == 0){
         out           = data.frame(bye=numeric(length(x)))
         rownames(out) = names(x)
         out           = out[,-1]
      }else{
         #----- Make sure that the names all match. ---------------------------------------#
         names.ok  = apply(X=sapply(X=x,FUN=names),MARGIN=1,FUN=length.unique)
         if (any(names.ok) != 1){
            stop("All list must be the same (same names, and same order)")
         }#end if
         #---------------------------------------------------------------------------------#




         #----- Keep the type. ------------------------------------------------------------#
         xwhat = sapply(X = x[[1]], FUN = typeof, simplify = FALSE)
         xchron        = ( sapply(X= x[[1]], FUN = is.chron) 
                         | sapply(X= x[[1]], FUN = is.dates)
                         )#end 
         xwhat[xchron] = "chron"
         #---------------------------------------------------------------------------------#



         #----- Bind the lists and create data drame. -------------------------------------#
         out        = rbind( apply(X=rbind(sapply(X=x,FUN=c)),MARGIN=1,FUN=unlist)
                           , deparse.level = 1
                           )#end rbind
         names.out  = colnames(out)
         out        = split(x=out,f=col(out))
         names(out) = names.out
         #---------------------------------------------------------------------------------#


         #----- Make sure that variable types are preserved. ------------------------------#
         out = as.data.frame( mapply(FUN=as,object=out,Class=xwhat,SIMPLIFY=FALSE)
                            , stringsAsFactors = FALSE
                            )#end as.data.frame
         rownames(out) = names(x)
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#

   return(out)
}#end function
#==========================================================================================#
#==========================================================================================#
