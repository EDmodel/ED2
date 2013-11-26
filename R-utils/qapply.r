#==========================================================================================#
#==========================================================================================#
# Function qapply.                                                                         #
# Developed by Marcos Longo - EPS/Harvard University                                       #
#                                                                                          #
#      This function is a combination of apply and tapply, so you can use tapply-like      #
# commands in matrices and arrays.                                                         #
#------------------------------------------------------------------------------------------#
qapply <<- function(X,INDEX,DIM,FUN,...){

   #---------------------------------------------------------------------------------------#
   #     Find the number of dimensions.                                                    #
   #---------------------------------------------------------------------------------------#
   dimexp = dim(X)
   ndims  = length(dimexp)
   if (is.null(dimexp)){
      dimexp = length(X)
      ndims  = 1
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Make sure that DIM <= ndims.                                                     #
   #---------------------------------------------------------------------------------------#
   if (length(DIM) != 1){
      stop(paste(" DIM must be a scalar! Yours has length ",length(DIM),"...",sep=""))
   }else if (DIM > ndims){
      cat (" - # of dimensions of X: ",ndims,"...","\n")
      cat (" - DIM: ",DIM,"...","\n")
      stop(" DIM must be less than or equal to the # of dimensions of X")
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Determine whether X is a matrix or an array of dimension 2.  If ndims=1, then    #
   # use tapply and return.                                                                #
   #---------------------------------------------------------------------------------------#
   if (ndims == 1){
      eout = tapply(X=X,INDEX=INDEX,FUN=FUN,...)
   }else{
      #------------------------------------------------------------------------------------#
      #   We apply "apply" to all margins but DIM.                                         #
      #------------------------------------------------------------------------------------#
      margin = sequence(ndims)
      margin = margin[-DIM]
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Define the "zapply" function, which is the tapply with different arguments.    #
      #------------------------------------------------------------------------------------#
      zapply = function(zx,zindex,zfunc,...){
         zout = tapply(X=zx,INDEX=zindex,FUN=zfunc,...)
         return(zout)
      } #end function zapply
      #------------------------------------------------------------------------------------#


      #----- Call zapply by the apply function. -------------------------------------------#
      eout  = apply(X=X,MARGIN=margin,FUN=zapply,zindex=INDEX,zfunc=FUN,...)
      if (is.list(INDEX) && length(INDEX) > 1){
         uniqlist       = sapply(X=lapply(X=INDEX,FUN=sort),FUN=unique)
         dimuniq        = sapply(X=uniqlist,FUN=length,simplify=TRUE)
         eout           = array (data=eout,dim=c(dimuniq,dim(X)[margin]))
         dimnames(eout) = uniqlist
         off            = length(dimuniq) - 1
      }else{
         off            = 0
      }#end if
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Check whether the original data was an array or a data frame.                 #
      #------------------------------------------------------------------------------------#
      if (DIM > 1){
         if(DIM < ndims){
            perm  = c(seq(from =     2+off,to = DIM+off,by=1)
                     ,seq(from =         1,to =   1+off,by=1)
                     ,seq(from = DIM+off+1,to =   ndims,by=1)
                     )#end c
            eout  = aperm(a=eout,perm=perm)
         }else{
            perm  = c(seq(from=2+off,to=DIM+off,by=1)
                     ,seq(from=    1,to=  1+off,by=1)
                     )#end c
            eout  = aperm(a=eout,perm=perm)
         }#end if (DIM < ndims)
      }#end if (DIM > 1)
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#

   return(eout)
}#end function qapply
#==========================================================================================#
#==========================================================================================#
