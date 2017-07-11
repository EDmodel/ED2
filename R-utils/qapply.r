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


      #----- Call zapply by the apply function. -------------------------------------------#
      eout  = apply( X       = X
                   , MARGIN  = margin
                   , FUN     = function(z,zIND,zFUN,...) tapply(X=z,INDEX=zIND,FUN=zFUN,...)
                   , zIND    = INDEX
                   , zFUN    = FUN
                   ,...
                   )#end apply
      if (! is.list(INDEX)) INDEX = list(INDEX)
      uniqlist       = lapply(X=lapply(X=INDEX,FUN=sort),FUN=unique)
      dimuniq        = sapply(X=uniqlist,FUN=length,simplify=TRUE)
      eout           = array (data=eout,dim=c(dimuniq,dim(X)[margin]))
      dimnames(eout) = uniqlist
      off            = length(dimuniq) - 1
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
         }else{
            perm  = c(seq(from=2+off,to=DIM+off,by=1)
                     ,seq(from=    1,to=  1+off,by=1)
                     )#end c
         }#end if (DIM < ndims)
         eout  = aperm(a=eout,perm=perm)
      }#end if (DIM > 1)
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#

   return(eout)
}#end function qapply
#==========================================================================================#
#==========================================================================================#
