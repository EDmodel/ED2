#==========================================================================================#
#==========================================================================================#
#      Find the maximum defined rectangle that is entirely defined.                        #
#------------------------------------------------------------------------------------------#
rect.maxdef <<- function(mat,exhaustive=FALSE,stpfac = 1/10,verbose=FALSE){
   #----- Make sure mat is a matrix. ------------------------------------------------------#
   err = FALSE
   if ( is.array(mat)){
      if (length(dim(mat)) != 2){
         err = TRUE
      }#end if (length(dim(mat)) != 2)
   }else if(! (is.matrix(mat) || is.data.frame(mat))){
      err = TRUE
   }#end if
   if (err) stop("mat must be a matrix, a data frame, or an array with dimension 2")
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Check step size.                                                                  #
   #---------------------------------------------------------------------------------------#
   if (! is.finite(stpfac)){
      stop("stpfac must be finite.")
   }else if (stpfac <= 0 || stpfac >= 1){
      stop("stpfac must be greater than 0 and less than 1.")
   }#end if (! is.finite(stpfac))
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Count the maximum number of contiguous measurements by row and by column.        #
   #---------------------------------------------------------------------------------------#
   fincol = t(apply(X = mat, MARGIN = 1, FUN = finite.counter))
   finrow = apply(X = mat, MARGIN = 2, FUN = finite.counter)
   finsiz = arrayInd(ind=which.max(finrow * fincol),.dim=dim(mat))
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #    Finsiz is a good candidate, but it may be overestimated.                           #
   #---------------------------------------------------------------------------------------#
   rowcheck = seq(from=1+finsiz[1]-finrow[finsiz],to=finsiz[1],by=1)
   colcheck = seq(from=1+finsiz[2]-fincol[finsiz],to=finsiz[2],by=1)
   nr.try   = sort(unique(finrow[finsiz[1],colcheck]))
   nc.try   = sort(unique(fincol[rowcheck,finsiz[2]]))
   nrc.try  = expand.grid(nr.try,nc.try)
   nsiz.try = nrc.try[,1]*nrc.try[,2]
   if (nsiz.try[1] > 1){
      nrc.try = rbind(c(1,1),nrc.try)
      nsiz.try = c(1,nsiz.try)
   }#end if (nsiz.try[1] > 1)
   o        = order(nsiz.try)
   nrc.try  = nrc.try[o,]
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Seek the largest finite array.                                                   #
   #---------------------------------------------------------------------------------------#
   row.sel = finsiz[1]
   col.sel = finsiz[2]
   fine    = all(is.finite(mat[row.sel,col.sel]))
   if (fine && exhaustive){
      #------------------------------------------------------------------------------------#
      #     At least a single cell case works.  Go through a "bisection" to find the       #
      # largest array.                                                                     #
      #------------------------------------------------------------------------------------#
      nmax   = nrow(nrc.try)
      loop.n = rev(sequence(nmax))
      fine   = FALSE
      n      = nmax + 1
      it     = 0
      while ((! fine) && (n > 1)){
         it      = it + 1
         n       = n  - 1

         #----- New guess is half way between the lower and the upper bound. --------------#
         row.try = seq(from=1+finsiz[1]-nrc.try[n,1],to=finsiz[1])
         col.try = seq(from=1+finsiz[2]-nrc.try[n,2],to=finsiz[2])
         #---------------------------------------------------------------------------------#
         fine    = all(is.finite(mat[row.try,col.try]))

         if (verbose){
            cat0(" -  Iteration: ",it,";   N      = ",n
                                     ,";   FINE   = ",fine
                )#end cat0
         }#end if (verbose)
         if (fine){
            #----- Guess worked, exit the loop. -------------------------------------------#
            row.sel = row.try
            col.sel = col.try
            #------------------------------------------------------------------------------#
         }#end if (fine)
         #---------------------------------------------------------------------------------#
      }#end while (fine & (n > 1))
      #------------------------------------------------------------------------------------#
   }else if (fine){
      #------------------------------------------------------------------------------------#
      #     At least a single cell case works.  Go through a "bisection" to find the       #
      # largest array.                                                                     #
      #------------------------------------------------------------------------------------#
      nmax   = nrow(nrc.try)
      n      = 1
      nlwr   = 1
      deltan = nmax - nlwr
      it     = 0
      while (deltan >= 1){
         it      = it + 1
         #----- New guess is half way between the lower and the upper bound. --------------#
         row.try = seq(from=1+finsiz[1]-nrc.try[n,1],to=finsiz[1])
         col.try = seq(from=1+finsiz[2]-nrc.try[n,2],to=finsiz[2])
         #---------------------------------------------------------------------------------#
         fine    = all(is.finite(mat[row.try,col.try]))

         if (verbose){
            cat0(" -  Iteration: ",it,";   N      = ",n
                                     ,";   DELTAN = ",deltan
                                     ,";   NLWR   = ",nlwr
                                     ,";   FINE   = ",fine
                )#end cat0
         }#end if
         if (fine){
            #----- Guess worked, try a larger step. ---------------------------------------#
            row.sel = row.try
            col.sel = col.try
            nlwr    = n
            deltan  = nmax - nlwr
            #------------------------------------------------------------------------------#
         }else{
            #----- Guess didn't work, shrink step size. -----------------------------------#
            deltan  = min(nmax - nlwr, floor(deltan * (1 - stpfac)))
            #------------------------------------------------------------------------------#
         }#end if
         n = nlwr + deltan
         #---------------------------------------------------------------------------------#
      }#end while
      #------------------------------------------------------------------------------------#
   }else{
      #----- Not even the single size worked.  Probably all points are undefined... -------#
      row.sel = integer(0)
      col.sel = integer(0)
      #------------------------------------------------------------------------------------#
   }#end if (fine)
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #      Answer is a list with indices for the largest defined matrix.                    #
   #---------------------------------------------------------------------------------------#
   ans = list(rows = row.sel, cols = col.sel)
   return(ans)
   #---------------------------------------------------------------------------------------#

}#end rect.maxdef
#==========================================================================================#
#==========================================================================================#



#==========================================================================================#
#==========================================================================================#
#     This function finds the number of contiguous finite numbers.  It resets every time   #
# it finds a zero.                                                                         #
#------------------------------------------------------------------------------------------#
finite.counter <<- function(x){
   zo = as.numeric(is.finite(x))

   #----- idx is used to define the continuous sequences with valid numbers. --------------#
   idx = 1 + cumsum(zo == 0)
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #    Use the indices to split the vector into continuous segments, then use cumsum to   #
   # count them.                                                                           #
   #---------------------------------------------------------------------------------------#
   ans = unlist(tapply(X=zo,INDEX=idx,FUN=cumsum))
   names(ans) = names(x)
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end finite.counter
#==========================================================================================#
#==========================================================================================#
