#==========================================================================================#
#==========================================================================================#
#     Find out whether the point is a good candidate for local maximum.                    #
#------------------------------------------------------------------------------------------#
is.peak <<- function( pt.cloud
                    , dwin         = 1.5
                    , retn.use     = 1
                    , clss.use     = c(4,5)
                    , zmin         = 1
                    ){ #end function

   #----- Check whether to use only first returns or all returns. -------------------------#
   idx.retn = with(pt.cloud, which(retn.number %in% retn.use))
   idx.clss = with(pt.cloud, which(pt.class    %in% clss.use))
   idx.agnd = with(pt.cloud, which(z %>=% zmin))
   idx.keep = intersect(idx.clss,intersect(idx.retn,idx.agnd))
   #---------------------------------------------------------------------------------------#

   #----- The default answer is a Boolean vector with no peak. ----------------------------#
   ans = rep(x=FALSE,times=nrow(pt.cloud))
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Make sure that at least one point has been selected.                              #
   #---------------------------------------------------------------------------------------#
   if (length(idx.keep) != 0){
      xy.keep   = pt.cloud[idx.keep,]

      #---- Find the distance matrix. -----------------------------------------------------#
      dist.mat       = sqrt( outer(X=xy.keep$x,Y=xy.keep$x,FUN = "-")^2 
                           + outer(X=xy.keep$y,Y=xy.keep$y,FUN = "-")^2
                           )#end sqrt
      diag(dist.mat) = Inf
      #------------------------------------------------------------------------------------#


      #---- Select only the elements that are within the window distance. -----------------#
      closeby = dist.mat <= dwin
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Split the elevation by columns, then evaluate each potential candidate .       #
      #------------------------------------------------------------------------------------#
      z.split  = as.list(xy.keep$z)
      z.list   = replicate(n=nrow(xy.keep),list(xy.keep$z))
      sel.list = split(x=closeby,f=col(closeby))

      pinnacle = mapply( FUN      = function(znow,z,sel){
                                       if (length(z[sel]) == 0){
                                          ans = TRUE
                                       }else{
                                          ans = all(znow > z[sel],na.rm=TRUE)
                                       }#end if
                                       return(ans)
                                    }#end function
                       , znow     = z.split
                       , z        = z.list
                       , sel      = sel.list
                       , SIMPLIFY = TRUE
                       )#end mapply
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Correct the information in the answer vector.                                  #
      #------------------------------------------------------------------------------------#
      ans[idx.keep[pinnacle]] = TRUE
      #------------------------------------------------------------------------------------#


      #----- Free memory. -----------------------------------------------------------------#
      rm(xy.keep,dist.mat,closeby,z.split,z.list,sel.list)
      #------------------------------------------------------------------------------------#
   }#end if (length(idx.keep) != 0)
   #---------------------------------------------------------------------------------------#

   #----- Send result back to the parent environment. -------------------------------------#
   rm(idx.clss,idx.retn,idx.agnd,idx.keep)
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function is.peak
#==========================================================================================#
#==========================================================================================#
