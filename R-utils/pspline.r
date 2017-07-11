#==========================================================================================#
#==========================================================================================#
#     This function calculates a parametric spline or linear interpolation, so it          #
# interpolates both x and y as a function of the index.  All methods from splinefun are    #
# allowed, and "linear" and "approx" are also options, but in these cases approxfun is     #
# used instead.                                                                            #
#------------------------------------------------------------------------------------------#
pspline <<- function(x,y=NULL,nfac=3,absol=FALSE,method="monoH.FC",ties=mean){

   #----- Find the method. ----------------------------------------------------------------#
   method = match.arg( arg        = method
                     , choices    = c("monoH.FC","fmm","periodic","natural","hyman"
                                     ,"linear","approx")
                     , several.ok = FALSE
                     )#end match.arg
   #---------------------------------------------------------------------------------------#


   #----- Find out whether x and y were given or if they are both in x. -------------------#
   if (is.null(y)){
      names(x) = tolower(names(x))
      y        = x$y
      x        = x$x
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Make sure that NAs are preserved.                                                #
   #---------------------------------------------------------------------------------------#
   breaks  = which(is.na(x) | is.na(y))
   breaks  = breaks[! (breaks %in% c(1,length(x)))]
   #---------------------------------------------------------------------------------------#


   #----- Find the break points. ----------------------------------------------------------#
   block.a = c(1,breaks+1)
   block.z = c(breaks-1,length(x))
   nblocks = length(block.a)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Initialise output vectors.                                                       #
   #---------------------------------------------------------------------------------------#
   xout = NULL
   yout = NULL
   for (b in sequence(nblocks)){

      #----- Find extremes of this block. -------------------------------------------------#
      ua   = block.a[b]
      uz   = block.z[b]
      nu   = uz-ua+1
      u    = seq(from=ua,to=uz,by=1)
      uout = seq(from=ua,to=uz,length.out=ifelse(absol,nfac,nfac*nu))
      #------------------------------------------------------------------------------------#



      #---- Find the interpolation functions. ---------------------------------------------#
      if (method %in% c("linear","approx")){
         xfun = approxfun(x=u,y=x[u],method="linear",ties=ties)
         yfun = approxfun(x=u,y=y[u],method="linear",ties=ties)
      }else{
         xfun = splinefun(x=u,y=x[u],method=method,ties=ties)
         yfun = splinefun(x=u,y=y[u],method=method,ties=ties)
      }#end if
      #------------------------------------------------------------------------------------#



      #----- Append output to the spline. -------------------------------------------------#
      if (b == 1){
         xout = xfun(uout)
         yout = yfun(uout)
      }else{
         xout = c(xout,NA,xfun(uout))
         yout = c(yout,NA,yfun(uout))
      }#end if (b ==1)
      #------------------------------------------------------------------------------------#
   }#end for (b in sequence(nblocks))
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Create output data frame.                                                         #
   #---------------------------------------------------------------------------------------#
   ans = data.frame(x=xout,y=yout)
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function pspline
#==========================================================================================#
#==========================================================================================#
