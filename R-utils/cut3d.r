#------------------------------------------------------------------------------------------#
#     Create factors for a three-dimension volime.  This works exactly the same way as the #
# cut function, except that the levels consider 3 variables.                               #
#------------------------------------------------------------------------------------------#
cut3d <<- function(x,y=NULL,z=NULL,xbreaks=10,ybreaks=10,zbreaks=10
                ,xlabels=NULL,ylabels=NULL,zlabels=NULL,...){

   #---------------------------------------------------------------------------------------#
   #     First thing, we check whether the input variables make sense.                     #
   #---------------------------------------------------------------------------------------#
   if (is.null(y) && is.null(z)){
      #----- Point. -----------------------------------------------------------------------#
      if (is.data.frame(x)){
         if ("x" %in% names(x) & "y" %in% names(x) & "z" %in% names(x) ){
            xyz = data.frame(x=x$x,y=x$y,z=x$z)
         }else if(dim(x)[2] == 3){
            xyz = data.frame(x=x[,1],y=x[,2],z=x[,3])
         }else{
            dum=bad.configuration("x")
         } 
      }else if (is.matrix(x)){
         numcol = dim(x)[2]
         if (numcol != 3) dum=bad.configuration("x")
         xyz = data.frame(x=x[,1],y=x[,2],z=x[,3])
      }else if (is.vector(x)){
         numcol = length(x)
         if (numcol != 3) dum=bad.configuration("x")

         #----- Make point a matrix. ------------------------------------------------------#
         xyz   = data.frame(x=x[1],y=x[2],z=x[3])
      }else{
         dum = bad.configuration("x")
      }#end if
   }else{
      xyz=as.data.frame(cbind(x,y,z))
   }#end if

   #---------------------------------------------------------------------------------------# 
   #    Split X, Y, and Z.                                                                 #
   #---------------------------------------------------------------------------------------#
   #----- X. ------------------------------------------------------------------------------#
   xcut  = cut(xyz$x,breaks=xbreaks,labels=xlabels,...)
   xlevs = levels(xcut)
   nxlev = length(xlevs)
   #----- Y. ------------------------------------------------------------------------------#
   ycut  = cut(xyz$y,breaks=ybreaks,labels=ylabels,...)
   ylevs = levels(ycut)
   nylev = length(ylevs)
   #----- Z. ------------------------------------------------------------------------------#
   zcut  = cut(xyz$z,breaks=zbreaks,labels=zlabels,...)
   zlevs = levels(zcut)
   nzlev = length(zlevs)
   #----- Combine X and Y. ----------------------------------------------------------------#
   xyzcut  = paste(xcut,ycut,zcut,sep=",")
   xyzlevs = paste(rep(xlevs,times=nylev*nzlev),rep(ylevs,each=nxlev,times=nzlev),
                   rep(zlevs,each =nxlev*nylev),sep=",")

   #----- Convert xylevs back to levels. --------------------------------------------------#
   xyzcut = factor(x=xyzcut,levels=xyzlevs)
   return(xyzcut)
}#end function cut3d
#------------------------------------------------------------------------------------------#






#------------------------------------------------------------------------------------------#
#    This function is called to tell that cut2d had an invalid argument for point...       #
#------------------------------------------------------------------------------------------#
bad.configuration <<- function(varname){
      print(paste("In function in.poly:",varname,"is not a valid variable.",sep=" "))
      print("It must be:")
      print(" a. A two-element vector (longitude and latitude)")
      print(" b. A two-column matrix or data frame (1st. column lon and 2nd. column lat")
      stop (paste("Incorrect variable type for ",varname,"...",sep=""))
}#end function bad.configuration
#------------------------------------------------------------------------------------------#
