#------------------------------------------------------------------------------------------#
#     Create factors for a two-dimension surface.  This works exactly the same way as the  #
# cut function, except that the levels consider 2 variables.                               #
#------------------------------------------------------------------------------------------#
cut2d <<- function(x,y=NULL,xbreaks=10,ybreaks=10,xlabels=NULL,ylabels=NULL,...){

   #---------------------------------------------------------------------------------------#
   #     First thing, we check whether the input variables make sense.                     #
   #---------------------------------------------------------------------------------------#
   if (is.null(y)){
      #----- Point. -----------------------------------------------------------------------#
      if (is.data.frame(x)){
         if (("x" %in% names(x) & "y" %in% names(x))){
            xy = with(x,data.frame(x=x,y=y))
         }else if(dim(x)[2] == 2){
            xy = data.frame(x=x[,1],y=x[,2])
         }else{
            dum=bad.input.syntax("x")
         } 
      }else if (is.matrix(x)){
         numcol = dim(x)[2]
         if (numcol != 2) dum=bad.input.syntax("x")
         xy = data.frame(x=x[,1],y=x[,2])
      }else if (is.vector(x)){
         numcol = length(x)
         if (numcol != 2) dum=bad.input.syntax("x")

         #----- Make point a matrix. ------------------------------------------------------#
         xy    = data.frame(x=x[1],y=x[2])
      }else{
         dum = bad.input.syntax("x")
      }#end if
   }else{
      xy=data.frame(x=x,y=y)
   }#end if

   #---------------------------------------------------------------------------------------# 
   #    Split both X and Y.                                                                #
   #---------------------------------------------------------------------------------------#
   #----- X. ------------------------------------------------------------------------------#
   xcut  = cut(xy$x,breaks=xbreaks,labels=xlabels,...)
   xlevs = levels(xcut)
   nxlev = length(xlevs)
   #----- Y. ------------------------------------------------------------------------------#
   ycut  = cut(xy$y,breaks=ybreaks,labels=ylabels,...)
   ylevs = levels(ycut)
   nylev = length(ylevs)
   #----- Combine X and Y. ----------------------------------------------------------------#
   xycut  = paste(xcut,ycut,sep=",")
   xylevs = paste(rep(xlevs,times=nylev),rep(ylevs,each=nxlev),sep=",")

   #----- Convert xylevs back to levels. --------------------------------------------------#
   xycut = factor(x=xycut,levels=xylevs)
   return(xycut)
}#end function cut2d
#------------------------------------------------------------------------------------------#






#------------------------------------------------------------------------------------------#
#    This function is called to tell that cut2d had an invalid argument for point...       #
#------------------------------------------------------------------------------------------#
bad.input.syntax <<- function(varname){
      print(paste("In function in.poly:",varname,"is not a valid variable.",sep=" "))
      print("It must be:")
      print(" a. A two-element vector (longitude and latitude)")
      print(" b. A two-column matrix or data frame (1st. column lon and 2nd. column lat")
      stop (paste("Incorrect variable type for ",varname,"...",sep=""))
}#end function
#------------------------------------------------------------------------------------------#
