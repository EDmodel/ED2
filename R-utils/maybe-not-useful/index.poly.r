#------------------------------------------------------------------------------------------#
#     Check whether a number is in a given interval.  Input data must contain three        #
# 2-column structures, one with the point coordinates (or a matrix with many points), and  #
# the other two with the southwestern and north-eastern edges.  The result will depend on  #
# the input data.  If the point is actually a matrix, then the result will be a matrix of  #
# Boolean results, with each row corresponding to the point, and each column corresponding #
# to the edge.  If the point is just a single point, the result will be a vector with the  #
# check.                                                                                   #
#------------------------------------------------------------------------------------------#
index.poly = function(point,swedge,needge){

   #---------------------------------------------------------------------------------------#
   #     First thing, we check whether the variables make sense.                           #
   #---------------------------------------------------------------------------------------#
   #----- Point. --------------------------------------------------------------------------#
   if (is.matrix(point) | is.data.frame(point)){
      numcol = dim(point)[2]
      if (numcol != 2) dum=invalid.variable("point")
      ptrow = dim(point)[1]
   }else if (is.vector(point)){
      numcol = length(point)
      if (numcol != 2) dum=invalid.variable("point")

      #----- Make point a matrix. ---------------------------------------------------------#
      point = matrix(point,ncol=2,nrow=1)
      ptrow = dim(point)[1]
   }else{
      dum = invalid.variable("point")
   }#end if
   #----- SW edge. ------------------------------------------------------------------------#
   if (is.matrix(swedge) | is.data.frame(swedge)){
      numcol = dim(swedge)[2]
      if (numcol != 2) dum=invalid.variable("swedge")
      swrow  = dim(swedge)[1]
   }else if (is.vector(swedge)){
      numcol = length(swedge)
      if (numcol != 2) dum=invalid.variable("swedge")
      swedge = matrix(swedge,nrow=1,ncol=2)
      swrow  = dim(swedge)[1]
   }else{
      dum = invalid.variable("swedge")
   }#end if
   #----- NE edge. ------------------------------------------------------------------------#
   if (is.matrix(needge) | is.data.frame(needge)){
      numcol = dim(needge)[2]
      if (numcol != 2) dum=invalid.variable("needge")
      nerow  = dim(needge)[1]
   }else if (is.vector(needge)){
      numcol = length(needge)
      if (numcol != 2) dum=invalid.variable("needge")
      needge = matrix(needge,nrow=1,ncol=2)
      nerow  = dim(needge)[1]
   }else{
      dum = invalid.variable("needge")
   }#end if
   #----- Now we check whether the SW and NE edges have the same size. --------------------#
   if (swrow != nerow) stop ("swedge and needge must have the same number of points.")


   #---------------------------------------------------------------------------------------#
   #    If we reached this point, we have valid data.                                      #
   #---------------------------------------------------------------------------------------#
   matmat    = matrix(NA,nrow=ptrow,ncol=nerow)
   polyindex = matrix(rep(1:nerow,each=ptrow),nrow=ptrow,ncol=nerow)

   #---- Errr... Here I will need a for loop... I promise this will be the only one :) .---#
   for (pp in 1:ptrow){
       matmat[pp,] = point[pp,1] > swedge[,1] & point[pp,1] <= needge[,1] &
                     point[pp,2] > swedge[,2] & point[pp,2] <= needge[,2]
   }#end (for pp in 1:ptrow)

   #---------------------------------------------------------------------------------------#
   #     Check to how many polygons each point belongs. No point can belong to two points  #
   # are inside each polygon. Currently no point can belong to more than one polygon, al-  #
   # though it is fine that each polygon has many points inside it.                        #
   #---------------------------------------------------------------------------------------#
   ninsi = apply(X=matmat,MARGIN=1,FUN=sum)
   if (any(ninsi) > 1) stop("Your polygons overlap and index.poly doesn't accept that!")
   polyindex[! matmat] = NA
   isinside = rowSums(polyindex,na.rm=TRUE)
   isinside[isinside == 0] = NA
   return(isinside)
}#end function index.poly
#------------------------------------------------------------------------------------------#






#------------------------------------------------------------------------------------------#
#    This function is called to tell that index.poly had an invalid argument for point...  #
#------------------------------------------------------------------------------------------#
invalid.variable = function(varname){
      print(paste("In function index.poly:",varname,"is not a valid variable.",sep=" "))
      print("It must be:")
      print(" a. A two-element vector (longitude and latitude)")
      print(" b. A two-column matrix or data frame (1st. column lon and 2nd. column lat")
      stop (paste("Incorrect variable type for ",varname,"...",sep=""))
}#end function
#------------------------------------------------------------------------------------------#
