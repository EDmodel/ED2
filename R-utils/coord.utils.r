#==========================================================================================#
#==========================================================================================#
#     This function finds the corners of a list of coordinates.                            #
#------------------------------------------------------------------------------------------#
four.corners <<- function(x,y=NULL,delta=0.005,verbose=FALSE){

   #----- Check whether y has been provided.  In case not, check the x object. ------------#
   if (is.null(y)){
      xy = x
      #----- Coerce data to a data frame. -------------------------------------------------#
      if (! is.data.frame(xy)){
         xy = try(as.data.frame(xy),silent=TRUE)
         if ("try-error" %in% xy){
            stop("x cannot be coerced to a data frame and y is not provided.")
         }#end if
      }#end if
      #------------------------------------------------------------------------------------#



      #----- Grab the columns with the x and y coordinates. -------------------------------#
      if (ncol(xy) < 2){
         stop("xy must have at least two columns!")
      }else{
         names(xy) = tolower(names(xy))
         if ("x" %in% names(xy)){x=xy$x}else{x=xy[[1]]}
         if ("y" %in% names(xy)){y=xy$y}else{y=xy[[2]]}
         
      }#end if
      #------------------------------------------------------------------------------------#
   }else if (length(x) != length(y)){
      stop("x and y must have the same length!")
   }#end if
   #----- Turn data set to a data frame. --------------------------------------------------#
   xy = data.frame(x = x, y = y)
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #   Find the coordinate range.                                                          #
   #---------------------------------------------------------------------------------------#
   idx = union( which(xy$y %in% min(xy$y,na.rm=TRUE))
              , union( which(xy$y %in% max(xy$y,na.rm=TRUE))
                     , union( which(xy$x %in% min(xy$x,na.rm=TRUE))
                            , which(xy$x %in% max(xy$x,na.rm=TRUE))
                            )#end union
                     )#end union
              )#end union
   corner = xy[idx,]
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Find the "ideal" corners, i.e., the corners of a rectangle that contains all       #
   # points and is perfectly aligned with the x and y axes.                                #
   #---------------------------------------------------------------------------------------#
   xa               = min(corner$x)
   xz               = max(corner$x)
   ya               = min(corner$y)
   yz               = max(corner$y)
   domain           = data.frame(x=c(xa,xa,xz,xz),y=c(ya,yz,yz,ya))
   rownames(domain) = c("sw","nw","ne","se")
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Decide the corner labels based on the closest distance to the "ideal" corners,    #
   # and label them according to the minimum distance to the perfect corners.              #
   #---------------------------------------------------------------------------------------#
   idx.sw = which.min(sqrt((corner$x-domain$x[1])^2+(corner$y-domain$y[1])^2))
   idx.nw = which.min(sqrt((corner$x-domain$x[2])^2+(corner$y-domain$y[2])^2))
   idx.ne = which.min(sqrt((corner$x-domain$x[3])^2+(corner$y-domain$y[3])^2))
   idx.se = which.min(sqrt((corner$x-domain$x[4])^2+(corner$y-domain$y[4])^2))
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Re-order the corners so the corners always go clockwise, starting from the SW     #
   # corner.                                                                               #
   #---------------------------------------------------------------------------------------#
   corner           = corner[c(idx.sw,idx.nw,idx.ne,idx.se),]
   rownames(corner) = c("sw","nw","ne","se")
   #---------------------------------------------------------------------------------------#



   #=======================================================================================#
   #=======================================================================================#
   #    If the point cloud has bulges, then some points may be outside the rectangle.      #
   # Adjust each corner until to make sure everyone is included.                           #
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Make sure that the current corners include all points.                             #
   #---------------------------------------------------------------------------------------#
   out.now  = sum(! inout(pts=xy,poly=rbind(corner,corner[1,])))
   if (out.now > 0){
      #----- Print message. ---------------------------------------------------------------#
      if (verbose){
         cat(" + Expanding the rectangle to include ",out.now," outsiders.","\n",sep="")
      }#end if
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Iterate until all points are included.  We go over each corner and expand the  #
      # vertices, one at a time, until all points are included.                            #
      #------------------------------------------------------------------------------------#
      ans     = corner
      n       = 0
      s       = 4
      while (out.now > 0 && n < 100){
         #---- Update counters and outsider count. ----------------------------------------#
         n        = n + (s == 4)
         s        = s %% 4 + 1
         out.prev = out.now
         #---------------------------------------------------------------------------------#

         #----- Expand corner. ------------------------------------------------------------#
         op        = ( (s-3) %% 4 ) + 1
         long      = ans
         long[s,]  = long[s,] + n * delta * (long[s,] - long[op,])
         #---------------------------------------------------------------------------------#


         #----- Update outsiders. ---------------------------------------------------------#
         out.now = sum(! inout(pts=xy,poly=rbind(long,long[1,])))
         #---------------------------------------------------------------------------------#


         #----- Check whether this expansion is helping. ----------------------------------#
         if (out.now < out.prev) ans[s,] = long[s,]
         #---------------------------------------------------------------------------------#


         #----- Entertain user. -----------------------------------------------------------#
         if (verbose){
            cat(" + Iteration: ",n,"; side: ",s,".  # outsiders: ",out.now,".","\n",sep="")
         }#end if
         #---------------------------------------------------------------------------------#
      }#end while (out.now && n < 100)
      #------------------------------------------------------------------------------------#
   }else{
      ans = corner
   }#end if (out.now > 0)
   #---------------------------------------------------------------------------------------#


   #---- Return the corners. --------------------------------------------------------------#
   return(ans)
   #---------------------------------------------------------------------------------------#

}#end function four.corners
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     Function that creates the grid mesh.                                                 #
#------------------------------------------------------------------------------------------#
grid.mesh <<- function(corners,nx,ny){


   #---- Find the normalised mesh. --------------------------------------------------------#
   x  = ( sequence(nx+1)-1 ) / nx
   y  = ( sequence(ny+1)-1 ) / ny
   xy = expand.grid(x,y)
   xx = xy[,1]
   yy = xy[,2]
   #---------------------------------------------------------------------------------------#


   #----- Project the mesh onto the grid. -------------------------------------------------#
   ee = ( (1 - xx) * (1 - yy) * corners$x[1] +      xx  * (1 - yy) * corners$x[2]
        +      xx  *      yy  * corners$x[3] + (1 - xx) *      yy  * corners$x[4] )
   nn = ( (1 - xx) * (1 - yy) * corners$y[1] +      xx  * (1 - yy) * corners$y[2]
        +      xx  *      yy  * corners$y[3] + (1 - xx) *      yy  * corners$y[4] )
   ans      = array(data=NA,dim=c(nx+1,ny+1,2),dimnames=list(NULL,NULL,c("x","y")))
   ans[,,1] = ee
   ans[,,2] = nn
   return(ans)
   #---------------------------------------------------------------------------------------#


}#end function grid.mesh
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function transforms normalised X;Y coordinates into native coordinates, given   #
# the four corners of the domain (always in the SW;NW;NE;SE order).  This function doesn't #
# assume that the shape is a rectangle.                                                    #
#------------------------------------------------------------------------------------------#
norm.to.coord <<- function(x,y=NULL,corners){

   #----- Check whether y has been provided.  In case not, check the x object. ------------#
   if (is.null(y)){
      xy = x
      #----- Coerce data to a data frame. -------------------------------------------------#
      if (! is.data.frame(xy)){
         xy = try(as.data.frame(xy),silent=TRUE)
         if ("try-error" %in% xy){
            stop("x cannot be coerced to a data frame and y is not provided.")
         }#end if
      }#end if
      #------------------------------------------------------------------------------------#



      #----- Grab the columns with the x and y coordinates. -------------------------------#
      if (ncol(xy) < 2){
         stop("xy must have at least two columns!")
      }else{
         names(xy) = tolower(names(xy))
         if ("x" %in% names(xy)){x=xy$x}else{x=xy[[1]]}
         if ("y" %in% names(xy)){y=xy$y}else{y=xy[[2]]}
         
      }#end if
      #------------------------------------------------------------------------------------#
   }else if (length(x) != length(y)){
      stop("x and y must have the same length!")
   }#end if
   #----- Turn data set to a data frame. --------------------------------------------------#
   xy = data.frame(x = x, y = y)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Make sure that 'corners' is a 4x2 matrix.
   #---------------------------------------------------------------------------------------#
   #----- Coerce data to a data frame. ----------------------------------------------------#
   if (! is.data.frame(corners)){
      corners = try(as.data.frame(corners),silent=TRUE)
      if ("try-error" %in% corners){
         stop("corners cannot be coerced to a data frame.")
      }#end if
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Grab the columns with the x and y coordinates. ----------------------------------#
   if (any(dim(corners) != c(4,2))){
      stop("corners must have four rows and two columns!")
   }else{
      names(corners) = tolower(names(corners))
      if ("x" %in% names(corners)) names(corners)[1] = "x"
      if ("y" %in% names(corners)) names(corners)[2] = "y"
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Find some auxiliary variables.                                                    #
   #---------------------------------------------------------------------------------------#
   xx = xy$x
   yy = xy$y
   #---------------------------------------------------------------------------------------#


   #----- Project the points onto the normalised grid. ------------------------------------#
   ee = ( (1 - xx) * (1 - yy) * corners$x[1] +      xx  * (1 - yy) * corners$x[2]
        +      xx  *      yy  * corners$x[3] + (1 - xx) *      yy  * corners$x[4] )
   nn = ( (1 - xx) * (1 - yy) * corners$y[1] +      xx  * (1 - yy) * corners$y[2]
        +      xx  *      yy  * corners$y[3] + (1 - xx) *      yy  * corners$y[4] )
   #---------------------------------------------------------------------------------------#


   #---- Return a data frame with the normalised variables. -------------------------------#
   ans = data.frame(x=ee,y=nn)
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end norm.to.coord
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function normalises the X and Y coordinates of quadrilaters.  It doesn't        #
# assume that the shape is a rectangle.                                                    #
#     This function normalises native coordinates (X 0-1; Y 0-1), given the four corners   #
# of the domain (always in the SW;NW;NE;SE order).  This function doesn't assume that the  #
# shape is a rectangle.                                                                    #
#------------------------------------------------------------------------------------------#
coord.to.norm <<- function(x,y=NULL,corners){

   #----- Check whether y has been provided.  In case not, check the x object. ------------#
   if (is.null(y)){
      xy = x
      #----- Coerce data to a data frame. -------------------------------------------------#
      if (! is.data.frame(xy)){
         xy = try(as.data.frame(xy),silent=TRUE)
         if ("try-error" %in% xy){
            stop("x cannot be coerced to a data frame and y is not provided.")
         }#end if
      }#end if
      #------------------------------------------------------------------------------------#



      #----- Grab the columns with the x and y coordinates. -------------------------------#
      if (ncol(xy) < 2){
         stop("xy must have at least two columns!")
      }else{
         names(xy) = tolower(names(xy))
         if ("x" %in% names(xy)){x=xy$x}else{x=xy[[1]]}
         if ("y" %in% names(xy)){y=xy$y}else{y=xy[[2]]}
         
      }#end if
      #------------------------------------------------------------------------------------#
   }else if (length(x) != length(y)){
      stop("x and y must have the same length!")
   }#end if
   #----- Turn data set to a data frame. --------------------------------------------------#
   xy = data.frame(x = x, y = y)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Make sure that 'corners' is a 4x2 matrix.
   #---------------------------------------------------------------------------------------#
   #----- Coerce data to a data frame. ----------------------------------------------------#
   if (! is.data.frame(corners)){
      corners = try(as.data.frame(corners),silent=TRUE)
      if ("try-error" %in% corners){
         stop("corners cannot be coerced to a data frame.")
      }#end if
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Grab the columns with the x and y coordinates. ----------------------------------#
   if (any(dim(corners) != c(4,2))){
      stop("corners must have four rows and two columns!")
   }else{
      names(corners) = tolower(names(corners))
      if ("x" %in% names(corners)) names(corners)[1] = "x"
      if ("y" %in% names(corners)) names(corners)[2] = "y"
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Find some auxiliary variables.                                                    #
   #---------------------------------------------------------------------------------------#
   et = xy$x         - corners$x[1]
   eb = corners$x[4] - corners$x[1]
   ec = corners$x[3] - corners$x[1]
   ed = corners$x[2] - corners$x[1]
   e3 = ec - eb - ed
   nt = xy$y         - corners$y[1]
   nb = corners$y[4] - corners$y[1]
   nc = corners$y[3] - corners$y[1]
   nd = corners$y[2] - corners$y[1]
   n3 = nc - nb - nd
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Find coefficients for quadratic formula.                                           #
   #---------------------------------------------------------------------------------------#
   aa = (nd * e3 - ed * n3)
   bb = ( nd * eb - nb * ed + n3 * et - nt * e3 ) / aa
   cc = ( nb * et - nt * eb )  / aa
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Find the y coordinates.                                                           #
   #---------------------------------------------------------------------------------------#
   dd = bb * bb - 4 * cc
   y1 = 0.5 * ( - bb - sqrt(dd) )
   y2 = 0.5 * ( - bb + sqrt(dd) )
   yy = ifelse( y1 %wr% c(0,1), y1, ifelse( y2 %wr% c(0,1), y2, NA) )
   xx = (et - yy * ed) / (eb + yy * e3)
   #---------------------------------------------------------------------------------------#


   #---- Return a data frame with the normalised variables. -------------------------------#
   ans = data.frame(x=xx,y=yy)
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end coord.to.norm
#==========================================================================================#
#==========================================================================================#
