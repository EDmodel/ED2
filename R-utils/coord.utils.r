#==========================================================================================#
#==========================================================================================#
#     This function finds the vertices of a list of coordinates.                           #
#------------------------------------------------------------------------------------------#
four.vertices <<- function(x,y=NULL,delta=0.005,verbose=FALSE){

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
   vertex = xy[idx,]
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Find the "ideal" vertices, i.e., the vertices of a rectangle that contains all     #
   # points and is perfectly aligned with the x and y axes.                                #
   #---------------------------------------------------------------------------------------#
   xa               = min(vertex$x)
   xz               = max(vertex$x)
   ya               = min(vertex$y)
   yz               = max(vertex$y)
   domain           = data.frame(x=c(xa,xa,xz,xz),y=c(ya,yz,yz,ya))
   rownames(domain) = c("sw","nw","ne","se")
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Decide the vertex labels based on the closest distance to the "ideal" vertices,   #
   # and label them according to the minimum distance to the perfect vertices.             #
   #---------------------------------------------------------------------------------------#
   idx.sw = which.min(sqrt((vertex$x-domain$x[1])^2+(vertex$y-domain$y[1])^2))
   idx.nw = which.min(sqrt((vertex$x-domain$x[2])^2+(vertex$y-domain$y[2])^2))
   idx.ne = which.min(sqrt((vertex$x-domain$x[3])^2+(vertex$y-domain$y[3])^2))
   idx.se = which.min(sqrt((vertex$x-domain$x[4])^2+(vertex$y-domain$y[4])^2))
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Re-order the vertices so the vertices always go clockwise, starting from the SW   #
   # vertex.                                                                               #
   #---------------------------------------------------------------------------------------#
   vertex           = vertex[c(idx.sw,idx.nw,idx.ne,idx.se),]
   rownames(vertex) = c("sw","nw","ne","se")
   #---------------------------------------------------------------------------------------#



   #=======================================================================================#
   #=======================================================================================#
   #    If the point cloud has bulges, then some points may be outside the rectangle.      #
   # Adjust each vertex until to make sure everyone is included.                           #
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Make sure that the current vertices include all points.                            #
   #---------------------------------------------------------------------------------------#
   out.now  = sum(! inout(pts=xy,poly=rbind(vertex,vertex[1,])))
   if (out.now > 0){
      #----- Print message. ---------------------------------------------------------------#
      if (verbose){
         cat(" + Expanding the rectangle to include ",out.now," outsiders.","\n",sep="")
      }#end if
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Iterate until all points are included.  We go over each vertex and expand the  #
      # vertices, one at a time, until all points are included.                            #
      #------------------------------------------------------------------------------------#
      ans     = vertex
      n       = 0
      s       = 4
      while (out.now > 0 && n < 100){
         #---- Update counters and outsider count. ----------------------------------------#
         n        = n + (s == 4)
         s        = s %% 4 + 1
         out.prev = out.now
         #---------------------------------------------------------------------------------#

         #----- Expand vertex. ------------------------------------------------------------#
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
      ans = vertex
   }#end if (out.now > 0)
   #---------------------------------------------------------------------------------------#


   #---- Return the vertices. -------------------------------------------------------------#
   return(ans)
   #---------------------------------------------------------------------------------------#

}#end function four.vertices
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     Function that creates the grid mesh.                                                 #
#------------------------------------------------------------------------------------------#
grid.mesh <<- function(vertex,nx,ny){


   #---- Find the normalised mesh. --------------------------------------------------------#
   x  = ( sequence(nx+1)-1 ) / nx
   y  = ( sequence(ny+1)-1 ) / ny
   xy = expand.grid(x,y)
   xx = xy[,1]
   yy = xy[,2]
   #---------------------------------------------------------------------------------------#


   #----- Project the mesh onto the grid. -------------------------------------------------#
   ee = ( (1 - xx) * (1 - yy) * vertex$x[1] +      xx  * (1 - yy) * vertex$x[2]
        +      xx  *      yy  * vertex$x[3] + (1 - xx) *      yy  * vertex$x[4] )
   nn = ( (1 - xx) * (1 - yy) * vertex$y[1] +      xx  * (1 - yy) * vertex$y[2]
        +      xx  *      yy  * vertex$y[3] + (1 - xx) *      yy  * vertex$y[4] )
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
# the four vertices of the domain (always in the SW;NW;NE;SE order).  This function        #
# doesn't assume that the shape is a rectangle.                                            #
#------------------------------------------------------------------------------------------#
norm.to.coord <<- function(x,y=NULL,vertex,xscl=1,yscl=1){

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


      #------------------------------------------------------------------------------------#
      #    Save row names for later.                                                       #
      #------------------------------------------------------------------------------------#
      rnmout = rownames(xy)
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
   }else{
      #------------------------------------------------------------------------------------#
      #    Save row names for later.                                                       #
      #------------------------------------------------------------------------------------#
      rnmout = names(x)
      #------------------------------------------------------------------------------------#
   }#end if
   #----- Turn data set to a data frame and normalise to 0-1. -----------------------------#
   xy = data.frame(x = x, y = y)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Make sure that 'vertex' is a 4x2 matrix.                                           #
   #---------------------------------------------------------------------------------------#
   #----- Coerce data to a data frame. ----------------------------------------------------#
   if (! is.data.frame(vertex)){
      vertex = try(as.data.frame(vertex),silent=TRUE)
      if ("try-error" %in% vertex){
         stop("vertex cannot be coerced to a data frame.")
      }#end if
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Grab the columns with the x and y coordinates. ----------------------------------#
   if (any(dim(vertex) != c(4,2))){
      stop("vertex must have four rows and two columns!")
   }else{
      names(vertex) = tolower(names(vertex))
      if ("x" %in% names(vertex)){vx=vertex$x}else{vx=vertex[[1]]}
      if ("y" %in% names(vertex)){vy=vertex$y}else{vy=vertex[[2]]}
      vertex = data.frame(x=vx,y=vy)
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Find some auxiliary variables.                                                    #
   #---------------------------------------------------------------------------------------#
   xx = xy$x / xscl
   yy = xy$y / yscl
   #---------------------------------------------------------------------------------------#


   #----- Project the points onto the normalised grid. ------------------------------------#
   ee = ( (1 - xx) * (1 - yy) * vertex$x[1] +      xx  * (1 - yy) * vertex$x[2]
        +      xx  *      yy  * vertex$x[3] + (1 - xx) *      yy  * vertex$x[4] )
   nn = ( (1 - xx) * (1 - yy) * vertex$y[1] +      xx  * (1 - yy) * vertex$y[2]
        +      xx  *      yy  * vertex$y[3] + (1 - xx) *      yy  * vertex$y[4] )
   #---------------------------------------------------------------------------------------#


   #---- Return a data frame with the normalised variables. -------------------------------#
   ans = data.frame(x=ee,y=nn)
   rownames(ans) = rnmout
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end norm.to.coord
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function normalises the X and Y coordinates of quadrilaters.  It doesn't        #
# assume that the shape is a rectangle.                                                    #
#     This function normalises native coordinates (X 0-1; Y 0-1), given the four vertices  #
# of the domain (always in the SW;NW;NE;SE order).  This function doesn't assume that the  #
# shape is a rectangle.                                                                    #
#------------------------------------------------------------------------------------------#
coord.to.norm <<- function(x,y=NULL,vertex,xscl=1,yscl=1){

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


      #------------------------------------------------------------------------------------#
      #    Save row names for later.                                                       #
      #------------------------------------------------------------------------------------#
      rnmout = rownames(xy)
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
   }else{
      #------------------------------------------------------------------------------------#
      #    Save row names for later.                                                       #
      #------------------------------------------------------------------------------------#
      rnmout = names(x)
      #------------------------------------------------------------------------------------#
   }#end if
   #----- Turn data set to a data frame. --------------------------------------------------#
   xy = data.frame(x = x, y = y)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Make sure that 'vertex' is a 4x2 matrix.
   #---------------------------------------------------------------------------------------#
   #----- Coerce data to a data frame. ----------------------------------------------------#
   if (! is.data.frame(vertex)){
      vertex = try(as.data.frame(vertex),silent=TRUE)
      if ("try-error" %in% vertex){
         stop("vertex cannot be coerced to a data frame.")
      }#end if
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Grab the columns with the x and y coordinates. ----------------------------------#
   if (any(dim(vertex) != c(4,2))){
      stop("vertex must have four rows and two columns!")
   }else{
      names(vertex) = tolower(names(vertex))
      if ("x" %in% names(vertex)){vx=vertex$x}else{vx=vertex[[1]]}
      if ("y" %in% names(vertex)){vy=vertex$y}else{vy=vertex[[2]]}
      vertex = data.frame(x=vx,y=vy)
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Find some auxiliary variables.                                                    #
   #---------------------------------------------------------------------------------------#
   et = xy$x         - vertex$x[1]
   eb = vertex$x[4] - vertex$x[1]
   ec = vertex$x[3] - vertex$x[1]
   ed = vertex$x[2] - vertex$x[1]
   e3 = ec - eb - ed
   nt = xy$y         - vertex$y[1]
   nb = vertex$y[4] - vertex$y[1]
   nc = vertex$y[3] - vertex$y[1]
   nd = vertex$y[2] - vertex$y[1]
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
   ans = data.frame( x = xx * xscl
                   , y = yy * yscl
                   )#end data.frame
   rownames(ans) = rnmout
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end coord.to.norm
#==========================================================================================#
#==========================================================================================#
