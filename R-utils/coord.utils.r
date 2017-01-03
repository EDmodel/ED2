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
   rownames(domain) = c("00","0Y","XY","X0")
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Decide the vertex labels based on the closest distance to the "ideal" vertices,   #
   # and label them according to the minimum distance to the perfect vertices.             #
   #---------------------------------------------------------------------------------------#
   idx.00 = which.min(sqrt((vertex$x-domain$x[1])^2+(vertex$y-domain$y[1])^2))
   idx.0Y = which.min(sqrt((vertex$x-domain$x[2])^2+(vertex$y-domain$y[2])^2))
   idx.XY = which.min(sqrt((vertex$x-domain$x[3])^2+(vertex$y-domain$y[3])^2))
   idx.X0 = which.min(sqrt((vertex$x-domain$x[4])^2+(vertex$y-domain$y[4])^2))
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Re-order the vertices so the vertices always go clockwise, starting from the      #
   # X=0,Y=0 vertex.                                                                       #
   #---------------------------------------------------------------------------------------#
   vertex           = vertex[c(idx.00,idx.0Y,idx.XY,idx.X0),]
   rownames(vertex) = c("00","Y0","XY","X0")
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
   ee = ( (1 - xx) * (1 - yy) * vertex$x[1] + (1 - xx) *      yy  * vertex$x[2]
        +      xx  *      yy  * vertex$x[3] +      xx  * (1 - yy) * vertex$x[4] )
   nn = ( (1 - xx) * (1 - yy) * vertex$y[1] + (1 - xx) *      yy  * vertex$y[2]
        +      xx  *      yy  * vertex$y[3] +      xx  * (1 - yy) * vertex$y[4] )
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
# the four vertices of the domain (always clockwise, first vertex corresponding to the     #
# X=0;Y=0 vertex).  This function doesn't assume that the shape is a rectangle.            #
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
   ee = ( (1 - xx) * (1 - yy) * vertex$x[1] + (1 - xx) *      yy  * vertex$x[2]
        +      xx  *      yy  * vertex$x[3] +      xx  * (1 - yy) * vertex$x[4] )
   nn = ( (1 - xx) * (1 - yy) * vertex$y[1] + (1 - xx) *      yy  * vertex$y[2]
        +      xx  *      yy  * vertex$y[3] +      xx  * (1 - yy) * vertex$y[4] )
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
# of the domain (always in the clockwise order, first vertex must be the X=0,Y=0 one).     #
# This function doesn't assume that the shape is a rectangle.                              #
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
   et = xy$x        - vertex$x[1]
   e2 = vertex$x[2] - vertex$x[1]
   e3 = vertex$x[3] - vertex$x[2] - vertex$x[4] + vertex$x[1]
   e4 = vertex$x[4] - vertex$x[1]
   nt = xy$y        - vertex$y[1]
   n2 = vertex$y[2] - vertex$y[1]
   n3 = vertex$y[3] - vertex$y[2] - vertex$y[4] + vertex$y[1]
   n4 = vertex$y[4] - vertex$y[1]
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    First, find out whether this equation is quadratic or not.                         #
   #---------------------------------------------------------------------------------------#
   aa         = (e3 * n2 - e2 * n3)
   if (aa == 0.){
      #----- Not quadratic. ---------------------------------------------------------------#
      yy = - ( et * n4 - e4 * nt ) / ( e4 * n2 - e2 * n4 + et * n3 - e3 * nt )
      xx = (et - yy * e2) / (e4 + yy * e3)
      #------------------------------------------------------------------------------------#
   }else{
      #----- Quadratic, search for the two solutions. -------------------------------------#
      bb         = ( e4 * n2 - e2 * n4 + et * n3 - e3 * nt ) / aa
      cc         = ( et * n4 - e4 * nt ) / aa
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Find the two possible y coordinates and the associated x coordinates (in case  #
      # this solution is quadratic).                                                       #
      #------------------------------------------------------------------------------------#
      dd = bb * bb - 4 * cc
      y1 = 0.5 * ( - bb - sqrt(dd) )
      y2 = 0.5 * ( - bb + sqrt(dd) )
      x1 = (et - y1 * e2) / (e4 + y1 * e3)
      x2 = (et - y2 * e2) / (e4 + y2 * e3)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Pick the one that makes the most sense.                                        #
      #------------------------------------------------------------------------------------#
      one   = (x1 %wr% c(0.0,1.0)) & (y1 %wr% c(0.0,1.0))
      two   = (x2 %wr% c(0.0,1.0)) & (y2 %wr% c(0.0,1.0))
      three = (x1^2+y1^2) < (x2^2+y2^2)
      xx  = ifelse(one,x1,ifelse(two,x2,ifelse(three,x1,x2)))
      yy  = ifelse(one,y1,ifelse(two,y2,ifelse(three,y1,y2)))
      #------------------------------------------------------------------------------------#
   }#end if (aa == 0.)
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






#==========================================================================================#
#==========================================================================================#
#      This function converts decimal degrees to degrees,minutes,seconds.                  #
#------------------------------------------------------------------------------------------#
dec2dms <<- function(lon=NULL,lat=NULL){

   #----- Transform longitude into degrees, minutes, seconds. -----------------------------#
   if (! is.null(lon)){
      lon    = (lon + 180.) %% 360 - 180.
      degree = sprintf("%3i"  ,floor(abs(lon))                      )
      minute = sprintf("%2.2i",floor(abs(lon) %% 1 * 60)            )
      second = sprintf("%2.2i",floor((abs(lon) %% 1 * 60) %% 1 * 60))
      hemisf = ifelse(lon %>=% 0,"E","W")
      olon   = paste0(degree,"-",minute,"\'",second,"\"",hemisf)
      olon   = ifelse(is.finite(lon),olon,NA_character_)
   }else{
      olon = NULL
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Transform latitude into degrees, minutes, seconds. ------------------------------#
   if (! is.null(lat)){
      degree = sprintf("%2i"  ,floor(abs(lat))                      )
      minute = sprintf("%2.2i",floor(abs(lat) %% 1 * 60)            )
      second = sprintf("%2.2i",floor((abs(lat) %% 1 * 60) %% 1 * 60))
      hemisf = ifelse(lat %>=% 0,"N","S")
      olat   = paste0(degree,"-",minute,"\'",second,"\"",hemisf)
      olat   = ifelse(is.finite(lat),olat,NA_character_)
   }else{
      olat = NULL
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Return result.  We must check whether lon or lat have been provided.              #
   #---------------------------------------------------------------------------------------#
   if (is.null(lon) && is.null(lat)){
      ans           = NULL
   }else if (is.null(lon)){
      ans           = olat
      names(ans)    = names(lat)
   }else if (is.null(lat)){
      ans           = olon
      names(ans)    = names(lon)
   }else{
      ans           = data.frame(lon=olon,lat=olat,stringsAsFactors=FALSE)
      rownames(ans) = names(lon)
   }#end if (is.null(lon) && is.null(lat))
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end dec2dms
#------------------------------------------------------------------------------------------#

