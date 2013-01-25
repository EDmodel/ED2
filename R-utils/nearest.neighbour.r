#------------------------------------------------------------------------------------------#
#     This function determines the gridded nearest neighbour for each point of the         #
# network.                                                                                 #
#------------------------------------------------------------------------------------------#
nearest.neighbour = function(netwcoord,gridcoord,g2n=TRUE){
    isok = require(fields)
    if (!isok) stop("The nearest.neighbour function requires the package fields...")

    #--------------------------------------------------------------------------------------#
    #     "Distmat" is the matrix with the distances between the each network point        #
    # (unstructured) and each grid point (structured).                                     #
    #--------------------------------------------------------------------------------------#
    distmat = rdist.earth(x1=netwcoord,x2=gridcoord,miles=FALSE)

    #--------------------------------------------------------------------------------------#
    #    "Nearest" is the index of the closest BRAMS grid point for each point of the      #
    # network.                                                                             #
    #--------------------------------------------------------------------------------------#
    if (g2n){
       nearest = apply(X=distmat,MARGIN=1,FUN=which.min)
    }else{
       nearest = apply(X=distmax,MARGIN=2,FUN=which.min)
    }#end if

    return(nearest)
}#end function
