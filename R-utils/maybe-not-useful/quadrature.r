#==========================================================================================#
#==========================================================================================#
#     This function evaluates the integral of a function using simple quadrature methods.  #
# Instead of finding the full integral, it finds the integral along all from the first to  #
# the last point.                                                                          #
#                                                                                          #
#   Input:                                                                                 #
#   ~ y      -- function evaluations at each x point                                       #
#   ~ x      -- the x points                                                               #
#   ~ method -- how to estimate the integral (case insensitive and partial match).         #
#               Possible values are:                                                       #
#               "trapezoid" - linear interpolation.                                        #
#               "spline"    - use spline then integrate.                                   #
#   ~ cum    -- return cumulative sum?                                                     #
#------------------------------------------------------------------------------------------#
quadrature <<- function(y,x,method="trapezoid",subdivisions=10000,cum=FALSE){

    #----- Decide which method to use. ----------------------------------------------------#
    if (length(y) != length(x)){
       cat (" Length of y: ",length(y),"\n")
       cat (" Length of x: ",length(x),"\n")
       stop(" X and Y must have the same size...")
    }else if (any(is.na(x))){
       stop(" X cannot contain NA...")
    }else if (any(is.na(y))){
       if (cum){
          cumint = NA * x
       }else{
          cumint = NA
       }#end if
       return(cumint)
    }#end if
    #--------------------------------------------------------------------------------------#


    #--------------------------------------------------------------------------------------#
    #      Find the vector length.                                                         #
    #--------------------------------------------------------------------------------------#
    nxy = length(x)
    #--------------------------------------------------------------------------------------#


    #----- Decide which method to use. ----------------------------------------------------#
    if (substr(tolower(method),1,1) == "t"){
       #----- Trapezoid method. -----------------------------------------------------------#
       ybar   = c(0,0.5*(y[-nxy]+y[-1]))
       dx     = c(0,diff(x))
       if (cum){
          cumint = cumsum(ybar*dx)
       }else{
          cumint = sum(ybar*dx)
       }#end if
       #-----------------------------------------------------------------------------------#
    }else if (substr(tolower(method),1,1) == "s"){
       smoothy = splinefun(x,y,method="monoH.FC")
       if (cum){
          xuse = x
       }else{
          xuse = x[nxy]
       }#end if
       intans  = mapply( FUN      = integrate
                       , upper    = xuse
                       , MoreArgs = list( f             = smoothy
                                        , lower         = x[1]
                                        , subdivisions  = subdivisions
                                        , stop.on.error = FALSE
                                        )#end list
                       )#end mapply
       if ( all(unlist(intans["message",]) == "OK") ){
          cumint  = unlist(intans["value",] )
       }else{
          warning("Spline method failed, using trapezoid instead")
          cumint  = quadrature(y,m,method="trapezoid",subdivisions=subdivisions,cum=cum)
       }#end if
    }else{
       stop(paste("Cannot recognise the method : ",method,"...",sep=" "))
    }#end if
    #--------------------------------------------------------------------------------------#
    return(cumint)
}#end function quadrature
#==========================================================================================#
#==========================================================================================#
