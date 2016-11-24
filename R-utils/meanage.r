#==========================================================================================#
#==========================================================================================#
#      This routine computes the mean age for a distribution of ages at steady state with  #
# constant disturbance rate.                                                               #
#------------------------------------------------------------------------------------------#
meanage <<- function(aa=0,az=Inf,lambda){
    if ( aa < 0 ){
       stop (paste( "Variable aa must be non-negative.  You set it to ",aa,"!",sep=""))
    }else if (aa > az){
       stop (paste( "aa (",aa,") cannot be greater than az (",az,")",sep=""))
    }else if (aa == az){
       ans         = aa
    }else if (is.infinite(az)){
       ans         = aa + 1 / lambda
    }else{
       numerator   = ( exp( - lambda * aa ) * ( lambda * aa + 1 )
                     - exp( - lambda * az ) * ( lambda * az + 1 ) )
       denominator = lambda * ( exp( - lambda * aa ) - exp( - lambda * az ) )
       ans         = numerator / denominator
    }#end if
    return(ans)
}#end function meanage
#==========================================================================================#
#==========================================================================================#
