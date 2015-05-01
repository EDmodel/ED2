#==========================================================================================#
#==========================================================================================#
#      This function appends the patch level either at the beginning or the end of each    #
# cohort.                                                                                  #
#------------------------------------------------------------------------------------------#
append.patch <<- function(ipaco,xpa,xco,left=TRUE){
   xpa = mapply(FUN=list,xpa,SIMPLIFY=TRUE)
   xco = split (x=xco,f=ipaco)
   if (left){
      xpaco = unlist(mapply(FUN=c,xpa,xco,SIMPLIFY=FALSE))
   }else{
      xpaco = unlist(mapply(FUN=c,xco,xpa,SIMPLIFY=FALSE))
   }#end if
   return(xpaco)
}#end append.patch
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      This function finds the total absorbed by each layer.                               #
#------------------------------------------------------------------------------------------#
layer.absorption <<- function(ipaco,use,down,up){
   #----- Break the patches. --------------------------------------------------------------#
   down = split(x=down,f=ipaco)
   up   = split(x=up  ,f=ipaco)
   use  = split(x=use ,f=ipaco)
   ans  = unlist( mapply( FUN      = function(down,up,use){
                                        #----- Find the absorption for 1 layer. -----------#
                                        x      = which(use)
                                        xm1    = x[-1]
                                        xmN    = x[-length(use)]
                                        ans    = 0 * use
                                        ans[x] = c(0,down[xmN]-down[xm1]+up[xm1]-up[xmN])
                                        return(ans)
                                        #--------------------------------------------------#
                                     }#end function
                        , down     = split(x=down,f=ipaco)
                        , up       = split(x=up  ,f=ipaco)
                        , use      = split(x=use ,f=ipaco)
                        , SIMPLIFY = FALSE
                        )#end mapply
                )#end unlist
   #---------------------------------------------------------------------------------------#


   return(ans)
}#end layer.absorption
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      This function finds the total absorbed by each layer.                               #
#------------------------------------------------------------------------------------------#
layer.absorption.one.cohort <<- function(use,down,up){
   #----- We must skip "invisible" cohorts. -----------------------------------------------#
}#end layer.absorption
#==========================================================================================#
#==========================================================================================#
