#==========================================================================================#
#==========================================================================================#
#      This function appends the patch level either at the beginning or the end of each    #
# cohort.                                                                                  #
#------------------------------------------------------------------------------------------#
append.patch <<- function(ipa,ipaco,xpa,xco,left=TRUE){
   if (length(ipaco) == 0){
      xpaco = xpa
   }else{
      xpa            = mapply(FUN=list,xpa,SIMPLIFY=TRUE)
      names(xpa)     = ipa
      xcofull        = replicate(n=length(ipa),list(vector(length=0)))
      names(xcofull) = ipa
      xco            = split(x=xco,f=ipaco)
      idx            = match(names(xco),names(xcofull))
      xcofull[idx]   = xco
      if (left){
         xpaco = unlist(mapply(FUN=c,xpa,xcofull,SIMPLIFY=FALSE))
      }else{
         xpaco = unlist(mapply(FUN=c,xcofull,xpa,SIMPLIFY=FALSE))
      }#end if
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
