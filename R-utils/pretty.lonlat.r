#==========================================================================================#
#==========================================================================================#
#     Function that creates a pretty name for longitude and latitude.                      #
#------------------------------------------------------------------------------------------#
pretty.lonlat <<- function(x,type="longitude",zero.ne=TRUE,...){
   s.type = tolower(substring(type,1,2))
   if (! s.type %in% c("lo","la")){
      cat(" type = ",type,"\n")
      stop(" Variable 'type' must be either 'longitude' or 'latitude'...")
   }#end if

   #----- Use pretty to find the typical tick marks. --------------------------------------#
   px  = pretty(range(x,na.rm=TRUE),...)
   npx = length(px)
   #---------------------------------------------------------------------------------------#



   #----- Decide whether to make 0 part of North (East) or South (West). ------------------#
   if (zero.ne){
      ne = px >= 0
   }else{
      ne = px > 0
   }#end if
   if (s.type == "lo"){
      flag = ifelse(ne,"E","W")
   }else{
      flag = ifelse(ne,"N","S")
   }#end if
   
   canvas = rep("xxx*degree*yyy",times=npx)
   canvas = mapply(FUN=sub,replacement=abs(px),x=canvas,MoreArgs=list(pattern="xxx"))
   canvas = mapply(FUN=sub,replacement=flag   ,x=canvas,MoreArgs=list(pattern="yyy"))
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Save the output list.                                                              #
   #---------------------------------------------------------------------------------------#
   ans = list(at=px,labels=parse(text=canvas))
   #---------------------------------------------------------------------------------------#
   return(ans)
}#end pretty.lonlat
#==========================================================================================#
#==========================================================================================#
