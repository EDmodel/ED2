#==========================================================================================#
#==========================================================================================#
#     This file contains some useful additional operators.                                 #
#------------------------------------------------------------------------------------------#


#----- Safe logical operators.  These will always return FALSE if x or y are not finite. --#
'%==%' <<- function(x,y){
   if (any(c(FALSE,is.numeric(x) & is.numeric(y)),na.rm=TRUE)){
      ans = is.finite(unlist(x)) & is.finite(unlist(y)) & x == y
   }else{
      ans = ! is.na(x) & ! is.na(y) & ! is.nan(x) & ! is.nan(y) & x == y
   }#end if
   return(ans)
}#end function
'%!=%' <<- function(x,y){
   if (any(c(FALSE,is.numeric(x) & is.numeric(y)),na.rm=TRUE)){
      ans = is.finite(unlist(x)) & is.finite(unlist(y)) & x != y
   }else{
      ans = ! is.na(x) & ! is.na(y) & ! is.nan(x) & ! is.nan(y) & x != y
   }#end if
   return(ans)
}#end function
'%>%'  <<- function(x,y) is.finite(unlist(x)) & is.finite(unlist(y)) & x  > y
'%<%'  <<- function(x,y) is.finite(unlist(x)) & is.finite(unlist(y)) & x  < y
'%>=%' <<- function(x,y) is.finite(unlist(x)) & is.finite(unlist(y)) & x >= y
'%<=%' <<- function(x,y) is.finite(unlist(x)) & is.finite(unlist(y)) & x <= y
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Operators to test whether the values are within range or outside range.              #
#                                                                                          #
# x %wr% y -- TRUE when x is within range of y     (x exactly at bounds returns TRUE )     #
# x %ir% y -- TRUE when x is inside range of y     (x exactly at bounds returns FALSE)     #
# x %or% y -- TRUE when x is outside range of y    (x exactly at bounds returns FALSE)     #
# x %nr% y -- TRUE when x is not inside range of y (x exactly at bounds returns TRUE )     #
#------------------------------------------------------------------------------------------#
'%wr%' <<- function(x,y){
   ans = (! (is.na(x) | is.nan(x))) & ( x >= min(y,na.rm=TRUE) & x <= max(y,na.rm=TRUE) )
   return(ans)
}#end function
'%ir%' <<- function(x,y){
   ans = (! (is.na(x) | is.nan(x))) & ( x >  min(y,na.rm=TRUE) & x <  max(y,na.rm=TRUE) )
   return(ans)
}#end function
'%or%' <<- function(x,y){
   ans = (! (is.na(x) | is.nan(x))) & ( x <  min(y,na.rm=TRUE) | x >  max(y,na.rm=TRUE) )
   return(ans)
}#end function
'%nr%' <<- function(x,y){
   ans = (! (is.na(x) | is.nan(x))) & ( x <= min(y,na.rm=TRUE) | x >= max(y,na.rm=TRUE) )
   return(ans)
}#end function
#------------------------------------------------------------------------------------------#
