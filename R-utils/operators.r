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


#----- wr is short for "within range", ir is short for inside range (excluding edges). ----#
'%wr%' <<- function(x,y) is.finite(x) & x %>=% min(y,na.rm=TRUE) & x %<=% max(y,na.rm=TRUE)
'%ir%' <<- function(x,y) is.finite(x) & x %>%  min(y,na.rm=TRUE) & x %<%  max(y,na.rm=TRUE)
#------------------------------------------------------------------------------------------#
