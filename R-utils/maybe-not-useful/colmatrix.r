#==========================================================================================#
#==========================================================================================#
#     This subroutine generates a colour matrix with the equally spaced hues and           #
# saturation for columns, and equally spaced vibrancy for rows.                            #
#------------------------------------------------------------------------------------------#
colmatrix <<- function(n.row,n.col){
   HH = rep( x = seq(from=0.0,to=n.col/(n.col+1),length.out=n.col),each=n.row)
   SS = rep( x = seq(from=1.0,to=0.05           ,length.out=n.row),times=n.col)
   VV = replicate(n=n.col,seq(from=runif(n=1,min=0.5,max=0.8),to=1,length.out=n.row))

   HH = matrix(HH,ncol=n.col,nrow=n.row)
   SS = matrix(SS,ncol=n.col,nrow=n.row)
   VV = matrix(VV,ncol=n.col,nrow=n.row)
   
   ans = matrix(hsv(h=HH,s=SS,v=VV),ncol=n.col,nrow=n.row)
   return(ans)
}#end function
#==========================================================================================#
#==========================================================================================#
