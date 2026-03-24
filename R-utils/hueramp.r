#------------------------------------------------------------------------------------------#
#   This routine creates a colour ramp based on different hues.                            #
#------------------------------------------------------------------------------------------#
hueramp <<- function(n,hue=0,alpha=1.0){
   sss   = c( 10,30,50,70,90) / 100.
   vvv   = c( 90,70,70,50,30) / 100.
   HHH   = rep(hue,times=length(sss))
   hhh   = HHH / 360.
   pivot = round(seq(from=1,to=n,by=(n-1)/(length(hhh)-1)),digits=0)

   hout  = pmax(0,pmin(1,spline(x=pivot,y=hhh,n=n)$y))
   sout  = pmax(0,pmin(1,spline(x=pivot,y=sss,n=n)$y))
   vout  = pmax(0,pmin(1,spline(x=pivot,y=vvv,n=n)$y))

   ans   = hsv(h=hout,s=sout,v=vout,alpha=alpha)
}#end function hueramp
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#   Inverse hue ramp #
#------------------------------------------------------------------------------------------#
ihueramp <<- function(...) rev(hueramp(...))
#------------------------------------------------------------------------------------------#
