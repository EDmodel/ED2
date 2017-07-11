#==========================================================================================#
#==========================================================================================#
#     Find elliptical sector radius function.                                              #
#------------------------------------------------------------------------------------------#
elliptical.radius <<- function(x,y,theta,degrees=FALSE){
   if (degrees) theta = theta * pio180

   both.zero = x %==% 0 & y %==% 0
   x.zero    = x %==% 0 & y %!=% 0
   y.zero    = x %!=% 0 & y %==% 0

   x[x.zero] = sqrt(.Machine$double.eps) * y[x.zero]
   y[y.zero] = sqrt(.Machine$double.eps) * x[y.zero]

   ans = ifelse( both.zero
               , 0
               , x*y / sqrt(y*y*(cos(theta))^2+x*x*(sin(theta))^2)
               )#end ifelse

   return(ans)
}#end ellsec.radius
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     Find elliptical sector area function.                                                #
#------------------------------------------------------------------------------------------#
elliptical.area <<- function(x,y,theta0=0,theta1=2*pi,degrees=FALSE){
   if (degrees){
      theta0 = theta0 * pio180
      theta1 = theta1 * pio180
   }#end if

   both.zero = x %==% 0 & y %==% 0
   x.zero    = x %==% 0 & y %!=% 0
   y.zero    = x %!=% 0 & y %==% 0

   x[x.zero] = sqrt(.Machine$double.eps) * y[x.zero]
   y[y.zero] = sqrt(.Machine$double.eps) * x[y.zero]

   fth0 = ifelse( both.zero
                , 0
                , 0.5*x*y*(theta0 - atan((y-x)*sin(2*theta0)/(y+x+(y-x)*cos(2*theta0))))
                )#end ifelse
   fth1 = ifelse( both.zero
                , 0
                , 0.5*x*y*(theta1 - atan((y-x)*sin(2*theta1)/(y+x+(y-x)*cos(2*theta1))))
                )#end ifelse
   ans  = fth1 - fth0
   return(ans)
}#end elliptical.area
#==========================================================================================#
#==========================================================================================#
