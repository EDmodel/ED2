#==========================================================================================#
#==========================================================================================#
#      This function creates a local polynomial regression fitting for a given population  #
# using R function loess.                                                                  #
#
# p     -- vector with population for each height class
# z     -- vector with height classes
# pfill -- vector with minimum population in case no entries exist. 
# ...   -- additional arguments to be passed to loess
#------------------------------------------------------------------------------------------#
pop.loess <<- function(p,z=seq_along(p),pfill=1e-6,zinf=1.25*max(z),...){

   #----- Data frame used to generate output. ---------------------------------------------#
   nz      = length(z)
   pz.pred = data.frame(z = z)
   #---------------------------------------------------------------------------------------#

   #----- For the prediction step, we only keep entries with actual population. -----------#
   keep = p %gt% 0
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Decide what to do depending on the number of finite points.                       #
   #---------------------------------------------------------------------------------------#
   nkeep = sum(keep)
   if (sum(keep) == 0){
      #----- Impossible to interpolate.  Assume very sparse population at extremes. -------#
      z.fit = c(0.99*min(z),mean(z),zinf)
      p.fit = c(pfill,pfill,pfill)
      #------------------------------------------------------------------------------------#
   }else if (sum(keep) == 1){
      #----- Impossible to interpolate.  Assume very sparse population at extremes. -------#
      z.fit = c(0.99*min(z),z[keep],zinf)
      p.fit = c(p[keep],p[keep],pfill)
      #------------------------------------------------------------------------------------#
   }else{
      #----- Create temporary data frame. -------------------------------------------------#
      z.fit = c(z[keep],zinf)
      p.fit = c(p[keep],pfill)
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Fit the curve and fill remaining NA data with extension. ------------------------#
   pz.fit = data.frame(z=z.fit,p=p.fit,ln.p=log(p.fit))
   lfit = try(loess(ln.p ~ z,data=pz.fit,...),silent=TRUE)
   if ("try-error" %in% is(lfit)) browser()
   ans  = try(predict(lfit,newdata=pz.pred),silent=TRUE)
   if ("try-error" %in% is(ans)) browser()
   ans  = na.fill(exp(ans),fill="extend")
   #---------------------------------------------------------------------------------------#


   #----- Return answer. ------------------------------------------------------------------#
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end pop.loess
#==========================================================================================#
#==========================================================================================#
