#==========================================================================================#
#==========================================================================================#
#      This function creates a local polynomial regression fitting for a given population  #
# using either loess (local fitting) or piecewise cubic spline.                            #
#                                                                                          #
# p      -- vector with population for each height class                                   #
# z      -- vector with height classes                                                     #
# pfill  -- vector with minimum population in case no entries exist.                       #
# method -- method to interpolate (loess filter or spline.                                 #
# ...    -- additional arguments to be passed to loess or splinefun                        #
#------------------------------------------------------------------------------------------#
pop.interp <<- function( p
                       , z      = seq_along(p)
                       , zout   = z
                       , zzero  = 0.
                       , pzero  = sqrt(0.1)
                       , pinf   = 1.e-6
                       , zinf   = 1.25*max(z)
                       , pfill  = sqrt(pinf*pzero)
                       , interp = c("loess","spline")
                       , ...
                       ){

   #----- Standardise method to interpolate. ----------------------------------------------#
   interp = match.arg(interp)
   #---------------------------------------------------------------------------------------#



   #----- Data frame used to generate output. ---------------------------------------------#
   nz      = length(z)
   pz.pred = data.frame(z = zout)
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
   }else{
      #------------------------------------------------------------------------------------#
      #     Create temporary data frame.  Add some boundary conditions at extremes to      #
      # avoid weird population profiles near the edges.                                    #
      #------------------------------------------------------------------------------------#
      z.fit = c(zzero,z[keep],zinf)
      p.fit = c(pzero,p[keep],pinf)
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Fit the curve and fill remaining NA data with extension.                         #
   #---------------------------------------------------------------------------------------#
   pz.fit = data.frame(z=z.fit,p=p.fit,ln.p=log(p.fit))
   if (interp %in% "loess"){
      lfit = try(loess(ln.p ~ z,data=pz.fit,...),silent=TRUE)
      if ("try-error" %in% is(lfit)) browser()
      ans  = try(predict(lfit,newdata=pz.pred),silent=TRUE)
      if ("try-error" %in% is(ans)) browser()
   }else{
      lspl = try(splinefun(x=pz.fit$z,y=pz.fit$ln.p,...))
      if ("try-error" %in% is(lspl)) browser()
      ans  = try(lspl(x=zout),silent=TRUE)
      if ("try-error" %in% is(ans)) browser()
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Fill in remaining NAs. ----------------------------------------------------------#
   ans  = na.fill(exp(ans),fill="extend")
   #---------------------------------------------------------------------------------------#


   #----- Return answer. ------------------------------------------------------------------#
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end pop.interp
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      This function creates a local polynomial regression fitting for a given population  #
# fraction using either loess (local fitting) or piecewise cubic spline.                   #
#                                                                                          #
# p      -- vector with population for each height class                                   #
# z      -- vector with height classes                                                     #
# pfill  -- vector with minimum population in case no entries exist.                       #
# method -- method to interpolate (loess filter or spline.                                 #
# ...    -- additional arguments to be passed to loess or splinefun                        #
#------------------------------------------------------------------------------------------#
frac.interp <<- function( f
                        , ffill
                        , z      = seq_along(p)
                        , zout   = z
                        , zzero  = 0.
                        , zinf   = 1.25*max(z)
                        , interp = c("loess","spline")
                        , ...
                        ){

   #----- Standardise method to interpolate. ----------------------------------------------#
   interp = match.arg(interp)
   #---------------------------------------------------------------------------------------#


   #----- Make sure required arguments are provided. --------------------------------------#
   if (missing(f) || missing(ffill)){
      cat0(" Required arguments missing!")
      cat0(" - Variable 'f' is missing:     ",missing(f    ))
      cat0(" - Variable 'ffill' is missing: ",missing(ffill))
      stop(" Please provide missing arguments")
   }#end if (missing(f) || missing(ffill))
   #---------------------------------------------------------------------------------------#



   #----- Data frame used to generate output. ---------------------------------------------#
   nz      = length(z)
   fz.pred = data.frame(z = zout)
   #---------------------------------------------------------------------------------------#

   #----- For the prediction step, we only keep entries with actual population. -----------#
   f       = ifelse(is.finite(f),f,0)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Decide what to do depending on the number of finite points.                       #
   #---------------------------------------------------------------------------------------#
   nkeep = sum(keep)
   if (sum(keep) == 0){
      #----- Impossible to interpolate.  Assume very sparse population at extremes. -------#
      z.fit = c(0.99*min(z),mean(z),zinf)
      f.fit = c(ffill,ffill,ffill)
      #------------------------------------------------------------------------------------#
   }else{
      #------------------------------------------------------------------------------------#
      #     Create temporary data frame.  Add some boundary conditions at extremes to      #
      # avoid weird population profiles near the edges.                                    #
      #------------------------------------------------------------------------------------#
      z.fit = c(zzero,z[keep],zinf )
      f.fit = c(ffill,f[keep],ffill)
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Fit the curve and fill remaining NA data with extension.                         #
   #---------------------------------------------------------------------------------------#
   fz.fit = data.frame(z=z.fit,f=f.fit)
   if (interp %in% "loess"){
      lfit = try(loess(f ~ z,data=fz.fit,...),silent=TRUE)
      if ("try-error" %in% is(lfit)) browser()
      ans  = try(predict(lfit,newdata=fz.pred),silent=TRUE)
      if ("try-error" %in% is(ans)) browser()
   }else{
      lspl = try(splinefun(x=fz.fit$z,y=fz.fit$f,...))
      if ("try-error" %in% is(lspl)) browser()
      ans  = try(lspl(x=zout),silent=TRUE)
      if ("try-error" %in% is(ans)) browser()
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Fill in remaining NAs. ----------------------------------------------------------#
   ans  = na.fill(ans,fill="extend")
   ans  = pmin(1,pmax(0,ans)) + 0. * ans
   #---------------------------------------------------------------------------------------#


   #----- Return answer. ------------------------------------------------------------------#
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end pop.interp
#==========================================================================================#
#==========================================================================================#
