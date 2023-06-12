#==========================================================================================#
#==========================================================================================#
#     Function time.to.edge.                                                               #
#                                                                                          #
#     This function computes the time it takes for air to reach the edge, assuming that    #
# air moves along the dominant wind.                                                       #
#                                                                                          #
#------------------------------------------------------------------------------------------#
time.to.edge <<- function( datum
                         , emask
                         , u.vnam
                         , v.vnam
                         , xy.m.fac = 1000.
                         , wmin     = 0.05
                         , tmax     = 5000./wmin
                         , cosmin   = cos(15.*pio180)
                         ){

   #------ Fortran code. ------------------------------------------------------------------#
   fortran.so = file.path(srcdir,"time_to_edge.so")
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #       Extract wind data from data table.                                              #
   #---------------------------------------------------------------------------------------#
   uvsel  = is.finite(datum[[u.vnam]]) & is.finite(datum[[v.vnam]])
   dsel   = datum[uvsel]
   nwind  = nrow(dsel)
   xwind  = dsel$x * xy.m.fac
   ywind  = dsel$y * xy.m.fac
   uwind  = dsel[[u.vnam]]
   vwind  = dsel[[v.vnam]]
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #       Extract wind data from mask table.                                              #
   #---------------------------------------------------------------------------------------#
   edsel = emask$mask & (! is.na(emask$mask) )
   edge  = emask[edsel]
   nedge = nrow(edge)
   xedge = edge$x * xy.m.fac
   yedge = edge$y * xy.m.fac
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #     Find the time to edge.                                                            #
   #---------------------------------------------------------------------------------------#
   if ( (nedge*nwind) > 0){


      #------ Load Fortran and run code. --------------------------------------------------#
      # if (! is.loaded("timetoedge")) dummy = dyn.load(fortran.so)
      # tedge = .Fortran( "timetoedge"
      #                 , nwind = as.integer(nwind)
      #                 , xwind = as.double (xwind)
      #                 , ywind = as.double (ywind)
      #                 , uwind = as.double (uwind)
      #                 , vwind = as.double (vwind)
      #                 , nedge = as.integer(nedge)
      #                 , xedge = as.double (xedge)
      #                 , yedge = as.double (yedge)
      #                 , wmin  = as.double (wmin)
      #                 , time  = as.double (tedge)
      #                )#end .Fortran
      #------------------------------------------------------------------------------------#




      #------ Populate answer. ------------------------------------------------------------#
      # ans        = rep(NA_real_,times=nrow(datum))
      # ans[uvsel] = pmin(tmax,tedge$time)
      # rm(tedge)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Calculate wind speed and direction.                                           #
      #------------------------------------------------------------------------------------#
      SPEED = matrix(data=sqrt(uwind*uwind+vwind*vwind),nrow=nwind,ncol=nedge)
      OMEGA = matrix(data=atan2(y=vwind,x=uwind)       ,nrow=nwind,ncol=nedge)
      #------------------------------------------------------------------------------------#


      #------ Calculate distance from point to edge. --------------------------------------#
      DX    = -outer(X=xwind,Y=xedge,FUN=`-`)
      DY    = -outer(X=ywind,Y=yedge,FUN=`-`)
      DIST  = sqrt(DX*DX+DY*DY)
      DELTA = atan2(y=DY,x=DX)
      #------------------------------------------------------------------------------------#


      #----- Find the baseline time to edge. ----------------------------------------------#
      WEDGE = matrix(data=pmax(wmin,SPEED * cos(OMEGA - DELTA)),nrow=nwind,ncol=nedge)
      #------------------------------------------------------------------------------------#


      #----- Select the minimum time as the time to edge. ---------------------------------#
      TEDGE = DIST / WEDGE
      tedge = apply(X=TEDGE,MARGIN=1,FUN=min,na.rm=TRUE)
      #------------------------------------------------------------------------------------#




      #------ Populate answer. ------------------------------------------------------------#
      ans        = rep(NA_real_,times=nrow(datum))
      ans[uvsel] = pmin(tmax,tedge)
      rm(SPEED,OMEGA,DX,DY,DIST,WEDGE,TEDGE,tedge)
      #------------------------------------------------------------------------------------#

   }else{
      #------ Either wind or edge were completely missing. Return dummy answer. -----------#
      ans = ifelse( test = uvsel
                  , yes  = tmax
                  , no   = NA_real_
                  )#end ifelse
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#



   return(ans)
}#end time.to.edge
#==========================================================================================#
#==========================================================================================#
