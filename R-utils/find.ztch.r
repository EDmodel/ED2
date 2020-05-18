#==========================================================================================#
#==========================================================================================#
#    This function finds the top canopy height and flags which trees are at the canopy     #
# layer.                                                                                   #
#------------------------------------------------------------------------------------------#
find.ztch <<- function(x0,y0,dbh,height,carea,cdepth,ipft,xsw=0,ysw=0,nx,ny,xyres,nretn=30){


   #---------------------------------------------------------------------------------------#
   #    Load Fortran.                                                                      #
   #---------------------------------------------------------------------------------------#
   dyn.load(file.path(srcdir,"find_ztch.so"))
   #---------------------------------------------------------------------------------------#



   #----- Check the input has everything that one needs. ----------------------------------#
   if (missing(x0)   ) stop("Variable \"x0\"  must be provided!"  )
   if (missing(y0)   ) stop("Variable \"y0\"  must be provided!"  )
   if (missing(dbh)  ) stop("Variable \"dbh\" must be provided!"  )
   if (missing(nx)   ) stop("Variable \"nx\"  must be provided!"   )
   if (missing(ny)   ) stop("Variable \"ny\"  must be provided!"   )
   if (missing(xyres)) stop("Variable \"xyres\" must be provided!")
   if (  ( missing(height) || missing(carea) || missing(cdepth)) 
      && ( missing(dbh)    || missing(ipft) )
      ){
      cat0(" Variable \"height\" is missing: ",missing(height),".")
      cat0(" Variable \"carea\"  is missing: ",missing(carea ),".")
      cat0(" Variable \"cdepth\" is missing: ",missing(cdepth),".")
      cat0(" Variable \"dbh\"    is missing: ",missing(dbh)   ,".")
      cat0(" Variable \"ipft\"   is missing: ",missing(ipft)  ,".")
      stop( paste0( "Variables \"dbh\" and \"ipft\" must be provided if any of the "
                  , "following variables are missing: \"height\", \"carea\", \"cdepth\"."
                  )#end paste0
          )#end stop
   }#end if ( ( missing(height) || missing(carea) || missing(cdepth)) && missing(ipft))
   #---------------------------------------------------------------------------------------#



   #----- Fill in missing information. ----------------------------------------------------#
   if (missing(height)) height = dbh2h (ipft=ipft,dbh=dbh)
   if (missing(carea )) carea  = dbh2ca(dbh=dbh,ipft=ipft)
   if (missing(cdepth)) cdepth = (1. - pft$b1Mh[ipft]) * height
   #---------------------------------------------------------------------------------------#



   #----- Find horizontal and vertical radii. ---------------------------------------------#
   rhoriz = sqrt(carea / pi)
   rvert  = 0.5 * cdepth
   #---------------------------------------------------------------------------------------#


   #----- Find dimensions. ----------------------------------------------------------------#
   nind = length(dbh)
   #---------------------------------------------------------------------------------------#



   #----- Find dimensions. ----------------------------------------------------------------#
   nseed = 64
   seed  = sample.int(n=.Machine$integer.max,size=nseed,replace=TRUE)
   #---------------------------------------------------------------------------------------#



   #----- Remove offset so minimum x and y coordinates are always zero. -------------------#
   xrel  = x0 - xsw
   yrel  = y0 - ysw
   #---------------------------------------------------------------------------------------#



   #----- Initialise output variables. ----------------------------------------------------#
   canopy = rep(0, times=nind)
   ztch   = matrix(data=0,nrow=nx,ncol=ny)
   fortout = .Fortran( "find_ztch"
                     , nind   = as.integer(nind  )
                     , nxtch  = as.integer(nx    )
                     , nytch  = as.integer(ny    )
                     , nseed  = as.integer(nseed )
                     , nretn  = as.integer(nretn )
                     , xyres  = as.double (xyres )
                     , seed   = as.integer(seed  )
                     , height = as.double (height)
                     , rh     = as.double (rhoriz)
                     , rv     = as.double (rvert )
                     , x0     = as.double (xrel  )
                     , y0     = as.double (yrel  )
                     , ztch   = as.double (ztch  )
                     , canopy = as.integer(canopy)
                     )#end .Fortran
   #---------------------------------------------------------------------------------------#


   #----- Save variables for output. ------------------------------------------------------#
   ans = list( ztch   = matrix(data=fortout$ztch,nrow=nx,ncol=ny)
             , canopy = as.logical(fortout$canopy)
             )#end list
   #---------------------------------------------------------------------------------------#
}#end function find.ztch
#==========================================================================================#
#==========================================================================================#
