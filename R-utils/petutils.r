#==========================================================================================#
#==========================================================================================#
#     Potential evapotranspiration, following Thornthwaite (1948).  This method assumes    #
# monthly data as inputs.                                                                  #
#                                                                                          #
# Reference:                                                                               #
#                                                                                          #
# Thornthwaite CW. 1948. An approach toward a rational classification of climate.          #
#    Geogr. Rev. 38: 55-94. doi:10.2307/210739.                                            #
#------------------------------------------------------------------------------------------#
pet.thornth <<- function(tempk,lon,lat,mon,year){

   #----- Find the mean declination and day length for each month. ------------------------#
   when0   = ( chron("12/31/2009")                    # Zero
             + rep(sequence(365)  ,each =288)         # Day of year
             + rep(sequence(288)-1,times=365) / 288 ) # Five-minute
   mon0    = nummonths(when0)
   year0   = numyears (when0)
   zen0    = ed.zen(lon,lat,when0)
   declin0 = tapply( X     = zen0$declin
                   , INDEX = mon0
                   , FUN   = mean
                   , na.rm = TRUE
                   )#end tapply
   sun0    = 24. * tapply( X     = zen0$day
                         , INDEX = mon0
                         , FUN   = function(x) sum(x,na.rm=TRUE) / length(x)
                         )#end tapply
   when    = chron(paste(mon,15,year,sep="/"))
   #---------------------------------------------------------------------------------------#



   #----- Expand the declination and day length for the entire time series. ---------------#
   declin  = declin0[mon]
   sun     = sun0   [mon]
   dmax    = daymax(month=mon,year=year)
   #---------------------------------------------------------------------------------------#



   #----- Calculate the heat index for each month. ----------------------------------------#
   tempc = 0. * tempk + pmax(0.,tempk - t00)
   hinst = (tempc/5)^1.514
   havg  = tapply(X=hinst,INDEX=mon,FUN=mean)
   #---------------------------------------------------------------------------------------#


   #----- The running average heat index. -------------------------------------------------#
   hindex    = rowSums(embed(x=c(havg[-1],hinst),12))
   #---------------------------------------------------------------------------------------#


   #----- Correction coefficient. ---------------------------------------------------------#
   kcoeff = sun / 12. * dmax / 30.
   #---------------------------------------------------------------------------------------#


   #----- Exponent coefficient. -----------------------------------------------------------#
   mcoeff = 4.92e-1 + hindex * ( 1.79e-2 + hindex * (-7.71e-5 + hindex * 6.75e-7 ) )
   #---------------------------------------------------------------------------------------#


   #----- Exponent coefficient. -----------------------------------------------------------#
   ans = 16. * kcoeff * (10. * tempc / hindex) ^ mcoeff
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end function pet.thornth
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     Potential evapotranspiration, following Penman-Monteith approach (Monteith, 1965).   #
# This method assumes monthly data as inputs.  It has the option to use the tabulated      #
# values for stomatal conductance by Maes et al. (2019).                                   #
#                                                                                          #
# Reference:                                                                               #
#                                                                                          #
# Monteith JL. 1965. Evaporation and environment. Symp. Soc. Exp. Biol. 19: 205-234.       #
#                                                                                          #
# Maes WH, Gentine P, Verhoest NEC , Miralles DG. 2019. Potential evaporation at eddy-     #
#    covariance sites across the globe. Hydrol. Earth Syst. Sci., 23: 925-948.             #
#    doi:10.5194/hess-23-925-2019.                                                         #
#------------------------------------------------------------------------------------------#
pet.penmon <<- function( rnet
                       , fgnd = 0.
                       , prss
                       , temp
                       , shv
                       , vels
                       , href
                       , hgt
                       , rc    = NA_real_
                       , biome = c("ebf","dbf","enf","mxf","csh","wsa"
                                  ,"osh","sav","wet","cro","gra")
                       ){

   #----- Estimate aerodynamic resistance. ------------------------------------------------#
   zd  = 0.63 * hgt
   z0m = 0.13 * hgt
   z0h = 0.10 * z0m
   ra  = log((href-zd)/z0m) * log((href-zd)/z0h) / (vonk * vonk * vels)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Canopy resistance: in case it is not provided, use the tabulated values from Maes  #
   # et al. (2019)                                                                         #
   #---------------------------------------------------------------------------------------#
   if (is.na(rc)){
      biome = match.arg(biome)
      rc    = switch( EXPR = biome
                    , ebf  =  23.80952
                    , dbf  =  30.67485
                    , enf  =  35.21127
                    , mxf  = 100.00000
                    , csh  = 117.64706
                    , wsa  = 119.04762
                    , osh  = 128.20513
                    , sav  = 232.55814
                    , wet  =  50.00000
                    , cro  =  26.10966
                    , gra  =  32.78689
                    , 81.4591
                    )#end switch
   }#end if (is.na(rc))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Calculate air density and vapour pressure deficit.                                #
   #---------------------------------------------------------------------------------------#
   rhos  = idealdenssh(pres=prss,temp=temp,qvpr=shv,qtot=shv)
   vpdef = eslif(temp) - prss * shv / (ep + (1.-ep)*shv)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Calculate the specific latent heat of vaporisation, the psychometric constant and #
   # the slope of saturation vapour pressure curve.                                        #
   #---------------------------------------------------------------------------------------#
   splat  = alvli(temp)
   psycho = cpdry * prss / (ep * splat)
   eslope = eslifp(temp)
   #---------------------------------------------------------------------------------------#


   #----- Find the potential ET (kgW/m2/s). -----------------------------------------------#
   anom = ( eslope * (rnet - fgnd) + rhos * cpdry * vpdef / ra )
   aden = splat * (eslope + psycho * (1. + rc / ra))
   ans  = ifelse( test = aden %!=% 0., yes = anom / aden, no = NA_real_)
   #---------------------------------------------------------------------------------------#
   return(ans)
}#end function pet.penmon
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     Potential evapotranspiration, following Priestley and Taylor (1972).  This method    #
# assumes monthly data as inputs.  It has the option to use tabulated values for the alpha #
# parameter from Maes et al. (2019).                                                       #
#                                                                                          #
# References:                                                                              #
#                                                                                          #
# Priestley CHB , Taylor RJ. 1972. On the assessment of surface heat flux and evaporation  #
#    using large-scale parameters. Mon. Wea. Rev., 100: 81-92.                             #
#    doi:10.1175/1520-0493(1972)100<0081:OTAOSH>2.3.CO;2.                                  #
#                                                                                          #
# Maes WH, Gentine P, Verhoest NEC , Miralles DG. 2019. Potential evaporation at eddy-     #
#    covariance sites across the globe. Hydrol. Earth Syst. Sci., 23: 925-948.             #
#    doi:10.5194/hess-23-925-2019.                                                         #
#------------------------------------------------------------------------------------------#
pet.prtay  <<- function( rnet
                       , fgnd = 0.
                       , prss
                       , temp
                       , alpha = NA_real_
                       , biome = c("ebf","dbf","enf","mxf","csh","wsa"
                                  ,"osh","sav","wet","cro","gra")
                       ){



   #---------------------------------------------------------------------------------------#
   #    The alpha term. In case it is not provided, use the tabulated values from Maes     #
   # et al. (2019)                                                                         #
   #---------------------------------------------------------------------------------------#
   if (is.na(alpha)){
      biome = match.arg(biome)
      alpha = switch( EXPR = biome
                    , ebf  = 1.09
                    , dbf  = 1.09
                    , enf  = 0.89
                    , mxf  = 0.88
                    , csh  = 0.90
                    , wsa  = 0.95
                    , osh  = 0.87
                    , sav  = 0.79
                    , wet  = 1.03
                    , cro  = 1.15
                    , gra  = 1.02
                    , 0.97
                    )#end switch
   }#end if (is.na(rc))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Calculate the specific latent heat of vaporisation, the psychometric constant and #
   # the slope of saturation vapour pressure curve.                                        #
   #---------------------------------------------------------------------------------------#
   splat  = alvli(temp)
   psycho = cpdry * prss / (ep * splat)
   eslope = eslifp(temp)
   #---------------------------------------------------------------------------------------#


   #----- Find the potential ET (kgW/m2/s). -----------------------------------------------#
   anom = alpha * eslope * ( rnet - fgnd )
   aden = splat * (eslope + psycho)
   ans  = ifelse( test = aden %!=% 0., yes = anom / aden, no = NA_real_)
   #---------------------------------------------------------------------------------------#
   return(ans)
}#end function pet.prtay
#==========================================================================================#
#==========================================================================================#
