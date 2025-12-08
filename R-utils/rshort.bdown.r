 #==========================================================================================#
#==========================================================================================#
#      This subroutine computes the split between direct and diffuse radiation, and        #
# between visible and near-infrared radiation.  Three methods are provided:                #
#                                                                                          #
# 1.  Default is Weiss-Norman:                                                             #
#     Weiss, A., J. M. Norman, 1985: Partitioning solar radiation into direct and diffuse, #
#        visible and near-infrared components.  Agric. For. Meteorol., 34, 205-213. (WN85) #
#                                                                                          #
# 2.  Clearness Index (clearidx).  This is a combination of methods suggested by:          #
#     Boland, J., B. Ridley, B. Brown, 2008: Models of diffuse solar radiation.  Renew.    #
#        Energy, 33, 575-584. doi:10.1016/j.renene.2007.04.012 (BD08).                     #
#                                                                                          #
#     Tsubo, M., S. Walker, 2005: Relationships between photosynthetically active          #
#        radiation and clearness index at Bloemfontein, South Africa.  Theor. Appl.        #
#         Climatol., 80, 17-25.  doi:10.1007/s00704-004-0080-5 (TW05).                     #
#                                                                                          #
#     Bendix, J., B. Silva, K. Roos, D. O. Gottlicher, R. Rollenbeck, T. Nauss, E. Beck,   #
#        2010: Model parameterization to simulate and compare the PAR absorption potential #
#        of two competing plant species.  Int. J. Biometeorol. 54, 283-295.                #
#        doi:10.1007/s00484-009-0279-3 (BX10).                                             #
#                                                                                          #
# Input variables:                                                                         #
#                                                                                          #
#    * rad.in   - The incoming radiation at surface, in W/m2.  This can be either PAR,     #
#                 NIR, or the total shortwave radiation, but it must be in W/m2 in any of  #
#                 the cases.                                                               #
#    * atm.prss - The atmospheric pressure at the surface, in Pa.  An actual measurement   #
#                 is better, but if you don't have one just use some typical value given   #
#                 the place altitude (higher elevation sites get more radiation).          #
#    * cosz     - The cosine of zenith angle.  This can be estimated using function ed.zen #
#                 in file zen.r                                                            #
#    * rad.type - The type of radiation provided in rad.in.  Default is total shortwave    #
#                 radiation, but the function also accepts PAR or NIR.  The value is case  #
#                 insensitive and only the first letter is checked.  "p" means PAR, "n"    #
#                 means NIR, and any other letter will be assumed shortwave.               #
#------------------------------------------------------------------------------------------#
rshort.bdown <<- function(rad.in,atm.prss,cosz,rad.type=c("rshort","par","nir")
                         ,rad.method=c("wn85","clearidx","sib"),apply.bx10.corr=FALSE){

   #----- Standard method. ----------------------------------------------------------------#
   rad.method = match.arg(rad.method)
   rad.type   = match.arg(rad.type  )
   #---------------------------------------------------------------------------------------#


   if (length(atm.prss) == 1) atm.prss = rep(atm.prss,times=length(rad.in))

   #----- Prevent clearness index method to be used when rad.type is not rshort. ----------#
   if ( (rad.method %in% "sib") && (! (rad.type %in% "rshort")) ){
      cat0(" - Radiation method: ",rad.method,".")
      cat0(" - Input radiation type: ",rad.type,".")
      stop(" Radiation input type must be \"rshort\" if using method \"sib\".)")
   }#end if
   #---------------------------------------------------------------------------------------#


   ans = switch( EXPR     = rad.method
               , wn85     = rshort.bdown.weissnorman(rad.in,atm.prss,cosz,rad.type)
               , clearidx = rshort.bdown.clearidx   (rad.in,atm.prss,cosz,apply.bx10.corr
                                                    ,rad.type)
               , sib      = rshort.bdown.sib        (rad.in,atm.prss,cosz)
               )#end switch
   return(ans)
}#end rshort.bdown
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#      This subroutine computes the split between direct and diffuse radiation, and        #
# between visible and near-infrared radiation, using the WN85 method.                      #
#------------------------------------------------------------------------------------------#
rshort.bdown.weissnorman <<- function(rad.in,atm.prss,cosz
                                     ,rad.type=c("rshort","par","nir")){

   #----- Standard method. ----------------------------------------------------------------#
   rad.type   = match.arg(rad.type  )
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #    Local constants.                                                                   #
   #---------------------------------------------------------------------------------------#
   #----- Extinction coefficient. (equations 1 and 4 of WN85) -----------------------------#
   par.beam.expext  = -0.185
   nir.beam.expext  = -0.060
   #----- This is the typical conversion of diffuse radiation in sunny days. --------------#
   par2diff.sun = 0.400
   nir2diff.sun = 0.600
   #----- Coefficients for various equations in WN85. -------------------------------------#
   wn85.06 = c( -1.1950, 0.4459, -0.0345 )
   wn85.11 = c(    0.90, 0.70  )
   wn85.12 = c(    0.88, 0.68  )
   #---------------------------------------------------------------------------------------#


   #------ Initialise the radiation with NAs. ---------------------------------------------#
   par.beam    = NA_real_ * rad.in
   nir.beam    = NA_real_ * rad.in
   par.diff    = NA_real_ * rad.in
   nir.diff    = NA_real_ * rad.in
   par.full    = NA_real_ * rad.in
   nir.full    = NA_real_ * rad.in
   rshort.beam = NA_real_ * rad.in
   rshort.diff = NA_real_ * rad.in
   rshort.full = NA_real_ * rad.in
   par.max     = NA_real_ * rad.in
   nir.max     = NA_real_ * rad.in
   rshort.max  = NA_real_ * rad.in
  #---------------------------------------------------------------------------------------#


   #------ Make day and night flags. ------------------------------------------------------#
   ntimes   = length(cosz)
   day      = cosz %gt% cosz.min
   night    = cosz %lt% cosz.twilight
   twilight = (! day) & (! night)
   #---------------------------------------------------------------------------------------#


   #----- Save the Chapman function to account for the earth's curvature. -----------------#
   chapman      = pmin(lnexp.max,huestis.fun(cosz=cosz))
   eff.cosz     = ifelse( test = night, yes = 0., no = 1./chapman)
   log10chapman = log10(chapman)
   #---------------------------------------------------------------------------------------#


   #----- Total radiation at the top of the atmosphere [  W/m2], using ED defaults. -------#
   rshort.beam.toa = solar * eff.cosz
   par.beam.toa    = fvis.beam.def * rshort.beam.toa
   nir.beam.toa    = fnir.beam.def * rshort.beam.toa
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Find the potential PAR components (beam, diffuse, total), using equations 1, 3,    #
   # and 9 of WN85.                                                                        #
   #---------------------------------------------------------------------------------------#
   par.beam.pot = ( par.beam.toa
                  * exp ( par.beam.expext * (atm.prss / prefsea) * chapman) )
   par.diff.pot = par2diff.sun * (par.beam.toa - par.beam.pot)
   par.full.pot = par.beam.pot + par.diff.pot
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the NIR absorption of 10 mm of precipitable water, using WN85 equation 6.    #
   #---------------------------------------------------------------------------------------#
   w10 = ( rshort.beam.toa
         * 10 ** ((wn85.06[1]) + log10chapman * (wn85.06[2] + wn85.06[3] * log10chapman))
         )#end w10
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the potential direct and diffuse near-infrared radiation, using equations    #
   # 4, 5, and 10 of WN85.                                                                 #
   #---------------------------------------------------------------------------------------#
   nir.beam.pot = ( ( nir.beam.toa
                    * exp ( nir.beam.expext * (atm.prss / prefsea) * chapman) - w10 ) )
   nir.diff.pot = nir2diff.sun * ( nir.beam.toa - nir.beam.pot - w10 )
   nir.full.pot = nir.beam.pot + nir.diff.pot
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    In case of twilight, set beam radiation to zero and diffuse radiation to           #
   # potential.                                                                            #
   #---------------------------------------------------------------------------------------#
   par.beam.pot[twilight] = 0.0
   par.diff.pot[twilight] = par.full.pot[twilight]
   nir.beam.pot[twilight] = 0.0
   nir.diff.pot[twilight] = nir.full.pot[twilight]
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Total maximum radiation.                                                          #
   #---------------------------------------------------------------------------------------#
   par.max    = par.full.pot
   nir.max    = nir.full.pot
   rshort.max = par.full.pot + nir.full.pot
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the actual total for PAR and NIR, using equations 7 and 8.                   #
   #---------------------------------------------------------------------------------------#
   if (rad.type %in% "par"){
      ratio = ifelse(test = night, yes = 0., no = rad.in / par.full.pot)
   }else if (rad.type %in% "nir"){
      ratio = ifelse(test = night, yes = 0., no = rad.in / nir.full.pot)
   }else{
      ratio = ifelse(test = night, yes = 0., no = rad.in / (par.full.pot + nir.full.pot))
   }#end if
   par.full = ratio * par.full.pot
   nir.full = ratio * nir.full.pot
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the fraction of PAR and NIR that stays as beam, using equations 11 and 12    #
   # of WN85.                                                                              #
   #---------------------------------------------------------------------------------------#
   #----- Make sure that the ratio is bounded. --------------------------------------------#
   aux.par  = pmin(wn85.11[1],pmax(0.,ratio))
   aux.nir  = pmin(wn85.12[1],pmax(0.,ratio))

   fvis.beam.act = ( par.beam.pot 
                   * (1. - ((wn85.11[1] - aux.par)/wn85.11[2]) ^ twothirds)
                   / par.full.pot )
   fvis.beam.act = pmin(1.,pmax(0.,fvis.beam.act))

   fnir.beam.act = ( nir.beam.pot 
                   * (1. - ((wn85.12[1] - aux.nir)/wn85.12[2]) ^ twothirds)
                   / nir.full.pot )
   fnir.beam.act = pmin(1.,pmax(0.,fvis.beam.act))

   fvis.diff.act = 1. - fvis.beam.act
   fnir.diff.act = 1. - fnir.beam.act
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the radiation components.                                                    #
   #---------------------------------------------------------------------------------------#
   par.beam    = fvis.beam.act * par.full
   par.diff    = fvis.diff.act * par.full
   nir.beam    = fnir.beam.act * nir.full
   nir.diff    = fnir.diff.act * nir.full
   rshort.beam = par.beam      + nir.beam
   rshort.diff = par.diff      + nir.diff
   rshort.full = rshort.beam   + rshort.diff
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Last thing to check is whether this is night time.  If the zenith angle is more   #
   # than the typical angle for civil twilight, assume nighttime and set radiation to zero #
   # (true radiation will still be slightly more than zero, but this avoids numerical      #
   # errors.                                                                               #
   #---------------------------------------------------------------------------------------#
   par.beam   [night] = 0.0
   nir.beam   [night] = 0.0
   par.diff   [night] = 0.0
   nir.diff   [night] = 0.0
   par.full   [night] = 0.0
   nir.full   [night] = 0.0
   rshort.beam[night] = 0.0
   rshort.diff[night] = 0.0
   rshort.full[night] = 0.0
   par.max    [night] = 0.0
   nir.max    [night] = 0.0
   rshort.max [night] = 0.0
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Save output data frame.                                                           #
   #---------------------------------------------------------------------------------------#
   rshort.bdown = data.frame( chapman     = chapman
                            , eff.cosz    = eff.cosz
                            , par.beam    = par.beam
                            , par.diff    = par.diff
                            , par.full    = par.full
                            , nir.beam    = nir.beam
                            , nir.diff    = nir.diff
                            , nir.full    = nir.full
                            , rshort.beam = rshort.beam
                            , rshort.diff = rshort.diff
                            , rshort.full = rshort.full
                            , par.max     = par.max
                            , nir.max     = nir.max
                            , rshort.max  = rshort.max
                            )#end list

   return(rshort.bdown)
}#end function rshort.bdown.weissnorman
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#      This subroutine computes the split between direct and diffuse radiation, and        #
# between visible and near-infrared radiation, using the clearness index method (BD08,     #
# TW05, BX10).                                                                             #
#------------------------------------------------------------------------------------------#
rshort.bdown.clearidx <<- function(rad.in,atm.prss,cosz,apply.bx10.corr
                                  ,rad.type=c("rshort","par","nir")){
   #---------------------------------------------------------------------------------------#
   #    Local constants.                                                                   #
   #---------------------------------------------------------------------------------------#
   #----- Extinction coefficient. (equations 1 and 4 of WN85) -----------------------------#
   par.beam.expext  = -0.185
   nir.beam.expext  = -0.060
   #----- Coefficients for various equations in WN85. -------------------------------------#
   wn85.06 = c( -1.1950, 0.4459, -0.0345 )
   #----- This is the typical conversion of diffuse radiation in sunny days. --------------#
   par2diff.sun = 0.400
   nir2diff.sun = 0.600
   #----- Coefficients for Diffuse/total fraction. (equations 24,26, 32 of BD10) ----------#
   bd08.eqn32      = c(-5.0033, 8.6025)
   #----- PAR/SW fraction coefficients (equation 2 of TW05). ------------------------------#
   tw05.eqn02      = c( 0.613, -0.334,  0.121)
   #---------------------------------------------------------------------------------------#



   #------ Initialise the radiation with NAs (except for rshort.full). --------------------#
   rshort.full = NA_real_ * rad.in
   rshort.beam = NA_real_ * rad.in
   rshort.diff = NA_real_ * rad.in
   par.beam    = NA_real_ * rad.in
   nir.beam    = NA_real_ * rad.in
   par.diff    = NA_real_ * rad.in
   nir.diff    = NA_real_ * rad.in
   par.full    = NA_real_ * rad.in
   nir.full    = NA_real_ * rad.in
   par.max     = NA_real_ * rad.in
   nir.max     = NA_real_ * rad.in
   rshort.max  = NA_real_ * rad.in
   #---------------------------------------------------------------------------------------#


   #------ Make day and night flags. ------------------------------------------------------#
   ntimes   = length(cosz)
   day      = cosz %ge% cosz.min
   night    = cosz %lt% cosz.twilight
   twilight = (! day) & (! night)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Save the Chapman function to account for the earth's curvature and find the        #
   # effective cosine of the zenith angle.                                                 #
   #---------------------------------------------------------------------------------------#
   chapman      = pmin(lnexp.max,huestis.fun(cosz=cosz))
   eff.cosz     = ifelse( test = night, yes = 0., no = 1./chapman)
   log10chapman = log10(chapman)
   #---------------------------------------------------------------------------------------#


   #----- Total radiation at the top of the atmosphere [  W/m2], using ED defaults. -------#
   rshort.beam.toa = solar * eff.cosz
   par.beam.toa    =       tw05.eqn02[1]  * rshort.beam.toa
   nir.beam.toa    = (1. - tw05.eqn02[1]) * rshort.beam.toa
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Find the potential PAR components (beam, diffuse, total), using equations 1, 3,    #
   # and 9 of WN85.                                                                        #
   #---------------------------------------------------------------------------------------#
   par.beam.pot = ( par.beam.toa
                  * exp ( par.beam.expext * (atm.prss / prefsea) * chapman) )
   par.diff.pot = par2diff.sun * (par.beam.toa - par.beam.pot)
   par.full.pot = par.beam.pot + par.diff.pot
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the NIR absorption of 10 mm of precipitable water, using WN85 equation 6.    #
   #---------------------------------------------------------------------------------------#
   w10 = ( rshort.beam.toa
         * 10 ** ((wn85.06[1]) + log10chapman * (wn85.06[2] + wn85.06[3] * log10chapman))
         )#end w10
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the potential direct and diffuse near-infrared radiation, using equations    #
   # 4, 5, and 10 of WN85.                                                                 #
   #---------------------------------------------------------------------------------------#
   nir.beam.pot = ( ( nir.beam.toa
                    * exp ( nir.beam.expext * (atm.prss / prefsea) * chapman) - w10 ) )
   nir.diff.pot = nir2diff.sun * ( nir.beam.toa - nir.beam.pot - w10 )
   nir.full.pot = nir.beam.pot + nir.diff.pot
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    In case of twilight, set beam radiation to zero and diffuse radiation to           #
   # potential.                                                                            #
   #---------------------------------------------------------------------------------------#
   par.beam.pot[twilight] = 0.0
   par.diff.pot[twilight] = par.full.pot[twilight]
   nir.beam.pot[twilight] = 0.0
   nir.diff.pot[twilight] = nir.full.pot[twilight]
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Find the clearness index based on total shortwave radiation, then find the         #
   # fraction of PAR radiation as a function of total radiation and clearness index,       #
   # following TW05, eqn. 2.                                                               #
   #---------------------------------------------------------------------------------------#
   if (rad.type %in% "par"){
      par.full    = rad.in
      fkt         = ifelse( test = night
                          , yes  = 0.
                          , no   = pmax(0.,pmin(1., par.full / par.beam.toa))
                          )#end ifelse
      fpar        = tw05.eqn02[1] + fkt * ( tw05.eqn02[2] + fkt * tw05.eqn02[3] )
      rshort.full = par.full / fpar
      nir.full    = rshort.full - par.full
   }else if (rad.type %in% "nir"){
      nir.full    = rad.in
      fkt         = ifelse( test = night
                          , yes  = 0.
                          , no   = pmax(0.,pmin(1., nir.full / nir.beam.toa))
                          )#end ifelse
      fpar        = tw05.eqn02[1] + fkt * ( tw05.eqn02[2] + fkt * tw05.eqn02[3] )
      rshort.full = nir.full / (1.-fpar)
      par.full    = rshort.full - nir.full

   }else{
      rshort.full = rad.in
      fkt         = ifelse( test = night
                          , yes  = 0.
                          , no   = pmax(0.,pmin(1., rshort.full / rshort.beam.toa))
                          )#end ifelse
      fpar        = tw05.eqn02[1] + fkt * ( tw05.eqn02[2] + fkt * tw05.eqn02[3] )
      par.full    = fpar * rshort.full
      nir.full    = rshort.full - par.full
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the uncorrected diffuse fraction based on BD08.  We use this equation        #
   # instead of the method proposed by BX10 because it is a continuous function and        #
   # results are very similar to BX10 curve except at high fkt, where BD08 predicts higher #
   # fraction of direct -- a reasonable assumption because total radiation should be       #
   # entirely direct if the it is the same as the TOA radiation.                           #
   #---------------------------------------------------------------------------------------#
   fdiff.aux = pmax(lnexp.min, pmin(lnexp.max, bd08.eqn32[1] + bd08.eqn32[2] * fkt ) )
   fdiff.1st = 1.0 / (1.0 + exp(fdiff.aux) )
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     One option is to apply the correction term by BX10.  I'm still not sure if this   #
   # term is really necessary, because the higher proportion of diffuse radiation at low   #
   # sun angles should be due to thicker optical depth, and a lower ratio between incident #
   # radiation at the ground and extraterrestrial radiation should be naturally lower.  If #
   # direct radiation is too high, then we may include this term.                          #
   #---------------------------------------------------------------------------------------#
   if (apply.bx10.corr){
      fdiff = fdiff.1st / ( (1.0-fdiff.1st) * cosz + fdiff.1st)
   }else{
      fdiff = fdiff.1st
   }#end if (apply_bx10_corr)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the radiation components.                                                    #
   #---------------------------------------------------------------------------------------#
   par.diff    = fdiff    * par.full
   par.beam    = par.full - par.diff
   nir.diff    = fdiff    * nir.full
   nir.beam    = nir.full - nir.diff
   rshort.diff = par.diff + nir.diff
   rshort.beam = par.beam + nir.beam
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Total maximum radiation.                                                          #
   #---------------------------------------------------------------------------------------#
   par.max    = par.full.pot
   nir.max    = nir.full.pot
   rshort.max = par.full.pot + nir.full.pot
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Last thing to check is whether this is night time.  If the zenith angle is more   #
   # than the typical angle for civil twilight, assume nighttime and set radiation to zero #
   # (true radiation will still be slightly more than zero, but this avoids numerical      #
   # errors.                                                                               #
   #---------------------------------------------------------------------------------------#
   par.beam   [night] = 0.0
   nir.beam   [night] = 0.0
   par.diff   [night] = 0.0
   nir.diff   [night] = 0.0
   par.full   [night] = 0.0
   nir.full   [night] = 0.0
   rshort.beam[night] = 0.0
   rshort.diff[night] = 0.0
   rshort.full[night] = 0.0
   par.max    [night] = 0.0
   nir.max    [night] = 0.0
   rshort.max [night] = 0.0
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Total maximum radiation.                                                          #
   #---------------------------------------------------------------------------------------#
   rshort.bdown = data.frame( chapman     = chapman
                            , eff.cosz    = eff.cosz
                            , par.beam    = par.beam
                            , par.diff    = par.diff
                            , par.full    = par.full
                            , nir.beam    = nir.beam
                            , nir.diff    = nir.diff
                            , nir.full    = nir.full
                            , rshort.beam = rshort.beam
                            , rshort.diff = rshort.diff
                            , rshort.full = rshort.full
                            , par.max     = par.max
                            , nir.max     = nir.max
                            , rshort.max  = rshort.max
                            )#end list

   return(rshort.bdown)
   #---------------------------------------------------------------------------------------#
}#end function rshort.bdown.clearidx
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#      This subroutine computes the split between direct and diffuse radiation, and        #
# between visible and near-infrared radiation, using the SiB method.                       #
#------------------------------------------------------------------------------------------#
rshort.bdown.sib <<- function(rad.in,atm.prss,cosz,apply.bx10.corr){
   #---------------------------------------------------------------------------------------#
   #    Local constants.                                                                   #
   #---------------------------------------------------------------------------------------#
   #----- Extinction coefficient. (equations 1 and 4 of WN85) -----------------------------#
   par.beam.expext  = -0.185
   nir.beam.expext  = -0.060
   #----- Coefficients for various equations in WN85. -------------------------------------#
   wn85.06 = c( -1.1950, 0.4459, -0.0345 )
   #----- This is the typical conversion of diffuse radiation in sunny days. --------------#
   par2diff.sun = 0.400
   nir2diff.sun = 0.600
   #----- Local constants. ----------------------------------------------------------------#
   csib         = c(580.,464.,499.,963.,1160.)
   #---------------------------------------------------------------------------------------#


   #------ Initialise the radiation with NAs (except for rshort.full). --------------------#
   rshort.full = rad.in
   rshort.beam = NA_real_ * rad.in
   rshort.diff = NA_real_ * rad.in
   par.beam    = NA_real_ * rad.in
   nir.beam    = NA_real_ * rad.in
   par.diff    = NA_real_ * rad.in
   nir.diff    = NA_real_ * rad.in
   par.full    = NA_real_ * rad.in
   nir.full    = NA_real_ * rad.in
   par.max     = NA_real_ * rad.in
   nir.max     = NA_real_ * rad.in
   rshort.max  = NA_real_ * rad.in
   #---------------------------------------------------------------------------------------#


   #------ Make day and night flags. ------------------------------------------------------#
   ntimes   = length(cosz)
   day      = cosz %gt% cosz.min
   night    = cosz %lt% cosz.twilight
   twilight = (! day) & (! night)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Save the Chapman function to account for the earth's curvature and find the        #
   # effective cosine of the zenith angle.                                                 #
   #---------------------------------------------------------------------------------------#
   chapman      = pmin(lnexp.max,huestis.fun(cosz=cosz))
   eff.cosz     = ifelse( test = night, yes = 0., no = 1./chapman)
   log10chapman = log10(chapman)
   #---------------------------------------------------------------------------------------#



   #----- Total radiation at the top of the atmosphere [  W/m2], using ED defaults. -------#
   rshort.beam.toa = solar * eff.cosz
   par.beam.toa    = fvis.beam.def * rshort.beam.toa
   nir.beam.toa    = fnir.beam.def * rshort.beam.toa
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Find the potential PAR components (beam, diffuse, total), using equations 1, 3,    #
   # and 9 of WN85.                                                                        #
   #---------------------------------------------------------------------------------------#
   par.beam.pot = ( par.beam.toa
                  * exp ( par.beam.expext * (atm.prss / prefsea) * chapman) )
   par.diff.pot = par2diff.sun * (par.beam.toa - par.beam.pot) 
   par.full.pot = par.beam.pot + par.diff.pot
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the NIR absorption of 10 mm of precipitable water, using WN85 equation 6.    #
   #---------------------------------------------------------------------------------------#
   w10 = ( rshort.beam.toa
         * 10 ** ((wn85.06[1]) + log10chapman * (wn85.06[2] + wn85.06[3] * log10chapman))
         )#end w10
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the potential direct and diffuse near-infrared radiation, using equations    #
   # 4, 5, and 10 of WN85.                                                                 #
   #---------------------------------------------------------------------------------------#
   nir.beam.pot = ( ( nir.beam.toa
                    * exp ( nir.beam.expext * (atm.prss / prefsea) * chapman) - w10 )  )
   nir.diff.pot = nir2diff.sun * ( nir.beam.toa - nir.beam.pot - w10 )
   nir.full.pot = nir.beam.pot + nir.diff.pot
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    In case of twilight, set beam radiation to zero and diffuse radiation to           #
   # potential.                                                                            #
   #---------------------------------------------------------------------------------------#
   par.beam.pot[twilight] = 0.0
   par.diff.pot[twilight] = par.full.pot[twilight]
   nir.beam.pot[twilight] = 0.0
   nir.diff.pot[twilight] = nir.full.pot[twilight]
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Find the cloud cover estimate and the fraction of diffuse radiation.               #
   #---------------------------------------------------------------------------------------#
   cloud  = pmin(1.,pmax(0.,(csib[5] * eff.cosz - rshort.full) / (csib[4] * eff.cosz)))
   difrat = pmin(1.,pmax(0.,0.0604 / ( eff.cosz -0.0223 ) + 0.0683))
   difrat = difrat + ( 1. - difrat ) * cloud
   vnrat  = ( ( csib[1] - cloud*csib[2] ) 
            / ( ( csib[1] - cloud*csib[3] ) + ( csib[1] - cloud*csib[2] ))
            )#end vnrat

   rshort.diff = ifelse(test=day,yes=difrat * vnrat * rshort.full,no=rshort.full)
   rshort.beam = rshort.full - rshort.diff
   par.diff    = fvis.diff.def * rshort.diff
   nir.diff    = fnir.diff.def * rshort.diff
   par.beam    = fvis.beam.def * rshort.beam
   nir.beam    = fnir.beam.def * rshort.beam
   par.full    = par.diff + par.beam
   nir.full    = nir.diff + nir.beam
   #---------------------------------------------------------------------------------------#
 


   #---------------------------------------------------------------------------------------#
   #     Total maximum radiation.                                                          #
   #---------------------------------------------------------------------------------------#
   par.max    = par.full.pot
   nir.max    = nir.full.pot
   rshort.max = par.full.pot + nir.full.pot
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Last thing to check is whether this is night time.  If the zenith angle is more   #
   # than the typical angle for civil twilight, assume nighttime and set radiation to zero #
   # (true radiation will still be slightly more than zero, but this avoids numerical      #
   # errors.                                                                               #
   #---------------------------------------------------------------------------------------#
   par.beam   [night] = 0.0
   nir.beam   [night] = 0.0
   par.diff   [night] = 0.0
   nir.diff   [night] = 0.0
   par.full   [night] = 0.0
   nir.full   [night] = 0.0
   rshort.beam[night] = 0.0
   rshort.diff[night] = 0.0
   rshort.full[night] = 0.0
   par.max    [night] = 0.0
   nir.max    [night] = 0.0
   rshort.max [night] = 0.0
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Total maximum radiation.                                                          #
   #---------------------------------------------------------------------------------------#
   rshort.bdown = data.frame( par.beam    = par.beam
                            , par.diff    = par.diff
                            , par.full    = par.full
                            , nir.beam    = nir.beam
                            , nir.diff    = nir.diff
                            , nir.full    = nir.full
                            , rshort.beam = rshort.beam
                            , rshort.diff = rshort.diff
                            , rshort.full = rshort.full
                            , par.max     = par.max
                            , nir.max     = nir.max
                            , rshort.max  = rshort.max
                            )#end list

   return(rshort.bdown)
   #---------------------------------------------------------------------------------------#
}#end function rshort.bdown.sib
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#    This function calculates the effect of sun angle increasing the optical depth in a    #
# sphere.  This allows accounting for the diffuse radiation at twilight, using the         #
# modified Chapman function following:                                                     #
#                                                                                          #
#    Reference:                                                                            #
#                                                                                          #
#    Huestis DL. 2001. Accurate evaluation of the chapman function for atmospheric         #
#       attenuation. J. Quant. Spectrosc. Radiat. Transf., 69: 709-721.                    #
#       doi:10.1016/S0022-4073(00)00107-2                                                  #
#                                                                                          #
#                                                                                          #
#    Input variables:                                                                      #
#    - tvir. Virtual temperature, in Kelvin.                                               #
#    - cosz. Cosine of the sun's zenith angle.                                             #
#    - dzen. Bin width of the zenith angle for the numerical integration                   #
#------------------------------------------------------------------------------------------#
huestis.fun <<- function(cosz,alt=0.,dzen=0.05){
   #----- Dimensionless curvature ratio. --------------------------------------------------#
   xx  = ( erad + alt ) / ehgt
   #---------------------------------------------------------------------------------------#


   #---- Create the look-up table. --------------------------------------------------------#
   zend.ref    = seq(from=0,to=180-dzen,by=dzen)
   lambda.ref  = mid.points(zend.ref)
   nzen        = length(zend.ref)
   ZEND.REF    = matrix(data=zend.ref,nrow=nzen,ncol=nzen,byrow=FALSE)
   LAMBDA.REF  = matrix(data=zend.ref,nrow=nzen,ncol=nzen,byrow=TRUE )
   LAMBDA.REF  = pmin(ZEND.REF,LAMBDA.REF)
   DLAMBDA.REF = t(apply(X=LAMBDA.REF,MARGIN=1,FUN=diff))
   LAMBDA.REF  = t(apply(X=LAMBDA.REF,MARGIN=1,FUN=mid.points))
   ZEND.REF    = ZEND.REF[,-nzen]
   SINZ.REF    = sin(ZEND.REF*pio180)
   SINL.REF    = sin(LAMBDA.REF*pio180)
   COSL.REF    = cos(LAMBDA.REF*pio180)
   KERNEL      = ifelse( test = SINL.REF %eq% 0.
                       , yes  = 0.
                       , no   = exp(xx*(1.-SINZ.REF/SINL.REF))/(1.+COSL.REF)
                       )#end ifelse
   huestis.ref = 1+xx*sin(zend.ref*pio180) * rowSums(KERNEL*DLAMBDA.REF*pio180)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Create matrices with the reference.                                               #
   #---------------------------------------------------------------------------------------#
   zend = dzen * round(acos(cosz)/pio180/dzen)
   idx  = match(zend,zend.ref)
   ans  = 0.*cosz + c(huestis.ref[idx])
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function huestis.fun
#==========================================================================================#
#==========================================================================================#
