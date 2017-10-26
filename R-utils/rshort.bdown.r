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
   if ( (! rad.method %in% "wn85") && (! (rad.type %in% "rshort")) ){
      cat0(" - Radiation method: ",rad.method,".")
      cat0(" - Input radiation type: ",rad.type,".")
      stop(" Radiation input type must be \"rshort\" unless using method \"wn85\".)")
   }#end if
   #---------------------------------------------------------------------------------------#


   ans = switch( EXPR     = rad.method
               , wn85     = rshort.bdown.weissnorman(rad.in,atm.prss,cosz,rad.type)
               , clearidx = rshort.bdown.clearidx   (rad.in,atm.prss,cosz,apply.bx10.corr)
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
   par.beam    = NA * rad.in
   nir.beam    = NA * rad.in
   par.diff    = NA * rad.in
   nir.diff    = NA * rad.in
   par.full    = NA * rad.in
   nir.full    = NA * rad.in
   rshort.beam = NA * rad.in
   rshort.diff = NA * rad.in
   rshort.full = NA * rad.in
   par.max     = NA * rad.in
   nir.max     = NA * rad.in
   rshort.max  = NA * rad.in


   #------ Make day and night flags. ------------------------------------------------------#
   ntimes = length(cosz)
   day    = cosz %>% cosz.min
   night  = ! day
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     First thing to check is whether this is daytime or "night-time".  If the zenith   #
   # angle is too close to horizon, we assume it's dawn/dusk and all radiation goes to     #
   # diffuse.                                                                              #
   #---------------------------------------------------------------------------------------#
   par.beam    [night] = 0.0
   nir.beam    [night] = 0.0
   if (rad.type %in% "par"){
      par.diff    [night] = rad.in[night]
      nir.diff    [night] = fnir.diff.def * rad.in[night] / fvis.diff.def
   }else if(rad.type %in% "nir"){
      par.diff    [night] = fvis.diff.def * rad.in[night] / fnir.diff.def
      nir.diff    [night] = rad.in[night]
   }else{
      par.diff    [night] = fvis.diff.def * rad.in[night]
      nir.diff    [night] = fnir.diff.def * rad.in[night]
   }#end if
   par.full    [night] = par.beam   [night] + par.diff   [night]
   nir.full    [night] = nir.beam   [night] + nir.diff   [night]
   rshort.beam [night] = par.beam   [night] + nir.beam   [night]
   rshort.diff [night] = par.diff   [night] + nir.diff   [night]
   rshort.full [night] = rshort.beam[night] + rshort.diff[night]
   par.max     [night] = 0.0
   nir.max     [night] = 0.0
   rshort.max  [night] = 0.0
   #---------------------------------------------------------------------------------------#



   #----- Save 1/cos(zen), which is the secant.  We will use this several times. ----------#
   secz      = ifelse(day, 1. / cosz, 0)
   log10secz = log10(secz)
   #---------------------------------------------------------------------------------------#


   #----- Total radiation at the top of the atmosphere [  W/m2], using ED defaults. -------#
   rshort.beam.toa = solar * ifelse(day,cosz,0.)
   par.beam.toa    = fvis.beam.def * rshort.beam.toa
   nir.beam.toa    = fnir.beam.def * rshort.beam.toa
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Find the potential PAR components (beam, diffuse, total), using equations 1, 3,    #
   # and 9 of WN85.                                                                        #
   #---------------------------------------------------------------------------------------#
   par.beam.pot = ( par.beam.toa
                  * exp ( par.beam.expext * (atm.prss / prefsea) * secz) )
   par.diff.pot = par2diff.sun * (par.beam.toa - par.beam.pot)
   par.full.pot = par.beam.pot + par.diff.pot
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the NIR absorption of 10 mm of precipitable water, using WN85 equation 6.    #
   #---------------------------------------------------------------------------------------#
   w10 = ( rshort.beam.toa
         * 10 ** ((wn85.06[1]) + log10secz * (wn85.06[2] + wn85.06[3] * log10secz))
         )#end w10
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the potential direct and diffuse near-infrared radiation, using equations    #
   # 4, 5, and 10 of WN85.                                                                 #
   #---------------------------------------------------------------------------------------#
   nir.beam.pot = ( ( nir.beam.toa
                    * exp ( nir.beam.expext * (atm.prss / prefsea) * secz) - w10 ) )
   nir.diff.pot = nir2diff.sun * ( nir.beam.toa - nir.beam.pot - w10 )
   nir.full.pot = nir.beam.pot + nir.diff.pot
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Total maximum radiation.                                                          #
   #---------------------------------------------------------------------------------------#
   par.max   [day] = par.full.pot[day]
   nir.max   [day] = nir.full.pot[day]
   rshort.max[day] = par.full.pot[day] + nir.full.pot[day]
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the actual total for PAR and NIR, using equations 7 and 8.                   #
   #---------------------------------------------------------------------------------------#
   if (rad.type %in% "par"){
      ratio = ifelse(day, rad.in / par.full.pot, 0.)
   }else if (rad.type %in% "nir"){
      ratio = ifelse(day, rad.in / nir.full.pot, 0.)
   }else{
      ratio = ifelse(day, rad.in / (par.full.pot + nir.full.pot), 0.)
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
   par.beam    [day] = fvis.beam.act[day] * par.full[day]
   par.diff    [day] = fvis.diff.act[day] * par.full[day]
   nir.beam    [day] = fnir.beam.act[day] * nir.full[day]
   nir.diff    [day] = fnir.diff.act[day] * nir.full[day]
   rshort.beam [day] = par.beam   [day] + nir.beam   [day]
   rshort.diff [day] = par.diff   [day] + nir.diff   [day]
   rshort.full [day] = rshort.beam[day] + rshort.diff[day]
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
}#end function rshort.bdown.weissnorman
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#      This subroutine computes the split between direct and diffuse radiation, and        #
# between visible and near-infrared radiation, using the clearness index method (BD08,     #
# TW05, BX10).                                                                             #
#------------------------------------------------------------------------------------------#
rshort.bdown.clearidx <<- function(rad.in,atm.prss,cosz,apply.bx10.corr){
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
   rshort.full = rad.in
   rshort.beam = NA * rad.in
   rshort.diff = NA * rad.in
   par.beam    = NA * rad.in
   nir.beam    = NA * rad.in
   par.diff    = NA * rad.in
   nir.diff    = NA * rad.in
   par.full    = NA * rad.in
   nir.full    = NA * rad.in
   par.max     = NA * rad.in
   nir.max     = NA * rad.in
   rshort.max  = NA * rad.in
   #---------------------------------------------------------------------------------------#


   #------ Make day and night flags. ------------------------------------------------------#
   ntimes = length(cosz)
   day    = cosz %>% cosz.min
   night  = ! day
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     First thing to check is whether this is daytime or "night-time".  If the zenith   #
   # angle is too close to horizon, we assume it's dawn/dusk and all radiation goes to     #
   # diffuse.                                                                              #
   #---------------------------------------------------------------------------------------#
   par.beam    [night] = 0.0
   nir.beam    [night] = 0.0
   par.diff    [night] =       tw05.eqn02[1]  * rshort.full[night]
   nir.diff    [night] = (1. - tw05.eqn02[1]) * rshort.full[night]
   par.full    [night] = par.beam   [night] + par.diff   [night]
   nir.full    [night] = nir.beam   [night] + nir.diff   [night]
   rshort.beam [night] = par.beam   [night] + nir.beam   [night]
   rshort.diff [night] = par.diff   [night] + nir.diff   [night]
   par.max     [night] = 0.0
   nir.max     [night] = 0.0
   rshort.max  [night] = 0.0
   #---------------------------------------------------------------------------------------#



   #----- Save 1/cos(zen), which is the secant.  We will use this several times. ----------#
   secz      = ifelse(day, 1. / cosz, 0.)
   log10secz = log10(secz)
   #---------------------------------------------------------------------------------------#


   #----- Total radiation at the top of the atmosphere [  W/m2], using ED defaults. -------#
   rshort.beam.toa = solar * ifelse(day,cosz,0.)
   par.beam.toa    =       tw05.eqn02[1]  * rshort.beam.toa
   nir.beam.toa    = (1. - tw05.eqn02[1]) * rshort.beam.toa
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Find the potential PAR components (beam, diffuse, total), using equations 1, 3,    #
   # and 9 of WN85.                                                                        #
   #---------------------------------------------------------------------------------------#
   par.beam.pot = ( par.beam.toa
                  * exp ( par.beam.expext * (atm.prss / prefsea) * secz) )
   par.diff.pot = par2diff.sun * (par.beam.toa - par.beam.pot)
   par.full.pot = par.beam.pot + par.diff.pot
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the NIR absorption of 10 mm of precipitable water, using WN85 equation 6.    #
   #---------------------------------------------------------------------------------------#
   w10 = ( rshort.beam.toa
         * 10 ** ((wn85.06[1]) + log10secz * (wn85.06[2] + wn85.06[3] * log10secz))
         )#end w10
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the potential direct and diffuse near-infrared radiation, using equations    #
   # 4, 5, and 10 of WN85.                                                                 #
   #---------------------------------------------------------------------------------------#
   nir.beam.pot = ( ( nir.beam.toa
                    * exp ( nir.beam.expext * (atm.prss / prefsea) * secz) - w10 ) )
   nir.diff.pot = nir2diff.sun * ( nir.beam.toa - nir.beam.pot - w10 )
   nir.full.pot = nir.beam.pot + nir.diff.pot
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Find the clearness index based on total shortwave radiation, then find the         #
   # fraction of PAR radiation as a function of total radiation and clearness index,       #
   # following TW05, eqn. 2.                                                               #
   #---------------------------------------------------------------------------------------#
   fkt           = ifelse(day,pmax(0.,pmin(1.,rshort.full / rshort.beam.toa)),0.)
   fpar          = tw05.eqn02[1] + fkt * ( tw05.eqn02[2] + fkt * tw05.eqn02[3] )
   par.full[day] = fpar[day] * rshort.full[day]
   nir.full[day] = rshort.full[day] - par.full[day]
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
   par.diff   [day] = fdiff   [day] * par.full[day]
   par.beam   [day] = par.full[day] - par.diff[day]
   nir.diff   [day] = fdiff   [day] * nir.full[day]
   nir.beam   [day] = nir.full[day] - nir.diff[day]
   rshort.diff[day] = par.diff[day] + nir.diff[day]
   rshort.beam[day] = par.beam[day] + nir.beam[day]
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Total maximum radiation.                                                          #
   #---------------------------------------------------------------------------------------#
   par.max   [day] = par.full.pot[day]
   nir.max   [day] = nir.full.pot[day]
   rshort.max[day] = par.full.pot[day] + nir.full.pot[day]
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
   rshort.beam = NA * rad.in
   rshort.diff = NA * rad.in
   par.beam    = NA * rad.in
   nir.beam    = NA * rad.in
   par.diff    = NA * rad.in
   nir.diff    = NA * rad.in
   par.full    = NA * rad.in
   nir.full    = NA * rad.in
   par.max     = NA * rad.in
   nir.max     = NA * rad.in
   rshort.max  = NA * rad.in
   #---------------------------------------------------------------------------------------#


   #------ Make day and night flags. ------------------------------------------------------#
   ntimes = length(cosz)
   day    = cosz %>% cosz.min
   night  = ! day
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     First thing to check is whether this is daytime or "night-time".  If the zenith   #
   # angle is too close to horizon, we assume it's dawn/dusk and all radiation goes to     #
   # diffuse.                                                                              #
   #---------------------------------------------------------------------------------------#
   par.beam    [night] = 0.0
   nir.beam    [night] = 0.0
   par.diff    [night] = fvis.diff.def * rad.in[night]
   nir.diff    [night] = fnir.diff.def * rad.in[night]
   par.full    [night] = par.beam   [night] + par.diff   [night]
   nir.full    [night] = nir.beam   [night] + nir.diff   [night]
   rshort.beam [night] = par.beam   [night] + nir.beam   [night]
   rshort.diff [night] = par.diff   [night] + nir.diff   [night]
   par.max     [night] = 0.0
   nir.max     [night] = 0.0
   rshort.max  [night] = 0.0
   #---------------------------------------------------------------------------------------#



   #----- Save 1/cos(zen), which is the secant.  We will use this several times. ----------#
   secz      = ifelse(day, 1. / cosz, 0.)
   log10secz = log10(secz)
   #---------------------------------------------------------------------------------------#


   #----- Total radiation at the top of the atmosphere [  W/m2], using ED defaults. -------#
   rshort.beam.toa = solar * ifelse(day,cosz,0.)
   par.beam.toa    = fvis.beam.def * rshort.beam.toa
   nir.beam.toa    = fnir.beam.def * rshort.beam.toa
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Find the potential PAR components (beam, diffuse, total), using equations 1, 3,    #
   # and 9 of WN85.                                                                        #
   #---------------------------------------------------------------------------------------#
   par.beam.pot = ( par.beam.toa
                  * exp ( par.beam.expext * (atm.prss / prefsea) * secz) )
   par.diff.pot = par2diff.sun * (par.beam.toa - par.beam.pot) 
   par.full.pot = par.beam.pot + par.diff.pot
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the NIR absorption of 10 mm of precipitable water, using WN85 equation 6.    #
   #---------------------------------------------------------------------------------------#
   w10 = ( rshort.beam.toa
         * 10 ** ((wn85.06[1]) + log10secz * (wn85.06[2] + wn85.06[3] * log10secz))
         )#end w10
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the potential direct and diffuse near-infrared radiation, using equations    #
   # 4, 5, and 10 of WN85.                                                                 #
   #---------------------------------------------------------------------------------------#
   nir.beam.pot = ( ( nir.beam.toa
                    * exp ( nir.beam.expext * (atm.prss / prefsea) * secz) - w10 )  )
   nir.diff.pot = nir2diff.sun * ( nir.beam.toa - nir.beam.pot - w10 )
   nir.full.pot = nir.beam.pot + nir.diff.pot
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Find the cloud cover estimate and the fraction of diffuse radiation.               #
   #---------------------------------------------------------------------------------------#
   cloud  = ifelse( day
                  , pmin(1.,pmax(0.,(csib[5] * cosz - rshort.full) / (csib[4] * cosz)))
                  , 0.5
                  )#end ifelse
   difrat = ifelse(day,pmin(1.,pmax(0.,0.0604 / ( cosz -0.0223 ) + 0.0683)),1.0)
   difrat = difrat + ( 1. - difrat ) * cloud
   vnrat  = ( ( csib[1] - cloud*csib[2] ) 
            / ( ( csib[1] - cloud*csib[3] ) + ( csib[1] - cloud*csib[2] ))
            )#end vnrat

   rshort.diff[day] = difrat[day] * vnrat[day] * rshort.full[day]
   rshort.beam[day] = rshort.full[day] - rshort.diff[day]
   par.diff   [day] = fvis.diff.def * rshort.diff[day]
   nir.diff   [day] = fnir.diff.def * rshort.diff[day]
   par.beam   [day] = fvis.beam.def * rshort.beam[day]
   nir.beam   [day] = fnir.beam.def * rshort.beam[day]
   par.full   [day] = par.diff[day] + par.beam[day]
   nir.full   [day] = nir.diff[day] + nir.beam[day]
   #---------------------------------------------------------------------------------------#
 


   #---------------------------------------------------------------------------------------#
   #     Total maximum radiation.                                                          #
   #---------------------------------------------------------------------------------------#
   par.max   [day] = par.full.pot[day]
   nir.max   [day] = nir.full.pot[day]
   rshort.max[day] = par.full.pot[day] + nir.full.pot[day]
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
