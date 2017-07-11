#==========================================================================================#
#==========================================================================================#
#     This function splits the photosynthetically active radiation into direct and         #
# diffuse.                                                                                 #
#                                                                                          #
# Weiss, A., J. M. Norman, 1985: Partitioning solar radiation into direct and diffuse,     #
#     visible and near-infrared components.  Agric. For. Meteorol., 34, 205-213. (WN85)    #
#------------------------------------------------------------------------------------------#
par.split = function(cosz,partop,atm.prss){
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
   wn85.11 = c(  0.90, 0.70 )
   wn85.12 = c(  0.88, 0.68 )

   #----- Initialise terms. ---------------------------------------------------------------#
   par.diff = cosz * NA
   par.beam = cosz * NA

   if (cosz < cosz.min){
      par.diff = partop
      par.beam = 0.
      par.out = list(beam=par.beam,diff=par.diff)
      return(par.out)
   }#end if

   #----- Save 1/cos(zen), which is the secant.  We will use this several times. ----------#
   secz      = 1. / cosz
   log10secz = log10(secz)
   #---------------------------------------------------------------------------------------#

   #----- Total radiation at the top [  W/m2], using ED defaults. -------------------------#
   par.beam.top = fvis.beam.def * solar

   #---------------------------------------------------------------------------------------#
   #    Find the potential PAR components (beam, diffuse, total), using equations 1, 3,    #
   # and 9 of WN85.                                                                        #
   #---------------------------------------------------------------------------------------#
   par.beam.pot = ( par.beam.top * exp ( par.beam.expext * (atm.prss / prefsea) * secz)
                  * cosz )
   par.diff.pot = par2diff.sun * (par.beam.top - par.beam.pot) * cosz
   par.full.pot = par.beam.pot + par.diff.pot
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the actual total for PAR and NIR, using equations 7 and 8.                   #
   #---------------------------------------------------------------------------------------#
   ratio    = partop / par.full.pot
   par.full = ratio * par.full.pot
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the fraction of PAR and NIR that stays as beam, using equations 11 and 12    #
   # of WN85.                                                                              #
   #---------------------------------------------------------------------------------------#
   #----- Make sure that the ratio is bounded. --------------------------------------------#
   aux.par       = min(wn85.11[1], max(0., ratio))
   fvis.beam.act = min(1., max(0., par.beam.pot
                              * (1. - ((wn85.11[1] - aux.par)/wn85.11[2]) ^ twothirds)
                              / par.full.pot ) )
   fvis.diff.act = 1. - fvis.beam.act
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the radiation components.                                                    #
   #---------------------------------------------------------------------------------------#
   par.beam    = fvis.beam.act * par.full
   par.diff    = fvis.diff.act * par.full
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   par.out = list(beam=par.beam,diff=par.diff)
   return(par.out)
}#end function

