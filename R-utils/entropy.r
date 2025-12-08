#==========================================================================================#
#==========================================================================================#
#      Functions to help computing entropy fluxes.  Most of the theory can be found in the #
# following papers.                                                                        #
#                                                                                          #
# Quijano JC , Lin H. 2015. Is spatially integrated entropy production useful to predict   #
#    the dynamics of ecosystems? Ecol. Model., 313: 341-354.                               #
#    doi:10.1016/j.ecolmodel.2015.06.012 (QJ15).                                           #
#                                                                                          #
# Wright SE, Scott DS, Haddow JB , Rosen MA. 2001. On the entropy of radiative heat        #
#    transfer in engineering thermodynamics. Int. J. Eng. Sci., 39: 1691-1706.             #
#    doi:10.1016/S0020-7225(01)00024-6 (W01).                                              #
#                                                                                          #
# Wu W , Liu Y. 2010. Radiation entropy flux and entropy production of the Earth system.   #
#    Rev. Geophys., 48: RG2003. doi:10.1029/2008RG000275 (WL10).                           #
#                                                                                          #
#------------------------------------------------------------------------------------------#


#------ Handy constants. ------------------------------------------------------------------#
w01.rf.sbeam <<- 2.31e-4 * tsun
w01.c0       <<- -45 / (4. * pi^4)
w01.c1       <<- +2.336
w01.c2       <<- -0.260
w01.rpfac    <<- 4./3.
emis.atm     <<- 0.85             # From Brunsell et al. (2011)
emis.sfc     <<- 0.98             # From ED2
#------------------------------------------------------------------------------------------#



#==========================================================================================#
#==========================================================================================#
#     This function finds the entropy flux due to shortwave radiation, following Q15,      #
# who base their formulation on WL10 (which is done for the top of the atmosphere).  This  #
# assumes that albedo is the same for direct and diffuse radiation, and that all reflected #
# radiation is diffuse.                                                                    #
#------------------------------------------------------------------------------------------#
fxentropy.rshort <<- function(rsdnwd,rswnet,rsbeam){

   #----- Find derived radiation properties. ----------------------------------------------#
   rsdiff = 0. * rsdnwd + pmax(0.,rsdnwd - rsbeam)        # Diffuse radiation as residual
   rsupwd = pmax(0.,rsdnwd - rswnet)                      # Upward radiation
   albedo = 0. * rsdnwd + ifelse( test = rsdnwd %gt% 0.   # Albedo
                                , yes  = rsupwd / rsdnwd  #
                                , no   = 0.1              #
                                )#end ifelse              #
   #---------------------------------------------------------------------------------------#



   #----- Entropy due to direct radiation. ------------------------------------------------#
   fs.rsbeam = w01.rf.sbeam * (1. - albedo) * rsbeam / tsun
   #---------------------------------------------------------------------------------------#


   #----- Entropy due to diffuse radiation. -----------------------------------------------#
   delta      = rsdiff / pi / radsol
   rf.rsdiff  = ifelse( test = delta %gt% 0.
                      , yes  = ( w01.c0 * ( w01.c1 + w01.c2 * delta) * log(delta) + 1. )
                             * w01.rpfac
                      , no   = 0.
                      )#end ifelse
   fs.rsdiff  = rf.rsdiff * (1. - albedo) * rsdiff / tsun
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Total entropy flux due to solar radiation.                                       #
   #---------------------------------------------------------------------------------------#
   ans = fs.rsbeam + fs.rsdiff
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end sflux.rshort
#==========================================================================================#
#==========================================================================================#




#==========================================================================================#
#==========================================================================================#
#     This function finds the entropy flux due to downwelling longwave radiation, follow-  #
# ing Q15, who base their formulation on WL10 and W01.  Two differences from their         #
# approach:                                                                                #
# 1. We don't use the air temperature, but the equivalent atmospheric temperature, similar #
#    to H10 (unless the user provides atmospheric temperature, in which case their values  #
#    take priority).                                                                       #
# 2. We also account for the small amount of downwelling LW irradiance that is not         #
#    absorbed, using Kirchhoff's law.                                                      #
#------------------------------------------------------------------------------------------#
fxentropy.rldnwd <<- function(rldnwd,tatm.k){

   #----- Find derived radiation properties. ----------------------------------------------#
   if (missing(tatm.k)){
      tatm.k = sqrt(sqrt(rldnwd / (emis.atm * stefan)))
   }#end if (missing(tsfc.k))
   #---------------------------------------------------------------------------------------#


   #----- Entropy due to diffuse radiation (note we add the "LW albedo"). -----------------#
   epsil     = emis.atm
   rf.rldnwd  = ( w01.c0 * ( w01.c1 + w01.c2 * epsil) * log(epsil) + 1. ) * w01.rpfac
   fs.rldnwd  = rf.rldnwd * emis.sfc * rldnwd / tatm.k
   #---------------------------------------------------------------------------------------#

   return(fs.rldnwd)
   #---------------------------------------------------------------------------------------#
}#end sflux.rldnwd
#==========================================================================================#
#==========================================================================================#




#==========================================================================================#
#==========================================================================================#
#     This function finds the entropy flux due to upwelling longwave radiation, following  #
# Q15, who base their formulation on WL10 and W01.  Because skin temperature may be        #
# available, we use that instead of calculating in here (unless skin temperature is not    #
# provided).                                                                               #
#------------------------------------------------------------------------------------------#
fxentropy.rlupwd <<- function(rlupwd,tsfc.k){

   #----- Find derived radiation properties. ----------------------------------------------#
   if (missing(tsfc.k)){
      tsfc.k = sqrt(sqrt(rlupwd / (emis.sfc * stefan)))
   }#end if (missing(tsfc.k))
   #---------------------------------------------------------------------------------------#


   #----- Entropy due to diffuse radiation (note we add the "LW albedo"). -----------------#
   epsil     = emis.sfc
   rf.rlupwd  = ( w01.c0 * ( w01.c1 + w01.c2 * epsil) * log(epsil) + 1. ) * w01.rpfac
   fs.rlupwd  = rf.rlupwd * rlupwd / tsfc.k
   #---------------------------------------------------------------------------------------#

   return(fs.rlupwd)
   #---------------------------------------------------------------------------------------#
}#end sflux.rlupwd
#==========================================================================================#
#==========================================================================================#
