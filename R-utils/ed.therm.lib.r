#==========================================================================================#
#==========================================================================================#
# This function calculates the individual heat capacity in (J/K/pl) of the cohort          #
# leaf biomass.  This function is primarily used to calculate leaf temperatures            #
# based on leaf energy.                                                                    #
# These methods follow the ways of Gu et al. 2007, with the only difference that for       #
# non-green biomass we dropped the temperature dependence and assumed T=T3ple, just to     #
# make it simpler.                                                                         #
# Reference:                                                                               #
#                                                                                          #
# Gu, L., T. Meyers, S. G. Pallardy, 2007: Influences of biomass heat and biochemical      #
#      energy storages on the land surface fluxes and radiative temperature.               #
#      J. Geophys. Res., v. 112, doi: 10.1029/2006JD007425.                                #
#------------------------------------------------------------------------------------------#
calc.ind.veg.hcap <<- function(dbh,ipft,hout=c("leaf","wood","veg"),use.crit=TRUE){

   #---- Standardise heat capacity output. ------------------------------------------------#
   hout = match.arg(hout)
   #---------------------------------------------------------------------------------------#


   #----- Find biomass of all pools. ------------------------------------------------------#
   hgt    = dbh2h  (dbh=dbh,ipft=ipft)
   bleaf  = size2bl(dbh=dbh,hgt=hgt,ipft=ipft)
   if (hout %in% c("wood","veg")){
      bdeada = pft$brf.wd[ipft] * pft$agf.bs[ipft] * size2bd(dbh=dbh,hgt=hgt,ipft=ipft)
      bsapa  = pft$brf.wd[ipft] * pft$agf.bs[ipft] * bleaf * pft$qsw  [ipft] * hgt
      bbarka = pft$brf.wd[ipft] * pft$agf.bs[ipft] * bleaf * pft$qbark[ipft] * hgt
      bwooda = bdeada + bsapa
   }#end if (hout %in% c("wood","veg"))
   #---------------------------------------------------------------------------------------#


   #----- Decide which heat capacity to return. -------------------------------------------#
   if (hout %in% "leaf"){
      ans = C2B *   bleaf  * pft$c.leaf[ipft] * (1. + pft$qwatdry.leaf[ipft])
   }else if (hout %in% "wood"){
      ans = C2B * ( bwooda * pft$c.wood[ipft] * (1. + pft$qwatdry.wood[ipft])
                  + bbarka * pft$c.bark[ipft] * (1. + pft$qwatdry.bark[ipft]) )
   }else if (hout %in% "veg" ){
      ans = C2B * ( bleaf  * pft$c.leaf[ipft] * (1. + pft$qwatdry.leaf[ipft])
                  + bwooda * pft$c.wood[ipft] * (1. + pft$qwatdry.wood[ipft])
                  + bbarka * pft$c.bark[ipft] * (1. + pft$qwatdry.bark[ipft]) )
   }#end if (hout %in% "leaf")
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end calc.ind.veg.hcap
#==========================================================================================#
#==========================================================================================#
