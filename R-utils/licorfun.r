#------------------------------------------------------------------------------------------#
#      This function finds the assimilation rate parameters for a given plant functional   #
# type, meteorological conditions, and limitation.                                         #
#------------------------------------------------------------------------------------------#
assim.params = function(ipft,met,limit){
   aparms=list()
   if (limit == "CLOSED"){
      #----- Closed stomata case, or night time.  These are the same for C3 and C4. -------#
      aparms$rho   = 0.0
      aparms$sigma = 0.0
      aparms$xi    = 0.0
      aparms$tau   = 1.0
      aparms$nu    = - met$leaf.resp
      #---------------------------------------------------------------------------------------#
      #     Open stomata case, so now we distinguish between C3 and C4 as their functional    #
      # forms are different.                                                                  #
      #---------------------------------------------------------------------------------------#
   }else if (pft$pathway[ipft] == 3){
      #------------------------------------------------------------------------------------#
      #     C3 case.  Decide whether this is the light- or Rubisco-limited case.           #
      #------------------------------------------------------------------------------------#
      if (limit == "LIGHT"){
         #---- Light-limited case. --------------------------------------------------------#
         aparms$rho   =  met$alpha * met$par
         aparms$sigma = -met$alpha * met$par * met$compp
         aparms$xi    = 1.0
         aparms$tau   = 2.0 * met$compp
         aparms$nu    = - met$leaf.resp

      }else if (limit == "RUBISCO"){
         #----- Rubisco-limited rate of photosynthesis case. -----------------------------#
         aparms$rho   =  met$vm
         aparms$sigma = -met$vm * met$compp
         aparms$xi    = 1.0
         aparms$tau   = met$kco2 * (1.0 + met$o2 / met$ko2)
         aparms$nu    = - met$leaf.resp

      }else if (limit == "CO2"){
         #----- CO2-limited for low CO2 concentration case. -------------------------------#
         aparms$rho   = 0.0
         aparms$sigma = 0.5 * met$vm
         aparms$xi    = 0.0
         aparms$tau   = 1.0
         aparms$nu    = - met$leaf.resp
      }#end if limit
      #------------------------------------------------------------------------------------#
   }else if (pft$pathway[ipft] == 4){
      #------------------------------------------------------------------------------------#
      #     C4 case.  There are three possibilities, the light-limited, the Rubisco-       #
      # limited, and the CO2-limited cases.                                                #
      #------------------------------------------------------------------------------------#
      if (limit == "LIGHT"){
         #----- Light-limited case. -------------------------------------------------------#
         aparms$rho   = 0.0
         aparms$sigma = met$alpha * met$par
         aparms$xi    = 0.0
         aparms$tau   = 1.0
         aparms$nu    = - met$leaf.resp

      }else if (limit == "RUBISCO"){
         #----- Rubisco-limited rate of photosynthesis case. ------------------------------#
         aparms$rho   = 0.0
         aparms$sigma = met$vm
         aparms$xi    = 0.0
         aparms$tau   = 1.0
         aparms$nu    = - met$leaf.resp

      }else if (limit == "CO2"){
         #----- CO2-limited for low CO2 concentration case. -------------------------------#
         aparms$rho   = klowco2 * met$vm
         aparms$sigma = 0.0
         aparms$xi    = 0.0
         aparms$tau   = 1.0
         aparms$nu    = - met$leaf.resp

      }#end if limit
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Find the minimum radiation that makes Ao = 0. (light compensation point) --------#
   aparms$par.min = ( pft$gamma[ipft] * met$vm * (met$co2 + 2. * met$compp)
                    / (met$co2 - met$compp) ) / met$alpha
   #---------------------------------------------------------------------------------------#

   #------ Find the minimum and the maximum possible Ci. ----------------------------------#
   if ( limit == "CLOSED"){
      aparms$cimin   = 0.
      aparms$cimax   = met$co2
      aparms$ok      = TRUE
   }else{
      aparms$cimin  = 0.
      aparms$cimax  = met$co2
      aparms$ok     = met$par > aparms$par.min
   }#end if
   #---------------------------------------------------------------------------------------#

   return(aparms)
}#end function
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#      This function finds the assimilation rate parameters given the parameters and the   #
# intercellular CO2.                                                                       #
#------------------------------------------------------------------------------------------#
assim.rate = function(ci,aparms){
   #------ Find the assimilation rate. ----------------------------------------------------#
   ao = (aparms$rho * ci + aparms$sigma) / (aparms$xi * ci + aparms$tau) + aparms$nu
   return(ao)
   #---------------------------------------------------------------------------------------#
}#end function
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#      This function finds the assimilation rate parameters given the parameters and the   #
# intercellular CO2.                                                                       #
#------------------------------------------------------------------------------------------#
stom.cond = function(ao,ci,met){
   #------ Find the assimilation rate. ----------------------------------------------------#
   if (ao < 0.){
      gsw = pft$b[ipft]
   }else{
      gsw = pft$m[ipft] * ao / ( (met$co2 - met$compp) * met$rterm) + pft$b[ipft]
   }#end if
   return(gsw)
   #---------------------------------------------------------------------------------------#
}#end function
#------------------------------------------------------------------------------------------#





#------------------------------------------------------------------------------------------#
#     Function whose root gives Ci.                                                        #
#------------------------------------------------------------------------------------------#
cifun = function(ci,met,ipft,aparms){
   ao  = assim.rate(ci,aparms)
   gsw = stom.cond(ao,ci,met)
   fun = met$co2 - ao / (gsw.2.gsc * gsw) - ci
   return(fun)
}#end function
#------------------------------------------------------------------------------------------#





#------------------------------------------------------------------------------------------#
#     Function whose root gives Ci.                                                        #
#------------------------------------------------------------------------------------------#
cisolver = function(met,ipft,aparms){
   #----- Find some auxiliary variables. --------------------------------------------------#
   bbb   = gsw.2.gsc * pft$b[ipft]
   mu    = gsw.2.gsc * pft$m[ipft] / ( (met$co2 - met$compp) * met$rterm)
   rho   = aparms$rho
   sigma = aparms$sigma
   xi    = aparms$xi
   tau   = aparms$tau
   nu    = aparms$nu
   cs    = met$co2

   #----- Find the quadratic terms. -------------------------------------------------------#
   aquad = mu * rho + mu * nu * xi + bbb * xi
   bquad = ( mu * sigma + mu * nu * tau + bbb * tau - mu * rho * cs - mu * nu * xi * cs
           - bbb * xi * cs + rho + nu * xi )
   cquad = nu * tau + sigma - mu * sigma * cs - mu * nu * tau * cs - tau * cs * bbb


   #----- Check whether this is really quadratic or not. ----------------------------------#
   if (aquad == 0){
      ci    = - cquad / bquad
   }else{
      delta = bquad * bquad - 4. * aquad * cquad
      if (delta > 0.){
         ci1 = - 0.5 * ( bquad + sqrt(delta)) / aquad
         ci2 = - 0.5 * ( bquad - sqrt(delta)) / aquad
         ci  = max(ci1,ci2)
      }else if(delta == 0.){
         ci = - 0.5 * bquad / aquad
      }else{
         #----- Delta is negative, crash! -------------------------------------------------#
         print(paste("MU    =",mu    ))
         print(paste("BETA  =",bbb   ))
         print(paste("RHO   =",rho   ))
         print(paste("SIGMA =",sigma ))
         print(paste("XI    =",xi    ))
         print(paste("TAU   =",tau   ))
         print(paste("NU    =",nu    ))
         print(paste("CS    =",cs    ))
         print(paste("AQUAD =",aquad ))
         print(paste("BQUAD =",bquad ))
         print(paste("CQUAD =",cquad ))
         print(paste("DELTA =",delta ))
         stop("Failed finding a reasonable solution")
      }#end if
   }#end if
   return(ci)
}#end function
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#      This function finds five parameters, namely:                                        #
#  - Vm        : the photosynthetic capacity.                                              #
#  - leaf.resp : the leaf respiration rate.                                                #
#  - compp     : the CO2 compensation point excluding respiration                          #
#  - kco2      : Michaelis-Mentel constant for CO2                                         #
#  - ko2       : Michaelis-Mentel constant for O2                                          #
#------------------------------------------------------------------------------------------#
photo.params = function(met,ipft,iphoto){

   #---------------------------------------------------------------------------------------#
   #     Decide which method to apply.  Currently we only solve differently if we are      #
   # using Collatz, otherwise we use the same Vm.                                          #
   #---------------------------------------------------------------------------------------#
   if (iphoto %in% 0:1){
      #------------------------------------------------------------------------------------#
      #     Use ED approach.  First we copy the parameters that will be used for scratch   #
      # local variables, then apply the same functional form.                              #
      #------------------------------------------------------------------------------------#
      if (iphoto == 0){
         #----- Foley et al. (1996) as is. ------------------------------------------------#
         kco2.ref  = kco2.ref.ibis
         kco2.hor  = kco2.hor.ibis
         ko2.ref   = ko2.ref.ibis
         ko2.hor   = ko2.hor.ibis
         compp.ref = compp.ref.ibis
         compp.hor = compp.hor.ibis

      }else if (iphoto == 1){
         #----- Foley et al. (1996) for KCO2 and KO2, find Gamma* like in CLM. ------------#
         kco2.ref  = kco2.ref.ibis
         kco2.hor  = kco2.hor.ibis
         ko2.ref   = ko2.ref.ibis
         ko2.hor   = ko2.hor.ibis

         compp.ref = kookc * met$o2 * kco2.ref / ko2.ref
         compp.hor = kco2.hor - ko2.hor

      }#end if

      lnexplow  = pft$vm.decay.e[ipft] * (pft$vm.low.temp[ipft]  - met$temp)
      tlow.fun  = 1.0 +  exp(lnexplow)
      lnexphigh = pft$vm.decay.e[ipft] * (met$temp - pft$vm.high.temp[ipft])
      thigh.fun = 1.0 + exp(lnexphigh)
      vm        = ( arrhenius(met$temp,pft$vm0[ipft],pft$vm.hor[ipft]) 
                  / (tlow.fun * thigh.fun) )

      lnexplow  = pft$lr.decay.e[ipft] * (pft$lr.low.temp[ipft]  - met$temp)
      tlow.fun  = 1.0 +  exp(lnexplow)
      lnexphigh = pft$lr.decay.e[ipft] * (met$temp - pft$lr.high.temp[ipft])
      thigh.fun = 1.0 + exp(lnexphigh)
      leaf.resp = ( arrhenius(met$temp,pft$gamma[ipft]*pft$vm0[ipft],pft$lr.hor[ipft]) 
                  / (tlow.fun * thigh.fun) )

      compp     = arrhenius(met$temp,compp.ref,compp.hor)
      kco2      = arrhenius(met$temp,kco2.ref,kco2.hor)
      ko2       = arrhenius(met$temp,ko2.ref,ko2.hor)

   }else if(iphoto %in% 2:3){
      #------------------------------------------------------------------------------------#
      #     Use the Collatz et al. scheme.  Here the only difference is that in option 5   #
      # we use the compensation point as defined in one of the Collatz papers, whereas in  #
      # option 6 we use the CLM approach, which is derived from Farquhar et al. (1980)     #
      # assuming that the ratio between the tunover for oxylase and carboxylase is         #
      # constant with temperature.                                                         #
      #------------------------------------------------------------------------------------#
      if (iphoto == 2){
         kco2.ref   = kco2.ref.coll
         kco2.base  = kco2.base.coll
         ko2.ref    = ko2.ref.coll
         ko2.base   = ko2.base.coll
         vm.ref     = pft$vm0[ipft]
         vm.base    = pft$vm.base[ipft]
         compp.ref  = compp.ref.coll
         compp.base = compp.base.coll
         lr.ref     = pft$gamma[ipft] * pft$vm0[ipft]
         lr.base    = pft$lr.base[ipft]
      }else if (iphoto == 3){
         kco2.ref   = kco2.ref.coll
         kco2.base  = kco2.base.coll
         ko2.ref    = ko2.ref.coll
         ko2.base   = ko2.base.coll
         compp.ref  = kookc * met$o2 * kco2.ref / (2.0 * ko2.ref)
         compp.base = kco2.base / ko2.base
         vm.ref     = pft$vm0[ipft]
         vm.base    = pft$vm.base[ipft]
         lr.ref     = pft$gamma[ipft] * pft$vm0[ipft]
         lr.base    = pft$lr.base[ipft]
      }#end if

      vm.nocorr  = collatz(met$temp,vm.ref,vm.base)
      lr.nocorr  = collatz(met$temp,lr.ref,lr.base)

      #------------------------------------------------------------------------------------#
      #      Apply the low/high temperature corrections for both Vm and Leaf respiration.  #
      # Here we must check whether we are solving C3 or C4, because the correction has a   #
      # different functional form for Vm.                                                  #
      #------------------------------------------------------------------------------------#
      if (c91vmcorr){
         lnexplow  = pft$vm.decay.e[ipft] * (pft$vm.low.temp[ipft]  - met$temp)
         tlow.fun  = 1.0 +  exp(lnexplow)
         thigh.fun = 1. + exp( (-pft$vm.decay.a[ipft] + pft$vm.decay.b[ipft] * met$temp)
                                  / (rmol * met$temp) )
         vm        = vm.nocorr / (tlow.fun * thigh.fun)

      }else{
         lnexplow  = pft$vm.decay.e[ipft] * (pft$vm.low.temp[ipft]  - met$temp)
         tlow.fun  = 1.0 +  exp(lnexplow)
         lnexphigh = pft$vm.decay.e[ipft] * (met$temp - pft$vm.high.temp[ipft])
         thigh.fun = 1.0 + exp(lnexphigh)
         vm        = vm.nocorr / (tlow.fun * thigh.fun)
      }#end if
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    Leaf respiration.  The correction is the same for both C3 and C4.               #
      #------------------------------------------------------------------------------------#
      lnexplow  = pft$lr.decay.e[ipft] * (pft$lr.low.temp[ipft]  - met$temp)
      tlow.fun  = 1.0 +  exp(lnexplow)
      lnexphigh = pft$lr.decay.e[ipft] * (met$temp - pft$lr.high.temp[ipft])
      thigh.fun = 1.0 + exp(lnexphigh)
      leaf.resp = lr.nocorr / (tlow.fun * thigh.fun)
      #------------------------------------------------------------------------------------#


      #----- Find the compensation point and the Michaelis-Mentel constants. --------------#
      compp     = collatz(met$temp,compp.ref,compp.base)
      kco2      = collatz(met$temp,kco2.ref,kco2.base)
      ko2       = collatz(met$temp,ko2.ref,ko2.base)

   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Before we leave, just make sure to set the CO2 compensation point and the         #
   # Michaelis-Mentel constants to zero in case this is a C4 plant.                        #
   #---------------------------------------------------------------------------------------#
   if (pft$pathway[ipft] == 4){
      compp = 0. * compp
      kco2  = 0. * kco2
      ko2   = 0. * ko2
   }#end if
   #---------------------------------------------------------------------------------------#


   photo = list (vm = vm, leaf.resp=leaf.resp,compp=compp,kco2=kco2,ko2=ko2)
   return(photo)
}#end function
#------------------------------------------------------------------------------------------#
