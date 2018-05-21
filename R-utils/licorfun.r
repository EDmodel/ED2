#==========================================================================================#
#==========================================================================================#
#      Some parameters for the solver.                                                     #
#------------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #    This is the tolerance for iterative methods.  Some of the functions are quite flat #
   # so it is a good idea to use a somewhat more strict tolerance than the ones used in    #
   # therm_lib8.                                                                           #
   #---------------------------------------------------------------------------------------#
   tolerfl          <<- 1.e-10
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Since this is more strict than therm_lib, the maximum number of attempts for the   #
   # Regula Falsi method should be also increased.                                         #
   #---------------------------------------------------------------------------------------#
   maxfpofl         <<- 320
   #---------------------------------------------------------------------------------------#


#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     This subroutine computes the maximum capacity of Rubisco to perform the carboxylase  #
# function at a certain temperature, given the PFT and phenology properties, and, in case  #
# of light controlled phenology, the leaf life span and the average reference value of the #
# Vm function (Vm0).  The output variables are stored in the aparms structure, as they are #
# going to be used in other sub-routines. Both Vm and the reference value Vm0 have units   #
# of umol/m2/s.                                                                            #
#     Compute the photosynthesis and leaf respiration parameters that depend on temper-    #
# ature, according to the parameters defined in the "thispft" structure and the functional #
# form chosen by the user.  The variables that are defined there are:                      #
# - alpha     - the quantum yield, which may be a function of temperature.                 #
# - Vm        - the maximum capacity of Rubisco to perform the carboxylase function.       #
# - leaf_resp - the leaf respiration.                                                      #
# - compp     - the CO2 compensation point for gross photosynthesis (Gamma*)               #
# - kco2      - Michaelis-Mentel coefficient for CO2                                       #
# - ko2       - Michaelis-Mentel coefficient for O2.                                       #
#------------------------------------------------------------------------------------------#
comp.photo.tempfun <<- function(thispft,met,quantum.t){

   #------ Initialise aparms data frame. --------------------------------------------------#
   aparms    = data.frame(matrix(nrow=1,ncol=0))
   #---------------------------------------------------------------------------------------#

   #----- Assume leaves are fully flushed. ------------------------------------------------#
   greenness = 1.0
   #---------------------------------------------------------------------------------------#


   #----- Quantum yield.  Check which method to use. --------------------------------------#
   if (quantum.t && thispft$photo.pathway == 3){
      tempc        = max(0,met$leaf.temp - t00)
      aparms$alpha = 0.081 - 0.000053*(met$temp-t00) - 0.000019*(met$temp-t00)^2
   }else{
      aparms$alpha = thispft$alpha0
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Find Vm using the Collatz equation, with no correction. -------------------------#
   vm.nocorr = greenness * collatz(met$leaf.temp,thispft$vm0,thispft$vm.q10)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Compute the functions that will control the Vm function for high temper-           #
   # ature.  In order to avoid floating point exceptions, we check whether the             #
   # temperature will make the exponential too small or too large.                         #
   #---------------------------------------------------------------------------------------#
   #----- Low temperature. ----------------------------------------------------------------#
   lnexplow  = thispft$vm.decay.ecold * (thispft$vm.low.temp - met$leaf.temp)
   lnexplow  = max(lnexp.min,min(lnexp.max,lnexplow))
   tlow.fun  = 1.0 +  exp(lnexplow)
   #----- High temperature. ---------------------------------------------------------------#
   lnexphigh = thispft$vm.decay.ehot  * (met$leaf.temp - thispft$vm.high.temp)
   lnexphigh = max(lnexp.min,min(lnexp.max,lnexphigh))
   thigh.fun = 1.0 + exp(lnexphigh)
   #---------------------------------------------------------------------------------------#

   #------ Correct Vm. --------------------------------------------------------------------#
   aparms$vm = vm.nocorr / (tlow.fun * thigh.fun)
   #---------------------------------------------------------------------------------------#


   #----- Find Rd using the Collatz equation, with no correction. -------------------------#
   rd.nocorr = greenness * collatz(met$leaf.temp,thispft$rd0,thispft$rd.q10)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Compute the functions that will control the rd function for low and high           #
   # temperature.  In order to avoid floating point exceptions, we check whether the       #
   # temperature will make the exponential too large or too small.                         #
   #---------------------------------------------------------------------------------------#
   #----- Low temperature. ----------------------------------------------------------------#
   lnexplow  = thispft$rd.decay.ecold * (thispft$rd.low.temp  - met$leaf.temp)
   lnexplow  = max(lnexp.min,min(lnexp.max,lnexplow))
   tlow.fun  = 1.0 +  exp(lnexplow)
   #----- High temperature. ---------------------------------------------------------------#
   lnexphigh = thispft$rd.decay.ehot  * (met$leaf.temp - thispft$rd.high.temp)
   lnexphigh = max(lnexp.min,min(lnexp.max,lnexphigh))
   thigh.fun = 1.0 + exp(lnexphigh)
   #---------------------------------------------------------------------------------------#

   #------ Correct Rd. --------------------------------------------------------------------#
   aparms$leaf.resp = rd.nocorr / (tlow.fun * thigh.fun)
   #---------------------------------------------------------------------------------------#


   #----- Find Jm using the Collatz equation, with no correction. -------------------------#
   jm.nocorr = greenness * collatz(met$leaf.temp,thispft$jm0,thispft$jm.q10)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Compute the functions that will control the rd function for low and high           #
   # temperature.  In order to avoid floating point exceptions, we check whether the       #
   # temperature will make the exponential too large or too small.                         #
   #---------------------------------------------------------------------------------------#
   #----- Low temperature. ----------------------------------------------------------------#
   lnexplow  = thispft$jm.decay.ecold * (thispft$jm.low.temp  - met$leaf.temp)
   lnexplow  = max(lnexp.min,min(lnexp.max,lnexplow))
   tlow.fun  = 1.0 +  exp(lnexplow)
   #----- High temperature. ---------------------------------------------------------------#
   lnexphigh = thispft$jm.decay.ehot  * (met$leaf.temp - thispft$jm.high.temp)
   lnexphigh = max(lnexp.min,min(lnexp.max,lnexphigh))
   thigh.fun = 1.0 + exp(lnexphigh)
   #---------------------------------------------------------------------------------------#

   #------ Correct Rd. --------------------------------------------------------------------#
   aparms$jm = jm.nocorr / (tlow.fun * thigh.fun)
   #---------------------------------------------------------------------------------------#



   #------ Find Jmax and Tpmax using Vcmax and the Jm0:Vm0 and Tp0:Vm0 ratios. ------------#
   aparms$tpm = aparms$vm * thispft$tpm0 / thispft$vm0
   #---------------------------------------------------------------------------------------#



   #------ Find the maximum electron transport ratio. -------------------------------------#
   if (thispft$curvpar > 1.){
      ipsII.max = ( aparms$jm  * ( 2.0 * thispft$curvpar - 1.
                                 - 2.0 * sqrt(thispft$curvpar * (thispft$curvpar-1.0) ) ) )
   }else{
      ipsII.max = discard
   }#end if (thispft$curvpar > 1.)
   ipsII = min(0.5 * thispft$phi.psII * met$par,ipsII.max)
   aterm  = thispft$curvpar
   bterm  = - ( ipsII + aparms$jm)
   cterm  = ipsII * aparms$jm
   if (aterm == 0.){
      aparms$j = - cterm / bterm
   }else{
      #------------------------------------------------------------------------------------#
      #     Quadratic function.                                                            #
      #------------------------------------------------------------------------------------#
      discr  = bterm*bterm - 4.0*aterm*cterm
      if (discr == 0.){
         aparms$j = - bterm / (2.0 * aterm)
      }else if (discr > 0.){
         jroot1 = ( - bterm - sqrt(discr) ) / (2. * aterm)
         jroot2 = ( - bterm + sqrt(discr) ) / (2. * aterm)
         if (jroot1 <= jroot2){
            aparms$j = jroot1
         }else{
            aparms$j = jroot2
         }#end if (jroot1 <= jroot2)
      }else{
         #------ Negative discriminant, assume J = Jcrit. ---------------------------------#
         cat0(" Negative discriminant")
         browser()
         #---------------------------------------------------------------------------------#
      }#end if
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    If this is a C3 plant, find the compensation point, and the Michaelis-Mentel       #
   # constants for CO2 and O2.  Otherwise, assign zeroes for compensation point and the    #
   # Michaelis-Mentel constant for CO2.  The oxygen one should have no impact,         #
   # but we always assign it to avoid divisions by zero.                                   #
   #---------------------------------------------------------------------------------------#
   if (thispft$photo.pathway == 3){
      aparms$compp = collatz(met$leaf.temp,compp.refval,compp.q10)
      aparms$kco2  = collatz(met$leaf.temp,kco2.refval ,kco2.q10 )
   }else{
      aparms$compp = 0.0
      aparms$kco2  = 0.0
   }#end if (thispft$photo.pathway == 3)
   aparms$ko2   = collatz(met$leaf.temp,ko2.refval,ko2.q10)
   #---------------------------------------------------------------------------------------#



   #----- Define the true CO2 compensation point. -----------------------------------------#
   grsp             = thispft$rd0 / thispft$vm0
   aparms$compptrue = ( ( aparms$compp + aparms$kco2 * (1. + met$can.o2/aparms$ko2) * grsp)
                      / (1. - grsp) )
   #---------------------------------------------------------------------------------------#

   return(aparms)
}#end comp.photo.tempfun
#==========================================================================================#
#==========================================================================================#




#==========================================================================================#
#==========================================================================================#
#      This function finds the assimilation rate parameters for a given plant functional   #
# type, meteorological conditions, and limitation.                                         #
#------------------------------------------------------------------------------------------#
set.co2.demand.params <<- function(thispft,aparms,met,limit){

   if (limit %in% "CLOSED"){
      #----- Closed stomata case, or night time.  These are the same for C3 and C4. -------#
      aparms$rho   = 0.0
      aparms$sigma = 0.0
      aparms$xi    = 0.0
      aparms$tau   = 1.0
      aparms$nu    = - aparms$leaf.resp
      #---------------------------------------------------------------------------------------#
      #     Open stomata case, so now we distinguish between C3 and C4 as their functional    #
      # forms are different.                                                                  #
      #---------------------------------------------------------------------------------------#
   }else if (thispft$photo.pathway %in% 3){
      #------------------------------------------------------------------------------------#
      #     C3 case.  Decide whether this is the light- or RuBP-saturated case.            #
      #------------------------------------------------------------------------------------#
      if (limit %in% "LIGHT"){
         #---- Light-limited case. --------------------------------------------------------#
         aparms$rho   =  aparms$j
         aparms$sigma = -aparms$j * aparms$compp
         aparms$xi    = 4.0
         aparms$tau   = 8.0 * aparms$compp
         aparms$nu    = - aparms$leaf.resp

      }else if (limit %in% "RUBP_SAT" ){
         #----- RuBP-saturated rate of photosynthesis case. -------------------------------#
         aparms$rho   =  aparms$vm
         aparms$sigma = -aparms$vm * aparms$compp
         aparms$xi    = 1.0
         aparms$tau   = aparms$kco2 * (1.0 + met$can.o2 / aparms$ko2)
         aparms$nu    = - aparms$leaf.resp

      }else if (limit %in% "TPU"){
         #----- Triose phosphate utilisation rate limitation case. ------------------------#
         aparms$rho   =  0.0
         aparms$sigma =  3.0 * aparms$tp 
         aparms$xi    =  0.0
         aparms$tau   =  1.0
         aparms$nu    = -aparms$leaf.resp
      }else if (limit %in% "CO2"){
         stop(" CO2 limitation is not allowed for C3 plants.")
      }#end if limit
      #------------------------------------------------------------------------------------#
   }else if (thispft$photo.pathway %in% 4){
      #------------------------------------------------------------------------------------#
      #     C4 case.  There are three possibilities, the light-limited, the RuBP-          #
      # saturated, and the CO2-limited cases.                                              #
      #------------------------------------------------------------------------------------#
      if (limit %in% "LIGHT"){
         #----- Light-limited case. -------------------------------------------------------#
         aparms$rho   = 0.0
         aparms$sigma = aparms$j
         aparms$xi    = 0.0
         aparms$tau   = 4.0
         aparms$nu    = - aparms$leaf.resp

      }else if (limit %in% "RUBP_SAT"){
         #----- RuBP-saturated rate of photosynthesis case. -------------------------------#
         aparms$rho   = 0.0
         aparms$sigma = aparms$vm
         aparms$xi    = 0.0
         aparms$tau   = 1.0
         aparms$nu    = - aparms$leaf.resp
      }else if (limit %in% "TPU"){
         stop(" Triose Phosphate limitation is not allowed for C4 plants.")
      }else if (limit %in% "CO2"){
         #----- CO2-limited for low CO2 concentration case. -------------------------------#
         aparms$rho   = klowco2 * aparms$vm
         aparms$sigma = 0.0
         aparms$xi    = 0.0
         aparms$tau   = 1.0
         aparms$nu    = - aparms$leaf.resp
      }#end if limit
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#

   return(aparms)
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#    Find the minimum amount of radiation for which we will consider daytime.  This        #
# cannot be a constant because it depends on the PFT and the environmental conditions.     #
#------------------------------------------------------------------------------------------#
find.twilight.min <<- function(aparms,thispft,met){

   #----- Find the electron transport that makes Ao = 0. ----------------------------------#
   j0 = ( 4.0 * aparms$leaf.resp * (met$can.co2 + 2. * aparms$compp)
                                 / (met$can.co2 -      aparms$compp) )
   #---------------------------------------------------------------------------------------#


   #----- Find the radiation associated with  that makes Ao = 0. (light compensation point) --------#
   if (j0 < aparms$jm){
      ans  = ( ( 2.0 * (aparms$jm * j0 - thispft$curvpar * j0 * j0) )
             / ( thispft$phi.psII * ( aparms$jm - j0 )              ) )
   }else{
      ans  = max( 0.4 * solar * Watts.2.Ein,met$par)
   }#end if
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end function find.twilight.min
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This sub-routine computes the minimum and maximum intercellular carbon dioxide       #
# concentration that we should seek the solution.  Both the function from which we seek    #
# the root and the stomatal conductance function have singularities, so we don't want the  #
# intercellular carbon dioxide to cross these singularities because the root-finding       #
# method assumes continuity.                                                               #
#------------------------------------------------------------------------------------------#
find.lint.co2.bounds <<- function(met,thispft,aparms){
   #---------------------------------------------------------------------------------------#
   # First case: This check will find when Aopen goes to 0., which causes a singularity    #
   # in the function of which we are looking for a root.                                   #
   #---------------------------------------------------------------------------------------#
   if (aparms$rho == 0. && aparms$xi == 0.){
      ciAo = c34smin.lint.co2
   }else{
      ciAo = ( - (aparms$tau * aparms$nu + aparms$sigma)
               / (aparms$xi  * aparms$nu + aparms$rho  ) )
   }#end if
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   # Second case: This check finds which ci causes the terms [gbc (ca -ci) - Aopen] to     #
   # be 0.   This will cause a singularity in the gsw function.                            #
   #---------------------------------------------------------------------------------------#
   #----- 1. Define the coefficients for the quadratic equation. --------------------------#
   aquad = met$blyr.cond.co2 * aparms$xi
   bquad = ( aparms$xi * (aparms$nu - met$blyr.cond.co2 * met$can.co2 )
           + met$blyr.cond.co2 * aparms$tau + aparms$rho )
   cquad = aparms$tau * (aparms$nu - met$blyr.cond.co2 * met$can.co2 ) + aparms$sigma
   #----- 2. Decide whether this is a true quadratic case or not. -------------------------#
   if (aquad != 0.0){
      #------------------------------------------------------------------------------------#
      #     This is a true quadratic case, the first step is to find the discriminant.     #
      #------------------------------------------------------------------------------------#
      discr = bquad * bquad - 4.0 * aquad * cquad
      if (discr == 0.0){
         #---------------------------------------------------------------------------------#
         #      Discriminant is zero, both roots are the same.  We save only one, and      #
         # make the other negative, which will make the guess discarded.                   #
         #---------------------------------------------------------------------------------#
         ciroot1      = - bquad / (2.0 * aquad)
         ciroot2      = - discard
      }else if (discr > 0.0){
         ciroot1 = (- bquad + sqrt(discr)) / (2.0 * aquad)
         ciroot2 = (- bquad - sqrt(discr)) / (2.0 * aquad)
      }else{
         #----- Discriminant is negative.  Impossible to solve. ---------------------------#
         ciroot1      = - discard
         ciroot2      = - discard
      }# end if
   }else{
      #------------------------------------------------------------------------------------#
      #    This is a linear case, the xi term is zero.  There is only one number that      #
      # works for this case.                                                               #
      #------------------------------------------------------------------------------------#
      ciroot1      = - cquad / bquad
      ciroot2      = - discard
      #----- Not used, just for the debugging process. ------------------------------------#
      discr        = bquad * bquad
   }# end if
   #---------------------------------------------------------------------------------------#
   #     Save the largest of the values.  In case both were discarded, we switch it to     #
   # the positive discard so this will never be chosen.                                    #
   #---------------------------------------------------------------------------------------#
   cigsw=max(ciroot1, ciroot2)
   if (cigsw == -discard) cigsw = discard
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   # Third case: This will find the intercellular CO2 mixing ratio that causes the qi-qs   #
   # term to be equal to -D0, which also creates a singularity in the stomatal             #
   # conductance.                                                                          #
   #---------------------------------------------------------------------------------------#
   #----- 1. Find some auxiliary variables. -----------------------------------------------#
   xtmp = met$blyr.cond.h2o * ( met$can.shv - met$lint.shv - thispft$d0 ) / thispft$d0
   ytmp = met$blyr.cond.co2 + xtmp * gbw.2.gbc
   ztmp = xtmp * gbw.2.gbc * met$blyr.cond.co2
   wtmp = ztmp * met$can.co2 - ytmp * aparms$nu
   #----- 2. Define the coefficients for the quadratic equation. --------------------------#
   aquad = ztmp * aparms$xi
   bquad = ytmp * aparms$rho - aparms$xi * wtmp + ztmp * aparms$tau
   cquad = - aparms$tau * wtmp + ytmp * aparms$sigma
   #----- 3. Decide whether this is a true quadratic case or not. -------------------------#
   if (aquad != 0.0){
      #------------------------------------------------------------------------------------#
      #     This is a true quadratic case, the first step is to find the discriminant.     #
      #------------------------------------------------------------------------------------#
      discr = bquad * bquad - 4.0 * aquad * cquad
      if (discr == 0.0){
         #---------------------------------------------------------------------------------#
         #      Discriminant is zero, both roots are the same.  We save only one, and      #
         # make the other negative, which will make the guess discarded.                   #
         #---------------------------------------------------------------------------------#
         ciroot1 = - bquad / (2.0 * aquad)
         ciroot2 = -discard
      }else if (discr > 0.0){
         ciroot1 = (- bquad + sqrt(discr)) / (2.0 * aquad)
         ciroot2 = (- bquad - sqrt(discr)) / (2.0 * aquad)
      }else{
         #----- Discriminant is negative.  Impossible to solve. ---------------------------#
         ciroot1      = -discard
         ciroot2      = -discard
      }#end if
   }else{
      #------------------------------------------------------------------------------------#
      #    This is a linear case, the xi term is zero.  There is only one number           #
      # that works for this case.                                                          #
      #------------------------------------------------------------------------------------#
      ciroot1 = - cquad / bquad
      ciroot2 = -discard
      #----- Not used, just for the debugging process. ------------------------------------#
      discr   = bquad * bquad
   }# end if
   #---------------------------------------------------------------------------------------#
   #     Save the largest of the values.  In case both were discarded, we switch it to     #
   # the positive discard so this will never be chosen.                                    #
   #---------------------------------------------------------------------------------------#
   ciQ=max(ciroot1, ciroot2)
   if (ciQ == -discard) ciQ = discard
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Make sure that there is no singularity within Cimin and Cimax.  From previous     #
   # tests, we know that cimin is the one associated with the case in which Aopen goes     #
   # to zero, and the maximum is the minimum between gsw case, q case, and the canopy      #
   # air CO2.                                                                              #
   #---------------------------------------------------------------------------------------#
   ciroot1 = max(c34smin.lint.co2,ciAo)
   ciroot2 = min(ciQ,cigsw,met$can.co2)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     The actual bounds are slightly squeezed so the edges will not be at the the       #
   # singularities.                                                                        #
   #---------------------------------------------------------------------------------------#
   cimin = ciroot1 + 1.e-5 * (ciroot2 - ciroot1)
   cimax = ciroot2 - 1.e-5 * (ciroot2 - ciroot1)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    In the last part we make sure that the two guesses are positively defined and      #
   # make sense.  If both bounds are non-positive, or if the minimum ci is greater than    #
   # the maximum cimax, the bounds make no sense, and the solution is not bounded.         #
   #---------------------------------------------------------------------------------------#
   if (cimin >= cimax){
      cimin   = c34smin.lint.co2
      cimax   = c34smin.lint.co2
      bounded = FALSE
   }else if (cimin <= 0.0 && cimax <= c34smin.lint.co2){
      cimin   = c34smin.lint.co2
      cimax   = c34smin.lint.co2
      bounded = FALSE
   }else if (cimin <= 0.0){
      cimin   = c34smin.lint.co2
      bounded = TRUE
   }else{
      bounded = TRUE
   }# end if

   ans = data.frame( cimin    = cimin
                   , cimax    = cimax
                   , bounded  = bounded
                   , ci.shv   = ciQ
                   , ci.assim = ciAo
                   , ci.gsw   = cigsw
                   )#end list
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This sub-routine will solve the case in which both the carbon demand and the         #
# stomatal conductance of water are functions of the intercellular CO2.  This function     #
# can be used for the simpler cases too, but it's not advisable because this is based      #
# on iterative methods, which makes the solution slower.                                   #
#     The iterative method is designed to use Newton's method as the default, and this     #
# should take care of most cases.  In case Newton's method fails, it will fall back to     #
# the modified Regula Falsi method (Illinois) and look for guesses with opposite sign.     #
# If the method fails finding the pair, it means that there is no viable solution          #
# within this range, so the method quits and return the error message.                     #
#------------------------------------------------------------------------------------------#
solve.iterative.case <<- function(met,thispft,aparms){



   #---------------------------------------------------------------------------------------#
   #   1. Initialise the success flag as true.  In case we have trouble solving this case, #
   #      we switch the flag to false before we quit.                                      #
   #---------------------------------------------------------------------------------------#
   answer               = aparms
   answer$lsfc.shv      = NA_real_
   answer$lsfc.co2      = NA_real_
   answer$lint.co2      = NA_real_
   answer$co2.demand    = NA_real_
   answer$stom.cond.h2o = NA_real_
   answer$stom.cond.co2 = NA_real_
   answer$transp        = NA_real_
   answer$cimin         = NA_real_
   answer$cimax         = NA_real_
   answer$ci.assim      = NA_real_
   answer$ci.gsw        = NA_real_
   answer$ci.shv        = NA_real_
   answer$success       = FALSE
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Initialise the convergence flag.  Here we start realistic, ops, I mean,          #
   # pessimistic, and assume that we are failing.  This will be switched to true only if a #
   # real answer is found, so in case we quit the sub-routine due to impossible solution,  #
   # the result will be failure.                                                           #
   #---------------------------------------------------------------------------------------#
   converged = FALSE
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Determine the minimum and maximum intercellular CO2 that will still produce  a     #
   # positive and bounded conductance, which should be above the cuticular conductance     #
   # (otherwise there is no reason to keep stomata opened.                                 #
   #---------------------------------------------------------------------------------------#
   co2.bounds = find.lint.co2.bounds(met,thispft,aparms)
   bounded    = co2.bounds$bounded
   cimin      = co2.bounds$cimin
   cimax      = co2.bounds$cimax
   ci.assim   = co2.bounds$ci.assim
   ci.gsw     = co2.bounds$ci.gsw
   ci.shv     = co2.bounds$ci.shv
   if (bounded){
      #------------------------------------------------------------------------------------#
      #     We have found bounds, find put the first guess in the middle and find the      #
      # function evaluation and derivative of this guess.                                  #
      #------------------------------------------------------------------------------------#
      cia     = sqrt(cimin*cimax)
      iter    = iter.solver.step(met,thispft,aparms,cia,newton=TRUE)
      funa    = iter$fun
      dfundci = iter$dfundci
   }else{
      #----- No reasonable bound was found, we don't even bother solving this case. -------#
      return(answer)
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Copy to the new guess just in case it fails at the first iteration --------------#
   ciz = cia
   fun = funa
   #---------------------------------------------------------------------------------------#


   if (fun == 0.0){
      #----- We have actually hit the jackpot, the answer is ciz. -------------------------#
      ci        = ciz
      converged = TRUE
   }#end if (fun == 0.0)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Enter Newton's method loop in case we haven't found the answer.                   #
   #---------------------------------------------------------------------------------------#
   if (! converged){
      iterate = TRUE
      itn     = 0
      while ((itn < floor(maxfpofl/6)) && iterate && (abs(dfundci) < tolerfl)){
         #---------------------------------------------------------------------------------#
         #    In case the derivative is bad, we give up on Newton's and go with Regula     #
         # Falsi.                                                                          #
         #---------------------------------------------------------------------------------#
         itn = itn + 1

         #----- Copy the previous guess. --------------------------------------------------#
         cia   = ciz
         funa  = fun
         #---------------------------------------------------------------------------------#


         #----- Update guess. -------------------------------------------------------------#
         ciz      = cia - fun/dfundci
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Check whether the method converged.                                         #
         #---------------------------------------------------------------------------------#
         converged = (2. * abs(cia-ciz)) < (tolerfl * (abs(cia)+abs(ciz)))
         #---------------------------------------------------------------------------------#

 
         #---------------------------------------------------------------------------------#
         #    At this point we test the current status of the guess.                       #
         #---------------------------------------------------------------------------------#
         if ((ciz < cimin) || (ciz > cimax)){
            #----- This guess went off-bounds, we give up on Newton's. --------------------#
            converged = FALSE
            iterate   = FALSE
         }else if (converged){
            #----- Converged, find the root as the mid-point. -----------------------------#
            ci      = 0.5 * (cia+ciz)
            iterate = FALSE
         }else{
            #----- Not there yet, update the function evaluation and the derivative. ------#
            iter    = iter.solver.step(met,thispft,aparms,ciz,newton=TRUE)
            fun     = iter$fun
            dfundci = iter$dfundci

            if (fun == 0.){
               #----- We have actually hit the jackpot, the answer is ciz. ----------------#
               ci        = ciz
               converged = TRUE
               iterate   = FALSE
            }# end if
         }#end if
         #---------------------------------------------------------------------------------#
      }#while ((itn < floor(maxfpofl/6)) && iterate && (abs(dfundci) < tolerfl))
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   if (! converged){

      #------------------------------------------------------------------------------------#
      #     If we have reached this point it means that Newton's method failed.  We        #
      # switch to the Regula Falsi instead.  The first step is to find two solutions of    #
      # the opposite side.                                                                 #
      #------------------------------------------------------------------------------------#
      if ((ciz < cimin) || (ciz > cimax)){
         #----- The guess is outside the range, discard it and start over. ----------------#
         cia      = sqrt(cimin*cimax)
         funa     = iter.solver.step(met,thispft,aparms,cia,newton=FALSE)
         zside    = FALSE
      }else if (funa * fun < 0.){
         funz     = fun
         zside    = TRUE
      }else{
         zside    = FALSE
      }#end if
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     If we still don't have a good second guess, look for one.                      #
      #------------------------------------------------------------------------------------#
      if (! zside){
         #----- Find the extrapolation term to try to hit the other side. -----------------#
         delta = 0.01 * min(abs(cimax-cia),abs(cimin-cia),abs(cimax-cimin))
         #---------------------------------------------------------------------------------#



         #----- First attempt. ------------------------------------------------------------#
         ciz   = cia + delta
         zside = FALSE
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Second guess seeker loop.                                                   #
         #---------------------------------------------------------------------------------#
         hitmin = FALSE
         hitmax = FALSE
         itb    = 0
         while (! zside){
            itb = itb + 1
            if (hitmin && hitmax){
               #---------------------------------------------------------------------------#
               #     We searched through the entire range of ci, and we couldn't find      #
               # any pair of roots of the opposite sign, it's likely that there is no      #
               # solution, so we give up.                                                  #
               #---------------------------------------------------------------------------#
               return(answer)
            }#end if

            ciz = cia + ((-1)**itb * (itb+3)/2) * delta
            if (ciz < cimin){
                #--------------------------------------------------------------------------#
                #    We have hit the minimum.  Force it to be the minimum, and make the    #
                # hitmin flag true.                                                        #
                #--------------------------------------------------------------------------#
                ciz    = cimin
                hitmin = TRUE
            }else if (ciz > cimax){
                ciz    = cimax
                hitmax = TRUE
            }#end if
            #------------------------------------------------------------------------------#
            #     Compute the function evaluate and check signs.                           #
            #------------------------------------------------------------------------------#
            funz  = iter.solver.step(met,thispft,aparms,ciz,newton=FALSE)
            zside = funa*funz < 0.
         }#end while
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     If we have reached this points, it means that there is a root that hasn't      #
      # been found yet, but at least we know that it is between cia and ciz.  Use the      #
      # modified Regula Falsi (Illinois) method to find it.                                #
      #------------------------------------------------------------------------------------#
      itb = itn - 1
      while (! converged && (itb < maxfpofl)){
         ci = (funz * cia - funa * ciz) / (funz-funa)

         #---------------------------------------------------------------------------------#
         #     Now that we updated the guess, check whether they are really close.         #
         # In case they are, this means that it converged, we can use this as our root.    #
         #---------------------------------------------------------------------------------#
         converged = 2. * abs(ci - cia) < tolerfl * (abs(cia)+abs(ciz))
         if (! converged){ 
            #----- Find the new function evaluation. --------------------------------------#
            fun  = iter.solver.step(met,thispft,aparms,ci,newton=FALSE)


            #------ Define the new interval based on the intermediate value theorem. ------#
            if (fun == 0.){
               #----- We have actually hit the jackpot, the answer is ciz. ----------------#
               converged = TRUE
            }else if (fun*funa < 0. ){
               ciz   = ci
               funz  = fun
               #----- If we are updating zside again, modify aside (Illinois method) ------#
               if (zside) funa=funa * 0.5
               #----- We just updated zside, set zside to true. ---------------------------#
               zside = TRUE
            }else{
               cia    = ci
               funa   = fun
               #----- If we are updating aside again, modify aside (Illinois method) ------#
               if (! zside) funz=funz * 0.5
               #----- We just updated aside, set zside to false. --------------------------#
               zside = FALSE
            }#end if
            #------------------------------------------------------------------------------#
         }#end if (! converged)
         #---------------------------------------------------------------------------------#
      }#end while (! converged && (itb < maxfpofl))
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #     In case it converged, we can compute the answer.                                  #
   #---------------------------------------------------------------------------------------#
   if (converged){
      #----- 1. Intercellular CO2, we just utilise the answer we've just got. -------------#
      answer$lint.co2      = ci
      #----- 3. Compute the CO2 demand. ---------------------------------------------------#
      answer$co2.demand    = calc.co2.demand(answer$lint.co2,aparms)
      #----- 4. Compute the actual stomatal conductance of water and CO2. -----------------#
      answer$stom.cond.h2o = calc.stom.cond.h2o(met,answer$lint.co2,answer$co2.demand)
      answer$stom.cond.co2 = gsw.2.gsc * answer$stom.cond.h2o
      #----- 5. Compute the leaf surface CO2. ---------------------------------------------#
      answer$lsfc.co2  = met$can.co2 - answer$co2.demand / met$blyr.cond.co2
      #----- 6. Compute the leaf surface vapour mixing ratio. -----------------------------#
      answer$lsfc.shv  = ( ( answer$stom.cond.h2o * met$lint.shv
                           + met$blyr.cond.h2o    * met$can.shv  )
                         / ( met$blyr.cond.h2o + answer$stom.cond.h2o) )
      #----- 7. Compute transpiration. ----------------------------------------------------#
      answer$transp    = met$blyr.cond.h2o * ( answer$lsfc.shv - met$can.shv )
      #----- 8. Export ci boundaries for model assessment. --------------------------------#
      answer$cimin     = cimin
      answer$cimax     = cimax
      answer$ci.assim  = ci.assim
      answer$ci.gsw    = ci.gsw
      answer$ci.shv    = ci.shv
      #----- 9. Update success flag. ------------------------------------------------------#
      answer$success   = TRUE
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     If we have hit this point, then we found an actual answer.                        #
   #---------------------------------------------------------------------------------------#
   return(answer)
   #---------------------------------------------------------------------------------------#
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function solves the model for the case where the carbon demand doesn't depend   #
# on the internal carbon.  This is simpler than the iterative case because we can solve    #
# through a quadratic equation for the stomatal conductance for water vapour.              #
#------------------------------------------------------------------------------------------#
solve.aofixed.case <<- function(met,thispft,aparms){
   #---------------------------------------------------------------------------------------#
   #   1. Initialise the success flag as false, and switch to true only when we obtain a   #
   #      valid solution.                                                                  #
   #---------------------------------------------------------------------------------------#
   answer               = aparms
   answer$lsfc.shv      = NA_real_
   answer$lsfc.co2      = NA_real_
   answer$lint.co2      = NA_real_
   answer$co2.demand    = NA_real_
   answer$stom.cond.h2o = NA_real_
   answer$stom.cond.co2 = NA_real_
   answer$transp        = NA_real_
   answer$cimin         = NA_real_
   answer$cimax         = NA_real_
   answer$ci.assim      = NA_real_
   answer$ci.gsw        = NA_real_
   answer$ci.shv        = NA_real_
   answer$success       = FALSE
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Determine the minimum and maximum intercellular CO2 that will still produce  a     #
   # positive and bounded conductance, which should be above the cuticular conductance     #
   # (otherwise there is no reason to keep stomata opened.                                 #
   #---------------------------------------------------------------------------------------#
   co2.bounds = find.lint.co2.bounds(met,thispft,aparms)
   bounded    = co2.bounds$bounded
   cimin      = co2.bounds$cimin
   cimax      = co2.bounds$cimax
   #----- Export ci boundaries for model assessment. --------------------------------------#
   answer$cimin     = co2.bounds$cimin
   answer$cimax     = co2.bounds$cimax
   answer$ci.assim  = co2.bounds$ci.assim
   answer$ci.gsw    = co2.bounds$ci.gsw
   answer$ci.shv    = co2.bounds$ci.shv
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Since the carbon demand doesn't depend on the intercellular CO2, compute it      #
   # using the first guess.                                                                #
   #---------------------------------------------------------------------------------------#
   answer$co2.demand = calc.co2.demand(met$can.co2,aparms)
   #---------------------------------------------------------------------------------------#



   #----- Compute the leaf surface CO2. ---------------------------------------------------#
   answer$lsfc.co2    = met$can.co2 - answer$co2.demand / met$blyr.cond.co2
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #       In case the leaf surface CO2 is unrealistic (i.e. negative), don't bother       #
   # looking for a solution.                                                               #
   #---------------------------------------------------------------------------------------#
   if (answer$lsfc.co2 < c34smin.lint.co2) return(answer)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #   Check the sign of the carbon demand.                                                #
   #---------------------------------------------------------------------------------------#
   if (answer$co2.demand <= 0.0){
      #------------------------------------------------------------------------------------#
      #     If carbon demand is zero or negative, this means that light is below the       #
      # light compensation point, so all stomata should remain closed.                     #
      #------------------------------------------------------------------------------------#
      answer$stom.cond.h2o = thispft$b
      answer$stom.cond.co2 = gsw.2.gsc * answer$stom.cond.h2o
      answer$lint.co2      = answer$lsfc.co2 - answer$co2.demand / answer$stom.cond.co2
   }else{
      #------------------------------------------------------------------------------------#
      #     Carbon demand is positive, look for a solution.                                #
      #------------------------------------------------------------------------------------#
      #----- Find auxiliary coefficients to compute the quadratic terms. ------------------#
      qterm1 = (met$can.co2 - aparms$compp) * met$blyr.cond.co2 - answer$co2.demand
      qterm2 = (thispft$d0 + met$lint.shv - met$can.shv)  * met$blyr.cond.h2o
      qterm3 = thispft$m * answer$co2.demand * thispft$d0 * met$blyr.cond.co2
      #----- Find the coefficients for the quadratic equation. ----------------------------#
      aquad = qterm1 * thispft$d0
      bquad = qterm1 * qterm2 - aquad * thispft$b - qterm3
      cquad = - qterm1 * qterm2 * thispft$b - qterm3 * met$blyr.cond.h2o
      #----- Solve the quadratic equation for gsw. ----------------------------------------#
      if (aquad == 0.0){
         #----- Not really a quadratic equation. ------------------------------------------#
         gswroot1 = -cquad / bquad
         gswroot2 = discard
      }else{
         #----- A quadratic equation, find the discriminant. ------------------------------#
         discr = bquad * bquad - 4.0 * aquad * cquad
         #----- Decide what to do based on the discriminant. ------------------------------#
         if (discr == 0.0){
            #----- Double root. -----------------------------------------------------------#
            gswroot1 = - bquad / (2.0 * aquad)
            gswroot2 = discard
         }else if (discr > 0.0){
            #----- Two distinct roots. ----------------------------------------------------#
            gswroot1 = (- bquad - sqrt(discr)) / (2.0 * aquad)
            gswroot2 = (- bquad + sqrt(discr)) / (2.0 * aquad)
         }else{
            #----- None of the roots are real, this solution failed. ----------------------#
            return(answer)
         }#end if
      }#end if
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Check both solutions, and decide which one makes sense.  Once the right        #
      # one is determined, compute the stomatal resistance for CO2 and the inter-          #
      # cellular CO2 concentration.  In case both make solutions make sense (unlikely),    #
      # we decide the root based on the intercellular CO2.                                 #
      #------------------------------------------------------------------------------------#
      bounded1 = gswroot1 >= thispft$b && gswroot1 <= c34smax.gsw
      bounded2 = gswroot2 >= thispft$b && gswroot2 <= c34smax.gsw
      if (bounded1 && bounded2){
         #----- Both solutions are valid, warn the user as this should never happen. ------#
         ciroot1 = answer$lsfc.co2 - answer$co2.demand / (gsw.2.gsc * gswroot1)
         ciroot2 = answer$lsfc.co2 - answer$co2.demand / (gsw.2.gsc * gswroot2)

         bounded1 = (  ciroot1 >= c34smin.lint.co2
                    && ciroot1 <= min(c34smax.lint.co28,met$can.co2)
                    )#end bounded1
         bounded2 = (  ciroot2 >= c34smin.lint.co2
                    && ciroot2 <= min(c34smax.lint.co28,met$can.co2)
                    )#end bounded2

         if (bounded1 && bounded2){
            #----- Both intercellular CO2 work, pick the highest and warn the user. -------#
            if (ciroot1 >= ciroot2){
               answer$stom.cond.h2o = gswroot1
               answer$stom.cond.co2 = gsw.2.gsc * answer$stom.cond.h2o
               answer$lint.co2      = ciroot1
            }else{
               answer$stom.cond.h2o = gswroot2
               answer$stom.cond.co2 = gsw.2.gsc * answer$stom.cond.h2o
               answer$lint.co2      = ciroot2
            }#end if
         }else if (bounded1){
            answer$stom.cond.h2o = gswroot1
            answer$stom.cond.co2 = gsw.2.gsc * answer$stom.cond.h2o
            answer$lint.co2      = ciroot1
         }else if (bounded2){
            answer$stom.cond.h2o = gswroot2
            answer$stom.cond.co2 = gsw.2.gsc * answer$stom.cond.h2o
            answer$lint.co2      = ciroot2
         }else{
            return(answer)
         }#end if
      }else if (bounded1){
         #----- First root is the only one that makes sense. ------------------------------#
         answer$stom.cond.h2o = gswroot1
         answer$stom.cond.co2 = gsw.2.gsc * answer$stom.cond.h2o
         answer$lint.co2      = answer$lsfc.co2 - answer$co2.demand / answer$stom.cond.co2

      }else if (bounded2){
         #----- Second root is the only one that makes sense. -----------------------------#
         answer$stom.cond.h2o = gswroot2
         answer$stom.cond.co2 = gsw.2.gsc * answer$stom.cond.h2o
         answer$lint.co2      = answer$lsfc.co2 - answer$co2.demand / answer$stom.cond.co2
      }else{
         #----- None of the solutions are bounded.  This solution failed. -----------------#
         return(answer)
      }#end if
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #   Find the surface water specific humidity and the transpiration.                     #
   #---------------------------------------------------------------------------------------#
   answer$lsfc.shv = ( ( answer$stom.cond.h2o * met$lint.shv
                       + met$blyr.cond.h2o    * met$can.shv  )
                     / ( met$blyr.cond.h2o + answer$stom.cond.h2o) )
   answer$transp   = answer$stom.cond.h2o * ( answer$lsfc.shv - met$can.shv )
   #---------------------------------------------------------------------------------------#

   #----- In case we reach this point, we successfully found the solution. ----------------#
   answer$success = TRUE
   #---------------------------------------------------------------------------------------#


   #----- Return the data frame. ----------------------------------------------------------#
   return(answer)
   #---------------------------------------------------------------------------------------#
}#end function solve.aofixed.case
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This subroutine finds the updated function evaluation and derivative for the         #
# iterative step for Newton's (or Regula Falsi) method.  The "function" here is a          #
# combination of the definition of the water stomatal conductance for open stomata         #
# (F96's equation 13), after substituting Ds by a combination of M09's equations B13       #
# and B16 so water demand is eliminated(Psi_open), and incorporating F96's equation 14     #
# to eliminate the surface carbon.  This function has the property of having one root      #
# that corresponds to the intercellular CO2 concentration.                                 #
#  The logical variable is used to decide whether to compute the derivatives or not.       #
# They are necessary only when it is a Newton's method call.                               #
#------------------------------------------------------------------------------------------#
iter.solver.step <<- function(met,thispft,aparms,lint.co2,newton=TRUE){

   #----- Find the CO2 demand. ------------------------------------------------------------#
   co2.demand       = calc.co2.demand(lint.co2,aparms)
   #----- Find the stomatal conductance of water. -----------------------------------------#
   stom.cond.h2o    = calc.stom.cond.h2o(met,lint.co2,co2.demand)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the function components, then the function evaluation.                       #
   #---------------------------------------------------------------------------------------#
   efun1 = (stom.cond.h2o - thispft$b) / (thispft$m * co2.demand)
   efun2 = (met$can.co2 - aparms$compp - co2.demand/met$blyr.cond.co2)
   efun3 = 1.0 + ( met$blyr.cond.h2o * (met$lint.shv - met$can.shv)
                 / (thispft$d0 * (met$blyr.cond.h2o + stom.cond.h2o) ) )
   fun   = efun1 * efun2 * efun3 - 1.0
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     We also compute the derivative.                                                   #
   #---------------------------------------------------------------------------------------#
   if (newton){
      #----- CO2 demand. ------------------------------------------------------------------#
      co2.demand.prime    = calc.co2.demand.prime(lint.co2,co2.demand,aparms)
      #----- stomatal conductance of water. -----------------------------------------------#
      stom.cond.h2o.prime = calc.stom.cond.h2o.prime(met,lint.co2,stom.cond.h2o
                                                    ,co2.demand,co2.demand.prime)
      #----- Function components. ---------------------------------------------------------#
      eprime1 = ( ( stom.cond.h2o.prime * co2.demand
                  - co2.demand.prime * (stom.cond.h2o - thispft$b) )
                / (thispft$m * co2.demand * co2.demand) )

      eprime2 = - co2.demand.prime / met$blyr.cond.co2

      eprime3 = ( - (efun3 - 1.0) *  stom.cond.h2o.prime
                / ( met$blyr.cond.h2o + stom.cond.h2o ) )

      dfundci = ( eprime1 * efun2   * efun3
                + efun1   * eprime2 * efun3
                + efun1   * efun2   * eprime3 )
   }else{
      dfundci = 0.0
   }#end if
   #------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #      Change the result depending on whether the user wants the derivative or not.     #
   #---------------------------------------------------------------------------------------#
   if (newton){
      ans=list(fun=fun,dfundci=dfundci)
   }else{
      ans=fun
   }#end if
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function iter.solver.step 
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This subroutine finds the updated function evaluation and derivative for the         #
# fixed case.  This is intended for diagnostics only.                                      #
#------------------------------------------------------------------------------------------#
fixed.solver.step <<- function(met,thispft,aparms,lint.co2){

   #----- Handy aliases. ------------------------------------------------------------------#
   ao      = aparms$co2.demand
   gsc     = aparms$stom.cond.co2
   gbc     = met$blyr.cond.co2
   can.co2 = met$can.co2
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the function components, then the function evaluation.                       #
   #---------------------------------------------------------------------------------------#
   fun   = gbc * gsc * (can.co2 - lint.co2) / (gbc + gsc) / ao - 1.
   return(fun)
   #---------------------------------------------------------------------------------------#
}#end function fixed.solver.step 
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      This function computes the CO2 demand given the intercellular CO2 concentration.    #
#------------------------------------------------------------------------------------------#
calc.co2.demand <<- function(lint.co2,aparms){

   ans = ( ( aparms$rho * lint.co2 + aparms$sigma)
         / ( aparms$xi  * lint.co2 + aparms$tau  ) ) + aparms$nu

   return(ans)
}#end calc.co2.demand
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function computes the derivative of the CO2 regarding the intercellular CO2     #
# concentration.                                                                           #
#------------------------------------------------------------------------------------------#
calc.co2.demand.prime <<- function(lint.co2,co2.demand,aparms){

   ans = ( ( aparms$rho - aparms$xi * (co2.demand - aparms$nu) )
         / ( aparms$xi  * lint.co2 + aparms$tau  ) )
   return(ans)
}#end function calc.co2.demand.prime
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function computes the stomatal conductance of water given the intercellular     #
# CO2 concentration and the CO2 demand.                                                    #
#------------------------------------------------------------------------------------------#
calc.stom.cond.h2o <<- function(met,lint.co2,co2.demand){

   ans = ( met$blyr.cond.co2 * co2.demand
         / (gsw.2.gsc * ( (met$can.co2 - lint.co2) * met$blyr.cond.co2  - co2.demand ) ) )

   return(ans)
}#end function calc_stom_cond_h2o
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function computes the stomatal conductance of water derivative given the inter- #
# cellular CO2 concentration, the water stomatal conductance, and the CO2 demand and its   #
# derivative.                                                                              #
#------------------------------------------------------------------------------------------#
calc.stom.cond.h2o.prime <<- function( met
                                     , lint.co2
                                     , stom.cond.h2o
                                     , co2.demand
                                     , co2.demand.prime
                                     ){

   ans = stom.cond.h2o * ( co2.demand.prime / co2.demand
                         + (met$blyr.cond.co2 + co2.demand.prime)
                         / (gsw.2.gsc * ( (met$can.co2 - lint.co2)
                                        * met$blyr.cond.co2 - co2.demand) ) )

   return(ans)
}#end end function calc.stom.cond.h2o.prime
#==========================================================================================#
#==========================================================================================#
