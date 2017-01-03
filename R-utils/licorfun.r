#==========================================================================================#
#==========================================================================================#
#      Some parameters for the solver.                                                     #
#------------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #     This is a flag used in various sub-routines and functions and denote that we      #
   # should ignore the result.                                                             #
   #---------------------------------------------------------------------------------------#
   discard          <<- Inf
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    This is the tolerance for iterative methods.  Some of the functions are quite flat #
   # so it is a good idea to use a somewhat more strict tolerance than the ones used in    #
   # therm_lib8.                                                                           #
   #---------------------------------------------------------------------------------------#
   tolerfl          <<- 1.e-8
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Since this is more strict than therm_lib, the maximum number of attempts for the   #
   # Regula Falsi method should be also increased.                                         #
   #---------------------------------------------------------------------------------------#
   maxfpofl         <<- 320
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Minimum and maximum acceptable intercellular CO2.  Don't be strict, but don't let  #
   # the boundaries be physically impossible.                                              #
   #---------------------------------------------------------------------------------------#
   c34smin.lint.co2 <<- 0.5   * umol.2.mol
   c34smax.lint.co2 <<- 1200. * umol.2.mol
   #---------------------------------------------------------------------------------------#


#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      This function finds the assimilation rate parameters for a given plant functional   #
# type, meteorological conditions, and limitation.                                         #
#------------------------------------------------------------------------------------------#
assim.params <<- function(ipft,env,limit){
   aparms=list()
   if (limit %in% "CLOSED"){
      #----- Closed stomata case, or night time.  These are the same for C3 and C4. -------#
      aparms$rho   = 0.0
      aparms$sigma = 0.0
      aparms$xi    = 0.0
      aparms$tau   = 1.0
      aparms$nu    = - env$leaf.resp
      #---------------------------------------------------------------------------------------#
      #     Open stomata case, so now we distinguish between C3 and C4 as their functional    #
      # forms are different.                                                                  #
      #---------------------------------------------------------------------------------------#
   }else if (pft$pathway[ipft] %in% 3){
      #------------------------------------------------------------------------------------#
      #     C3 case.  Decide whether this is the light- or RuBP-saturated case.            #
      #------------------------------------------------------------------------------------#
      if (limit %in% "LIGHT"){
         #---- Light-limited case. --------------------------------------------------------#
         aparms$rho   =  env$alpha * env$par
         aparms$sigma = -env$alpha * env$par * env$compp
         aparms$xi    = 1.0
         aparms$tau   = 2.0 * env$compp
         aparms$nu    = - env$leaf.resp

      }else if (limit %in% "RUBP" || limit %in% "CO2"){
         #----- RuBP-saturated rate of photosynthesis case. -------------------------------#
         aparms$rho   =  env$vm
         aparms$sigma = -env$vm * env$compp
         aparms$xi    = 1.0
         aparms$tau   = env$kco2 * (1.0 + env$o2 / env$ko2)
         aparms$nu    = - env$leaf.resp
      }#end if limit
      #------------------------------------------------------------------------------------#
   }else if (pft$pathway[ipft] %in% 4){
      #------------------------------------------------------------------------------------#
      #     C4 case.  There are three possibilities, the light-limited, the RuBP-          #
      # saturated, and the CO2-limited cases.                                              #
      #------------------------------------------------------------------------------------#
      if (limit %in% "LIGHT"){
         #----- Light-limited case. -------------------------------------------------------#
         aparms$rho   = 0.0
         aparms$sigma = env$alpha * env$par
         aparms$xi    = 0.0
         aparms$tau   = 1.0
         aparms$nu    = - env$leaf.resp

      }else if (limit %in% "RUBP"){
         #----- RuBP-saturated rate of photosynthesis case. -------------------------------#
         aparms$rho   = 0.0
         aparms$sigma = env$vm
         aparms$xi    = 0.0
         aparms$tau   = 1.0
         aparms$nu    = - env$leaf.resp

      }else if (limit %in% "CO2"){
         #----- CO2-limited for low CO2 concentration case. -------------------------------#
         aparms$rho   = klowco2 * env$vm
         aparms$sigma = 0.0
         aparms$xi    = 0.0
         aparms$tau   = 1.0
         aparms$nu    = - env$leaf.resp

      }#end if limit
      #------------------------------------------------------------------------------------#
   }#end if
   aparms$iterative = aparms$rho != 0 || aparms$xi != 0
   #---------------------------------------------------------------------------------------#


   #----- Find the minimum radiation that makes Ao = 0. (light compensation point) --------#
   aparms$par.min  = ( ( env$leaf.resp * (env$co2 + 2. * env$compp) )
                     / ( env$alpha     * (env$co2 -      env$compp) ) ) 

   aparms$success  = (limit %in% "CLOSED") || (env$par > aparms$par.min)
   aparms$cimin    = c34smin.lint.co2
   aparms$cimax    = c34smax.lint.co2
   aparms$ci.gsw   = c34smax.lint.co2
   aparms$ci.shv   = c34smax.lint.co2
   aparms$ci.assim = c34smin.lint.co2
   #---------------------------------------------------------------------------------------#

   return(aparms)
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      This function finds the assimilation rate parameters given the parameters and the   #
# intercellular CO2.                                                                       #
#------------------------------------------------------------------------------------------#
assim.rate <<- function(ci,aparms){
   #------ Find the assimilation rate. ----------------------------------------------------#
   ao = (aparms$rho * ci + aparms$sigma) / (aparms$xi * ci + aparms$tau) + aparms$nu
   return(ao)
   #---------------------------------------------------------------------------------------#
}#end function
#------------------------------------------------------------------------------------------#






#==========================================================================================#
#==========================================================================================#
#     This function computes the derivative of the CO2 regarding the intercellular CO2     #
# concentration.                                                                           #
#------------------------------------------------------------------------------------------#
assim.rate.prime <<- function(ci,assim,aparms){
   
   ao.prime = ( (aparms$rho - aparms$xi * (assim - aparms$nu) ) 
              / (aparms$xi  * ci + aparms$tau  ) )
   return(ao.prime)
}#end function assim.rate.prime
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      This function finds the assimilation rate parameters given the parameters and the   #
# intercellular CO2.                                                                       #
#------------------------------------------------------------------------------------------#
stom.cond <<- function(ci,assim,env){
   #------ Find the assimilation rate. ----------------------------------------------------#
   gsw = ( env$gbc * assim ) / ( gsw.2.gsc * ( (env$co2 - ci) * env$gbc - assim ) )
   return(gsw)
   #---------------------------------------------------------------------------------------#
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function computes the stomatal conductance of water derivative given the        #
# intercellular CO2 concentration, the water stomatal conductance, and the CO2 demand      #
# and its derivative.                                                                      #
#------------------------------------------------------------------------------------------#
stom.cond.prime <<- function(ci,gsw,assim,assim.prime,env){

   gsw.prime = ( gsw * ( assim.prime / assim
                       + (env$gbc + assim.prime)
                       / (gsw.2.gsc * ( (env$co2 - ci) * env$gbc - assim) ) ) )

   return(gsw.prime)
}#end function stom.cond.prime
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     Function whose root gives Ci.                                                        #
#------------------------------------------------------------------------------------------#
cibounds <<- function(env,ipft,aparms){
   #---------------------------------------------------------------------------------------#
   # First case: This check will find when Aopen goes to 0., which causes a singularity    #
   # in the function of which we are looking for a root.                                   #
   #---------------------------------------------------------------------------------------#
   ciAo = - (aparms$tau * aparms$nu + aparms$sigma) / (aparms$xi  * aparms$nu + aparms$rho)
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   # Second case: This check finds which ci causes the terms [gbc (ca -ci) - Aopen] to     #
   # be 0.   This will cause a singularity in the gsw function.                            #
   #---------------------------------------------------------------------------------------#
   #----- 1. Define the coefficients for the quadratic equation. --------------------------#
   aquad = env$gbc * aparms$xi
   bquad = ( aparms$xi * (aparms$nu - env$gbc * env$co2 )
           + env$gbc * aparms$tau + aparms$rho )
   cquad = aparms$tau * (aparms$nu - env$gbc * env$co2 ) + aparms$sigma
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
   # Third case: This is the upper threshold when the CO2 assimilation is negligible       #
   # because the intercellular CO2 is above the ambient.                                   #
   #---------------------------------------------------------------------------------------#
   epsil = gsc.2.gsw - gbc.2.gbw
   #----- Define the coefficients for the quadratic equation. -----------------------------#
   aquad = env$gbw * aparms$xi
   bquad = - ( epsil * ( aparms$rho + aparms$nu * aparms$xi)
             + env$gbw * aparms$xi * env$co2 )
   cquad = - epsil * ( aparms$sigma + aparms$tau * aparms$nu )
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
   ciroot2 = min(ciQ,cigsw,env$co2)
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
   if (! is.finite(cimin) || ! is.finite(cimax)){
      cimin   = c34smin.lint.co2
      cimax   = c34smin.lint.co2
      bounded = FALSE
   }else if (cimin > cimax){
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

   ans = list( cimin    = cimin
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
#     Function that finds Ci.                                                              #
#------------------------------------------------------------------------------------------#
cisolver <<- function(env,ipft,aparms){


   #----- Set up the anwer with the default (failure). ------------------------------------#
   ans = list( ci       =  discard
             , cs       =  discard
             , assim    =  discard
             , transp   =  discard
             , gsw      =  discard
             , gsc      =  discard
             , ws       =  discard
             , cimin    =  aparms$cimin
             , cimax    =  aparms$cimax
             , ci.shv   =  aparms$ci.shv
             , ci.assim =  aparms$ci.assim
             , ci.gsw   =  aparms$ci.gsw
             , success  = aparms$success
             )#end list
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #       Check whether there is a chance for solution...                                 #
   #---------------------------------------------------------------------------------------#
   if (! ans$success){ 
      return(ans)

   }else if (aparms$iterative){
      #----- Find the boundaries for ci. --------------------------------------------------#
      if (ans$success){
         cilim        = cibounds(env,ipft,aparms)
         ans$cimin    = cilim$cimin
         ans$cimax    = cilim$cimax
         ans$ci.shv   = cilim$ci.shv
         ans$ci.assim = cilim$ci.assim
         ans$ci.gsw   = cilim$ci.gsw
         #----- Quit if the function failed to find reasonable bounds. --------------------#
         if (! cilim$bounded){
            ans$success  = FALSE
            return(ans)
         }#end if
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #      Initialise the convergence flag.  Here we start realistic, ops, I mean,       #
      # pessimistic, and assume that we are failing.  This will be switched to true only   #
      # if a real answer is found, so in case we quit the sub-routine due to impossible    #
      # solution, the result will be failure.                                              #
      #------------------------------------------------------------------------------------#
      #------------------------------------------------------------------------------------#
      cia   = sqrt(cilim$cimin*cilim$cimax)
      ita   = cifun.eval(cia,env,ipft,aparms)
      #------------------------------------------------------------------------------------#



      #----- Copy to the new guess just in case it fails at the first iteration -----------#
      ciz = cia
      it  = ita
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Check whether we have hit the jackpot, if so, quit because we have found the  #
      # answer.                                                                            #
      #------------------------------------------------------------------------------------#
      if (it$fun == 0.0){ 
         ci        = ciz
         converged = TRUE
      }else{
         converged = FALSE
      }#end if
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Enter Newton's method loop in case we haven't found the answer.                #
      #------------------------------------------------------------------------------------#
      if (! converged){
         again = TRUE
         itn   = 0
         while (again){
            itn = itn + 1

            #------------------------------------------------------------------------------#
            #    In case the derivative is bad, we give up on Newton's and go with Regula  #
            # Falsi.                                                                       #
            #------------------------------------------------------------------------------#
            again = itn <= floor(maxfpofl/6) && abs(it$deriv) >= tolerfl
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Newton's step.                                                           #
            #------------------------------------------------------------------------------#
            if (again){

               #----- Copy the previous guess. --------------------------------------------#
               cia   = ciz
               ita   = it
               #---------------------------------------------------------------------------#


               #----- Update guess. -------------------------------------------------------#
               ciz      = cia - it$fun/it$deriv
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #     Check whether the method converged.                                   #
               #---------------------------------------------------------------------------#
               converged = 2.0 * abs(cia-ciz) < tolerfl * (abs(cia)+abs(ciz))
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #    At this point we test the current status of the guess.                 #
               #---------------------------------------------------------------------------#
               if (ciz < cilim$cimin || ciz > cilim$cimax){
                  #----- This guess went off-bounds, we give up on Newton's. --------------#
                  converged = FALSE
                  again     = FALSE
               }else if (converged){
                  #----- Converged, find the root as the mid-point. -----------------------#
                  ci    = 0.5 * (cia+ciz)
                  again = FALSE
               }else{
                  #----- Not there yet, update the function evaluation and the derivative. #
                  it = cifun.eval(ciz,env,ipft,aparms)
                  if (it$fun == 0.0){ 
                     #----- We have actually hit the jackpot, the answer is ciz. ----------#
                     ci        = ciz
                     converged = TRUE
                     again     = FALSE
                  }#end if (it$fun == 0.0)
                  #------------------------------------------------------------------------#
               }#end if (ciz < cilim$cimin || ciz > cilim$cimax)
               #---------------------------------------------------------------------------#
            }#end if (again)
            #------------------------------------------------------------------------------#
         }#end while (again)
         #---------------------------------------------------------------------------------#
      }#end if (! converged)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #    Check whether we have an answer or not.                                         #
      #------------------------------------------------------------------------------------#
      if (! converged){
         #---------------------------------------------------------------------------------#
         #     If we have reached this point it means that Newton's method failed.  We     #
         # switch to the Regula Falsi instead.  The first step is to find two solutions of #
         # the opposite side.                                                              #
         #---------------------------------------------------------------------------------#
         if (ciz < cilim$cimin || ciz > cilim$cimax){
            #----- The guess is outside the range, discard it and start over. -------------#
            cia   = sqrt(cilim$cimin*cilim$cimax)
            ita   = cifun.eval(cia,env,ipft,aparms)
            zside = FALSE
         }else if (ita$fun * it$fun < 0.0){
            itz   = it
            zside = TRUE
         }else{
            zside = FALSE
         }#end if
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #     If we still don't have a good second guess, look for one.                   #
         #---------------------------------------------------------------------------------#
         if (! zside){
            #----- Find the extrapolation term to try to hit the other side. --------------#
            delta = 0.01 * min( abs(cilim$cimax-cia        )
                              , abs(cilim$cimin-cia        )
                              , abs(cilim$cimax-cilim$cimin)
                              )#end min
            #------------------------------------------------------------------------------#



            #----- First attempt. ---------------------------------------------------------#
            ciz   = cia + delta
            zside = FALSE
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Second guess seeker loop.                                                #
            #------------------------------------------------------------------------------#
            hitmin = FALSE
            hitmax = FALSE
            itb    = 0
            again  = TRUE
            while (! zside){
               itb = itb + 1

               #---------------------------------------------------------------------------#
               #      Make sure that the function is behaving well...                      #
               #---------------------------------------------------------------------------#
               if (hitmin && hitmax){
                  #------------------------------------------------------------------------#
                  #     We searched through the entire range of ci, and we couldn't find   #
                  # any pair of roots of the opposite sign, it's likely that there is no   #
                  # solution, so we give up.                                               #
                  #------------------------------------------------------------------------#
                  ans$success = FALSE
                  return(ans)
               }#end if (hitmin && hitmax)
               #---------------------------------------------------------------------------#


               #----- Update the guess for ciz. -------------------------------------------#
               ciz = cia + ((-1)^itb * (itb+3)/2) * delta
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #     Flag any attempt to cross the reasonable limits.                      #
               #---------------------------------------------------------------------------#
               if (ciz < cilim$cimin){
                   #-----------------------------------------------------------------------#
                   #    We have hit the minimum.  Force it to be the minimum, and make the #
                   # hitmin flag true.                                                     #
                   #-----------------------------------------------------------------------#
                   ciz    = cilim$cimin
                   hitmin = TRUE
               }else if (ciz > cilim$cimax){
                   #-----------------------------------------------------------------------#
                   #    We have hit the maximum.  Force it to be the maximum, and make the #
                   # hitmax flag true.                                                     #
                   #-----------------------------------------------------------------------#
                   ciz    = cilim$cimax
                   hitmax = TRUE
               }#end if
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #     Compute the function evaluate and check signs.                        #
               #---------------------------------------------------------------------------#
               itz   = cifun.eval(ciz,env,ipft,aparms)
               zside = ita$fun * itz$fun < 0.0
               #---------------------------------------------------------------------------#
            }# end while (! zside)
            #------------------------------------------------------------------------------#
         }#end if (! zside)
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #     If we have reached this point, it means that there is a root that hasn't    #
         # been found yet, but at least we know that it is between cia and ciz.  Use the   #
         # modified Regula Falsi (Illinois) method to find it.                             #
         #---------------------------------------------------------------------------------#
         again = TRUE
         itb   = itn - 1
         while (again){
            #------ Update the step count, and update the flag for continuing. ------------#
            itb   = itb + 1
            again = again && itb <= maxfpofl
            #------------------------------------------------------------------------------#

            if (again){

               #----- Update the guess. ---------------------------------------------------#
               ci  = ( itz$fun * cia - ita$fun * ciz ) / ( itz$fun - ita$fun )
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #     Now that we updated the guess, check whether they are really close.   #
               # In case they are, this means that it converged, we can use this as our    #
               # root.                                                                     #
               #---------------------------------------------------------------------------#
               converged = 2.0 * abs(ci - cia) < tolerfl * (abs(cia)+abs(ciz))
               if (converged) again = FALSE
               #---------------------------------------------------------------------------#


               #----- Find the new function evaluation. -----------------------------------#
               it = cifun.eval(ci,env,ipft,aparms)
               #---------------------------------------------------------------------------#


               #------ Define the new interval based on the intermediate value theorem. ---#
               if (it$fun == 0.0){
                  #----- We have actually hit the jackpot, the answer is ciz. -------------#
                  converged = TRUE
                  again     = FALSE
               }else if (it$fun*ita$fun < 0.0 ){
                  ciz   = ci
                  itz   = it
                  #----- If we are updating zside again, modify aside (Illinois method) ---#
                  if (zside) ita$fun=ita$fun * 0.5
                  #----- We just updated zside, set zside to true. ------------------------#
                  zside = TRUE
               }else{
                  cia    = ci
                  ita    = it
                  #----- If we are updating aside again, modify aside (Illinois method) ---#
                  if (! zside) itz$fun=itz$fun * 0.5
                  #----- We just updated aside, set zside to false. -----------------------#
                  zside = FALSE
               }#end if (it$fun == 0.0)
               #---------------------------------------------------------------------------#
            }#end if (again)
            #------------------------------------------------------------------------------#
         }#end while (again)
         #---------------------------------------------------------------------------------#
      }#end if (! converged)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     We update the results only if the solver has found an answer.                  #
      #------------------------------------------------------------------------------------#
      if (converged){
         assim       = assim.rate(ci,aparms)
         gsw         = stom.cond(ci,assim,env)
         gsc         = gsw.2.gsc * gsw
         transp      = env$gbw * gsw * ( env$wi - env$h2o) / (env$gbw + gsw)
         cs          = env$co2 - assim  / env$gbc 
         ws          = env$h2o + transp / env$gbw
         ans$success = TRUE
      }else{
         ans$success = FALSE
         return(ans)
      }#end if
      #------------------------------------------------------------------------------------#

   }else{
      #------------------------------------------------------------------------------------#
      #     This sub-block solves the model for the case where the carbon demand doesn't   #
      # depend on the internal carbon.  This is simpler than the iterative case because we #
      # can solve through a quadratic equation for the stomatal conductance for water      #
      # vapour.                                                                            #
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Since the carbon demand doesn't depend on the intercellular CO2, compute it   #
      # using the first guess.                                                             #
      #------------------------------------------------------------------------------------#
      assim = assim.rate(env$co2,aparms)
      #------------------------------------------------------------------------------------#



      #----- Compute the leaf surface CO2. ------------------------------------------------#
      cs    = env$co2 - assim / env$gbc
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #   Check the sign of the carbon demand.                                             #
      #------------------------------------------------------------------------------------#
      if (assim <= 0.0){
         #---------------------------------------------------------------------------------#
         #     If carbon demand is zero or negative, this means that light is below the    #
         # light compensation point, so all stomata should remain closed.                  #
         #---------------------------------------------------------------------------------#
         gsw = pft$b[ipft]
         gsc = gsw.2.gsc * gsw
         ci  = cs - assim / gsc
      }else{
         #---------------------------------------------------------------------------------#
         #     Carbon demand is positive, look for a solution.                             #
         #---------------------------------------------------------------------------------#
         #----- Find auxiliary coefficients to compute the quadratic terms. ---------------#
         qterm1 = (env$co2 - env$compp) * env$gbc - assim
         qterm2 = (pft$d0[ipft] + env$wi - env$h2o) * env$gbw
         qterm3 = pft$m[ipft] * assim * pft$d0[ipft] * env$gbc
         #----- Find the coefficients for the quadratic equation. -------------------------#
         aquad = qterm1 * pft$d0[ipft]
         bquad = qterm1 * qterm2 - aquad * pft$b[ipft] - qterm3
         cquad = - qterm1 * qterm2 * pft$b[ipft] - qterm3 * env$gbw
         #----- Solve the quadratic equation for gsw. -------------------------------------#
         if (aquad == 0.0){
            #----- Not really a quadratic equation. ---------------------------------------#
            gswroot1 = -cquad / bquad
            gswroot2 = discard
         }else{
            #----- A quadratic equation, find the discriminant. ---------------------------#
            discr = bquad * bquad - 4.0 * aquad * cquad
            #----- Decide what to do based on the discriminant. ---------------------------#
            if (discr == 0.0){
               #----- Double root. --------------------------------------------------------#
               gswroot1 = - bquad / (2.0 * aquad)
               gswroot2 = discard
            }else if (discr > 0.0){
               #----- Two distinct roots. -------------------------------------------------#
               gswroot1 = (- bquad - sqrt(discr)) / (2.0 * aquad)
               gswroot2 = (- bquad + sqrt(discr)) / (2.0 * aquad)
            }else{
               #----- None of the roots are real, this solution failed. -------------------#
               ans$success = FALSE
               return(ans)
            }#end if
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Check both solutions, and decide which one makes sense.  Once the right     #
         # one is determined, compute the stomatal resistance for CO2 and the inter-       #
         # cellular CO2 concentration.  In case both make solutions make sense (unlikely), #
         # we decide the root based on the intercellular CO2.                              #
         #---------------------------------------------------------------------------------#
         bounded1 = gswroot1 >= pft$b[ipft] && gswroot1 <= c34smax.gsw
         bounded2 = gswroot2 >= pft$b[ipft] && gswroot2 <= c34smax.gsw
         if (bounded1 && bounded2){
            #----- Both solutions are valid, warn the user as this should never happen. ---#
            ciroot1 = cs - assim / (gsw.2.gsc * gswroot1)
            ciroot2 = cs - assim / (gsw.2.gsc * gswroot2)

            bounded1 = (  ciroot1 >= c34smin.lint.co2
                       && ciroot1 <= min(c34smax.lint.co2,env$co2) )
            bounded2 = (  ciroot2 >= c34smin.lint.co28
                       && ciroot2 <= min(c34smax.lint.co2,env$co2) )
             
            if (bounded1 && bounded2){
               #----- Both intercellular CO2 work, pick the highest and warn the user. ----#
               if (ciroot1 >= ciroot2){
                  gsw = gswroot1
                  gsc = gsw.2.gsc * gsw
                  ci  = ciroot1
               }else{
                  gsw = gswroot2
                  gsc = gsw.2.gsc * gsw
                  ci  = ciroot2
               }#end if
            }else if (bounded1){
               gsw = gswroot1
               gsc = gsw.2.gsc * gsw
               ci  = ciroot1
            }else if (bounded2){
               gsw = gswroot2
               gsc = gsw.2.gsc * gsw
               ci  = ciroot2
            }else{
               ans$success = FALSE
               return(ans)
            }#end if
         }else if (bounded1){
            #----- First root is the only one that makes sense. ---------------------------#
            gsw = gswroot1
            gsc = gsw.2.gsc * gsw
            ci  = cs - assim / gsc

         }else if (bounded2){
            #----- Second root is the only one that makes sense. --------------------------#
            gsw = gswroot2
            gsc = gsw.2.gsc * gsw
            ci  = cs - assim / gsc
         }else{
            #----- None of the solutions are bounded.  This solution failed. --------------#
            ans$success = FALSE
            return(ans)
         }#end if
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #   Find the surface water specific humidity and the transpiration.                  #
      #------------------------------------------------------------------------------------#
      ws          = ( gsw * env$wi + env$gbw * env$h2o) / (env$gbw + gsw )
      transp      = env$gbw * ( ws - env$h2o )
      ans$success = TRUE
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     If we have hit this point, then we found an actual answer.                        #
   #---------------------------------------------------------------------------------------#
   ans = modifyList(x=ans,val=list(ci=ci,cs=cs,assim=assim,transp=transp,gsw=gsw,gsc=gsc))
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      This function finds five parameters, namely:                                        #
#  - Vm        : the photosynthetic capacity.                                              #
#  - leaf.resp : the leaf respiration rate.                                                #
#  - compp     : the CO2 compensation point excluding respiration                          #
#  - kco2      : Michaelis-Mentel constant for CO2                                         #
#  - ko2       : Michaelis-Mentel constant for O2                                          #
#------------------------------------------------------------------------------------------#
photo.params <<- function(env,ipft,iphoto){
   #---------------------------------------------------------------------------------------#
   #     Use the Collatz et al. scheme to find Vm without temperature correction, then     #
   # apply the correction.                                                                 #
   #---------------------------------------------------------------------------------------#
   vm.nocorr    = collatz(env$temp,pft$vm0[ipft],pft$vm.base[ipft])
   vm.lnexplow  = pft$vm.decay.e.low [ipft] * (pft$vm.low.temp[ipft]  - env$temp)
   vm.tlow.fun  = 1.0 + exp(vm.lnexplow )
   vm.lnexphigh = pft$vm.decay.e.high[ipft] * (env$temp - pft$vm.high.temp[ipft])
   vm.thigh.fun = 1.0 + exp(vm.lnexphigh)
   vm           = vm.nocorr / (vm.tlow.fun * vm.thigh.fun)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Use the Collatz et al. scheme to find Lr without temperature correction, then     #
   # apply the correction.                                                                 #
   #---------------------------------------------------------------------------------------#
   lr.nocorr    = collatz(env$temp,pft$lr0[ipft],pft$lr.base[ipft])
   lr.lnexplow  = pft$lr.decay.e.low [ipft] * (pft$lr.low.temp[ipft]  - env$temp)
   lr.tlow.fun  = 1.0 + exp(lr.lnexplow )
   lr.lnexphigh = pft$lr.decay.e.high[ipft] * (env$temp - pft$lr.high.temp[ipft])
   lr.thigh.fun = 1.0 + exp(lr.lnexphigh)
   lr           = lr.nocorr / (lr.tlow.fun * lr.thigh.fun)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Find the compensation point and the Michaelis constants.  Make sure to set the    #
   # CO2 compensation point and the Michaelis constant of CO2 zero in case this is a C4    #
   # plant.                                                                                #
   #---------------------------------------------------------------------------------------#
   compp = ifelse( test = pft$pathway[ipft] == 4
                 , yes  = 0.
                 , no   = collatz(env$temp,compp.ref.coll,compp.base.coll)
                 )#end ifelse
   kco2  = ifelse( test = pft$pathway[ipft] == 4
                 , yes  = 0.
                 , no   = collatz(env$temp,kco2.ref.coll,kco2.base.coll)
                 )#end ifelse
   ko2   = collatz(env$temp,ko2.ref.coll,ko2.base.coll)
   #---------------------------------------------------------------------------------------#
   photo = list (vm = vm,leaf.resp=lr,compp=compp,kco2=kco2,ko2=ko2)
   return(photo)
}#end function
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
cifun.eval <<- function(ci,env,ipft,aparms,out.deriv=TRUE){

   #----- Find the CO2 demand. ------------------------------------------------------------#
   assim = assim.rate(ci,aparms)
   #----- Find the stomatal conductance of water. -----------------------------------------#
   gsw   = stom.cond(ci,assim,env)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the function components, then the function evaluation.                       #
   #---------------------------------------------------------------------------------------#
   efun1 = (gsw - pft$b[ipft]) / (pft$m[ipft] * assim)
   efun2 = (env$co2 - env$compp - assim/ env$gbc)
   efun3 = 1.0 + ( env$gbw * (env$wi - env$h2o) / (pft$d0[ipft] * (env$gbw + gsw)))
   fun   = efun1 * efun2 * efun3 - 1.0
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     We also compute the derivative.                                                   #
   #---------------------------------------------------------------------------------------#
   #----- CO2 demand. ---------------------------------------------------------------------#
   assim.prime = assim.rate.prime(ci,assim,aparms)
   #----- stomatal conductance of water. --------------------------------------------------#
   gsw.prime   = stom.cond.prime(ci,gsw,assim,assim.prime,env)
   #----- Function components. ------------------------------------------------------------#
   eprime1 = ( ( gsw.prime * assim - assim.prime * (gsw - pft$b[ipft]) )
             / ( pft$m[ipft] * assim * assim) )
   eprime2 = - assim.prime / env$gbc
   eprime3 = - (efun3 - 1.0) *  gsw.prime / ( env$gbw + gsw )
   dfundci = eprime1 * efun2 * efun3 + efun1 * eprime2 * efun3 + efun1 * efun2 * eprime3
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #      Change the result depending on whether the user wants the derivative or not.     #
   #---------------------------------------------------------------------------------------#
   if (out.deriv){
      ans=list(fun=fun,deriv=dfundci)
   }else{
      ans=fun
   }#end if
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function cifun.eval
#==========================================================================================#
#==========================================================================================#
