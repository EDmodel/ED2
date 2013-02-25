#==========================================================================================#
#==========================================================================================#
#    Some global dimensions.                                                               #
#------------------------------------------------------------------------------------------#

#----- Conversion for vegetation height. --------------------------------------------------#
vh2vr      <<- 0.13    # Vegetation height => roughness
vh2dh      <<- 0.63    # Vegetation height => 0-plane displacement height
iphim.clm  <<- 2       # Functional form of phi_m for the very unstable case:
                       # 1. As in CLM:      phi_m = chi_m * cbrt(- zeta)
                       # 2. My alternative: phi_m = chi_m * (-zeta)^(-1/6)
#------------------------------------------------------------------------------------------#



#----- Global constants. ------------------------------------------------------------------#
tprandtl   <<-  0.74
gamm       <<-  13.
gamh       <<-  13.
beta.s     <<-  5.0
#------------------------------------------------------------------------------------------#




#----- Parameters for Louis (1979) --------------------------------------------------------#
bb         <<- 5.0
csm        <<- 7.5
csh        <<- 5.0
dd         <<- 5.0
#------------------------------------------------------------------------------------------#





#----- Global constants. ------------------------------------------------------------------#
z0moz0h    <<-  1.
z0hoz0m    <<-  1. / z0moz0h 
#------------------------------------------------------------------------------------------#





#----- Constants for the Beljaars-Holtslag model. -----------------------------------------#
abh        <<- -1.0
bbh        <<- -2./3.
cbh        <<- 5.0
dbh        <<- 0.35
ebh        <<- -2./3.
fbh        <<- 1.5
#------------------------------------------------------------------------------------------#





#----- Assign the critical value for zeta. ------------------------------------------------#
zetac.um  <<- -1.5
zetac.uh  <<- -0.5
zetac.sm  <<-  1.0
zetac.sh  <<-  zetac.sm

#----- Define the coefficients for the very unstable and very stable cases. ---------------#
if (iphim.clm == 1){
   chim   <<- 1. / (cbrt(-zetac.um) * sqrt(sqrt(1.0 - gamm * zetac.um)))
}else if(iphim.clm == 2){
   chim   <<- sqrt(cbrt(-zetac.um)) / sqrt(sqrt(1.0 - gamm * zetac.um))
}#end if
chih      <<- cbrt(-zetac.uh) / sqrt(1.0 - gamh * zetac.uh)
beta.vs   <<- 1.0 - (1.0 - beta.s) * zetac.sm


#------------------------------------------------------------------------------------------#
#     Initialise these values with dummies, it will be updated after we define the         #
# functions.                                                                               #
#------------------------------------------------------------------------------------------#
psimc.um  <<- 0.
psihc.uh  <<- 0.

#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function computes the correction function for non-neutral conditions for        #
# momentum.                                                                                #
#------------------------------------------------------------------------------------------#
psim =  function(method,zeta){

   if (method == 2 | method == 3){
      stab  = zeta >= 0.
      unst  = zeta <  0.
   }else if (method == 4){
      vuns  = zeta <  zetac.um
      unst  = zeta <  0        & zeta >= zetac.um
      stab  = zeta >= 0        & zeta <= zetac.sm
      vsta  = zeta >  zetac.sm
   }#end if
   
   x        =  0. * zeta
   psimloc  =  0. * zeta


   x[stab]  =  zeta[stab]
   x[unst]  =  (1.0-gamm*zeta[unst])^0.25

   if (method == 2){
      psimloc[stab] = - beta.s * zeta[stab]
      psimloc[unst] = ( log(0.125*(1.0+x[unst])*(1.0+x[unst])*(1.0+x[unst]*x[unst])) 
                      - 2.0 * atan(x[unst]) + 0.5 * pi)

   }else if (method == 3){
      psimloc[stab] = ( abh * zeta[stab]
                      + bbh * (zeta[stab] - cbh/dbh) * exp(-dbh*zeta[stab]) 
                      + bbh * cbh / dbh )
      psimloc[unst] = ( log(0.125*(1.0+x[unst])*(1.0+x[unst])*(1.0+x[unst]*x[unst])) 
                      - 2.0 * atan(x[unst]) + 0.5 * pi )
   }else if (method == 4){
      psimloc[vsta] = ( (1.0 - beta.vs) * log(zeta[vsta]/zetac.sm)
                      + (1.0 - beta.s ) * zetac.sm - zeta[vsta] )
      psimloc[stab] = - beta.s * zeta[stab]
      psimloc[unst] = ( log(0.125*(1.0+x[unst])*(1.0+x[unst])*(1.0+x[unst]*x[unst])) 
                      - 2.0 * atan(x[unst]) + 0.5 * pi)
      if (iphim.clm == 1){
         psimloc[vuns] = ( log(zeta[vuns]/zetac.um)
                         + 3.0 * chim * (cbrt(zeta[vuns]) - cbrt(zetac.um)) + psimc.um)
      }else if (iphim.clm == 2){
         psimloc[vuns] = ( log(zeta[vuns]/zetac.um)
                         + 6.0 * chim  
                         * (1.0/sqrt(cbrt(-zeta[vuns])) - 1.0/sqrt(cbrt(-zetac.um))) 
                         + psimc.um )
      }#end if
   }#end if

   return(psimloc)
} # end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function computes the correction function for non-neutral conditions for        #
# heat.                                                                                    #
#------------------------------------------------------------------------------------------#
psih =  function(method,zeta){

   if (method == 2 | method == 3){
      stab  = zeta >= 0.
      unst  = zeta <  0.
   }else if (method == 4){
      vuns  = zeta <  zetac.uh
      unst  = zeta <  0        & zeta >= zetac.uh
      stab  = zeta >= 0        & zeta <= zetac.sh
      vsta  = zeta >  zetac.sh
   }#end if

   y        =  0. * zeta
   psihloc  =  0. * zeta

   y[stab]  =  zeta[stab]
   y[unst]  =  (1.0-gamh*zeta[unst])^0.50

   if (method == 2){
      psihloc[stab] =  - beta.s * zeta[stab]
      psihloc[unst] =  log(0.25*(1.0+y[unst])*(1.0+y[unst]))

   }else if(method == 3){
      psihloc[stab] =  ( 1. - (1. + abh * ebh * zeta[stab])^fbh
                       + bbh *(zeta[stab] - cbh / dbh) * exp(-dbh*zeta[stab])
                       + bbh * cbh / dbh )
      psihloc[unst]  =  log(0.25*(1.0+y[unst])*(1.0+y[unst]))

   }else if (method == 4){
      psihloc[vsta] = ( (1.0 - beta.vs) * log(zeta[vsta]/zetac.sh)
                      + (1.0 - beta.s ) * zetac.sh - zeta[vsta] )
      psihloc[stab] =  - beta.s * zeta[stab]
      psihloc[unst] =  log(0.25*(1.0+y[unst])*(1.0+y[unst]))
      psihloc[vuns] = ( log(zeta[vuns]/zetac.uh)
                      + 3.0 * chih  * (1./cbrt(zetac.uh) - 1.0/cbrt(zeta[vuns])) 
                      + psihc.uh )

   }#end if
   return(psihloc)
} # end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function computes the derivative of the correction function for non-neutral     #
# conditions for momentum.                                                                 #
#------------------------------------------------------------------------------------------#
dpsimdzeta =  function(method,zeta){

   if (method == 2 | method == 3){
      stab  = zeta >= 0.
      unst  = zeta <  0.
   }else if (method == 4){
      vuns  = zeta <  zetac.uh
      unst  = zeta <  0        & zeta >= zetac.uh
      stab  = zeta >= 0        & zeta <= zetac.sh
      vsta  = zeta >  zetac.sh
   }#end if

   xx          =  0. * zeta
   dpsimdzloc  =  0. * zeta


   xx[stab]  =  zeta[stab]
   xx[unst]  =  sqrt(sqrt(1.0-gamm*zeta[unst]))

   if (method == 2){
      dpsimdzloc[stab] = - beta.s + 0.*zeta[stab]
      # dpsimdzloc[unst] = - gamm / (xx[unst]*(1.0+xx[unst])*(1.0+xx[unst]*xx[unst]))
      dpsimdzloc[unst] = (1.0 - 1.0/xx[unst]) / zeta[unst]

   }else if (method == 3){
      dpsimdzloc[stab] = ( abh 
                         + bbh*(1. - dbh*zeta[stab] + cbh)*exp(-dbh*zeta[stab]) )
      # dpsimdzloc[unst] =  - gamm / (xx[unst]*(1.0+xx[unst])*(1.0+xx[unst]*xx[unst]))
      dpsimdzloc[unst] = (1.0 - 1.0/xx[unst]) / zeta[unst]

   }else if (method == 4){
      dpsimdzloc[vsta] = (1.0 - beta.vs) / zeta[vsta] - 1.0
      dpsimdzloc[stab] = - beta.s + 0.*zeta[stab]
      # dpsimdzloc[unst] = - gamm / (xx[unst]*(1.0+xx[unst])*(1.0+xx[unst]*xx[unst]))
      dpsimdzloc[unst] = (1.0 - 1.0/xx[unst]) / zeta[unst]
      if (iphim.clm == 1){
         dpsimdzloc[vuns] = (1.0 - chim * cbrt(-zeta[vuns])) / zeta[vuns]
      }else if (iphim.clm == 2){
         dpsimdzloc[vuns] = (1.0 - chim / (-zeta[vuns])^onesixth) / zeta[vuns]
      }#end if
   }
   return(dpsimdzloc)
} # end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function computes the derivative of the correction function for non-neutral     #
# conditions for heat.                                                                     #
#------------------------------------------------------------------------------------------#
dpsihdzeta =  function(method,zeta){

   if (method == 2 | method == 3){
      stab  = zeta >= 0.
      unst  = zeta <  0.
   }else if (method == 4){
      vuns  = zeta <  zetac.uh
      unst  = zeta <  0        & zeta >= zetac.uh
      stab  = zeta >= 0        & zeta <= zetac.sh
      vsta  = zeta >  zetac.sh
   }#end if

   yy          =  0. * zeta
   dpsihdzloc  =  0. * zeta

   yy[stab]  =  zeta[stab]
   yy[unst]  =  sqrt(1.0-gamh*zeta[unst])

   if (method == 2){
      dpsihdzloc[stab] =  - beta.s + 0.*zeta[stab]
      # dpsihdzloc[unst] =  - gamh / (yy[unst] * (1.0 + yy[unst]))
      dpsihdzloc[unst]   = (1.0 - 1.0/yy[unst]) / zeta[unst]
   }else if(method == 3){
      dpsihdzloc[stab] = ( -fbh * abh * ebh * ((1. + abh * ebh * zeta[stab])^(fbh-1))
                         + bbh * (1. - dbh * zeta[stab] + cbh) * exp(-dbh*zeta[stab]) )
      # dpsihdzloc[unst] = - gamh / (yy[unst] * (1.0 + yy[unst]))
      dpsihdzloc[unst]   = (1.0 - 1.0/yy[unst]) / zeta[unst]
   }else if(method == 4){
      dpsihdzloc[vsta] = (1.0 - beta.vs) / zeta[vsta] - 1.0
      dpsihdzloc[stab] =  - beta.s + 0.*zeta[stab]
      # dpsihdzloc[unst] =  - gamh / (yy[unst] * (1.0 + yy[unst]))
      dpsihdzloc[unst]   = (1.0 - 1.0/yy[unst]) / zeta[unst]
      dpsihdzloc[vuns] = (1 - chih / cbrt(zeta[vuns]))/zeta[vuns]
   }
   return(dpsihdzloc)
} # end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This subroutine combutes the residual of z/L, which is used to find the normalised   #
# height zeta.                                                                             #
#------------------------------------------------------------------------------------------#
zolmzeta =  function(zeta,method,rib,zstar,z0m,z0h){
  zeta0m     =  z0m*zeta/zstar
  zeta0h     =  z0h*zeta/zstar
  fm         =  log(zref/z0m) - psim(method,zeta) + psim(method,zeta0m)
  fh         =  log(zref/z0h) - psih(method,zeta) + psih(method,zeta0h)



   #---------------------------------------------------------------------------------------!
   #     Define the coefficient Ri * zstar / [Pr * (zstar-z0)]                             !
   #---------------------------------------------------------------------------------------!
   coeff = rib * zstar / (tprandtl * (zstar - z0m))
   #---------------------------------------------------------------------------------------!

  fff        =  coeff * fm * fm / fh - zeta
  if (! is.finite(fff)) browser()
  return(fff)
} #end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function finds the value of zeta for a given Richardson number, reference       #
# height and the roughness scale.  This is solved by using the definition of Obukhov       #
# length scale as stated in Louis (1979) equation (10), modified to define z/L rather      #
# than L.  The solution is found  iteratively since it's not a simple function to          #
# invert.  It tries to use Newton's method, which should take care of most cases.  In      #
# the unlikely case in which Newton's method fails, switch back to modified Regula         #
# Falsi method (Illinois).                                                                 #
#------------------------------------------------------------------------------------------#
zoobukhov = function (method,rib,zstar,rough,zoz0m,lnzoz0m,zoz0h,lnzoz0h){


   #---------------------------------------------------------------------------------------#
   #      Save the sign of the Richardson number.                                          #
   #---------------------------------------------------------------------------------------#
   stable = rib >= 0.
   #---------------------------------------------------------------------------------------#


   #----- Define some values that won't change during the iterative method. ---------------#
   z0moz = 1. / zoz0m
   z0hoz = 1. / zoz0h
   #---------------------------------------------------------------------------------------#

   #----- Find the maximum acceptable Richardson number for methods 2 and 4. --------------#
   if (method == 3){
      ribuse = rib
   }else{
      ribuse = min(rib,(1.0-toler) * tprandtl / beta.s )
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------!
   #     Define the coefficient Ri * zstar / [Pr * (zstar-z0)]                             !
   #---------------------------------------------------------------------------------------!
   coeff = ribuse * zstar / (tprandtl * (zstar - rough))
   #---------------------------------------------------------------------------------------!


   #---------------------------------------------------------------------------------------#
   #     If the bulk Richardson number is zero or almost zero, then we rather just         #
   # assign z/L to be the one similar to Oncley and Dudhia (1995).  This saves time and    #
   # also avoids the risk of having zeta with the opposite sign.                           #
   #---------------------------------------------------------------------------------------#
   zetasmall = rib * min(lnzoz0m,lnzoz0h) / tprandtl
   if (ribuse <= 0. && zetasmall > - z0moz0h * toler){
      zeta = zetasmall
      return(zeta)
   }else if (ribuse > 0. && zetasmall < z0moz0h * toler){
      zeta = zetasmall / (1.0 - beta.s * ribuse / tprandtl)
      return(zeta)
   }else{
      zetamin    =  toler
      zetamax    = -toler
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     First guess, using Oncley and Dudhia (1995) approximation for unstable case.      #
   # We won't use the stable case to avoid FPE or zeta with opposite sign when Ri > 0.20.  #
   #---------------------------------------------------------------------------------------#
   zetaa = ribuse * lnzoz0m / tprandtl
   #---------------------------------------------------------------------------------------#



   #----- Find the function and its derivative. -------------------------------------------#
   zeta0m   = zetaa * z0moz
   zeta0h   = zetaa * z0hoz
   fm       = lnzoz0m - psim(method,zetaa) + psim(method,zeta0m)
   fh       = lnzoz0h - psih(method,zetaa) + psih(method,zeta0h)
   dfmdzeta = z0moz * dpsimdzeta(method,zeta0m) - dpsimdzeta(method,zetaa)
   dfhdzeta = z0hoz * dpsihdzeta(method,zeta0h) - dpsihdzeta(method,zetaa)
   funa     = coeff * fm * fm / fh - zetaa
   myderiv  = coeff * (2. * fm * dfmdzeta * fh - fm * fm * dfhdzeta) / (fh * fh) - 1.
   #---------------------------------------------------------------------------------------#


   #----- Copy just in case it fails at the first iteration. ------------------------------#
   zetaz = zetaa
   fun   = funa
   #---------------------------------------------------------------------------------------#



   #----- Enter Newton's method loop. -----------------------------------------------------#
   converged = FALSE
   diverged  = FALSE
   itn = 0
   while ((! converged) && (! diverged) && itn < floor(maxfpo/6)){
      itn = itn + 1

      #------------------------------------------------------------------------------------#
      #     Newton's method converges fast when it's on the right track, but there are     #
      # cases in which it becomes ill-behaved.  Two situations are known to cause trouble: #
      # 1.  If the derivative is tiny, the next guess can be too far from the actual       #
      #     answer;                                                                        #
      # 2.  For this specific problem, when zeta is too close to zero.  In this case the   #
      #     derivative will tend to infinity at this point and Newton's method is not go-  #
      #     ing to perform well and can potentially enter in a weird behaviour or lead to  #
      #     the wrong answer.  In any case, so we rather go with bisection.                #
      #------------------------------------------------------------------------------------#
      flat      = abs(myderiv) < toler
      crossing0 = ( ( stable && (zetaz - fun/myderiv < zetamin))    ||
                    ((! stable) && (zetaz - fun/myderiv > zetamax)) )
      diverged  = flat || crossing0
      if (! diverged){

         #----- Copying the previous guess ------------------------------------------------#
         zetaa    = zetaz
         funa     = fun
         #----- New guess, its function and derivative evaluation -------------------------#
         zetaz    = zetaa - fun/myderiv
         zeta0m   = zetaz * z0moz
         zeta0h   = zetaz * z0hoz
         fm       = lnzoz0m - psim(method,zetaz) + psim(method,zeta0m)
         fh       = lnzoz0h - psih(method,zetaz) + psih(method,zeta0h)
         dfmdzeta = z0moz * dpsimdzeta(method,zeta0m) - dpsimdzeta(method,zetaz)
         dfhdzeta = z0hoz * dpsihdzeta(method,zeta0h) - dpsihdzeta(method,zetaz)
         fun      = coeff * fm * fm / fh - zetaz
         myderiv  = coeff * (2. * fm * dfmdzeta * fh - fm * fm * dfhdzeta) / ( fh * fh) - 1.

         converged = abs(zetaz-zetaa) < toler * abs(zetaz)

         if (converged){
            zeta = 0.5 * (zetaa+zetaz)
            return(zeta)
         }else if (fun == 0.0){ #---- Converged by luck. ----------------------------------#
            zeta = zetaz
            return(zeta)
         }#end if (converged)
      }#end if (! diverged)
   }#end while
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     If we reached this point then it's because Newton's method failed or it has       #
   # become too dangerous.  We use the Regula Falsi (Illinois) method, which is just a     #
   # fancier bisection.  For this we need two guesses, and the guesses must have opposite  #
   #signs.                                                                                 #
   #---------------------------------------------------------------------------------------#
   if (funa * fun < 0.0){
      funz  = fun
      zside = TRUE
   }else{
      if (abs(fun-funa) < 100. * toler * abs(zetaa)){
         if (stable){
            delta = max(0.5 * abs(zetaa-zetamin),100. * toler * abs(zetaa))
         }else{
            delta = max(0.5 * abs(zetaa-zetamax),100. * toler * abs(zetaa))
         }#end if
      }else{
         if (stable){
            delta = max(abs(funa * (zetaz-zetaa)/(fun-funa))
                       ,100. * toler * abs(zetaa)
                       ,0.5 * abs(zetaa-zetamin))
         }else{
            delta = max(abs(funa * (zetaz-zetaa)/(fun-funa))
                       ,100. * toler * abs(zetaa)
                       ,0.5 * abs(zetaa-zetamax))
         }#end if
      }#end if
      if (stable){
         zetaz = max(zetamin,zetaa + delta)
      }else{
         zetaz = min(zetamax,zetaa + delta)
      }#end if
      zside = FALSE
      itp   = 0
      while ((! zside) && itp < maxfpo){
         itp = itp + 1

         if (stable){
            zetaz    = max(zetamin,zetaa + ((-1)**itp * (itp+3)/2) * delta)
         }else{
            zetaz    = min(zetamax,zetaa + ((-1)**itp * (itp+3)/2) * delta)
         }#end if
         zeta0m   = zetaz * z0moz
         zeta0h   = zetaz * z0hoz
         fm       = lnzoz0m - psim(method,zetaz) + psim(method,zeta0m)
         fh       = lnzoz0h - psih(method,zetaz) + psih(method,zeta0h)
         funz     = coeff * fm * fm / fh - zetaz
         zside    = funa * funz < 0.0
      }#end while

      if (! zside){
         print(paste("=================================================="))
         print(paste("            No second guess for you...            "))
         print(paste("=================================================="))
         print(paste(" ZSTAR   =",zstar                                  ))
         print(paste(" ROUGH   =",rough                                  ))
         print(paste(" LNZOZ0M =",lnzoz0m                                ))
         print(paste(" LNZOZ0H =",lnzoz0h                                ))
         print(paste(" RIB     =",rib                                    ))
         print(paste(" RIBUSE  =",ribuse                                 ))
         print(paste(" STABLE  =",stable                                 ))
         print(paste(" FUN     =",fun                                    ))
         print(paste(" DELTA   =",delta                                  ))
         print(paste(" ZETAA   =",zetaa                                  ))
         print(paste(" FUNA    =",funa                                   ))
         print(paste(" ZETAZ   =",zetaz                                  ))
         print(paste(" FUNZ    =",funz                                   ))
         stop("Failed finding a second guess in zoobukhov, sorry...")
      }#end if
   }#end if
   #---------------------------------------------------------------------------------------#




   #----- Now we are ready to start the regula falsi method. ------------------------------#
   itb = itn-1
   while((! converged) & itb < maxfpo){
      zeta = (funz*zetaa-funa*zetaz)/(funz-funa)

      #------------------------------------------------------------------------------------#
      #     Now that we updated the guess, check whether they are really close. If so, it  #
      # converged, I can use this as my guess.                                             #
      #------------------------------------------------------------------------------------#
      converged = abs(zeta-zetaa) < toler * abs(zeta)
      if (! converged){

         #------ Find the new function ----------------------------------------------------#
         zeta0m   = zeta * z0moz
         zeta0h   = zeta * z0hoz
         fm       = lnzoz0m - psim(method,zeta) + psim(method,zeta0m)
         fh       = lnzoz0h - psih(method,zeta) + psih(method,zeta0h)
         fun      = coeff * fm * fm / fh - zeta

         #------ Define the new interval based on the intermediate value theorem. ---------#
         if (fun*funa < 0. ){
            zetaz = zeta
            funz  = fun
            #----- If we are updating zside again, modify aside (Illinois method) ---------#
            if (zside) funa = funa * 0.5
            #----- We just updated zside, setting zside to true. --------------------------#
            zside = TRUE
         }else{
            zetaa = zeta
            funa  = fun
            #----- If we are updating aside again, modify aside (Illinois method) ---------#
            if (! zside) funz = funz * 0.5
            #----- We just updated aside, setting aside to true. --------------------------#
            zside = FALSE
         }#end if
      }#end if
   }#end while
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   if (converged){
      return(zeta)
   }else{
      print(paste("-------------------------------------------------------"))
      print(paste(" Zeta finding didn't converge!!!"                       ))
      print(paste(" I gave up, after",maxfpo,"iterations..."               ))
      print(paste(" "                                                      ))
      print(paste(" Input values."                                         ))
      print(paste(" "                                                      ))
      print(paste("RIB             [   ---] =",rib                         ))
      print(paste("RIBUSE          [   ---] =",ribuse                      ))
      print(paste("ZSTAR           [     m] =",zstar                       ))
      print(paste("ROUGH           [     m] =",rough                       ))
      print(paste("ZOZ0M           [   ---] =",zoz0m                       ))
      print(paste("LNZOZ0M         [   ---] =",lnzoz0m                     ))
      print(paste("ZOZ0H           [   ---] =",zoz0h                       ))
      print(paste("LNZOZ0H         [   ---] =",lnzoz0h                     ))
      print(paste("STABLE          [   T|F] =",stable                      ))
      print(paste(" "                                                      ))
      print(paste(" Last iteration outcome (downdraft values)."            ))
      print(paste("ZETAA           [   ---] =",zetaa                       ))
      print(paste("ZETAZ           [   ---] =",zetaz                       ))
      print(paste("FUN             [   ---] =",fun                         ))
      print(paste("FM              [   ---] =",fm                          ))
      print(paste("FH              [   ---] =",fh                          ))
      print(paste("FUNA            [   ---] =",funa                        ))
      print(paste("FUNZ            [   ---] =",funz                        ))
      print(paste("MYDERIV         [   ---] =",myderiv                     ))
      print(paste("TOLER           [   ---] =",toler                       ))
      print(paste("ERROR           [   ---] =",abs(zetaz-zetaa)/abs(zetaz) ))
      print(paste("ZETA            [   ---] =",zeta                        ))
      print(paste("-------------------------------------------------------"))

      stop("Zeta didn't converge, giving up!")
   }# end if
   #---------------------------------------------------------------------------------------#

}#end function zoobukhov
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This subroutine combutes the residual of z/L, which is used to find the normalised   #
# height zeta.                                                                             #
#------------------------------------------------------------------------------------------#
ed.stars = function(rib,uspd,dthetav,thetav,dens,zref,z0m,z0h,ustmin,ribmin,ribmax,method){

   #---------------------------------------------------------------------------------------#
   #     Make sure all wind speed, and theta v and the gradient have the same dimension.   #
   #---------------------------------------------------------------------------------------#
   uuse    = uspd    + 0. * rib
   thvuse  = thetav  + 0. * rib
   dthvuse = dthetav + 0. * rib
   ribuse  = rib
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Change uspd for those cases in which Rib > Ribmax.                                #
   #---------------------------------------------------------------------------------------#
   unst      = rib < 0.
   stab      = rib >= 0.
   if (method == 2){
      ribmax.use = min(ribmax,0.3)
   }else{
      ribmax.use = ribmax
   }#end if
   if (method == 4 && iphim.clm == 1){
      ribmin.use = max(ribmin,-4.0)
   }else{
      ribmin.use = ribmin
   }#end if

   toostable = ribuse > ribmax.use
   uuse   [  toostable] = sqrt(ribuse[  toostable] / ribmax.use) * uuse[  toostable]
   ribuse [  toostable] = ribmax.use

   toounstable = ribuse < ribmin.use
   uuse   [toounstable] = sqrt(ribuse[toounstable] / ribmin.use) * uuse[toounstable]
   ribuse [toounstable] = ribmin.use
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the ratios between heights.                                                  #
   #---------------------------------------------------------------------------------------#
   zoz0m   =  zref / z0m
   zoz0h   =  zref / z0h
   lnzoz0m =  log(zref / z0m)
   lnzoz0h =  log(zref / z0h)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Select which method to use.                                                       #
   #---------------------------------------------------------------------------------------#
   if (method == 1){
      #------------------------------------------------------------------------------------#
      #    Louis (1979).                                                                   #
      #------------------------------------------------------------------------------------#
      #----- Find some useful parameters. -------------------------------------------------#
      a2     =  vonk*vonk/(lnzoz0m*lnzoz0m)
      c1     =  a2 * uuse
      ee     =  zoz0m^(1./3.) - 1.
      c2     =  bb * a2 * ee * sqrt(ee * abs(ribuse))
      cmo    =  csm * c2
      che    =  csh * c2
      #----- Determine the ones that should be stable and unstable. -----------------------#
      stab   = rib >  +1.e-5
      unst   = rib <  -1.e-5
      #----- Find the momentum flux corrections. ------------------------------------------#
      fm         =  1.0 + 0. * rib
      fm[unst]   =  (1.0 - 2.0 * bb * ribuse[unst] / (1.0 + 2.0 * cmo[unst]))
      fm[stab]   =  1.0 / (1.0 + (2.0 * bb * ribuse[stab] / sqrt(1.0 + dd * ribuse[stab])))
      #----- Find the heat flux corrections. ----------------------------------------------#
      fh         =  1.0 + 0. * rib
      fh[unst]   =  (1.0 - 3.0 * bb * ribuse[unst] / (1.0 + 3.0 * che[unst]))
      fh[stab]   =  1.0 / (1.0 + (3.0 * bb * ribuse[stab] * sqrt(1.0 + dd * ribuse[stab])))
      #----- Find u*, and impose minimum. -------------------------------------------------#
      ustar                 = c1 * uuse * fm
      ustar[ustar < ustmin] = ustmin
      #----- Find theta*. -----------------------------------------------------------------#
      tstar = c1 * fh * dthvuse / ustar
      #----- Find a diagnostic zeta. ------------------------------------------------------#
      zeta  = vonk * grav * tstar / (thvuse * ustar * ustar)
      #----- Sensible heat flux. ----------------------------------------------------------#
      sensible = - dens * cp * ustar * tstar
      #------------------------------------------------------------------------------------#

   }else{
      #------------------------------------------------------------------------------------#
      #     Any other method.  All the others are solved with the same general method, the #
      # difference is only in the functions psim and psih.                                 #
      #------------------------------------------------------------------------------------#
      #----- Find a first guess interval. -------------------------------------------------#
      zeta1st        =  0. * ribuse
      zeta1st[unst]  =  ribuse[unst] * lnzoz0m
      stab1 = stab & ribuse[stab] <= 0.19
      stab2 = stab & ribuse[stab] >  0.19
      zeta1st[stab1]  =  ribuse[stab1] * lnzoz0m/(1.1 - beta.s * ribuse[stab1])
      zeta1st[stab2]  =  ribuse[stab2] * lnzoz0m/(1.1 - beta.s * 0.19)


      zeta   =  NA + 0.*ribuse
      for (n in 1:length(ribuse)){
         if (ribuse[n] == 0.){
            zeta[n] = 0.
         }else{
            #------ Find two guesses with opposite signs. ---------------------------------#
            # delta0   = 1.1
            # ntry     = 0
            # found    = FALSE
            # while (! found){
            #    ntry     = ntry + 1
            #    delta    = delta0 ^ ntry
            #    zleft    = min(zeta1st[n]/delta,zeta1st[n]*delta)
            #    zright   = max(zeta1st[n]/delta,zeta1st[n]*delta)
            #    funleft  = zolmzeta(zleft,method,ribuse[n],zref,z0m,z0h)
            #    funright = zolmzeta(zright,method,ribuse[n],zref,z0m,z0h)
            #    found    = funleft * funright < 0.
            # }#end while

            # myroot  =  uniroot(f=zolmzeta,interval=c(zleft,zright),method=method
            #                   ,rib=ribuse[n],zref=zref,z0m=z0m,z0h=z0h)
            # zeta[n] = myroot$root

            zeta[n] = zoobukhov(method=method,rib=ribuse[n],zref=zref,rough=z0m
                               ,zoz0m=zoz0m,lnzoz0m=lnzoz0m,zoz0h=zoz0h,lnzoz0h=lnzoz0h)
         }#end if
      }#end for

      #----- Now that we have solved zeta, find the actual correction functions. ----------#
      zeta0m = z0m * zeta / zref
      zeta0h = z0h * zeta / zref
      psim1  = psim(method,zeta  )
      psim0  = psim(method,zeta0m)
      psih1  = psih(method,zeta  )
      psih0  = psih(method,zeta0h)
      #------------------------------------------------------------------------------------#



      #----- Determine u*. ----------------------------------------------------------------#
      ustar                 = vonk * uuse    / (lnzoz0m - psim1 + psim0)
      ustar[ustar < ustmin] = ustmin 
      #------------------------------------------------------------------------------------#



      #----- Determine theta*. ------------------------------------------------------------#
      tstar                 = vonk * dthvuse / (tprandtl * (lnzoz0h - psih1 + psih0))
      #------------------------------------------------------------------------------------#



      #----- Sensible heat flux. ----------------------------------------------------------#
      sensible = - dens * cp * ustar * tstar
      #------------------------------------------------------------------------------------#
   }#end if

   ans   = list(ustar=ustar,tstar=tstar,zeta=zeta,sensible=sensible)
   return(ans)
} #end function
#==========================================================================================#
#==========================================================================================#





#----- Find the psi_critical. -------------------------------------------------------------#
psimc.um  <<- psim(method=4,zeta=zetac.um)
psihc.uh  <<- psih(method=4,zeta=zetac.uh)
#------------------------------------------------------------------------------------------#
