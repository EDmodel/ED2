#==========================================================================================#
#==========================================================================================#
#     This function finds the LAI compensation point for any given PFT and PAR.            #
#------------------------------------------------------------------------------------------#
lai.compp = function(partop,parmin,prss,mu,ipft){
   orient = pft$orient.factor  [ipft]
   clump  = pft$clumping.factor[ipft]

   lailow   = 1.e-5
   laihigh  = 30.
   testlow  = mylightroot(x=lailow,orient,clump,partop,parmin,prss,mu)
   testhigh = mylightroot(x=laihigh,orient,clump,partop,parmin,prss,mu)
   while (testlow*testhigh >= 0.){
      lailow = lailow / 5.
      testlow  = mylightroot(x=lailow,orient,clump,partop,parmin,prss,mu)
   }#end while

   lai = uniroot(f=mylightroot,interval=c(lailow,laihigh),orient=orient,clump=clump
                ,partop=partop,parmin=parmin,prss=prss,mu=mu,tol=1.e-7)$root
   return(lai)
}#end for
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function finds the layer transmittance coefficient of diffuse radiation for a   #
# given LAI, clumping factor, and orientation factor.                                      #
#------------------------------------------------------------------------------------------#
comp.tau.diff = function(lai,orient,clump){

    phi1    = 0.5 - orient * (0.633 + 0.33 * orient)
    phi2    = 0.877 - (1 - 2 * phi1)
    elai    = clump * lai
    ediff1  = phi1 * elai
    ediff2  = phi2 * elai
    ediff12 = ediff1 + ediff2

    
    tau    = - exp(-ediff12) * (ediff1^2 * exp(ediff1) * eifun(-ediff1) + (ediff1 - 1))
    return(tau)
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function finds the layer transmittance coefficient of diffuse radiation for a   #
# given LAI, clumping factor, and orientation factor.                                      #
#------------------------------------------------------------------------------------------#
comp.tau.beam = function(lai,orient,clump,mu){

    phi1    = 0.5 - orient * (0.633 + 0.33 * orient)
    phi2    = 0.877 - (1 - 2 * phi1)
    elai    = clump * lai

    
    tau    = exp(- (phi1 + phi2 * mu) * lai / mu)
    return(tau)
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function is zero when LAI makes the light available after extinction to be      #
# exactly the light compensation point.                                                    #
#------------------------------------------------------------------------------------------#
mylightroot = function(x,orient,clump,partop,parmin,prss,mu){

   tau.diff = comp.tau.diff(lai=x,orient=orient,clump=clump)
   tau.beam = comp.tau.beam(lai=x,orient=orient,clump=clump,mu=mu)
   mypar = par.split(cosz=mu,partop=partop,atm.prss=prss)

   myfun = tau.diff * mypar$diff + tau.beam * mypar$beam - parmin
   return(myfun)
}#end if
#==========================================================================================#
#==========================================================================================#
