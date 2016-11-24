#------------------------------------------------------------------------------------------#
#      Parameters for the aerodynamic resistance between the leaf (flat surface) and       #
# wood (kind of cylinder surface), and the canopy air space.  These are the A, B, n,       #
# and m parameters that define the Nusselt number for forced and free convection, at       #
# equations 10.7 and 10.9.  The parameters are found at the appendix A, table A.5(a)       #
# and A.5(b).                                                                              #
#                                                                                          #
# M08 - Monteith, J. L., M. H. Unsworth, 2008. Principles of Environmental Physics,        #
#       3rd. edition, Academic Press, Amsterdam, 418pp.  (Mostly Chapter 10).              #
#                                                                                          #
# The coefficient B for flat plates under turbulent flow was changed to 0.19 so the        #
#     transition from laminar to turbulent regime will happen at Gr ~ 100,000, the         #
#     number suggested by M08.                                                             #
#------------------------------------------------------------------------------------------#
gbhmos.min <<- 1.e-9    # Minimum conductance (m/s)
aflat.lami <<- 0.600    # A (forced convection), laminar   flow
nflat.lami <<- 0.500    # n (forced convection), laminar   flow
aflat.turb <<- 0.032    # A (forced convection), turbulent flow
nflat.turb <<- 0.800    # n (forced convection), turbulent flow
bflat.lami <<- 0.500    # B (free   convection), laminar   flow
mflat.lami <<- 0.250    # m (free   convection), laminar   flow
bflat.turb <<- 0.190    # B (free   convection), turbulent flow
mflat.turb <<- onethird # m (free   convection), turbulent flow
ocyli.lami <<- 0.320    # intercept (forced convection), laminar   flow
acyli.lami <<- 0.510    # A (forced convection), laminar   flow
ncyli.lami <<- 0.520    # n (forced convection), laminar   flow
ocyli.turb <<- 0.000    # intercept (forced convection), turbulent flow
acyli.turb <<- 0.240    # A (forced convection), turbulent flow
ncyli.turb <<- 0.600    # n (forced convection), turbulent flow
bcyli.lami <<- 0.480    # B (free   convection), laminar   flow
mcyli.lami <<- 0.250    # m (free   convection), laminar   flow
bcyli.turb <<- 0.090    # B (free   convection), turbulent flow
mcyli.turb <<- onethird # m (free   convection), turbulent flow
#------------------------------------------------------------------------------------------#


#==========================================================================================#
#==========================================================================================#
#     This sub-routine computes the aerodynamic conductance between leaf and canopy        #
# air space for both heat and water vapour, based on:                                      #
#                                                                                          #
# L95 - Leuning, R., F. M. Kelliher, D. G. G. de Pury, E. D. Schulze, 1995: Leaf           #
#       nitrogen, photosynthesis, conductance and transpiration: scaling from leaves to    #
#       canopies.  Plant, Cell and Environ., 18, 1183-1200.                                #
# M08 - Monteith, J. L., M. H. Unsworth, 2008. Principles of Environmental Physics,        #
#       3rd. edition, Academic Press, Amsterdam, 418pp.  (Mostly Chapter 10).              #
#                                                                                          #
# Notice that the units are somewhat different from L95.                                   #
# - gbw is in kg_H2O/m2/s.                                                                 #
#------------------------------------------------------------------------------------------#
leaf.gbw.fun <<- function(wind,can.rhos,can.temp,leaf.temp,ipft){

   #----- Save the leaf width of this PFT. ------------------------------------------------#
   lwidth = pft$leaf.width[ipft]
   #---------------------------------------------------------------------------------------#

   #----- Define both the Reynolds and the Grashof numbers for the given variables. -------#
   th.expan = 1.0 / can.temp
   kin.visc = kin.visc.0 * ( 1.0 + dkin.visc * (can.temp - t00) )
   th.diff  = th.diff.0  * ( 1.0 + dth.diff  * (can.temp - t00) )
   gr.coeff = th.expan * grav / kin.visc^2
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Find the conductance, in m/s, associated with forced convection.                  #
   #---------------------------------------------------------------------------------------#
   #----- 1. Compute the Reynolds number. -------------------------------------------------#
   reynolds        = wind * lwidth / th.diff
   #----- 2. Compute the Nusselt number for both the laminar and turbulent case. ----------#
   nusselt.lami    = aflat.lami * reynolds ^ nflat.lami
   nusselt.turb    = aflat.turb * reynolds ^ nflat.turb
   #----- 3. The right Nusselt number is the largest. -------------------------------------#
   nusselt.forced  = pmax(nusselt.lami,nusselt.turb)
   #----- 4. The conductance is given by MU08 - equation 10.4 -----------------------------#
   forced.gbh.mos  = th.diff * nusselt.forced / lwidth
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the conductance, in m/s,  associated with free convection.                   #
   #---------------------------------------------------------------------------------------#
   #----- 1. Find the Grashof number. -----------------------------------------------------#
   grashof         = gr.coeff  * abs(leaf.temp - can.temp) * lwidth ^ 3
   #----- 2. Compute the Nusselt number for both the laminar and turbulent case. ----------#
   nusselt.lami    = bflat.lami * grashof ^ mflat.lami
   nusselt.turb    = bflat.turb * grashof ^ mflat.turb
   #----- 3. The right Nusselt number is the largest. -------------------------------------#
   nusselt.free    = pmax(nusselt.lami,nusselt.turb)
   #----- 4. The conductance is given by MU08 - equation 10.4 -----------------------------#
   free.gbh.mos    = th.diff * nusselt.free / lwidth
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     The heat conductance for the thermodynamic budget is the sum of conductances,     #
   # because we assume both forms of convection happen parallelly.  The conversion from    #
   # heat to water conductance (in m/s) can be found in L95, page 1198, after equation     #
   # E5.  For the ED purposes, the output variables are converted to the units of          #
   # entropy and water fluxes [J/K/m2/s and kg/m2/s, respectively].                        #
   #---------------------------------------------------------------------------------------#
   gbh.mos  = pmax(gbhmos.min, free.gbh.mos + forced.gbh.mos)
   leaf.gbw = gbh.2.gbw * gbh.mos * can.rhos
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   return(leaf.gbw)
   #---------------------------------------------------------------------------------------#
}#end function
#==========================================================================================#
#==========================================================================================#
