#==========================================================================================#
#==========================================================================================#
#     This module contains several parameters used in the Farquar Leuning photosynthesis   #
# solver.  The following references are used as a departing point:                         #
# - M09 - Medvigy, D.M., S. C. Wofsy, J. W. Munger, D. Y. Hollinger, P. R. Moorcroft,      #
#         2009: Mechanistic scaling of ecosystem function and dynamics in space and time:  #
#         Ecosystem Demography model version 2.  J. Geophys. Res., 114, G01002,            #
#         doi:10.1029/2008JG000812.                                                        #
# - M06 - Medvigy, D.M., 2006: The State of the Regional Carbon Cycle: results from a      #
#         constrained coupled ecosystem-atmosphere model, 2006.  Ph.D. dissertation,       #
#         Harvard University, Cambridge, MA, 322pp.                                        #
# - M01 - Moorcroft, P. R., G. C. Hurtt, S. W. Pacala, 2001: A method for scaling          #
#         vegetation dynamics: the ecosystem demography model, Ecological Monographs, 71,  #
#         557-586.                                                                         #
# - F96 - Foley, J. A., I. Colin Prentice, N. Ramankutty, S. Levis, D. Pollard, S. Sitch,  #
#         A. Haxeltime, 1996: An integrated biosphere model of land surface processes,     #
#         terrestrial carbon balance, and vegetation dynamics. Glob. Biogeochem. Cycles,   #
#         10, 603-602.                                                                     #
# - L95 - Leuning, R., F. M. Kelliher, D. G. G. de Pury, E. D. Schulze, 1995: Leaf         #
#         nitrogen, photosynthesis, conductance, and transpiration: scaling from leaves to #
#         canopies. Plant, Cell, and Environ., 18, 1183-1200.                              #
# - C91 - Collatz, G. J., J. T. Ball, C. Grivet, J. A. Berry, 1991: Physiological and      #
#         environmental regulation of stomatal conductance, photosynthesis and             #
#         transpiration: A model that includes a laminar boundary layer. Agric. For.       #
#         Meteor., 53, 107-136.                                                            #
# - B01 - Bernacchi, C. J., E. L. Singsaas, C. Pimentel, A. R. Portis Jr., S. P. Long,     #
#         2001: Improved temperature response functions for models of Rubisco-limited      #
#         photosynthesis. Plant, Cell and Environ., 24, 253-259.                           #
# - B02 - Bernacchi, C. J., A. R. Portis, H. Nakano, S. von Caemmerer, S. P. Long,         #
#         2002: Temperature response of mesophyll conductance.  Implications for the       #
#         determination of Rubisco enzyme kinetics and for limitations to photosynthesis   #
#         in vivo.  Plant physiol., 130, 1992-1998.                                        #
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Bounds for internal carbon and water stomatal conductance.                           #
#------------------------------------------------------------------------------------------#
c34smin.lint.co2 <<- 10.   * umol.2.mol # Minimum carbon dioxide concentration  [  mol/mol]
c34smax.lint.co2 <<- 1200. * umol.2.mol # Maximum carbon dioxide concentration  [  mol/mol]
c34smax.gsw      <<- 1.e+2              # Max. stomatal conductance (water)     [ mol/m2/s]
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Many parameters in the model are temperature-dependent and utilise a modified        #
# Arrhenius function to determine this dependence.  For that to work, reference values at  #
# a given temperature (tarrh, in Kelvin).                                                  #
#------------------------------------------------------------------------------------------#
tarrh        <<- 15.0+t00
tarrhi       <<- 1./tarrh
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     This is an alternative way to express the temperature dependence, which is used by   #
# C91.  (tcollatz, in Kelvin).  fcollatz is the factor that multiply the temperature       #
# departure from reference.                                                                #
#------------------------------------------------------------------------------------------#
tcollatz        <<- 15.0+t00
fcollatz        <<- 0.1
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     The next two variables are the parameters for the compensation point as in IBIS.     #
#------------------------------------------------------------------------------------------#
compp.ref.ibis  <<-  o2.ref / (2 * 4500.) # Gamma* reference                     [ mol/mol]
compp.hor.ibis  <<- 5000.                 # Gamma* "Activation energy"           [       K]
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     The next two variables are the parameters for the compensation point.  Notice that   #
# we give the compensation point rather than tau.                                          #
#------------------------------------------------------------------------------------------#
compp.ref.coll  <<-  o2.ref * 0.57 / (2 * 2600.) # Ref. value at 15C (not 25C)   [ mol/mol]
compp.base.coll <<-  1/0.57                      # Comp. point base parameter    [    ----]
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     The next three variables are the parameters for the compensation point using         #
# Bernacchi data.  Their equation is the Arrhenius formula, just expressed in a different  #
# way, so here we give the converted values.  They have 2 papers, in which the results are #
# quite different, so I will add both here.                                                #
#------------------------------------------------------------------------------------------#
compp.ref.b01 <<-  25.28064 * umol.2.mol # Reference compensation point at 25 C  [ mol/mol]
compp.hor.b01 <<-  4549.877              # Compensation point "c" parameter      [    ----]
compp.ref.b02 <<-  26.59113 * umol.2.mol # Reference compensation point at 25 C  [ mol/mol]
compp.hor.b02 <<-  2941.845              # Compensation point "c" parameter      [    ----]
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      The following terms are used to find the Michaelis-Menten coefficient for CO2.      #
#------------------------------------------------------------------------------------------#
kco2.ref.ibis  <<- 150. * umol.2.mol # Reference CO2 concentration               [ mol/mol]
kco2.hor.ibis  <<- 6000.             # Reference exponential coefficient         [       K]
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     The next two variables are the parameters for the Michaelis-Menten coefficients when #
# using Collatz equation.  Collatz equation uses partial pressure, so we must convert them #
# to mixing ratio.  Also, we use the reference at 15 C rather than 25 C, so Vm is the      #
# same.                                                                                    #
#------------------------------------------------------------------------------------------#
kco2.ref.coll  <<-  30. * mmco2 * mmdryi / prefsea / 2.1 # Ref. CO2 conc.        [ mol/mol]
kco2.base.coll <<-  2.1                                  # Power base            [    ----]
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     The next two variables are the parameters for the Michaelis-Menten coefficients when #
# using Bernacchi equations.                                                               #
#------------------------------------------------------------------------------------------#
kco2.ref.b01  <<-  133.38 * umol.2.mol            # Reference CO2 concentration  [ mol/mol]
kco2.hor.b01  <<-  9553.179                       # "Activation energy           [       K]
kco2.ref.b02  <<-  87.828 * umol.2.mol            # Reference CO2 concentration  [ mol/mol]
kco2.hor.b02  <<-  9740.803                       # "Activation energy           [       K]
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     These terms are used to find the Michaelis-Menten coefficient for O2.                #
#------------------------------------------------------------------------------------------#
ko2.ref.ibis  <<- 0.250     # Reference O2 concentration.                        [ mol/mol]
ko2.hor.ibis  <<-  1400.    # Reference exponential coefficient                  [       K]
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     These terms are used to find the Michaelis-Menten coefficient for O2.  Notice that   #
# the reference must be given at 15 C, not 25 C as usual.                                  #
#------------------------------------------------------------------------------------------#
ko2.ref.coll   <<- 3.e4 * mmo2 * mmdryi / prefsea / 1.2 # Ref. O2 concentration.  [ mol/mol]
ko2.base.coll  <<-  1.2                                 # Power base              [    ----]
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     These terms are used to find the Michaelis-Menten coefficient for O2 using the       #
# Bernacchi et al. (2002) method.                                                          #
#------------------------------------------------------------------------------------------#
ko2.ref.b01  <<-  0.1665438                       # Reference O2 concentration   [ mol/mol]
ko2.hor.b01  <<-  4376.483                        # "Activation energy           [       K]
ko2.ref.b02  <<-  0.1190386                       # Reference O2 concentration   [ mol/mol]
ko2.hor.b02  <<-  2852.844                        # "Activation energy           [       K]
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     These are constants obtained in Leuning et al. (1995) and Collatz et al. (1991) to   #
# convert different conductivities.                                                        #
#------------------------------------------------------------------------------------------#
gbh.2.gbw  <<- 1.075               # heat   to water  - leaf boundary layer
gsc.2.gsw  <<- 1.60                # carbon to water  - stomata
gbc.2.gbw  <<- gsc.2.gsw^twothirds # carbon to water  - leaf boundary layer
gsw.2.gsc  <<- 1.0 / gsc.2.gsw     # water  to carbon - stomata
gbw.2.gbc  <<- 1.0 / gbc.2.gbw     # water  to carbon - leaf boundary layer
gsc.2.gsw  <<- 1.0 / gsw.2.gsc     # carbon to water  - stomata
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     This is the minimum threshold for the photosynthetically active radiation, in        #
# umol/m2/s to consider non-night time conditions (day time or twilight).                  #
#------------------------------------------------------------------------------------------#
par.twilight.min <<- 0.5 * Watts.2.Ein # Minimum non-nocturnal PAR.              [mol/m2/s]
#------------------------------------------------------------------------------------------#

#----- This is a flag for the maximum representable number in R. --------------------------#
discard <<- 2^1023
#------------------------------------------------------------------------------------------#



#----- Fudging parameters to try to tune photosynthesis. ----------------------------------#
if (! "alpha.c3" %in% ls()){
   alpha.c3  <<-   0.08
}else{
   alpha.c3  <<-   alpha.c3
}#end if
if (! "alpha.c4" %in% ls()){
   alpha.c4  <<-   0.053
}else{
   alpha.c4  <<-   alpha.c4
}#end if
if (! "vmfact.c3" %in% ls()){
   vmfact.c3  <<-   1.0   
}else{
   vmfact.c3  <<-   vmfact.c3
}#end if
if (! "vmfact.c4" %in% ls()){
   vmfact.c4  <<-   1.0   
}else{
   vmfact.c4  <<-   vmfact.c4
}#end if
if (! "mphoto.c3" %in% ls()){
   mphoto.c3  <<-   8.0   
}else{
   mphoto.c3  <<-   mphoto.c3
}#end if
if (! "mphoto.aa" %in% ls()){
   mphoto.aa  <<-   6.4   
}else{
   mphoto.aa  <<-   mphoto.aa
}#end if
if (! "mphoto.c4" %in% ls()){
   mphoto.c4  <<-   4.0
}else{
   mphoto.c4  <<-   mphoto.c4
}#end if
if (! "gamma.c3" %in% ls()){
   gamma.c3  <<-   0.020
}else{
   gamma.c3  <<-   gamma.c3
}#end if
if (! "gamma.aa" %in% ls()){
   gamma.aa  <<-   0.028
}else{
   gamma.aa  <<-   gamma.aa
}#end if
if (! "gamma.c4" %in% ls()){
   gamma.c4  <<-   0.040
}else{
   gamma.c4  <<-   gamma.c4
}#end if
if (! "d0.grass" %in% ls()){
   d0.grass      <<-   0.01
}else{
   d0.grass      <<-   d0.grass
}#end if
if (! "d0.tree" %in% ls()){
   d0.tree      <<-   0.01
}else{
   d0.tree      <<-   d0.tree
}#end if


if (! "klowin" %in% ls()){
   klowin      <<-   18000.
}else{
   klowin    <<-   klowin
}#end if
if (! "base.c3" %in% ls()){
   base.c3       <<-   2.4
}else{
   base.c3       <<-   base.c3
}#end if
if (! "base.c4" %in% ls()){
   base.c4       <<-   2.0
}else{
   base.c4       <<-   base.c4
}#end if
if (! "b.c3" %in% ls()){
   b.c3       <<-   10000.
}else{
   b.c3       <<-   b.c3
}#end if
if (! "b.aa" %in% ls()){
   b.aa       <<-   1000.
}else{
   b.aa       <<-   b.aa
}#end if
if (! "b.c4" %in% ls()){
   b.c4       <<-    8000.
}else{
   b.c4       <<-   b.c4
}#end if
if (! "orient.tree" %in% ls()){
   orient.tree     <<-   0.00
}else{
   orient.tree     <<-   orient.tree
}#end if
if (! "orient.aa" %in% ls()){
   orient.aa       <<-   0.00
}else{
   orient.aa       <<-   orient.aa
}#end if
if (! "orient.grass" %in% ls()){
   orient.grass    <<-  0.00
}else{
   orient.grass    <<-   orient.grass
}#end if
if (! "clumping.tree" %in% ls()){
   clumping.tree     <<-   0.735
}else{
   clumping.tree     <<-   clumping.tree
}#end if
if (! "clumping.aa" %in% ls()){
   clumping.aa       <<-   0.735
}else{
   clumping.aa       <<-   clumping.aa
}#end if
if (! "clumping.grass" %in% ls()){
   clumping.grass    <<-  1.0
}else{
   clumping.grass    <<-   clumping.grass
}#end if
if (! "lwidth.grass" %in% ls()){
   lwidth.grass      <<-   0.05
}else{
   lwidth.grass      <<-   lwidth.grass
}#end if
if (! "lwidth.bltree" %in% ls()){
   lwidth.bltree     <<-   0.10
}else{
   lwidth.bltree     <<-   lwidth.bltree
}#end if
if (! "lwidth.nltree" %in% ls()){
   lwidth.nltree     <<-   0.05
}else{
   lwidth.nltree     <<-   lwidth.nltree
}#end if
if (! "vm.tcold.c3temp" %in% ls()){
   vm.tcold.c3temp   <<- 4.7137
}else{
   vm.tcold.c3temp   <<- vm.tcold.c3temp
}#end if
if (! "vm.tcold.c3trop" %in% ls()){
   vm.tcold.c3trop   <<- 8.0
}else{
   vm.tcold.c3trop   <<- vm.tcold.c3trop
}#end if
if (! "vm.tcold.aa" %in% ls()){
   vm.tcold.aa       <<- 4.7137
}else{
   vm.tcold.aa       <<- vm.tcold.aa
}#end if
if (! "vm.tcold.c4" %in% ls()){
   vm.tcold.c4       <<- 8.0
}else{
   vm.tcold.c4       <<- vm.tcold.aa
}#end if
if (! "vm.thot.c3temp" %in% ls()){
   vm.thot.c3temp       <<- 45.0
}else{
   vm.thot.c3temp       <<- vm.thot.c3temp
}#end if
if (! "vm.thot.c3trop" %in% ls()){
   vm.thot.c3trop       <<- 45.0
}else{
   vm.thot.c3trop       <<- vm.thot.c3trop
}#end if
if (! "vm.thot.aa" %in% ls()){
   vm.thot.aa           <<- 45.0
}else{
   vm.thot.aa           <<- vm.thot.aa
}#end if
if (! "vm.thot.c4" %in% ls()){
   vm.thot.c4           <<- 45.0
}else{
   vm.thot.c4           <<- vm.thot.c4
}#end if
if (! "vm.decay.ecold.c3" %in% ls()){
   vm.decay.ecold.c3    <<- 0.4
}else{
   vm.decay.ecold.c3    <<- vm.decay.ecold.c3
}#end if
if (! "vm.decay.ecold.c4" %in% ls()){
   vm.decay.ecold.c4    <<- 0.4
}else{
   vm.decay.ecold.c4    <<- vm.decay.ecold.c4
}#end if
if (! "vm.decay.ehot.c3" %in% ls()){
   vm.decay.ehot.c3    <<- 0.4
}else{
   vm.decay.ehot.c3    <<- vm.decay.ehot.c3
}#end if
if (! "vm.decay.ehot.c4" %in% ls()){
   vm.decay.ehot.c4    <<- 0.4
}else{
   vm.decay.ehot.c4    <<- vm.decay.ehot.c4
}#end if
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#    The following parameter is the k coefficient in Foley et al. (1996) that is used to   #
# determine the CO2-limited photosynthesis for C4 grasses.                                 #
#------------------------------------------------------------------------------------------#
klowco2      <<- klowin * mmco2 /mmdry  # Coefficient for low CO2      [ mol/mol]
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#    These parameters will be assigned to the PFT structure.  Notice that depending on the #
# physiology method (Foley or Collatz), some of them will be used, and some of them will   #
# not.                                                                                     #
#------------------------------------------------------------------------------------------#
vm.hor            <<- 3000.    # Ref. exp. coeff. for carboxylase fctn.          [       K]
vm.base.c3        <<- base.c3  # Base Vm for C3 - Collatz et al. (1991)          [mol/m2/s]
vm.base.c4        <<- base.c4  # Base Vm for C4 - Collatz et al. (1992)          [mol/m2/s]
vm.decay.a        <<- 220000.  # Decay function for warm temperatures -- A       [   J/mol]
vm.decay.b        <<- 695.     # Decay function for warm temperatures -- B       [ J/mol/K]

lr.hor            <<- 3000.      # Ref. exp. coeff. for carboxylase fctn.        [       K]
lr.base.c3        <<- base.c3
lr.base.c4        <<- base.c4

lr.tcold.c3temp   <<- vm.tcold.c3temp
lr.tcold.c3trop   <<- vm.tcold.c3trop
lr.tcold.aa       <<- vm.tcold.aa
lr.tcold.c4       <<- vm.tcold.c4

lr.thot.c3temp    <<- vm.thot.c3temp
lr.thot.c3trop    <<- vm.thot.c3trop
lr.thot.aa        <<- vm.thot.aa
lr.thot.c4        <<- vm.thot.c4

lr.decay.ecold.c3 <<- vm.decay.ecold.c3
lr.decay.ecold.c4 <<- vm.decay.ecold.c4
lr.decay.ehot.c3  <<- vm.decay.ehot.c3 
lr.decay.ehot.c4  <<- vm.decay.ehot.c4 

#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#    Define some common values to make it easier to change the list below...               #
#------------------------------------------------------------------------------------------#
C2B    <<- 2.0
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#    These variables define the allometry to use                                           #
# 0 - ED-2.1                                                                               #
# 1 - a.  AGB is based on Baker et al. (2004), Chave formula, keeping the same Bleaf as    #
#         ED-2.1 and Bdead to make up the difference.                                      #
#     b.  Use the crown area as in Poorter et al. (2006) for tropical PFTs                 #
#     c.  Use a simple rooting depth ranging from 0.5 to 5.0                               #
# 2 - Same as 1, but using the height as in Poorter et al. (2006) for tropical PFTs, and   #
#     finding a fit similar to Baker et al. (2004) that doesn't depend on height, using    #
#     Saldarriaga et al. (1988) as a starting point.                                       #
# 3 - New set of traits and parameters (mostly) based on literature and available data.    #
#------------------------------------------------------------------------------------------#
if ("iallom" %in% ls()){
   iallom <<- iallom
}else{
   iallom <<- 3
}#end if
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Factors for leaf:sapwood biomass ratio.  The original ED-1 number is incorrect, and  #
# we keep it incorrect unless the PFT is tropical and allometry is set to 4, in which case #
# we combine the pipe model with the data from Calvo-Alvarado et al. (2008) (and shape     #
# parameter from Falster et al. 2016) to derive the ratio.                                 #
#                                                                                          #
# References:                                                                              #
#                                                                                          #
# Calvo-Alvarado, J. C., N. G. McDowell, and R. H. Waring. Allometric relationships        #
#    predicting foliar biomass and leaf area:sapwood area ratio from tree height in five   #
#    Costa Rican rain forest species. Tree Physiol., 28 (11):1601-1608, Sep 2008.          #
#    doi:10.1093/treephys/28.11.1601. (CA08)                                               #
#                                                                                          #
# Falster, D. S., R. G. FitzJohn, A. Brannstrom, U. Dieckmann, and M. Westoby.  plant: A   #
#    package for modelling forest trait ecology and evolution.  Methods Ecol. Evol., 7(2): #
#    136-146, Feb 2016. doi:10.1111/2041-210X.12525. (F16)                                 #
#                                                                                          #
# McDowell, N., H. Barnard, B. Bond, T. Hinckley, R. Hubbard, H. Ishii, B. Kostner,        #
#    F. Magnani, J. Marshall, F. Meinzer, N. Phillips, M. Ryan, and D. Whitehead. The      #
#    relationship between tree height and leaf area: sapwood area ratio. Oecologia,        #
#    132(1):12-20, Jun 2002. doi:10.1007/s00442-002-0904-x. (MD02)                         #
#                                                                                          #
# Rosell, J. A., S. Gleason, R. Mendez-Alonzo, Y. Chang, and M. Westoby. Bark functional   #
#    ecology: evidence for tradeoffs, functional coordination, and environment producing   #
#    bark diversity. New Phytol., 201(2): 486-497, Jan 2014. doi:10.1111/nph.12541. (R14)  #
#                                                                                          #
# Yokozawa M., and T. Hara. Foliage profile, size structure and stem diameter-plant height #
#    relationship in crowded plant populations. Ann. Bot.-London, 76(3):271-285, Sep 1995. #
#    doi:10.1006/anbo.1995.1096. (YH95)                                                    #
#                                                                                          #
#      Leaf-to-sapwood area ratio (Al:As, or As/Al) for conifers was estimated from the    #
# average of all conifers listed in Table 1 of MD02, weighted by number of individuals.    #
# Broadleaf is currently tropical-only, and was obtained from the average of all points    #
# from Fig. 1 of CA08 (obtained from extracting data from the figure itself).              #
#                                                                                          #
#       Eta is a parameter that described the distribution of leaf area within the crown.  #
# According to F16 and the original referennce (YH95), eta=1 is typical of conifers,       #
# whereas eta=12 is closer to the profiles observed in angiosperms.  Araucarias somewhat   #
# intermediate, so we set eta=5.                                                           #
#------------------------------------------------------------------------------------------#
eta.f16            <<- c(bl=12.,cf=1.0,aa=5.0)
eta.c.f16          <<- 1. - 2.0 / (1.0+eta.f16) + 1.0 / (1 + 2.0 * eta.f16)
asal.bar           <<- c(bl=7.400139e-05,cf=3.709184e-05,aa=3.709184e-05)
sapwood.ratio.orig <<- 3900.
sapwood.factor     <<- 1.0 / (eta.c.f16 * asal.bar * 1000. / C2B)
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     These constants will help defining the allometric parameters for IALLOM 1 and 2.     #
#------------------------------------------------------------------------------------------#
ndead.small = c( -1.26395300, 2.43236100,  1.80180100)
ndead.large = c( -0.83468050, 2.42557360,  2.68228050)
nleaf       = c(  0.01925119, 0.97494935,  2.58585087)
uleaf       = c( -1.09254800, 1.28505100,  3.199019  )
ncrown.area = c(  0.11842950, 1.05211970)
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Additional parameters for determining specific heat.                                 #
#                                                                                          #
# spht.tref - Reference temperature (in Kelvin) for specific heat properties.              #
# wdr.fs    - Moisture content at fiber saturation. Excess moisture is treated as free     #
#             water (FPL10).                                                               #
#------------------------------------------------------------------------------------------#
spht.tref  <<- t00 + 15.
wdr.fs     <<- 0.30
#------------------------------------------------------------------------------------------#


#----- Define reference height and coefficients for tropical allometry. -------------------#
if (iallom %in% c(0,1)){
   hgt.ref.trop = NA
   b1Ht.trop    = 0.37 * log(10)
   b2Ht.trop    = 0.64
   hgt.max.trop = 35.0
}else if (iallom %in% c(2)){
   #---------------------------------------------------------------------------------------#
   #     Use the allometry proposed by:                                                    #
   #                                                                                       #
   # Poorter, L., L. Bongers, F. Bongers, 2006: Architecture of 54 moist-forest tree       #
   #    species: traits, trade-offs, and functional groups.  Ecology, 87, 1289-1301.       #
   #                                                                                       #
   #---------------------------------------------------------------------------------------#
   hgt.ref.trop = 61.7
   b1Ht.trop    = 0.0352
   b2Ht.trop    = 0.694
   hgt.max.trop = 35.0
   #---------------------------------------------------------------------------------------#
}else if (iallom %in% c(3)){
   #---------------------------------------------------------------------------------------#
   #     Allometric equation based on the Sustainable Landscapes data.                     #
   #                                                                                       #
   #    Longo, M. et al. Carbon Debt and Recovery time of degraded forests in              #
   #       the Amazon., in prep.                                                           #
   #                                                                                       #
   #    Equation was derived from multiple forest inventories carried out at multiple      #
   # locations in the Brazilian Amazon, and fitted using a heteroscedastic least           #
   # squares approach (though results converged to a homoscedastic fit).  This equation    #
   # is very similar to Feldpausch et al (2012) equation for South America.                #
   #                                                                                       #
   # Total number of trees: 17010                                                          #
   # hgt_ref = 47.2    (95% CI: [  44.8;   48.8])                                          #
   # b1Ht    = 0.0440  (95% CI: [0.0427; 0.0454])                                          #
   # b2Ht    = 0.802   (95% CI: [ 0.788;  0.822])                                          #
   # R2      = 0.677                                                                       #
   # RMSE    = 5.4                                                                         #
   #---------------------------------------------------------------------------------------#
   hgt.ref.trop = 47.173
   b1Ht.trop    = 0.044037
   b2Ht.trop    = 0.80248
   hgt.max.trop = 42.0
   #---------------------------------------------------------------------------------------#
}#end if
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#   Coefficients for DBH -> Bleaf allometry (iallom = 3).  We use the allometric           #
# parameters from X16 for evergreens because the biomass is similar to L83 for large trees #
# (based on > 1400 trees with DBH>1cm in French Guiana) and the small size is similar to   #
# values from K09 (based on ~100 trees with DBH < 32cm in Malaysia).                       #
#                                                                                          #
#  References:                                                                             #
#                                                                                          #
#   Lescure, J.-P., H. Puig, B. Riera, D. Leclerc, A. Beekman, and A. Beneteau. La         #
#      phytomasse epigee d'une foret dense en Guyane Francaise.  Acta Ecol.-Oec. Gen.,     #
#      4(3), 237--251, 1983. http://www.documentation.ird.fr/hor/fdi:010005089 (L83)       #
#                                                                                          #
#   Kenzo, T., T. Ichie, D. Hattori, T. Itioka, C. Handa, T. Ohkubo, J. J. Kendawang,      #
#      M. Nakamura, M. Sak- aguchi, N. Takahashi, M. Okamoto, A. Tanaka-Oda, K. Sakurai,   #
#      and I. Ninomiya. Development of al- lometric relationships for accurate estimation  #
#      of above- and below-ground biomass in tropical secondary forests in Sarawak,        #
#      Malaysia. J. Trop. Ecol., 25(4):371-386, Jul 2009.                                  #
#      doi:10.1017/S0266467409006129. (K09)                                                #
#                                                                                          #
#   Xu, X., D. Medvigy, J. S. Powers, J. M. Becknell, and K. Guan. Diversity in plant      #
#      hydraulic traits explains seasonal and inter-annual variations of vegetation        #
#      dynamics in seasonally dry tropical forests. New Phytol., 212(1):80-95, Oct 2016.   #
#      doi:10.1111/nph.14009. (X16).                                                       #
#------------------------------------------------------------------------------------------#
l83.l1 = 0.00873
l83.l2 = 2.1360
x16.l1 = c(0.006,0.024,0.016,0.046)[4]
x16.l2 = c( 2.04, 1.75, 2.12, 1.93)[4]
k09.l1 = 0.0180
k09.l2 = 1.83
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     The following variables will define the PFT characteristics regarding the            #
# photosynthesis.                                                                          #
#------------------------------------------------------------------------------------------#
pft01 = list( name               = "C4 grass"
            , key                = "C4G"
            , colour             = "#DDCC77"
            , tropical           = TRUE
            , savannah           = FALSE
            , conifer            = FALSE
            , grass              = TRUE
            , liana              = FALSE
            , pathway            = 4
            , d0                 = d0.grass
            , vm.hor             = vm.hor
            , vm.base            = vm.base.c4
            , vm.decay.a         = vm.decay.a
            , vm.decay.b         = vm.decay.b
            , vm.low.temp        = vm.tcold.c4   + t00
            , vm.high.temp       = vm.thot.c4    + t00
            , vm.decay.e.low     = vm.decay.ecold.c4
            , vm.decay.e.high    = vm.decay.ehot.c4
            , lr.hor             = lr.hor
            , lr.base            = lr.base.c4
            , lr.low.temp        = lr.tcold.c4   + t00
            , lr.high.temp       = lr.thot.c4    + t00
            , lr.decay.e.low     = lr.decay.ecold.c4
            , lr.decay.e.high    = lr.decay.ehot.c4
            , vm0                = 12.500 * vmfact.c4 * umol.2.mol
            , m                  = mphoto.c4
            , alpha              = alpha.c4
            , b                  = b.c4 * umol.2.mol
            , gamma.resp         = gamma.c4
            , effarea.transp     = 1.0
            , rho                = 0.20
            , leaf.turnover.rate = ifelse(iallom %in% 3, 3.0, 2.0)
            , root.turnover.rate = NA
            , bark.turnover.rate = 0.0
            , SLA                = 22.7
            , hgt.ref            = hgt.ref.trop
            , b1Ht               = b1Ht.trop
            , b2Ht               = b2Ht.trop
            , b1Bl.small         = NA
            , b2Bl.small         = NA
            , b1Bl.large         = NA
            , b2Bl.large         = NA
            , b1Bs.small         = NA
            , b2Bs.small         = NA
            , b1Bs.large         = NA
            , b2Bs.large         = NA
            , b1Ca               = 2.490154
            , b2Ca               = 0.8068806
            , b1Cl               = 0.99
            , b2Cl               = 1.00
            , b1Mh               = 0.495
            , b1Vol              = NA
            , b2Vol              = NA
            , b1Xb               = NA
            , b1Xs               = NA
            , hgt.min            = NA
            , hgt.max            = NA
            , qroot              = 1.0
            , qsw                = NA
            , qbark              = NA
            , qwai               = NA
            , qrhob              = NA
            , agf.bs             = 0.7
            , orient.factor      = orient.grass
            , clumping.factor    = clumping.grass
            , leaf.width         = lwidth.grass
            , init.density       = 0.1
            , veg.hcap.min       = 7.30807E+00
            )

pft02 = list( name               = "Early tropical"
            , key                = "ETR"
            , colour             = "#83CCC0"
            , tropical           = TRUE
            , savannah           = FALSE
            , conifer            = FALSE
            , grass              = FALSE
            , liana              = FALSE
            , pathway            = 3
            , d0                 = d0.tree
            , vm.hor             = vm.hor
            , vm.base            = vm.base.c3
            , vm.decay.a         = vm.decay.a
            , vm.decay.b         = vm.decay.b
            , vm.low.temp        = vm.tcold.c3trop + t00
            , vm.high.temp       = vm.thot.c3trop  + t00
            , vm.decay.e.low     = vm.decay.ecold.c3
            , vm.decay.e.high    = vm.decay.ehot.c3
            , lr.hor             = lr.hor
            , lr.base            = lr.base.c3
            , lr.low.temp        = lr.tcold.c3trop + t00
            , lr.high.temp       = lr.thot.c3trop  + t00
            , lr.decay.e.low     = lr.decay.ecold.c3
            , lr.decay.e.high    = lr.decay.ehot.c3
            , vm0                = 18.750 * vmfact.c3 * umol.2.mol
            , m                  = mphoto.c3
            , alpha              = alpha.c3
            , b                  = b.c3 * umol.2.mol
            , gamma.resp         = gamma.c3
            , effarea.transp     = 1.0
            , rho                = 0.53
            , leaf.turnover.rate = ifelse(iallom %in% 3, 1.25, 1.00)
            , root.turnover.rate = NA
            , bark.turnover.rate = NA
            , SLA                = NA
            , hgt.ref            = hgt.ref.trop
            , b1Ht               = b1Ht.trop
            , b2Ht               = b2Ht.trop
            , b1Bl.small         = NA
            , b2Bl.small         = NA
            , b1Bl.large         = NA
            , b2Bl.large         = NA
            , b1Bs.small         = NA
            , b2Bs.small         = NA
            , b1Bs.large         = NA
            , b2Bs.large         = NA
            , b1Ca               = 2.490154
            , b2Ca               = 0.8068806
            , b1Cl               = 0.3106775
            , b2Cl               = 1.098
            , b1Mh               = 0.8370557
            , b1Vol              = NA
            , b2Vol              = NA
            , b1Xb               = NA
            , b1Xs               = NA
            , hgt.min            = NA
            , hgt.max            = NA
            , qroot              = 1.0
            , qsw                = NA
            , qbark              = NA
            , qwai               = NA
            , qrhob              = NA
            , agf.bs             = 0.7
            , orient.factor      = orient.tree
            , clumping.factor    = clumping.tree
            , leaf.width         = lwidth.bltree
            , init.density       = 0.1
            , veg.hcap.min       = 9.53226E+00
            )

pft03 = list( name               = "Mid tropical"
            , key                = "MTR"
            , colour             = "#44AA99"
            , tropical           = TRUE
            , savannah           = FALSE
            , conifer            = FALSE
            , grass              = FALSE
            , liana              = FALSE
            , pathway            = 3
            , d0                 = d0.tree
            , vm.hor             = vm.hor
            , vm.base            = vm.base.c3
            , vm.decay.a         = vm.decay.a
            , vm.decay.b         = vm.decay.b
            , vm.low.temp        = vm.tcold.c3trop + t00
            , vm.high.temp       = vm.thot.c3trop  + t00
            , vm.decay.e.low     = vm.decay.ecold.c3
            , vm.decay.e.high    = vm.decay.ehot.c3
            , lr.hor             = lr.hor
            , lr.base            = lr.base.c3
            , lr.low.temp        = lr.tcold.c3trop + t00
            , lr.high.temp       = lr.thot.c3trop  + t00
            , lr.decay.e.low     = lr.decay.ecold.c3
            , lr.decay.e.high    = lr.decay.ehot.c3
            , vm0                = 12.500 * vmfact.c3 * umol.2.mol
            , m                  = mphoto.c3
            , alpha              = alpha.c3
            , b                  = b.c3 * umol.2.mol
            , gamma.resp         = gamma.c3
            , effarea.transp     = 1.0
            , rho                = 0.71
            , leaf.turnover.rate = ifelse(iallom %in% 3, 0.60, 0.50)
            , root.turnover.rate = NA
            , bark.turnover.rate = NA
            , SLA                = 16.0
            , hgt.ref            = hgt.ref.trop
            , b1Ht               = b1Ht.trop
            , b2Ht               = b2Ht.trop
            , b1Bl.small         = NA
            , b2Bl.small         = NA
            , b1Bl.large         = NA
            , b2Bl.large         = NA
            , b1Bs.small         = NA
            , b2Bs.small         = NA
            , b1Bs.large         = NA
            , b2Bs.large         = NA
            , b1Ca               = 2.490154
            , b2Ca               = 0.8068806
            , b1Cl               = 0.3106775
            , b2Cl               = 1.098
            , b1Mh               = 0.8370557
            , b1Vol              = NA
            , b2Vol              = NA
            , b1Xb               = NA
            , b1Xs               = NA
            , hgt.min            = NA
            , hgt.max            = NA
            , qroot              = 1.0
            , qsw                = NA
            , qbark              = NA
            , qwai               = NA
            , qrhob              = NA
            , agf.bs             = 0.7
            , orient.factor      = orient.tree
            , clumping.factor    = clumping.tree
            , leaf.width         = lwidth.bltree
            , init.density       = 0.1
            , veg.hcap.min       = 1.46161E+01
            )

pft04 = list( name               = "Late tropical"
            , key                = "LTR"
            , colour             = "#186659"
            , tropical           = TRUE
            , savannah           = FALSE
            , conifer            = FALSE
            , grass              = FALSE
            , liana              = FALSE
            , pathway            = 3
            , d0                 = d0.tree
            , vm.hor             = vm.hor
            , vm.base            = vm.base.c3
            , vm.decay.a         = vm.decay.a
            , vm.decay.b         = vm.decay.b
            , vm.low.temp        = vm.tcold.c3trop + t00
            , vm.high.temp       = vm.thot.c3trop  + t00
            , vm.decay.e.low     = vm.decay.ecold.c3
            , vm.decay.e.high    = vm.decay.ehot.c3
            , lr.hor             = lr.hor
            , lr.base            = lr.base.c3
            , lr.low.temp        = lr.tcold.c3trop + t00
            , lr.high.temp       = lr.thot.c3trop  + t00
            , lr.decay.e.low     = lr.decay.ecold.c3
            , lr.decay.e.high    = lr.decay.ehot.c3
            , vm0                = 6.250 * vmfact.c3 * umol.2.mol
            , m                  = mphoto.c3
            , alpha              = alpha.c3
            , b                  = b.c3 * umol.2.mol
            , gamma.resp         = gamma.c3
            , effarea.transp     = 1.0
            , rho                = 0.90
            , leaf.turnover.rate = ifelse(iallom %in% 3, 0.25, 1./3.)
            , root.turnover.rate = NA
            , bark.turnover.rate = NA
            , SLA                = NA
            , hgt.ref            = hgt.ref.trop
            , b1Ht               = b1Ht.trop
            , b2Ht               = b2Ht.trop
            , b1Bl.small         = NA
            , b2Bl.small         = NA
            , b1Bl.large         = NA
            , b2Bl.large         = NA
            , b1Bs.small         = NA
            , b2Bs.small         = NA
            , b1Bs.large         = NA
            , b2Bs.large         = NA
            , b1Ca               = 2.490154
            , b2Ca               = 0.8068806
            , b1Cl               = 0.3106775
            , b2Cl               = 1.098
            , b1Mh               = 0.8370557
            , b1Vol              = NA
            , b2Vol              = NA
            , b1Xb               = NA
            , b1Xs               = NA
            , hgt.min            = NA
            , hgt.max            = NA
            , qroot              = 1.0
            , qsw                = NA
            , qbark              = NA
            , qwai               = NA
            , qrhob              = NA
            , agf.bs             = 0.7
            , orient.factor      = orient.tree
            , clumping.factor    = clumping.tree
            , leaf.width         = lwidth.bltree
            , init.density       = 0.1
            , veg.hcap.min       = 2.43602E+01
            )

pft05 = list( name               = "Temperate C3 Grass"
            , key                = "TTG"
            , colour             = "#EBE0AA"
            , tropical           = FALSE
            , savannah           = FALSE
            , conifer            = FALSE
            , grass              = TRUE
            , liana              = FALSE
            , pathway            = 3
            , d0                 = d0.tree
            , vm.hor             = vm.hor
            , vm.base            = vm.base.c3
            , vm.decay.a         = vm.decay.a
            , vm.decay.b         = vm.decay.b
            , vm.low.temp        = vm.tcold.c3temp + t00
            , vm.high.temp       = vm.thot.c3temp  + t00
            , vm.decay.e.low     = vm.decay.ecold.c3
            , vm.decay.e.high    = vm.decay.ehot.c3
            , lr.hor             = lr.hor
            , lr.base            = lr.base.c3
            , lr.low.temp        = lr.tcold.c3temp + t00
            , lr.high.temp       = lr.thot.c3temp  + t00
            , lr.decay.e.low     = lr.decay.ecold.c3
            , lr.decay.e.high    = lr.decay.ehot.c3
            , vm0                = 18.3   * umol.2.mol
            , m                  = mphoto.c3
            , alpha              = alpha.c3
            , b                  = b.c3 * umol.2.mol
            , gamma.resp         = gamma.c3
            , effarea.transp     = 1.0
            , rho                = 0.32
            , leaf.turnover.rate = 2.0
            , root.turnover.rate = NA
            , bark.turnover.rate = 0.0
            , SLA                = 22.0
            , hgt.ref            = 0.0
            , b1Ht               = 0.4778
            , b2Ht               = -0.750
            , b1Bl.small         = 0.08
            , b2Bl.small         = 1.00
            , b1Bl.large         = 0.08
            , b2Bl.large         = 1.00
            , b1Bs.small         = 1.e-5
            , b2Bs.small         = 1.0
            , b1Bs.large         = 1.e-5
            , b2Bs.large         = 1.0
            , b1Ca               = 2.490154
            , b2Ca               = 0.8068806
            , b1Cl               = 0.99
            , b2Cl               = 1.0
            , b1Mh               = 0.495
            , b1Vol              = NA
            , b2Vol              = NA
            , b1Xb               = NA
            , b1Xs               = NA
            , hgt.min            = NA
            , hgt.max            = NA
            , qroot              = 1.0
            , qsw                = NA
            , qbark              = NA
            , qwai               = NA
            , qrhob              = NA
            , agf.bs             = 0.7
            , orient.factor      = -0.30
            , clumping.factor    =  1.00
            , leaf.width         = lwidth.grass
            , init.density       = 0.1
            , veg.hcap.min       = 9.16551E+00
            )

pft06 = list( name               = "North Pine"
            , key                = "NPN"
            , colour             = "#88CCEE"
            , tropical           = FALSE
            , savannah           = FALSE
            , conifer            = TRUE
            , grass              = FALSE
            , liana              = FALSE
            , pathway            = 3
            , d0                 = d0.tree
            , vm.hor             = vm.hor
            , vm.base            = vm.base.c3
            , vm.decay.a         = vm.decay.a
            , vm.decay.b         = vm.decay.b
            , vm.low.temp        = vm.tcold.c3temp + t00
            , vm.high.temp       = vm.thot.c3temp  + t00
            , vm.decay.e.low     = vm.decay.ecold.c3
            , vm.decay.e.high    = vm.decay.ehot.c3
            , lr.hor             = lr.hor
            , lr.base            = lr.base.c3
            , lr.low.temp        = lr.tcold.c3temp + t00
            , lr.high.temp       = lr.thot.c3temp  + t00
            , lr.decay.e.low     = lr.decay.ecold.c3
            , lr.decay.e.high    = lr.decay.ehot.c3
            , vm0                = 11.350 * umol.2.mol
            , m                  = mphoto.c3 * 6.3949 / 8.0
            , alpha              = alpha.c3
            , b                  = 1000. * umol.2.mol
            , gamma.resp         = gamma.c3
            , effarea.transp     = 2.0
            , rho                = NA
            , leaf.turnover.rate = 1./3.
            , root.turnover.rate = 3.927218
            , bark.turnover.rate = 0.0
            , SLA                = 6.0
            , hgt.ref            = 1.3
            , b1Ht               = 27.14
            , b2Ht               = -0.03884
            , b1Bl.small         = 0.024
            , b2Bl.small         = 1.899
            , b1Bl.large         = 0.024
            , b2Bl.large         = 1.899
            , b1Bs.small         = 0.147
            , b2Bs.small         = 2.238
            , b1Bs.large         = 0.147
            , b2Bs.large         = 2.238
            , b1Ca               = 2.490154
            , b2Ca               = 0.8068806
            , b1Cl               = 0.3106775
            , b2Cl               = 1.098
            , b1Mh               = 0.8370557
            , b1Vol              = NA
            , b2Vol              = NA
            , b1Xb               = NA
            , b1Xs               = NA
            , hgt.min            = NA
            , hgt.max            = NA
            , qroot              = 0.3463
            , qsw                = NA
            , qbark              = NA
            , qwai               = NA
            , qrhob              = NA
            , agf.bs             = 0.7
            , orient.factor      = 0.01
            , clumping.factor    = 0.735
            , leaf.width         = lwidth.nltree
            , init.density       = 0.1
            , veg.hcap.min       = 2.34683E-01
            )

pft07 = list( name               = "South Pine"
            , key                = "SPN"
            , colour             = "#B6E0F5"
            , tropical           = FALSE
            , savannah           = FALSE
            , conifer            = TRUE
            , grass              = FALSE
            , liana              = FALSE
            , pathway            = 3
            , d0                 = d0.tree
            , vm.hor             = vm.hor
            , vm.base            = vm.base.c3
            , vm.decay.a         = vm.decay.a
            , vm.decay.b         = vm.decay.b
            , vm.low.temp        = vm.tcold.c3temp + t00
            , vm.high.temp       = vm.thot.c3temp  + t00
            , vm.decay.e.low     = vm.decay.ecold.c3
            , vm.decay.e.high    = vm.decay.ehot.c3
            , lr.hor             = lr.hor
            , lr.base            = lr.base.c3
            , lr.low.temp        = lr.tcold.c3temp + t00
            , lr.high.temp       = lr.thot.c3temp  + t00
            , lr.decay.e.low     = lr.decay.ecold.c3
            , lr.decay.e.high    = lr.decay.ehot.c3
            , vm0                = 11.350 * umol.2.mol
            , m                  = mphoto.c3 * 6.3949 / 8.0
            , alpha              = alpha.c3
            , b                  = 1000. * umol.2.mol
            , gamma.resp         = gamma.c3
            , effarea.transp     = 2.0
            , rho                = NA
            , leaf.turnover.rate = 1./3.
            , root.turnover.rate = 4.117847
            , bark.turnover.rate = 0.4
            , SLA                = 9.0
            , hgt.ref            = 1.3
            , b1Ht               = 27.14
            , b2Ht               = -0.03884
            , b1Bl.small         = 0.024
            , b2Bl.small         = 1.899
            , b1Bl.large         = 0.024
            , b2Bl.large         = 1.899
            , b1Bs.small         = 0.147
            , b2Bs.small         = 2.238
            , b1Bs.large         = 0.147
            , b2Bs.large         = 2.238
            , b1Ca               = 2.490154
            , b2Ca               = 0.8068806
            , b1Cl               = 0.3106775
            , b2Cl               = 1.098
            , b1Mh               = 0.8370557
            , b1Vol              = NA
            , b2Vol              = NA
            , b1Xb               = NA
            , b1Xs               = NA
            , hgt.min            = NA
            , hgt.max            = NA
            , qroot              = 0.3463
            , qsw                = NA
            , qbark              = NA
            , qwai               = NA
            , qrhob              = NA
            , agf.bs             = 0.7
            , orient.factor      = 0.01
            , clumping.factor    = 0.735
            , leaf.width         = lwidth.nltree
            , init.density       = 0.1
            , veg.hcap.min       = 2.34683E-01
            )

pft08 = list( name               = "Late conifer"
            , key                = "LCN"
            , colour             = "#31708F"
            , tropical           = FALSE
            , savannah           = FALSE
            , conifer            = TRUE
            , grass              = FALSE
            , liana              = FALSE
            , pathway            = 3
            , d0                 = d0.tree
            , vm.hor             = vm.hor
            , vm.base            = vm.base.c3
            , vm.decay.a         = vm.decay.a
            , vm.decay.b         = vm.decay.b
            , vm.low.temp        = vm.tcold.c3temp + t00
            , vm.high.temp       = vm.thot.c3temp  + t00
            , vm.decay.e.low     = vm.decay.ecold.c3
            , vm.decay.e.high    = vm.decay.ehot.c3
            , lr.hor             = lr.hor
            , lr.base            = lr.base.c3
            , lr.low.temp        = lr.tcold.c3temp + t00
            , lr.high.temp       = lr.thot.c3temp  + t00
            , lr.decay.e.low     = lr.decay.ecold.c3
            , lr.decay.e.high    = lr.decay.ehot.c3
            , vm0                = 4.540 * umol.2.mol
            , m                  = mphoto.c3 * 6.3949 / 8.0
            , alpha              = alpha.c3
            , b                  = 1000. * umol.2.mol
            , gamma.resp         = gamma.c3
            , effarea.transp     = 2.0
            , rho                = NA
            , leaf.turnover.rate = 1./3.
            , root.turnover.rate = 3.800132
            , bark.turnover.rate = 0.0
            , SLA                = 10.0
            , hgt.ref            = 1.3
            , b1Ht               = 22.79
            , b2Ht               = -0.04445
            , b1Bl.small         = 0.0454
            , b2Bl.small         = 1.6829
            , b1Bl.large         = 0.0454
            , b2Bl.large         = 1.6829
            , b1Bs.small         = 0.1617
            , b2Bs.small         = 2.1536
            , b1Bs.large         = 0.1617
            , b2Bs.large         = 2.1536
            , b1Ca               = 2.490154
            , b2Ca               = 0.8068806
            , b1Cl               = 0.3106775
            , b2Cl               = 1.098
            , b1Mh               = 0.8370557
            , b1Vol              = NA
            , b2Vol              = NA
            , b1Xb               = NA
            , b1Xs               = NA
            , hgt.min            = NA
            , hgt.max            = NA
            , qroot              = 0.3463
            , qsw                = NA
            , qbark              = NA
            , qwai               = NA
            , qrhob              = NA
            , agf.bs             = 0.7
            , orient.factor      = 0.01
            , clumping.factor    = 0.735
            , leaf.width         = lwidth.nltree
            , init.density       = 0.1
            , veg.hcap.min       = 6.80074E-01
            )

pft09 = list( name               = "Early hardwood"
            , key                = "EHW"
            , colour             = "#CC83C0"
            , tropical           = FALSE
            , savannah           = FALSE
            , conifer            = FALSE
            , grass              = FALSE
            , liana              = FALSE
            , pathway            = 3
            , d0                 = d0.tree
            , vm.hor             = vm.hor
            , vm.base            = vm.base.c3
            , vm.decay.a         = vm.decay.a
            , vm.decay.b         = vm.decay.b
            , vm.low.temp        = vm.tcold.c3temp + t00
            , vm.high.temp       = vm.thot.c3temp  + t00
            , vm.decay.e.low     = vm.decay.ecold.c3
            , vm.decay.e.high    = vm.decay.ehot.c3
            , lr.hor             = lr.hor
            , lr.base            = lr.base.c3
            , lr.low.temp        = lr.tcold.c3temp + t00
            , lr.high.temp       = lr.thot.c3temp  + t00
            , lr.decay.e.low     = lr.decay.ecold.c3
            , lr.decay.e.high    = lr.decay.ehot.c3
            , vm0                = 20.387 * umol.2.mol
            , m                  = mphoto.c3 * 6.3949 / 8.0
            , alpha              = alpha.c3
            , b                  = 20000. * umol.2.mol
            , gamma.resp         = gamma.c3
            , effarea.transp     = 1.0
            , rho                = NA
            , leaf.turnover.rate = 0.0
            , root.turnover.rate = 5.772506
            , bark.turnover.rate = 0.0
            , SLA                = 30.0
            , hgt.ref            = 1.3
            , b1Ht               = 22.6799
            , b2Ht               = -0.06534
            , b1Bl.small         = 0.0129
            , b2Bl.small         = 1.7477
            , b1Bl.large         = 0.0129
            , b2Bl.large         = 1.7477
            , b1Bs.small         = 0.02648
            , b2Bs.small         = 2.95954
            , b1Bs.large         = 0.02648
            , b2Bs.large         = 2.95954
            , b1Ca               = 2.490154
            , b2Ca               = 0.8068806
            , b1Cl               = 0.3106775
            , b2Cl               = 1.098
            , b1Mh               = 0.8370557
            , b1Vol              = NA
            , b2Vol              = NA
            , b1Xb               = NA
            , b1Xs               = NA
            , hgt.min            = NA
            , hgt.max            = NA
            , qroot              = 1.1274
            , qsw                = NA
            , qbark              = NA
            , qwai               = NA
            , qrhob              = NA
            , agf.bs             = 0.7
            , orient.factor      = 0.25
            , clumping.factor    = 0.84
            , leaf.width         = lwidth.bltree
            , init.density       = 0.1
            , veg.hcap.min       = 8.95049E-02
            )

pft10 = list( name               = "Mid hardwood"
            , key                = "MHW"
            , colour             = "#AA4499"
            , tropical           = FALSE
            , savannah           = FALSE
            , conifer            = FALSE
            , grass              = FALSE
            , liana              = FALSE
            , pathway            = 3
            , d0                 = d0.tree
            , vm.hor             = vm.hor
            , vm.base            = vm.base.c3
            , vm.decay.a         = vm.decay.a
            , vm.decay.b         = vm.decay.b
            , vm.low.temp        = vm.tcold.c3temp + t00
            , vm.high.temp       = vm.thot.c3temp  + t00
            , vm.decay.e.low     = vm.decay.ecold.c3
            , vm.decay.e.high    = vm.decay.ehot.c3
            , lr.hor             = lr.hor
            , lr.base            = lr.base.c3
            , lr.low.temp        = lr.tcold.c3temp + t00
            , lr.high.temp       = lr.thot.c3temp  + t00
            , lr.decay.e.low     = lr.decay.ecold.c3
            , lr.decay.e.high    = lr.decay.ehot.c3
            , vm0                = 17.455 * umol.2.mol
            , m                  = mphoto.c3 * 6.3949 / 8.0
            , alpha              = alpha.c3
            , b                  = 20000. * umol.2.mol
            , gamma.resp         = gamma.c3
            , effarea.transp     = 1.0
            , rho                = NA
            , leaf.turnover.rate = 0.0
            , root.turnover.rate = 5.083700
            , bark.turnover.rate = 0.0
            , SLA                = 24.2
            , hgt.ref            = 1.3
            , b1Ht               = 25.18
            , b2Ht               = -0.04964
            , b1Bl.small         = 0.048
            , b2Bl.small         = 1.455
            , b1Bl.large         = 0.048
            , b2Bl.large         = 1.455
            , b1Bs.small         = 0.1617
            , b2Bs.small         = 2.4572
            , b1Bs.large         = 0.1617
            , b2Bs.large         = 2.4572
            , b1Ca               = 2.490154
            , b2Ca               = 0.8068806
            , b1Cl               = 0.3106775
            , b2Cl               = 1.098
            , b1Mh               = 0.8370557
            , b1Vol              = NA
            , b2Vol              = NA
            , b1Xb               = NA
            , b1Xs               = NA
            , hgt.min            = NA
            , hgt.max            = NA
            , qroot              = 1.1274
            , qsw                = NA
            , qbark              = NA
            , qwai               = NA
            , qrhob              = NA
            , agf.bs             = 0.7
            , orient.factor      = 0.25
            , clumping.factor    = 0.84
            , leaf.width         = lwidth.bltree
            , init.density       = 0.1
            , veg.hcap.min       = 7.65271E-01
            )

pft11 = list( name               = "Late hardwood"
            , key                = "LHW"
            , colour             = "#661859"
            , tropical           = FALSE
            , savannah           = FALSE
            , conifer            = FALSE
            , grass              = FALSE
            , liana              = FALSE
            , pathway            = 3
            , d0                 = d0.tree
            , vm.hor             = vm.hor
            , vm.base            = vm.base.c3
            , vm.decay.a         = vm.decay.a
            , vm.decay.b         = vm.decay.b
            , vm.low.temp        = vm.tcold.c3temp + t00
            , vm.high.temp       = vm.thot.c3temp  + t00
            , vm.decay.e.low     = vm.decay.ecold.c3
            , vm.decay.e.high    = vm.decay.ehot.c3
            , lr.hor             = lr.hor
            , lr.base            = lr.base.c3
            , lr.low.temp        = lr.tcold.c3temp + t00
            , lr.high.temp       = lr.thot.c3temp  + t00
            , lr.decay.e.low     = lr.decay.ecold.c3
            , lr.decay.e.high    = lr.decay.ehot.c3
            , vm0                = 6.982 * umol.2.mol
            , m                  = mphoto.c3 * 6.3949 / 8.0
            , alpha              = alpha.c3
            , b                  = 20000. * umol.2.mol
            , gamma.resp         = gamma.c3
            , effarea.transp     = 1.0
            , rho                = NA
            , leaf.turnover.rate = 0.0
            , root.turnover.rate = 5.070992
            , bark.turnover.rate = 0.4
            , SLA                = 60.0
            , hgt.ref            = 1.3
            , b1Ht               = 23.3874
            , b2Ht               = -0.05404
            , b1Bl.small         = 0.017
            , b2Bl.small         = 1.731
            , b1Bl.large         = 0.017
            , b2Bl.large         = 1.731
            , b1Bs.small         = 0.235
            , b2Bs.small         = 2.2518
            , b1Bs.large         = 0.235
            , b2Bs.large         = 2.2518
            , b1Ca               = 2.490154
            , b2Ca               = 0.8068806
            , b1Cl               = 0.3106775
            , b2Cl               = 1.098
            , b1Mh               = 0.8370557
            , b1Vol              = NA
            , b2Vol              = NA
            , b1Xb               = NA
            , b1Xs               = NA
            , hgt.min            = NA
            , hgt.max            = NA
            , qroot              = 1.1274
            , qsw                = NA
            , qbark              = NA
            , qwai               = NA
            , qrhob              = NA
            , agf.bs             = 0.7
            , orient.factor      = 0.25
            , clumping.factor    = 0.84
            , leaf.width         = lwidth.bltree
            , init.density       = 0.1
            , veg.hcap.min       = 1.60601E-01
            )

pft15 = list( name               = "Araucaria"
            , key                = "ARC"
            , colour             = "#7365B8"
            , tropical           = TRUE
            , savannah           = FALSE
            , conifer            = TRUE
            , grass              = FALSE
            , liana              = FALSE
            , pathway            = 3
            , d0                 = d0.tree
            , vm.hor             = vm.hor
            , vm.base            = vm.base.c3
            , vm.decay.a         = vm.decay.a
            , vm.decay.b         = vm.decay.b
            , vm.low.temp        = vm.tcold.aa + t00
            , vm.high.temp       = vm.thot.aa  + t00
            , vm.decay.e.low     = vm.decay.ecold.c3
            , vm.decay.e.high    = vm.decay.ehot.c3
            , lr.hor             = lr.hor
            , lr.base            = lr.base.c3
            , lr.low.temp        = lr.tcold.aa + t00
            , lr.high.temp       = lr.thot.aa  + t00
            , lr.decay.e.low     = lr.decay.ecold.c3
            , lr.decay.e.high    = lr.decay.ehot.c3
            , vm0                = 15.625  * vmfact.c3 * umol.2.mol
            , m                  = mphoto.aa
            , alpha              = alpha.c3
            , b                  = b.aa * umol.2.mol
            , gamma.resp         = gamma.aa
            , effarea.transp     = 2.0
            , rho                = 0.54
            , leaf.turnover.rate = 1./6.
            , root.turnover.rate = NA
            , bark.turnover.rate = NA
            , SLA                = 10.0
            , hgt.ref            = hgt.ref.trop
            , b1Ht               = b1Ht.trop
            , b2Ht               = b2Ht.trop
            , b1Bl.small         = NA
            , b2Bl.small         = NA
            , b1Bl.large         = NA
            , b2Bl.large         = NA
            , b1Bs.small         = NA
            , b2Bs.small         = NA
            , b1Bs.large         = NA
            , b2Bs.large         = NA
            , b1Ca               = 2.490154
            , b2Ca               = 0.8068806
            , b1Cl               = 0.3106775
            , b2Cl               = 1.098
            , b1Mh               = 0.8370557
            , b1Vol              = NA
            , b2Vol              = NA
            , b1Xb               = NA
            , b1Xs               = NA
            , hgt.min            = NA
            , hgt.max            = NA
            , qroot              = 1.0
            , qsw                = NA
            , qbark              = NA
            , qwai               = NA
            , qrhob              = NA
            , agf.bs             = 0.7
            , orient.factor      = orient.aa
            , clumping.factor    = clumping.aa
            , leaf.width         = lwidth.nltree
            , init.density       = 0.1
            , veg.hcap.min       = 2.19242E+01
            )

pft16 = list( name               = "C3 grass"
            , key                = "C3G"
            , colour             = "#85762B"
            , tropical           = TRUE
            , savannah           = FALSE
            , conifer            = FALSE
            , grass              = TRUE
            , liana              = FALSE
            , pathway            = 3
            , d0                 = d0.grass
            , vm.hor             = vm.hor
            , vm.base            = vm.base.c3
            , vm.decay.a         = vm.decay.a
            , vm.decay.b         = vm.decay.b
            , vm.low.temp        = vm.tcold.c3trop + t00
            , vm.high.temp       = vm.thot.c3trop  + t00
            , vm.decay.e.low     = vm.decay.ecold.c3
            , vm.decay.e.high    = vm.decay.ehot.c3
            , lr.hor             = lr.hor
            , lr.base            = lr.base.c3
            , lr.low.temp        = lr.tcold.c3trop + t00
            , lr.high.temp       = lr.thot.c3trop  + t00
            , lr.decay.e.low     = lr.decay.ecold.c3
            , lr.decay.e.high    = lr.decay.ehot.c3
            , vm0                = 18.750 * vmfact.c3 * umol.2.mol
            , m                  = mphoto.c3
            , alpha              = alpha.c3
            , b                  = b.c3 * umol.2.mol
            , gamma.resp         = gamma.c3
            , effarea.transp     = 1.0
            , rho                = 0.20
            , leaf.turnover.rate = ifelse(iallom %in% 3, 3.0, 2.0)
            , root.turnover.rate = NA
            , bark.turnover.rate = NA
            , SLA                = 22.7
            , hgt.ref            = hgt.ref.trop
            , b1Ht               = b1Ht.trop
            , b2Ht               = b2Ht.trop
            , b1Bl.small         = NA
            , b2Bl.small         = NA
            , b1Bl.large         = NA
            , b2Bl.large         = NA
            , b1Bs.small         = NA
            , b2Bs.small         = NA
            , b1Bs.large         = NA
            , b2Bs.large         = NA
            , b1Ca               = 2.490154
            , b2Ca               = 0.8068806
            , b1Cl               = 0.99
            , b2Cl               = 1.00
            , b1Mh               = 0.495
            , b1Vol              = NA
            , b2Vol              = NA
            , b1Xb               = NA
            , b1Xs               = NA
            , hgt.min            = NA
            , hgt.max            = NA
            , qroot              = 1.0
            , qsw                = NA
            , qbark              = NA
            , qwai               = NA
            , qrhob              = NA
            , agf.bs             = 0.7
            , orient.factor      = orient.grass
            , clumping.factor    = clumping.grass
            , leaf.width         = lwidth.grass
            , init.density       = 0.1
            , veg.hcap.min       = 7.30807E+00
            )

pft17 = list( name               = "Liana"
            , key                = "LNA"
            , colour             = "#332288"
            , tropical           = TRUE
            , savannah           = FALSE
            , conifer            = FALSE
            , grass              = FALSE
            , liana              = TRUE
            , pathway            = 3
            , d0                 = d0.tree
            , vm.hor             = vm.hor
            , vm.base            = vm.base.c3
            , vm.decay.a         = vm.decay.a
            , vm.decay.b         = vm.decay.b
            , vm.low.temp        = vm.tcold.aa + t00
            , vm.high.temp       = vm.thot.aa  + t00
            , vm.decay.e.low     = vm.decay.ecold.c3
            , vm.decay.e.high    = vm.decay.ehot.c3
            , lr.hor             = lr.hor
            , lr.base            = lr.base.c3
            , lr.low.temp        = lr.tcold.aa + t00
            , lr.high.temp       = lr.thot.aa  + t00
            , lr.decay.e.low     = lr.decay.ecold.c3
            , lr.decay.e.high    = lr.decay.ehot.c3
            , vm0                = 9.09  * vmfact.c3 * umol.2.mol
            , m                  = mphoto.aa
            , alpha              = alpha.c3
            , b                  = b.aa * umol.2.mol
            , gamma.resp         = gamma.aa
            , effarea.transp     = 2.0
            , rho                = 0.59
            , leaf.turnover.rate = 1.27
            , root.turnover.rate = NA
            , bark.turnover.rate = 0.0
            , SLA                = NA
            , hgt.ref            = hgt.ref.trop
            , b1Ht               = b1Ht.trop
            , b2Ht               = b2Ht.trop
            , b1Bl.small         = NA
            , b2Bl.small         = NA
            , b1Bl.large         = NA
            , b2Bl.large         = NA
            , b1Bs.small         = NA
            , b2Bs.small         = NA
            , b1Bs.large         = NA
            , b2Bs.large         = NA
            , b1Ca               = 2.490154
            , b2Ca               = 0.8068806
            , b1Cl               = 0.3106775
            , b2Cl               = 1.098
            , b1Mh               = 0.8370557
            , b1Vol              = NA
            , b2Vol              = NA
            , b1Xb               = NA
            , b1Xs               = NA
            , hgt.min            = NA
            , hgt.max            = NA
            , qroot              = 1.0
            , qsw                = NA
            , qbark              = NA
            , qwai               = NA
            , qrhob              = NA
            , agf.bs             = 0.7
            , orient.factor      = orient.aa
            , clumping.factor    = clumping.aa
            , leaf.width         = lwidth.nltree
            , init.density       = 0.1
            , veg.hcap.min       = 2.19242E+01
            )

#----- Derived PFTs. ----------------------------------------------------------------------#
pft12 = modifyList( x   = pft02
                  , val = list( name     = "Early savannah"
                              , key      = "ESV"
                              , colour   = "#CC839B"
                              , savannah = TRUE
                              )#end list
                  )#end modifyList

pft13 = modifyList( x   = pft03
                  , val = list( name     = "Mid savannah"
                              , key      = "MSV"
                              , colour   = "#AA4466"
                              , savannah = TRUE
                              )#end list
                  )#end modifyList

pft14 = modifyList( x   = pft04
                  , val = list( name     = "Late savannah"
                              , key      = "LSV"
                              , colour   = "#661832"
                              , savannah = TRUE
                              )#end list
                  )#end modifyList

pft18 = modifyList( x   = pft15
                  , val = list( name     = "Total"
                              , key      = "ALL"
                              , colour   = "#404040"
                              )#end list
                  )#end modifyList
#------------------------------------------------------------------------------------------#



#----- Build the structure of photosynthesis parameters by PFT. ---------------------------#
pft = list()
for (p in sequence(npft+1)){
  ppp  = sprintf("%2.2i",p)
  phph = paste0("pft",ppp)
  if (p == 1){
     pft = get(phph)
  }else{
     phmerge = get(phph)
     for (n in union(names(pft),names(phmerge))){
        pft[[n]] = c(pft[[n]],phmerge[[n]])
     } #end for
  }# end if
} #end for
pft = as.data.frame(pft,stringsAsFactors=FALSE)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Define minimum and maximum height based on life form and allometry.                  #
#------------------------------------------------------------------------------------------#
if (iallom %in% 3){
   pft$hgt.min = ifelse( test = pft$tropical
                       , yes  = ifelse(test=pft$grass,yes=0.15,no=0.5)
                       , no   = ifelse(test=pft$grass,yes=0.15,no=0.2+pft$hgt.ref)
                       )#end ifelse
   pft$hgt.max = ifelse( test = pft$tropical
                       , yes  = ifelse(test=pft$grass,yes=1.5         ,no=hgt.max.trop )
                       , no   = ifelse(test=pft$grass,yes=.95*pft$b1Ht,no=.999*pft$b1Ht)
                       )#end ifelse
}else{
   pft$hgt.min = ifelse( test = pft$tropical
                       , yes  = ifelse(test=pft$grass,yes=0.5,no=0.5+0.2*(iallom %in% 3))
                       , no   = ifelse(test=pft$grass,yes=0.15,no=0.2+pft$hgt.ref)
                       )#end ifelse
   pft$hgt.max = ifelse( test = pft$tropical
                       , yes  = ifelse(test=pft$grass,yes=1.5         ,no=hgt.max.trop )
                       , no   = ifelse(test=pft$grass,yes=.95*pft$b1Ht,no=.999*pft$b1Ht)
                       )#end ifelse
}#end if
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Set specific leaf area (SLA, m2leaf/kgC), and turnover rates for root and bark for   #
# those PFTs that have NA.  The curve relating SLA and leaf turnover rate came from        #
# fitting a model to trait data base (GLOPNET, W04).  The old model fitting was included   #
# in K12, whereas the new model is based on SMA fitting.                                   #
#                                                                                          #
# Wright, I. J., P. B. Reich, M. Westoby, et al., The worldwide leaf economics spectrum.   #
#    Nature, 428(6985):821-827, Apr 2004. doi:10.1038/nature02403 (W04).                   #
#                                                                                          #
# Kim, Y., R. G. Knox, M. Longo, D. Medvigy, L. R. Hutyra, E. H. Pyle, S. C. Wofsy,        #
#    R. L. Bras, and P. R. Moorcroft. Seasonal carbon dynamics and water fluxes in an      #
#    Amazon rainforest. Glob. Change Biol., 18 (4):1322 (K12).                             #
#                                                                                          #
# Chave, J., D. Coomes, S. Jansen, S. L. Lewis, N. G. Swenson, and A. E. Zanne. Towards a  #
#    worldwide wood economics spectrum. Ecol. Lett., 12(4):351-366, Apr 2009.              #
#    doi:10.1111/j.1461-0248.2009.01285.x (C09).                                           #
#------------------------------------------------------------------------------------------#
if (iallom %in% 3){
  sla.s0 = 19.41946059
  sla.s1 = 0.60550954
}else{
  sla.s0 = exp(log(0.1*C2B)+2.4*log(10.)-0.46*log(12.))
  sla.s1 = 0.46
}#end if
pft$SLA                = ifelse( test = is.finite(pft$SLA)
                               , yes  = pft$SLA
                               , no   = sla.s0*pft$leaf.turnover.rate^sla.s1
                               )#end ifelse
pft$root.turnover.rate = ifelse( test = is.finite(pft$root.turnover.rate)
                               , yes  = pft$root.turnover.rate
                               , no   = pft$leaf.turnover.rate
                               )#end ifelse
pft$bark.turnover.rate = ifelse( test = is.finite(pft$bark.turnover.rate)
                               , yes  = pft$bark.turnover.rate
                               , no   = 0.4
                               )#end ifelse
#------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------#
#     Set qsw according to the allometry and the PFT.                                      #
#------------------------------------------------------------------------------------------#
for (ipft in sequence(npft)){
   #---- Check PFT and allometry. ---------------------------------------------------------#
   if (pft$tropical[ipft] && pft$conifer[ipft] && iallom %in% 3){
      pft$qsw[ipft] = pft$SLA[ipft] * pft$rho[ipft] / sapwood.factor["aa"]
   }else if (pft$tropical[ipft] && pft$grass[ipft] && iallom %in% 3){
      pft$qsw[ipft] = 1.0e-5
   }else if (pft$tropical[ipft] && iallom %in% 3){
      pft$qsw[ipft] = pft$SLA[ipft] * pft$rho[ipft] / sapwood.factor["bl"]
   }else{
      pft$qsw[ipft] = pft$SLA[ipft] / sapwood.ratio.orig
   }#end if (pft$tropical[ipft] && is.finite(pft$rho[ipft]) && iallom %in% 3)
   #---------------------------------------------------------------------------------------#
}#end for (ipft in sequence(npft))
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Set the bark:wood density ratio (qrhob), and the water:biomass ratio for leaves      #
# (qwatdry.leaf), wood (qwatdry.wood), and bark (qwatdry.bark).  qwatdry.leaf is defined   #
# after G07.  The other variables were defined based on the PFT. For tropical broadleaf    #
# trees, we used a SMA analysis, similar to what was done for Vm0, SLA, and turnover, but  #
# using the data available in Table S1 (P14).  For other PFTs we use the default values    #
# from R14 and FPL10.                                                                      #
#                                                                                          #
# References:                                                                              #
#                                                                                          #
# Forest Products Laboratory. Wood handbook - wood as an engineering material. General     #
#    Technical Report FPL-GTR-190, U.S. Department of Agriculture, Madison, WI, 2010.      #
#    doi:10.2737/FPL-GTR-190 (FPL10)                                                       #
#                                                                                          #
# Gu, L., T. Meyers, S. G. Pallardy, P. J. Hanson, B. Yang, M. Heuer, K. P. Hosman,        #
#    Q. Liu, J. S. Riggs, D. Sluss, and S. D. Wullschleger. Influences of biomass heat and #
#    biochemical energy storages on the land surface fluxes and radiative temperature.     #
#    J. Geophys. Res., 112(D2):D02107, Jan 2007. doi:10.1029/2006JD007425 (G07)            #
#                                                                                          #
# Poorter, L., A. McNeil, V.-H. Hurtado, H. H. T. Prins, and F. E. Putz. Bark traits and   #
#    life-history strategies of tropical dry- and moist forest trees. Funct. Ecol.,        #
#    28(1):232-242, Feb 2014. doi:10.1111/1365-2435.12158 (P14)                            #
#                                                                                          #
# Rosell, J. A., S. Gleason, R. Mendez-Alonzo, Y. Chang, and M. Westoby. Bark              #
#    functional ecology: evidence for tradeoffs, functional coordination, and environ-     #
#    ment producing bark diversity. New Phytol., 201(2): 486-497, Jan 2014.                #
#    doi:10.1111/nph.12541. (R14)                                                          #
#------------------------------------------------------------------------------------------#
pft$qrhob        = rep(NA,times=npft+1)
pft$qwatdry.leaf = ifelse(pft$tropical,1.85,2.50)
pft$qwatdry.wood = rep(NA,times=npft+1)
pft$qwatdry.bark = rep(NA,times=npft+1)
pft$c.leaf.dry   = c(rep(3218.,times=npft),NA)
pft$c.wood.dry   = c(rep(103.1+3.867*spht.tref,times=npft),NA)
pft$c.bark.dry   = pft$c.wood.dry
for (ipft in sequence(npft)){
   usedef = with(pft,grass[ipft] || liana[ipft] || conifer[ipft] || (! tropical[ipft]))
   if (usedef){
      #----- Default values. --------------------------------------------------------------#
      pft$qrhob       [ipft] = 0.49/0.61
      pft$qwatdry.wood[ipft] = 0.7
      pft$qwatdry.bark[ipft] = 0.7
      #------------------------------------------------------------------------------------#
   }else{
      #----- Default values. --------------------------------------------------------------#
      pft$qrhob       [ipft] = exp(0.6966550 - 1.602123 * pft$rho[ipft])
      pft$qwatdry.wood[ipft] = exp(1.5018230 - 3.137476 * pft$rho[ipft])
      pft$qwatdry.bark[ipft] = exp(1.9892840 - 3.174365 * pft$rho[ipft])
      #------------------------------------------------------------------------------------#
   }#end if (usedef)
   #---------------------------------------------------------------------------------------#
}#end for (ipft in sequence(npft))
#------------------------------------------------------------------------------------------#

#----- Correction term for water-wood bond. -----------------------------------------------#
pft$delta.c.wood = ( 1.e5 * pmin(wdr.fs,pft$qwatdry.wood)
                   * ( - 0.06191 + 2.36e-4 * spht.tref 
                     - 1.33e-2 * pmin(wdr.fs,pft$qwatdry.wood)
                     )#end
                   )#end pft$delta.c.wood
pft$delta.c.bark = ( 1.e5 * pmin(wdr.fs,pft$qwatdry.bark)
                   * ( - 0.06191 + 2.36e-4 * spht.tref
                     - 1.33e-2 * pmin(wdr.fs,pft$qwatdry.bark)
                     )#end
                   )#end pft$delta.c.wood
pft$c.leaf = (pft$c.leaf.dry + pft$qwatdry.leaf * cliq) / (1. + pft$qwatdry.leaf)
pft$c.wood = ( (pft$c.wood.dry + pft$qwatdry.wood * cliq) / (1. + pft$qwatdry.wood)
             + pft$delta.c.wood
             )#end pft$c.wood
pft$c.bark = ( (pft$c.bark.dry + pft$qwatdry.bark * cliq) / (1. + pft$qwatdry.bark)
             + pft$delta.c.bark
             )#end pft$c.bark
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Set bark thickness and carbon allocation to bark.  This is currently done only for   #
# tropical trees when IALLOM=3, because all biomass pools must be corrected to ensure that #
# total aboveground biomass is consistent with the allometric equations.  This may and     #
# should be changed in the future.                                                         #
#                                                                                          #
# References:                                                                              #
#                                                                                          #
# Meinzer, F. C., G. Goldstein, and J. L. Andrade. Regulation of water flux through        #
#    tropical forest canopy trees: Do universal rules apply? Tree Physiol., 21(1):19-26,   #
#    Jan 2001. doi:10.1093/treephys/21.1.19. (M01)                                         #
#                                                                                          #
# de Mattos, P. P., A. T. dos Santos, H. Rivera, Y. M. M. de Oliveira, M. A. D. Rosot,     #
#    and M. C. Garrastazu. Growth of Araucaria angustifolia in the Embrapa/Epagri forest   #
#    reserve, Cacador, SC, Brazil. Pesq. Flor. Bras., 55(2):107-114, Jul 2007.             #
#    URL http://pfb.cnpf.embrapa.br/pfb/index.php/pfb/ article/view/124. In Portuguese.    #
#    (M07)                                                                                 #
#                                                                                          #
# Lawes, M. J. , J. J. Midgley, and P. J. Clarke. Costs and benefits of relative bark      #
#    thickness in relation to fire damage: a savanna/forest contrast. J. Ecol.,            #
#    101(2):517-524, Dec 2013. doi:10.1111/1365-2745.12035. (L13)                          #
#                                                                                          #
# Falster, D. S., R. G. FitzJohn, A. Brannstrom, U. Dieckmann, and M. Westoby.  plant: A   #
#    package for modelling forest trait ecology and evolution.  Methods Ecol. Evol., 7(2): #
#    136-146, Feb 2016. doi:10.1111/2041-210X.12525. (F16)                                 #
#------------------------------------------------------------------------------------------#
pft$b1Xs[sequence(npft)] = 0.315769481
for (ipft in sequence(npft)){
   skip = pft$grass[ipft] || pft$liana[ipft] || (! pft$tropical[ipft]) || (iallom != 3)
   if (skip){
      #------------------------------------------------------------------------------------#
      #   Set all bark variables to zero. in case this is not a tropical tree and in case  #
      # this is an old allometry.                                                          #
      #------------------------------------------------------------------------------------#
      pft$b1Xs [ipft] = 0.
      pft$b1Xb [ipft] = 0.
      pft$qbark[ipft] = 0.
      #------------------------------------------------------------------------------------#
   }else{
      #------------------------------------------------------------------------------------#
      #     Variable b1Xs is the ratio between sapwood thickness and DBH.  It is currently #
      # set to 0.316, based on a model fitting using M01 data.  This number is currently   #
      # used only to define carbon allocation to bark.                                     #
      #------------------------------------------------------------------------------------#
      pft$b1Xs[ipft] = 0.315769481
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Bark thickness slope depends on the life strategy.  Tropical broadleaf trees   #
      # use the meta-analysis by L13, and the tropical needleleaf trees use the slope      #
      # derived from M07's table 1.                                                        #
      #------------------------------------------------------------------------------------#
      if (pft$conifer[ipft]){
         pft$b1Xb[ipft] = 0.03936468
      }else if (pft$savannah[ipft]){
         pft$b1Xb[ipft] = 0.128
      }else{
         pft$b1Xb[ipft] = 0.019
      }#end (pft$conifer[ipft])
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     qbark is the factor that relates leaf biomass (and height) with bark biomass.  #
      #------------------------------------------------------------------------------------#
      abas            = ( pft$b1Xb[ipft] * (1.0 - pft$b1Xb[ipft]) 
                        / ( pft$b1Xs[ipft] * (1. + pft$b1Xs[ipft] - 2. * pft$b1Xb[ipft]) )
                        )#end abas
      pft$qbark[ipft] = pft$qrhob[ipft] * abas * pft$qsw[ipft]
      #------------------------------------------------------------------------------------#
   }#end if (skip)
   #---------------------------------------------------------------------------------------#
}#end for (ipft in sequence(npft))
#------------------------------------------------------------------------------------------#



#----- Minimum and Maximum DBH. -----------------------------------------------------------#
pft$dbh.min   = rep(NA,times=npft+1)
pft$dbh.crit  = rep(NA,times=npft+1)
for (ipft in sequence(npft)){
   if (pft$tropical[ipft]){
      if (iallom %in% c(0,1)){
         pft$dbh.min  [ipft] = exp((log(pft$hgt.min[ipft])-pft$b1Ht[ipft])/pft$b2Ht[ipft])
         pft$dbh.crit [ipft] = exp((log(pft$hgt.max[ipft])-pft$b1Ht[ipft])/pft$b2Ht[ipft])
      }else if (iallom %in% c(2,3)){
         pft$dbh.min [ipft] = ( log(   pft$hgt.ref[ipft]
                                   / ( pft$hgt.ref[ipft] - pft$hgt.min[ipft]) )
                              / pft$b1Ht[ipft] ) ^ (1.0 / pft$b2Ht[ipft])
         pft$dbh.crit[ipft] = ( log(   pft$hgt.ref[ipft]
                                   / ( pft$hgt.ref[ipft] - pft$hgt.max[ipft]) )
                              / pft$b1Ht[ipft] ) ^ (1.0 / pft$b2Ht[ipft])
     }#end if
   }else{
      pft$dbh.min [ipft] = ( log(1.0 - (pft$hgt.min[ipft]-pft$hgt.ref[ipft])
                                       / pft$b1Ht[ipft] ) / pft$b2Ht[ipft] )
      pft$dbh.crit[ipft] = ( log(1.0 - (pft$hgt.max[ipft]-pft$hgt.ref[ipft])
                                       / pft$b1Ht[ipft]) / pft$b2Ht[ipft] )
   }#end if
}#end for
#----- Initialise dbh.adult to be the same as dbh.crit.  This will change for iallom=3. ---#
pft$dbh.adult = pft$dbh.crit
#------------------------------------------------------------------------------------------#



#----- Constants shared by both bdead and bleaf -------------------------------------------#
a1    =  -1.981
b1    =   1.047
dcrit = 100.0
#----- Constants used by bdead only -------------------------------------------------------#
c1d   =   0.572
d1d   =   0.931
a2d   =  -1.086
b2d   =   0.876
c2d   =   0.604
d2d   =   0.871
#----- Constants used by bleaf only -------------------------------------------------------#
c1l   =  -0.584
d1l   =   0.550
a2l   =  -4.111
b2l   =   0.605
c2l   =   0.848
d2l   =   0.438
#------------------------------------------------------------------------------------------#
#pft$b1Bl       = rep(NA,times=npft+1)
#pft$b2Bl       = rep(NA,times=npft+1)
#pft$b1Bs.small = rep(NA,times=npft+1)
#pft$b2Bs.small = rep(NA,times=npft+1)
#pft$b1Bs.large = rep(NA,times=npft+1)
#pft$b2Bs.large = rep(NA,times=npft+1)
#pft$b1Ca       = rep(NA,times=npft+1)
#pft$b2Ca       = rep(NA,times=npft+1)
for (ipft in sequence(npft)){
   if (pft$tropical[ipft]){
      #------------------------------------------------------------------------------------#
      #      Fill in the structural biomass parameters.                                    #
      #------------------------------------------------------------------------------------#
      if (iallom %in% c(0)){
         #---- ED-2.1 allometry. ----------------------------------------------------------#
         pft$b1Bs.small[ipft] = exp(a1 + c1d * pft$b1Ht[ipft] + d1d * log(pft$rho[ipft]))
         pft$b1Bs.large[ipft] = exp(a1 + c1d * log(pft$hgt.max[ipft]) 
                                       + d1d * log(pft$rho[ipft]))
         aux                  = ( (a2d - a1) + pft$b1Ht[ipft] * (c2d - c1d)
                                + log(pft$rho[ipft]) * (d2d - d1d) ) * (1.0/log(dcrit))
         pft$b2Bs.small[ipft] = C2B * b2d + c2d * pft$b2Ht[ipft] + aux

         aux                  = ( (a2d - a1) + log(pft$hgt.max[ipft]) * (c2d - c1d)
                                + log(pft$rho[ipft]) * (d2d - d1d)) * (1.0/log(dcrit))
         pft$b2Bs.large[ipft] = C2B * b2d + aux

      }else if (iallom %in% c(1,2)){
         #---- Based an alternative modification of Chave et al. (2001) allometry. --------#
         pft$b1Bs.small[ipft] = C2B * exp(ndead.small[1]) * pft$rho[ipft] / ndead.small[3]
         pft$b2Bs.small[ipft] = ndead.small[2]
         pft$b1Bs.large[ipft] = C2B * exp(ndead.large[1]) * pft$rho[ipft] / ndead.large[3]
         pft$b2Bs.large[ipft] = ndead.large[2]
      }else if (iallom %in% c(3)){
         if (pft$grass[ipft]){
            #----- Small value just in case we are using old grass scheme. ----------------#
            pft$b1Bs.small[ipft] = 1.e-5
            pft$b2Bs.small[ipft] = 1.0
            pft$b1Bs.large[ipft] = pft$b1Bs.small[ipft]
            pft$b2Bs.large[ipft] = pft$b2Bs.small[ipft]
            #------------------------------------------------------------------------------#
         }else{
            #---- Based on a re-fit of the Chave et al. (2014) allometry. -----------------#
            pft$b1Bs.small[ipft] = C2B * 0.1684545 * pft$rho[ipft]
            pft$b2Bs.small[ipft] = 2.4629283
            pft$b2Bs.large[ipft] = 1.9670750
            pft$b1Bs.large[ipft] = ( pft$b1Bs.small[ipft] * pft$dbh.crit[ipft]
                                   ** (pft$b2Bs.small[ipft] - pft$b2Bs.large[ipft]) )
            #------------------------------------------------------------------------------#
         }#end if
      }#end if
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Replace the coefficients if we are going to use Poorter et al. (2006)         #
      # parameters for crown area.                                                         #
      #------------------------------------------------------------------------------------#
      if (iallom %in% c(0,1)){
         pft$b1Ca[ipft] = exp(-1.853) * exp(pft$b1Ht[ipft]) ^ 1.888
         pft$b2Ca[ipft] = pft$b2Ht[ipft] * 1.888
      }else if (iallom %in% c(2)){
         pft$b1Ca[ipft] = exp(ncrown.area[1])
         pft$b2Ca[ipft] = ncrown.area[2]
      }else if (iallom %in% c(3)){
         #---------------------------------------------------------------------------------#
         #     Allometry using the Sustainable Landscapes data.                            #
         #---------------------------------------------------------------------------------#
         #                                                                                 #
         #    Longo, M. et al.  Carbon Debt and Recovery time of degraded forests in       #
         #       the Amazon. Environ. Res. Lett., in prep.                                 #
         #                                                                                 #
         #    Equation was derived from forest inventory measurements carried out at       #
         # multiple locations in the Brazilian Amazon, and fitted using a                  #
         # heteroscedastic least squares approach.                                         #
         #                                                                                 #
         # Total number of trees: 17072                                                    #
         # b1Ca    = 0.582 (95% CI: [0.543;0.628])                                         #
         # b2Ca    = 1.224 (95% CI: [1.201;1.245])                                         #
         # R2      = 0.501                                                                 #
         # RMSE    = 29.89                                                                 #
         #---------------------------------------------------------------------------------#
         pft$b1Ca[ipft] = 0.582
         pft$b2Ca[ipft] = 1.224
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Replace b1Cl/b2Cl coefficients by those calculated by:                         #
      #                                                                                    #
      #    Longo, M. et al. 2016.  Carbon Debt and Recovery time of degraded forests in    #
      #       the Amazon. Biogeosciences, in prep.                                         #
      #                                                                                    #
      #    Equation was derived from forest inventory measurements carried out at multiple #
      # locations in the Brazilian Amazon, and fitted using a heteroscedastic least        #
      # squares approach.                                                                  #
      #                                                                                    #
      # Total number of trees: 16064                                                       #
      # b1Cl    = 0.298 (95% CI: [0.288;0.306])                                            #
      # b2Cl    = 1.032 (95% CI: [1.022;1.044])                                            #
      # R2      = 0.673                                                                    #
      # RMSE    = 2.29                                                                     #
      #------------------------------------------------------------------------------------#
      if (iallom %in% c(3) && (! pft$grass[ipft])){
         pft$b1Cl[ipft] = 0.29754
         pft$b2Cl[ipft] = 1.0324
      }#end if
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Fill in the leaf biomass parameters.                                          #
      #------------------------------------------------------------------------------------#
      if (iallom %in% c(0,1)){
         #---- ED-2.1 allometry. ----------------------------------------------------------#
         pft$b1Bl.small[ipft] = exp(a1 + c1l * pft$b1Ht[ipft] + d1l * log(pft$rho[ipft]))
         aux                  = ( (a2l - a1) + pft$b1Ht[ipft] * (c2l - c1l) 
                                 + log(pft$rho[ipft]) * (d2l - d1l)) * (1.0/log(dcrit))
         pft$b2Bl.small[ipft] = C2B * b2l + c2l * pft$b2Ht[ipft] + aux
         pft$b1Bl.large[ipft] = pft$b1Bl.small[ipft]
         pft$b2Bl.large[ipft] = pft$b2Bl.small[ipft]
      }else if(iallom %in% c(2)){
         pft$b1Bl.small[ipft] = C2B * exp(nleaf[1]) * pft$rho[ipft] / nleaf[3]
         pft$b2Bl.small[ipft] = nleaf[2]
         pft$b1Bl.large[ipft] = pft$b1Bl.small[ipft]
         pft$b2Bl.large[ipft] = pft$b2Bl.small[ipft]
      }else if(iallom %in% c(3)){
         #---------------------------------------------------------------------------------#
         #    Use L83 DBH-only equation for large trees, and a more modest exponential     #
         # slope for seedlings to avoid tiny leaf biomass at shorter classes.  For 
         # grasses, we apply a parameterised ratio for grasses that makes the individual-
         # level LAI to go from 0.5 at minimum height to 3.5 at maximum height, which is 
         # similar to the numbers for temperate grasses.                                   #
         #                                                                                 #
         # As a side note, these allometric equations diverge considerably at larger DBH   #
         # values.                                                                         #
         #                                                                                 #
         # References:                                                                     #
         #                                                                                 #
         # Lescure, J.-P., H. Puig, B. Riera, D. Leclerc, A. Beekman, and A. Beneteau.     #
         #    La phytomasse epigee d'une foret dense en Guyane Francaise.                  #
         #    Acta Ecol.-Oec. Gen., 4(3), 237--251, 1983.                                  #
         #    http://www.documentation.ird.fr/hor/fdi:010005089 (L83).                     #
         #                                                                                 #
         # Kenzo, T., R. Furutani, D. Hattori, J. J. Kendawang, S. Tanaka, K. Sakurai, and #
         #    I. Ninomiya. Allometric equations for accurate estimation of above-ground    #
         #    biomass in logged-over tropical rainforests in Sarawak, Malaysia.            #
         #    J. For. Res., 14(6):365, Sep 2009. doi:10.1007/s10310-009-0149-1 (K09).      #
         #---------------------------------------------------------------------------------#
         if (pft$grass[ipft]){
            lclai.max            = 3.0
            lclai.min            = 1.0
            pft$b2Bl.small[ipft] = ( pft$b2Ca[ipft]
                                    + log(lclai.max/lclai.min)
                                    / log(pft$dbh.crit[ipft] / pft$dbh.min[ipft])
                                    )#end pft$b2Bl.small
            pft$b1Bl.small[ipft] = ( lclai.min * C2B * pft$b1Ca[ipft] /pft$SLA[ipft]
                                   * pft$dbh.min[ipft]
                                   ^ (pft$b2Ca[ipft]-pft$b2Bl.small[ipft])
                                   )#end pft$b1Bl
            pft$b1Bl.large[ipft] = pft$b1Bl.small[ipft]
            pft$b2Bl.large[ipft] = pft$b2Bl.small[ipft]
         }else{
            pft$b1Bl.large[ipft] = 0.02439842 * pft$SLA[3] / pft$SLA[ipft]
            pft$b2Bl.large[ipft] = 1.86467176
            pft$dbh.adult [ipft] = ( ( C2B * pft$b1Ca[ipft] 
                                     / ( pft$SLA[ipft] * pft$b1Bl.large[ipft] ) )
                                     ^ ( 1.0 / ( pft$b2Bl.large[ipft] - pft$b2Ca[ipft] ) ) )
            pft$b1Bl.small[ipft] = C2B * pft$b1Ca[ipft] / pft$SLA[ipft]
            pft$b2Bl.small[ipft] = pft$b2Ca[ipft]
         }#end if
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#
   }#end if (pft$tropical[ipft]
   #---------------------------------------------------------------------------------------#
}#end for (ipft in sequence(npft))
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#    Ratio between WAI and maximum LAI, based on Olivas et al. (2013).  Currently this is  #
# applied to all PFTs.                                                                     #
#                                                                                          #
# Olivas, P. C., S. F. Oberbauer, D. B. Clark, D. A. Clark, M. G. Ryan, J. J. O'Brien, and #
#    H. Ordonez. Comparison of direct and indirect methods for assessing leaf area index   #
#    across a tropical rain forest landscape. Agric. For. Meteorol., 177:110-116,          #
#    Aug 2013. doi:10.1016/j.agrformet.2013.04.010.                                        #
#------------------------------------------------------------------------------------------#
pft$qwai = c(ifelse(pft$grass[-(npft+1)],0.00,0.11),NA)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Commercial volume of trees (stem/bole volume, in m3).  The current equation is a     #
# re-fit from Nogueira et al. (2008) so a single set of parameters can be used.  Their     #
# equation is for tropical trees only, so temperate and boreal forests may need a          #
# different set of parameters.  Grasses are assumed to have no commercial volume.          #
#                                                                                          #
# Nogueira, E. M., et al. Estimates of forest biomass in the Brazilian Amazon: new         #
#    allometric equations and adjustments to biomass from wood-volume inventories.         #
#    Forest Ecol. Manag., 256(11), 1853-1867, Nov. 2008, doi:10.1016/j.foreco.2008.07.022. #
#------------------------------------------------------------------------------------------#
pft$b1Vol = rep(NA,times=npft+1)
pft$b2Vol = rep(NA,times=npft+1)
for (ipft in sequence(npft)){
   if (pft$grass[ipft]){
      #----- Grasses have no commercial volume. -------------------------------------------#
      pft$b1Vol[ipft] = 0.0
      pft$b2Vol[ipft] = 1.0
      #------------------------------------------------------------------------------------#
   }else{
      #----- Nogueira et al. (2008) allometry. --------------------------------------------#
      pft$b1Vol[ipft] = 3.528e-5
      pft$b2Vol[ipft] = 0.976
      #------------------------------------------------------------------------------------#
   }#end if (pft$grass[ipft])
   #---------------------------------------------------------------------------------------#
}#end for (ipft %in% sequence(npft))
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Rooting depth coefficients.  They are updated to account for the volume.  The        #
# original volume equation did not have any source and it looks way too small (probably    #
# not SI units?), so I am replacing with commercial volume and making rooting depth self-  #
# -contained.                                                                              #
#------------------------------------------------------------------------------------------#
pft$b1Rd = rep(NA,times=npft+1)
pft$b2Rd = rep(NA,times=npft+1)
if (iallom %in% c(0)){
   #----- Original ED-2.1 scheme, based on standing volume. -------------------------------#
   for (ipft in sequence(npft)){
      if (pft$grass[ipft]){
         #----- Grasses have fixed rooting depth (70 cm). ---------------------------------#
         pft$b2Rd[ipft] = 0.0
         pft$b1Rd[ipft] = -0.700
         #---------------------------------------------------------------------------------#
      }else{
         #----- Volume-based allometry for trees. -----------------------------------------#
         pft$b2Rd[ipft] = 0.277
         pft$b1Rd[ipft] = -exp(0.545 *log(10.)) * (0.65*pi*0.11*0.11)^pft$b2Rd[ipft]
         #---------------------------------------------------------------------------------#
      }#end if (pft$grass[ipft])
      #------------------------------------------------------------------------------------#
   }#end for (ipft %in% sequence(npft))
   #---------------------------------------------------------------------------------------#
}else{
   #----- Simple allometry (0.5 m for seedlings, 5.0m for 35-m trees. ---------------------#
   pft$b1Rd[1:17]  = -1.1140580
   pft$b2Rd[1:17]  =  0.4223014
   #---------------------------------------------------------------------------------------#
}#end if
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#    Minimum bleaf and leaf area index that is resolvable.                                 #
#------------------------------------------------------------------------------------------#
pft$bleaf.min = c(dbh2bl(dbh=pft$dbh.min[1:npft],ipft=sequence(npft)),NA)
pft$lai.min   = onesixth * pft$init.dens * pft$bleaf.min * pft$SLA
#------------------------------------------------------------------------------------------#


#----- Reference leaf respiration. --------------------------------------------------------#
pft$lr0 = pft$gamma * pft$vm0
#------------------------------------------------------------------------------------------#


#----- Make it global. --------------------------------------------------------------------#
pft <<- pft
#------------------------------------------------------------------------------------------#
