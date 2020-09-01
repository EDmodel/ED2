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
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Bounds for internal carbon and water stomatal conductance.                           #
#------------------------------------------------------------------------------------------#
c34smin.lint.co2 <<- 0.001 * umol.2.mol # Minimum carbon dioxide concentration  [  mol/mol]
c34smax.lint.co2 <<- 3000. * umol.2.mol # Maximum carbon dioxide concentration  [  mol/mol]
c34smax.gsw      <<- 1.e+4              # Max. stomatal conductance (water)     [ mol/m2/s]
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Flag to decide whether to use effective functional rooting depth (sensu Brum et al.  #
# 2018).                                                                                   #
#------------------------------------------------------------------------------------------#
if (! "use.efrd.tree" %in% ls()) use.efrd.trtree <<- TRUE
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     This is an alternative way to express the temperature dependence, which is used by   #
# C91.  (tcollatz, in Kelvin).  fcollatz is the factor that multiply the temperature       #
# departure from reference.                                                                #
#------------------------------------------------------------------------------------------#
tphysref  <<- 15.0+t00
fcoll     <<- 0.1
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     The next two variables are the parameters for the Michaelis-Menten coefficients when #
# using Collatz equation.  Collatz equation uses partial pressure, so we must convert them #
# to mixing ratio.  Also, we use the reference at 15 C rather than 25 C, so Vm is the      #
# same.   The values themselves were obtained from von Caemmerer (2000).                   #
#------------------------------------------------------------------------------------------#
kco2.q10    <<- 2.24
kco2.refval <<- 40.4   / prefsea / kco2.q10
ko2.q10     <<- 1.63
ko2.refval  <<- 24800. / prefsea / ko2.q10
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Find the compensation point that is consistent with KCO2 and KO2, assuming that      #
# Vomax = 0.25 Vcmax (von Caemmerer 2000).                                                 #
#------------------------------------------------------------------------------------------#
compp.q10    <<- kco2.q10 / ko2.q10
compp.refval <<- o2.ref * kco2.refval / (8. * ko2.refval)
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


#----- Critical DBH for lianas. -----------------------------------------------------------#
liana.dbh.crit <<- 26.0
#------------------------------------------------------------------------------------------#


#----- Fudging parameters to try to tune photosynthesis. ----------------------------------#
if (! "alpha.c3"          %in% ls()) alpha.c3          <<- 0.08
if (! "alpha.c4"          %in% ls()) alpha.c4          <<- 0.040
if (! "phi.psII.c3"       %in% ls()) phi.psII.c3       <<- 0.85
if (! "phi.psII.c4"       %in% ls()) phi.psII.c4       <<- 0.85
if (! "curvpar.c3"        %in% ls()) curvpar.c3        <<- 0.70
if (! "curvpar.c4"        %in% ls()) curvpar.c4        <<- 0.70
if (! "vmfact.c3"         %in% ls()) vmfact.c3         <<- 1.0
if (! "vmfact.c4"         %in% ls()) vmfact.c4         <<- 1.0
if (! "jmfact.c3"         %in% ls()) jmfact.c3         <<- 1.767
if (! "jmfact.c4"         %in% ls()) jmfact.c4         <<- 1.767
if (! "tpmfact.c3"        %in% ls()) tpmfact.c3        <<- 0.110
if (! "mphoto.c3"         %in% ls()) mphoto.c3         <<- 9.0
if (! "mphoto.aa"         %in% ls()) mphoto.aa         <<- 7.2
if (! "mphoto.c4"         %in% ls()) mphoto.c4         <<- 4.0
if (! "gamma.c3"          %in% ls()) gamma.c3          <<- 0.015
if (! "gamma.aa"          %in% ls()) gamma.aa          <<- 0.015
if (! "gamma.c4"          %in% ls()) gamma.c4          <<- 0.025
if (! "d0.grass"          %in% ls()) d0.grass          <<- 0.016
if (! "d0.tree"           %in% ls()) d0.tree           <<- 0.016
if (! "klowin"            %in% ls()) klowin            <<- 0.7 / 39 * 1.e6
if (! "b.c3"              %in% ls()) b.c3              <<- 10000.
if (! "b.aa"              %in% ls()) b.aa              <<- 1000.
if (! "b.c4"              %in% ls()) b.c4              <<- 10000.
if (! "orient.tree"       %in% ls()) orient.tree       <<- +0.10
if (! "orient.aa"         %in% ls()) orient.aa         <<- +0.00
if (! "orient.grass"      %in% ls()) orient.grass      <<- -0.30
if (! "clumping.tree"     %in% ls()) clumping.tree     <<- 0.66
if (! "clumping.aa"       %in% ls()) clumping.aa       <<- 0.735
if (! "clumping.grass"    %in% ls()) clumping.grass    <<- 0.75
if (! "lwidth.grass"      %in% ls()) lwidth.grass      <<- 0.05
if (! "lwidth.bltree"     %in% ls()) lwidth.bltree     <<- 0.05
if (! "lwidth.nltree"     %in% ls()) lwidth.nltree     <<- 0.05
if (! "vm.q10.c3"         %in% ls()) vm.q10.c3         <<- 2.43 # 2.21
if (! "vm.q10.c4"         %in% ls()) vm.q10.c4         <<- 2.43 # 2.21
if (! "vm.tcold.c3temp"   %in% ls()) vm.tcold.c3temp   <<- 4.7137
if (! "vm.tcold.c3trop"   %in% ls()) vm.tcold.c3trop   <<- 10.0
if (! "vm.tcold.aa"       %in% ls()) vm.tcold.aa       <<- 4.7137
if (! "vm.tcold.c4"       %in% ls()) vm.tcold.c4       <<- 15.0
if (! "vm.thot.c3temp"    %in% ls()) vm.thot.c3temp    <<- 37.5
if (! "vm.thot.c3trop"    %in% ls()) vm.thot.c3trop    <<- 37.5
if (! "vm.thot.aa"        %in% ls()) vm.thot.aa        <<- 37.5
if (! "vm.thot.c4"        %in% ls()) vm.thot.c4        <<- 40.0
if (! "vm.decay.ecold.c3" %in% ls()) vm.decay.ecold.c3 <<- 0.3
if (! "vm.decay.ehot.c3"  %in% ls()) vm.decay.ehot.c3  <<- 0.6
if (! "vm.decay.ecold.c4" %in% ls()) vm.decay.ecold.c4 <<- 0.3
if (! "vm.decay.ehot.c4"  %in% ls()) vm.decay.ehot.c4  <<- 0.6
if (! "jm.q10.c3"         %in% ls()) jm.q10.c3         <<- 1.81 # 1.65
if (! "jm.q10.c4"         %in% ls()) jm.q10.c4         <<- 1.81 # 1.65
if (! "jm.tcold.c3temp"   %in% ls()) jm.tcold.c3temp   <<- vm.tcold.c3temp
if (! "jm.tcold.c3trop"   %in% ls()) jm.tcold.c3trop   <<- vm.tcold.c3trop
if (! "jm.tcold.aa"       %in% ls()) jm.tcold.aa       <<- vm.tcold.aa
if (! "jm.tcold.c4"       %in% ls()) jm.tcold.c4       <<- vm.tcold.c4
if (! "jm.thot.c3temp"    %in% ls()) jm.thot.c3temp    <<- vm.thot.c3temp
if (! "jm.thot.c3trop"    %in% ls()) jm.thot.c3trop    <<- vm.thot.c3trop
if (! "jm.thot.aa"        %in% ls()) jm.thot.aa        <<- vm.thot.aa
if (! "jm.thot.c4"        %in% ls()) jm.thot.c4        <<- vm.thot.c4
if (! "jm.decay.ecold.c3" %in% ls()) jm.decay.ecold.c3 <<- vm.decay.ecold.c3
if (! "jm.decay.ehot.c3"  %in% ls()) jm.decay.ehot.c3  <<- vm.decay.ehot.c3
if (! "jm.decay.ecold.c4" %in% ls()) jm.decay.ecold.c4 <<- vm.decay.ecold.c4
if (! "jm.decay.ehot.c4"  %in% ls()) jm.decay.ehot.c4  <<- vm.decay.ehot.c4
if (! "lr.q10.c3"         %in% ls()) lr.q10.c3         <<- 1.88 # 2.46
if (! "lr.q10.c4"         %in% ls()) lr.q10.c4         <<- 1.88 # 2.46
if (! "lr.tcold.c3temp"   %in% ls()) lr.tcold.c3temp   <<- vm.tcold.c3temp
if (! "lr.tcold.c3trop"   %in% ls()) lr.tcold.c3trop   <<- vm.tcold.c3trop
if (! "lr.tcold.aa"       %in% ls()) lr.tcold.aa       <<- vm.tcold.aa
if (! "lr.tcold.c4"       %in% ls()) lr.tcold.c4       <<- vm.tcold.c4
if (! "lr.thot.c3temp"    %in% ls()) lr.thot.c3temp    <<- vm.thot.c3temp
if (! "lr.thot.c3trop"    %in% ls()) lr.thot.c3trop    <<- vm.thot.c3trop
if (! "lr.thot.aa"        %in% ls()) lr.thot.aa        <<- vm.thot.aa
if (! "lr.thot.c4"        %in% ls()) lr.thot.c4        <<- vm.thot.c4
if (! "lr.decay.ecold.c3" %in% ls()) lr.decay.ecold.c3 <<- vm.decay.ecold.c3
if (! "lr.decay.ehot.c3"  %in% ls()) lr.decay.ehot.c3  <<- vm.decay.ehot.c3
if (! "lr.decay.ecold.c4" %in% ls()) lr.decay.ecold.c4 <<- vm.decay.ecold.c4
if (! "lr.decay.ehot.c4"  %in% ls()) lr.decay.ehot.c4  <<- vm.decay.ehot.c4
if (! "rr.q10.c3"         %in% ls()) rr.q10.c3         <<- lr.q10.c3 # 2.46
if (! "rr.q10.c4"         %in% ls()) rr.q10.c4         <<- lr.q10.c4 # 2.46
if (! "rr.tcold.c3temp"   %in% ls()) rr.tcold.c3temp   <<- lr.tcold.c3temp
if (! "rr.tcold.c3trop"   %in% ls()) rr.tcold.c3trop   <<- lr.tcold.c3trop
if (! "rr.tcold.aa"       %in% ls()) rr.tcold.aa       <<- lr.tcold.aa
if (! "rr.tcold.c4"       %in% ls()) rr.tcold.c4       <<- lr.tcold.c4
if (! "rr.thot.c3temp"    %in% ls()) rr.thot.c3temp    <<- lr.thot.c3temp
if (! "rr.thot.c3trop"    %in% ls()) rr.thot.c3trop    <<- lr.thot.c3trop
if (! "rr.thot.aa"        %in% ls()) rr.thot.aa        <<- lr.thot.aa
if (! "rr.thot.c4"        %in% ls()) rr.thot.c4        <<- lr.thot.c4
if (! "rr.decay.ecold.c3" %in% ls()) rr.decay.ecold.c3 <<- lr.decay.ecold.c3
if (! "rr.decay.ehot.c3"  %in% ls()) rr.decay.ehot.c3  <<- lr.decay.ehot.c3
if (! "rr.decay.ecold.c4" %in% ls()) rr.decay.ecold.c4 <<- lr.decay.ecold.c4
if (! "rr.decay.ehot.c4"  %in% ls()) rr.decay.ehot.c4  <<- lr.decay.ehot.c4
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#    The following parameter is the k coefficient in Foley et al. (1996) that is used to   #
# determine the CO2-limited photosynthesis for C4 grasses.                                 #
#------------------------------------------------------------------------------------------#
klowco2      <<- klowin 
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
# 4 - Similar to 3, but with a few changes by Xiangtao Xu.                                 #
#------------------------------------------------------------------------------------------#
if ("iallom" %in% ls()){
   iallom <<- iallom
}else{
   iallom <<- 3
}#end if
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#    These variables define the set of traits to use                                       #
# 0 - ED-2.1, mostly based on Reich et al. (1997), with a few updates using Wright et al.  #
#     (2004)                                                                               #
# 1 - Mostly TRY database, with additional data from NGEE-Tropics, RAINFOR, and GLOPNET.   #
#------------------------------------------------------------------------------------------#
if ("iecon" %in% ls()){
   iecon <<- iecon
}else{
   iecon <<- 1
}#end if
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Factors for leaf:sapwood biomass ratio.  The original ED-1 number is incorrect, and  #
# we keep it incorrect unless the PFT is tropical and allometry is set to 3, in which case #
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



#------------------------------------------------------------------------------------------#
#      Define reference heights for tropical allometry.                                    #
#------------------------------------------------------------------------------------------#
# IALLOM ........       0,       1,       2,       3,       4..............................#
hgt.ref.trop = c(NA_real_,NA_real_,    61.7,NA_real_,NA_real_)[iallom+1]
hgt.max.trop = c(     35.,     35.,     35.,     46.,     46.)[iallom+1]
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#   Coefficients for leaf and structural biomass (iallom = 3).  For adult individuals,     #
# we use the pantropical allometric equation from C14 that estimates AGB and the leaf      #
# biomass from an allometric equation derived from F15 data (tropical forest, wild         #
# flowering trees only), and the size- and site-dependent stratified sampling and          #
# aggregation (J17).  Total individual leaf area was fitted, so to get biomass we must     #
# divide by SLA.  The C2B term is added here but is removed when the coefficients are set. #
#                                                                                          #
#  References:                                                                             #
#                                                                                          #
#   Chave, J, Rejou-Mechain M, Burquez A, Chidumayo E, Colgan MS, Delitti WB, Duque A,     #
#      Eid T, Fearnside PM, Goodman RC et al. 2014. Improved allometric models to estimate #
#      the aboveground biomass of tropical trees. Glob. Change Biol., 20(10), 3177-3190.   #
#      doi:10.1111/gcb.12629 (C14).                                                        #
#                                                                                          #
#   Falster DS, Duursma RA, Ishihara MI, Barneche DR, FitzJohn RG, Vahammar A, Aiba M,     #
#      Ando M, Anten N, Aspinwall MJ et al. 2015. BAAD: a biomass and allometry database   #
#      for woody plants. Ecology, 96 (5):1445-1445. doi:10.1890/14-1889.1 (F15).           #
#                                                                                          #
#   Jucker T, Caspersen J, Chave J, Antin C, Barbier N, Bongers F, Dalponte M,             #
#      van Ewijk KY, Forrester DI, Haeni M et al. 2017. Allometric equations for           #
#      integrating remote sensing imagery into forest monitoring programmes.               #
#      Glob. Change Biol., 23(1):177-190. doi:10.1111/gcb.13388 (J17).                     #
#                                                                                          #
#------------------------------------------------------------------------------------------#
c14f15.bl.xx  = c(0.46769540,0.6410495) # c(0.09747026,0.7587000)
c14f15.bs.tf  = c(0.06080334,1.0044785) # c(0.08204475,0.9814422)
c14f15.bs.sv  = c(0.05602791,1.0093501) # c(0.08013971,0.9818603)
c14f15.bs.gr  = c(1.e-5,1.0) * c14f15.bl.xx
c14f15.ht.xx  = c( 0.5709,-0.1007, 0.6734)
c14f15.la.wd  = c(-0.5874, 0.5679, 0.5476)
SLA.ref       = 17.419
rho.ref       = 0.610
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
            , rho                = if(iecon == 1){0.08}else{0.20}
            , SLA                = if(iecon == 1){NA_real_}else{22.7}
            , c2n.leaf           = NA_real_
            , leaf.width         = NA_real_
            , vm0                = if(iecon == 1){NA_real_}else{12.500}
            , mort3              = if(iecon == 1){0.124}else{0.066}
            , leaf.turnover.rate = 2.0
            , root.turnover.rate = NA_real_
            , hgt.ref            = hgt.ref.trop
            , b1Ht               = NA_real_
            , b2Ht               = NA_real_
            , b1Bl               = NA_real_
            , b2Bl               = NA_real_
            , b1Bs.small         = NA_real_
            , b2Bs.small         = NA_real_
            , b1Bs.large         = NA_real_
            , b2Bs.large         = NA_real_
            , d1DBH.small        = NA_real_
            , d2DBH.small        = NA_real_
            , d1DBH.large        = NA_real_
            , d2DBH.large        = NA_real_
            , l1DBH              = NA_real_
            , l2DBH              = NA_real_
            , b1Ca               = 2.490154
            , b2Ca               = 0.8068806
            , b1Cl               = 0.99
            , b2Cl               = 1.00
            , veg.hcap.min       = 7.30807E+00
            , plant.min.temp     = 2.5
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
            , rho                = if(iecon == 1){0.450}else{0.53}
            , SLA                = NA_real_
            , c2n.leaf           = NA_real_
            , leaf.width         = NA_real_
            , vm0                = if(iecon == 1){NA_real_}else{18.75}
            , mort3              = NA_real_
            , leaf.turnover.rate = if(iecon == 1){NA_real_}else{1.0}
            , root.turnover.rate = NA_real_
            , hgt.ref            = hgt.ref.trop
            , b1Ht               = NA_real_
            , b2Ht               = NA_real_
            , b1Bl               = NA_real_
            , b2Bl               = NA_real_
            , b1Bs.small         = NA_real_
            , b2Bs.small         = NA_real_
            , b1Bs.large         = NA_real_
            , b2Bs.large         = NA_real_
            , d1DBH.small        = NA_real_
            , d2DBH.small        = NA_real_
            , d1DBH.large        = NA_real_
            , d2DBH.large        = NA_real_
            , l1DBH              = NA_real_
            , l2DBH              = NA_real_
            , b1Ca               = 2.490154
            , b2Ca               = 0.8068806
            , b1Cl               = 0.3106775
            , b2Cl               = 1.098
            , veg.hcap.min       = 9.53226E+00
            , plant.min.temp     = 2.5
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
            , rho                = if(iecon == 1){0.615}else{0.71}
            , SLA                = NA_real_
            , c2n.leaf           = NA_real_
            , leaf.width         = NA_real_
            , vm0                = if(iecon == 1){NA_real_}else{12.50}
            , mort3              = NA_real_
            , leaf.turnover.rate = if(iecon == 1){NA_real_}else{0.5}
            , root.turnover.rate = NA_real_
            , hgt.ref            = hgt.ref.trop
            , b1Ht               = NA_real_
            , b2Ht               = NA_real_
            , b1Bl               = NA_real_
            , b2Bl               = NA_real_
            , b1Bs.small         = NA_real_
            , b2Bs.small         = NA_real_
            , b1Bs.large         = NA_real_
            , b2Bs.large         = NA_real_
            , d1DBH.small        = NA_real_
            , d2DBH.small        = NA_real_
            , d1DBH.large        = NA_real_
            , d2DBH.large        = NA_real_
            , l1DBH              = NA_real_
            , l2DBH              = NA_real_
            , b1Ca               = 2.490154
            , b2Ca               = 0.8068806
            , b1Cl               = 0.3106775
            , b2Cl               = 1.098
            , veg.hcap.min       = 1.46161E+01
            , plant.min.temp     = 2.5
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
            , rho                = if(iecon == 1){0.790}else{0.90}
            , SLA                = NA_real_
            , c2n.leaf           = NA_real_
            , leaf.width         = NA_real_
            , vm0                = if(iecon == 1){NA_real_}else{6.25}
            , mort3              = NA_real_
            , leaf.turnover.rate = if(iecon == 1){NA_real_}else{1./3.}
            , root.turnover.rate = NA_real_
            , hgt.ref            = hgt.ref.trop
            , b1Ht               = NA_real_
            , b2Ht               = NA_real_
            , b1Bl               = NA_real_
            , b2Bl               = NA_real_
            , b1Bs.small         = NA_real_
            , b2Bs.small         = NA_real_
            , b1Bs.large         = NA_real_
            , b2Bs.large         = NA_real_
            , d1DBH.small        = NA_real_
            , d2DBH.small        = NA_real_
            , d1DBH.large        = NA_real_
            , d2DBH.large        = NA_real_
            , l1DBH              = NA_real_
            , l2DBH              = NA_real_
            , b1Ca               = 2.490154
            , b2Ca               = 0.8068806
            , b1Cl               = 0.3106775
            , b2Cl               = 1.098
            , veg.hcap.min       = 2.43602E+01
            , plant.min.temp     = 2.5
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
            , rho                = 0.08
            , SLA                = 22.0
            , c2n.leaf           = NA_real_
            , leaf.width         = NA_real_
            , leaf.turnover.rate = 2.0
            , root.turnover.rate = NA_real_
            , vm0                = 18.3
            , mort3              = 0.066
            , hgt.ref            = 0.0
            , b1Ht               = 0.4778
            , b2Ht               = -0.750
            , b1Bl               = 0.08
            , b2Bl               = 1.00
            , b1Bs.small         = 1.e-5
            , b2Bs.small         = 1.0
            , b1Bs.large         = 1.e-5
            , b2Bs.large         = 1.0
            , d1DBH.small        = NA_real_
            , d2DBH.small        = NA_real_
            , d1DBH.large        = NA_real_
            , d2DBH.large        = NA_real_
            , l1DBH              = NA_real_
            , l2DBH              = NA_real_
            , b1Ca               = 2.490154
            , b2Ca               = 0.8068806
            , b1Cl               = 0.99
            , b2Cl               = 1.0
            , veg.hcap.min       = 9.16551E+00
            , plant.min.temp     = -80
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
            , rho                = NA_real_
            , SLA                = 6.0
            , c2n.leaf           = NA_real_
            , leaf.width         = NA_real_
            , leaf.turnover.rate = 1./3.
            , root.turnover.rate = 3.927218
            , vm0                = 11.350
            , mort3              = 0.0033928
            , hgt.ref            = 1.3
            , b1Ht               = 27.14
            , b2Ht               = -0.03884
            , b1Bl               = 0.024
            , b2Bl               = 1.899
            , b1Bs.small         = 0.147
            , b2Bs.small         = 2.238
            , b1Bs.large         = 0.147
            , b2Bs.large         = 2.238
            , d1DBH.small        = NA_real_
            , d2DBH.small        = NA_real_
            , d1DBH.large        = NA_real_
            , d2DBH.large        = NA_real_
            , l1DBH              = NA_real_
            , l2DBH              = NA_real_
            , b1Ca               = 2.490154
            , b2Ca               = 0.8068806
            , b1Cl               = 0.3106775
            , b2Cl               = 1.098
            , veg.hcap.min       = 2.34683E-01
            , plant.min.temp     = -80
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
            , rho                = NA_real_
            , SLA                = 9.0
            , c2n.leaf           = NA_real_
            , leaf.width         = NA_real_
            , leaf.turnover.rate = 1./3.
            , root.turnover.rate = 4.117847
            , vm0                = 11.350
            , mort3              = 0.0043
            , hgt.ref            = 1.3
            , b1Ht               = 27.14
            , b2Ht               = -0.03884
            , b1Bl               = 0.024
            , b2Bl               = 1.899
            , b1Bs.small         = 0.147
            , b2Bs.small         = 2.238
            , b1Bs.large         = 0.147
            , b2Bs.large         = 2.238
            , d1DBH.small        = NA_real_
            , d2DBH.small        = NA_real_
            , d1DBH.large        = NA_real_
            , d2DBH.large        = NA_real_
            , l1DBH              = NA_real_
            , l2DBH              = NA_real_
            , b1Ca               = 2.490154
            , b2Ca               = 0.8068806
            , b1Cl               = 0.3106775
            , b2Cl               = 1.098
            , veg.hcap.min       = 2.34683E-01
            , plant.min.temp     = -10
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
            , rho                = NA_real_
            , SLA                = 10.0
            , c2n.leaf           = NA_real_
            , leaf.width         = NA_real_
            , leaf.turnover.rate = 1./3.
            , root.turnover.rate = 3.800132
            , vm0                = 4.540
            , mort3              = 0.0023568
            , hgt.ref            = 1.3
            , b1Ht               = 22.79
            , b2Ht               = -0.04445
            , b1Bl               = 0.0454
            , b2Bl               = 1.6829
            , b1Bs.small         = 0.1617
            , b2Bs.small         = 2.1536
            , b1Bs.large         = 0.1617
            , b2Bs.large         = 2.1536
            , d1DBH.small        = NA_real_
            , d2DBH.small        = NA_real_
            , d1DBH.large        = NA_real_
            , d2DBH.large        = NA_real_
            , l1DBH              = NA_real_
            , l2DBH              = NA_real_
            , b1Ca               = 2.490154
            , b2Ca               = 0.8068806
            , b1Cl               = 0.3106775
            , b2Cl               = 1.098
            , veg.hcap.min       = 6.80074E-01
            , plant.min.temp     = -60
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
            , rho                = NA_real_
            , SLA                = 30.0
            , c2n.leaf           = NA_real_
            , leaf.width         = NA_real_
            , leaf.turnover.rate = 0.0
            , root.turnover.rate = 5.772506
            , vm0                = 20.387
            , mort3              = 0.006144
            , hgt.ref            = 1.3
            , b1Ht               = 22.6799
            , b2Ht               = -0.06534
            , b1Bl               = 0.0129
            , b2Bl               = 1.7477
            , b1Bs.small         = 0.02648
            , b2Bs.small         = 2.95954
            , b1Bs.large         = 0.02648
            , b2Bs.large         = 2.95954
            , d1DBH.small        = NA_real_
            , d2DBH.small        = NA_real_
            , d1DBH.large        = NA_real_
            , d2DBH.large        = NA_real_
            , l1DBH              = NA_real_
            , l2DBH              = NA_real_
            , b1Ca               = 2.490154
            , b2Ca               = 0.8068806
            , b1Cl               = 0.3106775
            , b2Cl               = 1.098
            , veg.hcap.min       = 8.95049E-02
            , plant.min.temp     = -80
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
            , rho                = NA_real_
            , SLA                = 24.2
            , c2n.leaf           = NA_real_
            , leaf.width         = NA_real_
            , leaf.turnover.rate = 0.0
            , root.turnover.rate = 5.083700
            , vm0                = 17.455
            , mort3              = 0.003808
            , hgt.ref            = 1.3
            , b1Ht               = 25.18
            , b2Ht               = -0.04964
            , b1Bl               = 0.048
            , b2Bl               = 1.455
            , b1Bs.small         = 0.1617
            , b2Bs.small         = 2.4572
            , b1Bs.large         = 0.1617
            , b2Bs.large         = 2.4572
            , d1DBH.small        = NA_real_
            , d2DBH.small        = NA_real_
            , d1DBH.large        = NA_real_
            , d2DBH.large        = NA_real_
            , l1DBH              = NA_real_
            , l2DBH              = NA_real_
            , b1Ca               = 2.490154
            , b2Ca               = 0.8068806
            , b1Cl               = 0.3106775
            , b2Cl               = 1.098
            , veg.hcap.min       = 7.65271E-01
            , plant.min.temp     = -20
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
            , rho                = NA_real_
            , SLA                = 60.0
            , c2n.leaf           = NA_real_
            , leaf.width         = NA_real_
            , leaf.turnover.rate = 0.0
            , root.turnover.rate = 5.070992
            , vm0                = 6.982
            , mort3              = 0.00428
            , hgt.ref            = 1.3
            , b1Ht               = 23.3874
            , b2Ht               = -0.05404
            , b1Bl               = 0.017
            , b2Bl               = 1.731
            , b1Bs.small         = 0.235
            , b2Bs.small         = 2.2518
            , b1Bs.large         = 0.235
            , b2Bs.large         = 2.2518
            , d1DBH.small        = NA_real_
            , d2DBH.small        = NA_real_
            , d1DBH.large        = NA_real_
            , d2DBH.large        = NA_real_
            , l1DBH              = NA_real_
            , l2DBH              = NA_real_
            , b1Ca               = 2.490154
            , b2Ca               = 0.8068806
            , b1Cl               = 0.3106775
            , b2Cl               = 1.098
            , veg.hcap.min       = 1.60601E-01
            , plant.min.temp     = -20
            )

pft12 = list( name               = "Early savannah"
            , key                = "ESV"
            , colour             = "#CC839B"
            , tropical           = TRUE
            , savannah           = TRUE
            , conifer            = FALSE
            , grass              = FALSE
            , liana              = FALSE
            , pathway            = 3
            , rho                = if(iecon == 1){0.461}else{0.53}
            , SLA                = NA_real_
            , c2n.leaf           = NA_real_
            , leaf.width         = NA_real_
            , vm0                = if(iecon == 1){NA_real_}else{18.75}
            , mort3              = NA_real_
            , leaf.turnover.rate = if(iecon == 1){NA_real_}else{1.0}
            , root.turnover.rate = NA_real_
            , hgt.ref            = hgt.ref.trop
            , b1Ht               = NA_real_
            , b2Ht               = NA_real_
            , b1Bl               = NA_real_
            , b2Bl               = NA_real_
            , b1Bs.small         = NA_real_
            , b2Bs.small         = NA_real_
            , b1Bs.large         = NA_real_
            , b2Bs.large         = NA_real_
            , d1DBH.small        = NA_real_
            , d2DBH.small        = NA_real_
            , d1DBH.large        = NA_real_
            , d2DBH.large        = NA_real_
            , l1DBH              = NA_real_
            , l2DBH              = NA_real_
            , b1Ca               = 2.490154
            , b2Ca               = 0.8068806
            , b1Cl               = 0.3106775
            , b2Cl               = 1.098
            , veg.hcap.min       = 9.53226E+00
            , plant.min.temp     = 2.5
            )

pft13 = list( name               = "Mid savannah"
            , key                = "MSV"
            , colour             = "#AA4466"
            , tropical           = TRUE
            , savannah           = TRUE
            , conifer            = FALSE
            , grass              = FALSE
            , liana              = FALSE
            , pathway            = 3
            , rho                = if(iecon == 1){0.612}else{0.71}
            , SLA                = NA_real_
            , c2n.leaf           = NA_real_
            , leaf.width         = NA_real_
            , vm0                = if(iecon == 1){NA_real_}else{12.50}
            , mort3              = NA_real_
            , leaf.turnover.rate = if(iecon == 1){NA_real_}else{0.5}
            , root.turnover.rate = NA_real_
            , hgt.ref            = hgt.ref.trop
            , b1Ht               = NA_real_
            , b2Ht               = NA_real_
            , b1Bl               = NA_real_
            , b2Bl               = NA_real_
            , b1Bs.small         = NA_real_
            , b2Bs.small         = NA_real_
            , b1Bs.large         = NA_real_
            , b2Bs.large         = NA_real_
            , d1DBH.small        = NA_real_
            , d2DBH.small        = NA_real_
            , d1DBH.large        = NA_real_
            , d2DBH.large        = NA_real_
            , l1DBH              = NA_real_
            , l2DBH              = NA_real_
            , b1Ca               = 2.490154
            , b2Ca               = 0.8068806
            , b1Cl               = 0.3106775
            , b2Cl               = 1.098
            , veg.hcap.min       = 9.53226E+00
            , plant.min.temp     = 2.5
            )

pft14 = list( name               = "Late Savannah"
            , key                = "LSV"
            , colour             = "#661832"
            , tropical           = TRUE
            , savannah           = TRUE
            , conifer            = FALSE
            , grass              = FALSE
            , liana              = FALSE
            , pathway            = 3
            , rho                = if(iecon == 1){0.756}else{0.90}
            , SLA                = NA_real_
            , c2n.leaf           = NA_real_
            , leaf.width         = NA_real_
            , vm0                = if(iecon == 1){NA_real_}else{6.25}
            , mort3              = NA_real_
            , leaf.turnover.rate = if(iecon == 1){NA_real_}else{1./3.}
            , root.turnover.rate = NA_real_
            , hgt.ref            = hgt.ref.trop
            , b1Ht               = NA_real_
            , b2Ht               = NA_real_
            , b1Bl               = NA_real_
            , b2Bl               = NA_real_
            , b1Bs.small         = NA_real_
            , b2Bs.small         = NA_real_
            , b1Bs.large         = NA_real_
            , b2Bs.large         = NA_real_
            , d1DBH.small        = NA_real_
            , d2DBH.small        = NA_real_
            , d1DBH.large        = NA_real_
            , d2DBH.large        = NA_real_
            , l1DBH              = NA_real_
            , l2DBH              = NA_real_
            , b1Ca               = 2.490154
            , b2Ca               = 0.8068806
            , b1Cl               = 0.3106775
            , b2Cl               = 1.098
            , veg.hcap.min       = 9.53226E+00
            , plant.min.temp     = 2.5
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
            , rho                = 0.52
            , SLA                = 6.324555
            , c2n.leaf           = 86.29189
            , leaf.width         = NA_real_
            , leaf.turnover.rate = 0.04160842
            , root.turnover.rate = NA_real_
            , vm0                = 15.625
            , mort3              = 0.001111
            , hgt.ref            = hgt.ref.trop
            , b1Ht               = NA_real_
            , b2Ht               = NA_real_
            , b1Bl               = NA_real_
            , b2Bl               = NA_real_
            , b1Bs.small         = NA_real_
            , b2Bs.small         = NA_real_
            , b1Bs.large         = NA_real_
            , b2Bs.large         = NA_real_
            , d1DBH.small        = NA_real_
            , d2DBH.small        = NA_real_
            , d1DBH.large        = NA_real_
            , d2DBH.large        = NA_real_
            , l1DBH              = NA_real_
            , l2DBH              = NA_real_
            , b1Ca               = 2.490154
            , b2Ca               = 0.8068806
            , b1Cl               = 0.3106775
            , b2Cl               = 1.098
            , veg.hcap.min       = 2.19242E+01
            , plant.min.temp     = -15
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
            , rho                = if(iecon == 1){0.08}else{0.20}
            , SLA                = if(iecon == 1){NA_real_}else{22.7}
            , c2n.leaf           = NA_real_
            , leaf.width         = NA_real_
            , leaf.turnover.rate = 2.0
            , root.turnover.rate = NA_real_
            , vm0                = if(iecon == 1){NA_real_}else{18.75}
            , mort3              = if(iecon == 1){0.124}else{0.066}
            , hgt.ref            = hgt.ref.trop
            , b1Ht               = NA_real_
            , b2Ht               = NA_real_
            , b1Bl               = NA_real_
            , b2Bl               = NA_real_
            , b1Bs.small         = NA_real_
            , b2Bs.small         = NA_real_
            , b1Bs.large         = NA_real_
            , b2Bs.large         = NA_real_
            , d1DBH.small        = NA_real_
            , d2DBH.small        = NA_real_
            , d1DBH.large        = NA_real_
            , d2DBH.large        = NA_real_
            , l1DBH              = NA_real_
            , l2DBH              = NA_real_
            , b1Ca               = 2.490154
            , b2Ca               = 0.8068806
            , b1Cl               = 0.99
            , b2Cl               = 1.00
            , veg.hcap.min       = 7.30807E+00
            , plant.min.temp     = 2.5
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
            , rho                = 0.59
            , SLA                = NA_real_
            , c2n.leaf           = NA_real_
            , leaf.width         = NA_real_
            , leaf.turnover.rate = 1.27
            , root.turnover.rate = NA_real_
            , vm0                = 9.09
            , mort3              = 0.06311576
            , hgt.ref            = hgt.ref.trop
            , b1Ht               = NA_real_
            , b2Ht               = NA_real_
            , b1Bl               = NA_real_
            , b2Bl               = NA_real_
            , b1Bs.small         = NA_real_
            , b2Bs.small         = NA_real_
            , b1Bs.large         = NA_real_
            , b2Bs.large         = NA_real_
            , d1DBH.small        = NA_real_
            , d2DBH.small        = NA_real_
            , d1DBH.large        = NA_real_
            , d2DBH.large        = NA_real_
            , l1DBH              = NA_real_
            , l2DBH              = NA_real_
            , b1Ca               = 2.490154
            , b2Ca               = 0.8068806
            , b1Cl               = 0.3106775
            , b2Cl               = 1.098
            , veg.hcap.min       = 2.19242E+01
            , plant.min.temp     = 2.5
            )

#----- PFT 18 has the median wood density and SLA for tropical PFTs. ----------------------#
pft18 = modifyList( x   = pft03
                  , val = list( name     = "Total"
                              , key      = "ALL"
                              , colour   = "#404040"
                              , rho      = rho.ref
                              , SLA      = SLA.ref
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



#----- Adjust Vm0 based on settings and to make units consistent. -------------------------#
pft$vm0 = ifelse( test = pft$pathway == 4
                , yes  = pft$vm0 * vmfact.c4
                , no   = ifelse( test = pft$tropical & (! pft$liana)
                               , yes  = pft$vm0 * vmfact.c3
                               , no   = pft$vm0
                               )#end ifelse
                )#end ifelse
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Set specific leaf area (SLA, m2leaf/kgC), and turnover rates for root and bark for   #
# those PFTs that have NA.  The curve relating SLA and leaf turnover rate came from W04,   #
# after accounting for unit changes.  The old model fitting was included in K12, whereas   #
# the new model is based on SMA fitting.                                                   #
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
bltree  = pft$tropical & (! pft$conifer) & (! pft$liana) & (! pft$grass)
ltor.s0 = ifelse( test = rep(iecon,npft+1) == 1
                , yes  = 2.2349
                , no   = NA_real_
                )#end ifelse
ltor.s1 = ifelse( test = rep(iecon,npft+1) == 1
                , yes  = -3.98423
                , no   = NA_real_
                )#end ifelse
sla.s0  = ifelse( test = rep(iecon,npft+1) == 1
                , yes  = ifelse( test = pft$grass
                               , yes  = 15.15900
                               , no   = 21.66777
                               )#end if
                , no   = ifelse( test = pft$grass
                               , yes  = 22.7
                               , no   = exp(log(.1*C2B) + 2.4 * log(10.) - 0.46 * log(12.))
                               )#end if
                )#end ifelse
sla.s1  = ifelse( test = rep(iecon,npft+1) == 1
                , yes  = ifelse( test = pft$grass
                               , yes  = 0.8637294
                               , no   = 0.4314902
                               )#end if
                , no   = ifelse( test = pft$grass
                               , yes  = 0.0
                               , no   = 0.46
                               )#end if
                )#end ifelse
vm0.s0  = ifelse( test = pft$grass & iecon == 1
                , yes  = ifelse( test = pft$pathway == 4
                               , yes  = 51.6/vm.q10.c4
                               , no   = 78.2/vm.q10.c3
                               )#end ifelse
                , no   = ifelse( test = bltree & iecon == 1
                               , yes  = 7.82655 / vm.q10.c3
                               , no   = NA_real_
                               )#end ifelse
                )#end ltor.s0
vm0.s1  = ifelse( test = pft$grass & iecon == 1
                , yes  = ifelse( test = pft$pathway == 4
                               , yes  = 0.0
                               , no   = 0.0
                               )#end ifelse
                , no   = ifelse( test = bltree & iecon == 1
                               , yes  = 1.56402 - 1.0 # -1 is to convert to area based.
                               , no   = NA_real_
                               )#end ifelse
                )#end ltor.s0
vm0.ref = ifelse( test = pft$tropical | pft$grass
                , yes  = pft$vm0
                , no   = mapply( FUN      = switch
                               , EXPR     = pft$key
                               , MoreArgs = list(LCN=6.25,EHW=18.25,LHW=6.25,15.625)
                               )#end mapply
                )#end ifelse
#---- Tree and grasses had almost the same curve for SLA and C:N ratio. Fit one curve. ----#
c2nl.s0                = ifelse( test = (pft$grass | bltree) & iecon == 1
                               , yes  = 337.959
                               , no   = 1000./ ( 0.11289 + 0.129470 * vm0.ref )
                               )#end ltor.s0
c2nl.s1                = ifelse( test = (pft$grass | bltree) & iecon == 1
                               , yes  = -0.834527
                               , no   = -1.0
                               )#end ltor.s0
#------------------------------------------------------------------------------------------#



#----- Define traits that are scaled using standard major axis. ---------------------------#
pft$leaf.turnover.rate = ifelse( test = is.finite(pft$leaf.turnover.rate)
                               , yes  = pft$leaf.turnover.rate
                               , no   = ifelse( test = is.finite(pft$SLA)
                                              , yes  = (pft$SLA/sla.s0)^(1./sla.s1)
                                              , no   = exp(ltor.s0 + ltor.s1 * pft$rho)
                                              )#end ifelse
                               )#end ifelse
pft$SLA                = ifelse( test = is.finite(pft$SLA)
                               , yes  = pft$SLA
                               , no   = sla.s0 * pft$leaf.turnover.rate ^ sla.s1
                               )#end ifelse
pft$c2n.leaf           = ifelse( test = is.finite(pft$c2n.leaf)
                               , yes  = pft$c2n.leaf
                               , no   = c2nl.s0 * pft$SLA ^ c2nl.s1
                               )#end ifelse
pft$vm0                = ifelse( test = is.finite(pft$vm0)
                               , yes  = pft$vm0
                               , no   = vm0.s0 * pft$SLA ^ vm0.s1
                               )#end ifelse
pft$root.turnover.rate = ifelse( test = is.finite(pft$root.turnover.rate)
                               , yes  = pft$root.turnover.rate
                               , no   = ifelse(iecon == 1,0.9,1.0) * pft$leaf.turnover.rate
                               )#end ifelse
pft$bark.turnover.rate = ifelse( test = bltree & iecon == 1, yes = 0.4,no = 0.0)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Fill in ageing mortality.                                                            #
#------------------------------------------------------------------------------------------#

pft$mort3  = ifelse( test = is.finite(pft$mort3)
                   , yes  = pft$mort3
                   , no   = if (iecon == 1){
                               exp(-0.396-5.20*pft$rho) - 1/126.
                            }else{
                               0.15 * (1. - pft$rho/0.9)
                            }#end if (iecon == 1)
                   )#end ifelse
#------------------------------------------------------------------------------------------#


#----- Fill some allometric terms. --------------------------------------------------------#
pft$qroot  = ifelse( test = pft$tropical | pft$grass
                   , yes  = 1.0
                   , no   = ifelse(test=pft$conifer,yes=0.3463,no=1.1274)
                   )#end ifelse
pft$b1Mh   = ifelse( test = pft$grass
                   , yes  = 0.495
                   , no   = 0.8370557
                   )#end ifelse
pft$agf.bs = 0.7
pft$init.density = 0.1
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Fill in allometric coefficients.                                                     #
#------------------------------------------------------------------------------------------#
pft$qroot  = ifelse( test = pft$tropical | pft$grass
                   , yes  = 1.0
                   , no   = ifelse(test=pft$conifer,yes=0.3463,no=1.1274)
                   )#end ifelse
pft$b1Mh   = ifelse( test = pft$grass
                   , yes  = 0.495
                   , no   = 0.8370557
                   )#end ifelse
pft$agf.bs = 0.7
pft$init.density = 0.1
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Fill in optical and biophysical properties.                                          #
#------------------------------------------------------------------------------------------#
pft$orient.factor = ifelse( test = pft$grass
                          , yes  = orient.grass
                          , no   = ifelse( test = pft$conifer
                                         , yes  = orient.aa
                                         , no   = orient.tree
                                         )#end ifelse
                          )#end ifelse
pft$clumping.factor = ifelse( test = pft$grass
                            , yes  = clumping.grass
                            , no   = ifelse( test = pft$conifer
                                           , yes  = clumping.aa
                                           , no   = clumping.tree
                                           )#end ifelse
                            )#end ifelse
pft$leaf.width      = ifelse( test = is.finite(pft$leaf.width)
                            , yes  = pft$leaf.width
                            , no   = ifelse( test = pft$grass
                                           , yes  = lwidth.grass
                                           , no   = ifelse( test = pft$conifer
                                                          , yes  = lwidth.nltree
                                                          , no   = lwidth.bltree
                                                          )#end ifelse
                                           )#end ifelse
                            )#end ifelse
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Fill in the photosynthesis parameters.                                               #
#------------------------------------------------------------------------------------------#
pft$d0             = ifelse( test = pft$grass
                           , yes  = d0.grass
                           , no   = d0.tree
                           )#end ifelse
pft$vm.q10         = ifelse( test = pft$pathway == 4
                           , yes  = vm.q10.c4
                           , no   = vm.q10.c3
                           )#end ifelse
pft$vm.low.temp    = ifelse( test = pft$pathway == 4
                           , yes  = vm.tcold.c4
                           , no   = ifelse( test = pft$tropical & pft$conifer
                                          , yes  = vm.tcold.aa
                                          , no   = ifelse( test = pft$tropical
                                                         , yes  = vm.tcold.c3trop
                                                         , no   = vm.tcold.c3temp
                                                         )#end ifelse
                                          )#end ifelse
                           )#end ifelse
pft$vm.high.temp   = ifelse( test = pft$pathway == 4
                           , yes  = vm.thot.c4
                           , no   = ifelse( test = pft$tropical & pft$conifer
                                          , yes  = vm.thot.aa
                                          , no   = ifelse( test = pft$tropical
                                                         , yes  = vm.thot.c3trop
                                                         , no   = vm.thot.c3temp
                                                         )#end ifelse
                                          )#end ifelse
                           )#end ifelse
pft$vm.decay.ecold = ifelse( test = pft$pathway == 4
                           , yes  = vm.decay.ecold.c4
                           , no   = vm.decay.ecold.c3
                           )#end ifelse
pft$vm.decay.ehot = ifelse( test = pft$pathway == 4
                           , yes  = vm.decay.ehot.c4
                           , no   = vm.decay.ehot.c3
                           )#end ifelse
pft$jm.q10         = ifelse( test = pft$pathway == 4
                           , yes  = jm.q10.c4
                           , no   = jm.q10.c3
                           )#end ifelse
pft$jm.low.temp    = ifelse( test = pft$pathway == 4
                           , yes  = jm.tcold.c4
                           , no   = ifelse( test = pft$tropical & pft$conifer
                                          , yes  = jm.tcold.aa
                                          , no   = ifelse( test = pft$tropical
                                                         , yes  = jm.tcold.c3trop
                                                         , no   = jm.tcold.c3temp
                                                         )#end ifelse
                                          )#end ifelse
                           )#end ifelse
pft$jm.high.temp   = ifelse( test = pft$pathway == 4
                           , yes  = jm.thot.c4
                           , no   = ifelse( test = pft$tropical & pft$conifer
                                          , yes  = jm.thot.aa
                                          , no   = ifelse( test = pft$tropical
                                                         , yes  = jm.thot.c3trop
                                                         , no   = jm.thot.c3temp
                                                         )#end ifelse
                                          )#end ifelse
                           )#end ifelse
pft$jm.decay.ecold = ifelse( test = pft$pathway == 4
                           , yes  = jm.decay.ecold.c4
                           , no   = jm.decay.ecold.c3
                           )#end ifelse
pft$jm.decay.ehot = ifelse( test = pft$pathway == 4
                           , yes  = jm.decay.ehot.c4
                           , no   = jm.decay.ehot.c3
                           )#end ifelse
pft$lr.q10         = ifelse( test = pft$pathway == 4
                           , yes  = lr.q10.c4
                           , no   = lr.q10.c3
                           )#end ifelse
pft$lr.low.temp    = ifelse( test = pft$pathway == 4
                           , yes  = lr.tcold.c4
                           , no   = ifelse( test = pft$tropical & pft$conifer
                                          , yes  = lr.tcold.aa
                                          , no   = ifelse( test = pft$tropical
                                                         , yes  = lr.tcold.c3trop
                                                         , no   = lr.tcold.c3temp
                                                         )#end ifelse
                                          )#end ifelse
                           )#end ifelse
pft$lr.high.temp   = ifelse( test = pft$pathway == 4
                           , yes  = lr.thot.c4
                           , no   = ifelse( test = pft$tropical & pft$conifer
                                          , yes  = lr.thot.aa
                                          , no   = ifelse( test = pft$tropical
                                                         , yes  = lr.thot.c3trop
                                                         , no   = lr.thot.c3temp
                                                         )#end ifelse
                                          )#end ifelse
                           )#end ifelse
pft$lr.decay.ecold = ifelse( test = pft$pathway == 4
                           , yes  = lr.decay.ecold.c4
                           , no   = lr.decay.ecold.c3
                           )#end ifelse
pft$lr.decay.ehot  = ifelse( test = pft$pathway == 4
                           , yes  = lr.decay.ehot.c4
                           , no   = lr.decay.ehot.c3
                           )#end ifelse
pft$rr.q10         = ifelse( test = pft$pathway == 4
                           , yes  = rr.q10.c4
                           , no   = rr.q10.c3
                           )#end ifelse
pft$rr.low.temp    = ifelse( test = pft$pathway == 4
                           , yes  = rr.tcold.c4
                           , no   = ifelse( test = pft$tropical & pft$conifer
                                          , yes  = rr.tcold.aa
                                          , no   = ifelse( test = pft$tropical
                                                         , yes  = rr.tcold.c3trop
                                                         , no   = rr.tcold.c3temp
                                                         )#end ifelse
                                          )#end ifelse
                           )#end ifelse
pft$rr.high.temp   = ifelse( test = pft$pathway == 4
                           , yes  = rr.thot.c4
                           , no   = ifelse( test = pft$tropical & pft$conifer
                                          , yes  = rr.thot.aa
                                          , no   = ifelse( test = pft$tropical
                                                         , yes  = rr.thot.c3trop
                                                         , no   = rr.thot.c3temp
                                                         )#end ifelse
                                          )#end ifelse
                           )#end ifelse
pft$rr.decay.ecold = ifelse( test = pft$pathway == 4
                           , yes  = rr.decay.ecold.c4
                           , no   = rr.decay.ecold.c3
                           )#end ifelse
pft$rr.decay.ehot  = ifelse( test = pft$pathway == 4
                           , yes  = rr.decay.ehot.c4
                           , no   = rr.decay.ehot.c3
                           )#end ifelse
pft$m              = ifelse( test = pft$pathway == 4
                           , yes  = mphoto.c4
                           , no   = ifelse( test = pft$tropical & pft$conifer
                                          , yes  = mphoto.aa
                                          , no   = ifelse( test = pft$tropical
                                                         , yes  = mphoto.c3
                                                         , no   = mphoto.c3 * 6.3949 / 8.0
                                                         )#end ifelse
                                          )#end ifelse
                           )#end ifelse
pft$alpha          = ifelse( test = pft$pathway == 4
                           , yes  = alpha.c4
                           , no   = alpha.c3
                           )#end ifelse
pft$effarea.transp = ifelse( test = pft$conifer
                           , yes  = 2.0
                           , no   = 1.0
                           )#end ifelse
pft$b              = ifelse( test = pft$pathway == 4
                           , yes  = b.c4
                           , no   = ifelse( test = pft$tropical & pft$conifer
                                          , yes  = b.aa
                                          , no   = ifelse( test = pft$tropical
                                                         , yes  = b.c3
                                                         , no   = ifelse( test = pft$conifer
                                                                        , yes  = 1000.
                                                                        , no   = 20000.
                                                                        )#end ifelse
                                                         )#end ifelse
                                          )#end ifelse
                           )#end ifelse
pft$curvpar        = ifelse( test = pft$pathway == 4
                           , yes  = curvpar.c4
                           , no   = curvpar.c3
                           )#end ifelse
pft$phi.psII       = ifelse( test = pft$pathway == 4
                           , yes  = phi.psII.c4
                           , no   = phi.psII.c3
                           )#end ifelse
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#      Fill in the electron transport and respiration at reference temperature.  We must   #
# account for differences in Q10 to obtain the correct ratio at 15degC.                    #
#------------------------------------------------------------------------------------------#
pft$jm0            = ifelse( test = pft$pathway == 4
                           , yes  = jmfact.c4 * pft$vm0 * pft$vm.q10 / pft$jm.q10
                           , no   = jmfact.c3 * pft$vm0 * pft$vm.q10 / pft$jm.q10
                           )#end ifelse
pft$lr0            = ifelse( test = pft$pathway == 4
                           , yes  = gamma.c4 * pft$vm0 * pft$vm.q10 / pft$lr.q10
                           , no   = gamma.c3 * pft$vm0 * pft$vm.q10 / pft$lr.q10
                           )#end ifelse
pft$rr0            = 0.2455212 + 0. * pft$vm0
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Fill in the triose phosphate utilisation rate.  The C4 value is left undefined      #
# becausei its never used.  Instead, we use the PEP carboxylase-limited rate for them.     #
#------------------------------------------------------------------------------------------#
pft$tpm0           = ifelse( test = pft$pathway == 4
                           , yes  = 0.
                           , no   = tpmfact.c3 * pft$vm0
                           )#end ifelse
#------------------------------------------------------------------------------------------#


#----- Define reference height and coefficients for tropical allometry. -------------------#
tropical = is.na(pft$b1Ht)
if (iallom %in% c(0,1)){
   pft$b1Ht[tropical] = 0.37 * log(10)
   pft$b2Ht[tropical] = 0.64
}else if (iallom %in% c(2)){
   #---------------------------------------------------------------------------------------#
   #     Use the allometry proposed by:                                                    #
   #                                                                                       #
   # Poorter, L., L. Bongers, F. Bongers, 2006: Architecture of 54 moist-forest tree       #
   #    species: traits, trade-offs, and functional groups.  Ecology, 87, 1289-1301.       #
   #                                                                                       #
   #---------------------------------------------------------------------------------------#
   pft$b1Ht[tropical] = 0.0352
   pft$b2Ht[tropical] = 0.694
   #---------------------------------------------------------------------------------------#
}else if (iallom %in% c(3)){
   #---------------------------------------------------------------------------------------#
   #     Allometric equation based on the fitted curve using the Sustainable Landscapes    #
   # data set (L16) and the size- and site-dependent stratified sampling and aggregation   #
   # (J17).  This relationship is fitted using Standardised Major Axis (SMA) so the same   #
   # parameters can be used for y=f(x) and x=f(y).  This is particularly useful when       #
   # initialising the model with airborne lidar data.  Because it would be extremely       #
   # cumbersome to derive a SMA-based regression based on Weibull function, we use a log-  #
   # linear relationship.  The maximum height is based on the 99% quantile of all trees    #
   # measured by the SL team.                                                              #
   #                                                                                       #
   # References:                                                                           #
   #                                                                                       #
   # Jucker T, Caspersen J, Chave J, Antin C, Barbier N, Bongers F, Dalponte M,            #
   #    van Ewijk KY, Forrester DI, Haeni M et al. 2017. Allometric equations for          #
   #    integrating remote sensing imagery into forest monitoring programmes.              #
   #    Glob. Change Biol., 23(1):177-190. doi:10.1111/gcb.13388 (J17).                    #
   #                                                                                       #
   # Longo M, Keller M, dos-Santos MN, Leitold V, Pinage ER, Baccini A, Saatchi S,         #
   #    Nogueira EM, Batistella M , Morton DC. 2016. Aboveground biomass variability       #
   #    across intact and degraded forests in the Brazilian Amazon.  Global Biogeochem.    #
   #    Cycles, 30(11):1639-1660. doi:10.1002/2016GB005465 (L16).                          #
   #---------------------------------------------------------------------------------------#
   pft$b1Ht[tropical] = 1.139963 #   0.0400   # 0.0482
   pft$b2Ht[tropical] = 0.564899 #   0.8136   # 0.8307
   #---------------------------------------------------------------------------------------#
}else if (iallom %in% c(4)){
   #---------------------------------------------------------------------------------------#
   #     Allometric equation based on C14 and F15.                                         #
   #                                                                                       #
   # References:                                                                           #
   #                                                                                       #
   #   Chave, J, Rejou-Mechain M, Burquez A, Chidumayo E, Colgan MS, Delitti WB, Duque A,  #
   #      Eid T, Fearnside PM, Goodman RC et al. 2014. Improved allometric models to       #
   #      estimate the aboveground biomass of tropical trees. Glob. Change Biol., 20(10),  #
   #      3177-3190. doi:10.1111/gcb.12629 (C14).                                          #
   #                                                                                       #
   #   Falster DS, Duursma RA, Ishihara MI, Barneche DR, FitzJohn RG, Vahammar A, Aiba M,  #
   #      Ando M, Anten N, Aspinwall MJ et al. 2015. BAAD: a biomass and allometry         #
   #      database for woody plants. Ecology, 96 (5):1445-1445. doi:10.1890/14-1889.1      #
   #      (F15).                                                                           #
   #                                                                                       #
   #---------------------------------------------------------------------------------------#
   pft$b1Ht[tropical] = c14f15.ht.xx[1] + c14f15.ht.xx[2] * log(pft$rho[tropical])
   pft$b2Ht[tropical] = c14f15.ht.xx[3]
   #---------------------------------------------------------------------------------------#
}#end if
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Define minimum and maximum height based on life form and allometry.                  #
#------------------------------------------------------------------------------------------#
pft$hgt.min = ifelse( test = pft$tropical
                    , yes  = ifelse(test=pft$grass,yes=0.50,no=0.50           )
                    , no   = ifelse(test=pft$grass,yes=0.15,no=0.2+pft$hgt.ref)
                    )#end ifelse
pft$hgt.max = ifelse( test = pft$tropical
                    , yes  = ifelse(test=pft$grass,yes=1.5         ,no=hgt.max.trop )
                    , no   = ifelse(test=pft$grass,yes=.95*pft$b1Ht,no=.999*pft$b1Ht)
                    )#end ifelse
pft$hgt.show = pft$hgt.min + 0.015 * (pft$hgt.max - pft$hgt.min)
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Set qsw according to the allometry and the PFT.                                      #
#------------------------------------------------------------------------------------------#
pft$qsw = pft$SLA / sapwood.ratio.orig
for (ipft in sequence(npft+1)){
   #---- Check PFT and allometry. ---------------------------------------------------------#
   if (pft$tropical[ipft] && pft$conifer[ipft] && iallom %in% c(3)){
      pft$qsw[ipft] = pft$SLA[ipft] * pft$rho[ipft] / sapwood.factor["aa"]
   }else if (pft$tropical[ipft] && pft$grass[ipft] && iallom %in% c(3)){
      pft$qsw[ipft] = 1.0e-5
   }else if (pft$tropical[ipft] && iallom %in% c(3)){
      pft$qsw[ipft] = pft$SLA[ipft] * pft$rho[ipft] / sapwood.factor["bl"]
   }else{
      pft$qsw[ipft] = pft$SLA[ipft] / sapwood.ratio.orig
   }#end if (pft$tropical[ipft] && is.finite(pft$rho[ipft]) && iallom %in% c(3,4))
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
pft$qrhob        = rep(NA_real_,times=npft+1)
pft$qwatdry.leaf = ifelse(pft$tropical,1.85,2.50)
pft$qwatdry.wood = rep(NA_real_,times=npft+1)
pft$qwatdry.bark = rep(NA_real_,times=npft+1)
pft$c.leaf.dry   = rep(3218.,times=npft+1)
pft$c.wood.dry   = rep(103.1+3.867*spht.tref,times=npft+1)
pft$c.bark.dry   = pft$c.wood.dry
pft$brf.wd       = ifelse(pft$grass,0.0,0.16)
for (ipft in sequence(npft+1)){
   usedef = with(pft,grass[ipft] || liana[ipft] || conifer[ipft] || (! tropical[ipft]))
   usedef = usedef || (iecon == 0)
   if (usedef){
      #----- Default values. --------------------------------------------------------------#
      pft$qrhob       [ipft] = 0.49/0.61
      pft$qwatdry.wood[ipft] = 0.7
      pft$qwatdry.bark[ipft] = 0.7
      #------------------------------------------------------------------------------------#
   }else{
      #----- Wood-density dependent values. -----------------------------------------------#
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
#    136-146, Feb 2016. doi:10.1111/2041-210X.12525. (F15)                                 #
#------------------------------------------------------------------------------------------#
pft$b1Xs  = 0.315769481
pft$b1Xb  = 0.
pft$qbark = 0.
for (ipft in sequence(npft+1)){
   skip = pft$grass[ipft] || pft$liana[ipft] || (! pft$tropical[ipft]) || ( iallom != 3 )
   if (! skip){
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
for (ipft in sequence(npft+1)){
   if (pft$tropical[ipft]){
      if (iallom %in% c(0,1,3,4)){
         pft$dbh.min  [ipft] = exp((log(pft$hgt.min[ipft])-pft$b1Ht[ipft])/pft$b2Ht[ipft])
         pft$dbh.crit [ipft] = exp((log(pft$hgt.max[ipft])-pft$b1Ht[ipft])/pft$b2Ht[ipft])
      }else if (iallom %in% c(2)){
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
for (ipft in sequence(npft+1)){
   if (pft$tropical[ipft]){
      #------------------------------------------------------------------------------------#
      #      Fill in the structural biomass parameters.                                    #
      #------------------------------------------------------------------------------------#
      if (pft$liana[ipft]){
         #----- Lianas.  Ignore iallom. ---------------------------------------------------#
         pft$b1Bs.small[ipft] = 0.2749
         pft$b2Bs.small[ipft] = 2.69373
         pft$b1Bs.large[ipft] = pft$b1Bs.small[ipft]
         pft$b2Bs.large[ipft] = pft$b2Bs.small[ipft]
         #---------------------------------------------------------------------------------#
      }else if (iallom %in% c(0)){
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
      }else if (iallom %in% c(3,4)){
         if (pft$grass[ipft]){
            c14f15.bs.xx = c14f15.bs.gr
         }else if (pft$savannah[ipft]){
            c14f15.bs.xx = c14f15.bs.sv
         }else{
            c14f15.bs.xx = c14f15.bs.tf
         }#end if (pft$grass[ipft])
         #---- Based on Chave/Lescure based allometry, with single coefficients. ----------#
         pft$b1Bs.small[ipft] = c14f15.bs.xx[1] * pft$rho[ipft]^c14f15.bs.xx[2]
         pft$b2Bs.small[ipft] = c14f15.bs.xx[2]
         pft$b1Bs.large[ipft] = pft$b1Bs.small[ipft]
         pft$b2Bs.large[ipft] = pft$b2Bs.small[ipft]
         #---------------------------------------------------------------------------------#
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
      }else if (iallom %in% c(3,4)){
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
         #     The functional form uses both DBH and Height.                               #
         #                                                                                 #
         #                       CA = b1Ca * (DBH^2 * Hgt)^b2Ca                            #
         #                       m2            cm      m                                   #
         #                                                                                 #
         # Total number of trees: 17072                                                    #
         # b1Ca    = 0.370 (95% CI: [0.346;0.398])                                         #
         # b2Ca    = 0.464 (95% CI: [0.457;0.472])                                         #
         # R2      = 0.521                                                                 #
         # RMSE    = 29.78                                                                 #
         #---------------------------------------------------------------------------------#
         pft$b1Ca[ipft] = 0.370
         pft$b2Ca[ipft] = 0.464
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Replace b1Cl/b2Cl coefficients by those calculated by:                         #
      #                                                                                    #
      #    Longo, M. et al. Carbon Debt and Recovery time of degraded forests in           #
      #       the Amazon, in prep.                                                         #
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
      if (iallom %in% c(3,4) && (! pft$grass[ipft])){
         pft$b1Cl[ipft] = 0.29754
         pft$b2Cl[ipft] = 1.0324
      }#end if
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Fill in the leaf biomass parameters.                                          #
      #------------------------------------------------------------------------------------#
      if (pft$liana[ipft]){
         #----- Lianas.  Ignore iallom. ---------------------------------------------------#
         pft$b1Bl[ipft] = 0.0856
         pft$b2Bl[ipft] = 2.0
         #---------------------------------------------------------------------------------#
      }else if (iallom %in% c(0,1)){
         #---- ED-2.1 allometry. ----------------------------------------------------------#
         pft$b1Bl[ipft] = exp(a1 + c1l * pft$b1Ht[ipft] + d1l * log(pft$rho[ipft]))
         aux            = ( (a2l - a1) + pft$b1Ht[ipft] * (c2l - c1l) 
                           + log(pft$rho[ipft]) * (d2l - d1l)) * (1.0/log(dcrit))
         pft$b2Bl[ipft] = C2B * b2l + c2l * pft$b2Ht[ipft] + aux
      }else if(iallom %in% c(2)){
         pft$b1Bl[ipft] = C2B * exp(nleaf[1]) * pft$rho[ipft] / nleaf[3]
         pft$b2Bl[ipft] = nleaf[2]
      }else if(iallom %in% c(3)){
         #---------------------------------------------------------------------------------#
         #    Leaf allometry, use the individual leaf area allometry derived from the BAAD #
         # data base (F15), scaled by the specific leaf area.                              #
         #                                                                                 #
         # Reference:                                                                      #
         #                                                                                 #
         # Falster DS, Duursma RA, Ishihara MI, Barneche DR, FitzJohn RG, Vahammar A, Aiba #
         #      M, Ando M, Anten N, Aspinwall MJ. 2015. BAAD: a biomass and allometry      #
         #      database for woody plants. Ecology, 96 (5):1445-1445.                      #
         #      doi:10.1890/14-1889.1 (F15).                                               #
         #                                                                                 #
         #---------------------------------------------------------------------------------#
         pft$b1Bl[ipft] = c14f15.bl.xx[1] / pft$SLA[ipft] 
         pft$b2Bl[ipft] = c14f15.bl.xx[2]
         #---------------------------------------------------------------------------------#
      }else{
         #---------------------------------------------------------------------------------#
         #    Allometry based on the BAAD data based (F15).  We only used leaves from      #
         # wild tropical, note that b1Bl has the unit of m2 leaf under this scenario       #
         # and will be converted to leaf carbon using SLA in size2bl.                      #
         #---------------------------------------------------------------------------------#
         pft$b1Bl[ipft] = exp( c14f15.la.wd[1] + c14f15.la.wd[2] * log(pft$rho[ipft]))
         pft$b2Bl[ipft] = c14f15.la.wd[3]
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#
   }#end if (pft$tropical[ipft]
   #---------------------------------------------------------------------------------------#
}#end for (ipft in sequence(npft))
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#    Sapwood area allometry for plants                                                     #
#    - Grasses: 100% of sapwood area over basal area but sapwood area (not used)           #
#    - Temperate angiosperms: 30% of sapwood (oaks)                                        #
#    - Tropical angiosperms: allometry-dependent, see below.                               #
#    - Conifers: 30% of sapwood. (loblolly pine)                                           #
#    Users are welcome to update those numbuers based on literatures or measurements       #
#------------------------------------------------------------------------------------------#
if (iallom %in% c(4)){
   #---------------------------------------------------------------------------------------#
   #     Tropical angiosperms derived from:                                                #
   #                                                                                       #
   # Falster DS, Duursma RA, Ishihara MI, Barneche DR, FitzJohn RG, Vahammar A,            #
   #    Aiba M,  Ando M, Anten N, Aspinwall MJ et al. 2015. BAAD: a biomass and            #
   #    allometry database for woody plants. Ecology, 96 (5):1445-1445.                    #
   #    doi:10.1890/14-1889.1 (F15).                                                       #
   #                                                                                       #
   # Christoffersen BO, Gloor M, Fauset S, Fyllas NM, Galbraith DR, Baker TR, Kruijt B,    #
   #    Rowland L, Fisher RA, Binks OJ et al. 2016. Linking hydraulic traits to            #
   #    tropical forest function in a size- structured and trait-driven model (TFS         #
   #    v.1-Hydro). Geosci. Model Dev., 9: 4227-4255. doi:10.5194/gmd-9- 4227-2016         #
   #    (C16).                                                                             #
   #---------------------------------------------------------------------------------------#
   pft$b1SA = ifelse( test = pft$grass
                    , yes  = 1.0
                    , no   = ifelse( test = pft$conifer
                                   , yes  = 0.30
                                   , no   = ifelse(test=pft$tropical,yes=0.6572,no=0.3000)
                                   )#end ifelse
                    )#end ifelse
   pft$b2SA = ifelse( test = pft$tropical & (! pft$grass)
                    , yes  = 1.8530
                    , no   = 2.0000
                    )#end ifelse
   #---------------------------------------------------------------------------------------#
}else{
   #---------------------------------------------------------------------------------------#
   #     Tropical angiosperms derived from:                                                #
   #                                                                                       #
   # Meinzer, F. C., G. Goldstein, and J. L. Andrade (2001), Regulation of water flux      #
   #    through tropical forest canopy trees: Do universal rules apply?, Tree Physiol.,    #
   #    21(1), 19-26, doi:10.1093/treephys/21.1.19 (M01).                                  #
   #---------------------------------------------------------------------------------------#
   pft$b1SA = ifelse( test = pft$grass
                    , yes  = 1.0
                    , no   = ifelse( test = pft$conifer
                                   , yes  = 0.30
                                   , no   = ifelse(test=pft$tropical,yes=1.5820,no=0.3000)
                                   )#end ifelse
                    )#end ifelse
   pft$b2SA = ifelse( test = pft$tropical & (! pft$grass)
                    , yes  = 1.7640
                    , no   = 2.0000
                    )#end ifelse
   #---------------------------------------------------------------------------------------#
}#end if (iallom %in% c(4))
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
for (ipft in sequence(npft+1)){
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
pft$b1Rd     = rep(NA,times=npft+1)
pft$b2Rd     = rep(NA,times=npft+1)
if (iallom %in% c(0)){
   #----- Original ED-2.1 scheme, based on standing volume. -------------------------------#
   for (ipft in sequence(npft+1)){
      if (pft$grass[ipft]){
         #----- Grasses have fixed rooting depth (70 cm). ---------------------------------#
         pft$b2Rd   [ipft] = 0.0
         pft$b1Rd   [ipft] = -0.700
         #---------------------------------------------------------------------------------#
      }else{
         #----- Volume-based allometry for trees. -----------------------------------------#
         pft$b2Rd   [ipft] = 0.277
         pft$b1Rd   [ipft] = -exp(0.545 *log(10.)) * (0.65*pi*0.11*0.11)^pft$b2Rd[ipft]
         #---------------------------------------------------------------------------------#
      }#end if (pft$grass[ipft])
      #------------------------------------------------------------------------------------#
   }#end for (ipft %in% sequence(npft))
   #---------------------------------------------------------------------------------------#
}else if (iallom %in% c(1,2)){
   #----- Simple allometry (0.5 m for seedlings, 5.0m for 35-m trees. ---------------------#
   pft$b1Rd  [sequence(npft+1)] = -1.1140580
   pft$b2Rd  [sequence(npft+1)] =  0.4223014
   #---------------------------------------------------------------------------------------#
}else if (iallom %in% c(3)){
   #----- Experimental allometry combining Xu's data from Panama and Brum et al. (2018). --#
   pft$b1Rd = ifelse( test = pft$tropical & (! pft$liana)
                    , yes  = ifelse(test=pft$grass,yes=-2.572,no=-0.947)
                    , no   = -1.1140580
                    )#end ifelse
   pft$b2Rd = ifelse( test = pft$tropical & (! pft$liana)
                    , yes  = ifelse(test=pft$grass,yes=0.246,no=0.148)
                    , no   = 0.4223014
                    )#end ifelse
   #---------------------------------------------------------------------------------------#
}else if (iallom %in% c(3)){
   #----- Test allometry based on excavation data in Panama, using height as predictor. ---#
   pft$b1Rd  [sequence(npft+1)] = -0.609
   pft$b2Rd  [sequence(npft+1)] =  0.580
   #---------------------------------------------------------------------------------------#
}#end if
#------------------------------------------------------------------------------------------#






#------------------------------------------------------------------------------------------#
#     BDead -> DBH parameters.  These are calculated a posteriori because they are just    #
# the inverse of size->BDead allometry.                                                    #
#------------------------------------------------------------------------------------------#
pft$d2DBH.small = ifelse( test = pft$tropical & (! pft$liana) & (iallom %in% c(3,4))
                        , yes  = 1.  / ( ( 2. + pft$b2Ht ) * pft$b2Bs.small )
                        , no   = 1.  / pft$b2Bs.small
                        )#end ifelse
pft$d1DBH.small = ifelse( test = pft$tropical & (! pft$liana) & (iallom %in% c(3,4))
                        , yes  = ( C2B / ( pft$b1Bs.small * exp(pft$b1Ht*pft$b2Bs.small) ) )
                               ^ pft$d2DBH.small
                        , no   = ( C2B / pft$b1Bs.small ) ^ pft$d2DBH.small
                        )#end ifelse
pft$d2DBH.large = ifelse( test = pft$tropical & (! pft$liana) & (iallom %in% c(3,4))
                        , yes  = 1.  / ( ( 2. + pft$b2Ht ) * pft$b2Bs.large )
                        , no   = 1.  / pft$b2Bs.large
                        )#end ifelse
pft$d1DBH.large = ifelse( test = pft$tropical & (! pft$liana) & (iallom %in% c(3,4))
                        , yes  = ( C2B / ( pft$b1Bs.large * exp(pft$b1Ht*pft$b2Bs.large) ) )
                               ^ pft$d2DBH.large
                        , no   = ( C2B / pft$b1Bs.large ) ^ pft$d2DBH.large
                        )#end ifelse
pft$l2DBH       = ifelse( test = pft$tropical & (! pft$liana) & (iallom %in% c(3,4))
                        , yes  = 1.  / ( ( 2. + pft$b2Ht ) * pft$b2Bl )
                        , no   = 1.  / pft$b2Bl
                        )#end ifelse
pft$l1DBH       = ifelse( test = pft$tropical & (! pft$liana) & (iallom %in% c(3,4))
                        , yes  = ( C2B / ( pft$b1Bl * exp(pft$b1Ht*pft$b2Bl) ) )
                               ^ pft$l2DBH
                        , no   = ( C2B / pft$b1Bl ) ^ pft$l2DBH
                        )#end ifelse
pft$bdead.crit  = ifelse( test = pft$tropical & (! pft$liana) & (iallom %in% c(3,4))
                        , yes  = pft$b1Bs.small / C2B
                               * ( pft$dbh.crit * pft$dbh.crit * pft$hgt.max)
                               ^ pft$b2Bs.small
                        , no   = pft$b1Bs.small / C2B * pft$dbh.crit ^ pft$b2Bs.small
                        )#end ifelse
#------------------------------------------------------------------------------------------#






#------------------------------------------------------------------------------------------#
#    Alternative rooting depth parameter for tropical trees.  Use the allometric model     #
# to obtain the Effective Functional Rooting Depth based on B18.  We made a slight         #
# modification in their equation relating delta 18O and depth:                             #
#                                                                                          #
#    depth = exp(a + b * d18O^2)                                                           #
#                                                                                          #
# because it fits the data better than the original equation without the square, and it    #
# avoids extremely shallow soils for small trees.  We also use a heteroscedastic least     #
# squares, using the algorithm developed by L16.                                           #
#                                                                                          #
# References:                                                                              #
#                                                                                          #
# Brum M, Vadeboncoeur MA, Ivanov V, Asbjornsen H, Saleska S, Alves LF, Penha D,           #
#    Dias JD, Aragao LEOC, Barros F, Bittencourt P, Pereira L, Oliveira RS, 2018.          #
#    Hydrological niche segregation defines forest structure and drought tolerance         #
#    strategies in a seasonal Amazonian forest. J. Ecol., in press.                        #
#    doi:10.1111/1365-2745.13022 (B18).                                                    #
#                                                                                          #
# Longo M, Keller M, dos-Santos MN, Leitold V, Pinage ER, Baccini A, Saatchi S,            #
#    Nogueira EM, Batistella M , Morton DC. 2016. Aboveground biomass variability          #
#    across intact and degraded forests in the Brazilian Amazon.  Global Biogeochem.       #
#    Cycles, 30(11):1639-1660. doi:10.1002/2016GB005465 (L16).                             #
#------------------------------------------------------------------------------------------#
pft$d18O.ref = rep(x=-5.356,times=npft+1)
pft$b1d18O   = rep(x=0.0516,times=npft+1)
pft$b2d18O   = rep(x=1.0   ,times=npft+1)
pft$b1Efrd   = rep(x=-2.436,times=npft+1)
pft$b2Efrd   = rep(x=0.1822,times=npft+1)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Parameters to estimate wood area index.  Use the equation by:                       #
#                                                                                          #
# Hormann, G., S. Irrgan, H. Jochheim, M. Lukes, H. Meesenburg, J. Muller, B. Scheler,     #
#    J. Scherzer, G. Schuler, B. Schultze, B. Strohbach, F. Suckow, M. Wegehenkel,         #
#    and G. Wessolek.   Wasserhaushalt von waldokosystemen: methodenleitfaden zur          #
#    bestimmung der wasserhaushaltskomponenten auf level II-flachen. Technical note,       #
#    Bundesministerium fur Verbraucherschutz, Ernahrung und Landwirtschaft (BMVEL),        #
#    Bonn, Germany, 2003.                                                                  #
#    URL http://www.wasklim.de/download/Methodenband.pdf.                                  #
#                                                                                          #
# In the case of allometry 3 and tropical trees, use coefficients that scale WAI with      #
# dbh*dbh*height.                                                                          #
#------------------------------------------------------------------------------------------#
pft$b1WAI = ifelse( test = pft$grass
                  , yes  = 0.0
                  , no   = ifelse( test = pft$tropical & (! pft$liana)
                                                       & (iallom %in% c(3,4))
                                 , yes  = ifelse( test = pft$conifer
                                                , yes  = 0.01148449
                                                , no   = 0.00378399
                                                )#end ifelse
                                 , no   = ifelse( test = pft$conifer
                                                , yes  = 0.0553 * 0.5
                                                , no   = 0.0192 * 0.5
                                                )#end ifelse
                                 )#end ifelse
                  )#end ifelse
pft$b2WAI = ifelse( test = pft$grass
                  , yes  = 0.0
                  , no   = ifelse( test = pft$tropical & (! pft$liana)
                                                       & (iallom %in% c(3,4))
                                 , yes  = ifelse( test = pft$conifer
                                                , yes  = 0.77075160
                                                , no   = 0.81667933
                                                )#end ifelse
                                 , no   = ifelse( test = pft$conifer
                                                , yes  = 1.9769
                                                , no   = 2.0947
                                                )#end ifelse
                                 )#end ifelse
                  )#end ifelse
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#    Minimum bleaf and leaf area index that is resolvable.                                 #
#------------------------------------------------------------------------------------------#
pft$bleaf.min = size2bl( dbh  = pft$dbh.min[sequence(npft+1)]
                       , hgt  = pft$hgt.min[sequence(npft+1)]
                       , sla  = pft$SLA    [sequence(npft+1)]
                       , ipft = sequence(npft+1)
                       )#end size2bl
pft$lai.min   = onesixth * pft$init.dens * pft$bleaf.min * pft$SLA
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#    Minimum and maximum rooting depth.                                                    #
#------------------------------------------------------------------------------------------#
pft$rdp.min = size2rd( dbh  = pft$dbh.min[sequence(npft+1)]
                     , hgt  = pft$hgt.min[sequence(npft+1)]
                     , ipft = sequence(npft+1)
                     )#end size2rd
pft$rdp.max = size2rd( dbh  = pft$dbh.crit[sequence(npft+1)]
                     , hgt  = pft$hgt.max [sequence(npft+1)]
                     , ipft = sequence(npft+1)
                     )#end size2rd
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#    Find Vcmax at 25deg C, which is used by the plasticity terms.                         #
#------------------------------------------------------------------------------------------#
t25C           = 25. + t00
Vcmax25.nocorr = collatz(temp=t25C,refval=pft$vm0,q10=pft$vm.q10)
#------------------------------------------------------------------------------------------#
#    Find the temperature correction factors to shut down Vm at extreme temperatures.      #
# This should have minimal effect at 25degC, though.  It is unlikely, but to avoid float-  #
# ing point exceptions, we also check whether the temperature will make the exponential    #
# too large or too small.                                                                  #
#------------------------------------------------------------------------------------------#
#----- Low temperature. -------------------------------------------------------------------#
lnexplow   = pft$vm.decay.ecold  * (pft$vm.low.temp + t00 - t25C)
lnexplow   = pmax(lnexp.min,pmin(lnexp.max,lnexplow))
tlow.fun   = 1.0 +  exp(lnexplow)
#----- High temperature. ------------------------------------------------------------------#
lnexphigh = pft$vm.decay.ehot * (t25C - pft$vm.high.temp - t00)
lnexphigh = pmax(lnexp.min,pmin(lnexp.max,lnexphigh))
thigh.fun = 1.0 + exp(lnexphigh)
#----- Correct Vcmax. ---------------------------------------------------------------------#
pft$vcmax25 = Vcmax25.nocorr / (tlow.fun * thigh.fun)
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#    Find the trait plasticity terms.                                                      #
#                                                                                          #
# kplastic.vm0:  use the function from L10.  This is not the same fit, instead it was      #
#                found using robust standardised major axis to reduce leverage.            #
# kplastic.sla:  use an empirical function obtained using the data from KN2016, assuming   #
#                Beer's law to find the LAI that would cause the drop of light conditions  #
#                from 40 mol/m2/d to 3 mol/m2/d (their standard light conditions), and     #
#                assuming that changes in SLA  would follow a exponential decay/expansion, #
#                akin to Vm0.  Only tropical trees were used for the standardised major    #
#                axis regression, and we used a linear-log function because kplastic.SLA   #
#                may switch sign: in fact KN2016 observations had negative values for very #
#                high-canopy-level SLA values.                                             #
# eplastic.xxx:  these variables control the expansion/reduction exponent for leaf turn-   #
#                over rate when using trait plasticity.  Currently we start from Eq. 1 of  #
#                X17 (originally from K91), by substituting b with their SMA regression    #
#                against mass-based Vcmax25, and we fit a curve relating Aa with           #
#                Vcmax25_m, which was obtained by simulating ED-2 for average diurnal      #
#                cycles pooled from multiple tower sites in South America (mean diurnal    #
#                cycle of a month/year randomly selected), using random values of          #
#                Vcmax25_m within the range of values found in multiple trait data bases   #
#                (TRY, GLOPNET, RAINFOR, and NGEE-Tropics).                                #
# kplastic.ltor: the expansion/extinction rate derived from kplastic and eplastic terms    #
#                above.                                                                    #
#                                                                                          #
# References:                                                                              #
#                                                                                          #
# Keenan TF, Niinemets U. 2016. Global leaf trait estimates biased due to plasticity in    #
#    the shade. Nat. Plants 3:16201, Dec 2016. doi:10.1038/nplants.2016.201 (KN16).        #
#                                                                                          #
# Kikuzawa K. 1991. A cost-benefit analysis of leaf habit and leaf longevity of trees and  #
#    their geographical pattern. Am. Nat. 138(5): 1250-1263. doi:10.1086/285281 (K91).     #
#                                                                                          #
# Lloyd J, Patino S, Paiva RQ, Nardoto GB, Quesada CA, Santos AJB, Baker TR, Brand WA,     #
#    Hilke I, Gielmann H, Raessler M, Luizao FJ, Martinelli LA, Mercado LM. 2009.          #
#    Optimisation of photosynthetic carbon gain and within-canopy gradients of associated  #
#    foliar traits for Amazon forest trees. Biogeosciences 7(6): 1833-1859,                #
#    doi:10.5194/bg-7-1833-2010 (L10).                                                     #
#                                                                                          #
# Xu X, Medvigy D, Wright SJ, Kitajima K, Wu J, Albert LP, Martins GA, Saleska SR,         #
#    Pacala SW. 2017. Variations of leaf longevity in tropical moist forests predicted     #
#    by a trait-driven carbon optimality model. Ecol. Lett. 20(9): 1097-1106.              #
#    doi:10.1111/ele.12804 (X17).                                                          #
#------------------------------------------------------------------------------------------#
pft$kplastic.vm0  = - exp(-2.788 + 0.01439 * pft$vcmax25)
pft$kplastic.sla  = 0.462 - 0.1239 * log(pft$SLA)
pft$eplastic.vm0  = 0.5 * (-1.36 - 0.8399) + 0. * pft$kplastic.vm0
pft$eplastic.sla  = 0.5 * (-1.36 - 1.0000) + 0. * pft$kplastic.vm0
pft$kplastic.ltor = -1 * with(pft,kplastic.vm0*eplastic.vm0+kplastic.sla*eplastic.sla)
#------------------------------------------------------------------------------------------#



#----- Make it global. --------------------------------------------------------------------#
pft <<- pft
#------------------------------------------------------------------------------------------#




#==========================================================================================#
#==========================================================================================#
#     Global variables to define plant properties.                                         #
#------------------------------------------------------------------------------------------#
bltree = ( pft$tropical 
         & (! (pft$grass | pft$savannah | pft$conifer | pft$liana | (pft$key %in% "ALL")))
         )#end bltree
if (! "mypfts"     %in% ls()) mypfts     <<- which(bltree)
if (! "n.pfts"     %in% ls()) n.pfts     <<- length(mypfts)
if (! "pft.dens"   %in% ls()) pft.dens   <<- pft$rho[mypfts]
if (! "pft.breaks" %in% ls()) pft.breaks <<- c(-Inf,mid.points(pft.dens),Inf)
if (! "pft.cols"   %in% ls()) pft.cols   <<- pft$colour[mypfts]
if (! "pft.pch"    %in% ls()) pft.pch    <<- pch.stock(n=n.pfts)
if (! "pft.lwd"    %in% ls()) pft.lwd    <<- rep(1.5,n.pfts)
if (! "pft.names"  %in% ls()) pft.names  <<- pft$name[mypfts]
#------------------------------------------------------------------------------------------#
