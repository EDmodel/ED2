#------------------------------------------------------------------------------------------#
# Trigonometric constants                                                                  #
#------------------------------------------------------------------------------------------#
pii       <<- 1./pi             # 1/Pi                                          [      ---]
halfpi    <<- pi/2              # Pi/2                                          [      ---]
sqrtpii   <<- 0.564189583547756 # 1/(pi**0.5)                                   [      ---]
twopi     <<- pi* 2.            # 2 Pi                                          [      ---]
pio180    <<- pi/ 180.          # Pi/180 (deg -> rad)                           [      ---]
onerad    <<- 180. / pi         # 180/pi (rad -> deg)                           [      ---]
pi4       <<- pi  * 4.          # 4 Pi                                          [      ---]
pio4      <<- pi  / 4.          # Pi/4                                          [      ---]
pio6      <<- pi  / 6.          # Pi/6                                          [      ---]
pio6i     <<- 6.  / pi          # 6/Pi                                          [      ---]
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
# Algebraic shortcuts                                                                      #
#------------------------------------------------------------------------------------------#
srtwo     <<- sqrt(2.)          # Square root of 2.                            [       ---]
srthree   <<- sqrt(3.)          # Square root of 3.                            [       ---]
sqrt2o2   <<- 0.5 * srtwo       # Half of srtwo                                [       ---]
srtwoi    <<- 1./srtwo          # 1./ Square root of 2.                        [       ---]
srthreei  <<- 1./srthree        # 1./ Square root of 3.                        [       ---]
onethird  <<- 1./3.             # 1/3                                          [       ---]
twothirds <<- 2./3.             # 2/3                                          [       ---]
onesixth  <<- 1./6.             # 1/6                                          [       ---]
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
# Universal constants                                                                      #
#------------------------------------------------------------------------------------------#
stefan    <<- 5.6696e-8         # Stefan-Boltzmann constant                    [   W/m2/K4]
boltzmann <<- 1.3806503e-23     # Boltzmann constant                           [m2 kg/s2/K]
t00       <<- 273.15            # 0 degC                                       [      degC]
rmol      <<- 8.314510          # Molar gas constant                           [   J/mol/K]
volmol    <<- 0.022710980       # Molar volume at STP                          [        m3]
volmoll   <<- volmol*1.e3       # Molar volume at STP                          [         L]
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
# Molar masses and derived variables                                                       #
#------------------------------------------------------------------------------------------#
mmdry       <<- 0.02897        # Mean dry air molar mass                       [    kg/mol]
mmo2        <<- 0.03199880     # Mean water molar mass                         [    kg/mol]
mmh2o       <<- 0.01801505     # Mean water molar mass                         [    kg/mol]
mmco2       <<- 0.0440095      # Mean CO2 molar mass                           [    kg/mol]
mmdoc       <<- mmdry/mmco2    # mmdry/mmco2                                   [      ----]
mmcod       <<- mmco2/mmdry    # mmco2/mmdry                                   [      ----]
mmdry1000   <<- 1000.*mmdry    # Mean dry air molar mass                       [    kg/mol]
mmcod1em6   <<- mmcod * 1.e-6  # Convert ppm to kgCO2/kgair                    [      ----]
mmdryi      <<- 1./mmdry       # 1./mmdry                                      [    mol/kg]
mmco2i      <<- 1./mmco2       # 1./mmco2                                      [    mol/kg]
mmh2oi      <<- 1./mmh2o       # 1./mmh2o                                      [    mol/kg]
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
# Time conversion units                                                                    #
#------------------------------------------------------------------------------------------#
yr.day   <<- 365.2425         # # of days in a year                            [    day/yr]
yr.mon   <<- 12               # # of months in a year                          [    mon/yr]
day.sec  <<- 86400.           # # of seconds in a day                          [     s/day]
day.sec2 <<- day.sec^2        # # Square of day.sec                            [   s2/day2]
day.mon  <<- yr.day/yr.mon    # # of days in a month                           [   day/mon]
day.min  <<- 1440.            # # of minutes in a day                          [   min/day]
day.hr   <<- 24.              # # of hours in a day                            [    hr/day]
hr.sec   <<- 3600.            # # of seconds in an hour                        [      s/hr]
hr.min   <<- 60.              # # of minutes in an hour                        [    min/hr]
min.sec  <<- 60.              # # of seconds in a minute                       [     s/min]
yr.sec   <<- yr.day * day.sec # # of seconds in a year                         [      s/yr]
kg2g     <<- 1000.            # # of grams in a kilogram                       [      g/kg]
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
# General Earth properties                                                                 #
#------------------------------------------------------------------------------------------#
vonk      <<- 0.40            # Von Karman constant                             [      ---]
grav      <<- 9.80665         # Gravity acceleration                            [     m/s2]
gg        <<- .5 * grav       # Half of grav                                    [     m/s2]
erad      <<- 6370997.        # Earth radius                                    [        m]
spcon     <<- pio180*erad     # One degree of latitude                          [        m]
spconkm   <<- spcon*0.001     # One degree of latitude                          [       km]
eradi     <<- 1./erad         # Inverse of Earth radius                         [      1/m]
erad2     <<- erad*2          # Earth diameter                                  [        m]
ss60      <<- 1.8663          # Polar stereo conversion to 60 deg               [         ]
omega     <<- 7.292e-5        # Earth's rotation speed                          [    rad/s]
viscos    <<- .15e-4          # Viscosity coefficient                           [         ]
solar     <<- 1.3533e3        # Solar constant                                  [     W/m2]
p00       <<- 1.e5            # Reference pressure                              [       Pa]
prefsea   <<- 101325.         # Reference sea level pressure                    [       Pa]
p00i      <<- 1. / p00        # 1/p00                                           [     1/Pa]
o2.ref    <<- 0.209           # Nominal O2 concentration                        [  mol/mol]
capri     <<- -23.44 * pio180 # Tropic of Capricornium latitude                 [      rad]
shsummer  <<- -10             # Day of year of S.Hemisphere summer solstice     [      day]
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
# Reference for this block:                                                                #
# MU08 - Monteith, J. L., M. H. Unsworth, 2008. Principles of Environmental Physics,       #
#        third edition, Academic Press, Amsterdam, 418pp.  (Chapters 3 and 10).            #
#                                                                                          #
#     Air diffusion properties. These properties are temperature-dependent in reality,     #
# but for simplicity we assume them constants, using the value at 20C.                     #
#                                                                                          #
# Thermal diffusivity - Computed from equation on page 32 of MU08;                         #
# Kinematic viscosity - Computed from equation on page 32 of MU08;                         #
# Thermal expansion coefficient - determined by inverting the coefficient at equation      #
#                                 10.11 (MU08).                                            #
# These terms could be easily made function of temperature in the future if needed be.     #
#------------------------------------------------------------------------------------------#
th.diff.0  <<- 1.890e-5     # Air thermal diffusivity                           [     m2/s]
dth.diff   <<- 0.007        # Temperature dependency slop                       [      1/K]
kin.visc.0 <<- 1.330e-5     # Kinematic viscosity                               [     m2/s]
dkin.visc  <<- 0.007        # Temperature dependency slop                       [      1/K]
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
# Dry air properties                                                                       #
#------------------------------------------------------------------------------------------#
rdry    <<- rmol/mmdry       # Gas constant for dry air (Ra)                    [   J/kg/K]
rdryi   <<- mmdry/rmol       # 1./Gas constant for dry air (Ra)                 [   kg K/J]
cpdry   <<- 1004.            # Specific heat at constant pressure               [   J/kg/K]
cvdry   <<- 717.             # Specific heat at constant volume                 [   J/kg/K]
cpog    <<- cpdry /grav      # cp/g                                             [      m/K]
rocp    <<- rdry / cpdry     # Ra/cp                                            [     ----]
cpor    <<- cpdry / rdry     # Cp/Ra                                            [     ----]
rocv    <<- rdry / cvdry     # Ra/Cv                                            [     ----]
gocp    <<- grav / cpdry     # g/Cp, dry adiabatic lapse rate                   [      K/m]
gordry  <<- grav / rdry      # g/Ra                                             [      K/m]
cpdryi  <<- 1. / cpdry       # 1/Cp                                             [   kg K/J]
cpdryi4 <<- 4. * cpdryi      # 4/Cp                                             [   kg K/J]
p00k    <<- p00^(rdry/cpdry) # p0 ** (Ra/Cp)                                    [ Pa^0.286]
p00ki   <<- 1. / p00k        # p0 ** (-Ra/Cp)                                   [Pa^-0.286]
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Typical lapse rates for tropics and temperate zone.                                 #
#      The value for tropics came from:                                                    #
#      Gaffen, D. J.; B. D. Santer, J. S. Boyle, J. R. Christy, N. E. Graham, R. J. Ross,  #
#           2000: Multidecadal changes in the vertical temperature structure of the        #
#           tropical troposphere.  Science, 287, 1242-1245.                                #
#------------------------------------------------------------------------------------------#
lapse.trop <<- 5.5e-3    # g/Cp, dry adiabatic lapse rate                       [      K/m]
lapse.temp <<- 6.5e-3    # g/Cp, dry adiabatic lapse rate                       [      K/m]
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
# Water vapour properties                                                                  #
#------------------------------------------------------------------------------------------#
rh2o    <<- rmol/mmh2o   # Gas const. for water vapour (Rv)                     [   J/kg/K]
cph2o   <<- 1859.        # Heat capacity at const. pres.                        [   J/kg/K]
cph2oi  <<- 1. / cph2o   # Inverse of heat capacity                             [   kg K/J]
cvh2o   <<- cph2o-rh2o   # Heat capacity at const. volume                       [   J/kg/K]
gorh2o  <<- grav / rh2o  # g/Rv                                                 [      K/m]
ep      <<- mmh2o/mmdry  # or Ra/Rv, epsilon                                    [    kg/kg]
epi     <<- mmdry/mmh2o  # or Rv/Ra, 1/epsilon                                  [    kg/kg]
epim1   <<- epi-1.       # that 0.61 term of virtual temp.                      [    kg/kg]
rh2oocp <<- rh2o / cpdry # Rv/cp                                                [     ----]
toodry  <<- 1.e-8        # Minimum acceptable mixing ratio.                     [    kg/kg]
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
# Liquid water properties                                                                  #
#------------------------------------------------------------------------------------------#
wdns     <<- 1.000e3    # Liquid water density                                  [    kg/m3]
wdnsi    <<- 1./wdns    # Inverse of liquid water density                       [    m3/kg]
cliq     <<- 4.186e3    # Liquid water specific heat (Cl)                       [   J/kg/K]
cliqi    <<- 1./cliq    # Inverse of water heat capacity                        [   kg K/J]
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
# Ice properties                                                                           #
#------------------------------------------------------------------------------------------#
idns     <<- 9.167e2      # "Hard" ice density                                  [    kg/m3]
idnsi    <<- 1./idns      # Inverse of ice density                              [    m3/kg]
cice     <<- 2.093e3      # Ice specific heat (Ci)                              [   J/kg/K]
cicei    <<- 1. / cice    # Inverse of ice heat capacity                        [   kg K/J]
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
# Phase change properties                                                                  #
#------------------------------------------------------------------------------------------#
t3ple     <<- 273.16         # Water triple point temp. (T3)                    [        K]
t3plei    <<- 1./t3ple       # 1./T3                                            [      1/K]
es3ple    <<- 611.65685464   # Vapour pressure at T3 (es3)                      [       Pa]
es3plei   <<- 1./es3ple      # 1./es3                                           [     1/Pa]
epes3ple  <<- ep * es3ple    # epsilon * es3                                    [ Pa kg/kg]
rh2ot3ple <<- rh2o * t3ple   # Rv * T3                                          [     J/kg]
alli      <<- 3.34e5         # Lat. heat - fusion       (Lf)                    [     J/kg]
alvl3     <<- 2.50e6         # Lat. heat - vaporisation (Lv) at T = T3          [     J/kg]
alvi3     <<- alli + alvl3   # Lat. heat - sublimation  (Ls) at T = T3          [     J/kg]
allii     <<- 1.   / alli    # 1./Lf                                            [     kg/J]
aklv      <<- alvl3 / cpdry  # Lv/Cp                                            [        K]
akiv      <<- alvi3 / cpdry  # Ls/Cp                                            [        K]
lvordry   <<- alvl3 / rdry   # Lv/Ra                                            [        K]
lvorvap   <<- alvl3 / rh2o   # Lv/Rv                                            [        K]
lsorvap   <<- alvi3 / rh2o   # Ls/Rv                                            [        K]
lvt3ple   <<- alvl3 * t3ple  # Lv * T3                                          [   K J/kg]
lst3ple   <<- alvi3 * t3ple  # Ls * T3                                          [   K J/kg]
uiicet3   <<- cice  * t3ple  # internal energy at triple point, only ice        [     J/kg]
uiliqt3   <<- uiicet3 + alli # internal energy at triple point, only liq.       [     J/kg]
dcpvl     <<- cph2o - cliq   # Difference of sp. heat                           [   J/kg/K]
dcpvi     <<- cph2o - cice   # Difference of sp. heat                           [   J/kg/K]
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     The following variables are useful when defining the derivatives of theta_il.  They  #
# correspond to L?(T) - L?' T.                                                             #
#------------------------------------------------------------------------------------------#
del.alvl3 <<- alvl3 - dcpvl * t3ple
del.alvi3 <<- alvi3 - dcpvi * t3ple
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#    Tsupercool are defined as temperatures of supercooled liquid water (water vapour)     #
# that will cause the internal energy (enthalpy) to be the same as ice at 0K.  It can be   #
# used as an offset for temperature when defining internal energy (enthalpy).  The next    #
# two methods of defining the internal energy for the liquid part:                         #
#                                                                                          #
#   Uliq = Mliq [ Cice T3 + Cliq (T - T3) + Lf]                                            #
#   Uliq = Mliq Cliq (T - Tsupercool_liq)                                                  #
#                                                                                          #
#   H    = Mliq [ Cice T3 + Cliq (Ts - T3) + Lv3 + (Cpv - Cliq) (Ts-T3) + Cpv (T-T3) ]     #
#   H    = Mliq Cpv (T - Tsupercool_vap) ]                                                 #
#                                                                                          #
#     You may be asking yourself why would we have the ice term in the internal energy     #
# definition. The reason is that we can think that internal energy is the amount of energy #
# a parcel received to leave the 0K state to reach the current state (or if you prefer the #
# inverse way, Uliq is the amount of energy the parcel would need to lose to become solid  #
# at 0K.)                                                                                  #
#------------------------------------------------------------------------------------------#
tsupercool.liq <<- t3ple - (uiicet3+alli) * cliqi
tsupercool.vap <<- cph2oi * ( (cph2o - cice) * t3ple - alvi3 )
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Minimum temperature for computing the condensation effect of temperature on          #
# theta_il, thetae_iv, and associates. Below this temperature, assuming the latent         #
# heats as constants becomes a really bad assumption. See :                                #
#                                                                                          #
# Tripoli, J. T.; and Cotton, W.R., 1981: The use of ice-liquid water potential temper-    #
#    ature as a thermodynamic variable in deep atmospheric models. Mon. Wea. Rev.,         #
#    v. 109, 1094-1102.                                                                    #
#------------------------------------------------------------------------------------------#
ttripoli  <<- 253.           # "Tripoli-Cotton" temp. (Ttr)                     [        K]
htripoli  <<- cpdry*ttripoli # Sensible enthalpy at T=Ttr                       [     J/kg]
htripolii <<- 1./htripoli    # 1./htripoli                                      [     kg/J]
#------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------#
# Unit conversion, it must be defined locally even for coupled runs.                       #
#------------------------------------------------------------------------------------------#
mol.2.umol     <<- 1.e6                 # mol         => umol
umol.2.mol     <<- 1.e-6                # umol        => mol
umol.2.kgC     <<- 1.20107e-8           # umol(CO2)   => kg(C)
Watts.2.Ein    <<- 4.6e-6               # W/m2        => mol/m2/s
Ein.2.Watts    <<- 1./Watts.2.Ein       # mol/m2/s    => W/m2
kgC.2.umol     <<- 1. / umol.2.kgC      # kg(C)       => umol(CO2)
kgom2.2.tonoha <<- 10.                  # kg(C)/m2    => ton(C)/ha
tonoha.2.kgom2 <<- 0.1                  # ton(C)/ha   => kg(C)/m2
umols.2.kgCyr  <<- umol.2.kgC * yr.sec  # umol(CO2)/s => kg(C)/yr
kgCyr.2.umols  <<- 1. / umols.2.kgCyr   # kg(C)/yr    => umol(CO2)/s
kgCday.2.umols <<- kgC.2.umol / day.sec # kg(C)/day   => umol(CO2)/s
Torr.2.Pa      <<- prefsea / 760.       # Torr        => Pa
Pa.2.Torr      <<- 1. / Torr.2.Pa       # Pa          => Torr
kt.2.mos       <<- 1852 / hr.sec        # knots       => m/s
mos.2.kt       <<- 1. / kt.2.mos        # m/s         => knots
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Guess for ET from :                                                                  #
#   Malhi, Y., et al. 2009: Exploring the likelihood and mechanism of a climate-change-    #
#      -induced dieback of the Amazon rainforest. Proc. Natl. Ac. Sci., 106, 20610--20615. #
#------------------------------------------------------------------------------------------#
et.malhi <<- 100  # [mm/month]
#------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------#
#     These are the lower and upper bounds in which we compute exponentials.  This is to   #
# avoid overflows and/or underflows when we compute exponentials.                          #
#------------------------------------------------------------------------------------------#
lnexp.min <<- -38.
lnexp.max <<-  38.
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     These are the just default huge and tiny numbers that are not the actual huge or     #
# tiny values from Fortran intrinsic functions, so if you do any numerical operations you  #
# will still be fine.                                                                      #
#------------------------------------------------------------------------------------------#
huge.num  <<- 1.e+19
tiny.num  <<- 1.e-19
#------------------------------------------------------------------------------------------#
