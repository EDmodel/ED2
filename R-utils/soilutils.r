#==========================================================================================#
#==========================================================================================#
#    General settings for soil parametrisation.                                            #
#------------------------------------------------------------------------------------------#
#----- Soil hydraulic scheme. -------------------------------------------------------------#
if (! "isoil.hydro" %in% ls()) isoil.hydro = 0
#----- Minimum hydraulic conductivity. ----------------------------------------------------#
hydcond.min = 1.e-5 / wdns / day.sec
#------------------------------------------------------------------------------------------#
#==========================================================================================#
#==========================================================================================#


#==========================================================================================#
#==========================================================================================#
#      This variable has the "edges" of all soil types.                                    #
#------------------------------------------------------------------------------------------#
stext.lines       = list()
stext.lines[[ 1]] = list(sand=c(0.900,0.850),clay=c(0.100,0.000))
stext.lines[[ 2]] = list(sand=c(0.850,0.700),clay=c(0.150,0.000))
stext.lines[[ 3]] = list(sand=c(0.800,0.525),clay=c(0.200,0.200))
stext.lines[[ 4]] = list(sand=c(0.520,0.525),clay=c(0.200,0.075))
stext.lines[[ 5]] = list(sand=c(0.425,0.525),clay=c(0.075,0.075))
stext.lines[[ 6]] = list(sand=c(0.225,0.500),clay=c(0.275,0.000))
stext.lines[[ 7]] = list(sand=c(0.200,0.075),clay=c(0.000,0.125))
stext.lines[[ 8]] = list(sand=c(0.075,0.000),clay=c(0.125,0.125))
stext.lines[[ 9]] = list(sand=c(0.525,0.450),clay=c(0.200,0.275))
stext.lines[[10]] = list(sand=c(0.450,0.000),clay=c(0.275,0.275))
stext.lines[[11]] = list(sand=c(0.200,0.200),clay=c(0.275,0.400))
stext.lines[[12]] = list(sand=c(0.650,0.450),clay=c(0.350,0.350))
stext.lines[[13]] = list(sand=c(0.450,0.450),clay=c(0.275,0.550))
stext.lines[[14]] = list(sand=c(0.450,0.000),clay=c(0.400,0.400))
stext.lines[[15]] = list(sand=c(0.200,0.000),clay=c(0.400,0.600))
stext.lines[[16]] = list(sand=c(0.300,0.000),clay=c(0.400,0.700))
stext.lines[[17]] = list(sand=c(0.300,0.300),clay=c(0.400,0.700))
stext.lines[[18]] = list(sand=c(0.300,0.000),clay=c(0.700,0.700))
nstext.lines      = length(stext.lines)
for(n in 1:nstext.lines){
   stext.lines[[n]]$silt = pmax(0,pmin(1,1.-stext.lines[[n]]$sand-stext.lines[[n]]$clay))
}#end for
stext.lines   <<- stext.lines
nstext.lines  <<- nstext.lines
#==========================================================================================#
#==========================================================================================#




#==========================================================================================#
#==========================================================================================#
#      This variable has the "polygons" for all soil types.                                #
#------------------------------------------------------------------------------------------#
stext.polygon       = list()
stext.polygon[[ 1]] = list(sand = c(1.000,0.900,0.850)
                          ,clay = c(0.000,0.100,0.000)
                          )#end list
stext.polygon[[ 2]] = list(sand = c(0.900,0.850,0.700,0.850)
                          ,clay = c(0.100,0.150,0.000,0.000)
                          )#end list
stext.polygon[[ 3]] = list(sand = c(0.850,0.800,0.525,0.525,0.425,0.500,0.700)
                          ,clay = c(0.150,0.200,0.200,0.075,0.075,0.000,0.000)
                          )#end list
stext.polygon[[ 4]] = list(sand = c(0.500,0.225,0.000,0.000,0.075,0.200)
                          ,clay = c(0.000,0.275,0.275,0.125,0.125,0.000)
                          )#end list
stext.polygon[[ 5]] = list(sand = c(0.525,0.450,0.225,0.425,0.525)
                          ,clay = c(0.200,0.275,0.275,0.075,0.075)
                          )#end list
stext.polygon[[ 6]] = list(sand = c(0.800,0.650,0.450,0.450,0.525)
                          ,clay = c(0.200,0.350,0.350,0.275,0.200)
                          )#end list
stext.polygon[[ 7]] = list(sand = c(0.200,0.000,0.000,0.200)
                          ,clay = c(0.400,0.400,0.275,0.275)
                          )#end list
stext.polygon[[ 8]] = list(sand = c(0.450,0.200,0.200,0.450)
                          ,clay = c(0.400,0.400,0.275,0.275)
                          )#end list
stext.polygon[[ 9]] = list(sand = c(0.650,0.450,0.450)
                          ,clay = c(0.350,0.550,0.350)
                          )#end list
stext.polygon[[10]] = list(sand = c(0.200,0.000,0.000)
                          ,clay = c(0.400,0.600,0.400)
                          )#end list
stext.polygon[[11]] = list(sand = c(0.300,0.300,0.000)
                          ,clay = c(0.400,0.700,0.700)
                          )#end list
stext.polygon[[12]] = list(sand = c(NA,NA)
                          ,clay = c(NA,NA)
                          )#end list
stext.polygon[[13]] = list(sand = c(NA,NA)
                          ,clay = c(NA,NA)
                          )#end list
stext.polygon[[14]] = list(sand = c(0.200,0.075,0.000,0.000)
                          ,clay = c(0.000,0.125,0.125,0.000)
                          )#end list
stext.polygon[[15]] = list(sand = c(0.300,0.000,0.000)
                          ,clay = c(0.700,1.000,0.700)
                          )#end list
stext.polygon[[16]] = list(sand = c(0.450,0.300,0.300,0.450)
                          ,clay = c(0.550,0.700,0.400,0.400)
                          )#end list
stext.polygon[[17]] = list(sand = c(0.300,0.000,0.000,0.200)
                          ,clay = c(0.400,0.700,0.600,0.400)
                          )#end list
nstext.polygon      = length(stext.polygon)

for(n in sequence(nstext.polygon)){
   sand.now = stext.polygon[[n]]$sand
   clay.now = stext.polygon[[n]]$clay
   stext.polygon[[n]]$silt = pmax(0,pmin(1,1.-sand.now-clay.now))
}#end for
stext.polygon  <<- stext.polygon
nstext.polygon <<- nstext.polygon
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     This function finds the soil parameters.                                             #
#------------------------------------------------------------------------------------------#
soil.params <<- function(ntext,isoilflg,slxsand,slxclay,slsoc,slph,slcec,sldbd
                        ,out.dfr=FALSE){
   #----- Define some prescribed fractions. -----------------------------------------------#
   xsand.def = c( 0.920, 0.825, 0.660, 0.200, 0.410, 0.590
                , 0.100, 0.320, 0.520, 0.060, 0.200, 0.200
                , 0.333, 0.075, 0.100, 0.375, 0.125)
   xclay.def = c( 0.030, 0.060, 0.110, 0.160, 0.170, 0.270
                , 0.340, 0.340, 0.420, 0.470, 0.600, 0.200
                , 0.333, 0.050, 0.800, 0.525, 0.525)
                
   soil.name = c("Sand","Loamy sand","Sandy loam","Silt loam","Loam","Sandy clay loam"
                ,"Silty clay loam","Clayey loam","Sandy clay","Silty clay","Clay"
                ,"Peat","Bedrock","Silt","Heavy clay","Clayey sand","Clayey silt")
   soil.key  = c("Sa","LSa","SaL","SiL","L","SaCL","SiCL","CL","SaC","SiC","C","P","BR"
                ,"Si","CC","CSa","CSi")
   #----- Define some constants. ----------------------------------------------------------#
   fieldcp.K  = 0.1    # hydraulic conduct. at field capacity                     [ mm/day]
   slpotst.MPa = 5.e-4 # "saturated" soil approximation for van Genuchten         [    MPa]
   slpotth.MPa = 0.01  # soil-water potential: field capacity (Tomasella-Hodnett) [    MPa]
   slpotsr.MPa = 0.033 # soil-water potential: field capacity (Saxton-Rawls)      [    MPa]
   slpotcp.MPa = 3.1   # soil-water potential: air dry soil                       [    MPa]
   slpotwp.MPa = 1.5   # soil-water potential: wilting point                      [    MPa]
   slpotld.MPa = 0.75  # soil-water potential that plants start dropping leaves   [    MPa]
   slpotfr.MPa = 1.40  # soil-water potential that triggers fires                 [    MPa]
   theta.crit  = 0.11  # frac. soil moisture that plants start dropping leaves    [  m3/m3]
   sm.fire     = 0.08  # frac. soil moisture that triggers fires                  [  m3/m3]
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   # Soil heat capacity.  Didn't find silt values, using average between sand and clay     #
   #---------------------------------------------------------------------------------------#
   sand.hcap = 2.128e6
   clay.hcap = 2.385e6
   silt.hcap = .5 * (sand.hcap + clay.hcap)
   air.hcap  = 1212
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   # Soil heat capacity.  Didn't find silt values, using average between sand and clay     #
   #---------------------------------------------------------------------------------------#
   sand.cond = 8.80
   clay.cond = 2.92
   silt.cond = .5 * (sand.cond + clay.cond)
   air.cond  = 0.025
   h2o.cond  = 0.57
   #---------------------------------------------------------------------------------------#



   #----- Initialise the list with meaningless parameters. --------------------------------#
   mysoil = list( ntext      = NA_integer_
                , key        = NA_character_
                , name       = NA_character_
                , method     = NA_character_
                , xsand      = NA_real_
                , xsilt      = NA_real_
                , xclay      = NA_real_
                , slsoc      = NA_real_
                , slph       = NA_real_
                , slcec      = NA_real_
                , sldbd      = NA_real_
                , soilre     = 0.
                , soilcp     = 0.
                , soilwp     = 0.
                , soilfr     = 0.
                , soilld     = 0.
                , sfldcap    = 0.
                , soilbp     = 0.
                , slmsts     = 0.
                , soilpo     = 0.
                , slpotcp    = 0.
                , slpotwp    = 0.
                , slpotfr    = 0.
                , slpotld    = 0.
                , slpotfc    = 0.
                , slpotbp    = 0.
                , slpots     = 0.
                , slpotpo    = 0.
                , sltt       = 0.
                , slnm       = 0.
                , slbs       = 0.
                , slmm       = 0.
                , slmu       = 0.
                , malpha     = 0.
                , slcons     = 0.
                , fhydraul   = 0.
                , slcpd      = 0.
                , thcond0    = 0.
                , thcond1    = 0.
                , thcond2    = 0.
                , thcond3    = 0.
                , slas       = 0.
                , thcr.orig  = 0.
                , thcr.spot  = 0.
                , fire.orig  = 0.
                , fire.spot  = 0.
                )#end list
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find soil class and sand, silt, and clay fractions.                               #
   #---------------------------------------------------------------------------------------#
   if (missing(slxsand) || missing(slxclay) || missing(isoilflg)){
      mysoil$ntext = ntext
      mysoil$name  = soil.name[mysoil$ntext]
      mysoil$key   = soil.key [mysoil$ntext]
      mysoil$xsand = xsand.def[ntext]
      mysoil$xclay = xclay.def[ntext]
      mysoil$xsilt = 1. - mysoil$xsand - mysoil$xclay
   }else if (isoilflg == 2 & slxsand > 0 & slxclay > 0 & (slxsand + slxclay) < 1 ){
      mysoil$ntext = sclass(slxsand,slxclay)
      mysoil$name  = soil.name[mysoil$ntext]
      mysoil$key   = soil.key [mysoil$ntext]
      mysoil$xsand = slxsand
      mysoil$xclay = slxclay
      mysoil$xsilt = 1. - mysoil$xsand - mysoil$xclay
   }else{
      mysoil$ntext = ntext
      mysoil$name  = soil.name[mysoil$ntext]
      mysoil$key   = soil.key [mysoil$ntext]
      mysoil$xsand = xsand.def[ntext]
      mysoil$xclay = xclay.def[ntext]
      mysoil$xsilt = 1. - mysoil$xsand - mysoil$xclay
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Decide which method to use. -----------------------------------------------------#
   if ( mysoil$ntext %in% 13){
      mysoil$method = "BDRK"
   }else if ( (mysoil$ntext %in% c(12)) || isoil.hydro %in% c(0,1)){
      mysoil$method = "BC64"
   }else if (isoil.hydro %in% 2){
      mysoil$method = "vG80"
   }else{
      mysoil$method = "SR06"
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Additional soil parameters must be provided.                                     #
   #---------------------------------------------------------------------------------------#
   if (mysoil$method %in% "vG80"){
      #----- Make sure all required variables are provided. -------------------------------#
      fine = ! (missing(slsoc) || missing(slph) || missing(slcec) || missing(sldbd))
      if (fine){
         #----- Copy values from inputs. --------------------------------------------------#
         mysoil$slsoc = slsoc
         mysoil$slph  = slph
         mysoil$slcec = slcec
         mysoil$sldbd = sldbd
         #---------------------------------------------------------------------------------#
      }else{
         #----- Missing data, issue an error message. -------------------------------------#
         cat0("---------------------------------------------------------------------------")
         cat0(" Missing parameters for isoil.hydro = 2 (Hodnett and Tomasella 2002)! ")
         cat0(" Missing \"slsoc\" = ",missing(slsoc),".")
         cat0(" Missing \"slph\"  = ",missing(slph) ,".")
         cat0(" Missing \"slcec\" = ",missing(slcec),".")
         cat0(" Missing \"sldbd\" = ",missing(sldbd),".")
         cat0("---------------------------------------------------------------------------")
         stop(" Invalid settings for soil.params.")
         #---------------------------------------------------------------------------------#
      }#end if (fine)
      #------------------------------------------------------------------------------------#
   }else if (mysoil$method %in% "SR06"){
      #----- Make sure all required variables are provided. -------------------------------#
      fine = ! (missing(slsoc) || missing(slph) || missing(slcec) || missing(sldbd))
      if (fine){
         #----- Copy values from inputs. --------------------------------------------------#
         mysoil$slsoc = slsoc
         #---------------------------------------------------------------------------------#
      }else{
         #----- Missing data, issue an error message. -------------------------------------#
         cat0("----------------------------------------------------------------------")
         cat0(" Missing \"slsoc\" for isoil.hydro = 3 (Saxton and Rawls 2006)! ")
         cat0("----------------------------------------------------------------------")
         stop(" Invalid settings for soil.params.")
         #---------------------------------------------------------------------------------#
      }#end if (fine)
      #------------------------------------------------------------------------------------#

      #----- Variables are dummy but we honour the input variables if provided. -----------#
      if (! missing(slph )) mysoil$slph  = slph
      if (! missing(slcec)) mysoil$slcec = slcec
      if (! missing(sldbd)) mysoil$sldbd = sldbd
      #------------------------------------------------------------------------------------#
   }else{
      #----- Variables are dummy but we honour the input variables if provided. -----------#
      if (! missing(slsoc)) mysoil$slsoc = slsoc
      if (! missing(slph )) mysoil$slph  = slph
      if (! missing(slcec)) mysoil$slcec = slcec
      if (! missing(sldbd)) mysoil$sldbd = sldbd
      #------------------------------------------------------------------------------------#
   }#end if (mysoil$method %in% "vG1980")
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #       Set up primary properties.                                                      #
   #---------------------------------------------------------------------------------------#
   if (mysoil$ntext == 13){
      #----- Bedrock.  Most things are zero, because it is an impermeable soil. -----------#
      mysoil$slcpd     =  2130000.
      #------------------------------------------------------------------------------------#
   }else if (mysoil$ntext == 12){
      #------------------------------------------------------------------------------------#
      #      Peat.  High concentration of organic matter.  Mineral soil equations don't    #
      # apply here.                                                                        #
      #------------------------------------------------------------------------------------#
      mysoil$sltt    = 1.0
      mysoil$slbs    =  7.75
      mysoil$slnm    = 1. / mysoil$slbs
      mysoil$slmm    = 2. + mysoil$sltt + 2. * mysoil$slbs
      mysoil$slmu    = -1. / mysoil$slmm
      mysoil$slpots  = -0.356
      mysoil$slpotpo = mysoil$slpots
      mysoil$slpotbp = mysoil$slpots
      mysoil$slcons  =  2.357930e-6
      mysoil$slmsts  =  0.863
      mysoil$soilpo  =  mysoil$slmsts
      mysoil$soilbp  =  mysoil$slmsts
      mysoil$sfldcap =  0.535
      mysoil$slcpd   =  874000.
      mysoil$slpotfc = smoist2mpot(mysoil$sfldcap, mysoil)
      mysoil$soilre  = 0.
      mysoil$slpotwp = - slpotwp.MPa * 1.e6 / ( grav * wdns )
      mysoil$soilwp  = mpot2smoist(mysoil$slpotwp, mysoil)
      #------------------------------------------------------------------------------------#
   }else if (isoil.hydro == 0){
      #------------------------------------------------------------------------------------#
      #      Mineral soil.  Use the standard ED-2.2 equations.                             #
      #------------------------------------------------------------------------------------#
      mysoil$sltt    = 1.
      mysoil$slnm    = 1. / (3.10 + 15.7*mysoil$xclay - 0.3*mysoil$xsand)
      mysoil$slbs    = 1. / mysoil$slnm
      mysoil$slmm    = 2. + mysoil$sltt + 2. * mysoil$slbs
      mysoil$slmu    = -1. / mysoil$slmm
      mysoil$slpots  = -0.01 * (10.^(2.17 - 0.63*mysoil$xclay - 1.58*mysoil$xsand))
      mysoil$slpotpo = mysoil$slpots
      mysoil$slpotbp = mysoil$slpots
      mysoil$slcons  = ( (10.^(-0.60 + 1.26*mysoil$xsand - 0.64*mysoil$xclay)) 
                       * 0.0254/hr.sec )
      mysoil$slmsts  = (50.5 - 14.2*mysoil$xsand - 3.7*mysoil$xclay) / 100.
      mysoil$soilpo  =  mysoil$slmsts
      mysoil$soilbp  =  mysoil$slmsts
      mysoil$sfldcap = ( mysoil$slmsts 
                       * ( mysoil$slcons / (fieldcp.K /(wdns*day.sec))) ^ mysoil$slmu
                       )#end mysoil$sfldcap
      mysoil$slcpd   = ( (1. - mysoil$slmsts)
                       * ( mysoil$xsand * sand.hcap + mysoil$xsilt * silt.hcap
                         + mysoil$xclay * clay.hcap )
                       + 0.5 * (mysoil$slmsts - mysoil$soilcp) * air.hcap 
                       )#end mysoil$slcpd
      mysoil$slpotfc = smoist2mpot(mysoil$sfldcap, mysoil)
      mysoil$slpotwp = - slpotwp.MPa * 1.e6 / ( grav * wdns )
      mysoil$soilwp  = mpot2smoist(mysoil$slpotwp, mysoil)
   }else if (isoil.hydro == 1){
      #------------------------------------------------------------------------------------#
      #      Use Tomasella and Hodnett (1998).                                             #
      #------------------------------------------------------------------------------------#
      mysoil$sltt    = 0.5
      mysoil$slnm    = with( mysoil
                           , exp( - 1.197 - 0.417 * xsilt + 0.450 * xclay
                                  -  8.94 * xsilt * xclay +  10.0 * xsilt * xsilt * xclay
                                )#end exp
                           )#end with
      mysoil$slbs    = 1. / mysoil$slnm
      mysoil$slmm    = 2. + mysoil$sltt + 2. * mysoil$slbs
      mysoil$slmu    = -1. / mysoil$slmm
      mysoil$slpots  = with( mysoil
                           , -1. / grav * ( 0.285 + 7.33 * xsilt * xsilt 
                                          - 1.30 * xsilt * xclay 
                                          + 3.60 * xsilt * xsilt * xclay
                                          )#end
                           )#end with
      mysoil$slpotpo = mysoil$slpots
      mysoil$slpotbp = mysoil$slpots
      mysoil$slmsts  = with( mysoil
                           , 0.4061  + 0.165 * xsilt + 0.162 * xclay 
                           + 1.37e-3 * xsilt * xsilt + 1.80e-5 * xsilt * xsilt * xclay
                           )#end with
      mysoil$soilpo  =  mysoil$slmsts
      mysoil$soilbp  =  mysoil$slmsts
      mysoil$soilre  = with( mysoil
                           , 0. * slmsts + pmax( 0.01
                                               , - 0.02095 + 0.047 * xsilt + 0.431 * xclay
                                                 - 0.00827 * xsilt * xclay
                                               )#end pmax
                           )#end with
      mysoil$slpotfc = - slpotth.MPa * 1.e6 / ( grav * wdns )
      mysoil$sfldcap = mpot2smoist(mysoil$slpotfc, mysoil)

      mysoil$slcons  = ( (10.^(-0.60 + 1.26*mysoil$xsand - 0.64*mysoil$xclay)) 
                       * 0.0254/hr.sec )
      mysoil$slcpd   = ( (1. - mysoil$slmsts)
                       * ( mysoil$xsand * sand.hcap + mysoil$xsilt * silt.hcap
                         + mysoil$xclay * clay.hcap )
                       + 0.5 * (mysoil$slmsts - mysoil$soilcp) * air.hcap )
      mysoil$slpotwp = - slpotwp.MPa * 1.e6 / ( grav * wdns )
      mysoil$soilwp  = mpot2smoist(mysoil$slpotwp, mysoil)
      #------------------------------------------------------------------------------------#
   }else if (isoil.hydro == 2){
      #------------------------------------------------------------------------------------#
      #     Use Hodnett and Tomasella (2002).                                              #
      #------------------------------------------------------------------------------------#
      mysoil$sltt = -1.0
      #------------------------------------------------------------------------------------#

      mysoil$slnm    = with( mysoil
                           , exp( 0.62986 - 0.833 * xclay - 0.529 * slsoc + 0.00593 * slph
                                          + 0.700 * xclay * xclay - 1.400 * xsand * xsilt
                                )#end exp
                           )#end with
      mysoil$slmm    = 1. - 1. / mysoil$slnm
      mysoil$slmu    = mysoil$slnm / (1. - mysoil$slnm)
      mysoil$slpotbp = with( mysoil
                           , -1. / grav 
                           * exp(  0.02294 + 3.526 * xsilt - 2.440 * slsoc + 0.076 * slcec
                                  +0.11331 * slph - 1.900*xsilt*xsilt
                                )#end exp
                           )#end with
      mysoil$malpha  = 1. / mysoil$slpotbp
      mysoil$soilpo  = with( mysoil
                           , 0.81799 + 0.099 * xclay - 3.142e-4 * sldbd
                           + 0.018 * slcec + 0.00451 * slph - 0.050 * xsand * xclay
                           )#end with
      mysoil$soilre  = with( mysoil
                           , 0.22733 - 0.164 * xsand + 0.235 * slcec - 0.00831 * slph
                           + 0.18 * xclay * xclay + 0.26 * xsand * xclay
                           )#end with
      mysoil$soilbp  = with( mysoil
                           , soilre + 2^(1./slmu) * (soilpo - soilre)
                           )#end with
      mysoil$slpotpo = smoist2mpot(mysoil$soilpo,mysoil)
      mysoil$slpots  = - slpotst.MPa * 1.e6 / ( grav * wdns )
      mysoil$slmsts  = mpot2smoist(mysoil$slpots, mysoil)

      #----- Saturated hydraulic conductivity (Ottoni et al. 2019). -----------------------#
      slpot33       = - slpotsr.MPa * 1.e6 / ( grav * wdns )
      effpo         = mysoil$slmsts - mpot2smoist(slpot33, mysoil)
      mysoil$slcons = 19.31 / day.sec * effpo^1.948
      #------------------------------------------------------------------------------------#

      mysoil$slpotfc = - slpotth.MPa * 1.e6 / ( grav * wdns )
      mysoil$sfldcap = mpot2smoist(mysoil$slpotfc, mysoil)

      mysoil$slcpd   = ( (1. - mysoil$slmsts)
                       * ( mysoil$xsand * sand.hcap + mysoil$xsilt * silt.hcap
                         + mysoil$xclay * clay.hcap )
                       + 0.5 * (mysoil$slmsts - mysoil$soilcp) * air.hcap )
      mysoil$slpotwp = - slpotwp.MPa * 1.e6 / ( grav * wdns )
      mysoil$soilwp  = mpot2smoist(mysoil$slpotwp, mysoil)
      #------------------------------------------------------------------------------------#
   }else if (isoil.hydro == 3){
      #------------------------------------------------------------------------------------#
      #     Use Saxton and Rawls (2006).                                                   #
      #------------------------------------------------------------------------------------#
      #----- First solution (mind that psi is in kPa). ------------------------------------#
      th1500t = with( mysoil
                    , - 0.024 * xsand + 0.4870 * xclay + 0.006 * slsoc
                      + 0.005 * xsand * slsoc  - 0.013 * xclay * slsoc
                      + 0.068 * xsand * xclay  + 0.031
                    )#end with
      th33t   = with( mysoil
                    , - 0.251 * xsand + 0.1950 * xclay + 0.011 * slsoc
                      + 0.006 * xsand * slsoc  - 0.027 * xclay * slsoc
                      + 0.452 * xsand * xclay + 0.299
                    )#end with
      thsm33t = with( mysoil
                    , + 0.278 * xsand + 0.0340 * xclay + 0.022 * slsoc
                      - 0.018 * xsand * slsoc  - 0.027 * xclay * slsoc
                      - 0.584 * xsand * xclay + 0.078
                    )#end with
      psiet   = with( mysoil
                    , - 21.67 * xsand - 27.93 * xclay - 81.97  * thsm33t
                      + 71.12 * xsand * thsm33t + 8.29 * xclay * thsm33t
                      + 14.05 * xsand * xclay + 27.16
                    )#end with
      #---- Second solution. --------------------------------------------------------------#
      mysoil$soilwp  = 1.14 * th1500t  - 0.02
      mysoil$sfldcap = 1.283 * th33t^2 + 0.626 * th33t - 0.015
      thsm33         = 1.636 * thsm33t - 0.107
      mysoil$slmsts  = mysoil$sfldcap + thsm33 - 0.097 * mysoil$xsand + 0.043
      mysoil$soilpo  = mysoil$slmsts
      mysoil$soilbp  = mysoil$slmsts
      mysoil$slpots  = - ( 0.02 * psiet^2 + 0.887 * psiet - 0.70 ) / grav
      mysoil$slpotpo = mysoil$slpots
      mysoil$slpotbp = mysoil$slpots
      mysoil$soilre  = 0.
      #---- Derived parameters. -----------------------------------------------------------#
      mysoil$slpotfc = - slpotsr.MPa * 1.e6 / ( grav * wdns )
      mysoil$slpotwp = - slpotwp.MPa * 1.e6 / ( grav * wdns )
      mysoil$sltt    = 1.0
      mysoil$slnm    = with(mysoil, log(sfldcap/soilwp) / log(slpotwp/slpotfc))
      mysoil$slbs    = 1. / mysoil$slnm
      mysoil$slmm    = 2. + mysoil$sltt + 2. * mysoil$slbs
      mysoil$slmu    = -1. / mysoil$slmm
      mysoil$slas    = with(mysoil,(slpotpo-slpotfc)/(slmsts-sfldcap))
      mysoil$slcons  = with(mysoil,1.93/hr.sec * (slmsts-sfldcap)^(3-slnm))
      mysoil$slcpd   = ( (1. - mysoil$slmsts)
                       * ( mysoil$xsand * sand.hcap + mysoil$xsilt * silt.hcap
                         + mysoil$xclay * clay.hcap )
                       + 0.5 * (mysoil$slmsts - mysoil$soilcp) * air.hcap )
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Calculate the derived properties in case this is not bedrock.                    #
   #---------------------------------------------------------------------------------------#
   if (mysoil$ntext != 13){
      mysoil$fhydraul  = 2.0
      mysoil$slpotcp   = - slpotcp.MPa * 1.e6 / ( wdns * grav )
      mysoil$slpotld   = - slpotld.MPa * 1.e6 / ( wdns * grav )
      mysoil$slpotfr   = - slpotfr.MPa * 1.e6 / ( wdns * grav )
      mysoil$soilcp    = mpot2smoist(mysoil$slpotcp, mysoil)
      mysoil$soilfr    = mpot2smoist(mysoil$slpotfr, mysoil)
      mysoil$soilld    = mpot2smoist(mysoil$slpotld, mysoil)
      mysoil$thcr.orig = mysoil$soilwp + theta.crit * (mysoil$slmsts - mysoil$soilwp)
      mysoil$fire.orig = mysoil$soilcp + sm.fire    * (mysoil$slmsts - mysoil$soilcp)
      mysoil$thcr.spot = mpot2smoist(mysoil$slpotld, mysoil)
      mysoil$fire.spot = mpot2smoist(mysoil$slpotfr, mysoil)
      #------------------------------------------------------------------------------------#
   }#end if (mysoil$ntext != 13)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Soil conductivity.                                                               #
   #---------------------------------------------------------------------------------------#
   ksand = 3. * h2o.cond / ( 2. * h2o.cond + sand.cond )
   ksilt = 3. * h2o.cond / ( 2. * h2o.cond + silt.cond )
   kclay = 3. * h2o.cond / ( 2. * h2o.cond + clay.cond )
   kair  = 3. * h2o.cond / ( 2. * h2o.cond +  air.cond )

   mysoil$thcond0 = (   ksand * mysoil$xsand  * ( 1. - mysoil$slmsts ) * sand.cond 
                      + ksilt * mysoil$xsilt  * ( 1. - mysoil$slmsts ) * silt.cond
                      + kclay * mysoil$xclay  * ( 1. - mysoil$slmsts ) * clay.cond
                      + kair                  *       mysoil$slmsts    *  air.cond  )
   mysoil$thcond1 = h2o.cond - kair * air.cond
   mysoil$thcond2 = (   ksand * mysoil$xsand  * ( 1. - mysoil$slmsts )
                      + ksilt * mysoil$xsilt  * ( 1. - mysoil$slmsts )
                      + kclay * mysoil$xclay  * ( 1. - mysoil$slmsts )
                      + kair                  *        mysoil$slmsts   )
   mysoil$thcond3 = 1. - kair
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Make sure everything makes sense.  In case it doesn't, check.                     #
   #---------------------------------------------------------------------------------------#
   cp.fine = all(mysoil$soilre  %<=% mysoil$soilcp )
   wp.fine = all(mysoil$soilcp  %<=% mysoil$soilwp )
   fc.fine = all(mysoil$soilwp  %<=% mysoil$sfldcap)
   po.fine = all(mysoil$sfldcap %<=% mysoil$slmsts )
   if (! (cp.fine && wp.fine && fc.fine && po.fine) ) browser()
   #---------------------------------------------------------------------------------------#


   #----- Select the output format according to the users' choice. ------------------------#
   if (out.dfr) mysoil = as.data.frame(mysoil,stringsAsFactors=FALSE)
   #---------------------------------------------------------------------------------------#

   return(mysoil)
}#end function
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     This function determines the soil class based on the fraction of sand, clay, and     #
# silt separates.                                                                          #
#------------------------------------------------------------------------------------------#
sclass <<- function(sandfrac,clayfrac){
    
   #----- Define the percentage of sand, clay, and silt. ----------------------------------#
   sand = 100. * sandfrac
   clay = 100. * clayfrac
   silt = 100. - sand - clay
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #     Here there is not much we can do other than explore where in the triangle space   #
   # we are.                                                                               #
   #---------------------------------------------------------------------------------------#

   if (silt > 100. | silt < 0. | sand > 100. | sand < 0. | clay > 100. | clay < 0. ) {
      print("---------------------------------------------------")
      print(" At least one of your percentages is screwy...")
      print(paste0("SAND = ",sprintf("%.2f",sand),"%"))
      print(paste0("CLAY = ",sprintf("%.2f",clay),"%"))
      print(paste0("SILT = ",sprintf("%.2f",silt),"%"))
      print("---------------------------------------------------")
      stop ("This soil doesn''t fit into any category...")
      
   }else if(sand > (85.0 + 0.5 * clay)) {
      mysoil =  1 #----- Sand. ------------------------------------------------------------#
   }else if(sand > (70.0 + clay)) {
      mysoil =  2 #----- Loamy sand. ------------------------------------------------------#
   }else if((clay <= 20.0 & sand >= 52.5) | (clay <= 7.5 & silt <= 50.0)) {
      mysoil =  3 #----- Sandy loam. ------------------------------------------------------#
   }else if((clay <= 27.5 & silt > 50.0 & silt <= 80.0) | (silt >  80.0 & clay > 12.5)) {
      mysoil =  4 #----- Silt loam. -------------------------------------------------------#
   }else if(clay > 7.5 & clay <= 27.5 & silt > 27.5 & silt <= 50.0 & sand <= 52.5) {
      mysoil =  5 #----- Loam. ------------------------------------------------------------#
   }else if(clay > 20.0 & clay <= 35.0 & silt <= 27.5 & sand > 45.0) {
      mysoil =  6 #----- Sandy clay loam. -------------------------------------------------#
   }else if(clay > 27.5 & clay <= 40.0 & sand <= 20.0) {
      mysoil =  7 #----- Silty clay loam. -------------------------------------------------#
   }else if(clay > 27.5 & clay <= 40.0 & sand > 20.0 & sand <= 45.0) {
      mysoil =  8 #----- Clayey loam. -----------------------------------------------------#
   }else if(clay > 35.0 & sand > 45.0) {
      mysoil =  9 #----- Sandy clay. ------------------------------------------------------#
   }else if(clay > 40.0 & silt > 40.0) {
      mysoil = 10 #----- Silty clay. ------------------------------------------------------#
   }else if(clay <= 70.0 & sand <= 30.0 & silt <= 30.0) {
      mysoil = 11 #----- Clay. ------------------------------------------------------------#
   }else if( silt > 80.0 & clay <= 12.5) {
      mysoil = 14 #----- Silt. ------------------------------------------------------------#
   }else if( clay > 70.0) {
      mysoil = 15 #----- Heavy clay. ------------------------------------------------------#
   }else if( clay > 40.0 & sand > 30.0 & sand <= 45.0) {
      mysoil = 16 #----- Clayey sand. -----------------------------------------------------#
   }else if( clay > 40.0 & silt > 30.0 & silt <= 40.0) {
      mysoil = 17 #----- Clayey silt. -----------------------------------------------------#
  }else{
      print("---------------------------------------------------")
      print(paste0("SAND = ",sprintf("%.2f",sand),"%"))
      print(paste0("CLAY = ",sprintf("%.2f",clay),"%"))
      print(paste0("SILT = ",sprintf("%.2f",silt),"%"))
      print("---------------------------------------------------")
      stop ("This soil doesn''t fit into any category...")
  }#end if

  return(mysoil)
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function normalises the soil moisture as a function of some of the most         #
# important references.  This is not a linear scale, but it helps deciding how dry or wet  #
# the soil is, and how easy or hard is to evaporate.                                       #
# The scale is linear between the important points:                                        #
# -1. : Residual water (dry air soil moisture)                                             #
#  0. : Wilting point                                                                      #
#  1. : Field capacity                                                                     #
#  2. : Porosity (saturation)                                                              #
#------------------------------------------------------------------------------------------#
soil.scale <<- function(soil.water,soil){
   low = soil.water <  soil$soilcp
   dry = soil.water <= soil$soilwp 
   mid = soil.water >  soil$soilwp & soil.water <= soil$sfldcap
   wet = soil.water >  soil$sfldcap
   sat = soil.water >  soil$slmsts

   sindex      = soil.water
   sindex[low] = - 1.
   sindex[dry] = - 1. + (soil.water[dry] - soil$soilcp ) / (soil$soilwp  - soil$soilcp )
   sindex[mid] =   0. + (soil.water[mid] - soil$soilwp ) / (soil$sfldcap - soil$soilwp )
   sindex[wet] = + 1. + (soil.water[wet] - soil$sfldcap) / (soil$slmsts  - soil$sfldcap)
   sindex[sat] = + 1.

   return(sindex)
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function normalises the soil moisture as a function of some of the most         #
# important references.  This is not a linear scale, but it helps deciding how dry or wet  #
# the soil is, and how easy or hard is to evaporate.                                       #
# The scale is linear between the important points:                                        #
# -1. : Residual water (dry air soil moisture)                                             #
#  0. : Wilting point                                                                      #
#  1. : Field capacity                                                                     #
#  2. : Porosity (saturation)                                                              #
#------------------------------------------------------------------------------------------#
soil.idx2water <<- function(soil.index,soil){
   soil.index[soil.index < -1.0] = -1.0
   soil.index[soil.index > +2.0] = +2.0
   dry = soil.index <= 0.0
   mid = soil.index >  0.0 & soil.index <= 1.0
   wet = soil.index >  1.0

   swater      = soil.index

   swater[dry] = soil$soilcp  + (soil.index[dry] + 1.0) * (soil$soilwp  - soil$soilcp )
   swater[mid] = soil$soilwp  + (soil.index[mid] + 0.0) * (soil$sfldcap - soil$soilwp )
   swater[wet] = soil$sfldcap + (soil.index[wet] - 1.0) * (soil$slmsts  - soil$sfldcap)

   return(swater)
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     Function that converts matric potential (m) into soil moisture (m3/m3).              #
#------------------------------------------------------------------------------------------#
smoist2mpot <<- function(smoist,mysoil){
   #----- Pick method. --------------------------------------------------------------------#
   if (mysoil$method %in% "BDRK"){
      #----- Bedrock, ignore. -------------------------------------------------------------#
      mpot = 0.
      #------------------------------------------------------------------------------------#
   }else if (mysoil$method %in% "vG80"){
      #----- van Genuchten (1980). --------------------------------------------------------#
      smfrac = (smoist - mysoil$soilre) / (mysoil$soilpo - mysoil$soilre)
      smfrac = 0. * smfrac + pmax(0.,pmin(1.,smfrac))
      mpot   = mysoil$slpotbp * (smfrac^mysoil$slmu - 1.)^(1.-mysoil$slmm)
      #------------------------------------------------------------------------------------#
   }else if (mysoil$method %in% "SR06"){
      #----- Saxton and Rawls (2006). -----------------------------------------------------#
      smbelow = pmin(1.,smoist/mysoil$sfldcap)
      smabove = pmax(0.,smoist - mysoil$sfldcap)
      mpot    = with(mysoil, slpotfc / smbelow^slbs + slas * smabove)
      #------------------------------------------------------------------------------------#
   }else{
      #----- Brooks and Corey (1964). -----------------------------------------------------#
      smfrac = (smoist - mysoil$soilre) / (mysoil$slmsts - mysoil$soilre)
      smfrac = 0. * smfrac + pmax(0.,pmin(1.,smfrac))
      mpot   = mysoil$slpots / smfrac ^ mysoil$slbs
      #------------------------------------------------------------------------------------#
   }#end if (mysoil$method %in% "vG80")
   #---------------------------------------------------------------------------------------#

   return(mpot)
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     Function that converts soil moisture (m3/m3) into matric potential (m).              #
#------------------------------------------------------------------------------------------#
mpot2smoist <<- function(mpot,mysoil){
   #----- Pick method. --------------------------------------------------------------------#
   if (mysoil$method %in% "BDRK"){
      #----- Bedrock, ignore. -------------------------------------------------------------#
      smoist = 0.
      #------------------------------------------------------------------------------------#
   }else if (mysoil$method %in% "vG80"){
      #----- van Genuchten (1980). --------------------------------------------------------#
      smfrac = ( 1. / ( 1. + (mysoil$malpha*mpot)^mysoil$slnm) )^mysoil$slmm
      smfrac = 0. * smfrac + pmax(0.,pmin(1.,smfrac))
      smoist = (1. - smfrac) * mysoil$soilre + smfrac * mysoil$soilpo
      #------------------------------------------------------------------------------------#
   }else if (mysoil$method %in% "SR06"){
      #----- Saxton and Rawls (2006). -----------------------------------------------------#
      smoist = with( mysoil
                   , ifelse( test = mpot %<% slpotfc
                           , yes  = sfldcap * ( slpotfc / mpot  ) ^ (1. / slbs )
                           , no   = sfldcap + ( slpotfc - slpots) / slas
                           )#end ifelse
                   )#end with
      #------------------------------------------------------------------------------------#
   }else{
      #----- Brooks and Corey (1964). -----------------------------------------------------#
      smfrac = ( mysoil$slpots / mpot ) ^ mysoil$slnm
      smfrac = 0. * smfrac + pmax(0.,pmin(1.,smfrac))
      smoist = (1. - smfrac) * mysoil$soilre + smfrac * mysoil$slmsts
      #------------------------------------------------------------------------------------#
   }#end if (mysoil$method %in% "vG80")
   #---------------------------------------------------------------------------------------#



   return(smoist)
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     Function that converts soil moisture (m3/m3) into hydraulic conductivity (kg/m2/s).  #
#------------------------------------------------------------------------------------------#
smoist2hydcond <<- function(smoist,mysoil){


   #----- Pick method. --------------------------------------------------------------------#
   if (mysoil$method %in% "vG80"){
      #----- van Genuchten (1980). --------------------------------------------------------#
      smfrac  = (smoist - mysoil$soilre) / (mysoil$soilpo - mysoil$soilre)
      smfrac  = 0. * smfrac + pmax(0.,pmin(1.,smfrac))
      fibf    = 1. - (1. - smfrac^(1./mysoil$slmm))^mysoil$slmm
      hydcond = mysoil$slcons * smfrac ^ mysoil$sltt * fibf * fibf
      #------------------------------------------------------------------------------------#
   }else{
      #----- Saxton and Rawls (2006) or Brooks and Corey (1964). --------------------------#
      smfrac  = (smoist - mysoil$soilre) / (mysoil$slmsts - mysoil$soilre)
      smfrac  = 0. * smfrac + pmax(0.,pmin(1.,smfrac))
      hydcond = mysoil$slcons * smfrac ^ mysoil$slmm
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#

   hydcond = 0. * hydcond + pmax(hydcond.min,hydcond)
   return(hydcond)
}#end function
#==========================================================================================#
#==========================================================================================#



#==========================================================================================#
#==========================================================================================#
#     This function creates soil layers for ED-2 that monotonically increase in depth.     #
#------------------------------------------------------------------------------------------#
make.slz <<- function(from,to=-0.02,length.out=16,digits=3,out.char=FALSE){
   fmt  = paste0("%.",digits,"f")
   ehgt = log(from/to) / log(length.out)
   ilyr = rev(sequence(length.out))
   zlyr = to * ilyr ^ ehgt
   if (out.char){
      ans  = sprintf(fmt,zlyr)
   }else{
      ans  = round(zlyr,digits)
   }#end if (out.char)
   return(ans)
}#end make.slz
#==========================================================================================#
#==========================================================================================#
