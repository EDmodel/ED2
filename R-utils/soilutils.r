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
#==========================================================================================#
#==========================================================================================#


#==========================================================================================#
#==========================================================================================#
#     This function finds the soil parameters.                                             #
#------------------------------------------------------------------------------------------#
soil.params = function(ntext,isoilflg,slxsand,slxclay){
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
   fieldcp.K  = 0.1  # hydraulic conduct. at field capacity                       [ mm/day]
   soilcp.MPa = 3.1  # soil-water potential for air dry soil                      [    MPa]
   soilwp.MPa = 1.5  # soil-water potential at wilting point                      [    MPa]
   soilld.MPa = 0.75 # soil-water potential that plants start dropping leaves     [    MPa]
   soilfr.MPa = 1.40 # soil-water potential that triggers fires                   [    MPa]
   theta.crit = 0.11 # fractional soil moisture that plants start dropping leaves [  m3/m3]
   sm.fire    = 0.08 # fractional soil moisture that triggers fires               [  m3/m3]
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
   mysoil = list( ntext      = NA
                , name       = NA
                , key        = NA
                , xsand      = NA
                , xclay      = NA
                , xsilt      = NA
                , slbs       = 0.
                , slpots     = 0.
                , slcons     = 0.
                , slmsts     = 0.
                , sfldcap    = 0.
                , soilcp     = 0.
                , soilwp     = 0.
                , slcpd      = 0.
                , thcr.orig  = 0.
                , thcr.spot  = 0.
                , fire.orig  = 0.
                , fire.spot  = 0.
                , slpotcp    = 0.
                , slpotwp    = 0.
                , slpotfc    = 0.
                , slpotld    = 0.
                , slpotfr    = 0.
                , thcond0    = 0.
                , thcond1    = 0.
                , thcond2    = 0.
                , thcond3    = 0.
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




   #---------------------------------------------------------------------------------------#
   #       Set up primary properties.                                                      #
   #---------------------------------------------------------------------------------------#
   if (mysoil$ntext == 13){
      #----- Bedrock.  Most things are zero, because it is an impermeable soil. -----------#
      mysoil$slbs      =  0.
      mysoil$slpots    =  0.
      mysoil$slcons    =  0.
      mysoil$slmsts    =  0.
      mysoil$sfldcap   =  0.
      mysoil$soilcp    =  0.
      mysoil$soilwp    =  0.
      mysoil$slcpd     =  2130000.
      #------------------------------------------------------------------------------------#
   }else if (mysoil$ntext == 12){
      #------------------------------------------------------------------------------------#
      #      Peat.  High concentration of organic matter.  Mineral soil equations don't    #
      # apply here.                                                                        #
      #------------------------------------------------------------------------------------#
      mysoil$slbs    =  6.180000
      mysoil$slpots  = -0.534564359
      mysoil$slcons  =  2.357930e-6
      mysoil$slmsts  =  0.469200
      mysoil$sfldcap =  0.285709966
      mysoil$slcpd   =  874000.
      #------------------------------------------------------------------------------------#
   }else{
      #------------------------------------------------------------------------------------#
      #      Mineral soil.  Use the standard ED-2.2 equations.                             #
      #------------------------------------------------------------------------------------#
      mysoil$slbs    = 3.10 + 15.7*mysoil$xclay - 0.3*mysoil$xsand
      mysoil$slpots  = -0.01 * (10.^(2.17 - 0.63*mysoil$xclay - 1.58*mysoil$xsand))
      mysoil$slcons  = ( (10.^(-0.60 + 1.26*mysoil$xsand - 0.64*mysoil$xclay)) 
                       * 0.0254/hr.sec )
      mysoil$slmsts  = (50.5 - 14.2*mysoil$xsand - 3.7*mysoil$xclay) / 100.
      mysoil$sfldcap = ( mysoil$slmsts * ( (fieldcp.K/1000./day.sec)/mysoil$slcons)
                                        ^ (1. / (2.*mysoil$slbs+3.)) )
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
      mysoil$slpotcp   = - soilcp.MPa * 1000. / grav
      mysoil$slpotwp   = - soilwp.MPa * 1000. / grav
      mysoil$slpotld   = - soilld.MPa * 1000. / grav
      mysoil$slpotfr   = - soilfr.MPa * 1000. / grav
      mysoil$soilcp    = mpot2smoist(mysoil$slpotcp, mysoil)
      mysoil$soilwp    = mpot2smoist(mysoil$slpotwp, mysoil)
      mysoil$thcr.orig = mysoil$soilwp + theta.crit * (mysoil$slmsts - mysoil$soilwp)
      mysoil$fire.orig = mysoil$soilcp + sm.fire    * (mysoil$slmsts - mysoil$soilcp)
      mysoil$thcr.spot = mpot2smoist(mysoil$slpotld, mysoil)
      mysoil$fire.spot = mpot2smoist(mysoil$slpotfr, mysoil)
      mysoil$slpotfc   = smoist2mpot(mysoil$sfldcap, mysoil)
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Soil conductivity.                                                               #
   #---------------------------------------------------------------------------------------#
   ksand = 3. * h2o.cond / ( 2. * h2o.cond + sand.cond )
   ksilt = 3. * h2o.cond / ( 2. * h2o.cond + silt.cond )
   kclay = 3. * h2o.cond / ( 2. * h2o.cond + clay.cond )
   kair  = 3. * h2o.cond / ( 2. * h2o.cond +  air.cond )
     
   mysoil$thcond0 = ( ksand * mysoil$xsand  * ( 1. - mysoil$slmsts ) * sand.cond 
                      + ksilt * mysoil$xsilt  * ( 1. - mysoil$slmsts ) * silt.cond
                      + kclay * mysoil$xclay  * ( 1. - mysoil$slmsts ) * clay.cond
                      + kair                  *       mysoil$slmsts    *  air.cond  )
   mysoil$thcond1 = h2o.cond - kair * air.cond
   mysoil$thcond2 = ( ksand * mysoil$xsand  * ( 1. - mysoil$slmsts )
                      + ksilt * mysoil$xsilt  * ( 1. - mysoil$slmsts )
                      + kclay * mysoil$xclay  * ( 1. - mysoil$slmsts )
                      + kair                  *        mysoil$slmsts   )
   mysoil$thcond3 = 1. - kair
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
sclass = function(sandfrac,clayfrac){
    
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
      print(paste("SAND = ",sprintf("%.2f",sand),"%",sep=""))
      print(paste("CLAY = ",sprintf("%.2f",clay),"%",sep=""))
      print(paste("SILT = ",sprintf("%.2f",silt),"%",sep=""))
      print("---------------------------------------------------")
      stop ("This soil doesn''t fit into any category...")
      
   }else if(sand > 85.0 + 0.5 * clay) {
      mysoil =  1 #----- Sand. ------------------------------------------------------------#
   }else if(sand > 70.0 + clay) {
      mysoil =  2 #----- Loamy sand. ------------------------------------------------------#
   }else if((clay <= 20.0 & sand > 52.5) | (clay <= 7.5 & silt <= 50.0)) {
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
      print(paste("SAND = ",sprintf("%.2f",sand),"%",sep=""))
      print(paste("CLAY = ",sprintf("%.2f",clay),"%",sep=""))
      print(paste("SILT = ",sprintf("%.2f",silt),"%",sep=""))
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
soil.scale = function(soil.water,soil){
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
soil.idx2water = function(soil.index,soil){
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
smoist2mpot = function(smoist,mysoil){
   smfrac = smoist / mysoil$slmsts
   mpot   = mysoil$slpots / smfrac ^ mysoil$slbs
   return(mpot)
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
mpot2smoist = function(mpot,mysoil){
   smfrac = ( mpot / mysoil$slpots) ^ (-1. / mysoil$slbs)
   smoist = smfrac * mysoil$slmsts
   return(smoist)
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
smoist2hydcond = function(smoist,mysoil){
   smfrac  = smoist / mysoil$slmsts
   hydcond = mysoil$slcons * smfrac ^ (2. * mysoil$slbs + 3.)
   return(hydcond)
}#end function
#==========================================================================================#
#==========================================================================================#
