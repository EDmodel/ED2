#----- Some dimensions based on ED-2.2 default. -------------------------------------------#
npft       <<- 17 # Number of plant functional types
nlu        <<-  6 # Number of land use types.
nstyp      <<- 17 # Number of default soil types
#------------------------------------------------------------------------------------------#



#----- Radiation thresholds. --------------------------------------------------------------#
cosz.min      <<- 0.03 # cos(89*pi/180) # Minimum cosine of zenith angle
cosz.highsun  <<- cos(84*pi/180)        # Zenith angle to not be called sunrise or sunset
cosz.twilight <<- cos(96*pi/180)        # Cosine of the end of civil twilight
fvis.beam.def <<- 0.43
fnir.beam.def <<- 1.0 - fvis.beam.def
fvis.diff.def <<- 0.52
fnir.diff.def <<- 1.0 - fvis.diff.def
phap.min      <<- 25                    # Minimum incoming radiation to be considered 
                                        # daytime
#------------------------------------------------------------------------------------------#


#----- Minimum R2 that we consider meaningful. --------------------------------------------#
r2.min        <<- 0.36
#----- Typical p-value below which we reject the NULL hypothesis. -------------------------#
pval.max      <<- 0.05
#----- p-value below which we reject the NULL hypothesis for u* filters. ------------------#
pval.ustar    <<- 0.05
#------------------------------------------------------------------------------------------#



#----- Minimum pixel resolution which would allow MacArthur-Horn corrections. -------------#
min.mh.pixres  <<- 5 # metres
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Tolerance for root-finding methods.                                                  #
#------------------------------------------------------------------------------------------#
toler  <<- 1.e-6  # Tolerance for the root-finder algorithms
maxfpo <<- 60     # Maximum number of Regula-Falsi iterations
maxit  <<- 150    # Maximum number of iterations in general
#------------------------------------------------------------------------------------------#


#----- Colours for months. ----------------------------------------------------------------#
month.at   <<- (sequence(12) - 0.5 ) / 12
month.off  <<- -1.5/12
month.cols <<- hsv( h = month.at
                  , s = 0.75 + 0.25 * cos((month.at - month.off) * 2 * pi)
                  , v = 0.75 + 0.25 * sin((month.at - month.off) * 2 * pi)
                  )#end hsv
#------------------------------------------------------------------------------------------#

#------ Define the colours and labels for IGBP land use maps. -----------------------------#
igbp.col <<- c( H2O = RGB(   0,  20,  82)
              , ENF = RGB(   0, 100, 164)
              , EBF = RGB(   0,  63,   0)
              , DNF = RGB(  85, 192, 255)
              , DBF = RGB(  97, 255,  96)
              , MXF = RGB(   0, 156,  76)
              , CSH = RGB( 255, 105,   0)
              , OSH = RGB( 207, 151, 114)
              , WSV = RGB(  72, 144,   0)
              , SAV = RGB( 150, 255,  40)
              , GSL = RGB( 233, 198,   9)
              , PWL = RGB(  69,  35, 232)
              , CRL = RGB( 166, 131,   0)
              , URB = RGB(  96,  96,  96)
              , CNV = RGB( 115, 158,  80)
              , ICE = RGB(   0, 222, 255)
              , BRN = RGB( 144,   3,   2)
              , UND = RGB( 222, 222, 222)
              )#end igbp.col
igbp.leg <<- names(igbp.col)
igbp.val <<- seq_along(igbp.col)-1
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Colours, labels, and conversion to simplified TerraClass categories.                 #
#------------------------------------------------------------------------------------------#
tc2012.col  <<- c( H2O = RGB(  12,  84, 170) #  0 - Water
                 , FOR = RGB(   0,  79,   0) #  1 - Forest
                 , SCD = RGB(   0, 177,   0) #  2 - Secondary growth
                 , RFR = RGB(  58, 204,   0) #  3 - Reforestation
                 , MXL = RGB( 155, 155, 255) #  4 - Mixed landscape
                 , NFR = RGB( 255, 175,   0) #  5 - Non-forest
                 , PST = RGB( 243, 240,  96) #  6 - Pasture
                 , CRL = RGB( 176, 147,  81) #  7 - Cropland
                 , D12 = RGB( 205,  85,   0) #  8 - Deforestation 2012
                 , MIN = RGB( 148,  33, 147) #  9 - Mining
                 , URB = RGB( 113,  40, 197) # 10 - Urban
                 , UND = RGB( 222, 222, 222) # 11 - Undefined
                 )#end tc2012.col
tc2012.leg  <<- names(tc2012.col)
tc2012.desc <<- c("Water","Forest","Secondary Growth","Reforestation"
                 ,"Mixed Landscape","Non-forest","Pasture","Cropland"
                 ,"Deforestation 2012","Mining","Urban","Unclassified")
tc2012.val  <<- seq_along(tc2012.col)-1
tc2012.t2c  <<- c(7,11,10,8,1,0,9,4,5,5,6,6,6,3,2,2)
#------------------------------------------------------------------------------------------#




#------ Define the colours and labels for ASPRS (LAS/LiDAR) classes. ----------------------#
asprs.col <<- c( CNC = RGB( 222, 222, 222)
               , UND = RGB( 202, 202, 202)
               , GND = RGB( 144,   3,   2)
               , LVG = RGB( 166, 131,   0)
               , MVG = RGB(   0, 100, 164)
               , HVG = RGB(   0,  63,   0)
               , BLD = RGB(  96,  96,  96)
               , MKP = RGB( 233, 198,   9)
               , H2O = RGB(   0,  20,  82)
               , RSV = RGB( 192, 192, 192)
               , R11 = RGB( 192, 192, 192)
               , OLP = RGB(   0, 222, 255)
               )#end asprs.col
asprs.leg <<- names(asprs.col)
asprs.val <<- seq_along(asprs.col)-1
#------------------------------------------------------------------------------------------#




#----- Define some default legend colours and names. --------------------------------------#
lunames   <<- c("Agricultural","Plantation","Tree fall"
               ,"Burnt","Abandoned","Logged","Total")
lucols    <<- c("#E65C17","#2996CC","#306614","#990F0F","#2996CC","#A3CC52",all.colour)

distnames <<- c("AGR->AGR" ,"FPL->AGR" ,"TFL->AGR","BRN->AGR","ABN->AGR","LOG->AGR"
               ,"AGR->FPL" ,"FPL->FPL" ,"TFL->FPL","BRN->FPL","ABN->FPL","LOG->FPL"
               ,"AGR->TFL" ,"FPL->TFL" ,"TFL->TFL","BRN->TFL","ABN->TFL","LOG->TFL"
               ,"AGR->BRN" ,"FPL->BRN" ,"TFL->BRN","BRN->BRN","ABN->BRN","LOG->BRN"
               ,"AGR->ABN" ,"FPL->ABN" ,"TFL->ABN","BRN->ABN","ABN->ABN","LOG->ABN"
               ,"AGR->LOG" ,"FPL->LOG" ,"TFL->LOG","BRN->LOG","ABN->LOG","LOG->LOG")
distcols  <<- c("#66005A","#A80095","#E500CB","#FF4CEA","#FF998F","#FFCCF9"
               ,"#005A66","#00889A","#00B3CB","#00D3EF","#96F3FF","#CBF9FF"
               ,"#DEDEDE","#DEDEDE","#596500","#A0B310","#D1ED00","#F2FF91"
               ,"#DEDEDE","#DEDEDE","#983A2E","#D64F3D","#F29C91","#FFD1CB"
               ,"#007F0F","#00BE16","#DEDEDE","#DEDEDE","#DEDEDE","#DEDEDE"
               ,"#DEDEDE","#DEDEDE","#4C3FB2","#6E5AFF","#B0A6FF","#DED9FF")
#------------------------------------------------------------------------------------------#


#----- Growth respiration factor (to estimate when the actual variable is not there). -----#
growth.resp.fac <<- c(rep(0.333,times=5),rep(0.4503,times=3),rep(0,times=3)
                     ,rep(0.333,times=4))
#------------------------------------------------------------------------------------------#



#----- fswh is the FSW that plants experience and they are happy (wilting point = 0). -----#
fswh <<- 0.99
#------------------------------------------------------------------------------------------#


#----- Years for which we have eddy flux tower. -------------------------------------------#
eft.year <<- c(     1998,     1999,     2000,     2001,     2002,     2003,     2004
              ,     2005,     2006,     2007,     2008,     2009,     2010,     2011
              ,     2012,     2013)
eft.pch  <<- c(        0,        1,        6,        3,        4,        5,       13
              ,       15,       16,       17,       18,        9,        8,        7
              ,       14,        2)
eft.col  <<- c("#7D6E93","#B49ED2","#39025D","#520485","#042E88","#0742C3","#00480E"
              ,"#006715","#31B223","#46FF32","#AB8C3D","#F5C858","#B23C00","#FF5700"
              ,"#70000E","#A00014")
#------------------------------------------------------------------------------------------#



#----- Standard colours and names for soil classes. ---------------------------------------#
stext.cols  <<- c("gold","chartreuse","limegreen","darkgreen","purple3"
                 ,"deepskyblue","aquamarine","slateblue2","darkorange3","sienna"
                 ,"firebrick","grey61","grey29","orchid","olivedrab","goldenrod"
                 ,"steelblue")
stext.names <<- c("Sand","Loamy Sand","Sandy loam","Silt loam","Loam","Sandy clay loam"
                 ,"Silty clay loam","Clayey loam","Sandy clay","Silty clay","Clay","Peat"
                 ,"Bedrock","Silt","Heavy clay","Clayey sand","Clayey silt")
stext.acron <<- c("Sa","LSa","SaL","SiL","L","SaCL","SiCL","CL"
                 ,"SaC","SiC","C","Pe","BR","Si","CC","CSa","CSi")
#------------------------------------------------------------------------------------------#




#==========================================================================================#
#==========================================================================================#
#     Patch information for POV-Ray.                                                       #
#------------------------------------------------------------------------------------------#
pov.dbh.min    <<-    5    # Minimum DBH to be plotted
pov.patch.xmax <<-   16    # Size of each sub-plot in the x direction [m]
pov.patch.ymax <<-   16    # Size of each sub-plot in the y direction [m]
pov.nx.patch   <<-   25    # Number of sub-plots in each x transect
pov.ny.patch   <<-   25    # Number of sub-plots in each y transect
pov.nxy.patch  <<- pov.nx.patch  * pov.ny.patch
pov.total.area <<- pov.nxy.patch * pov.patch.xmax * pov.patch.ymax
pov.x0         <<- rep( x     = -200 + seq(from=0,to=pov.nx.patch-1) * pov.patch.xmax
                      , times = pov.ny.patch
                      )#end rep
pov.y0         <<- rep( x     = -200 + seq(from=0,to=pov.ny.patch-1) * pov.patch.ymax
                      , each  = pov.nx.patch
                      )#end rep
#------------------------------------------------------------------------------------------#




#==========================================================================================#
#==========================================================================================#
#     Census-related thresholds.  If the script has different values, ignore these and     #
# use whatever the main script says.                                                       #
#------------------------------------------------------------------------------------------#
#----- Minimum height to be considered for "ground-based observations". -------------------#
if ( "census.height.min" %in% ls()){
   census.height.min <<- census.height.min
}else{
   census.height.min <<- 1.5 
}#end if
#----- Minimum DBH to be considered for "ground-based observations". ----------------------#
if ( "census.dbh.min" %in% ls()){
   census.dbh.min <<- census.dbh.min
}else{
   census.dbh.min <<- 10.
}#end if
#----- Minimum DBH to be considered for "ground-based observations". ----------------------#
if ( "recruit.dbh.min" %in% ls()){
   recruit.dbh.min <<- recruit.dbh.min
}else{
   recruit.dbh.min <<- 0.16
}#end if
#------------------------------------------------------------------------------------------#


#==========================================================================================#
#==========================================================================================#
#     Define which DBH classes to use based on the DBH flag.                               #
#------------------------------------------------------------------------------------------#
if ("idbh.type" %in% ls()){
   idbh.type <<- idbh.type
}else{
   idbh.type <<- 1   
}#end if
if (idbh.type == 1){
   ndbh       <<- 11
   ddbh       <<- 10
   classdbh   <<- seq(from=0,to=(ndbh-1)*ddbh,by=ddbh)
   breakdbh   <<- c(-Inf,classdbh[-1],Inf)
   dbhlabel   <<- "11_szclss"
   dbhkeys    <<- paste(classdbh,"-",c(classdbh[-1],Inf),sep="")
   dbhnames   <<- paste( c("<",paste(classdbh[-c(1,ndbh)],"-",sep=""),">")
                       , c(classdbh[-1],classdbh[ndbh]),"cm"
                       , sep=""
                       )#end paste
   dbhcols    <<- c(         "purple3",   "mediumpurple1",      "royalblue4"
                   ,      "steelblue3",     "deepskyblue",         "#004E00"
                   ,     "chartreuse3",      "olivedrab3", "lightgoldenrod3"
                   ,         "yellow3",     "darkorange1",       "firebrick"
                   ,        all.colour
                   )#end c

}else if (idbh.type == 2){
   ndbh       <<-  5
   classdbh   <<- c(0,10,20,35,55)
   breakdbh   <<- c(-Inf,classdbh[-1],Inf)
   dbhlabel   <<- "05_szclss"
   dbhkeys    <<- paste(classdbh,"-",c(classdbh[-1],Inf),sep="")
   dbhnames   <<- paste( c("<",paste(classdbh[-c(1,ndbh)],"-",sep=""),">")
                       , c(classdbh[-1],classdbh[ndbh]),"cm"
                       , sep=""
                       )#end paste
   dbhcols    <<- c(      "royalblue3",     "chartreuse3" ,         "yellow3"
                   ,     "darkorange1",       "firebrick" ,        all.colour
                   )#end c
}else if (idbh.type == 3){
   ndbh       <<-  4
   classdbh   <<- c(0,10,35,55)
   dbhlabel   <<- "04_szclss"
   breakdbh   <<- c(-Inf,classdbh[-1],Inf)
   dbhkeys    <<- paste(classdbh,"-",c(classdbh[-1],Inf),sep="")
   dbhnames   <<- paste( c("<",paste(classdbh[-c(1,ndbh)],"-",sep=""),">")
                       , c(classdbh[-1],classdbh[ndbh]),"cm"
                       , sep=""
                       )#end paste
   dbhcols    <<- c(         "#3B24B3",     "#2996CC"
                   ,         "#FF9466",     "#990F0F",         all.colour
                   )#end c
}else if (idbh.type == 4){
   ndbh       <<-  6
   classdbh   <<- c(0,2,10,20,45,70)
   breakdbh   <<- c(-Inf,classdbh[-1],Inf)
   dbhlabel   <<- "06_szclss"
   dbhkeys    <<- paste(classdbh,"-",c(classdbh[-1],Inf),sep="")
   dbhnames   <<- paste( c("<",paste(classdbh[-c(1,ndbh)],"-",sep=""),">")
                       , c(classdbh[-1],classdbh[ndbh]),"cm"
                       , sep=""
                       )#end paste
   dbhcols    <<- c(         "purple3",       "royalblue3",     "chartreuse3"
                   ,         "yellow3",      "darkorange1",       "firebrick"
                   ,        all.colour
                   )#end c
}else if (idbh.type == 5){
  ndbh       <<-  17
  classdbh   <<- c(0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,20)
  breakdbh   <<- c(-Inf,classdbh[-1],Inf)
  dbhlabel   <<- "07_szclss"
  dbhkeys    <<- paste(classdbh,"-",c(classdbh[-1],Inf),sep="")
  dbhnames   <<- paste( c("<",paste(classdbh[-c(1,ndbh)],"-",sep=""),">")
                        , c(classdbh[-1],classdbh[ndbh]),"cm"
                        , sep=""
  )#end paste
  dbhcols    <<- c(         "purple3",       "royalblue3",     "chartreuse3"
                  ,         "yellow3",      "darkorange1",       "firebrick"
                  ,        all.colour
  )#end c
}else{
   cat(" In globdims.r:","\n")
   cat(" IDBH.TYPE = ",idbh.type,"\n")
   stop(" Invalid IDBH.TYPE, it must be between 1 and 5 (feel free to add more options.")
}#end if
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     Define which height classes to use based on the DBH flag.                            #
#------------------------------------------------------------------------------------------#
if ("ihgt.type" %in% ls()){
   ihgt.type <<- ihgt.type
}else{
   ihgt.type <<- 1   
}#end if
if (ihgt.type == 1){
   nhgt       <<- 8
   classhgt   <<- c(0,5,10,15,20,25,30,35)-0.001
   roundhgt   <<- round(classhgt,0)
   breakhgt   <<- c(-Inf,classhgt[-1],Inf)
   hgtlabel   <<- "08_htclss"
   hgtkeys    <<- paste(classhgt,"-",c(roundhgt[-1],Inf),sep="")
   hgtnames   <<- paste( c("<",paste(roundhgt[-c(1,nhgt)],"-",sep=""),">")
                       , c(roundhgt[-1],roundhgt[nhgt]),"m"
                       , sep=""
                       )#end paste
   hgtcols    <<- c(      "slateblue4",      "steelblue3",     "deepskyblue"
                   ,         "#004E00",     "chartreuse3", "lightgoldenrod3"
                   ,     "darkorange1",       "firebrick",        all.colour
                   )#end c
}else if (ihgt.type == 2){
   nhgt       <<-  13
   classhgt   <<- c(0,1,4,7,10,13,16,19,22,25,28,31,34)
   hgtlabel   <<- "13_szclss"
   breakhgt   <<- c(-Inf,classhgt[-1],Inf)
   hgtkeys    <<- paste(classhgt,"-",c(classhgt[-1],Inf),sep="")
   hgtnames   <<- paste( c("<",paste(classhgt[-c(1,nhgt)],"-",sep=""),">")
                       , c(classhgt[-1],classhgt[nhgt]),"cm"
                       , sep=""
                       )#end paste
   hgtcols    <<- c(         "purple3",   "mediumpurple1",      "royalblue4"
                   ,      "steelblue3",     "deepskyblue",         "#004E00"
                   ,     "chartreuse3",      "olivedrab3", "lightgoldenrod3"
                   ,         "yellow3",     "darkorange1",            "red3"
                   ,      "firebrick4",        all.colour
                   )#end c
}else{
   cat(" In globdims.r:","\n")
   cat(" IHGT.TYPE = ",ihgt.type,"\n")
   stop(" Invalid IHGT.TYPE, it must be between 1 and 2 (feel free to add more options.")
}#end if
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     Define which height classes to use based on the DBH flag.                            #
#------------------------------------------------------------------------------------------#
if ("isld.type" %in% ls()){
   isld.type <<- isld.type
}else{
   isld.type <<- 1
}#end if
if (isld.type == 1){
   nsld       <<- 12
   belowsld   <<- c(-6.0,-5.5,-5.0,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5)
   abovesld   <<- c(belowsld[-1],0)
   thicksld  <<- diff(c(belowsld,0))
   roundsld   <<- round(belowsld,2)
   breaksld   <<- c(-Inf,belowsld[-1],Inf)
   sldlabel   <<- "cc_sdclss"
   sldkeys    <<- paste(roundsld,"-",c(roundsld[-1],0),sep="")
   sldnames   <<- paste( c("<",paste(roundsld[-c(1,nsld)],"-",sep=""),">")
                       , c(roundsld[-1],roundsld[nhgt]),"m"
                       , sep=""
                       )#end paste
   sldcols    <<- c(         "purple3",   "mediumpurple1",      "royalblue4"
                   ,      "steelblue3",     "deepskyblue",         "#004E00"
                   ,     "chartreuse3",      "olivedrab3", "lightgoldenrod3"
                   ,         "yellow3",     "darkorange1",            "red3"
                   ,      "firebrick4",        all.colour
                   )#end c
}else{
   cat(" In globdims.r:","\n")
   cat(" ISLD.TYPE = ",isld.type,"\n")
   stop(" Invalid ISLD.TYPE, it must be C or H (feel free to add more options.")
}#end if
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     Define default weighting factor for carbon balance.                                  #
#------------------------------------------------------------------------------------------#
if ("klight" %in% ls()){
   klight <<- klight
}else{
   klight <<- 0.8
}#end if 
#==========================================================================================#
#==========================================================================================#
