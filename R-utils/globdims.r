#----- Some dimensions based on ED-2.2 default. -------------------------------------------#
npft       <<- 17 # Number of plant functional types
nlu        <<-  8 # Number of land use types.
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


#----- It should be defined by now. -------------------------------------------------------#
if (! "all.colour" %in% ls()) all.colour    = "grey22"
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
lunames   <<- c("Pasture","Plantation","Tree fall","Burnt","Abandoned","Logging"
               ,"Skid trail","Cropland","Total")
lucols    <<- c("#E65C17","#2996CC","#306614","#990F0F","#2996CC","#A3CC52"
               ,"#B49ED2","#F5C858",all.colour)
distnames <<- c("PST->PST" ,"FPL->PST" ,"TFL->PST","BRN->PST","ABN->PST","LOG->PST","SKD->PST","CPL->PST"
               ,"PST->FPL" ,"FPL->FPL" ,"TFL->FPL","BRN->FPL","ABN->FPL","LOG->FPL","SKD->FPL","CPL->FPL"
               ,"PST->TFL" ,"FPL->TFL" ,"TFL->TFL","BRN->TFL","ABN->TFL","LOG->TFL","SKD->TFL","CPL->TFL"
               ,"PST->BRN" ,"FPL->BRN" ,"TFL->BRN","BRN->BRN","ABN->BRN","LOG->BRN","SKD->BRN","CPL->BRN"
               ,"PST->ABN" ,"FPL->ABN" ,"TFL->ABN","BRN->ABN","ABN->ABN","LOG->ABN","SKD->ABN","CPL->ABN"
               ,"PST->LOG" ,"FPL->LOG" ,"TFL->LOG","BRN->LOG","ABN->LOG","LOG->LOG","SKD->LOG","CPL->LOG"
               ,"PST->SKD" ,"FPL->SKD" ,"TFL->SKD","BRN->SKD","ABN->SKD","LOG->SKD","SKD->SKD","CPL->SKD"
               ,"PST->CPL" ,"FPL->CPL" ,"TFL->CPL","BRN->CPL","ABN->CPL","LOG->CPL","SKD->CPL","CPL->CPL"
               )
distcols  <<- c("#500C22","#FF99A3","#6B0020","#FF7780","#8C001B","#FF7168","#B30010","#CD2A0A"
               ,"#5B1500","#FF9F7F","#FF7F4F","#AD4000","#FFD4B6","#904E00","#BA7100","#683D00"
               ,"#FFB055","#FFC27C","#FFB330","#F3B700","#907400","#4F3F00","#F5E074","#B4BA00"
               ,"#E3E39E","#89A900","#607800","#B9EDA7","#007218","#35CE4D","#8BF790","#002F09"
               ,"#01B259","#C5E8CE","#009357","#00C377","#006E45","#00C597","#008072","#005654"
               ,"#71F2F5","#018CD1","#00466D","#53A1FF","#00467E","#9CAAFF","#0E3EC4","#00207B"
               ,"#D7C5FF","#311B53","#5B008D","#46006D","#E27DFF","#F5B6FF","#F49AFF","#381E3A"
               ,"#52004E","#FF99E2","#4E053E","#DC1BA7","#92005E","#FFA3CD","#FF47A8","#FFAEC6")
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
#     Define necromass pools.                                                              #
#------------------------------------------------------------------------------------------#
nnecro     = 7
necrokeys  = c("PSC","SSC","STSC","MSC","FSC","STGC","FGC","ALL")
necronames = c("Passive","Humified","BG woody debris","Microbial"
              ,"BG litter","AG woody debris","AG litter","Total")
necrocols  = c("#811F9E","#1BA2F7","#880D32","#CCCA3D","#107C92","#F87856","#2BD2DB"
              ,all.colour)
necroltys  = c("dotted","dotdash","twodash","dashed","longdash","dotdash","twodash","solid")
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
   dbhkeys    <<- paste0(classdbh,"-",c(classdbh[-1],Inf))
   dbhnames   <<- paste0( c("<",paste(classdbh[-c(1,ndbh)],"-",sep=""),">")
                        , c(classdbh[-1],classdbh[ndbh]),"cm"
                        )#end paste0
   dbhcols    <<- c(         "purple3",   "mediumpurple1",      "royalblue4"
                   ,      "steelblue3",     "deepskyblue",     "chartreuse3"
                   , "lightgoldenrod3",         "yellow3",     "darkorange1"
                   ,       "firebrick",        all.colour
                   )#end c
   dbhltys    <<- c("twodash","dashed","longdash"
                   ,"dotdash","twodash","dashed"
                   ,"longdash","dotdash","twodash"
                   ,"dashed","solid")
}else if (idbh.type == 2){
   ndbh       <<-  6
   classdbh   <<- c(0,10,20,35,55,80)
   breakdbh   <<- c(-Inf,classdbh[-1],Inf)
   dbhlabel   <<- "06_szclss"
   dbhkeys    <<- paste0(classdbh,"-",c(classdbh[-1],Inf))
   dbhnames   <<- paste0( c("<",paste(classdbh[-c(1,ndbh)],"-",sep=""),">")
                        , c(classdbh[-1],classdbh[ndbh]),"cm"
                        )#end paste0
   dbhcols    <<- c("#811F9E","#107C92","#1BA2F7","#2BD2DB","#CCCA3D","#F87856",all.colour)
   dbhltys    <<- c("twodash","dashed","longdash","dotdash","twodash","longdash","solid")
}else if (idbh.type == 3){
   ndbh       <<-  4
   classdbh   <<- c(0,10,35,55)
   dbhlabel   <<- "04_szclss"
   breakdbh   <<- c(-Inf,classdbh[-1],Inf)
   dbhkeys    <<- paste0(classdbh,"-",c(classdbh[-1],Inf))
   dbhnames   <<- paste0( c("<",paste(classdbh[-c(1,ndbh)],"-",sep=""),">")
                        , c(classdbh[-1],classdbh[ndbh]),"cm"
                        )#end paste0
   dbhcols    <<- c("#811F9E","#1BA2F7","#CCCA3D","#F87856",all.colour)
   dbhltys    <<- c("twodash","dashed","longdash","dotdash","solid")
}else if (idbh.type == 4){
   ndbh       <<-  5
   classdbh   <<- c(0,10,30,50,80)
   breakdbh   <<- c(-Inf,classdbh[-1],Inf)
   dbhlabel   <<- "05_szclss"
   dbhkeys    <<- paste0(classdbh,"-",c(classdbh[-1],Inf))
   dbhnames   <<- paste0( c("<",paste(classdbh[-c(1,ndbh)],"-",sep=""),">")
                        , c(classdbh[-1],classdbh[ndbh]),"cm"
                        )#end paste0
#   dbhcols    <<- c("#811F9E","#1BA2F7","#2BD2DB","#CCCA3D","#F87856",all.colour)
   dbhcols    <<- c("#10002A","#3F046B","#A12663","#F87D00","#FCC200",all.colour)
   dbhltys    <<- c("twodash","dashed","longdash","dotdash","twodash","solid")
}else if (idbh.type == 5){
   ndbh       <<-  8
   classdbh   <<- c(0,20,40,60,80,100,120,140)
   breakdbh   <<- c(-Inf,classdbh[-1],Inf)
   dbhlabel   <<- "08_szclss"
   dbhkeys    <<- paste0(classdbh,"-",c(classdbh[-1],Inf))
   dbhnames   <<- paste0( c("<",paste(classdbh[-c(1,ndbh)],"-",sep=""),">")
                        , c(classdbh[-1],classdbh[ndbh]),"cm"
                        )#end paste0
   dbhcols    <<- c("#811F9E","#107C92","#1BA2F7","#2BD2DB","#F9E5C0"
                   ,"#CCCA3D","#F87856","#880D32",all.colour)
   dbhltys    <<- c("twodash","dashed"  ,"longdash","dotdash","twodash"
                   ,"dashed" ,"longdash","dotdash" ,"solid")
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
   classhgt   <<- c(0,2,4,7,10,14,18,22,26,30,34,38,42)
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
}else if (ihgt.type == 3){
   nhgt       <<-  12
   classhgt   <<- c(0,5,10,15,20,24,28,31,34,36,38,40)
   hgtlabel   <<- "12_szclss"
   breakhgt   <<- c(-Inf,classhgt[-1],Inf)
   hgtkeys    <<- paste(classhgt,"-",c(classhgt[-1],Inf),sep="")
   hgtnames   <<- paste( c("<",paste(classhgt[-c(1,nhgt)],"-",sep=""),">")
                       , c(classhgt[-1],classhgt[nhgt]),"cm"
                       , sep=""
                       )#end paste
   hgtcols    <<- c(         "purple3",   "mediumpurple1",      "royalblue4"
                   ,      "steelblue3",     "deepskyblue",     "chartreuse3"
                   ,      "olivedrab3", "lightgoldenrod3",         "yellow3"
                   ,     "darkorange1",            "red3",      "firebrick4"
                   ,        all.colour
                   )#end c
}else{
   cat(" In globdims.r:","\n")
   cat(" IHGT.TYPE = ",ihgt.type,"\n")
   stop(" Invalid IHGT.TYPE, it must be between 1 and 3 (feel free to add more options.")
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





#==========================================================================================#
#==========================================================================================#
#     Default above ground fraction and C:N ratios for necromass pools.                    #
#------------------------------------------------------------------------------------------#
agf.fsc        <<- 0.5
agf.stsc       <<- 0.7
c2n.fast       <<- 29.24
c2n.structural <<- 150.
c2n.slow       <<- 10.0
litter.2.fsc   <<- 0.32   # Brechet et al. (2017) Ecosystems
fwd.2.cwd      <<- 0.125  # Fine:coarse woody debris ratio, based on 
                          #    Berenguer et al. (2015) PLoS ONE
#==========================================================================================#
#==========================================================================================#
