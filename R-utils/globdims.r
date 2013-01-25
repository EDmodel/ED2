#----- Some dimensions based on ED-2.2 default. -------------------------------------------#
npft       <<- 17 # Number of plant functional types
nlu        <<-  3 # Number of land use types.
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
#------------------------------------------------------------------------------------------#


#----- Minimum R2 that we consider meaningful. --------------------------------------------#
r2.min        <<- 0.36
#----- Typical p-value below which we reject the NULL hypothesis. -------------------------#
pval.max      <<- 0.05
#----- p-value below which we reject the NULL hypothesis for u* filters. ------------------#
pval.ustar    <<- 0.05
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Tolerance for root-finding methods.                                                  #
#------------------------------------------------------------------------------------------#
toler  <<- 1.e-6  # Tolerance for the root-finder algorithms
maxfpo <<- 60     # Maximum number of Regula-Falsi iterations
maxit  <<- 150    # Maximum number of iterations in general
#------------------------------------------------------------------------------------------#




#----- Define some default legend colours and names. --------------------------------------#
lunames   <<- c("Agricultural","Secondary","Primary","Total")
lucols    <<- c("goldenrod","chartreuse","darkgreen",all.colour)

distnames <<- c("Agr->Agr" ,"2nd->Agr" ,"Prim->Agr"
               ,"Agr->2nd" ,"2nd->2nd" ,"Prim->2nd"
               ,"Agr->Prim","2nd->Prim","Prim->Prim")
distcols  <<- c("gold","darkorange2","firebrick"
               ,"lightskyblue","turquoise","steelblue"
               ,"palegreen","chartreuse","forestgreen")
#------------------------------------------------------------------------------------------#


#----- Growth respiration factor (to estimate when the actual variable is not there). -----#
growth.resp.fac <<- c(rep(0.333,times=5),rep(0.4503,times=3),rep(0,times=3)
                     ,rep(0.333,times=4))
#------------------------------------------------------------------------------------------#



#----- fswh is the FSW that plants experience and they are happy (wilting point = 0). -----#
fswh <<- 0.99
#------------------------------------------------------------------------------------------#


#----- Years for which we have eddy flux tower. -------------------------------------------#
eft.year <<- c(1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012)
eft.pch  <<- c(   0,   1,   6,   3,   4,   5,  13,  15,  16,  17,  18,   9,   8,   7,  14)
eft.col  <<- c("purple4","purple1","mediumpurple1","royalblue4","steelblue","deepskyblue"
              ,"#004000","chartreuse4","lawngreen","red3","#FF6000","#FFA000"
              ,"chocolate4","lightgoldenrod3","yellow2")
#------------------------------------------------------------------------------------------#



#----- Standard colours and names for soil classes. ---------------------------------------#
stext.cols  <<- c("gold","chartreuse","limegreen","darkgreen","purple3"
                 ,"deepskyblue","aquamarine","midnightblue","darkorange3","sienna"
                 ,"firebrick","grey61","grey29","orchid","olivedrab","goldenrod"
                 ,"steelblue")
stext.names <<- c("Sand","Loamy Sand","Sandy loam","Silt loam","Loam","Sandy clay loam"
                 ,"Silty clay loam","Clayey loam","Sandy clay","Silty clay","Clay","Peat"
                 ,"Bedrock","Silt","Heavy clay","Clayey sand","Clayey silt")
stext.acron <<- c("Sa","LoSa","SaLo","SiLo","Lo","SaClLo","SiClLo","ClLo"
                 ,"SaCl","SiCl","Cl","Pe","Br","Si","HCl","ClSa","ClSi")
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
   dbhkeys    <<- paste(classdbh,"-",c(classdbh[-1],Inf),sep="")
   dbhnames   <<- c(paste("> ",sprintf("%.1f",classdbh),"cm",sep=""),"Total")
   dbhcols    <<- c(         "purple3",   "mediumpurple1",      "royalblue4"
                   ,      "steelblue3",     "deepskyblue",         "#004E00"
                   ,     "chartreuse3",      "olivedrab3", "lightgoldenrod3"
                   ,         "yellow3",     "darkorange1",       "firebrick"
                   ,        all.colour
                   )#end c

}else if (idbh.type == 2){
   ndbh       <<-  6
   classdbh   <<- c(0,10,20,35,50,70)
   breakdbh   <<- c(-Inf,classdbh[-1],Inf)
   dbhkeys    <<- paste(classdbh,"-",c(classdbh[-1],Inf),sep="")
   dbhnames   <<- c(paste("> ",sprintf("%.1f",classdbh),"cm",sep=""),"Total")
   dbhcols    <<- c(         "purple3",      "royalblue3",     "chartreuse3"
                  ,          "yellow3",     "darkorange1",       "firebrick"
                  ,         all.colour
                  )#end c
}else if (idbh.type == 3){
   ndbh       <<-  4
   classdbh   <<- c(0,10,35,55)
   breakdbh   <<- c(-Inf,classdbh[-1],Inf)
   dbhkeys    <<- paste(classdbh,"-",c(classdbh[-1],Inf),sep="")
   dbhnames   <<- c(paste("> ",sprintf("%.1f",classdbh),"cm",sep=""),"Total")
   dbhcols    <<- c(      "royalblue3",     "chartreuse3"
                  ,          "yellow3",     "darkorange1",         all.colour
                  )#end c
}else{
   cat(" In globdims.r:","\n")
   cat(" IDBH.TYPE = ",idbh.type,"\n")
   stop(" Invalid IDBH.TYPE, it must be either 1, 2, or 3 (feel free to add more options.")
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
