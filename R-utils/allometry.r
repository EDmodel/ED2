#==========================================================================================#
#==========================================================================================#
h2dbh = function(h,ipft){

   if (length(ipft) == 1){
     zpft = rep(ipft,times=length(h))
   }else{
     zpft = ipft
   }#end if

   tropo = pft$tropical[zpft] & (iallom == 0 | iallom == 1)
   tropn = pft$tropical[zpft] & iallom == 2
   tempe = ! pft$tropical[zpft]

   dbh = NA * h
   dbh[tropo] = exp((log(h[tropo])-pft$b1Ht[zpft[tropo]])/pft$b2Ht[zpft[tropo]])
   dbh[tropn] = ( log( pft$hgt.ref[zpft[tropn]] / ( pft$hgt.ref[zpft[tropn]] - h[tropn] ) )
                / pft$b1Ht[zpft[tropn]] ) ^ ( 1. / pft$b2Ht[zpft[tropn]])
   dbh[tempe] = log( 1.0 - ( h[tempe] - pft$hgt.ref[zpft[tempe]])
                   / pft$b1Ht[zpft[tempe]] ) / pft$b2Ht[zpft[tempe]]

   return(dbh)
}#end function h2dbh
#==========================================================================================#
#==========================================================================================#






#==========================================================================================!
#==========================================================================================!
dbh2h = function(ipft,dbh){

   if (length(ipft) == 1){
     zpft = rep(ipft,times=length(dbh))
   }else{
     zpft = ipft
   }#end if

   dbhuse        = dbh
   large         = is.finite(dbh) & dbh > pft$dbh.crit[zpft]
   dbhuse[large] = pft$dbh.crit[zpft[large]]

   tropo         = pft$tropical[zpft] & (iallom == 0 | iallom == 1)
   tropn         = pft$tropical[zpft] & iallom == 2
   tempe         = ! pft$tropical[zpft]

   h         = NA * dbh
   h[tropo]  = exp(pft$b1Ht[zpft[tropo]] + pft$b2Ht[zpft[tropo]] * log(dbhuse[tropo]) )
   h[tropn]  = ( pft$hgt.ref[zpft[tropn]] 
               * (1.0 - exp( -pft$b1Ht[zpft[tropn]] * dbhuse^pft$b2Ht[zpft[tropn]] ) ) )
   h[tempe]  = ( pft$hgt.ref[zpft[tempe]] + pft$b1Ht[zpft[tempe]] 
               * (1.0 - exp(pft$b2Ht[zpft[tempe]] * dbhuse[tempe] ) ) )

   return(h)
}#end function dbh2h
#==========================================================================================!
#==========================================================================================!






#==========================================================================================#
#==========================================================================================#
dbh2bl = function(dbh,ipft){

   if (length(ipft) == 1){
     zpft = rep(ipft,times=length(dbh))
   }else{
     zpft = ipft
   }#end if

   dbhuse       = dbh
   huge         = is.finite(dbh) & dbh > pft$dbh.crit[zpft]
   dbhuse[huge] = pft$dbh.crit[zpft[huge]]

   bleaf = pft$b1Bl [zpft] /C2B * dbhuse ^ pft$b2Bl [zpft] 

   return(bleaf)
}# end function dbh2bl
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
dbh2bd = function(dbh,ipft){
   if (length(ipft) == 1){
     zpft = rep(ipft,times=length(dbh))
   }else{
     zpft = ipft
   }#end if

   small = is.finite(dbh) & dbh <= pft$dbh.crit[zpft]
   large = is.finite(dbh) & dbh >  pft$dbh.crit[zpft]

   bdead = NA * dbh
   bdead[small] = ( pft$b1Bs.small[zpft[small]] / C2B * dbh[small] 
                  ^ pft$b2Bs.small[zpft[small]] )
   bdead[large] = ( pft$b1Bs.large[zpft[large]] / C2B * dbh[large] 
                  ^ pft$b2Bs.large[zpft[large]] )
   return(bdead)
}# end function dbh2bl
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#    Canopy Area allometry from Dietze and Clark (2008).                                   #
#------------------------------------------------------------------------------------------#
dbh2ca = function(dbh,ipft){
   if (length(ipft) == 1){
     zpft = rep(ipft,times=length(dbh))
   }else{
     zpft = ipft
   }#end if

   #----- Find local LAI, the minimum size for a crown area. ------------------------------#
   bleaf  = dbh2bl(dbh,ipft)
   loclai = pft$SLA[zpft] * bleaf
   
   dbhuse        = dbh
   large         = is.finite(dbh) & dbh > pft$dbh.crit[zpft]
   dbhuse[large] = pft$dbh.crit[zpft[large]]
   crown         = pft$b1Ca[zpft] * dbhuse ^ pft$b2Ca[zpft]

   #----- Local LAI / Crown area should never be less than one. ---------------------------#
   small        = crown < loclai
   crown[small] = loclai[small]

   return(crown)
}#end function dbh2ca
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#    Wood area index from Ahrends et al. (2010).                                           #
#------------------------------------------------------------------------------------------#
dbh2wai = function(dbh,ipft){
   if (length(ipft) == 1){
     zpft = rep(ipft,times=length(dbh))
   }else{
     zpft = ipft
   }#end if

   dbhuse        = dbh
   large         = is.finite(dbh) & dbh > pft$dbh.crit[zpft]
   dbhuse[large] = pft$dbh.crit[zpft[large]]
   wai           = pft$b1WAI[zpft] * dbhuse ^ pft$b2WAI[zpft]

   return(wai)
}#end function dbh2ca
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#    Standing volume of a tree.                                                            #
#------------------------------------------------------------------------------------------#
dbh2vol = function(hgt,dbh,ipft){
   vol  = pft$b1Vol[ipft] * hgt * dbh ^ pft$b2Vol[ipft]
   return(vol)
}#end function dbh2ca
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#    Rooting depth.                                                                        #
#------------------------------------------------------------------------------------------#
dbh2rd = function(hgt,dbh,ipft){
   if (iallom == 0){
      #------------------------------------------------------------------------------------#
      #    Original ED-2.1 (I don't know the source for this equation, though).            #
      #------------------------------------------------------------------------------------#
      vol  = dbh2vol(hgt,dbh,ipft)
      rd   = pft$b1Rd[ipft] * vol ^ pft$b2Rd[ipft]
   }else if (iallom == 1 || iallom == 2){
       #-----------------------------------------------------------------------------------#
       #    This is just a test allometry, that imposes root depth to be 0.5 m for         #
       # plants that are 0.15-m tall, and 5.0 m for plants that are 35-m tall.             #
       #-----------------------------------------------------------------------------------#
       rd = pft$b1Rd[ipft] * hgt ^ pft$b2Rd[ipft]
       #-----------------------------------------------------------------------------------#
   }#end if
   return(rd)
}#end function dbh2rd
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#    This function finds the trunk height.  Currently this is based on the following       #
# reference, which is for a site in Bolivia:                                               #
#                                                                                          #
# Poorter L., L. Bongers, F. Bongers, 2006: Architecture of 54 moist-forest tree           #
#     species: traits, trade-offs, and functional groups. Ecology, 87, 1289-1301.          #
#------------------------------------------------------------------------------------------#
h2crownbh = function (height,ipft){
   crown_length = pft$b1Cl[ipft] * height ^ pft$b2Cl[ipft]
   h2crownbh    = height - crown_length

   h2crownbh[is.finite(h2crowndbh) & h2crownbh < 0.05] = 0.05

   return(h2crownbh)
}#end function h2crownbh
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#    This function finds the leaf biomass for different plants.  This is based on the      #
# following papers:                                                                        #
#                                                                                          #
#  Cole, T. J., J. J. Ewel, 2006:  Allometric equations for four valuable tropical tree    #
#      species.  Forest Ecol. Management, 229, 351-360.                                    #
#                                                                                          #
#  Calvo-Alvarado, J. C., N. G. McDowell, R. H. Waring, 2008:  Allometric relationships    #
#      predicting foliar biomass and leaf area:sapwood area ratio from tree height in five #
#      Costa Rican rain forest species.  Tree Physiol. 28, 1601-1608.                      #
#------------------------------------------------------------------------------------------#
dbh2bl.alt = function (dbh,genus){
   #----- Make genus case insensitive. ----------------------------------------------------#
   genushere = tolower(genus)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Decide which equation to use based on the genus, or if it is to call ED-2.1, just #
   # call the good old dbh2bl...                                                           #
   #---------------------------------------------------------------------------------------#
   if (genushere == "cedrela"){
      h = dbh2h(3,dbh)
      x = dbh^2 * h
      
      small = is.finite(dbh) & dbh <= 10.
      large = is.finite(dbh) & dbh >  10.
      
      bleaf = NA * dbh
      bleaf[small] = 0.1265 / C2B * x[small] ^ 0.2787
      bleaf[large] = 0.0013 / C2B * x[large] ^ 0.9218

   }else if(genushere == "cordia"){
      h = dbh2h(2,dbh)
      x = dbh^2 * h
      
      small = is.finite(dbh) & dbh <= 5.
      large = is.finite(dbh) & dbh >  5.
      
      bleaf = NA * dbh
      bleaf[small] = 0.3041 / C2B * x[small] ^ 0.1082
      bleaf[large] = 0.0391 / C2B * x[large] ^ 0.5151


   }else if(genushere == "hyeronima"){
      h = dbh2h(4,dbh)
      x = dbh^2 * h
      
      small = is.finite(dbh) & dbh <= 10.
      large = is.finite(dbh) & dbh >  10.
      
      bleaf = NA * dbh
      bleaf[small] = 0.2144 / C2B * x[small] ^ 0.2852
      bleaf[large] = 0.0094 / C2B * x[large] ^ 0.6910

   }else if(genushere == "tetragastris"){
      bleaf = 0.040 / C2B * dbh ^ 1.737

   }else if(genushere == "virola"){
      bleaf = 0.002 / C2B * dbh ^ 2.468

   }else if(genushere == "carapa"){
      bleaf = 0.012 / C2B * dbh ^ 2.089

   }else if(genushere == "vochysia"){
      bleaf = 0.673 / C2B * dbh ^ 1.058

   }else if(genushere == "pentaclethra"){
      bleaf = 0.958 / C2B * dbh ^ 0.757

   }else if(genushere == "grass"){
      bleaf = dbh2bl(dbh,1)

   }else if(genushere == "early"){
      bleaf = dbh2bl(dbh,2)

   }else if(genushere == "mid"){
      bleaf = dbh2bl(dbh,3)

   }else if(genushere == "late"){
      bleaf = dbh2bl(dbh,4)

   }else{
      stop (paste("Genus ",genus," wasn't found.  ",
                 ,"Sorry, I can't find bleaf for this one...",sep=""))
   }#end if
   return(bleaf)
}#end function h2crownbh
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     We find the above-ground biomass based on dbh and wood density, following the        #
# allometry proposed by:                                                                   #
#                                                                                          #
# Baker, T. R., and co-authors, 2004: Variation in wood density determines spatial         #
#     patterns in Amazonian forest biomass. Global Change Biol., 10, 545-562.              #
#                                                                                          #
# Chave, J., B. Riera, M.A. Dubois, 2001: Estimation of biomass in a neotropical forest of #
#     French Guiana: spatial and temporal variability.  J. Trop. Ecol., 17, 79-96.         #
#------------------------------------------------------------------------------------------#
dbh2agb.baker <<- function(dbh,wdens,allom="baker.chave"){


   ln.dbh = log(dbh)


   if (allom == "baker.chave"){
      #------ Use Chave's based function (equation 2, Baker et al., 2004). ----------------#
      agb = wdens / 0.58 / C2B * exp(2.42 * ln.dbh - 2.00)
      #------------------------------------------------------------------------------------#
   }else if (allom == "baker.chambers"){
      #------ Use Chambers' function. -----------------------------------------------------#
      agb = ( wdens / 0.67 / C2B 
            * exp(0.33 * ln.dbh + 0.933 * ln.dbh^2 - 0.122 * ln.dbh^3 - 0.37) )
      #------------------------------------------------------------------------------------#
   }else if (allom == "chave.2006"){
      #------ Use Chambers' function. -----------------------------------------------------#
      agb = ( wdens / C2B 
            * exp( -1.499 + 2.1481 * ln.dbh + 0.207 * ln.dbh^2 - 0.0281 * ln.dbh^3 ) )
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#
   return(agb)
}#end if
#==========================================================================================#
#==========================================================================================#
