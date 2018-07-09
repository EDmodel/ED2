#==========================================================================================#
#==========================================================================================#
h2dbh <<- function(h,ipft){

   if (length(ipft) == 1){
     zpft = rep(ipft,times=length(h))
   }else{
     zpft = ipft
   }#end if

   tropo = pft$tropical[zpft] & iallom %in% c(0,1)
   tropn = pft$tropical[zpft] & iallom %in% c(2,3)
   tempe = ! pft$tropical[zpft]


   huse    = pmin(pft$hgt.max[zpft],h)

   hgt.ref = pft$hgt.ref[zpft]
   b1Ht    = pft$b1Ht   [zpft]
   b2Ht    = pft$b2Ht   [zpft]

   dbh = ifelse( test = tempe
               , yes  = log( 1.0 - ( huse - hgt.ref) / b1Ht ) / b2Ht
               , no   = ifelse ( test = tropn
                               , yes  = (log(hgt.ref / (hgt.ref-huse))/b1Ht ) ^ (1./b2Ht)
                               , no   = exp((log(huse)-b1Ht)/b2Ht)
                               )#end ifelse
               )#end ifelse

   return(dbh)
}#end function h2dbh
#==========================================================================================#
#==========================================================================================#






#==========================================================================================!
#==========================================================================================!
dbh2h <<- function(dbh,ipft,use.crit=TRUE){

   if (length(ipft) == 1){
     zpft = rep(ipft,times=length(dbh))
   }else{
     zpft = ipft
   }#end if


   if (use.crit){
      dbhuse = pmin(pft$dbh.crit[zpft],dbh) + 0. * dbh
   }else{
      dbhuse = dbh
   }#end if (use.crit)

   tropo         = pft$tropical[zpft] & iallom %in% c(0,1)
   tropn         = pft$tropical[zpft] & iallom %in% c(2,3)
   tempe         = ! pft$tropical[zpft]

   hgt.ref = pft$hgt.ref[zpft]
   b1Ht    = pft$b1Ht   [zpft]
   b2Ht    = pft$b2Ht   [zpft]

   h = ifelse( test = tempe
             , yes  = hgt.ref + b1Ht * (1. - exp(b2Ht * dbhuse) )
             , no   = ifelse( test = tropn
                            , yes  = hgt.ref * (1.-exp(-b1Ht * dbhuse^b2Ht ) )
                            , no   = exp(b1Ht + b2Ht * log(dbhuse) )
                            )#end ifelse
             )#end ifelse

   return(h)
}#end function dbh2h
#==========================================================================================!
#==========================================================================================!






#==========================================================================================#
#==========================================================================================#
size2bl <<- function(dbh,hgt,ipft,use.crit=TRUE){
   #----- Make sure that the PFT variable has the same length as dbh. ---------------------#
   if (length(ipft) == 1){
     zpft = rep(ipft,times=length(dbh))
   }else{
     zpft = ipft
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Limit dbh to dbh.crit. ----------------------------------------------------------#
   if (use.crit){
      dbhuse  = pmin(dbh,pft$dbh.crit[zpft]) + 0. * dbh
   }else{
      dbhuse  = dbh
   }#end if (use.crit)
   #---------------------------------------------------------------------------------------#


   #----- Decide which variable to use as dependent variable (DBH or DBH^2*Hgt). ----------#
   size     = ifelse( test = pft$tropical[zpft] & (! pft$liana[zpft]) & (iallom %in% 3)
                    , yes  = dbhuse * dbhuse * hgt
                    , no   = dbhuse
                    )#end ifelse
   #---------------------------------------------------------------------------------------#



   #----- Find leaf biomass. --------------------------------------------------------------#
   bleaf = pft$b1Bl[zpft] / C2B * size ^ pft$b2Bl[zpft]
   #---------------------------------------------------------------------------------------#


   return(bleaf)
}# end function size2bl
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
size2bd <<- function(dbh,hgt,ipft){
   #----- Make sure that the PFT variable has the same length as dbh. ---------------------#
   if (length(ipft) == 1){
     zpft = rep(ipft,times=length(dbh))
   }else{
     zpft = ipft
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Decide which variable to use as dependent variable (DBH or DBH^2*Hgt). ----------#
   size     = ifelse( test = pft$tropical[zpft] & (! pft$liana[zpft]) & (iallom %in% 3)
                    , yes  = dbh * dbh * dbh2h(dbh,ipft=zpft)
                    , no   = dbh
                    )#end ifelse
   #---------------------------------------------------------------------------------------#



   #----- Select allometric parameters based on the size. ---------------------------------#
   bdead = ifelse( test = dbh %<% pft$dbh.crit[zpft]
                 , yes  = pft$b1Bs.small[zpft] / C2B * size ^ pft$b2Bs.small[zpft]
                 , no   = pft$b1Bs.large[zpft] / C2B * size ^ pft$b2Bs.large[zpft]
                 )#end ifelse
   #---------------------------------------------------------------------------------------#

   return(bdead)
}# end function size2bd
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
size2bw <<- function(dbh,hgt,ipft,use.crit=TRUE){
   if (length(ipft) == 1){
     zpft = rep(ipft,times=length(dbh))
   }else{
     zpft = ipft
   }#end if


   bdead  = size2bd(dbh=dbh,hgt=hgt,ipft=zpft)
   bleaf  = size2bl(dbh=dbh,hgt=hgt,ipft=zpft,use.crit=use.crit)
   bsapw  = pft$qsw  [ipft] * hgt * bleaf
   bbark  = pft$qbark[ipft] * hgt * bleaf
   bwood  = bdead + bsapw + bbark

   return(bwood)
}# end function size2bw
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#    This function finds the equivalent on-allometry DBH given observed height and dbh.    #
#------------------------------------------------------------------------------------------#
size2de <<- function(dbh,hgt,ipft,dbh.by=0.1,...){
   #----- Make sure that the PFT variable has the same length as dbh. ---------------------#
   if (length(ipft) == 1){
     zpft = rep(ipft,times=length(dbh))
   }else{
     zpft = ipft
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Identify cohorts with size that outside resolvable height range. ----------------#
   large = (dbh*dbh*hgt) %>=% (pft$dbh.crit[zpft]*pft$dbh.crit[zpft]*pft$hgt.max[zpft])
   small = (dbh*dbh*hgt) %<=% (pft$dbh.min [zpft]*pft$dbh.min [zpft]*pft$hgt.min[zpft])
   heq   = ifelse( test = large
                 , yes  = pft$hgt.max[zpft]
                 , no   = ifelse(test = small,yes=pft$hgt.min[zpft],no=NA_real_)
                 )#end ifelse
   #---------------------------------------------------------------------------------------#


   #----- First solve the equivalent dbh for the simplest cases. --------------------------#
   dbheq = sqrt(hgt/heq) * dbh
   #---------------------------------------------------------------------------------------#


   #----- Find "size" (dbh^2 * hgt). ------------------------------------------------------#
   size = dbh * dbh * hgt
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Go through each PFT type that must be filled. , then use interpolation to find    #
   # the best match.                                                                       #
   #---------------------------------------------------------------------------------------#
   loop.pft = sort(unique(zpft[is.na(dbheq)]))
   for (wpft in loop.pft){
       dbh.lut  = seq(from=pft$dbh.min[wpft],to=pft$dbh.crit[wpft],by=dbh.by)
       hgt.lut  = dbh2h(dbh=dbh.lut,ipft=wpft)
       size.lut = dbh.lut*dbh.lut*hgt.lut

       sel        = is.na(dbheq) & (zpft %in% wpft)
       idx        = mapply(FUN=which.closest,x=size[sel], MoreArgs=list(A=size.lut))
       dbheq[sel] = dbh.lut[idx]
   }#end for (wpft in loop.pft)
   #---------------------------------------------------------------------------------------#


   return(dbheq)
}# end function size2de
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#    Canopy Area allometry from Dietze and Clark (2008).                                   #
#------------------------------------------------------------------------------------------#
size2ca <<- function(dbh,hgt,ipft,use.crit=TRUE){
   if (length(ipft) == 1){
     zpft = rep(ipft,times=length(dbh))
   }else{
     zpft = ipft
   }#end if

   #----- Find local LAI, the minimum size for a crown area. ------------------------------#
   loclai = pft$SLA[zpft] * size2bl(dbh,hgt,ipft)
   #---------------------------------------------------------------------------------------#


   #----- Limit dbh to dbh.crit. ----------------------------------------------------------#
   if (use.crit){
      dbhuse  = pmin(dbh,pft$dbh.crit[zpft]) + 0. * dbh
   }else{
      dbhuse  = dbh
   }#end if (use.crit)
   #---------------------------------------------------------------------------------------#


   #----- Decide which variable to use as dependent variable (DBH or DBH^2*Hgt). ----------#
   size     = ifelse( test = pft$tropical[zpft] & (! pft$liana[zpft]) & (iallom %in% 3)
                    , yes  = dbhuse * dbhuse * hgt
                    , no   = dbhuse
                    )#end ifelse
   #---------------------------------------------------------------------------------------#



   #----- Select allometric parameters based on the size. ---------------------------------#
   crown = pmin(loclai,pft$b1Ca[zpft] * size ^ pft$b2Ca[zpft])
   #---------------------------------------------------------------------------------------#

   return(crown)
}#end function size2ca
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#    Wood area index.                                                                      #
#------------------------------------------------------------------------------------------#
size2wai <<- function(dbh,hgt,ipft,use.crit=TRUE){
   #----- Make sure the size of variable ipft matches the number of dbh entries. ----------#
   if (length(ipft) == 1){
     zpft = rep(ipft,times=length(dbh))
   }else{
     zpft = ipft
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Cap dbh size to not exceed dbh.crit. --------------------------------------------#
   if (use.crit){
      dbhuse  = pmin(dbh,pft$dbh.crit[zpft]) + 0. * dbh
   }else{
      dbhuse  = dbh
   }#end if (use.crit)
   #---------------------------------------------------------------------------------------#



   #----- Decide which variable to use as dependent variable (DBH or DBH^2*Hgt). ----------#
   size     = ifelse( test = pft$tropical[zpft] & (! pft$liana[zpft]) & (iallom %in% 3)
                    , yes  = dbhuse * dbhuse * hgt
                    , no   = dbhuse
                    )#end ifelse
   #---------------------------------------------------------------------------------------#



   #----- Find the wood area index. -------------------------------------------------------#
   wai      = pft$b1WAI[zpft] * size ^ pft$b2WAI[zpft]
   #---------------------------------------------------------------------------------------#
   

   return(wai)
}#end function dbh2wai
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#    Standing volume of a tree.                                                            #
#------------------------------------------------------------------------------------------#
size2vol <<- function(hgt,dbh,ipft){
   vol  = pft$b1Vol[ipft] * ( hgt * dbh * dbh ) ^ pft$b2Vol[ipft]
   return(vol)
}#end function dbh2ca
#==========================================================================================#
#==========================================================================================#







#==========================================================================================#
#==========================================================================================#
#    Above-ground biomass.                                                                 #
#------------------------------------------------------------------------------------------#
ed.biomass <<- function(dbh,ipft,use.crit=TRUE){
   if (length(ipft) == 1){
     zpft = rep(ipft,times=length(dbh))
   }else{
     zpft = ipft
   }#end if

   hgt    = dbh2h(dbh=dbh,ipft=zpft)
   bleaf  = size2bl(dbh=dbh,hgt=hgt,ipft=zpft,use.crit=use.crit)
   bsapa  = pft$agf.bs[zpft] * pft$qsw  [zpft] * hgt * bleaf
   bbarka = pft$agf.bs[zpft] * pft$qbark[zpft] * hgt * bleaf
   bdeada = pft$agf.bs[zpft] * size2bd(dbh=dbh,hgt=hgt,ipft=zpft)
   agb    = bleaf + bsapa + bbarka + bdeada
   return(agb)
}#end function ed.biomass
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#    This function computes the commercial timber biomass for different PFTs.              #
#------------------------------------------------------------------------------------------#
size2bt <<- function(dbh,hgt,ipft){
   if (length(ipft) == 1){
     zpft = rep(ipft,times=length(dbh))
   }else{
     zpft = ipft
   }#end if

   #----- Find bdead and bsapwooda. -------------------------------------------------------#
   bdead    = size2bd(dbh=dbh,hgt=hgt,ipft=zpft)
   bleaf    = size2bl(dbh=dbh,hgt=hgt,ipft=zpft)
   bsapwood = pft$agf.bs[zpft] * pft$qsw  [zpft] * hgt * bleaf
   bbark    = pft$agf.bs[zpft] * pft$qbark[zpft] * hgt * bleaf
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Find tropical timber and above-ground biomass.                                   #
   #---------------------------------------------------------------------------------------#
   btimber = 1000. * pft$rho[zpft] * size2vol(dbh=dbh,hgt=hgt,ipft=zpft) / C2B
   bagwood = pft$agf.bs[zpft] * (bdead + bsapwood + bbark)
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #      Calculate commercial timber biomass based on the life form and habitat.          #
   #---------------------------------------------------------------------------------------#
   ans = ifelse( test = pft$grass[zpft]
               , yes  = 0.
               , no   = ifelse( test = pft$tropical[zpft]
                              , yes  = pmin(btimber,bagwood)
                              , no   = bagwood
                              )#end ifelse
               )#end ifelse
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end function size2bt
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#    Rooting depth.                                                                        #
#------------------------------------------------------------------------------------------#
size2rd <<- function(hgt,dbh,ipft){
   if (length(ipft) == 1){
     zpft = rep(ipft,times=length(dbh))
   }else{
     zpft = ipft
   }#end if




   if (iallom %in% c(0)){
      #------------------------------------------------------------------------------------#
      #    Original ED-2.1 (I don't know the source for this equation, though).            #
      #------------------------------------------------------------------------------------#
      vol  = size2vol(hgt=hgt,dbh=dbh,ipft=zpft)
      rd   = pft$b1Rd[zpft] * (hgt * dbh * dbh) ^ pft$b2Rd[zpft]
   }else if (iallom %in% c(1,2)){
      #------------------------------------------------------------------------------------#
      #    This is just a test allometry, that imposes root depth to be 0.5 m for          #
      # plants that are 0.15-m tall, and 5.0 m for plants that are 35-m tall.              #
      #------------------------------------------------------------------------------------#
      rd = pft$b1Rd[zpft] * hgt ^ pft$b2Rd[zpft]
      #------------------------------------------------------------------------------------#
   }else{
      #------------------------------------------------------------------------------------#
      #    For tropical trees, use the allometric model to obtain the Effective            #
      # Functional Rooting Depth based on B18.  We made a slight modification in their     #
      # equation relating delta 18O and depth:                                             #
      #                                                                                    #
      #    depth = exp(a + b * d18O^2)                                                     #
      #                                                                                    #
      # because it fits the data better than the original equation without the square,     #
      # and it avoids extremely shallow soils for small trees.  We also use a              #
      # heteroscedastic least squares, using the algorithm developed by L16.               #
      #                                                                                    #
      # References:                                                                        #
      #                                                                                    #
      # Brum M, Vadeboncoeur MA, Ivanov V, Asbjornsen H, Saleska S, Alves LF, Penha D,     #
      #    Dias JD, Aragao LEOC, Barros F, Bittencourt P, Pereira L, Oliveira RS, 2018.    #
      #    Hydrological niche segregation defines forest structure and drought             #
      #    tolerance strategies in a seasonal Amazonian forest. J. Ecol., in press.        #
      #    doi:10.1111/1365-2745.13022 (B18).                                              #
      #                                                                                    #
      # Longo M, Keller M, dos-Santos MN, Leitold V, Pinage ER, Baccini A, Saatchi S,      #
      #    Nogueira EM, Batistella M , Morton DC. 2016. Aboveground biomass variability    #
      #    across intact and degraded forests in the Brazilian Amazon.                     #
      #    Global Biogeochem. Cycles, 30(11):1639-1660. doi:10.1002/2016GB005465 (L16).    #
      #------------------------------------------------------------------------------------#
      dbhuse = pmin(dbh,pft$dbh.crit[zpft])
      d18O   = pft$d18O.ref[zpft] * (1. - exp(-pft$b1d18O[zpft] * dbhuse^pft$b2d18O[zpft]))
      trtree = pft$tropical[zpft] & (! pft$liana[zpft]) & (! pft$grass[zpft])


      #------------------------------------------------------------------------------------#
      #    For tropical trees, use an allometry loosely based on B18.  The original        #
      # approach by B18 yields extremely shallow roots for small- and medium-sized trees,  #
      # causing excessive water stress and poor agreement with GPP estimated from towers.  #
      #------------------------------------------------------------------------------------#
      size   = ifelse( test = pft$tropical[zpft] & (! pft$liana[zpft])
                     , yes  = dbhuse * dbhuse * hgt
                     , no   = hgt
                     )#end ifelse
      rd     = ifelse( test = trtree & use.efrd.trtree
                     , yes  = -exp(pft$b1Efrd[zpft] + pft$b2Efrd[zpft] * d18O * d18O)
                     , no   = pft$b1Rd[zpft] * size ^ pft$b2Rd[zpft]
                     )#end ifelse
      #------------------------------------------------------------------------------------#
   }#end if
   return(rd)
}#end function size2rd
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#    This function finds the trunk height.                                                 #
#------------------------------------------------------------------------------------------#
h2crownbh <<- function (height,ipft){
   crown.length = pft$b1Cl[ipft] * height ^ pft$b2Cl[ipft]
   ans          = pmax(0.05,height - crown.length)
   return(ans)
}#end function h2crownbh
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#    This function finds the crown length.                                                 #
#------------------------------------------------------------------------------------------#
h2cl <<- function (height,ipft){
   ans = pft$b1Cl[ipft] * height ^ pft$b2Cl[ipft]
   return(ans)
}#end function h2cl
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
dbh2bl.alt <<- function (dbh,genus){
   #----- Make genus case insensitive. ----------------------------------------------------#
   genushere = tolower(genus)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Decide which equation to use based on the genus, or if it is to call ED-2.1, just #
   # call size2bl...                                                                       #
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


   }else if(genushere %in% c("hieronima","hyeronima","hieronyma")){
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
      hgt   = dbh2h  (dbh=dbh,ipft=1)
      bleaf = size2bl(dbh=dbh,hgt=hgt,ipft=1)

   }else if(genushere == "early"){
      hgt   = dbh2h  (dbh=dbh,ipft=2)
      bleaf = size2bl(dbh=dbh,hgt=hgt,ipft=2)

   }else if(genushere == "mid"){
      hgt   = dbh2h  (dbh=dbh,ipft=3)
      bleaf = size2bl(dbh=dbh,hgt=hgt,ipft=3)

   }else if(genushere == "late"){
      hgt   = dbh2h  (dbh=dbh,ipft=4)
      bleaf = size2bl(dbh=dbh,hgt=hgt,ipft=4)

   }else{
      stop (paste("Genus ",genus," wasn't found.  ",
                 ,"Sorry, I can't find bleaf for this one...",sep=""))
   }#end if
   return(bleaf)
}#end function dbh2bl.alt
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
}#end function dbh2agb.baker
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     Biomass allometry based on Chave et al. (2014).  Results are always in kgC/plant.    #
#                                                                                          #
# References:                                                                              #
#                                                                                          #
# Chave, J., and co-authors, 2014: Improved allometric models to estimate tha aboveground  #
#     biomass of tropical trees.  Glob. Change Biol., 20, 3177-3190                        #
#     doi:10.1111/gcb.12629                                                                #
#                                                                                          #
#                                                                                          #
#                                                                                          #
# Input:                                                                                   #
# ---------------------------------------------------------------------------------------- #
# dbh        --- Diameter at breast height [cm]                                            #
# height     --- Height [m]                                                                #
# wdens      --- Wood density [g/cm3]                                                      #
# ---------------------------------------------------------------------------------------- #
#------------------------------------------------------------------------------------------#
size2agb.chave <<- function(dbh,height,wdens){
   #---------------------------------------------------------------------------------------#
   #     Make sure all terms have the same length.                                         #
   #---------------------------------------------------------------------------------------#
   lens = unique(c(length(dbh),length(height),length(wdens)))
   if ( length(lens) != 1 ){
      cat0("-----------------------------------------------------------")
      cat0("   Variables don't have the same length."                   )
      cat0("   DBH    = ",length(dbh)                                   )
      cat0("   HEIGHT = ",length(height)                                )
      cat0("   WDENS  = ",length(wdens)                                 )
      cat0("-----------------------------------------------------------")
      stop(" Incorrect input data.")
   }else{
      fine.dbh    = is.numeric  (dbh)    || all(is.na(dbh   ))
      fine.height = is.numeric  (height) || all(is.na(height))
      fine.wdens  = is.numeric  (wdens)  || all(is.na(wdens ))
      if (! all(c(fine.dbh,fine.height,fine.wdens))){
         cat0("-----------------------------------------------------------")
         cat0("   Not all variables have the correct type."                )
         cat0("   DBH    (numeric)   = ",fine.dbh                          )
         cat0("   HEIGHT (numeric)   = ",fine.height                       )
         cat0("   WDENS  (numeric)   = ",fine.wdens                        )
         cat0("-----------------------------------------------------------")
         stop(" Incorrect data types.")
      }#end if (! all(c(fine.dbh,fine.height,fine.wdens,fine.type,fine.dead)))
   }#end if ( length(lens) != 1)
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #     Find the allometry based on Chave et al. (2014).                                  #
   #---------------------------------------------------------------------------------------#
   ans = 0.0673 * (wdens*dbh^2*height)^0.976 / C2B
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end function size2agb.chave
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     We find the dbh based on above-ground biomass and wood density, following the        #
# allometry proposed by:                                                                   #
#                                                                                          #
# Baker, T. R., and co-authors, 2004: Variation in wood density determines spatial         #
#     patterns in Amazonian forest biomass. Global Change Biol., 10, 545-562.              #
#                                                                                          #
# Chave, J., B. Riera, M.A. Dubois, 2001: Estimation of biomass in a neotropical forest of #
#     French Guiana: spatial and temporal variability.  J. Trop. Ecol., 17, 79-96.         #
#------------------------------------------------------------------------------------------#
agb2dbh.baker <<- function(agb,wdens,allom="baker.chave"){


   ln.agb  = log(agb)
   ln.rhon = log(wdens / 0.58 / C2B)


   if (allom == "baker.chave"){
      #------ Use Chave's based function (equation 2, Baker et al., 2004). ----------------#
      dbh = exp( 1 / 2.42 * ( ln.agb - ln.rhon + 2.00 ) )
      #------------------------------------------------------------------------------------#
   }else if (allom == "baker.chambers"){
      #------ Use Chambers' function. -----------------------------------------------------#
      stop("Cannot invert Chambers' equation...")
      #------------------------------------------------------------------------------------#
   }else if (allom == "chave.2006"){
      #------ Use Chambers' function. -----------------------------------------------------#
      stop ("Cannot invert Chave 2006 equation...")
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#
   return(dbh)
}#end function agb2dbh.baker
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     Biomass allometry that is used by Sustainable Landscapes.  Results are always in     #
# kgC/plant.                                                                               #
#                                                                                          #
# References:                                                                              #
#                                                                                          #
# Chave, J., and co-authors, 2014: Improved allometric models to estimate tha aboveground  #
#     biomass of tropical trees.  Glob. Change Biol., 20, 3177-3190                        #
#     doi:10.1111/gcb.12629                                                                #
#                                                                                          #
# Goodman, R., and co-authors, 2013: Amazon palm biomass and allometry.  Forest Ecol.      #
#     Manag., 310, 994-1004. doi:10.1016/j.foreco.2013.09.045                              #
#                                                                                          #
# Palace, M., and co-authors, 2007: Necromass in undisturbed ad logged forests in the      #
#     Brazilian Amazon.  Forest Ecol. Manag., 238, 309-318.                                #
#     doi:10.1016/j.foreco.2006.10.026                                                     #
#                                                                                          #
# Schnitzer, S. A., and co-authors, 2006: Censusing and measuring lianas: a quantitative   #
#     comparison of the common methods.  Biotropica, 38, 581-591                           #
#     doi:10.1111/j.1744-7429.2006.00187.x                                                 #
#                                                                                          #
#                                                                                          #
#                                                                                          #
# Input:                                                                                   #
# ---------------------------------------------------------------------------------------- #
# dbh        --- Diameter at breast height [cm]                                            #
# height     --- Height [m]                                                                #
# wdens      --- Wood density [g/cm3]                                                      #
# type       --- Plant type:                                                               #
#                L - liana                                                                 #
#                P - palm                                                                  #
#                O - others (aka trees)                                                    #
#                In case type = NULL, all plants are assumed to be trees.                  #
# dead       --- Life status:                                                              #
#                TRUE  - plant is dead                                                     #
#                FALSE - plant is alive                                                    #
#                In case dead = NULL, all plants are assumed to be alive.                  #
# eps.dbh    --- Relative uncertainty for DBH [1 means 100%]                               #
# eps.height --- Relative uncertainty for height [1 means 100%]                            #
# eps.wdens  --- Relative uncertainty for wood density [1 means 100%]                      #
# out.err    --- Output error in addition to the estimates of biomass/necromass?           #
# ---------------------------------------------------------------------------------------- #
#
#
#
# ---------------------------------------------------------------------------------------- #
# Output:
# ---------------------------------------------------------------------------------------- #
# - In case out.err is FALSE, the function returns a vector with biomass for each entry.   #
# - In case out.err is TRUE, the output is a data frame with the following vectors         #
#   with the same length as the entries:                                                   #
#   * agb      -- biomass (necromass)                       [kgC]                          #
#   * ae.agb   -- uncertainty in biomass due to allometry   [kgC, not relative]            #
#   * me.agb   -- uncertainty in biomass due to measurement [kgC, not relative]            #
#   * lnagb    -- log(biomass), used for error propagation.                                #
#   * sd.lnagb -- standard error of log-biomass                                            #
#------------------------------------------------------------------------------------------#
agb.SL <<- function( dbh
                   , height
                   , wdens
                   , type       = NULL
                   , dead       = NULL
                   , eps.dbh    = 0.02
                   , eps.height = 0.167
                   , eps.wdens  = 0.10
                   , out.err    = FALSE
                   ){
   #---------------------------------------------------------------------------------------#
   #     "type" and "dead" may not be present, in which case we use dummy values.          #
   #---------------------------------------------------------------------------------------#
   if (is.null(type)) type = rep(  "O",times=length(dbh))
   if (is.null(dead)) dead = rep(FALSE,times=length(dbh))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Make sure all terms have the same length.                                         #
   #---------------------------------------------------------------------------------------#
   lens = unique(c(length(dbh),length(height),length(wdens),length(type),length(dead)))
   if ( length(lens) != 1 ){
      cat0("-----------------------------------------------------------")
      cat0("   Variables don't have the same length."                   )
      cat0("   DBH    = ",length(dbh)                                   )
      cat0("   HEIGHT = ",length(height)                                )
      cat0("   WDENS  = ",length(wdens)                                 )
      cat0("   TYPE   = ",length(type)                                  )
      cat0("   DEAD   = ",length(dead)                                  )
      cat0("-----------------------------------------------------------")
      stop(" Incorrect input data.")
   }else{
      fine.dbh    = is.numeric  (dbh)    || all(is.na(dbh   ))
      fine.height = is.numeric  (height) || all(is.na(height))
      fine.wdens  = is.numeric  (wdens)  || all(is.na(wdens ))
      fine.type   = is.character(type)   || all(is.na(type  ))
      fine.dead   = is.logical  (dead)   || all(is.na(dead  ))
      if (! all(c(fine.dbh,fine.height,fine.wdens,fine.type,fine.dead))){
         cat0("-----------------------------------------------------------")
         cat0("   Not all variables have the correct type."                )
         cat0("   DBH    (numeric)   = ",fine.dbh                          )
         cat0("   HEIGHT (numeric)   = ",fine.height                       )
         cat0("   WDENS  (numeric)   = ",fine.wdens                        )
         cat0("   TYPE   (character) = ",fine.type                         )
         cat0("   DEAD   (logical)   = ",fine.dead                         )
         cat0("-----------------------------------------------------------")
         stop(" Incorrect data types.")
      }#end if (! all(c(fine.dbh,fine.height,fine.wdens,fine.type,fine.dead)))
   }#end if ( length(lens) != 1)
   #---------------------------------------------------------------------------------------#

   #----- Initialise the output. ----------------------------------------------------------#
   agb = NA * dbh
   #---------------------------------------------------------------------------------------#
   
   #---------------------------------------------------------------------------------------#
   #     Find the possible statuses, then choose the best equation.                        #
   #---------------------------------------------------------------------------------------#
   tree  = type %in% "O" & (! dead)
   palm  = type %in% "P"
   liana = type %in% "L"
   dead  = type %in% "O" & dead
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #     AGB by type.                                                                      #
   #---------------------------------------------------------------------------------------#
   #----- Living tree: Chave et al. (2014). -----------------------------------------------#
   agb[tree ] = 0.0673 * (wdens[tree]*dbh[tree]^2*height[tree])^0.976 / C2B
   #----- Palm: Goodman et al. (2013). ----------------------------------------------------#
   agb[palm ] = exp(-3.448+0.588^2/2) * dbh[palm]^2.7483 / C2B
   #----- Liana: Schnitzer et al. (2006). -------------------------------------------------#
   agb[liana] = exp(-0.968) * dbh[liana]^2.657 / C2B
   #----- Dead trees: Palace et al. (2007). -----------------------------------------------#
   v1         = 0.091
   v0         = 0.01 / (1.3^-v1)
   a0         = 0.25 * pi * v0^2 / (1. - 2*v1)
   a1         = 1 - 2*v1
   agb[dead ] = 1000. * wdens[dead] * a0 * dbh[dead]^2 * height[dead]^a1 / C2B
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Check whether to estimate associated errors (measurement and allometry),         #
   # following:                                                                            #
   #                                                                                       #
   #   Chave, J., and co-authors, 2004: Error propagation and scaling for tropical forest  #
   #      biomass estimates. Phil. Trans. R. Soc. Lond. B., 359, 409-420.                  #
   #      doi:10.1098/rstb.2003.1425                                                       #
   #---------------------------------------------------------------------------------------#
   if (out.err){
      #------------------------------------------------------------------------------------#
      #       Find error associated with allometry.                                        #
      #------------------------------------------------------------------------------------#
      ae.agb        = NA * agb
      #----- Living tree: Chave et al. (2014). --------------------------------------------#
      ae.agb[tree ] = sqrt(exp(0.357^2)-1)*agb[tree ]
      #----- Palm: Goodman et al. (2013). -------------------------------------------------#
      ae.agb[palm ] = sqrt(exp(0.588^2)-1)*agb[palm ]
      #----- Liana: Schnitzer et al. (2006). ----------------------------------------------#
      ae.agb[liana] = sqrt(exp(1.016^2)-1)*agb[liana]
      #----- Standing dead: assumed the same as living trees. -----------------------------#
      ae.agb[dead ] = sqrt(exp(0.357^2)-1)*agb[dead ]
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #       Find error associated with measurements.                                     #
      #------------------------------------------------------------------------------------#
      me.agb  = NA * agb
      eps.dh  = sqrt( cov(x=dbh[tree|dead],y=height[tree|dead])
                    / ( mean(dbh[tree|dead])*mean(height[tree|dead]) )
                    )#end sqrt
      #----- Living tree: Chave et al. (2014). --------------------------------------------#
      me.agb[tree ] = agb[tree ] * sqrt( ( 2.0 * 0.976 * eps.dbh    )^2
                                       + (       0.976 * eps.height )^2
                                       + (       0.976 * eps.wdens  )^2
                                       + 2.0 * (2.0 * 0.976) * 0.976 * eps.dh * eps.dh
                                       )#end sqrt
      #----- Palm: Goodman et al. (2013). -------------------------------------------------#
      me.agb[palm ] = agb[palm ] * 2.7483 * eps.dbh
      #----- Liana: Schnitzer et al. (2006). ----------------------------------------------#
      me.agb[liana] = agb[liana] * 2.657  * eps.dbh
      #----- Standing dead: assumed the same as living trees. -----------------------------#
      me.agb[dead ] = agb[dead ] * sqrt( ( 2.0 * eps.dbh    )^2
                                       + ( a1  * eps.height )^2
                                       + (       eps.wdens  )^2
                                       + 2.0 * 2.0 * a1 * eps.dh
                                       )#end sqrt
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Save the standard error of the log scale: it will be useful for error analysis #
      #------------------------------------------------------------------------------------#
      sd.lnagb        = NA * agb
      #----- Living tree: Chave et al. (2014). --------------------------------------------#
      sd.lnagb[tree ] = 0.357
      #----- Palm: Goodman et al. (2013). -------------------------------------------------#
      sd.lnagb[palm ] = 0.588
      #----- Liana: Schnitzer et al. (2006). ----------------------------------------------#
      sd.lnagb[liana] = 1.016
      #----- Standing dead: unavailable, assume the same as living trees. -----------------#
      sd.lnagb[dead ] = 0.357
      #------------------------------------------------------------------------------------#

      #------------------------------------------------------------------------------------#
      #      Combine estimates and errors in a data frame.                                 #
      #------------------------------------------------------------------------------------#
      ans = data.frame( agb      = agb
                      , ae.agb   = ae.agb
                      , me.agb   = me.agb
                      , lnagb    = log(agb) - 0.5 * sd.lnagb^2
                      , sd.lnagb = sd.lnagb
                      )#end data.frame
      #------------------------------------------------------------------------------------#

   }else{
      #----- No error needed.  Return estimate only. --------------------------------------#
      ans = agb
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end function agb.SL
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     Volume allometry: this is literally the biomass equation above divided by wood       #
# density.                                                                                 #
#                                                                                          #
# Palace, M., and co-authors, 2007: Necromass in undisturbed ad logged forests in the      #
#     Brazilian Amazon.  Forest Ecol. Manag., 238, 309-318.                                #
#     doi:10.1016/j.foreco.2006.10.026                                                     #
#                                                                                          #
#------------------------------------------------------------------------------------------#
vol.SL <<- function(dbh,height,wdens,type=NULL,dead=NULL){

   vol = 0.002 * agb.SL(dbh,height,wdens,type=type,dead=dead) / wdens
   return(vol)
}#end function vol.SL
#==========================================================================================#
#==========================================================================================#
