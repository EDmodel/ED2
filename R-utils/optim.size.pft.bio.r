#==========================================================================================#
#==========================================================================================#
#      This function finds the relative error of both the size structure and the basal     #
# area by species.                                                                         #
#------------------------------------------------------------------------------------------#
optim.size.pft.bio <<- function( logdbh.minus.min            # Log of DBH above minimum
                               , obs.dbh.cnt                 # Frequency of each DBH class
                                                             #    (dbh.cnt.bks)
                               , dbh.cnt.bks                 # DBH breaks for size
                                                             #    evaluation
                               , dbh.table.obs               # DBH table by species
                               , dbh.table.bks = dbh.cnt.bks # DBH breaks by species
                               , dbh.max       = NULL        # Maximum DBH by species
                               , bsa.summ      = NULL        # Basal area summary
                               , dbh.bsa.bks   = dbh.cnt.bks #  DBH breaks for basal area
                               , cnt.log       = TRUE        # Use log for count?
                               , bsa.log       = FALSE       # Use log for basal area?
                               , survey.area                 # Total surveyed area (m2)
                               , census                      # Table with the full census
                               , datum                       # Table with the sought 
                                                             #    properties
                               , pft.var   = "pft"           # Variable to determine the PFT
                               , dbh.min   = 10.             # Minimum DBH
                               , tiny.off  = 1.e-7           # Tiny offset
                               ){


   #----- Revert back to DBH. -------------------------------------------------------------#
   dbh.now     = dbh.min + exp(logdbh.minus.min)
   pft.now     = rownames(dbh.table.obs)
   #---------------------------------------------------------------------------------------#



   #----- Find the error due to size distribution. ----------------------------------------#
   mod.dbh.cut = cut(dbh.now,dbh.cnt.bks)
   mod.dbh.cnt = table(mod.dbh.cut)
   obs.dbh.tot = sum(obs.dbh.cnt)
   mod.dbh.cnt = mod.dbh.cnt / obs.dbh.tot
   obs.dbh.cnt = obs.dbh.cnt / obs.dbh.tot
   if (cnt.log){
      ln.obs.dbh.cnt = log(tiny.off + obs.dbh.cnt)
      ln.mod.dbh.cnt = log(tiny.off + mod.dbh.cnt)
      mse.size       = ( sum( ( ln.obs.dbh.cnt - ln.mod.dbh.cnt )^2 )
                       / length(ln.obs.dbh.cnt) )
   }else{
      mse.size       = sum( ( mod.dbh.cnt - obs.dbh.cnt)^2 ) / length(obs.dbh.cnt)
   }#end if
   wgt.size = length(mse.size)
   #---------------------------------------------------------------------------------------#




   #----- Find the error due to size distribution. ----------------------------------------#
   if (is.list(dbh.table.bks)){
      #----- Get dbh thresholds for each species. -----------------------------------------#
      o              = order(names(dbh.table.bks))
      dbh.table.cut  = dbh.table.bks[o]
      #------------------------------------------------------------------------------------#


      #----- Get size count using species-specific thresholds. ----------------------------#
      dbh.list       = split(x=dbh.now,f=census[[pft.var]])
      dbh.table.tmp  = t( mapply( FUN  = function(dbh,dcut){
                                            dbh.cut  = cut(dbh,dcut)
                                            dbh.tbl  = table(dbh.cut)
                                            idx      = match(names(dbh.tbl),levels(dbh.cut))
                                            ans      = rep(0L,times=nlevels(dbh.cut))
                                            ans[idx] = dbh.tbl
                                            return(ans)
                                         }#end function
                                , dbh  = dbh.list
                                , dcut = dbh.table.cut
                                )#end mapply
                        )#end t
      #------------------------------------------------------------------------------------#
   }else{
      #---- Get size count using global thresholds. ---------------------------------------#
      dbh.table.cut       = cut(dbh.now,dbh.table.bks)
      dbh.table.tmp       = table(census[[pft.var]],dbh.table.cut)
      #------------------------------------------------------------------------------------#
   }#end if (is.null(dbh.table.bks))
   #---------------------------------------------------------------------------------------#



   #----- Align results so it matches the reference table. --------------------------------#
   dbh.table.mod       = array( data     = 0
                              , dim      = dim(dbh.table.obs)
                              , dimnames = dimnames(dbh.table.obs)
                              )#end array
   idx                 = match(rownames(dbh.table.tmp),rownames(dbh.table.mod))
   dbh.table.mod[idx,] = dbh.table.tmp
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Check for valid data, as sometimes PFT-specific distribution exists for a few     #
   # species only.  Then normalise data for both observations and estimates                #
   #---------------------------------------------------------------------------------------#
   dbh.table.use = ! is.na(dbh.table.obs)
   dbh.table.tot = sum(dbh.table.obs[dbh.table.use])
   dbh.table.obs = dbh.table.obs[dbh.table.use] / dbh.table.tot
   dbh.table.mod = dbh.table.mod[dbh.table.use] / dbh.table.tot
   #---------------------------------------------------------------------------------------#

   if (cnt.log){
      ln.dbh.table.mod = log(tiny.off + dbh.table.mod)
      ln.dbh.table.obs = log(tiny.off + dbh.table.obs)
      mse.class        = ( sum( ( ln.dbh.table.mod - ln.dbh.table.obs)^2 )
                         / length(dbh.table.obs) )
   }else{
      mse.class        = sum( (dbh.table.mod - dbh.table.obs)^2 ) / length(dbh.table.obs)
   }#end if
   wgt.class = length(mse.class)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Not every survey provides the basal area by species. In case it is missing, make a #
   # dummy error and dummy length.                                                         #
   #---------------------------------------------------------------------------------------#
   if ("bsa" %in% names(datum)){


      #----- Find modelled basal area by PFT. ---------------------------------------------#
      bsa.now            = pio4 * dbh.now^2 / survey.area
      bsa.pft.obs        = datum$bsa
      names(bsa.pft.obs) = datum[[pft.var]]
      bsa.pft.mod        = tapply(X=bsa.now,INDEX=census[[pft.var]],FUN=sum,na.rm=TRUE)
      idx                = match(names(bsa.pft.obs),names(bsa.pft.mod))
      bsa.pft.mod        = bsa.pft.mod[idx]
      #------------------------------------------------------------------------------------#


      #----- Normalise basal area, taking observed total as the reference. ----------------#
      bsa.pft.tot = sum(bsa.pft.obs)
      bsa.pft.obs = bsa.pft.obs / bsa.pft.tot
      bsa.pft.mod = bsa.pft.mod / bsa.pft.mod
      #------------------------------------------------------------------------------------#


      #----- Calculate error. -------------------------------------------------------------#
      if (bsa.log){
         ln.bsa.pft.obs = log(tiny.off + bsa.pft.obs)
         ln.bsa.pft.mod = log(tiny.off + bsa.pft.mod)
         mse.bsa = sum( ( ln.bsa.pft.mod - ln.bsa.pft.obs )^2 ) / length(ln.bsa.pft.obs)
      }else{
         mse.bsa = sum( ( bsa.pft.mod - bsa.pft.obs)^2 ) / length(bsa.pft.obs)
      }#end if (bsa.log)
      wgt.bsa  = length(mse.bsa)
      #------------------------------------------------------------------------------------#
   }else{
      #----- Make dummy error. ------------------------------------------------------------#
      mse.bsa = 0.
      wgt.bsa = 0.
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    In case a general basal area by size class.                                        #
   #---------------------------------------------------------------------------------------#
   if (is.null(bsa.summ)){
      #----- Make dummy error. ------------------------------------------------------------#
      mse.bss = 0.
      wgt.bss  = 0.
      #------------------------------------------------------------------------------------#
   }else{
      #----- Find the total basal area by DBH class. --------------------------------------#
      dbh.bsa.cut = cut(dbh.now,dbh.bsa.bks)
      bsa.now     = pio4 * dbh.now^2 / survey.area
      dbh.bsa.mod = tapply(X=bsa.now,INDEX=dbh.bsa.cut,FUN=sum,na.rm=TRUE)
      #------------------------------------------------------------------------------------#


      #----- Normalise basal area distribution, using observations as reference. ----------#
      dbh.bsa.tot = sum(bsa.summ)
      dbh.bsa.obs = bsa.summ    / dbh.bsa.tot
      dbh.bsa.mod = dbh.bsa.mod / dbh.bsa.tot
      #------------------------------------------------------------------------------------#



      #----- Find error associated with basal area size distribution. ---------------------#
      if (bsa.log){
         ln.dbh.bsa.obs = log(tiny.off + dbh.bsa.obs)
         ln.dbh.bsa.mod = log(tiny.off + dbh.bsa.mod)
         mse.bss        = ( sum( ( ln.dbh.bsa.mod - ln.dbh.bsa.obs )^2 )
                          / length(ln.dbh.bsa.obs) )
      }else{
         mse.bss        = sum( ( dbh.bsa.mod - dbh.bsa.obs )^2 ) / length(dbh.bsa.obs)
      }#end if (bsa.log)
      wgt.bss     = length(mse.bss)
      #------------------------------------------------------------------------------------#
   }#end if (is.null(bsa.summ))
   #---------------------------------------------------------------------------------------#

   #----- Find total error. ---------------------------------------------------------------#
   rmse = sqrt( ( wgt.size  * mse.size + wgt.class * mse.class
                + wgt.bsa   * mse.bsa  + wgt.bss   * mse.bss   )
              / ( wgt.size + wgt.class + wgt.bsa + wgt.bss ) )
   #---------------------------------------------------------------------------------------#

   return(rmse)
}#end function
#==========================================================================================#
#==========================================================================================#
