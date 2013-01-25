#==========================================================================================#
#==========================================================================================#
#      This function finds the relative error of both the size structure and the basal     #
# area by species.                                                                         #
#------------------------------------------------------------------------------------------#
optim.size.pft.bio <<- function( logdbh.minus.min   # Log of DBH above minimum
                               , obs.dbh.cnt        # Frequency of each DBH class (dbh.bks)
                               , dbh.bks            # DBH breaks for size evaluation
                               , census             # Table with the full census
                               , datum              # Table with the sought properties
                               , pft.var   = "pft"  # Variable to determine the PFT
                               , dbh.min   = 10.    # Minimum DBH
                               , tiny.off  = 1.e-20 # Tiny offset
                               ){

   #----- Revert back to DBH. -------------------------------------------------------------#
   dbh.now     = dbh.min + exp(logdbh.minus.min)
   #---------------------------------------------------------------------------------------#



   #----- Find the error due to size distribution. ----------------------------------------#
   dbh.cut     = cut(dbh.now,dbh.bks)
   dbh.cnt     = tiny.off + table(dbh.cut)
   if (tiny.off == 0){
      rmse.size   = ( sum((dbh.cnt/sum(dbh.cnt)-obs.dbh.cnt/sum(obs.dbh.cnt))^2) 
                    / length(obs.dbh.cnt) )
   }else{
      rmse.size   = sum((log(dbh.cnt)-log(obs.dbh.cnt))^2) / length(obs.dbh.cnt)
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Find the error due to basal area by PFT. ----------------------------------------#
   bsa.now     = 0.25 * pi * dbh.now^2
   bsa.pft     = tapply(X=bsa.now,INDEX=census[[pft.var]],FUN=sum,na.rm=TRUE)
   idx         = match(datum[[pft.var]],names(bsa.pft))
   rmse.bsa    = sum((log(datum$bsa) - log(bsa.pft[idx]))^2) / length(datum$bsa)
   #---------------------------------------------------------------------------------------#

   rmse        = sqrt(rmse.size^2+rmse.bsa^2)

   return(rmse)
}#end function
#==========================================================================================#
#==========================================================================================#
