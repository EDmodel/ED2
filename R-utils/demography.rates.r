#==========================================================================================#
#==========================================================================================#
#      This function computes the recruitment rates of each group, the confidence          #
# interval, and the community-wide rate along with the confidence interval using the       #
# binomial distribution or bootstrap.  Recruitment rates are given in the interest rate    #
# format that accounts for varying time scales as in:                                      #
#                                                                                          #
# Nakagawa, M.; Tanaka, K.; Nakashizuka, T.; Ohkubo, T.; Kato, T.; Maeda, T.; Sato, K.;    #
#    Miguchi, H.; Nagamasu, H.; Ogino, K.; Teo, S.; Hamid, A. A.; Seng, L. H., 2000:       #
#    Impact of severe drought associated with the 1997-1998 El Nino in a tropical forest   #
#    in Sarawak. J. Trop. Ecol., 16, 355-367.                                              #
#                                                                                          #
#------------------------------------------------------------------------------------------#
recruitment.rate <<- function( property
                             , count
                             , global        = count
                             , p.use
                             , p.established
                             , taxon
                             , dtime
                             , R             = 1000
                             ){


   #---------------------------------------------------------------------------------------#
   #     Check whether the mandatory variables are given.                                  #
   #---------------------------------------------------------------------------------------#
   if (  missing(p.use) || missing(p.established) || missing(count) || missing(taxon) 
      || missing(dtime) || missing(property) ){
      cat("  At least one required variable is missing:  "     ,"\n")
      cat("  - Missing property:      ",missing(property)      ,"\n")
      cat("  - Missing count:         ",missing(count)         ,"\n")
      cat("  - Missing p.use:         ",missing(p.use)         ,"\n")
      cat("  - Missing p.established: ",missing(p.established) ,"\n")
      cat("  - Missing taxon:         ",missing(taxon)         ,"\n")
      cat("  - Missing dtime:         ",missing(dtime)         ,"\n")
      stop("   Required variables not provided...")
   }#end if
   #---------------------------------------------------------------------------------------#



   #---- Find the inverse of the time step (and check whether it is an scalar). -----------#
   if (length(dtime) == 1){
      dtime = rep(dtime,times=length(taxon))
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the total rates for community and taxon.                                     #
   #---------------------------------------------------------------------------------------#
   N.tx = tapply(X = count  * p.use                , INDEX = taxon, FUN = sum, na.rm = TRUE)
   E.tx = tapply(X = count  * p.use * p.established, INDEX = taxon, FUN = sum, na.rm = TRUE)
   N.gb = sum   (x = global * p.use                                          , na.rm = TRUE)
   E.gb = sum   (x = global * p.use * p.established                          , na.rm = TRUE)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Convert N.tx and E.tx to lists, and find the median.  In case none of the        #
   # trees were established or all of them were established, we add two trees with         #
   # median size and assume that one was recruited and the other was not, so bootstrap     #
   # can do something.                                                                     #
   #---------------------------------------------------------------------------------------#
   property.tx            = split (x = property     , f = taxon)
   p.use.tx               = split (x = p.use        , f = taxon)
   p.established.tx       = split (x = p.established, f = taxon)
   dtime.tx               = split (x = dtime        , f = taxon)
   median.tx              = sapply(X = property.tx, FUN = median, na.rm = TRUE)
   dtbar.tx               = sapply(X = dtime.tx   , FUN = mean  , na.rm = TRUE)
   zero.append            = N.tx == 0
   median.tx[zero.append] = 1.
   dont.append            = N.tx > 0 & ( N.tx != E.tx & E.tx != 0 )
   median.tx[dont.append] = NA
   dtbar.tx [dont.append] = NA
   property.tx            = lapply( X   = mapply( FUN = c
                                                , mapply(FUN=c,property.tx,median.tx)
                                                , median.tx
                                                )#end mapply
                                  , FUN = na.omit
                                  )#end lapply
   p.established.tx       = lapply( X   = mapply( FUN = c
                                                , mapply( FUN = c
                                                        , p.established.tx
                                                        , 1. + 0. * median.tx
                                                        )#end mapply
                                                , 0. * median.tx
                                                )#end mapply
                                  , FUN = na.omit
                                  )#end lapply
   p.use.tx               = lapply( X   = mapply( FUN = c
                                                , mapply( FUN = c
                                                        , p.use.tx
                                                        , 1. + 0. * median.tx
                                                        )#end mapply
                                                , 1. + 0. * median.tx
                                                )#end mapply
                                  , FUN = na.omit
                                  )#end lapply
   dtime.tx               = lapply( X   = mapply( FUN = c
                                                , mapply( FUN = c, dtime.tx, dtbar.tx)
                                                , dtbar.tx
                                                )#end mapply
                                  , FUN = na.omit
                                  )#end lapply
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Collapse the data into lists, and run bootstrap for each of the groups.          #
   #---------------------------------------------------------------------------------------#
   datum.tx    = lapply( X        = mapply( FUN           = list
                                          , property      = property.tx
                                          , p.established = p.established.tx
                                          , p.use         = p.use.tx
                                          , dtime         = dtime.tx
                                          , SIMPLIFY      = FALSE
                                          )#end mapply
                       , FUN      = data.frame
                       , MoreArgs = list(stringsAsFactors = FALSE)
                       )#end lapply
   boot.tx     = try(lapply(X= datum.tx,FUN=boot,statistic=boot.recruit,R=R))
   if ("try-error" %in% is(boot.tx)) browser()
   expected.tx = unlist(sapply(X=boot.tx,FUN=c)["t0",])
   q025.tx     = sapply(X= boot.tx ,FUN=boot.ci.lower,conf=0.95,type="perc")
   q975.tx     = sapply(X= boot.tx ,FUN=boot.ci.upper,conf=0.95,type="perc")
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Find the statistics for the global properties.  In case of extreme               #
   # probability (no recruits or all trees are recruits), add two trees, one that          #
   # recruited and another that did not, both with the median value of the property.       #
   #---------------------------------------------------------------------------------------#
   property.gb         = property
   p.established.gb    = p.established
   p.use.gb            = p.use
   dtime.gb            = dtime
   if ( N.gb > 0 && ( N.gb == E.gb || E.gb == 0 ) ){
      median.gb        = median(x = property, na.rm = TRUE)
      dtbar.gb         = mean  (x = dtime   , na.rm = TRUE)
      property.gb      = c(property.gb     ,median.gb,median.gb)
      p.established.gb = c(p.established.gb,       1.,       0.)
      p.use.gb         = c(p.use.gb        ,       1.,       1.)
      dtime.gb         = c(dtime.gb        , dtbar.gb, dtbar.gb)
   }else if (N.gb == 0){
      dtbar.gb         = mean  (x = dtime   , na.rm = TRUE)
      property.gb      = c(1.,1.)
      p.established.gb = c(1.,0.)
      p.use.gb         = c(1.,1.)
      dtime.gb         = c(dtbar.gb,dtbar.gb)
   }#end if
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #      Collapse the census to a data frame and run bootstrap.                           #
   #---------------------------------------------------------------------------------------#
   datum.gb    = data.frame   ( property         = property.gb
                              , p.established    = p.established.gb
                              , p.use            = p.use.gb
                              , dtime            = dtime.gb
                              , stringsAsFactors = FALSE
                              )#end data.frame
   boot.gb     = boot         (data=datum.gb,statistic=boot.recruit,R=R)
   expected.gb = boot.gb$t0
   q025.gb     = boot.ci.lower(boot.out=boot.gb,conf=0.95,type="perc")
   q975.gb     = boot.ci.upper(boot.out=boot.gb,conf=0.95,type="perc")
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #       Make a list with the answers.                                                   #
   #---------------------------------------------------------------------------------------#
   ans = list( taxon = data.frame( expected         = expected.tx
                                 , q025             = q025.tx
                                 , q975             = q975.tx
                                 , stringsAsFactors = FALSE
                                 )
             , comm  = data.frame( expected         = expected.gb
                                 , q025             = q025.gb
                                 , q975             = q975.gb
                                 , stringsAsFactors = FALSE
                                 )
             )#end list
   #---------------------------------------------------------------------------------------#
   return(ans)
}#end function recruitment.rate
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      This function computes the mortality rates of each group, the confidence interval,  #
# and the community-wide rate along with the confidence interval using the binomial        #
# distribution or bootstrap.  The mortality rates are given in the interest rate format,   #
# which accounts for varying time scales as in:                                            #
#                                                                                          #
# Sheil, D.; May, R. M., 1996: Motality and recruitment rate evaluations in heterogeneous  #
#     tropical forests.  J. Ecology, 84, 91-100.                                           #
#------------------------------------------------------------------------------------------#
mortality.rate <<- function( property
                           , count
                           , global     = count
                           , p.use
                           , p.survivor
                           , taxon
                           , dtime
                           , R          = 1000
                           ){


   #---------------------------------------------------------------------------------------#
   #     Check whether the mandatory variables are given.                                  #
   #---------------------------------------------------------------------------------------#
   if (  missing(p.use) || missing(p.survivor) || missing(count) || missing(taxon) 
      || missing(dtime) || missing(property) ){
      cat ("  At least one required variable is missing: ","\n")
      cat ("  - Missing property:   ",missing(property)   ,"\n")
      cat ("  - Missing count:      ",missing(count)      ,"\n")
      cat ("  - Missing p.use:      ",missing(p.use)      ,"\n")
      cat ("  - Missing p.survivor: ",missing(p.survivor) ,"\n")
      cat ("  - Missing taxon:      ",missing(taxon)      ,"\n")
      cat ("  - Missing dtime:      ",missing(dtime)      ,"\n")
      stop("   Required variables not provided...")
   }#end if
   #---------------------------------------------------------------------------------------#



   #---- Find the inverse of the time step (and check whether it is an scalar). -----------#
   if (length(dtime) == 1){
      dtime = rep(dtime,times=length(taxon))
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the total rates for community and taxon.                                     #
   #---------------------------------------------------------------------------------------#
   N.tx = tapply(X = count  * p.use             , INDEX = taxon, FUN = sum, na.rm = TRUE )
   S.tx = tapply(X = count  * p.use * p.survivor, INDEX = taxon, FUN = sum, na.rm = TRUE )
   N.gb = sum   (x = global * p.use                                       , na.rm = TRUE )
   S.gb = sum   (x = global * p.use * p.survivor                          , na.rm = TRUE )
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     For the mean mortality rates and 95% confidence interval, we must check whether   #
   # this distribution is at the edge (i.e., none of the trees died, or all trees died).   #
   # In case we are at the limit, we add two trees with median size, and assume that one   #
   # is still alive whereas the other died, so bootstrap can still be of some use.  The    #
   # extreme cases are likely to happen when the sample size is very small, so this should #
   # make the error bars large.                                                            #      
   #---------------------------------------------------------------------------------------#
   property.tx            = split (x = property  , f = taxon)
   p.survivor.tx          = split (x = p.survivor, f = taxon)
   p.use.tx               = split (x = p.use     , f = taxon)
   dtime.tx               = split (x = dtime        , f = taxon)
   median.tx              = sapply(X = property.tx, FUN = median, na.rm = TRUE)
   dtbar.tx               = sapply(X = dtime.tx   , FUN = median, na.rm = TRUE)
   zero.append            = N.tx == 0
   median.tx[zero.append] = 1.
   dont.append            = N.tx > 0 & ( N.tx != S.tx & S.tx != 0)
   median.tx[dont.append] = NA
   dtbar.tx [dont.append] = NA
   property.tx            = lapply( X   = mapply( FUN = c
                                                , mapply( FUN = c
                                                        , property.tx
                                                        , median.tx
                                                        , SIMPLIFY = FALSE
                                                        )#end mapply
                                                , median.tx
                                                , SIMPLIFY = FALSE
                                                )#end mapply
                                  , FUN = na.omit
                                  )#end lapply
   p.survivor.tx          = lapply( X   = mapply( FUN = c
                                                , mapply( FUN     = c
                                                        , p.survivor.tx
                                                        , 1. + 0. * median.tx
                                                        , SIMPLIFY = FALSE
                                                        )#end mapply
                                                , 0. * median.tx
                                                , SIMPLIFY = FALSE
                                                )#end mapply
                                  , FUN = na.omit
                                  )#end lapply
   p.use.tx               = lapply( X   = mapply( FUN = c
                                                , mapply( FUN     = c
                                                        , p.use.tx
                                                        , 1. + 0. * median.tx
                                                        , SIMPLIFY = FALSE
                                                        )#end mapply
                                                , 1. + 0. * median.tx
                                                , SIMPLIFY = FALSE
                                                )#end mapply
                                  , FUN = na.omit
                                  )#end lapply
   dtime.tx               = lapply( X   = mapply( FUN = c
                                                , mapply( FUN = c
                                                        , dtime.tx
                                                        , dtbar.tx
                                                        )#end mapply
                                                , dtbar.tx
                                                )#end mapply
                                  , FUN = na.omit
                                  )#end lapply
   #---------------------------------------------------------------------------------------#







   #---------------------------------------------------------------------------------------#
   #      Collapse the data into lists, and run bootstrap for each of the groups.          #
   #---------------------------------------------------------------------------------------#
   datum.tx = lapply( X        = mapply( FUN         = list
                                       , property    = property.tx
                                       , p.survivor  = p.survivor.tx
                                       , p.use       = p.use.tx
                                       , dtime       = dtime.tx
                                       , SIMPLIFY    = FALSE
                                       )#end mapply
                    , FUN      = data.frame
                    , MoreArgs = list(stringsAsFactors = FALSE)
                    )#end lapply
   boot.tx     = lapply(X= datum.tx,FUN=boot,statistic=boot.mortality,R=R)
   expected.tx = unlist(sapply(X=boot.tx,FUN=c)["t0",])
   q025.tx     = sapply(X= boot.tx ,FUN=boot.ci.lower,conf=0.95,type="perc")
   q975.tx     = sapply(X= boot.tx ,FUN=boot.ci.upper,conf=0.95,type="perc")
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Find the statistics for the global properties.  In case of extreme probability   #
   # (no survivors or all trees survived), add two trees, one that survived and another    #
   # that did not, both with the median value of the property.                             #
   #---------------------------------------------------------------------------------------#
   property.gb   = property
   p.survivor.gb = p.survivor
   p.use.gb      = p.use
   dtime.gb      = dtime
   if ( N.gb > 0 && ( N.gb == S.gb || S.gb == 0 ) ){
      median.gb      = median(x = property, na.rm = TRUE)
      dtbar.gb       = mean  (x = dtime   , na.rm = TRUE)
      property.gb    = c(property.gb  ,median.gb,median.gb)
      p.survivor.gb  = c(p.survivor.gb,       1.,       0.)
      p.use.gb       = c(p.use.gb     ,       1.,       1.)
      dtime.gb       = c(dtime.gb     , dtbar.gb, dtbar.gb)
   }else if (N.gb == 0){
      dtbar.gb       = mean  (x = dtime   , na.rm = TRUE)
      property.gb    = c(1.,1.)
      p.survivor.gb  = c(1.,0.)
      p.use.gb       = c(1.,1.)
      dtime.gb       = c(dtbar.gb,dtbar.gb)
   }#end if
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #      Collapse the census to a data frame and run bootstrap.                           #
   #---------------------------------------------------------------------------------------#
   datum.gb    = data.frame   ( property         = property.gb
                              , p.survivor       = p.survivor.gb
                              , p.use            = p.use.gb
                              , dtime            = dtime.gb
                              , stringsAsFactors = FALSE
                              )#end data.frame
   boot.gb     = boot         (data=datum.gb,statistic=boot.mortality,R=R)
   expected.gb = boot.gb$t0
   q025.gb     = boot.ci.lower(boot.out=boot.gb,conf=0.95,type="perc")
   q975.gb     = boot.ci.upper(boot.out=boot.gb,conf=0.95,type="perc")
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #       Make a list with the answers.                                                   #
   #---------------------------------------------------------------------------------------#
   ans = list( taxon = data.frame( expected         = expected.tx
                                 , q025             = q025.tx
                                 , q975             = q975.tx 
                                 , stringsAsFactors = FALSE
                                 )
             , comm  = data.frame( expected = expected.gb
                                 , q025     = q025.gb
                                 , q975     = q975.gb
                                 , stringsAsFactors = FALSE
                                 )
             )#end list
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end function mortality.rate
#==========================================================================================#
#==========================================================================================#







#==========================================================================================#
#==========================================================================================#
#      This function computes the expected growth rates and the confidence interval, using #
# bootstrap so we do not need to assume any distribution.                                  #
#------------------------------------------------------------------------------------------#
growth.rate <<- function(growth,count,taxon,R=1000){

   #------ Split the data according to the class. -----------------------------------------#
   growth.tx = split(x = growth, f = taxon)
   count.tx  = split(x = count , f = taxon)
   #---------------------------------------------------------------------------------------#



   #------ Pair the growth and the weights of each taxon into data frames. ----------------#
   datum.tx = lapply( X                = mapply( FUN      = list
                                               , growth   = growth.tx
                                               , count    = count.tx
                                               , SIMPLIFY = FALSE
                                               )#end mapply
                    , FUN              = data.frame
                    , stringsAsFactors = FALSE
                    )#end lapply
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Find the statistics.                                                              #
   #---------------------------------------------------------------------------------------#
   boot.tx     = lapply(X= datum.tx,FUN=boot,statistic=boot.growth,R=R)
   expected.tx = unlist(sapply(X=boot.tx,FUN=c)["t0",])
   q025.tx     = sapply(X= boot.tx ,FUN=boot.ci.lower,conf=0.95,type="perc")
   q975.tx     = sapply(X= boot.tx ,FUN=boot.ci.upper,conf=0.95,type="perc")
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Global data frame.                                                                #
   #---------------------------------------------------------------------------------------#
   datum.gb = data.frame( growth = growth, count = count , stringsAsFactors = FALSE)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the global rates.                                                            #
   #---------------------------------------------------------------------------------------#
   boot.gb     = boot(data=datum.gb,statistic=boot.growth,R=R)
   expected.gb = boot.gb$t0
   q025.gb     = boot.ci.lower(boot.out=boot.gb,conf=0.95,type="perc")
   q975.gb     = boot.ci.upper(boot.out=boot.gb,conf=0.95,type="perc")
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #       Make a list with the answers.                                                   #
   #---------------------------------------------------------------------------------------#
   ans = list( taxon = data.frame( expected         = expected.tx
                                 , q025             = q025.tx
                                 , q975             = q975.tx 
                                 , stringsAsFactors = FALSE
                                 )
             , comm  = data.frame( expected = expected.gb
                                 , q025     = q025.gb
                                 , q975     = q975.gb
                                 , stringsAsFactors = FALSE
                                 )
             )#end list
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end function growth.rate
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#      This function computes the accumulated recruitment rates of each group, the         #
# confidence interval, and the community-wide rate along with the confidence interval.     #
#------------------------------------------------------------------------------------------#
acc.recruitment.rate <<- function( property
                                 , count
                                 , global        = count
                                 , p.use
                                 , p.established
                                 , taxon
                                 , dtime
                                 , R             = 1000
                                 ){


   #---------------------------------------------------------------------------------------#
   #     Check whether the mandatory variables are given.                                  #
   #---------------------------------------------------------------------------------------#
   if (  missing(p.use) || missing(p.established) || missing(count) || missing(taxon) 
      || missing(dtime) || missing(property) ){
      cat ("  At least one required variable is missing:  "     ,"\n")
      cat ("  - Missing property:      ",missing(property)      ,"\n")
      cat ("  - Missing count:         ",missing(count)         ,"\n")
      cat ("  - Missing p.use:         ",missing(p.use)         ,"\n")
      cat ("  - Missing p.established: ",missing(p.established) ,"\n")
      cat ("  - Missing taxon:         ",missing(taxon)         ,"\n")
      cat ("  - Missing dtime:         ",missing(dtime)         ,"\n")
      stop("   Required variables not provided...")
   }#end if
   #---------------------------------------------------------------------------------------#



   #---- Find the inverse of the time step (and check whether it is an scalar). -----------#
   if (length(dtime) == 1){
      dtime = rep(dtime,times=length(taxon))
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the total rates for community and taxon.                                     #
   #---------------------------------------------------------------------------------------#
   N.tx = tapply(X = count  * p.use                , INDEX = taxon, FUN = sum, na.rm = TRUE)
   E.tx = tapply(X = count  * p.use * p.established, INDEX = taxon, FUN = sum, na.rm = TRUE)
   N.gb = sum   (x = global * p.use                                          , na.rm = TRUE)
   E.gb = sum   (x = global * p.use * p.established                          , na.rm = TRUE)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Convert N.tx and E.tx to lists, and find the median.  In case none of the        #
   # trees were established or all of them were established, we add two trees with         #
   # median size and assume that one was recruited and the other was not, so bootstrap     #
   # can do something.                                                                     #
   #---------------------------------------------------------------------------------------#
   property.tx            = split (x = count * property, f = taxon)
   p.use.tx               = split (x = p.use           , f = taxon)
   p.established.tx       = split (x = p.established   , f = taxon)
   dtime.tx               = split (x = dtime           , f = taxon)
   median.tx              = sapply(X = property.tx     , FUN = median, na.rm = TRUE)
   dtbar.tx               = sapply(X = dtime.tx        , FUN = median, na.rm = TRUE)
   zero.append            = N.tx == 0
   median.tx[zero.append] = 1.
   dont.append            = N.tx > 0 & ( N.tx != E.tx & E.tx != 0 )
   median.tx[dont.append] = NA
   dtbar.tx [dont.append] = NA
   property.tx            = lapply( X   = mapply( FUN = c
                                                , mapply(FUN=c,property.tx,median.tx)
                                                , median.tx
                                                )#end mapply
                                  , FUN = na.omit
                                  )#end lapply
   p.established.tx       = lapply( X   = mapply( FUN = c
                                                , mapply( FUN = c
                                                        , p.established.tx
                                                        , 1. + 0. * median.tx
                                                        )#end mapply
                                                , 0. * median.tx
                                                )#end mapply
                                  , FUN = na.omit
                                  )#end lapply
   p.use.tx               = lapply( X   = mapply( FUN = c
                                                , mapply( FUN = c
                                                        , p.use.tx
                                                        , 1. + 0. * median.tx
                                                        )#end mapply
                                                , 1. + 0. * median.tx
                                                )#end mapply
                                  , FUN = na.omit
                                  )#end lapply
   dtime.tx               = lapply( X   = mapply( FUN = c
                                                , mapply( FUN = c
                                                        , dtime.tx
                                                        , dtbar.tx
                                                        )#end mapply
                                                , dtbar.tx
                                                )#end mapply
                                  , FUN = na.omit
                                  )#end lapply
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #      Collapse the data into lists, and run bootstrap for each of the groups.          #
   #---------------------------------------------------------------------------------------#
   datum.tx    = lapply( X        = mapply( FUN           = list
                                          , property      = property.tx
                                          , p.established = p.established.tx
                                          , p.use         = p.use.tx
                                          , dtime         = dtime.tx
                                          , SIMPLIFY      = FALSE
                                          )#end mapply
                       , FUN      = data.frame
                       , MoreArgs = list(stringsAsFactors = FALSE)
                       )#end lapply
   boot.tx     = lapply(X= datum.tx,FUN=boot,statistic=boot.acc.recruit,R=R)
   expected.tx = unlist(sapply(X=boot.tx,FUN=c)["t0",])
   q025.tx     = sapply(X= boot.tx ,FUN=boot.ci.lower,conf=0.95,type="perc")
   q975.tx     = sapply(X= boot.tx ,FUN=boot.ci.upper,conf=0.95,type="perc")
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Find the statistics for the global properties.  In case of extreme               #
   # probability (no recruits or all trees are recruits), add two trees, one that          #
   # recruited and another that did not, both with the median value of the property.       #
   #---------------------------------------------------------------------------------------#
   property.gb         = global * property
   p.established.gb    = p.established
   p.use.gb            = p.use
   dtime.gb            = dtime
   if ( N.gb > 0 && ( N.gb == E.gb || E.gb == 0 ) ){
      median.gb        = median(x = property.gb, na.rm = TRUE)
      dtbar.gb         = mean  (x = dtime      , na.rm = TRUE)
      property.gb      = c(property.gb     ,median.gb,median.gb)
      p.established.gb = c(p.established.gb,       1.,       0.)
      p.use.gb         = c(p.use.gb        ,       1.,       1.)
      dtime.gb         = c(dtime.gb        , dtbar.gb, dtbar.gb)
   }else if (N.gb == 0){
      dtbar.gb         = mean  (x = dtime   , na.rm = TRUE)
      property.gb      = c(1.,1.)
      p.established.gb = c(1.,0.)
      p.use.gb         = c(1.,1.)
      dtime.gb         = c(dtbar.gb,dtbar.gb)
   }#end if
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #      Collapse the census to a data frame and run bootstrap.                           #
   #---------------------------------------------------------------------------------------#
   datum.gb    = data.frame   ( property         = property.gb
                              , p.established    = p.established.gb
                              , p.use            = p.use.gb
                              , dtime            = dtime.gb
                              , stringsAsFactors = FALSE
                              )#end data.frame
   boot.gb     = boot         (data=datum.gb,statistic=boot.acc.recruit,R=R)
   expected.gb = boot.gb$t0
   q025.gb     = boot.ci.lower(boot.out=boot.gb,conf=0.95,type="perc")
   q975.gb     = boot.ci.upper(boot.out=boot.gb,conf=0.95,type="perc")
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #       Make a list with the answers.                                                   #
   #---------------------------------------------------------------------------------------#
   ans = list( taxon = data.frame( expected         = expected.tx
                                 , q025             = q025.tx
                                 , q975             = q975.tx 
                                 , stringsAsFactors = FALSE
                                 )
             , comm  = data.frame( expected = expected.gb
                                 , q025     = q025.gb
                                 , q975     = q975.gb
                                 , stringsAsFactors = FALSE
                                 )
             )#end list
   #---------------------------------------------------------------------------------------#
   return(ans)
}#end function acc.recruitment.rate
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      This function computes the accumulated mortality rates of each group, the confi-    #
# dence interval, and the community-wide rate along with the confidence interval.          #
#------------------------------------------------------------------------------------------#
acc.mortality.rate <<- function( property
                               , count
                               , global     = count
                               , p.use
                               , p.survivor
                               , taxon
                               , dtime
                               , R          = 1000
                               ){


   #---------------------------------------------------------------------------------------#
   #     Check whether the mandatory variables are given.                                  #
   #---------------------------------------------------------------------------------------#
   if (  missing(p.use) || missing(p.survivor) || missing(count) || missing(taxon) 
      || missing(dtime) || missing(property) ){
      cat ("  At least one required variable is missing: ","\n")
      cat ("  - Missing property:   ",missing(property)   ,"\n")
      cat ("  - Missing count:      ",missing(count)      ,"\n")
      cat ("  - Missing p.use:      ",missing(p.use)      ,"\n")
      cat ("  - Missing p.survivor: ",missing(p.survivor) ,"\n")
      cat ("  - Missing taxon:      ",missing(taxon)      ,"\n")
      cat ("  - Missing dtime:      ",missing(dtime)      ,"\n")
      stop("   Required variables not provided...")
   }#end if
   #---------------------------------------------------------------------------------------#



   #---- Find the inverse of the time step (and check whether it is an scalar). -----------#
   if (length(dtime) == 1){
      dtime = rep(dtime,times=length(taxon))
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the total rates for community and taxon.                                     #
   #---------------------------------------------------------------------------------------#
   N.tx = tapply(X = count  * p.use             , INDEX = taxon, FUN = sum, na.rm = TRUE )
   S.tx = tapply(X = count  * p.use * p.survivor, INDEX = taxon, FUN = sum, na.rm = TRUE )
   N.gb = sum   (x = global * p.use                                       , na.rm = TRUE )
   S.gb = sum   (x = global * p.use * p.survivor                          , na.rm = TRUE )
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     For the mean mortality rates and 95% confidence interval, we must check whether   #
   # this distribution is at the edge (i.e., none of the trees died, or all trees died).   #
   # In case we are at the limit, we add two trees with median size, and assume that one   #
   # is still alive whereas the other died, so bootstrap can still be of some use.  The    #
   # extreme cases are likely to happen when the sample size is very small, so this should #
   # make the error bars large.                                                            #      
   #---------------------------------------------------------------------------------------#
   property.tx            = split (x = count*property, f = taxon)
   p.survivor.tx          = split (x = p.survivor    , f = taxon)
   p.use.tx               = split (x = p.use         , f = taxon)
   dtime.tx               = split (x = dtime         , f = taxon)
   median.tx              = sapply(X = property.tx   , FUN = median, na.rm = TRUE)
   dtbar.tx               = sapply(X = dtime.tx   , FUN = median, na.rm = TRUE)
   zero.append            = N.tx == 0
   median.tx[zero.append] = 1.
   dont.append            = N.tx > 0 & ( N.tx != S.tx & S.tx != 0)
   median.tx[dont.append] = NA
   dtbar.tx [dont.append] = NA
   property.tx            = lapply( X   = mapply( FUN = c
                                                , mapply( FUN = c
                                                        , property.tx
                                                        , median.tx
                                                        , SIMPLIFY = FALSE
                                                        )#end mapply
                                                , median.tx
                                                , SIMPLIFY = FALSE
                                                )#end mapply
                                  , FUN = na.omit
                                  )#end lapply
   p.survivor.tx          = lapply( X   = mapply( FUN = c
                                                , mapply( FUN     = c
                                                        , p.survivor.tx
                                                        , 1. + 0. * median.tx
                                                        , SIMPLIFY = FALSE
                                                        )#end mapply
                                                , 0. * median.tx
                                                , SIMPLIFY = FALSE
                                                )#end mapply
                                  , FUN = na.omit
                                  )#end lapply
   p.use.tx               = lapply( X   = mapply( FUN = c
                                                , mapply( FUN     = c
                                                        , p.use.tx
                                                        , 1. + 0. * median.tx
                                                        , SIMPLIFY = FALSE
                                                        )#end mapply
                                                , 1. + 0. * median.tx
                                                , SIMPLIFY = FALSE
                                                )#end mapply
                                  , FUN = na.omit
                                  )#end lapply
   dtime.tx               = lapply( X   = mapply( FUN = c
                                                , mapply( FUN = c
                                                        , dtime.tx
                                                        , dtbar.tx
                                                        )#end mapply
                                                , dtbar.tx
                                                )#end mapply
                                  , FUN = na.omit
                                  )#end lapply
   #---------------------------------------------------------------------------------------#







   #---------------------------------------------------------------------------------------#
   #      Collapse the data into lists, and run bootstrap for each of the groups.          #
   #---------------------------------------------------------------------------------------#
   datum.tx = lapply( X        = mapply( FUN         = list
                                       , property    = property.tx
                                       , p.survivor  = p.survivor.tx
                                       , p.use       = p.use.tx
                                       , dtime       = dtime.tx
                                       , SIMPLIFY    = FALSE
                                       )#end mapply
                    , FUN      = data.frame
                    , MoreArgs = list(stringsAsFactors = FALSE)
                    )#end lapply
   boot.tx     = lapply(X= datum.tx,FUN=boot,statistic=boot.acc.mortality,R=R)
   expected.tx = unlist(sapply(X=boot.tx,FUN=c)["t0",])
   q025.tx     = sapply(X= boot.tx ,FUN=boot.ci.lower,conf=0.95,type="perc")
   q975.tx     = sapply(X= boot.tx ,FUN=boot.ci.upper,conf=0.95,type="perc")
   if (any(expected.tx %>% q975.tx | expected.tx %<% q025.tx,na.rm=TRUE)) browser()
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Find the statistics for the global properties.  In case of extreme probability   #
   # (no survivors or all trees survived), add two trees, one that survived and another    #
   # that did not, both with the median value of the property.                             #
   #---------------------------------------------------------------------------------------#
   property.gb   = global * property
   p.survivor.gb = p.survivor
   p.use.gb      = p.use
   dtime.gb      = dtime
   if ( N.gb > 0 && ( N.gb == S.gb || S.gb == 0 ) ){
      median.gb      = median(x = property * global, na.rm = TRUE)
      dtbar.gb       = mean  (x = dtime   , na.rm = TRUE)
      property.gb    = c(property.gb  ,median.gb,median.gb)
      p.survivor.gb  = c(p.survivor.gb,       1.,       0.)
      p.use.gb       = c(p.use.gb     ,       1.,       1.)
      dtime.gb       = c(dtime.gb     , dtbar.gb, dtbar.gb)
   }else if (N.gb == 0){
      dtbar.gb       = mean  (x = dtime   , na.rm = TRUE)
      property.gb    = c(1.,1.)
      p.survivor.gb  = c(1.,0.)
      p.use.gb       = c(1.,1.)
      dtime.gb       = c(dtbar.gb,dtbar.gb)
   }#end if
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #      Collapse the census to a data frame and run bootstrap.                           #
   #---------------------------------------------------------------------------------------#
   datum.gb    = data.frame   ( property         = property.gb
                              , p.survivor       = p.survivor.gb
                              , p.use            = p.use.gb
                              , dtime            = dtime.gb
                              , stringsAsFactors = FALSE
                              )#end data.frame
   boot.gb     = boot         (data=datum.gb,statistic=boot.acc.mortality,R=R)
   expected.gb = boot.gb$t0
   q025.gb     = boot.ci.lower(boot.out=boot.gb,conf=0.95,type="perc")
   q975.gb     = boot.ci.upper(boot.out=boot.gb,conf=0.95,type="perc")
   if (any(expected.gb %>% q975.gb | expected.gb %<% q025.gb,na.rm=TRUE)) browser()
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #       Make a list with the answers.                                                   #
   #---------------------------------------------------------------------------------------#
   ans = list( taxon = data.frame( expected         = expected.tx
                                 , q025             = q025.tx
                                 , q975             = q975.tx 
                                 , stringsAsFactors = FALSE
                                 )
             , comm  = data.frame( expected = expected.gb
                                 , q025     = q025.gb
                                 , q975     = q975.gb
                                 , stringsAsFactors = FALSE
                                 )
             )#end list
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end function acc.mortality.rate
#==========================================================================================#
#==========================================================================================#







#==========================================================================================#
#==========================================================================================#
#      This function computes the expected accumulated growth rate and the confidence      #
# interval, using bootstrap so we do not need to assume any distribution.                  #
#------------------------------------------------------------------------------------------#
acc.growth.rate <<- function(nok,lok,pop,gpop=pop,dtime,taxon,R=1000){

   #------ Split the data according to the class. -----------------------------------------#
   nok.tx    = split(x = nok   , f = taxon)
   lok.tx    = split(x = lok   , f = taxon)
   pop.tx    = split(x = pop   , f = taxon)
   dtime.tx  = split(x = dtime , f = taxon)
   #---------------------------------------------------------------------------------------#



   #------ Pair the growth and the weights of each taxon into data frames. ----------------#
   datum.tx = lapply( X                = mapply( FUN      = list
                                               , nok      = nok.tx
                                               , lok      = lok.tx
                                               , pop      = pop.tx
                                               , dtime    = dtime.tx
                                               , SIMPLIFY = FALSE
                                               )#end mapply
                    , FUN              = data.frame
                    , stringsAsFactors = FALSE
                    )#end lapply
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Find the statistics.                                                              #
   #---------------------------------------------------------------------------------------#
   boot.tx     = lapply(X= datum.tx,FUN=boot,statistic=boot.acc.growth,R=R)
   expected.tx = unlist(sapply(X=boot.tx,FUN=c)["t0",])
   q025.tx     = sapply(X= boot.tx ,FUN=boot.ci.lower,conf=0.95,type="perc")
   q975.tx     = sapply(X= boot.tx ,FUN=boot.ci.upper,conf=0.95,type="perc")
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Global data frame.                                                                #
   #---------------------------------------------------------------------------------------#
   datum.gb = data.frame( nok              = nok
                        , lok              = lok
                        , pop              = gpop
                        , dtime            = dtime
                        , stringsAsFactors = FALSE)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the global rates.                                                            #
   #---------------------------------------------------------------------------------------#
   boot.gb     = boot(data=datum.gb,statistic=boot.acc.growth,R=R)
   expected.gb = boot.gb$t0
   q025.gb     = boot.ci.lower(boot.out=boot.gb,conf=0.95,type="perc")
   q975.gb     = boot.ci.upper(boot.out=boot.gb,conf=0.95,type="perc")
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #       Make a list with the answers.                                                   #
   #---------------------------------------------------------------------------------------#
   ans = list( taxon = data.frame( expected         = expected.tx
                                 , q025             = q025.tx
                                 , q975             = q975.tx 
                                 , stringsAsFactors = FALSE
                                 )
             , comm  = data.frame( expected = expected.gb
                                 , q025     = q025.gb
                                 , q975     = q975.gb
                                 , stringsAsFactors = FALSE
                                 )
             )#end list
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end function acc.growth.rate
#==========================================================================================#
#==========================================================================================#
