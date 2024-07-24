#==========================================================================================#
#==========================================================================================#
#     This function assigns a PFT given the genus and the full data set.  This will        #
# use imputed data for those genera that do not have all traits.                           #
#------------------------------------------------------------------------------------------#
assign.pft <<- function(datum,path=srcdir,approach=c("rho","cluster"),verbose=FALSE,...){
   #----- Standardise the approach for defining the PFT. ----------------------------------#
   approach = match.arg(approach)
   #---------------------------------------------------------------------------------------#



   if (approach %in% "rho"){
      #---- Find variable that represents wood density. -----------------------------------#
      idx.rho = which(names(datum) %in% c("rho","wood.dens","WSD"))
      if (length(idx.rho) == 1){
         #---- Make sure that the column is named wood.dens. ------------------------------#
         nm.datum          = names(datum)
         nm.datum[idx.rho] = "wood.dens"
         names(datum)      = nm.datum
         #---------------------------------------------------------------------------------#
      }else if(length(idx.rho) == 0){
         datum$wood.dens = rep(NA_real_,times=nrow(datum))
      }else{
         cat0("----------------------------------------------")
         cat0(" The following columns have names that look like wood density.")
         for (idx in idx.rho) cat0("   - ",idx," (",names(datum)[idx],")")
         cat0("----------------------------------------------")
         stop(" Ambiguous data frame, it can't have more than one wood density column.")
      }#end if
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Make sure wood density is entirely defined.                                    #
      #------------------------------------------------------------------------------------#
      fill = is.na(datum$wood.dens)
      if (any(fill)){
         warning(" Wood density has NA values, filling the data.")
         dfill                 = datum[fill,,drop=FALSE]
         dfill                 = find.trait( datum   = dfill
                                           , trait   = "wood.dens"
                                           , verbose = verbose
                                           ,...
                                           )#end find.trait
         datum$wood.dens[fill] = dfill$wood.dens
      }#end if (any(fill))
      #------------------------------------------------------------------------------------#


      #------ Assign PFT values based on most similar wood density. -----------------------#
      ans = mypfts[as.integer(cut(datum$wood.dens,pft.breaks))]
      #------------------------------------------------------------------------------------#
   }else{


      #----- Load the cluster analysis with the medoid matrix. ----------------------------#
      dummy = load(file.path(srcdir,"cluster_traits_species.RData"))
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Initialise genus with the actual one, but replace those not find in the        #
      # look-up table.                                                                     #
      #------------------------------------------------------------------------------------#
      if ("scientific.name" %in% names(datum)) datum$scientific = datum$scientific.name
      if ("genus.name"      %in% names(datum)) datum$genus      = datum$genus.name
      if ("family.name"     %in% names(datum)) datum$family     = datum$family.name
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #    Initialise use.scientific with the actual scientific name.  We will replace     #
      # those that are not found in the look-up table for traits.                          #
      #------------------------------------------------------------------------------------#
      use.scientific = datum$scientific
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #    Retrieve scientific name.  In case the scientific name was not used in the      #
      # cluster analysis, randomly select one name from the same genus, or from the same   #
      # family in case the genus is not available.  In case nothing is known, make a       #
      # completely random selection.                                                       #
      #------------------------------------------------------------------------------------#
      yes.sci         = datum$scientific %in% taxon.summ$scientific
      yes.gen         = datum$genus      %in% taxon.summ$genus
      yes.fam         = datum$family     %in% taxon.summ$family
      nsciygen        = (! yes.sci) & yes.gen
      ngenyfam        = (! yes.gen) & yes.fam
      nothing         = (! yes.fam)
      oddity          = ( (yes.sci & (! (yes.gen & yes.fam)) )
                        | (yes.gen & (! yes.fam) )
                        )#end oddity
      if (any(oddity)) browser()
      sci.genus.mean  = datum[nsciygen,c("scientific","genus","family")]
      sci.family.mean = datum[ngenyfam,c("scientific","genus","family")]
      sci.loose.mean  = datum[nothing ,c("scientific","genus","family")]
      genus.keep      = ! duplicated(sci.genus.mean$scientific )
      family.keep     = ! duplicated(sci.family.mean$scientific)
      loose.keep      = ! duplicated(sci.loose.mean$scientific )
      sci.genus.mean  = sci.genus.mean [genus.keep ,,drop=FALSE]
      sci.family.mean = sci.family.mean[family.keep,,drop=FALSE]
      sci.loose.mean  = sci.loose.mean [loose.keep ,,drop=FALSE]
      oge             = order(sci.genus.mean$scientific )
      ofm             = order(sci.family.mean$scientific)
      ols             = order(sci.loose.mean$scientific )
      sci.genus.mean  = sci.genus.mean [oge,,drop=FALSE]
      sci.family.mean = sci.family.mean[ofm,,drop=FALSE]
      sci.loose.mean  = sci.loose.mean [ols,,drop=FALSE]
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    Loop through scientific names, and fill in with sampling from the same genus.   #
      #------------------------------------------------------------------------------------#
      for (s in seq_along(sci.genus.mean$scientific)){
         #----- Handy aliases for current genus. ------------------------------------------#
         s.scientific = sci.genus.mean$scientific[s]
         s.genus      = sci.genus.mean$genus     [s]
         s.family     = sci.genus.mean$family    [s]
         s.sel        = datum$scientific %in% s.scientific
         ns.sel       = sum(s.sel)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      List of scientific names from the same genus that are available in the     #
         # look up table.  We also include scientific names that exist in the source data  #
         # set, belong to the same genus, and are found in the look-up table, to increase  #
         # chances of selecting a scientific name observed in the source data itself.      #
         #---------------------------------------------------------------------------------#
         d.sel    = yes.sci & (datum$genus %in% s.genus)
         t.sel    = taxon.summ$genus %in% s.genus
         sci.pool = c(datum$scientific[d.sel],taxon.summ$scientific[t.sel])
         #---------------------------------------------------------------------------------#


         #----- Overwrite the names with those that are found in the look-up table. -------#
         use.scientific[s.sel] = lit.sample(x=sci.pool,size=ns.sel,replace=TRUE)
         #---------------------------------------------------------------------------------#
      }#end or (s in seq_along(sci.genus.mean$scientific))
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    Loop through scientific names, and fill in with sampling from the same family.  #
      #------------------------------------------------------------------------------------#
      for (s in seq_along(sci.family.mean$scientific)){
         #----- Handy aliases for current scientific name. --------------------------------#
         s.scientific = sci.family.mean$scientific[s]
         s.genus      = sci.family.mean$genus     [s]
         s.family     = sci.family.mean$family    [s]
         s.sel        = datum$scientific %in% s.scientific
         ns.sel       = sum(s.sel)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      List of scientific names from the same family that are available in the    #
         # look up table.  We also include genera that exist in the source data set,       #
         # belong to the same family, and are found in the look-up table, to increase      #
         # chances of selecting a genus observed in the source data itself.                #
         #---------------------------------------------------------------------------------#
         d.sel    = yes.sci & (datum$family %in% s.family)
         t.sel    = taxon.summ$family %in% s.family
         sci.pool = c(datum$scientific[d.sel],taxon.summ$scientific[t.sel])
         #---------------------------------------------------------------------------------#


         #----- Overwrite the species with those that are found in the look-up table. -----#
         use.scientific[s.sel] = lit.sample(x=sci.pool,size=ns.sel,replace=TRUE)
         #---------------------------------------------------------------------------------#
      }#end for (s in seq_along(sci.family.mean$scientific))
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #    Loop through genera, and fill in with sampling from any family.                 #
      #------------------------------------------------------------------------------------#
      for (s in seq_along(sci.loose.mean$scientific)){
         #----- Handy aliases for current scientific name. --------------------------------#
         s.scientific = sci.loose.mean$scientific[s]
         s.genus      = sci.loose.mean$genus     [s]
         s.family     = sci.loose.mean$family    [s]
         s.sel        = datum$scientific %in% s.scientific
         ns.sel       = sum(s.sel)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      List of scientific names that are available in the look up table.  We also #
         # include species that exist in the source data set are found in the look-up      #
         # table, to increase chances of selecting a species observed in the source data   #
         # itself.                                                                         #
         #---------------------------------------------------------------------------------#
         d.sel    = yes.sci
         sci.pool = c(datum$scientific[d.sel],taxon.summ$scientific)
         #---------------------------------------------------------------------------------#


         #----- Overwrite the species with those that are found in the look-up table. -----#
         use.scientific[s.sel] = lit.sample(x=sci.pool,size=ns.sel,replace=TRUE)
         #---------------------------------------------------------------------------------#
      }#end for (s in seq_along(sci.loose.mean$scientific))
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     List all species that were not filled with species information.                #
      #------------------------------------------------------------------------------------#
      if (verbose && (nrow(sci.genus.mean) > 0)){
         cat0("")
         cat0("---------------------------------------------------------------------------")
         cat0(" Found species that were not found in trait data base (check for synonyms)!")
         print (sci.genus.mean,quote=FALSE)
         cat0("---------------------------------------------------------------------------")
         cat0("")
      }#end if (verbose && nrow(sci.genus.mean[genus.only,,drop=FALSE]) > 0)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     List all species that were not filled with species information.                #
      #------------------------------------------------------------------------------------#
      if (verbose && (nrow(sci.family.mean) > 0)){
         cat0("")
         cat0("--------------------------------------------------------------------------")
         cat0(" Found genera that were not found in trait data base (check for synonyms)!")
         print (sci.family.mean,quote=FALSE)
         cat0("--------------------------------------------------------------------------")
         cat0("")
      }#end if (verbose && nrow(sci.genus.mean[genus.only,,drop=FALSE]) > 0)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     List all genera that didn't belong to any known family.                        #
      #------------------------------------------------------------------------------------#
      if (verbose && (nrow(sci.loose.mean) > 0)){
         cat0("")
         cat0("-----------------------------------------------------------------------")
         cat0(" Found species from families with no valid data in the trait data base!")
         print(sci.loose.mean,quote=FALSE)
         cat0("-----------------------------------------------------------------------")
         cat0("")
      }#end if
      #------------------------------------------------------------------------------------#


      #----- List of genera to obtain data. -----------------------------------------------#
      un.scientific  = sort(unique(use.scientific))
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Create data frame with traits that are going to be used for the cluster       #
      # assessment.                                                                        #
      #------------------------------------------------------------------------------------#
      tcluster           = data.frame(matrix(nrow=nrow(datum),ncol=0))
      for (tt in which(trait$cluster)){
         #----- Handy aliases. ------------------------------------------------------------#
         t.vnam    = trait$vnam   [tt]
         t.desc    = trait$desc   [tt]
         t.type    = trait$type   [tt]
         t.units   = trait$units  [tt]
         t.impute  = trait$impute [tt]
         t.cluster = trait$cluster[tt]
         t.sma     = trait$sma    [tt]
         t.vlog    = trait$vlog   [tt]
         #---------------------------------------------------------------------------------#


         #----- Initialise data set. ------------------------------------------------------#
         tcluster[[t.vnam]] = as(rep(NA,times=nrow(datum)),t.type)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Loop through all genera, and sample from the imputation matrix.             #
         #---------------------------------------------------------------------------------#
         for (s in seq_along(un.scientific)){
            #---- Select lines for sampling and replacement. ------------------------------#
            sel          = use.scientific %in% un.scientific[s]
            nsel         = sum(sel)
            r.scientific = which(rownames(mimpute[[t.vnam]]) %in% un.scientific[s])
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Obtain sample.  This ensures that imputed data will not be a unique      #
            # value, whereas well-constrained genera will have no variability (so all      #
            # plants in the same genus will have the same number).                         #
            #------------------------------------------------------------------------------#
            tcluster[[t.vnam]][sel] = lit.sample( x       = mimpute[[t.vnam]][r.scientific,]
                                                , size    = nsel
                                                , replace = TRUE
                                                )#end lit.sample
            #------------------------------------------------------------------------------#
         }#end for (s in seq_along(un.scientific))
         #---------------------------------------------------------------------------------#
      }#end for (tt in which(trait$cluster))
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      For each row, we compute the dissimilarity matrix with the medoid, then pick  #
      # the least dissimilar one.                                                          #
      #------------------------------------------------------------------------------------#
      ans        = which.pft(datum=tcluster,medoid=medoid,trait=trait,wcluster=wcluster)
      names(ans) = NULL
      return(ans)
      #------------------------------------------------------------------------------------#
   }#end if (approach %in% "rho")
   #---------------------------------------------------------------------------------------#
}#end function assign.pft
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#    This function finds the most similar group to the data input.
#------------------------------------------------------------------------------------------#
which.pft <<- function(datum,medoid,trait,wcluster){


   #---------------------------------------------------------------------------------------#
   #    Make sure dat is a data frame, or at least it can be coerced to a data frame.      #
   #---------------------------------------------------------------------------------------#
   if (! is.data.frame(datum)){
      #----- Try turning dat into a data frame. -------------------------------------------#
      datum.try = try(as.data.frame(datum))
      if ("try-error" %in% is(datum.try)){
         stop("Argument \"datum\" must be coercible to a data frame.")
      }else{
         datum = datum.try
      }#end if
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Check dimensions of dat.   In case it has multiple rows, call it recursively, to   #
   # avoid unnecessary similarity comparisons amongst input rows (we only need to compare  #
   # with the medoids.                                                                     #
   #---------------------------------------------------------------------------------------#
   if (nrow(datum) > 1){
      ans = mapply( FUN = which.pft
                  , dat = split(x=datum,f=sequence(nrow(datum)))
                  , MoreArgs = list( medoid   = medoid
                                   , trait    = trait
                                   , wcluster = wcluster
                                   )#end list
                  )#end mapply
   }else{
      #----- Run the cluster analysis with medoids and the current datum. -----------------#
      tcluster = rbind(medoid[,trait$vnam[trait$cluster],drop=FALSE],datum)
      ndatum   = nrow(tcluster)
      #------------------------------------------------------------------------------------#


      #----- Turn categorical variables into logical. -------------------------------------#
      matmod   = as.formula(paste0("~ 0 + ",paste(names(tcluster),collapse=" + ")))
      mcluster = model.matrix(frml=matmod,data=tcluster)
      #------------------------------------------------------------------------------------#


      #----- Prepare data types. ----------------------------------------------------------#
      ordvnam  = trait$vnam[ (! trait$vlog) & trait$type %in% "numeric"   ]
      logvnam  = trait$vnam[ trait$vlog     & trait$type %in% "numeric"   ]
      iordtr   = which(names(mcluster) %in% ordvnam)
      ilogtr   = which(names(mcluster) %in% logvnam)
      iasytr   = seq_along(names(mcluster))[-c(iordtr,ilogtr)]
      #------------------------------------------------------------------------------------#




      #----- Calculate dissimilarity matrix. ----------------------------------------------#
      adissmat = daisy( x       = mcluster
                      , metric  = "gower"
                      , type    = list(ordratio=iordtr,logratio=ilogtr,asymm=iasytr)
                      , weights = wcluster
                      )#end daisy
      #------------------------------------------------------------------------------------#



      #----- Return the PFT of the most similar medoid. -----------------------------------#
      idx = which.min(as.matrix(adissmat)[ndatum,-ndatum])
      ans = medoid$pft[idx]
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#

   #----- Return answer. ------------------------------------------------------------------#
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function which.pft
#==========================================================================================#
#==========================================================================================#
