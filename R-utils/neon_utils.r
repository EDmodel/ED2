#----- List of flags for undetermined species. --------------------------------------------#
unknown_wildcard     <<- c("aff","cf","deleteme","ind","indet","na","ni","sp","spp"
                          ,"spnov","unknown","unidentified"
                          ,paste0("sp"  ,sequence(99))
                          ,paste0("spb" ,sequence(99))
                          ,paste0("spp" ,sequence(99))
                          ,paste0("sp"  ,sequence(99),"cay-atdn")
                          ,paste0("sp"  ,sequence(99),"-cay"    )
                          ,paste0("sp"  ,sequence(99),"guyafor" )
                          ,paste0("spfg",sequence(99),"-holst"  )
                          )#end c
unknown_common       <<- "unknown"
unk_liana_common     <<- "liana"
unknown_phylum       <<- "Ignotophyta"
unknown_class        <<- "Ignotopsida"
unknown_order        <<- "Ignotales"
unknown_family       <<- "Ignotaceae"
unknown_genus        <<- "Ignotum"
unknown_epithet      <<- "indet"
unk_liana_phylum     <<- "Lianophyta"
unk_liana_class      <<- "Lianopsida"
unk_liana_order      <<- "Lianales"
unk_liana_family     <<- "Lianaceae"
unk_liana_genus      <<- "Liana"
unknown_scientific   <<- paste(unknown_genus,unknown_epithet)
unk_liana_scientific <<- paste(unk_liana_genus,unknown_epithet)
#------------------------------------------------------------------------------------------#



#==========================================================================================
#==========================================================================================
#      This function standardises the spelling of common names of trees.
#------------------------------------------------------------------------------------------
standard_common_name_NEON <<- function(x){
   #---~---
   #   Make sure common names are lower case.
   #---~---
   x   = tolower(x)
   #---~---


   #---~---
   #   Remove underscore marks.
   #---~---
   x   = gsub(pattern="_",replacement=" ",x=x)
   #---~---



   #---~---
   #     General substitutions.  These are very aggressive, so don't use it too much.
   # Good things to put here are names that are often misspelt in a way that cannot occur
   # in other words. For example, abiuarana instead of abiurana is a good case; replacing
   # abiu with abiurana is a bad idea because abiurana would become abiuranarana. It's 
   # wise to use regexpr rules such as ^ and $ to make sure gsub won't substitute more
   # than it is supposed to.
   #---~---
   # x = gsub(pattern="^abiuarana"      ,replacement="abiurana"    ,x=x)
   #---~---



   #---~---
   #   Specific substitutions. Most cases should come here.
   #---~---
   sel = (x %in% "dense sedge"                              ); x[sel] = NA_character_
   #---~---

   return(x)
}#end function standard.common.name
#==========================================================================================
#==========================================================================================






#==========================================================================================
#==========================================================================================
#      This function corrects scientific names and families that are not correctly typed,
# are synonyms or have become obsolete.
#------------------------------------------------------------------------------------------#






#==========================================================================================
#==========================================================================================
#      This attributes scientific names based on common names for NEON surveys.
# It is NOT a good idea to use this anywhere else because common names may mean completely
# diferent things...
#------------------------------------------------------------------------------------------
scientific_lookup_NEON <<- function(datum,lookup_path,site){

   #---~---
   #   Read in the look-up table.
   #---~---
   lookup_file        = file.path(lookup_path,"NEON_taxon_lookup.csv")
   look_up            = as.data.table(read.csv(file=lookup_file,stringsAsFactors=FALSE))
   look_up$common     = tolower(trim(look_up$common    ))
   look_up$scientific = trim(look_up$scientific)
   look_up$family     = trim(look_up$family    )
   #---~---


   #---~---
   #   Keep only entries for this NEON site.
   #---~---
   if (site %in% names(look_up)){
      look_up = look_up[look_up[[site]],,drop=FALSE]
   }else{
      cat0("-------------------------------------------------------------------")
      cat0("  ERROR! Site is not available at the look-up data base."           )
      cat0("-------------------------------------------------------------------")
      cat0(" "                                                                  )
      cat0(" Site:           ",site                                             )
      cat0(" Data base file: ",lookup_file                                      )
      cat0(" "                                                                  )
      cat0(" Edit data base file, adding all the common names occurring at"     )
      cat0(" your site, and add a column named site with the entries that"      )
      cat0(" should be applicable to your site. Look for existing nearby"       )
      cat0(" sites, it is fine to have the same common name valid for multiple" )
      cat0(" locations.  If the common name already exists in the data base"    )
      cat0(" but refers to a different species in your site, duplicate the"     )
      cat0(" common name, make sure to set the new entry to TRUE for the new"   )
      cat0(" site, and FALSE to all the other sites. Also set the new site to"  )
      cat0(" FALSE in the existing entry."                                      )
      cat0(" "                                                                  )
      cat0("-------------------------------------------------------------------")
      stop("Invalid site.")
   }#end if (site %in% look_up)
   #---~---



   #---~---
   #   Break into genus and epithet.
   #---~---
   ge_list                  = sapply(X = tolower(look_up$scientific),FUN=strsplit,split=" ")
   ge_length                = sapply(X = ge_list, FUN = length)
   ge_mat                   = cbind( mapply(FUN="[",ge_list,MoreArgs=list(1))
                                   , mapply(FUN="[",ge_list,MoreArgs=list(2))
                                   )#end cbind
   g                        = capwords(ge_mat[,1],strict=TRUE)
   e                        = tolower(ge_mat[,2])
   g_e                      = paste(g,e,sep=" ")
   g_e[is.na(g) & is.na(e)] = NA_character_
   look_up$scientific       = g_e
   look_up$genus            = g
   #---~---



   #---~---
   #   Trim the common names, and simplify/replace some names.
   #---~---
   datum$common                      = tolower(trim(datum$common))
   datum$common[is.na(datum$common)] = unknown_common
   #---~---


   #---~---
   #   Initialise column gf.scientific in case it is not available.
   #---~---
   if (! "gf.scientific" %in% names(datum)){
      datum$gf.scientific = rep(NA_integer_,times=nrow(datum))
   }#end (! "gf.scientific" %in% names(datum))
   #---~---


   #---~---
   #   Find all unique common names.
   #---~---
   unique_common = unique(datum$common)
   n_common      = length(unique_common)
   not_found     = character(0L)
   for (n in sequence(n_common)){
      #---~---
      #   Find the trees that have the same common name as this one.
      #---~---
      cat0(" - ",n,"/",n_common," -- ",unique_common[n],".")
      w_dat  = which(datum$common %in% unique_common[n])
      n_dat  = length(w_dat)
      #---~---


      #---~---
      #   Find the trees in the look-up table with the same common name.
      #---~---
      w_look = which(look_up$common %in% unique_common[n])
      n_look = length(w_look)
      #---~---



      #---~---
      #   Check how many trees have the same common name in the look-up table.
      #---~---
      if (n_look == 1){
         #---~---
         #   Only one.  Use it.
         #---~---
         datum$scientific   [w_dat] = look_up$scientific[w_look]
         datum$genus        [w_dat] = look_up$genus     [w_look]
         datum$gf.scientific[w_dat] = 1L
      }else if (n_look > 1){
         datum$scientific   [w_dat] = sample( x       = look_up$scientific[w_look]
                                            , size    = n_dat
                                            , replace = TRUE
                                            )#end sample
         datum$genus        [w_dat] = look_up$genus     [w_look]
         datum$gf.scientific[w_dat] = 1
      }else{
         not_found                  = c(not_found,unique_common[n])
         datum$scientific   [w_dat] = unknown_scientific
         datum$genus        [w_dat] = unknown_genus
         datum$gf.scientific[w_dat] = 0
      }#end if
      #---~---
   }#end for (n in sequence(n_common))
   #---~---


   #---~---
   #   In case species were not found, list them so the user is aware.
   #---~---
   not_found = not_found[! not_found %in% unknown_common]
   if (length(not_found) %gt% 0L){
      #---~---
      #   Report the unmatched names.
      #---~---
      cat0(" "                                                                  )
      cat0(" "                                                                  )
      cat0("-------------------------------------------------------------------")
      cat0("  WARNING! Some names were not found in the common name data base." )
      cat0("-------------------------------------------------------------------")
      cat0(" "                                                                  )
      cat0(" Site:           ",site                                             )
      cat0(" Data base file: ",lookup_file                                      )
      cat0(" "                                                                  )
      cat0(" Species not found:"                                                )
      #---~---

      #---~---
      #   Loop through unmatched names.
      #---~---
      for (n in seq_along(not_found)){
         cat0(" - ",not_found[n],"."                                            )
      }#end for (n in seq_along(not_found))
      #---~---

      #---~---
      #   End report
      #---~---
      cat0(" "                                                                  )
      cat0("-------------------------------------------------------------------")
      cat0(" "                                                                  )
      cat0(" "                                                                  )
      #---~---
   }#end if (length(not_found) %gt% 0L)
   #---~---

   #---~---
   #   Return standardised data.
   #---~---
   return(datum)
   #---~---
}#end function scientific_lookup_NEON
#==========================================================================================#
#==========================================================================================#






#==========================================================================================
#==========================================================================================
#     Fill in the traits for all individuals.  We do this in three stages, and save how 
# the trait was determined.  The numbers represent the flag given for each trait.
#  0 -- Individuals have full identification and the species is listed in the
#       database; we used the reported density for that species.
#  1 -- The species/genus is identified but it isn't listed in the database, we use an
#       average of all other species of that genus that exist in database.
#  2 -- We could not identify the individual to the genus level, but we know the
#       family.  We use the average wood density of all individuals of this census that
#       belong to that family.
#  3 -- No taxonomic information could be retrieved for this individual, so we filled
#       with random sampling.
#
#  We also add 10 when the scientific name was gap filled.
#
#  INPUT variables:
#
#  - datum       -- data frame with data.  This is going to be the output as well
#  - trait       -- trait to fill
#  - tdb.csv     -- csv file containing the trait data base
#  - country     -- restrict attribution to these countries (NULL uses all countries)
#  - continents  -- restrict attribution to these continents (NULL uses all continents)
#  - fsample     -- use random sampling to fill unidentified individuals, or individuals
#                   with genus that is not available at the trait data base?
#                   If FALSE then it uses averages
#  - weight      -- A weighting factor to give either probability or to weight
#                   the average.  This could be a vector with weights, or a character with
#                   the name of the variable in datum to use as the weight, or an integer
#                   with the column to be used as the weighting factor
#  - verbose     -- Flag to control the amount of information
#------------------------------------------------------------------------------------------
find_trait_NEON <<- function( datum
                            , trait     = c( "wood.dens", "SLA", "leaf.phen"
                                           , "growth.form", "woodland")
                            , tdb.csv   = file.path(srcdir,"NEON_trait_table.csv")
                            , country   = NULL
                            , continent = NULL
                            , fsample   = TRUE
                            , weight    = NULL
                            , verbose   = FALSE
                            ){

   #---~---
   #   Check that trait is valid.
   #---~---
   trait    = match.arg(trait)
   n.trait  = paste0("n." ,trait)
   gf.trait = paste0("gf.",trait)
   #---~---


   #---~---
   #   Initialise gap filling flag and region in case they are not there.
   #---~---
   if (! trait %in% names(datum)){
      if ( trait %in% c("leaf.phen","growth.form") ){
         datum[[trait]]    = rep(NA_character_,times=nrow(datum))
      }else if ( trait %in% "woodland" ){
         datum[[trait]]    = rep(NA_logical_  ,times=nrow(datum))
      }else{
         datum[[trait]]    = rep(NA_real_     ,times=nrow(datum))
      }#end if ( trait %in% c("leaf.phen","growth.form") )
   }#end (! trait %in% names(datum))
   if (! gf.trait %in% names(datum)){
      datum[[gf.trait]] = rep(NA_integer_,times=nrow(datum))
   }#end (! gf.trait %in% names(datum))
   if (! "country" %in% names(datum)){
      datum$country = rep(NA_character_,times=nrow(datum))
   }#end (! "country" %in% names(datum))
   if (! "continent" %in% names(datum)){
      datum$continent = rep(NA_character_,times=nrow(datum))
   }#end (! "continent" %in% names(datum))
   #---~---


   #---~---
   #   Check whether weight is a vector, or a character.  Make them a vector here.
   #---~---
   if (is.null(weight)){
      #---~---
      #   No weight provided, use equal weights.
      #---~---
      wgtfac = rep(x=1/nrow(datum),times=nrow(datum))
      #---~---
   }else if (is.character(weight) && (length(weight) == 1)){
      #---~---
      #   Character with column name was provided.
      #---~---
      if (weight %in% names(datum)){
         wgtfac = ifelse(datum[[weight]] %gt% 0, datum[[weight]], 0)
      }else{
         stop(paste0(" Weight Variable name (",weight,") not found in datum!"))
      }#end if
      #---~---
   }else if (is.numeric(weight) && (length(weight) == 1)){
      #---~---
      #   Column index provided.
      #---~---
      if (! (weight %wr% c(1,ncol(datum)))){
         stop(paste0(" Weight column index (",weight,") doesn't make sense"))
      }else if (is.numeric(datum[,weight])){
         wgtfac = ifelse(datum[,weight] %gt% 0, datum[,weight], 0)
      }else{
         stop(paste0(" Column ",weight," of data frame is not numeric!"))
      }#end if
      #---~---
   }else if (is.numeric(weight) && (length(weight) == nrow(datum))){
      wgtfac = ifelse(weight %gt% 0, weight, 0)
   }else{
      stop("Variable weight is not properly set!")
   }#end if (is.null(weight))
   #---~---



   #---~---
   #   Check that the weighting factor makes sense.
   #---~---
   if (sum(wgtfac) > 0){
      wgtfac = wgtfac / sum(wgtfac)
   }else{
      stop(" Invalid weighting variable. Most numbers should be positive...")
   }#end if
   #---~---



   #---~---
   #   Read data base to fill traits.
   #---~---
   tdb  = read.csv(file=tdb.csv,header=TRUE,stringsAsFactors=FALSE)
   keep = ! is.na(tdb[[trait]])
   tdb  = tdb[keep,,drop=FALSE]
   #---~---



   #---~---
   #   Restrict data base to valid traits in selected regions, and make sure the data base
   # can be used.
   #---~---
   if (! (trait %in% names(tdb))){
      stop(paste0(" Trait ",trait," is not available in trait file \""
                 ,basename(tdb.csv),"\"."))
   }else if (! all(c("country","continent") %in% names(tdb))){
      stop(paste0(" \"country\" and/or \"continent\" missing in trait file \""
                 ,basename(tdb.csv),"\"."))
   }else{
      #---~---
      #   Restrict species to regions of interest.
      #---~---
      if (is.null(country)){
         sel.country = rep(x=TRUE,times=nrow(tdb))
      }else{
         sel.country = tdb$country %in% country
      }#end if (is.null(country))
      if (is.null(continent)){
         sel.continent = rep(x=TRUE,times=nrow(tdb))
      }else{
         sel.continent = tdb$continent %in% continent
      }#end if (is.null(continent))
      keep = (! is.na(tdb[[trait]]) ) & sel.country & sel.continent
      ntdb = sum(keep)
      if (ntdb == 0){
         cat0(" - Selected trait: ",trait)
         if (! is.null(country)){
            cat0(" - Selected countries: " ,paste(country  ,collapse="; "))
         }#end if (! is.null(country))
         if (! is.null(continent)){
            cat0(" - Selected continents: ",paste(continent,collapse="; "))
         }#end if (! is.null(continent))
         stop(paste0("No data available for trait ",trait," in file \""
                    ,basename(tdb.csv),"\"."))
      }else{
         tdb  = tdb[keep,,drop=FALSE]
      }#end if (ntdb == 0)
      #---~---
   }#end if (! (trait %in% names(tdb)))
   #---~---


   #=======================================================================================
   #=======================================================================================
   #     First loop, fill in information to all individuals for which the genus is
   # known.
   #---------------------------------------------------------------------------------------
   #---~---
   #   Separate all species.
   #---~---
   species          = unique(datum$scientific)
   nspecies         = length(species)
   sci.genus.mean   = matrix(nrow=0,ncol=3
                            ,dimnames=list(NULL,c("scientific","genus","family")))
   sci.loose.mean   = matrix(nrow=0,ncol=3
                            ,dimnames=list(NULL,c("scientific","genus","family")))
   #---~---


   #---~---
   #   Loop through species
   #---~---
   for (s in sequence(nspecies)){
      if ((! is.na(species[s])) && length(grep(unknown_genus,species[s])) == 0){
         if (verbose) cat0("   - ",s,"/",nspecies,"  -  ",species[s],".") 

         #---~---
         #   Get the genus and family of this species, in case we need it.
         #---~---
         igen        = which (datum$scientific %in% species[s])
         this.genus  = unique(datum$genus[igen])
         this.family = unique(capwords(datum$family[igen],strict=TRUE))
         if (length(this.genus) != 1 || length(this.family) != 1){
            cat0(" - Perhaps a bastard genus?")
            browser()
         }#end if
         #---~---



         #---~---
         #     Check whether we have biomass for this species.
         #---~---
         if (species[s] %in% tdb$scientific){
            #---~---
            #     Find the index of this species in the database, and assign the trait
            # from there.
            #---~---
            itrait = intersect( which(tdb$scientific %in% species[s] )
                              , intersect( which(tdb$genus  %in% this.genus )
                                         , which(tdb$family %in% this.family) ) )

            if (length(itrait) != 1){
               #---~---
               #   Intersection was zero! Misidentification, maybe?
               #---~---
               loose      = TRUE
               cat0("Weird, length(itrait)=",length(itrait),"!")
               browser()
               #---~---
            }else{
               if (any(c(FALSE,! is.na(tdb[[trait]][itrait])),na.rm=TRUE)){
                  sel                    = datum$scientific %in% species[s]
                  datum[[trait]]   [sel] = tdb[[trait]][itrait]
                  datum[[gf.trait]][sel] = 0L
                  datum$country    [sel] = tdb$country  [itrait]
                  datum$continent  [sel] = tdb$continent[itrait]
                  fill.genus              = FALSE
                  loose                   = FALSE
               }else{
                  loose                   = TRUE
                  fill.genus              = this.genus %in% tdb$genus
               }#end if
               #---~---
            }#end if
            #---~---
         }else{
            loose      = TRUE
            fill.genus = this.genus %in% tdb$genus
         }#end if
         #---~---



         #---~---
         #   If this species is to be filled with genus average.
         #---~---
         if (fill.genus){
            #---~---
            #     Find all the plants from this genus, take the average, and attribute
            # to this species.  In this case, warn the user about the species, it may
            # be a typo in the scientific name that is easy to fix.
            #---~---
            itrait  = intersect( which( tdb$genus  %in% this.genus )
                               , which( tdb$family %in% this.family)
                              )#end intersect

            #---~---
            #   Find individuals that belong to this species.
            #---~---
            sel  = datum$scientific %in% species[s]
            nsel = sum(sel)
            #---~---

            if (length(itrait) == 0){
               #---~---
               #   Intersection was zero! Misidentification, maybe?
               #---~---
               loose      = TRUE
               #---~---
            }else if (nsel > 0){
               #---~---
               #   Select data that can be used for weighting.
               #---~---
               t.value     = tdb[[trait  ]][itrait]
               t.weight    = tdb[[n.trait]][itrait] / sum(tdb[[n.trait]][itrait])
               t.continent = tdb$continent [itrait]
               t.country   = tdb$country   [itrait]
               t.idx       = seq_along(t.value)
               if (any(is.na(t.weight))) browser()
               #---~---


               #---~---
               #   Decide whether to use sample (non-numeric) or average (numeric).
               #---~---
               if (fsample || ( trait %in% c("leaf.phen","growth.form","woodland") )){
                  smp.idx       = lit.sample(x=t.idx,size=nsel,replace=TRUE,prob=t.weight)
                  smp.trait     = t.value    [smp.idx]
                  smp.country   = t.country  [smp.idx]
                  smp.continent = t.continent[smp.idx]
               }else{
                  smp.trait     = weighted.mean     (x=t.value    ,w=t.weight,na.rm=TRUE)
                  smp.country   = weighted.commonest(x=t.country  ,w=t.weight,na.rm=TRUE)
                  smp.continent = weighted.commonest(x=t.continent,w=t.weight,na.rm=TRUE)
               }#end if (fsample || ( trait %in% c("leaf.phen") ))
               #---~---


               #---~---
               #   Fill in missing entries.
               #---~---
               datum[[trait]]   [sel] = smp.trait
               datum[[gf.trait]][sel] = 1L
               datum$country    [sel] = smp.country
               datum$continent  [sel] = smp.continent
               sci.genus.mean         = rbind(sci.genus.mean
                                             ,c(species[s],this.genus,this.family))
               if (verbose){
                  cat0("     * Species ",species[s]," not found."
                                        ,"  Use genus average instead.")
               }#end if
               loose = FALSE
               #---~---
            }#end if
            #---~---
         }#end if
         #---~---


         #---~---
         #    Append plants that have no family together.
         #---~---
         if (loose){
            #---~---
            #   The plant probably doesn't have any family.
            #---~---
            sci.loose.mean = rbind(sci.loose.mean
                                  ,c(species[s],this.genus,this.family))
            #---~---
         }#end if
         #---~---
      }#end if
      #---~---
   }#end for
   #---~---



   #---~---
   #   List all genera that were not filled with species information.
   #---~---
   sci.genus.mean = sci.genus.mean[order(sci.genus.mean[,"scientific"]),,drop=FALSE]
   genus.only     = grepl( pattern     = " NA$"
                         , x           = sci.genus.mean[,"scientific",drop=FALSE]
                         , ignore.case = TRUE
                         )#end grepl
   if (verbose && (nrow(sci.genus.mean[!genus.only,,drop=FALSE]) > 0)){
      cat0("")
      cat0("-----------------------------------------------------------------------")
      cat0(" Found species that were not in trait data base (check for synonyms)!")
      print(sci.genus.mean[! genus.only,,drop=FALSE],quote=FALSE)
      cat0("-----------------------------------------------------------------------")
      cat0("")
   }#end if (verbose && nrow(sci.genus.mean[!genus.only,,drop=FALSE]) > 0)
   if (verbose && (nrow(sci.genus.mean[genus.only,,drop=FALSE]) > 0)){
      cat0("")
      cat0("-----------------------------------------------------------------------")
      cat0("Only genus was provided: filled with average genus value:")
      print (sci.genus.mean[genus.only,,drop=FALSE],quote=FALSE)
      cat0("-----------------------------------------------------------------------")
      cat0("")
   }#end if (verbose && nrow(sci.genus.mean[genus.only,,drop=FALSE]) > 0)
   #---~---


   #---~---
   #   List all genera that didn't belong to any known family.
   #---~---
   if (verbose && (nrow(sci.loose.mean) > 0)){
      cat0("")
      cat0("-----------------------------------------------------------------------")
      cat0(" Found genera from families with no valid data in the trait data base!")
      print(sci.loose.mean,quote=FALSE)
      cat0("-----------------------------------------------------------------------")
      cat0("")
   }#end if
   #---~---
   #=======================================================================================
   #=======================================================================================





   #=======================================================================================
   #=======================================================================================
   #     Second loop: we list all families, and look for individuals that have no genus
   # associated. We compute the mean trait of the known individuals for that family and
   # use that as an estimate of the trait.
   #---------------------------------------------------------------------------------------
   #---~---
   families          = unique(datum$family)
   nfamilies         = length(families)
   sci.family.sample = matrix(nrow=0,ncol=3
                             ,dimnames=list(NULL,c("scientific","family",trait)))
   #---~---


   #---~---
   #   Loop over all families
   #---~---
   for (f in sequence(nfamilies)){
      if (! (families[f] %in% unknown_family)){
         if (verbose) cat0("   - ",f,"/",nfamilies,"  -  ",families[f],".")

         #---~---
         #   Get the individuals that belong to this family.
         #---~---
         ifam  = which (datum$family %in% families[f] & (! is.na(datum[[trait]])) )
         imiss = which (datum$family %in% families[f] &    is.na(datum[[trait]])  )
         nmiss = length(imiss)
         #---~---

         if (length(imiss) > 0 && length(ifam) > 0){
            #---~---
            #   Select data that can be used for weighting.
            #---~---
            t.value     = datum[[trait]] [ifam]
            t.weight    = wgtfac         [ifam]
            t.continent = datum$continent[ifam]
            t.country   = datum$country  [ifam]
            t.idx       = seq_along(t.value)
            #---~---



            #---~---
            #   Decide whether to use sample or average.
            #---~---
            if (fsample || ( trait %in% c("leaf.phen","growth.form","woodland") )){
               smp.idx       = lit.sample(x=t.idx,size=nmiss,replace=TRUE,prob=t.weight)
               smp.trait     = t.value    [smp.idx]
               smp.country   = t.country  [smp.idx]
               smp.continent = t.continent[smp.idx]
            }else{
               smp.trait     = weighted.mean     (x=t.value    ,w=t.weight,na.rm=TRUE)
               smp.country   = weighted.commonest(x=t.country  ,w=t.weight,na.rm=TRUE)
               smp.continent = weighted.commonest(x=t.continent,w=t.weight,na.rm=TRUE)
            }#end if (fsample || ( trait %in% c("leaf.phen") ))
            #---~---



            #---~---
            #   Fill in with the sample/mean.
            #---~---
            datum[[trait]]   [imiss] = smp.trait
            datum[[gf.trait]][imiss] = 2L
            datum$country    [imiss] = smp.country
            datum$continent  [imiss] = smp.continent
            gf2                      = cbind(datum$scientific[imiss]
                                            ,datum$family    [imiss]
                                            ,sprintf("%6.3f",smp.trait)
                                            )#end cbind
            sci.family.sample        = rbind(sci.family.sample,gf2)
            #---~---
         }#end if
         #---~---
      }#end if
      #---~---
   }#end for
   #---~---


   #---~---
   #     Stop if there was any genus that didn't belong to any known family.
   #---~---
   if (nrow(sci.family.sample) > 0){
      if (verbose){
         cat0(" Found families with unidentified genera!")
         print(sci.family.sample,quote=FALSE)
      }#end if
   }#end if
   #---~---


   #---~---
   #    Final block.  We fill in the trait for unknown individuals, by randomly sampling
   # from the individuals we know the density.
   #---~---
   imiss  = which(is.na(datum[[trait]]))
   nmiss  = length(imiss)
   nvalid = nrow(datum) - nmiss

   if (nmiss > 0){
      if (verbose){
         cat0(" The following families are filled with global sampling: ")
         fam.global.sampling = t(t(sort(unique(datum$family[imiss]))))
         print(fam.global.sampling,quote=FALSE)
      }#end if
      #---~---


      #---~---
      #   Check whether this is going to be a complete guess or something slightly more
      # elegant.
      #---~---
      if (nvalid == 0){
         #---~---
         #   Use any data from the trait data base.
         #---~---
         warning(" None of the trees are known!  Trait-filling is going to be very crude!")
         itrait      = which(tdb$life.type %in% "T")
         t.value     = tdb[[trait  ]][itrait]
         t.weight    = tdb[[n.trait]][itrait] / sum(tdb[[n.trait]][itrait])
         t.continent = tdb$continent [itrait]
         t.country   = tdb$country   [itrait]
         t.idx       = seq_along(t.value)
         #---~---
      }else{

         #---~---
         #   Select data that can be used for weighting.
         #---~---
         t.value     = datum[[trait]] [-imiss]
         t.weight    = wgtfac         [-imiss]
         t.continent = datum$continent[-imiss]
         t.country   = datum$country  [-imiss]
         t.idx       = seq_along(t.value)
         #---~---
      }#end if (nvalid == 0)
      #---~---


      #---~---
      #   Decide whether to use sample or averaged values.
      #---~---
      if (fsample || ( trait %in% c("leaf.phen") )){
         smp.idx       = lit.sample(x=t.idx,size=nmiss,replace=TRUE,prob=t.weight)
         smp.trait     = t.value    [smp.idx]
         smp.country   = t.country  [smp.idx]
         smp.continent = t.continent[smp.idx]
      }else{
         smp.trait     = weighted.mean     (x=t.value    ,w=t.weight,na.rm=TRUE)
         smp.country   = weighted.commonest(x=t.country  ,w=t.weight,na.rm=TRUE)
         smp.continent = weighted.commonest(x=t.continent,w=t.weight,na.rm=TRUE)
      }#end if (fsample || ( trait %in% c("leaf.phen") ))
      #---~---


      #---~---
      #   Fill the remaining gaps.
      #---~---
      datum[[trait]]   [imiss] = smp.trait
      datum$country    [imiss] = smp.country
      datum$continent  [imiss] = smp.continent
      datum[[gf.trait]][imiss] = 3L
      #---~---
   }#end if
   #---~---


   #---~---
   #   Adjust the gap-filling flag for trait by adding whether the scientific name
   # itself was gap-filled.
   #---~---
   if ("gf.scientific" %in% names(datum)){
      datum[[gf.trait]] = datum[[gf.trait]] + as.integer(10 * datum$gf.scientific)
   }#end if ("gf.scientific" %in% names(datum))
   #---~---


   #---~---
   #   Assign a plant functional type based on the wood density.
   #---~---
   if (trait %in% "wood.dens"){
      pft.cut    = cut(datum[[trait]],breaks=pft.breaks)
      pft.levels = levels(pft.cut)
      pft.idx    = match(pft.cut,pft.levels)
      datum$pft  = mypfts[pft.idx]
   }#end if (trait %in% "wood.dens")
   #---~---
   return(datum)
}#end function find_trait_NEON
#==========================================================================================
#==========================================================================================






#==========================================================================================
#==========================================================================================
#     This function assigns a PFT given the genus and the full data set.  This will
# use imputed data for those genera that do not have all traits.
#------------------------------------------------------------------------------------------
assign_pft_NEON <<- function(datum,path=srcdir,refsites,verbose=FALSE,...){
   #---~---
   #   Load the cluster analysis with the medoid matrix.
   #---~---
   neon_lut = read.csv( file             = file.path(srcdir,"NEON_taxon_lookup.csv")
                      , header           = TRUE
                      , stringsAsFactors = FALSE
                      )#end read.csv
   #---~---


   #---~---
   #   Select species from the sites of interest (it is fine to use more than one 
   # reference site, just mind that they should be somewhat similar.
   #---~---
   sel = rep(FALSE,times=nrow(neon_lut))
   for (r in seq_along(refsites)){
      site = refsites[r]
      sel  = sel | neon_lut[[site]]
   }#end for (r in seq_along(refsites))
   neon_lut = neon_lut[sel,,drop=FALSE]
   #---~---



   #---~---
   #     Initialise genus with the actual one, but replace those not find in the
   # look-up table.
   #---~---
   if ("scientific.name" %in% names(datum)) datum$scientific = datum$scientific.name
   if ("genus.name"      %in% names(datum)) datum$genus      = datum$genus.name
   if ("family.name"     %in% names(datum)) datum$family     = datum$family.name
   #---~---



   #---~---
   #    Initialise use.scientific with the actual scientific name.  We will replace
   # those that are not found in the look-up table for traits.
   #---~---
   ans = rep(NA_integer_,times=nrow(datum))
   #---~---


   #---~---
   #    First attempt, fill in the species that have an exact match.
   #---~---
   idx      = match(datum$scientific,neon_lut$scientific)
   sel      = ! is.na(idx)
   ans[sel] = neon_lut$pft[idx[sel]]
   #---~---


   #---~---
   #   Second attempt, find genera that can be filled with existing data.
   #---~---
   un_genus = sort(unique(datum$genus[is.na(ans)]))
   un_genus = un_genus[un_genus %in% neon_lut$genus]
   for (u in seq_along(un_genus)){
      #---~---
      #   Handy alias.
      #---~---
      genus_now = un_genus[u]
      #---~---


      #---~---
      #   Retrieve PFTs associated with this genus.
      #---~---
      pft_pool  = neon_lut$pft[neon_lut$genus %in% genus_now]
      #---~---


      #---~---
      #   Fill in missing data with random sampling.
      #---~---
      sel      = ( datum$genus %in% genus_now ) & is.na(ans)
      ans[sel] = lit.sample(x=pft_pool,size=sum(sel),replace=TRUE)
      #---~---
   }#end for (u in seq_along(un_genus))
   #---~---


   #---~---
   #   Third attempt, find families that can be filled with existing data.
   #---~---
   un_family = sort(unique(datum$family[is.na(ans)]))
   un_family = un_genus[un_family %in% neon_lut$family]
   for (u in seq_along(un_family)){
      #---~---
      #   Handy alias.
      #---~---
      family_now = un_family[u]
      #---~---

      #---~---
      #   Retrieve PFTs associated with this family. We include any existing data from 
      # this family that has been already filled, as families can be broad and it may
      # be better to assume unidentified species are close to the identified species.
      #---~---
      pft_pool = c( neon_lut$pft[neon_lut$family %in% family_now]
                  , ans[ ( datum$family %in% family_now ) & (! is.na(ans))]
                  )#end c
      #---~---


      #---~---
      #   Fill in missing data with random sampling.
      #---~---
      sel      = ( datum$family %in% family_now ) & is.na(ans)
      ans[sel] = lit.sample(x=pft_pool,size=sum(sel),replace=TRUE)
      #---~---
   }#end for (u in seq_along(un_family))
   #---~---


   #---~---
   #   Last attempt, sample from within the entries that have been identified.
   #---~---
   sel      = is.na(ans)
   ans[sel] = lit.sample(x=ans[! sel],size=sum(sel),replace=TRUE)
   
   un_family = sort(unique(datum$family[is.na(ans)]))
   un_family = un_genus[un_family %in% neon_lut$family]
   for (u in seq_along(un_family)){
      #---~---
      #   Handy alias.
      #---~---
      family_now = un_family[u]
      #---~---

      #---~---
      #   Retrieve PFTs associated with this genus.
      #---~---
      pft_pool  = neon_lut$pft[neon_lut$family %in% family_now]
      #---~---


      #---~---
      #   Fill in missing data with random sampling.
      #---~---
      sel      = ( datum$family %in% family_now ) & is.na(ans)
      ans[sel] = lit.sample(pft_pool,size=sum(sel),replace=TRUE)
      #---~---
   }#end for (u in seq_along(un_family))
   #---~---




   #---~---
   #   Assign PFT
   #---~---
   return(ans)
   #---~---
}#end function assign_pft_NEON
#==========================================================================================
#==========================================================================================






#==========================================================================================
#==========================================================================================
#     Biomass allometry that is used for the NEON plots.  Results are always in kgC/plant.
#
# References:
#
# Chojnacky DC, Heath LS , Jenkins JC. 2014. Updated generalized biomass equations for 
#    North American tree species. Forestry (Lond), 87: 129-151. 
#    doi:10.1093/forestry/cpt053.
#
# Lutz JA, Furniss TJ, Germain SJ, Becker KML, Blomdahl EM, Jeronimo SMA, Cansler CA, 
#    Freund JA, Swanson ME , Larson AJ. 2017. Shrub communities, spatial patterns, and
#    shrub- mediated tree mortality following reintroduced fire in Yosemite National Park,
#    California, USA. Fire Ecol., 13: 104-126. doi:10.4996/fireecology.1301104.
#
# Input:
# ----------------------------------------------------------------------------------------
# dxh        --- Diameter of reference (either DBH or basal diameter) [cm]
# wdens      --- Wood density [g/cm3]
# scientific --- Species
# genus      --- Genus
# family     --- Family
# type       --- Plant type:
#                B - broadleaf tree
#                N - needleleaf tree
#                S - shrubs
# dead       --- Life status:
#                TRUE  - plant is dead
#                FALSE - plant is alive
#                In case dead = NULL, all plants are assumed to be alive.
# eps.dbh    --- Relative uncertainty for DBH [1 means 100%]
# eps.height --- Relative uncertainty for height [1 means 100%]
# eps.wdens  --- Relative uncertainty for wood density [1 means 100%]
# out.err    --- Output error in addition to the estimates of biomass/necromass?
# ----------------------------------------------------------------------------------------
#
#
#
# ----------------------------------------------------------------------------------------
# Output:
# ----------------------------------------------------------------------------------------
# - In case out.err is FALSE, the function returns a vector with biomass for each entry.
# - In case out.err is TRUE, the output is a data frame with the following vectors
#   with the same length as the entries:
#   * agb      -- biomass (necromass)                       [kgC]
#   * ae.agb   -- uncertainty in biomass due to allometry   [kgC, not relative]
#   * me.agb   -- uncertainty in biomass due to measurement [kgC, not relative]
#   * lnagb    -- log(biomass), used for error propagation.
#   * sd.lnagb -- standard error of log-biomass
#------------------------------------------------------------------------------------------
agb_NEON <<- function( dbh
                     , ddh
                     , wdens
                     , scientific
                     , genus
                     , family
                     , type
                     , phenology
                     , woodland
                     , dead       = NULL
                     , eps.dxh    = 0.02
                     , eps.height = 0.167
                     , eps.wdens  = 0.10
                     , out.err    = FALSE
                     ){
   #---~---
   #   If variable "dead" is missing, assume all individuals are alive
   #---~---
   if (is.null(dead)) dead = rep(FALSE,times=length(dbh))
   #---~---



   #---~---
   #   Make sure all terms have the same length and correct type.
   #---~---
   lens = unique( c( length(dbh)     , length(ddh)   , length(wdens), length(scientific)
                   , length(genus)   , length(family), length(type) , length(phenology)
                   , length(woodland), length(dead)  ) )
   if ( length(lens) != 1 ){
      #---~---
      #   Stop if any variable has a different length.
      #---~---
      cat0("-----------------------------------------------------------")
      cat0("   Variables don't have the same length."                   )
      cat0("   DBH         = ",length(dbh)                              )
      cat0("   DDH         = ",length(ddh)                              )
      cat0("   WDENS       = ",length(wdens)                            )
      cat0("   SCIENTIFIC  = ",length(scientific)                       )
      cat0("   GENUS       = ",length(genus)                            )
      cat0("   FAMILY      = ",length(family)                           )
      cat0("   TYPE        = ",length(type)                             )
      cat0("   PHENOLOGY   = ",length(phenology)                        )
      cat0("   WOODLAND    = ",length(woodland)                         )
      cat0("   DEAD        = ",length(dead)                             )
      cat0("-----------------------------------------------------------")
      stop(" Incorrect input data.")
      #---~---
   }else{
      fine.dbh        = is.numeric  (dbh)        || all(is.na(dbh       ))
      fine.ddh        = is.numeric  (ddh)        || all(is.na(ddh       ))
      fine.wdens      = is.numeric  (wdens)      || all(is.na(wdens     ))
      fine.scientific = is.character(scientific) || all(is.na(scientific))
      fine.genus      = is.character(genus)      || all(is.na(genus     ))
      fine.family     = is.character(family)     || all(is.na(family    ))
      fine.type       = is.character(type)       || all(is.na(type      ))
      fine.phenology  = is.character(phenology)  || all(is.na(phenology ))
      fine.woodland   = is.logical  (woodland)   || all(is.na(woodland  ))
      fine.dead       = is.logical  (dead)       || all(is.na(dead      ))
      all.fine        = (  fine.dbh      && fine.ddh    && fine.wdens && fine.scientific
                        && fine.genus    && fine.family && fine.type  && fine.phenology
                        && fine.woodland && fine.dead
                        )#end all.fine
      #---~---
      #   Stop if anything is unexpected.
      #---~---
      if (! all.fine){
         cat0("-----------------------------------------------------------")
         cat0("   Not all variables have the correct type."                )
         cat0("   DBH        (numeric)   = ",fine.dbh                      )
         cat0("   DDH        (numeric)   = ",fine.ddh                      )
         cat0("   WDENS      (numeric)   = ",fine.wdens                    )
         cat0("   SCIENTIFIC (character) = ",fine.scientific               )
         cat0("   GENUS      (character) = ",fine.genus                    )
         cat0("   FAMILY     (character) = ",fine.family                   )
         cat0("   TYPE       (character) = ",fine.type                     )
         cat0("   PHENOLOGY  (character) = ",fine.phenology                )
         cat0("   WOODLAND   (logical)   = ",fine.woodland                 )
         cat0("   DEAD       (logical)   = ",fine.dead                     )
         cat0("-----------------------------------------------------------")
         stop(" Incorrect data types.")
      }#end if (! all.fine)
      #---~---
   }#end if ( length(lens) != 1)
   #---~---



   #---~---
   #   Initialise the output.
   #---~---
   a0.agb = NA_real_ * dbh
   a1.agb = NA_real_ * dbh
   #---~---



   #---~---
   #   Simplify phenology labels.
   #---~---
   phenology = toupper(substr(phenology,1,1))
   #---~---


   #---~---
   #   Link names to the equation sets.
   #---~---
   bleaf         = ( type %in% "B" ) & (! dead)
   nleaf         = ( type %in% "N" ) & (! dead)
   shrub         = ( type %in% "S" ) & (! dead)
   tree          = ( bleaf | nleaf ) & (! shrub)
   wlnd          = tree & woodland
   multigroup.01 = family %in% c( "Cornaceae"  , "Ericaceae", "Lauraceae"
                                , "Platanaceae", "Rosaceae" , "Ulmaceae" )
   multigroup.02 = family %in% c("Fabaceae","Juglandaceae")
   multigroup.03 = family %in% c("Hippocastanaceae","Tiliaceae")
   multigroup.04 = family %in% c("Fabaceae","Rosaceae")
   #---~---



   #---~---
   #   Define generic diameter (useful for shrubs).
   #---~---
   dxh = ifelse( test = shrub | wlnd
               , yes  = ifelse( test = is.finite(ddh), yes = ddh, no = dbh )
               , no   = ifelse( test = is.finite(dbh), yes = dbh, no = ddh )
               )#end ifelse
   #---~---


   #---~---
   #   Handles for the tree equations
   #---~---
   t.abies.lwr        = (genus  %in% "Abies"         ) & (wdens %lt% 0.35        ) & tree
   t.abies.upr        = (genus  %in% "Abies"         ) & (wdens %ge% 0.35        ) & tree
   t.cupressaceae.lwr = (family %in% "Cupressaceae"  ) & (wdens %lt% 0.30        ) & tree
   t.cupressaceae.mid = (family %in% "Cupressaceae"  ) & (wdens %wl% c(0.30,0.40)) & tree
   t.cupressaceae.upr = (family %in% "Cupressaceae"  ) & (wdens %ge% 0.40        ) & tree
   t.larix            = (genus  %in% "Larix"         )                             & tree
   t.picea.lwr        = (genus  %in% "Picea"         ) & (wdens %lt% 0.35        ) & tree
   t.picea.upr        = (genus  %in% "Picea"         ) & (wdens %ge% 0.35        ) & tree
   t.pinus.lwr        = (genus  %in% "Pinus"         ) & (wdens %lt% 0.45        ) & tree
   t.pinus.upr        = (genus  %in% "Pinus"         ) & (wdens %ge% 0.45        ) & tree
   t.pseudotsuga      = (genus  %in% "Pseudotsuga"   )                             & tree
   t.tsuga.lwr        = (genus  %in% "Tsuga"         ) & (wdens %lt% 0.40        ) & tree
   t.tsuga.upr        = (genus  %in% "Tsuga"         ) & (wdens %ge% 0.40        ) & tree
   t.aceraceae.lwr    = (family %in% "Aceraceae"     ) & (wdens %lt% 0.50        ) & tree
   t.aceraceae.upr    = (family %in% "Aceraceae"     ) & (wdens %ge% 0.50        ) & tree
   t.betulaceae.lwr   = (family %in% "Betulaceae"    ) & (wdens %lt% 0.40        ) & tree
   t.betulaceae.lmd   = (family %in% "Betulaceae"    ) & (wdens %wl% c(0.40,0.50)) & tree
   t.betulaceae.umd   = (family %in% "Betulaceae"    ) & (wdens %wl% c(0.50,0.60)) & tree
   t.betulaceae.upr   = (family %in% "Betulaceae"    ) & (wdens %ge% 0.60        ) & tree
   t.multigroup.01    = (family %in% multigroup.01   )                             & tree
   t.carya            = (genus  %in% "Carya"         )                             & tree
   t.multigroup.02    = (family %in% multigroup.02   ) & ( ! genus %in% "Carya"  ) & tree
   t.fagaceae.dcd     = (family %in% "Fagaceae"      ) & (phenology %in% "D"     ) & tree
   t.fagaceae.evg     = (family %in% "Fagaceae"      ) & (phenology %in% "E"     ) & tree
   t.hamamelidaceae   = (family %in% "Hamamelidaceae")                             & tree
   t.multigroup.03    = (family %in% multigroup.03   )                             & tree
   t.magnoliaceae     = (family %in% "Magnoliaceae"  )                             & tree
   t.oleaceae.lwr     = (family %in% "Oleaceae"      ) & (wdens %lt% 0.55        ) & tree
   t.oleaceae.upr     = (family %in% "Oleaceae"      ) & (wdens %ge% 0.55        ) & tree
   t.salicaceae.lwr   = (family %in% "Salicaceae"    ) & (wdens %lt% 0.35        ) & tree
   t.salicaceae.upr   = (family %in% "Salicaceae"    ) & (wdens %ge% 0.35        ) & tree
   w.cupressaceae     = (family %in% "Cupressaceae"  )                             & wlnd
   w.multigroup.04    = (family %in% multigroup.04   )                             & wlnd
   w.fagaceae         = (family %in% "Fagaceae"      )                             & wlnd
   w.pinaceae         = (family %in% "Pinaceae"      )                             & wlnd
   #---~---


   #---~---
   #   Handles for the shrub equations.
   #---~---
   s.arctostaphylos   = (genus %in% "Arctostaphylos") & shrub
   s.ceanothus        = (genus %in% "Ceanothus"     ) & shrub
   s.chrysolepis      = (genus %in% "Chrysolepis"   ) & shrub
   s.corylus          = (genus %in% "Corylus"       ) & shrub
   s.cornus           = (genus %in% "Cornus"        ) & shrub
   s.leucothoe        = (genus %in% "Leucothoe"     ) & shrub
   s.rhododendron     = (genus %in% "Rhododendron"  ) & shrub
   s.ribes            = (genus %in% "Ribes"         ) & shrub
   s.rosa             = (genus %in% "Rosa"          ) & shrub
   s.rubus            = (genus %in% "Rubus"         ) & shrub
   s.sambucus         = (genus %in% "Sambucus"      ) & shrub
   s.symphoricarpos   = (genus %in% "Symphoricarpos") & shrub
   s.vaccinium        = (genus %in% "Vaccinium"     ) & shrub
   #---~---



   #---~---
   #     Set coefficients. We initialise them with a generic equation that only 
   # distinguishes whether they are shrubs, broadleaf trees, or needleleaf trees. For 
   # trees, we currently used some common families, but this can be refined later.
   # After the initial assignment, we update the coefficients if we find a better match
   # for the tree/shrub.
   #---~---
   #--- General coefficients
   a0.agb[bleaf             ] = -2.2118 ; a1.agb[bleaf             ] = 2.4133
   a0.agb[nleaf             ] = -2.6177 ; a1.agb[nleaf             ] = 2.4638
   a0.agb[shrub             ] = -3.1478 ; a1.agb[shrub             ] = 2.3750
   #--- Specific coefficients
   a0.agb[t.abies.lwr       ] = -2.3123 ; a1.agb[t.abies.lwr       ] = 2.3482
   a0.agb[t.abies.upr       ] = -3.1774 ; a1.agb[t.abies.upr       ] = 2.6426
   a0.agb[t.cupressaceae.lwr] = -1.9615 ; a1.agb[t.cupressaceae.lwr] = 2.1063
   a0.agb[t.cupressaceae.mid] = -2.7765 ; a1.agb[t.cupressaceae.mid] = 2.4195
   a0.agb[t.cupressaceae.upr] = -2.6327 ; a1.agb[t.cupressaceae.upr] = 2.4757
   a0.agb[t.larix           ] = -2.3012 ; a1.agb[t.larix           ] = 2.3853
   a0.agb[t.picea.lwr       ] = -3.0300 ; a1.agb[t.picea.lwr       ] = 2.5567
   a0.agb[t.picea.upr       ] = -2.1364 ; a1.agb[t.picea.upr       ] = 2.3233
   a0.agb[t.pinus.lwr       ] = -2.6177 ; a1.agb[t.pinus.lwr       ] = 2.4638
   a0.agb[t.pinus.upr       ] = -3.0506 ; a1.agb[t.pinus.upr       ] = 2.6465
   a0.agb[t.pseudotsuga     ] = -2.4623 ; a1.agb[t.pseudotsuga     ] = 2.4852
   a0.agb[t.tsuga.lwr       ] = -2.3480 ; a1.agb[t.tsuga.lwr       ] = 2.3876
   a0.agb[t.tsuga.upr       ] = -2.9208 ; a1.agb[t.tsuga.upr       ] = 2.5697
   a0.agb[t.aceraceae.lwr   ] = -2.0470 ; a1.agb[t.aceraceae.lwr   ] = 2.3852
   a0.agb[t.aceraceae.upr   ] = -1.8011 ; a1.agb[t.aceraceae.upr   ] = 2.3852
   a0.agb[t.betulaceae.lwr  ] = -2.5932 ; a1.agb[t.betulaceae.lwr  ] = 2.5349
   a0.agb[t.betulaceae.lmd  ] = -2.2271 ; a1.agb[t.betulaceae.lmd  ] = 2.4513
   a0.agb[t.betulaceae.umd  ] = -1.8096 ; a1.agb[t.betulaceae.umd  ] = 2.3480
   a0.agb[t.betulaceae.upr  ] = -2.2652 ; a1.agb[t.betulaceae.upr  ] = 2.5349
   a0.agb[t.multigroup.01   ] = -2.2118 ; a1.agb[t.multigroup.01   ] = 2.4133
   a0.agb[t.carya           ] = -2.5095 ; a1.agb[t.carya           ] = 2.6175
   a0.agb[t.multigroup.02   ] = -2.5095 ; a1.agb[t.multigroup.02   ] = 2.5437
   a0.agb[t.fagaceae.dcd    ] = -2.0705 ; a1.agb[t.fagaceae.dcd    ] = 2.4410
   a0.agb[t.fagaceae.evg    ] = -2.2198 ; a1.agb[t.fagaceae.evg    ] = 2.4410
   a0.agb[t.hamamelidaceae  ] = -2.6390 ; a1.agb[t.hamamelidaceae  ] = 2.5466
   a0.agb[t.multigroup.03   ] = -2.4108 ; a1.agb[t.multigroup.03   ] = 2.4177
   a0.agb[t.magnoliaceae    ] = -2.5497 ; a1.agb[t.magnoliaceae    ] = 2.5011
   a0.agb[t.oleaceae.lwr    ] = -2.0314 ; a1.agb[t.oleaceae.lwr    ] = 2.3524
   a0.agb[t.oleaceae.upr    ] = -1.8384 ; a1.agb[t.oleaceae.upr    ] = 2.3524
   a0.agb[t.salicaceae.lwr  ] = -2.6863 ; a1.agb[t.salicaceae.lwr  ] = 2.4561
   a0.agb[t.salicaceae.upr  ] = -2.4441 ; a1.agb[t.salicaceae.upr  ] = 2.4561
   a0.agb[w.cupressaceae    ] = -2.7096 ; a1.agb[w.cupressaceae    ] = 2.1942
   a0.agb[w.multigroup.04   ] = -2.9255 ; a1.agb[w.multigroup.04   ] = 2.4109
   a0.agb[w.fagaceae        ] = -3.0304 ; a1.agb[w.fagaceae        ] = 2.4982
   a0.agb[w.pinaceae        ] = -3.2007 ; a1.agb[w.pinaceae        ] = 2.5339
   a0.agb[s.arctostaphylos  ] = -3.5892 ; a1.agb[s.arctostaphylos  ] = 2.6846
   a0.agb[s.ceanothus       ] = -3.2406 ; a1.agb[s.ceanothus       ] = 2.6502
   a0.agb[s.chrysolepis     ] = -3.0198 ; a1.agb[s.chrysolepis     ] = 2.3110
   a0.agb[s.corylus         ] = -3.3776 ; a1.agb[s.corylus         ] = 2.3720
   a0.agb[s.cornus          ] = -3.5928 ; a1.agb[s.cornus          ] = 2.6470
   a0.agb[s.leucothoe       ] = -4.1588 ; a1.agb[s.leucothoe       ] = 2.3060
   a0.agb[s.rhododendron    ] = -3.1478 ; a1.agb[s.rhododendron    ] = 2.3750
   a0.agb[s.ribes           ] = -3.1478 ; a1.agb[s.ribes           ] = 2.3750
   a0.agb[s.rosa            ] = -3.1478 ; a1.agb[s.rosa            ] = 2.3750
   a0.agb[s.rubus           ] = -3.1478 ; a1.agb[s.rubus           ] = 2.3750
   a0.agb[s.sambucus        ] = -3.3378 ; a1.agb[s.sambucus        ] = 2.3720
   a0.agb[s.symphoricarpos  ] = -3.1478 ; a1.agb[s.symphoricarpos  ] = 2.3750
   a0.agb[s.vaccinium       ] = -3.1478 ; a1.agb[s.vaccinium       ] = 2.3750
   #---~---


   #---~---
   #   Find biomass.
   #---~---
   agb    = exp(a0.agb + a1.agb * log(dxh)) / C2B
   #---~---


   #---~---
   #      Check whether to estimate associated errors (measurement and allometry),
   # following:
   #
   #   Chave, J., and co-authors, 2004: Error propagation and scaling for tropical forest
   #      biomass estimates. Phil. Trans. R. Soc. Lond. B., 359, 409-420.
   #      doi:10.1098/rstb.2003.1425
   #---~---
   if (out.err){
      #---~---
      #   Find error associated with allometry. The papers do not present the actual error
      # of each equation. For the time being, we use the same error as reported by 
      # Chave et al. (2014) as a zero-order guess.
      #
      # Chave, J., and co-authors, 2014: Improved allometric models to estimate the
      #     aboveground biomass of tropical trees.  Glob. Change Biol., 20, 3177-3190
      #     doi:10.1111/gcb.12629
      #---~---
      ae.agb = sqrt(exp(0.357^2)-1)*agb
      #---~---



      #---~---
      #   Find error associated with measurements.
      #---~---
      me.agb = agb * a1.agb * eps.dxh
      #---~---


      #---~---
      #   Save the standard error of the log scale: it will be useful for error analysis.
      # The 0*agb term will ensure that sd.lnagb will be NA when agb is NA.
      #---~---
      sd.lnagb = 0.357 + 0. * agb
      #---~---

      #---~---
      #   Combine estimates and errors in a data frame.
      #---~---
      ans = data.frame( agb      = agb
                      , a0       = a0.agb
                      , a1       = a1.agb
                      , ae.agb   = ae.agb
                      , me.agb   = me.agb
                      , lnagb    = log(agb) - 0.5 * sd.lnagb^2
                      , sd.lnagb = sd.lnagb
                      )#end data.frame
      #---~---

   }else{
      #---~---
      #   No error needed.  Return estimate only.
      #---~---
      ans = agb
      #---~---
   }#end if
   #---~---

   return(ans)
}#end function agb_NEON
#==========================================================================================
#==========================================================================================





#==========================================================================================
#==========================================================================================
#      This function determines the plot aggregated error due to measurement and
# allometry for NEON.
#
# Input:
# ----------------------------------------------------------------------------------------
# datum       - a data frame with all trees in this plot (note, you can't run this for all
#               plots at once, either use a for loop or mapply). The data frame must
#               contain the following variables:
#               * nplant   - 1/area sampled if tree is alive, zero if tree is dead [1/m2]
#               * ntotal   - 1/area sampled                                        [1/m2]
#                 For nplant/ntotal, if the sampled area was 2500 m2, the number should be
#                 0.0004.  In case a differential sampling effort was used, then nplant
#                 varies.  For example, if trees with 10 <= DBH < 35 cm were measured in a
#                 50x5 subplot and trees with DBH >= 35 cm were measured in the entire
#                 2500m2 plots, then nplant should be 0.004 for the trees with DBH < 35cm
#                 and 0.0004 for trees with DBH >= 35cm.
#               * AGC      - above-ground carbon                                   [ kgC]
#               * ME.AGC   - measurement uncertainty of above-ground carbon        [ kgC]
#               * LNAGC    - log of above-ground carbon
#               * SD.LNAGC - allometry uncertainty, using the log-scale
#               * X        - x position in the plot (used only if epsilon.smp is NULL)
#               * Y        - y position in the plot (used only if epsilon.smp is NULL)
# xmax        - maximum size along the x axis (used only if epsilon.smp is NULL)
# ymax        - maximum size along the y axis (used only if epsilon.smp is NULL)
# subalong    - which axis has the subplot
# epsilon.smp - in case epsilon.smp is null, the function will try to estimate within-plot
#               sampling uncertainty.  This is unlikely to work unless you have large
#               plots or at the very least have plots without sub-sampling.
#               Alternatively, you may provide the number from previous studies.
# n.sub       - number of subplot samples
# n.real      - number of replicates for estimating uncertainty.  Large numbers
#               (10000 or more) are needed for stable results.
# ----------------------------------------------------------------------------------------
#
#
#
# ----------------------------------------------------------------------------------------
# Output:
# ----------------------------------------------------------------------------------------
# A vector with 8 numbers:
# se.agb.xxxxx - uncertainties for plot estimate of biomass (no standing dead) [kgC/m2]
# se.acd.xxxxx - uncertainties for plot estimate of biomass+necromass          [kgC/m2]
# se.xxx.measurement - contribution of measurement uncertainty to plot uncertainty
# se.xxx.allometry   - contribution of allometry uncertainty to plot uncertainty
# se.xxx.sampling    - contribution of sampling uncertainty to plot uncertainty
# se.xxx.census      - total uncertainty (combining the three terms above).
#------------------------------------------------------------------------------------------
find_agb_error_NEON <<- function(datum,xmax,ymax,subalong=c("x","y"),epsilon.smp=NULL
                           ,n.sub=25,n.real=10000){

   #---~---
   #   Standardise subalong.
   #---~---
   subalong = match.arg(subalong)
   #---~---


   #---~---
   #   Number of data points.
   #---~---
   ndatum = nrow(datum)
   #---~---


   #---~---
   #   First, find the combined error.
   #---~---
   se.use = data.frame( measurement = sqrt(log(1+(datum$ME.AGC / datum$AGC)^2))
                      , allometry   = datum$SD.LNAGC
                      )#end datum
   #---~---


   #---~---
   #   Initialise answer.
   #---~---
   se.agb        = rep(x=NA_real_,length(se.use)+2)
   names(se.agb) = c(names(se.use),"sampling","census")
   #---~---



   #---~---
   #   Create a population matrix.
   #---~---
   NPLANT = matrix( data = rep(datum$nplant,times=n.real), nrow=ndatum,ncol=n.real)
   #---~---


   #---~---
   #   Initialise answer.  We will create the vector before and loop through the errors.
   #---~---
   for (e in seq_along(se.use)){
      AGC       = matrix( data = rlnorm( n = n.real*ndatum
                                       , meanlog = rep(datum$LNAGC,times=n.real)
                                       , sdlog   = rep(se.use[[e]],times=n.real)
                                       )#end rlnorm
                        , nrow = ndatum
                        , ncol = n.real
                        )#end matrix
      AGB       = colSums(NPLANT * AGC)
      se.agb[e] = sd(AGB)
   }#end for (e in seq_along(se.use))
   #---~---


   #---~---
   #   Find total biomass, to be used by the within plot sampling error.
   #---~---
   agb.bar = sum(datum$nplant * datum$AGC)
   #---~---


   #---~---
   #   Unless epsilon.smp is provided, we estimate the within-plot sampling error.
   #---~---
   if (is.null(epsilon.smp)){
      #---~---
      #   Split domain into smaller subsplots.
      #---~---
      if (subalong %in% "x"){
         xbreaks = seq(from=0,to=xmax,length.out=n.sub+1)
         xlwr    =       sqrt(.Machine$double.eps)  * xmax
         xupr    = (1. - sqrt(.Machine$double.eps)) * xmax
         xnow    = pmax(xlwr,pmin(xupr,datum$X))
         isub    = as.integer(cut(xnow,breaks=xbreaks))
      }else{
         ybreaks = seq(from=0,to=ymax,length.out=n.sub+1)
         ylwr    =       sqrt(.Machine$double.eps)  * ymax
         yupr    = (1. - sqrt(.Machine$double.eps)) * ymax
         ynow    = pmax(ylwr,pmin(yupr,datum$Y))
         isub    = as.integer(cut(ynow,breaks=ybreaks))
      }#end if
      #---~---


      #---~---
      #   Find the AGB/ACD for each subplot.
      #---~---
      agb.sub          = rep(0.,times=n.sub)
      acd.sub          = rep(0.,times=n.sub)
      agb.tmp          = tapply(X=n.sub*datum$nplant*datum$AGC,INDEX=isub,FUN=sum)
      acd.tmp          = tapply(X=n.sub*datum$ntotal*datum$AGC,INDEX=isub,FUN=sum)
      agb.idx          = as.integer(names(agb.tmp))
      acd.idx          = as.integer(names(acd.tmp))
      agb.sub[agb.idx] = agb.tmp
      acd.sub[acd.idx] = acd.tmp
      #---~---




      #---~---
      #   Create replicates using bootstrap with replacement.
      #---~---
      AGB = colMeans( matrix( data = sample(agb.sub,size=n.real*n.sub,replace=TRUE)
                            , nrow = n.sub
                            , ncol = n.real
                            )#end matrix
                    )#end colMeans
      se.agb["sampling"] = sd(AGB)
      #---~---
   }else{
      #---~---
      #   Use pre-defined sampling error.
      #---~---
      se.agb["sampling"] = epsilon.smp * agb.bar
      #---~---
   }#end if (is.null(epsilon.smp))
   #---~---



   #---~---
   #   Find combined source of errors.
   #---~---
   se.agb["census"] = sqrt(sum(se.agb[c("measurement","allometry","sampling")]^2))
   #---~---


   #---~---
   #   Error is the standard deviation of all realisations.
   #---~---
   ans        = se.agb
   names(ans) = paste0("se.agb.",names(se.agb))
   #---~---


   #---~---
   #   Return answer.
   #---~---
   return(ans)
   #---~---
}#end find_agb_error_NEON
#==========================================================================================
#==========================================================================================
