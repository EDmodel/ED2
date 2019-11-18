#==========================================================================================#
#==========================================================================================#
#     Function fuse.trees                                                                  #
#                                                                                          #
#     This function fuse one group of trees to the other, once the fusion indices have     #
# been already set.                                                                        #
#                                                                                          #
#     receptor        -- The data frame that will receive the data                         #
#     donor           -- The data frame that will send the data                            #
#     don.2.rec       -- The indices of the receptor that will receive the data.           #
#     update.yr.notes -- Should we use notes.orig to make year notes? (TRUE/FALSE).        #
#------------------------------------------------------------------------------------------#
fuse.trees <<- function(receptor,donor,don.2.rec,survey.years,sci.strict=TRUE){

   #---------------------------------------------------------------------------------------#
   #       Get the list of all variables.                                                  #
   #---------------------------------------------------------------------------------------#
   variables   = names(receptor)
   n.variables = length(variables)
   n.receptor  = nrow(receptor)
   empty       = rep(NA,n.receptor)
   try.genus   = "genus" %in% names(donor)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Copy the scientific name and family.                                              #
   #---------------------------------------------------------------------------------------#
   receptor.s.ignotum = receptor$scientific %in% c(unknown.scientific,NA_character_)
   receptor.o.ignotum = receptor$old.common %in% c(unknown.common    ,NA_character_)
   receptor.c.ignotum = receptor$common     %in% c(unknown.common    ,NA_character_)
   receptor.f.ignotum = receptor$common.1st %in% c(unknown.common    ,NA_character_)
   this.1st           = empty; this.1st[don.2.rec] = donor$common.1st
   this.sci           = empty; this.sci[don.2.rec] = donor$scientific
   this.fam           = empty; this.fam[don.2.rec] = donor$family
   this.s.ignotum     = this.sci            %in% c(unknown.scientific,NA_character_)
   this.f.ignotum     = this.1st            %in% c(unknown.common    ,NA_character_)
   if (try.genus){
      this.gen  = empty; this.gen[don.2.rec] = donor$genus
      this.full = ifelse(this.s.ignotum,FALSE, this.sci != this.gen)
   }else{
      this.gen  = empty
      this.full = ! this.s.ignotum
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Copy new data to the receptor dataset, one variable at a time as they have       #
   # different rules.                                                                      #
   #---------------------------------------------------------------------------------------#
   for (nv in sequence(n.variables)){
      #------------------------------------------------------------------------------------#
      #       Grab the variable name, and make some auxiliary names.                       #
      #------------------------------------------------------------------------------------#
      vname     = variables[nv]
      old       = paste("old",vname,sep=".")
      first     = paste(vname,"1st",sep=".")
      new.input = paste(substring(vname,1,nchar(vname)-5),survey.years,sep=".")
      #------------------------------------------------------------------------------------#

      #----- Save variable to the template. -----------------------------------------------#
      this.val            = empty
      this.val[don.2.rec] = donor[[vname]]
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Decide what to do based on the variable.                                      #
      #------------------------------------------------------------------------------------#
      if ( length(grep(pattern="trans"     ,x=vname)) > 0 ||
           length(grep(pattern="full.tag"  ,x=vname)) > 0 ||
           length(grep(pattern="old."      ,x=vname)) > 0 ||
           length(grep(pattern=".1st"      ,x=vname)) > 0 ||
           length(grep(pattern="scientific",x=vname)) > 0 ||
           length(grep(pattern="genus"     ,x=vname)) > 0 ||
           length(grep(pattern="family"    ,x=vname)) > 0 ||
           length(grep(pattern="conflict." ,x=vname)) > 0){
         #---- We never copy these variables. ---------------------------------------------#
         this.val = this.val
         #---------------------------------------------------------------------------------#




      }else if (vname %in% c("year.last")){
         #---- We always copy these variables. --------------------------------------------#
         receptor[[vname]][don.2.rec] = this.val[don.2.rec]
         #---------------------------------------------------------------------------------#




      }else if (vname %in% c("year.added","year.recruit")){
         #---- We always keep the lowest value. -------------------------------------------#
         receptor[[vname]][don.2.rec] = pmin( this.val[don.2.rec]
                                            , receptor[[vname]][don.2.rec]
                                            , na.rm = TRUE )
         #---------------------------------------------------------------------------------#




      }else if (vname %in% c("year.death")){
         #---- We always keep the highest value. ------------------------------------------#
         receptor[[vname]][don.2.rec] = pmax( this.val[don.2.rec]
                                            , receptor[[vname]][don.2.rec]
                                            , na.rm = TRUE )
         #---------------------------------------------------------------------------------#


      }else if (vname %in% c("tag","x","y")){
         #---- Variables that may change and we update using later censuses. --------------#
         sel.1st                      = ( ( ! is.na(this.val)) 
                                        &     is.na(receptor[[first]]) )
         old.sel                      = ( (! is.na(this.val))
                                        & (! is.na(receptor[[vname]])) )
         receptor[[first]][sel.1st]   = this.val[sel.1st]
         receptor[[old  ]][old.sel]   = receptor[[vname]][old.sel]
         receptor[[vname]][don.2.rec] = this.val[don.2.rec]
         #---------------------------------------------------------------------------------#


      }else if (vname %in% c("common")){
         #---------------------------------------------------------------------------------#
         #     We must check whether the donor or receptor have scientific name and        #
         # family.  We always keep the one that has scientific name.  If both have         #
         # scientific name and they don't match, we stop.                                  #
         #---------------------------------------------------------------------------------#
         this.c.ignotum     = this.val            %in% c(unknown.common,NA_character_)



         #---------------------------------------------------------------------------------#
         #  1. Receptor has scientific name, donor doesn't; we copy the info from donor to #
         #     old.common but don't update common and scientific.                          #
         #---------------------------------------------------------------------------------#
         rec.sci                      = ( ! receptor.s.ignotum) & this.s.ignotum
         old.sel                      = ( rec.sci & (! this.c.ignotum)
                                        & receptor$common != this.val  )
         receptor$old.common[old.sel] = this.val[old.sel]
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #  2. Donor has scientific name, receptor doesn't; we copy common name (current   #
         #     and first), scientific name, and family info from donor to receptor, and    #
         #     move receptor common name to old.common.                                    #
         #---------------------------------------------------------------------------------#
         don.sci                      = receptor.s.ignotum & (! this.s.ignotum)
         old.sel                      = don.sci & (! receptor.c.ignotum)
         receptor$old.common[old.sel] = receptor$common[old.sel]
         receptor$common    [don.sci] = this.val       [don.sci]
         receptor$common.1st[don.sci] = this.1st       [don.sci]
         receptor$scientific[don.sci] = this.sci       [don.sci]
         receptor$genus     [don.sci] = this.sci       [don.sci]
         receptor$family    [don.sci] = this.fam       [don.sci]
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #  3. Neither donor nor receptor have scientific name.  Replace those that are NA #
         #     by the donor info.                                                          #
         #---------------------------------------------------------------------------------#
         no.sci  = receptor.s.ignotum    & this.s.ignotum
         sel.1st = no.sci & ( ! this.c.ignotum) & receptor.f.ignotum
         old.sel = no.sci & ( ! this.c.ignotum) & (! receptor.c.ignotum)
         new.sel = no.sci & ( ! this.c.ignotum)
         receptor[[first]][sel.1st] = this.val[sel.1st]
         receptor[[old  ]][old.sel] = receptor[[vname]][old.sel]
         receptor[[vname]][new.sel] = this.val[new.sel]
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #  4. Both receptor and donor have scientific name; and they are the same.  Do    #
         #     nothing, they are consistent already.                                       #
         #---------------------------------------------------------------------------------#
         same.sci = ( (! receptor.s.ignotum) & (! this.s.ignotum)
                    & receptor$scientific == this.sci )
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #  5. Receptor and donor have scientific name; and they are not the same.  This   #
         #     should not happen; crash!                                                   #
         #---------------------------------------------------------------------------------#
         diff.sci = ( (! receptor.s.ignotum) & (! this.s.ignotum)
                    & receptor$scientific != this.sci )
         if (any(diff.sci) && sci.strict){
            stop ("Mismatch between scientific names!")
         }else if (any(diff.sci)){
            update.scientific   = this.full | receptor.s.ignotum
            receptor$scientific = ifelse(update.scientific,this.sci,receptor$scientific)
            if (try.genus){
               receptor$genus   = ifelse(update.scientific,this.gen,receptor$genus     )
            }#end if
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Make sure common.1st is initialised.                                       #
         #---------------------------------------------------------------------------------#
         sel                      = (! receptor.c.ignotum ) & receptor.f.ignotum
         receptor$common.1st[sel] = receptor$common[sel]
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Flag conflicts.                                                            #
         #---------------------------------------------------------------------------------#
         mismatch = ( (! receptor.c.ignotum)
                    & (! receptor.o.ignotum)
                    & receptor[[vname]] != receptor[[old]] )
         receptor$conflict.common[mismatch] = TRUE
         #---------------------------------------------------------------------------------#


      }else if ( length(grep(pattern="dbh."  ,x=vname)) > 0 ||
                 length(grep(pattern="gbh."  ,x=vname)) > 0 ||
                 ( length(grep(pattern="notes.",x=vname)) > 0 && 
                   ! vname %in% c("notes.orig","notes.qaqc") ) ){
         #---------------------------------------------------------------------------------#
         #      Variables that we never update unless this is the actual year.             #
         #---------------------------------------------------------------------------------#
         if (vname %in% new.input){
            #----- Copy only the cells that had no info before. ---------------------------#
            sel = (! is.na(this.val)) & (is.na(receptor[[vname]]))
            receptor[[vname]][sel] = this.val[sel]
            #------------------------------------------------------------------------------#

         }else if (length(grep(pattern="dbh.",x=vname)) > 0){
            #------------------------------------------------------------------------------#
            #    If this is DBH from previous surveys, we check whether they               #
            # match, and flag any mismatch in any given year.                              #
            #------------------------------------------------------------------------------#
            mismatch.na           = is.na(receptor[[vname]]) != is.na(this.val)
            mismatch.val          = receptor[[vname]] != this.val
            receptor$conflict.dbh = ( receptor$conflict.dbh
                                    | mismatch.na | mismatch.val )
            #------------------------------------------------------------------------------#
         }#end if
         #---------------------------------------------------------------------------------#


      }else if (  length(grep(pattern="dead."     ,x=vname)) > 0
               || length(grep(pattern="canopy."   ,x=vname)) > 0
               || length(grep(pattern="height."   ,x=vname)) > 0
               || length(grep(pattern="band."     ,x=vname)) > 0
               || length(grep(pattern="pom."      ,x=vname)) > 0
               || length(grep(pattern="ladder"    ,x=vname)) > 0
               || length(grep(pattern="cnpj"      ,x=vname)) > 0
               || length(grep(pattern="soil"      ,x=vname)) > 0
               || length(grep(pattern="code"      ,x=vname)) > 0
               ){
         #---- Variables that we update only if they have NA. -----------------------------#
         sel                      = ( is.na(receptor[[vname]])
                                    & ( ! is.na(this.val)) )
         receptor[[vname]][sel]   = this.val[sel]
         #---------------------------------------------------------------------------------#


      }else if ( vname == "notes.orig"){
         for (yr in survey.years){
            notes.year = paste("notes",yr,sep=".")
            #------------------------------------------------------------------------------#
            #     If original notes were NA, then the comments are purely for              #
            # this year.                                                                   #
            #------------------------------------------------------------------------------#
            goto.yr                         = is.na(receptor[[notes.year]])
            goto.yr[is.na(goto.yr)]         = FALSE
            receptor[[notes.year]][goto.yr] = this.val[goto.yr]
            #----- Find lines where comments were added. ----------------------------------#
            idx        = which( (! is.na(receptor$notes.orig))
                              & (! is.na(this.val           )) )
            added      = rep(FALSE,times=n.receptor)
            added[idx] = ( mapply( FUN      = regexpr
                                 , pattern  = receptor$notes.orig[idx]
                                 , text     = this.val         [idx]
                                 , MoreArgs = list(ignore.case=TRUE)
                                 , SIMPLIFY = TRUE) > 0
                         & receptor$notes.orig[idx] != this.val[idx] )
            added      = which(added)
            #----- Find lines where comments were different. ------------------------------#
            other      = rep(FALSE,times=n.receptor)
            other[idx] = mapply( FUN      = regexpr
                               , pattern  = receptor$notes.orig[idx]
                               , text     = this.val         [idx]
                               , MoreArgs = list(ignore.case=TRUE)
                               , SIMPLIFY = TRUE) == -1
            other      = which(other)
            #----- 3. Remove old comments from the new notes before we copy. --------------#
            bye.old = mapply( FUN      = sub
                            , pattern  = receptor$notes.orig
                            , x        = this.val
                            , MoreArgs = list(replacement="")
                            )#end mapply
            bye.old = paste("AAAA",bye.old,"EEEE",sep="")
            bye.old = sub(pattern="AAAA " ,replacement="" ,x=bye.old)
            bye.old = sub(pattern="AAAA"  ,replacement="" ,x=bye.old)
            bye.old = sub(pattern="AAAA; ",replacement="" ,x=bye.old)
            bye.old = sub(pattern="AAAA;" ,replacement="" ,x=bye.old)
            bye.old = sub(pattern="AAAA. ",replacement="" ,x=bye.old)
            bye.old = sub(pattern="AAAA." ,replacement="" ,x=bye.old)
            bye.old = sub(pattern=";EEEE" ,replacement="" ,x=bye.old)
            bye.old = sub(pattern="; EEEE",replacement="" ,x=bye.old)
            bye.old = sub(pattern="EEEE"  ,replacement="" ,x=bye.old)
            receptor[[notes.year]] [added] = bye.old [added]
            receptor[[notes.year]] [other] = this.val[other]
            receptor$conflict.notes[added] = TRUE
            receptor$conflict.notes[other] = TRUE
            #----- 4. Replace original notes. ---------------------------------------------#
            receptor$notes.orig[don.2.rec] = this.val[don.2.rec]
         }#end for
         #---------------------------------------------------------------------------------#

      }else if ( vname == "notes.qaqc"){

            #------------------------------------------------------------------------------#
            #      This variable is likely to go away, for the time being we just          #
            # concatenate the information.                                                 #
            #------------------------------------------------------------------------------#
            receptor[[vname]] = concatenate.message(receptor[[vname]],this.val)
            #------------------------------------------------------------------------------#

      }else{
         #----- Missing instructions, quit... ---------------------------------------------#
         cat ("     * VARIABLE: ",vname," has no instrucion!","\n")
         stop("Could not process the variable, sorry!")
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#
   return(receptor)
}#end fuse.trees
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     Function blend.trees                                                                 #
#                                                                                          #
#     This function blends one tree to the other.  This is slightly different from the     #
# fusion routine above because blending adds information but compares the periods for      #
# which information is to be added.  The other difference is that IT CANNOT BE USED IN     #
# VECTOR MODE: it works between one donor and one receptor only.                           #
#                                                                                          #
#     receptor        -- The data frame that will receive the data                         #
#     donor           -- The data frame that will send the data                            #
#     don.2.rec       -- The indices of the receptor that will receive the data.           #
#     update.yr.notes -- Should we use notes.orig to make year notes? (TRUE/FALSE).        #
#------------------------------------------------------------------------------------------#
blend.trees = function(receptor,donor,years.blend){

   #---------------------------------------------------------------------------------------#
   #       Get the list of all variables.                                                  #
   #---------------------------------------------------------------------------------------#
   variables   = names(receptor)
   n.variables = length(variables)
   n.receptor  = nrow(receptor)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      We will create a vector of DBH variables.                                        #
   #---------------------------------------------------------------------------------------#
   dbh.vars = NULL
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Copy new data to the receptor dataset, one variable at a time as they have       #
   # different rules.                                                                      #
   #---------------------------------------------------------------------------------------#
   for (nv in sequence(n.variables)){
      #------------------------------------------------------------------------------------#
      #       Grab the variable name, and make some auxiliary names.                       #
      #------------------------------------------------------------------------------------#
      vname       = variables[nv]
      old         = paste("old",vname,sep=".")
      first       = paste(vname,"1st",sep=".")
      blend.input = paste(substring(vname,1,nchar(vname)-5),years.blend,sep=".")
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Decide what to do based on the variable.                                      #
      #------------------------------------------------------------------------------------#
      if ( length(grep(pattern="trans"     ,x=vname)) > 0 ||
           length(grep(pattern="full.tag"  ,x=vname)) > 0 ||
           length(grep(pattern="old."      ,x=vname)) > 0 ||
           length(grep(pattern=".1st"      ,x=vname)) > 0 ||
           length(grep(pattern="scientific",x=vname)) > 0 ||
           length(grep(pattern="genus"     ,x=vname)) > 0 ||
           length(grep(pattern="family"    ,x=vname)) > 0 ||
           length(grep(pattern="year.last" ,x=vname)) > 0 ||
           length(grep(pattern="cnpj"      ,x=vname)) > 0 ||
           length(grep(pattern="conflict." ,x=vname)) > 0){
         #---------------------------------------------------------------------------------#
         #     We skip these variables.                                                    #
         #---------------------------------------------------------------------------------#
         dummy = NA
         #---------------------------------------------------------------------------------#



      }else if (vname %in% c("year.added")){
         #---- We always keep the lowest value. -------------------------------------------#
         receptor$year.added = min(receptor$year.added,min(years.blend))
         #---------------------------------------------------------------------------------#



      }else if (vname %in% c("year.recruit")){
         #---- We always keep the lowest value. -------------------------------------------#
         receptor$year.recruit = min(receptor$year.recruit,donor$year.recruit)
         #---------------------------------------------------------------------------------#



      }else if (vname %in% c("year.death")){
         #---- We always keep the highest value. ------------------------------------------#
         receptor$year.death = max(receptor$year.death,donor$year.death)
         #---------------------------------------------------------------------------------#


      }else if (vname %in% c("tag","x","y")){
         #---------------------------------------------------------------------------------#
         #     We will consider changing the values only if the donor has information.     #
         #---------------------------------------------------------------------------------#
         if (! is.na(donor[[vname]])){
            #------------------------------------------------------------------------------#
            #     We update coordinates only if the last year of the original dataset was  #
            # less than the years we are blending.                                         #
            #------------------------------------------------------------------------------#
            if (  (  ( max(years.blend)  >  receptor$year.last )
                  && ( receptor[[vname]] != donor[[vname]]     ) )
               || is.na(receptor[[vname]]) ){
               receptor[[old  ]] = receptor[[vname]]
               receptor[[vname]] = donor   [[vname]]
            }#end if
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     We don't change the "1st" variables unless the data we are blending is   #
            # older than the first year the receptor has data.                             #
            #------------------------------------------------------------------------------#
            if (  (  ( min(years.blend)  <  receptor$year.1st  )
                  && ( receptor[[vname]] != donor[[vname]]     ) )
               || is.na(receptor[[vname]]) ){
               receptor[[first]] = donor[[first]]
            }#end if
            #------------------------------------------------------------------------------#
         }#end if
         #---------------------------------------------------------------------------------#


      }else if (vname %in% c("common")){
         #---------------------------------------------------------------------------------#
         #     We must check whether the donor or receptor have scientific name and        #
         # family.  We always keep the one that has scientific name.  If both have         #
         # scientific name and they don't match, we stop.                                  #
         #---------------------------------------------------------------------------------#
         receptor.s.ignotum = receptor$scientific %in% c(unknown.scientific,NA_character_)
         receptor.c.ignotum = receptor$common     %in% c(unknown.common    ,NA_character_)
         receptor.o.ignotum = receptor$old.common %in% c(unknown.common    ,NA_character_)
         receptor.f.ignotum = receptor$common.1st %in% c(unknown.common    ,NA_character_)
         donor.s.ignotum    = donor$scientific    %in% c(unknown.scientific,NA_character_)
         donor.c.ignotum    = donor$common        %in% c(unknown.common    ,NA_character_)
         donor.f.ignotum    = donor$common.1st    %in% c(unknown.common    ,NA_character_)



         #---------------------------------------------------------------------------------#
         #  1. Receptor has scientific name, donor doesn't; we copy the info from donor to #
         #     old.common but don't update common and scientific.                          #
         #---------------------------------------------------------------------------------#
         if (  ( ! receptor.s.ignotum) && donor.s.ignotum
            && ( ! donor.c.ignotum   ) && ( receptor$common != donor$common) ){
            receptor$old.common = donor$common
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #  2. Donor has scientific name, receptor doesn't; we copy common name (current   #
         #     and first), scientific name, and family info from donor to receptor, and    #
         #     move receptor common name to old.common.                                    #
         #---------------------------------------------------------------------------------#
         if ( receptor.s.ignotum & (! donor.s.ignotum )){

            if ( ( ! receptor.c.ignotum) & receptor$common != donor$common ){
               receptor$old.common = receptor$common
            }#end if

            receptor$common     = donor$common
            receptor$common.1st = donor$common.1st
            receptor$scientific = donor$scientific
            receptor$family     = donor$family
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #  3. Neither donor nor receptor have scientific name.  Check the years and       #
         #     decide which common name has preference (the most recent one).              #
         #---------------------------------------------------------------------------------#
         if ( receptor.s.ignotum && donor.s.ignotum ){
            #------------------------------------------------------------------------------#
            #     Check whether we update the common name.                                 #
            #------------------------------------------------------------------------------#
            if ( receptor$year.last < max(years.blend) ){
               copy.common = ! donor.c.ignotum
               update.old  = ( copy.common && (! receptor.c.ignotum)
                                           && receptor$common != donor$common )
               if (update.old ) receptor$old.common = receptor$common
               if (copy.common) receptor$common     = donor$common
            }#end if
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Check whether the donor first name is older than the receptor one.       #
            #------------------------------------------------------------------------------#
            if ( receptor$year.1st > min(years.blend) ){
               if ( ! donor.f.ignotum ){
                  receptor$common.1st = donor$common.1st
               }#end if
            }#end if
            #------------------------------------------------------------------------------#
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #  4. Both receptor and donor have scientific name; and they are the same.  Do    #
         #     nothing, they are consistent already.                                       #
         #---------------------------------------------------------------------------------#
         same.sci = ( (! receptor.s.ignotum) && (! donor.s.ignotum)
                    & receptor$scientific == donor$scientific )
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #  5. Receptor and donor have scientific name; and they are not the same.  This   #
         #     should not happen; crash!                                                   #
         #---------------------------------------------------------------------------------#
         diff.sci = ( (! receptor.s.ignotum) & (! donor.s.ignotum)
                    & receptor$scientific != donor$scientific )
         if (diff.sci){
            cat ("Mismatch between scientific names!","\n")
            browser()
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Make sure common.1st is initialised.                                       #
         #---------------------------------------------------------------------------------#
         sel                      = (! receptor.c.ignotum ) & receptor.f.ignotum
         receptor$common.1st[sel] = receptor$common[sel]
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Flag conflicts.                                                            #
         #---------------------------------------------------------------------------------#
         if( (! receptor.c.ignotum) & (! receptor.o.ignotum)
           & receptor[[vname]] != receptor[[old]] ){
            receptor$conflict.common = TRUE
         }#end if
         #---------------------------------------------------------------------------------#


      }else if (    length(grep(pattern="dbh."  ,x=vname)) > 0
               ||   length(grep(pattern="dead." ,x=vname)) > 0
               || ( length(grep(pattern="notes.",x=vname)) > 0
                  &&  ! vname %in% c("notes.orig","notes.qaqc") ) ){
         #---------------------------------------------------------------------------------#
         #      Variables that we update only if this is the blending year.                #
         #---------------------------------------------------------------------------------#
         if (vname %in% blend.input && (! is.na(donor[[vname]]))){
            #----- Copy only the cells that had no info before. ---------------------------#
            receptor[[vname]] = donor[[vname]]
            #------------------------------------------------------------------------------#
         }#end if
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     For DBH, we also check whether they match in case of overlap, and flag any  #
         # mismatch for blending years.                                                    #
         #---------------------------------------------------------------------------------#
         if ( length(grep(pattern="dbh."  ,x=vname)) > 0){
            dbh.vars              = c(dbh.vars,vname)
            mismatch.na           = is.na(receptor[[vname]]) != is.na(donor[[vname]])
            mismatch.val          = receptor[[vname]] != donor[[vname]]
            receptor$conflict.dbh = receptor$conflict.dbh | mismatch.na | mismatch.val
            #------------------------------------------------------------------------------#
         }#end if
         #---------------------------------------------------------------------------------#


      }else if (  length(grep(pattern="canopy."   ,x=vname)) > 0
               || length(grep(pattern="height."   ,x=vname)) > 0
               || length(grep(pattern="band."     ,x=vname)) > 0
               || length(grep(pattern="pom."      ,x=vname)) > 0
               || length(grep(pattern="ladder"    ,x=vname)) > 0
               ){
         #---- Variables that we update only if they have NA. -----------------------------#
         sel                      = ( is.na(receptor[[vname]])
                                    & ( ! is.na(donor[[vname]])) )
         receptor[[vname]][sel]   = donor[[vname]][sel]
         #---------------------------------------------------------------------------------#

      }else if ( vname %in% c("notes.qaqc","notes.orig")){
            #------------------------------------------------------------------------------#
            #      We concatenate the information.                                         #
            #------------------------------------------------------------------------------#
            receptor[[vname]] = concatenate.message(receptor[[vname]],donor[[vname]])
            #------------------------------------------------------------------------------#

      }else{
         #----- Missing instructions, quit... ---------------------------------------------#
         cat ("     * VARIABLE: ",vname," has no instrucion!","\n")
         stop("Could not process the variable, sorry!")
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Update year.1st and year.last with the updated info.                              #
   #---------------------------------------------------------------------------------------#
   dbh.years = sort(as.numeric(substring(dbh.vars,5)))
   n.years   = length(dbh.years)
   receptor$year.1st  = NA
   receptor$year.last = NA
   for (y in sequence(n.years)){
      yr      = dbh.years[y]
      dbh.now = receptor[[paste("dbh",yr,sep=".")]]
      #------------------------------------------------------------------------------------#
      #      Check whether to update the first and last year with valid DBH.               #
      #------------------------------------------------------------------------------------#
      if (is.finite(dbh.now)){
         receptor$year.last = yr
         if (is.na(receptor$year.1st)) receptor$year.1st = yr
      }#end if
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#

   return(receptor)
}#end fuse.trees
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function partially removes information from the line.                           #
#------------------------------------------------------------------------------------------#
delete.tree.period <<- function(tree,period){
   #---- Find the names of the variables we will clean. -----------------------------------#
   dbh.names    = paste("dbh"   ,period,sep=".")
   notes.names  = paste("notes" ,period,sep=".")
   canopy.names = paste("canopy",period,sep=".")
   dead.names   = paste("dead"  ,period,sep=".")
   height.names = paste("height",period,sep=".")
   band.names   = paste("band"  ,period,sep=".")
   reset.names  = c(dbh.names,notes.names,canopy.names,dead.names,height.names,band.names)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the columns we will discard.                                                 #
   #---------------------------------------------------------------------------------------#
   reset        = which(names(tree) %in% reset.names)
   tree[reset]  = NA
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Update year.1st and year.last with the updated info.                              #
   #---------------------------------------------------------------------------------------#
   dbh.vars       = names(tree)[grep("dbh.",names(tree))]
   dbh.years      = sort(as.numeric(substring(dbh.vars,5)))
   n.years        = length(dbh.years)
   tree$year.1st  = NA
   tree$year.last = NA
   for (y in sequence(n.years)){
      yr      = dbh.years[y]
      dbh.now = tree[[paste("dbh",yr,sep=".")]]
      #------------------------------------------------------------------------------------#
      #      Check whether to update the first and last year with valid DBH.               #
      #------------------------------------------------------------------------------------#
      if (is.finite(dbh.now)){
         tree$year.last = yr
         if (is.na(tree$year.1st)) tree$year.1st = yr
      }#end if
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#

   return(tree)
}#end if
#==========================================================================================#
#==========================================================================================#
