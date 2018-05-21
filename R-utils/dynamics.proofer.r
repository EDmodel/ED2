#==========================================================================================#
#==========================================================================================#
#      Common constants to control the proofs.                                             #
#------------------------------------------------------------------------------------------#
low.wrong     <<- -1.0   # Growth rates below this threshold [cm/yr] are considered wrong
high.wrong    <<- +5.0   # Growth rates above this threshold [cm/yr] are considered wrong
low.weird     <<- 0.025  # Bottom threshold for non-suspicious quantiles
high.weird    <<- 0.975  # Top threshold for non-suspicious quantiles
low.kink      <<- 0.10   # Threshold for considering a suspicious kink (trough)
high.kink     <<- 0.90   # Threshold for considering a suspicious kink (ridge)
lngrowth.min  <<- 0.002  # Minimum growth rate (relative)
dbh.min.sub   <<- 10.    # Minimum dbh to be included in the subplot
dbh.min.trans <<- 35.    # Minimum dbh to be included anywhere
dbh.min.toler <<- 0.5    # Tolerance for minimum DBH (so we are not too strict)
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      This function goes through the notes and spot death comments, and make sure the     #
# death flag is consistent with the notes and with the DBH.  This must be done in 4 steps: #
# 1.  There are lots of false alarms (e.g. not dead, 1/2 dead, previously dead now alive,  #
#     etc.) that must be replaced by something else so we don't flag them as death.  The   #
#     current list is good for census up to 2011; you may need to add more cases if you    #
#     are going to run for post 2011.                                                      #
# 2.  Flag trees that have notes that suggest they died.  We set the flags to 1 for the    #
#     first survey they are flagged.                                                       #
# 3.  Look for trees that turned out not to be dead, or underwent fusion, junction, or     #
#     blending.  Remove the death flags for these trees.                                   #
# 4.  Once the flags are consistent, move DBH for the death year and subsequent years to   #
#     the notes.                                                                           #
#                                                                                          #
#     We never flag missing trees as dead; death must be certified in order to be flagged. #
# The last time a recruited stop reporting data is saved as 'year.goodbye', and that can   #
# be also used as death year if you want to assume permanently missing trees are also dead #
#------------------------------------------------------------------------------------------#
death.proofer <<- function( datum                # Forest inventory data set
                          , year4                # Inventory years
                          , use.flags  = FALSE   # Use dead flags? (TRUE|FALSE)
                          , use.notes  = TRUE    # Use notes? (TRUE|FALSE)
                          , dbh.stdead = FALSE   # Standing dead may have DBH? (TRUE|FALSE)
                          , dead.pref  = "dead"  # Prefix of dead flag variables
                          , dbh.pref   = "dbh"   # Prefix of DBH variables
                          , notes.pref = "notes" # Prefix of note variables
                          , vsep       = "."     # Separator between prefixes and year
                          ){

   #----- Find data and years. ------------------------------------------------------------#
   n.years = length(year4)
   n.datum = nrow(datum)
   next.year = c(year4[-1],Inf)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Initialise the death year.                                                        #
   #---------------------------------------------------------------------------------------#
   datum$year.death = rep(Inf,times=n.datum)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Initialise notes in case the column does not exist.                               #
   #---------------------------------------------------------------------------------------#
   for (y in sequence(n.years)){
      this.notes = paste(notes.pref,year4[y],sep=vsep)
      if (! (this.notes %in% names(datum))){
         datum[[this.notes]] = rep(NA_character_,times=n.datum)
      }#end if (! (this.notes %in% names(datum))){
      #------------------------------------------------------------------------------------#
   }#end for (y in sequence(n.years))
   #---------------------------------------------------------------------------------------#

   
   #---------------------------------------------------------------------------------------#
   #     If the user wants to use the flags, loop over them and check for dead trees.      #
   #---------------------------------------------------------------------------------------#
   if (use.flags){
      for (y in sequence(n.years)){
         this.dead             = paste(dead.pref,year4[y],sep=vsep)
         if (this.dead %in% names(datum)){
            sel                   = as.numeric(datum[[this.dead]]) %==% 1
            datum$year.death[sel] = pmin(datum$year.death[sel],year4[y])
         }#end if (this.dead %in% names(datum))
         #---------------------------------------------------------------------------------#
      }#end for (y in sequence(n.years))
      #------------------------------------------------------------------------------------#
   }#end if (use.flags)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     In case standing dead were not supposed to have DBH measurements, shift death     #
   # year of dead trees with DBH measurements to the following survey.                     #
   #---------------------------------------------------------------------------------------#
   if (! dbh.stdead){
      for (y in sequence(n.years)){
         this.dbh                = paste(dbh.pref ,year4[y],sep=vsep)
         alive                   = ( is.finite(datum[[this.dbh]]) 
                                   & datum$year.death %<=% year4[y]
                                   )#end alive
         datum$year.death[alive] = next.year[y]
         #---------------------------------------------------------------------------------#
      }#end for (y in sequence(n.years))
      #------------------------------------------------------------------------------------#
   }#end if (use.flags)
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     If the user wants to use notes, then we first look for misleading messages        #
   # throughout the years; and replace words so dead means dead and alive means alive.     #
   # Then we assign the death year.                                                        #
   #---------------------------------------------------------------------------------------#
   if (use.notes){
      for (y in sequence(n.years)){
         this.notes = paste(notes.pref,year4[y],sep=vsep)
         x          = datum[[this.notes]]

         #----- Standardise life and death comments. --------------------------------------#
         x = sub( pattern     = "not dead"
                , replacement = "alive"
                , x           = x )
         x = sub( pattern     = "looks almost dead"
                , replacement = "dying"
                , x           = x )
         x = sub( pattern     = "almost dead"
                , replacement = "dying"
                , x           = x )
         x = sub( pattern     = "nearly dead"
                , replacement = "dying"
                , x           = x )
         x = sub( pattern     = "thought to be dead"
                , replacement = "thought to have perished"
                , x           = x )
         x = sub( pattern     = "thought as dead"
                , replacement = "thought to have perished"
                , x           = x )
         x = sub( pattern     = "half-dead"
                , replacement = "dying"
                , x           = x )
         x = sub( pattern     = "not all dead"
                , replacement = "still alive"
                , x           = x )
         x = sub( pattern     = "looks dead"
                , replacement = "dying"
                , x           = x )
         x = sub( pattern     = "1/2 dead"
                , replacement = "dying"
                , x           = x )
         x = sub( pattern     = "half dead"
                , replacement = "dying"
                , x           = x )
         x = sub( pattern     = "dead in 2003 but actually alive"
                , replacement = "alive"
                , x           = x )
         x = sub( pattern     = "not found -- assume dead"
                , replacement = "missing"
                , x           = x )
         x = sub( pattern     = "not  -- assume dead"
                , replacement = "missing"
                , x           = x )
         x = sub( pattern     = "not found\\.assume dead"
                , replacement = "missing"
                , x           = x )
         x = sub( pattern     = "couldn\\.t find??? dead??"
                , replacement = "missing"
                , x           = x )
         x = sub( pattern     = "couldn\\.t find; probably dead; in a large tree fall\\."
                , replacement = "missing; area has tree fall."
                , x           = x )
         x = sub( pattern     = "partially dead"
                , replacement = "dying"
                , x           = x )
         x = sub( pattern     = "99% dead"
                , replacement = "dying"
                , x           = x )
         x = sub( pattern     = "dead; it\\.s alive now"
                , replacement = "alive"
                , x           = x )
         x = sub( pattern     = "dead; alive now"
                , replacement = "alive"
                , x           = x )
         x = sub( pattern     = "half of tree it\\.s dead"
                , replacement = "dying"
                , x           = x )
         x = sub( pattern     = "not found\\.  assume dead"
                , replacement = "missing"
                , x           = x )
         x = sub( pattern     = "half of the tree is dead"
                , replacement = "dying"
                , x           = x )
         x = sub( pattern     = "dead; on the ground in a large tree fall; alive; broken"
                , replacement = "on the ground in a large tree fall; alive; broken"
                , x           = x )
         x = sub( pattern     = "dead; broken; alive"
                , replacement = "broken; alive"
                , x           = x )
         x = sub( pattern     = "deading"
                , replacement = "dying"
                , x           = x )
         x = sub( pattern     = "probably it is dead"
                , replacement = "likely dead"
                , x           = x )
         x = sub( pattern     = "one of the bifurcations is dead"
                , replacement = "one of the bifurcations perished"
                , x           = x )
         x = sub( pattern     = "not found; presumed dead"
                , replacement = "missing"
                , x           = x )
         x = sub( pattern     = "not found; probably dead"
                , replacement = "missing"
                , x           = x )
         x = sub( pattern     = "dead? or not measured?"
                , replacement = "missing"
                , x           = x )
         x = sub( pattern     = "alive\\. previously reported as dead"
                , replacement = "alive"
                , x           = x )
         x = sub( pattern     = "grew inside of other dead tree"
                , replacement = "grew inside of tree that had perished"
                , x           = x )
         x = sub( pattern     = "meia morto"
                , replacement = "dying"
                , x           = x )
         x = sub( pattern     = "nao foi encontrado\\.\\.\\. perdido ou morto?"
                , replacement = "missing"
                , x           = x )
         x = sub( pattern     = "nao morto"
                , replacement = "alive"
                , x           = x )
         x = sub( pattern     = "parece como um arvore caido e morto mais nao sei com certeza"
                , replacement = "it may be a fallen and dead tree but I am not sure"
                , x           = x )
         x = sub( pattern     = "decom 1 - morta?"
                , replacement = "decomp 1 - dead?"
                , x           = x )
         x = sub( pattern     = "nao foi encontrada; provavalmente caiu; tem caida grande aqui"
                , replacement = "missing; likely fallen; there was a large treefall there"
                , x           = x )
         x = sub( pattern     = "morta 2010"
                , replacement = "dead in 2010"
                , x           = x )
         x = sub( pattern     = "quase morta"
                , replacement = "dying"
                , x           = x )
         x = sub( pattern     = "morta? desaparecida"
                , replacement = "missing"
                , x           = x )
         x = sub( pattern     = "morta?"
                , replacement = "dead?"
                , x           = x )
         x = sub( pattern     = "living"
                , replacement = "alive"
                , x           = x )
         x = sub( pattern     = "vivo"
                , replacement = "alive"
                , x           = x )
         x = sub( pattern     = "viva"
                , replacement = "alive"
                , x           = x )
         x = sub( pattern     = "morto"
                , replacement = "dead"
                , x           = x )
         x = sub( pattern     = "morta"
                , replacement = "dead"
                , x           = x )
         #---------------------------------------------------------------------------------#

         datum[[this.notes]] = x
      }#end for
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Now we loop through all years and find out which year is the most appropriate  #
      # for death.  Because of the possible misidentifications, we keep updating it every  #
      # year.                                                                              #
      #------------------------------------------------------------------------------------#
      for (y in sequence(n.years)){
         this.notes = paste(notes.pref,year4[y],sep=vsep)
         x          = datum[[this.notes]]

         #----- Look for trees that were alive but died this year. ------------------------#
         dead       = grep("dead",x)
         prev.alive = which(datum$year.death > year4[y])
         new.dead   = intersect(dead,prev.alive)
         if (length(new.dead) > 0) datum$year.death[new.dead] = year4[y]
         #---------------------------------------------------------------------------------#



         #----- Look for trees that were flagged as dead but are actually alive. ----------#
         prev.dead  = which(datum$year.death <= year4[y])
         alive      = grep("alive",x)
         new.alive  = intersect(alive,prev.dead)
         if (length(new.alive) > 0) datum$year.death[new.alive] = Inf
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Look for trees that were flagged as dead but have undergone fusion at this  #
         # year.                                                                           #
         #---------------------------------------------------------------------------------#
         action       = unique(c(grep("Fusion"  ,datum$notes.qaqc)
                                ,grep("Junction",datum$notes.qaqc)
                                ,grep("Blend"   ,datum$notes.qaqc)))
         post.death   = which(datum$year.last > year4[y] & datum$year.death <= year4[y])
         concatenated = intersect(action,post.death)
         if (length(concatenated) > 0) datum$year.death[concatenated] = Inf
         #---------------------------------------------------------------------------------#
      }#end for
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Now we loop over all trees, correct the death flags, and move the DBH of dead     #
   # trees to the notes.                                                                   #
   #---------------------------------------------------------------------------------------#
   for (y in sequence(n.years)){
      #------------------------------------------------------------------------------------#
      #      Retrieve information from 
      #------------------------------------------------------------------------------------#
      this.dbh   = paste(dbh.pref,year4[y],sep=vsep)
      this.notes = paste(notes.pref,year4[y],sep=vsep)
      message    = ifelse( test = is.na(datum[[this.dbh]])
                         , yes  = NA_character_
                         , no   = paste0(dbh.pref,vsep,year4[y],"=",datum[[this.dbh]])
                         )#end ifelse
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Find the trees we must modify stuff.                                           #
      #------------------------------------------------------------------------------------#
      dead.new    = which( datum$year.death == year4[y] )
      dead.old    = which( datum$year.death <  year4[y] )
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Update the death flag and notes for trees that have recently died.             #
      #------------------------------------------------------------------------------------#
      if (length(dead.new) > 0 && (! dbh.stdead)){
         datum[[this.notes]][dead.new] = concatenate.message( datum[[this.notes]][dead.new]
                                                            , message            [dead.new]
                                                            )#end concatenate.message
         datum[[this.dbh  ]][dead.new] = NA
      }#end if
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Update the death flag and notes for trees that have died before but may still  #
      # have DBH data.                                                                     #
      #------------------------------------------------------------------------------------#
      if (length(dead.old) > 0 && (! dbh.stdead)){
         datum[[this.notes]][dead.old] = concatenate.message( datum[[this.notes]][dead.old]
                                                            , message            [dead.old]
                                                            )#end concatenate.message
         datum[[this.dbh  ]][dead.old] = NA
      }#end if
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #      Update the first and last year information.                                      #
   #---------------------------------------------------------------------------------------#
   datum$year.1st  = rep(-Inf,times=n.datum)
   datum$year.last = rep(-Inf,times=n.datum)
   for (y in sequence(n.years)){
      yr                       = year4[y]
      dbh.now                  = datum[[paste(dbh.pref,yr,sep=vsep)]]
      sel.dbh                  = is.finite(dbh.now)
      sel.living               = datum$year.death > yr
      sel.1st                  = sel.dbh & sel.living & (! is.finite(datum$year.1st))
      sel.last                 = sel.dbh & sel.living
      datum$year.1st [sel.1st] = yr
      datum$year.last[sel.dbh] = yr
   }#end for
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #      Set up the "goodbye" year (first year in which plant DBH is never reported       #
   # again.  In case DBH is reported for standing dead trees, we make sure that            #
   # year.goodbye is never after year.death.                                               #
   #---------------------------------------------------------------------------------------#
   datum$year.goodbye = rep(Inf,times=n.datum)
   for (y in sequence(n.years)){
      yr                             = year4[y]
      dbh.now                        = datum[[paste(dbh.pref,yr,sep=vsep)]]
      sel.dbh                        = is.finite(dbh.now)
      sel.bye                        = datum$year.1st < yr & is.na(dbh.now)
      datum$year.goodbye[sel.dbh]    = Inf
      sel.bye                        = ifelse(is.na(sel.bye),TRUE,sel.bye)
      if (any(sel.bye)){
         datum$year.goodbye[sel.bye] = pmin(datum$year.goodbye[sel.bye],yr)
      }#end if
      #------------------------------------------------------------------------------------#
   }#end for
   datum$year.goodbye = pmin(datum$year.goodbye,datum$year.death)
   #---------------------------------------------------------------------------------------#

   
   return(datum)
}#end function death.proofer
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      This function goes through the DBH data and spot data that are too strange to be    #
# true. We create a growth matrix and compare with the growth distribution for that year,  #
# and the growth history for that tree.  If both are weird, we move the datum to the notes #
# and temporarily replace by the negative number (to make it easier for recruitment).      #
#------------------------------------------------------------------------------------------#
growth.proofer <<- function( datum                   # Forest inventory data set
                           , month2                  # Months
                           , year4                   # Years
                           , indiv.dates = FALSE     # Use individual dates?
                           , miss.pref   = "missing" # Prefix of dead flag variables
                           , dbh.pref    = "dbh"     # Prefix of DBH variables
                           , date.pref   = "date"    # Prefix of date variables
                           , notes.pref  = "notes"   # Prefix of note variables
                           , vsep        = "."       # Separator between prefixes and year
                           ){


   #----- Define some auxiliary variables. ------------------------------------------------#
   n.years  = length(year4)
   n.datum  = nrow(datum)
   when     = year4 + (month2-0.5)/12
   #---------------------------------------------------------------------------------------#



   #----- Create a table with dbh, so it is easier to find growth rates. ------------------#
   mat.names    = list(datum$full.tag,year4)
   dbh.table    = matrix(NA   ,ncol=n.years,nrow=n.datum,dimnames=mat.names)
   date.table   = matrix(NA   ,ncol=n.years,nrow=n.datum,dimnames=mat.names)
   growth.bef   = matrix(NA   ,ncol=n.years,nrow=n.datum,dimnames=mat.names)
   growth.aft   = matrix(NA   ,ncol=n.years,nrow=n.datum,dimnames=mat.names)
   lngrowth.bef = matrix(NA   ,ncol=n.years,nrow=n.datum,dimnames=mat.names)
   lngrowth.aft = matrix(NA   ,ncol=n.years,nrow=n.datum,dimnames=mat.names)
   dtime.bef    = matrix(NA   ,ncol=n.years,nrow=n.datum,dimnames=mat.names)
   dtime.aft    = matrix(NA   ,ncol=n.years,nrow=n.datum,dimnames=mat.names)
   quantile.bef = matrix(NA   ,ncol=n.years,nrow=n.datum,dimnames=mat.names)
   quantile.aft = matrix(NA   ,ncol=n.years,nrow=n.datum,dimnames=mat.names)
   remove.data  = matrix(FALSE,ncol=n.years,nrow=n.datum,dimnames=mat.names)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Loop over years and fill in the DBH table.                                       # 
   #---------------------------------------------------------------------------------------#
   for (y in sequence(n.years)){
      this.dbh         = paste(dbh.pref  ,year4[y],sep=vsep)
      this.date        = paste(date.pref ,year4[y],sep=vsep)
      this.missing     = paste(miss.pref ,year4[y],sep=vsep)
      this.notes       = paste(notes.pref,year4[y],sep=vsep)
      dbh.table   [,y] = datum[[this.dbh]]
      if (indiv.dates){
         this.when      = as.numeric(datum[[this.date]])
         this.year      = numyears(datum[[this.date]])
         this.when0     = as.numeric(chron(paste(12,31,this.year-1,sep="/")))
         this.yr.day    = ifelse(test=is.leap(this.year),yes=366,no=365)
         date.table[,y] = this.year + (this.when - this.when0) / this.yr.day
      }else{
         date.table[,y] = when[y]
      }#end if

      if (this.missing %in% names(datum)){
         del              = as.numeric(datum[[this.missing]]) == 1
         dbh.table[del,y] = NA
         message          = paste0(this.dbh,"=",datum[[this.dbh]]
                                  ," was gap filled thus removed")
         datum[[this.notes]][del] = concatenate.message( datum[[this.notes]][del]
                                                       , message[del])
         datum[[this.dbh  ]][del] = NA
      }#end if (this.missing %in% names(datum))
      #------------------------------------------------------------------------------------#
   }#end for (y in sequence(n.years))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the growth rates between surveys and previous one with valid data, and       #
   # organise them by quantiles.                                                           #
   #---------------------------------------------------------------------------------------#
   for (y in sequence(n.years)[-1]){
      dbh.bef  = rep(NA,times=n.datum)
      when.bef = rep(NA,times=n.datum)
      for (b in sequence(y-1)){
         sel           = is.finite(dbh.table[,b])
         dbh.bef [sel] = dbh.table [sel,b]
         when.bef[sel] = date.table[sel,b]
      }#end for (b in seq(from=1,y-1,+1))
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Define the growth rates between the previous measurement and the current.      #
      #------------------------------------------------------------------------------------#
      dtime.bef   [,y] = date.table[,y]-when.bef
      lngrowth.bef[,y] = log(dbh.table[,y]/dbh.bef)/dtime.bef[,y]
      growth.bef  [,y] = (dbh.table[,y]-dbh.bef)/dtime.bef[,y]
      if (any(is.finite(lngrowth.bef[,y]))){
         efun             = try(ecdf(lngrowth.bef[,y]))
         quantile.bef[,y] = efun(lngrowth.bef[,y])
      }#end if (any(is.finite(lngrowth.bef[,y])))
      #------------------------------------------------------------------------------------#
   }#end for (y in sequence(n.years)[-1])
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the growth rates between surveys and following one with valid data, and      #
   # organise them by quantiles.                                                           #
   #---------------------------------------------------------------------------------------#
   for (y in sequence(n.years-1)){
      dbh.aft  = rep(NA,times=n.datum)
      when.aft = rep(NA,times=n.datum)
      for (b in seq(from=n.years,to=y+1,by=-1)){
         sel           = is.finite(dbh.table[,b])
         dbh.aft [sel] = dbh.table [sel,b]
         when.aft[sel] = date.table[sel,b]
      }#end for (b in seq(from=n.years,y+1,-1))
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Define the growth rates between the previous measurement and the current.      #
      #------------------------------------------------------------------------------------#
      dtime.aft   [,y] = when.aft-date.table[,y]
      lngrowth.aft[,y] = log(dbh.aft/dbh.table[,y])/dtime.aft[,y]
      growth.aft  [,y] = (dbh.aft-dbh.table[,y])/dtime.aft[,y]
      if (any(is.finite(lngrowth.aft[,y]))){
         efun             = ecdf(lngrowth.aft[,y])
         quantile.aft[,y] = efun(lngrowth.aft[,y])
      }#end if (any(is.finite(lngrowth.aft[,y])))
      #------------------------------------------------------------------------------------#
   }#end for (y in sequence(n.years-1))
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Find the "leapfrog" growth rates at any given year.                               #
   #---------------------------------------------------------------------------------------#
   lngrowth.leap = ( (lngrowth.bef/dtime.bef + lngrowth.aft/dtime.aft )
                   * (dtime.bef + dtime.aft) )
   growth.leap   = ( (growth.bef/dtime.bef + growth.aft/dtime.aft )
                   * (dtime.bef + dtime.aft) )
   quantile.leap = NA * lngrowth.leap
   for (y in sequence(n.years-1)[-1]){
      if (any(is.finite(lngrowth.leap[,y]))){
         efun              = ecdf(lngrowth.leap[,y])
         quantile.leap[,y] = efun(lngrowth.leap[,y])
      }#end if (any(lngrowth.leap[,y]))
      #------------------------------------------------------------------------------------#
   }#end for (y in sequence(n.years-1)[-1])
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Look for data that is unnaceptable in absolute numbers.  The threshold is some-   #
   # what arbitrary, but they can be easily adjusted.  We first look at the middle ones,   #
   # then we search at the edges.                                                          #
   #---------------------------------------------------------------------------------------#
   for (y in sequence(n.years)[-1]){

      #------------------------------------------------------------------------------------#
      #     Detect data that "shrunk" too much.  Then we decide which data point was bad,  #
      # the current dbh or the previous one.                                               #
      #------------------------------------------------------------------------------------#
      shrink  = ( ! remove.data[,max(1,y-1)] & ! remove.data[,min(n.years,y+1)]
                & is.finite(quantile.bef [,y])  & quantile.bef [,y] < low.weird
                & is.finite(growth.bef   [,y])  & growth.bef   [,y] < low.wrong )
      #----- Find out whether this is a kink. ---------------------------------------------#
      kink    = ( shrink
                & is.finite(quantile.aft [,y])  & quantile.aft [,y] > high.kink )
      #----- Find out whether the previous one is a bad data. -----------------------------#
      plateau = ( shrink & ! kink 
                & is.finite(quantile.leap[,y])  & quantile.leap[,y] < low.weird
                & is.finite(growth.leap  [,y])  & growth.leap  [,y] < low.wrong )
      #------------------------------------------------------------------------------------#


      #----- Flag data that is suspicious and should be removed. --------------------------#
      remove.data   [kink   ,  y] = TRUE
      remove.data   [plateau,y-1] = TRUE
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Detect trees that grew too much.  Then we decide which data point was bad,     #
      # the current dbh or the previous one.                                               #
      #------------------------------------------------------------------------------------#
      boom    = ( ! remove.data[,max(1,y-1)] & ! remove.data[,min(n.years,y+1)]
                & is.finite(quantile.bef [,y])  & quantile.bef [,y] > high.weird
                & is.finite(growth.bef   [,y])  & growth.bef   [,y] > high.wrong )
      #----- Find out whether this is a kink. ---------------------------------------------#
      kink    = ( boom
                & is.finite(quantile.aft [,y])  & quantile.aft [,y] < low.kink )
      #----- Find out whether the previous one is a bad data. -----------------------------#
      plateau = ( boom & ! kink 
                & is.finite(quantile.leap[,y])  & quantile.leap[,y] > high.weird
                & is.finite(growth.leap  [,y])  & growth.leap  [,y] > high.wrong )
      #------------------------------------------------------------------------------------#


      #----- Flag data that is suspicious and should be removed. --------------------------#
      remove.data   [kink   ,  y] = TRUE
      remove.data   [plateau,y-1] = TRUE
      #------------------------------------------------------------------------------------#
   }#end for (y in sequence(n.years)[-1])
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Similar to the first loop, but we look at the after derivatives so we can flag    #
   # the last census in case it is strange.                                                #
   #---------------------------------------------------------------------------------------#
   for (y in rev(sequence(n.years-1))){

      #------------------------------------------------------------------------------------#
      #     Detect data that "shrunk" too much.  Then we decide which data point was bad,  #
      # the current dbh or the previous one.                                               #
      #------------------------------------------------------------------------------------#
      shrink  = ( ! remove.data[,max(1,y-1)] & ! remove.data[,min(n.years,y+1)]
                & is.finite(quantile.aft [,y])  & quantile.aft [,y] < low.weird 
                & is.finite(growth.aft   [,y])  & growth.aft   [,y] < low.wrong )
      #----- Find out whether this is a kink. ---------------------------------------------#
      kink    = ( shrink
                & is.finite(quantile.bef [,y])  & quantile.bef [,y] > high.kink )
      #----- Find out whether the previous one is a bad data. -----------------------------#
      plateau = ( shrink & ! kink
                & is.finite(quantile.leap[,y])  & quantile.leap[,y] < low.weird
                & is.finite(growth.leap  [,y])  & growth.leap  [,y] < low.wrong )
      #------------------------------------------------------------------------------------#


      #----- Flag data that is suspicious and should be removed. --------------------------#
      remove.data   [kink   ,  y] = TRUE
      remove.data   [plateau,y+1] = TRUE
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Detect trees that grew too much.  Then we decide which data point was bad,     #
      # the current dbh or the previous one.                                               #
      #------------------------------------------------------------------------------------#
      boom    = ( ! remove.data[,max(1,y-1)] & ! remove.data[,min(n.years,y+1)]
                & is.finite(quantile.aft [,y])  & quantile.aft [,y] > high.weird 
                & is.finite(growth.aft   [,y])  & growth.aft   [,y] > high.wrong )
      #----- Find out whether this is a kink. ---------------------------------------------#
      kink    = ( boom
                & is.finite(quantile.bef [,y])  & quantile.bef [,y] < low.kink )
      #----- Find out whether the previous one is a bad data. -----------------------------#
      plateau = ( boom & ! kink
                & is.finite(quantile.leap[,y])  & quantile.leap[,y] > high.weird
                & is.finite(growth.leap  [,y])  &  growth.leap [,y] > high.wrong )
      #------------------------------------------------------------------------------------#


      #----- Flag data that is suspicious and should be removed. --------------------------#
      remove.data   [kink   ,  y] = TRUE
      remove.data   [plateau,y+1] = TRUE
      #------------------------------------------------------------------------------------#
   }#end for (y in rev(sequence(n.years-1)))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Loop over all years, discard the flagged data, and send the measurement to the   #
   # notes.                                                                                #
   #---------------------------------------------------------------------------------------#
   for (y in sequence(n.years)){
      this.dbh   = paste (dbh.pref  ,year4[y],sep=vsep)
      this.notes = paste (notes.pref,year4[y],sep=vsep)
      message    = paste0(this.dbh,"=",datum[[this.dbh]]
                         ," looks suspicious thus removed"
                         )#end paste0

      sel        = remove.data[,y]
      datum[[this.notes]][sel] = concatenate.message(datum[[this.notes]][sel],message[sel])

      #------------------------------------------------------------------------------------#
      #     Make the data to be removed temporarily negative, so they are finite but we    #
      # know that they are bad data that shall be purged.                                  #
      #------------------------------------------------------------------------------------#
      datum[[this.dbh]][sel] = -datum[[this.dbh]][sel]
      dbh.table[sel,y]       = -dbh.table[sel,y]
   }#end for (y in sequence(n.years))
   #---------------------------------------------------------------------------------------#

   return(datum)
}#end function growth.proofer
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      This function finds missing DBH data and fill them.  DBH is filled in 3 cases:      #
# 1.  Spurious DBH had been discarded (wrong tree, bad point of measurements, tight vines, #
#     measurement on buttress).                                                            #
# 2.  Tumbleweed effect: tree had bad coordinates and was not found in some surveys, but   #
#     it was reported before and after                                                     #
# 3.  Saint Exupery's effect: recruit appeared out of nowhere with a DBH that is too big   #
#     to be a true recruit.  We extrapolate the data based on growth rates after the tree  #
#     is established.                                                                      #
#------------------------------------------------------------------------------------------#
dbh.gap.filler <<- function( datum
                           , month2
                           , year4
                           , abs.y.sub     = Inf
                           , indiv.dates   = FALSE
                           , gf.individual = FALSE
                           , dbh.brks      = c(10,20,35,55,Inf)
                           , dbh.pref      = "dbh"
                           , gf.dbh.pref   = "gf.dbh"
                           , date.pref     = "date"
                           , notes.pref    = "notes"
                           , vsep          = "."
                           , dbhmin.small  = dbh.min.sub
                           , dbhmin.large  = dbh.min.trans
                           ){

   #----- Define some auxiliary variables. ------------------------------------------------#
   n.years      = length(year4)
   n.datum      = nrow(datum)
   when         = year4 + (month2-0.5)/12
   #---------------------------------------------------------------------------------------#
 
   #---------------------------------------------------------------------------------------#
   #     Loop over years and fill in the DBH table.                                        #
   #---------------------------------------------------------------------------------------#
   for (y in sequence(n.years)){
      this.gf.dbh      = paste(gf.dbh.pref,year4[y],sep=vsep)
      if (! this.gf.dbh %in% names(datum)){
         datum[[this.gf.dbh]] = rep(0,times=n.datum)
      }#end if (! this.gf.dbh %in% names(datum))
      #------------------------------------------------------------------------------------#
   }#end for (y in sequence(n.years))
   #---------------------------------------------------------------------------------------#



   #----- Create a table with dbh, so it is easier to find growth rates. ------------------#
   mat.names               = list(datum$full.tag,year4)
   dbh.table               = matrix(NA   ,ncol=n.years,nrow=n.datum,dimnames=mat.names)
   date.table              = matrix(NA   ,ncol=n.years,nrow=n.datum,dimnames=mat.names)
   lngrowth                = matrix(NA   ,ncol=n.years,nrow=n.datum,dimnames=mat.names)
   dtime                   = matrix(NA   ,ncol=n.years,nrow=n.datum,dimnames=mat.names)
   if (! ("subplot" %in% names(datum))){
      datum$subplot = (abs(datum$y) <= abs.y.sub)
   }#end if (! ("subsample" %in% names(datum)))
   dbh.min = ifelse(test = datum$subplot, yes=dbhmin.small,no=dbhmin.large)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Loop over years and fill in the DBH table.                                       # 
   #---------------------------------------------------------------------------------------#
   y.1st   = rep(NA,times=n.datum)
   dbh.1st = rep(NA,times=n.datum)
   for (y in sequence(n.years)){
      this.dbh         = paste(dbh.pref ,year4[y],sep=vsep)
      this.date        = paste(date.pref,year4[y],sep=vsep)
      dbh.ok           = is.finite(datum[[this.dbh]]) & datum$year.goodbye > year4[y]

      #------------------------------------------------------------------------------------#
      #     Some trees don't have coordinates, fill in the minimum DBH based on the first  #
      # instance.                                                                          #
      #------------------------------------------------------------------------------------#
      fill.min          = is.na(dbh.min) & dbh.ok & datum[[this.dbh]] >= dbhmin.small
      dbh.min[fill.min] = dbhmin.small
      fill.min          = is.na(dbh.min) & dbh.ok & datum[[this.dbh]] >= dbhmin.large
      dbh.min[fill.min] = dbhmin.large
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Check whether this is a valid data set for the neat table.                     #
      #------------------------------------------------------------------------------------#
      sel              = ( dbh.ok 
                         & datum[[this.dbh]] > 0
                         & ( datum[[this.dbh]] >= dbh.min
                           | ( ! is.na(y.1st) & y.1st <= y) )
                         )
      dbh.table[sel,y] = datum[[this.dbh]][sel]
      #------------------------------------------------------------------------------------#

      #----- Fill in date table. ----------------------------------------------------------#
      if (indiv.dates){
         this.when      = datum[[this.date]]
         this.year      = numyears(this.when)
         this.when0     = chron(paste(12,31,this.year-1,sep="/"))
         this.yr.day    = ifelse(test=is.leap(this.year),yes=366,no=365)
         date.table[,y] = this.year + as.numeric((this.when - this.when0)) / this.yr.day
      }else{
         date.table[,y] = when[y]
      }#end if
      #------------------------------------------------------------------------------------#

      #----- Check whether this is the first time this tree has good data. ----------------#
      sel              = sel & is.na(y.1st)
      y.1st    [sel  ] = y
      #------------------------------------------------------------------------------------#
   }#end for (y in sequence(n.years))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #       Get the mean DBH for each tree, then find to which size class it belongs.       #
   #---------------------------------------------------------------------------------------#
   dbh.mean = rowMeans(dbh.table,na.rm=TRUE)
   dbh.size = as.integer(cut(x=dbh.mean,breaks=dbh.brks,right=FALSE))
   dbh.size = ifelse(is.finite(dbh.size),dbh.size,max(dbh.size,na.rm=TRUE)+1)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Find the growth rates between surveys and previous one with valid data, and       #
   # organise them by quantiles.                                                           #
   #---------------------------------------------------------------------------------------#
   for (y in sequence(n.years)[-1]){
      dbh.bef  = rep(NA,times=n.datum)
      when.bef = rep(NA,times=n.datum)
      for (b in seq(from=1,y-1,+1)){
         sel           = is.finite(dbh.table[,b]) & dbh.table[,b] > 0
         dbh.bef [sel] = dbh.table [sel,b]
         when.bef[sel] = date.table[sel,b]
      }#end for
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Define the growth rates between the previous measurement and the current.      #
      #------------------------------------------------------------------------------------#
      dtime    [,y] = date.table[,y]-when.bef
      lngrowth [,y] = log(dbh.table[,y] / dbh.bef) / dtime[,y]
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the mean growth rates for each tree.  First we find the median rates for     #
   # each genus.                                                                           #
   #---------------------------------------------------------------------------------------#
   mean.tree.lngrowth    = apply(X=lngrowth,MARGIN=1,FUN=mean,na.rm=TRUE)
   mean.genus.size       = tapply( X     = mean.tree.lngrowth
                                 , INDEX = list(datum$genus,dbh.size)
                                 , FUN   = mean
                                 , na.rm = TRUE
                                 )#end tapply
   mean.genus            = matrix( data     = tapply( X     = mean.tree.lngrowth
                                                    , INDEX = datum$genus
                                                    , FUN   = mean
                                                    , na.rm = TRUE
                                                    )#end tapply
                                 , ncol     = ncol(mean.genus.size)
                                 , nrow     = nrow(mean.genus.size)
                                 , dimnames = dimnames(mean.genus.size)
                                 , byrow    = FALSE
                                 )#end matrix
   mean.size             = matrix( data     = tapply( X     = mean.tree.lngrowth
                                                    , INDEX = dbh.size
                                                    , FUN   = mean
                                                    , na.rm = TRUE
                                                    )#end tapply
                                 , ncol     = ncol(mean.genus.size)
                                 , nrow     = nrow(mean.genus.size)
                                 , dimnames = dimnames(mean.genus.size)
                                 , byrow    = TRUE
                                 )#end matrix
   mean.genus.size       = ifelse( is.finite(mean.genus.size)
                                 , mean.genus.size
                                 , ifelse(is.finite(mean.size),mean.size,mean.genus)
                                 )#end ifelse
   mean.genus.size       = ifelse(mean.genus.size<lngrowth.min,lngrowth.min,mean.genus.size)
   idx.genus             = match(datum$genus,dimnames(mean.genus.size)[[1]])
   idx.size              = match(dbh.size   ,dimnames(mean.genus.size)[[2]])
   typical.tree.lngrowth = mean.genus.size[cbind(idx.genus,idx.size)]
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     If the user prefer to fill using the individual typical growth rate, overwrite    #
   # the data (unless the tree had only one valid measurement.                             #
   #---------------------------------------------------------------------------------------#
   if (gf.individual){
      typical.tree.growth = ifelse( is.finite(mean.tree.lngrowth)
                                  , mean.tree.lngrowth
                                  , typical.tree.lngrowth
                                  )#end ifelse
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     We loop again and detect the years that shall be filled due to bad measurements   #
   # or tumbleweed effect.                                                                 #
   #---------------------------------------------------------------------------------------#
   datum$year.recruit = rep(Inf  , times=n.datum)
   removed            = rep(FALSE,times=n.datum)
   for (y in sequence(n.years)){
      dbh.label    = paste(dbh.pref   ,year4[y],sep=vsep)
      gf.dbh.label = paste(gf.dbh.pref,year4[y],sep=vsep)
      notes.label  = paste(notes.pref ,year4[y],sep=vsep)



      #----- Find out whether to fill the data. -------------------------------------------#
      measured                     = ( ! is.na(dbh.table[,y]))
      if (any(measured)){
         datum$year.recruit[measured] = pmin(year4[y],datum$year.recruit[measured])
      }#end if
      #------------------------------------------------------------------------------------#



      #----- Find out whether to fill the data. -------------------------------------------#
      dbh.miss      = is.na(datum[[dbh.label]]) | datum[[dbh.label]] %<% 0
      if (any(is.na(dbh.miss))){
         cat0(" dbh.miss has NA!")
      }#end if
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Trash information from years before the recruitment that had data.  This had   #
      # been taken care of before, but after discarding some bad measurements a few other  #
      # pre-recruitment points may appear.                                                 #
      #------------------------------------------------------------------------------------#
      bye     = datum[[dbh.label]] %>% 0 & is.na(dbh.table[,y])
      message = paste0(dbh.label,"=",abs(datum[[dbh.label]])
                      ," is pre-recruitment thus removed",sep="")
      datum[[notes.label]][bye] = concatenate.message(datum[[notes.label]][bye]
                                                     ,message[bye])
      datum[[dbh.label  ]][bye] = NA
      datum$year.recruit  [bye] = Inf
      removed = removed | bye
      #---------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Find the previous DBH.                                                         #
      #------------------------------------------------------------------------------------#
      dbh.bef  = rep(NA,times=n.datum)
      when.bef = rep(NA,times=n.datum)
      if (y != 1){
         for (b in seq(from=1,y-1,+1)){
            sel           = is.finite(dbh.table[,b])
            dbh.bef [sel] = dbh.table [sel,b]
            when.bef[sel] = date.table[sel,b]
         }#end for
      }#end if
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Find the next DBH.                                                             #
      #------------------------------------------------------------------------------------#
      dbh.aft  = rep(NA,times=n.datum)
      when.aft = rep(NA,times=n.datum)
      if (y != n.years){
         for (b in seq(from=n.years,y+1,-1)){
            sel           = is.finite(dbh.table[,b])
            dbh.aft [sel] = dbh.table [sel,b]
            when.aft[sel] = date.table[sel,b]
         }#end for
      }#end if
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Estimate the DBH by using the median growth rate starting from before and      #
      # after.  If both are available, then use the mean value, otherwise, take the one    #
      # that is available.                                                                 #
      #------------------------------------------------------------------------------------#
      dbh.guess = data.frame( before = dbh.bef
                                     * exp( typical.tree.lngrowth 
                                          * (date.table[,y] - when.bef) )
                            , after  = dbh.aft
                                     * exp( typical.tree.lngrowth
                                          * (date.table[,y] - when.aft))
                            )#end data.frame
      gfflg     = 1. + rowSums(is.na(dbh.guess),na.rm=TRUE)
      dbh.guess = round(rowMeans(dbh.guess,na.rm=TRUE),1)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Choose the years to update.  It has to be NA and it has to be post-            #
      # recruitment to be accepted.                                                        #
      #------------------------------------------------------------------------------------#
      is.recruited = ( dbh.guess %>% ( dbh.min + dbh.min.toler )
                     | datum[[dbh.label]] %>=% dbh.min
                     | datum$year.recruit <= year4[y] )
      is.alive     = datum$year.death > year4[y]
      update.year  = dbh.miss & is.recruited & is.alive
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Fill in the dbh and gap filling flag.                                          #
      #------------------------------------------------------------------------------------#
      if (any(update.year,na.rm=TRUE)){
         datum[[dbh.label   ]][update.year] = dbh.guess[update.year]
         datum[[gf.dbh.label]][update.year] = gfflg    [update.year]
         datum$year.recruit   [update.year] = pmin(datum$year.recruit[update.year]
                                                  ,year4[y]
                                                  )#end pmin
      }#end if
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#

   datum$year.recruit[datum$year.recruit == year4[1]] = -Inf
   return(datum)
}#end function dbh.gap.filler
#==========================================================================================#
#==========================================================================================#
