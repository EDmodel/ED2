#==========================================================================================#
#==========================================================================================#
#     This function is supposed to work for Santarem census only, feel free to adapt it    #
# elsewhere.  It creates a collection of unique tags based on various variables such as    #
# tags, coordinates, common name, and the older values that these variables once had.      #
#------------------------------------------------------------------------------------------#
census.tagger.s67 <<- function(merged,this,survey.years){
   #=======================================================================================#
   #=======================================================================================#
   #      Create all possible tags.  We must declare them as global variables because they #
   # must be accessed from sapply.                                                         #
   #---------------------------------------------------------------------------------------#


   #----- Transects. ----------------------------------------------------------------------#
   n.trans.merged                      <<- merged$trans
   n.trans.this                        <<- this$trans
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Scientific name flag, so it never fuses 2 trees with non-NA scientific names     #
   # unless they are the same.                                                             #
   #---------------------------------------------------------------------------------------#
   merged.ignotum = merged$scientific %in% c(unknown.scientific,NA_character_)
   this.ignotum   = this$scientific   %in% c(unknown.scientific,NA_character_)
   n.scientific.merged                 <<- rep("M" ,times=nrow(merged))
   n.scientific.this                   <<- rep("T" ,times=nrow(this)  )
   n.scientific.merged[merged.ignotum] <<- "T"
   n.scientific.this  [this.ignotum  ] <<- "M"
   o.scientific.merged                 <<- rep("M" ,times=nrow(merged))
   o.scientific.this                   <<- rep("T" ,times=nrow(this)  )
   o.scientific.merged[merged.ignotum] <<- "O"
   o.scientific.this  [this.ignotum  ] <<- "O"
   #---------------------------------------------------------------------------------------#


   #----- Coordinates. --------------------------------------------------------------------#
   n.x.merged                     = as.integer(10 * merged$x)
   o.x.merged                     = as.integer(10 * merged$old.x)
   f.x.merged                     = as.integer(10 * merged$x.1st)
   n.x.this                       = as.integer(10 * this$x)
   s.x.this                       = as.integer(10 * (1000 - this$x))
   n.y.merged                     = as.integer(10 * merged$y)
   o.y.merged                     = as.integer(10 * merged$old.y)
   f.y.merged                     = as.integer(10 * merged$y.1st)
   n.y.this                       = as.integer(10 * this$y)
   s.y.this                       = as.integer(10 * (- this$y))
   f.x.merged [is.na(f.x.merged)] = n.x.merged[is.na(f.x.merged)     ]
   f.y.merged [is.na(f.y.merged)] = n.y.merged[is.na(f.y.merged)     ]
   o.x.merged [is.na(o.x.merged)] = f.x.merged[is.na(o.x.merged)     ]
   o.y.merged [is.na(o.y.merged)] = f.y.merged[is.na(o.y.merged)     ]
   n.xy.merged <<- paste0("x",n.x.merged,"+y",n.y.merged)
   o.xy.merged <<- paste0("x",o.x.merged,"+y",o.y.merged)
   f.xy.merged <<- paste0("x",f.x.merged,"+y",f.y.merged)
   n.xy.this   <<- paste0("x",n.x.this  ,"+y",n.y.this  )
   s.xy.this   <<- paste0("x",s.x.this  ,"+y",n.y.this  )
   c.xy.this   <<- paste0("x",n.x.this  ,"+y",s.y.this  )
   z.xy.this   <<- paste0("x",s.x.this  ,"+y",s.y.this  )
   #---------------------------------------------------------------------------------------#




   #------ Tags. --------------------------------------------------------------------------#
   n.tag.merged                      <<- merged$tag
   o.tag.merged                      <<- merged$tag.1st
   f.tag.merged                      <<- merged$tag.1st
   n.tag.this                        <<- this$tag
   o.tag.this                        <<- this$old.tag
   f.tag.merged[is.na(f.tag.merged)] <<- n.tag.merged   [is.na(f.tag.merged)]
   o.tag.merged[is.na(o.tag.merged)] <<- f.tag.merged   [is.na(o.tag.merged)]
   o.tag.this  [is.na(o.tag.this)  ] <<- n.tag.this     [is.na(o.tag.this)  ]
   #---------------------------------------------------------------------------------------#




   #------ Common name. -------------------------------------------------------------------#
   o.merged.ignotum = merged$old.common %in% c(unknown.common,NA_character_)
   f.merged.ignotum = merged$common.1st %in% c(unknown.common,NA_character_)
   o.this.ignotum   = this$old.common   %in% c(unknown.common,NA_character_)
   f.this.ignotum   = this$common.1st   %in% c(unknown.common,NA_character_)
   n.common.merged                   <<- merged$common
   o.common.merged                   <<- merged$old.common
   f.common.merged                   <<- merged$common.1st
   n.common.this                     <<- this$common
   f.common.merged[f.merged.ignotum] <<- n.common.merged[f.merged.ignotum]
   o.common.merged[o.merged.ignotum] <<- f.common.merged[o.merged.ignotum]
   #---------------------------------------------------------------------------------------#



   #------- DBH. --------------------------------------------------------------------------#
   f.dbh.merged                      <<- as.integer(10 * round(merged$dbh.1999,1))
   f.dbh.this                        <<- as.integer(10 * round(this$dbh.1999  ,1))
   #----- The DBH tag changes according to the census. ------------------------------------#
   if (survey.years[1] < 2001){
      o.dbh.merged                   <<- as.integer(10 * round(merged$dbh.1999,1))
      o.dbh.this                     <<- as.integer(10 * round(this$dbh.1999  ,1))
      n.dbh.merged                   <<- as.integer(10 * round(merged$dbh.1999,1))
      n.dbh.this                     <<- as.integer(10 * round(this$dbh.1999  ,1))
   }else if (survey.years[1] < 2005){
      o.dbh.merged                   <<- as.integer(10 * round(merged$dbh.1999,1))
      o.dbh.this                     <<- as.integer(10 * round(this$dbh.1999  ,1))
      n.dbh.merged                   <<- as.integer(10 * round(merged$dbh.2001,1))
      n.dbh.this                     <<- as.integer(10 * round(this$dbh.2001  ,1))
   }else if (survey.years[1] < 2007){
      o.dbh.merged                   <<- as.integer(10 * round(merged$dbh.2001,1))
      o.dbh.this                     <<- as.integer(10 * round(this$dbh.2001  ,1))
      n.dbh.merged                   <<- as.integer(10 * round(merged$dbh.2005,1))
      n.dbh.this                     <<- as.integer(10 * round(this$dbh.2005  ,1))
   }else if (survey.years[1] < 2013){
      o.dbh.merged                   <<- as.integer(10 * round(merged$dbh.2005,1))
      o.dbh.this                     <<- as.integer(10 * round(this$dbh.2005  ,1))
      n.dbh.merged                   <<- as.integer(10 * round(merged$dbh.2007,1))
      n.dbh.this                     <<- as.integer(10 * round(this$dbh.2007  ,1))
   }else{
      o.dbh.merged                   <<- as.integer(10 * round(merged$dbh.2010,1))
      o.dbh.this                     <<- as.integer(10 * round(this$dbh.2010  ,1))
      n.dbh.merged                   <<- as.integer(10 * round(merged$dbh.2012,1))
      n.dbh.this                     <<- as.integer(10 * round(this$dbh.2012  ,1))
   }#end if
   f.dbh.merged[is.na(f.dbh.merged)] <<- "MM"
   o.dbh.merged[is.na(o.dbh.merged)] <<- "MM"
   n.dbh.merged[is.na(n.dbh.merged)] <<- "MM"
   f.dbh.this  [is.na(f.dbh.this  )] <<- "TT"
   o.dbh.this  [is.na(o.dbh.this  )] <<- "TT"
   n.dbh.this  [is.na(n.dbh.this  )] <<- "TT"
   #---------------------------------------------------------------------------------------#



   #------- DBH (alternative). ------------------------------------------------------------#
   if (survey.years[1] == 2013){
      x.dbh.merged                      <<- as.integer(10 * round(merged$dbh.2010,1))
      x.dbh.this                        <<- as.integer(10 * round(this$dbh.2010  ,1))
      y.dbh.merged                      <<- as.integer(10 * round(merged$dbh.2011,1))
      y.dbh.this                        <<- as.integer(10 * round(this$dbh.2011  ,1))
      z.dbh.merged                      <<- as.integer(10 * round(merged$dbh.2012,1))
      z.dbh.this                        <<- as.integer(10 * round(this$dbh.2012  ,1))
      x.dbh.merged[is.na(x.dbh.merged)] <<- "MM"
      x.dbh.this  [is.na(x.dbh.this  )] <<- "MM"
      y.dbh.merged[is.na(y.dbh.merged)] <<- "MM"
      y.dbh.this  [is.na(y.dbh.this  )] <<- "TT"
      z.dbh.merged[is.na(z.dbh.merged)] <<- "TT"
      z.dbh.this  [is.na(z.dbh.this  )] <<- "TT"
   }else{
      #----- Don't use alternative dbh tag for any other survey. --------------------------#
      x.dbh.merged <<- rep("MM",times=nrow(merged))
      x.dbh.this   <<- rep("MM",times=nrow(merged))
      y.dbh.merged <<- rep("MM",times=nrow(merged))
      y.dbh.this   <<- rep("TT",times=nrow(this  ))
      z.dbh.merged <<- rep("TT",times=nrow(this  ))
      z.dbh.this   <<- rep("TT",times=nrow(this  ))
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#



   #------- DBH (alternative). ------------------------------------------------------------#
   if (survey.years[1] == 2013){
      i.dbh.merged                      <<- as.integer(10 * round(merged$dbh.1999,1))
      i.dbh.this                        <<- as.integer(10 * round(this$dbh.2005  ,1))
      j.dbh.merged                      <<- as.integer(10 * round(merged$dbh.2010,1))
      j.dbh.this                        <<- as.integer(10 * round(this$dbh.1999  ,1))
      k.dbh.merged                      <<- as.integer(10 * round(merged$dbh.2010,1))
      k.dbh.this                        <<- as.integer(10 * round(this$dbh.2010  ,1))
      i.dbh.merged[is.na(i.dbh.merged)] <<- "MM"
      i.dbh.this  [is.na(i.dbh.this  )] <<- "MM"
      j.dbh.merged[is.na(j.dbh.merged)] <<- "MM"
      j.dbh.this  [is.na(j.dbh.this  )] <<- "TT"
      k.dbh.merged[is.na(k.dbh.merged)] <<- "TT"
      k.dbh.this  [is.na(k.dbh.this  )] <<- "TT"
   }else{
      #----- Don't use alternative dbh tag for any other survey. --------------------------#
      i.dbh.merged <<- rep("MM",times=nrow(merged))
      i.dbh.this   <<- rep("MM",times=nrow(merged))
      j.dbh.merged <<- rep("MM",times=nrow(merged))
      j.dbh.this   <<- rep("TT",times=nrow(this  ))
      k.dbh.merged <<- rep("TT",times=nrow(this  ))
      k.dbh.this   <<- rep("TT",times=nrow(this  ))
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     List of flags.  Here the order matters, the first is the preferred method and the #
   # last is the desperate method.                                                         #
   #---------------------------------------------------------------------------------------#
   cat0("     * Create tags for matching.")
   a.trans.merged        = c("n.trans.merged"     )
   a.trans.this          = c("n.trans.this"       )
   a.scientific.merged   = c("n.scientific.merged","o.scientific.merged")
   a.scientific.this     = c("n.scientific.this"  ,"o.scientific.this"  )
   #---- Check whether there is any tag change. -------------------------------------------#
   if ( any(o.tag.merged != f.tag.merged,na.rm=TRUE) 
      & any(o.tag.merged != n.tag.merged,na.rm=TRUE)){
      a.tag.merged     = c("n.tag.merged"    , "f.tag.merged"    , "o.tag.merged"   )
   }else if ( any(f.tag.merged != n.tag.merged,na.rm=TRUE)){
      a.tag.merged     = c("n.tag.merged"    , "f.tag.merged"    )
   }else{
      a.tag.merged     = c("n.tag.merged"    )
   }#end if
   if ( any(o.tag.this   != n.tag.this,na.rm=TRUE)){
      a.tag.this       = c("n.tag.this", "o.tag.this")
   }else{
      a.tag.this       = c("n.tag.this"      )
   }#end if
   #---- Check whether there is any coordinate change. ------------------------------------#
   if ( any(o.xy.merged != f.xy.merged,na.rm=TRUE) 
      & any(o.xy.merged != n.xy.merged,na.rm=TRUE)){
      a.xy.merged     = c("n.xy.merged"    , "f.xy.merged"    , "o.xy.merged"   )
   }else if ( any(f.xy.merged != n.xy.merged,na.rm=TRUE)){
      a.xy.merged     = c("n.xy.merged"    , "f.xy.merged"    )
   }else{
      a.xy.merged     = c("n.xy.merged"    )
   }#end if
   a.xy.this        = c("n.xy.this"       )
   b.xy.this        = c("s.xy.this"       , "c.xy.this"       , "z.xy.this"      )
   #---- Check whether there is any common name change. -----------------------------------#
   if ( any(o.common.merged != f.common.merged,na.rm=TRUE) 
      & any(o.common.merged != n.common.merged,na.rm=TRUE)){
      a.common.merged     = c("n.common.merged", "f.common.merged", "o.common.merged")
   }else if ( any(f.common.merged != n.common.merged,na.rm=TRUE)){
      a.common.merged     = c("n.common.merged", "f.common.merged")
   }else{
      a.common.merged     = c("n.common.merged"    )
   }#end if
   a.common.this    = c("n.common.this"   )
   #---- Check whether there are various DBH to compare. ----------------------------------#
   if (survey.years[1] < 2001){
      a.dbh.merged  = c("f.dbh.merged"    )
      a.dbh.this    = c("f.dbh.this"      )
   }else if (survey.years[1] < 2005){
      a.dbh.merged  = c("n.dbh.merged"    , "f.dbh.merged"    )
      a.dbh.this    = c("n.dbh.this"      , "f.dbh.this"      )
   }else{
      a.dbh.merged  = c("n.dbh.merged"    , "o.dbh.merged"    , "f.dbh.merged"   )
      a.dbh.this    = c("n.dbh.this"      , "o.dbh.this"      , "f.dbh.this"     )
   }#end if
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Initialise the tag.                                                               #
   #---------------------------------------------------------------------------------------#
   uni.tag          = list()
   u = 0
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     First we create the long tags that explore all variables.                         #
   #---------------------------------------------------------------------------------------#
   cat0("     * Full tags.")
   eg.merged = expand.grid( trans            = a.trans.merged
                          , scientific       = a.scientific.merged
                          , tag              = a.tag.merged
                          , xy               = a.xy.merged
                          , common           = a.common.merged
                          , dbh              = a.dbh.merged
                          , stringsAsFactors = FALSE
                          )#end expand.grid
   eg.this   = expand.grid( trans            = a.trans.this
                          , scientific       = a.scientific.this
                          , tag              = a.tag.this
                          , xy               = a.xy.this
                          , common           = a.common.this
                          , dbh              = a.dbh.this
                          , stringsAsFactors = FALSE
                          )#end expand.grid
   if (is.null(dim(eg.merged))) eg.merged = matrix(eg.merged,nrow=1)
   if (is.null(dim(eg.this  ))) eg.this   = matrix(eg.this  ,nrow=1)

   combo     = expand.grid( merged = sequence(nrow(eg.merged))
                          , this   = sequence(nrow(eg.this  ))
                          )#end expand.grid
   #----- Loop over all combinations. -----------------------------------------------------#
   for (co in sequence(nrow(combo))){
      u                 = u + 1
      em                = combo[co,1]
      et                = combo[co,2]
      retrieved.merged  = t(t(sapply(X=eg.merged[em,],FUN=get)))
      retrieved.this    = t(t(sapply(X=eg.this  [et,],FUN=get)))
      t.merged          = apply(X=retrieved.merged,MARGIN=1,FUN=paste,collapse="+")
      t.this            = apply(X=retrieved.this  ,MARGIN=1,FUN=paste,collapse="+")
      #----- Save the tags to the output list. --------------------------------------------#
      uni.tag[[u]]      = list(t.merged = t.merged, t.this = t.this)
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Now we only check the transect, the tag, and the coordinates.                     #
   #---------------------------------------------------------------------------------------#
   cat0("     * Transect, tags, xy.")
   eg.merged = expand.grid( trans            = a.trans.merged
                          , scientific       = a.scientific.merged
                          , tag              = a.tag.merged
                          , xy               = a.xy.merged
                          , stringsAsFactors = FALSE
                          )#end expand.grid
   eg.this   = expand.grid( trans            = a.trans.this
                          , scientific       = a.scientific.this
                          , tag              = a.tag.this
                          , xy               = a.xy.this
                          , stringsAsFactors = FALSE
                          )#end expand.grid
   if (is.null(dim(eg.merged))) eg.merged = matrix(eg.merged,nrow=1)
   if (is.null(dim(eg.this  ))) eg.this   = matrix(eg.this  ,nrow=1)
   combo     = expand.grid( merged = sequence(nrow(eg.merged))
                          , this   = sequence(nrow(eg.this  ))
                          )#end expand.grid
   #----- Loop over all combinations. -----------------------------------------------------#
   for (co in sequence(nrow(combo))){
      u                 = u + 1
      em                = combo[co,1]
      et                = combo[co,2]
      retrieved.merged  = t(t(sapply(X=eg.merged[em,],FUN=get)))
      retrieved.this    = t(t(sapply(X=eg.this  [et,],FUN=get)))
      t.merged          = apply(X=retrieved.merged,MARGIN=1,FUN=paste,collapse="+")
      t.this            = apply(X=retrieved.this  ,MARGIN=1,FUN=paste,collapse="+")
      #----- Save the tags to the output list. --------------------------------------------#
      uni.tag[[u]]      = list(t.merged = t.merged, t.this = t.this)
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Now we check only transect, tag, common name, and DBH (no XY).                    #
   #---------------------------------------------------------------------------------------#
   cat0("     * Transect, tags, common, dbh.")
   eg.merged = expand.grid( trans            = a.trans.merged
                          , scientific       = a.scientific.merged
                          , tag              = a.tag.merged
                          , common           = a.common.merged
                          , dbh              = a.dbh.merged
                          , stringsAsFactors = FALSE
                          )#end expand.grid
   eg.this   = expand.grid( trans            = a.trans.this
                          , scientific       = a.scientific.this
                          , tag              = a.tag.this  
                          , common           = a.common.this  
                          , dbh              = a.dbh.this  
                          , stringsAsFactors = FALSE
                          )#end expand.grid
   if (is.null(dim(eg.merged))) eg.merged = matrix(eg.merged,nrow=1)
   if (is.null(dim(eg.this  ))) eg.this   = matrix(eg.this  ,nrow=1)
   combo     = expand.grid( merged = sequence(nrow(eg.merged))
                          , this   = sequence(nrow(eg.this  ))
                          )#end expand.grid
   #----- Loop over all combinations. -----------------------------------------------------#
   for (co in sequence(nrow(combo))){
      u                 = u + 1
      em                = combo[co,1]
      et                = combo[co,2]
      retrieved.merged  = t(t(sapply(X=eg.merged[em,],FUN=get)))
      retrieved.this    = t(t(sapply(X=eg.this  [et,],FUN=get)))
      t.merged          = apply(X=retrieved.merged,MARGIN=1,FUN=paste,collapse="+")
      t.this            = apply(X=retrieved.this  ,MARGIN=1,FUN=paste,collapse="+")
      #----- Save the tags to the output list. --------------------------------------------#
      uni.tag[[u]]      = list(t.merged = t.merged, t.this = t.this)
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Last we check for individuals with swapped coordinates but otherwise same tag and #
   # same dbh.                                                                             #
   #---------------------------------------------------------------------------------------#
   cat0("     * Transect, tags, swapped coordinates, common.")
   eg.merged = expand.grid( trans            = a.trans.merged
                          , scientific       = a.scientific.merged
                          , tag              = a.tag.merged
                          , xy               = a.xy.merged
                          , common           = a.common.merged
                          , stringsAsFactors = FALSE
                          )#end expand.grid
   eg.this   = expand.grid( trans            = a.trans.this
                          , scientific       = a.scientific.this
                          , tag              = a.tag.this
                          , xy               = b.xy.this
                          , common           = a.common.this
                          , stringsAsFactors = FALSE
                          )#end expand.grid
   if (is.null(dim(eg.merged))) eg.merged = matrix(eg.merged,nrow=1)
   if (is.null(dim(eg.this  ))) eg.this   = matrix(eg.this  ,nrow=1)
   combo     = expand.grid( merged = sequence(nrow(eg.merged))
                          , this   = sequence(nrow(eg.this  ))
                          )#end expand.grid
   #----- Loop over all combinations. -----------------------------------------------------#
   for (co in sequence(nrow(combo))){
      u                 = u + 1
      em                = combo[co,1]
      et                = combo[co,2]
      retrieved.merged  = t(t(sapply(X=eg.merged[em,],FUN=get)))
      retrieved.this    = t(t(sapply(X=eg.this  [et,],FUN=get)))
      t.merged          = apply(X=retrieved.merged,MARGIN=1,FUN=paste,collapse="+")
      t.this            = apply(X=retrieved.this  ,MARGIN=1,FUN=paste,collapse="+")
      #----- Save the tags to the output list. --------------------------------------------#
      uni.tag[[u]]      = list(t.merged = t.merged, t.this = t.this)
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Now we only check the transect, the tag, and three DBHs.                          #
   #---------------------------------------------------------------------------------------#
   cat0("     * Transect, tags, three DBHs (1999, 2005, 2010).")
   eg.merged = expand.grid( trans            = a.trans.merged
                          , scientific       = a.scientific.merged
                          , tag              = a.tag.merged
                          , idbh             = "i.dbh.merged"
                          , jdbh             = "j.dbh.merged"
                          , kdbh             = "k.dbh.merged"
                          , stringsAsFactors = FALSE
                          )#end expand.grid
   eg.this   = expand.grid( trans            = a.trans.this
                          , scientific       = a.scientific.this
                          , tag              = a.tag.this
                          , idbh             = "i.dbh.this"
                          , jdbh             = "j.dbh.this"
                          , kdbh             = "k.dbh.this"
                          , stringsAsFactors = FALSE
                          )#end expand.grid
   if (is.null(dim(eg.merged))) eg.merged = matrix(eg.merged,nrow=1)
   if (is.null(dim(eg.this  ))) eg.this   = matrix(eg.this  ,nrow=1)
   combo     = expand.grid( merged = sequence(nrow(eg.merged))
                          , this   = sequence(nrow(eg.this  ))
                          )#end expand.grid
   #----- Loop over all combinations. -----------------------------------------------------#
   for (co in sequence(nrow(combo))){
      u                 = u + 1
      em                = combo[co,1]
      et                = combo[co,2]
      retrieved.merged  = t(t(sapply(X=eg.merged[em,],FUN=get)))
      retrieved.this    = t(t(sapply(X=eg.this  [et,],FUN=get)))
      t.merged          = apply(X=retrieved.merged,MARGIN=1,FUN=paste,collapse="+")
      t.this            = apply(X=retrieved.this  ,MARGIN=1,FUN=paste,collapse="+")
      #----- Save the tags to the output list. --------------------------------------------#
      uni.tag[[u]]      = list(t.merged = t.merged, t.this = t.this)
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Now we only check the transect, the tag, and three DBHs.                          #
   #---------------------------------------------------------------------------------------#
   cat0("     * Transect, tags, three DBHs (2010, 2011, 2012).")
   eg.merged = expand.grid( trans            = a.trans.merged
                          , scientific       = a.scientific.merged
                          , tag              = a.tag.merged
                          , xdbh             = "x.dbh.merged"
                          , ydbh             = "y.dbh.merged"
                          , zdbh             = "z.dbh.merged"
                          , stringsAsFactors = FALSE
                          )#end expand.grid
   eg.this   = expand.grid( trans            = a.trans.this
                          , scientific       = a.scientific.this
                          , tag              = a.tag.this
                          , xdbh             = "x.dbh.this"
                          , ydbh             = "y.dbh.this"
                          , zdbh             = "z.dbh.this"
                          , stringsAsFactors = FALSE
                          )#end expand.grid
   if (is.null(dim(eg.merged))) eg.merged = matrix(eg.merged,nrow=1)
   if (is.null(dim(eg.this  ))) eg.this   = matrix(eg.this  ,nrow=1)
   combo     = expand.grid( merged = sequence(nrow(eg.merged))
                          , this   = sequence(nrow(eg.this  ))
                          )#end expand.grid
   #----- Loop over all combinations. -----------------------------------------------------#
   for (co in sequence(nrow(combo))){
      u                 = u + 1
      em                = combo[co,1]
      et                = combo[co,2]
      retrieved.merged  = t(t(sapply(X=eg.merged[em,],FUN=get)))
      retrieved.this    = t(t(sapply(X=eg.this  [et,],FUN=get)))
      t.merged          = apply(X=retrieved.merged,MARGIN=1,FUN=paste,collapse="+")
      t.this            = apply(X=retrieved.this  ,MARGIN=1,FUN=paste,collapse="+")
      #----- Save the tags to the output list. --------------------------------------------#
      uni.tag[[u]]      = list(t.merged = t.merged, t.this = t.this)
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Now we only check the transect, the tag, and two DBHs.                            #
   #---------------------------------------------------------------------------------------#
   cat0("     * Transect, tags, two DBHs (2010,2012).")
   eg.merged = expand.grid( trans            = a.trans.merged
                          , scientific       = a.scientific.merged
                          , tag              = a.tag.merged
                          , xdbh             = "x.dbh.merged"
                          , zdbh             = "z.dbh.merged"
                          , stringsAsFactors = FALSE
                          )#end expand.grid
   eg.this   = expand.grid( trans            = a.trans.this
                          , scientific       = a.scientific.this
                          , tag              = a.tag.this
                          , xdbh             = "x.dbh.this"
                          , zdbh             = "z.dbh.this"
                          , stringsAsFactors = FALSE
                          )#end expand.grid
   if (is.null(dim(eg.merged))) eg.merged = matrix(eg.merged,nrow=1)
   if (is.null(dim(eg.this  ))) eg.this   = matrix(eg.this  ,nrow=1)
   combo     = expand.grid( merged = sequence(nrow(eg.merged))
                          , this   = sequence(nrow(eg.this  ))
                          )#end expand.grid
   #----- Loop over all combinations. -----------------------------------------------------#
   for (co in sequence(nrow(combo))){
      u                 = u + 1
      em                = combo[co,1]
      et                = combo[co,2]
      retrieved.merged  = t(t(sapply(X=eg.merged[em,],FUN=get)))
      retrieved.this    = t(t(sapply(X=eg.this  [et,],FUN=get)))
      t.merged          = apply(X=retrieved.merged,MARGIN=1,FUN=paste,collapse="+")
      t.this            = apply(X=retrieved.this  ,MARGIN=1,FUN=paste,collapse="+")
      #----- Save the tags to the output list. --------------------------------------------#
      uni.tag[[u]]      = list(t.merged = t.merged, t.this = t.this)
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Return the tags.                                                                  #
   #---------------------------------------------------------------------------------------#
   return (uni.tag)
   #---------------------------------------------------------------------------------------#
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function is supposed to work for Paracou census only, feel free to adapt it     #
# elsewhere.  It creates a collection of unique tags based on various variables such as    #
# tags, coordinates, common name, and the older values that these variables once had.      #
#------------------------------------------------------------------------------------------#
census.tagger.gyf <<- function(merged,this,survey.years){
   #=======================================================================================#
   #=======================================================================================#
   #      Create all possible tags.  We must declare them as global variables because they #
   # must be accessed from sapply.                                                         #
   #---------------------------------------------------------------------------------------#


   #----- Transects. ----------------------------------------------------------------------#
   n.trans.merged        <<- merged$trans
   n.trans.this          <<- this$trans
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Scientific name flag, use the scientific name.                                   #
   #---------------------------------------------------------------------------------------#
   n.scientific.merged   <<- merged$scientific
   n.scientific.this     <<- this$scientific
   #---------------------------------------------------------------------------------------#


   #----- Coordinates. --------------------------------------------------------------------#
   n.x.merged                     = as.integer(10 * merged$x)
   o.x.merged                     = as.integer(10 * merged$old.x)
   f.x.merged                     = as.integer(10 * merged$x.1st)
   n.x.this                       = as.integer(10 * this$x)
   n.y.merged                     = as.integer(10 * merged$y)
   o.y.merged                     = as.integer(10 * merged$old.y)
   f.y.merged                     = as.integer(10 * merged$y.1st)
   n.y.this                       = as.integer(10 * this$y)
   f.x.merged [is.na(f.x.merged)] = n.x.merged[is.na(f.x.merged)     ]
   f.y.merged [is.na(f.y.merged)] = n.y.merged[is.na(f.y.merged)     ]
   o.x.merged [is.na(o.x.merged)] = f.x.merged[is.na(o.x.merged)     ]
   o.y.merged [is.na(o.y.merged)] = f.y.merged[is.na(o.y.merged)     ]
   n.xy.merged <<- paste("x",n.x.merged,"+y",n.y.merged,sep="")
   o.xy.merged <<- paste("x",o.x.merged,"+y",o.y.merged,sep="")
   f.xy.merged <<- paste("x",f.x.merged,"+y",f.y.merged,sep="")
   n.xy.this   <<- paste("x",n.x.this  ,"+y",n.y.this  ,sep="")
   #---------------------------------------------------------------------------------------#




   #------ Tags. --------------------------------------------------------------------------#
   n.tag.merged                      <<- merged$tag
   o.tag.merged                      <<- merged$tag.1st
   f.tag.merged                      <<- merged$tag.1st
   n.tag.this                        <<- this$tag
   o.tag.this                        <<- this$old.tag
   f.tag.merged[is.na(f.tag.merged)] <<- n.tag.merged   [is.na(f.tag.merged)]
   o.tag.merged[is.na(o.tag.merged)] <<- f.tag.merged   [is.na(o.tag.merged)]
   o.tag.this  [is.na(o.tag.this)  ] <<- n.tag.this     [is.na(o.tag.this)  ]
   #---------------------------------------------------------------------------------------#



   #------- gbh. --------------------------------------------------------------------------#
   f.gbh.merged                      <<- as.integer(10 * round(merged$gbh.2004,1))
   f.gbh.this                        <<- as.integer(10 * round(this$gbh.2004  ,1))
   #----- The gbh tag changes according to the census. ------------------------------------#
   if (survey.years[1] < 2006){
      o.gbh.merged                   <<- as.integer(10 * round(merged$gbh.2004,1))
      o.gbh.this                     <<- as.integer(10 * round(this$gbh.2004  ,1))
      n.gbh.merged                   <<- as.integer(10 * round(merged$gbh.2004,1))
      n.gbh.this                     <<- as.integer(10 * round(this$gbh.2004  ,1))
   }else if (survey.years[1] < 2008){
      o.gbh.merged                   <<- as.integer(10 * round(merged$gbh.2004,1))
      o.gbh.this                     <<- as.integer(10 * round(this$gbh.2004  ,1))
      n.gbh.merged                   <<- as.integer(10 * round(merged$gbh.2006,1))
      n.gbh.this                     <<- as.integer(10 * round(this$gbh.2006  ,1))
   }else if (survey.years[1] < 2010){
      o.gbh.merged                   <<- as.integer(10 * round(merged$gbh.2006,1))
      o.gbh.this                     <<- as.integer(10 * round(this$gbh.2006  ,1))
      n.gbh.merged                   <<- as.integer(10 * round(merged$gbh.2008,1))
      n.gbh.this                     <<- as.integer(10 * round(this$gbh.2008  ,1))
   }else{
      o.gbh.merged                   <<- as.integer(10 * round(merged$gbh.2008,1))
      o.gbh.this                     <<- as.integer(10 * round(this$gbh.2008  ,1))
      n.gbh.merged                   <<- as.integer(10 * round(merged$gbh.2010,1))
      n.gbh.this                     <<- as.integer(10 * round(this$gbh.2010  ,1))
   }#end if
   f.gbh.merged[is.na(f.gbh.merged)] <<- "MM"
   o.gbh.merged[is.na(o.gbh.merged)] <<- "MM"
   n.gbh.merged[is.na(n.gbh.merged)] <<- "MM"
   f.gbh.this  [is.na(f.gbh.this  )] <<- "TT"
   o.gbh.this  [is.na(o.gbh.this  )] <<- "TT"
   n.gbh.this  [is.na(n.gbh.this  )] <<- "TT"
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     List of flags.  Here the order matters, the first is the preferred method and the #
   # last is the desperate method.                                                         #
   #---------------------------------------------------------------------------------------#
   cat("     * Create tags for matching...","\n")
   a.trans.merged        = c("n.trans.merged"     )
   a.trans.this          = c("n.trans.this"       )
   a.scientific.merged   = c("n.scientific.merged")
   a.scientific.this     = c("n.scientific.this"  )
   #---- Check whether there is any tag change. -------------------------------------------#
   if ( any(o.tag.merged != f.tag.merged,na.rm=TRUE) 
      & any(o.tag.merged != n.tag.merged,na.rm=TRUE)){
      a.tag.merged     = c("n.tag.merged"    , "f.tag.merged"    , "o.tag.merged"   )
   }else if ( any(f.tag.merged != n.tag.merged,na.rm=TRUE)){
      a.tag.merged     = c("n.tag.merged"    , "f.tag.merged"    )
   }else{
      a.tag.merged     = c("n.tag.merged"    )
   }#end if
   if ( any(o.tag.this   != n.tag.this,na.rm=TRUE)){
      a.tag.this       = c("n.tag.this", "o.tag.this")
   }else{
      a.tag.this       = c("n.tag.this"      )
   }#end if
   #---- Check whether there is any coordinate change. ------------------------------------#
   if ( any(o.xy.merged != f.xy.merged,na.rm=TRUE) 
      & any(o.xy.merged != n.xy.merged,na.rm=TRUE)){
      a.xy.merged     = c("n.xy.merged"    , "f.xy.merged"    , "o.xy.merged"   )
   }else if ( any(f.xy.merged != n.xy.merged,na.rm=TRUE)){
      a.xy.merged     = c("n.xy.merged"    , "f.xy.merged"    )
   }else{
      a.xy.merged     = c("n.xy.merged"    )
   }#end if
   a.xy.this          = c("n.xy.this"       )
   #---- Check whether there are various gbh to compare. ----------------------------------#
   if (survey.years[1] < 2006){
      a.gbh.merged  = c("f.gbh.merged"    )
      a.gbh.this    = c("f.gbh.this"      )
   }else if (survey.years[1] < 2008){
      a.gbh.merged  = c("n.gbh.merged"    , "f.gbh.merged"    )
      a.gbh.this    = c("n.gbh.this"      , "f.gbh.this"      )
   }else{
      a.gbh.merged  = c("n.gbh.merged"    , "o.gbh.merged"    , "f.gbh.merged"   )
      a.gbh.this    = c("n.gbh.this"      , "o.gbh.this"      , "f.gbh.this"     )
   }#end if
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Initialise the tag.                                                               #
   #---------------------------------------------------------------------------------------#
   uni.tag          = list()
   u = 0
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     First we create the long tags that explore all variables.                         #
   #---------------------------------------------------------------------------------------#
   cat("     * Full tags...","\n")
   eg.merged = expand.grid( trans       = a.trans.merged
                          , scientific  = a.scientific.merged
                          , tag         = a.tag.merged
                          , xy          = a.xy.merged
                          , gbh         = a.gbh.merged
                          )#end expand.grid
   eg.this   = expand.grid( trans       = a.trans.this
                          , scientific  = a.scientific.this
                          , tag         = a.tag.this  
                          , xy          = a.xy.this  
                          , gbh         = a.gbh.this  
                          )#end expand.grid
   eg.merged = sapply(eg.merged,as.character)
   eg.this   = sapply(eg.this  ,as.character)
   if (is.null(dim(eg.merged))) eg.merged = matrix(eg.merged,nrow=1)
   if (is.null(dim(eg.this  ))) eg.this   = matrix(eg.this  ,nrow=1)

   combo     = expand.grid( merged = sequence(nrow(eg.merged))
                          , this   = sequence(nrow(eg.this  ))
                          )#end expand.grid
   #----- Loop over all combinations. -----------------------------------------------------#
   for (co in 1:nrow(combo)){
      u                 = u + 1
      em                = combo[co,1]
      et                = combo[co,2]
      retrieved.merged  = t(t(sapply(X=eg.merged[em,],FUN=get)))
      retrieved.this    = t(t(sapply(X=eg.this  [et,],FUN=get)))
      t.merged          = apply(X=retrieved.merged,MARGIN=1,FUN=paste,collapse="+")
      t.this            = apply(X=retrieved.this  ,MARGIN=1,FUN=paste,collapse="+")
      #----- Save the tags to the output list. --------------------------------------------#
      uni.tag[[u]]      = list(t.merged = t.merged, t.this = t.this)
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Now we only check the transect, the tag, and the coordinates.                     #
   #---------------------------------------------------------------------------------------#
   cat("     * Trans, scientific, tag, xy...","\n")
   eg.merged = expand.grid( trans       = a.trans.merged
                          , scientific  = a.scientific.merged
                          , tag         = a.tag.merged
                          , xy          = a.xy.merged
                          )#end expand.grid
   eg.this   = expand.grid( trans       = a.trans.this
                          , scientific  = a.scientific.this
                          , tag         = a.tag.this
                          , xy          = a.xy.this
                          )#end expand.grid
   eg.merged = sapply(eg.merged,as.character)
   eg.this   = sapply(eg.this  ,as.character)
   if (is.null(dim(eg.merged))) eg.merged = matrix(eg.merged,nrow=1)
   if (is.null(dim(eg.this  ))) eg.this   = matrix(eg.this  ,nrow=1)
   combo     = expand.grid( merged = sequence(nrow(eg.merged))
                          , this   = sequence(nrow(eg.this  ))
                          )#end expand.grid
   #----- Loop over all combinations. -----------------------------------------------------#
   for (co in 1:nrow(combo)){
      u                 = u + 1
      em                = combo[co,1]
      et                = combo[co,2]
      retrieved.merged  = t(t(sapply(X=eg.merged[em,],FUN=get)))
      retrieved.this    = t(t(sapply(X=eg.this  [et,],FUN=get)))
      t.merged          = apply(X=retrieved.merged,MARGIN=1,FUN=paste,collapse="+")
      t.this            = apply(X=retrieved.this  ,MARGIN=1,FUN=paste,collapse="+")
      #----- Save the tags to the output list. --------------------------------------------#
      uni.tag[[u]]      = list(t.merged = t.merged, t.this = t.this)
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Now we check trans, scientific, xy, and gbh (no tags).                            #
   #---------------------------------------------------------------------------------------#
   cat("     * Trans, scientific, xy, gbh...","\n")
   eg.merged = expand.grid( trans       = a.trans.merged
                          , scientific  = a.scientific.merged
                          , xy          = a.xy.merged
                          , gbh         = a.gbh.merged
                          )#end expand.grid
   eg.this   = expand.grid( trans       = a.trans.this
                          , scientific  = a.scientific.this
                          , xy          = a.xy.this  
                          , gbh         = a.gbh.this  
                          )#end expand.grid
   eg.merged = sapply(eg.merged,as.character)
   eg.this   = sapply(eg.this  ,as.character)
   if (is.null(dim(eg.merged))) eg.merged = matrix(eg.merged,nrow=1)
   if (is.null(dim(eg.this  ))) eg.this   = matrix(eg.this  ,nrow=1)

   combo     = expand.grid( merged = sequence(nrow(eg.merged))
                          , this   = sequence(nrow(eg.this  ))
                          )#end expand.grid
   #----- Loop over all combinations. -----------------------------------------------------#
   for (co in 1:nrow(combo)){
      u                 = u + 1
      em                = combo[co,1]
      et                = combo[co,2]
      retrieved.merged  = t(t(sapply(X=eg.merged[em,],FUN=get)))
      retrieved.this    = t(t(sapply(X=eg.this  [et,],FUN=get)))
      t.merged          = apply(X=retrieved.merged,MARGIN=1,FUN=paste,collapse="+")
      t.this            = apply(X=retrieved.this  ,MARGIN=1,FUN=paste,collapse="+")
      #----- Save the tags to the output list. --------------------------------------------#
      uni.tag[[u]]      = list(t.merged = t.merged, t.this = t.this)
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Now we only check the transect, the tag, and the coordinates.                     #
   #---------------------------------------------------------------------------------------#
   cat("     * Transect, tags, xy...","\n")
   eg.merged = expand.grid( trans       = a.trans.merged
                          , tag         = a.tag.merged
                          , xy          = a.xy.merged
                          )#end expand.grid
   eg.this   = expand.grid( trans       = a.trans.this
                          , tag         = a.tag.this
                          , xy          = a.xy.this
                          )#end expand.grid
   eg.merged = sapply(eg.merged,as.character)
   eg.this   = sapply(eg.this  ,as.character)
   if (is.null(dim(eg.merged))) eg.merged = matrix(eg.merged,nrow=1)
   if (is.null(dim(eg.this  ))) eg.this   = matrix(eg.this  ,nrow=1)
   combo     = expand.grid( merged = sequence(nrow(eg.merged))
                          , this   = sequence(nrow(eg.this  ))
                          )#end expand.grid
   #----- Loop over all combinations. -----------------------------------------------------#
   for (co in 1:nrow(combo)){
      u                 = u + 1
      em                = combo[co,1]
      et                = combo[co,2]
      retrieved.merged  = t(t(sapply(X=eg.merged[em,],FUN=get)))
      retrieved.this    = t(t(sapply(X=eg.this  [et,],FUN=get)))
      t.merged          = apply(X=retrieved.merged,MARGIN=1,FUN=paste,collapse="+")
      t.this            = apply(X=retrieved.this  ,MARGIN=1,FUN=paste,collapse="+")
      #----- Save the tags to the output list. --------------------------------------------#
      uni.tag[[u]]      = list(t.merged = t.merged, t.this = t.this)
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Now we check only transect, tag, name, and gbh (no XY).                           #
   #---------------------------------------------------------------------------------------#
   cat("     * Transect, tags, gbh...","\n")
   eg.merged = expand.grid( trans       = a.trans.merged
                          , tag         = a.tag.merged
                          , gbh         = a.gbh.merged
                          )#end expand.grid
   eg.this   = expand.grid( trans       = a.trans.this
                          , tag         = a.tag.this  
                          , gbh         = a.gbh.this  
                          )#end expand.grid
   eg.merged = sapply(eg.merged,as.character)
   eg.this   = sapply(eg.this  ,as.character)
   if (is.null(dim(eg.merged))) eg.merged = matrix(eg.merged,nrow=1)
   if (is.null(dim(eg.this  ))) eg.this   = matrix(eg.this  ,nrow=1)
   combo     = expand.grid( merged = sequence(nrow(eg.merged))
                          , this   = sequence(nrow(eg.this  ))
                          )#end expand.grid
   #----- Loop over all combinations. -----------------------------------------------------#
   for (co in 1:nrow(combo)){
      u                 = u + 1
      em                = combo[co,1]
      et                = combo[co,2]
      retrieved.merged  = t(t(sapply(X=eg.merged[em,],FUN=get)))
      retrieved.this    = t(t(sapply(X=eg.this  [et,],FUN=get)))
      t.merged          = apply(X=retrieved.merged,MARGIN=1,FUN=paste,collapse="+")
      t.this            = apply(X=retrieved.this  ,MARGIN=1,FUN=paste,collapse="+")
      #----- Save the tags to the output list. --------------------------------------------#
      uni.tag[[u]]      = list(t.merged = t.merged, t.this = t.this)
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Return the tags.                                                                  #
   #---------------------------------------------------------------------------------------#
   return (uni.tag)
   #---------------------------------------------------------------------------------------#
}#end function
#------------------------------------------------------------------------------------------#
