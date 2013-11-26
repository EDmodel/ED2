#==========================================================================================#
#==========================================================================================#
#     This function is supposed to work for Santarem census only, feel free to adapt it    #
# elsewhere.  It creates a collection of unique tags based on various variables such as    #
# tags, coordinates, common name, and the older values that these variables once had.      #
#------------------------------------------------------------------------------------------#
census.tagger.s67 = function(merged,this,survey.years){
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
   n.scientific.merged                           <<- rep("M" ,times=nrow(merged))
   n.scientific.this                             <<- rep("T" ,times=nrow(this)  )
   n.scientific.merged[is.na(merged$scientific)] <<- "T"
   n.scientific.this  [is.na(this$scientific)  ] <<- "M"
   o.scientific.merged                           <<- rep("M" ,times=nrow(merged))
   o.scientific.this                             <<- rep("T" ,times=nrow(this)  )
   o.scientific.merged[is.na(merged$scientific)] <<- "O"
   o.scientific.this  [is.na(this$scientific)  ] <<- "O"
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
   n.xy.merged <<- paste("x",n.x.merged,"+y",n.y.merged,sep="")
   o.xy.merged <<- paste("x",o.x.merged,"+y",o.y.merged,sep="")
   f.xy.merged <<- paste("x",f.x.merged,"+y",f.y.merged,sep="")
   n.xy.this   <<- paste("x",n.x.this  ,"+y",n.y.this  ,sep="")
   s.xy.this   <<- paste("x",s.x.this  ,"+y",n.y.this  ,sep="")
   c.xy.this   <<- paste("x",n.x.this  ,"+y",s.y.this  ,sep="")
   z.xy.this   <<- paste("x",s.x.this  ,"+y",s.y.this  ,sep="")
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
   n.common.merged                         <<- merged$common
   n.common.this                           <<- this$common
   f.common.merged                         <<- merged$common.1st
   o.common.merged                         <<- merged$old.common
   f.common.merged[is.na(f.common.merged)] <<- n.common.merged[is.na(f.common.merged)]
   o.common.merged[is.na(o.common.merged)] <<- f.common.merged[is.na(o.common.merged)]
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
   }else{
      o.dbh.merged                   <<- as.integer(10 * round(merged$dbh.2005,1))
      o.dbh.this                     <<- as.integer(10 * round(this$dbh.2005  ,1))
      n.dbh.merged                   <<- as.integer(10 * round(merged$dbh.2007,1))
      n.dbh.this                     <<- as.integer(10 * round(this$dbh.2007  ,1))
   }#end if
   f.dbh.merged[is.na(f.dbh.merged)] <<- "MM"
   o.dbh.merged[is.na(o.dbh.merged)] <<- "MM"
   n.dbh.merged[is.na(n.dbh.merged)] <<- "MM"
   f.dbh.this  [is.na(f.dbh.this  )] <<- "TT"
   o.dbh.this  [is.na(o.dbh.this  )] <<- "TT"
   n.dbh.this  [is.na(n.dbh.this  )] <<- "TT"
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     List of flags.  Here the order matters, the first is the preferred method and the #
   # last is the desperate method.                                                         #
   #---------------------------------------------------------------------------------------#
   cat("     * Create tags for matching...","\n")
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
   cat("     * Full tags...","\n")
   eg.merged = expand.grid( trans       = a.trans.merged
                          , scientific  = a.scientific.merged
                          , tag         = a.tag.merged
                          , xy          = a.xy.merged
                          , common      = a.common.merged
                          , dbh         = a.dbh.merged
                          )#end expand.grid
   eg.this   = expand.grid( trans       = a.trans.this
                          , scientific  = a.scientific.this
                          , tag         = a.tag.this  
                          , xy          = a.xy.this  
                          , common      = a.common.this  
                          , dbh         = a.dbh.this  
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
   #     Now we check only transect, tag, common name, and DBH (no XY).                    #
   #---------------------------------------------------------------------------------------#
   cat("     * Transect, tags, common, dbh...","\n")
   eg.merged = expand.grid( trans       = a.trans.merged
                          , scientific  = a.scientific.merged
                          , tag         = a.tag.merged
                          , common      = a.common.merged
                          , dbh         = a.dbh.merged
                          )#end expand.grid
   eg.this   = expand.grid( trans       = a.trans.this
                          , scientific  = a.scientific.this
                          , tag         = a.tag.this  
                          , common      = a.common.this  
                          , dbh         = a.dbh.this  
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
   #     Last we check for individuals with swapped coordinates but otherwise same tag and #
   # same dbh.                                                                             #
   #---------------------------------------------------------------------------------------#
   cat("     * Transect, tags, swapped coordinates, common...","\n")
   eg.merged = expand.grid( trans       = a.trans.merged
                          , scientific  = a.scientific.merged
                          , tag         = a.tag.merged
                          , xy          = a.xy.merged
                          , common      = a.common.merged
                          )#end expand.grid
   eg.this   = expand.grid( trans       = a.trans.this
                          , scientific  = a.scientific.this
                          , tag         = a.tag.this
                          , xy          = b.xy.this
                          , common      = a.common.this
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
