#==========================================================================================#
#==========================================================================================#
#     Function read.ptable.r                                                               #
#                                                                                          #
#     This function reads in partial tables, by skipping additional columns.  Unlike       #
# read.table, here you may also specify the number of columns.  If you don't, then it is   #
# the same thing as read.table.                                                            #
#                                                                                          #
#------------------------------------------------------------------------------------------#
read.ptable <<- function( file
                        , header           = FALSE
                        , sep              = ""
                        , quote            = "\"'"
                        , comment.char     = "#"
                        , dec              = "."
                        , row.names
                        , col.names
                        , as.is            = !stringsAsFactors
                        , na.strings       = "NA"
                        , colClasses       = NA
                        , nrows            = -1
                        , ncols            = if (  all(is.na(colClasses))
                                                && missing(col.names)){
                                                -1
                                             }else if(all(is.na(colClasses))){
                                                length(col.names)
                                             }else{
                                                length(colClasses)
                                             }#end if
                        , skip             =  0
                        , check.names      = TRUE
                        , fill             = !blank.lines.skip
                        , strip.white      = FALSE
                        , blank.lines.skip = TRUE
                        , allowEscapes     = FALSE
                        , flush = FALSE
                        , stringsAsFactors = default.stringsAsFactors()
                        , fileEncoding     = ""
                        , encoding         = "unknown"
                        , text
                        ){

   #------ If ncols is -1, then this should fall back to traditional read.table. ----------#
   if ((! is.null(ncols)) & is.finite(ncols) && ncols > 0){
      #---- Find out how many columns there are in each column. ---------------------------#
      cidx = count.fields( file             = file
                         , sep              = sep
                         , quote            = quote
                         , skip             = skip
                         , blank.lines.skip = blank.lines.skip
                         , comment.char     = comment.char
                         )#end count.fields

      #------------------------------------------------------------------------------------#
      #     Find out how many rows and columns exist.  If the data is going to be cropped  #
      # (nrows is valid), then we reduce the number of rows.                               #
      #------------------------------------------------------------------------------------#
      nc   = max     (cidx)
      nr   = length  (cidx)
      if ((! is.null(nrows)) & is.finite(nrows) && nrows > 0) nr = min(nr,nrows)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Initialise the matrix, and define the indices for rows.                       #
      #------------------------------------------------------------------------------------#
      ans  = matrix(nrow=nr,ncol=nc)
      ridx = sequence(nr)
      #------------------------------------------------------------------------------------#



      #-----Find the indices for a nice table, the remainder will be NA. ------------------#
      idx  = cbind( row = c(unlist(mapply(FUN=rep,x=ridx,times=cidx[ridx])))
                  , col = c(unlist(mapply(FUN=sequence,cidx[ridx])))
                  )#end cbind
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Read the header using "scan".                                                  #
      #------------------------------------------------------------------------------------#
      if (header && missing(col.names)){
         col.names = scan( file             = file
                         , what             = character()
                         , nmax             = min(nc,ncols)
                         , sep              = sep
                         , quote            = if(identical(sep, "\n")){""}else{quote}
                         , dec              = dec
                         , skip             = skip
                         , nlines           = 1
                         , na.strings       = na.strings
                         , flush            = flush
                         , fill             = fill
                         , strip.white      = strip.white
                         , quiet            = TRUE
                         , blank.lines.skip = blank.lines.skip
                         , multi.line       = TRUE
                         , comment.char     = comment.char
                         , allowEscapes     = FALSE
                         , fileEncoding     = fileEncoding
                         , encoding         = encoding
                         , text             = text
                         )#end scan
      }else if (missing(col.names)){
         col.names = paste("V",sequence(min(nc,ncols)),sep=".")
      }#end if
      skip.tot = skip + header
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Read the data using "scan".                                                    #
      #------------------------------------------------------------------------------------#
      ans[idx] = scan( file             = file
                     , what             = character()
                     , nmax             = -1
                     , sep              = sep
                     , quote            = if(identical(sep, "\n")){""}else{quote}
                     , dec              = dec
                     , skip             = skip.tot
                     , nlines           = nr
                     , na.strings       = na.strings
                     , flush            = flush
                     , fill             = fill
                     , strip.white      = strip.white
                     , quiet            = TRUE
                     , blank.lines.skip = blank.lines.skip
                     , multi.line       = TRUE
                     , comment.char     = comment.char
                     , allowEscapes     = FALSE
                     , fileEncoding     = fileEncoding
                     , encoding         = encoding
                     , text             = text
                     )#end scan
      ans[ans == character(0L)] = NA_character_
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Trim the number of columns to the ones sought in the output.                   #
      #------------------------------------------------------------------------------------#
      nckeep = sequence(min(nc,ncols))
      ans    = ans[,nckeep]
      nc     = ncol(ans)
      #------------------------------------------------------------------------------------#


      #----- Transform matrix into data.frame, and use mapply to correct the types. -------#
      if (all(is.na(colClasses)) && length(colClasses) == 1) colClasses = rep(NA,times=nc)

      #----- Discard any column with class "NULL". ----------------------------------------#
      keep       = ! colClasses %in% "NULL"
      colClasses = colClasses[keep]
      ans        = ans[,keep]
      #------------------------------------------------------------------------------------#



      #----- Convert data for each column, then re-organise into a data frame. ------------#
      ans        = split(x=ans,f=col(ans))
      colClasses = split(x=colClasses,f=seq_along(colClasses))
      ans        = data.frame( mapply( FUN      = as.or.guess
                                     , x        = ans
                                     , Class    = colClasses
                                     , SIMPLIFY = FALSE
                                     , MoreArgs = list( na.strings = character(0L)
                                                      , as.is      = as.is
                                                      , dec        = dec
                                                      )#end list
                                     )#end mapply
                             , stringsAsFactors = FALSE
                             )#end data.frame
      names(ans) = col.names
      if (! missing(row.names)) rownames(ans) = row.names
      #------------------------------------------------------------------------------------#
   }else{
      #----- Use common read.table. -------------------------------------------------------#
      ans = read.table( file             = file
                      , header           = header
                      , sep              = sep
                      , quote            = quote
                      , comment.char     = comment.char
                      , dec              = dec
                      , row.names        = row.names
                      , col.names        = col.names
                      , as.is            = as.is
                      , na.strings       = na.strings
                      , colClasses       = colClasses
                      , nrows            = nrows
                      , skip             = skip
                      , check.names      = check.names
                      , fill             = fill
                      , strip.white      = strip.white
                      , blank.lines.skip = blank.lines.skip
                      , allowEscapes     = allowEscapes
                      , flush            = flush
                      , stringsAsFactors = stringsAsFactors
                      , fileEncoding     = fileEncoding
                      , encoding         = encoding
                      , text             = text
                      )#end read.table
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      This is just an if block really.  If Class is known, then it converts data using    #
# as, otherwise, uses type.convert to give a try.                                          #
#------------------------------------------------------------------------------------------#
as.or.guess <<- function(x,Class=NA,strict=TRUE,na.strings="NA",as.is=FALSE,dec="."){
   if (is.na(Class)){
      ans = type.convert(x=x,na.strings=na.strings,as.is=as.is,dec=dec)
   }else{
      ans = as(object=x,Class=Class,strict=strict)
   }#end if
   return(ans)
}#end function as.or.guess
#==========================================================================================#
#==========================================================================================#
