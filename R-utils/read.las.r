#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#     Functions to read in LAS files in R.                                                 #
#                                                                                          #
#     The function now uses package rlas, with a few modifications.  The rlas function     #
# readlas is preferrable because it handles LAZ files.  The output is written similarly to #
# the previous read.las for back-compatibility.                                            #
#------------------------------------------------------------------------------------------#
read.las <<- function( lasfile
                     , skip          = 0
                     , nrows         = NULL
                     , return.sp     = FALSE
                     , return.header = FALSE
                     , int.nbits     = 12L   # number of bits for intensity 
                     , select        = c("xyzitrndecskwupo")
                     ){ 


   #----- Stop if return.sp is true but package sp can't be loaded. -----------------------#
   if (return.sp && ! return.header){
      if (! "package:sp" %in% search()){
         isok = require(sp)
         if (! isok) stop("Function read.las requires package sp if return.sp = TRUE!")
      }#end if
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Check whether the las is compressed.  In case so, make a temporary file.         #
   #---------------------------------------------------------------------------------------#
   is.las = any(sapply(X=c("\\.las$"      ,"\\.laz$"      ),FUN=grepl,x=lasfile))
   is.gz  = any(sapply(X=c("\\.las\\.gz$" ,"\\.laz\\.gz$" ),FUN=grepl,x=lasfile))
   is.bz2 = any(sapply(X=c("\\.las\\.bz2$","\\.laz\\.bz2$"),FUN=grepl,x=lasfile))
   if (! file.exists(lasfile)){
      stop(paste0(" File ",lasfile," doesn't exist!"))
   }else if (is.las){
      temp.las = lasfile
   }else if (is.bz2){
      temp.las = file.path( tempdir()
                          , gsub(pattern="\\.bz2$",replacement="",x=basename(lasfile))
                          )#end file.path
      if (file.exists(temp.las)) file.remove(temp.las)
      dummy    = bunzip2(filename=lasfile,destname=temp.las,remove=FALSE)
   }else if (is.gz){
      temp.las = file.path( tempdir()
                          , gsub(pattern="\\.gz$",replacement="",x=basename(lasfile))
                          )#end file.path
      if (file.exists(temp.las)) file.remove(temp.las)
      dummy    = gunzip(filename=lasfile,destname=temp.las,remove=FALSE)
   }else{
      cat0("Unrecognised format for file ",basename(lasfile),".")
      cat0("Acceptable formats are las, laz, las.gz, laz.gz, las.bz2, or laz.bz2.")
      stop("Invalid lidar file!")
   }#end if
   #---------------------------------------------------------------------------------------#


   #------ Read the header. ---------------------------------------------------------------#
   pheader               = read.lasheader(file=temp.las)
   numberPointRecords    = pheader[["Number of point records" ]]
   offsetToPointData     = pheader[["Offset to point data"    ]]
   pointDataRecordLength = pheader[["Point Data Record Length"]]
   x.fac                 = pheader[["X scale factor"          ]]
   x.off                 = pheader[["X offset"                ]]
   y.fac                 = pheader[["Y scale factor"          ]]
   y.off                 = pheader[["Y offset"                ]]
   z.fac                 = pheader[["Z scale factor"          ]]
   z.off                 = pheader[["Z offset"                ]]
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Decide whether to return the header or the full data set.                        #
   #---------------------------------------------------------------------------------------#
   if (return.header){
      #----- Copy header to answer. -------------------------------------------------------#
      ans = pheader
      #------------------------------------------------------------------------------------#
   }else{
      #----- Make translation table. ------------------------------------------------------#
      eqname = rbind( c("X"                 ,"x"               )
                    , c("Y"                 ,"y"               )
                    , c("Z"                 ,"z"               )
                    , c("gpstime"           ,"gpstime"         )
                    , c("Intensity"         ,"intensity"       )
                    , c("ReturnNumber"      ,"retn.number"     )
                    , c("NumberOfReturns"   ,"number.retn.gp"  )
                    , c("ScanDirectionFlag" ,"scan.dir.flag"   )
                    , c("EdgeOfFlightline"  ,"edge.flight.line")
                    , c("Classification"    ,"pt.class"        )
                    , c("Synthetic_flag"    ,"synthetic"       )
                    , c("Keypoint_flag"     ,"key.point"       )
                    , c("Withheld_flag"     ,"withheld"        )
                    , c("ScanAngleRank"     ,"scan.angle.rank" )
                    , c("UserData"          ,"user.data"       )
                    , c("PointSourceID"     ,"pt.source.ID"    )
                    )#end cbind
      #------------------------------------------------------------------------------------#



      #----- Keep only the columns that are 
      ans        = rlas::read.las(files=lasfile,select=select)
      both       = intersect(eqname[,1],names(ans))
      ans        = as.data.frame(ans[,..both])
      keep       = eqname[,1] %in% both
      names(ans) = eqname[keep,2]
      #------------------------------------------------------------------------------------#



      #---- Estimate the pulse number. ----------------------------------------------------#
      if ("retn.number" %in% names(ans)){
         retn.before      = c(Inf,ans$retn.number[-nrow(ans)])
         ans$pulse.number = cumsum(ans$retn.number == 1 | ans$retn.number <= retn.before)
      }#end if ("retn.number" %in% names(ans))
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Decide the output format.  Default is a data frame, but it can be also Spatial #
      # Points (package sp).                                                               #
      #------------------------------------------------------------------------------------#
      if (return.sp) ans = SpatialPoints(ans)
      #------------------------------------------------------------------------------------#
   }#end if (return.header)
   #---------------------------------------------------------------------------------------#
   


   #----- Return answer. ------------------------------------------------------------------#
   return (ans)
   #---------------------------------------------------------------------------------------#
}#end function read.las
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     This function reads in the LAS files.                                                #
#                                                                                          #
# Originally developed by:                                                                 #
#    Michael Sumner                                                                        #
#    Antarctic Climate & Ecosystems Cooperative Research Centre                            #
#    Hobart, Tasmania, Australia                                                           #
#                                                                                          #
# Minor modifications by Marcos Longo.                                                     #
#                                                                                          #
# To do:                                                                                   #
#   - Generalise to any version                                                            #
#   - Figure out what this gpstime is to convert to POSIXct...                             #
#   - How do we get coordinate system?                                                     #
#   - Bits after Intensity                                                                 #
#   - Can we be more efficient to "seek" without actually reading bytes on windows?        #
#                                                                                          #
# Done:                                                                                    #
#   - Ensure the entire header is read, from Header Size                                   #
#   - Parse header                                                                         #
#   - Provide chunked read                                                                 #
#------------------------------------------------------------------------------------------#
old.read.las <<- function( lasfile
                         , skip          = 0
                         , nrows         = NULL
                         , return.sp     = FALSE
                         , return.header = FALSE
                         , int.nbits     = 12L   # number of bits for intensity 
                         ){ 


   #----- Stop if return.sp is true but package sp can't be loaded. -----------------------#
   if (return.sp && ! return.header){
      if (! "package:sp" %in% search()){
         isok = require(sp)
         if (! isok) stop("Function read.las requires package sp if return.sp = TRUE!")
      }#end if
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Check whether the las is compressed.  In case so, make a temporary file.         #
   #---------------------------------------------------------------------------------------#
   is.las = grepl(pattern="\\.las$"      ,x=lasfile)
   is.gz  = grepl(pattern="\\.las\\.gz$" ,x=lasfile)
   is.bz2 = grepl(pattern="\\.las\\.bz2$",x=lasfile)
   if (! file.exists(lasfile)){
      stop(paste0(" File ",lasfile," doesn't exist!"))
   }else if (is.las){
      temp.las = lasfile
   }else if (is.bz2){
      temp.las = file.path( tempdir()
                          , gsub(pattern="\\.bz2$",replacement="",x=basename(lasfile))
                          )#end file.path
      if (file.exists(temp.las)) file.remove(temp.las)
      dummy    = bunzip2(filename=lasfile,destname=temp.las,remove=FALSE)
   }else if (is.gz){
      temp.las = file.path( tempdir()
                          , gsub(pattern="\\.gz$",replacement="",x=basename(lasfile))
                          )#end file.path
      if (file.exists(temp.las)) file.remove(temp.las)
      dummy    = gunzip(filename=lasfile,destname=temp.las,remove=FALSE)
   }else{
      stop("Unrecognised format: point cloud must be las, las.gz, or las.bz2!")
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Get the header. -----------------------------------------------------------------#
   hd             = las.header()
   nhd            = nrow(hd)
   pheader        = vector("list", nhd)
   names(pheader) = hd$Item
   #---------------------------------------------------------------------------------------#

   #---- Open connection, and read the LAS File bytes. ------------------------------------#
   con                   = file( description = temp.las, open = "rb")
   isLASFbytes           = readBin( con    = con
                                  , what   = "raw"
                                  , size   = 1
                                  , n      = 4L
                                  , endian = "little"
                                  )#end readBin
   #---------------------------------------------------------------------------------------#

   #---- Read first header item, and make sure that it understands the LASF header. -------#
   pheader[[hd$Item[1]]] = readBin( con    = isLASFbytes
                                  , what   = "character"
                                  , size   = 4
                                  , n      = 1L
                                  , endian = "little"
                                  )#end readBin
   if (! pheader[[hd$Item[1]]] %in% "LASF") {
      stop(paste("File ",basename(temp.las)," is not a valid LAS file!",sep=""))
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Go through the additional header items. -----------------------------------------#
   warn.orig = getOption("warn")
   options(warn=-1)
   for (i in sequence(nhd)[-1]){
      now            = hd$Item[i]
      pheader[[now]] = readBin( con    = con
                              , what   = hd$what  [i]
                              , signed = hd$signed[i]
                              , size   = hd$Rsize [i]
                              , n      = hd$n     [i]
                              , endian = "little"
                              )#end readBin
   }#end for (i in sequence(nhd)[-1])
   options(warn=warn.orig)
   #---------------------------------------------------------------------------------------#



   #----- Close the file. -----------------------------------------------------------------#
   close(con=con)
   #---------------------------------------------------------------------------------------#




   #----- Read the data. ------------------------------------------------------------------#
   numberPointRecords    = pheader[["Number of point records" ]]
   offsetToPointData     = pheader[["Offset to point data"    ]]
   pointDataRecordLength = pheader[["Point Data Record Length"]]
   x.fac                 = pheader[["X scale factor"          ]]
   x.off                 = pheader[["X offset"                ]]
   y.fac                 = pheader[["Y scale factor"          ]]
   y.off                 = pheader[["Y offset"                ]]
   z.fac                 = pheader[["Z scale factor"          ]]
   z.off                 = pheader[["Z offset"                ]]
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Check what is to be returned.                                                     #
   #---------------------------------------------------------------------------------------#
   if (return.header){
      #----- Return the header if the user only wants to know the header. -----------------#
      ans = pheader
      #------------------------------------------------------------------------------------#
   }else{ 

      #------------------------------------------------------------------------------------#
      #     Return the actual data.                                                        #
      #------------------------------------------------------------------------------------#


      #----- Read in the actual data. -----------------------------------------------------#
      con  = file   (description=temp.las, open = "rb")
      junk = readBin(con=con,what="raw",size=1L,n=offsetToPointData)
      #------------------------------------------------------------------------------------#



      #----- Deal with rows to skip. ------------------------------------------------------#
      if (skip > 0) {

         #----- Junk the bytes to skip. ---------------------------------------------------#
         junk               = readBin( con  = con
                                     , what = "raw"
                                     , size = 1L
                                     , n    = pointDataRecordLength * skip
                                     )#end readBin
         numberPointRecords = numberPointRecords - skip
         #---------------------------------------------------------------------------------#
      }#end if (skip > 0)
      #------------------------------------------------------------------------------------#




      #---- Deal with maximum number of rows to be read. ----------------------------------#
      if (! is.null(nrows)) {
         if (numberPointRecords > nrows) numberPointRecords = nrows
      }#end if (! is.null(nrows))
      #------------------------------------------------------------------------------------#



      #----- Check that there are points left to be read. ---------------------------------#
      if (numberPointRecords < 1) stop("No records left to read!")
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Read the data.                                                                 #
      #------------------------------------------------------------------------------------#
      allbytes = readBin( con    = con
                        , what   = "raw"
                        , n      = pointDataRecordLength * numberPointRecords
                        , size   = 1L
                        , endian = "little"
                        )#end readBin
      allbytes = matrix( data  = allbytes
                       , ncol  = pointDataRecordLength
                       , nrow  = numberPointRecords
                       , byrow = TRUE
                       )#end matrix
      close(con = con)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Initialise the output.                                                         #
      #------------------------------------------------------------------------------------#
      empty = rep(x=NA,times=numberPointRecords)
      ans   = data.frame( x                = empty
                        , y                = empty
                        , z                = empty
                        , intensity        = empty
                        , pulse.number     = empty
                        , retn.number      = empty
                        , number.retn.gp   = empty
                        , scan.dir.flag    = empty
                        , edge.flight.line = empty
                        , pt.class         = empty
                        , synthetic        = empty
                        , key.point        = empty
                        , withheld         = empty
                        , scan.angle.rank  = empty
                        , user.data        = empty
                        , pt.source.ID     = empty
                        , gpstime          = empty
                        )#end data.frame
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Read the X coordinates.                                                        #
      #------------------------------------------------------------------------------------#
      ans$x = readBin( t(allbytes[,1:4])
                     , what   = "integer"
                     , size   = 4L
                     , n      = numberPointRecords
                     , endian = "little"
                     )#end readBin
      ans$x = ans$x * x.fac + x.off
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Read the Y coordinates.                                                        #
      #------------------------------------------------------------------------------------#
      ans$y = readBin( t(allbytes[,5:8])
                     , what   = "integer"
                     , size   = 4L
                     , n      = numberPointRecords
                     , endian = "little"
                     )#end readBin
      ans$y = ans$y * y.fac + y.off
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Read the Z coordinates.                                                        #
      #------------------------------------------------------------------------------------#
      ans$z = readBin( t(allbytes[,9:12])
                     , what   = "integer"
                     , size   = 4L
                     , n      = numberPointRecords
                     , endian = "little"
                     )#end readBin
      ans$z = ans$z * z.fac + z.off
      #------------------------------------------------------------------------------------#




      #----- Read intensity. --------------------------------------------------------------#
      ans$intensity = readBin( t(allbytes[, 13:14])
                             , what   = "integer"
                             , size   = 2L
                             , n      = numberPointRecords
                             , signed = FALSE
                             , endian = "little"
                             )#end readBin
      ans$intensity = ans$intensity * 2^(16L-int.nbits)
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #      Read return byte.                                                             #
      #------------------------------------------------------------------------------------#
      retn.byte              = readBin( t(allbytes[, 15])
                                      , what   = "integer"
                                      , size   = 1L
                                      , n      = numberPointRecords
                                      , signed = FALSE
                                      , endian = "little"
                                      )#end readBin
      #----- Convert retn.byte to raw. ----------------------------------------------------#
      retn.byte              = matrix(intToBits(retn.byte),nrow=32,ncol=numberPointRecords)
      retn.number            = matrix(raw(),nrow=8,ncol=numberPointRecords)
      number.retn.gp         = matrix(raw(),nrow=8,ncol=numberPointRecords)
      scan.dir.flag          = matrix(raw(),nrow=8,ncol=numberPointRecords)
      edge.flight.line       = matrix(raw(),nrow=8,ncol=numberPointRecords)
      #----- Copy bits. -------------------------------------------------------------------#
      retn.number     [1:3,] = retn.byte[1:3,]
      number.retn.gp  [1:3,] = retn.byte[4:6,]
      scan.dir.flag   [  1,] = retn.byte[7  ,]
      edge.flight.line[  1,] = retn.byte[8  ,]
      #----- Convert to numeric or logical. -----------------------------------------------#
      ans$retn.number        = as.numeric(packBits(unlist(retn.number     )))
      ans$number.retn.gp     = as.numeric(packBits(unlist(number.retn.gp  )))
      ans$scan.dir.flag      = as.logical(packBits(unlist(scan.dir.flag   )))
      ans$edge.flight.line   = as.logical(packBits(unlist(edge.flight.line)))
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #      Read point classification byte.                                               #
      #------------------------------------------------------------------------------------#
      class.byte = readBin( t(allbytes[, 16])
                          , what   = "integer"
                          , size   = 1L
                          , n      = numberPointRecords
                          , signed = FALSE
                          , endian = "little"
                          )#end readBin
      #----- Convert class.byte to raw. ---------------------------------------------------#
      class.byte       = matrix(intToBits(class.byte),nrow=32,ncol=numberPointRecords)
      pt.class         = matrix(raw(),nrow=8,ncol=numberPointRecords)
      synthetic        = matrix(raw(),nrow=8,ncol=numberPointRecords)
      key.point        = matrix(raw(),nrow=8,ncol=numberPointRecords)
      withheld         = matrix(raw(),nrow=8,ncol=numberPointRecords)
      #----- Copy bits, and convert back to integer. --------------------------------------#
      pt.class  [1:5,] = class.byte[1:5,]
      synthetic [  1,] = class.byte[  6,]
      key.point [  1,] = class.byte[  7,]
      withheld  [  1,] = class.byte[  8,]
      ans$pt.class     = as.numeric(packBits(unlist(pt.class      )))
      ans$synthetic    = as.logical(packBits(unlist(synthetic     )))
      ans$key.point    = as.numeric(packBits(unlist(key.point     )))
      ans$withheld     = as.logical(packBits(unlist(withheld      )))
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #      Read the scan angle rank (left side).                                         #
      #------------------------------------------------------------------------------------#
      ans$scan.angle.rank = readBin( t(allbytes[, 17])
                                   , what   = "integer"
                                   , size   = 1L
                                   , n      = numberPointRecords
                                   , signed = TRUE
                                   , endian = "little"
                                   )#end readBin
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #      Read the scan angle rank (left side).                                         #
      #------------------------------------------------------------------------------------#
      ans$user.data = readBin( t(allbytes[, 18])
                             , what   = "integer"
                             , size   = 1L
                             , n      = numberPointRecords
                             , signed = FALSE
                             , endian = "little"
                             )#end readBin
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #      Read the point source ID.                                                     #
      #------------------------------------------------------------------------------------#
      ans$pt.source.ID = readBin( t(allbytes[, 19:20])
                                , what   = "integer"
                                , size   = 2L
                                , n      = numberPointRecords
                                , signed = FALSE
                                , endian = "little"
                                )#end readBin
      #------------------------------------------------------------------------------------#





      #----- Read the GPS time (if available). --------------------------------------------#
      if (ncol(allbytes) >= 28){
         ans$gpstime = readBin( t(allbytes[ , 21:28])
                              , what   = "numeric"
                              , size   = 8L
                              , n      = numberPointRecords
                              , endian = "little"
                              )#end readBin
      }#end if
      #------------------------------------------------------------------------------------#



      #---- Estimate the pulse number. ----------------------------------------------------#
      retn.before      = c(Inf,ans$retn.number[-nrow(ans)])
      ans$pulse.number = cumsum(ans$retn.number == 1 | ans$retn.number <= retn.before)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Decide the output format.  Default is a data frame, but it can be also Spatial #
      # Points (package sp).                                                               #
      #------------------------------------------------------------------------------------#
      if (return.sp) ans = SpatialPoints(ans)
      #------------------------------------------------------------------------------------#
   }#end if (return.header)
   #---------------------------------------------------------------------------------------#

   #----- Remove temporary file. ----------------------------------------------------------#
   if (is.gz || is.bz2) file.remove(temp.las)
   #---------------------------------------------------------------------------------------#


   #----- Return answer. ------------------------------------------------------------------#
   return (ans)
   #---------------------------------------------------------------------------------------#
}#end function read.las
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function sets the LAS header.                                                   #
#------------------------------------------------------------------------------------------#
las.header <<- function() {
   #----- List with variable information. -------------------------------------------------#
   n       = 0
   hd      = list()
   n       = n + 1
   hd[[n]] = list( Item     = "File Signature (``LASF'')"
                 , Format   = "char[4]"
                 , Size     = "4 bytes"
                 , Required = "*"
                 )#end list
   n       = n + 1
   hd[[n]] = list( Item     = "(1.1) File Source ID"
                 , Format   = "unsigned short"
                 , Size     = "2 bytes"
                 , Required = "*"
                 )#end list
   n       = n + 1
   hd[[n]] = list( Item     = "(1.1) Global Encoding"
                 , Format   = "unsigned short"
                 , Size     = "2 bytes"
                 , Required = "*"
                 )#end list
   n       = n + 1
   hd[[n]] = list( Item     = "(1.1) Project ID - GUID data 1"
                 , Format   = "unsigned long"
                 , Size     = "4 bytes"
                 , Required = "*"
                 )#end list
   n       = n + 1
   hd[[n]] = list( Item     = "(1.1) Project ID - GUID data 2"
                 , Format   = "unsigned short"
                 , Size     = "2 bytes"
                 , Required = "*"
                 )#end list
   n       = n + 1
   hd[[n]] = list( Item     = "(1.1) Project ID - GUID data 3"
                 , Format   = "unsigned short"
                 , Size     = "2 bytes"
                 , Required = "*"
                 )#end list
   n       = n + 1
   hd[[n]] = list( Item     = "(1.1) Project ID - GUID data 4"
                 , Format   = "unsigned char[8]"
                 , Size     = "8 bytes"
                 , Required = "*"
                 )#end list
   n       = n + 1
   hd[[n]] = list( Item     = "Version Major"
                 , Format   = "unsigned char"
                 , Size     = "1 byte"
                 , Required = "*"
                 )#end list
   n       = n + 1
   hd[[n]] = list( Item     = "Version Minor"
                 , Format   = "unsigned char"
                 , Size     = "1 byte"
                 , Required = "*"
                 )#end list
   n       = n + 1
   hd[[n]] = list( Item     = "(1.1) System Identifier"
                 , Format   = "char[32]"
                 , Size     = "32 bytes"
                 , Required = "*"
                 )#end list
   n       = n + 1
   hd[[n]] = list( Item     = "Generating Software"
                 , Format   = "char[32]"
                 , Size     = "32 bytes"
                 , Required = "*"
                 )#end list
   n       = n + 1
   hd[[n]] = list( Item     = "(1.1) File Creation Day of Year"
                 , Format   = "unsigned short"
                 , Size     = "2 bytes"
                 , Required = "*"
                 )#end list
   n       = n + 1
   hd[[n]] = list( Item     = "(1.1) File Creation Year"
                 , Format   = "unsigned short"
                 , Size     = "2 bytes"
                 , Required = "*"
                 )#end list
   n       = n + 1
   hd[[n]] = list( Item     = "Header Size"
                 , Format   = "unsigned short"
                 , Size     = "2 bytes"
                 , Required = "*"
                 )#end list
   n       = n + 1
   hd[[n]] = list( Item     = "Offset to point data"
                 , Format   = "unsigned long"
                 , Size     = "4 bytes"
                 , Required = "*"
                 )#end list
   n       = n + 1
   hd[[n]] = list( Item     = "Number of variable length records"
                 , Format   = "unsigned long"
                 , Size     = "4 bytes"
                 , Required = "*"
                 )#end list
   n       = n + 1
   hd[[n]] = list( Item     = "Point Data Format ID (0-99 for spec)"
                 , Format   = "unsigned char"
                 , Size     = "1 byte"
                 , Required = "*"
                 )#end list
   n       = n + 1
   hd[[n]] = list( Item     = "Point Data Record Length"
                 , Format   = "unsigned short"
                 , Size     = "2 bytes"
                 , Required = "*"
                 )#end list
   n       = n + 1
   hd[[n]] = list( Item     = "Number of point records"
                 , Format   = "unsigned long"
                 , Size     = "4 bytes"
                 , Required = "*"
                 )#end list
   n       = n + 1
   hd[[n]] = list( Item     = "Number of points by return"
                 , Format   = "unsigned long[5]"
                 , Size     = "20 bytes"
                 , Required = "*"
                 )#end list
   n       = n + 1
   hd[[n]] = list( Item     = "X scale factor"
                 , Format   = "double"
                 , Size     = "8 bytes"
                 , Required = "*"
                 )#end list
   n       = n + 1
   hd[[n]] = list( Item     = "Y scale factor"
                 , Format   = "double"
                 , Size     = "8 bytes"
                 , Required = "*"
                 )#end list
   n       = n + 1
   hd[[n]] = list( Item     = "Z scale factor"
                 , Format   = "double"
                 , Size     = "8 bytes"
                 , Required = "*"
                 )#end list
   n       = n + 1
   hd[[n]] = list( Item     = "X offset"
                 , Format   = "double"
                 , Size     = "8 bytes"
                 , Required = "*"
                 )#end list
   n       = n + 1
   hd[[n]] = list( Item     = "Y offset"
                 , Format   = "double"
                 , Size     = "8 bytes"
                 , Required = "*"
                 )#end list
   n       = n + 1
   hd[[n]] = list( Item     = "Z offset"
                 , Format   = "double"
                 , Size     = "8 bytes"
                 , Required = "*"
                 )#end list
   n       = n + 1
   hd[[n]] = list( Item     = "Max X"
                 , Format   = "double"
                 , Size     = "8 bytes"
                 , Required = "*"
                 )#end list
   n       = n + 1
   hd[[n]] = list( Item     = "Min X"
                 , Format   = "double"
                 , Size     = "8 bytes"
                 , Required = "*"
                 )#end list
   n       = n + 1
   hd[[n]] = list( Item     = "Max Y"
                 , Format   = "double"
                 , Size     = "8 bytes"
                 , Required = "*"
                 )#end list
   n       = n + 1
   hd[[n]] = list( Item     = "Min Y"
                 , Format   = "double"
                 , Size     = "8 bytes"
                 , Required = "*"
                 )#end list
   n       = n + 1
   hd[[n]] = list( Item     = "Max Z"
                 , Format   = "double"
                 , Size     = "8 bytes"
                 , Required = "*"
                 )#end list
   n       = n + 1
   hd[[n]] = list( Item     = "Min Z"
                 , Format   = "double"
                 , Size     = "8 bytes"
                 , Required = "*"
                 )#end list
   #---------------------------------------------------------------------------------------#


   #----- Convert list to data frame. -----------------------------------------------------#
   hd           = list.2.data.frame(hd)
   nhd          = nrow(hd)
   rownames(hd) = sequence(nhd)+1
   #---------------------------------------------------------------------------------------#


   #----- Initialise the variable type. ---------------------------------------------------#
   hd$what = character(length=nhd)
   hd$what[grep("unsigned", hd$Format)] = "integer"
   hd$what[grep("char"    , hd$Format)] = "raw"
   hd$what[grep("short"   , hd$Format)] = "integer"
   hd$what[grep("long"    , hd$Format)] = "integer"
   hd$what[grep("double"  , hd$Format)] = "numeric"
   #---------------------------------------------------------------------------------------#


   #----- Set Boolean flag for signed variable. -------------------------------------------#
   hd$signed  = ! grepl("unsigned",hd$Format)
   #---------------------------------------------------------------------------------------#


   #----- Number of values in record. -----------------------------------------------------#
   hd$n                         = as.numeric(gsub("[[:alpha:][:punct:]]", "", hd$Format))
   hd$n[hd$what == "character"] = 1
   hd$n[is.na(hd$n)]            = 1
   #---------------------------------------------------------------------------------------#


   #----- Size of record. -----------------------------------------------------------------#
   hd$Hsize = as.numeric(gsub("[[:alpha:]]", "", hd$Size))
   #---------------------------------------------------------------------------------------#


   #----- Size of each value in record. ---------------------------------------------------#
   hd$Rsize         = hd$Hsize / hd$n
   is.raw           = hd$what %in% "raw"
   hd$Rsize[is.raw] = 1
   hd$n    [is.raw] = hd$Hsize[is.raw]
   #---------------------------------------------------------------------------------------#


   return(hd)
}#end function las.header
#==========================================================================================#
#==========================================================================================#

