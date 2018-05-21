#==========================================================================================#
#==========================================================================================#
#    General function to open plots.                                                       #
#------------------------------------------------------------------------------------------#
open.plot <<- function( fichier
                      , outform  = switch( EXPR    = get.os()
                                         , macos   = "quartz"
                                         , linux   = "x11"
                                         , windows = "windows"
                                         , "x11"
                                         )#end switch
                      , size     = plotsize(proje=FALSE,paper="square")
                      , ptsz     = 16
                      , depth    = 300
                      , bg       = "transparent"
                      , compress = "lzw"
                      , ...
                      ){

   #----- Fichier must be provided unless we are plotting on screen. ----------------------#
   if (missing(fichier) && (! outform %in% c("x11","quartz","windows"))){
      stop("Output file (\"fichier\") must be provided when plotting on file!")
   }#end if (missing(fichier) && (! outform %in% c("x11","quartz","windows")))
   #---------------------------------------------------------------------------------------#



   #----- Open the file or the plot window. -----------------------------------------------#
   if (outform %in% "x11"){
      X11       ( width     = size$width
                , height    = size$height
                , pointsize = ptsz
                , ...
                )#end X11
   }else if (outform %in% "quartz"){
      quartz    ( width     = size$width
                , height    = size$height
                , pointsize = ptsz
                , ...
                )#end quartz
   }else if (outform %in% "windows"){
      windows   ( width     = size$width
                , height    = size$height
                , pointsize = ptsz
                , ...
                )#end quartz
   }else if (outform %in% "png"){
      png       ( filename  = fichier
                , width     = size$width*depth
                , height    = size$height*depth
                , pointsize = ptsz
                , res       = depth
                , bg        = bg
                , ...
                )#end png
   }else if (outform %in% "tif"){
      tiff      ( filename  = fichier
                , width     = size$width*depth
                , height    = size$height*depth
                , pointsize = ptsz
                , res       = depth
                , bg        = bg
                , compress  = compress
                , ...
                )#end tiff
   }else if (outform %in% "eps"){
      postscript( file      = fichier
                , width     = size$width
                , height    = size$height
                , pointsize = ptsz
                , paper     = size$paper
                , ...
                )#end postscript
   }else if (outform[o] %in% "pdf"){
      pdf       ( file      = fichier
                , onefile   = FALSE
                , width     = size$width
                , height    = size$height
                , pointsize = ptsz
                , paper     = size$paper
                )#end pdf
   }#end if
   #---------------------------------------------------------------------------------------#

   invisible()
}#end open.plot
#==========================================================================================#
#==========================================================================================#



#==========================================================================================#
#==========================================================================================#
#    General function to close plots.                                                      #
#------------------------------------------------------------------------------------------#
close.plot <<- function( outform = switch( EXPR    = get.os()
                                         , macos   = "quartz"
                                         , linux   = "x11"
                                         , windows = "windows"
                                         , "x11"
                                         )#end switch
                       , n       = 1
                       ){

   if (outform %in% c("x11","quartz","windows")){
      locator(n=n)
      dev.off()
   }else{
      dev.off()
   }#end if
   dummy = clean.tmp()


   invisible()
}#end close.plot
#==========================================================================================#
#==========================================================================================#
