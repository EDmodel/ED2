#==========================================================================================#
#==========================================================================================#
#  This function plots the map for cities.                                                 #
#------------------------------------------------------------------------------------------#
citymap <<- function( mapdir   = NULL
                    , cities   = NULL
                    , skip.off = TRUE
                    , xlim     = NULL
                    , ylim     = NULL
                    , fun      = c("lines","points","polygon","panel.points")
                    , ...){

   #----- Decide which type we should use. ------------------------------------------------#
   fun = match.arg(fun)
   #---------------------------------------------------------------------------------------#


   #----- Make sure package lattice is properly loaded in case of panel.points. -----------#
   if (fun %in% "panel.points"){
      isok = require("lattice")
      if (! isok){
         stop("Package 'lattice' must be loaded if using panel.points to draw the map!")
      }#end if
   }#end if
   #---------------------------------------------------------------------------------------#


   #---- Get the limits. ------------------------------------------------------------------#
   if (is.null(xlim)) xlim = par("usr")[c(1,2)]
   if (is.null(ylim)) ylim = par("usr")[c(3,4)]
   #---------------------------------------------------------------------------------------#


   #----- Make the default path to search in case none is given. --------------------------#
   if (is.null(mapdir)) mapdir = file.path(srcdir,"citymap")
   #---------------------------------------------------------------------------------------#


   #----- Figure all the countries and country subdivisions available. --------------------#
   all.files       = dir(mapdir)
   mapped.cities   = substring(basename(all.files),6,8)
   n.mapped.cities = length(mapped.cities)
   #---------------------------------------------------------------------------------------#



   #----- Make the list with all cities to plot if it isn't given. ------------------------#
   if (is.null(cities)){
      cities    = mapped.cities
   }else{
      #----- Match cities with standardised "iata" code. ----------------------------------#
      use.short = poilist$short    %in% tolower(cities)
      use.full  = poilist$longname %in% tolower(cities)
      use.iata  = poilist$iata     %in% tolower(cities)
      use       = use.short | use.full | use.iata
      cities    = poilist$iata[use]
      #------------------------------------------------------------------------------------#

      #----- Keep only the cities that have been mapped. ----------------------------------#
      cities    = cities[cities %in% mapped.cities]
      #------------------------------------------------------------------------------------#
   }#end if
   ct.file  = file.path(mapdir,paste("city_",cities,".txt",sep=""))
   n.cities = length(cities)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Grab options for plot.                                                             #
   #---------------------------------------------------------------------------------------#
   dotdotdot = list(...)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Loop over all cities to plot.                                                      #
   #---------------------------------------------------------------------------------------#
   for (ct in sequence(n.cities)){
      dat        = read.table(ct.file[ct],na.string=c("NA","999"))
      names(dat) = c("x","y")
      #------------------------------------------------------------------------------------#
      #     Check whether any of the point is within range.  Plot the entire group if it   #
      # is, but replace outside numbers by NA.                                             #
      #------------------------------------------------------------------------------------#
      if ((! skip.off) || any(dat$x %wr% xlim & dat$y %wr% ylim)){
         if (! fun %in% "polygon"){
            keep        = dat$x %wr% xlim & dat$y %wr% ylim
            dat[!keep,] = NA
         }#end if

         #----- Append X and Y to the argument list and call the function. ----------------#
         argnow      = modifyList(x=dotdotdot,val=list(x=dat$x,y=dat$y))
         dummy       = do.call(what=fun,args=argnow)
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#  This function plots the map for cities in 3D plots.                                     #
#------------------------------------------------------------------------------------------#
citymap3d <<- function( height      = 0
                      , pmat
                      , mapdir   = NULL
                      , cities   = NULL
                      , skip.off = TRUE
                      , xlim     = NULL
                      , ylim     = NULL
                      , ...){


   #---- Get the limits. ------------------------------------------------------------------#
   if (is.null(xlim)) xlim = par("usr")[c(1,2)]
   if (is.null(ylim)) ylim = par("usr")[c(3,4)]
   #---------------------------------------------------------------------------------------#


   #----- Make the default path to search in case none is given. --------------------------#
   if (is.null(mapdir)) mapdir = file.path(srcdir,"citymap")
   #---------------------------------------------------------------------------------------#


   #----- Figure all the countries and country subdivisions available. --------------------#
   all.files       = dir(mapdir)
   mapped.cities   = substring(basename(all.files),6,8)
   n.mapped.cities = length(mapped.cities)
   #---------------------------------------------------------------------------------------#



   #----- Make the list with all cities to plot if it isn't given. ------------------------#
   if (is.null(cities)){
      cities    = mapped.cities
   }else{
      #----- Match cities with standardised "iata" code. ----------------------------------#
      use.short = poilist$short    %in% tolower(cities)
      use.full  = poilist$longname %in% tolower(cities)
      use.iata  = poilist$iata     %in% tolower(cities)
      use       = use.short | use.full | use.iata
      cities    = poilist$iata[use]
      #------------------------------------------------------------------------------------#

      #----- Keep only the cities that have been mapped. ----------------------------------#
      cities    = cities[cities %in% mapped.cities]
      #------------------------------------------------------------------------------------#
   }#end if
   ct.file  = file.path(mapdir,paste("city_",cities,".txt",sep=""))
   n.cities = length(cities)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Loop over all cities to plot.                                                      #
   #---------------------------------------------------------------------------------------#
   for (ct in sequence(n.cities)){
      dat        = read.table(ct.file[ct],na.string=c("NA","999"))
      names(dat) = c("x","y")
      #------------------------------------------------------------------------------------#
      #     Check whether any of the point is within range.  Plot the entire group if it   #
      # is, but replace outside numbers by NA.                                             #
      #------------------------------------------------------------------------------------#
      if ((! skip.off) || any(dat$x %wr% xlim & dat$y %wr% ylim)){
         keep        = dat$x %wr% xlim & dat$y %wr% ylim
         dat[!keep,] = NA
         dummy       = lines(trans3d(x=dat$x,y=dat$y,z=height,pmat),...)
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#
}#end function citymap3d
#==========================================================================================#
#==========================================================================================#
