#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#  This function plots the map for South America, including Brazilian states.              #
#------------------------------------------------------------------------------------------#
southammap <<- function( mapdir         = NULL
                       , countries      = NULL
                       , states         = TRUE
                       , skip.countries = FALSE
                       , skip.states    = FALSE
                       , skip.off       = TRUE
                       , xlim           = NULL
                       , ylim           = NULL
                       , fun            = c("lines","points","polygon","panel.points")
                       ,...){


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
   if (is.null(mapdir)) mapdir = file.path(srcdir,"samap")
   #---------------------------------------------------------------------------------------#



   #----- Figure all the countries and country subdivisions available. --------------------#
   all.files = dir(mapdir)
   n.files   = length(all.files)
   idx.ct    = grep("_zz.txt",all.files)
   idx.dv    = sequence(n.files)[-idx.ct]
   #---------------------------------------------------------------------------------------#



   #----- Check which countries to plot. --------------------------------------------------#
   if (is.null(countries) && is.null(states)){
      countries   = tolower(country.list$iata)
      st.only     = FALSE
      def.country = TRUE
   }else if (is.null(countries) && (! is.logical(states))){
      countries   = character(0)
      st.only     = TRUE
      def.country = TRUE
   }else{
      if (is.null(countries)){
         countries   = tolower(country.list$iata)
         def.country = TRUE
      }else{
         countries   = tolower(countries)
         def.country = FALSE
      }#end if
      st.only    = FALSE
      use.iata   = country.list$iata %in% countries
      use.full   = country.list$full %in% countries
      use.alte   = country.list$alte %in% countries
      use        = use.iata | use.full | use.alte
      countries  = country.list$iata[use]
   }#end if
   ct.file     = file.path(mapdir,paste(countries,"_zz.txt",sep=""))
   n.countries = length(countries)
   #---------------------------------------------------------------------------------------#



   #----- Check which states to plot. -----------------------------------------------------#
   if (is.null(states)){
      use       = state.list$country %in% countries & state.list$std
      states    = state.list$iata   [use]
      ct.states = state.list$country[use]
   }else if (is.logical(states)){
      if (states && def.country){
         use       = ( state.list$country %in% countries & state.list$std ) | st.only
         states    = state.list$iata   [use]
         ct.states = state.list$country[use]
      }else if (states){
         use       = state.list$country %in% countries | st.only
         states    = state.list$iata   [use]
         ct.states = state.list$country[use]
      }else{
         states    = character(0)
         ct.states = character(0)
      }#end if
   }else{
      states    = tolower(states)
      ct.use    = state.list$country %in% countries | st.only
      st.use    = ( state.list$iata %in% states 
                  | state.list$full %in% states )
      use       = ct.use & st.use
      states    = state.list$iata   [use]
      ct.states = state.list$country[use]
   }#end if
   st.file  = file.path(mapdir,paste(ct.states,"_",states,".txt",sep=""))
   n.states = length(states)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Grab options for plot.                                                             #
   #---------------------------------------------------------------------------------------#
   dotdotdot = list(...)
   #---------------------------------------------------------------------------------------#


   #----- Set loops for countries and states. ---------------------------------------------#
   if (skip.countries){
      loop.countries = numeric(0)
   }else{
      loop.countries = sequence(n.countries)
   }#end if (skip.countries)
   if (skip.states){
      loop.states    = numeric(0)
   }else{
      loop.states    = sequence(n.states)
   }#end if (skip.states)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Loop over all countries.                                                           #
   #---------------------------------------------------------------------------------------#
   for (ct in loop.countries){
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



   #---------------------------------------------------------------------------------------#
   #    Loop over all states.                                                              #
   #---------------------------------------------------------------------------------------#
   for (st in loop.states){
      dat        = read.table(st.file[st],na.string=c("NA","999"))
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
#------------------------------------------------------------------------------------------#





#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#  This function plots the map for South America, including Brazilian states.              #
#------------------------------------------------------------------------------------------#
southammap3d <<- function( height         = 0
                         , pmat
                         , mapdir         = NULL
                         , countries      = NULL
                         , states         = TRUE
                         , skip.countries = FALSE
                         , skip.states    = FALSE
                         , skip.off       = TRUE
                         , xlim           = NULL
                         , ylim           = NULL
                         ,...){
   #---- Get the limits. ------------------------------------------------------------------#
   if (is.null(xlim)) xlim = attr(pmat,"ranges")$x
   if (is.null(ylim)) ylim = attr(pmat,"ranges")$y
   #---------------------------------------------------------------------------------------#


   #----- Make the default path to search in case none is given. --------------------------#
   if (is.null(mapdir)) mapdir = file.path(srcdir,"samap")
   #---------------------------------------------------------------------------------------#



   #----- Figure all the countries and country subdivisions available. --------------------#
   all.files = dir(mapdir)
   n.files   = length(all.files)
   idx.ct    = grep("_zz.txt",all.files)
   idx.dv    = sequence(n.files)[-idx.ct]
   #---------------------------------------------------------------------------------------#



   #----- Check which countries to plot. --------------------------------------------------#
   if (is.null(countries) && is.null(states)){
      countries  = tolower(country.list$iata)
      st.only    = FALSE
   }else if (is.null(countries) && (! is.logical(states))){
      countries  = character(0)
      st.only    = TRUE
   }else{
      if (is.null(countries)){
         countries  = tolower(country.list$iata)
      }else{
         countries  = tolower(countries)
      }#end if
      st.only    = FALSE
      use.iata   = country.list$iata %in% countries
      use.full   = country.list$full %in% countries
      use.alte   = country.list$alte %in% countries
      use        = use.iata | use.full | use.alte
      countries  = country.list$iata[use]
   }#end if
   ct.file     = file.path(mapdir,paste(countries,"_zz.txt",sep=""))
   n.countries = length(countries)
   #---------------------------------------------------------------------------------------#



   #----- Check which states to plot. -----------------------------------------------------#
   if (is.null(states)){
      use       = state.list$country %in% countries & state.list$std
      states    = state.list$iata   [use]
      ct.states = state.list$country[use]
   }else if (is.logical(states)){
      if (states){
         use       = ( state.list$country %in% countries & state.list$std ) | st.only
         states    = state.list$iata   [use]
         ct.states = state.list$country[use]
      }else{
         states    = character(0)
         ct.states = character(0)
      }#end if
   }else{
      states    = tolower(states)
      ct.use    = state.list$country %in% countries | st.only
      st.use    = ( state.list$iata %in% states 
                  | state.list$full %in% states )
      use       = ct.use & st.use
      states    = state.list$iata   [use]
      ct.states = state.list$country[use]
   }#end if
   st.file  = file.path(mapdir,paste(ct.states,"_",states,".txt",sep=""))
   n.states = length(states)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Grab options for plot.                                                             #
   #---------------------------------------------------------------------------------------#
   dotdotdot = list(...)
   #---------------------------------------------------------------------------------------#


   #----- Set loops for countries and states. ---------------------------------------------#
   if (skip.countries){
      loop.countries = numeric(0)
   }else{
      loop.countries = sequence(n.countries)
   }#end if (skip.countries)
   if (skip.states){
      loop.states    = numeric(0)
   }else{
      loop.states    = sequence(n.states)
   }#end if (skip.states)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Loop over all countries.                                                           #
   #---------------------------------------------------------------------------------------#
   for (ct in loop.countries){
      dat        = read.table(ct.file[ct],na.string=c("NA","999"))
      names(dat) = c("x","y")
      #------------------------------------------------------------------------------------#
      #     Check whether any of the point is within range.  Plot the entire group if it   #
      # is, but replace outside numbers by NA.                                             #
      #------------------------------------------------------------------------------------#
      if (skip.off){
         keep        = dat$x %wr% xlim & dat$y %wr% ylim
         dat[!keep,] = NA
      }#end if
      #------------------------------------------------------------------------------------#

      #----- Append X and Y to the argument list and call the function. -------------------#
      u3d    = trans3d(x=dat$x,y=dat$y,z=height,pmat=pmat)
      argnow = modifyList(x=dotdotdot,val=list(x=u3d$x,y=u3d$y))
      dummy  = do.call(what="lines",args=argnow)
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Loop over all states.                                                              #
   #---------------------------------------------------------------------------------------#
   for (st in loop.states){
      dat        = read.table(st.file[st],na.string=c("NA","999"))
      names(dat) = c("x","y")
      #------------------------------------------------------------------------------------#
      #     Check whether any of the point is within range.  Plot the entire group if it   #
      # is, but replace outside numbers by NA.                                             #
      #------------------------------------------------------------------------------------#
      if (skip.off){
         keep        = dat$x %wr% xlim & dat$y %wr% ylim
         dat[!keep,] = NA
      }#end if
      #------------------------------------------------------------------------------------#

      #----- Append X and Y to the argument list and call the function. -------------------#
      u3d    = trans3d(x=dat$x,y=dat$y,z=height,pmat=pmat)
      argnow = modifyList(x=dotdotdot,val=list(x=u3d$x,y=u3d$y))
      dummy  = do.call(what="lines",args=argnow)
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#
}#end function
#------------------------------------------------------------------------------------------#
