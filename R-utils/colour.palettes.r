#==========================================================================================#
#==========================================================================================#
#     List of functions that generate colour palettes.  They all call a single inter-      #
# polation function.  Most functions are derived from QGIS colour palettes.                #
#------------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Blues.                                                                           #
   #---------------------------------------------------------------------------------------#
   blues <<- function(n,alpha=1.0){
      #----- Color entries. ---------------------------------------------------------------#
      nodes = c( "#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6"
               , "#2171B5", "#08519C", "#08306B")
      #------------------------------------------------------------------------------------#

      #----- Call the interpolator. -------------------------------------------------------#
      ans   = colour.interpol(nodes=nodes,n=n)
      ans   = scales::alpha(colour=ans,alpha=alpha)
      return(ans)
      #------------------------------------------------------------------------------------#
   }#end function blues
   #----- Inverse colour palette. ---------------------------------------------------------#
   iblues <- function(n,alpha=1.0){ rev(blues(n=n,alpha=alpha))}
   #---------------------------------------------------------------------------------------#






   #---------------------------------------------------------------------------------------#
   #      BrBG.                                                                            #
   #---------------------------------------------------------------------------------------#
   brbg <<- function(n,alpha=1.0){
      #----- Color entries. ---------------------------------------------------------------#
      nodes = c( "#A6611A", "#DFC27D", "#F5F5F5", "#80CDC1", "#018571")
      #------------------------------------------------------------------------------------#

      #----- Call the interpolator. -------------------------------------------------------#
      ans   = colour.interpol(nodes=nodes,n=n)
      ans   = scales::alpha(colour=ans,alpha=alpha)
      return(ans)
      #------------------------------------------------------------------------------------#
   }#end function brbg
   #----- Inverse colour palette. ---------------------------------------------------------#
   ibrbg <- function(n,alpha=1.0){ rev(brbg(n=n,alpha=alpha))}
   #---------------------------------------------------------------------------------------#







   #---------------------------------------------------------------------------------------#
   #      BuGn.                                                                            #
   #---------------------------------------------------------------------------------------#
   bugn <<- function(n,alpha=1.0){
      #----- Color entries. ---------------------------------------------------------------#
      nodes = c( "#C2D2F2", "#7AB1CC", "#63A6A6", "#338066", "#115929", "#003300")
      #------------------------------------------------------------------------------------#

      #----- Call the interpolator. -------------------------------------------------------#
      ans   = colour.interpol(nodes=nodes,n=n)
      ans   = scales::alpha(colour=ans,alpha=alpha)
      return(ans)
      #------------------------------------------------------------------------------------#
   }#end function bugn
   #----- Inverse colour palette. ---------------------------------------------------------#
   ibugn <- function(n,alpha=1.0){ rev(bugn(n=n,alpha=alpha))}
   #---------------------------------------------------------------------------------------#







   #---------------------------------------------------------------------------------------#
   #      BuGy.                                                                            #
   #---------------------------------------------------------------------------------------#
   bugy <<- function(n,alpha=1.0){
      #----- Color entries. ---------------------------------------------------------------#
      nodes = c( "#0050CE", "#1E64FF", "#2896FF", "#78D2FF", "#B4E6FF", "#FDFDFD"
               , "#DEDEDE", "#B6B6B6", "#8E8E8E", "#666666", "#3E3E3E")
      #------------------------------------------------------------------------------------#

      #----- Call the interpolator. -------------------------------------------------------#
      ans   = colour.interpol(nodes=nodes,n=n)
      ans   = scales::alpha(colour=ans,alpha=alpha)
      return(ans)
      #------------------------------------------------------------------------------------#
   }#end function bugy
   #----- Inverse colour palette. ---------------------------------------------------------#
   ibugy <- function(n,alpha=1.0){ rev(bugy(n=n,alpha=alpha))}
   #---------------------------------------------------------------------------------------#







   #---------------------------------------------------------------------------------------#
   #      BuPu.                                                                            #
   #---------------------------------------------------------------------------------------#
   bupu <<- function(n,alpha=1.0){
      #----- Color entries. ---------------------------------------------------------------#
      nodes     = c( "#EDF8FB", "#CCE0EE", "#AEC5DF", "#97A6CF", "#8B84BD", "#895FAC"
                   , "#853795", "#810F7C")
      #------------------------------------------------------------------------------------#

      #----- Call the interpolator. -------------------------------------------------------#
      ans   = colour.interpol(nodes=nodes,n=n)
      ans   = scales::alpha(colour=ans,alpha=alpha)
      return(ans)
      #------------------------------------------------------------------------------------#
   }#end function bupu
   #----- Inverse colour palette. ---------------------------------------------------------#
   ibupu <- function(n,alpha=1.0){ rev(bupu(n=n,alpha=alpha))}
   #---------------------------------------------------------------------------------------#







   #---------------------------------------------------------------------------------------#
   #      Greys.                                                                          #
   #---------------------------------------------------------------------------------------#
   greys <<- function(n,alpha=1.0){
      #----- Color entries. ---------------------------------------------------------------#
      nodes     = c("#fafafa","#050505")
      #------------------------------------------------------------------------------------#

      #----- Call the interpolator. -------------------------------------------------------#
      ans   = colour.interpol(nodes=nodes,n=n)
      ans   = scales::alpha(colour=ans,alpha=alpha)
      return(ans)
      #------------------------------------------------------------------------------------#
   }#end function greys
   #----- Inverse colour palette. ---------------------------------------------------------#
   igreys <- function(n,alpha=1.0){ rev(greys(n=n,alpha=alpha))}
   #---------------------------------------------------------------------------------------#







   #---------------------------------------------------------------------------------------#
   #      OrRd.                                                                            #
   #---------------------------------------------------------------------------------------#
   orrd <<- function(n,alpha=1.0){
      #----- Color entries. ---------------------------------------------------------------#
      nodes     = c( "#FEFAE4", "#FEE7B7", "#F2BB83", "#CC7549", "#A6341E", "#800000")
      #------------------------------------------------------------------------------------#

      #----- Call the interpolator. -------------------------------------------------------#
      ans   = colour.interpol(nodes=nodes,n=n)
      ans   = scales::alpha(colour=ans,alpha=alpha)
      return(ans)
      #------------------------------------------------------------------------------------#
   }#end function orrd
   #----- Inverse colour palette. ---------------------------------------------------------#
   iorrd <- function(n,alpha=1.0){ rev(orrd(n=n,alpha=alpha))}
   #---------------------------------------------------------------------------------------#







   #---------------------------------------------------------------------------------------#
   #      PiYG.                                                                            #
   #---------------------------------------------------------------------------------------#
   piyg <<- function(n,alpha=1.0){
      #----- Color entries. ---------------------------------------------------------------#
      nodes     = c( "#D01C8B", "#F1B6DA", "#F7F7F7", "#B8E186", "#4DAC26")
      #------------------------------------------------------------------------------------#

      #----- Call the interpolator. -------------------------------------------------------#
      ans   = colour.interpol(nodes=nodes,n=n)
      ans   = scales::alpha(colour=ans,alpha=alpha)
      return(ans)
      #------------------------------------------------------------------------------------#
   }#end function piyg
   #----- Inverse colour palette. ---------------------------------------------------------#
   ipiyg <- function(n,alpha=1.0){ rev(piyg(n=n,alpha=alpha))}
   #---------------------------------------------------------------------------------------#







   #---------------------------------------------------------------------------------------#
   #      PRGn.                                                                            #
   #---------------------------------------------------------------------------------------#
   prgn <<- function(n,alpha=1.0){
      #----- Color entries. ---------------------------------------------------------------#
      nodes     = c( "#7B3294", "#C2A5CF", "#F7F7F7", "#A6DBA0", "#008837")
      #------------------------------------------------------------------------------------#

      #----- Call the interpolator. -------------------------------------------------------#
      ans   = colour.interpol(nodes=nodes,n=n)
      ans   = scales::alpha(colour=ans,alpha=alpha)
      return(ans)
      #------------------------------------------------------------------------------------#
   }#end function prgn
   #----- Inverse colour palette. ---------------------------------------------------------#
   iprgn <- function(n,alpha=1.0){ rev(prgn(n=n,alpha=alpha))}
   #---------------------------------------------------------------------------------------#







   #---------------------------------------------------------------------------------------#
   #      PuOr.                                                                            #
   #---------------------------------------------------------------------------------------#
   puor <<- function(n,alpha=1.0){
      #----- Color entries. ---------------------------------------------------------------#
      nodes     = c("#5E3C99", "#B2ABD2", "#F7F7F7", "#FDB863", "#E66101")
      #------------------------------------------------------------------------------------#

      #----- Call the interpolator. -------------------------------------------------------#
      ans   = colour.interpol(nodes=nodes,n=n)
      ans   = scales::alpha(colour=ans,alpha=alpha)
      return(ans)
      #------------------------------------------------------------------------------------#
   }#end function puor
   #----- Inverse colour palette. ---------------------------------------------------------#
   ipuor <- function(n,alpha=1.0){ rev(puor(n=n,alpha=alpha))}
   #---------------------------------------------------------------------------------------#







   #---------------------------------------------------------------------------------------#
   #      PuRd.                                                                            #
   #---------------------------------------------------------------------------------------#
   purd <<- function(n,alpha=1.0){
      #----- Color entries. ---------------------------------------------------------------#
      nodes     = c( "#F1EEF6", "#D7B5D8", "#DF65B0", "#DD1C77", "#980043")
      #------------------------------------------------------------------------------------#

      #----- Call the interpolator. -------------------------------------------------------#
      ans   = colour.interpol(nodes=nodes,n=n)
      ans   = scales::alpha(colour=ans,alpha=alpha)
      return(ans)
      #------------------------------------------------------------------------------------#
   }#end function purd
   #----- Inverse colour palette. ---------------------------------------------------------#
   ipurd <- function(n,alpha=1.0){ rev(purd(n=n,alpha=alpha))}
   #---------------------------------------------------------------------------------------#







   #---------------------------------------------------------------------------------------#
   #      PuBuGn.                                                                          #
   #---------------------------------------------------------------------------------------#
   pubugn <<- function(n,alpha=1.0){
      #----- Color entries. ---------------------------------------------------------------#
      nodes     = c( "#F6EFF7", "#D6DAEB", "#B1C5DF", "#80B2D4", "#51A2C0", "#2694A1"
                   , "#10817E", "#016C59")
      #------------------------------------------------------------------------------------#

      #----- Call the interpolator. -------------------------------------------------------#
      ans   = colour.interpol(nodes=nodes,n=n)
      ans   = scales::alpha(colour=ans,alpha=alpha)
      return(ans)
      #------------------------------------------------------------------------------------#
   }#end function pubugn
   #----- Inverse colour palette. ---------------------------------------------------------#
   ipubugn <- function(n,alpha=1.0){ rev(pubugn(n=n,alpha=alpha))}
   #---------------------------------------------------------------------------------------#







   #---------------------------------------------------------------------------------------#
   #      RdBu.                                                                            #
   #---------------------------------------------------------------------------------------#
   rdbu <<- function(n,alpha=1.0){
      #----- Color entries. ---------------------------------------------------------------#
      nodes = c("#CA0020", "#F4A582", "#F7F7F7", "#92C5DE", "#0571B0")
      #------------------------------------------------------------------------------------#

      #----- Call the interpolator. -------------------------------------------------------#
      ans   = colour.interpol(nodes=nodes,n=n)
      ans   = scales::alpha(colour=ans,alpha=alpha)
      return(ans)
      #------------------------------------------------------------------------------------#
   }#end function rdbu
   #----- Inverse colour palette. ---------------------------------------------------------#
   irdbu <- function(n,alpha=1.0){ rev(rdbu(n=n,alpha=alpha))}
   #---------------------------------------------------------------------------------------#







   #---------------------------------------------------------------------------------------#
   #      RdGy.                                                                            #
   #---------------------------------------------------------------------------------------#
   rdgy <<- function(n,alpha=1.0){
      #----- Color entries. ---------------------------------------------------------------#
      nodes = C( "#CA0020", "#F4A582", "#FFFFFF", "#BABABA", "#404040")
      #------------------------------------------------------------------------------------#

      #----- Call the interpolator. -------------------------------------------------------#
      ans   = colour.interpol(nodes=nodes,n=n)
      ans   = scales::alpha(colour=ans,alpha=alpha)
      return(ans)
      #------------------------------------------------------------------------------------#
   }#end function rdgy
   #----- Inverse colour palette. ---------------------------------------------------------#
   irdgy <- function(n,alpha=1.0){ rev(rdgy(n=n,alpha=alpha))}
   #---------------------------------------------------------------------------------------#







   #---------------------------------------------------------------------------------------#
   #      YlGnBu.                                                                          #
   #---------------------------------------------------------------------------------------#
   ylgnbu <<- function(n,alpha=1.0){
      #----- Color entries. ---------------------------------------------------------------#
      nodes     = c( "#FFFFCC", "#CAEABF", "#93D5B6", "#5CC1C0", "#3BA6C1", "#2F87BA"
                   , "#295FA9", "#253494")
      #------------------------------------------------------------------------------------#

      #----- Call the interpolator. -------------------------------------------------------#
      ans   = colour.interpol(nodes=nodes,n=n)
      ans   = scales::alpha(colour=ans,alpha=alpha)
      return(ans)
      #------------------------------------------------------------------------------------#
   }#end function ylgnbu
   #----- Inverse colour palette. ---------------------------------------------------------#
   iylgnbu <- function(n,alpha=1.0){ rev(ylgnbu(n=n,alpha=alpha))}
   #---------------------------------------------------------------------------------------#







   #---------------------------------------------------------------------------------------#
   #      YlOrRd.                                                                          #
   #---------------------------------------------------------------------------------------#
   ylorrd <<- function(n,alpha=1.0){
      #----- Color entries. ---------------------------------------------------------------#
      nodes     = c( "#FFFFB2", "#FECC5C", "#FD8D3C", "#F03B20", "#BD0026")
      #------------------------------------------------------------------------------------#

      #----- Call the interpolator. -------------------------------------------------------#
      ans   = colour.interpol(nodes=nodes,n=n)
      ans   = scales::alpha(colour=ans,alpha=alpha)
      return(ans)
      #------------------------------------------------------------------------------------#
   }#end function ylorrd
   #----- Inverse colour palette. ---------------------------------------------------------#
   iylorrd <- function(n,alpha=1.0){ rev(ylorrd(n=n,alpha=alpha))}
   #---------------------------------------------------------------------------------------#







   #---------------------------------------------------------------------------------------#
   #      For the following palettes, we use the default from package "viridis".  We only  #
   # define the inverse functions for convenience.                                         #
   #---------------------------------------------------------------------------------------#
   iinferno <- function(n,alpha=1.0){ inferno(n=n,alpha=alpha,direction=-1)}
   imagma   <- function(n,alpha=1.0){ magma  (n=n,alpha=alpha,direction=-1)}
   iviridis <- function(n,alpha=1.0){ viridis(n=n,alpha=alpha,direction=-1)}
   #---------------------------------------------------------------------------------------#





#==========================================================================================#
#==========================================================================================#
