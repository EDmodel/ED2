#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#  This function plots the map for the Amazon forest.                                      #
#------------------------------------------------------------------------------------------#
amazonmap = function(mapdir=NULL,legal=FALSE,...){
  if (is.null(mapdir)) mapdir = file.path(srcdir,"amzmap")
  if (legal){
     fichier = file.path(mapdir,"brazilian_legal_amazon.txt")
  }else{
     fichier = file.path(mapdir,"amazon.ecoregion.txt")
  }#end if
  the.amazon = read.table(fichier,na.string="NA")
  names(the.amazon) = c("x","y")
  lines(x=the.amazon,...)
} #end function
#------------------------------------------------------------------------------------------#





#------------------------------------------------------------------------------------------#
#     This function plots the map for the Amazon forest.  This is useful for when running  #
# lattice package.                                                                         #
#------------------------------------------------------------------------------------------#
panel.amazonmap = function(mapdir=NULL,legal=FALSE,...){
  if (is.null(mapdir)) mapdir = file.path(srcdir,"amzmap")
  if (legal){
     fichier = file.path(mapdir,"amazon.forest.txt")
  }else{
     fichier = file.path(mapdir,"amazon.ecoregion.txt")
  }#end if
  the.amazon = read.table(fichier,na.string="NA")
  names(the.amazon) = c("x","y")
  panel.points(x=the.amazon,type="l",...)
} #end function
#------------------------------------------------------------------------------------------#





#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#  This function plots the map for the Amazon forest in 3D maps.                           #
#------------------------------------------------------------------------------------------#
amazonmap3d = function(height=0,pmat,xlim=NULL,ylim=NULL,mapdir=NULL,legal=FALSE,...){

  if (is.null(mapdir)) mapdir = file.path(srcdir,"amzmap")
  if (legal){
     fichier = file.path(mapdir,"amazon.forest.txt")
  }else{
     fichier = file.path(mapdir,"amazon.ecoregion.txt")
  }#end if
  the.amazon = read.table(fichier,na.string="NA")
  names(the.amazon) = c("x","y")
  if (! is.null(xlim)){
     sel = the.amazon$x >= xlim[1] & the.amazon$x <= xlim[2]
     the.amazon = the.amazon[sel,]
  }#end if
  if (! is.null(ylim)){
     sel = the.amazon$y >= ylim[1] & the.amazon$y <= ylim[2]
     the.amazon = the.amazon[sel,]
  }#end if
  lines(trans3d(x=the.amazon$x,y=the.amazon$y,z=height,pmat),...)
} #end function
#------------------------------------------------------------------------------------------#
