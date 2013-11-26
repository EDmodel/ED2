#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#  This function plots the map for South America, including Brazilian states.              #
#------------------------------------------------------------------------------------------#
southammap <<- function(mapdir=NULL,...){
  if (is.null(mapdir)){
     arquivo = file.path(srcdir,"samap","americadosul")
  }else{
     arquivo = file.path(mapdir,"americadosul")
  }#end if
  americadosul = read.table(arquivo,na.string="999")
  names(americadosul) = c("lon","lat")
  lines(x=americadosul,...)
} #end function
#------------------------------------------------------------------------------------------#





#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#     This function plots the map for South America, including Brazilian states.  This is  #
# useful for when running lattice package.                                                 #
#------------------------------------------------------------------------------------------#
panel.southammap <<- function(mapdir=NULL,...){
  if (is.null(mapdir)){
     arquivo = file.path(srcdir,"samap","americadosul")
  }else{
     arquivo = file.path(mapdir,"americadosul")
  }#end if
  americadosul = read.table(arquivo,na.string="999")
  names(americadosul) = c("lon","lat")
  panel.points(x=americadosul,type="l",...)
} #end function
#------------------------------------------------------------------------------------------#





#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#  This function plots the map for South America, including Brazilian states.              #
#------------------------------------------------------------------------------------------#
southammap3d <<- function(height=0,pmat,xlim=NULL,ylim=NULL,mapdir=NULL,...){
  if (is.null(mapdir)){
     arquivo = file.path(srcdir,"samap","americadosul")
  }else{
     arquivo = file.path(mapdir,"americadosul")
  }#end if
  americadosul = read.table(arquivo,na.string="999")
  names(americadosul) = c("lon","lat")
  if (! is.null(xlim)){
     sel = americadosul$lon >= xlim[1] & americadosul$lon <= xlim[2]
     americadosul = americadosul[sel,]
  }#end if
  if (! is.null(ylim)){
     sel = americadosul$lat >= ylim[1] & americadosul$lat <= ylim[2]
     americadosul = americadosul[sel,]
  }#end if
  lines(trans3d(x=americadosul$lon,y=americadosul$lat,z=height,pmat),...)
} #end function
#------------------------------------------------------------------------------------------#
