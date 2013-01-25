#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#     This function plots roads, given the road name and the path.  Additional options for #
# the "lines" command are welcome!                                                         #
#------------------------------------------------------------------------------------------#
roadmap <<- function(rpath=NULL,roads,...){

  if (is.null(rpath)) rpath = file.path(srcdir,"roads")
  if (missing(roads)) roads = sub(pattern=".csv",replacement="",x=basename(dir(rpath)))

  roadfiles = paste(rpath,"/",roads,".csv",sep="")
  nroads    = length(roads)

  for (r in 1:nroads){
     datum=read.csv(file=roadfiles[r],header=TRUE,na.strings="NA")
     names(datum) = c("lon","lat","hgt")
     lines(x=datum$lon,y=datum$lat,...)
  }#end for
} #end function
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#     This function plots roads, given the road name and the path, when using the lattice  #
# plots.  Additional options for the "llines" command are welcome!                         #
#------------------------------------------------------------------------------------------#
panel.roadmap <<- function(rpath=NULL,roads,...){

  if (is.null(rpath)) rpath = file.path(srcdir,"roads")
  if (missing(roads)) roads = sub(pattern=".csv",replacement="",x=basename(dir(rpath)))

  roadfiles = paste(rpath,"/",roads,".csv",sep="")
  nroads    = length(roads)

  for (r in 1:nroads){
     datum=read.csv(file=roadfiles[r],header=TRUE,na.strings="NA")
     names(datum) = c("lon","lat","hgt")
     panel.lines(x=datum$lon,y=datum$lat,...)
  }#end for
} #end function
#------------------------------------------------------------------------------------------#
