#==========================================================================================#
#==========================================================================================#
#     This script loads all other scripts in this path, and also loads all the necessary   #
# packages.                                                                                #
#------------------------------------------------------------------------------------------#
if ("srcdir" %in% ls()){
   srcdir <<- srcdir
}else{
   srcdir <<- getwd()
}#end if
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Find which major version of R is calling this script.                               #
#------------------------------------------------------------------------------------------#
R.major <<- as.numeric(R.version$major)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Make the screen output as wide as the screen permits.                               #
#------------------------------------------------------------------------------------------#
ncstring = as.integer(Sys.getenv("COLUMNS"))
if (! is.na(ncstring)){
   if (ncstring > 80 & ncstring < 500) options(width=ncstring)
}#end if (! is.na(ncstring))
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Define the default colours for a white background and black foreground.             #
#------------------------------------------------------------------------------------------#
grid.colour   = "grey66"
map.colour    = "black"
miss.colour   = "grey94"
miss.colour.2 = "grey86"
all.colour    = "grey22"
washed.3g     = c("grey80"              ,"grey87"              ,"grey94"              )
grey.3g       = c("grey25"              ,"grey50"              ,"grey75"              )
purple.3g     = c("purple4"             ,"purple1"             ,"mediumpurple1"       )
indigo.3g     = c("slateblue4"          ,"slateblue"           ,"#A0A0FF"             )
blue.3g       = c("midnightblue"        ,"dodgerblue3"         ,"dodgerblue"          )
royalblue.3g  = c("royalblue4"          ,"royalblue"           ,"#6496FF"             )
steelblue.3g  = c("steelblue4"          ,"steelblue3"          ,"steelblue1"          )
sky.3g        = c("deepskyblue3"        ,"deepskyblue"         ,"skyblue1"            )
pink.3g       = c("deeppink3"           ,"hotpink"             ,"orchid1"             )
green.3g      = c("#004000"             ,"green4"              ,"green2"              )
chartreuse.3g = c("#408000"             ,"#60C000"             ,"#80FF00"             )
olive.3g      = c("darkolivegreen"      ,"olivedrab4"          ,"olivedrab2"          )
khaki.3g      = c("burlywood4"          ,"burlywood1"          ,"bisque3"             )
yellow.3g     = c("saddlebrown"         ,"yellow3"             ,"yellow"              )
gold.3g       = c("darkgoldenrod"       ,"goldenrod"           ,"lightgoldenrod"      )
orange.3g     = c("orangered"           ,"orange3"             ,"gold"                )
red.3g        = c("red3"                ,"red"                 ,"lightcoral"          )
firebrick.3g  = c("firebrick4"          ,"firebrick2"          ,"#FFA0A0"             )
brown.3g      = c("#422D11"             ,"#5A3D17"             ,"#7E5123"             )
grey.7g       = c("grey25"              ,"grey33"              ,"grey41"
                 ,"grey49"              ,"grey57"              ,"grey65"
                 ,"grey73"
                 )#end c
purple.7g     = c("purple4"             ,"purple2"             ,"mediumpurple3"
                 ,"mediumpurple1"       ,"slateblue3"          ,"slateblue1"  
                 ,"#B0B0FF"
                 )#end c
blue.7g       = c("midnightblue"        ,"royalblue4"          ,"royalblue2"
                 ,"dodgerblue"          ,"deepskyblue"         ,"cadetblue3"
                 ,"darkslategray1"
                 )#end c
green.7g      = c("#004000"             ,"chartreuse4"         ,"chartreuse"
                 ,"olivedrab3"          ,"olivedrab1"          ,"darkseagreen3"
                 ,"darkseagreen1"
                 )#end c
orange.7g     = c("#502800"             ,"sienna"              ,"darkorange1"
                 ,"orange"              ,"gold"                ,"lightgoldenrod2"
                 ,"khaki1"
                 )#end c
red.7g        = c("red4"                ,"firebrick3"          ,"tomato"
                 ,"lightcoral"          ,"palevioletred"       ,"lightpink3"
                 ,"lightpink"
                 )#end c
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Fix the colours according to the current background.                                 #
#------------------------------------------------------------------------------------------#
if (! "ibackground" %in% ls()) ibackground = 0
if (ibackground == 0){
  foreground    = "black"
  background    = "white"
}else if (ibackground == 1){
   foreground    <<- "white"
   background    <<- "black"
}else if (ibackground == 2){
   foreground    <<- "white"
   background    <<- "#282828"
}else{
   stop(paste0(" Invalid ibackground value (",ibackground,")"))
}#end if
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Correct all colours.                                                                 #
#------------------------------------------------------------------------------------------#
source(file.path(srcdir,"switch.col.r"))
grid.colour   <<- negative.col(x=grid.colour  ,fg=foreground,bg=background)
map.colour    <<- negative.col(x=map.colour   ,fg=foreground,bg=background)
miss.colour   <<- negative.col(x=miss.colour  ,fg=foreground,bg=background)
miss.colour.2 <<- negative.col(x=miss.colour.2,fg=foreground,bg=background)
all.colour    <<- negative.col(x=all.colour   ,fg=foreground,bg=background)
gcol          <<- c(switch.col(x="lightblue"  ,fg=foreground,bg=background),background)
washed.3g     <<- switch.col  (x=washed.3g    ,fg=foreground,bg=background)
grey.3g       <<- switch.col  (x=grey.3g      ,fg=foreground,bg=background)
purple.3g     <<- switch.col  (x=purple.3g    ,fg=foreground,bg=background)
indigo.3g     <<- switch.col  (x=indigo.3g    ,fg=foreground,bg=background)
blue.3g       <<- switch.col  (x=blue.3g      ,fg=foreground,bg=background)
royalblue.3g  <<- switch.col  (x=royalblue.3g ,fg=foreground,bg=background)
steelblue.3g  <<- switch.col  (x=steelblue.3g ,fg=foreground,bg=background)
sky.3g        <<- switch.col  (x=sky.3g       ,fg=foreground,bg=background)
pink.3g       <<- switch.col  (x=pink.3g      ,fg=foreground,bg=background)
green.3g      <<- switch.col  (x=green.3g     ,fg=foreground,bg=background)
chartreuse.3g <<- switch.col  (x=chartreuse.3g,fg=foreground,bg=background)
olive.3g      <<- switch.col  (x=olive.3g     ,fg=foreground,bg=background)
khaki.3g      <<- switch.col  (x=khaki.3g     ,fg=foreground,bg=background)
yellow.3g     <<- switch.col  (x=yellow.3g    ,fg=foreground,bg=background)
gold.3g       <<- switch.col  (x=gold.3g      ,fg=foreground,bg=background)
orange.3g     <<- switch.col  (x=orange.3g    ,fg=foreground,bg=background)
red.3g        <<- switch.col  (x=red.3g       ,fg=foreground,bg=background)
firebrick.3g  <<- switch.col  (x=firebrick.3g ,fg=foreground,bg=background)
brown.3g      <<- switch.col  (x=brown.3g     ,fg=foreground,bg=background)
grey.rbow     <<- switch.col  (x=grey.7g      ,fg=foreground,bg=background,may.rev=FALSE)
purple.rbow   <<- switch.col  (x=purple.7g    ,fg=foreground,bg=background,may.rev=FALSE)
blue.rbow     <<- switch.col  (x=blue.7g      ,fg=foreground,bg=background,may.rev=FALSE)
green.rbow    <<- switch.col  (x=green.7g     ,fg=foreground,bg=background,may.rev=FALSE)
orange.rbow   <<- switch.col  (x=orange.7g    ,fg=foreground,bg=background,may.rev=FALSE)
red.rbow      <<- switch.col  (x=red.7g       ,fg=foreground,bg=background,may.rev=FALSE)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Define the colours for different hues.                                               #
#------------------------------------------------------------------------------------------#
fmg.cols  <<- c("washed","grey","purple","indigo","blue","royalblue","steelblue","sky"
               ,"pink","green","chartreuse","olive","khaki","yellow","gold","orange","red"
               ,"firebrick","brown")
for (fmg.col in fmg.cols){
   threeg = get(paste(fmg.col,"3g",sep="."))
   assign(x=paste(fmg.col,"fg",sep="."),value=threeg[1],envir=.GlobalEnv)
   assign(x=paste(fmg.col,"mg",sep="."),value=threeg[2],envir=.GlobalEnv)
   assign(x=paste(fmg.col,"bg",sep="."),value=threeg[3],envir=.GlobalEnv)
}#end for
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#      Define the size of the titles and axes.                                             #
#------------------------------------------------------------------------------------------#
if (! "ptsz" %in% ls()){
   ptsz <<- 16
}else{
   ptsz <<- ptsz
}#end if
cex.ptsz <<- 1.0 * min(1.0,ptsz / 15)
cex.main <<- 1.1 * min(1.0,ptsz / 14)
cex.lab  <<- 1.0 * min(1.0,ptsz / 14)
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#      Define the integration period for photosynthesis-related variables.                 #
#------------------------------------------------------------------------------------------#
if (! "iint.photo" %in% ls()){
   iint.photo <<- 0
}else{
   iint.photo <<- iint.photo
}#end if
#------------------------------------------------------------------------------------------#






#----- Create the default plotting settings for R. ----------------------------------------#
par.user <<- list( bg       = "transparent"
                 , col      = foreground
                 , col.axis = foreground
                 , col.lab  = foreground
                 , col.main = foreground
                 , col.sub  = foreground
                 , fg       = foreground
                 , cex.main = cex.main
                 , cex.lab  = cex.lab
                 , family   = "Helvetica"
                 , mar      = c(5.1,4.4,4.1,2.1)
                 )#end list
#------------------------------------------------------------------------------------------#


#----- Wrapper for loading packages without pop ups. --------------------------------------#
discreet.require <<- function(...){
   dummy = suppressPackageStartupMessages(suppressWarnings(require(...)))
   return(dummy)
}#end discreet.require
#------------------------------------------------------------------------------------------#

#-----------------------------------------------------------------------------------------#
#     Load all packages needed.                                                            #
#------------------------------------------------------------------------------------------#
loaded.package = list()
loaded.package[["abind"       ]] = discreet.require(abind       )
loaded.package[["agricolae"   ]] = discreet.require(agricolae   )
loaded.package[["akima"       ]] = discreet.require(akima       )
loaded.package[["beanplot"    ]] = discreet.require(beanplot    )
loaded.package[["boot"        ]] = discreet.require(boot        )
loaded.package[["car"         ]] = discreet.require(car         )
loaded.package[["caTools"     ]] = discreet.require(caTools     )
loaded.package[["chron"       ]] = discreet.require(chron       )
loaded.package[["cluster"     ]] = discreet.require(cluster     )
loaded.package[["compiler"    ]] = discreet.require(compiler    )
loaded.package[["devtools"    ]] = discreet.require(devtools    )
loaded.package[["fields"      ]] = discreet.require(fields      )
loaded.package[["gbm"         ]] = discreet.require(gbm         )
loaded.package[["gdalUtils"   ]] = discreet.require(gdalUtils   )
loaded.package[["geoR"        ]] = discreet.require(geoR        )
loaded.package[["gpclib"      ]] = discreet.require(gpclib      )
loaded.package[["grDevices"   ]] = discreet.require(grDevices   )
loaded.package[["gstat"       ]] = discreet.require(gstat       )
loaded.package[["hdf5"        ]] = discreet.require(hdf5        )
loaded.package[["Hmisc"       ]] = discreet.require(Hmisc       )
loaded.package[["klaR"        ]] = discreet.require(klaR        )
loaded.package[["kriging"     ]] = discreet.require(kriging     )
loaded.package[["leaps"       ]] = discreet.require(leaps       )
loaded.package[["maps"        ]] = discreet.require(maps        )
loaded.package[["mapdata "    ]] = discreet.require(mapdata     )
loaded.package[["MASS"        ]] = discreet.require(MASS        )
loaded.package[["MCMCpack"    ]] = discreet.require(MCMCpack    )
loaded.package[["nlme"        ]] = discreet.require(nlme        )
loaded.package[["numDeriv"    ]] = discreet.require(numDeriv    )
loaded.package[["onls"        ]] = discreet.require(onls        )
loaded.package[["PBSmapping"  ]] = discreet.require(PBSmapping  )
loaded.package[["plotrix"     ]] = discreet.require(plotrix     )
loaded.package[["proto"       ]] = discreet.require(proto       )
loaded.package[["randomForest"]] = discreet.require(randomForest)
loaded.package[["raster"      ]] = discreet.require(raster      )
loaded.package[["rgdal"       ]] = discreet.require(rgdal       )
loaded.package[["rgeos"       ]] = discreet.require(rgeos       )
loaded.package[["rJava"       ]] = discreet.require(rJava       )
loaded.package[["robustbase"  ]] = discreet.require(robustbase  )
loaded.package[["rworldmap"   ]] = discreet.require(rworldmap   )
loaded.package[["RSEIS"       ]] = discreet.require(RSEIS       )
loaded.package[["R.utils"     ]] = discreet.require(R.utils     )
loaded.package[["shapefiles"  ]] = discreet.require(shapefiles  )
loaded.package[["splancs"     ]] = discreet.require(splancs     )
loaded.package[["smatr"       ]] = discreet.require(smatr       )
loaded.package[["sn"          ]] = discreet.require(sn          )
loaded.package[["sp"          ]] = discreet.require(sp          )
loaded.package[["stats4"      ]] = discreet.require(stats4      )
loaded.package[["vioplot"     ]] = discreet.require(vioplot     )
loaded.package[["VoxR"        ]] = discreet.require(VoxR        )
loaded.package[["zoo"         ]] = discreet.require(zoo         )
loaded.package = unlist(loaded.package)
if (! all(loaded.package)){
   miss = which(! loaded.package)
   cat(" You must install the following packages before using the scripts:","\n")
   for (m in miss) cat(" -> ",names(loaded.package)[m],"\n",sep="")
   risky = readline(" Are you sure you want to proceed [y|N]? ")
   risky = tolower(risky)
   if (! risky %in% c("y","yes")) stop("Missing packages!!!")
}#end if
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#  SHADY BUSINESS...  We must unlock grav from package boot and replace by our good        #
#                     old value from rconstants.r.                                         #
#------------------------------------------------------------------------------------------#
envir = as.environment("package:boot")
try(unlockBinding("grav",envir),silent=TRUE)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#  SHADY BUSINESS...  We must unlock trim from package R.oo and raster and replace by our  #
#                     function that has more options than the package one.                 #
#------------------------------------------------------------------------------------------#
envir = as.environment("package:R.oo")
try(unlockBinding("trim",envir),silent=TRUE)
envir = as.environment("package:raster")
try(unlockBinding("trim",envir),silent=TRUE)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#  SHADY BUSINESS...  We must unlock RGB from package raster and replace by our function.  #
#------------------------------------------------------------------------------------------#
envir = as.environment("package:raster")
try(unlockBinding("RGB",envir),silent=TRUE)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#  SHADY BUSINESS...  We must unlock tobin from package survival and replace by our        #
# function.                                                                                #
#------------------------------------------------------------------------------------------#
envir = as.environment("package:survival")
try(unlockBinding("tobin",envir),silent=TRUE)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#  SHADY BUSINESS...  We must unlock pspline from package survival and replace by our      #
# function.                                                                                #
#------------------------------------------------------------------------------------------#
envir = as.environment("package:survival")
try(unlockBinding("pspline",envir),silent=TRUE)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#  SHADY BUSINESS...  We must unlock theme from package ggplot2 and replace by our         #
# function.                                                                                #
#------------------------------------------------------------------------------------------#
envir = as.environment("package:ggplot2")
try(unlockBinding("theme",envir),silent=TRUE)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Organise the files so we load them in the right order.                               #
#------------------------------------------------------------------------------------------#
at.first      = c("colour.utils.r","rconstants.r","globdims.r","unitlist.r")
at.end        = c("iso3166.r","pft.coms.r","pmonthly_varlist.r","pcomp_varlist.r")
myself        = c("load.everything.r")
all.scripts   = sort(list.files(path=srcdir,pattern="\\.[RrSsQq]$"))
back.up       = sort(list.files(path=srcdir,pattern="^[~]"))
keep          = ! ( all.scripts %in% at.first
                  | all.scripts %in% at.end
                  | all.scripts %in% myself
                  | all.scripts %in% back.up
                  )#end
middle        = all.scripts[keep]
order.scripts = c(at.first,middle,at.end)
nscripts      = length(order.scripts)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Load all files, in order.  Here we replace the warnings by errors, just to make sure #
# that all the functions are clean.                                                        #
#------------------------------------------------------------------------------------------#
warn.orig = getOption("warn")
options(warn=2)
cat(" + Load scripts from ",srcdir,".","\n",sep="")
for (iscript in sequence(nscripts)){
   script.now  = order.scripts[iscript]
   full        = file.path(srcdir,script.now)
   isok        = try(source(full),silent=TRUE)
   if ("try-error" %in% is(isok)){
      options(warn=warn.orig)
      cat("   - Script ",script.now," has bugs!  Check the errors/warnings: ","\n",sep="")
      source(full)
      stop("Source code problem")
   }#end if
}#end for
options(warn=warn.orig)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Check for fortran code to be loaded.                                                #
#------------------------------------------------------------------------------------------#
all.f90  = sort( c( list.files(path=srcdir,pattern="\\.[Ff]90$")
                  , list.files(path=srcdir,pattern="\\.[Ff]$")
                  )#end c
               )#end sort
nall.f90 = length(all.f90)
for (if90 in sequence(nall.f90)){
   fnow    = file.path(srcdir,all.f90[if90])
   flib.o  = fnow
   flib.o  = gsub(pattern = "\\.[Ff]90$",replacement=".o",x=flib.o)
   flib.o  = gsub(pattern = "\\.[Ff]$"  ,replacement=".o",x=flib.o)
   flib.so = fnow
   flib.so = gsub(pattern = "\\.[Ff]90$",replacement=".so",x=flib.so)
   flib.so = gsub(pattern = "\\.[Ff]$"  ,replacement=".so",x=flib.so)
   flib.sl = fnow
   flib.sl = gsub(pattern = "\\.[Ff]90$",replacement=".sl",x=flib.sl)
   flib.sl = gsub(pattern = "\\.[Ff]$"  ,replacement=".sl",x=flib.sl)

   #----- Select library. -----------------------------------------------------------------#
   if (file.exists(flib.so)){
      flib.sx = flib.so
   }else if (file.exists(flib.sl)){
      flib.sx = flib.sl
   }else{
      #----- This is guaranteed to fail, so it will force recompilation. ------------------#
      flib.sx = flib.o
      #------------------------------------------------------------------------------------#
   }#end if (file.exists(flib.so))
   #---------------------------------------------------------------------------------------#



   #----- Check whether dynamic library can be loaded.  In case not, recompile. -----------#
   dummy = try(dyn.load(flib.sx),silent=TRUE)
   if ("try-error" %in% is(dummy)){
      dummy = file.remove(flib.so)
      dummy = file.remove(flib.sl)
      dummy = file.remove(flib.o )
      dummy = RCMD(cmd="SHLIB",options=fnow,path=srcdir)
   }#end if ("try-error" %in% is(dummy))
   #---------------------------------------------------------------------------------------#
}#end for (if90 in sequence(nall.f90))
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      CBal had unit issues in some versions, check whether to fix it in the post-         #
# processor.                                                                               #
#------------------------------------------------------------------------------------------#
if (! "kludgecbal" %in% ls()){
   kludgecbal <<- FALSE
}else{
   kludgecbal <<- kludgecbal
}#end if
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Growth, storage, and vleaf respiration had wrong units in the output in previous    #
# versions, check whether to correct them.                                                 #
#------------------------------------------------------------------------------------------#
if (! "corr.growth.storage" %in% ls()){
   corr.growth.storage <<- 1.
}else{
   corr.growth.storage <<- corr.growth.storage
}#end if
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#       Get rid of the extremely annoying and unnecessary bell.  Also, force the system to #
# use Helvetica as the default font family.                                                #
#------------------------------------------------------------------------------------------#
options(locatorBell=FALSE,family="Helvetica")
#------------------------------------------------------------------------------------------#

#==========================================================================================#
#==========================================================================================#
