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
   stop(paste(" Invalid ibackground value (",ibackground,")",sep=""))
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


#------------------------------------------------------------------------------------------#
#     Load all packages needed.                                                            #
#------------------------------------------------------------------------------------------#
isok = list()
isok[["abind"   ]] = require(abind   )
isok[["akima"   ]] = require(akima   )
isok[["boot"    ]] = require(boot    )
isok[["chron"   ]] = require(chron   )
isok[["fields"  ]] = require(fields  )
isok[["hdf5"    ]] = require(hdf5    )
isok[["klaR"    ]] = require(klaR    )
isok[["maps"    ]] = require(maps    )
isok[["mapdata "]] = require(mapdata )
isok[["MASS"    ]] = require(MASS    )
isok[["MCMCpack"]] = require(MCMCpack)
isok[["ncdf"    ]] = require(ncdf    )
isok[["numDeriv"]] = require(numDeriv)
isok[["plotrix" ]] = require(plotrix )
isok[["R.utils" ]] = require(R.utils )
isok[["RSEIS"   ]] = require(RSEIS   )
isok[["sn"      ]] = require(sn      )
isok[["zoo"     ]] = require(zoo     )
isok = unlist(isok)
if (! all(isok)){
   miss = which(! isok)
   cat(" You must install the following packages before using the scripts:","\n")
   for (m in miss) cat(" -> ",names(isok)[m],"\n")
   stop("Missing packages!!!")
}#end if
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#  SHADY BUSINESS...  We must unlock grav from package boot and replace by our good        #
#                     old value from rconstants.r.                                         #
#------------------------------------------------------------------------------------------#
envir = as.environment("package:boot")
unlockBinding("grav",envir)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#  SHADY BUSINESS...  We must unlock trim from package R.utils and replace by our function #
#                     that has more options than the package one.                          #
#------------------------------------------------------------------------------------------#
envir = as.environment("package:R.oo")
unlockBinding("trim",envir)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Organise the files so we load them in the right order.                               #
#------------------------------------------------------------------------------------------#
at.first      = c("rconstants.r","globdims.r","unitlist.r")
at.end        = c("pft.coms.r","pmonthly_varlist.r")
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
cat(" + Loading scripts from ",srcdir,"...","\n")
for (iscript in 1:nscripts){
   script.now  = order.scripts[iscript]
   full        = file.path(srcdir,script.now)
   isok        = try(source(full),silent=TRUE)
   if ("try-error" %in% is(isok)){
      options(warn=warn.orig)
      cat("   - Script ",script.now," has bugs!  Check the errors/warnings: ","\n")
      source(full)
      stop("Source code problem")
   }#end if
}#end for
options(warn=warn.orig)
#------------------------------------------------------------------------------------------#


if (! "kludgecbal" %in% ls()){
   kludgecbal <<- FALSE
}else{
   kludgecbal <<- kludgecbal
}#end if



#------------------------------------------------------------------------------------------#
#       Get rid of the extremely annoying and unnecessary bell.  Also, force the system to #
# use Helvetica as the default font family.                                                #
#------------------------------------------------------------------------------------------#
options(locatorBell=FALSE,family="Helvetica")
#------------------------------------------------------------------------------------------#

#==========================================================================================#
#==========================================================================================#
