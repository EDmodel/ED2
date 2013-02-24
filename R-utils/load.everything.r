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
#      Define the background and foreground colour.                                        #
#------------------------------------------------------------------------------------------#
if (! "ibackground" %in% ls()) ibackground = 0
if (ibackground == 0){
   foreground    <<- "black"
   background    <<- "white"
   grid.colour   <<- "grey66"
   map.colour    <<- "black"
   miss.colour   <<- "grey94"
   all.colour    <<- "grey22"
}else if (ibackground == 1){
   foreground    <<- "white"
   background    <<- "black"
   grid.colour   <<- "grey33"
   map.colour    <<- "white"
   miss.colour   <<- "grey10"
   all.colour    <<- "grey78"
   gcol          <<- c("lightblue","white")
}else if (ibackground == 2){
   foreground    <<- "white"
   background    <<- "#282828"
   grid.colour   <<- "grey38"
   map.colour    <<- "white"
   miss.colour   <<- "grey20"
   all.colour    <<- "grey73"
   gcol          <<- c("lightblue","#282828")
}else{
   stop(paste(" Invalid ibackground value (",ibackground,")",sep=""))
}#end if
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#      Define the size of the titles and axes.                                             #
#------------------------------------------------------------------------------------------#
if (! "ptsz" %in% ls()){
   ptsz <<- 16
}else{
   ptsz <<- ptsz
}#end if
if (ptsz <= 14){
   cex.main      <<- 1.1
   cex.lab       <<- 1.0
}else if (ptsz <= 16){
   cex.main      <<- 0.9
   cex.lab       <<- 1.0
}else{
   cex.main      <<- 0.7
   cex.lab       <<- 0.9
}#end if
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Define the colours for different hues.                                               #
#------------------------------------------------------------------------------------------#
if (ibackground == 0){
   grey.fg   <<- "grey25"        ;grey.mg   <<- "grey50"     ;grey.bg   <<- "grey75"
   purple.fg <<- "purple4"       ;purple.mg <<- "purple1"    ;purple.bg <<- "mediumpurple1"
   indigo.fg <<- "slateblue4"    ;indigo.mg <<- "slateblue"  ;indigo.bg <<- "slateblue1"
   blue.fg   <<- "royalblue4"    ;blue.mg   <<- "royalblue"  ;blue.bg   <<- "steelblue2"
   sky.fg    <<- "dodgerblue3"   ;sky.mg    <<- "deepskyblue";sky.bg    <<- "lightskyblue1"
   pink.fg   <<- "deeppink3"     ;pink.mg   <<- "hotpink"    ;pink.bg   <<- "orchid1"
   green.fg  <<- "darkgreen"     ;green.mg  <<- "chartreuse4";green.bg  <<- "chartreuse"
   olive.fg  <<- "darkolivegreen";olive.mg  <<- "olivedrab4" ;olive.bg  <<- "olivedrab2"
   khaki.fg  <<- "burlywood4"    ;khaki.mg  <<- "burlywood1" ;khaki.bg  <<- "bisque3"
   yellow.fg <<- "saddlebrown"   ;yellow.mg <<- "yellow3"    ;yellow.bg <<- "yellow"
   orange.fg <<- "orangered"     ;orange.mg <<- "orange3"    ;orange.bg <<- "gold"
   red.fg    <<- "red3"          ;red.mg    <<- "firebrick1" ;red.bg    <<- "indianred1"
}else{
   grey.bg   <<- "grey75"        ;grey.mg   <<- "grey50"     ;grey.fg   <<- "grey25"
   purple.bg <<- "purple4"       ;purple.mg <<- "purple1"    ;purple.fg <<- "mediumpurple1"
   indigo.bg <<- "slateblue4"    ;indigo.mg <<- "slateblue"  ;indigo.fg <<- "slateblue1"
   blue.bg   <<- "royalblue4"    ;blue.mg   <<- "royalblue"  ;blue.fg   <<- "steelblue2"
   sky.bg    <<- "dodgerblue3"   ;sky.mg    <<- "deepskyblue";sky.fg    <<- "lightskyblue1"
   pink.bg   <<- "deeppink3"     ;pink.mg   <<- "hotpink"    ;pink.fg   <<- "orchid1"
   green.bg  <<- "darkgreen"     ;green.mg  <<- "chartreuse4";green.fg  <<- "chartreuse"
   olive.bg  <<- "darkolivegreen";olive.mg  <<- "olivedrab4" ;olive.fg  <<- "olivedrab2"
   khaki.bg  <<- "burlywood4"    ;khaki.mg  <<- "burlywood1" ;khaki.fg  <<- "bisque3"
   yellow.bg <<- "saddlebrown"   ;yellow.mg <<- "yellow3"    ;yellow.fg <<- "yellow"
   orange.bg <<- "orangered"     ;orange.mg <<- "orange3"    ;orange.fg <<- "gold"
   red.bg    <<- "red3"          ;red.mg    <<- "firebrick1" ;red.fg    <<- "indianred1"
}#end if
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Make some rainbows with 7 colours, always from darker to lighter.                    #
#------------------------------------------------------------------------------------------#
grey.rbow    <<- c("grey25","grey33","grey41","grey49","grey57","grey65","grey73")
purple.rbow  <<- c("purple4","purple2","mediumpurple3","mediumpurple1"
                  ,"slateblue3","slateblue1","#B0B0FF")
blue.rbow    <<- c("midnightblue","royalblue4","royalblue2","dodgerblue"
                   ,"deepskyblue","cadetblue3","darkslategray1")
green.rbow   <<- c("#004000","chartreuse4","chartreuse","olivedrab3","olivedrab1"
                  ,"darkseagreen3","darkseagreen1")
orange.rbow  <<- c("#502800","sienna","darkorange1","orange","gold"
                  ,"lightgoldenrod2","khaki1")
red.rbow     <<- c("red4","firebrick3","tomato","lightcoral","palevioletred"
                  ,"lightpink3","lightpink")
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
at.first      = c("rconstants.r","globdims.r")
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
for (n in 1:nscripts){
   full = file.path(srcdir,order.scripts[n])
   isok = try(source(full),silent=TRUE)
   if ("try-error" %in% is(isok)){
      options(warn=warn.orig)
      cat("   - Script ",order.scripts[n]," has bugs!  Check the errors/warnings: ","\n")
      source(full)
      stop("Source code problem")
   }#end if
}#end for
options(warn=warn.orig)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#       Get rid of the extremely annoying and unnecessary bell.  Also, force the system to #
# use Helvetica as the default font family.                                                #
#------------------------------------------------------------------------------------------#
options(locatorBell=FALSE,family="Helvetica")
#------------------------------------------------------------------------------------------#

#==========================================================================================#
#==========================================================================================#
