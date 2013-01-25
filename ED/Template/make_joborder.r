#==========================================================================================#
#==========================================================================================#
#     Reset session.                                                                       #
#------------------------------------------------------------------------------------------#
rm(list=ls())
graphics.off()
#------------------------------------------------------------------------------------------#





#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#      Here is the user defined variable section.                                          #
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
here    = getwd()                                  #   Current directory
srcdir  = "/n/moorcroft_data/mlongo/util/Rsc"      #   Script directory
outfile = file.path(here,"newjoborder.txt")        #   Job order
#------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------#
#      Variables that will vary (check the list below for all possibilities.               #
#------------------------------------------------------------------------------------------#
# rscen    = c("r+000","r-020","r-040","r-060","r-080","r-100")
# tscen    = c("t+000","t+100","t+200","t+300")
# escen    = paste("real",sprintf("%2.2i",seq(from=0,to=9,by=1)),sep="-")
# 
# scenario = apply( X        = expand.grid(list(rscen,tscen,escen), stringsAsFactors = FALSE)
#                 , MARGIN   = 1
#                 , FUN      = paste
#                 , collapse = "_"
#                 )#end apply


varrun  = list( iata      = c("gyf","s67","rja")
              , iscenario = scenario
              , iphen     = c(-1,2)
              , yeara     = 1962
              , yearz     = 2012
              , isoilflg  = 2
              , istext    = c(16,6,2,8,11)
              )#end list
varlabel = list( iata      = varrun$iata
               , iscenario = scenario
               , iphen     = paste("iphen",sprintf("%+2.2i",varrun$iphen) ,sep="")
               , yeara     = "yr1962"
               , yearz     = "yr2012"
               , isoilflg  = "isoil02"
               , istext    = paste("stext",sprintf("%2.2i",varrun$istext),sep="")
               )#end list
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#      NO NEED TO CHANGE ANYTHING BEYOND THIS POINT UNLESS YOU ARE DEVELOPING THE CODE...  #
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#




#----- Load some packages and scripts. ----------------------------------------------------#
source(file.path(srcdir,"load.everything.r"))
#------------------------------------------------------------------------------------------#




#----- Check that varrun and varlabel have the same length and dimensions. ----------------#
if (length(varrun) != length(varlabel)){
   stop(" Varrun and varlabel must have the same number of variables!")
}else if (any(names(varrun) != names(varlabel))){
   stop(" Variable names of varrun and varlabel must match an be in the same order!")
}else if (any(sapply(X=varrun,FUN=length) != sapply(X=varlabel,FUN=length))){
   length.matrix = cbind( run   = sapply(X=varrun  ,FUN=length)
                        , label = sapply(X=varlabel,FUN=length)
                        )#end cbind
   print(length.matrix)
   stop(" Length of all variables in varrun and varlabel must match!")
}#end if
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Default properties.                                                                  #
#------------------------------------------------------------------------------------------#
default = list( run           = "unnamed"
              , iata          = "xxx"
              , lon           = 0.00
              , lat           = 0.00
              , yeara         = "1967"
              , montha        = "01"
              , daya          = "01"
              , timea         = "0000"
              , yearz         = "2013"
              , monthz        = "01"
              , dayz          = "01"
              , timez         = "0000"
              , init.mode     = 6
              , iscenario     = "default"
              , isizepft      = 0
              , iage          = 20
              , isoilflg      = 1
              , istext        = 1
              , sand          = -1.0
              , clay          = -1.0
              , depth         = "F"
              , isoilbc       = 1
              , sldrain       = 90.
              , scolour       = 14
              , slzres        = 0
              , queue         = "long_serial"
              , met.driver    = "tower"
              , dtlsm         = 600.
              , vmfact.c3     = 1.00
              , vmfact.c4     = 1.00
              , mphoto.trc3   = 9.0
              , mphoto.tec3   = 7.2
              , mphoto.c4     = 5.2
              , bphoto.blc3   = 10000.
              , bphoto.nlc3   = 1000.
              , bphoto.c4     = 10000.
              , kw.grass      = 900.
              , kw.tree       = 600.
              , gamma.c3      = 0.0145
              , gamma.c4      = 0.035
              , d0.grass      = 0.016
              , d0.tree       = 0.016
              , alpha.c3      = 0.080
              , alpha.c4      = 0.055
              , klowco2       = 4000.
              , decomp.scheme = 2
              , rrffact       = 1.000
              , growthresp    = 0.333
              , lwidth.grass  = 0.05
              , lwidth.bltree = 0.10
              , lwidth.nltree = 0.05
              , q10.c3        = 2.4
              , q10.c4        = 2.4
              , h2o.limit     = 2
              , imort.scheme  = 1
              , ddmort.const  = 0.8
              , isfclyrm      = 3
              , icanturb      = 2
              , ubmin         = 0.65
              , ugbmin        = 0.25
              , ustmin        = 0.05
              , gamm          = 13.0
              , gamh          = 13.0
              , tprandtl      = 0.74
              , ribmax        = 0.50
              , atmco2        = 378.
              , thcrit        = -1.20
              , sm.fire       = -1.40
              , ifire         = 0
              , fire.parm     = 0.5
              , ipercol       = 0
              , runoff.time   = 3600.
              , imetrad       = 2
              , ibranch       = 1
              , icanrad       = 1
              , crown.mod     = 0
              , ltrans.vis    = 0.050
              , lreflect.vis  = 0.100
              , ltrans.nir    = 0.230
              , lreflect.nir  = 0.460
              , orient.tree   = 0.100
              , orient.grass  = 0.000
              , clump.tree    = 0.800
              , clump.grass   = 1.000
              , ivegtdyn      = 1
              , igndvap       = 0
              , iphen         = -1
              , iallom        = 2
              , ibigleaf      = 0
              , irepro        = 2
              , treefall      = 0.0111
              ) #end list
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Create the full combination of runs.                                                 #
#------------------------------------------------------------------------------------------#
myruns = expand.grid(varrun,stringsAsFactors=FALSE)
nvars  = ncol(myruns)
nruns  = nrow(myruns)
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Create a table with room for all simulations.                                        #
#------------------------------------------------------------------------------------------#
joborder = data.frame(sapply(X=default,FUN=rep,times=nruns),stringsAsFactors=FALSE)
for (n in 1:nvars) joborder[[names(myruns)[n]]] = myruns[[names(myruns)[n]]]
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Create a table with room for all simulations.                                        #
#------------------------------------------------------------------------------------------#
poidata = poilist[match(joborder$iata,poilist$iata),]
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Replace met driver when it is supposed to be tower data.                             #
#------------------------------------------------------------------------------------------#
is.tower                      = joborder$met.driver == "tower"
joborder$met.driver[is.tower] = poidata$met.driver[is.tower]
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Replace other POI-specific variables as long as they are not to be specified by the  #
# user settings.                                                                           #
#------------------------------------------------------------------------------------------#
keep    = ( names(poidata) %in% names(joborder) 
          & ( ! names(poidata) %in% names(myruns) )
          & ( ! names(poidata) %in% c("iata","met.driver") )
          )#end keep
poidata = poidata[,keep]
npois   = ncol(poidata)
if (npois > 0){
   for (p in 1:npois) joborder[[names(poidata)[p]]] = poidata[[names(poidata)[p]]]
}#end if
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Make sure that soil texture is standard when the user provides soil texture and the  #
# value is not zero.                                                                       #
#------------------------------------------------------------------------------------------#
if ("istext" %in% names(varrun)){
   sel                    = joborder$istext != 0
   joborder$sand  [  sel] = -1.0
   joborder$clay  [  sel] = -1.0
   joborder$istext[! sel] = poidata$istext[! sel]
}#end if
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#    Build the name of the simulations.  The polygon name always stays, even if the run is #
# for one site only.                                                                       #
#------------------------------------------------------------------------------------------#
stay = which(names(varlabel) %in% c("iata"))
bye  = which(sapply(X=varlabel,FUN=length) == 1)
bye  = bye[! bye %in% stay]
if (length(bye) > 0) for (b in sort(bye,decreasing=TRUE)) varlabel[[b]] = NULL
runname  = apply( X        = expand.grid(varlabel, stringsAsFactors = FALSE)
                , MARGIN   = 1
                , FUN      = paste
                , collapse = "_"
                )#end apply
metname           = "s"
metname[is.tower] = "t"
joborder$run      = paste(metname,runname,sep="")
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Write the joborder table.                                                            #
#------------------------------------------------------------------------------------------#
dash       = sapply( FUN      = paste
                   , X        = mapply( FUN      = rep
                                      , times    = nchar(names(joborder))
                                      , MoreArgs = list(x="-")
                                      )#end mapply
                   , collapse = "")
jobneat    = format(rbind(dash,toupper(names(joborder)),dash,joborder),justify="right"
                   ,width=nchar(names(joborder)))
jobneat    = apply(X=jobneat,MARGIN=1,FUN=paste,collapse="  ")
jobneat[1] = gsub(pattern=" ",replacement="-",x=jobneat[1])
jobneat[3] = gsub(pattern=" ",replacement="-",x=jobneat[3])
dum        = write.table(x=jobneat,file=outfile,quote=FALSE,row.names=FALSE,col.names=FALSE)
#------------------------------------------------------------------------------------------#
