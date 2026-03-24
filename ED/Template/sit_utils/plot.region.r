#==========================================================================================#
#==========================================================================================#
#     Leave these commands at the beginning.  They will refresh the session.               #
#------------------------------------------------------------------------------------------#
rm(list=ls())
graphics.off()
gc()
options(warn=0)
#==========================================================================================#
#==========================================================================================#



#==========================================================================================#
#==========================================================================================#
#      Here is the user defined variable section.                                          #
#------------------------------------------------------------------------------------------#


#----- Paths. -----------------------------------------------------------------------------#
main    = "mypath"                                      # Main path.
here    = file.path(main,"sit_utils")                   # This path
srcdir  = c( "/home/mlongo/Util/Rsc"                    # Possible paths with libraries
           , "/Users/mlongo/Util/Rsc"                   #    R will select the first
           , "/prj/prjidfca/marcosl/Util/Rsc"           #    one that is found.
           , "/prj/bramsolam/marcos.longo/Util/Rsc"     #
           , "/scratch/bramsolam/marcos.longo/Util/Rsc" #
           , "/n/home00/mlongo/Util/Rsc"                #
           )#end c                                      #
outroot = here                                          # Directory for figures
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#      List files to be read.  The current and the previous check, and the path with the   #
# job order.                                                                               #
#------------------------------------------------------------------------------------------#
joborder  = file.path(main,"joborder.txt" )
lastcheck = file.path(here,"lastcheck.txt")
mycheck   = file.path(here,"mycheck.txt"  )
#------------------------------------------------------------------------------------------#




#----- Plot options. ----------------------------------------------------------------------#
outform  = c("png")             # Formats for output file.  Supported formats are:
                                #   - "X11" - for printing on screen
                                #   - "eps" - for postscript printing
                                #   - "png" - for PNG printing
                                #   - "pdf" - for PDF printing
depth    = 96                   # PNG resolution, in pixels per inch
paper    = "letter"             # Paper size, to define the plot shape
ptsz     = 16                   # Font size.
lwidth   = 2.5                  # Line width
inset    = 0.01                 # inset between legend and edge of plot region.
fracexp  = 0.40                 # Expand the y axis by this amount.
ncolours = 20                   # Number of colours to split the real variables
st.leg   = 1./5.
ey.leg   = 1./6.
#------------------------------------------------------------------------------------------#


#------ List of variables to make the maps. -----------------------------------------------#
n            = 0
outvars      = list()
n            = n + 1
outvars[[n]] = list( vnam = "lai"
                   , desc = "Leaf area index"
                   , unit = "m2lom2"
                   , csch = "pubugn"
                   )#end list
n            = n + 1
outvars[[n]] = list( vnam = "bsa"
                   , desc = "Basal area"
                   , unit = "cm2om2"
                   , csch = "bupu"
                   )#end list
n            = n + 1
outvars[[n]] = list( vnam = "agb"
                   , desc = "Aboveground carbon"
                   , unit = "kgcom2"
                   , csch = "ylgnbu"
                   )#end list
n            = n + 1
outvars[[n]] = list( vnam = "scb"
                   , desc = "Soil carbon"
                   , unit = "kgcom2"
                   , csch = "orrd"
                   )#end list
n            = n + 1
outvars[[n]] = list( vnam = "npa"
                   , desc = "Patch count"
                   , unit = "empty"
                   , csch = "magma"
                   )#end list
#------------------------------------------------------------------------------------------#



#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#      NO NEED TO CHANGE ANYTHING BEYOND THIS POINT UNLESS YOU ARE DEVELOPING THE CODE...  #
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#

#----- Check that directory main has been set. --------------------------------------------#
if (! dir.exists(main)){
   cat (" Main: ",main,"\n")
   stop(" Directory main is incorrect or has not been set!!!")
}#end if
#------------------------------------------------------------------------------------------#




#----- Loading some packages and scripts. -------------------------------------------------#
srcdir = (srcdir[file.exists(srcdir)])[1]
source(file.path(srcdir,"load.everything.r"))
#------------------------------------------------------------------------------------------#



#----- Set how many formats we must output. -----------------------------------------------#
outform = tolower(outform)
nout    = length (outform)
#------------------------------------------------------------------------------------------#



#----- Convert output var list to data table. ---------------------------------------------#
outvars  = list.2.data.table(outvars)
noutvars = nrow(outvars)
#------------------------------------------------------------------------------------------#



#---- Create the main output directory in case there is none. -----------------------------#
dummy = dir.create(outroot,recursive=TRUE,showWarnings=FALSE)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Read job order.  We always use this file to build the array.                         #
#------------------------------------------------------------------------------------------#
cat0(" + Read ",basename(joborder),".")
names.jobs  = scan(file=joborder,skip=1,nlines=1,what="character",quiet=TRUE)
names.jobs  = gsub(pattern="_",replacement=".",x=tolower(names.jobs))
jobs        = read.table(file             = joborder
                        ,skip             = 3
                        ,header           = FALSE
                        ,comment.char     = ""
                        ,stringsAsFactors = FALSE
                        )#end read.table
names(jobs) = names.jobs
jobs        = data.table(jobs)
njobs       = nrow(jobs)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Read the last and the current check.  For the current check, we normally skip the    #
# last line to avoid trouble, unless the file is complete.  Sometimes one of the files     #
# (typically last) may be missing.  Account for that too.                                  #
#------------------------------------------------------------------------------------------#
names.check   = c("run","lon","lat","year","month","day","hhmm","stall","runt"
                 ,"agb","bsa","lai","scb","npa")
if (file.exists(lastcheck) && file.exists(mycheck)){
   #------ Both files exist, read both. ---------------------------------------------------#
   cat0(" + Read ",basename(lastcheck),".")
   last = read.table( file             = lastcheck
                    , skip             = 0
                    , header           = FALSE
                    , comment.char     = ""
                    , col.names        = names.check
                    , stringsAsFactors = FALSE
                    )#end read.table
   cat0(" + Read ",basename(mycheck),".")
   ncurr = length(readLines(mycheck))
   if (ncurr == njobs){
      curr = read.table( file             = mycheck
                       , skip             = 0
                       , header           = FALSE
                       , comment.char     = ""
                       , col.names        = names.check
                       , stringsAsFactors = FALSE
                       )#end read.table
   }else if (ncurr > 0){
      curr = read.table( file             = mycheck
                       , skip             = 0
                       , nrows            = ncurr-1
                       , header           = FALSE
                       , comment.char     = ""
                       , col.names        = names.check
                       , stringsAsFactors = FALSE
                       )#end read.table
   }else{
      curr = last[ 1,,drop=FALSE]
      curr = curr[-1,,drop=FALSE]
   }#end if (ncurr == njobs)
   #---------------------------------------------------------------------------------------#
}else if(file.exists(lastcheck)){
   #------ Only lastcheck exists, duplicate it. -------------------------------------------#
   cat0(" + Read ",basename(lastcheck),".")
   last = read.table( file             = lastcheck
                    , skip             = 0
                    , header           = FALSE
                    , comment.char     = ""
                    , col.names        = names.check
                    , stringsAsFactors = FALSE
                    )#end read.table
   #---------------------------------------------------------------------------------------#

   #---- Duplicate last. ------------------------------------------------------------------#
   curr = last
   #---------------------------------------------------------------------------------------#
}else if(file.exists(mycheck  )){
   #------ Only mycheck exists. -----------------------------------------------------------#
   cat0(" + Read ",basename(mycheck),".")
   ncurr = length(readLines(mycheck))
   if (ncurr == njobs){
      curr = read.table( file             = mycheck
                       , skip             = 0
                       , header           = FALSE
                       , comment.char     = ""
                       , col.names        = names.check
                       , stringsAsFactors = FALSE
                       )#end read.table
   }else if (ncurr > 0){
      curr = read.table( file             = mycheck
                       , skip             = 0
                       , nrows            = ncurr-1
                       , header           = FALSE
                       , comment.char     = ""
                       , col.names        = names.check
                       , stringsAsFactors = FALSE
                       )#end read.table
   }else{
      #------ Nothing to read, stop the run. ----------------------------------------------#
      cat0("-----------------------------------------------------------------------------")
      cat0("    Last check file not found and current check file is empty/almost empty!")
      cat0(" Lastcheck: "              ,lastcheck,".")
      cat0(" Mycheck:   "              ,mycheck  ,".")
      cat0(" Line count in mycheck:   ",ncurr    ,".")
      cat0("-----------------------------------------------------------------------------")
      stop(" At least one of these files must exist (ideally both of them).")
      #------------------------------------------------------------------------------------#
   }#end if (ncurr == njobs)
   #---------------------------------------------------------------------------------------#

   #---- Duplicate curr. ------------------------------------------------------------------#
   last = curr
   #---------------------------------------------------------------------------------------#
}else{
   #------ Nothing to read, stop the run. -------------------------------------------------#
   cat0("-----------------------------------------------------------------------------")
   cat0("    None of the check files was found!")
   cat0(" Lastcheck: ",lastcheck,".")
   cat0(" Mycheck:   ",mycheck  ,".")
   cat0("-----------------------------------------------------------------------------")
   stop(" At least one of these files must exist (ideally both of them).")
   #---------------------------------------------------------------------------------------#
}#end if (file.exists(lastcheck) && file.exists(mycheck))
#------------------------------------------------------------------------------------------#


#----- Turn data frames into data tables. -------------------------------------------------#
last = data.table(last)
curr = data.table(curr)
#------------------------------------------------------------------------------------------#





#------------------------------------------------------------------------------------------#
#     Find the dimensions.                                                                 #
#------------------------------------------------------------------------------------------#
lon             = unique(sort(jobs$lon))
lat             = unique(sort(jobs$lat))
dlon            = median(diff(lon))
dlat            = median(diff(lat))
n.lon           = length(lon)
n.lat           = length(lat)
key.lon         = sprintf("%+06.2f",lon)
key.lat         = sprintf("%+06.2f",lat)
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Initialise the variables.                                                            #
#------------------------------------------------------------------------------------------#
r.template  = array( data     = NA_real_
                  , dim      = c(n.lon,n.lat)
                  , dimnames = list(key.lon,key.lat)
                  )#end array
i.template  = array( data     = NA_integer_
                  , dim      = c(n.lon,n.lat)
                  , dimnames = list(key.lon,key.lat)
                  )#end array
datum      = list ( agb    = r.template
                  , lai    = r.template
                  , bsa    = r.template
                  , scb    = r.template
                  , npa    = i.template
                  , status = i.template
                  , yearn  = i.template
                  )#end list
#------------------------------------------------------------------------------------------#



#----- Define plot window size ------------------------------------------------------------#
limlon  = range(lon)
limlat  = range(lat)
st.ext  = st.leg / (1. - st.leg)
st.size = plotsize( proje     = TRUE
                  , limlon    = limlon
                  , limlat    = limlat
                  , extendfc  = "lat"
                  , extfactor = st.ext
                  , paper     = paper
                  )
ey.ext  = ey.leg / (1. - ey.leg)
ey.size = plotsize( proje     = TRUE
                  , limlon    = limlon
                  , limlat    = limlat
                  , extendfc  = "lat"
                  , extfactor = ey.ext
                  , paper     = paper
                  )
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Copy the information from last check to joborder.                                    #
#------------------------------------------------------------------------------------------#
jobs$status = rep("INITIAL",times=njobs)
jobs$yearn  = jobs$yeara
jobs$agb    = rep(NA_real_   ,times=njobs)
jobs$bsa    = rep(NA_real_   ,times=njobs)
jobs$lai    = rep(NA_real_   ,times=njobs)
jobs$scb    = rep(NA_real_   ,times=njobs)
jobs$npa    = rep(NA_integer_,times=njobs)

il = match(last$run,jobs$run); l.sel = ! is.na(il)
ic = match(curr$run,jobs$run); c.sel = ! is.na(ic)
jobs$status[il[l.sel]] = last$runt[l.sel]  ; jobs$status[ic[c.sel]] = curr$runt[c.sel]
jobs$agb   [il[l.sel]] = last$agb [l.sel]  ; jobs$agb   [ic[c.sel]] = curr$agb [c.sel]
jobs$lai   [il[l.sel]] = last$lai [l.sel]  ; jobs$lai   [ic[c.sel]] = curr$lai [c.sel]
jobs$bsa   [il[l.sel]] = last$bsa [l.sel]  ; jobs$bsa   [ic[c.sel]] = curr$bsa [c.sel]
jobs$scb   [il[l.sel]] = last$scb [l.sel]  ; jobs$scb   [ic[c.sel]] = curr$scb [c.sel]
jobs$npa   [il[l.sel]] = last$npa [l.sel]  ; jobs$npa   [ic[c.sel]] = curr$npa [c.sel]
jobs$yearn [il[l.sel]] = last$year[l.sel]  ; jobs$yearn [ic[c.sel]] = curr$year[c.sel]
#------------------------------------------------------------------------------------------#




#----- Select layers that are used for plotting. ------------------------------------------#
keep  = names(jobs) %in% c("lon","lat","iata","yeara","yearz","yearn","status"
                           ,"agb","bsa","lai","scb","npa")
jobs  = jobs[,..keep,drop=FALSE]
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Find the indices to map the data to the arrays.                                      #
#------------------------------------------------------------------------------------------#
i.lon               = match(jobs$lon,lon)
i.lat               = match(jobs$lat,lat)
index               = cbind(i.lon,i.lat)
datum$agb   [index] = jobs$agb
datum$lai   [index] = jobs$lai
datum$bsa   [index] = jobs$bsa
datum$scb   [index] = jobs$scb
datum$npa   [index] = jobs$npa
datum$status[index] = jobs$status
datum$yearn [index] = jobs$yearn
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Run the matrices.                                                                    #
#------------------------------------------------------------------------------------------#
yeara    = min(jobs$yeara)
yearz    = max(jobs$yearz)
yr.range = pretty.xylim(u=c(yeara,yearz))
yr.brks  = pretty(yr.range,n=10)
yr.keep  = (yr.brks %wr% yr.range) & (! yr.brks %in% c(yeara,yearz))
yr.brks  = unique(c(-Inf,yr.brks[yr.keep],Inf))
n.cut    = length(yr.brks)-1
yr.cut   = cut(datum$yearn,yr.brks)
yr.level = levels(yr.cut)
n.level  = length(yr.level)
#------------------------------------------------------------------------------------------#



#------ Find the year colours. ------------------------------------------------------------#
datum$yr.idx          = match(yr.cut,yr.level) + 0 * datum$yearn
initial               = datum$status == "INITIAL"
crashed               = datum$status == "CRASHED"
hydfail               = datum$status == "HYDFAIL"
bad.met               = datum$status == "BAD_MET"
metmiss               = datum$status == "METMISS"
stopped               = datum$status == "STOPPED"
extinct               = datum$status == "EXTINCT"
ststate               = datum$status == "STSTATE"
the.end               = datum$status == "THE_END"
datum$yr.idx[initial] = 0
datum$yr.idx[the.end] = n.level + 1
datum$yr.idx[ststate] = n.level + 2
datum$yr.idx[extinct] = n.level + 3
datum$yr.idx[stopped] = n.level + 4
datum$yr.idx[metmiss] = n.level + 5
datum$yr.idx[bad.met] = n.level + 6
datum$yr.idx[hydfail] = n.level + 7
datum$yr.idx[crashed] = n.level + 8
yr.cscheme            = c("grey89",atlas(n=n.level),"royalblue4","steelblue3","purple3"
                         ,"mediumpurple1","deepskyblue","hotpink","red3","firebrick4")
ybottom               = rep(0,times=n.level+9)
ytop                  = rep(1,times=n.level+9)
xleft                 = sequence(n.level+9) - 1
xright                = sequence(n.level+9)
xat                   = c(sequence(n.level+2)-1,n.level+sequence(8)+0.5)
xbrks                 = sequence(n.level+10)-1-0.5
xlabel                = c("Initial",yeara,yr.brks[-c(1,n.level+1)],yearz,"Finish","StState"
                         ,"Extinct","Stopped","MetMiss","Bad Met","HydFail","Crashed")
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Limits and labels for longitude and latitude.                                       #
#------------------------------------------------------------------------------------------#
limlon  = c(min(lon)-0.5*dlon,max(lon)+0.5*dlon)
limlat  = c(min(lat)-0.5*dlat,max(lat)+0.5*dlat)
lonplot = pretty.lonlat(x=limlon,n=6,type="lon")
latplot = pretty.lonlat(x=limlat,n=6,type="lat")
#------------------------------------------------------------------------------------------#




#==========================================================================================#
#==========================================================================================#
#      Create the status map for all sites.                                                #
#------------------------------------------------------------------------------------------#
cat0(" Plot the current status.")

   #------ Make plot annotation. ----------------------------------------------------------#
   letitre      = "Polygon status"
   #---------------------------------------------------------------------------------------#

   for (o in sequence(nout)){
      #----- Open file. -------------------------------------------------------------------#
      fichier = file.path(here,paste0("stt_region.",outform[o]))
      dummy   = open.plot( fichier = fichier
                         , outform = outform[o]
                         , size    = st.size
                         , ptsz    = ptsz
                         , depth   = depth
                         )#end open.plot
      #------------------------------------------------------------------------------------#


      #----- Load settings and split window. ----------------------------------------------#
      par.orig = par(par.user)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Split the plotting window.                                                    #
      #------------------------------------------------------------------------------------#
      layout(mat = rbind(2,1), heights = c(1.-st.leg,st.leg))
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     First, let's plot the legend.                                                  #
      #------------------------------------------------------------------------------------#
      par(mar=c(3.1,3.1,1.1,1.1))
      plot.new()
      plot.window( xlim = range(xleft,xright)
                 , ylim = range(ybottom,ytop)
                 , xaxs = "i"
                 , yaxs = "i"
                 )#end plot.window
      rect( xleft   = xleft
          , ybottom = ybottom
          , xright  = xright
          , ytop    = ytop
          , col     = yr.cscheme
          , border  = "grey50"
          )#end rect
      box()
      axis(side=1,at=xat,srt=45,labels=FALSE)
      text(x=xat,y=par("usr")[3]-0.2,labels=xlabel,srt=45,adj=1,xpd=TRUE,cex=0.95)
      title(main="Status",cex.main=1.0)
      #------------------------------------------------------------------------------------#





      #----- Set the window. --------------------------------------------------------------#
      par(mar = c(2.1,3.1,1.1,1.1))
      plot.new()
      plot.window(xlim=limlon,ylim=limlat,xaxs="i",yaxs="i")
      axis(side=1,las=1,at=lonplot$at,labels=lonplot$labels)
      axis(side=2,las=1,at=latplot$at,labels=latplot$labels)
      box()
      image( x      = lon
           , y      = lat
           , z      = datum$yr.idx
           , col    = yr.cscheme
           , breaks = xbrks
           , add    = TRUE
           )#end image
      southammap(col="grey30",lwd=1)
      amazonmap (col="black" ,lwd=2)
      #------------------------------------------------------------------------------------#



      #----- Close the device. ------------------------------------------------------------#
      dummy = close.plot(outform = outform[o])
      #------------------------------------------------------------------------------------#
   }#end for (o in sequence(nout))
   #---------------------------------------------------------------------------------------#
#==========================================================================================#
#==========================================================================================#




#------------------------------------------------------------------------------------------#
#      Create parameter space maps for all other variables.                                #
#------------------------------------------------------------------------------------------#
cat0(" + Plot the current properties.")
for (v in sequence(noutvars)){
   #----- Handy aliases. ------------------------------------------------------------------#
   v.vnam = outvars$vnam[v]
   v.desc = outvars$desc[v]
   v.unit = untab[[outvars$unit[v]]]
   v.csch = match.fun(outvars$csch[v])
   cat0("   - ",v.desc,".")
   #---------------------------------------------------------------------------------------#


   #----- Collapse realisations using the median. -----------------------------------------#
   v.value      = datum[[v.vnam]]
   del          = ! is.finite(v.value)
   v.value[del] = NA_real_
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Break the data into bins.                                                        #
   #---------------------------------------------------------------------------------------#
   v.limit   = pretty.xylim(v.value)
   v.brks    = pretty(v.value,n=ncolours)
   n.brks    = length(v.brks)
   v.cut     = cut(as.numeric(v.value),breaks=v.brks)
   v.lev     = levels(v.cut)
   v.idx     = match(v.cut,v.lev) + 0 * v.value
   v.cscheme = v.csch(n=n.brks-1)
   #---------------------------------------------------------------------------------------#


   #----- Make the edges. -----------------------------------------------------------------#
   xleft       = v.brks[-n.brks]
   xright      = v.brks[     -1]
   ybottom     = rep(0,times=n.brks)
   ytop        = rep(1,times=n.brks)
   xat         = v.brks
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Loop through formats.                                                             #
   #---------------------------------------------------------------------------------------#
   for (o in sequence(nout)){
      #----- Open file. -------------------------------------------------------------------#
      fichier = file.path(here,paste0(v.vnam,"_region.",outform[o]))
      dummy   = open.plot( fichier = fichier
                         , outform = outform[o]
                         , size    = ey.size
                         , ptsz    = ptsz
                         , depth   = depth
                         )#end open.plot
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Split the plotting window.                                                    #
      #------------------------------------------------------------------------------------#
      par(par.user)
      layout(mat=rbind(2,1),heights=c(1.-ey.leg,ey.leg))
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     First, let's plot the legend.                                                  #
      #------------------------------------------------------------------------------------#
      par(mar=c(2.1,3.1,1.6,1.1))
      plot.new()
      plot.window(xlim=range(xleft,xright),ylim=range(ybottom,ytop),xaxs="i",yaxs="i")
      rect( xleft   = xleft
          , ybottom = ybottom
          , xright  = xright
          , ytop    = ytop
          , col     = v.cscheme
          , border  = "transparent"
          )#end rect
      box()
      axis(side=1,at=xat)
      title(main=desc.unit(desc=v.desc,unit=v.unit),cex.main=1.0)
      #------------------------------------------------------------------------------------#



      #----- Set the window. --------------------------------------------------------------#
      par(mar = c(2.1,3.1,1.1,1.1))
      plot.new()
      plot.window(xlim=limlon,ylim=limlat,xaxs="i",yaxs="i")
      axis(side=1,las=1,at=lonplot$at,labels=lonplot$labels)
      axis(side=2,las=1,at=latplot$at,labels=latplot$labels)
      box()
      image(x=lon,y=lat,z=v.value,col=v.cscheme,breaks=v.brks,add=TRUE)
      southammap(col="grey30",lwd=1)
      amazonmap (col="black" ,lwd=2)
      #------------------------------------------------------------------------------------#



      #----- Close the device. ------------------------------------------------------------#
      dummy = close.plot(outform=outform[o])
      #------------------------------------------------------------------------------------#
   }#end for (o in sequence(nout))
   #---------------------------------------------------------------------------------------#
}#end for (v in sequence(noutvars))
#------------------------------------------------------------------------------------------#
