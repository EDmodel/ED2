#==========================================================================================#
#==========================================================================================#
#     Leave these commands at the beginning.  They will refresh the session.               #
#------------------------------------------------------------------------------------------#
rm(list=ls())
graphics.off()
#==========================================================================================#
#==========================================================================================#



#==========================================================================================#
#==========================================================================================#
#      Here is the user defined variable section.                                          #
#------------------------------------------------------------------------------------------#


#----- Paths. -----------------------------------------------------------------------------#
main    = "mypath"                                      # Main path (one up sit_utils).
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
outform        = c("png","eps","pdf") # Formats for output file.  Supported formats are:
                                      #   - "X11" - for printing on screen
                                      #   - "eps" - for postscript printing
                                      #   - "png" - for PNG printing
                                      #   - "pdf" - for PDF printing
depth          = 96                   # PNG resolution, in pixels per inch
paper          = "letter"             # Paper size, to define the plot shape
ptsz           = 14                   # Font size.
lwidth         = 2.5                  # Line width
inset          = 0.01                 # inset between legend and edge of plot region.
fracexp        = 0.40                 # Expand the y axis by this amount.
ncolours       = 20                   # Number of colours to split the real variables
mtext.xoff     = -8.50                # Offset for the x label
mtext.yoff     = -1.00                # Offset for the y label
mtext.xadj     =  0.50                # Offset for the x label
mtext.yadj     =  0.65                # Offset for the y label
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



#----- Convert output var list to data table. ---------------------------------------------#
outvars  = list.2.data.table(outvars)
noutvars = nrow(outvars)
#------------------------------------------------------------------------------------------#



#----- Set how many formats we must output. -----------------------------------------------#
outform = tolower(outform)
nout    = length (outform)
#------------------------------------------------------------------------------------------#



#----- Define plot window size ------------------------------------------------------------#
size = plotsize(proje=FALSE,paper=paper)
#------------------------------------------------------------------------------------------#



#---- Create the main output directory in case there is none. -----------------------------#
dummy = dir.create(outroot,recursive=TRUE,showWarnings=FALSE)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Read job order.  We always use this file to build the array.                         #
#------------------------------------------------------------------------------------------#
cat0(" + Read ",basename(joborder),".")
names.jobs       = scan(file=joborder,skip=1,nlines=1,what="character",quiet=TRUE)
names.jobs       = gsub(pattern="_",replacement=".",x=tolower(names.jobs))
jobs             = read.table(file=joborder,skip=3,header=FALSE,comment.char=""
                             ,stringsAsFactors=FALSE)
names(jobs)      = names.jobs
njobs            = nrow(jobs)
jobs$drain       = as.numeric(substring(jobs$run, 7,10)) / 100
jobs$dtemp       = as.numeric(substring(jobs$run,13,16)) / 100
jobs$realisation = as.numeric(substring(jobs$run,23,24))
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
iata            = unique(sort(jobs$iata       ))
drain           = unique(sort(jobs$drain      ))
dtemp           = unique(sort(jobs$dtemp      ))
realisation     = unique(sort(jobs$realisation))
iphen           = unique(sort(jobs$iphen      ))
stext           = unique(sort(jobs$istext     ))
n.iata          = length(iata       )
n.drain         = length(drain      )
n.dtemp         = length(dtemp      )
n.realisation   = length(realisation)
n.stext         = length(stext      )
n.iphen         = length(iphen      )
key.iata        = toupper(iata      )
key.drain       = paste0("r"    ,sprintf("%+3.3i",sort(100*unique(jobs$drain      ))))
key.dtemp       = paste0("t"    ,sprintf("%+3.3i",sort(100*unique(jobs$dtemp      ))))
key.realisation = paste0("real" ,sprintf("%2.2i" ,sort(    unique(jobs$realisation))))
key.stext       = paste0("stext",sprintf("%2.2i" ,sort(    unique(jobs$istext     ))))
key.iphen       = c("Evergreen","Deciduous")
desc.iata       = poilist$longname[match(iata,poilist$iata)]
desc.stext      = stext.names[stext]
desc.iphen      = c("Evergreen","Drought deciduous")
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Initialise the variables.                                                            #
#------------------------------------------------------------------------------------------#
template   = array( data     = NA
                  , dim      = c(n.drain,n.dtemp,n.realisation,n.iphen,n.iata,n.stext)
                  , dimnames = list(key.drain,key.dtemp,key.realisation,key.iphen
                                   ,key.iata,key.stext)
                  )#end array
datum      = list ( agb    = template
                  , lai    = template
                  , bsa    = template
                  , scb    = template
                  , npa    = template
                  , status = template
                  , yearn  = template
                  )#end list
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Copy the information from last check to joborder.                                    #
#------------------------------------------------------------------------------------------#
jobs$status = rep("INITIAL",times=njobs)
jobs$yearn  = jobs$yeara
jobs$agb    = rep(NA,times=njobs)
jobs$bsa    = rep(NA,times=njobs)
jobs$lai    = rep(NA,times=njobs)
jobs$scb    = rep(NA,times=njobs)
jobs$npa    = rep(NA,times=njobs)

il = match(last$run,jobs$run); l.sel = ! is.na(il)
ic = match(curr$run,jobs$run); c.sel = ! is.na(ic)
jobs$status[il[l.sel]] = last$runt[l.sel]  ; jobs$status[ic[c.sel]] = curr$runt[c.sel]
jobs$agb   [il[l.sel]] = last$agb [l.sel]  ; jobs$agb   [ic[c.sel]] = curr$agb [c.sel]
jobs$lai   [il[l.sel]] = last$lai [l.sel]  ; jobs$lai   [ic[c.sel]] = curr$lai [c.sel]
jobs$bsa   [il[l.sel]] = last$bsa [l.sel]  ; jobs$bsa   [ic[c.sel]] = curr$bsa [c.sel]
jobs$scb   [il[l.sel]] = last$scb [l.sel]  ; jobs$scb   [ic[c.sel]] = curr$scb [c.sel]
jobs$npa   [il[l.sel]] = last$npa [l.sel]  ; jobs$npa   [ic[c.sel]] = curr$npa [c.sel]
jobs$yearn [il[l.sel]] = last$year[l.sel]  ; jobs$yearn [ic[c.sel]] = curr$year[c.sel]





keep            = names(jobs) %in% c("drain","dtemp","realisation","iphen","iata","istext"
                                    ,"yeara","yearz","yearn","status"
                                    ,"agb","bsa","lai","scb","npa")
jobs            = jobs[,keep]
weird           = is.finite(jobs$lai) & abs(jobs$lai) > 20
jobs$lai[weird] = NA
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Find the indices to map the data to the arrays.                                      #
#------------------------------------------------------------------------------------------#
i.drain             = match(jobs$drain      ,drain      )
i.dtemp             = match(jobs$dtemp      ,dtemp      )
i.realisation       = match(jobs$realisation,realisation)
i.iphen             = match(jobs$iphen      ,iphen      )
i.iata              = match(jobs$iata       ,iata       )
i.stext             = match(jobs$istext     ,stext      )
index               = cbind(i.drain,i.dtemp,i.realisation,i.iphen,i.iata,i.stext)
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
#    Pick a random realisation to send to the user.                                        #
#------------------------------------------------------------------------------------------#
use.yr.idx = apply( X      = datum$yr.idx
                  , FUN    = max
                  , MARGIN = c(1,2,4,5,6)
                  , na.rm  = TRUE
                  )#end use.idx
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Create the status map for all sites.                                                #
#------------------------------------------------------------------------------------------#
cat0(" + lot the current status.")
lo.box = pretty.box(n.iphen*n.iata)
for (st in sequence(n.stext)){
   #------ Find the soil texture key and description. -------------------------------------#
   this.st.key  = key.stext [st]
   letitre      = paste("Polygon status   --   Soil type: ",desc.stext[st]
                       ,sep="")
   lex          = paste("Rainfall change [scale parameter]")
   ley          = paste("Temperature change [K]")
   #---------------------------------------------------------------------------------------#

   for (o in sequence(nout)){
      #----- Open file. -------------------------------------------------------------------#
      fichier = file.path(here,paste0("stt_",this.st.key,".",outform[o]))
      dummy   = open.plot( fichier = fichier
                         , outform = outform[o]
                         , size    = size
                         , ptsz    = ptsz
                         , depth   = depth
                         )#end open.plot
      #------------------------------------------------------------------------------------#


      #----- Save the margins to avoid losing the data. -----------------------------------#
      par(par.user)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Split the plotting window.                                                    #
      #------------------------------------------------------------------------------------#
      par(oma = c(0.2,3,4,0))
      layout( mat     = rbind(lo.box$mat+1,rep(1,times=lo.box$ncol))
            , heights = c(rep(5/lo.box$nrow,lo.box$nrow),1)
            )#end layout
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     First, let's plot the legend.                                                  #
      #------------------------------------------------------------------------------------#
      par(mar=c(3,3,2,3)+0.1)
      plot.new()
      plot.window(xlim=range(xleft,xright),ylim=range(ybottom,ytop),xaxs="i",yaxs="i")
      rect(xleft=xleft,ybottom=ybottom,xright=xright,ytop=ytop,col=yr.cscheme)
      box()
      axis(side=1,at=xat,srt=45,labels=FALSE)
      text(x=xat,y=par("usr")[3]-0.6,labels=xlabel,srt=45,adj=1,xpd=TRUE,cex=0.95)
      title(main="Status",ylab="",xlab="")
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Now we loop over sites and pheonologies.                                       #
      #------------------------------------------------------------------------------------#
      k = 0
      for (ph in sequence(n.iphen)){
         for (pl in sequence(n.iata)){
            #----- Make the sub-title. ----------------------------------------------------#
            lesub=paste(desc.iata[pl],desc.iphen[ph],sep=" - ")
            #------------------------------------------------------------------------------#

            #----- Find out where is this box going, and set up axes and margins. ---------#
            k       = k + 1
            left    = (k %% lo.box$ncol) == 1 || lo.box$ncol == 1
            right   = (k %% lo.box$ncol) == 0
            top     = k <= lo.box$ncol
            bottom  = k > (lo.box$nrow - 1) * lo.box$ncol
            mar.now = c(3 + 1 * bottom,1 + 1 * left,1 + 1 * top,1 + 1 * right) + 0.1
            if (bottom){
               xaxt = "s"
            }else{
               xaxt = "n"
            }#end if
            if (left){
               yaxt = "s"
            }else{
               yaxt = "n"
            }#end if
            #------------------------------------------------------------------------------#



            #----- Set the window. --------------------------------------------------------#
            par(mar = mar.now)
            image(x=drain,y=dtemp,z=use.yr.idx[,,ph,pl,st],col=yr.cscheme,breaks=xbrks
                 ,xaxt=xaxt,yaxt=yaxt,main=lesub,xlab="",ylab="")
            #------------------------------------------------------------------------------#
         }#end for
         #---------------------------------------------------------------------------------#
      }#end for
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Plot the global title.                                                         #
      #------------------------------------------------------------------------------------#
      par(las=0)
      mtext(side=1,text=lex    ,outer=TRUE,adj=mtext.xadj,padj=mtext.xoff)
      mtext(side=2,text=ley    ,outer=TRUE,adj=mtext.yadj,padj=mtext.yoff)
      mtext(side=3,text=letitre,outer=TRUE,font=2)
      #------------------------------------------------------------------------------------#



      #----- Close the device. ------------------------------------------------------------#
      dummy = close.plot(outform=outform[o])
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#
}#end for
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#      Create parameter space maps for all other variables.                                #
#------------------------------------------------------------------------------------------#
cat0(" Plot the current properties.")
lo.box = pretty.box(n.iphen*n.iata)
for (v in sequence(noutvars)){
   #----- Handy aliases. ------------------------------------------------------------------#
   v.vnam = outvars$vnam[v]
   v.desc = outvars$desc[v]
   v.unit = untab[[outvars$unit[v]]]
   v.csch = match.fun(outvars$csch[v])
   cat0("   - ",v.desc,".")
   #---------------------------------------------------------------------------------------#


   #----- Collapse realisations using the median. -----------------------------------------#
   v.value       = apply( X      = datum[[v.vnam]]
                        , MARGIN = c(1,2,4,5,6)
                        , FUN    = median
                        , na.rm  = TRUE
                        )#end apply
   rien          = ! is.finite(v.value)
   v.value[rien] = NA
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
   xleft       = var.brks[-n.brks]
   xright      = var.brks[     -1]
   ybottom     = rep(0,times=n.brks)
   ytop        = rep(1,times=n.brks)
   xat         = v.brks
   #---------------------------------------------------------------------------------------#


   for (st in sequence(n.stext)){
      #------ Find the soil texture key and description. ----------------------------------#
      this.st.key  = key.stext [st]
      letitre      = paste0(desc.var[v],"   --   Soil type: ",desc.stext[st])
      lex          = paste0("Rainfall change [scale parameter]")
      ley          = paste0("Temperature change [K]")
      #------------------------------------------------------------------------------------#

      for (o in sequence(nout)){
         #----- Open file. ----------------------------------------------------------------#
         fichier = file.path(here,paste0(key.var[v],"_",this.st.key,".",outform[o]))
         dummy   = open.plot( fichier = fichier
                            , outform = outform[o]
                            , size    = size
                            , ptsz    = ptsz
                            , depth   = depth
                            )#end open.plot
         #---------------------------------------------------------------------------------#


         #----- Save the margins to avoid losing the data. --------------------------------#
         par.orig = par(no.readonly = TRUE)
         mar.orig = par.orig$mar
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Split the plotting window.                                                 #
         #---------------------------------------------------------------------------------#
         par(oma = c(0.2,3,4,0))
         layout( mat     = rbind(lo.box$mat+1,rep(1,times=lo.box$ncol))
               , heights = c(rep(5/lo.box$nrow,lo.box$nrow),1)
               )#end layout
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     First, let's plot the legend.                                               #
         #---------------------------------------------------------------------------------#
         par(mar=c(3,3,2,3)+0.1)
         plot.new()
         plot.window(xlim=range(xleft,xright),ylim=range(ybottom,ytop),xaxs="i",yaxs="i")
         rect(xleft=xleft,ybottom=ybottom,xright=xright,ytop=ytop,col=v.cscheme)
         box()
         axis(side=1,at=xat)
         title(main=desc.unit(desc=v.desc,unit=v.unit),cex.main=1.0)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Now we loop over sites and pheonologies.                                    #
         #---------------------------------------------------------------------------------#
         k = 0
         for (ph in sequence(n.iphen)){
            for (pl in sequence(n.iata)){
               #----- Make the sub-title. -------------------------------------------------#
               lesub=paste(desc.iata[pl],desc.iphen[ph],sep=" - ")
               #---------------------------------------------------------------------------#

               #----- Find out where is this box going, and set up axes and margins. ------#
               k       = k + 1
               left    = (k %% lo.box$ncol) == 1 || lo.box$ncol == 1
               right   = (k %% lo.box$ncol) == 0
               top     = k <= lo.box$ncol
               bottom  = k > (lo.box$nrow - 1) * lo.box$ncol
               mar.now = c(3 + 1 * bottom,1 + 1 * left,1 + 1 * top,1 + 1 * right) + 0.1
               if (bottom){
                  xaxt = "s"
               }else{
                  xaxt = "n"
               }#end if
               if (left){
                  yaxt = "s"
               }else{
                  yaxt = "n"
               }#end if
               #---------------------------------------------------------------------------#



               #----- Set the window. -----------------------------------------------------#
               par(mar = mar.now)
               image(x=drain,y=dtemp,z=v.value[,,ph,pl,st],col=v.cscheme
                    ,breaks=var.brks,xaxt=xaxt,yaxt=yaxt,main=lesub,xlab="",ylab="")
               #---------------------------------------------------------------------------#
            }#end for
            #------------------------------------------------------------------------------#
         }#end for
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Plot the global title.                                                      #
         #---------------------------------------------------------------------------------#
         par(las=0)
         mtext(side=1,text=lex    ,outer=TRUE,adj=mtext.xadj,padj=mtext.xoff)
         mtext(side=2,text=ley    ,outer=TRUE,adj=mtext.yadj,padj=mtext.yoff)
         mtext(side=3,text=letitre,outer=TRUE,font=2)
         #---------------------------------------------------------------------------------#



         #----- Close the device. ---------------------------------------------------------#
         dummy = close.plot(outform=outform[o])
         #---------------------------------------------------------------------------------#
      }#end for (o in sequence(nout))
      #------------------------------------------------------------------------------------#
   }#end for (st in sequence(n.stext))
   #---------------------------------------------------------------------------------------#
}#end for (v in sequence(noutvars))
#------------------------------------------------------------------------------------------#
