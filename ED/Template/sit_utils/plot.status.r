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
main    = "/x/xxxxxxxxxxxx/xxxxxx/xxxxxxx/xxxxxxxx" # Main simulation directory.
here    = file.path(main,"sit_utils")               # This directory.

srcdir  = "/n/home00/mlongo/util/Rsc"               # Source  directory.
outroot = here                                      # Directory for figures
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
if ( main == "/x/xxxxxxxxxxxx/xxxxxx/xxxxxxx/xxxxxxxx"){
   cat (" Main: ",main,"\n")
   stop(" Directory main has not been set!!!")
}#end if
#------------------------------------------------------------------------------------------#



#----- Loading some packages and scripts. -------------------------------------------------#
source(file.path(srcdir,"load.everything.r"))
#------------------------------------------------------------------------------------------#



#----- Set how many formats we must output. -----------------------------------------------#
outform = tolower(outform)
nout    = length (outform)
#------------------------------------------------------------------------------------------#


#----- Avoid unecessary and extremely annoying beeps. -------------------------------------#
options(locatorBell=FALSE)
#------------------------------------------------------------------------------------------#



#----- Define plot window size ------------------------------------------------------------#
size = plotsize(proje=FALSE,paper=paper)
#------------------------------------------------------------------------------------------#



#---- Create the main output directory in case there is none. -----------------------------#
if (! file.exists(outroot)) dir.create(outroot)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Read job order.  We always use this file to build the array.                         #
#------------------------------------------------------------------------------------------#
cat (" + Reading ",basename(joborder),"...","\n")
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
# last line to avoid trouble, unless the file is complete.                                 #
#------------------------------------------------------------------------------------------#
names.check   = c("run","lon","lat","year","month","day","hhmm","runt"
                 ,"agb","bsa","lai","vel")
cat (" + Reading ",basename(lastcheck),"...","\n")
last          = read.table(file=lastcheck,skip=0,header=FALSE,comment.char=""
                          ,col.names=names.check,stringsAsFactors=FALSE)
cat (" + Reading ",basename(mycheck),"...","\n")
ncurr         = length(readLines(mycheck))
if (ncurr == njobs){
   curr = read.table(file=mycheck,skip=0,header=FALSE,comment.char=""
                    ,col.names=names.check,stringsAsFactors=FALSE)
}else if (ncurr > 0){
   curr = read.table(file=mycheck,skip=0,nrows=ncurr-1,header=FALSE,comment.char=""
                    ,col.names=names.check,stringsAsFactors=FALSE)
}else{
   curr = data.frame(rep(NA,times=length(names.check)),names=names.check)
   curr = curr[-1,]
}#end f
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
key.drain       = paste("r"    ,sprintf("%+3.3i",sort(100*unique(jobs$drain      ))),sep="")
key.dtemp       = paste("t"    ,sprintf("%+3.3i",sort(100*unique(jobs$dtemp      ))),sep="")
key.realisation = paste("real" ,sprintf("%2.2i" ,sort(    unique(jobs$realisation))),sep="")
key.stext       = paste("stext",sprintf("%2.2i" ,sort(    unique(jobs$istext     ))),sep="")
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

il = match(last$run,jobs$run); l.sel = ! is.na(il)
ic = match(curr$run,jobs$run); c.sel = ! is.na(ic)
jobs$status[il[l.sel]] = last$runt[l.sel]  ; jobs$status[ic[c.sel]] = curr$runt[c.sel]
jobs$agb   [il[l.sel]] = last$agb [l.sel]  ; jobs$agb   [ic[c.sel]] = curr$agb [c.sel]
jobs$lai   [il[l.sel]] = last$lai [l.sel]  ; jobs$lai   [ic[c.sel]] = curr$lai [c.sel]
jobs$bsa   [il[l.sel]] = last$bsa [l.sel]  ; jobs$bsa   [ic[c.sel]] = curr$bsa [c.sel]
jobs$yearn [il[l.sel]] = last$year[l.sel]  ; jobs$yearn [ic[c.sel]] = curr$year[c.sel]





keep            = names(jobs) %in% c("drain","dtemp","realisation","iphen","iata","istext"
                                    ,"yeara","yearz","yearn","status","agb","bsa","lai")
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
datum$status[index] = jobs$status
datum$yearn [index] = jobs$yearn
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Run the matrices.                                                                    #
#------------------------------------------------------------------------------------------#
yr.range              = range(c(jobs$yeara,jobs$yearz))
yr.cut                = pretty(yr.range,n=10)
yr.keep               = yr.cut > yr.range[1] & yr.cut < yr.range[2]
yr.brks               = c(-Inf,c(yr.range[1],yr.cut[yr.keep],yr.range[2]))
n.cut                 = length(yr.cut)-1
yr.cut                = cut(datum$yearn,yr.brks)
yr.level              = levels(yr.cut)
n.level               = length(yr.level)
#------------------------------------------------------------------------------------------#



#------ Find the year colours. ------------------------------------------------------------#
datum$yr.idx          = match(yr.cut,yr.level) + 0 * datum$yearn
initial               = datum$status == "INITIAL"
crashed               = datum$status == "CRASHED"
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
datum$yr.idx[crashed] = n.level + 6
yr.cscheme            = c("grey89",iatlas(n=n.level),"royalblue4","steelblue3","purple3"
                         ,"mediumpurple1","deepskyblue","hotpink")
ybottom               = rep(0,times=n.level+7)
ytop                  = rep(1,times=n.level+7)
xleft                 = seq(from=-1,to=n.level+5)
xright                = seq(from= 0,to=n.level+6)
xat                   = seq(from=-1,to=n.level+5)+0.5
xbrks                 = seq(from=-1,to=n.level+6)+0.5
xlabel                = c("Initial",yr.brks[-1],"Finish","StState","Extinct"
                         ,"Stopped","MetMiss","Crashed")
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
cat(" Plot the current status...","\n")
lo.box = pretty.box(n.iphen*n.iata)
for (st in 1:n.stext){
   #------ Find the soil texture key and description. -------------------------------------#
   this.st.key  = key.stext [st]
   letitre      = paste("Polygon status   --   Soil type: ",desc.stext[st]
                       ,sep="")
   lex          = paste("Rainfall change [scale parameter]")
   ley          = paste("Temperature change [K]")
   #---------------------------------------------------------------------------------------#

   for (o in 1:nout){
      #----- Open file. -------------------------------------------------------------------#
      fichier = paste(here,"/stt_",this.st.key,".",outform[o],sep="")
      if(outform[o] == "x11"){
         X11(width=size$width,height=size$height,pointsize=ptsz)
      }else if(outform[o] == "png"){
         png(filename=fichier,width=size$width*depth,height=size$height*depth
            ,pointsize=ptsz,res=depth)
      }else if(outform[o] == "eps"){
         postscript(file=fichier,width=size$width,height=size$height
                   ,pointsize=ptsz,paper=size$paper)
      }else if(outform[o] == "pdf"){
         pdf(file=fichier,onefile=FALSE,width=size$width,height=size$height
            ,pointsize=ptsz,paper=size$paper)
      }#end if
      #------------------------------------------------------------------------------------#


      #----- Save the margins to avoid losing the data. -----------------------------------#
      par.orig = par(no.readonly = TRUE)
      mar.orig = par.orig$mar
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
      text(x=xat,y=par("usr")[3]-0.6,labels=xlabel,srt=30,adj=1,xpd=TRUE,cex=1.15)
      title(main="Status",ylab="",xlab="")
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Now we loop over sites and pheonologies.                                       #
      #------------------------------------------------------------------------------------#
      k = 0
      for (ph in 1:n.iphen){
         for (pl in 1:n.iata){
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
      if (outform[o] == "x11"){
         locator(n=1)
         dev.off()
      }else{
         dev.off()
      }#end if
      bye = clean.tmp()
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#
}#end for
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#      Create parameter space maps for all other variables.                                #
#------------------------------------------------------------------------------------------#
cat(" Plot the current properties...","\n")
key.var  = c("lai","bsa","agb")
desc.var = c("Leaf area index [m2/m2]","Basal area [cm2/m2]"
            ,"Above-ground biomass [kgC/m2]")
n.var    = length(key.var)
lo.box = pretty.box(n.iphen*n.iata)
for (v in 1:n.var){
   cat("   - ",desc.var[v],"...","\n")
   #----- Collapse realisations using the median. -----------------------------------------#
   this.var       = apply( X      = datum[[key.var[v]]]
                         , MARGIN = c(1,2,4,5,6)
                         , FUN    = median
                         , na.rm  = TRUE
                         )#end apply
   rien           = ! is.finite(this.var)
   this.var[rien] = NA
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Break the data into bins.                                                        #
   #---------------------------------------------------------------------------------------#
   if (all(is.na(this.var))){
     var.brks    = c(-1,0,1)
     n.brks      = length(var.brks)
     var.cut     = cut(as.numeric(this.var),breaks=var.brks)
   }else{
     var.brks    = pretty(this.var,n=ncolours)
     n.brks      = length(var.brks)
     var.cut     = cut(this.var,breaks=var.brks)
   }#end if
   n.brks      = length(var.brks)
   var.cut     = cut(this.var,breaks=var.brks)
   var.lev     = levels(var.cut)
   var.idx     = match(var.cut,var.lev) + 0 * this.var
   var.cscheme = iatlas(n=n.brks-1)
   #---------------------------------------------------------------------------------------#


   #----- Make the edges. -----------------------------------------------------------------#
   xleft       = var.brks[-n.brks]
   xright      = var.brks[     -1]
   ybottom     = rep(0,times=n.brks)
   ytop        = rep(1,times=n.brks)
   xat         = var.brks
   #---------------------------------------------------------------------------------------#


   for (st in 1:n.stext){
      #------ Find the soil texture key and description. ----------------------------------#
      this.st.key  = key.stext [st]
      letitre      = paste(desc.var[v],"   --   Soil type: ",desc.stext[st],sep="")
      lex          = paste("Rainfall change [scale parameter]")
      ley          = paste("Temperature change [K]")
      #------------------------------------------------------------------------------------#

      for (o in 1:nout){
         #----- Open file. ----------------------------------------------------------------#
         fichier = paste(here,"/",key.var[v],"_",this.st.key,".",outform[o],sep="")
         if(outform[o] == "x11"){
            X11(width=size$width,height=size$height,pointsize=ptsz)
         }else if(outform[o] == "png"){
            png(filename=fichier,width=size$width*depth,height=size$height*depth
               ,pointsize=ptsz,res=depth)
         }else if(outform[o] == "eps"){
            postscript(file=fichier,width=size$width,height=size$height
                      ,pointsize=ptsz,paper=size$paper)
         }else if(outform[o] == "pdf"){
            pdf(file=fichier,onefile=FALSE,width=size$width,height=size$height
               ,pointsize=ptsz,paper=size$paper)
         }#end if
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
         rect(xleft=xleft,ybottom=ybottom,xright=xright,ytop=ytop,col=var.cscheme)
         box()
         axis(side=1,at=xat)
         title(main=desc.var[v],xlab="",ylab="")
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Now we loop over sites and pheonologies.                                    #
         #---------------------------------------------------------------------------------#
         k = 0
         for (ph in 1:n.iphen){
            for (pl in 1:n.iata){
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
               image(x=drain,y=dtemp,z=this.var[,,ph,pl,st],col=var.cscheme
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
         if (outform[o] == "x11"){
            locator(n=1)
            dev.off()
         }else{
            dev.off()
         }#end if
         bye = clean.tmp()
         #---------------------------------------------------------------------------------#
      }#end for
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#
}#end for
#------------------------------------------------------------------------------------------#
