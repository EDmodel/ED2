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
here    = getwd()                               #   Current directory
srcdir  = "/prj/prjidfca/marcosl/Util/Rsc"      #   Script directory
ibackground    = 0                              # Make figures compatible to background
                                                # 0 -- white
                                                # 1 -- black
                                                # 2 -- dark grey
#----- Output directory -------------------------------------------------------------------#
outroot = file.path(here,paste("longterm_comp_ibg",sprintf("%2.2i",ibackground),sep=""))
#------------------------------------------------------------------------------------------#



#----- Info on hourly data. ---------------------------------------------------------------#
reload.hour  = TRUE
reload.range = TRUE
rdata.path   = file.path(here,"RData_longterm")
rdata.suffix = "equilibrium_ed22.RData"
#------------------------------------------------------------------------------------------#



#----- Default phenology and fire model. --------------------------------------------------#
default.iphen = "phen+02"
default.ifire = "fire03"
emean.yeara   = 1972
emean.yearz   = 2011
#------------------------------------------------------------------------------------------#



#----- Number of "bootstrap" realisations. ------------------------------------------------#
nboot = 1000
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Site settings:                                                                       #
#                                                                                          #
# eort      -- first letter ("e" or "t")                                                   #
#                                                                                          #
# sites     -- a list of the sites that could be processed (see use.sites). Each list      #
#              element is also a list with the following information:                      #
#              iata      -- 3-letter code to refer to this site                            #
#              desc      -- Full name of the site (for titles).                            #
#              pch       -- site symbols                                                   #
#              col       -- background colour                                              #
#              fg        -- foreground colour                                              #
#              drya      -- beginning of the dry season (MM/DD) -- NA skips dry season     #
#              dryz      -- end of the dry season       (MM/DD) -- NA skips dry season     #
#                           dryz can be less than drya, in which case two rectangles are   #
#                           plotted.                                                       #
#                                                                                          #
# use.sites -- controls which sites are actually run.  It does different things depending  #
#              on the type.                                                                #
#              Logical   -- If TRUE, it uses all sites; FALSE doesn't make sense...        #
#              Character -- vector with the 3-letter code ("iata") of selected sites       #
#              Numeric   -- Site index.                                                    #
#                                                                                          #
#------------------------------------------------------------------------------------------#
eort  = "t"
n     = 0
sites = list()
n          = n + 1
sites[[n]] = list( iata    = "gyf"
                 , desc    = "Paracou"
                 , pch     = 2
                 , col     = "#3B24B3"
                 , fg      = "#160959"
                 , dbh.min = 10.
                 , ord     = 1
                 , bsa     = 31.0
                 , agb     = 14.9
                 , lai     = 6.15
                 )#end list
n          = n + 1
sites[[n]] = list( iata    = "s67"
                 , desc    = "Santarem km 67"
                 , pch     =  5
                 , col     = "#A3CC52"
                 , fg      = "#4B6614"
                 , dbh.min = 10.
                 , ord     = 3
                 , bsa     = 27.6
                 , agb     = 14.5
                 , lai     = 5.3
                 )#end list
n          = n + 1
sites[[n]] = list( iata    = "pdg"
                 , desc    = "Pe-de-Gigante"
                 , pch     = 13
                 , col     = "#990F0F"
                 , fg      = "#4D0404"
                 , dbh.min = 3.
                 , ord     = 4
                 , bsa     = 12.5
                 , agb     = 4.8
                 , lai     = 2.19
                 )#end list
n          = n + 1
sites[[n]] = list( iata    = "rja"
                 , desc    = "Rebio Jaru"
                 , pch     =  1
                 , col     = "#306614"
                 , fg      = "#143305"
                 , dbh.min = 5.
                 , ord     = 2
                 , bsa     = 32.6
                 , agb     = 16.9
                 , lai     = 4.9
                 )#end list
n          = n + 1
sites[[n]] = list( iata    = "m34"
                 , desc    = "Manaus K34"
                 , pch     =  6
                 , col     = "#2996CC"
                 , fg      = "#0A4766"
                 , dbh.min = 10.
                 , ord     = 9
                 , bsa     = 27.3
                 , agb     = 13.2
                 , lai     = 5.7
                 )#end list
n          = n + 1
sites[[n]] = list( iata    = "pnz"
                 , desc    = "Petrolina"
                 , pch     =  4
                 , col     = "#B49ED2"
                 , fg      = "#7D6E93"
                 , dbh.min = 3.
                 , ord     = 6
                 , bsa     = 19.3
                 , agb     = NA
                 , lai     = NA
                 )#end list
n          = n + 1
sites[[n]] = list( iata    = "ban"
                 , desc    = "Bananal"
                 , pch     =  8
                 , col     = "#F5C858"
                 , fg      = "#AB8C3D"
                 , dbh.min = 5.
                 , ord     = 8
                 , bsa     = 25.9
                 , agb     = NA
                 , lai     = NA
                 )#end list
n          = n + 1
sites[[n]] = list( iata    = "bsb"
                 , desc    = "Brasilia"
                 , pch     =  9
                 , col     = "#E65C17"
                 , fg      = "#732A06"
                 , dbh.min = 5.
                 , ord     = 5
                 , bsa     = 7.7
                 , agb     = NA
                 , lai     = 0.86
                 )#end list
n          = n + 1
sites[[n]] = list( iata    = "nat"
                 , desc    = "Natal"
                 , pch     =  0
                 , col     = "#00F3FB"
                 , fg      = "#00AAAF"
                 , dbh.min = 3.
                 , ord     = 7
                 , bsa     = 19.4
                 , agb     = NA
                 , lai     = NA
                 )#end list
use.sites  = TRUE
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#    Simulation settings:                                                                  #
# name -- the suffix of the simulations (list all combinations.                            #
# desc -- description (for legends)                                                        #
# verbose -- long description (for titles)                                                 #
# colour  -- colour to represent this simulation                                           #
#------------------------------------------------------------------------------------------#
sim.struct  = list( name        = c("ble_age30_pft02"  ,"ble_age30_pft05"
                                   ,"sas_age01_pft02"  ,"sas_age01_pft05"
                                   ,"sas_age30_pft02"  ,"sas_age30_pft05"  )
                  , desc        = c("Size 01 + Age 01 + PFT 02"
                                   ,"Size 01 + Age 01 + PFT 05"
                                   ,"Size 80 + Age 01 + PFT 02"
                                   ,"Size 80 + Age 01 + PFT 05"
                                   ,"Size 80 + Age 30 + PFT 02"
                                   ,"Size 80 + Age 30 + PFT 05"
                                   )#end c
                  , verbose     = c("2 PFTs"             ,"5 PFTs"
                                   ,"Size + 2 PFTs"      ,"Size + 5 PFTs"
                                   ,"Size + Age + 2 PFTs","Size + Age + 5 PFTs"
                                   )#end c
                   , colour     = c("#3B24B3","#2996CC"
                                   ,"#990F0F","#E65C17"
                                   ,"#306614","#A3CC52"
                                   )#end c
                   , fgcol      = c("#160959","#0A4766"
                                   ,"#4D0404","#732A06"
                                   ,"#143305","#4B6614"
                                   )#end c
                   , age.interp = c(    FALSE,    FALSE
                                   ,    FALSE,    FALSE
                                   ,     TRUE,     TRUE
                                   )#end c
                   )#end list
#----- List the default simulation. -------------------------------------------------------#
sim.default = 6
#------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------#
#       Plot options.                                                                      #
#------------------------------------------------------------------------------------------#
outform           = c("pdf")  # Formats for output file.  Supported formats are:
                              #   - "X11" - for printing on screen
                              #   - "eps" - for postscript printing
                              #   - "png" - for PNG printing
                              #   - "pdf" - for PDF printing

byeold            = TRUE      # Remove old files of the given format?

depth             = 96        # PNG resolution, in pixels per inch
paper             = "letter"  # Paper size, to define the plot shape
wpaper            = "legal"   # Wide paper size, to define the plot shape
ptsz              = 22        # Font size.
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Switch controls to plot only the needed ones.                                       #
#------------------------------------------------------------------------------------------#
plot.spider            = c(FALSE,TRUE)[2]
col.dryseason          = "papayawhip"
col.ust.altern         = "firebrick4"
col.ust.default        = "deeppink"
slz.cscheme            = "visible"
#------------------------------------------------------------------------------------------#






#------------------------------------------------------------------------------------------#
#      Colours for model comparison.                                                       #
#------------------------------------------------------------------------------------------#
obs.col   = "#6D6D6D"
obs.fg    = "#131313"
ed22.col  = "#99FF02"
ed22.fg   = "#548901"
alt.col   = "#9D57F8"
alt.fg    = "#2B0071"
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



#----- Load some packages. ----------------------------------------------------------------#
source(file.path(srcdir,"load.everything.r"))
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Eddy flux comparisons.                                                               #
#------------------------------------------------------------------------------------------#
n = 0
compvar       = list()
n             = n + 1
compvar[[ n]] = list( vnam         = "agb"
                    , desc         = "Above-ground biomass"
                    , unit         = untab$kgcom2
                    , cscheme.mean = "clife"
                    , hue.low      = "purple"
                    , hue.high     = "green"
                    , szpftvar     = TRUE
                    , patchvar     = TRUE
                    , cohortvar    = TRUE
                    , scalevar     = "nplant"
                    , zlog         = TRUE
                    , spider       = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "bsa"
                    , desc         = "Basal area"
                    , unit         = untab$m2om2
                    , cscheme.mean = "clife"
                    , hue.low      = "purple"
                    , hue.high     = "green"
                    , szpftvar     = TRUE
                    , patchvar     = TRUE
                    , cohortvar    = TRUE
                    , scalevar     = "nplant"
                    , zlog         = TRUE
                    , spider       = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "lai"
                    , desc         = "Leaf area index"
                    , unit         = untab$m2lom2
                    , cscheme.mean = "clife"
                    , hue.low      = "purple"
                    , hue.high     = "green"
                    , szpftvar     = TRUE
                    , patchvar     = TRUE
                    , cohortvar    = TRUE
                    , scalevar     = NA
                    , zlog         = FALSE
                    , spider       = TRUE
                    )#end list
n             = n + 1
compvar[[ n]] = list( vnam         = "tai"
                    , desc         = "Tree area index"
                    , unit         = untab$m2lom2
                    , cscheme.mean = "clife"
                    , hue.low      = "purple"
                    , hue.high     = "green"
                    , szpftvar     = TRUE
                    , patchvar     = TRUE
                    , cohortvar    = TRUE
                    , scalevar     = NA
                    , zlog         = FALSE
                    , spider       = TRUE
                    )#end list
#------------------------------------------------------------------------------------------#


#----- Set how many formats we must output. -----------------------------------------------#
outform = tolower(outform)
nout    = length(outform)
#------------------------------------------------------------------------------------------#



#----- Define plot window size ------------------------------------------------------------#
zsize   = plotsize( proje     = FALSE
                  , stdheight = 8.5 * 7/6
                  , stdwidth  = 11
                  )#end plotsize
#------------------------------------------------------------------------------------------#



#----- Create directory with RData. -------------------------------------------------------#
if (! file.exists(rdata.path)) dir.create(rdata.path)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Create all output directories, separated by format.                                 #
#------------------------------------------------------------------------------------------#
if (! file.exists(outroot)) dir.create(outroot)
out = list()
for (o in sequence(nout)){
   is.figure   = ! outform[o] %in% c("quartz","x11")
   this.form   = outform[o]


   #---- Main path for this output format. ------------------------------------------------#
   o.form = list()
   o.form$main = file.path(outroot,this.form)
   if (is.figure && ! file.exists(o.form$main)) dir.create(o.form$main)
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Create paths for the spider web plots.                                            #
   #---------------------------------------------------------------------------------------#
   o.main   = file.path(o.form$main,"spider")
   o.spider = list( main      = o.main
                  , variables = file.path(o.main,"variables")
                  )#end list
   if (is.figure){
      if (! file.exists(o.spider$main     )) dir.create(o.spider$main     )
      if (! file.exists(o.spider$variables)) dir.create(o.spider$variables)
   }#end if
   o.form$spider = o.spider
   #---------------------------------------------------------------------------------------#


   #----- Save the full list to the main path list. ---------------------------------------#
   out[[this.form]] = o.form
   #---------------------------------------------------------------------------------------#
}#end for (o in 1:nout)
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Combine all structures into a consistent list.                                       #
#------------------------------------------------------------------------------------------#
n.sim       = length(sim.struct$name)
#----- Simulation keys. -------------------------------------------------------------------#
simul.key   = sim.struct$name
#----- Description. -----------------------------------------------------------------------#
simleg.key  = sim.struct$desc
#---- Create the colours and line type for legend. ----------------------------------------#
simcol.key     = sim.struct$colour
simfgc.key     = sim.struct$fgcol
simlty.key     = rep("solid",times=n.sim)
simcex.key     = rep(2.0    ,times=n.sim)
simain.key     = sim.struct$age.interp
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Dump the information to a list.                                                      #
#------------------------------------------------------------------------------------------#
simul       = data.frame( name             = simul.key
                        , desc             = simleg.key
                        , colour           = simcol.key
                        , fgcol            = simfgc.key
                        , lty              = simlty.key
                        , cex              = simcex.key
                        , age.interp       = simain.key
                        , stringsAsFactors = FALSE
                        )#end data.frame
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#      Convert sites to data frame, and keep only those to be run.                         #
#------------------------------------------------------------------------------------------#
sites        = list.2.data.frame(sites)
if (is.logical(use.sites)){
   use.sites = which(rep(use.sites,times=nrow(sites))[sequence(nrow(sites))])
}else if (is.character(use.sites)){
   use.sites = match(use.sites,sites$iata)
   use.sites = use.sites[! is.na(use.sites)]
}#end if
if (length(use.sites) == 0){
   stop ("You must specify at least one valid site")
}else{
   sites = sites[use.sites,]
}#end if
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#      List the keys for all dimensions.                                                   #
#------------------------------------------------------------------------------------------#
compvar.key  = list.2.data.frame(compvar)$vnam
compvar.pch  = list.2.data.frame(compvar)$pch
#------------------------------------------------------------------------------------------#




#----- Set the various dimensions associated with variables, simulations, and sites. ------#
nsites     = length(sites$iata )
nsimul     = length(simul.key  )
ncompvar   = length(compvar.key)
#------------------------------------------------------------------------------------------#



#----- Initialise summary array. ----------------------------------------------------------#
summ = array( data     = NA
            , dim      = c(ncompvar,nsites,nsimul,3)
            , dimnames = list( compvar.key
                             , sites$iata
                             , simul.key
                             , c("expected","std.err","cw95")
                             )#end list
            )#end array
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#      Loop over all sites, variables, and simulations, to prepare the data for the        #
# model comparison.  We want to use only times with actual measurements, so we will        #
# discard model results from times with no observation so all derived quantities have      #
# the same number of defined points (so if measurements are biased towards daytime, the    #
# model will also be equally biased).                                                      #
#------------------------------------------------------------------------------------------#
for (p in sequence(nsites)){
   #----- Get the basic information. ------------------------------------------------------#
   iata          = sites$iata[p]
   im            = match(iata,poilist$iata)
   this          = list()
   this$short    = poilist$short   [im]
   this$longname = sites$desc      [ p]
   this$iata     = poilist$iata    [im]
   this$lon      = poilist$lon     [im]
   this$lat      = poilist$lat     [im]
   this$dbh.min  = sites$dbh.min   [ p]
   cat("   - Site :",this$longname,"...","\n")
   #---------------------------------------------------------------------------------------#





   #---------------------------------------------------------------------------------------#
   #     Get all the statistics and actual values for every simulation.                    #
   #---------------------------------------------------------------------------------------#
   for (s in sequence(nsimul)){
      cat("    # Simulation: ",simul$desc[s],"...","\n")

      #----- Load hourly averages. --------------------------------------------------------#
      ans.name = paste("e",iata,"_wmo_",simul$name[s],"_",default.iphen,"_",default.ifire
                      ,sep="")
      ans.path = file.path(here,ans.name)
      ans.file = file.path(ans.path,"rdata_month"
                          ,paste("40_years_",ans.name,".RData",sep=""))
      load(ans.file)
      #------------------------------------------------------------------------------------#





      #----- Create some variables to describe season and time of the day. ----------------#
      model = list()
      #----- Create time stamp for annual and monthly means. ------------------------------#
      esel     = numyears(datum$when) %wr% c(emean.yeara,emean.yearz)
      tomonth  = datum$when[esel]
      nemean   = length(tomonth)
      patch    = datum$patch
      cohort   = datum$cohort
      #------------------------------------------------------------------------------------#





      #------------------------------------------------------------------------------------#
      #      Initialise all variables.                                                     #
      #------------------------------------------------------------------------------------#
      template  = array( data     = NA
                       , dim      = c(nemean,2)
                       , dimnames = list(paste(tomonth),c("expected","std.err"))
                      )#end dimnames
      agb.table = template
      bsa.table = template
      lai.table = template
      tai.table = template
      #------------------------------------------------------------------------------------#





      #------------------------------------------------------------------------------------#
      #     Aggregate the data to patch level, and estimate the expected value and         #
      # standard error.                                                                    #
      #------------------------------------------------------------------------------------#
      for (e in sequence(nemean)){
         now     = tomonth[e]
         mm      = nummonths(now)
         yyyy    = numyears (now)
         stamp   = paste("y",sprintf("%4.4i",yyyy),"m",sprintf("%2.2i",mm),sep="")
         #if (mm == 1) cat("        > ",yyyy,"...","\n",sep="")

         #----- Get cohort-level variables. -----------------------------------------------#
         ipaco    = cohort$ipa   [[stamp]]
         areaco   = cohort$area  [[stamp]]
         agbco    = cohort$agb   [[stamp]] * cohort$nplant[[stamp]] / cohort$area[[stamp]]
         bsaco    = cohort$ba    [[stamp]] / cohort$area  [[stamp]]
         laico    = cohort$lai   [[stamp]]
         taico    = cohort$tai   [[stamp]]
         dbhco    = cohort$dbh   [[stamp]]
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Eliminate cohorts that are too small for AGB and basal area estimates.     #
         #---------------------------------------------------------------------------------#
         agbco = ifelse(dbhco %>=% this$dbh.min, agbco, 0.)
         bsaco = ifelse(dbhco %>=% this$dbh.min, bsaco, 0.)
         #---------------------------------------------------------------------------------#


         #----- Find the patch averages. --------------------------------------------------#
         areapa   = tapply( X = areaco, INDEX = ipaco, FUN = mean, na.rm = TRUE )
         agbpa    = tapply( X = agbco , INDEX = ipaco, FUN = sum , na.rm = TRUE )
         bsapa    = tapply( X = bsaco , INDEX = ipaco, FUN = sum , na.rm = TRUE )
         laipa    = tapply( X = laico , INDEX = ipaco, FUN = sum , na.rm = TRUE )
         taipa    = tapply( X = taico , INDEX = ipaco, FUN = sum , na.rm = TRUE )
         #---------------------------------------------------------------------------------#


         #----- Find the sampling size that makes the rarest patch appear twice. ----------#
         nps = ceiling(2*sum(areapa/min(areapa)))
         nnn = nps * nboot
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #       Sample data with replacement.                                             #
         #---------------------------------------------------------------------------------#
         ptsmp    = function(x,n,p) lit.sample(x=x,size=n,replace=TRUE,prob=p)
         agb.boot = matrix( data = ptsmp(x=agbpa,n=nnn,p=areapa),ncol=nps,nrow=nboot)
         bsa.boot = matrix( data = ptsmp(x=bsapa,n=nnn,p=areapa),ncol=nps,nrow=nboot)
         lai.boot = matrix( data = ptsmp(x=laipa,n=nnn,p=areapa),ncol=nps,nrow=nboot)
         tai.boot = matrix( data = ptsmp(x=taipa,n=nnn,p=areapa),ncol=nps,nrow=nboot)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #       Find the mean value for each realisation.                                 #
         #---------------------------------------------------------------------------------#
         agb.real = rowMeans(agb.boot)
         bsa.real = rowMeans(bsa.boot)
         lai.real = rowMeans(lai.boot)
         tai.real = rowMeans(tai.boot)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #       Find the first two moments.                                               #
         #---------------------------------------------------------------------------------#
         agb.table[e,] = c(mean(agb.real),sd(agb.real))
         bsa.table[e,] = c(mean(bsa.real),sd(bsa.real))
         lai.table[e,] = c(mean(lai.real),sd(lai.real))
         tai.table[e,] = c(mean(tai.real),sd(tai.real))
         #---------------------------------------------------------------------------------#
      }#end for (e in emean.loop)
      #------------------------------------------------------------------------------------#


      #----- Add summary for this simulation and site. ------------------------------------#
      summ["agb",p,s,] = c(mean(agb.table[,1]),mean2(agb.table[,2])*c(1,qt(.975,nboot-1)))
      summ["bsa",p,s,] = c(mean(bsa.table[,1]),mean2(bsa.table[,2])*c(1,qt(.975,nboot-1)))
      summ["lai",p,s,] = c(mean(lai.table[,1]),mean2(lai.table[,2])*c(1,qt(.975,nboot-1)))
      summ["tai",p,s,] = c(mean(tai.table[,1]),mean2(tai.table[,2])*c(1,qt(.975,nboot-1)))
      #------------------------------------------------------------------------------------#
   }#end for (s in sequence(nsimul))
   #---------------------------------------------------------------------------------------#

}#end for (p in sequence(loop.nsites))
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#      Save processed data to RData.                                                       #
#------------------------------------------------------------------------------------------#
rdata.summ = file.path(rdata.path,rdata.suffix)
cat(" + Saving processed data to ",basename(rdata.summ),"...","\n")
dummy = save(list=c("summ"), file=rdata.summ)
#------------------------------------------------------------------------------------------#







#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#     Plot the spider-web plots at equilibrium.                                            #
#------------------------------------------------------------------------------------------#
if (plot.spider){
   cat(" + Plotting spider web plots for the last cycle...","\n")


   #---------------------------------------------------------------------------------------#
   #     Loop over variables.                                                              #
   #---------------------------------------------------------------------------------------#
   for (v in sequence(ncompvar)){
      #----- Load variable settings. ------------------------------------------------------#
      this.compvar    = compvar[[v]]
      this.vnam       = this.compvar$vnam
      this.desc       = this.compvar$desc
      this.unit       = this.compvar$unit
      this.spider     = this.compvar$spider
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #     Check that this variable is a spider-web variable.                             #
      #------------------------------------------------------------------------------------#
      if (this.spider){
         cat("   - ",this.desc,"...","\n")


         #----- Site order. ---------------------------------------------------------------#
         p.ord = order(sites$ord)
         #---------------------------------------------------------------------------------#
         web.range = c(0,max(summ[this.vnam,,,"expected"],na.rm=TRUE))
         web.lim   = pretty(web.range,n=4)
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #     Make one panel for each simulation, and one plot per site.                  #
         #---------------------------------------------------------------------------------#


         #------ Set title. ---------------------------------------------------------------#
         letitre = desc.unit(desc=this.desc,unit=this.unit)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Loop over formats.                                                         #
         #---------------------------------------------------------------------------------#
         for (o in sequence(nout)){
            #----- Make the file name. ----------------------------------------------------#
            out.now = out[[outform[o]]]$spider$variables
            fichier = file.path(out.now,paste("spider-",this.vnam,".",outform[o],sep=""))
            if (outform[o] == "x11"){
              X11(width=zsize$width,height=zsize$height,pointsize=col.use)
            }else if(outform[o] == "png"){
               png(filename=fichier,width=zsize$width*depth,height=zsize$height*depth
                  ,pointsize=ptsz,res=depth)
            }else if(outform[o] == "eps"){
               postscript(file=fichier,width=zsize$width,height=zsize$height
                         ,pointsize=ptsz,paper=zsize$paper)
            }else if(outform[o] == "pdf"){
               pdf(file=fichier,onefile=FALSE,width=zsize$width,height=zsize$height
                  ,pointsize=ptsz,paper=zsize$paper)
            }#end if
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Split the window into 2, and add simulation legend at the bottom.        #
            #------------------------------------------------------------------------------#
            par(par.user)
            par.orig = par(no.readonly = TRUE)
            mar.orig = par.orig$mar
            par(oma = c(0.2,0,3.0,0))
            layout(mat = rbind(2,1),height = c(6.0,1.0))
            #------------------------------------------------------------------------------#




            #----- Legend: the simulations. -----------------------------------------------#
            par(mar=c(0.1,0.1,0.1,0.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
            legend ( x       = "bottom"
                   , inset   = 0.0
                   , legend  = simul$desc
                   , fill    = simul$col
                   , border  = simul$col
                   , ncol    = 3
                   , title   = expression(bold("Structure"))
                   , pt.cex  = simul$cex
                   , cex     = 0.75 * cex.ptsz
                   )#end legend
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #     Plot the spider web.                                                     #
            #------------------------------------------------------------------------------#
            radial.flex( lengths          = t(summ[this.vnam,p.ord,,"expected"])
                       , labels           = toupper(sites$iata[p.ord])
                       , radlab           = FALSE
                       , start            = 90
                       , clockwise        = TRUE
                       , rp.type          = "p"
                       , main             = ""
                       , line.col         = simul$col
                       , lty              = simul$lty
                       , lwd              = 3.0
                       , show.grid        = TRUE
                       , show.grid.labels = 1
                       , show.radial.grid = TRUE
                       , grid.col         = grid.colour
                       , radial.lim       = web.lim
                       , poly.col         = NA
                       , mar              = c(2,1,2,1)+0.1
                       , cex.lab          = 0.5
                       )#end radial.plot
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #    Add the observed value in case there is a reference for at least one      #
            # site.                                                                        #
            #------------------------------------------------------------------------------#
            if (this.vnam %in% names(sites)){

               #---------------------------------------------------------------------------#
               #     Plot the spider web.                                                  #
               #---------------------------------------------------------------------------#
               radial.flex( lengths          = matrix(sites[[this.vnam]][p.ord],nrow=1)
                          , labels           = toupper(sites$iata[p.ord])
                          , radlab           = FALSE
                          , start            = 90
                          , clockwise        = TRUE
                          , rp.type          = "s"
                          , main             = ""
                          , point.symbols    = 16
                          , cex              = 0.5
                          , line.col         = grey.mg
                          , radial.lim       = web.lim
                          , lwd              = 3.0
                          , add              = TRUE
                          )#end radial.plot
               radial.flex( lengths          = matrix(sites[[this.vnam]][p.ord],nrow=1)
                          , labels           = toupper(sites$iata[p.ord])
                          , radlab           = FALSE
                          , start            = 90
                          , clockwise        = TRUE
                          , rp.type          = "s"
                          , main             = ""
                          , point.symbols    = 1
                          , line.col         = grey.mg
                          , radial.lim       = web.lim
                          , lwd              = 3.0
                          , add              = TRUE
                          )#end radial.plot
               #---------------------------------------------------------------------------#
            }#end if (this.vnam %in% names(obser))
            #------------------------------------------------------------------------------#

            #------------------------------------------------------------------------------#
            #     Plot the global title.                                                   #
            #------------------------------------------------------------------------------#
            par(las=0)
            mtext(text=letitre,side=3,outer=TRUE,cex=1.1,font=2)
            #------------------------------------------------------------------------------#




            #----- Close the device. ------------------------------------------------------#
            if (outform[o] == "x11"){
               locator(n=1)
               dev.off()
            }else{
               dev.off()
            }#end if
            dummy = clean.tmp()
            #------------------------------------------------------------------------------#
         }#end for (o in sequence(nout))
         #---------------------------------------------------------------------------------#
      }#end if (this.spider)
      #------------------------------------------------------------------------------------#
   }#end for (v in sequence(ncompvar))
   #---------------------------------------------------------------------------------------#
}#end if (plot.spider)
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
