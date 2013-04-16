#==========================================================================================#
#==========================================================================================#
#     List of possible time series plots. Pay attention to this list because all the other #
# plots will depend on this list.                                                          #
#                                                                                          #
#     Variable is.sum tells whether to sum the total (TRUE), or to average (FALSE).        #
#------------------------------------------------------------------------------------------#
   n            = 0
   scen.ts      = list()
   n            = n + 1
   scen.ts[[n]] = list( vname    = "agb"
                      , desc     = "Above ground biomass"
                      , unit     = untab$kgcom2
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end lis
   n            = n + 1
   scen.ts[[n]] = list( vname    = "bgb"
                      , desc     = "Below ground biomass"
                      , unit     = untab$kgcom2
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "biomass"
                      , desc     = "Total biomass"
                      , unit     = untab$kgcom2
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end lis
   n            = n + 1
   scen.ts[[n]] = list( vname    = "lai"
                      , desc     = "Leaf area index"
                      , unit     = untab$m2lom2
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "nplant"
                      , desc     = "Plant density"
                      , unit     = untab$plom2
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end lis
   n            = n + 1
   scen.ts[[n]] = list( vname    = "gpp"
                      , desc     = "Gross primary productivity"
                      , unit     = untab$kgcom2oyr
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "npp"
                      , desc     = "Net primary productivity"
                      , unit     = untab$kgcom2oyr
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "mco"
                      , desc     = "Maintenance costs"
                      , unit     = untab$kgcom2oyr
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "cba"
                      , desc     = "Carbon balance"
                      , unit     = untab$kgcom2oyr
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "ldrop"
                      , desc     = "Leaf drop"
                      , unit     = untab$kgcom2oyr
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "bstorage"
                      , desc     = "Storage biomass"
                      , unit     = untab$kgcom2
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = FALSE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "bseeds"
                      , desc     = "Seed biomass"
                      , unit     = untab$kgcom2
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = FALSE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "ba"
                      , desc     = "Basal area"
                      , unit     = untab$cm2om2
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "fs.open"
                      , desc     = "Soil moisture stress factor"
                      , unit     = untab$empty
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "mort"
                      , desc     = "Mortality rate"
                      , unit     = untab$pcpopoyr
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = TRUE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = TRUE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "ncbmort"
                      , desc     = "Mortality rate - Neg. C balance"
                      , unit     = untab$pcpopoyr
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = TRUE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = TRUE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "dimort"
                      , desc     = "Mortality rate - Density-independent"
                      , unit     = untab$pcpopoyr
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = TRUE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = TRUE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "recr"
                      , desc     = "Recruitment rate"
                      , unit     = untab$pcpopoyr
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = TRUE
                      , plog     = FALSE
                      , plog.dbh = TRUE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "growth"
                      , desc     = "Growth rate (DBH)"
                      , unit     = untab$pcdbhoyr
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = TRUE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "agb.mort"
                      , desc     = "Mortality rate"
                      , unit     = untab$pcagboyr
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = TRUE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = TRUE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "agb.ncbmort"
                      , desc     = "Mortality rate - Neg. C balance"
                      , unit     = untab$pcagboyr
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = TRUE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = TRUE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "agb.dimort"
                      , desc     = "Mortality rate - Density-independent"
                      , unit     = untab$pcagboyr
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = TRUE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = TRUE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "agb.recr"
                      , desc     = "Recruitment rate"
                      , unit     = untab$pcagboyr
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = TRUE
                      , plog     = FALSE
                      , plog.dbh = TRUE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "agb.growth"
                      , desc     = "Growth rate (AGB)"
                      , unit     = untab$pcagboyr
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = TRUE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "cbamax"
                      , desc     = "Maximum carbon balance"
                      , unit     = untab$kgcom2oyr
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "cbarel"
                      , desc     = "Relative carbon balance"
                      , unit     = untab$empty
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "demand"
                      , desc     = "Water demand"
                      , unit     = untab$kgwom2oday
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "supply"
                      , desc     = "Water supply"
                      , unit     = untab$kgwom2oday
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "rain"
                      , desc     = "Precipitation"
                      , unit     = untab$mmoyr
                      , f.aggr   = "sum"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "runoff"
                      , desc     = "Total runoff"
                      , unit     = untab$mmoyr
                      , f.aggr   = "sum"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "intercepted"
                      , desc     = "Canopy interception"
                      , unit     = untab$mmoyr
                      , f.aggr   = "sum"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "last.1yr.rain"
                      , desc     = "Previous 12 months rainfall"
                      , unit     = untab$mmoyr
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "last.2yr.rain"
                      , desc     = "Previous 24 months rainfall"
                      , unit     = untab$mmoyr
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "last.3yr.rain"
                      , desc     = "Previous 36 months rainfall"
                      , unit     = untab$mmoyr
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "water.deficit"
                      , desc     = "Water deficit (ED-2.2)"
                      , unit     = untab$mmoyr
                      , f.aggr   = "max"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "malhi.deficit"
                      , desc     = "Water deficit (Malhi 2009)"
                      , unit     = untab$mmoyr
                      , f.aggr   = "max"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "atm.temp"
                      , desc     = "Air temperature"
                      , unit     = untab$degC
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "leaf.temp"
                      , desc     = "Leaf temperature"
                      , unit     = untab$degC
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "rshort"
                      , desc     = "Incoming shortwave radiation"
                      , unit     = untab$wom2
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "rlong"
                      , desc     = "Incoming longwave radiation"
                      , unit     = untab$wom2
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "atm.vpd"
                      , desc     = "Air vapour pressure deficit"
                      , unit     = untab$pa
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "leaf.vpd"
                      , desc     = "Leaf vapour pressure deficit"
                      , unit     = untab$pa
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "paw"
                      , desc     = "Potential Available Water"
                      , unit     = untab$pcsat
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "smpot"
                      , desc     = "Soil Matric Potential"
                      , unit     = untab$mpa
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = TRUE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "nep"
                      , desc     = "Net Ecosystem Productivity"
                      , unit     = untab$kgcom2oyr
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "reco"
                      , desc     = "Ecosystem Respiration"
                      , unit     = untab$kgcom2oyr
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "fast.soil.c"
                      , desc     = "Fast soil carbon"
                      , unit     = untab$kgcom2
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "struct.soil.c"
                      , desc     = "Structural soil carbon"
                      , unit     = untab$kgcom2
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "slow.soil.c"
                      , desc     = "Slow soil carbon"
                      , unit     = untab$kgcom2
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "het.resp"
                      , desc     = "Heterotrophic Respiration"
                      , unit     = untab$kgcom2oyr
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "plant.resp"
                      , desc     = "Plant Respiration"
                      , unit     = untab$kgcom2oyr
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "hflxlc"
                      , desc     = "Leaf sensible heat"
                      , unit     = untab$wom2
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "wflxlc"
                      , desc     = "Leaf Evaporation"
                      , unit     = untab$kgwom2oday
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "transp"
                      , desc     = "Leaf Transpiration"
                      , unit     = untab$kgwom2oday
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "i.gpp"
                      , desc     = "Mean Gross Primary Production"
                      , unit     = untab$kgcoployr
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "i.npp"
                      , desc     = "Mean Net Primary Production"
                      , unit     = untab$kgcoployr
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "i.plant.resp"
                      , desc     = "Mean Plant Respiration"
                      , unit     = untab$kgcoployr
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "i.cba"
                      , desc     = "Mean Carbon balance"
                      , unit     = untab$kgcoployr
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "i.cbamax"
                      , desc     = "Mean Maximum C balance"
                      , unit     = untab$kgcoployr
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "i.hflxlc"
                      , desc     = "Mean Leaf sensible heat flux"
                      , unit     = untab$wopl
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "i.wflxlc"
                      , desc     = "Mean Leaf evaporation"
                      , unit     = untab$kgwoploday
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "i.transp"
                      , desc     = "Mean Transpiration"
                      , unit     = untab$kgwoploday
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "hflxgc"
                      , desc     = "Sensible heat - Gnd->CAS"
                      , unit     = untab$wom2
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "hflxca"
                      , desc     = "Sensible heat - CAS->ATM"
                      , unit     = untab$wom2
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "wflxgc"
                      , desc     = "Water flux - Gnd->CAS"
                      , unit     = untab$wom2
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "wflxca"
                      , desc     = "Water flux - CAS->ATM"
                      , unit     = untab$kgwom2oday
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "wue"
                      , desc     = "Water use efficiency"
                      , unit     = untab$gcokgw
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "leaf.gbw"
                      , desc     = "Leaf Bnd. Lyr. Conductance"
                      , unit     = untab$kgwom2loday
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "leaf.gsw"
                      , desc     = "Leaf stomatal Conductance"
                      , unit     = untab$kgwom2loday
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "f.gpp"
                      , desc     = "Gross Primary Productivity"
                      , unit     = untab$pcbiooyr
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "f.plant.resp"
                      , desc     = "Plant respiration"
                      , unit     = untab$pcbiooyr
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "f.npp"
                      , desc     = "Net Primary Productivity"
                      , unit     = untab$pcbiooyr
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "f.cba"
                      , desc     = "Carbon balance"
                      , unit     = untab$pcbiooyr
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "f.bstorage"
                      , desc     = "Relative storage biomass"
                      , unit     = untab$gcokgcbio
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "f.bleaf"
                      , desc     = "Relative leaf biomass"
                      , unit     = untab$gcokgcbio
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "f.broot"
                      , desc     = "Relative root biomass"
                      , unit     = untab$gcokgcbio
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "f.bseeds"
                      , desc     = "Relative seed biomass"
                      , unit     = untab$gcokgcbio
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "leaf.par"
                      , desc     = "Absorbed PAR - Leaf"
                      , unit     = untab$umolom2los
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1.e6
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "leaf.rshort"
                      , desc     = "Absorbed SW - Leaf"
                      , unit     = untab$wom2l
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "leaf.rlong"
                      , desc     = "Absorbed LW - Leaf"
                      , unit     = untab$wom2l
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Replace the list by a data frame.                                                 #
   #---------------------------------------------------------------------------------------#
   scen.ts = data.frame( apply( X = sapply(X=scen.ts,FUN=c), MARGIN = 1, FUN = unlist )
                       , stringsAsFactors = FALSE
                       )#end data.frame
   for (nl in c("pftvar","dbhvar","mort","recr","plog","plog.dbh","plt")){
      scen.ts[[nl]] = as.logical(scen.ts[[nl]])
   }#end for
   for (nl in c("add","mult")){
      scen.ts[[nl]] = as.numeric(scen.ts[[nl]])
   }#end for
   #---------------------------------------------------------------------------------------#
#==========================================================================================#
#==========================================================================================#




#==========================================================================================#
#==========================================================================================#
#     This list tells which variables to do the plot of size and PFT-dependent box plot.   #
# Units, description, and log scale will be copied from scen.ts.                           #
#                                                                                          #
# IMPORTANT:  All variables here MUST come from one of the variables defined in scen.ts,   #
#             and they must be PFT-dependent or PFT- and DBH- dependent (in which case two #
#             plots will be generated).                                                    #
#------------------------------------------------------------------------------------------#
   #----- All that we need here is the variable name. -------------------------------------#
   scen.szpft = data.frame( vname = c(          "agb",          "lai",           "ba"
                                     ,         "recr",         "mort",      "ncbmort"
                                     ,       "dimort",       "growth",     "agb.recr"
                                     ,     "agb.mort",  "agb.ncbmort",   "agb.dimort"
                                     ,   "agb.growth",          "gpp",          "npp"
                                     ,   "plant.resp",          "cba",          "mco"
                                     ,     "bstorage",       "bseeds",      "fs.open"
                                     ,       "supply",       "demand",       "hflxlc"
                                     ,       "wflxlc",       "transp",        "i.gpp"
                                     , "i.plant.resp",        "i.npp",        "i.cba"
                                     ,     "i.transp",       "cbamax",     "i.hflxlc"
                                     ,     "i.wflxlc",          "wue",     "i.cbamax"
                                     ,     "leaf.gsw",     "leaf.gbw",        "f.gpp"
                                     ,        "f.npp",        "f.cba",   "f.bstorage"
                                     ,      "f.bleaf",     "f.bseeds",         "rain"
                                     ,"last.1yr.rain","last.2yr.rain","last.3yr.rain"
                                     ,       "runoff",  "intercepted","water.deficit"
                                     ,"malhi.deficit",      "atm.vpd",     "leaf.vpd"
                                     ,    "leaf.temp",     "atm.temp",     "leaf.par"
                                     ,       "nplant"
                                     )#end c
                          , stringsAsFactors = FALSE
                          )#end vname
   #---------------------------------------------------------------------------------------#






   #---------------------------------------------------------------------------------------#
   #     Fill in the box plot list with information brought from scen.ts.                  #
   #---------------------------------------------------------------------------------------#
   #----- Find the variable names to be added. --------------------------------------------#
   which.names = names(scen.ts)
   keep        = ! ( which.names %in% names(scen.szpft))
   which.names = which.names[keep]
   #----- Match the lists based on vname. -------------------------------------------------#
   sz.idx = match(scen.szpft$vname,scen.ts$vname)
   #----- Look for variables that weren't defined in scen.ts. -----------------------------#
   sz.miss = is.na(sz.idx)
   if (any(sz.miss)){
      cat(" - The following variables in scen.szpft are missing from scen.ts:","\n")
      cat(paste("   * ",scen.szpft$vname[sz.miss],sep=""),sep="\n")
      stop(" - All variables defined in scen.szpft must be defined in scen.ts!!!")
   }#end if(any(sz.miss))
   #---------------------------------------------------------------------------------------#



   #----- Append the information to the data frame. ---------------------------------------#
   for (wn in which.names) scen.szpft[[wn]] = scen.ts[[wn]][sz.idx]
   #---------------------------------------------------------------------------------------#
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This list tells which variables to do the bar plot of size and PFT-dependent, long-  #
# term means for each season.                                                              #
# Units, description, and log scale will be copied from scen.ts.                           #
#                                                                                          #
# IMPORTANT:  All variables here MUST come from one of the variables defined in scen.ts,   #
#             and they must be PFT- and DBH- dependent.  Avoid variables that can be       #
#             negative (e.g. npp), or that cannot be stacked (e.g. fs.open)                #
#------------------------------------------------------------------------------------------#
   #----- All that we need here is the variable name. -------------------------------------#
   scen.barplot = data.frame( vname = c(         "agb",         "lai",          "ba"
                                       ,         "gpp",       "ldrop",         "mco"
                                       ,      "supply",      "demand",      "transp"
                                       ,  "plant.resp",       "i.gpp","i.plant.resp"
                                       ,    "i.transp",      "bseeds",    "bstorage"
                                       )#end c
                            , stringsAsFactors = FALSE
                            )#end vname
   #---------------------------------------------------------------------------------------#






   #---------------------------------------------------------------------------------------#
   #     Fill in the box plot list with information brought from scen.ts.                  #
   #---------------------------------------------------------------------------------------#
   #----- Find the variable names to be added. --------------------------------------------#
   which.names = names(scen.ts)
   keep        = ! ( which.names %in% names(scen.barplot))
   which.names = which.names[keep]
   #----- Match the lists based on vname. -------------------------------------------------#
   bar.idx = match(scen.barplot$vname,scen.ts$vname)
   #----- Look for variables that weren't defined in scen.ts. -----------------------------#
   bar.miss = is.na(bar.idx)
   if (any(bar.miss)){
      cat(" - The following variables in scen.barplot are missing from scen.ts:","\n")
      cat(paste("   * ",scen.barplot$vname[bar.miss],sep=""),sep="\n")
      stop(" - All variables defined in scen.barplot must be defined in scen.ts!!!")
   }#end if(any(bar.miss))
   #---------------------------------------------------------------------------------------#



   #----- Append the information to the data frame. ---------------------------------------#
   for (wn in which.names) scen.barplot[[wn]] = scen.ts[[wn]][bar.idx]
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    All variables here must be PFT- and DBH- dependent, so we make sure that they are  #
   # indeed.                                                                               #
   #---------------------------------------------------------------------------------------#
   if (! all(scen.barplot$pftvar & scen.barplot$dbhvar)){
      not.pftdbh = ! ( scen.barplot$pftvar & scen.barplot$dbhvar)
      cat (" - The following variables in scen.barplot are not PFT+DBH variables:","\n")
      cat (paste("   * ",scen.barplot$vname[not.pftdbh],sep=""),sep="\n")
      stop(" - All variables defined in scen.barplot must be PFT+DBH variables!!!")
   }#end if
   #---------------------------------------------------------------------------------------#
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     XYZ plots, to explore the parameter space.  Units, description, and log scale will   #
# be copied from scen.ts.                                                                  #
#                                                                                          #
#   IMPORTANT: All variables here MUST come from one of the variables defined in scen.ts!! #
#------------------------------------------------------------------------------------------#
   #----- All that we need here is the variable name, legend position or colour scheme. ---#
   scen.xyz = list()
   scen.xyz$xvar      = list( list( vname = "rshort"       , leg        = "right"  )
                            , list( vname = "leaf.temp"    , leg        = "right"  )
                            , list( vname = "leaf.par"     , leg        = "right"  )
                            , list( vname = "leaf.vpd"     , leg        = "right"  )
                            , list( vname = "leaf.gsw"     , leg        = "right"  )
                            , list( vname = "rain"         , leg        = "right"  )
                            , list( vname = "intercepted"  , leg        = "right"  )
                            , list( vname = "runoff"       , leg        = "right"  )
                            , list( vname = "last.1yr.rain", leg        = "right"  )
                            , list( vname = "last.2yr.rain", leg        = "right"  )
                            , list( vname = "last.3yr.rain", leg        = "right"  )
                            , list( vname = "smpot"        , leg        = "right"  )
                            , list( vname = "water.deficit", leg        = "right"  )
                            , list( vname = "bstorage"     , leg        = "right"  )
                            , list( vname = "f.bstorage"   , leg        = "right"  )
                            , list( vname = "bseeds"       , leg        = "right"  )
                            , list( vname = "f.bseeds"     , leg        = "right"  )
                            , list( vname = "leaf.gsw"     , leg        = "right"  )
                            , list( vname = "cba"          , leg        = "top"    )
                            , list( vname = "f.cba"        , leg        = "top"    )
                            )#end list
   scen.xyz$yvar      = list( list( vname = "recr"         , leg        = "top"    )
                            , list( vname = "mort"         , leg        = "top"    )
                            , list( vname = "ncbmort"      , leg        = "top"    )
                            , list( vname = "growth"       , leg        = "top"    )
                            , list( vname = "agb.recr"     , leg        = "top"    )
                            , list( vname = "agb.mort"     , leg        = "top"    )
                            , list( vname = "agb.ncbmort"  , leg        = "top"    )
                            , list( vname = "agb.growth"   , leg        = "top"    )
                            , list( vname = "leaf.gsw"     , leg        = "top"    )
                            )#end list
   scen.xyz$zvar      = list( list( vname = "lai"          , col.scheme = "clife" )
                            , list( vname = "ba"           , col.scheme = "clife" )
                            )#end list
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Replace the list by a data frame.                                                 #
   #---------------------------------------------------------------------------------------#
   scen.xyz$xvar = data.frame( apply(X=sapply(X=scen.xyz$xvar,FUN=c),MARGIN=1,FUN=unlist)
                             , stringsAsFactors = FALSE
                             )#end data.frame
   scen.xyz$yvar = data.frame( apply(X=sapply(X=scen.xyz$yvar,FUN=c),MARGIN=1,FUN=unlist)
                             , stringsAsFactors = FALSE
                             )#end data.frame
   scen.xyz$zvar = data.frame( apply(X=sapply(X=scen.xyz$zvar,FUN=c),MARGIN=1,FUN=unlist)
                             , stringsAsFactors = FALSE
                             )#end data.frame
   #---------------------------------------------------------------------------------------#






   #---------------------------------------------------------------------------------------#
   #     Fill in the XYZ list with information brought from scen.ts.                       #
   #---------------------------------------------------------------------------------------#
   #----- Get all the names that shall be added. ------------------------------------------#
   which.names = names(scen.ts) 
   keep        = ( ! which.names %in% union( union( names(scen.xyz$xvar)
                                                  , names(scen.xyz$yvar)
                                                  )#end union
                                           , names(scen.xyz$zvar)
                                           )#end union
                 )#end keep
   which.names = which.names[keep]
   #----- Match the lists based on vname. -------------------------------------------------#
   x.idx = match(scen.xyz$xvar$vname,scen.ts$vname)
   y.idx = match(scen.xyz$yvar$vname,scen.ts$vname)
   z.idx = match(scen.xyz$zvar$vname,scen.ts$vname)
   #----- Look for variables that weren't defined in scen.ts. -----------------------------#
   x.sel = is.na(x.idx)
   y.sel = is.na(y.idx)
   z.sel = is.na(z.idx)
   if (any(x.sel) | any(y.sel) | any(z.sel)){
      if (any(x.sel)){
         cat(" - The following variables in scen.xyz$xvar are missing from scen.ts:","\n")
         cat(paste("   * ",scen.xyz$xvar$vname[x.sel],sep=""),sep="\n")
      }#end if(any(x.sel))
      if (any(y.sel)){
         cat(" - The following variables in scen.xyz$yvar are missing from scen.ts:","\n")
         cat(paste("   * ",scen.xyz$yvar$vname[y.sel],sep=""),sep="\n")
      }#end if(any(y.sel))
      if (any(z.sel)){
         cat(" - The following variables in scen.xyz$zvar are missing from scen.ts:","\n")
         cat(paste("   * ",scen.xyz$zvar$vname[z.sel],sep=""),sep="\n")
      }#end if(any(z.sel))
      stop(" - All variables defined in scen.xyz must be defined in scen.ts!!!")
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Append the information to the data frame. ---------------------------------------#
   for (wn in which.names){
      scen.xyz$xvar[[wn]] = scen.ts[[wn]][x.idx]
      scen.xyz$yvar[[wn]] = scen.ts[[wn]][y.idx]
      scen.xyz$zvar[[wn]] = scen.ts[[wn]][z.idx]
   }#end for
   #---------------------------------------------------------------------------------------#
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This list tells which variables to do the scenario comparison panel, with each panel #
# representing a collection of scenarios.  Only mean differences will be plotted.  Because #
# this works with differences, log variables will not be plotted as log variables.         #
#                                                                                          #
# IMPORTANT:  All variables here MUST come from one of the variables defined in scen.ts.   #
#------------------------------------------------------------------------------------------#
   #----- All that we need here is the variable name. -------------------------------------#
   scen.comp = list( list( vname =           "agb", low = "purple"   , high = "green"    )
                   , list( vname =           "lai", low = "purple"   , high = "green"    )
                   , list( vname =            "ba", low = "purple"   , high = "green"    )
                   , list( vname =        "nplant", low = "purple"   , high = "green"    )
                   , list( vname =           "gpp", low = "purple"   , high = "green"    )
                   , list( vname =           "npp", low = "purple"   , high = "green"    )
                   , list( vname =           "mco", low = "purple"   , high = "green"    )
                   , list( vname =           "cba", low = "purple"   , high = "green"    )
                   , list( vname =         "ldrop", low = "green"    , high = "purple"   )
                   , list( vname =      "bstorage", low = "purple"   , high = "green"    )
                   , list( vname =        "bseeds", low = "purple"   , high = "green"    )
                   , list( vname =       "fs.open", low = "orangered", high = "blue"     )
                   , list( vname =          "mort", low = "green"    , high = "purple"   )
                   , list( vname =          "recr", low = "purple"   , high = "green"    )
                   , list( vname =        "growth", low = "purple"   , high = "green"    )
                   , list( vname =       "ncbmort", low = "green"    , high = "purple"   )
                   , list( vname =        "dimort", low = "green"    , high = "purple"   )
                   , list( vname =      "agb.mort", low = "green"    , high = "purple"   )
                   , list( vname =      "agb.recr", low = "purple"   , high = "green"    )
                   , list( vname =    "agb.growth", low = "purple"   , high = "green"    )
                   , list( vname =   "agb.ncbmort", low = "green"    , high = "purple"   )
                   , list( vname =    "agb.dimort", low = "green"    , high = "purple"   )
                   , list( vname =        "cbarel", low = "purple"   , high = "green"    )
                   , list( vname =        "demand", low = "blue"     , high = "orangered")
                   , list( vname =        "supply", low = "orangered", high = "blue"     )
                   , list( vname =          "rain", low = "orangered", high = "blue"     )
                   , list( vname =        "runoff", low = "orangered", high = "blue"     )
                   , list( vname =   "intercepted", low = "orangered", high = "blue"     )
                   , list( vname = "water.deficit", low = "blue"     , high = "orangered")
                   , list( vname =      "atm.temp", low = "blue"     , high = "orangered")
                   , list( vname =     "leaf.temp", low = "blue"     , high = "orangered")
                   , list( vname =        "rshort", low = "grey"     , high = "blue"     )
                   , list( vname =         "rlong", low = "blue"     , high = "grey"     )
                   , list( vname =       "atm.vpd", low = "blue"     , high = "orangered")
                   , list( vname =      "leaf.vpd", low = "blue"     , high = "orangered")
                   , list( vname =           "paw", low = "orangered", high = "blue"     )
                   , list( vname =         "smpot", low = "blue"     , high = "orangered")
                   , list( vname =           "nep", low = "purple"   , high = "green"    )
                   , list( vname =          "reco", low = "blue"     , high = "orangered")
                   , list( vname =   "fast.soil.c", low = "blue"     , high = "orangered")
                   , list( vname = "struct.soil.c", low = "blue"     , high = "orangered")
                   , list( vname =   "slow.soil.c", low = "blue"     , high = "orangered")
                   , list( vname =      "het.resp", low = "blue"     , high = "orangered")
                   , list( vname =    "plant.resp", low = "blue"     , high = "orangered")
                   , list( vname =        "hflxlc", low = "blue"     , high = "orangered")
                   , list( vname =        "wflxlc", low = "orangered", high = "blue"     )
                   , list( vname =        "transp", low = "orangered", high = "blue"     )
                   , list( vname =        "hflxgc", low = "blue"     , high = "orangered")
                   , list( vname =        "hflxca", low = "blue"     , high = "orangered")
                   , list( vname =        "wflxgc", low = "orangered", high = "blue"     )
                   , list( vname =        "wflxca", low = "orangered", high = "blue"     )
                   , list( vname =      "i.hflxlc", low = "blue"     , high = "orangered")
                   , list( vname =      "i.wflxlc", low = "orangered", high = "blue"     )
                   , list( vname =      "i.transp", low = "orangered", high = "blue"     )
                   , list( vname =         "i.gpp", low = "purple"   , high = "green"    )
                   , list( vname =         "i.npp", low = "purple"   , high = "green"    )
                   , list( vname =  "i.plant.resp", low = "green"    , high = "purple"   )
                   , list( vname =         "i.cba", low = "purple"   , high = "green"    )
                   , list( vname =      "i.cbamax", low = "purple"   , high = "green"    )
                   , list( vname =           "wue", low = "purple"   , high = "green"    )
                   , list( vname =      "leaf.gbw", low = "orangered", high = "blue"     )
                   , list( vname =      "leaf.gsw", low = "orangered", high = "blue"     )
                   , list( vname =      "leaf.par", low = "blue"     , high = "orangered")
                   , list( vname =   "leaf.rshort", low = "blue"     , high = "orangered")
                   , list( vname =    "leaf.rlong", low = "blue"     , high = "orangered")
                   , list( vname =    "f.bstorage", low = "purple"   , high = "green"    )
                   , list( vname =      "f.bseeds", low = "purple"   , high = "green"    )
                   )#end list



      #------------------------------------------------------------------------------------#
      #     Replace the list by a data frame.                                              #
      #------------------------------------------------------------------------------------#
      scen.comp = data.frame( apply(X=sapply(X=scen.comp,FUN=c),MARGIN=1,FUN=unlist)
                            , stringsAsFactors = FALSE
                            )#end data.frame
      #------------------------------------------------------------------------------------#
   #---------------------------------------------------------------------------------------#






   #---------------------------------------------------------------------------------------#
   #     Fill in the box plot list with information brought from scen.ts.                  #
   #---------------------------------------------------------------------------------------#
   #----- Find the variable names to be added. --------------------------------------------#
   which.names = names(scen.ts)
   keep        = ! ( which.names %in% names(scen.comp))
   which.names = which.names[keep]
   #----- Match the lists based on vname. -------------------------------------------------#
   comp.idx = match(scen.comp$vname,scen.ts$vname)
   #----- Look for variables that weren't defined in scen.ts. -----------------------------#
   comp.miss = is.na(comp.idx)
   if (any(comp.miss)){
      cat(" - The following variables in scen.comp are missing from scen.ts:","\n")
      cat(paste("   * ",scen.comp$vname[comp.miss],sep=""),sep="\n")
      stop(" - All variables defined in scen.comp must be defined in scen.ts!!!")
   }#end if(any(x.sel))
   #---------------------------------------------------------------------------------------#



   #----- Append the information to the data frame. ---------------------------------------#
   for (wn in which.names) scen.comp[[wn]] = scen.ts[[wn]][comp.idx]
   #---------------------------------------------------------------------------------------#
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This list tells which variables to do the box plot comparing the scenarios, each     #
# plot having the scenario[[1]] in the X scale, and scenario 2 in the group of box plots   #
# Panels are panels.                                                                       #
#                                                                                          #
# IMPORTANT:  All variables here MUST come from one of the variables defined in scen.ts.   #
#------------------------------------------------------------------------------------------#
   #----- All that we need here is the variable name. -------------------------------------#
   panel.box = data.frame( vname = c(          "agb",          "lai",           "ba"
                                    ,          "gpp",          "npp",          "mco"
                                    ,          "cba",        "ldrop",     "bstorage"
                                    ,       "bseeds",      "fs.open",         "mort"
                                    ,      "ncbmort",       "dimort",         "recr"
                                    ,       "growth",     "agb.mort",  "agb.ncbmort"
                                    ,   "agb.dimort",     "agb.recr",   "agb.growth"
                                    ,       "cbarel",       "demand",       "supply"
                                    ,         "rain","water.deficit","malhi.deficit"
                                    ,       "runoff",     "atm.temp",    "leaf.temp"
                                    ,       "rshort",        "rlong",      "atm.vpd"
                                    ,     "leaf.vpd",          "paw",        "smpot"
                                    ,          "nep",         "reco",  "fast.soil.c"
                                    ,"struct.soil.c",  "slow.soil.c",     "het.resp"
                                    ,   "plant.resp",       "hflxlc",       "wflxlc"
                                    ,       "transp",       "hflxgc",       "hflxca"
                                    ,       "wflxgc",       "wflxca",     "i.hflxlc"
                                    ,     "i.wflxlc",     "i.transp",        "i.gpp"
                                    ,        "i.npp", "i.plant.resp",        "i.cba"
                                    ,     "i.cbamax",          "wue",     "leaf.gbw"
                                    ,     "leaf.gsw",     "leaf.par",  "leaf.rshort"
                                    ,        "f.gpp",        "f.npp",        "f.cba"
                                    ,   "f.bstorage",     "f.bseeds",      "f.bleaf"
                                    ,       "nplant"
                                    )#end vname
                         , stringsAsFactors = FALSE
                         )#end data.frame
   #---------------------------------------------------------------------------------------#






   #---------------------------------------------------------------------------------------#
   #     Fill in the box plot list with information brought from scen.ts.                  #
   #---------------------------------------------------------------------------------------#
   #----- Find the variable names to be added. --------------------------------------------#
   which.names = names(scen.ts)
   keep        = ! ( which.names %in% names(panel.box))
   which.names = which.names[keep]
   #----- Match the lists based on vname. -------------------------------------------------#
   comp.idx = match(panel.box$vname,scen.ts$vname)
   #----- Look for variables that weren't defined in scen.ts. -----------------------------#
   comp.miss = is.na(comp.idx)
   if (any(comp.miss)){
      cat(" - The following variables in panel.box are missing from scen.ts:","\n")
      cat(paste("   * ",panel.box$vname[comp.miss],sep=""),sep="\n")
      stop(" - All variables defined in panel.box must be defined in scen.ts!!!")
   }#end if(any(x.sel))
   #---------------------------------------------------------------------------------------#



   #----- Append the information to the data frame. ---------------------------------------#
   for (wn in which.names) panel.box[[wn]] = scen.ts[[wn]][comp.idx]
   #---------------------------------------------------------------------------------------#
#==========================================================================================#
#==========================================================================================#




#==========================================================================================#
#==========================================================================================#
#     Turn the variables that matter global.                                               #
#------------------------------------------------------------------------------------------#
scen.ts       <<- scen.ts
scen.szpft    <<- scen.szpft
scen.barplot  <<- scen.barplot
scen.xyz      <<- scen.xyz 
scen.comp     <<- scen.comp
panel.box     <<- panel.box
nscen.ts      <<- nrow(scen.ts)
nscen.szpft   <<- nrow(scen.szpft   )
nscen.barplot <<- nrow(scen.barplot )
nscen.xvar    <<- nrow(scen.xyz$xvar)
nscen.yvar    <<- nrow(scen.xyz$yvar)
nscen.zvar    <<- nrow(scen.xyz$zvar)
nscen.comp    <<- nrow(scen.comp    )
npanel.box    <<- nrow(panel.box    )
#==========================================================================================#
#==========================================================================================#
