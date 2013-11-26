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
   scen.ts[[n]] = list( vname  = "agb"
                      , desc   = "Above ground biomass"
                      , unit   = "kgC/m2"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end lis
   n            = n + 1
   scen.ts[[n]] = list( vname  = "bgb"
                      , desc   = "Below ground biomass"
                      , unit   = "kgC/m2"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end lis
   n            = n + 1
   scen.ts[[n]] = list( vname  = "biomass"
                      , desc   = "Total biomass"
                      , unit   = "kgC/m2"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end lis
   n            = n + 1
   scen.ts[[n]] = list( vname  = "lai"
                      , desc   = "Leaf area index"
                      , unit   = "m2/m2"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "nplant"
                      , desc   = "Plant density"
                      , unit   = "plant/m2"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end lis
   n            = n + 1
   scen.ts[[n]] = list( vname  = "gpp"
                      , desc   = "Gross primary productivity"
                      , unit   = "kgC/m2/yr"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "npp"
                      , desc   = "Net primary productivity"
                      , unit   = "kgC/m2/yr"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "mco"
                      , desc   = "Maintenance costs"
                      , unit   = "kgC/m2/yr"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "cba"
                      , desc   = "Carbon balance"
                      , unit   = "kgC/m2/yr"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "ldrop"
                      , desc   = "Leaf drop"
                      , unit   = "kgC/m2/yr"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "bstorage"
                      , desc   = "Storage biomass"
                      , unit   = "kgC/m2"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = FALSE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "ba"
                      , desc   = "Basal area"
                      , unit   = "cm2/m2"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "fs.open"
                      , desc   = "Soil moisture stress factor"
                      , unit   = "--"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "mort"
                      , desc   = "Mortality rate"
                      , unit   = "%pop/yr"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = TRUE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "recr"
                      , desc   = "Recruitment rate"
                      , unit   = "%pop/yr"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = TRUE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "growth"
                      , desc   = "Growth rate (DBH)"
                      , unit   = "%DBH/yr"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "ncbmort"
                      , desc   = "Mortality rate - Neg. C balance"
                      , unit   = "%pop/yr"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = TRUE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "dimort"
                      , desc   = "Mortality rate - Density-independent"
                      , unit   = "%pop/yr"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = TRUE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "cbamax"
                      , desc   = "Maximum carbon balance"
                      , unit   = "kgC/m2/yr"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "cbarel"
                      , desc   = "Relative carbon balance"
                      , unit   = "--"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "demand"
                      , desc   = "Water demand"
                      , unit   = "kg/m2/day"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "supply"
                      , desc   = "Water supply"
                      , unit   = "kg/m2/day"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "rain"
                      , desc   = "Precipitation"
                      , unit   = "mm"
                      , f.aggr = "sum"
                      , pftvar = FALSE
                      , dbhvar = FALSE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "water.deficit"
                      , desc   = "Water deficit (ED-2.2)"
                      , unit   = "mm"
                      , f.aggr = "max"
                      , pftvar = FALSE
                      , dbhvar = FALSE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "malhi.deficit"
                      , desc   = "Water deficit (Malhi 2009)"
                      , unit   = "mm"
                      , f.aggr = "max"
                      , pftvar = FALSE
                      , dbhvar = FALSE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "atm.temp"
                      , desc   = "Air temperature"
                      , unit   = "degC"
                      , f.aggr = "mean"
                      , pftvar = FALSE
                      , dbhvar = FALSE
                      , plog   = FALSE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "leaf.temp"
                      , desc   = "Leaf temperature"
                      , unit   = "degC"
                      , f.aggr = "mean"
                      , pftvar = FALSE
                      , dbhvar = FALSE
                      , plog   = FALSE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "rshort"
                      , desc   = "Incoming shortwave radiation"
                      , unit   = "W/m2"
                      , f.aggr = "mean"
                      , pftvar = FALSE
                      , dbhvar = FALSE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "rlong"
                      , desc   = "Incoming longwave radiation"
                      , unit   = "W/m2"
                      , f.aggr = "mean"
                      , pftvar = FALSE
                      , dbhvar = FALSE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "atm.vpd"
                      , desc   = "Air vapour pressure deficit"
                      , unit   = "Pa"
                      , f.aggr = "mean"
                      , pftvar = FALSE
                      , dbhvar = FALSE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "leaf.vpd"
                      , desc   = "Leaf vapour pressure deficit"
                      , unit   = "Pa"
                      , f.aggr = "mean"
                      , pftvar = FALSE
                      , dbhvar = FALSE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "paw"
                      , desc   = "Potential Available Water"
                      , unit   = "%"
                      , f.aggr = "mean"
                      , pftvar = FALSE
                      , dbhvar = FALSE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "smpot"
                      , desc   = "Soil Matric Potential"
                      , unit   = "MPa"
                      , f.aggr = "mean"
                      , pftvar = FALSE
                      , dbhvar = FALSE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "nep"
                      , desc   = "Net Ecosystem Productivity"
                      , unit   = "kgC/m2/yr"
                      , f.aggr = "mean"
                      , pftvar = FALSE
                      , dbhvar = FALSE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "reco"
                      , desc   = "Ecosystem Respiration"
                      , unit   = "kgC/m2/yr"
                      , f.aggr = "mean"
                      , pftvar = FALSE
                      , dbhvar = FALSE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "fast.soil.c"
                      , desc   = "Fast soil carbon"
                      , unit   = "kgC/m2"
                      , f.aggr = "mean"
                      , pftvar = FALSE
                      , dbhvar = FALSE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "struct.soil.c"
                      , desc   = "Structural soil carbon"
                      , unit   = "kgC/m2"
                      , f.aggr = "mean"
                      , pftvar = FALSE
                      , dbhvar = FALSE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "slow.soil.c"
                      , desc   = "Slow soil carbon"
                      , unit   = "kgC/m2"
                      , f.aggr = "mean"
                      , pftvar = FALSE
                      , dbhvar = FALSE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "het.resp"
                      , desc   = "Heterotrophic Respiration"
                      , unit   = "kgC/m2/yr"
                      , f.aggr = "mean"
                      , pftvar = FALSE
                      , dbhvar = FALSE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "plant.resp"
                      , desc   = "Plant Respiration"
                      , unit   = "kgC/m2/yr"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "hflxlc"
                      , desc   = "Leaf sensible heat"
                      , unit   = "W/m2"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "wflxlc"
                      , desc   = "Leaf Evaporation"
                      , unit   = "kg/m2/day"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "transp"
                      , desc   = "Leaf Transpiration"
                      , unit   = "kg/m2/day"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "i.gpp"
                      , desc   = "Mean Gross Primary Production"
                      , unit   = "kgC/plant/yr"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "i.npp"
                      , desc   = "Mean Net Primary Production"
                      , unit   = "kgC/plant/yr"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "i.plant.resp"
                      , desc   = "Mean Plant Respiration"
                      , unit   = "kgC/plant/yr"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "i.cba"
                      , desc   = "Mean Carbon balance"
                      , unit   = "kgC/plant/yr"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "i.cbamax"
                      , desc   = "Mean Maximum C balance"
                      , unit   = "kgC/plant/yr"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "i.hflxlc"
                      , desc   = "Mean Leaf sensible heat flux"
                      , unit   = "W/plant"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "i.wflxlc"
                      , desc   = "Mean Leaf evaporation"
                      , unit   = "kg/plant/day"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "i.transp"
                      , desc   = "Mean Transpiration"
                      , unit   = "kg/plant/day"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "hflxgc"
                      , desc   = "Sensible heat - Gnd->CAS"
                      , unit   = "W/m2"
                      , f.aggr = "mean"
                      , pftvar = FALSE
                      , dbhvar = FALSE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "hflxca"
                      , desc   = "Sensible heat - CAS->ATM"
                      , unit   = "W/m2"
                      , f.aggr = "mean"
                      , pftvar = FALSE
                      , dbhvar = FALSE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "wflxgc"
                      , desc   = "Water flux - Gnd->CAS"
                      , unit   = "kg/m2/day"
                      , f.aggr = "mean"
                      , pftvar = FALSE
                      , dbhvar = FALSE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "wflxca"
                      , desc   = "Water flux - CAS->ATM"
                      , unit   = "kg/m2/day"
                      , f.aggr = "mean"
                      , pftvar = FALSE
                      , dbhvar = FALSE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "wue"
                      , desc   = "Water use efficiency"
                      , unit   = "gC/kgH2O"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "leaf.gbw"
                      , desc   = "Leaf Bnd. Lyr. Conductance"
                      , unit   = "kg/m2/day"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "leaf.gsw"
                      , desc   = "Leaf stomatal Conductance"
                      , unit   = "kg/m2/day"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "f.gpp"
                      , desc   = "Gross Primary Productivity"
                      , unit   = "%biomass/yr"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "f.plant.resp"
                      , desc   = "Plant respiration"
                      , unit   = "%biomass/yr"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "f.npp"
                      , desc   = "Net Primary Productivity"
                      , unit   = "%biomass/yr"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "f.cba"
                      , desc   = "Carbon balance"
                      , unit   = "%biomass/yr"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "f.bstorage"
                      , desc   = "Relative storage biomass"
                      , unit   = "gC_st/kgC_bio"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "f.bleaf"
                      , desc   = "Relative leaf biomass"
                      , unit   = "gC_leaf/kgC_bio"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "f.broot"
                      , desc   = "Relative root biomass"
                      , unit   = "gC_root/kgC_bio"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "f.bseeds"
                      , desc   = "Relative seed biomass"
                      , unit   = "gC_seed/kgC_bio"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "leaf.par"
                      , desc   = "Absorbed PAR - Leaf"
                      , unit   = "umol/m2leaf/s"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "leaf.rshort"
                      , desc   = "Absorbed SW - Leaf"
                      , unit   = "W/m2leaf"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname  = "leaf.rlong"
                      , desc   = "Absorbed LW - Leaf"
                      , unit   = "W/m2leaf"
                      , f.aggr = "mean"
                      , pftvar = TRUE
                      , dbhvar = TRUE
                      , mort   = FALSE
                      , recr   = FALSE
                      , plog   = FALSE
                      , plt    = TRUE
                      )#end list
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Replace the list by a data frame.                                                 #
   #---------------------------------------------------------------------------------------#
   scen.ts = data.frame( apply( X = sapply(X=scen.ts,FUN=c), MARGIN = 1, FUN = unlist )
                       , stringsAsFactors = FALSE
                       )#end data.frame
   for (nl in c("pftvar","dbhvar","mort","recr","plog","plt")){
      scen.ts[[nl]] = as.logical(scen.ts[[nl]])
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
                                     ,       "dimort",       "growth",          "gpp"
                                     ,          "npp",   "plant.resp",          "cba"
                                     ,          "mco",     "bstorage",      "fs.open"
                                     ,       "supply",       "demand",       "hflxlc"
                                     ,       "wflxlc",       "transp",        "i.gpp"
                                     , "i.plant.resp",        "i.npp",        "i.cba"
                                     ,     "i.transp",       "cbamax",     "i.hflxlc"
                                     ,     "i.wflxlc",          "wue",     "i.cbamax"
                                     ,     "leaf.gsw",     "leaf.gbw",        "f.gpp"
                                     ,        "f.npp",        "f.cba",   "f.bstorage"
                                     ,      "f.bleaf",     "f.bseeds",         "rain"
                                     ,"water.deficit","malhi.deficit",      "atm.vpd"
                                     ,     "leaf.vpd",    "leaf.temp",     "atm.temp"
                                     ,     "leaf.par"
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
                                       ,    "i.transp"
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
                            , list( vname = "atm.temp"     , leg        = "right"  )
                            , list( vname = "leaf.temp"    , leg        = "right"  )
                            , list( vname = "leaf.par"     , leg        = "right"  )
                            , list( vname = "leaf.vpd"     , leg        = "right"  )
                            , list( vname = "rain"         , leg        = "right"  )
                            , list( vname = "smpot"        , leg        = "right"  )
                            , list( vname = "water.deficit", leg        = "right"  )
                            , list( vname = "malhi.deficit", leg        = "right"  )
                            , list( vname = "i.cbamax"     , leg        = "right"  )
                            , list( vname = "cbarel"       , leg        = "right"  )
                            , list( vname = "bstorage"     , leg        = "right"  )
                            , list( vname = "f.bstorage"   , leg        = "right"  )
                            , list( vname = "leaf.gsw"     , leg        = "right"  )
                            )#end list
   scen.xyz$yvar      = list( list( vname = "recr"         , leg        = "top"    )
                            , list( vname = "mort"         , leg        = "top"    )
                            , list( vname = "ncbmort"      , leg        = "top"    )
                            , list( vname = "dimort"       , leg        = "top"    )
                            , list( vname = "growth"       , leg        = "top"    )
                            , list( vname = "gpp"          , leg        = "top"    )
                            , list( vname = "npp"          , leg        = "top"    )
                            , list( vname = "i.gpp"        , leg        = "top"    )
                            , list( vname = "i.npp"        , leg        = "top"    )
                            , list( vname = "cba"          , leg        = "top"    )
                            , list( vname = "mco"          , leg        = "top"    )
                            , list( vname = "wue"          , leg        = "top"    )
                            , list( vname = "i.transp"     , leg        = "top"    )
                            , list( vname = "i.hflxlc"     , leg        = "top"    )
                            , list( vname = "leaf.gsw"     , leg        = "top"    )
                            )#end list
   scen.xyz$zvar      = list( list( vname = "lai"          , col.scheme = "clife" )
                            , list( vname = "ba"           , col.scheme = "clife" )
                            , list( vname = "agb"          , col.scheme = "clife" )
                            , list( vname = "nplant"       , col.scheme = "clife" )
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
                   , list( vname =           "gpp", low = "purple"   , high = "green"    )
                   , list( vname =           "npp", low = "purple"   , high = "green"    )
                   , list( vname =           "mco", low = "purple"   , high = "green"    )
                   , list( vname =           "cba", low = "purple"   , high = "green"    )
                   , list( vname =         "ldrop", low = "green"    , high = "purple"   )
                   , list( vname =      "bstorage", low = "purple"   , high = "green"    )
                   , list( vname =       "fs.open", low = "orangered", high = "blue"     )
                   , list( vname =          "mort", low = "green"    , high = "purple"   )
                   , list( vname =          "recr", low = "purple"   , high = "green"    )
                   , list( vname =        "growth", low = "purple"   , high = "green"    )
                   , list( vname =       "ncbmort", low = "green"    , high = "purple"   )
                   , list( vname =        "dimort", low = "green"    , high = "purple"   )
                   , list( vname =        "cbarel", low = "purple"   , high = "green"    )
                   , list( vname =        "demand", low = "blue"     , high = "orangered")
                   , list( vname =        "supply", low = "orangered", high = "blue"     )
                   , list( vname =          "rain", low = "orangered", high = "blue"     )
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
                                    ,      "fs.open",         "mort",         "recr"
                                    ,       "growth",      "ncbmort",       "dimort"
                                    ,       "cbarel",       "demand",       "supply"
                                    ,         "rain","water.deficit","malhi.deficit"
                                    ,     "atm.temp",    "leaf.temp",       "rshort"
                                    ,        "rlong",      "atm.vpd",     "leaf.vpd"
                                    ,          "paw",        "smpot",          "nep"
                                    ,         "reco",  "fast.soil.c","struct.soil.c"
                                    ,  "slow.soil.c",     "het.resp",   "plant.resp"
                                    ,       "hflxlc",       "wflxlc",       "transp"
                                    ,       "hflxgc",       "hflxca",       "wflxgc"
                                    ,       "wflxca",     "i.hflxlc",     "i.wflxlc"
                                    ,     "i.transp",        "i.gpp",        "i.npp"
                                    , "i.plant.resp",        "i.cba",     "i.cbamax"
                                    ,          "wue",     "leaf.gbw",     "leaf.gsw"
                                    ,     "leaf.par",  "leaf.rshort",        "f.gpp"
                                    ,        "f.npp",        "f.cba",   "f.bstorage"
                                    ,     "f.bseeds",      "f.bleaf"
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
