#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#     List of possible plots. In case you don't want some of them, simply switch plt to F. #
#------------------------------------------------------------------------------------------#
#----- Time series per PFT. ---------------------------------------------------------------#
n                 = 0
tspftdbh          = list()
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "agb"
                        , desc     = "Above ground biomass"
                        , e.unit   = "kgC/m2"
                        , i.unit   = "kgC/plant"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = TRUE
                        , bar.plot = TRUE
                        , stack    = TRUE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "bgb"
                        , desc     = "Below ground biomass"
                        , e.unit   = "kgC/m2"
                        , i.unit   = "kgC/plant"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = TRUE
                        , bar.plot = TRUE
                        , stack    = TRUE
                        , scsout   = FALSE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "biomass"
                        , desc     = "Total biomass"
                        , e.unit   = "kgC/m2"
                        , i.unit   = "kgC/plant"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = TRUE
                        , bar.plot = TRUE
                        , stack    = TRUE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "ba"
                        , desc     = "Basal area"
                        , e.unit   = "cm2/m2"
                        , i.unit   = "cm2/plant"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = TRUE
                        , bar.plot = TRUE
                        , stack    = TRUE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "lai"
                        , desc     = "Leaf area index"
                        , e.unit   = "m2/m2"
                        , i.unit   = "m2/m2"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = TRUE
                        , bar.plot = TRUE
                        , stack    = TRUE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "wai"
                        , desc     = "Wood area index"
                        , e.unit   = "m2/m2"
                        , i.unit   = "m2/m2"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = TRUE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "tai"
                        , desc     = "Total area index"
                        , e.unit   = "m2/m2"
                        , i.unit   = "m2/m2"
                        , plog     = FALSE
                        , pft      = FALSE
                        , pftdbh   = FALSE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = TRUE
                        , scsout   = FALSE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "gpp"
                        , desc     = "Gross primary productivity"
                        , e.unit   = "kgC/m2/yr"
                        , i.unit   = "kgC/plant/yr"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = TRUE
                        , bar.plot = TRUE
                        , stack    = TRUE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "npp"
                        , desc     = "Net primary productivity"
                        , e.unit   = "kgC/m2/yr"
                        , i.unit   = "kgC/plant/yr"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = FALSE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "mco"
                        , desc     = "Maintenance costs"
                        , e.unit   = "kgC/m2/yr"
                        , i.unit   = "kgC/plant/yr"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "cba"
                        , desc     = "Carbon balance"
                        , e.unit   = "kgC/m2/yr"
                        , i.unit   = "kgC/plant/yr"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "ldrop"
                        , desc     = "Leaf drop"
                        , e.unit   = "kgC/m2/yr"
                        , i.unit   = "kgC/plant/yr"
                        , plog     = FALSE
                        , pft      = TRUE 
                        , pftdbh   = TRUE 
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "balive"
                        , desc     = "Biomass of active tissues"
                        , e.unit   = "kgC/m2"
                        , i.unit   = "kgC/plant"
                        , plog     = FALSE
                        , pft      = FALSE
                        , pftdbh   = FALSE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = TRUE
                        , scsout   = FALSE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "bdead"
                        , desc     = "Structural biomass"
                        , e.unit   = "kgC/m2"
                        , i.unit   = "kgC/plant"
                        , plog     = FALSE
                        , pft      = FALSE
                        , pftdbh   = FALSE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = TRUE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "bleaf"
                        , desc     = "Leaf biomass"
                        , e.unit   = "kgC/m2"
                        , i.unit   = "kgC/plant"
                        , plog     = FALSE
                        , pft      = FALSE
                        , pftdbh   = FALSE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = TRUE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "broot"
                        , desc     = "Root biomass"
                        , e.unit   = "kgC/m2"
                        , i.unit   = "kgC/plant"
                        , plog     = FALSE
                        , pft      = FALSE
                        , pftdbh   = FALSE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = TRUE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "bsapwood"
                        , desc     = "Sapwood biomass"
                        , e.unit   = "kgC/m2"
                        , i.unit   = "kgC/plant"
                        , plog     = FALSE
                        , pft      = FALSE
                        , pftdbh   = FALSE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = TRUE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "bstorage"
                        , desc     = "Storage biomass"
                        , e.unit   = "kgC/m2"
                        , i.unit   = "kgC/plant"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = TRUE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "bseeds"
                        , desc     = "Seed biomass"
                        , e.unit   = "kgC/m2"
                        , i.unit   = "kgC/plant"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = TRUE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "leaf.resp"
                        , desc     = "Leaf respiration"
                        , e.unit   = "kgC/m2/yr"
                        , i.unit   = "kgC/plant/yr"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = TRUE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "root.resp"
                        , desc     = "Root respiration"
                        , e.unit   = "kgC/m2/yr"
                        , i.unit   = "kgC/plant/yr"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = TRUE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "growth.resp"
                        , desc     = "Growth respiration"
                        , e.unit   = "kgC/m2/yr"
                        , i.unit   = "kgC/plant/yr"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = TRUE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "nplant"
                        , desc     = "Plant density"
                        , e.unit   = "plant/m2"
                        , i.unit   = "plant/m2"
                        , plog     = TRUE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = TRUE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "fs.open"
                        , desc     = "Soil moisture stress factor"
                        , e.unit   = "--"
                        , i.unit   = "--"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "leaf.gsw"
                        , desc     = "Stomatal conductance"
                        , e.unit   = "kg/m2leaf/day"
                        , i.unit   = "kg/m2leaf/day"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = TRUE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "leaf.gbw"
                        , desc     = "Leaf boundary layer conductance"
                        , e.unit   = "kg/m2leaf/day"
                        , i.unit   = "kg/m2leaf/day"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = TRUE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "mort"
                        , desc     = "Mortality rate"
                        , e.unit   = "%pop/yr"
                        , i.unit   = "%pop/yr"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = TRUE
                        , stack    = FALSE
                        , scsout   = FALSE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "recr"
                        , desc     = "Recruitment rate"
                        , e.unit   = "%pop/yr"
                        , i.unit   = "%pop/yr"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = TRUE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "growth"
                        , desc     = "Growth rate"
                        , e.unit   = "%dbh/yr"
                        , i.unit   = "%dbh/yr"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = TRUE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "census.lai"
                        , desc     = "Leaf area index (H > 1.5m)"
                        , e.unit   = "m2/m2"
                        , i.unit   = "m2/m2"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = TRUE
                        , scsout   = FALSE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "census.agb"
                        , desc     = "Above-ground biomass (DBH > 10cm)"
                        , sas      = FALSE
                        , e.unit   = "m2/m2"
                        , i.unit   = "m2/m2"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = TRUE
                        , scsout   = FALSE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "census.ba"
                        , desc     = "Basal area (DBH > 10cm)"
                        , sas      = FALSE
                        , e.unit   = "cm2/m2"
                        , i.unit   = "cm2/m2"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = TRUE
                        , scsout   = FALSE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "ncbmort"
                        , desc     = "Mortality rate - Neg. C balance"
                        , e.unit   = "%pop/yr"
                        , i.unit   = "%pop/yr"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = TRUE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "dimort"
                        , desc     = "Mortality rate - Density independent"
                        , e.unit   = "%pop/yr"
                        , i.unit   = "%pop/yr"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = TRUE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "cbamax"
                        , desc     = "Maximum Carbon balance"
                        , e.unit   = "kgC/m2/yr"
                        , i.unit   = "kgC/plant/yr"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = FALSE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "cbalight"
                        , desc     = "Carbon balance - Maximum light"
                        , e.unit   = "kgC/m2/yr"
                        , i.unit   = "kgC/plant/yr"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = FALSE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "cbamoist"
                        , desc     = "Carbon balance - Maximum moisture"
                        , e.unit   = "kgC/m2/yr"
                        , i.unit   = "kgC/plant/yr"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = FALSE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "cbarel"
                        , desc     = "Relative carbon balance"
                        , e.unit   = "--"
                        , i.unit   = "--"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = TRUE
                        , bar.plot = TRUE
                        , stack    = FALSE
                        , scsout   = FALSE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "hflxlc"
                        , desc     = "Leaf sensible heat flux"
                        , e.unit   = "W/m2"
                        , i.unit   = "W/plant"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = FALSE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "wflxlc"
                        , desc     = "Leaf evaporation"
                        , e.unit   = "kg/m2/day"
                        , i.unit   = "kg/plant/day"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "transp"
                        , desc     = "Leaf transpiration"
                        , e.unit   = "kg/m2/day"
                        , i.unit   = "kg/plant/day"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = TRUE
                        , stack    = TRUE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "wue"
                        , desc     = "Water use efficiency"
                        , e.unit   = "gC/kgH2O"
                        , i.unit   = "gC/kgH2O"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = TRUE
                        , stack    = FALSE
                        , scsout   = FALSE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "i.gpp"
                        , desc     = "Mean gross primary production"
                        , e.unit   = "kgC/plant/yr"
                        , i.unit   = "kgC/plant/yr"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "i.npp"
                        , desc     = "Mean net primary production"
                        , e.unit   = "kgC/plant/yr"
                        , i.unit   = "kgC/plant/yr"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = FALSE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "i.plant.resp"
                        , desc     = "Mean plant respiration"
                        , e.unit   = "kgC/plant/yr"
                        , i.unit   = "kgC/plant/yr"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "i.mco"
                        , desc     = "Mean maintenance costs"
                        , e.unit   = "kgC/plant/yr"
                        , i.unit   = "kgC/plant/yr"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "i.cba"
                        , desc     = "Mean carbon balance"
                        , e.unit   = "kgC/plant/yr"
                        , i.unit   = "kgC/plant/yr"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "i.cbamax"
                        , desc     = "Maximum carbon balance"
                        , e.unit   = "kgC/plant/yr"
                        , i.unit   = "kgC/plant/yr"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = FALSE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "i.cbalight"
                        , desc     = "Mean carbon balance - maximum light"
                        , e.unit   = "kgC/plant/yr"
                        , i.unit   = "kgC/plant/yr"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = FALSE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "i.cbamoist"
                        , desc     = "Mean carbon balance - maximum moisture"
                        , e.unit   = "kgC/plant/yr"
                        , i.unit   = "kgC/plant/yr"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = FALSE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "i.hflxlc"
                        , desc     = "Mean leaf sensible heat flux"
                        , e.unit   = "W/plant"
                        , i.unit   = "W/plant"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = FALSE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "i.wflxlc"
                        , desc     = "Mean leaf evaporation"
                        , e.unit   = "kg/plant/day"
                        , i.unit   = "kg/plant/day"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "i.transp"
                        , desc     = "Mean leaf evaporation"
                        , e.unit   = "kg/plant/day"
                        , i.unit   = "kg/plant/day"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "height"
                        , desc     = "Height"
                        , e.unit   = "m"
                        , i.unit   = "m"
                        , plog     = FALSE
                        , pft      = FALSE
                        , pftdbh   = FALSE
                        , sas      = TRUE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = FALSE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "f.gpp"
                        , desc     = "Gross Primary Productivity"
                        , e.unit   = "%biomass/yr"
                        , i.unit   = "%biomass/yr"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "f.plant.resp"
                        , desc     = "Plant respiration"
                        , e.unit   = "%biomass/yr"
                        , i.unit   = "%biomass/yr"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "f.npp"
                        , desc     = "Net Primary Productivity"
                        , e.unit   = "%biomass/yr"
                        , i.unit   = "%biomass/yr"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = FALSE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "f.cba"
                        , desc     = "Carbon balance"
                        , e.unit   = "%biomass/yr"
                        , i.unit   = "%biomass/yr"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "f.bstorage"
                        , desc     = "Relative storage biomass"
                        , e.unit   = "gC_st/kgC_bio"
                        , i.unit   = "gC_st/kgC_bio"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "f.bleaf"
                        , desc     = "Relative leaf biomass"
                        , e.unit   = "gC_leaf/kgC_bio"
                        , i.unit   = "gC_leaf/kgC_bio"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "f.broot"
                        , desc     = "Relative root biomass"
                        , e.unit   = "gC_root/kgC_bio"
                        , i.unit   = "gC_root/kgC_bio"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "f.bseeds"
                        , desc     = "Relative seed biomass"
                        , e.unit   = "gC_seed/kgC_bio"
                        , i.unit   = "gC_seed/kgC_bio"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "leaf.par"
                        , desc     = "Absorbed PAR - Leaf"
                        , e.unit   = "umol/m2leaf/s"
                        , i.unit   = "umol/m2leaf/s"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "leaf.rshort"
                        , desc     = "Absorbed SW - Leaf"
                        , e.unit   = "W/m2leaf"
                        , i.unit   = "W/m2leaf"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "leaf.rlong"
                        , desc     = "Net absorbed LW - Leaf"
                        , e.unit   = "W/m2leaf"
                        , i.unit   = "W/m2leaf"
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
#------------------------------------------------------------------------------------------#






#----- Time series per Land use. ----------------------------------------------------------#
n         = 0
tslu      = list()
n         = n + 1
tslu[[n]] = list( vnam = "agb"
                , desc = "Above ground biomass"
                , unit = "kgC/m2"
                , plog = FALSE
                , plt  = FALSE)
n         = n + 1
tslu[[n]] = list( vnam = "bgb"
                , desc = "Below ground biomass"
                , unit = "kgC/m2"
                , plog = FALSE
                , plt  = FALSE)
n         = n + 1
tslu[[n]] = list( vnam = "biomass"
                , desc = "Total biomass"
                , unit = "kgC/m2"
                , plog = FALSE
                , plt  = FALSE)
n         = n + 1
tslu[[n]] = list( vnam = "lai"
                , desc = "Leaf area index"
                , unit = "m2/m2"
                , plog = FALSE
                , plt  = FALSE)
n         = n + 1
tslu[[n]] = list( vnam = "gpp"
                , desc = "Gross primary productivity"
                , unit = "kgC/m2/yr"
                , plog = FALSE
                , plt  = FALSE)
n         = n + 1
tslu[[n]] = list( vnam = "npp"
                , desc = "Net primary productivity"
                , unit = "kgC/m2/yr"
                , plog = FALSE
                , plt  = FALSE)
n         = n + 1
tslu[[n]] = list( vnam = "area"
                , desc = "Fraction of area"
                , unit = ""
                , plog = FALSE
                , plt  = FALSE)
n         = n + 1
tslu[[n]] = list( vnam = "basarea"
                , desc = "Basal area"
                , unit = "cm2/m2"
                , plog = FALSE
                , plt  = FALSE)
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#      This plots distributions of the properties over time, in three different ways.      #
# -- fco.mmean: filled contour, with month in the X axis and year in the y axis            #
# -- fco.qmean: filled contour, with month/year in the X axis and hour in the y axis       #
# -- box.plot: box plot with fixed months and year spread.                                 #
#                                                                                          #
#    Soil variables should not be placed here, because they have two dimensions...         #
# Integrated soil properties can come here, though.                                        #
#------------------------------------------------------------------------------------------#
n            = 0
squeeze      = list()
n            = n + 1
squeeze[[n]] = list( vnam       = "gpp"
                   , desc       = "Gross Primary productivity"
                   , unit       = "kgC/m2/yr"
                   , col.scheme = "atlas"
                   , fco.mmean  = TRUE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "npp"
                   , desc       = "Net Primary productivity"
                   , unit       = "kgC/m2/yr"
                   , col.scheme = "atlas"
                   , fco.mmean  = TRUE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "plant.resp"
                   , desc       = "Plant respiration"
                   , unit       = "kgC/m2/yr"
                   , col.scheme = "muitas"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "het.resp"
                   , desc       = "Heterotrophic respiration"
                   , unit       = "kgC/m2/yr"
                   , col.scheme = "muitas"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "nep"
                   , desc       = "Net ecosystem production"
                   , unit       = "kgC/m2/yr"
                   , col.scheme = "muitas"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "reco"
                   , desc       = "Ecosystem respiration"
                   , unit       = "kgC/m2/yr"
                   , col.scheme = "muitas"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "hflxca"
                   , desc       = "Sensible heat flux"
                   , unit       = "W/m2"
                   , col.scheme = "muitas"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "wflxca"
                   , desc       = "Water flux"
                   , unit       = "kg/m2/day"
                   , col.scheme = "imuitas"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "wflxgc"
                   , desc       = "Ground evaporation"
                   , unit       = "kg/m2/day"
                   , col.scheme = "imuitas"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "wflxlc"
                   , desc       = "Leaf evaporation"
                   , unit       = "kg/m2/day"
                   , col.scheme = "imuitas"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "wflxwc"
                   , desc       = "Wood evaporation"
                   , unit       = "kg/m2/day"
                   , col.scheme = "imuitas"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "evap"
                   , desc       = "Evaporation"
                   , unit       = "kg/m2/day"
                   , col.scheme = "imuitas"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "transp"
                   , desc       = "Transpiration"
                   , unit       = "kg/m2/day"
                   , col.scheme = "imuitas"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "nee"
                   , desc       = "Net ecosystem exchange"
                   , unit       = "umol/m2/s"
                   , col.scheme = "imuitas"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "cflxca"
                   , desc       = "CO2 flux"
                   , unit       = "umol/m2/s"
                   , col.scheme = "imuitas"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "cflxst"
                   , desc       = "CO2 flux"
                   , unit       = "umol/m2/s"
                   , col.scheme = "imuitas"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "can.temp"
                   , desc       = "Canopy air temperature"
                   , unit       = "degC"
                   , col.scheme = "muitas"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "atm.temp"
                   , desc       = "Atmospheric temperature"
                   , unit       = "degC"
                   , col.scheme = "muitas"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "leaf.temp"
                   , desc       = "Leaf temperature"
                   , unit       = "degC"
                   , col.scheme = "muitas"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "wood.temp"
                   , desc       = "Wood temperature"
                   , unit       = "degC"
                   , col.scheme = "muitas"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "gnd.temp"
                   , desc       = "Ground temperature"
                   , unit       = "degC"
                   , col.scheme = "muitas"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "atm.shv"
                   , desc       = "Atmospheric specific humidity" 
                   , unit       = "g/kg"
                   , col.scheme = "imuitas"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "can.shv"
                   , desc       = "Canopy air specific humidity"
                   , unit       = "g/kg"
                   , col.scheme = "imuitas"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "gnd.shv"
                   , desc       = "Ground specific humidity"
                   , unit       = "g/kg"
                   , col.scheme = "imuitas"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "atm.co2"
                   , desc       = "Atmospheric CO2 mixing ratio"
                   , unit       = "umol/mol"
                   , col.scheme = "muitas"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "can.co2"
                   , desc       = "Canopy air CO2 mixing ratio"
                   , unit       = "umol/mol"
                   , col.scheme = "muitas"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "rain"
                   , desc       = "Total monthly precipitation"
                   , unit       = "mm"
                   , col.scheme = "imuitas"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "fs.open"
                   , desc       = "Fraction of open stomata"
                   , unit       = "---"
                   , col.scheme = "imuitas"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "leaf.gbw"
                   , desc       = "Leaf boundary layer conductance"
                   , unit       = "kg/m2l/day"
                   , col.scheme = "imuitas"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "leaf.gsw"
                   , desc       = "Stomatal conductance"
                   , unit       = "kg/m2l/day"
                   , col.scheme = "imuitas"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "rshort"
                   , desc       = "Downward shortwave radiation"
                   , unit       = "W/m2"
                   , col.scheme = "icloudy"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "rshort.gnd"
                   , desc       = "Abs. gnd. shortwave radiation"
                   , unit       = "W/m2"
                   , col.scheme = "icloudy"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "rshortup"
                   , desc       = "Outgoing shortwave radiation"
                   , unit       = "W/m2"
                   , col.scheme = "icloudy"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "rlong"
                   , desc       = "Downward longwave radiation"
                   , unit       = "W/m2"
                   , col.scheme = "cloudy"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "rlong.gnd"
                   , desc       = "Abs. gnd. longwave radiation"
                   , unit       = "W/m2"
                   , col.scheme = "cloudy"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "rlongup"
                   , desc       = "Outgoing longwave radiation"
                   , unit       = "W/m2"
                   , col.scheme = "cloudy"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "albedo"
                   , desc       = "Shortwave albedo"
                   , unit       = "---"
                   , col.scheme = "muitas"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "ustar"
                   , desc       = "Friction velocity"
                   , unit       = "m/s"
                   , col.scheme = "muitas"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "paw"
                   , desc       = "Potential available water"
                   , unit       = "%"
                   , col.scheme = "imuitas"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "smpot"
                   , desc       = "Integrated matric potential"
                   , unit       = "MPa"
                   , col.scheme = "muitas"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#      Theme plots (multiple variables in the same plot).                                  #
#------------------------------------------------------------------------------------------#
n          = 0
theme      = list()
n          = n + 1
theme[[n]] = list( vnam      = c(        "gpp", "plant.resp",  "het.resp",      "reco"
                                ,        "npp",        "nep")
                 , desc      = c(        "GPP","Plant resp.","Het. resp.","Ecos. Resp."
                                ,        "NPP",        "NEP")
                 , colour    = c(  "darkgreen",       "gold",   "purple3",  "orangered"
                                ,"chartreuse3","dodgerblue3")
                 , lwd       = c(          2.5,          2.5,         2.5,          2.5
                                ,          2.5,          2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "ecoflux"
                 , title     = "Ecosystem fluxes"
                 , unit      = "kgC/m2/yr"
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = TRUE
                 , ymean     = TRUE
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(   "rshort",      "rlong","rshort.gnd",    "qwflxca"
                                ,   "hflxca")
                 , desc      = c(  "Down SW",    "Down LW", "Abs. Grnd",   "'Latent'"
                                , "Sensible")
                 , colour    = c("goldenrod","chartreuse4",   "purple4","dodgerblue3"
                                ,"firebrick")
                 , lwd       = c(        2.5,          2.5,         2.5,          2.5
                                        ,2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "eneflux"
                 , title     = "Energy fluxes"
                 , unit      = "W/m2"
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = TRUE
                 , ymean     = TRUE
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(        "wflxgc",         "wflxca",      "wflxlc"
                                ,        "wflxwc",         "transp")
                 , desc      = c("Ground->Canopy",    "Canopy->Atm","Leaf->Canopy"
                                ,  "Wood->Canopy",  "Transpiration")
                 , colour    = c(     "firebrick",     "royalblue4", "chartreuse4"
                                , "darkgoldenrod","darkolivegreen2")
                 , lwd       = c(             2.5,              2.5,           2.5
                                ,             2.5,              2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "h2oflux"
                 , title     = "Water fluxes"
                 , unit      = "kg/m2/day"
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = TRUE
                 , ymean     = TRUE
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(     "atm.temp",  "can.temp",  "leaf.temp"
                                ,    "wood.temp",  "gnd.temp")
                 , desc      = c(   "Atmosphere","Canopy air",       "Leaf"
                                ,         "Wood",    "Ground")
                 , colour    = c(  "deepskyblue",    "grey45","chartreuse4"
                                ,"darkgoldenrod", "orangered")
                 , lwd       = c(            2.5,         2.5,          2.5
                                ,            2.5,         2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "temperature"
                 , title     = "Temperature"
                 , unit      = "degC"
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = TRUE
                 , ymean     = TRUE
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(    "atm.shv",   "can.shv",      "gnd.shv")
                 , desc      = c( "Atmosphere","Canopy air",       "Ground")
                 , colour    = c("deepskyblue", "slateblue","darkgoldenrod")
                 , lwd       = c(          2.5,         2.5,            2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "h2ovapour"
                 , title     = "Water vapour mixing ratio"
                 , unit      = "g/kg"
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = TRUE
                 , ymean     = TRUE
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(    "atm.co2",    "can.co2")
                 , desc      = c( "Atmosphere", "Canopy air")
                 , colour    = c("deepskyblue",  "slateblue")
                 , lwd       = c(2.5,2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "co2"
                 , title     = "CO2 mixing ratio"
                 , unit      = "umol/mol"
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = TRUE
                 , ymean     = TRUE
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(         "rain")
                 , desc      = c("Precipitation")
                 , colour    = c(   "royalblue4")
                 , lwd       = c(            2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "prec"
                 , title     = "Precipitation"
                 , unit      = "mm"
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = TRUE
                 , ymean     = TRUE
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c("npat.global")
                 , desc      = c("Patch count")
                 , colour    = c(  "orangered")
                 , lwd       = c(          2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "npatch"
                 , title     = "Total number of patches"
                 , unit      = "---"
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = FALSE
                 , ymean     = TRUE
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c( "ncoh.global")
                 , desc      = c("Cohort count")
                 , colour    = c( "chartreuse4")
                 , lwd       = c(           2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "ncohort"
                 , title     = "Total number of cohorts"
                 , unit      = "---"
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = FALSE
                 , ymean     = TRUE
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c( "plant.resp",     "het.resp",   "cwd.resp", "leaf.resp"
                                ,  "root.resp",  "growth.resp",       "reco")
                 , desc      = c("Plant resp.",   "Het. resp.",  "CWD resp.","Leaf resp."
                                , "Root resp.", "Growth resp.","Ecos. resp.")
                 , colour    = c("chartreuse3",      "purple3","saddlebrown", "darkgreen"
                                ,"orangered"  ,"darkgoldenrod","dodgerblue3")
                 , lwd       = c(          2.5,            2.5,          2.5,         2.5
                                ,          2.5,            2.5,          2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "respiration"
                 , title     = "Respiration fluxes"
                 , unit      = "kgC/m2/yr"
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = TRUE
                 , ymean     = TRUE
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(  "atm.vels",            "ustar")
                 , desc      = c("Wind speed","Friction velocity")
                 , colour    = c("deepskyblue",       "slateblue")
                 , lwd       = c(          2.5,               2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "wind"
                 , title     = "Winds"
                 , unit      = "m/s"
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = TRUE
                 , ymean     = TRUE
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c("fast.soil.c","slow.soil.c","struct.soil.c")
                 , desc      = c(       "Fast",       "Slow",   "Structural")
                 , colour    = c("chartreuse3","dodgerblue3",    "orangered")
                 , lwd       = c(          2.5,          2.5,            2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "soil_carbon"
                 , title     = "Soil Carbon"
                 , unit      = "kgC/m2"
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = FALSE
                 , ymean     = TRUE
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(    "atm.vpd",    "can.vpd",   "leaf.vpd")
                 , desc      = c( "Atmosphere", "Canopy air",       "Leaf")
                 , colour    = c("deepskyblue","dodgerblue4","chartreuse3")
                 , lwd       = c(          2.5,          2.5,          2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "vpdef"
                 , title     = "Vapour pressure deficit"
                 , unit      = "hPa"
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = TRUE
                 , ymean     = TRUE
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(         "paw")
                 , desc      = c("Pot.Av.Water")
                 , colour    = c(   "steelblue")
                 , lwd       = c(           2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "paw"
                 , title     = "Average potential available water"
                 , unit      = "%"
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = FALSE
                 , ymean     = TRUE
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(         "smpot")
                 , desc      = c("Neg. Potential")
                 , colour    = c(    "royalblue4")
                 , lwd       = c(             2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "smpot"
                 , title     = "Average soil matric potential"
                 , unit      = "MPa"
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = FALSE
                 , ymean     = TRUE
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c("water.deficit",      "malhi.deficit")
                 , desc      = c(       "ED-2.2","Malhi et al. (2009)")
                 , colour    = c(    "orangered",               "gold")
                 , lwd       = c(            2.5,                  2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "water_deficit"
                 , title     = "Monthly Water deficit"
                 , unit      = "kg/m2/month"
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = FALSE
                 , ymean     = TRUE
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(        "nee",   "cflxca",     "cflxst")
                 , desc      = c(        "NEE", "CO2 Flux","CO2 Storage")
                 , colour    = c("chartreuse4","steelblue",  "orangered")
                 , lwd       = c(          2.5,        2.5,          2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "carbflux"
                 , title     = "CO2 fluxes"
                 , unit      = "umol/m2/s"
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = TRUE
                 , ymean     = TRUE
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(         "rshort",  "rshort.beam","rshort.diff"
                                ,     "rshort.gnd",     "rshortup")
                 , desc      = c("Down Top canopy",         "Beam",    "Diffuse"
                                ,    "Abs. Ground","Up Top canopy")
                 , colour    = c(    "deepskyblue","darkgoldenrod",     "grey45"
                                ,      "firebrick",   "royalblue3")
                 , lwd       = c(              2.5,            2.5,          2.5
                                ,              2.5,            2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "rshort"
                 , title     = "Short wave radiation"
                 , unit      = "W/m2"
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = TRUE
                 , ymean     = TRUE
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(          "rlong",  "rlongup",    "rlong.gnd")
                 , desc      = c("Down Top canopy","Upward LW",  "Abs. Ground")
                 , colour    = c(    "deepskyblue","orangered","darkgoldenrod")
                 , lwd       = c(              2.5,        2.5,            2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "rlong"
                 , title     = "Long wave radiation"
                 , unit      = "W/m2"
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = TRUE
                 , ymean     = TRUE
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(         "albedo", "albedo.par","albedo.nir")
                 , desc      = c("SW Albedo (Net)", "PAR Albedo","NIR Albedo")
                 , colour    = c(    "deepskyblue","chartreuse3", "orangered")
                 , lwd       = c(              2.5,          2.5,         2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "albedo"
                 , title     = "Albedo"
                 , unit      = "---"
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = TRUE
                 , ymean     = TRUE
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(        "par.tot",     "par.beam", "par.diff"
                                ,        "par.gnd",        "parup")
                 , desc      = c("Down Top canopy",         "Beam",  "Diffuse"
                                ,    "Abs. Ground","Up Top canopy")
                 , colour    = c(    "deepskyblue",  "darkorange3","slateblue"
                                ,  "darkgoldenrod",  "chartreuse3")
                 , lwd       = c(              2.5,            2.5,        2.5
                                ,              2.5,            2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "par"
                 , title     = "Photosynthetically Active Radiation"
                 , unit      = "umol/m2/s"
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = TRUE
                 , ymean     = TRUE
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(       "leaf.gsw",        "leaf.gbw",        "wood.gbw")
                 , desc      = c( "Leaf (Stomata)","Leaf (Bnd. Lyr.)","Wood (Bnd. Lyr.)")
                 , colour    = c(    "chartreuse4",       "steelblue",          "sienna")
                 , lwd       = c(              2.5,               2.5,               2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "conduct"
                 , title     = "Conductance"
                 , unit      = "kg/m2/day"
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = TRUE
                 , ymean     = TRUE
                 )#end list
#------------------------------------------------------------------------------------------#




#----- Comparison between observations and model averages. --------------------------------#
n              = 0
compmodel      = list()
n              = n + 1
compmodel[[n]] = list( vnam   = "nep"
                     , desc   = "Net Ecosystem Productivity"
                     , unit   = "kgC/m2/yr"
                     , plotsd = TRUE
                     , colour = c(green.fg,grey.fg)
                     , errcol = c(green.bg,grey.bg)
                     , angle  = c(           45,     -45)
                     , dens   = c(           40,      40)
                     , lwd    = c(          2.5,     2.5)
                     , shwd   = c(          1.0,     1.0)
                     , type   = "o"
                     , plog   = ""
                     , legpos = "topleft"
                     , mmean  = TRUE
                     , qmean  = TRUE
                     , emean  = TRUE
                     , scsout = FALSE
                     )#end list
n              = n + 1
compmodel[[n]] = list( vnam   = "gpp"
                     , desc   = "Gross Primary Productivity"
                     , unit   = "kgC/m2/yr"
                     , plotsd = TRUE
                     , colour = c(green.fg,grey.fg)
                     , errcol = c(green.bg,grey.bg)
                     , angle  = c(           45,     -45)
                     , dens   = c(           40,      40)
                     , lwd    = c(          2.5,     2.5)
                     , shwd   = c(          1.0,     1.0)
                     , type   = "o"
                     , plog   = ""
                     , legpos = "topleft"
                     , mmean  = TRUE
                     , qmean  = TRUE
                     , emean  = TRUE
                     , scsout = FALSE
                     )#end list
n              = n + 1
compmodel[[n]] = list( vnam   = "atm.co2"
                     , desc   = "Air CO2 mixing ratio"
                     , unit   = "g/kg"
                     , plotsd = TRUE
                     , colour = c(green.fg,grey.fg)
                     , errcol = c(green.bg,grey.bg)
                     , angle  = c(           45,     -45)
                     , dens   = c(           40,      40)
                     , lwd    = c(          2.5,     2.5)
                     , shwd   = c(          1.0,     1.0)
                     , type   = "o"
                     , plog   = ""
                     , legpos = "topleft"
                     , mmean  = TRUE
                     , qmean  = TRUE
                     , emean  = TRUE
                     , scsout = TRUE
                     )#end list
n              = n + 1
compmodel[[n]] = list( vnam   = "hflxca"
                     , desc   = "Sensible heat flux"
                     , unit   = "W/m2"
                     , plotsd = TRUE
                     , colour = c(orange.fg,grey.fg)
                     , errcol = c(orange.bg,grey.bg)
                     , angle  = c(         45,     -45)
                     , dens   = c(         40,      40)
                     , lwd    = c(        2.5,     2.5)
                     , shwd   = c(        1.0,     1.0)
                     , type   = "o"
                     , plog   = ""
                     , legpos = "topleft"
                     , mmean  = TRUE
                     , qmean  = TRUE
                     , emean  = TRUE
                     , scsout = TRUE
                     )#end list
n              = n + 1
compmodel[[n]] = list( vnam   = "wflxca"
                     , desc   = "Water flux"
                     , unit   = "kgH2O/m2/day"
                     , plotsd = TRUE
                     , colour = c(blue.fg,grey.fg)
                     , errcol = c(blue.bg,grey.bg)
                     , angle  = c(          45,     -45)
                     , dens   = c(          40,      40)
                     , lwd    = c(         2.5,     2.5)
                     , shwd   = c(         1.0,     1.0)
                     , type   = "o"
                     , plog   = ""
                     , legpos = "topleft"
                     , mmean  = TRUE
                     , qmean  = TRUE
                     , emean  = TRUE
                     , scsout = TRUE
                     )#end list
n              = n + 1
compmodel[[n]] = list( vnam   = "rshort"
                     , desc   = "Incoming shortwave radiation"
                     , unit   = "W/m2"
                     , plotsd = TRUE
                     , colour = c(orange.fg,grey.fg)
                     , errcol = c(orange.bg,grey.bg)
                     , angle  = c(         45,     -45)
                     , dens   = c(         40,      40)
                     , lwd    = c(        2.5,     2.5)
                     , shwd   = c(        1.0,     1.0)
                     , type   = "o"
                     , plog   = ""
                     , legpos = "topleft"
                     , mmean  = TRUE
                     , qmean  = TRUE
                     , emean  = TRUE
                     , scsout = TRUE
                     )#end list
n              = n + 1
compmodel[[n]] = list( vnam   = "par.tot"
                     , desc   = "Incoming PAR"
                     , unit   = "umol/m2/s"
                     , plotsd = TRUE
                     , colour = c(orange.fg,grey.fg)
                     , errcol = c(orange.bg,grey.bg)
                     , angle  = c(         45,     -45)
                     , dens   = c(         40,      40)
                     , lwd    = c(        2.5,     2.5)
                     , shwd   = c(        1.0,     1.0)
                     , type   = "o"
                     , plog   = ""
                     , legpos = "topleft"
                     , mmean  = FALSE
                     , qmean  = FALSE
                     , emean  = FALSE
                     , scsout = TRUE
                     )#end list
n              = n + 1
compmodel[[n]] = list( vnam   = "rlong"
                     , desc   = "Incoming longwave radiation"
                     , unit   = "W/m2"
                     , plotsd = TRUE
                     , colour = c(blue.fg,grey.fg)
                     , errcol = c(blue.bg,grey.bg)
                     , angle  = c(          45,     -45)
                     , dens   = c(          40,      40)
                     , lwd    = c(         2.5,     2.5)
                     , shwd   = c(         1.0,     1.0)
                     , type   = "o"
                     , plog   = ""
                     , legpos = "topleft"
                     , mmean  = TRUE
                     , qmean  = TRUE
                     , emean  = TRUE
                     , scsout = TRUE
                     )#end list
n              = n + 1
compmodel[[n]] = list( vnam   = "atm.temp"
                     , desc   = "Air temperature"
                     , unit   = "C"
                     , plotsd = TRUE
                     , colour = c(orange.fg,grey.fg)
                     , errcol = c(orange.bg,grey.bg)
                     , angle  = c(           45,     -45)
                     , dens   = c(           40,      40)
                     , lwd    = c(          2.5,     2.5)
                     , shwd   = c(          1.0,     1.0)
                     , type   = "o"
                     , plog   = ""
                     , legpos = "topleft"
                     , mmean  = TRUE
                     , qmean  = TRUE
                     , emean  = TRUE
                     , scsout = TRUE
                     )#end list
n              = n + 1
compmodel[[n]] = list( vnam   = "atm.shv"
                     , desc   = "Air specific humidity"
                     , unit   = "g/kg"
                     , plotsd = TRUE
                     , colour = c(blue.fg,grey.fg)
                     , errcol = c(blue.bg,grey.bg)
                     , angle  = c(          45,     -45)
                     , dens   = c(          40,      40)
                     , lwd    = c(         2.5,     2.5)
                     , shwd   = c(         1.0,     1.0)
                     , type   = "o"
                     , plog   = ""
                     , legpos = "topleft"
                     , mmean  = TRUE
                     , qmean  = TRUE
                     , emean  = TRUE
                     , scsout = TRUE
                     )#end list
n              = n + 1
compmodel[[n]] = list( vnam   = "rain"
                     , desc   = "Precipitation rate"
                     , unit   = "kg/m2/day"
                     , plotsd = FALSE
                     , colour = c(blue.fg,grey.fg)
                     , errcol = c(blue.bg,grey.bg)
                     , angle  = c(          45,     -45)
                     , dens   = c(          40,      40)
                     , lwd    = c(         2.5,     2.5)
                     , shwd   = c(         1.0,     1.0)
                     , type   = "o"
                     , plog   = ""
                     , legpos = "topleft"
                     , mmean  = TRUE
                     , qmean  = TRUE
                     , emean  = TRUE
                     , scsout = TRUE
                     )#end list
n              = n + 1
compmodel[[n]] = list( vnam   = "atm.vels"
                     , desc   = "Wind speed"
                     , unit   = "m/s"
                     , plotsd = TRUE
                     , colour = c(indigo.fg,grey.fg)
                     , errcol = c(indigo.bg,grey.bg)
                     , angle  = c(             45,     -45)
                     , dens   = c(             40,      40)
                     , lwd    = c(            2.5,     2.5)
                     , shwd   = c(            1.0,     1.0)
                     , type   = "o"
                     , plog   = ""
                     , legpos = "topleft"
                     , mmean  = TRUE
                     , qmean  = TRUE
                     , emean  = TRUE
                     , scsout = TRUE
                     )#end list
n              = n + 1
compmodel[[n]] = list( vnam   = "atm.prss"
                     , desc   = "Air pressure"
                     , unit   = "hPa"
                     , plotsd = TRUE
                     , colour = c(purple.fg,grey.fg)
                     , errcol = c(purple.bg,grey.bg)
                     , angle  = c(             45,     -45)
                     , dens   = c(             40,      40)
                     , lwd    = c(            2.5,     2.5)
                     , shwd   = c(            1.0,     1.0)
                     , type   = "o"
                     , plog   = ""
                     , legpos = "topleft"
                     , mmean  = TRUE
                     , qmean  = TRUE
                     , emean  = TRUE
                     , scsout = TRUE
                     )#end list
n              = n + 1
compmodel[[n]] = list( vnam   = "cflxca"
                     , desc   = "CO2 flux"
                     , unit   = "umol/m2/s"
                     , plotsd = TRUE
                     , colour = c(green.fg,grey.fg)
                     , errcol = c(green.bg,grey.bg)
                     , angle  = c(           45,     -45)
                     , dens   = c(           40,      40)
                     , lwd    = c(          2.5,     2.5)
                     , shwd   = c(          1.0,     1.0)
                     , type   = "o"
                     , plog   = ""
                     , legpos = "topleft"
                     , mmean  = TRUE
                     , qmean  = TRUE
                     , emean  = TRUE
                     , scsout = TRUE
                     )#end list
n              = n + 1
compmodel[[n]] = list( vnam   = "cflxst"
                     , desc   = "CO2 storage"
                     , unit   = "umol/m2/s"
                     , plotsd = TRUE
                     , colour = c(orange.fg,grey.fg)
                     , errcol = c(orange.bg,grey.bg)
                     , angle  = c(45,-45)
                     , dens   = c(40, 40)
                     , lwd    = c(2.5,2.5)
                     , shwd   = c(1.0,1.0)
                     , type   = "o"
                     , plog   = ""
                     , legpos = "topleft"
                     , mmean  = TRUE
                     , qmean  = TRUE
                     , emean  = TRUE
                     , scsout = TRUE
                     )#end list
n              = n + 1
compmodel[[n]] = list( vnam   = "nee"
                     , desc   = "Net Ecosystem Exchange"
                     , unit   = "umol/m2/s"
                     , plotsd = TRUE
                     , colour = c(green.fg,grey.fg)
                     , errcol = c(green.bg,grey.bg)
                     , angle  = c(           45,     -45)
                     , dens   = c(           40,      40)
                     , lwd    = c(          2.5,     2.5)
                     , shwd   = c(          1.0,     1.0)
                     , type   = "o"
                     , plog   = ""
                     , legpos = "topleft"
                     , mmean  = TRUE
                     , qmean  = TRUE
                     , emean  = TRUE
                     , scsout = FALSE
                     )#end list
n              = n + 1
compmodel[[n]] = list( vnam   = "ustar"
                     , desc   = "Friction velocity"
                     , unit   = "m/s"
                     , plotsd = TRUE
                     , colour = c(purple.fg,grey.fg)
                     , errcol = c(purple.bg,grey.bg)
                     , angle  = c(             45,     -45)
                     , dens   = c(             40,      40)
                     , lwd    = c(            2.5,     2.5)
                     , shwd   = c(            1.0,     1.0)
                     , type   = "o"
                     , plog   = ""
                     , legpos = "topleft"
                     , mmean  = TRUE
                     , qmean  = TRUE
                     , emean  = TRUE
                     , scsout = TRUE
                     )#end list
n              = n + 1
compmodel[[n]] = list( vnam   = "reco"
                     , desc   = "Ecosystem respiration"
                     , unit   = "kgC/m2/yr"
                     , plotsd = TRUE
                     , colour = c(yellow.bg,grey.fg)
                     , errcol = c(yellow.bg,grey.bg)
                     , angle  = c(           45,     -45)
                     , dens   = c(           40,      40)
                     , lwd    = c(          2.5,     2.5)
                     , shwd   = c(          1.0,     1.0)
                     , type   = "o"
                     , plog   = ""
                     , legpos = "topleft"
                     , mmean  = TRUE
                     , qmean  = TRUE
                     , emean  = TRUE
                     , scsout = FALSE
                     )#end list
n              = n + 1
compmodel[[n]] = list( vnam   = "albedo"
                     , desc   = "Albedo"
                     , unit   = "--"
                     , plotsd = TRUE
                     , colour = c(orange.fg,grey.fg)
                     , errcol = c(orange.bg,grey.bg)
                     , angle  = c(         45,     -45)
                     , dens   = c(         40,      40)
                     , lwd    = c(        2.5,     2.5)
                     , shwd   = c(        1.0,     1.0)
                     , type   = "o"
                     , plog   = ""
                     , legpos = "topleft"
                     , mmean  = TRUE
                     , qmean  = TRUE
                     , emean  = TRUE
                     , scsout = FALSE
                     )#end list
n              = n + 1
compmodel[[n]] = list( vnam   = "rshortup"
                     , desc   = "Outgoing SW Radiation"
                     , unit   = "W/m2"
                     , plotsd = TRUE
                     , colour = c(blue.fg,grey.fg)
                     , errcol = c(blue.bg,grey.bg)
                     , angle  = c(          45,     -45)
                     , dens   = c(          40,      40)
                     , lwd    = c(         2.5,     2.5)
                     , shwd   = c(         1.0,     1.0)
                     , type   = "o"
                     , plog   = ""
                     , legpos = "topleft"
                     , mmean  = TRUE
                     , qmean  = TRUE
                     , emean  = TRUE
                     , scsout = TRUE
                     )#end list
n              = n + 1
compmodel[[n]] = list( vnam   = "rlongup"
                     , desc   = "Outgoing LW Radiation"
                     , unit   = "W/m2"
                     , plotsd = TRUE
                     , colour = c(red.fg,grey.fg)
                     , errcol = c(red.bg,grey.bg)
                     , angle  = c(          45,     -45)
                     , dens   = c(          40,      40)
                     , lwd    = c(         2.5,     2.5)
                     , shwd   = c(         1.0,     1.0)
                     , type   = "o"
                     , plog   = ""
                     , legpos = "topleft"
                     , mmean  = TRUE
                     , qmean  = TRUE
                     , emean  = TRUE
                     , scsout = TRUE
                     )#end list
n              = n + 1
compmodel[[n]] = list( vnam   = "parup"
                     , desc   = "Outgoing PAR"
                     , unit   = "umol/m2/s"
                     , plotsd = TRUE
                     , colour = c(green.fg,grey.fg)
                     , errcol = c(green.bg,grey.bg)
                     , angle  = c(               45,     -45)
                     , dens   = c(               40,      40)
                     , lwd    = c(              2.5,     2.5)
                     , shwd   = c(              1.0,     1.0)
                     , type   = "o"
                     , plog   = ""
                     , legpos = "topleft"
                     , mmean  = TRUE
                     , qmean  = TRUE
                     , emean  = TRUE
                     , scsout = TRUE
                     )#end list
n              = n + 1
compmodel[[n]] = list( vnam   = "rnet"
                     , desc   = "Net radiation"
                     , unit   = "W/m2"
                     , plotsd = TRUE
                     , colour = c(yellow.bg,grey.fg)
                     , errcol = c(yellow.bg,grey.bg)
                     , angle  = c(           45,     -45)
                     , dens   = c(           40,      40)
                     , lwd    = c(          2.5,     2.5)
                     , shwd   = c(          1.0,     1.0)
                     , type   = "o"
                     , plog   = ""
                     , legpos = "topleft"
                     , mmean  = TRUE
                     , qmean  = TRUE
                     , emean  = TRUE
                     , scsout = TRUE
                     )#end list
#------------------------------------------------------------------------------------------#



#----- Annual mean. -----------------------------------------------------------------------#
n             = 0
soilplot      = list()
n             = n + 1
soilplot[[n]] = list( vnam   = "soil.water"
                    , desc   = "Soil moisture"
                    , unit   = "m3H2O/m3"
                    , csch   = "imuitas"
                    , pnlog  = FALSE
                    , mmean  = TRUE
                    , emean  = TRUE
                    , ymean  = TRUE
                    , scsout = TRUE
                    )#end list
n             = n + 1
soilplot[[n]] = list( vnam   = "soil.temp"
                    , desc   = "Soil temperature"
                    , unit   = "C"
                    , csch   = "muitas"
                    , pnlog  = FALSE
                    , mmean  = TRUE
                    , emean  = TRUE
                    , ymean  = TRUE
                    , scsout = TRUE
                    )#end list
n             = n + 1
soilplot[[n]] = list( vnam   = "soil.mstpot"
                    , desc   = "(Negative) Soil moisture potential"
                    , unit   = "m"
                    , csch   = "muitas"
                    , pnlog  = TRUE
                    , mmean  = TRUE
                    , emean  = TRUE
                    , ymean  = TRUE
                    , scsout = TRUE
                    )#end list
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#    List of variables to be displayed by patch.                                           #
#------------------------------------------------------------------------------------------#
n              = 0
plotpatch      = list()
n              = n + 1
plotpatch[[n]] = list( vnam       = "nep"
                     , desc       = "Net ecosystem production"
                     , unit       = "kgC/m2/yr"
                     , col.scheme = "muitas"
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "het.resp"
                     , desc       = "Heterotrophic respiration"
                     , unit       = "kgC/m2/yr"
                     , col.scheme = "muitas"
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "can.temp"
                     , desc       = "CAS temperature"
                     , unit       = "degC"
                     , col.scheme = "muitas"
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "gnd.temp"
                     , desc       = "Ground temperature"
                     , unit       = "degC"
                     , col.scheme = "muitas"
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "can.shv"
                     , desc       = "CAS specific humidity"
                     , unit       = "g/kg"
                     , col.scheme = "muitas"
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "gnd.shv"
                     , desc       = "Ground specific humidity"
                     , unit       = "g/kg"
                     , col.scheme = "muitas"
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "can.vpd"
                     , desc       = "CAS vapour pressure deficit"
                     , unit       = "hPa"
                     , col.scheme = "muitas"
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "can.co2"
                     , desc       = "CAS CO2 mixing ratio"
                     , unit       = "umol/mol"
                     , col.scheme = "muitas"
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "can.prss"
                     , desc       = "CAS pressure"
                     , unit       = "hPa"
                     , col.scheme = "muitas"
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "cflxca"
                     , desc       = "CO2 flux (CAS->Atm)"
                     , unit       = "umol/m2/s"
                     , col.scheme = "muitas"
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "cflxst"
                     , desc       = "CO2 storage at CAS"
                     , unit       = "umol/m2/s"
                     , col.scheme = "muitas"
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "nee"
                     , desc       = "Net ecosystem exchange"
                     , unit       = "umol/m2/s"
                     , col.scheme = "muitas"
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "qwflxca"
                     , desc       = "'Latent' heat flux (CAS->Atm)"
                     , unit       = "W/m2"
                     , col.scheme = "muitas"
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "hflxca"
                     , desc       = "Sensible heat flux (CAS->Atm)"
                     , unit       = "W/m2"
                     , col.scheme = "muitas"
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "hflxgc"
                     , desc       = "Sensible heat flux (Grnd->CAS)"
                     , unit       = "W/m2"
                     , col.scheme = "muitas"
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "hflxlc"
                     , desc       = "Sensible heat flux (Leaf->CAS)"
                     , unit       = "W/m2"
                     , col.scheme = "muitas"
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "hflxwc"
                     , desc       = "Sensible heat flux (Wood->CAS)"
                     , unit       = "W/m2"
                     , col.scheme = "muitas"
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "wflxca"
                     , desc       = "Water flux (CAS->Atm)"
                     , unit       = "kg/m2/day"
                     , col.scheme = "muitas"
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "wflxgc"
                     , desc       = "Water flux (Grnd->CAS)"
                     , unit       = "kg/m2/day"
                     , col.scheme = "muitas"
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "wflxlc"
                     , desc       = "Water flux (Leaf->CAS)"
                     , unit       = "kg/m2/day"
                     , col.scheme = "muitas"
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "wflxwc"
                     , desc       = "Water flux (Wood->CAS)"
                     , unit       = "kg/m2/day"
                     , col.scheme = "muitas"
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "transp"
                     , desc       = "Leaf transpiration"
                     , unit       = "kg/m2/day"
                     , col.scheme = "muitas"
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "ustar"
                     , desc       = "Friction velocity"
                     , unit       = "m/s"
                     , col.scheme = "muitas"
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "rshortup"
                     , desc       = "Outgoing SW radiation"
                     , unit       = "W/m2"
                     , col.scheme = "muitas"
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "rlongup"
                     , desc       = "Outgoing LW radiation"
                     , unit       = "W/m2"
                     , col.scheme = "muitas"
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "parup"
                     , desc       = "Outgoing PAR"
                     , unit       = "umol/m2/s"
                     , col.scheme = "muitas"
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "rnet"
                     , desc       = "Net radiation at ToC"
                     , unit       = "W/m2"
                     , col.scheme = "muitas"
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "lai"
                     , desc       = "Leaf area index"
                     , unit       = "m2/m2"
                     , col.scheme = "muitas"
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "leaf.temp"
                     , desc       = "Leaf temperature"
                     , unit       = "degC"
                     , col.scheme = "muitas"
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "leaf.vpd"
                     , desc       = "Leaf vapour pressure deficit"
                     , unit       = "hPa"
                     , col.scheme = "muitas"
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "wood.temp"
                     , desc       = "Wood temperature"
                     , unit       = "degC"
                     , col.scheme = "muitas"
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "gpp"
                     , desc       = "Gross primary productivity"
                     , unit       = "kgC/m2/yr"
                     , col.scheme = "muitas"
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "npp"
                     , desc       = "Net primary productivity"
                     , unit       = "kgC/m2/yr"
                     , col.scheme = "muitas"
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "plant.resp"
                     , desc       = "Plant respiration"
                     , unit       = "kgC/m2/yr"
                     , col.scheme = "muitas"
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "reco"
                     , desc       = "Ecosystem respiration"
                     , unit       = "kgC/m2/yr"
                     , col.scheme = "muitas"
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#    List for variables to be compared by year.                                            #
#------------------------------------------------------------------------------------------#
n            = 0
yc.tvar      = list()
n            = n + 1
yc.tvar[[n]] = list( vnam      = "rshort"
                   , desc      = "Incident SW"
                   , unit      = "W/m2"
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "fast.soil.c"
                   , desc      = "Fast soil carbon"
                   , unit      = "kgC/m2"
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "struct.soil.c"
                   , desc      = "Structural soil carbon"
                   , unit      = "kgC/m2"
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "nep"
                   , desc      = "Net Ecosystem Production"
                   , unit      = "kgC/m2"
                   , plt       = TRUE
                   , cumul     = TRUE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "nee"
                   , desc      = "Net Ecosystem Exchange"
                   , unit      = "umol/m2"
                   , plt       = TRUE
                   , cumul     = TRUE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "mco"
                   , desc      = "Maintenance Costs"
                   , unit      = "kgC/m2/yr"
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "npp"
                   , desc      = "Net primary productivity"
                   , unit      = "kgC/m2/yr"
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "cba"
                   , desc      = "Carbon balance"
                   , unit      = "kgC/m2"
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "ldrop"
                   , desc      = "Leaf drop"
                   , unit      = "kgC/m2/yr"
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "cflxca"
                   , desc      = "Carbon flux"
                   , unit      = "umol/m2/s"
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "cflxst"
                   , desc      = "Carbon storage"
                   , unit      = "umol/m2/s"
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "evap"
                   , desc      = "Evaporation"
                   , unit      = "kg/m2/day"
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "transp"
                   , desc      = "Transpiration"
                   , unit      = "kg/m2/day"
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "ustar"
                   , desc      = "Friction velocity"
                   , unit      = "m/s"
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "hflxca"
                   , desc      = "Sensible heat flux"
                   , unit      = "W/m2"
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "wflxca"
                   , desc      = "Water vapour flux"
                   , unit      = "kg/m2/day"
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "fs.open"
                   , desc      = "Water Stress index"
                   , unit      = "--"
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "gpp"
                   , desc      = "Gross Primary Productivity"
                   , unit      = "kgC/m2/yr"
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "reco"
                   , desc      = "Ecosystem respiration"
                   , unit      = "kgC/m2/yr"
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "rshortup"
                   , desc      = "Upward SW radiation"
                   , unit      = "W/m2"
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "parup"
                   , desc      = "Upward PAR"
                   , unit      = "umol/m2/s"
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "rnet"
                   , desc      = "Net radiation at the top of the canopy"
                   , unit      = "W/m2"
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "rlong"
                   , desc      = "Incoming longwave radiation"
                   , unit      = "W/m2"
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "rlongup"
                   , desc      = "Outgoing longwave radiation"
                   , unit      = "W/m2"
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "cbamax"
                   , desc      = "Maximum carbon balance"
                   , unit      = "kgC/m2"
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "cbarel"
                   , desc      = "Relative carbon balance"
                   , unit      = "--"
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "f.gpp"
                   , desc      = "Gross Primary Productivity"
                   , unit      = "%biomass/yr"
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "f.npp"
                   , desc      = "Net Primary Productivity"
                   , unit      = "%biomass/yr"
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "f.cba"
                   , desc      = "Carbon balance"
                   , unit      = "%biomass/yr"
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "f.bstorage"
                   , desc      = "Relative storage biomass"
                   , unit      = "gC_stor/kgC_bio"
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "f.bseeds"
                   , desc      = "Relative seed biomass"
                   , unit      = "gC_seed/kgC_bio"
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "f.bleaf"
                   , desc      = "Relative leaf biomass"
                   , unit      = "gC_leaf/kgC_bio"
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "agb"
                   , desc      = "Above-ground biomass"
                   , unit      = "kgC/m2"
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
#------------------------------------------------------------------------------------------#




#----- XYZ plots, to explore the parameter space. -----------------------------------------#
yc.xyzvar      = list()
yc.xyzvar$zvar = list( list ( vname = "gpp"
                            , desc  = "Gross Primary Productivity"
                            , key   = "GPP"
                            , unit  = "kgC/m2/yr"
                            )#end list
                     , list ( vname = "reco"
                            , desc  = "Ecosystem respiration"
                            , key   = "RE"
                            , unit  = "kgC/m2/yr"
                            )#end list
                     , list ( vname = "plant.resp"
                            , desc  = "Plant respiration"
                            , key   = "PR"
                            , unit  = "kgC/m2/yr"
                            )#end list
                     , list ( vname = "het.resp"
                            , desc  = "Heterotrophic respiration"
                            , key   = "RH"
                            , unit  = "kgC/m2/yr"
                            )#end list
                     , list ( vname = "nep"
                            , desc  = "Net Ecosystem Productivity"
                            , key   = "NEP"
                            , unit  = "kgC/m2/yr"
                            )#end list
                     , list ( vname = "cba"
                            , desc  = "Carbon balance"
                            , key   = "CBA"
                            , unit  = "kgC/m2/yr"
                            )#end list
                     )#end list
yc.xyzvar$xvar = list( list ( vname = "rain"
                            , desc  = "Rainfall"
                            , unit  = "mm/month"
                            , leg   = "right"
                            )#end list
                     , list ( vname = "demand"
                            , desc  = "Water demand"
                            , unit  = "kg/m2/day"
                            , leg   = "left"
                            )#end list
                     , list ( vname = "fs.open"
                            , desc  = "Soil water stress"
                            , unit  = "--"
                            , leg   = "left"
                            )#end list
                     )#end list
yc.xyzvar$yvar = list( list ( vname = "rshort"
                            , desc  = "Shortwave radiation"
                            , unit  = "W/m2"
                            , leg   = "top"
                            )#end list
                     , list ( vname = "leaf.temp"
                            , desc  = "Leaf temperature"
                            , unit  = "degC"
                            , leg   = "top"
                            )#end list
                     )#end list
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Define the number of plots of each kind, and make the lists global.                  #
#------------------------------------------------------------------------------------------#
tspftdbh    <<- tspftdbh
tserdist    <<- TRUE                  # Time series of disturbance rates
tslu        <<- tslu
squeeze     <<- squeeze
theme       <<- theme
compmodel   <<- compmodel
soilplot    <<- soilplot
plotpatch   <<- plotpatch
yc.tvar     <<- yc.tvar
yc.xyzvar   <<- yc.xyzvar
ntspftdbh   <<- length(tspftdbh      )
ntslu       <<- length(tslu          )
nsqueeze    <<- length(squeeze       )
ntheme      <<- length(theme         )
ncompmodel  <<- length(compmodel     )
nsoilplot   <<- length(soilplot      )
nplotpatch  <<- length(plotpatch     )
nyc.tvar    <<- length(yc.tvar       )
nyc.xvar    <<- length(yc.xyzvar$xvar)
nyc.yvar    <<- length(yc.xyzvar$yvar)
nyc.zvar    <<- length(yc.xyzvar$zvar)
#------------------------------------------------------------------------------------------#
