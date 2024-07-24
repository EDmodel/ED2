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
                        , e.unit   = untab$kgcom2
                        , i.unit   = untab$kgcopl
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = TRUE
                        , bar.plot = TRUE
                        , stack    = TRUE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "biomass"
                        , desc     = "Total biomass"
                        , e.unit   = untab$kgcom2
                        , i.unit   = untab$kgcopl
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
                        , e.unit   = untab$cm2om2
                        , i.unit   = untab$cm2opl
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
                        , e.unit   = untab$m2lom2
                        , i.unit   = untab$m2lom2
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
                        , e.unit   = untab$m2wom2
                        , i.unit   = untab$m2wom2
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
                        , e.unit   = untab$m2om2
                        , i.unit   = untab$m2om2
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
                        , e.unit   = untab$kgcom2oyr
                        , i.unit   = untab$kgcom2opl
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
                        , e.unit   = untab$kgcom2oyr
                        , i.unit   = untab$kgcoployr
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
                        , e.unit   = untab$kgcom2oyr
                        , i.unit   = untab$kgcoployr
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
                        , e.unit   = untab$kgcom2oyr
                        , i.unit   = untab$kgcoployr
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
                        , e.unit   = untab$kgcom2oyr
                        , i.unit   = untab$kgcoployr
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
                        , e.unit   = untab$kgcom2
                        , i.unit   = untab$kgcopl
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
                        , e.unit   = untab$kgcom2
                        , i.unit   = untab$kgcopl
                        , plog     = FALSE
                        , pft      = FALSE
                        , pftdbh   = FALSE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = TRUE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "btimber"
                        , desc     = "(Commercial) timber biomass"
                        , e.unit   = untab$kgcom2
                        , i.unit   = untab$kgcopl
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = TRUE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "bleaf"
                        , desc     = "Leaf biomass"
                        , e.unit   = untab$kgcom2
                        , i.unit   = untab$kgcopl
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
                        , e.unit   = untab$kgcom2
                        , i.unit   = untab$kgcopl
                        , plog     = FALSE
                        , pft      = FALSE
                        , pftdbh   = FALSE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = TRUE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "bstem"
                        , desc     = "Stem biomass"
                        , e.unit   = untab$kgcom2
                        , i.unit   = untab$kgcopl
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
                        , e.unit   = untab$kgcom2
                        , i.unit   = untab$kgcopl
                        , plog     = FALSE
                        , pft      = FALSE
                        , pftdbh   = FALSE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = TRUE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "bbark"
                        , desc     = "Bark biomass"
                        , e.unit   = untab$kgcom2
                        , i.unit   = untab$kgcopl
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = TRUE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "bstorage"
                        , desc     = "Storage biomass"
                        , e.unit   = untab$kgcom2
                        , i.unit   = untab$kgcopl
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
                        , e.unit   = untab$kgcom2
                        , i.unit   = untab$kgcopl
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = TRUE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "byield"
                        , desc     = "Yield (seed) biomass"
                        , e.unit   = untab$kgcom2
                        , i.unit   = untab$kgcopl
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = TRUE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "thbark"
                        , desc     = "Bark thickness"
                        , e.unit   = untab$cm
                        , i.unit   = untab$cm
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = TRUE
                        , bar.plot = FALSE
                        , stack    = TRUE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "leaf.resp"
                        , desc     = "Leaf respiration"
                        , e.unit   = untab$kgcom2oyr
                        , i.unit   = untab$kgcoployr
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = TRUE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "stem.resp"
                        , desc     = "Stem respiration"
                        , e.unit   = untab$kgcom2oyr
                        , i.unit   = untab$kgcoployr
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
                        , e.unit   = untab$kgcom2oyr
                        , i.unit   = untab$kgcoployr
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = TRUE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "aerobic.resp"
                        , desc     = "Aerobic respiration"
                        , e.unit   = untab$kgcom2oyr
                        , i.unit   = untab$kgcoployr
                        , plog     = FALSE
                        , pft      = FALSE
                        , pftdbh   = FALSE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = TRUE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "growth.resp"
                        , desc     = "Growth respiration"
                        , e.unit   = untab$kgcom2oyr
                        , i.unit   = untab$kgcoployr
                        , plog     = FALSE
                        , pft      = FALSE
                        , pftdbh   = FALSE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = TRUE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "storage.resp"
                        , desc     = "Storage respiration"
                        , e.unit   = untab$kgcom2oyr
                        , i.unit   = untab$kgcoployr
                        , plog     = FALSE
                        , pft      = FALSE
                        , pftdbh   = FALSE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = TRUE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "nplant"
                        , desc     = "Number density"
                        , e.unit   = untab$plom2
                        , i.unit   = untab$plom2
                        , plog     = TRUE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = TRUE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "sm.stress"
                        , desc     = "Soil moisture stress factor"
                        , e.unit   = untab$empty
                        , i.unit   = untab$empty
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
                        , e.unit   = untab$kgwom2loday
                        , i.unit   = untab$kgwom2loday
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "dmin.leaf.psi"
                        , desc     = "Midday leaf water potential"
                        , e.unit   = untab$mpa
                        , i.unit   = untab$mpa
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "dmax.leaf.psi"
                        , desc     = "Pre-dawn leaf water potential"
                        , e.unit   = untab$mpa
                        , i.unit   = untab$mpa
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "leaf.gbw"
                        , desc     = "Leaf boundary layer conductance"
                        , e.unit   = untab$kgwom2loday
                        , i.unit   = untab$kgwom2loday
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "census.nplant"
                        , desc     = "Number density (DBH > 10cm)"
                        , sas      = FALSE
                        , e.unit   = untab$plom2
                        , i.unit   = untab$plom2
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = TRUE
                        , scsout   = FALSE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "census.lai"
                        , desc     = "Leaf area index (H > 1.5m)"
                        , e.unit   = untab$m2lom2
                        , i.unit   = untab$m2lom2
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
                        , e.unit   = untab$kgcom2
                        , i.unit   = untab$kgcopl
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
                        , e.unit   = untab$cm2om2
                        , i.unit   = untab$cm2om2
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = TRUE
                        , scsout   = FALSE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "mort"
                        , desc     = "Mortality rate"
                        , e.unit   = untab$pcpopoyr
                        , i.unit   = untab$pcpopoyr
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = FALSE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "ncbmort"
                        , desc     = "Mortality rate - Neg. C balance"
                        , e.unit   = untab$pcpopoyr
                        , i.unit   = untab$pcpopoyr
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "hydmort"
                        , desc     = "Mortality rate - Hydr. failure"
                        , e.unit   = untab$pcpopoyr
                        , i.unit   = untab$pcpopoyr
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "dimort"
                        , desc     = "Mortality rate - Density independent"
                        , e.unit   = untab$pcpopoyr
                        , i.unit   = untab$pcpopoyr
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "fire.lethal"
                        , desc     = "Fire lethality rate"
                        , e.unit   = untab$pcpopoyr
                        , i.unit   = untab$pcpopoyr
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "recr"
                        , desc     = "Recruitment rate"
                        , e.unit   = untab$pcpopoyr
                        , i.unit   = untab$pcpopoyr
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "growth"
                        , desc     = "Growth rate"
                        , e.unit   = untab$pcdbhoyr
                        , i.unit   = untab$pcdbhoyr
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "agb.mort"
                        , desc     = "Mortality rate"
                        , e.unit   = untab$pcagboyr
                        , i.unit   = untab$pcagboyr
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = FALSE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "agb.ncbmort"
                        , desc     = "Mortality rate - Neg. C balance"
                        , e.unit   = untab$pcagboyr
                        , i.unit   = untab$pcagboyr
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "agb.hydmort"
                        , desc     = "Mortality rate - Hydr. failure"
                        , e.unit   = untab$pcagboyr
                        , i.unit   = untab$pcagboyr
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "agb.dimort"
                        , desc     = "Mortality rate - Density independent"
                        , e.unit   = untab$pcagboyr
                        , i.unit   = untab$pcagboyr
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "agb.recr"
                        , desc     = "Recruitment rate"
                        , e.unit   = untab$pcagboyr
                        , i.unit   = untab$pcagboyr
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "agb.growth"
                        , desc     = "Growth rate"
                        , e.unit   = untab$pcagboyr
                        , i.unit   = untab$pcagboyr
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "agb.change"
                        , desc     = "Change rate"
                        , e.unit   = untab$pcagboyr
                        , i.unit   = untab$pcagboyr
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = FALSE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "acc.change"
                        , desc     = "Change rate"
                        , e.unit   = untab$kgcom2oyr
                        , i.unit   = untab$kgcoployr
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = FALSE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "acc.mort"
                        , desc     = "Mortality rate"
                        , e.unit   = untab$kgcom2oyr
                        , i.unit   = untab$kgcoployr
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = FALSE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "acc.ncbmort"
                        , desc     = "Mortality rate - Neg. C balance"
                        , e.unit   = untab$kgcom2oyr
                        , i.unit   = untab$kgcoployr
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "acc.hydmort"
                        , desc     = "Mortality rate - Hydr. failure"
                        , e.unit   = untab$kgcom2oyr
                        , i.unit   = untab$kgcoployr
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "acc.dimort"
                        , desc     = "Mortality rate - Density independent"
                        , e.unit   = untab$kgcom2oyr
                        , i.unit   = untab$kgcoployr
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "acc.recr"
                        , desc     = "Recruitment rate"
                        , e.unit   = untab$kgcom2oyr
                        , i.unit   = untab$kgcoployr
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "acc.growth"
                        , desc     = "Growth rate"
                        , e.unit   = untab$kgcom2oyr
                        , i.unit   = untab$kgcoployr
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "cbamax"
                        , desc     = "Maximum Carbon balance"
                        , e.unit   = untab$kgcom2oyr
                        , i.unit   = untab$kgcoployr
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
                        , e.unit   = untab$kgcom2oyr
                        , i.unit   = untab$kgcoployr
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
                        , e.unit   = untab$kgcom2oyr
                        , i.unit   = untab$kgcoployr
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
                        , e.unit   = untab$empty
                        , i.unit   = untab$empty
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = TRUE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = FALSE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "hflxlc"
                        , desc     = "Leaf sensible heat flux"
                        , e.unit   = untab$wom2
                        , i.unit   = untab$wopl
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
                        , e.unit   = untab$kgwom2oday
                        , i.unit   = untab$kgwoploday
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
                        , e.unit   = untab$kgwom2oday
                        , i.unit   = untab$kgwoploday
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
                        , e.unit   = untab$gcokgw
                        , i.unit   = untab$gcokgw
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = FALSE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "cue"
                        , desc     = "Carbon use efficiency"
                        , e.unit   = untab$kgcokgc
                        , i.unit   = untab$kgcokgc
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = FALSE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "ecue"
                        , desc     = "Effective carbon use efficiency"
                        , e.unit   = untab$kgcokgc
                        , i.unit   = untab$kgcokgc
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = FALSE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "i.gpp"
                        , desc     = "Mean gross primary production"
                        , e.unit   = untab$kgcoployr
                        , i.unit   = untab$kgcoployr
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
                        , e.unit   = untab$kgcoployr
                        , i.unit   = untab$kgcoployr
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
                        , e.unit   = untab$kgcoployr
                        , i.unit   = untab$kgcoployr
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
                        , e.unit   = untab$kgcoployr
                        , i.unit   = untab$kgcoployr
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
                        , e.unit   = untab$kgcoployr
                        , i.unit   = untab$kgcoployr
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
                        , e.unit   = untab$kgcoployr
                        , i.unit   = untab$kgcoployr
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
                        , e.unit   = untab$kgcoployr
                        , i.unit   = untab$kgcoployr
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
                        , e.unit   = untab$kgcoployr
                        , i.unit   = untab$kgcoployr
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
                        , e.unit   = untab$wopl
                        , i.unit   = untab$wopl
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
                        , e.unit   = untab$kgwoploday
                        , i.unit   = untab$kgwoploday
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
                        , e.unit   = untab$kgwoploday
                        , i.unit   = untab$kgwoploday
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
                        , e.unit   = untab$m
                        , i.unit   = untab$m
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
                        , e.unit   = untab$pcbiooyr
                        , i.unit   = untab$pcbiooyr
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
                        , e.unit   = untab$pcbiooyr
                        , i.unit   = untab$pcbiooyr
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
                        , e.unit   = untab$pcbiooyr
                        , i.unit   = untab$pcbiooyr
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
                        , e.unit   = untab$pcbiooyr
                        , i.unit   = untab$pcbiooyr
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
                        , e.unit   = untab$gcokgcbio
                        , i.unit   = untab$gcokgcbio
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
                        , e.unit   = untab$gcokgcbio
                        , i.unit   = untab$gcokgcbio
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "f.bstem"
                        , desc     = "Relative stem biomass"
                        , e.unit   = untab$gcokgcbio
                        , i.unit   = untab$gcokgcbio
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
                        , e.unit   = untab$gcokgcbio
                        , i.unit   = untab$gcokgcbio
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
                        , e.unit   = untab$gcokgcbio
                        , i.unit   = untab$gcokgcbio
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "f.bbark"
                        , desc     = "Relative bark biomass"
                        , e.unit   = untab$gcokgcbio
                        , i.unit   = untab$gcokgcbio
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
                        , desc     = "Norm. Absorbed PAR - Leaf"
                        , e.unit   = untab$umolom2los
                        , i.unit   = untab$umolom2los
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "par.leaf"
                        , desc     = "Absolute Absorbed PAR - Leaf"
                        , e.unit   = untab$umolom2os
                        , i.unit   = untab$umolom2los
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "leaf.gpp"
                        , desc     = "Leaf-level GPP"
                        , e.unit   = untab$kgcom2loyr
                        , i.unit   = untab$kgcom2loyr
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
                        , desc     = "Norm. Absorbed SW - Leaf"
                        , e.unit   = untab$wom2l
                        , i.unit   = untab$wom2l
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
                        , desc     = "Norm. Net absorbed LW - Leaf"
                        , e.unit   = untab$wom2l
                        , i.unit   = untab$wom2l
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "assim.light"
                        , desc     = "Light-limited assimilation rate"
                        , e.unit   = untab$umolom2los
                        , i.unit   = untab$umolom2los
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = TRUE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "assim.rubp"
                        , desc     = "RuBP-limited assimilation rate"
                        , e.unit   = untab$umolom2los
                        , i.unit   = untab$umolom2los
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = TRUE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "assim.co2"
                        , desc     = "CO2-limited assimilation rate"
                        , e.unit   = untab$umolom2los
                        , i.unit   = untab$umolom2los
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = TRUE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "assim.ratio"
                        , desc     = "Light saturation"
                        , e.unit   = untab$empty
                        , i.unit   = untab$empty
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = TRUE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "vm0"
                        , desc     = "Maximum carboxylation rate (15C)"
                        , e.unit   = untab$umolom2los
                        , i.unit   = untab$umolom2los
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "llspan"
                        , desc     = "Leaf longevity"
                        , e.unit   = untab$mo
                        , i.unit   = untab$mo
                        , plog     = FALSE
                        , pft      = TRUE
                        , pftdbh   = TRUE
                        , sas      = FALSE
                        , bar.plot = FALSE
                        , stack    = FALSE
                        , scsout   = TRUE
                        )#end list
n                 = n + 1
tspftdbh[[n]]     = list( vnam     = "sla"
                        , desc     = "Specific leaf area"
                        , e.unit   = untab$m2lokgc
                        , i.unit   = untab$m2lokgc
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
                , unit = untab$kgcom2
                , plog = FALSE
                , plt  = TRUE)
n         = n + 1
tslu[[n]] = list( vnam = "biomass"
                , desc = "Total biomass"
                , unit = untab$kgcom2
                , plog = FALSE
                , plt  = TRUE)
n         = n + 1
tslu[[n]] = list( vnam = "btimber"
                , desc = "(Commercial) timber biomass"
                , unit = untab$kgcom2
                , plog = FALSE
                , plt  = TRUE)
n         = n + 1
tslu[[n]] = list( vnam = "byield"
                , desc = "Yield (seed) biomass"
                , unit = untab$kgcom2
                , plog = FALSE
                , plt  = TRUE)
n         = n + 1
tslu[[n]] = list( vnam = "thbark"
                , desc = "Bark thickness"
                , unit = untab$cm
                , plog = FALSE
                , plt  = TRUE)
n         = n + 1
tslu[[n]] = list( vnam = "lai"
                , desc = "Leaf area index"
                , unit = untab$m2lom2
                , plog = FALSE
                , plt  = TRUE)
n         = n + 1
tslu[[n]] = list( vnam = "gpp"
                , desc = "Gross primary productivity"
                , unit = untab$kgcom2oyr
                , plog = FALSE
                , plt  = TRUE)
n         = n + 1
tslu[[n]] = list( vnam = "npp"
                , desc = "Net primary productivity"
                , unit = untab$kgcom2oyr
                , plog = FALSE
                , plt  = TRUE)
n         = n + 1
tslu[[n]] = list( vnam = "area"
                , desc = "Fraction of area"
                , unit = untab$empty
                , plog = FALSE
                , plt  = TRUE)
n         = n + 1
tslu[[n]] = list( vnam = "ba"
                , desc = "Basal area"
                , unit = untab$cm2om2
                , plog = FALSE
                , plt  = TRUE)
n         = n + 1
tslu[[n]] = list( vnam = "f.agb"
                , desc = "Mean above ground biomass"
                , unit = untab$kgcom2
                , plog = FALSE
                , plt  = TRUE)
n         = n + 1
tslu[[n]] = list( vnam = "f.biomass"
                , desc = "Mean total biomass"
                , unit = untab$kgcom2
                , plog = FALSE
                , plt  = TRUE)
n         = n + 1
tslu[[n]] = list( vnam = "f.lai"
                , desc = "Mean leaf area index"
                , unit = untab$m2lom2
                , plog = FALSE
                , plt  = TRUE)
n         = n + 1
tslu[[n]] = list( vnam = "f.gpp"
                , desc = "Mean gross primary productivity"
                , unit = untab$kgcom2oyr
                , plog = FALSE
                , plt  = TRUE)
n         = n + 1
tslu[[n]] = list( vnam = "f.npp"
                , desc = "Mean net primary productivity"
                , unit = untab$kgcom2oyr
                , plog = FALSE
                , plt  = TRUE)
n         = n + 1
tslu[[n]] = list( vnam = "f.ba"
                , desc = "Mean basal area"
                , unit = untab$cm2om2
                , plog = FALSE
                , plt  = TRUE)
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
                   , unit       = untab$kgcom2oyr
                   , col.scheme = "brbg"
                   , fco.mmean  = TRUE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "npp"
                   , desc       = "Net Primary productivity"
                   , unit       = untab$kgcom2oyr
                   , col.scheme = "brbg"
                   , fco.mmean  = TRUE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "plant.resp"
                   , desc       = "Plant respiration"
                   , unit       = untab$kgcom2oyr
                   , col.scheme = "irdbu"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "het.resp"
                   , desc       = "Heterotrophic respiration"
                   , unit       = untab$kgcom2oyr
                   , col.scheme = "irdbu"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "fgc.resp"
                   , desc       = "Surface litter respiration"
                   , unit       = untab$kgcom2oyr
                   , col.scheme = "irdbu"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "fsc.resp"
                   , desc       = "Sub-surface litter respiration"
                   , unit       = untab$kgcom2oyr
                   , col.scheme = "irdbu"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "stgc.resp"
                   , desc       = "Surface woody debris respiration"
                   , unit       = untab$kgcom2oyr
                   , col.scheme = "irdbu"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "stsc.resp"
                   , desc       = "Sub-surface woody debris respiration"
                   , unit       = untab$kgcom2oyr
                   , col.scheme = "irdbu"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "msc.resp"
                   , desc       = "Microbial soil respiration"
                   , unit       = untab$kgcom2oyr
                   , col.scheme = "irdbu"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "ssc.resp"
                   , desc       = "Humified soil respiration"
                   , unit       = untab$kgcom2oyr
                   , col.scheme = "irdbu"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "psc.resp"
                   , desc       = "Passive soil respiration"
                   , unit       = untab$kgcom2oyr
                   , col.scheme = "irdbu"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "fgc.in"
                   , desc       = "Surface litter input"
                   , unit       = untab$kgcom2oyr
                   , col.scheme = "magma"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "fsc.in"
                   , desc       = "Sub-surface litter input"
                   , unit       = untab$kgcom2oyr
                   , col.scheme = "magma"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "stgc.in"
                   , desc       = "Surface woody debris input"
                   , unit       = untab$kgcom2oyr
                   , col.scheme = "magma"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "stsc.in"
                   , desc       = "Surface woody debris input"
                   , unit       = untab$kgcom2oyr
                   , col.scheme = "magma"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "nep"
                   , desc       = "Net ecosystem production"
                   , unit       = untab$kgcom2oyr
                   , col.scheme = "irdbu"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "reco"
                   , desc       = "Ecosystem respiration"
                   , unit       = untab$kgcom2oyr
                   , col.scheme = "irdbu"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "hflxca"
                   , desc       = "Sensible heat flux"
                   , unit       = untab$wom2
                   , col.scheme = "irdbu"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "wflxca"
                   , desc       = "Water flux"
                   , unit       = untab$kgwom2oday
                   , col.scheme = "rdbu"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "wflxgc"
                   , desc       = "Ground evaporation"
                   , unit       = untab$kgwom2oday
                   , col.scheme = "rdbu"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "wflxlc"
                   , desc       = "Leaf evaporation"
                   , unit       = untab$kgwom2oday
                   , col.scheme = "rdbu"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "wflxwc"
                   , desc       = "Wood evaporation"
                   , unit       = untab$kgwom2oday
                   , col.scheme = "rdbu"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "evap"
                   , desc       = "Evaporation"
                   , unit       = untab$kgwom2oday
                   , col.scheme = "rdbu"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "transp"
                   , desc       = "Transpiration"
                   , unit       = untab$kgwom2oday
                   , col.scheme = "rdbu"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "nee"
                   , desc       = "Net ecosystem exchange"
                   , unit       = untab$umolcom2os
                   , col.scheme = "rdbu"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "cflxca"
                   , desc       = "CO2 flux"
                   , unit       = untab$umolcom2os
                   , col.scheme = "rdbu"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "cflxst"
                   , desc       = "CO2 flux"
                   , unit       = untab$umolcom2os
                   , col.scheme = "rdbu"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "can.temp"
                   , desc       = "Canopy air temperature"
                   , unit       = untab$degC
                   , col.scheme = "puor"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "atm.temp"
                   , desc       = "Atmospheric temperature"
                   , unit       = untab$degC
                   , col.scheme = "puor"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "leaf.temp"
                   , desc       = "Leaf temperature"
                   , unit       = untab$degC
                   , col.scheme = "puor"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "wood.temp"
                   , desc       = "Wood temperature"
                   , unit       = untab$degC
                   , col.scheme = "puor"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "gnd.temp"
                   , desc       = "Ground temperature"
                   , unit       = untab$degC
                   , col.scheme = "puor"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "atm.shv"
                   , desc       = "Atmospheric specific humidity" 
                   , unit       = untab$gwokg
                   , col.scheme = "rdbu"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "can.shv"
                   , desc       = "Canopy air specific humidity"
                   , unit       = untab$gwokg
                   , col.scheme = "rdbu"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "gnd.shv"
                   , desc       = "Ground specific humidity"
                   , unit       = untab$gwokg
                   , col.scheme = "rdbu"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "atm.co2"
                   , desc       = "Atmospheric CO2 mixing ratio"
                   , unit       = untab$molcomol
                   , col.scheme = "puor"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "can.co2"
                   , desc       = "Canopy air CO2 mixing ratio"
                   , unit       = untab$molcomol
                   , col.scheme = "puor"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "rain"
                   , desc       = "Total monthly precipitation"
                   , unit       = untab$mm
                   , col.scheme = "rdbu"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "runoff"
                   , desc       = "Total monthly runoff"
                   , unit       = untab$mm
                   , col.scheme = "rdbu"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "intercepted"
                   , desc       = "Total monthly interception"
                   , unit       = untab$mm
                   , col.scheme = "rdbu"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "wshed"
                   , desc       = "Total monthly dripping"
                   , unit       = untab$mm
                   , col.scheme = "rdbu"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "sm.stress"
                   , desc       = "Fraction of open stomata"
                   , unit       = untab$empty
                   , col.scheme = "rdbu"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "leaf.gbw"
                   , desc       = "Leaf boundary layer conductance"
                   , unit       = untab$kgwom2loday
                   , col.scheme = "rdbu"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "leaf.gsw"
                   , desc       = "Stomatal conductance"
                   , unit       = untab$kgwom2loday
                   , col.scheme = "rdbu"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "rshort"
                   , desc       = "Downward shortwave radiation"
                   , unit       = untab$wom2
                   , col.scheme = "ibugy"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "rshort.gnd"
                   , desc       = "Abs. gnd. shortwave radiation"
                   , unit       = untab$wom2
                   , col.scheme = "ibugy"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "rshortup"
                   , desc       = "Outgoing shortwave radiation"
                   , unit       = untab$wom2
                   , col.scheme = "ibugy"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "rlong"
                   , desc       = "Downward longwave radiation"
                   , unit       = untab$wom2
                   , col.scheme = "bugy"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "rlong.gnd"
                   , desc       = "Abs. gnd. longwave radiation"
                   , unit       = untab$wom2
                   , col.scheme = "bugy"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "rlongup"
                   , desc       = "Outgoing longwave radiation"
                   , unit       = untab$wom2
                   , col.scheme = "bugy"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "albedo"
                   , desc       = "Shortwave albedo"
                   , unit       = untab$empty
                   , col.scheme = "irdbu"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "ustar"
                   , desc       = "Friction velocity"
                   , unit       = untab$mos
                   , col.scheme = "irdbu"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "paw"
                   , desc       = "Potential available water"
                   , unit       = untab$pcsat
                   , col.scheme = "rdbu"
                   , fco.mmean  = FALSE
                   , fco.qmean  = FALSE
                   , box.plot   = FALSE
                   )#end list
n            = n + 1
squeeze[[n]] = list( vnam       = "smpot"
                   , desc       = "Integrated matric potential"
                   , unit       = untab$mpa
                   , col.scheme = "irdbu"
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
                 , colour    = c(    "#009E73",    "#F0E442",   "#332288",    "#882255"
                                ,    "#56B4E9",    "#0072B2")
                 , lwd       = c(          2.5,          2.5,         2.5,          2.5
                                ,          2.5,          2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "ecoflux"
                 , title     = "Ecosystem fluxes"
                 , unit      = untab$kgcom2oyr
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = TRUE
                 , ymean     = TRUE
                 , stack     = FALSE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(NA,NA)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(NA,NA)
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(   "rshort",      "rlong","rshort.gnd",    "qwflxca"
                                ,   "hflxca")
                 , desc      = c(  "Down SW",    "Down LW", "Abs. Grnd",   "'Latent'"
                                , "Sensible")
                 , colour    = c(  "#E69F00",    "#56B4E9",   "#332288",    "#0072B2"
                                ,  "#882255")
                 , lwd       = c(        2.5,          2.5,         2.5,          2.5
                                        ,2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "eneflux"
                 , title     = "Energy fluxes"
                 , unit      = untab$wom2
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = TRUE
                 , ymean     = TRUE
                 , stack     = FALSE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(NA,NA)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(NA,NA)
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(        "wflxgc",         "wflxca",        "wflxlc"
                                ,        "wflxwc",         "transp")
                 , desc      = c("Ground->Canopy",    "Canopy->Atm",  "Leaf->Canopy"
                                ,  "Wood->Canopy",  "Transpiration")
                 , colour    = c(       "#332288",        "#0072B2",       "#56B4E9"
                                ,       "#882255",        "#009E73")
                 , lwd       = c(             2.5,              2.5,             2.5
                                ,             2.5,              2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "h2oflux"
                 , title     = "Water fluxes"
                 , unit      = untab$kgwom2oday
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = TRUE
                 , ymean     = TRUE
                 , stack     = FALSE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(NA,NA)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(NA,NA)
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(        "hflxgc",         "hflxca",        "hflxlc"
                                ,        "hflxwc")
                 , desc      = c("Ground->Canopy",    "Canopy->Atm",  "Leaf->Canopy"
                                ,  "Wood->Canopy")
                 , colour    = c(       "#332288",        "#0072B2",       "#56B4E9"
                                ,       "#882255",        "#009E73")
                 , lwd       = c(             2.5,              2.5,             2.5
                                ,             2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "sensflux"
                 , title     = "Sensible heat fluxes"
                 , unit      = untab$wom2
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = TRUE
                 , ymean     = TRUE
                 , stack     = FALSE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(NA,NA)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(NA,NA)
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(     "atm.temp",  "can.temp",  "leaf.temp"
                                ,    "wood.temp",  "gnd.temp")
                 , desc      = c(   "Atmosphere","Canopy air",       "Leaf"
                                ,         "Wood",    "Ground")
                 , colour    = c(      "#56B4E9",    "grey45",    "#009E73"
                                ,      "#882255",   "#332288")
                 , lwd       = c(            2.5,         2.5,          2.5
                                ,            2.5,         2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "temperature"
                 , title     = "Temperature"
                 , unit      = untab$degC
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = TRUE
                 , ymean     = TRUE
                 , stack     = FALSE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(NA,NA)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(NA,NA)
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(    "atm.shv",   "can.shv",      "gnd.shv")
                 , desc      = c( "Atmosphere","Canopy air",       "Ground")
                 , colour    = c(    "#56B4E9",   "#0072B2",      "#332288")
                 , lwd       = c(          2.5,         2.5,            2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "h2ovapour"
                 , title     = "Water vapour mixing ratio"
                 , unit      = untab$gwokg
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = TRUE
                 , ymean     = TRUE
                 , stack     = FALSE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(NA,NA)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(NA,NA)
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(    "atm.co2",    "can.co2")
                 , desc      = c( "Atmosphere", "Canopy air")
                 , colour    = c(    "#56B4E9",    "#0072B2")
                 , lwd       = c(2.5,2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "co2"
                 , title     = "CO2 mixing ratio"
                 , unit      = untab$umolcomol
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = TRUE
                 , ymean     = TRUE
                 , stack     = FALSE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(NA,NA)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(NA,NA)
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(         "rain",      "runoff", "intercepted",   "wshed")
                 , desc      = c("Precipitation","Total runoff","Interception","Dripping")
                 , colour    = c(      "#0072B2",     "#E69F00",     "#009E73", "#332288")
                 , lwd       = c(            2.5,           2.5,           2.5,       2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "prec"
                 , title     = "Precipitation distribution"
                 , unit      = untab$mm
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = TRUE
                 , ymean     = TRUE
                 , stack     = FALSE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(NA,NA)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(NA,NA)
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c("npat.global")
                 , desc      = c("Patch count")
                 , colour    = c(    "#D55E00")
                 , lwd       = c(          2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "npatch"
                 , title     = "Total number of patches"
                 , unit      = untab$empty
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = FALSE
                 , ymean     = TRUE
                 , stack     = TRUE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(NA,NA)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(NA,NA)
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c( "ncoh.global")
                 , desc      = c("Cohort count")
                 , colour    = c(     "#009E73")
                 , lwd       = c(           2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "ncohort"
                 , title     = "Total number of cohorts"
                 , unit      = untab$empty
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = FALSE
                 , ymean     = TRUE
                 , stack     = TRUE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(NA,NA)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(NA,NA)
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(         "workload",                "specwork")
                 , desc      = c("RK4 steps (Total)","Avg. RK4 steps per patch")
                 , colour    = c(          "#56B4E9",                 "#009E73")
                 , lwd       = c(                2.5,                       2.5)
                 , type      = "o"
                 , plog      = TRUE
                 , prefix    = "workload"
                 , title     = "Work Load"
                 , unit      = untab$empty
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = FALSE
                 , ymean     = TRUE
                 , stack     = FALSE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(NA,NA)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(NA,NA)
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(                "rk4step")
                 , desc      = c("Thermodynamic time step")
                 , colour    = c(                "#0072B2")
                 , lwd       = c(                      2.5)
                 , type      = "o"
                 , plog      = TRUE
                 , prefix    = "rk4step"
                 , title     = "Thermodynamic time step"
                 , unit      = untab$s
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = TRUE
                 , ymean     = TRUE
                 , stack     = FALSE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(NA,NA)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(NA,NA)
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(     "root.resp",     "stem.resp",     "leaf.resp")
                 , desc      = c(          "Root",          "Stem",          "Leaf")
                 , colour    = c(       "#332288",       "#E69F00",       "#009E73")
                 , lwd       = c(             2.5,             2.5,             2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "pltissueresp"
                 , title     = "Auto. respiration fluxes"
                 , unit      = untab$kgcom2oyr
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = FALSE
                 , ymean     = TRUE
                 , stack     = TRUE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(NA,NA)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(NA,NA)
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(  "aerobic.resp",   "growth.resp",  "storage.resp")
                 , desc      = c(       "Aerobic",        "Growth",       "Storage")
                 , colour    = c(       "#0072B2",       "#E69F00",       "#785EF0")
                 , lwd       = c(             2.5,             2.5,             2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "plprocresp"
                 , title     = "Auto. respiration fluxes"
                 , unit      = untab$kgcom2oyr
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = FALSE
                 , ymean     = TRUE
                 , stack     = TRUE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(0.,4.5)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(0.,4.5)
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(       "psc.resp",       "ssc.resp",      "stsc.resp"
                                ,       "msc.resp",       "fsc.resp",      "stgc.resp"
                                ,       "fgc.resp")
                 , desc      = c(   "Passive soil",  "Humified soil","BG Woody Debris"
                                ,      "Microbial",      "BG Litter","AG Woody Debris"
                                ,      "AG Litter")
                 , colour    = c(        "#811F9E",        "#1BA2F7",        "#880D32"
                                ,        "#CCCA3D",        "#107C92",        "#F87856"
                                ,        "#2BD2DB")
                 , lwd       = c(              2.5,              2.5,              2.5
                                ,              2.5,              2.5,              2.5
                                ,              2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "hetresp"
                 , title     = "Het. Respiration fluxes"
                 , unit      = untab$kgcom2oyr
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = FALSE
                 , ymean     = TRUE
                 , stack     = TRUE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(NA,NA)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(NA,NA)
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(        "stsc.in",         "fsc.in"
                                ,        "stgc.in",         "fgc.in")
                 , desc      = c("BG Woody Debris",      "BG Litter"
                                ,"AG Woody Debris",      "AG Litter")
                 , colour    = c(        "#880D32",        "#107C92"
                                ,        "#F87856",        "#2BD2DB")
                 , lwd       = c(              2.5,              2.5
                                ,              2.5,              2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "necro.input"
                 , title     = "Necromass inputs"
                 , unit      = untab$kgcom2oyr
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = FALSE
                 , ymean     = TRUE
                 , stack     = TRUE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(NA,NA)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(NA,NA)
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(  "atm.vels",            "ustar")
                 , desc      = c("Wind speed","Friction velocity")
                 , colour    = c(   "#56B4E9",          "#332288")
                 , lwd       = c(          2.5,               2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "wind"
                 , title     = "Winds"
                 , unit      = untab$mos
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = TRUE
                 , ymean     = TRUE
                 , stack     = FALSE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(NA,NA)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(NA,NA)
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c( "passive.soil.c",    "slow.soil.c", "microbe.soil.c"
                                ,  "struct.soil.c",    "fast.soil.c",  "struct.grnd.c"
                                ,    "fast.grnd.c")
                 , desc      = c(        "Passive",       "Humified",      "Microbial"
                                ,"BG Woody Debris",      "BG Litter","AG Woody Debris"
                                ,      "AG Litter")
                 , colour    = c(        "#811F9E",        "#1BA2F7",        "#880D32"
                                ,        "#CCCA3D",        "#107C92",        "#F87856"
                                ,        "#2BD2DB")
                 , lwd       = c(              2.5,              2.5,              2.5
                                ,              2.5,              2.5,              2.5
                                ,              2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "soil_carbon"
                 , title     = "Soil Carbon"
                 , unit      = untab$kgcom2
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = FALSE
                 , ymean     = TRUE
                 , stack     = TRUE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(NA,NA)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(NA,NA)
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(       "crop.yield",     "crop.harvest"
                                ,  "logging.harvest",   "combusted.fuel")
                 , desc      = c(     "Crop (Seeds)",     "Crop (Other)"
                                ,  "Logging Harvest","Combusted Biomass")
                 , colour    = c(          "#009E73",          "#0072B2"
                                ,          "#E69F00",          "#882255")
                 , lwd       = c(                2.5,                2.5
                                ,                2.5,                2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "output_carbon"
                 , title     = "Carbon output"
                 , unit      = untab$kgcom2
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = FALSE
                 , ymean     = TRUE
                 , stack     = TRUE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(NA,NA)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(NA,NA)
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(    "atm.vpd",    "can.vpd",   "leaf.vpd")
                 , desc      = c( "Atmosphere", "Canopy air",       "Leaf")
                 , colour    = c(    "#56B4E9",    "#0072B2",    "#009E73")
                 , lwd       = c(          2.5,          2.5,          2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "vpdef"
                 , title     = "Vapour pressure deficit"
                 , unit      = untab$hpa
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = TRUE
                 , ymean     = TRUE
                 , stack     = FALSE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(NA,NA)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(NA,NA)
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(         "paw")
                 , desc      = c("Pot.Av.Water")
                 , colour    = c(     "#0072B2")
                 , lwd       = c(           2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "paw"
                 , title     = "Average potential available water"
                 , unit      = untab$pcsat
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = FALSE
                 , ymean     = TRUE
                 , stack     = TRUE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(NA,NA)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(NA,NA)
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(         "smpot")
                 , desc      = c("Neg. Potential")
                 , colour    = c(       "#0072B2")
                 , lwd       = c(             2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "smpot"
                 , title     = "Average soil matric potential"
                 , unit      = untab$mpa
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = FALSE
                 , ymean     = TRUE
                 , stack     = FALSE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(NA,NA)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(NA,NA)
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c("water.deficit",      "malhi.deficit")
                 , desc      = c(       "ED-2.2","Malhi et al. (2009)")
                 , colour    = c(      "#882255",            "#E69F00")
                 , lwd       = c(            2.5,                  2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "water_deficit"
                 , title     = "Monthly Water deficit"
                 , unit      = untab$mmomo
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = FALSE
                 , ymean     = TRUE
                 , stack     = FALSE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(NA,NA)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(NA,NA)
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(        "nee",   "cflxca",     "cflxst")
                 , desc      = c(        "NEE", "CO2 Flux","CO2 Storage")
                 , colour    = c(    "#009E73",  "#56B4E9",    "#E69F00")
                 , lwd       = c(          2.5,        2.5,          2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "carbflux"
                 , title     = "CO2 fluxes"
                 , unit      = untab$umolcom2os
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = TRUE
                 , ymean     = TRUE
                 , stack     = FALSE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(NA,NA)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(NA,NA)
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(         "rshort",  "rshort.beam","rshort.diff"
                                ,     "rshort.gnd",     "rshortup")
                 , desc      = c("Down Top canopy",         "Beam",    "Diffuse"
                                ,    "Abs. Ground","Up Top canopy")
                 , colour    = c(        "#56B4E9",      "#E69F00",     "grey45"
                                ,        "#882255",      "#0072B2")
                 , lwd       = c(              2.5,            2.5,          2.5
                                ,              2.5,            2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "rshort"
                 , title     = "Short wave radiation"
                 , unit      = untab$wom2
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = TRUE
                 , ymean     = TRUE
                 , stack     = FALSE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(NA,NA)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(NA,NA)
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(          "rlong",  "rlongup",    "rlong.gnd")
                 , desc      = c("Down Top canopy","Upward LW",  "Abs. Ground")
                 , colour    = c(        "#56B4E9",  "#E69F00",      "#882255")
                 , lwd       = c(              2.5,        2.5,            2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "rlong"
                 , title     = "Long wave radiation"
                 , unit      = untab$wom2
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = TRUE
                 , ymean     = TRUE
                 , stack     = FALSE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(NA,NA)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(NA,NA)
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(         "albedo", "albedo.par","albedo.nir")
                 , desc      = c("SW Albedo (Net)", "PAR Albedo","NIR Albedo")
                 , colour    = c(        "#56B4E9",    "#009E73",   "#882255")
                 , lwd       = c(              2.5,          2.5,         2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "albedo"
                 , title     = "Albedo"
                 , unit      = untab$empty
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = TRUE
                 , ymean     = TRUE
                 , stack     = FALSE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(NA,NA)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(NA,NA)
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(        "par.tot",     "par.beam", "par.diff"
                                ,        "par.gnd",        "parup")
                 , desc      = c("Down Top canopy",         "Beam",  "Diffuse"
                                ,    "Abs. Ground","Up Top canopy")
                 , colour    = c(        "#56B4E9",      "#882255",  "#332288"
                                ,        "#E69F00",      "#009E73")
                 , lwd       = c(              2.5,            2.5,        2.5
                                ,              2.5,            2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "par"
                 , title     = "Photosynthetically Active Radiation"
                 , unit      = untab$umolom2os
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = TRUE
                 , ymean     = TRUE
                 , stack     = FALSE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(NA,NA)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(NA,NA)
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(       "leaf.gsw",        "leaf.gbw",        "wood.gbw")
                 , desc      = c( "Leaf (Stomata)","Leaf (Bnd. Lyr.)","Wood (Bnd. Lyr.)")
                 , colour    = c(        "#009E73",         "#56B4E9",         "#882255")
                 , lwd       = c(              2.5,               2.5,               2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "conduct"
                 , title     = "Conductance"
                 , unit      = untab$kgwom2oday
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = TRUE
                 , ymean     = TRUE
                 , stack     = FALSE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(NA,NA)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(NA,NA)
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(                "vm0")
                 , desc      = c( "Max. Carboxylation")
                 , colour    = c(            "#0072B2")
                 , lwd       = c(                  2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "vm0"
                 , title     = "Maximum carboxylation"
                 , unit      = untab$umolom2los
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = FALSE
                 , ymean     = TRUE
                 , stack     = FALSE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(NA,NA)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(NA,NA)
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(                "sla")
                 , desc      = c( "Specific leaf area")
                 , colour    = c(            "#009E73")
                 , lwd       = c(                  2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "sla"
                 , title     = "Specific leaf area"
                 , unit      = untab$m2lokgc
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = FALSE
                 , ymean     = TRUE
                 , stack     = FALSE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(NA,NA)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(NA,NA)
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(         "llspan")
                 , desc      = c( "Leaf longevity")
                 , colour    = c(        "#332288")
                 , lwd       = c(              2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "llspan"
                 , title     = "Leaf longevity"
                 , unit      = untab$mo
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = FALSE
                 , ymean     = TRUE
                 , stack     = FALSE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(NA,NA)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(NA,NA)
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(          "veg.height",           "can.depth"
                                ,        "veg.displace",           "veg.rough"
                                ,           "can.rough")
                 , desc      = c(   "Vegetation height",        "Canopy depth"
                                , "Displacement height","Vegetation Roughness"
                                ,       "Net roughness")
                 , colour    = c(             "#56B4E9",             "#882255"
                                ,             "#332288",             "#E69F00"
                                ,             "#009E73")
                 , lwd       = c(2.5,2.5,2.5,2.5,2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "can.struct"
                 , title     = "Canopy structure"
                 , unit      = untab$m
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = FALSE
                 , ymean     = TRUE
                 , stack     = FALSE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(NA,NA)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(NA,NA)
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(      "fire.intensity")
                 , desc      = c(      "Fire intensity")
                 , colour    = c(             "#882255")
                 , lwd       = c(2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "fire.intensity"
                 , title     = "Fire intensity"
                 , unit      = untab$kwom
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = FALSE
                 , ymean     = TRUE
                 , stack     = FALSE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(NA,NA)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(NA,NA)
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(      "fire.density")
                 , desc      = c("Fire count density")
                 , colour    = c(           "#882255")
                 , lwd       = c(2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "fire.density"
                 , title     = "Fire count density"
                 , unit      = untab$oneokm2
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = FALSE
                 , ymean     = TRUE
                 , stack     = FALSE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(NA,NA)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(NA,NA)
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(      "fire.ignition")
                 , desc      = c(      "Fire ignition")
                 , colour    = c(            "#882255")
                 , lwd       = c(2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "fire.ignition"
                 , title     = "Fire ignition rate"
                 , unit      = untab$oneokm2omo
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = FALSE
                 , ymean     = TRUE
                 , stack     = FALSE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(NA,NA)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(NA,NA)
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c("fire.extinction")
                 , desc      = c("Fire extinction")
                 , colour    = c(        "#0072B2")
                 , lwd       = c(2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "fire.extinction"
                 , title     = "Fire extinction rate"
                 , unit      = untab$pcoday
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = FALSE
                 , ymean     = TRUE
                 , stack     = FALSE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(NA,NA)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(NA,NA)
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(      "fire.spread")
                 , desc      = c(      "Fire spread")
                 , colour    = c(          "#332288")
                 , lwd       = c(2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "fire.spread"
                 , title     = "Fire spread rate"
                 , unit      = untab$momin
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = FALSE
                 , ymean     = TRUE
                 , stack     = FALSE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(NA,NA)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(NA,NA)
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(     "fire.f.bherb",    "fire.f.bwoody"
                                ,       "fire.f.fgc",      "fire.f.stgc")
                 , desc      = c(       "Herbaceous",     "Woody (live)"
                                ,      "Fine litter","Structural litter")
                 , colour    = c(          "#009E73",          "#E69F00"
                                ,          "#56B4E9",          "#332288")
                 , lwd       = c(2.5,2.5,2.5,2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "fire.consumption"
                 , title     = "Relative fuel consumption"
                 , unit      = untab$pc
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = FALSE
                 , ymean     = TRUE
                 , stack     = FALSE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(NA,NA)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(NA,NA)
                 )#end list
n          = n + 1
theme[[n]] = list( vnam      = c(          "burnt.area")
                 , desc      = c(          "Burnt area")
                 , colour    = c(             "#D55E00")
                 , lwd       = c(2.5)
                 , type      = "o"
                 , plog      = FALSE
                 , prefix    = "burnt.area"
                 , title     = "Burnt area"
                 , unit      = untab$pc
                 , legpos    = "topleft"
                 , emean     = TRUE
                 , mmean     = TRUE
                 , qmean     = FALSE
                 , ymean     = TRUE
                 , stack     = FALSE
                 , emean.lim = c(NA,NA)
                 , mmean.lim = c(NA,NA)
                 , qmean.lim = c(NA,NA)
                 , ymean.lim = c(NA,NA)
                 )#end list
#------------------------------------------------------------------------------------------#




#----- Comparison between observations and model averages. --------------------------------#
n              = 0
compmodel      = list()
n              = n + 1
compmodel[[n]] = list( vnam   = "nep"
                     , desc   = "Net Ecosystem Productivity"
                     , unit   = untab$kgcom2oyr
                     , plotsd = TRUE
                     , colour = c(green.mg,grey.fg)
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
                     , unit   = untab$kgcom2oyr
                     , plotsd = TRUE
                     , colour = c(green.mg,grey.fg)
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
                     , unit   = untab$umolcomol
                     , plotsd = TRUE
                     , colour = c(green.mg,grey.fg)
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
                     , unit   = untab$wom2
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
                     , unit   = untab$kgwom2oday
                     , plotsd = TRUE
                     , colour = c(sky.fg,grey.fg)
                     , errcol = c(sky.bg,grey.bg)
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
                     , unit   = untab$wom2
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
                     , unit   = untab$umolom2os
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
                     , unit   = untab$wom2
                     , plotsd = TRUE
                     , colour = c(sky.fg,grey.fg)
                     , errcol = c(sky.bg,grey.bg)
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
                     , unit   = untab$degC
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
                     , unit   = untab$gwokg
                     , plotsd = TRUE
                     , colour = c(sky.fg,grey.fg)
                     , errcol = c(sky.bg,grey.bg)
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
                     , unit   = untab$mmoday
                     , plotsd = FALSE
                     , colour = c(sky.fg,grey.fg)
                     , errcol = c(sky.bg,grey.bg)
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
                     , unit   = untab$mos
                     , plotsd = TRUE
                     , colour = c(indigo.mg,grey.fg)
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
                     , unit   = untab$hpa
                     , plotsd = TRUE
                     , colour = c(purple.mg,grey.fg)
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
                     , unit   = untab$umolcom2os
                     , plotsd = TRUE
                     , colour = c(green.mg,grey.fg)
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
                     , unit   = untab$umolcom2os
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
                     , unit   = untab$umolcom2os
                     , plotsd = TRUE
                     , colour = c(green.mg,grey.fg)
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
                     , unit   = untab$mos
                     , plotsd = TRUE
                     , colour = c(purple.mg,grey.fg)
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
                     , unit   = untab$kgcom2oyr
                     , plotsd = TRUE
                     , colour = c(orange.mg,grey.fg)
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
                     , scsout = FALSE
                     )#end list
n              = n + 1
compmodel[[n]] = list( vnam   = "albedo"
                     , desc   = "Albedo"
                     , unit   = untab$empty
                     , plotsd = TRUE
                     , colour = c(yellow.fg,grey.fg)
                     , errcol = c(yellow.bg,grey.bg)
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
                     , unit   = untab$wom2
                     , plotsd = TRUE
                     , colour = c(sky.fg,grey.fg)
                     , errcol = c(sky.bg,grey.bg)
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
                     , unit   = untab$wom2
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
                     , unit   = untab$umolom2os
                     , plotsd = TRUE
                     , colour = c(green.mg,grey.fg)
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
                     , unit   = untab$wom2
                     , plotsd = TRUE
                     , colour = c(yellow.mg,grey.fg)
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
n              = n + 1
compmodel[[n]] = list( vnam   = "soil.temp.top"
                     , desc   = "Soil temperature (0-20cm)"
                     , unit   = untab$degC
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
compmodel[[n]] = list( vnam   = "soil.water.top"
                     , desc   = "Soil water (0-50cm)"
                     , unit   = untab$kgwom2
                     , plotsd = TRUE
                     , colour = c(sky.fg,grey.fg)
                     , errcol = c(sky.bg,grey.bg)
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
compmodel[[n]] = list( vnam   = "soil.water.bot"
                     , desc   = "Soil water (50-200cm)"
                     , unit   = untab$kgwom2
                     , plotsd = TRUE
                     , colour = c(sky.fg,grey.fg)
                     , errcol = c(sky.bg,grey.bg)
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
compmodel[[n]] = list( vnam   = "soil.wetness.top"
                     , desc   = "Soil wetness (0-50cm)"
                     , unit   = untab$empty
                     , plotsd = TRUE
                     , colour = c(sky.fg,grey.fg)
                     , errcol = c(sky.bg,grey.bg)
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
compmodel[[n]] = list( vnam   = "soil.wetness.bot"
                     , desc   = "Soil wetness (50-200cm)"
                     , unit   = untab$empty
                     , plotsd = TRUE
                     , colour = c(sky.fg,grey.fg)
                     , errcol = c(sky.bg,grey.bg)
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
#------------------------------------------------------------------------------------------#



#----- Annual mean. -----------------------------------------------------------------------#
n             = 0
soilplot      = list()
n             = n + 1
soilplot[[n]] = list( vnam   = "soil.water"
                    , desc   = "Soil moisture"
                    , unit   = untab$m3wom3
                    , csch   = "rdbu"
                    , pnlog  = FALSE
                    , mmean  = TRUE
                    , emean  = TRUE
                    , ymean  = TRUE
                    , scsout = TRUE
                    )#end list
n             = n + 1
soilplot[[n]] = list( vnam   = "soil.temp"
                    , desc   = "Soil temperature"
                    , unit   = untab$degC
                    , csch   = "irdbu"
                    , pnlog  = FALSE
                    , mmean  = TRUE
                    , emean  = TRUE
                    , ymean  = TRUE
                    , scsout = TRUE
                    )#end list
n             = n + 1
soilplot[[n]] = list( vnam   = "soil.mstpot"
                    , desc   = "(Negative) Soil moisture potential"
                    , unit   = untab$mpa
                    , csch   = "irdbu"
                    , pnlog  = TRUE
                    , mmean  = TRUE
                    , emean  = TRUE
                    , ymean  = TRUE
                    , scsout = TRUE
                    )#end list
n             = n + 1
soilplot[[n]] = list( vnam   = "soil.extracted"
                    , desc   = "Water extraction by plants"
                    , unit   = untab$kgwom3oday
                    , csch   = "rdbu"
                    , pnlog  = FALSE
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
                     , unit       = untab$kgcom2oyr
                     , col.scheme = "inferno"
                     , vmin       = -Inf
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "het.resp"
                     , desc       = "Heterotrophic respiration"
                     , unit       = untab$kgcom2oyr
                     , col.scheme = "inferno"
                     , vmin       = 0
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "can.temp"
                     , desc       = "CAS temperature"
                     , unit       = untab$degC
                     , col.scheme = "inferno"
                     , vmin       = -t00
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "gnd.temp"
                     , desc       = "Ground temperature"
                     , unit       = untab$degC
                     , col.scheme = "inferno"
                     , vmin       = -t00
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "can.shv"
                     , desc       = "CAS specific humidity"
                     , unit       = untab$gwokg
                     , col.scheme = "inferno"
                     , vmin       = 0
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "gnd.shv"
                     , desc       = "Ground specific humidity"
                     , unit       = untab$gwokg
                     , col.scheme = "inferno"
                     , vmin       = 0
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "can.vpd"
                     , desc       = "CAS vapour pressure deficit"
                     , unit       = untab$hpa
                     , col.scheme = "inferno"
                     , vmin       = 0
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "can.co2"
                     , desc       = "CAS CO2 mixing ratio"
                     , unit       = untab$umolcomol
                     , col.scheme = "inferno"
                     , vmin       = 0
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "can.prss"
                     , desc       = "CAS pressure"
                     , unit       = untab$hpa
                     , col.scheme = "inferno"
                     , vmin       = 0
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "cflxca"
                     , desc       = "CO2 flux (CAS->Atm)"
                     , unit       = untab$umolcom2os
                     , col.scheme = "inferno"
                     , vmin       = -Inf
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "cflxst"
                     , desc       = "CO2 storage at CAS"
                     , unit       = untab$umolcom2os
                     , col.scheme = "inferno"
                     , vmin       = -Inf
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "nee"
                     , desc       = "Net ecosystem exchange"
                     , unit       = untab$umolcom2os
                     , col.scheme = "inferno"
                     , vmin       = -Inf
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "qwflxca"
                     , desc       = "'Latent' heat flux (CAS->Atm)"
                     , unit       = untab$wom2
                     , col.scheme = "inferno"
                     , vmin       = -Inf
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "hflxca"
                     , desc       = "Sensible heat flux (CAS->Atm)"
                     , unit       = untab$wom2
                     , col.scheme = "inferno"
                     , vmin       = -Inf
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "hflxgc"
                     , desc       = "Sensible heat flux (Grnd->CAS)"
                     , unit       = untab$wom2
                     , col.scheme = "inferno"
                     , vmin       = -Inf
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "hflxlc"
                     , desc       = "Sensible heat flux (Leaf->CAS)"
                     , unit       = untab$wom2
                     , col.scheme = "inferno"
                     , vmin       = -Inf
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "hflxwc"
                     , desc       = "Sensible heat flux (Wood->CAS)"
                     , unit       = untab$wom2
                     , col.scheme = "inferno"
                     , vmin       = -Inf
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "wflxca"
                     , desc       = "Water flux (CAS->Atm)"
                     , unit       = untab$kgwom2oday
                     , col.scheme = "inferno"
                     , vmin       = -Inf
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "wflxgc"
                     , desc       = "Water flux (Grnd->CAS)"
                     , unit       = untab$kgwom2oday
                     , col.scheme = "inferno"
                     , vmin       = -Inf
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "wflxlc"
                     , desc       = "Water flux (Leaf->CAS)"
                     , unit       = untab$kgwom2oday
                     , col.scheme = "inferno"
                     , vmin       = -Inf
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "wflxwc"
                     , desc       = "Water flux (Wood->CAS)"
                     , unit       = untab$kgwom2oday
                     , col.scheme = "inferno"
                     , vmin       = -Inf
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "transp"
                     , desc       = "Leaf transpiration"
                     , unit       = untab$kgwom2oday
                     , col.scheme = "inferno"
                     , vmin       = 0.
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "ustar"
                     , desc       = "Friction velocity"
                     , unit       = untab$mos
                     , col.scheme = "inferno"
                     , vmin       = 0.
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "rshortup"
                     , desc       = "Outgoing SW radiation"
                     , unit       = untab$wom2
                     , col.scheme = "inferno"
                     , vmin       = 0.
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "rlongup"
                     , desc       = "Outgoing LW radiation"
                     , unit       = untab$wom2
                     , col.scheme = "inferno"
                     , vmin       = 0.
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "parup"
                     , desc       = "Outgoing PAR"
                     , unit       = untab$umolom2os
                     , col.scheme = "inferno"
                     , vmin       = 0.
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "rnet"
                     , desc       = "Net radiation at ToC"
                     , unit       = untab$wom2
                     , col.scheme = "inferno"
                     , vmin       = -Inf
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "lai"
                     , desc       = "Leaf area index"
                     , unit       = untab$m2lom2
                     , col.scheme = "inferno"
                     , vmin       = 0.
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "leaf.temp"
                     , desc       = "Leaf temperature"
                     , unit       = untab$degC
                     , col.scheme = "inferno"
                     , vmin       = 0.
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "leaf.vpd"
                     , desc       = "Leaf vapour pressure deficit"
                     , unit       = untab$hpa
                     , col.scheme = "inferno"
                     , vmin       = 0.
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "wood.temp"
                     , desc       = "Wood temperature"
                     , unit       = untab$degC
                     , col.scheme = "inferno"
                     , vmin       = -t00
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "gpp"
                     , desc       = "Gross primary productivity"
                     , unit       = untab$kgcom2oyr
                     , col.scheme = "inferno"
                     , vmin       = 0.
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "npp"
                     , desc       = "Net primary productivity"
                     , unit       = untab$kgcom2oyr
                     , col.scheme = "inferno"
                     , vmin       = 0.
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "plant.resp"
                     , desc       = "Plant respiration"
                     , unit       = untab$kgcom2oyr
                     , col.scheme = "inferno"
                     , vmin       = 0.
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "reco"
                     , desc       = "Ecosystem respiration"
                     , unit       = untab$kgcom2oyr
                     , col.scheme = "inferno"
                     , vmin       = 0.
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "fast.grnd.c"
                     , desc       = "Aboveground litter"
                     , unit       = untab$kgcom2
                     , col.scheme = "inferno"
                     , vmin       = 0.
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "fast.soil.c"
                     , desc       = "Belowground litter"
                     , unit       = untab$kgcom2
                     , col.scheme = "inferno"
                     , vmin       = 0.
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "struct.grnd.c"
                     , desc       = "Aboveground woody debris"
                     , unit       = untab$kgcom2
                     , col.scheme = "inferno"
                     , vmin       = 0.
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "struct.soil.c"
                     , desc       = "Belowground woody debris"
                     , unit       = untab$kgcom2
                     , col.scheme = "inferno"
                     , vmin       = 0.
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "microbe.soil.c"
                     , desc       = "Microbial soil carbon"
                     , unit       = untab$kgcom2
                     , col.scheme = "inferno"
                     , vmin       = 0.
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "slow.soil.c"
                     , desc       = "Humified soil carbon"
                     , unit       = untab$kgcom2
                     , col.scheme = "inferno"
                     , vmin       = 0.
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "passive.soil.c"
                     , desc       = "Passive soil carbon"
                     , unit       = untab$kgcom2
                     , col.scheme = "inferno"
                     , vmin       = 0.
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "fgc.in"
                     , desc       = "Surface litter influx"
                     , unit       = untab$kgcom2oyr
                     , col.scheme = "magma"
                     , vmin       = 0.
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "fsc.in"
                     , desc       = "Sub-surface litter influx"
                     , unit       = untab$kgcom2oyr
                     , col.scheme = "magma"
                     , vmin       = 0.
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "stgc.in"
                     , desc       = "Surface woody debris influx"
                     , unit       = untab$kgcom2oyr
                     , col.scheme = "magma"
                     , vmin       = 0.
                     , vmax       = Inf
                     , plog       = TRUE
                     , emean      = TRUE
                     , mmean      = TRUE
                     , ymean      = TRUE
                     )#end list
n              = n + 1
plotpatch[[n]] = list( vnam       = "stsc.in"
                     , desc       = "Sub-surface woody debris influx"
                     , unit       = untab$kgcom2oyr
                     , col.scheme = "magma"
                     , vmin       = 0.
                     , vmax       = Inf
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
                   , unit      = untab$wom2
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "fast.soil.c"
                   , desc      = "Fast soil carbon"
                   , unit      = untab$kgcom2
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "struct.soil.c"
                   , desc      = "Structural soil carbon"
                   , unit      = untab$kgcom2
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "nep"
                   , desc      = "Net Ecosystem Production"
                   , unit      = untab$kgcom2
                   , plt       = TRUE
                   , cumul     = TRUE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "nee"
                   , desc      = "Net Ecosystem Exchange"
                   , unit      = untab$umolcom2os
                   , plt       = TRUE
                   , cumul     = TRUE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "mco"
                   , desc      = "Maintenance Costs"
                   , unit      = untab$kgcom2oyr
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "npp"
                   , desc      = "Net primary productivity"
                   , unit      = untab$kgcom2oyr
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "cba"
                   , desc      = "Carbon balance"
                   , unit      = untab$kgcom2
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "ldrop"
                   , desc      = "Leaf drop"
                   , unit      = untab$kgcom2oyr
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "cflxca"
                   , desc      = "Carbon flux"
                   , unit      = untab$umolcom2os
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "cflxst"
                   , desc      = "Carbon storage"
                   , unit      = untab$umolcom2os
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "evap"
                   , desc      = "Evaporation"
                   , unit      = untab$kgwom2oday
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "transp"
                   , desc      = "Transpiration"
                   , unit      = untab$kgwom2oday
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "ustar"
                   , desc      = "Friction velocity"
                   , unit      = untab$mos
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "hflxca"
                   , desc      = "Sensible heat flux"
                   , unit      = untab$wom2
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "wflxca"
                   , desc      = "Water vapour flux"
                   , unit      = untab$kgwom2oday
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "sm.stress"
                   , desc      = "Water Stress factor"
                   , unit      = untab$empty
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "gpp"
                   , desc      = "Gross Primary Productivity"
                   , unit      = untab$kgcom2oyr
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "reco"
                   , desc      = "Ecosystem respiration"
                   , unit      = untab$kgcom2oyr
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "rshortup"
                   , desc      = "Upward SW radiation"
                   , unit      = untab$wom2
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "parup"
                   , desc      = "Upward PAR"
                   , unit      = untab$umolom2os
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "rnet"
                   , desc      = "Net radiation at the top of the canopy"
                   , unit      = untab$wom2
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "rlong"
                   , desc      = "Incoming longwave radiation"
                   , unit      = untab$wom2
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "rlongup"
                   , desc      = "Outgoing longwave radiation"
                   , unit      = untab$wom2
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "cbamax"
                   , desc      = "Maximum carbon balance"
                   , unit      = untab$kgcom2
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "cbarel"
                   , desc      = "Relative carbon balance"
                   , unit      = untab$empty
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "f.gpp"
                   , desc      = "Gross Primary Productivity"
                   , unit      = untab$pcbiooyr
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "f.npp"
                   , desc      = "Net Primary Productivity"
                   , unit      = untab$pcbiooyr
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "f.cba"
                   , desc      = "Carbon balance"
                   , unit      = untab$pcbiooyr
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "f.bstorage"
                   , desc      = "Relative storage biomass"
                   , unit      = untab$gcokgcbio
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "f.bseeds"
                   , desc      = "Relative seed biomass"
                   , unit      = untab$gcokgcbio
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "f.bleaf"
                   , desc      = "Relative leaf biomass"
                   , unit      = untab$gcokgcbio
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
n            = n + 1
yc.tvar[[n]] = list( vnam      = "agb"
                   , desc      = "Above-ground biomass"
                   , unit      = untab$kgcom2
                   , plt       = TRUE
                   , cumul     = FALSE
                   )#end list
#------------------------------------------------------------------------------------------#




#----- XYZ plots, to explore the parameter space. -----------------------------------------#
yc.xyzvar      = list()
yc.xyzvar$zvar = list( list ( vname = "gpp"
                            , desc  = "Gross Primary Productivity"
                            , key   = "GPP"
                            , unit  = untab$kgcom2oyr
                            )#end list
                     , list ( vname = "reco"
                            , desc  = "Ecosystem respiration"
                            , key   = "RE"
                            , unit  = untab$kgcom2oyr
                            )#end list
                     , list ( vname = "plant.resp"
                            , desc  = "Plant respiration"
                            , key   = "PR"
                            , unit  = untab$kgcom2oyr
                            )#end list
                     , list ( vname = "het.resp"
                            , desc  = "Heterotrophic respiration"
                            , key   = "RH"
                            , unit  = untab$kgcom2oyr
                            )#end list
                     , list ( vname = "nep"
                            , desc  = "Net Ecosystem Productivity"
                            , key   = "NEP"
                            , unit  = untab$kgcom2oyr
                            )#end list
                     , list ( vname = "cba"
                            , desc  = "Carbon balance"
                            , key   = "CBA"
                            , unit  = untab$kgcom2oyr
                            )#end list
                     )#end list
yc.xyzvar$xvar = list( list ( vname = "rain"
                            , desc  = "Rainfall"
                            , unit  = untab$mmomo
                            , leg   = "right"
                            )#end list
                     , list ( vname = "demand"
                            , desc  = "Water demand"
                            , unit  = untab$kgwom2oday
                            , leg   = "left"
                            )#end list
                     , list ( vname = "sm.stress"
                            , desc  = "Soil moisture stress factor"
                            , unit  = untab$empty
                            , leg   = "left"
                            )#end list
                     )#end list
yc.xyzvar$yvar = list( list ( vname = "rshort"
                            , desc  = "Shortwave radiation"
                            , unit  = untab$wom2
                            , leg   = "top"
                            )#end list
                     , list ( vname = "leaf.temp"
                            , desc  = "Leaf temperature"
                            , unit  = untab$degC
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
