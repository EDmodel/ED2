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
   scen.ts[[n]] = list( vname    = "wood.dens"
                      , desc     = "Wood density"
                      , lname    = "Wood density"
                      , short    = "omega"
                      , unit     = untab$gocm3
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
   scen.ts[[n]] = list( vname    = "can.depth"
                      , desc     = "Mean height"
                      , lname    = "Height"
                      , short    = "z[c*a*n]"
                      , unit     = untab$m
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
   scen.ts[[n]] = list( vname    = "can.area"
                      , desc     = "Canopy fraction"
                      , lname    = "Canopy fraction"
                      , short    = "alpha[c*a*n]"
                      , unit     = untab$empty
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
   scen.ts[[n]] = list( vname    = "agb"
                      , desc     = "Above ground biomass"
                      , lname    = "AGB"
                      , short    = "A*G*B"
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
   scen.ts[[n]] = list( vname    = "bgb"
                      , desc     = "Below ground biomass"
                      , lname    = "BGB"
                      , short    = "B*G*B"
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
                      , lname    = "Total biomass"
                      , short    = "B[t*o*t]"
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
   scen.ts[[n]] = list( vname    = "lai"
                      , desc     = "Leaf area index"
                      , lname    = "LAI"
                      , short    = "L*A*I"
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
                      , lname    = "Plant density"
                      , short    = "n[p*l]"
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
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "gpp"
                      , desc     = "Gross primary productivity"
                      , lname    = "GPP"
                      , short    = "G*P*P"
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
   scen.ts[[n]] = list( vname    = "last.1yr.gpp"
                      , desc     = "GPP of the past 12 months"
                      , lname    = "GPP"
                      , short    = "G*P*P"
                      , unit     = untab$kgcom2oyr
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.2yr.gpp"
                      , desc     = "GPP of the past 24 months"
                      , lname    = "GPP"
                      , short    = "G*P*P"
                      , unit     = untab$kgcom2oyr
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.3yr.gpp"
                      , desc     = "GPP of the past 36 months"
                      , lname    = "GPP"
                      , short    = "G*P*P"
                      , unit     = untab$kgcom2oyr
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "mco"
                      , desc     = "Maintenance costs"
                      , lname    = "Maintenance"
                      , short    = "dot(B)[M*C*o]"
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
                      , lname    = "C Balance"
                      , short    = "C[B*a*l]"
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
   scen.ts[[n]] = list( vname    = "last.1yr.cba"
                      , desc     = "CBA of the past 12 months"
                      , lname    = "C Balance"
                      , short    = "C[B*a*l]"
                      , unit     = untab$kgcom2
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.2yr.cba"
                      , desc     = "CBA of the past 24 months"
                      , lname    = "C Balance"
                      , short    = "C[B*a*l]"
                      , unit     = untab$kgcom2
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.3yr.cba"
                      , desc     = "CBA of the past 36 months"
                      , lname    = "C Balance"
                      , short    = "C[B*a*l]"
                      , unit     = untab$kgcom2
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "ldrop"
                      , desc     = "Leaf drop"
                      , lname    = "Leaf drop"
                      , short    = "dot(B)[L*D*r]"
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
                      , lname    = "Storage"
                      , short    = "B[S*t*o*r]"
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
                      , lname    = "Seed"
                      , short    = "B[S*e*e*d]"
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
                      , lname    = "Basal area"
                      , short    = "B*A"
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
   scen.ts[[n]] = list( vname    = "sm.stress"
                      , desc     = "Soil moisture stress factor"
                      , lname    = "SM Stress"
                      , short    = "beta[S*W]"
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
   scen.ts[[n]] = list( vname    = "phap.sms"
                      , desc     = "PhAP Soil moisture stress"
                      , lname    = "SM Stress"
                      , short    = "beta[S*W]"
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
   scen.ts[[n]] = list( vname    = "last.1yr.sms"
                      , desc     = "SMS of the past 12 months"
                      , lname    = "SM Stress"
                      , short    = "beta[S*W]"
                      , unit     = untab$empty
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.2yr.sms"
                      , desc     = "SMS of the past 24 months"
                      , lname    = "SM Stress"
                      , short    = "beta[S*W]"
                      , unit     = untab$empty
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.3yr.sms"
                      , desc     = "SMS of the past 36 months"
                      , lname    = "SM Stress"
                      , short    = "beta[S*W]"
                      , unit     = untab$empty
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "mort"
                      , desc     = "Mortality rate"
                      , lname    = "Mortality"
                      , short    = "dot(mu)"
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
                      , lname    = "DD Mortality"
                      , short    = "dot(mu)[D*D]"
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
                      , lname    = "DI Mortality"
                      , short    = "dot(mu)[D*I]"
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
                      , lname    = "Recruitment"
                      , short    = "dot(rho)"
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
                      , lname    = "Growth"
                      , short    = "dot(gamma)"
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
   scen.ts[[n]] = list( vname    = "change"
                      , desc     = "Population change"
                      , lname    = "Growth"
                      , short    = "dot(n)[p*l]"
                      , unit     = untab$oneoyr
                      , f.aggr   = "mean"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = TRUE
                      , dbhvar   = TRUE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = TRUE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "agb.mort"
                      , desc     = "Mortality rate"
                      , lname    = "Mortality"
                      , short    = "dot(mu)"
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
   scen.ts[[n]] = list( vname    = "last.1yr.mort"
                      , desc     = "AGB mortality of the past 12 months"
                      , lname    = "Mortality"
                      , short    = "dot(mu)"
                      , unit     = untab$pcagboyr
                      , f.aggr   = "get.last"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = TRUE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "last.2yr.mort"
                      , desc     = "AGB mortality of the past 24 months"
                      , lname    = "Mortality"
                      , short    = "dot(mu)"
                      , unit     = untab$pcagboyr
                      , f.aggr   = "get.last"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = TRUE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "last.3yr.mort"
                      , desc     = "AGB mortality of the past 36 months"
                      , lname    = "Mortality"
                      , short    = "dot(mu)"
                      , unit     = untab$pcagboyr
                      , f.aggr   = "get.last"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = TRUE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "agb.ncbmort"
                      , desc     = "Mortality rate - Neg. C balance"
                      , lname    = "Mortality"
                      , short    = "dot(mu)[D*D]"
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
   scen.ts[[n]] = list( vname    = "last.1yr.ncbmort"
                      , desc     = "NCB mortality of the past 12 months"
                      , lname    = "DD Mortality"
                      , short    = "dot(mu)[D*D]"
                      , unit     = untab$pcagboyr
                      , f.aggr   = "get.last"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = TRUE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "last.2yr.ncbmort"
                      , desc     = "NCB mortality of the past 24 months"
                      , lname    = "DD Mortality"
                      , short    = "dot(mu)[D*D]"
                      , unit     = untab$pcagboyr
                      , f.aggr   = "get.last"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = TRUE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "last.3yr.ncbmort"
                      , desc     = "NCB mortality of the past 36 months"
                      , lname    = "DD Mortality"
                      , short    = "dot(mu)[D*D]"
                      , unit     = untab$pcagboyr
                      , f.aggr   = "get.last"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = TRUE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "agb.dimort"
                      , desc     = "Mortality rate - Density-independent"
                      , lname    = "DI Mortality"
                      , short    = "dot(mu)[D*I]"
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
   scen.ts[[n]] = list( vname    = "last.1yr.dimort"
                      , desc     = "DI mortality of the past 12 months"
                      , lname    = "DI Mortality"
                      , short    = "dot(mu)[D*I]"
                      , unit     = untab$pcagboyr
                      , f.aggr   = "get.last"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = TRUE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "last.2yr.dimort"
                      , desc     = "DI mortality of the past 24 months"
                      , lname    = "DI Mortality"
                      , short    = "dot(mu)[D*I]"
                      , unit     = untab$pcagboyr
                      , f.aggr   = "get.last"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = TRUE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "last.3yr.dimort"
                      , desc     = "DI mortality of the past 36 months"
                      , lname    = "DI Mortality"
                      , short    = "dot(mu)[D*I]"
                      , unit     = untab$pcagboyr
                      , f.aggr   = "get.last"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = TRUE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "agb.recr"
                      , desc     = "Recruitment rate"
                      , lname    = "Recruitment"
                      , short    = "dot(rho)"
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
   scen.ts[[n]] = list( vname    = "last.1yr.recr"
                      , desc     = "Recruitment of the past 12 months"
                      , lname    = "Recruitment"
                      , short    = "dot(rho)"
                      , unit     = untab$pcagboyr
                      , f.aggr   = "get.last"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = FALSE
                      , recr     = TRUE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "last.2yr.recr"
                      , desc     = "Recruitment of the past 24 months"
                      , lname    = "Recruitment"
                      , short    = "dot(rho)"
                      , unit     = untab$pcagboyr
                      , f.aggr   = "get.last"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = FALSE
                      , recr     = TRUE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "last.3yr.recr"
                      , desc     = "Recruitment of the past 36 months"
                      , lname    = "Recruitment"
                      , short    = "dot(rho)"
                      , unit     = untab$pcagboyr
                      , f.aggr   = "get.last"
                      , add      = 0
                      , mult     = 1
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = FALSE
                      , recr     = TRUE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "agb.growth"
                      , desc     = "Growth rate (AGB)"
                      , lname    = "Growth"
                      , short    = "dot(gamma)"
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
   scen.ts[[n]] = list( vname    = "last.1yr.growth"
                      , desc     = "Growth of the past 12 months"
                      , lname    = "Growth"
                      , short    = "dot(gamma)"
                      , unit     = untab$pcagboyr
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.2yr.growth"
                      , desc     = "Growth of the past 24 months"
                      , lname    = "Growth"
                      , short    = "dot(gamma)"
                      , unit     = untab$pcagboyr
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.3yr.growth"
                      , desc     = "Growth of the past 36 months"
                      , lname    = "Growth"
                      , short    = "dot(gamma)"
                      , unit     = untab$pcagboyr
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "bsa.mort"
                      , desc     = "Mortality rate"
                      , lname    = "Mortality"
                      , short    = "dot(mu)"
                      , unit     = untab$pcbaoyr
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
   scen.ts[[n]] = list( vname    = "bsa.ncbmort"
                      , desc     = "Mortality rate - Neg. C balance"
                      , lname    = "Mortality"
                      , short    = "dot(mu)[D*D]"
                      , unit     = untab$pcbaoyr
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
   scen.ts[[n]] = list( vname    = "bsa.dimort"
                      , desc     = "Mortality rate - Density-independent"
                      , lname    = "Mortality"
                      , short    = "dot(mu)[D*I]"
                      , unit     = untab$pcbaoyr
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
   scen.ts[[n]] = list( vname    = "bsa.recr"
                      , desc     = "Recruitment rate"
                      , lname    = "Recruitment"
                      , short    = "dot(rho)"
                      , unit     = untab$pcbaoyr
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
   scen.ts[[n]] = list( vname    = "bsa.growth"
                      , desc     = "Growth rate"
                      , lname    = "Growth"
                      , short    = "dot(gamma)"
                      , unit     = untab$pcbaoyr
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
                      , lname    = "Maximum CB"
                      , short    = "C*B[m*a*x]"
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
   scen.ts[[n]] = list( vname    = "cbarel"
                      , desc     = "Relative carbon balance"
                      , lname    = "Relative CB"
                      , short    = "C*B[r*e*l]"
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
                      , lname    = "Water demand"
                      , short    = "W[D*e*m]"
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
                      , lname    = "Water supply"
                      , short    = "W[S*u*p]"
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
                      , lname    = "Precipitation"
                      , short    = "dot(R)"
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
                      , lname    = "Runoff"
                      , short    = "dot(W)[R*O]"
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
                      , lname    = "Intercepted"
                      , short    = "dot(W)[I*n*t]"
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
                      , lname    = "Rainfall"
                      , short    = "dot(R)"
                      , unit     = untab$mmoyr
                      , f.aggr   = "get.last"
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
                      , lname    = "Rainfall"
                      , short    = "dot(R)"
                      , unit     = untab$mmoyr
                      , f.aggr   = "get.last"
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
                      , lname    = "Rainfall"
                      , short    = "dot(R)"
                      , unit     = untab$mmoyr
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "nmon.lt.090"
                      , desc     = "Drought length"
                      , lname    = "Drought length"
                      , short    = "t[d*r*g*t]"
                      , unit     = untab$nmo.090
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
   scen.ts[[n]] = list( vname    = "nmon.lt.100"
                      , desc     = "Drought length"
                      , lname    = "Drought length"
                      , short    = "t[d*r*g*t]"
                      , unit     = untab$nmo.100
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
   scen.ts[[n]] = list( vname    = "nmon.lt.110"
                      , desc     = "Drought length"
                      , lname    = "Drought length"
                      , short    = "t[d*r*g*t]"
                      , unit     = untab$nmo.110
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
   scen.ts[[n]] = list( vname    = "nmon.lt.120"
                      , desc     = "Drought length"
                      , lname    = "Drought length"
                      , short    = "t[d*r*g*t]"
                      , unit     = untab$nmo.120
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
   scen.ts[[n]] = list( vname    = "nmon.wdef"
                      , desc     = "Drought length (ET)"
                      , lname    = "Drought length"
                      , short    = "t[d*r*g*t]"
                      , unit     = untab$nmo.wdef
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
   scen.ts[[n]] = list( vname    = "nmon.mdef"
                      , desc     = "Drought length (Malhi)"
                      , lname    = "Drought length"
                      , short    = "t[d*r*g*t]"
                      , unit     = untab$nmo.wdef
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
   scen.ts[[n]] = list( vname    = "npp"
                      , desc     = "Net Primary Production"
                      , lname    = "NPP"
                      , short    = "N*P*P"
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
   scen.ts[[n]] = list( vname    = "last.1yr.npp"
                      , desc     = "NPP of the past 12 months"
                      , lname    = "NPP"
                      , short    = "N*P*P"
                      , unit     = untab$kgcom2oyr
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.2yr.npp"
                      , desc     = "NPP of the past 24 months"
                      , lname    = "NPP"
                      , short    = "N*P*P"
                      , unit     = untab$kgcom2oyr
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.3yr.npp"
                      , desc     = "NPP of the past 36 months"
                      , lname    = "NPP"
                      , short    = "N*P*P"
                      , unit     = untab$kgcom2oyr
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "dcbadt"
                      , desc     = "Change in C Balance"
                      , lname    = "Delta C Bal."
                      , short    = "Delta*C[B*a*l]"
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
   scen.ts[[n]] = list( vname    = "last.1yr.dcbadt"
                      , desc     = "Delta CB of the past 12 months"
                      , lname    = "Delta C Bal."
                      , short    = "Delta*C[B*a*l]"
                      , unit     = untab$kgcom2oyr
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.2yr.dcbadt"
                      , desc     = "Delta CB of the past 24 months"
                      , lname    = "Delta C Bal."
                      , short    = "Delta*C[B*a*l]"
                      , unit     = untab$kgcom2oyr
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.3yr.dcbadt"
                      , desc     = "Delta CB of the past 36 months"
                      , lname    = "Delta C Bal."
                      , short    = "Delta*C[B*a*l]"
                      , unit     = untab$kgcom2oyr
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "rue"
                      , desc     = "Rainfall Use Efficiency"
                      , lname    = "RUE"
                      , short    = "R*U*E"
                      , unit     = untab$gcokgw
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.1yr.rue"
                      , desc     = "RUE of the past 12 months"
                      , lname    = "RUE"
                      , short    = "R*U*E"
                      , unit     = untab$gcokgw
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.2yr.rue"
                      , desc     = "RUE of the past 24 months"
                      , lname    = "RUE"
                      , short    = "R*U*E"
                      , unit     = untab$gcokgw
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.3yr.rue"
                      , desc     = "RUE of the past 36 months"
                      , lname    = "RUE"
                      , short    = "R*U*E"
                      , unit     = untab$gcokgw
                      , f.aggr   = "get.last"
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
                      , lname    = "MWD"
                      , short    = "-Delta*W[m*a*x]"
                      , unit     = untab$mm
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
   scen.ts[[n]] = list( vname    = "last.1yr.mwd"
                      , desc     = "MWD of the past 12 months"
                      , lname    = "MWD"
                      , short    = "-Delta*W[m*a*x]"
                      , unit     = untab$mm
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.2yr.mwd"
                      , desc     = "MWD of the past 24 months"
                      , lname    = "MWD"
                      , short    = "-Delta*W[m*a*x]"
                      , unit     = untab$mm
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.3yr.mwd"
                      , desc     = "MWD of the past 36 months"
                      , lname    = "MWD"
                      , short    = "-Delta*W[m*a*x]"
                      , unit     = untab$mm
                      , f.aggr   = "get.last"
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
                      , lname    = "MWD"
                      , short    = "-Delta*W[m*a*x]"
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
                      , lname    = "Air Temp."
                      , short    = "T[a*i*r]"
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
                      , lname    = "Leaf Temp."
                      , short    = "T[l*e*a*f]"
                      , unit     = untab$degC
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
   scen.ts[[n]] = list( vname    = "phap.ltemp"
                      , desc     = "PhAP Leaf temperature"
                      , lname    = "Leaf Temp."
                      , short    = "T[l*e*a*f]"
                      , unit     = untab$degC
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
   scen.ts[[n]] = list( vname    = "last.1yr.ltemp"
                      , desc     = "Leaf temp. of past 12 months"
                      , lname    = "Leaf Temp."
                      , short    = "T[l*e*a*f]"
                      , unit     = untab$degC
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.2yr.ltemp"
                      , desc     = "Leaf temp. of past 24 months"
                      , lname    = "Leaf Temp."
                      , short    = "T[l*e*a*f]"
                      , unit     = untab$degC
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.3yr.ltemp"
                      , desc     = "Leaf temp. of past 36 months"
                      , lname    = "Leaf Temp."
                      , short    = "T[l*e*a*f]"
                      , unit     = untab$degC
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "leaf.water"
                      , desc     = "Leaf intercepted water"
                      , lname    = "Leaf Water"
                      , short    = "W[l*e*a*f]"
                      , unit     = untab$kgwom2l
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
   scen.ts[[n]] = list( vname    = "phap.lwater"
                      , desc     = "PhAP Leaf intercepted water"
                      , lname    = "Leaf Water"
                      , short    = "W[l*e*a*f]"
                      , unit     = untab$kgwom2l
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
   scen.ts[[n]] = list( vname    = "last.1yr.lwater"
                      , desc     = "Leaf water of past 12 months"
                      , lname    = "Leaf Water"
                      , short    = "W[l*e*a*f]"
                      , unit     = untab$kgwom2l
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.1yr.lwater"
                      , desc     = "Leaf water of past 24 months"
                      , lname    = "Leaf Water"
                      , short    = "W[l*e*a*f]"
                      , unit     = untab$kgwom2l
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.3yr.lwater"
                      , desc     = "Leaf water of past 36 months"
                      , lname    = "Leaf Water"
                      , short    = "W[l*e*a*f]"
                      , unit     = untab$kgwom2l
                      , f.aggr   = "get.last"
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
                      , lname    = "SW Rad."
                      , short    = "S*W*symbol(\"\335\")"
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
   scen.ts[[n]] = list( vname    = "last.1yr.rshort"
                      , desc     = "SW of the past 12 months"
                      , lname    = "SW Rad."
                      , short    = "S*W*symbol(\"\335\")"
                      , unit     = untab$wom2
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.2yr.rshort"
                      , desc     = "SW of the past 24 months"
                      , lname    = "SW Rad."
                      , short    = "S*W*symbol(\"\335\")"
                      , unit     = untab$wom2
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.3yr.rshort"
                      , desc     = "SW of the past 36 months"
                      , lname    = "SW Rad."
                      , short    = "S*W*symbol(\"\335\")"
                      , unit     = untab$wom2
                      , f.aggr   = "get.last"
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
                      , lname    = "LW Rad."
                      , short    = "L*W*symbol(\"\335\")"
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
                      , lname    = "AVPD"
                      , short    = "A*V*P*D"
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
                      , lname    = "LVPD"
                      , short    = "L*V*P*D"
                      , unit     = untab$hpa
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
   scen.ts[[n]] = list( vname    = "phap.lvpd"
                      , desc     = "PhAP Leaf VPD"
                      , lname    = "LVPD"
                      , short    = "L*V*P*D"
                      , unit     = untab$hpa
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
   scen.ts[[n]] = list( vname    = "last.1yr.lvpd"
                      , desc     = "Leaf VPD of the past 12 months"
                      , lname    = "LVPD"
                      , short    = "L*V*P*D"
                      , unit     = untab$hpa
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.2yr.lvpd"
                      , desc     = "Leaf VPD of the past 24 months"
                      , lname    = "LVPD"
                      , short    = "L*V*P*D"
                      , unit     = untab$hpa
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.3yr.lvpd"
                      , desc     = "Leaf VPD of the past 36 months"
                      , lname    = "LVPD"
                      , short    = "L*V*P*D"
                      , unit     = untab$hpa
                      , f.aggr   = "get.last"
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
                      , lname    = "PAW"
                      , short    = "P*A*W"
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
                      , lname    = "Matric Potl."
                      , short    = "Psi[m]"
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
   scen.ts[[n]] = list( vname    = "last.1yr.smpot"
                      , desc     = "SMPot of the past 12 months"
                      , lname    = "Matric Potl."
                      , short    = "Psi[m]"
                      , unit     = untab$mpa
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.2yr.smpot"
                      , desc     = "SMPot of the past 24 months"
                      , lname    = "Matric Potl."
                      , short    = "Psi[m]"
                      , unit     = untab$mpa
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.3yr.smpot"
                      , desc     = "SMPot of the past 36 months"
                      , lname    = "Matric Potl."
                      , short    = "Psi[m]"
                      , unit     = untab$mpa
                      , f.aggr   = "get.last"
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
                      , lname    = "NEP"
                      , short    = "N*E*P"
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
                      , lname    = "Ecos. Resp."
                      , short    = "R[E*c*o]"
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
                      , lname    = "FSC"
                      , short    = "C[f*a*s*t]"
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
                      , lname    = "StSC"
                      , short    = "C[s*t*r*u*c*t]"
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
                      , lname    = "SSC"
                      , short    = "C[s*l*o*w]"
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
                      , lname    = "Het. Resp."
                      , short    = "R[H*e*t]"
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
                      , lname    = "Auto. Resp."
                      , short    = "R[A*u*t*o]"
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
   scen.ts[[n]] = list( vname    = "last.1yr.plresp"
                      , desc     = "Pl. Resp. of the past 12 months"
                      , lname    = "Auto. Resp."
                      , short    = "R[A*u*t*o]"
                      , unit     = untab$kgcom2oyr
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.2yr.plresp"
                      , desc     = "Pl. Resp. of the past 24 months"
                      , lname    = "Auto. Resp."
                      , short    = "R[A*u*t*o]"
                      , unit     = untab$kgcom2oyr
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.3yr.plresp"
                      , desc     = "Pl. Resp. of the past 36 months"
                      , lname    = "Auto. Resp."
                      , short    = "R[A*u*t*o]"
                      , unit     = untab$kgcom2oyr
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "hflxlc"
                      , desc     = "Leaf sensible heat"
                      , lname    = "Sensible heat"
                      , short    = "dot(theta)[L*C]"
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
                      , lname    = "Leaf Evap."
                      , short    = "dot(epsilon)[L*e*a*f]"
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
                      , lname    = "Leaf Transp."
                      , short    = "dot(tau)"
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
   scen.ts[[n]] = list( vname    = "last.1yr.transp"
                      , desc     = "Transp. of the past 12 months"
                      , lname    = "Leaf Transp."
                      , short    = "dot(tau)"
                      , unit     = untab$kgwom2oday
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.2yr.transp"
                      , desc     = "Transp. of the past 24 months"
                      , lname    = "Leaf Transp."
                      , short    = "dot(tau)"
                      , unit     = untab$kgwom2oday
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.3yr.transp"
                      , desc     = "Transp. of the past 36 months"
                      , lname    = "Leaf Transp."
                      , short    = "dot(tau)"
                      , unit     = untab$kgwom2oday
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "et"
                      , desc     = "Evapotranspiration"
                      , lname    = "ET"
                      , short    = "dot(epsilon)+dot(tau)"
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
   scen.ts[[n]] = list( vname    = "last.1yr.et"
                      , desc     = "ET of the past 12 months"
                      , lname    = "ET"
                      , short    = "dot(epsilon)+dot(tau)"
                      , unit     = untab$kgwom2oday
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.2yr.et"
                      , desc     = "ET of the past 24 months"
                      , lname    = "ET"
                      , short    = "dot(epsilon)+dot(tau)"
                      , unit     = untab$kgwom2oday
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.3yr.et"
                      , desc     = "ET of the past 36 months"
                      , lname    = "ET"
                      , short    = "dot(epsilon)+dot(tau)"
                      , unit     = untab$kgwom2oday
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "evap"
                      , desc     = "Total evaporation"
                      , lname    = "ET"
                      , short    = "dot(epsilon)"
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
   scen.ts[[n]] = list( vname    = "last.1yr.evap"
                      , desc     = "Evap. of the past 12 months"
                      , lname    = "ET"
                      , short    = "dot(epsilon)"
                      , unit     = untab$kgwom2oday
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.2yr.evap"
                      , desc     = "Evap. of the past 24 months"
                      , lname    = "ET"
                      , short    = "dot(epsilon)"
                      , unit     = untab$kgwom2oday
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.3yr.evap"
                      , desc     = "Evap. of the past 36 months"
                      , lname    = "ET"
                      , short    = "dot(epsilon)"
                      , unit     = untab$kgwom2oday
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "i.gpp"
                      , desc     = "Mean Gross Primary Production"
                      , lname    = "GPP"
                      , short    = "G*P*P[I*n*d]"
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
                      , lname    = "NPP"
                      , short    = "bar(N*P*P)"
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
                      , lname    = "Auto. Resp."
                      , short    = "bar(R[A*u*t*o])"
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
                      , lname    = "C Balance"
                      , short    = "bar(C[B*a*l])"
                      , unit     = untab$kgcopl
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
   scen.ts[[n]] = list( vname    = "i.mco"
                      , desc     = "Mean Maintenance Costs"
                      , lname    = "Maintenance"
                      , short    = "bar(B[M*C*o])"
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
                      , lname    = "Sens. heat"
                      , short    = "dot(theta)[L*e*a*f]"
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
   scen.ts[[n]] = list( vname    = "i.wflxlc"
                      , desc     = "Mean Leaf evaporation"
                      , lname    = "Leaf evap."
                      , short    = "dot(theta)[L*e*a*f]"
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
   scen.ts[[n]] = list( vname    = "i.transp"
                      , desc     = "Mean Transpiration"
                      , lname    = "Leaf transp."
                      , short    = "dot(tau)"
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
   scen.ts[[n]] = list( vname    = "hflxgc"
                      , desc     = "Sensible heat - Gnd->CAS"
                      , lname    = "Gnd. Sens. heat"
                      , short    = "dot(theta)[G*n*d]"
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
                      , lname    = "EF - sensible"
                      , short    = "dot(theta)[E*F]"
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
                      , lname    = "Ground Evap."
                      , short    = "dot(epsilon)[G*n*d]"
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
                      , lname    = "EF - water"
                      , short    = "dot(epsilon)[E*F]"
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
   scen.ts[[n]] = list( vname    = "cue"
                      , desc     = "Carbon use efficiency"
                      , lname    = "CUE"
                      , short    = "C*U*E"
                      , unit     = untab$kgcokgc
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.1yr.cue"
                      , desc     = "CUE of the past 12 months"
                      , lname    = "CUE"
                      , short    = "C*U*E"
                      , unit     = untab$kgcokgc
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.2yr.cue"
                      , desc     = "CUE of the past 24 months"
                      , lname    = "CUE"
                      , short    = "C*U*E"
                      , unit     = untab$kgcokgc
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.3yr.cue"
                      , lname    = "CUE"
                      , desc     = "CUE of the past 36 months"
                      , short    = "C*U*E"
                      , unit     = untab$kgcokgc
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "ecue"
                      , desc     = "Effective Carbon use efficiency"
                      , lname    = "CUE"
                      , short    = "C*U*E"
                      , unit     = untab$kgcokgc
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.1yr.ecue"
                      , desc     = "ECUE of the past 12 months"
                      , lname    = "CUE"
                      , short    = "C*U*E"
                      , unit     = untab$kgcokgc
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.2yr.ecue"
                      , desc     = "ECUE of the past 24 months"
                      , lname    = "CUE"
                      , short    = "C*U*E"
                      , unit     = untab$kgcokgc
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.3yr.ecue"
                      , desc     = "ECUE of the past 36 months"
                      , lname    = "CUE"
                      , short    = "C*U*E"
                      , unit     = untab$kgcokgc
                      , f.aggr   = "get.last"
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
                      , desc     = "Actual Water use efficiency"
                      , lname    = "LWUE"
                      , short    = "L*W*U*E"
                      , unit     = untab$gcokgw
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.1yr.wue"
                      , desc     = "LWUE of the past 12 months"
                      , lname    = "LWUE"
                      , short    = "L*W*U*E"
                      , unit     = untab$gcokgw
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.2yr.wue"
                      , desc     = "LWUE of the past 24 months"
                      , lname    = "LWUE"
                      , short    = "L*W*U*E"
                      , unit     = untab$gcokgw
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.3yr.wue"
                      , desc     = "LWUE of the past 36 months"
                      , lname    = "LWUE"
                      , short    = "L*W*U*E"
                      , unit     = untab$gcokgw
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "etue"
                      , desc     = "Bulk water use efficiency"
                      , lname    = "BWUE"
                      , short    = "B*W*U*E"
                      , unit     = untab$gcokgw
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.1yr.etue"
                      , desc     = "BWUE of the past 12 months"
                      , lname    = "BWUE"
                      , short    = "B*W*U*E"
                      , unit     = untab$gcokgw
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.2yr.etue"
                      , desc     = "BWUE of the past 24 months"
                      , lname    = "BWUE"
                      , short    = "B*W*U*E"
                      , unit     = untab$gcokgw
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.3yr.etue"
                      , desc     = "BWUE of the past 36 months"
                      , lname    = "BWUE"
                      , short    = "B*W*U*E"
                      , unit     = untab$gcokgw
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "leaf.gbw"
                      , desc     = "Leaf Bnd. Lyr. Conductance"
                      , lname    = "LBL Condct."
                      , short    = "g[b*w]"
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
   scen.ts[[n]] = list( vname    = "phap.lgbw"
                      , desc     = "PhAP Leaf Bnd. Lyr. Conduct."
                      , lname    = "LBL Condct."
                      , short    = "g[b*w]"
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
                      , lname    = "Stom. Condct."
                      , short    = "g[s*w]"
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
   scen.ts[[n]] = list( vname    = "phap.lgsw"
                      , desc     = "PhAP Leaf stomatal Conductance"
                      , lname    = "Stom. Condct."
                      , short    = "g[s*w]"
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
   scen.ts[[n]] = list( vname    = "last.1yr.lgsw"
                      , desc     = "GSW of the past 12 months"
                      , lname    = "Stom. Condct."
                      , short    = "g[s*w]"
                      , unit     = untab$kgwom2loday
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.2yr.lgsw"
                      , desc     = "GSW of the past 24 months"
                      , lname    = "Stom. Condct."
                      , short    = "g[s*w]"
                      , unit     = untab$kgwom2loday
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.3yr.lgsw"
                      , desc     = "GSW of the past 36 months"
                      , lname    = "Stom. Condct."
                      , short    = "g[s*w]"
                      , unit     = untab$kgwom2loday
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "f.gpp"
                      , desc     = "Gross Primary Productivity"
                      , lname    = "GPP"
                      , short    = "G*P*P"
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
                      , lname    = "Auto. Resp."
                      , short    = "R[A*u*t*o]"
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
                      , lname    = "NPP"
                      , short    = "N*P*P"
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
   scen.ts[[n]] = list( vname    = "f.dcbadt"
                      , desc     = "Change in C Balance"
                      , lname    = "Delta C Bal."
                      , short    = "Delta*C[B*a*l]"
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
   scen.ts[[n]] = list( vname    = "f.mco"
                      , desc     = "Maintenance costs"
                      , lname    = "Maintenance"
                      , short    = "dot(B)[M*C*o]"
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
                      , lname    = "C Balance"
                      , short    = "C[B*a*l]"
                      , unit     = untab$kgcokgc
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
                      , lname    = "Storage"
                      , short    = "B[S*t*o*r]"
                      , unit     = untab$kgcokgc
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
                      , lname    = "Leaf"
                      , short    = "B[L*e*a*f]"
                      , unit     = untab$kgcokgc
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
                      , lname    = "Root"
                      , short    = "B[R*o*o*t]"
                      , unit     = untab$kgcokgc
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
                      , lname    = "Seed"
                      , short    = "B[S*e*e*d]"
                      , unit     = untab$kgcokgc
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
                      , lname    = "Leaf PAR"
                      , short    = "P*A*R[L*e*a*f]"
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
   scen.ts[[n]] = list( vname    = "phap.lpar"
                      , desc     = "PhAP Absorbed PAR - Leaf"
                      , lname    = "Leaf PAR"
                      , short    = "P*A*R[L*e*a*f]"
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
   scen.ts[[n]] = list( vname    = "last.1yr.lpar"
                      , desc     = "LPAR of the past 12 months"
                      , lname    = "Leaf PAR"
                      , short    = "P*A*R[L*e*a*f]"
                      , unit     = untab$umolom2los
                      , f.aggr   = "get.last"
                      , add      = 0
                      , mult     = 1.e6
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "last.2yr.lpar"
                      , desc     = "LPAR of the past 24 months"
                      , lname    = "Leaf PAR"
                      , short    = "P*A*R[L*e*a*f]"
                      , unit     = untab$umolom2los
                      , f.aggr   = "get.last"
                      , add      = 0
                      , mult     = 1.e6
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "last.3yr.lpar"
                      , desc     = "LPAR of the past 36 months"
                      , lname    = "Leaf PAR"
                      , short    = "P*A*R[L*e*a*f]"
                      , unit     = untab$umolom2los
                      , f.aggr   = "get.last"
                      , add      = 0
                      , mult     = 1.e6
                      , pftvar   = FALSE
                      , dbhvar   = FALSE
                      , mort     = FALSE
                      , recr     = FALSE
                      , plog     = FALSE
                      , plog.dbh = FALSE
                      , plt      = TRUE
                      )#end list
   n            = n + 1
   scen.ts[[n]] = list( vname    = "leaf.rshort"
                      , desc     = "Absorbed SW - Leaf"
                      , lname    = "Leaf SW"
                      , short    = "S*W[L*e*a*f]"
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
                      , lname    = "Leaf LW"
                      , short    = "L*W[L*e*a*f]"
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
   scen.ts[[n]] = list( vname    = "agb.change"
                      , desc     = "AGB change"
                      , lname    = "AGB change"
                      , short    = "delta[A*G*B]"
                      , unit     = untab$oneoyr
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
   scen.ts[[n]] = list( vname    = "last.1yr.change"
                      , desc     = "AGB change of the past 12 months"
                      , lname    = "AGB change"
                      , short    = "delta[A*G*B]"
                      , unit     = untab$oneoyr
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.2yr.change"
                      , desc     = "AGB change of the past 24 months"
                      , lname    = "AGB change"
                      , short    = "delta[A*G*B]"
                      , unit     = untab$oneoyr
                      , f.aggr   = "get.last"
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
   scen.ts[[n]] = list( vname    = "last.3yr.change"
                      , desc     = "AGB change of the past 36 months"
                      , lname    = "AGB change"
                      , short    = "delta[A*G*B]"
                      , unit     = untab$oneoyr
                      , f.aggr   = "mean.log"
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
   scen.szpft = data.frame( vname = c(            "agb",            "lai",             "ba"
                                     ,           "recr",           "mort",        "ncbmort"
                                     ,         "dimort",         "growth",         "change"
                                     ,       "agb.recr",       "agb.mort",    "agb.ncbmort"
                                     ,     "agb.dimort",     "agb.growth",     "agb.change"
                                     ,            "gpp",            "npp",     "plant.resp"
                                     ,         "dcbadt",            "cba",            "mco"
                                     ,       "bstorage",         "bseeds",      "sm.stress"
                                     ,       "phap.sms",         "wflxlc",         "transp"
                                     ,          "i.gpp",   "i.plant.resp",          "i.npp"
                                     ,          "i.cba",       "i.transp",         "cbamax"
                                     ,       "i.hflxlc",       "i.wflxlc",          "i.mco"
                                     ,       "leaf.gsw",      "phap.lgsw",          "f.gpp"
                                     ,          "f.npp",          "f.mco",       "f.dcbadt"
                                     ,          "f.cba",     "f.bstorage",        "f.bleaf"
                                     ,       "f.bseeds",           "rain",         "runoff"
                                     ,    "intercepted",        "atm.vpd",       "leaf.vpd"
                                     ,      "phap.lvpd",     "leaf.water",      "leaf.temp"
                                     ,     "phap.ltemp",       "atm.temp",       "leaf.par"
                                     ,      "phap.lpar",         "nplant",      "nmon.wdef"
                                     ,      "wood.dens",   "last.1yr.wue",  "last.1yr.etue"
                                     ,   "last.1yr.rue",   "last.1yr.cue",  "last.1yr.ecue"
                                     ,      "can.depth",       "can.area",            "wue"
                                     ,           "etue",            "rue",            "cue"
                                     ,           "ecue"
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
#             negative (e.g. npp), or that cannot be stacked (e.g. sm.stress)              #
#------------------------------------------------------------------------------------------#
   #----- All that we need here is the variable name. -------------------------------------#
   scen.barplot = data.frame( vname = c(         "agb",         "lai",          "ba"
                                       ,         "gpp",       "ldrop",         "mco"
                                       ,      "transp",  "plant.resp",      "bseeds"
                                       ,    "bstorage"
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
   scen.xyz      = list()
   scen.xyz$xvar = list( list( vname = "f.bseeds"        , leg        = "centre" )
                       , list( vname = "f.bstorage"      , leg        = "centre" )
                       , list( vname = "f.cba"           , leg        = "centre" )
                       , list( vname = "atm.temp"        , leg        = "centre" )
                       , list( vname = "atm.vpd"         , leg        = "centre" )
                       , list( vname = "last.1yr.change" , leg        = "centre" )
                       , list( vname = "last.1yr.et"     , leg        = "centre" )
                       , list( vname = "last.1yr.evap"   , leg        = "centre" )
                       , list( vname = "last.1yr.sms"    , leg        = "centre" )
                       , list( vname = "last.1yr.gpp"    , leg        = "centre" )
                       , list( vname = "last.1yr.growth" , leg        = "centre" )
                       , list( vname = "last.1yr.lgsw"   , leg        = "centre" )
                       , list( vname = "last.1yr.ltemp"  , leg        = "centre" )
                       , list( vname = "last.1yr.lwater" , leg        = "centre" )
                       , list( vname = "last.1yr.lvpd"   , leg        = "centre" )
                       , list( vname = "last.1yr.mwd"    , leg        = "centre" )
                       , list( vname = "last.1yr.ncbmort", leg        = "centre" )
                       , list( vname = "last.1yr.npp"    , leg        = "centre" )
                       , list( vname = "last.1yr.plresp" , leg        = "centre" )
                       , list( vname = "last.1yr.rshort" , leg        = "centre" )
                       , list( vname = "last.1yr.rain"   , leg        = "centre" )
                       , list( vname = "last.1yr.smpot"  , leg        = "centre" )
                       , list( vname = "last.1yr.transp" , leg        = "centre" )
                       , list( vname = "last.1yr.dcbadt" , leg        = "centre" )
                       , list( vname = "nmon.wdef"       , leg        = "centre" )
                       )#end list
   scen.xyz$yvar = scen.xyz$xvar
   scen.xyz$zvar = list( list( vname      = "agb"
                             , col.scheme = "clife"
                             , plog.xyz   = FALSE
                             )#end list
                       , list( vname      = "wood.dens"
                             , col.scheme = "iclife"
                             , plog.xyz   = TRUE
                             )#end list
                       , list( vname      = "last.1yr.growth"
                             , col.scheme = "clife"
                             , plog.xyz   = TRUE
                             )#end list
                       , list( vname      = "last.1yr.ncbmort"
                             , col.scheme = "iclife"
                             , plog.xyz   = TRUE
                             )#end list
                       , list( vname      = "last.1yr.change"
                             , col.scheme = "clife"
                             , plog.xyz   = FALSE
                             )#end list
                       , list( vname      = "last.1yr.cue"
                             , col.scheme = "clife"
                             , plog.xyz   = FALSE
                             )#end list
                       , list( vname      = "last.1yr.ecue"
                             , col.scheme = "clife"
                             , plog.xyz   = FALSE
                             )#end list
                       , list( vname      = "last.1yr.wue"
                             , col.scheme = "ivisible"
                             , plog.xyz   = FALSE
                             )#end list
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
   scen.comp = list( list( vname =     "wood.dens", low = "green"    , high = "purple"   )
                   , list( vname =           "agb", low = "purple"   , high = "green"    )
                   , list( vname =           "lai", low = "purple"   , high = "green"    )
                   , list( vname =            "ba", low = "purple"   , high = "green"    )
                   , list( vname =        "nplant", low = "purple"   , high = "green"    )
                   , list( vname =           "gpp", low = "purple"   , high = "green"    )
                   , list( vname =           "npp", low = "purple"   , high = "green"    )
                   , list( vname =           "mco", low = "purple"   , high = "green"    )
                   , list( vname =           "cba", low = "purple"   , high = "green"    )
                   , list( vname =        "dcbadt", low = "purple"   , high = "green"    )
                   , list( vname =         "ldrop", low = "green"    , high = "purple"   )
                   , list( vname =      "bstorage", low = "purple"   , high = "green"    )
                   , list( vname =        "bseeds", low = "purple"   , high = "green"    )
                   , list( vname =     "sm.stress", low = "blue"     , high = "orangered")
                   , list( vname =          "mort", low = "green"    , high = "purple"   )
                   , list( vname =          "recr", low = "purple"   , high = "green"    )
                   , list( vname =        "growth", low = "purple"   , high = "green"    )
                   , list( vname =       "ncbmort", low = "green"    , high = "purple"   )
                   , list( vname =        "dimort", low = "green"    , high = "purple"   )
                   , list( vname =        "change", low = "green"    , high = "purple"   )
                   , list( vname =      "agb.mort", low = "green"    , high = "purple"   )
                   , list( vname =      "agb.recr", low = "purple"   , high = "green"    )
                   , list( vname =    "agb.growth", low = "purple"   , high = "green"    )
                   , list( vname =   "agb.ncbmort", low = "green"    , high = "purple"   )
                   , list( vname =    "agb.dimort", low = "green"    , high = "purple"   )
                   , list( vname =    "agb.change", low = "green"    , high = "purple"   )
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
                   , list( vname =           "wue", low = "orangered", high = "blue"     )
                   , list( vname =           "cue", low = "purple"   , high = "green"    )
                   , list( vname =          "ecue", low = "purple"   , high = "green"    )
                   , list( vname =          "etue", low = "orangered", high = "blue"     )
                   , list( vname =      "leaf.gbw", low = "orangered", high = "blue"     )
                   , list( vname =      "leaf.gsw", low = "orangered", high = "blue"     )
                   , list( vname =      "leaf.par", low = "blue"     , high = "orangered")
                   , list( vname =   "leaf.rshort", low = "blue"     , high = "orangered")
                   , list( vname =    "leaf.rlong", low = "blue"     , high = "orangered")
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
                                    ,       "bseeds",    "sm.stress",     "phap.sms"
                                    ,         "mort",      "ncbmort",       "dimort"
                                    ,         "recr",       "growth",     "agb.mort"
                                    ,  "agb.ncbmort",   "agb.dimort",     "agb.recr"
                                    ,   "agb.growth",         "rain","water.deficit"
                                    ,       "runoff",     "atm.temp",    "leaf.temp"
                                    ,   "phap.ltemp",       "rshort",        "rlong"
                                    ,      "atm.vpd",     "leaf.vpd",    "phap.lvpd"
                                    ,        "smpot",          "nep",         "reco"
                                    ,  "fast.soil.c","struct.soil.c",  "slow.soil.c"
                                    ,     "het.resp",   "plant.resp",       "hflxlc"
                                    ,       "wflxlc",       "transp",       "hflxgc"
                                    ,       "hflxca",       "wflxgc",       "wflxca"
                                    ,          "wue",     "leaf.gbw",     "leaf.gsw"
                                    ,    "phap.lgsw",     "leaf.par",    "phap.lpar"
                                    ,  "leaf.rshort",       "nplant",          "cue"
                                    ,         "ecue",    "wood.dens",    "can.depth"
                                    ,     "can.area",       "dcbadt"
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
#     XYZ plots, to explore the parameter space.  Units, description, and log scale will   #
# be copied from scen.ts.                                                                  #
#                                                                                          #
#   IMPORTANT: All variables here MUST come from one of the variables defined in scen.ts!! #
#------------------------------------------------------------------------------------------#
   #----- All that we need here is the variable name, legend position or colour scheme. ---#
   panel.xyz = list()
   panel.xyz$xvar = list( list( vname = "f.bstorage"      , leg        = "centre" )
                        , list( vname = "f.dcbadt"        , leg        = "centre" )
                        , list( vname = "last.1yr.change" , leg        = "centre" )
                        , list( vname = "last.1yr.cue"    , leg        = "centre" )
                        , list( vname = "last.1yr.dcbadt" , leg        = "centre" )
                        , list( vname = "last.1yr.ecue"   , leg        = "centre" )
                        , list( vname = "last.1yr.et"     , leg        = "centre" )
                        , list( vname = "last.1yr.etue"   , leg        = "centre" )
                        , list( vname = "last.1yr.evap"   , leg        = "centre" )
                        , list( vname = "last.1yr.sms"    , leg        = "centre" )
                        , list( vname = "last.1yr.gpp"    , leg        = "centre" )
                        , list( vname = "last.1yr.growth" , leg        = "centre" )
                        , list( vname = "last.1yr.lpar"   , leg        = "centre" )
                        , list( vname = "last.1yr.lgsw"   , leg        = "centre" )
                        , list( vname = "last.1yr.ltemp"  , leg        = "centre" )
                        , list( vname = "last.1yr.lwater" , leg        = "centre" )
                        , list( vname = "last.1yr.lvpd"   , leg        = "centre" )
                        , list( vname = "last.1yr.mwd"    , leg        = "centre" )
                        , list( vname = "last.1yr.ncbmort", leg        = "centre" )
                        , list( vname = "last.1yr.npp"    , leg        = "centre" )
                        , list( vname = "last.1yr.plresp" , leg        = "centre" )
                        , list( vname = "last.1yr.rain"   , leg        = "centre" )
                        , list( vname = "last.1yr.rue"    , leg        = "centre" )
                        , list( vname = "last.1yr.smpot"  , leg        = "centre" )
                        , list( vname = "last.1yr.transp" , leg        = "centre" )
                        , list( vname = "last.1yr.wue"    , leg        = "centre" )
                        , list( vname = "nmon.wdef"       , leg        = "centre" )
                        )#end list
   panel.xyz$yvar = panel.xyz$xvar
   panel.xyz$zvar = list( list( vname      = "agb"
                              , col.scheme = "clife"
                              , plog.xyz   = FALSE
                              )#end list
                        , list( vname      = "f.bstorage"
                              , col.scheme = "clife"
                              , plog.xyz   = FALSE
                              )#end list
                        , list( vname      = "last.1yr.change"
                              , col.scheme = "clife"
                              , plog.xyz   = FALSE
                              )#end list
                        , list( vname      = "last.1yr.cue"
                              , col.scheme = "clife"
                              , plog.xyz   = FALSE
                              )#end list
                        , list( vname      = "last.1yr.dcbadt"
                              , col.scheme = "clife"
                              , plog.xyz   = FALSE
                              )#end list
                        , list( vname      = "last.1yr.ecue"
                              , col.scheme = "clife"
                              , plog.xyz   = FALSE
                              )#end list
                        , list( vname      = "last.1yr.etue"
                              , col.scheme = "ivisible"
                              , plog.xyz   = FALSE
                              )#end list
                        , list( vname      = "last.1yr.gpp"
                              , col.scheme = "clife"
                              , plog.xyz   = FALSE
                              )#end list
                        , list( vname      = "last.1yr.growth"
                              , col.scheme = "clife"
                              , plog.xyz   = TRUE
                              )#end list
                        , list( vname      = "last.1yr.ncbmort"
                              , col.scheme = "iclife"
                              , plog.xyz   = TRUE
                              )#end list
                        , list( vname      = "last.1yr.npp"
                              , col.scheme = "clife"
                              , plog.xyz   = FALSE
                              )#end list
                        , list( vname      = "last.1yr.plresp"
                              , col.scheme = "iclife"
                              , plog.xyz   = FALSE
                              )#end list
                        , list( vname      = "last.1yr.rue"
                              , col.scheme = "ivisible"
                              , plog.xyz   = FALSE
                              )#end list
                        , list( vname      = "last.1yr.transp"
                              , col.scheme = "ivisible"
                              , plog.xyz   = FALSE
                              )#end list
                        , list( vname      = "last.1yr.wue"
                              , col.scheme = "ivisible"
                              , plog.xyz   = FALSE
                              )#end list
                        , list( vname      = "wood.dens"
                              , col.scheme = "iclife"
                              , plog.xyz   = FALSE
                              )#end list
                        )#end list
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Replace the list by a data frame.                                                 #
   #---------------------------------------------------------------------------------------#
   panel.xyz$xvar = data.frame( apply(X=sapply(X=panel.xyz$xvar,FUN=c),MARGIN=1,FUN=unlist)
                             , stringsAsFactors = FALSE
                             )#end data.frame
   panel.xyz$yvar = data.frame( apply(X=sapply(X=panel.xyz$yvar,FUN=c),MARGIN=1,FUN=unlist)
                             , stringsAsFactors = FALSE
                             )#end data.frame
   panel.xyz$zvar = data.frame( apply(X=sapply(X=panel.xyz$zvar,FUN=c),MARGIN=1,FUN=unlist)
                             , stringsAsFactors = FALSE
                             )#end data.frame
   #---------------------------------------------------------------------------------------#






   #---------------------------------------------------------------------------------------#
   #     Fill in the XYZ list with information brought from scen.ts.                       #
   #---------------------------------------------------------------------------------------#
   #----- Get all the names that shall be added. ------------------------------------------#
   which.names = names(scen.ts) 
   keep        = ( ! which.names %in% union( union( names(panel.xyz$xvar)
                                                  , names(panel.xyz$yvar)
                                                  )#end union
                                           , names(panel.xyz$zvar)
                                           )#end union
                 )#end keep
   which.names = which.names[keep]
   #----- Match the lists based on vname. -------------------------------------------------#
   x.idx = match(panel.xyz$xvar$vname,scen.ts$vname)
   y.idx = match(panel.xyz$yvar$vname,scen.ts$vname)
   z.idx = match(panel.xyz$zvar$vname,scen.ts$vname)
   #----- Look for variables that weren't defined in scen.ts. -----------------------------#
   x.sel = is.na(x.idx)
   y.sel = is.na(y.idx)
   z.sel = is.na(z.idx)
   if (any(x.sel) | any(y.sel) | any(z.sel)){
      if (any(x.sel)){
         cat(" - The following variables in panel.xyz$xvar are missing from scen.ts:","\n")
         cat(paste("   * ",panel.xyz$xvar$vname[x.sel],sep=""),sep="\n")
      }#end if(any(x.sel))
      if (any(y.sel)){
         cat(" - The following variables in panel.xyz$yvar are missing from scen.ts:","\n")
         cat(paste("   * ",panel.xyz$yvar$vname[y.sel],sep=""),sep="\n")
      }#end if(any(y.sel))
      if (any(z.sel)){
         cat(" - The following variables in panel.xyz$zvar are missing from scen.ts:","\n")
         cat(paste("   * ",panel.xyz$zvar$vname[z.sel],sep=""),sep="\n")
      }#end if(any(z.sel))
      stop(" - All variables defined in panel.xyz must be defined in scen.ts!!!")
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Append the information to the data frame. ---------------------------------------#
   for (wn in which.names){
      panel.xyz$xvar[[wn]] = scen.ts[[wn]][x.idx]
      panel.xyz$yvar[[wn]] = scen.ts[[wn]][y.idx]
      panel.xyz$zvar[[wn]] = scen.ts[[wn]][z.idx]
   }#end for
   for (nl in c("plog.xyz")){
      panel.xyz$zvar[[nl]] = as.logical(panel.xyz$zvar[[nl]])
   }#end for
   #---------------------------------------------------------------------------------------#
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This list tells which variables to do the principal component analysis by DBH class. #
#                                                                                          #
# IMPORTANT:  All variables here MUST come from one of the variables defined in scen.ts.   #
#------------------------------------------------------------------------------------------#
   #----- All that we need here is the variable name. -------------------------------------#
   pca.explain = list( list( vname = "water.deficit" , colour = "#FF80FF"              )
                     , list( vname =      "leaf.gsw" , colour = "#C800C8"              )
                     , list( vname =    "f.bstorage" , colour = "#800080"              )
                     , list( vname =      "leaf.par" , colour = "#FF0000"              )
                     , list( vname =     "leaf.temp" , colour = "#BE0000"              )
                     , list( vname =      "leaf.vpd" , colour = "#800000"              )
                     , list( vname =    "agb.growth" , colour = "#F0F000"              )
                     , list( vname =     "wood.dens" , colour = "#A08240"              )
                     , list( vname =      "f.dcbadt" , colour = "#604020"              )
                     , list( vname =      "agb.recr" , colour = "#40FF00"              )
                     , list( vname =          "ecue" , colour = "#009000"              )
                     , list( vname =           "agb" , colour = "#004000"              )
                     , list( vname =           "wue" , colour = "#26D4EE"              )
                     , list( vname =    "leaf.water" , colour = "#1A9BE3"              )
                     , list( vname =          "rain" , colour = "#0C54AA"              )
                     , list( vname =     "nmon.wdef" , colour = "#9B9BFF"              )
                     , list( vname =         "smpot" , colour = "#694AFF"              )
                     , list( vname =      "agb.mort" , colour = "#5C00B8"              )
                     )#end list
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Replace the list by a data frame.                                                 #
   #---------------------------------------------------------------------------------------#
   pca.explain    = data.frame( apply(X=sapply(X=pca.explain,FUN=c),MARGIN=1,FUN=unlist)
                              , stringsAsFactors = FALSE
                              )#end data.frame
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Fill in the box plot list with information brought from scen.ts.                  #
   #---------------------------------------------------------------------------------------#
   #----- Find the variable names to be added. --------------------------------------------#
   which.names = names(scen.ts)
   keep        = ! ( which.names %in% names(pca.explain))
   which.names = which.names[keep]
   #----- Match the lists based on vname. -------------------------------------------------#
   comp.idx = match(pca.explain$vname,scen.ts$vname)
   #----- Look for variables that weren't defined in scen.ts. -----------------------------#
   comp.miss = is.na(comp.idx)
   if (any(comp.miss)){
      cat(" - The following variables in pca.explain are missing from scen.ts:","\n")
      cat(paste("   * ",pca.explain$vname[comp.miss],sep=""),sep="\n")
      stop(" - All variables defined in pca.explain must be defined in scen.ts!!!")
   }#end if(any(x.sel))
   #---------------------------------------------------------------------------------------#



   #----- Append the information to the data frame. ---------------------------------------#
   for (wn in which.names) pca.explain[[wn]] = scen.ts[[wn]][comp.idx]
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
pca.explain   <<- pca.explain
nscen.ts      <<- nrow(scen.ts)
nscen.szpft   <<- nrow(scen.szpft    )
nscen.barplot <<- nrow(scen.barplot  )
nscen.xvar    <<- nrow(scen.xyz$xvar )
nscen.yvar    <<- nrow(scen.xyz$yvar )
nscen.zvar    <<- nrow(scen.xyz$zvar )
nscen.comp    <<- nrow(scen.comp     )
npanel.box    <<- nrow(panel.box     )
npanel.xvar   <<- nrow(panel.xyz$xvar)
npanel.yvar   <<- nrow(panel.xyz$yvar)
npanel.zvar   <<- nrow(panel.xyz$zvar)
npca.explain  <<- nrow(pca.explain   )
#==========================================================================================#
#==========================================================================================#
