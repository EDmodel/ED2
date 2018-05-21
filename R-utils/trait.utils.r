trop.southam     <<- c("Bolivia","Brazil","Colombia","Ecuador","French Guiana","Guyana"
                      ,"Paraguay","Peru","Suriname","Venezuela")
trop.mesoam      <<- c("Belize","Costa Rica","El Salvador","Guatemala","Honduras"
                      ,"Nicaragua","Panama")
trop.caribbean   <<- c("Aruba","Anguilla","Antigua and Barbuda","Bahamas","Barbados"
                      ,"Bonaire","British Virgin Islands","Cayman Islands","Cuba","Curacao"
                      ,"Dominica","Dominican Republic","Grenada","Guadeloupe","Haiti"
                      ,"Jamaica","Martinique","Montserrat","Puerto Rico","Saba"
                      ,"Saint Barthelemy","Saint Eustatius","Saint Kitts and Nevis"
                      ,"Saint Lucia","Saint-Martin","Sint Maarten"
                      ,"Saint Vincent and the Grenadines","Trinidad and Tobago"
                      ,"United States Virgin Islands")
trop.northam     <<- c("Mexico")
trop.atlantic    <<- c("Ascencion Island","Cape Verde","Saint Helena")
trop.africa      <<- c("Angola","Benin","Botswana","Burkina-Faso","Burundi","Cameroon"
                      ,"Central African Republic","Chad","Congo","Cote d.Ivoire"
                      ,"Democratic Republic of the Congo","Djibouti","Equatorial Guinea"
                      ,"Eritrea","Ethiopia","Gabon","Ghana","Guinea","Guinea-Bissau"
                      ,"Kenya","Liberia","Malawi","Mali","Mauritania","Mozambique"
                      ,"Namibia","Niger","Nigeria","Rwanda","Senegal","Sierra Leone"
                      ,"Somalia","South Sudan","Sudan","Tanzania","The Gambia","Togo"
                      ,"Uganda","Zambia","Zimbabwe")
trop.indian      <<- c("British Indian Ocean Territory","Comoros","Madagascar"
                      ,"Maldives","Mauritius","Mayotte","Reunion","Rodrigues","Seychelles")
trop.asia        <<- c("Bahrain","Brunei","Cambodia","Cocos Islands","Christmas Island"
                      ,"East Timor","India","Indonesia","Laos","Malaysia","Myanmar"
                      ,"Oman","Qatar","Philippines","Saudi Arabia","Singapore","Sri Lanka"
                      ,"Thailand","Timor-Leste","United Arab Emirates","Vietnam","Viet Nam"
                      ,"Yemen")
trop.oceania     <<- c("American Samoa","Australia","Cook Islands"
                      ,"Federated States of Micronesia","Fiji","French Polynesia"
                      ,"Galapagos","Guam","Hawaii","Kiribati","Marshall Islands","Nauru"
                      ,"New Caledonia","Niue","Northern Mariana Islands","Palau"
                      ,"Papua New Guinea","Rotuma","Samoa","Solomon Islands","Tonga"
                      ,"Tuvalu","Vanuatu","Wallis and Futuna")
trop.countries   <<- c(trop.southam,trop.mesoam,trop.caribbean,trop.northam,trop.atlantic
                      ,trop.africa,trop.indian,trop.asia,trop.oceania)


#==========================================================================================#
#==========================================================================================#
#    Fix longitude from the TRY data base.                                                 #
#------------------------------------------------------------------------------------------#
fix.try.lon <<- function(x){
   ans = rep(NA_real_,times=length(x))

   #----- First step, replace bogus "0" entries. ------------------------------------------#
   miss    = x %in% "0"
   x[miss] = NA_character_
   #---------------------------------------------------------------------------------------#


   #----- Make western hemisphere negative. -----------------------------------------------#
   west    = grepl(pattern="W",x=x) | grepl(pattern="O",x=x)
   x[west] = paste0("-",x[west])
   #---------------------------------------------------------------------------------------#


   #----- Remove letters to denote longitude. ---------------------------------------------#
   x = gsub(pattern="W",replacement="",x=x)
   x = gsub(pattern="O",replacement="",x=x)
   x = gsub(pattern="E",replacement="",x=x)
   x = gsub(pattern="L",replacement="",x=x)
   #---------------------------------------------------------------------------------------#


   #------ Decide whether these entries are decimal or DMS. -------------------------------#
   with.x        = grepl(pattern="x",x=x)
   #---------------------------------------------------------------------------------------#


   #------ DMS, transform entries in decimal. ---------------------------------------------#
   xdms        = x[with.x]
   xdms        = gsub(pattern="\\." ,replacement="x",x=xdms)
   xdms        = gsub(pattern=" "   ,replacement="x",x=xdms)
   xdms        = gsub(pattern="x$"  ,replacement="" ,x=xdms)
   xdms        = gsub(pattern="xxxx",replacement="x",x=xdms)
   xdms        = gsub(pattern="xxx" ,replacement="x",x=xdms)
   xdms        = gsub(pattern="xx"  ,replacement="x",x=xdms)
   xdms.list   = strsplit(x=xdms,split="x")
   ans[with.x] = sapply(X=xdms.list,FUN=fix.try.dms.dec)
   #---------------------------------------------------------------------------------------#


   #------ Decimal, simply convert to numeric. --------------------------------------------#
   ans  = as.numeric(x)
   #---------------------------------------------------------------------------------------#

   #------ Make sure the coordinates are between 180W and 180E. ---------------------------#
   subtr      = (ans %>=% 180.) & (ans %<=% 360.)
   ans[subtr] = ans[subtr] - 360.
   #---------------------------------------------------------------------------------------#


   #------ Discard data outside range. ----------------------------------------------------#
   weird      = ! (ans %wr% c(-180.,180.))
   ans[weird] = NA
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end fix.try.lon
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#    Fix latitude from the TRY data base.                                                  #
#------------------------------------------------------------------------------------------#
fix.try.lat <<- function(x){
   ans = rep(NA_real_,times=length(x))

   #----- First step, replace bogus "0" entries. ------------------------------------------#
   miss   = x %in% "0"
   x[miss] = NA_character_
   #---------------------------------------------------------------------------------------#


   #----- Make southern hemisphere negative. ----------------------------------------------#
   south    = grepl(pattern="S",x=x) 
   x[south] = paste0("-",x[south])
   #---------------------------------------------------------------------------------------#


   #----- Remove letters to denote latitude. ----------------------------------------------#
   x = gsub(pattern="S",replacement="",x=x)
   x = gsub(pattern="N",replacement="",x=x)
   #---------------------------------------------------------------------------------------#


   #------ Decide whether these entries are decimal or DMS. -------------------------------#
   with.x        = grepl(pattern="x",x=x)
   #---------------------------------------------------------------------------------------#


   #------ DMS, transform entries in decimal. ---------------------------------------------#
   xdms        = x[with.x]
   xdms        = gsub(pattern="\\." ,replacement="x",x=xdms)
   xdms        = gsub(pattern=" "   ,replacement="x",x=xdms)
   xdms        = gsub(pattern="x$"  ,replacement="" ,x=xdms)
   xdms        = gsub(pattern="xxxx",replacement="x",x=xdms)
   xdms        = gsub(pattern="xxx" ,replacement="x",x=xdms)
   xdms        = gsub(pattern="xx"  ,replacement="x",x=xdms)
   xdms.list   = strsplit(x=xdms,split="x")
   ans[with.x] = sapply(X=xdms.list,FUN=fix.try.dms.dec)
   #---------------------------------------------------------------------------------------#

   #------ Decimal, simply convert to numeric. --------------------------------------------#
   ans = as.numeric(x)
   #---------------------------------------------------------------------------------------#


   #------ Discard data outside range. ----------------------------------------------------#
   weird      = ! (ans %wr% c(-90.,90.))
   ans[weird] = NA
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end fix.try.lat
#==========================================================================================#
#==========================================================================================#




#==========================================================================================#
#==========================================================================================#
#      DMS to decimal (TRY Data Base).                                                     #
#------------------------------------------------------------------------------------------#
fix.try.dms.dec <<- function(x){
   x = as.numeric(x)
   if (length(x) >= 2) x[2] = sign(x[1])*x[2]/60.
   if (length(x) >= 3) x[3] = sign(x[1])*x[3]/3600.
   if (length(x) >= 4) x[4] = sign(x[1])*x[4]/36000.
   if (length(x) >= 5) x = x[1:4]
   ans = sum(x)
   return(ans)
}#end fix.try.dms.dec
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#       Standardise countries, fix typos.                                                  #
#------------------------------------------------------------------------------------------#
fix.try.country <<- function(x){

   #----- Make sure all countries are spelt with capital first letter. --------------------#
   ans= ifelse(test = x %in% "", yes = NA_character_, no = capwords(x))
   #---------------------------------------------------------------------------------------#


   #----- Replace bogus countries. --------------------------------------------------------#
   ans[ans %in% "Alberta"               ] = "Canada"
   ans[ans %in% "Alpen"                 ] = "Czech Republic"
   ans[ans %in% "Alps"                  ] = "Czech Republic"
   ans[ans %in% "Argentian"             ] = "Argentina"
   ans[ans %in% "Arizona"               ] = "United States"
   ans[ans %in% "Austira"               ] = "Austria"
   ans[ans %in% "British Columbia"      ] = "Canada"
   ans[ans %in% "California"            ] = "United States"
   ans[ans %in% "Central Europa"        ] = "Czech Republic"
   ans[ans %in% "Central Europe"        ] = "Czech Republic"
   ans[ans %in% "Connecticut"           ] = "United States"
   ans[ans %in% "Czechia"               ] = "Czech Republic"
   ans[ans %in% "Demark"                ] = "Denmark"
   ans[ans %in% "Entral Europe"         ] = "Czech Republic"
   ans[ans %in% "ESP"                   ] = "Spain"
   ans[ans %in% "Ex Jugoslavia"         ] = "Macedonia"
   ans[ans %in% "Far East Of Russia"    ] = "Russia"
   ans[ans %in% "Florida"               ] = "United States"
   ans[ans %in% "Freance"               ] = "France"
   ans[ans %in% "Gemany"                ] = "Germany"
   ans[ans %in% "Germeny"               ] = "Germany"
   ans[ans %in% "Great Britain"         ] = "United Kingdom"
   ans[ans %in% "Great Britian"         ] = "United Kingdom"
   ans[ans %in% "Greatm Britain"        ] = "United Kingdom"
   ans[ans %in% "Greatn Britain"        ] = "United Kingdom"
   ans[ans %in% "Grenland"              ] = "Greenland"
   ans[ans %in% "Gret Britain"          ] = "United Kingdom"
   ans[ans %in% "Groenland"             ] = "Greenland"
   ans[ans %in% "Gronland"              ] = "Greenland"
   ans[ans %in% "Grxnland"              ] = "Greenland"
   ans[ans %in% "Holland"               ] = "Netherlands"
   ans[ans %in% "Icland"                ] = "Iceland"
   ans[ans %in% "Kansas"                ] = "United States"
   ans[ans %in% "Lituania"              ] = "Lithuania"
   ans[ans %in% "Macaronesia Portugal"  ] = "Macaronesia"
   ans[ans %in% "Macaronesia, Portugal" ] = "Macaronesia"
   ans[ans %in% "Macaronesia, Spain"    ] = "Macaronesia"
   ans[ans %in% "Macaronesia, Tenerife" ] = "Macaronesia"
   ans[ans %in% "Macaronesia,Portugal"  ] = "Macaronesia"
   ans[ans %in% "Macaronesia,Spain"     ] = "Macaronesia"
   ans[ans %in% "Marocco"               ] = "Morocco"
   ans[ans %in% "Massachusetts"         ] = "United States"
   ans[ans %in% "Michigan"              ] = "United States"
   ans[ans %in% "Montana"               ] = "United States"
   ans[ans %in% "Nebraska"              ] = "United States"
   ans[ans %in% "Netharland"            ] = "Netherlands"
   ans[ans %in% "Netherland"            ] = "Netherlands"
   ans[ans %in% "New York"              ] = "United States"
   ans[ans %in% "North America"         ] = "United States"
   ans[ans %in% "North Carolina"        ] = "United States"
   ans[ans %in% "Northwest Costa Rica:" ] = "Costa Rica"
   ans[ans %in% "Norway Svalbard"       ] = "Svalbard"
   ans[ans %in% "Oregon"                ] = "United States"
   ans[ans %in% "Pennsylvania"          ] = "United States"
   ans[ans %in% "PRT"                   ] = "Portugal"
   ans[ans %in% "Rusiia"                ] = "Russia"
   ans[ans %in% "Schwitzerland"         ] = "Switzerland"
   ans[ans %in% "Scotland"              ] = "United Kingdom"
   ans[ans %in% "South Sweden"          ] = "Sweden"
   ans[ans %in% "Sveden"                ] = "Sweden"
   ans[ans %in% "Switzerland, Alps"     ] = "Switzerland"
   ans[ans %in% "Tennessee"             ] = "United States"
   ans[ans %in% "Texas"                 ] = "United States"
   ans[ans %in% "The Netherland"        ] = "Netherlands"
   ans[ans %in% "The Netherlands"       ] = "Netherlands"
   ans[ans %in% "Ukraina"               ] = "Ukraine"
   ans[ans %in% "USA"                   ] = "United States"
   ans[ans %in% "Utah"                  ] = "United States"
   ans[ans %in% "Virginia"              ] = "United States"
   ans[ans %in% "Washington"            ] = "United States"
   ans[ans %in% "Wisconsin"             ] = "United States"
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end fix.try.country
#==========================================================================================#
#==========================================================================================#




#------------------------------------------------------------------------------------------#
#       Standardise countries, fix typos.                                                  #
#------------------------------------------------------------------------------------------#
fix.try.continent <<- function(x){

   #----- Make sure all countries are spelt with capital first letter. --------------------#
   ans= ifelse( test = x %in% "",yes=NA_character_,no=capwords(x))
   #---------------------------------------------------------------------------------------#


   #----- Replace bogus countries. --------------------------------------------------------#
   ans[ans %in% "Australasia"           ] = "Oceania"
   ans[ans %in% "Australia"             ] = "Oceania"
   ans[ans %in% "Austria"               ] = "Europe"
   ans[ans %in% "Central Asia"          ] = "Asia"
   ans[ans %in% "Europa; Russia;Georgia"] = "Europe"
   ans[ans %in% "Hawaii"                ] = "Oceania"
   ans[ans %in% "N=America"             ] = "North America"
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end fix.try.continent
#==========================================================================================#
#==========================================================================================#




#==========================================================================================#
#==========================================================================================#
#      Separate location information from Chave's data base into habitat, country and      #
# continent.                                                                               #
#------------------------------------------------------------------------------------------#
fix.chave.siteinfo <<- function(dat){

   #----- Assign habitat. -----------------------------------------------------------------#
   dat$habitat = rep(NA_character_,nrow(dat))
   dat$habitat  [dat$region %in% "Africa (extratropical)"       ] = "temperate"
   dat$habitat  [dat$region %in% "Africa (tropical)"            ] = "tropical"
   dat$habitat  [dat$region %in% "Australia"                    ] = "temperate"
   dat$habitat  [dat$region %in% "Australia/PNG (tropical)"     ] = "tropical"
   dat$habitat  [dat$region %in% "Central America (tropical)"   ] = "tropical"
   dat$habitat  [dat$region %in% "China"                        ] = "temperate"
   dat$habitat  [dat$region %in% "Europe"                       ] = "temperate"
   dat$habitat  [dat$region %in% "India"                        ] = "subtropical"
   dat$habitat  [dat$region %in% "Madagascar"                   ] = "tropical"
   dat$habitat  [dat$region %in% "Mexico"                       ] = "subtropical"
   dat$habitat  [dat$region %in% "NorthAmerica"                 ] = "temperate"
   dat$habitat  [dat$region %in% "Oceania"                      ] = "tropical"
   dat$habitat  [dat$region %in% "South America (extratropical)"] = "temperate"
   dat$habitat  [dat$region %in% "South America (tropical)"     ] = "tropical"
   dat$habitat  [dat$region %in% "South America (Tropical)"     ] = "tropical"
   dat$habitat  [dat$region %in% "South-East Asia"              ] = "tropical"
   dat$habitat  [dat$region %in% "South-East Asia (tropical)"   ] = "tropical"
   #---------------------------------------------------------------------------------------#

   #----- Assign continent. ---------------------------------------------------------------#
   dat$continent = rep(NA_character_,nrow(dat))
   dat$continent[dat$region %in% "Africa (extratropical)"       ] = "Africa"
   dat$continent[dat$region %in% "Africa (tropical)"            ] = "Africa"
   dat$continent[dat$region %in% "Australia"                    ] = "Oceania"
   dat$continent[dat$region %in% "Australia/PNG (tropical)"     ] = "Oceania"
   dat$continent[dat$region %in% "Central America (tropical)"   ] = "Central America"
   dat$continent[dat$region %in% "China"                        ] = "Asia"
   dat$continent[dat$region %in% "Europe"                       ] = "Europe"
   dat$continent[dat$region %in% "India"                        ] = "Asia"
   dat$continent[dat$region %in% "Madagascar"                   ] = "Africa"
   dat$continent[dat$region %in% "Mexico"                       ] = "North America"
   dat$continent[dat$region %in% "NorthAmerica"                 ] = "North America"
   dat$continent[dat$region %in% "Oceania"                      ] = "Oceania"
   dat$continent[dat$region %in% "South America (extratropical)"] = "South America"
   dat$continent[dat$region %in% "South America (tropical)"     ] = "South America"
   dat$continent[dat$region %in% "South America (Tropical)"     ] = "South America"
   dat$continent[dat$region %in% "South-East Asia"              ] = "Asia"
   dat$continent[dat$region %in% "South-East Asia (tropical)"   ] = "Asia"
   #---------------------------------------------------------------------------------------#

   #----- Assign country. -----------------------------------------------------------------#
   dat$country = rep(NA_character_,nrow(dat))
   dat$country  [dat$region %in% "Australia"                    ] = "Australia"
   dat$country  [dat$region %in% "China"                        ] = "China"
   dat$country  [dat$region %in% "India"                        ] = "India"
   dat$country  [dat$region %in% "Madagascar"                   ] = "Madagascar"
   dat$country  [dat$region %in% "Mexico"                       ] = "Mexico"
   #---------------------------------------------------------------------------------------#


   return(dat)
}#end fix.chave.siteinfo
#==========================================================================================#
#==========================================================================================#




#==========================================================================================#
#==========================================================================================#
#      Separate location information from GLOPNET data base into habitat, country and      #
# continent.                                                                               #
#------------------------------------------------------------------------------------------#
fix.glopnet.siteinfo <<- function(dat){


   #------ Handy aliases. -----------------------------------------------------------------#
   d.s = dat$site
   d.b = dat$biome
   #---------------------------------------------------------------------------------------#


   #----- Assign continent. ---------------------------------------------------------------#
   dat$country = rep(NA_character_,nrow(dat))
   dat$country[d.s %in% "Ackerly_Jasper"                          ] = "United States"
   dat$country[d.s %in% "Baruch&Goldstein_Hawaii_High_Old"        ] = "Hawaii"
   dat$country[d.s %in% "Baruch&Goldstein_Hawaii_High_Rec"        ] = "Hawaii"
   dat$country[d.s %in% "Baruch&Goldstein_Hawaii_Low_Old"         ] = "Hawaii"
   dat$country[d.s %in% "Baruch&Goldstein_Hawaii_Med_Old"         ] = "Hawaii"
   dat$country[d.s %in% "Baruch&Goldstein_Hawaii_Med_Rec"         ] = "Hawaii"
   dat$country[d.s %in% "Bassow&Bazzaz_PETERSHAM__MA"             ] = "United States"
   dat$country[d.s %in% "Bongers_et_al_Los_Tuxtlas"               ] = "Mexico"
   dat$country[d.s %in% "Cavendar_Bares_Florida"                  ] = "United States"
   dat$country[d.s %in% "Chapin_etc_Toolik_Lake"                  ] = "United States"
   dat$country[d.s %in% "Christodoulakis_Malakasa"                ] = "Greece"
   dat$country[d.s %in% "Chua_et_al_Malaysia"                     ] = "Malaysia"
   dat$country[d.s %in% "Coley_BCI"                               ] = "Panama"
   dat$country[d.s %in% "Cornelissen_UK_Sheffield"                ] = "United Kingdom"
   dat$country[d.s %in% "DeLucia91Ecol_RENO__NEVADA"              ] = "United States"
   dat$country[d.s %in% "DeLucia95_Okefenokee_Swamp"              ] = "United States"
   dat$country[d.s %in% "Diemer_Korner_Austria_high"              ] = "Austria"
   dat$country[d.s %in% "Diemer_Korner_Austria_low"               ] = "Austria"
   dat$country[d.s %in% "Diemer-Ecuador_high"                     ] = "Ecuador"
   dat$country[d.s %in% "Diemer-Ecuador_highest"                  ] = "Ecuador"
   dat$country[d.s %in% "Diemer-Ecuador_low"                      ] = "Ecuador"
   dat$country[d.s %in% "Diemer-Ecuador_lowest"                   ] = "Ecuador"
   dat$country[d.s %in% "Field_et_al_83_Jasper_Ridge"             ] = "United States"
   dat$country[d.s %in% "Franco&Luttge_Brasilia"                  ] = "Brazil"
   dat$country[d.s %in% "Garnier_etal_F/CR"                       ] = "France"
   dat$country[d.s %in% "Garnier_etal_Les_Agros"                  ] = "France"
   dat$country[d.s %in% "Garnier_etal_SM/C"                       ] = "France"
   dat$country[d.s %in% "Gulias_Binifaldo"                        ] = "Spain"
   dat$country[d.s %in% "Gulias_Puigpunyent"                      ] = "Spain"
   dat$country[d.s %in% "Gulias_Soller"                           ] = "Spain"
   dat$country[d.s %in% "Gulias_UIB"                              ] = "Spain"
   dat$country[d.s %in% "Hikosaka-Japan_Chiba_Japan"              ] = "Japan"
   dat$country[d.s %in% "Hikosaka-Malaysia_Mt_Kinabalu_high"      ] = "Malaysia"
   dat$country[d.s %in% "Hikosaka-Malaysia_Mt_Kinabalu_highest"   ] = "Malaysia"
   dat$country[d.s %in% "Hikosaka-Malaysia_Mt_Kinabalu_low"       ] = "Malaysia"
   dat$country[d.s %in% "Hikosaka-Malaysia_Mt_Kinabalu_med"       ] = "Malaysia"
   dat$country[d.s %in% "Hogan_etal_PNM_crane"                    ] = "Panama"
   dat$country[d.s %in% "Jayasekara_Sri_Lanka"                    ] = "Sri Lanka"
   dat$country[d.s %in% "Jose_Gillespie_Indiana"                  ] = "United States"
   dat$country[d.s %in% "Jurik86_Pellston_MI"                     ] = "United States"
   dat$country[d.s %in% "Kitajima_Panama"                         ] = "Panama"
   dat$country[d.s %in% "Koike_SAPPORO__JAPAN"                    ] = "Japan"
   dat$country[d.s %in% "Korner_et_al_86_Haast_Valley_NZ"         ] = "New Zealand"
   dat$country[d.s %in% "Kudo_Cornelissen_Abisko"                 ] = "Sweden"
   dat$country[d.s %in% "Kudo_Cornelissen_Latnjajaure"            ] = "Sweden"
   dat$country[d.s %in% "Kudo_Cornelissen_Svalbard"               ] = "Svalbard"
   dat$country[d.s %in% "Kudo96_high"                             ] = "Japan"
   dat$country[d.s %in% "Kudo96_low"                              ] = "Japan"
   dat$country[d.s %in% "Kuppers_Bayreuth"                        ] = "Germany"
   dat$country[d.s %in% "Lal_etal_Inceptisol"                     ] = "India"
   dat$country[d.s %in% "Lal_etal_Ultisol"                        ] = "India"
   dat$country[d.s %in% "Lamont_S_Africa_1_Citrusdal"             ] = "South Africa"
   dat$country[d.s %in% "Lamont_S_Africa_10_Kylemore"             ] = "South Africa"
   dat$country[d.s %in% "Lamont_S_Africa_11_Soetanysberg"         ] = "South Africa"
   dat$country[d.s %in% "Lamont_S_Africa_12_Salmonsdam"           ] = "South Africa"
   dat$country[d.s %in% "Lamont_S_Africa_13_Salmonsdam"           ] = "South Africa"
   dat$country[d.s %in% "Lamont_S_Africa_14_Herrmanus"            ] = "South Africa"
   dat$country[d.s %in% "Lamont_S_Africa_2_Stellenbosch"          ] = "South Africa"
   dat$country[d.s %in% "Lamont_S_Africa_3_Opdieberg"             ] = "South Africa"
   dat$country[d.s %in% "Lamont_S_Africa_4_Matjies_River"         ] = "South Africa"
   dat$country[d.s %in% "Lamont_S_Africa_5_Algeria"               ] = "South Africa"
   dat$country[d.s %in% "Lamont_S_Africa_6_Scarborough"           ] = "South Africa"
   dat$country[d.s %in% "Lamont_S_Africa_7_Hopefield "            ] = "South Africa"
   dat$country[d.s %in% "Lamont_S_Africa_9_Jonkershonk"           ] = "South Africa"
   dat$country[d.s %in% "Lamont_WA_Darling_Scarp"                 ] = "Australia"
   dat$country[d.s %in% "Lamont_WA_Eneabba"                       ] = "Australia"
   dat$country[d.s %in% "Lamont_WA_Esperance"                     ] = "Australia"
   dat$country[d.s %in% "Lamont_WA_Fitzgerald_River"              ] = "Australia"
   dat$country[d.s %in% "Lamont_WA_Kalbarri"                      ] = "Australia"
   dat$country[d.s %in% "Lamont_WA_Lake_King"                     ] = "Australia"
   dat$country[d.s %in% "Lamont_WA_Merridin_etc"                  ] = "Australia"
   dat$country[d.s %in% "Lamont_WA_Millbrook"                     ] = "Australia"
   dat$country[d.s %in% "Lamont_WA_Stirling_Ranges"               ] = "Australia"
   dat$country[d.s %in% "Lamont_WA_Walpole"                       ] = "Australia"
   dat$country[d.s %in% "Lamont_WA_Watheroo"                      ] = "Australia"
   dat$country[d.s %in% "Lee_Cedar_Creek2"                        ] = "United States"
   dat$country[d.s %in% "Lee_NZ_Murchison_Mtns"                   ] = "New Zealand"
   dat$country[d.s %in% "Lusk_saplings_Cordillera_Pelada"         ] = "Chile"
   dat$country[d.s %in% "Lusk-adults_Concepcion"                  ] = "Chile"
   dat$country[d.s %in% "Lusk-adults_Los Lleuques"                ] = "Chile"
   dat$country[d.s %in% "Lusk-adults_Puyehue"                     ] = "Chile"
   dat$country[d.s %in% "Marin_Medina_Piritu_Venezuela"           ] = "Venezuela"
   dat$country[d.s %in% "Martin_etal_Guanacaste"                  ] = "Costa Rica"
   dat$country[d.s %in% "McAllister_Konza"                        ] = "United States"
   dat$country[d.s %in% "Mediavilla_et_al_Salamanca"              ] = "Spain"
   dat$country[d.s %in% "MidgelySA_Alexandria"                    ] = "South Africa"
   dat$country[d.s %in% "MidgelySA_Amatolas"                      ] = "South Africa"
   dat$country[d.s %in% "MidgelySA_Dukuduku"                      ] = "South Africa"
   dat$country[d.s %in% "MidgelySA_Jonkershoek_Mtn"               ] = "South Africa"
   dat$country[d.s %in% "MidgelySA_Jonkershoek_Rip"               ] = "South Africa"
   dat$country[d.s %in% "MidgelySA_Knysna"                        ] = "South Africa"
   dat$country[d.s %in% "MidgelySA_Mapelane"                      ] = "South Africa"
   dat$country[d.s %in% "MidgelySA_Sand_Forest"                   ] = "South Africa"
   dat$country[d.s %in% "MidgelySA_Umtiza"                        ] = "South Africa"
   dat$country[d.s %in% "MitchellNC_Coweeta"                      ] = "United States"
   dat$country[d.s %in% "Miyazawa_Chiba_M"                        ] = "Japan"
   dat$country[d.s %in% "Mooney_etal_81_desert"                   ] = "United States"
   dat$country[d.s %in% "Mooney_etal_81_old-field"                ] = "United States"
   dat$country[d.s %in% "Mooney_etal83_Jonk_Mtn"                  ] = "South Africa"
   dat$country[d.s %in% "Mulkey9193_BCI_Panama"                   ] = "Panama"
   dat$country[d.s %in% "Nelson_etal_Texas"                       ] = "United States"
   dat$country[d.s %in% "Niinemets_Kull94_Estonia"                ] = "Estonia"
   dat$country[d.s %in% "Niinemets_Kull98_Tartu"                  ] = "Estonia"
   dat$country[d.s %in% "Nitta_Chiba_N"                           ] = "Japan"
   dat$country[d.s %in% "OleksynPol_Siemanice"                    ] = "Poland"
   dat$country[d.s %in% "Olivares_Caracas_Venezuela"              ] = "Venezuela"
   dat$country[d.s %in% "Osada_Thomas_Pasoh"                      ] = "Malaysia"
   dat$country[d.s %in% "Poorter_de_Jong_Along_ditch "            ] = "Netherlands"
   dat$country[d.s %in% "Poorter_de_Jong_Chalk_grassland"         ] = "Netherlands"
   dat$country[d.s %in% "Poorter_de_Jong_Dry_heath"               ] = "Netherlands"
   dat$country[d.s %in% "Poorter_de_Jong_Dry_open_grassland"      ] = "Netherlands"
   dat$country[d.s %in% "Poorter_de_Jong_Poor_hay_meadow"         ] = "Netherlands"
   dat$country[d.s %in% "Poorter_de_Jong_Quaking_Fen"             ] = "Netherlands"
   dat$country[d.s %in% "Poorter_de_Jong_Reed_marsh"              ] = "Netherlands"
   dat$country[d.s %in% "Poorter_de_Jong_Wet_heath"               ] = "Netherlands"
   dat$country[d.s %in% "Prado&DeMoraes1997_SAO_CARLOS"           ] = "Brazil"
   dat$country[d.s %in% "Prior_dry monsoon forest"                ] = "Australia"
   dat$country[d.s %in% "Prior_open forest"                       ] = "Australia"
   dat$country[d.s %in% "Prior_swamp"                             ] = "Australia"
   dat$country[d.s %in% "Prior_woodland"                          ] = "Australia"
   dat$country[d.s %in% "Pyankov_Tadjikistan_Tadjikistan_high"    ] = "Tajikistan"
   dat$country[d.s %in% "Pyankov_Tadjikistan_Tadjikistan_higher"  ] = "Tajikistan"
   dat$country[d.s %in% "Pyankov_Tadjikistan_Tadjikistan_highest" ] = "Tajikistan"
   dat$country[d.s %in% "Pyankov_Urals_Yekaterinburg"             ] = "Russia"
   dat$country[d.s %in% "Reichetal_Colorado"                      ] = "United States"
   dat$country[d.s %in% "Reichetal_N_Carolina"                    ] = "United States"
   dat$country[d.s %in% "Reichetal_New_Mexico"                    ] = "United States"
   dat$country[d.s %in% "Reichetal_S_Carolina"                    ] = "United States"
   dat$country[d.s %in% "Reichetal_Venezuela"                     ] = "Venezuela"
   dat$country[d.s %in% "Reichetal_Wisconsin"                     ] = "United States"
   dat$country[d.s %in% "Ricklefs_SE_Ontario"                     ] = "Canada"
   dat$country[d.s %in% "Schulze_Kapalga"                         ] = "Australia"
   dat$country[d.s %in% "Schulze_Katherine"                       ] = "Australia"
   dat$country[d.s %in% "Schulze_Kidman Springs"                  ] = "Australia"
   dat$country[d.s %in% "Schulze_Melville Island"                 ] = "Australia"
   dat$country[d.s %in% "Schulze_Mt_Sanford"                      ] = "Australia"
   dat$country[d.s %in% "Shipley_Sherbrooke"                      ] = "Canada"
   dat$country[d.s %in% "Small1972_OTTAWA"                        ] = "Canada"
   dat$country[d.s %in% "Sobrado_Charallave"                      ] = "Venezuela"
   dat$country[d.s %in% "Sobrado&Medina_SanCarlos_bana"           ] = "Venezuela"
   dat$country[d.s %in% "Specht_Rundel_Dark_Island_heath"         ] = "Australia"
   dat$country[d.s %in% "Specht_Rundel_Dark_Island_mallee"        ] = "Australia"
   dat$country[d.s %in% "Specht_Rundel_Mt_Lofty"                  ] = "Australia"
   dat$country[d.s %in% "Tan_et_al_adinandra_trema_belukar"       ] = "Singapore"
   dat$country[d.s %in% "Terashima_Nepal"                         ] = "Nepal"
   dat$country[d.s %in% "Tezara_etal98_Coro"                      ] = "Venezuela"
   dat$country[d.s %in% "Tjoelker_Cedar_Creek"                    ] = "United States"
   dat$country[d.s %in% "Tuohy_etal_Zimbabwe_CHID"                ] = "Zimbabwe"
   dat$country[d.s %in% "Tuohy_etal_Zimbabwe_CRST_MCLW"           ] = "Zimbabwe"
   dat$country[d.s %in% "Tuohy_etal_Zimbabwe_MAT"                 ] = "Zimbabwe"
   dat$country[d.s %in% "Turner_&_Tan_Adinandra_Belukar"          ] = "Singapore"
   dat$country[d.s %in% "Turner_&_Tan_Beach_forest"               ] = "Singapore"
   dat$country[d.s %in% "Turner_&_Tan_Mangroves"                  ] = "Singapore"
   dat$country[d.s %in% "Turner_&_Tan_Undegraded_secondary_forest"] = "Singapore"
   dat$country[d.s %in% "Veneklaas_W_Australia"                   ] = "Australia"
   dat$country[d.s %in% "Villar_Andalucia_mesic"                  ] = "Spain"
   dat$country[d.s %in% "Villar_Andalucia_xeric"                  ] = "Spain"
   dat$country[d.s %in% "Villar_California_chaparral"             ] = "United States"
   dat$country[d.s %in% "Villar_California_forest"                ] = "United States"
   dat$country[d.s %in% "Villar_Canary_Is_lauriphyll"             ] = "Canary Islands"
   dat$country[d.s %in% "Villar_Canary_Is_xeric"                  ] = "Canary Islands"
   dat$country[d.s %in% "Villar_Chihuahua"                        ] = "Mexico"
   dat$country[d.s %in% "Villar_Devon_Is_Canada"                  ] = "Canada"
   dat$country[d.s %in% "Villar_Douala-Edea Forest, Cameroon"     ] = "Cameroon"
   dat$country[d.s %in% "Villar_Kibale Forest, Uganda"            ] = "Uganda"
   dat$country[d.s %in% "Villar_N_Carolina_forest"                ] = "United States"
   dat$country[d.s %in% "Villar_Tierra_del_Fuego"                 ] = "Argentina"
   dat$country[d.s %in% "Villar_Toronto"                          ] = "Canada"
   dat$country[d.s %in% "Williams et al_LosTuxtlas2"              ] = "Mexico"
   dat$country[d.s %in% "Williams_Linera_Mexico"                  ] = "Mexico"
   dat$country[d.s %in% "Wright_Oz_syd_hiP"                       ] = "Australia"
   dat$country[d.s %in% "Wright_Oz_syd_loP"                       ] = "Australia"
   dat$country[d.s %in% "Wright_Oz_wnsw_hiP"                      ] = "Australia"
   dat$country[d.s %in% "Wright_Oz_wnsw_loP"                      ] = "Australia"
   dat$country[d.s %in% "Zotz_Fortuna_Panama"                     ] = "Panama"
   #---------------------------------------------------------------------------------------#

   #----- Assign country. -----------------------------------------------------------------#
   dat$continent = rep(NA_character_,nrow(dat))
   dat$continent[d.s %in% "Ackerly_Jasper"                          ] = "North America"
   dat$continent[d.s %in% "Baruch&Goldstein_Hawaii_High_Old"        ] = "Hawaii"
   dat$continent[d.s %in% "Baruch&Goldstein_Hawaii_High_Rec"        ] = "Hawaii"
   dat$continent[d.s %in% "Baruch&Goldstein_Hawaii_Low_Old"         ] = "Hawaii"
   dat$continent[d.s %in% "Baruch&Goldstein_Hawaii_Med_Old"         ] = "Hawaii"
   dat$continent[d.s %in% "Baruch&Goldstein_Hawaii_Med_Rec"         ] = "Hawaii"
   dat$continent[d.s %in% "Bassow&Bazzaz_PETERSHAM__MA"             ] = "North America"
   dat$continent[d.s %in% "Bongers_et_al_Los_Tuxtlas"               ] = "North America"
   dat$continent[d.s %in% "Cavendar_Bares_Florida"                  ] = "North America"
   dat$continent[d.s %in% "Chapin_etc_Toolik_Lake"                  ] = "North America"
   dat$continent[d.s %in% "Christodoulakis_Malakasa"                ] = "Europe"
   dat$continent[d.s %in% "Chua_et_al_Malaysia"                     ] = "Asia"
   dat$continent[d.s %in% "Coley_BCI"                               ] = "Central America"
   dat$continent[d.s %in% "Cornelissen_UK_Sheffield"                ] = "Europe"
   dat$continent[d.s %in% "DeLucia91Ecol_RENO__NEVADA"              ] = "North America"
   dat$continent[d.s %in% "DeLucia95_Okefenokee_Swamp"              ] = "North America"
   dat$continent[d.s %in% "Diemer_Korner_Austria_high"              ] = "Austria"
   dat$continent[d.s %in% "Diemer_Korner_Austria_low"               ] = "Austria"
   dat$continent[d.s %in% "Diemer-Ecuador_high"                     ] = "South America"
   dat$continent[d.s %in% "Diemer-Ecuador_highest"                  ] = "South America"
   dat$continent[d.s %in% "Diemer-Ecuador_low"                      ] = "South America"
   dat$continent[d.s %in% "Diemer-Ecuador_lowest"                   ] = "South America"
   dat$continent[d.s %in% "Field_et_al_83_Jasper_Ridge"             ] = "North America"
   dat$continent[d.s %in% "Franco&Luttge_Brasilia"                  ] = "South America"
   dat$continent[d.s %in% "Garnier_etal_F/CR"                       ] = "Europe"
   dat$continent[d.s %in% "Garnier_etal_Les_Agros"                  ] = "Europe"
   dat$continent[d.s %in% "Garnier_etal_SM/C"                       ] = "Europe"
   dat$continent[d.s %in% "Gulias_Binifaldo"                        ] = "Europe"
   dat$continent[d.s %in% "Gulias_Puigpunyent"                      ] = "Europe"
   dat$continent[d.s %in% "Gulias_Soller"                           ] = "Europe"
   dat$continent[d.s %in% "Gulias_UIB"                              ] = "Europe"
   dat$continent[d.s %in% "Hikosaka-Japan_Chiba_Japan"              ] = "Asia"
   dat$continent[d.s %in% "Hikosaka-Malaysia_Mt_Kinabalu_high"      ] = "Asia"
   dat$continent[d.s %in% "Hikosaka-Malaysia_Mt_Kinabalu_highest"   ] = "Asia"
   dat$continent[d.s %in% "Hikosaka-Malaysia_Mt_Kinabalu_low"       ] = "Asia"
   dat$continent[d.s %in% "Hikosaka-Malaysia_Mt_Kinabalu_med"       ] = "Asia"
   dat$continent[d.s %in% "Hogan_etal_PNM_crane"                    ] = "Central America"
   dat$continent[d.s %in% "Jayasekara_Sri_Lanka"                    ] = "Asia"
   dat$continent[d.s %in% "Jose_Gillespie_Indiana"                  ] = "North America"
   dat$continent[d.s %in% "Jurik86_Pellston_MI"                     ] = "North America"
   dat$continent[d.s %in% "Kitajima_Panama"                         ] = "Central America"
   dat$continent[d.s %in% "Koike_SAPPORO__JAPAN"                    ] = "Asia"
   dat$continent[d.s %in% "Korner_et_al_86_Haast_Valley_NZ"         ] = "Oceania"
   dat$continent[d.s %in% "Kudo_Cornelissen_Abisko"                 ] = "Europe"
   dat$continent[d.s %in% "Kudo_Cornelissen_Latnjajaure"            ] = "Europe"
   dat$continent[d.s %in% "Kudo_Cornelissen_Svalbard"               ] = "Europe"
   dat$continent[d.s %in% "Kudo96_high"                             ] = "Asia"
   dat$continent[d.s %in% "Kudo96_low"                              ] = "Asia"
   dat$continent[d.s %in% "Kuppers_Bayreuth"                        ] = "Europe"
   dat$continent[d.s %in% "Lal_etal_Inceptisol"                     ] = "Asia"
   dat$continent[d.s %in% "Lal_etal_Ultisol"                        ] = "Asia"
   dat$continent[d.s %in% "Lamont_S_Africa_1_Citrusdal"             ] = "Africa"
   dat$continent[d.s %in% "Lamont_S_Africa_10_Kylemore"             ] = "Africa"
   dat$continent[d.s %in% "Lamont_S_Africa_11_Soetanysberg"         ] = "Africa"
   dat$continent[d.s %in% "Lamont_S_Africa_12_Salmonsdam"           ] = "Africa"
   dat$continent[d.s %in% "Lamont_S_Africa_13_Salmonsdam"           ] = "Africa"
   dat$continent[d.s %in% "Lamont_S_Africa_14_Herrmanus"            ] = "Africa"
   dat$continent[d.s %in% "Lamont_S_Africa_2_Stellenbosch"          ] = "Africa"
   dat$continent[d.s %in% "Lamont_S_Africa_3_Opdieberg"             ] = "Africa"
   dat$continent[d.s %in% "Lamont_S_Africa_4_Matjies_River"         ] = "Africa"
   dat$continent[d.s %in% "Lamont_S_Africa_5_Algeria"               ] = "Africa"
   dat$continent[d.s %in% "Lamont_S_Africa_6_Scarborough"           ] = "Africa"
   dat$continent[d.s %in% "Lamont_S_Africa_7_Hopefield "            ] = "Africa"
   dat$continent[d.s %in% "Lamont_S_Africa_9_Jonkershonk"           ] = "Africa"
   dat$continent[d.s %in% "Lamont_WA_Darling_Scarp"                 ] = "Oceania"
   dat$continent[d.s %in% "Lamont_WA_Eneabba"                       ] = "Oceania"
   dat$continent[d.s %in% "Lamont_WA_Esperance"                     ] = "Oceania"
   dat$continent[d.s %in% "Lamont_WA_Fitzgerald_River"              ] = "Oceania"
   dat$continent[d.s %in% "Lamont_WA_Kalbarri"                      ] = "Oceania"
   dat$continent[d.s %in% "Lamont_WA_Lake_King"                     ] = "Oceania"
   dat$continent[d.s %in% "Lamont_WA_Merridin_etc"                  ] = "Oceania"
   dat$continent[d.s %in% "Lamont_WA_Millbrook"                     ] = "Oceania"
   dat$continent[d.s %in% "Lamont_WA_Stirling_Ranges"               ] = "Oceania"
   dat$continent[d.s %in% "Lamont_WA_Walpole"                       ] = "Oceania"
   dat$continent[d.s %in% "Lamont_WA_Watheroo"                      ] = "Oceania"
   dat$continent[d.s %in% "Lee_Cedar_Creek2"                        ] = "North America"
   dat$continent[d.s %in% "Lee_NZ_Murchison_Mtns"                   ] = "Oceania"
   dat$continent[d.s %in% "Lusk_saplings_Cordillera_Pelada"         ] = "South America"
   dat$continent[d.s %in% "Lusk-adults_Concepcion"                  ] = "South America"
   dat$continent[d.s %in% "Lusk-adults_Los Lleuques"                ] = "South America"
   dat$continent[d.s %in% "Lusk-adults_Puyehue"                     ] = "South America"
   dat$continent[d.s %in% "Marin_Medina_Piritu_Venezuela"           ] = "South America"
   dat$continent[d.s %in% "Martin_etal_Guanacaste"                  ] = "Central America"
   dat$continent[d.s %in% "McAllister_Konza"                        ] = "North America"
   dat$continent[d.s %in% "Mediavilla_et_al_Salamanca"              ] = "Europe"
   dat$continent[d.s %in% "MidgelySA_Alexandria"                    ] = "Africa"
   dat$continent[d.s %in% "MidgelySA_Amatolas"                      ] = "Africa"
   dat$continent[d.s %in% "MidgelySA_Dukuduku"                      ] = "Africa"
   dat$continent[d.s %in% "MidgelySA_Jonkershoek_Mtn"               ] = "Africa"
   dat$continent[d.s %in% "MidgelySA_Jonkershoek_Rip"               ] = "Africa"
   dat$continent[d.s %in% "MidgelySA_Knysna"                        ] = "Africa"
   dat$continent[d.s %in% "MidgelySA_Mapelane"                      ] = "Africa"
   dat$continent[d.s %in% "MidgelySA_Sand_Forest"                   ] = "Africa"
   dat$continent[d.s %in% "MidgelySA_Umtiza"                        ] = "Africa"
   dat$continent[d.s %in% "MitchellNC_Coweeta"                      ] = "North America"
   dat$continent[d.s %in% "Miyazawa_Chiba_M"                        ] = "Asia"
   dat$continent[d.s %in% "Mooney_etal_81_desert"                   ] = "North America"
   dat$continent[d.s %in% "Mooney_etal_81_old-field"                ] = "North America"
   dat$continent[d.s %in% "Mooney_etal83_Jonk_Mtn"                  ] = "Africa"
   dat$continent[d.s %in% "Mulkey9193_BCI_Panama"                   ] = "Central America"
   dat$continent[d.s %in% "Nelson_etal_Texas"                       ] = "North America"
   dat$continent[d.s %in% "Niinemets_Kull94_Estonia"                ] = "Europe"
   dat$continent[d.s %in% "Niinemets_Kull98_Tartu"                  ] = "Europe"
   dat$continent[d.s %in% "Nitta_Chiba_N"                           ] = "Asia"
   dat$continent[d.s %in% "OleksynPol_Siemanice"                    ] = "Europe"
   dat$continent[d.s %in% "Olivares_Caracas_Venezuela"              ] = "South America"
   dat$continent[d.s %in% "Osada_Thomas_Pasoh"                      ] = "Asia"
   dat$continent[d.s %in% "Poorter_de_Jong_Along_ditch "            ] = "Europe"
   dat$continent[d.s %in% "Poorter_de_Jong_Chalk_grassland"         ] = "Europe"
   dat$continent[d.s %in% "Poorter_de_Jong_Dry_heath"               ] = "Europe"
   dat$continent[d.s %in% "Poorter_de_Jong_Dry_open_grassland"      ] = "Europe"
   dat$continent[d.s %in% "Poorter_de_Jong_Poor_hay_meadow"         ] = "Europe"
   dat$continent[d.s %in% "Poorter_de_Jong_Quaking_Fen"             ] = "Europe"
   dat$continent[d.s %in% "Poorter_de_Jong_Reed_marsh"              ] = "Europe"
   dat$continent[d.s %in% "Poorter_de_Jong_Wet_heath"               ] = "Europe"
   dat$continent[d.s %in% "Prado&DeMoraes1997_SAO_CARLOS"           ] = "South America"
   dat$continent[d.s %in% "Prior_dry monsoon forest"                ] = "Oceania"
   dat$continent[d.s %in% "Prior_open forest"                       ] = "Oceania"
   dat$continent[d.s %in% "Prior_swamp"                             ] = "Oceania"
   dat$continent[d.s %in% "Prior_woodland"                          ] = "Oceania"
   dat$continent[d.s %in% "Pyankov_Tadjikistan_Tadjikistan_high"    ] = "Asia"
   dat$continent[d.s %in% "Pyankov_Tadjikistan_Tadjikistan_higher"  ] = "Asia"
   dat$continent[d.s %in% "Pyankov_Tadjikistan_Tadjikistan_highest" ] = "Asia"
   dat$continent[d.s %in% "Pyankov_Urals_Yekaterinburg"             ] = "Asia"
   dat$continent[d.s %in% "Reichetal_Colorado"                      ] = "North America"
   dat$continent[d.s %in% "Reichetal_N_Carolina"                    ] = "North America"
   dat$continent[d.s %in% "Reichetal_New_Mexico"                    ] = "North America"
   dat$continent[d.s %in% "Reichetal_S_Carolina"                    ] = "North America"
   dat$continent[d.s %in% "Reichetal_Venezuela"                     ] = "South America"
   dat$continent[d.s %in% "Reichetal_Wisconsin"                     ] = "North America"
   dat$continent[d.s %in% "Ricklefs_SE_Ontario"                     ] = "North America"
   dat$continent[d.s %in% "Schulze_Kapalga"                         ] = "Oceania"
   dat$continent[d.s %in% "Schulze_Katherine"                       ] = "Oceania"
   dat$continent[d.s %in% "Schulze_Kidman Springs"                  ] = "Oceania"
   dat$continent[d.s %in% "Schulze_Melville Island"                 ] = "Oceania"
   dat$continent[d.s %in% "Schulze_Mt_Sanford"                      ] = "Oceania"
   dat$continent[d.s %in% "Shipley_Sherbrooke"                      ] = "North America"
   dat$continent[d.s %in% "Small1972_OTTAWA"                        ] = "North America"
   dat$continent[d.s %in% "Sobrado_Charallave"                      ] = "South America"
   dat$continent[d.s %in% "Sobrado&Medina_SanCarlos_bana"           ] = "South America"
   dat$continent[d.s %in% "Specht_Rundel_Dark_Island_heath"         ] = "Oceania"
   dat$continent[d.s %in% "Specht_Rundel_Dark_Island_mallee"        ] = "Oceania"
   dat$continent[d.s %in% "Specht_Rundel_Mt_Lofty"                  ] = "Oceania"
   dat$continent[d.s %in% "Tan_et_al_adinandra_trema_belukar"       ] = "Asia"
   dat$continent[d.s %in% "Terashima_Nepal"                         ] = "Asia"
   dat$continent[d.s %in% "Tezara_etal98_Coro"                      ] = "South America"
   dat$continent[d.s %in% "Tjoelker_Cedar_Creek"                    ] = "North America"
   dat$continent[d.s %in% "Tuohy_etal_Zimbabwe_CHID"                ] = "Africa"
   dat$continent[d.s %in% "Tuohy_etal_Zimbabwe_CRST_MCLW"           ] = "Africa"
   dat$continent[d.s %in% "Tuohy_etal_Zimbabwe_MAT"                 ] = "Africa"
   dat$continent[d.s %in% "Turner_&_Tan_Adinandra_Belukar"          ] = "Asia"
   dat$continent[d.s %in% "Turner_&_Tan_Beach_forest"               ] = "Asia"
   dat$continent[d.s %in% "Turner_&_Tan_Mangroves"                  ] = "Asia"
   dat$continent[d.s %in% "Turner_&_Tan_Undegraded_secondary_forest"] = "Asia"
   dat$continent[d.s %in% "Veneklaas_W_Australia"                   ] = "Oceania"
   dat$continent[d.s %in% "Villar_Andalucia_mesic"                  ] = "Europe"
   dat$continent[d.s %in% "Villar_Andalucia_xeric"                  ] = "Europe"
   dat$continent[d.s %in% "Villar_California_chaparral"             ] = "North America"
   dat$continent[d.s %in% "Villar_California_forest"                ] = "North America"
   dat$continent[d.s %in% "Villar_Canary_Is_lauriphyll"             ] = "Africa"
   dat$continent[d.s %in% "Villar_Canary_Is_xeric"                  ] = "Africa"
   dat$continent[d.s %in% "Villar_Chihuahua"                        ] = "North America"
   dat$continent[d.s %in% "Villar_Devon_Is_Canada"                  ] = "North America"
   dat$continent[d.s %in% "Villar_Douala-Edea Forest, Cameroon"     ] = "Africa"
   dat$continent[d.s %in% "Villar_Kibale Forest, Uganda"            ] = "Africa"
   dat$continent[d.s %in% "Villar_N_Carolina_forest"                ] = "North America"
   dat$continent[d.s %in% "Villar_Tierra_del_Fuego"                 ] = "South America"
   dat$continent[d.s %in% "Villar_Toronto"                          ] = "North America"
   dat$continent[d.s %in% "Williams et al_LosTuxtlas2"              ] = "North America"
   dat$continent[d.s %in% "Williams_Linera_Mexico"                  ] = "North America"
   dat$continent[d.s %in% "Wright_Oz_syd_hiP"                       ] = "Oceania"
   dat$continent[d.s %in% "Wright_Oz_syd_loP"                       ] = "Oceania"
   dat$continent[d.s %in% "Wright_Oz_wnsw_hiP"                      ] = "Oceania"
   dat$continent[d.s %in% "Wright_Oz_wnsw_loP"                      ] = "Oceania"
   dat$continent[d.s %in% "Zotz_Fortuna_Panama"                     ] = "Central America"
   #---------------------------------------------------------------------------------------#


   #----- Assign habitat. -----------------------------------------------------------------#
   dat$habitat = rep(NA_character_,nrow(dat))
   dat$habitat[d.b %in% "ALPINE"  ] = "alpine"
   dat$habitat[d.b %in% "BOREAL"  ] = "boreal"
   dat$habitat[d.b %in% "DESERT"  ] = "desert"
   dat$habitat[d.b %in% "GRASS/M" ] = "grassland"
   dat$habitat[d.b %in% "TEMP_FOR"] = "temperate forest"
   dat$habitat[d.b %in% "TEMP_RF" ] = "temperate rainforest"
   dat$habitat[d.b %in% "TROP_FOR"] = "tropical forest"
   dat$habitat[d.b %in% "TROP_RF" ] = "tropical rainforest"
   dat$habitat[d.b %in% "TUNDRA"  ] = "tundra"
   dat$habitat[d.b %in% "WLAND"   ] = "woodland"
   #---------------------------------------------------------------------------------------#


   #------ Append "tropical" to some habitats in tropical countries. ----------------------#
   trop.pot    = c("desert","grassland","woodland")
   trop.append = ( ( dat$country %in% trop.countries )
                 & ( dat$habitat %in% c("desert","grassland","woodland") )
                 )#end trop.append
   dat$habitat[trop.append] = paste("tropical",dat$habitat[trop.append])
   #---------------------------------------------------------------------------------------#

   return(dat)
}#end fix.glopnet.siteinfo
#==========================================================================================#
#==========================================================================================#




#==========================================================================================#
#==========================================================================================#
#      Separate location information from RAINFOR data base into habitat, country and      #
# continent.                                                                               #
#------------------------------------------------------------------------------------------#
fix.rainfor.siteinfo <<- function(dat){


   #------ Handy aliases. -----------------------------------------------------------------#
   d.s = dat$site
   #---------------------------------------------------------------------------------------#


   #----- Assign country. -----------------------------------------------------------------#
   dat$country = mapply( FUN  = switch
                       , EXPR = d.s
                       , MoreArgs = list( SUC_05 = "Peru"
                                        , TAM_05 = "Peru"
                                        , JEN_11 = "Peru"
                                        , ALP_01 = "Peru"
                                        , SUC_01 = "Peru"
                                        , JEN_12 = "Peru"
                                        , ALP_30 = "Peru"
                                        , CUZ_03 = "Peru"
                                        , ALP_40 = "Peru"
                                        , TAM_09 = "Peru"
                                        , TAM_06 = "Peru"
                                        , SPD_02 = "Peru"
                                        , SPD_01 = "Peru"
                                        , TRU_08 = "Peru"
                                        , ESP_01 = "Peru"
                                        , TRU_03 = "Peru"
                                        , WAQ_01 = "Peru"
                                        , TRU_01 = "Peru"
                                        , NA_character_
                                        )#end list
                       )#end mapply
   #---------------------------------------------------------------------------------------#


   #----- Assign continent. ---------------------------------------------------------------#
   dat$continent = mapply( FUN  = switch
                         , EXPR = d.s
                         , MoreArgs = list( SUC_05 = "South America"
                                          , TAM_05 = "South America"
                                          , JEN_11 = "South America"
                                          , ALP_01 = "South America"
                                          , SUC_01 = "South America"
                                          , JEN_12 = "South America"
                                          , ALP_30 = "South America"
                                          , CUZ_03 = "South America"
                                          , ALP_40 = "South America"
                                          , TAM_09 = "South America"
                                          , TAM_06 = "South America"
                                          , SPD_02 = "South America"
                                          , SPD_01 = "South America"
                                          , TRU_08 = "South America"
                                          , ESP_01 = "South America"
                                          , TRU_03 = "South America"
                                          , WAQ_01 = "South America"
                                          , TRU_01 = "South America"
                                          , NA_character_
                                          )#end list
                         )#end mapply
   #---------------------------------------------------------------------------------------#


   #----- Assign longitude. ---------------------------------------------------------------#
   dat$lon = mapply( FUN  = switch
                   , EXPR = d.s
                   , MoreArgs = list( SUC_05 = -72.8942
                                    , TAM_05 = -69.2705
                                    , JEN_11 = -73.6295
                                    , ALP_01 = -73.4333
                                    , SUC_01 = -72.9078
                                    , JEN_12 = -73.6276
                                    , ALP_30 = -73.4267
                                    , CUZ_03 = -69.0539
                                    , ALP_40 = -73.4400
                                    , TAM_09 = -69.2843
                                    , TAM_06 = -69.2960
                                    , SPD_02 = -71.5365
                                    , SPD_01 = -71.5423
                                    , TRU_08 = -71.5559
                                    , ESP_01 = -71.5948
                                    , TRU_03 = -71.5995
                                    , WAQ_01 = -71.5874
                                    , TRU_01 = -71.6069
                                    , NA_real_
                                    )#end list
                   )#end mapply
   #---------------------------------------------------------------------------------------#


   #----- Assign latitude. ----------------------------------------------------------------#
   dat$lat = mapply( FUN  = switch
                   , EXPR = d.s
                   , MoreArgs = list( SUC_05 =  -3.2558
                                    , TAM_05 = -12.8309
                                    , JEN_11 =  -4.8781
                                    , ALP_01 =  -3.9500
                                    , SUC_01 =  -3.2519
                                    , JEN_12 =  -4.8990
                                    , ALP_30 =  -3.9543
                                    , CUZ_03 = -12.5344
                                    , ALP_40 =  -3.9410
                                    , TAM_09 = -12.8309
                                    , TAM_06 = -12.8385
                                    , SPD_02 = -13.0491
                                    , SPD_01 = -13.0475
                                    , TRU_08 = -13.0702
                                    , ESP_01 = -13.1751
                                    , TRU_03 = -13.1097
                                    , WAQ_01 = -13.1908
                                    , TRU_01 = -13.1136
                                    , NA_real_
                                    )#end list
                   )#end mapply
   #---------------------------------------------------------------------------------------#


   #----- Assign altitude. ----------------------------------------------------------------#
   dat$alt = mapply( FUN  = switch
                   , EXPR = d.s
                   , MoreArgs = list( SUC_05 = 132.
                                    , TAM_05 = 223.
                                    , JEN_11 = 131.
                                    , ALP_01 = 120.
                                    , SUC_01 = 117.
                                    , JEN_12 = 135.
                                    , ALP_30 = 150.
                                    , CUZ_03 = 205.
                                    , ALP_40 = 142.
                                    , TAM_09 = 219.
                                    , TAM_06 = 215.
                                    , SPD_02 = 1527.
                                    , SPD_01 = 1776.
                                    , TRU_08 = 1885.
                                    , ESP_01 = 2863.
                                    , TRU_03 = 3044.
                                    , WAQ_01 = 3045.
                                    , TRU_01 = 3379.
                                    , NA_real_
                                    )#end list
                   )#end mapply
   #---------------------------------------------------------------------------------------#


   #----- Assign habitat. Assume alpine above 2000m. --------------------------------------#
   dat$habitat = mapply( FUN  = switch
                       , EXPR = d.s
                       , MoreArgs = list( SUC_05 = "tropical rainforest"
                                        , TAM_05 = "tropical rainforest"
                                        , JEN_11 = "tropical rainforest"
                                        , ALP_01 = "tropical rainforest"
                                        , SUC_01 = "tropical rainforest"
                                        , JEN_12 = "tropical rainforest"
                                        , ALP_30 = "tropical rainforest"
                                        , CUZ_03 = "tropical rainforest"
                                        , ALP_40 = "tropical rainforest"
                                        , TAM_09 = "tropical rainforest"
                                        , TAM_06 = "tropical rainforest"
                                        , SPD_02 = "tropical rainforest"
                                        , SPD_01 = "tropical rainforest"
                                        , TRU_08 = "tropical rainforest"
                                        , ESP_01 = "tropical forest"
                                        , TRU_03 = "tropical forest"
                                        , WAQ_01 = "tropical forest"
                                        , TRU_01 = "tropical forest"
                                        , NA_real_
                                        )#end list
                       )#end mapply
   #---------------------------------------------------------------------------------------#

   return(dat)
}#end fix.rainfor.siteinfo
#==========================================================================================#
#==========================================================================================#




#==========================================================================================#
#==========================================================================================#
#      Separate location information from NGEE-Tropics data base into habitat, country and #
# continent.                                                                               #
#------------------------------------------------------------------------------------------#
fix.ngee.tropics.siteinfo <<- function(dat){


   #------ Handy aliases. -----------------------------------------------------------------#
   d.s = dat$site
   d.g = dat$growth.form
   #---------------------------------------------------------------------------------------#


   #----- Assign country. -----------------------------------------------------------------#
   dat$country = mapply( FUN  = switch
                       , EXPR = d.s
                       , MoreArgs = list( PNM = "Panama"
                                        , SLZ = "Panama"
                                        , NA_character_
                                        )#end list
                       )#end mapply
   #---------------------------------------------------------------------------------------#


   #----- Assign continent. ---------------------------------------------------------------#
   dat$continent = mapply( FUN  = switch
                         , EXPR = d.s
                         , MoreArgs = list( PNM = "Central America"
                                          , SLZ = "Central America"
                                          , NA_character_
                                          )#end list
                         )#end mapply
   #---------------------------------------------------------------------------------------#


   #----- Assign longitude. ---------------------------------------------------------------#
   dat$lon = mapply( FUN  = switch
                   , EXPR = d.s
                   , MoreArgs = list( PNM = -79.54330
                                    , SLZ = -79.97452
                                    )#end list
                   )#end mapply
   #---------------------------------------------------------------------------------------#


   #----- Assign latitude. ----------------------------------------------------------------#
   dat$lat = mapply( FUN  = switch
                   , EXPR = d.s
                   , MoreArgs = list( PNM = 8.99441
                                    , SLZ = 9.28103
                                    , NA_real_
                                    )#end list
                   )#end mapply
   #---------------------------------------------------------------------------------------#


   #----- Assign altitude. ----------------------------------------------------------------#
   dat$alt = mapply( FUN  = switch
                   , EXPR = d.s
                   , MoreArgs = list( PNM =  71.
                                    , SLZ = 172.
                                    , NA_real_
                                    )#end list
                   )#end mapply
   #---------------------------------------------------------------------------------------#


   #----- Assign habitat. Assume alpine above 2000m. --------------------------------------#
   dat$habitat = mapply( FUN  = switch
                       , EXPR = d.s
                       , MoreArgs = list( PNM = "tropical forest"
                                        , SLZ = "tropical rainforest"
                                        , NA_real_
                                        )#end list
                       )#end mapply
   #---------------------------------------------------------------------------------------#


   #----- Assign habitat. Assume alpine above 2000m. --------------------------------------#
   dat$life.type = mapply( FUN  = switch
                         , EXPR = d.g
                         , MoreArgs = list( epiphye        = "EP"
                                          , fern           = "F"
                                          , liana          = "V"
                                          , perennial_herb = "H"
                                          , tree           = "T"
                                          , NA_real_
                                          )#end list
                         )#end mapply
   #---------------------------------------------------------------------------------------#

   return(dat)
}#end fix.ngee.tropics.siteinfo
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      Simplify habitats.                                                                  #
#------------------------------------------------------------------------------------------#
fix.try.habitat <<- function(x){
   #-----  List all possible habitats, then simplify them. --------------------------------#
   x50 = ifelse(test = x %in% "",yes=NA_character_,no=substring(x,1,50))
   #---------------------------------------------------------------------------------------#


   #----- o.s is original, simple classification. -----------------------------------------#
   o.s = rbind( c( origin = "Abandoned railway site (urban)"
                 , simple = "ruderal"
                 )#end c
              , c( origin = "alpine"
                 , simple = "alpine"
                 )#end c
              , c( origin = "Amazonian terra firme forest"
                 , simple = "tropical forest"
                 )#end c
              , c( origin = "aquatic"
                 , simple = "ruderal"
                 )#end c
              , c( origin = "arable"
                 , simple = "ruderal"
                 )#end c
              , c( origin = "coast; saline site"
                 , simple = "ruderal"
                 )#end c
              , c( origin = "desert erg"
                 , simple = "dryland"
                 )#end c
              , c( origin = "desert shrub"
                 , simple = "dryland"
                 )#end c
              , c( origin = "desert wadi"
                 , simple = "dryland"
                 )#end c
              , c( origin = "dry site"
                 , simple = "dryland"
                 )#end c
              , c( origin = "Evergreen lowland tropical forest on a terra firme"
                 , simple = "tropical forest"
                 )#end c
              , c( origin = "In gaps and as epiphytes; wet forest"
                 , simple = "ruderal"
                 )#end c
              , c( origin = "meadow"
                 , simple = "ruderal"
                 )#end c
              , c( origin = "Oak forest"
                 , simple = "temperate forest"
                 )#end c
              , c( origin = "ruderal"
                 , simple = "ruderal"
                 )#end c
              , c( origin = "subtropical"
                 , simple = "tropical forest"
                 )#end c
              , c( origin = "temperate"
                 , simple = "temperate forest"
                 )#end c
              , c( origin = "temperate forest"
                 , simple = "temperate forest"
                 )#end c
              , c( origin = "temperate rainforest"
                 , simple = "temperate forest"
                 )#end c
              , c( origin = "tropical"
                 , simple = "tropical forest"
                 )#end c
              , c( origin = "tropical desert"
                 , simple = "desert"
                 )#end c
              , c( origin = "tropical forest"
                 , simple = "tropical forest"
                 )#end c
              , c( origin = "tropical rainforest"
                 , simple = "tropical forest"
                 )#end c
              , c( origin = "tropical woodland"
                 , simple = "woodland"
                 )#end c
              , c( origin = "Wet forest understories"
                 , simple = "tropical forest"
                 )#end c
              , c( origin = "young secondary seasonally dry tropical forest"
                 , simple = "tropical forest"
                 )#end c
              )#end rbind
   o.s = as.data.frame(o.s,stringsAsFactors=FALSE)
   #---------------------------------------------------------------------------------------#



   #----- Match habitats with the original data. ------------------------------------------#
   idx      = match(x50,o.s$origin)
   sel      = is.finite(idx)
   ans      = rep(NA_character_,times=length(x))
   ans[sel] = o.s$simple[idx[sel]]
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end fix.try.habitat
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     Fix photosynthetic path.                                                             #
#                                                                                          #
#  Look-up table was built by combining two data sets.                                     #
#                                                                                          #
#  Christin P-A, Osborne CP, Sage RF, Araraki M, Edwards EJ. 2011. C4 eudicots are not     #
#     younger than C4 monocots. J. Ex. Bot., 62(9), 3171-3181. doi:10.1093/jxb/err041      #
#                                                                                          #
#  Osborne CP, Salomaa A, Kluyer TA, Visser V, Kellogg EA, Morrone O, Vorontsova MS,       #
#     Clayton WD, Simpson DA. 2014. A global database of C4 photosynthesis in grasses.     #
#     New Phytol., 204(3), 441-446. doi:10.1111/nph.12942                                  #
#                                                                                          #
#  Sage RF. 2016. A portrait of the C4 photosynthetic family on the 50th anniversary of    #
#     its discovery: species number, evolutionary lineages, and Hall of Fame.              #
#     J. Exp. Bot., 67(14), 4039-4056. doi:10.1093/jxb/erw156.                             #
#     (Currently only the 'Hall of fame' is included as species were not provided for      #
#     other lineages)                                                                      #
#------------------------------------------------------------------------------------------#
fix.try.photo.path <<- function(dat){
   #----- Read Osborne/Sage list of photosynthetic pathway. -------------------------------#
   photodb = read.csv( file = file.path(srcdir,"PhotoPath_DataBase.csv")
                     , header = TRUE
                     , stringsAsFactors = FALSE
                     )#end read.csv
   #---------------------------------------------------------------------------------------#

   #----- Initialise the results with the original values. --------------------------------#
   ans = ifelse( test = dat$photo.path %in% "",yes=NA_character_,no=dat$photo.path)
   #---------------------------------------------------------------------------------------#

   #----- Find species in the original data base that can be filled. ----------------------#
   sel      = is.na(ans) & (dat$scientific %in% photodb$scientific)
   idx      = match(dat$scientific,photodb$scientific)
   ans[sel] = photodb$photo.path[idx[sel]]
   #---------------------------------------------------------------------------------------#



   #------ Assume that non-filled trees must be C3. ---------------------------------------#
   c3.brute = is.na(ans) & (dat$life.type %in% "T")
   ans[c3.brute] = "C3"
   #---------------------------------------------------------------------------------------#


   #------ Standardise notation.  Everything in uppercase, with "C" attached. -------------#
   ans[ans %in% "3"] = "C3"
   ans[ans %in% "4"] = "C4"
   ans               = toupper(ans)
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end fix.try.photo.path
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     Fix life type based on a few additional information.                                 #
#------------------------------------------------------------------------------------------#
fix.try.life.type <<- function(dat){

   #----- Initialise the results with the original values. --------------------------------#
   ans  = ifelse( test = dat$life.type %in% "",yes=NA_character_,no=dat$life.type)
   #---------------------------------------------------------------------------------------#


   #----- First we use well-established taxonomy to flag ferns and mosses. ----------------#
   ans[dat$phylum %in% "Bryophyta"                         ] = "M"
   ans[dat$phylum %in% c("Polypodiophyta","Lycopodiophyta")] = "F"
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #       Second, assume that the life type is generally consistent for plants in the     #
   # same genus.                                                                           #
   #---------------------------------------------------------------------------------------#
   use        = ! is.na(ans)
   life.gen   = tapply(X=dat$life.type[use],INDEX=dat$genus[use],FUN=commonest,na.rm=TRUE)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      For the time being, assume that genera are homogeneous regarding the life type.  #
   #---------------------------------------------------------------------------------------#
   sel      = is.na(ans) & (dat$genus %in% names(life.gen))
   idx      = match(dat$genus,names(life.gen))
   ans[sel] = life.gen[idx[sel]]
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #       To-do, include a data base with lianas to not misattribute trees.  For the time #
   # being, assume that individuals that have finite wood density and are not flagged as   #
   # trees are also trees.                                                                 #
   #---------------------------------------------------------------------------------------#
   sel      = is.na(ans) & is.finite(dat$wood.dens)
   ans[sel] = "T"
   #---------------------------------------------------------------------------------------#


   #----- Return and hope for the best. ---------------------------------------------------#
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end fix.try.photo.path
#==========================================================================================#
#==========================================================================================#
