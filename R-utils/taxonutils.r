#==========================================================================================#
#==========================================================================================#
#     Global variables to define plant properties.                                         #
#------------------------------------------------------------------------------------------#
mypfts     <<- c(2,3,4)
pft.dens   <<- c(0.53,0.71,0.90)
pft.names  <<- c("Early tropical","Mid tropical","Late tropical")
n.pfts     <<- length(mypfts)
pft.breaks <<- c(-Inf,pft.dens[2:n.pfts]-0.5*diff(pft.dens),Inf)
#------------------------------------------------------------------------------------------#





#==========================================================================================#
#==========================================================================================#
#      This function standardises the spelling of common names of trees.  Based on         #
# Santarem survey,  but feel free to add more.                                             #
#------------------------------------------------------------------------------------------#
standard.common.name <<- function(x){
   #----- Make sure common names are lower case. ------------------------------------------#
   x   = tolower(x)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     General substitutions.  These are very aggressive, so don't use it too much.      #
   # Good things to put here are names that are often misspelt.  It's wise to use regexpr  #
   # rules such as ^ and $ to make sure gsub won't substitute more than it is supposed to. #
   #---------------------------------------------------------------------------------------#
   x = gsub(pattern="^abiuarana"      ,replacement="abiurana"    ,x=x)
   x = gsub(pattern="^angelin"        ,replacement="angelim"     ,x=x)
   x = gsub(pattern="^axicha"         ,replacement="axixa"       ,x=x)
   x = gsub(pattern="^jaracatia"      ,replacement="jacaratia"   ,x=x)
   x = gsub(pattern="^mara mara"      ,replacement="maramara"    ,x=x)
   x = gsub(pattern="^mara-mara"      ,replacement="maramara"    ,x=x)
   x = gsub(pattern="^mata mata"      ,replacement="matamata"    ,x=x)
   x = gsub(pattern="^mata-mata"      ,replacement="matamata"    ,x=x)
   x = gsub(pattern="^moratinga"      ,replacement="muiratinga"  ,x=x)
   x = gsub(pattern="^muiracatiara"   ,replacement="maracatiara" ,x=x)
   x = gsub(pattern="^muiraitinga"    ,replacement="muiratinga"  ,x=x)
   x = gsub(pattern="^muiriatinga"    ,replacement="muiratinga"  ,x=x)
   x = gsub(pattern="^muracatiara"    ,replacement="maracatiara" ,x=x)
   x = gsub(pattern="^muratinga"      ,replacement="muiratinga"  ,x=x)
   x = gsub(pattern="^murici"         ,replacement="muruci"      ,x=x)
   x = gsub(pattern="^mutuxi"         ,replacement="mututi"      ,x=x)
   x = gsub(pattern="muida$"          ,replacement="miuda"       ,x=x)
   x = gsub(pattern="^ocooba"         ,replacement="ucuuba"      ,x=x)
   x = gsub(pattern="^ocuuba"         ,replacement="ucuuba"      ,x=x)
   x = gsub(pattern="^papa terra"     ,replacement="papa-terra"  ,x=x)
   x = gsub(pattern="^papaterra"      ,replacement="papa-terra"  ,x=x)
   x = gsub(pattern="^para para"      ,replacement="parapara"    ,x=x)
   x = gsub(pattern="^para-para"      ,replacement="parapara"    ,x=x)
   x = gsub(pattern="^quina quina"    ,replacement="quinaquina"  ,x=x)
   x = gsub(pattern="^quina-quina"    ,replacement="quinaquina"  ,x=x)
   x = gsub(pattern="^tamarino"       ,replacement="tamarindo"   ,x=x)
   x = gsub(pattern="^taxi"           ,replacement="tachi"       ,x=x)
   x = gsub(pattern="^uchi"           ,replacement="uxi"         ,x=x)
   x = gsub(pattern="verdadiera"      ,replacement="verdadeira"  ,x=x)
   x = gsub(pattern="^xixua"          ,replacement="chichua"     ,x=x)
   #---------------------------------------------------------------------------------------#



   #----- Full replacements. --------------------------------------------------------------#
   sel = (x %in% "?"                             ); x[sel] = NA_character_
   sel = (x %in% "abacaba"                       ); x[sel] = "bacaba"
   sel = (x %in% "abicuiba"                      ); x[sel] = "ucuuba"
   sel = (x %in% "abirana rosadinha"             ); x[sel] = "abiu rosadinho"
   sel = (x %in% "abiui"                         ); x[sel] = "abiu"
   sel = (x %in% "abiu acariquara"               ); x[sel] = "abiu-acariquara"
   sel = (x %in% "abiu cramuri"                  ); x[sel] = "abiu-cramuri"
   sel = (x %in% "abiu cutite"                   ); x[sel] = "abiu-cutite"
   sel = (x %in% "abiu cutite folha verde"       ); x[sel] = "abiu-cutite"
   sel = (x %in% "abiu cutiti"                   ); x[sel] = "abiu-cutite"
   sel = (x %in% "abiu-cutiti"                   ); x[sel] = "abiu-cutite"
   sel = (x %in% "abiu guajara"                  ); x[sel] = "abiu-guajara"
   sel = (x %in% "abiu goiabao"                  ); x[sel] = "abiu-goiabao"
   sel = (x %in% "abiu mangabinha"               ); x[sel] = "abiu-mangabinha"
   sel = (x %in% "abiu tauari"                   ); x[sel] = "tauari"
   sel = (x %in% "abiu vermelha"                 ); x[sel] = "abiu vermelho"
   sel = (x %in% "abiurana vermelho"             ); x[sel] = "abiurana vermelha"
   sel = (x %in% "abiurana vermlho"              ); x[sel] = "abiurana vermelha"
   sel = (x %in% "abiuarana"                     ); x[sel] = "abiurana"
   sel = (x %in% "abiurana rosadinha"            ); x[sel] = "abiu rosadinho"
   sel = (x %in% "acoita cavalo"                 ); x[sel] = "acoita-cavalo"
   sel = (x %in% "algodao brabo"                 ); x[sel] = "algodao-bravo"
   sel = (x %in% "algodao bravo"                 ); x[sel] = "algodao-bravo"
   sel = (x %in% "amalelinha"                    ); x[sel] = "amarelinho"
   sel = (x %in% "amapa tirana"                  ); x[sel] = "amapatirana"
   sel = (x %in% "amarelinha"                    ); x[sel] = "amarelinho"
   sel = (x %in% "ameixa"                        ); x[sel] = "ameixa-do-para"
   sel = (x %in% "amerelinho"                    ); x[sel] = "amarelinho"
   sel = (x %in% "amescla"                       ); x[sel] = "breu mescla"
   sel = (x %in% "amoninha"                      ); x[sel] = "mamoninha"
   sel = (x %in% "ananin"                        ); x[sel] = "anani"
   sel = (x %in% "angelim aroreira"              ); x[sel] = "angelim-aroeira"
   sel = (x %in% "angelim aroeira"               ); x[sel] = "angelim-aroeira"
   sel = (x %in% "angelim margoso"               ); x[sel] = "angelim amargoso"
   sel = (x %in% "angelim pedra"                 ); x[sel] = "angelim-pedra"
   sel = (x %in% "angelim pedro"                 ); x[sel] = "angelim-pedra"
   sel = (x %in% "angelim peroba"                ); x[sel] = "angelim-pedra"
   sel = (x %in% "apuii"                         ); x[sel] = "apui"
   sel = (x %in% "aquariquara"                   ); x[sel] = "acariquara" 
   sel = (x %in% "aquariquarana"                 ); x[sel] = "acariquarana"
   sel = (x %in% "araca nego"                    ); x[sel] = "araca"
   sel = (x %in% "araca de anta"                 ); x[sel] = "araca-de-anta"
   sel = (x %in% "aratacio"                      ); x[sel] = "arataciu"
   sel = (x %in% "araticu"                       ); x[sel] = "araticum"
   sel = (x %in% "ata-menju"                     ); x[sel] = "atameju"
   sel = (x %in% "bacaba de leque"               ); x[sel] = "bacaba-de-leque"
   sel = (x %in% "bacabinha"                     ); x[sel] = "bacabi"
   sel = (x %in% "bacupari"                      ); x[sel] = "bacuripari"
   sel = (x %in% "bacuri-pari"                   ); x[sel] = "bacuripari"
   sel = (x %in% "bacuri pari"                   ); x[sel] = "bacuripari"
   sel = (x %in% "barbatimao"                    ); x[sel] = "fava-barbatimao"
   sel = (x %in% "bate puta"                     ); x[sel] = "batiputa"
   sel = (x %in% "baubarana"                     ); x[sel] = "embaubarana"
   sel = (x %in% "breu/louro preto?"             ); x[sel] = NA_character_
   sel = (x %in% "babao"                         ); x[sel] = "macauba"
   sel = (x %in% "bolao"                         ); x[sel] = "fava-bolota"
   sel = (x %in% "bombeira"                      ); x[sel] = "pau-pombo"
   sel = (x %in% "brau"                          ); x[sel] = "breu"
   sel = (x %in% "brejauba"                      ); x[sel] = "brejauva"              
   sel = (x %in% "breu sucuuba"                  ); x[sel] = "breu sucuruba"         
   sel = (x %in% "cabeca de urubu"               ); x[sel] = "cabeca-de-urubu"       
   sel = (x %in% "cabela"                        ); x[sel] = "louro canela"          
   sel = (x %in% "cabriuva"                      ); x[sel] = "cabreuva"              
   sel = (x %in% "cabriuna"                      ); x[sel] = "cabreuva"              
   sel = (x %in% "cabreu"                        ); x[sel] = "cabreuva"              
   sel = (x %in% "cabreuba"                      ); x[sel] = "cabreuva"              
   sel = (x %in% "caca piolho"                   ); x[sel] = "mata-piolho"           
   sel = (x %in% "cacau bravo"                   ); x[sel] = "cacaui"                
   sel = (x %in% "cacau da mata"                 ); x[sel] = "cacau"  
   sel = (x %in% "cacaurana"                     ); x[sel] = "cacau"                 
   sel = (x %in% "cachudinha"                    ); x[sel] = "cascudinha"            
   sel = (x %in% "cagaca"                        ); x[sel] = "abiurana-cagaca"
   sel = (x %in% "cajarana"                      ); x[sel] = "jarana"                
   sel = (x %in% "caju acu"                      ); x[sel] = "cajuacu"
   sel = (x %in% "cajuba"                        ); x[sel] = "caju"
   sel = (x %in% "calcho"                        ); x[sel] = "caucho"                
   sel = (x %in% "canafistula"                   ); x[sel] = "fava-marimari"
   sel = (x %in% "canela brava"                  ); x[sel] = "catuaba"
   sel = (x %in% "canela de anta"                ); x[sel] = "canela-de-anta"        
   sel = (x %in% "canela de veado"               ); x[sel] = "canela-de-veado"       
   sel = (x %in% "canela de velho"               ); x[sel] = "canela-de-velho"       
   sel = (x %in% "canelha velha"                 ); x[sel] = "canela-de-velho"       
   sel = (x %in% "canela de jacamim"             ); x[sel] = "canela-de-jacamim"     
   sel = (x %in% "canella de jacami"             ); x[sel] = "canela-de-jacamim"     
   sel = (x %in% "canella ge jacami"             ); x[sel] = "canela-de-jacamim"     
   sel = (x %in% "canniela"                      ); x[sel] = "canela"                
   sel = (x %in% "capa bode"                     ); x[sel] = "capa-bode"
   sel = (x %in% "capoeiro"                      ); x[sel] = "capueiro"
   sel = (x %in% "capoeiro preto"                ); x[sel] = "capueiro preto"
   sel = (x %in% "capoeiro branco"               ); x[sel] = "capueiro"
   sel = (x %in% "captiurana"                    ); x[sel] = "capitiurana"
   sel = (x %in% "capueiro branco"               ); x[sel] = "capueiro"
   sel = (x %in% "caqui branco"                  ); x[sel] = "caqui"
   sel = (x %in% "caqui folha grande"            ); x[sel] = "caqui"
   sel = (x %in% "carobia"                       ); x[sel] = "caroba"                
   sel = (x %in% "cascudinho"                    ); x[sel] = "cascudinha"            
   sel = (x %in% "cascudo"                       ); x[sel] = "cascudinha"            
   sel = (x %in% "castanha do brasil"            ); x[sel] = "castanha-do-para"
   sel = (x %in% "castanha do para"              ); x[sel] = "castanha-do-para"
   sel = (x %in% "castanha"                      ); x[sel] = "castanha-do-para"
   sel = (x %in% "castanheiro"                   ); x[sel] = "castanha-do-para"
   sel = (x %in% "castanha de galinha"           ); x[sel] = "castanha-de-galinha"   
   sel = (x %in% "castanha de periquito"         ); x[sel] = "castanha-de-periquito" 
   sel = (x %in% "castanha de sapocaia"          ); x[sel] = "castanha-sapucaia"     
   sel = (x %in% "castanha de sapucaia"          ); x[sel] = "castanha-sapucaia"     
   sel = (x %in% "castanha sapocaia"             ); x[sel] = "castanha-sapucaia"     
   sel = (x %in% "castanheira"                   ); x[sel] = "castanha-do-para"      
   sel = (x %in% "cauba"                         ); x[sel] = "macacauba"
   sel = (x %in% "cauxo"                         ); x[sel] = "caucho"                
   sel = (x %in% "caxeta"                        ); x[sel] = "caixeta"               
   sel = (x %in% "caxicha"                       ); x[sel] = "caxixa"           
   sel = (x %in% "caximbeira"                    ); x[sel] = "cachimbeiro"           
   sel = (x %in% "caximbeiro"                    ); x[sel] = "cachimbeiro"           
   sel = (x %in% "caxudinha"                     ); x[sel] = "cascudinha"            
   sel = (x %in% "cedroarana"                    ); x[sel] = "cedrorana"
   sel = (x %in% "cega jumenta"                  ); x[sel] = "cega-jumento"
   sel = (x %in% "chamaecrista"                  ); x[sel] = "coracao-de-negro"
   sel = (x %in% "chapeu de sol"                 ); x[sel] = "chapeu-de-sol"
   sel = (x %in% "chocolate"                     ); x[sel] = "cacau"                 
   sel = (x %in% "cip"                           ); x[sel] = "cipo"
   sel = (x %in% "cipo(dbh a 0.9m do chao)"      ); x[sel] = "cipo"
   sel = (x %in% "cipo+arapo"                    ); x[sel] = "cipo"
   sel = (x %in% "coco pau"                      ); x[sel] = "coco-pau"
   sel = (x %in% "comida de jaboti"              ); x[sel] = "comida-de-jabuti"      
   sel = (x %in% "comida de jabuti"              ); x[sel] = "comida-de-jabuti"      
   sel = (x %in% "comida de pomba"               ); x[sel] = "comida-de-pombo"       
   sel = (x %in% "comida de pombo"               ); x[sel] = "comida-de-pombo"       
   sel = (x %in% "coracao  de negro"             ); x[sel] = "coracao-de-negro"      
   sel = (x %in% "coracao de nego"               ); x[sel] = "coracao-de-negro"      
   sel = (x %in% "coracao de negro"              ); x[sel] = "coracao-de-negro"      
   sel = (x %in% "corante de indio"              ); x[sel] = "urucum"                
   sel = (x %in% "coro preto"                    ); x[sel] = "louro preto"           
   sel = (x %in% "coussarea racemosa"            ); x[sel] = "caferana"              
   sel = (x %in% "crista de mutum"               ); x[sel] = "crista"  
   sel = (x %in% "cucutiteriba"                  ); x[sel] = "cucutitiriba"
   sel = (x %in% "cumari"                        ); x[sel] = "cumaru"                
   sel = (x %in% "cumaru/apui"                   ); x[sel] = "apui"                  
   sel = (x %in% "cumaru de ferro"               ); x[sel] = "cumaru-ferro"          
   sel = (x %in% "cumatezinho/goiabinha fm"      ); x[sel] = "cumatezinho"           
   sel = (x %in% "cupiu"                         ); x[sel] = "cupui"
   sel = (x %in% "cupuacu da mata"               ); x[sel] = "cupuacu-da-mata"
   sel = (x %in% "cutiti"                        ); x[sel] = "abiu cutite"
   sel = (x %in% "cutiti"                        ); x[sel] = "abiu cutite"
   sel = (x %in% "cutite"                        ); x[sel] = "abiu cutite"
   sel = (x %in% "cutiteriba"                    ); x[sel] = "cucutitiriba"
   sel = (x %in% "cutitiriba"                    ); x[sel] = "cucutitiriba"
   sel = (x %in% "cruzeiro"                      ); x[sel] = "quina-cruzeiro"
   sel = (x %in% "dulacia/cachaceira"            ); x[sel] = "cachaceiro"
   sel = (x %in% "dulacia/cachaceira"            ); x[sel] = "cachaceiro"
   sel = (x %in% "einvira preta"                 ); x[sel] = "envira preta" 
   sel = (x %in% "embauba branco"                ); x[sel] = "embauba branca"        
   sel = (x %in% "embauba vick"                  ); x[sel] = "embauba"               
   sel = (x %in% "embauba torem"                 ); x[sel] = "embauba toren"  
   sel = (x %in% "embirata"                      ); x[sel] = "envira ata"            
   sel = (x %in% "embireira branca"              ); x[sel] = "envira"                
   sel = (x %in% "embireira rosa"                ); x[sel] = "envira"                
   sel = (x %in% "envira biriba"                 ); x[sel] = "envira-biriba"         
   sel = (x %in% "envira cana"                   ); x[sel] = "envira-cana"           
   sel = (x %in% "envira cabo de rodo"           ); x[sel] = "envira cabo-de-rodo"   
   sel = (x %in% "envira caju"                   ); x[sel] = "envira-caju"           
   sel = (x %in% "envira conduru"                ); x[sel] = "envira-conduru"        
   sel = (x %in% "envira cunauaru"               ); x[sel] = "envira-cunauaru"       
   sel = (x %in% "envira-mao-de-onca"            ); x[sel] = "envira mao-de-onca"    
   sel = (x %in% "envira preto"                  ); x[sel] = "envira preta"          
   sel = (x %in% "envira quiabo"                 ); x[sel] = "axixa"
   sel = (x %in% "envira-quiabo"                 ); x[sel] = "axixa"
   sel = (x %in% "envira sombrera"               ); x[sel] = "envira-sombreiro"
   sel = (x %in% "envira sombreira"              ); x[sel] = "envira-sombreiro"
   sel = (x %in% "envira sombreiro"              ); x[sel] = "envira-sombreiro"
   sel = (x %in% "envira surucu"                 ); x[sel] = "envira-surucucu"
   sel = (x %in% "envira surucucu"               ); x[sel] = "envira-surucucu"
   sel = (x %in% "envira taia"                   ); x[sel] = "envira-taia"           
   sel = (x %in% "envira turi"                   ); x[sel] = "envira-turi"
   sel = (x %in% "envira turi duro"              ); x[sel] = "envira-turi"
   sel = (x %in% "envira vermelho"               ); x[sel] = "envira vermelha"       
   sel = (x %in% "envireira"                     ); x[sel] = "envira"       
   sel = (x %in% "escorrega macaco"              ); x[sel] = "escorrega-macaco"      
   sel = (x %in% "escurrega macaco"              ); x[sel] = "escorrega-macaco"      
   sel = (x %in% "espeturana f. g."              ); x[sel] = "espeturana"            
   sel = (x %in% "fava amarg"                    ); x[sel] = "fava amargosa"
   sel = (x %in% "fava arara"                    ); x[sel] = "fava-arara-tucupi"
   sel = (x %in% "fava arara tucupi"             ); x[sel] = "fava-arara-tucupi"     
   sel = (x %in% "fava atana"                    ); x[sel] = "fava-atana"
   sel = (x %in% "fiora preta"                   ); x[sel] = NA_character_
   sel = (x %in% "fava barbatimao"               ); x[sel] = "fava-barbatimao"       
   sel = (x %in% "fava bolacha"                  ); x[sel] = "fava-bolacha"          
   sel = (x %in% "fava bolota"                   ); x[sel] = "fava-bolota"           
   sel = (x %in% "fava cana"                     ); x[sel] = "fava-cana"           
   sel = (x %in% "fava core"                     ); x[sel] = "fava-core"           
   sel = (x %in% "fava de anta"                  ); x[sel] = "fava-de-anta"          
   sel = (x %in% "fava marimari"                 ); x[sel] = "fava-marimari"         
   sel = (x %in% "fava mapuxique"                ); x[sel] = "fava-mapuxiqui"        
   sel = (x %in% "fava mapuxiqui"                ); x[sel] = "fava-mapuxiqui"        
   sel = (x %in% "fava orelha de macaco"         ); x[sel] = "fava orelha-de-macaco" 
   sel = (x %in% "fava paricarana"               ); x[sel] = "fava-paricana"         
   sel = (x %in% "fava saboeira"                 ); x[sel] = "fava-saboeira"         
   sel = (x %in% "fava tambori"                  ); x[sel] = "fava-tamboril"         
   sel = (x %in% "fava tamburi"                  ); x[sel] = "fava-tamboril"         
   sel = (x %in% "fava tamboril"                 ); x[sel] = "fava-tamboril" 
   sel = (x %in% "fava tamboriu"                 ); x[sel] = "fava-tamboril"  
   sel = (x %in% "favera amargosa"               ); x[sel] = "fava amargosa"         
   sel = (x %in% "faveira branca"                ); x[sel] = "fava branca"           
   sel = (x %in% "feijo branco"                  ); x[sel] = "freijo branco"         
   sel = (x %in% "ferdinandusa elliptica"        ); x[sel] = "bacabinha quina"       
   sel = (x %in% "figado de preguisa"            ); x[sel] = "figado-de-preguica"
   sel = (x %in% "figueira brava"                ); x[sel] = "figueira"              
   sel = (x %in% "gameleiro"                     ); x[sel] = "gameleira"             
   sel = (x %in% "gapeba"                        ); x[sel] = "guapeva"               
   sel = (x %in% "guapeba"                       ); x[sel] = "guapeva"               
   sel = (x %in% "gema de ovo"                   ); x[sel] = "gema-de-ovo"           
   sel = (x %in% "geniparana"                    ); x[sel] = "jeniparana"            
   sel = (x %in% "genipapo"                      ); x[sel] = "jenipapo"              
   sel = (x %in% "goiaba"                        ); x[sel] = "araca"
   sel = (x %in% "goiaba de anta"                ); x[sel] = "goiaba-de-anta"
   sel = (x %in% "goibarana"                     ); x[sel] = "goiabarana"            
   sel = (x %in% "gombeira fg"                   ); x[sel] = "gombeira folha grande"
   sel = (x %in% "gombeira vermelho"             ); x[sel] = "gombeira vermelha"     
   sel = (x %in% "grao de galo"                  ); x[sel] = "grao-de-galo"          
   sel = (x %in% "grao de guariba"               ); x[sel] = "grao-de-guariba"       
   sel = (x %in% "grao de macaco"                ); x[sel] = "grao-de-guariba" 
   sel = (x %in% "gravilola brava"               ); x[sel] = "graviola-brava"        
   sel = (x %in% "graviola brava"                ); x[sel] = "graviola-brava"        
   sel = (x %in% "guaiba"                        ); x[sel] = "araca"                 
   sel = (x %in% "guaiaba"                       ); x[sel] = "araca"
   sel = (x %in% "guajara bolacha"               ); x[sel] = "guajara-bolacha"       
   sel = (x %in% "guajara ferro"                 ); x[sel] = "guajara-ferro"         
   sel = (x %in% "guajara pedra"                 ); x[sel] = "guajara-pedra"         
   sel = (x %in% "guajara mirim"                 ); x[sel] = "guajara-mirim"         
   sel = (x %in% "guariuva"                      ); x[sel] = "guariuba"              
   sel = (x %in% "guaruba"                       ); x[sel] = "quaruba"
   sel = (x %in% "ibirucu"                       ); x[sel] = "embirucu"              
   sel = (x %in% "imbauba torem"                 ); x[sel] = "embauba toren"         
   sel = (x %in% "imbirata"                      ); x[sel] = "envira ata"            
   sel = (x %in% "imbireira"                     ); x[sel] = "envira"                
   sel = (x %in% "imbireira rosa"                ); x[sel] = "envira"                
   sel = (x %in% "imbiricu"                      ); x[sel] = "embirucu"              
   sel = (x %in% "imbirucu"                      ); x[sel] = "embirucu"              
   sel = (x %in% "inga a"                        ); x[sel] = "inga"
   sel = (x %in% "inga amarela"                  ); x[sel] = "inga amarelo"
   sel = (x %in% "inga branca"                   ); x[sel] = "inga branco"           
   sel = (x %in% "inga chichica"                 ); x[sel] = "inga-xixica"           
   sel = (x %in% "inga de orelha"                ); x[sel] = "inga-de-orelha"        
   sel = (x %in% "inga de preguica"              ); x[sel] = "inga-de-preguica"      
   sel = (x %in% "inga de rodo"                  ); x[sel] = "inga-de-rodo"          
   sel = (x %in% "inga de rosca"                 ); x[sel] = "inga-de-rosca"         
   sel = (x %in% "inga f g"                      ); x[sel] = "inga folha grauda"
   sel = (x %in% "inga f.p."                     ); x[sel] = "inga folha peluda"     
   sel = (x %in% "inga folha grande"             ); x[sel] = "inga folha grauda"
   sel = (x %in% "inga folhao"                   ); x[sel] = "inga folha grauda"
   sel = (x %in% "inga peluda"                   ); x[sel] = "inga peludo"         
   sel = (x %in% "inga titica"                   ); x[sel] = "inga-xixica"           
   sel = (x %in% "inga vermelha"                 ); x[sel] = "inga vermelho"         
   sel = (x %in% "inga xixica"                   ); x[sel] = "inga-xixica"           
   sel = (x %in% "ipe amerelo"                   ); x[sel] = "ipe amarelo"           
   sel = (x %in% "jaboticaba"                    ); x[sel] = "jabuticaba"            
   sel = (x %in% "jacare"                        ); x[sel] = "pau-jacare"
   sel = (x %in% "jacariuba"                     ); x[sel] = "jacareuba"
   sel = (x %in% "jambo"                         ); x[sel] = "jambo-do-mato"
   sel = (x %in% "jara"                          ); x[sel] = "jarana"
   sel = (x %in% "jaruma"                        ); x[sel] = "taruma"
   sel = (x %in% "jauari"                        ); x[sel] = "tauari"
   sel = (x %in% "jenita"                        ); x[sel] = "janita"                
   sel = (x %in% "jito"                          ); x[sel] = "gito"
   sel = (x %in% "joao mole"                     ); x[sel] = "joao-mole"             
   sel = (x %in% "joao moleza"                   ); x[sel] = "joao-moleza"           
   sel = (x %in% "jotobazinho"                   ); x[sel] = "jatobazinho"           
   sel = (x %in% "jutai mirim"                   ); x[sel] = "jutai-mirim"           
   sel = (x %in% "jutai acu"                     ); x[sel] = "jutai-acu"             
   sel = (x %in% "jutai pororoca"                ); x[sel] = "jutai-pororoca"        
   sel = (x %in% "lacre da mata"                 ); x[sel] = "lacre-da-mata"         
   sel = (x %in% "laranginha"                    ); x[sel] = "laranjinha"            
   sel = (x %in% "leiteiro"                      ); x[sel] = "leiteira"              
   sel = (x %in% "leitera"                       ); x[sel] = "leiteira"              
   sel = (x %in% "loro amarelo"                  ); x[sel] = "louro amarelo"         
   sel = (x %in% "louro?"                        ); x[sel] = "louro"
   sel = (x %in% "louro abacate"                 ); x[sel] = "louro-abacate"         
   sel = (x %in% "louro aritu"                   ); x[sel] = "louro-aritu"           
   sel = (x %in% "louro branco"                  ); x[sel] = "louro"                 
   sel = (x %in% "louro bosta"                   ); x[sel] = "louro-bosta"           
   sel = (x %in% "louro canela"                  ); x[sel] = "louro canelado" 
   sel = (x %in% "louro chumbo"                  ); x[sel] = "louro-chumbo"          
   sel = (x %in% "louro p"                       ); x[sel] = "louro preto"
   sel = (x %in% "louro pimenta"                 ); x[sel] = "louro-pimenta"         
   sel = (x %in% "louro seda"                    ); x[sel] = "louro-seda"            
   sel = (x %in% "macucu de sangue"              ); x[sel] = "macucu-de-sangue"
   sel = (x %in% "mafim"                         ); x[sel] = "marfim"
   sel = (x %in% "mamao jacatia"                 ); x[sel] = "jacaratia"             
   sel = (x %in% "mamonini"                      ); x[sel] = "mamoninha"             
   sel = (x %in% "mandioqueiro"                  ); x[sel] = "mandioqueira"
   sel = (x %in% "mandioqueiro escamoso"         ); x[sel] = "mandioqueira"
   sel = (x %in% "mangueira"                     ); x[sel] = "manguerana"
   sel = (x %in% "manguerano"                    ); x[sel] = "manguerana"
   sel = (x %in% "maparajuba"                    ); x[sel] = "parajuba"
   sel = (x %in% "maprounea"                     ); x[sel] = "caxixa"
   sel = (x %in% "mapuxique"                     ); x[sel] = "fava-mapuxiqui"
   sel = (x %in% "mapuxiqui"                     ); x[sel] = "fava-mapuxiqui"
   sel = (x %in% "maquira"                       ); x[sel] = "muiratinga"
   sel = (x %in% "maracata"                      ); x[sel] = "marassacaca"
   sel = (x %in% "maracatia"                     ); x[sel] = "muiracatiara"
   sel = (x %in% "maracacaca"                    ); x[sel] = "marassacaca"
   sel = (x %in% "marasacaca"                    ); x[sel] = "marassacaca"
   sel = (x %in% "marimari"                      ); x[sel] = "fava-marimari"
   sel = (x %in% "massaranduba"                  ); x[sel] = "macaranduba"
   sel = (x %in% "matamata branca"               ); x[sel] = "matamata branco"
   sel = (x %in% "matamata ci"                   ); x[sel] = "matamata-ci"
   sel = (x %in% "matamata cinza"                ); x[sel] = "matamata-ci"
   sel = (x %in% "matamata jiboia"               ); x[sel] = "matamata-jiboia"
   sel = (x %in% "matamata vermelha"             ); x[sel] = "matamata vermelho"
   sel = (x %in% "mata pau+jito"                 ); x[sel] = "gito"
   sel = (x %in% "mata caldo"                    ); x[sel] = "mata-calado"
   sel = (x %in% "melanciera"                    ); x[sel] = "melancieira"
   sel = (x %in% "morto"                         ); x[sel] = "morta"
   sel = (x %in% "muiratinga folha grande/amapa" ); x[sel] = "muiratinga folha grande"
   sel = (x %in% "muiratinga fura fura"          ); x[sel] = "muiratinga fura-fura"
   sel = (x %in% "mulatero"                      ); x[sel] = "mulateiro"
   sel = (x %in% "mulugu"                        ); x[sel] = "mulungu"
   sel = (x %in% "murta da mata"                 ); x[sel] = "murta-da-mata"
   sel = (x %in% "murucidu mata"                 ); x[sel] = "muruci-da-mata"
   sel = (x %in% "muruci da mata"                ); x[sel] = "muruci-da-mata"
   sel = (x %in% "muruci fp"                     ); x[sel] = "muruci folha peluda"
   sel = (x %in% "mutama"                        ); x[sel] = "mutambo"
   sel = (x %in% "mutamba"                       ); x[sel] = "mutambo"
   sel = (x %in% "mututiassu"                    ); x[sel] = "mututi-acu"
   sel = (x %in% "orelha de burro"               ); x[sel] = "orelha-de-burro"
   sel = (x %in% "orelha de macaco"              ); x[sel] = "fava orelha-de-macaco"
   sel = (x %in% "olho de sapo"                  ); x[sel] = "olho-de-sapo"
   sel = (x %in% "olho de veado"                 ); x[sel] = "olho-de-veado"
   sel = (x %in% "olho de viado"                 ); x[sel] = "olho-de-veado"
   sel = (x %in% "ouro branco"                   ); x[sel] = "seringueira"
   sel = (x %in% "p bolacha"                     ); x[sel] = "guajara-bolacha"
   sel = (x %in% "paineira"                      ); x[sel] = "sumauma"
   sel = (x %in% "palmito"                       ); x[sel] = "acai"
   sel = (x %in% "palmito babosa"                ); x[sel] = "acai"
   sel = (x %in% "pao de sangue"                 ); x[sel] = "pau-sangue"
   sel = (x %in% "papo de mutum"                 ); x[sel] = "papo-de-mutum"
   sel = (x %in% "papo de mutum ff"              ); x[sel] = "papo-de-mutum folha fina"
   sel = (x %in% "papo de mutum  ff"             ); x[sel] = "papo-de-mutum folha fina"
   sel = (x %in% "papo-de-mutum ff"              ); x[sel] = "papo-de-mutum folha fina"
   sel = (x %in% "papo-de-mutum  ff"             ); x[sel] = "papo-de-mutum folha fina"
   sel = (x %in% "paricarana"                    ); x[sel] = "fava-paricana"
   sel = (x %in% "passarinhiera"                 ); x[sel] = "passarinheira"
   sel = (x %in% "pata de vaca"                  ); x[sel] = "pata-de-vaca"
   sel = (x %in% "pata preta"                    ); x[sel] = "ata preta"
   sel = (x %in% "patua"                         ); x[sel] = "pataua"
   sel = (x %in% "pau de arco"                   ); x[sel] = "pau-de-arco"
   sel = (x %in% "pau d.arco"                    ); x[sel] = "pau-de-arco"
   sel = (x %in% "pau de bicho"                  ); x[sel] = "pau-de-bicho"
   sel = (x %in% "pau de cobra"                  ); x[sel] = "pau-cobra"
   sel = (x %in% "pau colher"                    ); x[sel] = "pau-de-colher"
   sel = (x %in% "pau de colher"                 ); x[sel] = "pau-de-colher"
   sel = (x %in% "pau de jacare"                 ); x[sel] = "pau-jacare"
   sel = (x %in% "pau de macaco"                 ); x[sel] = "pau-de-macaco"
   sel = (x %in% "pau de rego"                   ); x[sel] = "pau-de-remo"
   sel = (x %in% "pau de remo"                   ); x[sel] = "pau-de-remo"
   sel = (x %in% "pau de sangue"                 ); x[sel] = "pau-sangue"
   sel = (x %in% "pau jacare"                    ); x[sel] = "pau-jacare"
   sel = (x %in% "pau marfim"                    ); x[sel] = "pau-marfim"
   sel = (x %in% "pau mulato"                    ); x[sel] = "pau-mulato"
   sel = (x %in% "pau para tudo"                 ); x[sel] = "pau-para-tudo"
   sel = (x %in% "pau pereira"                   ); x[sel] = "peroba mica"
   sel = (x %in% "pau-pra-tudo"                  ); x[sel] = "pau-para-tudo"
   sel = (x %in% "pau pra tudo"                  ); x[sel] = "pau-para-tudo"
   sel = (x %in% "paupratudo"                    ); x[sel] = "pau-para-tudo"
   sel = (x %in% "pau prego"                     ); x[sel] = "pau-de-remo"
   sel = (x %in% "pau purui"                     ); x[sel] = "purui"
   sel = (x %in% "pau rego"                      ); x[sel] = "pau-de-remo"
   sel = (x %in% "pau sangue"                    ); x[sel] = "pau-sangue"
   sel = (x %in% "pereauna"                      ); x[sel] = "perebuna"
   sel = (x %in% "pedra umi"                     ); x[sel] = "pedra ume-caa"
   sel = (x %in% "pelo de cutia"                 ); x[sel] = "pelo-de-cutia"
   sel = (x %in% "pente de macaco"               ); x[sel] = "pente-de-macaco"
   sel = (x %in% "pepino da mata"                ); x[sel] = "pepino-do-mato"
   sel = (x %in% "pepino-da-mata"                ); x[sel] = "pepino-do-mato"
   sel = (x %in% "pepino do mato"                ); x[sel] = "pepino-do-mato"
   sel = (x %in% "pepino-do-mato"                ); x[sel] = "pepino-do-mato"
   sel = (x %in% "perna de moca"                 ); x[sel] = "perna-de-moca"
   sel = (x %in% "piqui"                         ); x[sel] = "piquia"  
   sel = (x %in% "piqui rosa"                    ); x[sel] = "piquia"
   sel = (x %in% "piquiazeiro"                   ); x[sel] = "piquia"
   sel = (x %in% "pororoca"                      ); x[sel] = "jutai-pororoca"
   sel = (x %in% "prapara"                       ); x[sel] = "parapara"
   sel = (x %in% "pratudo"                       ); x[sel] = "pau-para-tudo"
   sel = (x %in% "puruirana/purui branco"        ); x[sel] = "purui branco"
   sel = (x %in% "quaiquara"                     ); x[sel] = "acariquara"
   sel = (x %in% "quariquara"                    ); x[sel] = "acariquara"
   sel = (x %in% "quariquarana"                  ); x[sel] = "acariquara"
   sel = (x %in% "quariquari"                    ); x[sel] = "acariquara"            
   sel = (x %in% "quari quari"                   ); x[sel] = "acariquara"            
   sel = (x %in% "quebrado"                      ); x[sel] = NA_character_
   sel = (x %in% "quina"                         ); x[sel] = "quinarana" 
   sel = (x %in% "quina cruzeiro"                ); x[sel] = "quina-cruzeiro"        
   sel = (x %in% "rim de paca"                   ); x[sel] = "rim-de-paca"           
   sel = (x %in% "ripeiro"                       ); x[sel] = "ripeira"               
   sel = (x %in% "roxao"                         ); x[sel] = "roxinho"               
   sel = (x %in% "roxinao"                       ); x[sel] = "roxinho"               
   sel = (x %in% "sangra de agua"                ); x[sel] = "sangra-de-agua"
   sel = (x %in% "sapucaia"                      ); x[sel] = "castanha-sapucaia"
   sel = (x %in% "saboeira"                      ); x[sel] = "fava-saboeiro"
   sel = (x %in% "saboeira amarela"              ); x[sel] = "fava-saboeiro amarela"
   sel = (x %in% "saboeiro"                      ); x[sel] = "fava-saboeiro"
   sel = (x %in% "saboiera"                      ); x[sel] = "fava-saboeiro"
   sel = (x %in% "sabueira"                      ); x[sel] = "fava-saboeiro"  
   sel = (x %in% "sabueiro"                      ); x[sel] = "fava-saboeiro"
   sel = (x %in% "sabuguero"                     ); x[sel] = "sabugueiro"
   sel = (x %in% "sajinera"                      ); x[sel] = NA_character_
   sel = (x %in% "samauma de terra firme"        ); x[sel] = "sumauma da terra firme" 
   sel = (x %in% "sangra d`agua"                 ); x[sel] = "sangra-de-agua" 
   sel = (x %in% "sardinheiro"                   ); x[sel] = "sardinheira"
   sel = (x %in% "segador"                       ); x[sel] = "cegador"               
   sel = (x %in% "seringa"                       ); x[sel] = "seringueira"
   sel = (x %in% "seringa branca"                ); x[sel] = "seringueira"           
   sel = (x %in% "seringa branco"                ); x[sel] = "seringueira"           
   sel = (x %in% "seringa verdadeira"            ); x[sel] = "seringueira"           
   sel = (x %in% "seringarana preta"             ); x[sel] = "seringarana"           
   sel = (x %in% "seritinga"                     ); x[sel] = "seringueira"           
   sel = (x %in% "sorveira"                      ); x[sel] = "sorva"
   sel = (x %in% "sorveira leite"                ); x[sel] = "sorva"                 
   sel = (x %in% "sorvo"                         ); x[sel] = "sorva"
   sel = (x %in% "sova"                          ); x[sel] = "sorva"
   sel = (x %in% "sucupira p sapo"               ); x[sel] = "sucupira pele-de-sapo" 
   sel = (x %in% "sucupira pele de sapo"         ); x[sel] = "sucupira pele-de-sapo" 
   sel = (x %in% "sucuuba preta"                 ); x[sel] = "sucuuba" 
   sel = (x %in% "tachi branca"                  ); x[sel] = "tachi branco"          
   sel = (x %in% "tachi preta"                   ); x[sel] = "tachi preto"           
   sel = (x %in% "tachi preto ???"               ); x[sel] = "tachi preto"
   sel = (x %in% "tachi preto folh"              ); x[sel] = "tachi preto"
   sel = (x %in% "tachi vermelha"                ); x[sel] = "tachi vermelho"        
   sel = (x %in% "talquari"                      ); x[sel] = "tauari"                
   sel = (x %in% "tamaquarao"                    ); x[sel] = "tamaquare"             
   sel = (x %in% "tamarindu"                     ); x[sel] = "tamarindo"              
   sel = (x %in% "tamauma"                       ); x[sel] = "sumauma"             
   sel = (x %in% "tamboril"                      ); x[sel] = "fava-tamboril"
   sel = (x %in% "tamboriul"                     ); x[sel] = "fava-tamboril"
   sel = (x %in% "tanari roxo"                   ); x[sel] = "tauari"
   sel = (x %in% "tangarana"                     ); x[sel] = "tangirana"
   sel = (x %in% "tanimbuca"                     ); x[sel] = "tanibuca"
   sel = (x %in% "tapiririca"                    ); x[sel] = "tatapiririca"
   sel = (x %in% "tatapirirca"                   ); x[sel] = "tatapiririca"
   sel = (x %in% "tatapiririca verm."            ); x[sel] = "tatapiririca vermelha"
   sel = (x %in% "taturana"                      ); x[sel] = "taturuba"
   sel = (x %in% "tauri"                         ); x[sel] = "tauari"                
   sel = (x %in% "tento folha"                   ); x[sel] = "tento"
   sel = (x %in% "tento foha grauda"             ); x[sel] = "tento folha grauda"    
   sel = (x %in% "tintero"                       ); x[sel] = "tinteiro"              
   sel = (x %in% "titiriba"                      ); x[sel] = "cucutitiriba"              
   sel = (x %in% "tucuma acu"                    ); x[sel] = "tucuma-acu"
   sel = (x %in% "ucuarana"                      ); x[sel] = "urucurana"             
   sel = (x %in% "ucuuba da varzea"              ); x[sel] = "ucuuba-da-varzea"    
   sel = (x %in% "ucuuba da terra firme"         ); x[sel] = "ucuuba terra-firme"    
   sel = (x %in% "ucuuba terra firme"            ); x[sel] = "ucuuba terra-firme"    
   sel = (x %in% "ucuuba tf"                     ); x[sel] = "ucuuba terra-firme"    
   sel = (x %in% "ucuuba vermelho"               ); x[sel] = "ucuuba vermelha"       
   sel = (x %in% "umbia"                         ); x[sel] = "goiabarana"
   sel = (x %in% "unha de vaca"                  ); x[sel] = "pata-de-vaca"          
   sel = (x %in% "uruci"                         ); x[sel] = "muruci"
   sel = (x %in% "urucu"                         ); x[sel] = "urucum"
   sel = (x %in% "urucuri"                       ); x[sel] = "urucum"
   sel = (x %in% "uruucurana"                    ); x[sel] = "urucurana"             
   sel = (x %in% "virola"                        ); x[sel] = "ucuuba"
   sel = (x %in% "xaonoquito"                    ); x[sel] = "pau vermelho"
   #---------------------------------------------------------------------------------------#

   return(x)
}#end function standard.common.name
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      This function corrects scientific names and families that are not correctly typed,  #
# are synonyms or have become obsolete.                                                    #
#------------------------------------------------------------------------------------------#
standard.scientific.name <<- function(dat){
   #---------------------------------------------------------------------------------------#
   #     First we make sure all scientific names and families have the genus and family    #
   # capitalised (e.g.  scientific name: Araucaria angustifolia; family: Araucariaceae.    #
   #---------------------------------------------------------------------------------------#
   nplants = nrow(dat)
   #----- Remove cf., assume that we know all species for sure. ---------------------------#
   dat$scientific = tolower(dat$scientific)
   #----- Break into genus and species. ---------------------------------------------------#
   g.s     = keep.gen.spe.only(dat$scientific)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Substitutions of scientific name.  We standardise these first so family becomes  #
   # easier.                                                                               #
   #---------------------------------------------------------------------------------------#
   g.s = sub("Aiouea densiflora"          ,"Aiouea laevis"                ,x=g.s)
   g.s = sub("Alibertia myrciifolia"      ,"Cordiera myrciifolia"         ,x=g.s)
   g.s = sub("Allophyllus floribunda"     ,"Allophylus floribundus"       ,x=g.s)
   g.s = sub("Amburana acreana"           ,"Amburana cearensis"           ,x=g.s)
   g.s = sub("Ampelocera endentula"       ,"Ampelocera edentula"          ,x=g.s)
   g.s = sub("Amphiodon effusus"          ,"Poecilanthe effusa"           ,x=g.s)
   g.s = sub("Amphirrhox longiflora"      ,"Amphirrhox longifolia"        ,x=g.s)
   g.s = sub("Amphirrhox surinamensis"    ,"Amphirrhox longifolia"        ,x=g.s)
   g.s = sub("Anadenanthera falcata"      ,"Anadenanthera peregrina"      ,x=g.s)
   g.s = sub("Anartia"                    ,"Tabernaemontana"              ,x=g.s)
   g.s = sub("Aniba roseodora"            ,"Aniba rosaeodora"             ,x=g.s)
   g.s = sub("Annona decicoma"            ,"Annona densicoma"             ,x=g.s)
   g.s = sub("Annona exsucca"             ,"Rollinia exsucca"             ,x=g.s)
   g.s = sub("Annona longifolia"          ,"Fusaea longifolia"            ,x=g.s)
   g.s = sub("Annona mucosa"              ,"Rollinia mucosa"              ,x=g.s)
   g.s = sub("Antirrhoea"                 ,"Antirhea"                     ,x=g.s)
   g.s = sub("Arthrophyllum"              ,"Polyscias"                    ,x=g.s)
   g.s = sub("Apeiba burchelii"           ,"Apeiba glabra"                ,x=g.s)
   g.s = sub("Apeiba echinata"            ,"Apeiba petoumo"               ,x=g.s)
   g.s = sub("Aspidosperma aracanga"      ,"Aspidosperma araracanga"      ,x=g.s)
   g.s = sub("Aspidosperma auriculata"    ,"Aspidosperma auriculatum"     ,x=g.s)
   g.s = sub("Aspidosperma cruentum"      ,"Aspidosperma desmanthum"      ,x=g.s)
   g.s = sub("Aspidosperma desmantum"     ,"Aspidosperma desmanthum"      ,x=g.s)
   g.s = sub("Aspidosperma desmathum"     ,"Aspidosperma desmanthum"      ,x=g.s)
   g.s = sub("Aspidosperma eteanun"       ,"Aspidosperma eteanum"         ,x=g.s)
   g.s = sub("Aspidosperma nitidum"       ,"Aspidosperma excelsum"        ,x=g.s)
   g.s = sub("Astrocaryum gynacant"       ,"Astrocaryum aculeatum"        ,x=g.s)
   g.s = sub("Astrocaryum gynacanthum"    ,"Astrocaryum aculeatum"        ,x=g.s)
   g.s = sub("Astronium le-cointei"       ,"Astronium lecointei"          ,x=g.s)
   g.s = sub("Austroplenckia populnea"    ,"Plenckia populnea"            ,x=g.s)
   g.s = sub("Balisia pedicelares"        ,"Albizia pedicellaris"         ,x=g.s)
   g.s = sub("Balizia pedicellaris"       ,"Albizia pedicellaris"         ,x=g.s)
   g.s = sub("Bauhinia jarensis"          ,"Bauhinia deleteme"            ,x=g.s)
   g.s = sub("Bellucia grossulariodis"    ,"Bellucia grossularioides"     ,x=g.s)
   g.s = sub("Beureria"                   ,"Calycanthus"                  ,x=g.s)
   g.s = sub("Bosqueia"                   ,"Trilepisium"                  ,x=g.s)
   g.s = sub("Bracteanthus glycycarpus"   ,"Siparuna glycycarpa"          ,x=g.s)
   g.s = sub("Brosimopsis obovata"        ,"Brosimum acutifolium"         ,x=g.s)
   g.s = sub("Brosimum autifolium"        ,"Brosimum acutifolium"         ,x=g.s)
   g.s = sub("Brosimum lactascens"        ,"Brosimum lactescens"          ,x=g.s)
   g.s = sub("Brosimum guianensis"        ,"Brosimum guianense"           ,x=g.s)
   g.s = sub("Brosimum obovata"           ,"Brosimum acutifolium"         ,x=g.s)
   g.s = sub("Byrsonima chrysophylla"     ,"Byrsonima spicata"            ,x=g.s)
   g.s = sub("Byrsonima estipulacea"      ,"Byrsonima stipulacea"         ,x=g.s)
   g.s = sub("Byrsonima schultesiana"     ,"Byrsonima arthropoda"         ,x=g.s)
   g.s = sub("Capirona ulei"              ,"Capirona decorticans"         ,x=g.s)
   g.s = sub("Casearia bracteifera"       ,"Casearia combaymensis"        ,x=g.s)
   g.s = sub("Capparis frondosa"          ,"Capparidastrum frondosum"     ,x=g.s)
   g.s = sub("Cecropia distachia"         ,"Cecropia distachya"           ,x=g.s)
   g.s = sub("Cedrela fistis"             ,"Cedrela fissilis"             ,x=g.s)
   g.s = sub("Cedrella odorata"           ,"Cedrela odorata"              ,x=g.s)
   g.s = sub("Caeseria"                   ,"Casearia"                     ,x=g.s)
   g.s = sub("Chanouchiton kapleri"       ,"Chaunochiton kappleri"        ,x=g.s)
   g.s = sub("Chaunochiton kapleri"       ,"Chaunochiton kappleri"        ,x=g.s)
   g.s = sub("Clarisia elicifolia"        ,"Clarisia ilicifolia"          ,x=g.s)
   g.s = sub("Compomanesia"               ,"Campomanesia"                 ,x=g.s)
   g.s = sub("Connarus perrottetes"       ,"Connarus perrottetii"         ,x=g.s)
   g.s = sub("Connarus spermattetii"      ,"Connarus perrotettii"         ,x=g.s)
   g.s = sub("Cordia scrabida"            ,"Cordia exaltata"              ,x=g.s)
   g.s = sub("Gaulettia racemosa"         ,"Couepia racemosa"             ,x=g.s)
   g.s = sub("Coussarea racemosa"         ,"Coussarea albescens"          ,x=g.s)
   g.s = sub("Crepidospermum gondotiano"  ,"Crepidospermum goudotianum"   ,x=g.s)
   g.s = sub("Crepidospermum goudotiano"  ,"Crepidospermum goudotianum"   ,x=g.s)
   g.s = sub("Cupania hirta"              ,"Cupania hirsuta"              ,x=g.s)
   g.s = sub("Cybistax antisiphyllitica"  ,"Cybistax antisyphilitica"     ,x=g.s)
   g.s = sub("Dendrobrangea boliviana"    ,"Dendrobangia boliviana"       ,x=g.s)
   g.s = sub("Dialium guianensis"         ,"Dialium guianense"            ,x=g.s)
   g.s = sub("Diclinanonna matogrossensis","Diclinanona matogrossensis"   ,x=g.s)
   g.s = sub("Didymopanax vinosum"        ,"Schefflera vinosa"            ,x=g.s)
   g.s = sub("Didymopanax"                ,"Schefflera"                   ,x=g.s)
   g.s = sub("Diospyros praetermissa?"    ,"Diospyros vestita"            ,x=g.s)
   g.s = sub("Diospyros praetermissa"     ,"Diospyros vestita"            ,x=g.s)
   g.s = sub("Diploon venezuelano"        ,"Diploon cuspidatum"           ,x=g.s)
   g.s = sub("Dipteriyx odorata"          ,"Dipteryx odorata"             ,x=g.s)
   g.s = sub("Endlecheria"                ,"Endlicheria"                  ,x=g.s)
   g.s = sub("Enterolobium schobburgkii"  ,"Enterolobium schomburgkii"    ,x=g.s)
   g.s = sub("Enterolobium schombburgkii" ,"Enterolobium schomburgkii"    ,x=g.s)
   g.s = sub("Ephedrantus parviflorus"    ,"Ephedranthus parviflorus"     ,x=g.s)
   g.s = sub("Ephedrantus parvifolius"    ,"Ephedranthus parviflorus"     ,x=g.s)
   g.s = sub("Eriotreca globosa"          ,"Eriotheca globosa"            ,x=g.s)
   g.s = sub("Eschweilera amazomica"      ,"Eschweilera amazonica"        ,x=g.s)
   g.s = sub("Eschweilera apiculatum"     ,"Eschweilera apiculata"        ,x=g.s)
   g.s = sub("Eschweilera grandifolia"    ,"Eschweilera coriacea"         ,x=g.s)
   g.s = sub("Eschweilera idatimom"       ,"Lecythis idatimon"            ,x=g.s)
   g.s = sub("Eschweilera observa"        ,"Eschweilera obversa"          ,x=g.s)
   g.s = sub("Eschweilera pedicelata"     ,"Eschweilera pedicellata"      ,x=g.s)
   g.s = sub("Eugenia pacumensis"         ,"Austromyrtus ploumensis"      ,x=g.s)
   g.s = sub("Eugenia schumburgkii"       ,"Eugenia lambertiana"          ,x=g.s)
   g.s = sub("Ganua"                      ,"Madhuca"                      ,x=g.s)
   g.s = sub("Garcinia gardneriana"       ,"Garcinia brasiliensis"        ,x=g.s)
   g.s = sub("Geissospermum velosii"      ,"Geissospermum laeve"          ,x=g.s)
   g.s = sub("Geissospermum velozii"      ,"Geissospermum laeve"          ,x=g.s)
   g.s = sub("Glicoxilum"                 ,"Pouteria oppositifolia"       ,x=g.s)
   g.s = sub("Glycydendrom amazonicus"    ,"Glycydendron amazonicum"      ,x=g.s)
   g.s = sub("Glycydendron amazonicus"    ,"Glycydendron amazonicum"      ,x=g.s)
   g.s = sub("Guarea guianensis"          ,"Guarea"                       ,x=g.s)
   g.s = sub("Guarea subsessiliflora"     ,"Guarea macrophylla"           ,x=g.s)
   g.s = sub("Guatteria cardoniana"       ,"Guatteria recurvisepala"      ,x=g.s)
   g.s = sub("Hieronyma alcheornoides"    ,"Hieronyma alchorneoides"      ,x=g.s)
   g.s = sub("Hyeronima alcheornoides"    ,"Hieronyma alchorneoides"      ,x=g.s)
   g.s = sub("Hyeronima"                  ,"Hieronyma"                    ,x=g.s)
   g.s = sub("Hymenaea parviflora"        ,"Hymenaea parvifolia"          ,x=g.s)
   g.s = sub("Hymenolobium flavium"       ,"Hymenolobium flavum"          ,x=g.s)
   g.s = sub("Ilex parviflora"            ,"Ilex petiolaris"              ,x=g.s)
   g.s = sub("Inga aff\\."                ,"Inga affinis"                 ,x=g.s)
   g.s = sub("Inga jenmanii"              ,"Inga sertulifera"             ,x=g.s)
   g.s = sub("Inga dibaldiana"            ,"Inga thibaudiana"             ,x=g.s)
   g.s = sub("Inga nitida"                ,"Inga pilosula"                ,x=g.s)
   g.s = sub("Inga paraenses"             ,"Inga paraensis"               ,x=g.s)
   g.s = sub("Inga poliphylla"            ,"Inga"                         ,x=g.s)
   g.s = sub("Jacaraitia espinhosa"       ,"Jacaratia spinosa"            ,x=g.s)
   g.s = sub("Joannesia hevioides"        ,"Joannesia heveoides"          ,x=g.s)
   g.s = sub("Lacmelea aculeata"          ,"Lacmellea aculeata"           ,x=g.s)
   g.s = sub("Laetia procerra"            ,"Laetia procera"               ,x=g.s)
   g.s = sub("Lecythis praeclara"         ,"Eschweilera praealta"         ,x=g.s)
   g.s = sub("Licania densa"              ,"Licania densiflora"           ,x=g.s)
   g.s = sub("Licaria heteromorpha"       ,"Licania heteromorpha"         ,x=g.s)
   g.s = sub("Lucanarea"                  ,"Lacunaria"                    ,x=g.s)
   g.s = sub("Luheopsis duckeana"         ,"Lueheopsis duckeana"          ,x=g.s)
   g.s = sub("Lustema pubescei"           ,"Lacistema pubescens"          ,x=g.s)
   g.s = sub("Manilkara amazonica"        ,"Manilkara bidentata"          ,x=g.s)
   g.s = sub("Maluria"                    ,"Marlierea umbraticola"        ,x=g.s)
   g.s = sub("Maquira callophylla"        ,"Maquira calophylla"           ,x=g.s)
   g.s = sub("Marmaroxylon racemosum"     ,"Zygia racemosa"               ,x=g.s)
   g.s = sub("Mauriri chamissoana"        ,"Mouriri chamissoana"          ,x=g.s)
   g.s = sub("Maytenos guianensis"        ,"Maytenus guyanensis"          ,x=g.s)
   g.s = sub("Maytenus guianensis"        ,"Maytenus guyanensis"          ,x=g.s)
   g.s = sub("Meia maderensis"            ,"Neea"                         ,x=g.s)
   g.s = sub("Memora"                     ,"Adenocalymma"                 ,x=g.s)
   g.s = sub("Meterosideros"              ,"Metrosideros"                 ,x=g.s)
   g.s = sub("Mezelaurus"                 ,"Mezilaurus"                   ,x=g.s)
   g.s = sub("Mezelaurus itauba"          ,"Mezilaurus itauba"            ,x=g.s)
   g.s = sub("Miconia chrysophyllum"      ,"Miconia chrysophylla"         ,x=g.s)
   g.s = sub("Miconia surinamensis"       ,"Miconia poeppigii"            ,x=g.s)
   g.s = sub("Michopholis venulosa"       ,"Micropholis venulosa"         ,x=g.s)
   g.s = sub("Microphilis"                ,"Micropholis"                  ,x=g.s)
   g.s = sub("Micropholis guianensis"     ,"Micropholis guyanensis"       ,x=g.s)
   g.s = sub("Micropholis cuneata"        ,"Micropholis crassipedicellata",x=g.s)
   g.s = sub("Microphylis acutangula"     ,"Micropholis acutangula"       ,x=g.s)
   g.s = sub("Milletia"                   ,"Millettia"                    ,x=g.s)
   g.s = sub("Mimosa hostilis"            ,"Mimosa tenuiflora"            ,x=g.s)
   g.s = sub("Mouriri abnormis"           ,"Votomita guianensis"          ,x=g.s)
   g.s = sub("Myrciaria reticulata"       ,"Myrcia reticulata"            ,x=g.s)
   g.s = sub("Myrcia rutipula"            ,"Myrcia rufipila"              ,x=g.s)
   g.s = sub("Myrcia velutina"            ,"Myrcia"                       ,x=g.s)
   g.s = sub("Nephrolepsis"               ,"Nephrolepis"                  ,x=g.s)
   g.s = sub("Newtonia psilostachya"      ,"Pseudopiptadenia psilostachya",x=g.s)
   g.s = sub("Ocotea baturitensis"        ,"Ocotea"                       ,x=g.s)
   g.s = sub("Ocotea caudata"             ,"Ocotea cernua"                ,x=g.s)
   g.s = sub("Ocotea rubra"               ,"Sextonia rubra"               ,x=g.s)
   g.s = sub("Omedia perebea"             ,"Perebea mollis"               ,x=g.s)
   g.s = sub("Onichiopetalum amazonico"   ,"Onychopetalum amazonicum"     ,x=g.s)
   g.s = sub("Orbignya phalerata"         ,"Attalea speciosa"             ,x=g.s)
   g.s = sub("Pancovia"                   ,"Eurhynchium"                  ,x=g.s)
   g.s = sub("Parkia oppositifolia"       ,"Parkia nitida"                ,x=g.s)
   g.s = sub("Pausandra densiflora"       ,"Pausandra trianae"            ,x=g.s)
   g.s = sub("Peltogyne leicointei"       ,"Peltogyne lecointei"          ,x=g.s)
   g.s = sub("Piptadenia cobi"            ,"Stryphnodendron pulcherrimum" ,x=g.s)
   g.s = sub("Pouroma guianensis"         ,"Pourouma guianensis"          ,x=g.s)
   g.s = sub("Pouruma guianensis"         ,"Pourouma guianensis"          ,x=g.s)
   g.s = sub("Pourouma vilosa"            ,"Pourouma villosa"             ,x=g.s)
   g.s = sub("Pouteria ambelanifolia"     ,"Pouteria ambelaniifolia"      ,x=g.s)
   g.s = sub("Pouteria biloculares"       ,"Pouteria bilocularis"         ,x=g.s)
   g.s = sub("Pouteria filipis"           ,"Pouteria filipes"             ,x=g.s)
   g.s = sub("Pouteria lasiocarpa"        ,"Pouteria caimito"             ,x=g.s)
   g.s = sub("Pouteria gonggrijpii"       ,"Pouteria gongrijpii"          ,x=g.s)
   g.s = sub("Pouteria heterosepala"      ,"Pouteria polysepala"          ,x=g.s)
   g.s = sub("Pouteria macrophilla"       ,"Pouteria macrophylla"         ,x=g.s)
   g.s = sub("Pouteria paraensis"         ,"Pouteria macrocarpa"          ,x=g.s)
   g.s = sub("Pradosia praealta"          ,"Pradosia cochlearia"          ,x=g.s)
   g.s = sub("Protiumuceanum"             ,"Protium spruceanum"           ,x=g.s)
   g.s = sub("Protium heptafilum"         ,"Protium heptaphyllum"         ,x=g.s)
   g.s = sub("Protium hepytaphillum"      ,"Protium heptaphyllum"         ,x=g.s)
   g.s = sub("Protium pernervatum"        ,"Protium"                      ,x=g.s)
   g.s = sub("Prumes myrtifoliu"          ,"Prunus myrtifolia"            ,x=g.s)
   g.s = sub("Prunus myrtifolius"         ,"Prunus myrtifolia"            ,x=g.s)
   g.s = sub("Pseudolmedia murure"        ,"Pseudolmedia macrophylla"     ,x=g.s)
   g.s = sub("Ptecelobuim jucumba"        ,"Abarema jupunba"              ,x=g.s)
   g.s = sub("Pterocarpos rorhii"         ,"Pterocarpus rohrii"           ,x=g.s)
   g.s = sub("Pterocarpus amazonium"      ,"Pterocarpus santalinoides"    ,x=g.s)
   g.s = sub("Pterocarpus rhoire"         ,"Pterocarpus rohrii"           ,x=g.s)
   g.s = sub("Pterocarpus rhoiri"         ,"Pterocarpus rohrii"           ,x=g.s)
   g.s = sub("Ragala guianensis"          ,"Chrysophyllum sanguinolentum" ,x=g.s)
   g.s = sub("Rapanea ferruginea"         ,"Myrsine coriacea"             ,x=g.s)
   g.s = sub("Rapanea guianensis"         ,"Myrsine guianensis"           ,x=g.s)
   g.s = sub("Rheedia acuminata"          ,"Garcinia madruno"             ,x=g.s)
   g.s = sub("Rheedia gardneriana"        ,"Garcinia gardneriana"         ,x=g.s)
   g.s = sub("Rhodognaphalopsis"          ,"Pachira"                      ,x=g.s)
   g.s = sub("Rim de"                     ,"Crudia"                       ,x=g.s)
   g.s = sub("Rinorea pectino-squamata"   ,"Rinorea pectinosquamata"      ,x=g.s)
   g.s = sub("Rinoria guianensis"         ,"Rinorea guianensis"           ,x=g.s)
   g.s = sub("Rinorea passoura"           ,"Rinorea pubiflora"            ,x=g.s)
   g.s = sub("Rinoria racenosa"           ,"Rinorea racemosa"             ,x=g.s)
   g.s = sub("Rolinha esxuccar"           ,"Rollinia exsucca"             ,x=g.s)
   g.s = sub("Rollinia esxucca"           ,"Rollinia exsucca"             ,x=g.s)
   g.s = sub("Sacrogrotis guianensis"     ,"Sacoglottis guianensis"       ,x=g.s)
   g.s = sub("Salacea"                    ,"Salacia"                      ,x=g.s)
   g.s = sub("Salacia imprissifolia"      ,"Salacia impressifolia"        ,x=g.s)
   g.s = sub("Schefflera morototonii"     ,"Schefflera morototoni"        ,x=g.s)
   g.s = sub("Sclerolobium chrysophylum"  ,"Tachigali chrysophylla"       ,x=g.s)
   g.s = sub("Sclerolobium chrysophyllum" ,"Tachigali chrysophylla"       ,x=g.s)
   g.s = sub("Sclerolobium guianensis"    ,"Tachigali guianensis"         ,x=g.s)
   g.s = sub("Sclerolobium guianense"     ,"Tachigali guianensis"         ,x=g.s)
   g.s = sub("Simaba guianenses"          ,"Simaba guianensis"            ,x=g.s)
   g.s = sub("Simabacedron planch\\."     ,"Simaba cedron"                ,x=g.s)
   g.s = sub("Simarouba armara"           ,"Simarouba amara"              ,x=g.s)
   g.s = sub("Simaruba armara"            ,"Simarouba amara"              ,x=g.s)
   g.s = sub("Sipararuna"                 ,"Siparuna"                     ,x=g.s)
   g.s = sub("Stryphnodendron pulchrrimum","Stryphnodendron pulcherrimum" ,x=g.s)
   g.s = sub("Styrax ferrugineum"         ,"Styrax ferrugineus"           ,x=g.s)
   g.s = sub("Swartzia arborensis"        ,"Swartzia arborescens"         ,x=g.s)
   g.s = sub("Swartzia flamingii"         ,"Swartzia flaemingii"          ,x=g.s)
   g.s = sub("Swartizia microcarpum"      ,"Swartzia microcarpa"          ,x=g.s)
   g.s = sub("Swartzia retusa"            ,"Swartzia recurva"             ,x=g.s)
   g.s = sub("Swartzia viridiflora"       ,"Bocoa viridiflora"            ,x=g.s)
   g.s = sub("Tabebuia avellanedae"       ,"Handroanthus impetiginosus"   ,x=g.s)
   g.s = sub("Tabebuia barbata"           ,"Handroanthus barbatus"        ,x=g.s)
   g.s = sub("Tabebuia billbergii"        ,"Handroanthus billbergii"      ,x=g.s)
   g.s = sub("Tabebuia capitata"          ,"Handroanthus capitatus"       ,x=g.s)
   g.s = sub("Tabebuia chrysantha"        ,"Handroanthus chrysanthus"     ,x=g.s)
   g.s = sub("Tabebuia chrysotricha"      ,"Handroanthus chrysotrichus"   ,x=g.s)
   g.s = sub("Tabebuia donnell-smithii"   ,"Roseodendron donnell-smithii" ,x=g.s)
   g.s = sub("Tabebuia guayacan"          ,"Handroanthus guayacan"        ,x=g.s)
   g.s = sub("Tabebuia heptaphylla"       ,"Handroanthus heptaphyllus"    ,x=g.s)
   g.s = sub("Tabebuia impetiginosa"      ,"Handroanthus impetiginosus"   ,x=g.s)
   g.s = sub("Tabebuia incana"            ,"Handroanthus incanus"         ,x=g.s)
   g.s = sub("Tabebuia lapacho"           ,"Handroanthus lapacho"         ,x=g.s)
   g.s = sub("Tabebuia obscura"           ,"Handroanthus obscurus"        ,x=g.s)
   g.s = sub("Tabebuia pedicellata"       ,"Handroanthus pedicellatus"    ,x=g.s)
   g.s = sub("Tabebuia serratifolia"      ,"Handroanthus serratifolius"   ,x=g.s)
   g.s = sub("Tabebuia vellosoi"          ,"Handroanthus vellosoi"        ,x=g.s)
   g.s = sub("Tachigalia alba"            ,"Tachigali paniculata"         ,x=g.s)
   g.s = sub("Tachigali alba"             ,"Tachigali paniculata"         ,x=g.s)
   g.s = sub("Tachigalia myrmecophila"    ,"Tachigali myrmecophila"       ,x=g.s)
   g.s = sub("Tachigalia paniculata"      ,"Tachigali paniculata"         ,x=g.s)
   g.s = sub("Tapirira marchandi"         ,"Tapirira obtusa"              ,x=g.s)
   g.s = sub("Tapirira marchandii"        ,"Tapirira obtusa"              ,x=g.s)
   g.s = sub("Tapirira myriantha"         ,"Tapirira guianensis"          ,x=g.s)
   g.s = sub("Tapirira peckoltiana"       ,"Tapirira obtusa"              ,x=g.s)
   g.s = sub("Terminalea argentea"        ,"Terminalia argentea"          ,x=g.s)
   g.s = sub("Terminalia amazonica"       ,"Terminalia amazonia"          ,x=g.s)
   g.s = sub("Theobroma subincaum"        ,"Theobroma subincanum"         ,x=g.s)
   g.s = sub("Thyrsodium paraense"        ,"Thyrsodium spruceanum"        ,x=g.s)
   g.s = sub("Trattinickia laurencei"     ,"Trattinnickia lawrancei"      ,x=g.s)
   g.s = sub("Trattinickia laurense"      ,"Trattinnickia lawrancei"      ,x=g.s)
   g.s = sub("Trattinickia rhoifolia"     ,"Trattinnickia rhoifolia"      ,x=g.s)
   g.s = sub("Trattinickia"               ,"Trattinnickia"                ,x=g.s)
   g.s = sub("Trattinnickia burseraefolia","Trattinnickia burserifolia"   ,x=g.s)
   g.s = sub("Trichilia lequente"         ,"Trichilia lecointei"          ,x=g.s)
   g.s = sub("Triquilia lequente"         ,"Trichilia lecointei"          ,x=g.s)
   g.s = sub("Trichilia sipo"             ,"Trichilia cipo"               ,x=g.s)
   g.s = sub("Trichillia"                 ,"Trichilia"                    ,x=g.s)
   g.s = sub("Trymatococcus amazonicum"   ,"Trymatococcus amazonicus"     ,x=g.s)
   g.s = sub("Trymatococcus paraensis"    ,"Trymatococcus amazonicus"     ,x=g.s)
   g.s = sub("Vataicreopesis speciosa"    ,"Vataireopsis speciosa"        ,x=g.s)
   g.s = sub("Vataireopesis speciosa"     ,"Vataireopsis speciosa"        ,x=g.s)
   g.s = sub("Vernonia diffusa"           ,"Vernonia erigeroides"         ,x=g.s)
   g.s = sub("Vimia guianensis"           ,"Vismia guianensis"            ,x=g.s)
   g.s = sub("Virola crebinervia"         ,"Virola crebrinervia"          ,x=g.s)
   g.s = sub("Virola melinoni$"           ,"Virola michelii"              ,x=g.s)
   g.s = sub("Virola melinonii"           ,"Virola michelii"              ,x=g.s)
   g.s = sub("Virola melionii"            ,"Virola michelii"              ,x=g.s)
   g.s = sub("Virola melioni"             ,"Virola michelii"              ,x=g.s)
   g.s = sub("Virola michelli"            ,"Virola michelii"              ,x=g.s)
   g.s = sub("Vochisia surinamensis"      ,"Vochysia surinamensis"        ,x=g.s)
   g.s = sub("Vochysia micrantha"         ,"Vochysia"                     ,x=g.s)
   g.s = sub("Zizyphus itacaiunensis"     ,"Ziziphus"                     ,x=g.s)
   #---------------------------------------------------------------------------------------#



   #----- Break again into genus and species. ---------------------------------------------#
   g.s = keep.gen.spe.only(x=g.s)
   #---------------------------------------------------------------------------------------#



   #----- Remove family names from scientific columns. ------------------------------------#
   bye                 = g.s %in% c("Anacardiaceae","Myrtaceae","Moraceae","Ind"
                                   ,"Lauraceae","Nyctaginaceae")
   g.s[bye]            = NA_character_
   liana               = g.s %in% c("Liane")
   g.s[liana]          = "Liana"
   #---------------------------------------------------------------------------------------#


   #------ Find genus. --------------------------------------------------------------------#
   g   = keep.gen.spe.only(x=g.s,out="genus")
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #      Copy the data back to the structure.                                             #
   #---------------------------------------------------------------------------------------#
   dat$scientific = g.s
   dat$genus      = g
   #---------------------------------------------------------------------------------------#

   return(dat)
}#end function standard.common.name
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      This is the most up-to-date table to find families given the genus.  This contains  #
# all genera that occur at TNF and Paracou, feel free to add more.  In case the family is  #
# known but the genus is not, we add a unique genus name that is non-informative.          #
#------------------------------------------------------------------------------------------#
standard.family.name <<- function(datum){
   #----- Make sure families are properly capitalised. ------------------------------------#
   if ("family" %in% datum){
      datum$family = capwords(datum$family,strict=TRUE)
   }else{
      datum$family = rep("Ignotaceae",times=nrow(datum))
   }#end if ("family" %in% datum)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Substitute obsolete/misspelt families that are straightforward.                   #
   #---------------------------------------------------------------------------------------#
   datum$family = sub("Bombacaceae"                 ,"Malvaceae"     ,x=datum$family)
   datum$family = sub("Caesalpinaceae"              ,"Fabaceae"      ,x=datum$family)
   datum$family = sub("Cecropiaceae"                ,"Urticaceae"    ,x=datum$family)
   datum$family = sub("Elacocarpaceae"              ,"Elaeocarpaceae",x=datum$family)
   datum$family = sub("Guttiferae"                  ,"Clusiaceae"    ,x=datum$family)
   datum$family = sub("Hippocrateaceae"             ,"Celastraceae"  ,x=datum$family)
   datum$family = sub("Leguminosae"                 ,"Fabaceae"      ,x=datum$family)
   datum$family = sub("Leguminosae-Mimosoideae"     ,"Fabaceae"      ,x=datum$family)
   datum$family = sub("Leguminosae-Caesalpinioideae","Fabaceae"      ,x=datum$family)
   datum$family = sub("Leguminosae-Papilionoideae"  ,"Fabaceae"      ,x=datum$family)
   datum$family = sub("Quiinaceae"                  ,"Ochnaceae"     ,x=datum$family)
   datum$family = sub("Ramnaceae"                   ,"Rhamnaceae"    ,x=datum$family)
   datum$family = sub("Mimosaceae"                  ,"Fabaceae"      ,x=datum$family)
   datum$family = sub("Papilionaceae"               ,"Fabaceae"      ,x=datum$family)
   datum$family = sub("Tiliaceae"                   ,"Malvaceae"     ,x=datum$family)
   datum$family = sub("Sterculiaceae"               ,"Malvaceae"     ,x=datum$family)
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Look-up table for all families.                                                   #
   #---------------------------------------------------------------------------------------#
   n=0  ; g2f      = list()
   n=n+1; g2f[[n]] = list( genus = "Abarema"           , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Abuta"             , family = "Menispermaceae"   )
   n=n+1; g2f[[n]] = list( genus = "Acacia"            , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Acalypha"          , family = "Euphorbiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Acer"              , family = "Sapindaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Acmanthera"        , family = "Malpighiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Acosmium"          , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Acrocomia"         , family = "Arecaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Adenanthos"        , family = "Proteaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Adenocalymma"      , family = "Bignoniaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Adenocarpus"       , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Adenophaedra"      , family = "Euphorbiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Adenostoma"        , family = "Rosaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Adinandra"         , family = "Pentaphylacaceae" )
   n=n+1; g2f[[n]] = list( genus = "Aegiphila"         , family = "Lamiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Aegle"             , family = "Rutaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Aesculus"          , family = "Sapindaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Afrocarpus"        , family = "Podocarpaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Agonandra"         , family = "Opiliaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Aiouea"            , family = "Lauraceae"        )
   n=n+1; g2f[[n]] = list( genus = "Albizia"           , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Alchornea"         , family = "Euphorbiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Alchorneopsis"     , family = "Euphorbiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Aldina"            , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Alexa"             , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Alibertia"         , family = "Rubiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Allantoma"         , family = "Lecythidaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Allocasuarina"     , family = "Casuarinaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Allophylus"        , family = "Sapindaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Alnus"             , family = "Betulaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Alseis"            , family = "Rubiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Alternanthera"     , family = "Amaranthaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Amaioua"           , family = "Rubiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Amaranthus"        , family = "Amaranthaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Ambelania"         , family = "Apocynaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Amburana"          , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Amelanchier"       , family = "Rosaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Ampelocera"        , family = "Ulmaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Amphiodon"         , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Amphirrhox"        , family = "Violaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Amphitecna"        , family = "Bignoniaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Anacardium"        , family = "Anacardiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Anadenanthera"     , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Anagyris"          , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Anaxagorea"        , family = "Annonaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Andersonia"        , family = "Ericaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Andira"            , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Aniba"             , family = "Lauraceae"        )
   n=n+1; g2f[[n]] = list( genus = "Anigozanthos"      , family = "Haemodoraceae"    )
   n=n+1; g2f[[n]] = list( genus = "Anisophyllea"      , family = "Anisophylleaceae" )
   n=n+1; g2f[[n]] = list( genus = "Annona"            , family = "Annonaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Anogeissus"        , family = "Combretaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Anomalocalyx"      , family = "Euphorbiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Anthonotha"        , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Antirhea"          , family = "Rubiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Antonia"           , family = "Loganiaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Aparisthmium"      , family = "Euphorbiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Apeiba"            , family = "Malvaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Aphanamixis"       , family = "Meliaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Apodytes"          , family = "Icacinaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Apollonias"        , family = "Lauraceae"        )
   n=n+1; g2f[[n]] = list( genus = "Aporosa"           , family = "Phyllanthaceae"   )
   n=n+1; g2f[[n]] = list( genus = "Aptandra"          , family = "Olacaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Apuleia"           , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Arbutus"           , family = "Ericaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Arctostaphylos"    , family = "Ericaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Artemisia"         , family = "Asteraceae"       )
   n=n+1; g2f[[n]] = list( genus = "Artocarpus"        , family = "Moraceae"         )
   n=n+1; g2f[[n]] = list( genus = "Aspidosperma"      , family = "Apocynaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Astrocaryum"       , family = "Arecaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Astroloma"         , family = "Ericaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Astronium"         , family = "Anacardiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Atriplex"          , family = "Amaranthaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Attalea"           , family = "Arecaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Aulax"             , family = "Proteaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Austromyrtus"      , family = "Myrtaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Avicennia"         , family = "Acanthaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Azadirachta"       , family = "Meliaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Baccaurea"         , family = "Phyllanthaceae"   )
   n=n+1; g2f[[n]] = list( genus = "Baccharis"         , family = "Asteraceae"       )
   n=n+1; g2f[[n]] = list( genus = "Bactris"           , family = "Arecaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Baeckea"           , family = "Myrtaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Bagassa"           , family = "Moraceae"         )
   n=n+1; g2f[[n]] = list( genus = "Balfourodendron"   , family = "Rutaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Balizia"           , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Bambusa"           , family = "Poaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Banara"            , family = "Salicaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Banksia"           , family = "Proteaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Barteria"          , family = "Passifloraceae"   )
   n=n+1; g2f[[n]] = list( genus = "Batesia"           , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Batocarpus"        , family = "Moraceae"         )
   n=n+1; g2f[[n]] = list( genus = "Bauhinia"          , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Bellucia"          , family = "Melastomataceae"  )
   n=n+1; g2f[[n]] = list( genus = "Berberis"          , family = "Berberidaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Berlinia"          , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Bertholletia"      , family = "Lecythidaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Bertya"            , family = "Euphorbiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Beta"              , family = "Amaranthaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Betula"            , family = "Betulaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Beyeria"           , family = "Euphorbiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Bixa"              , family = "Bixaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Blakea"            , family = "Melastomataceae"  )
   n=n+1; g2f[[n]] = list( genus = "Blepharocalyx"     , family = "Myrtaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Bocageopsis"       , family = "Annonaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Bocoa"             , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Boehmeria"         , family = "Urticaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Bolusanthus"       , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Bombax"            , family = "Malvaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Bonellia"          , family = "Primulaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Boronia"           , family = "Rutaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Bossiaea"          , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Boswellia"         , family = "Burseraceae"      )
   n=n+1; g2f[[n]] = list( genus = "Botryarrhena"      , family = "Rubiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Bougainvillea"     , family = "Nyctaginaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Bowdichia"         , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Brachychiton"      , family = "Malvaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Brachystegia"      , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Bridelia"          , family = "Phyllanthaceae"   )
   n=n+1; g2f[[n]] = list( genus = "Brosimum"          , family = "Moraceae"         )
   n=n+1; g2f[[n]] = list( genus = "Bruguiera"         , family = "Rhizophoraceae"   )
   n=n+1; g2f[[n]] = list( genus = "Buchanania"        , family = "Anacardiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Buchenavia"        , family = "Combretaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Buddleja"          , family = "Scrophulariaceae" )
   n=n+1; g2f[[n]] = list( genus = "Bulnesia"          , family = "Zygophyllaceae"   )
   n=n+1; g2f[[n]] = list( genus = "Bupleurum"         , family = "Apiaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Burkea"            , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Bursera"           , family = "Burseraceae"      )
   n=n+1; g2f[[n]] = list( genus = "Butea"             , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Buxus"             , family = "Buxaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Byrsonima"         , family = "Malpighiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Caesalpinia"       , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Calathea"          , family = "Marantaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Calliandra"        , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Callitris"         , family = "Cupressaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Calluna"           , family = "Ericaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Calophyllum"       , family = "Calophyllaceae"   )
   n=n+1; g2f[[n]] = list( genus = "Calothamnus"       , family = "Myrtaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Calycanthus"       , family = "Calycanthaceae"   )
   n=n+1; g2f[[n]] = list( genus = "Calycophyllum"     , family = "Rubiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Calyptranthes"     , family = "Myrtaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Calytrix"          , family = "Myrtaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Campomanesia"      , family = "Myrtaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Capirona"          , family = "Rubiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Capparidastrum"    , family = "Capparaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Capparis"          , family = "Capparaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Caraipa"           , family = "Calophyllaceae"   )
   n=n+1; g2f[[n]] = list( genus = "Carapa"            , family = "Meliaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Carex"             , family = "Cyperaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Cariniana"         , family = "Lecythidaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Carissa"           , family = "Apocynaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Carpinus"          , family = "Betulaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Caryocar"          , family = "Caryocaraceae"    )
   n=n+1; g2f[[n]] = list( genus = "Cascabela"         , family = "Apocynaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Casearia"          , family = "Salicaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Cassia"            , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Cassinia"          , family = "Asteraceae"       )
   n=n+1; g2f[[n]] = list( genus = "Cassipourea"       , family = "Rhizophoraceae"   )
   n=n+1; g2f[[n]] = list( genus = "Castanea"          , family = "Fagaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Castilla"          , family = "Moraceae"         )
   n=n+1; g2f[[n]] = list( genus = "Catostemma"        , family = "Malvaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Cavanillesia"      , family = "Malvaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Ceanothus"         , family = "Rhamnaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Cecropia"          , family = "Urticaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Cedrela"           , family = "Meliaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Cedrelinga"        , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Ceiba"             , family = "Malvaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Celtis"            , family = "Cannabaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Cephalaria"        , family = "Caprifoliaceae"   )
   n=n+1; g2f[[n]] = list( genus = "Ceratonia"         , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Cerbera"           , family = "Apocynaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Cercocarpus"       , family = "Rosaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Cereus"            , family = "Cactaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Ceriops"           , family = "Rhizophoraceae"   )
   n=n+1; g2f[[n]] = list( genus = "Chaetachme"        , family = "Ulmaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Chaetocarpus"      , family = "Euphorbiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Chaetochlamys"     , family = "Acanthaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Chamaecrista"      , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Chamaedorea"       , family = "Arecaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Chaunochiton"      , family = "Olacaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Cheiloclinium"     , family = "Celastraceae"     )
   n=n+1; g2f[[n]] = list( genus = "Chenopodium"       , family = "Amaranthaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Chimarrhis"        , family = "Rubiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Chionanthus"       , family = "Oleaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Chloroleucon"      , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Chromolucuma"      , family = "Sapotaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Chrysobalanus"     , family = "Chrysobalanaceae" )
   n=n+1; g2f[[n]] = list( genus = "Chrysocephalum"    , family = "Asteraceae"       )
   n=n+1; g2f[[n]] = list( genus = "Chrysophyllum"     , family = "Sapotaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Cichorium"         , family = "Asteraceae"       )
   n=n+1; g2f[[n]] = list( genus = "Cinnamomum"        , family = "Lauraceae"        )
   n=n+1; g2f[[n]] = list( genus = "Cirsium"           , family = "Asteraceae"       )
   n=n+1; g2f[[n]] = list( genus = "Cissus"            , family = "Vitaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Cistus"            , family = "Cistaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Citharexylum"      , family = "Verbenaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Citrus"            , family = "Rutaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Cladostemon"       , family = "Capparaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Claoxylon"         , family = "Euphorbiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Clarisia"          , family = "Moraceae"         )
   n=n+1; g2f[[n]] = list( genus = "Clausena"          , family = "Rutaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Cleistanthus"      , family = "Phyllanthaceae"   )
   n=n+1; g2f[[n]] = list( genus = "Clematis"          , family = "Ranunculaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Clethra"           , family = "Clethraceae"      )
   n=n+1; g2f[[n]] = list( genus = "Clidemia"          , family = "Melastomataceae"  )
   n=n+1; g2f[[n]] = list( genus = "Clitoria"          , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Clusia"            , family = "Clusiaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Cneorum"           , family = "Rutaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Cnidoscolus"       , family = "Euphorbiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Coccoloba"         , family = "Polygonaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Cochlospermum"     , family = "Bixaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Cola"              , family = "Malvaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Colubrina"         , family = "Rhamnaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Combretum"         , family = "Combretaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Commersonia"       , family = "Malvaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Commiphora"        , family = "Burseraceae"      )
   n=n+1; g2f[[n]] = list( genus = "Conceveiba"        , family = "Euphorbiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Condalia"          , family = "Rhamnaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Connarus"          , family = "Connaraceae"      )
   n=n+1; g2f[[n]] = list( genus = "Conostephium"      , family = "Ericaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Convolvulus"       , family = "Convolvulaceae"   )
   n=n+1; g2f[[n]] = list( genus = "Copaiba"           , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Copaifera"         , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Cordia"            , family = "Boraginaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Cordiera"          , family = "Rubiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Cornus"            , family = "Cornaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Corylus"           , family = "Betulaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Corymbia"          , family = "Myrtaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Corythophora"      , family = "Lecythidaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Couepia"           , family = "Chrysobalanaceae" )
   n=n+1; g2f[[n]] = list( genus = "Coula"             , family = "Olacaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Couma"             , family = "Apocynaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Couratari"         , family = "Lecythidaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Coursetia"         , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Coussapoa"         , family = "Urticaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Coussarea"         , family = "Rubiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Coutarea"          , family = "Rubiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Crataegus"         , family = "Rosaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Crepidospermum"    , family = "Burseraceae"      )
   n=n+1; g2f[[n]] = list( genus = "Crepis"            , family = "Asteraceae"       )
   n=n+1; g2f[[n]] = list( genus = "Crescentia"        , family = "Bignoniaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Croton"            , family = "Euphorbiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Crudia"            , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Cupania"           , family = "Sapindaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Curtisia"          , family = "Curtisiaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Cussonia"          , family = "Araliaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Cyathocalyx"       , family = "Annonaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Cybianthus"        , family = "Primulaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Cybistax"          , family = "Bignoniaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Cycas"             , family = "Cycadaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Cymbopetalum"      , family = "Annonaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Cynometra"         , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Cynophalla"        , family = "Capparaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Cyrilla"           , family = "Cyrillaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Cytisus"           , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Dacrydium"         , family = "Podocarpaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Dacryodes"         , family = "Burseraceae"      )
   n=n+1; g2f[[n]] = list( genus = "Dalbergia"         , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Daphne"            , family = "Thymelaeaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Datura"            , family = "Solanaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Davilla"           , family = "Dilleniaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Deidamia"          , family = "Passifloraceae"   )
   n=n+1; g2f[[n]] = list( genus = "Delphinium"        , family = "Ranunculaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Dendrobangia"      , family = "Cardiopteridaceae")
   n=n+1; g2f[[n]] = list( genus = "Dendropanax"       , family = "Araliaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Deschampsia"       , family = "Poaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Desmoncus"         , family = "Arecaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Desmopsis"         , family = "Annonaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Dialium"           , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Dianthus"          , family = "Caryophyllaceae"  )
   n=n+1; g2f[[n]] = list( genus = "Diastella"         , family = "Proteaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Dichostemma"       , family = "Euphorbiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Diclinanona"       , family = "Annonaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Dicorynia"         , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Dicranopteris"     , family = "Gleicheniaceae"   )
   n=n+1; g2f[[n]] = list( genus = "Dicranostyles"     , family = "Convolvulaceae"   )
   n=n+1; g2f[[n]] = list( genus = "Dictyocaryum"      , family = "Arecaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Dicypellium"       , family = "Lauraceae"        )
   n=n+1; g2f[[n]] = list( genus = "Dillenia"          , family = "Dilleniaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Dimorphandra"      , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Dinizia"           , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Diospyros"         , family = "Ebenaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Diplospora"        , family = "Rubiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Diplotropis"       , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Dipterocarpus"     , family = "Dipterocarpaceae" )
   n=n+1; g2f[[n]] = list( genus = "Dipteryx"          , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Dirca"             , family = "Thymelaeaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Discophora"        , family = "Stemonuraceae"    )
   n=n+1; g2f[[n]] = list( genus = "Dittrichia"        , family = "Asteraceae"       )
   n=n+1; g2f[[n]] = list( genus = "Dodecastigma"      , family = "Euphorbiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Dodonaea"          , family = "Sapindaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Doliocarpus"       , family = "Dilleniaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Dombeya"           , family = "Malvaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Dorycnium"         , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Dovyalis"          , family = "Salicaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Drimia"            , family = "Asparagaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Drimys"            , family = "Winteraceae"      )
   n=n+1; g2f[[n]] = list( genus = "Dryas"             , family = "Rosaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Drypetes"          , family = "Putranjivaceae"   )
   n=n+1; g2f[[n]] = list( genus = "Duckeodendron"     , family = "Solanaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Duckesia"          , family = "Humiriaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Duguetia"          , family = "Annonaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Dulacia"           , family = "Olacaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Duroia"            , family = "Rubiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Dussia"            , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Dystovomita"       , family = "Clusiaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Ecclinusa"         , family = "Sapotaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Efulensia"         , family = "Passifloraceae"   )
   n=n+1; g2f[[n]] = list( genus = "Ekebergia"         , family = "Meliaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Elaeagnus"         , family = "Elaeagnaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Elaeis"            , family = "Arecaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Elaeoluma"         , family = "Sapotaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Elateriospermum"   , family = "Euphorbiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Elizabetha"        , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Embothrium"        , family = "Proteaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Emmotum"           , family = "Emmotaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Empetrum"          , family = "Ericaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Endlicheria"       , family = "Lauraceae"        )
   n=n+1; g2f[[n]] = list( genus = "Endopleura"        , family = "Humiriaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Enterolobium"      , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Eperua"            , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Ephedranthus"      , family = "Annonaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Eremaea"           , family = "Myrtaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Eremophila"        , family = "Scrophulariaceae" )
   n=n+1; g2f[[n]] = list( genus = "Erica"             , family = "Ericaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Eriodictyon"       , family = "Boraginaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Eriostemon"        , family = "Rutaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Eriotheca"         , family = "Malvaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Erisma"            , family = "Vochysiaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Erythrina"         , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Erythrophleum"     , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Erythroxylum"      , family = "Erythroxylaceae"  )
   n=n+1; g2f[[n]] = list( genus = "Eschweilera"       , family = "Lecythidaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Esenbeckia"        , family = "Rutaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Eubotrys"          , family = "Ericaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Eucalyptus"        , family = "Myrtaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Euclea"            , family = "Ebenaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Eugenia"           , family = "Myrtaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Euonymus"          , family = "Celastraceae"     )
   n=n+1; g2f[[n]] = list( genus = "Euphorbia"         , family = "Euphorbiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Eurhynchium"       , family = "Brachytheciaceae" )
   n=n+1; g2f[[n]] = list( genus = "Eurycoma"          , family = "Simaroubaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Eutaxia"           , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Euterpe"           , family = "Arecaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Euxylophora"       , family = "Rutaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Fagraea"           , family = "Gentianaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Fagus"             , family = "Fagaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Faramea"           , family = "Rubiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Faurea"            , family = "Proteaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Ferdinandusa"      , family = "Rubiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Ficus"             , family = "Moraceae"         )
   n=n+1; g2f[[n]] = list( genus = "Firmiana"          , family = "Malvaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Frangula"          , family = "Rhamnaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Fraunhofera"       , family = "Celastraceae"     )
   n=n+1; g2f[[n]] = list( genus = "Fraxinus"          , family = "Oleaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Funtumia"          , family = "Apocynaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Fusaea"            , family = "Annonaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Galipea"           , family = "Rutaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Gallesia"          , family = "Phytolaccaceae"   )
   n=n+1; g2f[[n]] = list( genus = "Garcinia"          , family = "Clusiaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Gardenia"          , family = "Rubiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Gaulettia"         , family = "Chrysobalanaceae" )
   n=n+1; g2f[[n]] = list( genus = "Geijera"           , family = "Rutaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Geissospermum"     , family = "Apocynaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Genipa"            , family = "Rubiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Genista"           , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Gilbertiodendron"  , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Gladiolus"         , family = "Iridaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Globularia"        , family = "Plantaginaceae"   )
   n=n+1; g2f[[n]] = list( genus = "Glycydendron"      , family = "Euphorbiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Gochnatia"         , family = "Asteraceae"       )
   n=n+1; g2f[[n]] = list( genus = "Gompholobium"      , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Goupia"            , family = "Goupiaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Grevillea"         , family = "Proteaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Grewia"            , family = "Malvaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Guamia"            , family = "Annonaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Guapira"           , family = "Nyctaginaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Guarea"            , family = "Meliaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Guatteria"         , family = "Annonaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Guazuma"           , family = "Malvaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Gustavia"          , family = "Lecythidaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Gutierrezia"       , family = "Asteraceae"       )
   n=n+1; g2f[[n]] = list( genus = "Hakea"             , family = "Proteaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Haldina"           , family = "Rubiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Halimium"          , family = "Cistaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Hamelia"           , family = "Rubiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Hampea"            , family = "Malvaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Handroanthus"      , family = "Bignoniaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Hardwickia"        , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Harpephyllum"      , family = "Anacardiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Hebe"              , family = "Plantaginaceae"   )
   n=n+1; g2f[[n]] = list( genus = "Hebepetalum"       , family = "Linaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Heberdenia"        , family = "Primulaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Hedera"            , family = "Araliaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Hedyosmum"         , family = "Chloranthaceae"   )
   n=n+1; g2f[[n]] = list( genus = "Heisteria"         , family = "Olacaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Helianthemum"      , family = "Cistaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Helianthostylis"   , family = "Moraceae"         )
   n=n+1; g2f[[n]] = list( genus = "Helicostylis"      , family = "Moraceae"         )
   n=n+1; g2f[[n]] = list( genus = "Helictotrichon"    , family = "Poaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Heliocarpus"       , family = "Malvaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Helleborus"        , family = "Ranunculaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Henriettea"        , family = "Melastomataceae"  )
   n=n+1; g2f[[n]] = list( genus = "Henriettella"      , family = "Melastomataceae"  )
   n=n+1; g2f[[n]] = list( genus = "Heteropogon"       , family = "Poaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Hevea"             , family = "Euphorbiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Hibbertia"         , family = "Dilleniaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Hibiscus"          , family = "Malvaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Hieronyma"         , family = "Phyllanthaceae"   )
   n=n+1; g2f[[n]] = list( genus = "Himatanthus"       , family = "Apocynaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Hippobromus"       , family = "Sapindaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Hippocratea"       , family = "Celastraceae"     )
   n=n+1; g2f[[n]] = list( genus = "Hirtella"          , family = "Chrysobalanaceae" )
   n=n+1; g2f[[n]] = list( genus = "Holodiscus"        , family = "Rosaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Homalium"          , family = "Salicaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Homalium"          , family = "Salicaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Humiria"           , family = "Humiriaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Humiriastrum"      , family = "Humiriaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Hura"              , family = "Euphorbiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Hymenaea"          , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Hymenolobium"      , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Hyperacanthus"     , family = "Rubiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Hypericum"         , family = "Hypericaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Hypoestes"         , family = "Acanthaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Hypolaena"         , family = "Restionaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Ignotum"           , family = "Ignotaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Ilex"              , family = "Aquifoliaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Inga"              , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Inhambanella"      , family = "Sapotaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Ipomoea"           , family = "Convolvulaceae"   )
   n=n+1; g2f[[n]] = list( genus = "Iryanthera"        , family = "Myristicaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Isertia"           , family = "Rubiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Isoglossa"         , family = "Acanthaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Ixonanthes"        , family = "Ixonanthaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Ixora"             , family = "Rubiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Jacaranda"         , family = "Bignoniaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Jacaratia"         , family = "Caricaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Jacksonia"         , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Jacquinia"         , family = "Primulaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Jasminum"          , family = "Oleaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Jatropha"          , family = "Euphorbiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Joannesia"         , family = "Euphorbiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Juglans"           , family = "Juglandaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Julbernardia"      , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Juniperus"         , family = "Cupressaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Justicia"          , family = "Acanthaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Kielmeyera"        , family = "Calophyllaceae"   )
   n=n+1; g2f[[n]] = list( genus = "Kigelia"           , family = "Bignoniaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Laburnum"          , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Lacistema"         , family = "Lacistemataceae"  )
   n=n+1; g2f[[n]] = list( genus = "Lacmellea"         , family = "Apocynaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Lacunaria"         , family = "Ochnaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Ladenbergia"       , family = "Rubiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Laetia"            , family = "Salicaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Lafoensia"         , family = "Lythraceae"       )
   n=n+1; g2f[[n]] = list( genus = "Lagerstroemia"     , family = "Lythraceae"       )
   n=n+1; g2f[[n]] = list( genus = "Lannea"            , family = "Anacardiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Laportea"          , family = "Urticaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Larix"             , family = "Pinaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Larrea"            , family = "Zygophyllaceae"   )
   n=n+1; g2f[[n]] = list( genus = "Laurus"            , family = "Lauraceae"        )
   n=n+1; g2f[[n]] = list( genus = "Lavandula"         , family = "Lamiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Lecythis"          , family = "Lecythidaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Leonia"            , family = "Violaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Lepechinia"        , family = "Lamiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Lepidosperma"      , family = "Cyperaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Leptaulus"         , family = "Cardiopteridaceae")
   n=n+1; g2f[[n]] = list( genus = "Leptospermum"      , family = "Myrtaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Leucadendron"      , family = "Proteaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Leucospermum"      , family = "Proteaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Liana"             , family = "Lianaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Librevillea"       , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Licania"           , family = "Chrysobalanaceae" )
   n=n+1; g2f[[n]] = list( genus = "Licaria"           , family = "Lauraceae"        )
   n=n+1; g2f[[n]] = list( genus = "Ligusticum"        , family = "Apiaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Ligustrum"         , family = "Oleaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Limonium"          , family = "Plumbaginaceae"   )
   n=n+1; g2f[[n]] = list( genus = "Lindackeria"       , family = "Achariaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Lippia"            , family = "Verbenaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Liquidambar"       , family = "Altingiaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Litchi"            , family = "Sapindaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Livistona"         , family = "Arecaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Lonchocarpus"      , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Lonicera"          , family = "Caprifoliaceae"   )
   n=n+1; g2f[[n]] = list( genus = "Lophira"           , family = "Ochnaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Lophostemon"       , family = "Myrtaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Luehea"            , family = "Malvaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Lueheopsis"        , family = "Malvaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Lumnitzera"        , family = "Combretaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Lunania"           , family = "Salicaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Lyginia"           , family = "Anarthriaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Lyonia"            , family = "Ericaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Lysiphyllum"       , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Mabea"             , family = "Euphorbiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Macairea"          , family = "Melastomataceae"  )
   n=n+1; g2f[[n]] = list( genus = "Macaranga"         , family = "Euphorbiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Machaerium"        , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Macoubea"          , family = "Apocynaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Macrocnemum"       , family = "Rubiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Macrolobium"       , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Macrozamia"        , family = "Zamiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Madhuca"           , family = "Sapotaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Magnolia"          , family = "Magnoliaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Mallotus"          , family = "Euphorbiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Malouetia"         , family = "Apocynaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Malus"             , family = "Rosaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Malva"             , family = "Malvaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Mammea"            , family = "Calophyllaceae"   )
   n=n+1; g2f[[n]] = list( genus = "Mangifera"         , family = "Anacardiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Manicaria"         , family = "Arecaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Manihot"           , family = "Euphorbiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Manilkara"         , family = "Sapotaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Mansoa"            , family = "Bignoniaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Maprounea"         , family = "Euphorbiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Maquira"           , family = "Moraceae"         )
   n=n+1; g2f[[n]] = list( genus = "Marcgravia"        , family = "Marcgraviaceae"   )
   n=n+1; g2f[[n]] = list( genus = "Margaritaria"      , family = "Phyllanthaceae"   )
   n=n+1; g2f[[n]] = list( genus = "Markhamia"         , family = "Bignoniaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Marlierea"         , family = "Myrtaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Matayba"           , family = "Sapindaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Matisia"           , family = "Malvaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Mauritia"          , family = "Arecaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Mauritiella"       , family = "Arecaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Mayna"             , family = "Achariaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Maytenus"          , family = "Celastraceae"     )
   n=n+1; g2f[[n]] = list( genus = "Meiogyne"          , family = "Annonaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Melaleuca"         , family = "Myrtaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Melastoma"         , family = "Melastomataceae"  )
   n=n+1; g2f[[n]] = list( genus = "Melicoccus"        , family = "Sapindaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Metrodorea"        , family = "Rutaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Metrosideros"      , family = "Myrtaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Mezilaurus"        , family = "Lauraceae"        )
   n=n+1; g2f[[n]] = list( genus = "Miconia"           , family = "Melastomataceae"  )
   n=n+1; g2f[[n]] = list( genus = "Micrandra"         , family = "Euphorbiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Micrandropsis"     , family = "Euphorbiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Micromyrtus"       , family = "Myrtaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Micropholis"       , family = "Sapotaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Mikania"           , family = "Asteraceae"       )
   n=n+1; g2f[[n]] = list( genus = "Millettia"         , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Mimetes"           , family = "Proteaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Mimosa"            , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Mimulus"           , family = "Phrymaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Mimusops"          , family = "Sapotaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Minquartia"        , family = "Olacaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Misanteca"         , family = "Lauraceae"        )
   n=n+1; g2f[[n]] = list( genus = "Mitragyna"         , family = "Rubiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Molinia"           , family = "Poaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Monocarpia"        , family = "Annonaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Monopteryx"        , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Morisonia"         , family = "Capparaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Moronobea"         , family = "Clusiaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Mortoniodendron"   , family = "Malvaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Morus"             , family = "Moraceae"         )
   n=n+1; g2f[[n]] = list( genus = "Mouriri"           , family = "Melastomataceae"  )
   n=n+1; g2f[[n]] = list( genus = "Mucuna"            , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Myracrodruon"      , family = "Anacardiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Myrcia"            , family = "Myrtaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Myrciaria"         , family = "Myrtaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Myrica"            , family = "Myricaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Myriocarpa"        , family = "Urticaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Myrocarpus"        , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Myroxylon"         , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Myrsine"           , family = "Primulaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Myrtus"            , family = "Myrtaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Naucleopsis"       , family = "Moraceae"         )
   n=n+1; g2f[[n]] = list( genus = "Nectandra"         , family = "Lauraceae"        )
   n=n+1; g2f[[n]] = list( genus = "Neea"              , family = "Nyctaginaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Neolamarckia"      , family = "Rubiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Nephrolepis"       , family = "Davalliaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Nerium"            , family = "Apocynaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Newtonia"          , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Nothofagus"        , family = "Nothofagaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Nyctanthes"        , family = "Oleaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Ochroma"           , family = "Malvaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Ocotea"            , family = "Lauraceae"        )
   n=n+1; g2f[[n]] = list( genus = "Oenocarpus"        , family = "Arecaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Olea"              , family = "Oleaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Olearia"           , family = "Asteraceae"       )
   n=n+1; g2f[[n]] = list( genus = "Olinia"            , family = "Penaeaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Omphalea"          , family = "Euphorbiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Onychopetalum"     , family = "Annonaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Oreopanax"         , family = "Araliaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Ormosia"           , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Orthion"           , family = "Violaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Osteophloeum"      , family = "Myristicaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Ouratea"           , family = "Ochnaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Oxandra"           , family = "Annonaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Pachira"           , family = "Malvaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Paeonia"           , family = "Paeoniaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Palicourea"        , family = "Rubiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Parahancornia"     , family = "Apocynaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Paraia"            , family = "Lauraceae"        )
   n=n+1; g2f[[n]] = list( genus = "Paranomus"         , family = "Proteaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Parashorea"        , family = "Dipterocarpaceae" )
   n=n+1; g2f[[n]] = list( genus = "Parathesis"        , family = "Primulaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Parinari"          , family = "Chrysobalanaceae" )
   n=n+1; g2f[[n]] = list( genus = "Parkia"            , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Patersonia"        , family = "Iridaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Paullinia"         , family = "Sapindaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Pausandra"         , family = "Euphorbiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Paypayrola"        , family = "Violaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Peltogyne"         , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Pera"              , family = "Euphorbiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Perebea"           , family = "Moraceae"         )
   n=n+1; g2f[[n]] = list( genus = "Pereskia"          , family = "Cactaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Peridiscus"        , family = "Peridiscaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Persea"            , family = "Lauraceae"        )
   n=n+1; g2f[[n]] = list( genus = "Petalostigma"      , family = "Picrodendraceae"  )
   n=n+1; g2f[[n]] = list( genus = "Phenakospermum"    , family = "Strelitziaceae"   )
   n=n+1; g2f[[n]] = list( genus = "Phillyrea"         , family = "Oleaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Phlomis"           , family = "Lamiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Phyllota"          , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Picconia"          , family = "Oleaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Picea"             , family = "Pinaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Pickeringia"       , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Pimelea"           , family = "Thymelaeaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Pinus"             , family = "Pinaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Piper"             , family = "Piperaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Piptadenia"        , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Piptocarpha"       , family = "Asteraceae"       )
   n=n+1; g2f[[n]] = list( genus = "Pipturus"          , family = "Urticaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Pistacia"          , family = "Anacardiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Pithecellobium"    , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Pittoniotis"       , family = "Rubiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Planchonella"      , family = "Sapotaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Planchonia"        , family = "Lecythidaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Plathymenia"       , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Platonia"          , family = "Clusiaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Platymiscium"      , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Platypodium"       , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Pleiostachya"      , family = "Marantaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Plenckia"          , family = "Celastraceae"     )
   n=n+1; g2f[[n]] = list( genus = "Pleuranthodendron" , family = "Salicaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Podocarpus"        , family = "Podocarpaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Poecilanthe"       , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Poeppigia"         , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Pogonophora"       , family = "Euphorbiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Polyalthia"        , family = "Annonaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Polyscias"         , family = "Araliaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Pongamia"          , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Populus"           , family = "Salicaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Poraqueiba"        , family = "Icacinaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Posoqueria"        , family = "Rubiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Poulsenia"         , family = "Moraceae"         )
   n=n+1; g2f[[n]] = list( genus = "Pourouma"          , family = "Urticaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Pouteria"          , family = "Sapotaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Pradosia"          , family = "Sapotaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Prioria"           , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Prosopis"          , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Protea"            , family = "Proteaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Protium"           , family = "Burseraceae"      )
   n=n+1; g2f[[n]] = list( genus = "Protomegabaria"    , family = "Phyllanthaceae"   )
   n=n+1; g2f[[n]] = list( genus = "Prunus"            , family = "Rosaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Pseudima"          , family = "Sapindaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Pseudobombax"      , family = "Malvaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Pseudolmedia"      , family = "Moraceae"         )
   n=n+1; g2f[[n]] = list( genus = "Pseudopiptadenia"  , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Pseudoxandra"      , family = "Annonaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Psidium"           , family = "Myrtaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Psychotria"        , family = "Rubiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Ptaeroxylon"       , family = "Rutaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Pteleopsis"        , family = "Combretaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Pterandra"         , family = "Malpighiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Pteridium"         , family = "Dennstaedtiaceae" )
   n=n+1; g2f[[n]] = list( genus = "Pterocarpus"       , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Pterospermum"      , family = "Malvaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Pterygota"         , family = "Malvaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Ptychopetalum"     , family = "Olacaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Ptychopyxis"       , family = "Euphorbiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Purshia"           , family = "Rosaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Pyrus"             , family = "Rosaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Quadrella"         , family = "Capparaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Qualea"            , family = "Vochysiaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Quararibea"        , family = "Malvaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Quercus"           , family = "Fagaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Quiina"            , family = "Ochnaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Rapanea"           , family = "Primulaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Raputia"           , family = "Rutaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Rauvolfia"         , family = "Apocynaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Recordoxylon"      , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Regelia"           , family = "Myrtaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Remijia"           , family = "Rubiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Retiniphyllum"     , family = "Rubiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Rhabdodendron"     , family = "Rhabdodendraceae" )
   n=n+1; g2f[[n]] = list( genus = "Rhamnus"           , family = "Rhamnaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Rhizophora"        , family = "Rhizophoraceae"   )
   n=n+1; g2f[[n]] = list( genus = "Rhodamnia"         , family = "Myrtaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Rhododendron"      , family = "Ericaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Rhodostemonodaphne", family = "Lauraceae"        )
   n=n+1; g2f[[n]] = list( genus = "Ribes"             , family = "Grossulariaceae"  )
   n=n+1; g2f[[n]] = list( genus = "Richeria"          , family = "Phyllanthaceae"   )
   n=n+1; g2f[[n]] = list( genus = "Rinorea"           , family = "Violaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Robinia"           , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Robinsonella"      , family = "Malvaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Rollinia"          , family = "Annonaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Rosa"              , family = "Rosaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Roseodendron"      , family = "Bignoniaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Rosmarinus"        , family = "Lamiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Roucheria"         , family = "Linaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Roupala"           , family = "Proteaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Rubus"             , family = "Rosaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Ruizterania"       , family = "Vochysiaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Rumex"             , family = "Polygonaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Ryania"            , family = "Salicaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Sacoglottis"       , family = "Humiriaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Sagotia"           , family = "Euphorbiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Salacia"           , family = "Celastraceae"     )
   n=n+1; g2f[[n]] = list( genus = "Salix"             , family = "Salicaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Sambucus"          , family = "Adoxaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Santalum"          , family = "Santalaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Sapindus"          , family = "Sapindaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Sapium"            , family = "Euphorbiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Saraca"            , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Sarcaulus"         , family = "Sapotaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Scaevola"          , family = "Goodeniaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Schefflera"        , family = "Araliaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Schinopsis"        , family = "Anacardiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Schizolobium"      , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Schleichera"       , family = "Sapindaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Scholtzia"         , family = "Myrtaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Sclerolobium"      , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Scleronema"        , family = "Malvaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Scolopia"          , family = "Salicaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Scyphiphora"       , family = "Rubiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Senefeldera"       , family = "Euphorbiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Senegalia"         , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Senna"             , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Serruria"          , family = "Proteaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Sextonia"          , family = "Lauraceae"        )
   n=n+1; g2f[[n]] = list( genus = "Shorea"            , family = "Dipterocarpaceae" )
   n=n+1; g2f[[n]] = list( genus = "Sideroxylon"       , family = "Sapotaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Silene"            , family = "Caryophyllaceae"  )
   n=n+1; g2f[[n]] = list( genus = "Simaba"            , family = "Simaroubaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Simarouba"         , family = "Simaroubaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Siparuna"          , family = "Siparunaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Sloanea"           , family = "Elaeocarpaceae"   )
   n=n+1; g2f[[n]] = list( genus = "Smilax"            , family = "Smilacaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Socratea"          , family = "Arecaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Solanum"           , family = "Solanaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Sonneratia"        , family = "Lythraceae"       )
   n=n+1; g2f[[n]] = list( genus = "Sorbus"            , family = "Rosaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Sorocea"           , family = "Moraceae"         )
   n=n+1; g2f[[n]] = list( genus = "Soymida"           , family = "Meliaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Spartothamnella"   , family = "Lamiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Spathodea"         , family = "Bignoniaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Spondias"          , family = "Anacardiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Stemmadenia"       , family = "Apocynaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Sterculia"         , family = "Malvaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Sterigmapetalum"   , family = "Rhizophoraceae"   )
   n=n+1; g2f[[n]] = list( genus = "Strombosia"        , family = "Olacaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Strychnos"         , family = "Loganiaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Stryphnodendron"   , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Styrax"            , family = "Styracaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Succisa"           , family = "Caprifoliaceae"   )
   n=n+1; g2f[[n]] = list( genus = "Suregada"          , family = "Euphorbiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Swartzia"          , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Swietenia"         , family = "Meliaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Syagrus"           , family = "Arecaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Symphonia"         , family = "Clusiaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Syrmatium"         , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Syzygium"          , family = "Myrtaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Tabebuia"          , family = "Bignoniaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Tabernaemontana"   , family = "Apocynaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Tachigali"         , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Talisia"           , family = "Sapindaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Tamarindus"        , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Tamilnadia"        , family = "Rubiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Tapirira"          , family = "Anacardiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Tapura"            , family = "Dichapetalaceae"  )
   n=n+1; g2f[[n]] = list( genus = "Taralea"           , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Tarenna"           , family = "Rubiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Tasmannia"         , family = "Winteraceae"      )
   n=n+1; g2f[[n]] = list( genus = "Taxus"             , family = "Taxaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Teclea"            , family = "Rutaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Tectona"           , family = "Lamiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Terminalia"        , family = "Combretaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Tetragastris"      , family = "Burseraceae"      )
   n=n+1; g2f[[n]] = list( genus = "Theobroma"         , family = "Malvaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Thespesia"         , family = "Malvaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Thymus"            , family = "Lamiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Thyrsodium"        , family = "Anacardiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Tibouchina"        , family = "Melastomataceae"  )
   n=n+1; g2f[[n]] = list( genus = "Tocoyena"          , family = "Rubiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Toulicia"          , family = "Sapindaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Touroulia"         , family = "Ochnaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Tovomita"          , family = "Clusiaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Toxicodendron"     , family = "Anacardiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Trattinnickia"     , family = "Burseraceae"      )
   n=n+1; g2f[[n]] = list( genus = "Trema"             , family = "Cannabaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Triadica"          , family = "Euphorbiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Trichilia"         , family = "Meliaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Trichocladus"      , family = "Hamamelidaceae"   )
   n=n+1; g2f[[n]] = list( genus = "Trichoscypha"      , family = "Anacardiaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Trichospermum"     , family = "Malvaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Trigonobalanus"    , family = "Fagaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Trilepisium"       , family = "Moraceae"         )
   n=n+1; g2f[[n]] = list( genus = "Tristaniopsis"     , family = "Myrtaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Trophis"           , family = "Moraceae"         )
   n=n+1; g2f[[n]] = list( genus = "Trymatococcus"     , family = "Moraceae"         )
   n=n+1; g2f[[n]] = list( genus = "Turpinia"          , family = "Staphyleaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Uapaca"            , family = "Phyllanthaceae"   )
   n=n+1; g2f[[n]] = list( genus = "Ulex"              , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Ulmus"             , family = "Ulmaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Umbellularia"      , family = "Lauraceae"        )
   n=n+1; g2f[[n]] = list( genus = "Umtiza"            , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Unonopsis"         , family = "Annonaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Urera"             , family = "Urticaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Urtica"            , family = "Urticaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Uvaria"            , family = "Annonaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Vaccinium"         , family = "Ericaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Vantanea"          , family = "Humiriaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Vatairea"          , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Vataireopsis"      , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Vepris"            , family = "Rutaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Verbascum"         , family = "Scrophulariaceae" )
   n=n+1; g2f[[n]] = list( genus = "Vernonia"          , family = "Asteraceae"       )
   n=n+1; g2f[[n]] = list( genus = "Verticordia"       , family = "Myrtaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Viburnum"          , family = "Adoxaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Viguieranthus"     , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Viola"             , family = "Violaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Virola"            , family = "Myristicaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Vismia"            , family = "Hypericaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Vitex"             , family = "Lamiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Vitis"             , family = "Vitaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Vochysia"          , family = "Vochysiaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Votomita"          , family = "Melastomataceae"  )
   n=n+1; g2f[[n]] = list( genus = "Vouacapoua"        , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Warszewiczia"      , family = "Rubiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Woodfordia"        , family = "Lythraceae"       )
   n=n+1; g2f[[n]] = list( genus = "Wrightia"          , family = "Apocynaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Xanthophyllum"     , family = "Polygalaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Xanthorrhoea"      , family = "Xanthorrhoeaceae" )
   n=n+1; g2f[[n]] = list( genus = "Xanthostemon"      , family = "Myrtaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Ximenia"           , family = "Olacaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Xylopia"           , family = "Annonaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Zanthoxylum"       , family = "Rutaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Ziziphus"          , family = "Rhamnaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Zollernia"         , family = "Fabaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Zuelania"          , family = "Salicaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Zygia"             , family = "Fabaceae"         )
   #---------------------------------------------------------------------------------------#



   #----- Convert the g2f list into a data frame. -----------------------------------------#
   g2f = list.2.data.frame(g2f)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Check whether the family has all genera needed.                                   #
   #---------------------------------------------------------------------------------------#
   idx = which( ( ! is.na(datum$genus) ) 
              & ( ! datum$genus %in% g2f$genus )
              & ( regexpr(pattern="Ignotum",text=datum$genus,ignore.case=TRUE) %==% -1)
              )#end which
   if (length(idx) > 0){
     tofill.gen = t(t(sort(unique(datum$genus[idx]))))
     toshow     = which(names(datum) %in% c("trans","tag","x","y","genus","family"))
     toshow     = datum[idx,toshow]
     toshow     = toshow[order(toshow[,"genus"]),]
     toshow     = toshow[! duplicated(toshow[,"genus"]),]
     cat0("-----------------------------------------------------------")
     cat0(" You must add genera to g2f:")
     cat0(" ")
     print(toshow,quote=FALSE)
     cat0("-----------------------------------------------------------")
     cat0(" ")
     browser()
   }#end if
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     List the default name for all families if genera was unknown.                     #
   #---------------------------------------------------------------------------------------#
   n=0  ; f2nog      = list()
   n=n+1; f2nog[[n]] = list(family = "Acanthaceae"      , genus = "Ignotum.acanthus"      )
   n=n+1; f2nog[[n]] = list(family = "Achariaceae"      , genus = "Ignotum.acharia"       )
   n=n+1; f2nog[[n]] = list(family = "Adoxaceae"        , genus = "Ignotum.adoxa"         )
   n=n+1; f2nog[[n]] = list(family = "Altingiaceae"     , genus = "Ignotum.altingia"      )
   n=n+1; f2nog[[n]] = list(family = "Amaranthaceae"    , genus = "Ignotum.amaranthus"    )
   n=n+1; f2nog[[n]] = list(family = "Anacardiaceae"    , genus = "Ignotum.anacardium"    )
   n=n+1; f2nog[[n]] = list(family = "Anarthriaceae"    , genus = "Ignotum.anarthria"     )
   n=n+1; f2nog[[n]] = list(family = "Anisophylleaceae" , genus = "Ignotum.anisophyllea"  )
   n=n+1; f2nog[[n]] = list(family = "Annonaceae"       , genus = "Ignotum.annona"        )
   n=n+1; f2nog[[n]] = list(family = "Apiaceae"         , genus = "Ignotum.apium"         )
   n=n+1; f2nog[[n]] = list(family = "Apocynaceae"      , genus = "Ignotum.apocynum"      )
   n=n+1; f2nog[[n]] = list(family = "Aquifoliaceae"    , genus = "Ignotum.ilex"          )
   n=n+1; f2nog[[n]] = list(family = "Araliaceae"       , genus = "Ignotum.aralia"        )
   n=n+1; f2nog[[n]] = list(family = "Arecaceae"        , genus = "Ignotum.areca"         )
   n=n+1; f2nog[[n]] = list(family = "Asparagaceae"     , genus = "Ignotum.asparagus"     )
   n=n+1; f2nog[[n]] = list(family = "Asteraceae"       , genus = "Ignotum.aster"         )
   n=n+1; f2nog[[n]] = list(family = "Berberidaceae"    , genus = "Ignotum.berberis"      )
   n=n+1; f2nog[[n]] = list(family = "Betulaceae"       , genus = "Ignotum.betula"        )
   n=n+1; f2nog[[n]] = list(family = "Bignoniaceae"     , genus = "Ignotum.bignonia"      )
   n=n+1; f2nog[[n]] = list(family = "Bixaceae"         , genus = "Ignotum.bixa"          )
   n=n+1; f2nog[[n]] = list(family = "Boraginaceae"     , genus = "Ignotum.borago"        )
   n=n+1; f2nog[[n]] = list(family = "Burseraceae"      , genus = "Ignotum.bursera"       )
   n=n+1; f2nog[[n]] = list(family = "Buxaceae"         , genus = "Ignotum.buxus"         )
   n=n+1; f2nog[[n]] = list(family = "Brachytheciaceae" , genus = "Ignotum.brachythecium" )
   n=n+1; f2nog[[n]] = list(family = "Cactaceae"        , genus = "Ignotum.cactus"        )
   n=n+1; f2nog[[n]] = list(family = "Calophyllaceae"   , genus = "Ignotum.calophyllum"   )
   n=n+1; f2nog[[n]] = list(family = "Calycanthaceae"   , genus = "Ignotum.calycanthus"   )
   n=n+1; f2nog[[n]] = list(family = "Cannabaceae"      , genus = "Ignotum.cannabis"      )
   n=n+1; f2nog[[n]] = list(family = "Capparaceae"      , genus = "Ignotum.capparis"      )
   n=n+1; f2nog[[n]] = list(family = "Caprifoliaceae"   , genus = "Ignotum.caprifolium"   )
   n=n+1; f2nog[[n]] = list(family = "Cardiopteridaceae", genus = "Ignotum.cardiopteris"  )
   n=n+1; f2nog[[n]] = list(family = "Caricaceae"       , genus = "Ignotum.carica"        )
   n=n+1; f2nog[[n]] = list(family = "Caryocaraceae"    , genus = "Ignotum.caryocar"      )
   n=n+1; f2nog[[n]] = list(family = "Caryophyllaceae"  , genus = "Ignotum.charyophyllum" )
   n=n+1; f2nog[[n]] = list(family = "Casuarinaceae"    , genus = "Ignotum.casuarina"     )
   n=n+1; f2nog[[n]] = list(family = "Celastraceae"     , genus = "Ignotum.celastrus"     )
   n=n+1; f2nog[[n]] = list(family = "Chloranthaceae"   , genus = "Ignotum.chloranthus"   )
   n=n+1; f2nog[[n]] = list(family = "Chrysobalanaceae" , genus = "Ignotum.chrysobalanus" )
   n=n+1; f2nog[[n]] = list(family = "Cistaceae"        , genus = "Ignotum.cistus"        )
   n=n+1; f2nog[[n]] = list(family = "Clethraceae"      , genus = "Ignotum.clethra"       )
   n=n+1; f2nog[[n]] = list(family = "Clusiaceae"       , genus = "Ignotum.clusia"        )
   n=n+1; f2nog[[n]] = list(family = "Combretaceae"     , genus = "Ignotum.combretum"     )
   n=n+1; f2nog[[n]] = list(family = "Connaraceae"      , genus = "Ignotum.connarus"      )
   n=n+1; f2nog[[n]] = list(family = "Convolvulaceae"   , genus = "Ignotum.convolvulus"   )
   n=n+1; f2nog[[n]] = list(family = "Cornaceae"        , genus = "Ignotum.cornus"        )
   n=n+1; f2nog[[n]] = list(family = "Cupressaceae"     , genus = "Ignotum.cupressus"     )
   n=n+1; f2nog[[n]] = list(family = "Curtisiaceae"     , genus = "Curtisia"              )
   n=n+1; f2nog[[n]] = list(family = "Cycadaceae"       , genus = "Cycas"                 )
   n=n+1; f2nog[[n]] = list(family = "Cyperaceae"       , genus = "Ignotum.cyperus"       )
   n=n+1; f2nog[[n]] = list(family = "Cyrillaceae"      , genus = "Ignotum.cyrilla"       )
   n=n+1; f2nog[[n]] = list(family = "Davalliaceae"     , genus = "Ignotum.davallia"      )
   n=n+1; f2nog[[n]] = list(family = "Dennstaedtiaceae" , genus = "Ignotum.dennstaedtia"  )
   n=n+1; f2nog[[n]] = list(family = "Dichapetalaceae"  , genus = "Ignotum.dichapetalum"  )
   n=n+1; f2nog[[n]] = list(family = "Dilleniaceae"     , genus = "Ignotum.dillenia"      )
   n=n+1; f2nog[[n]] = list(family = "Dipterocarpaceae" , genus = "Ignotum.dipterocarpus" )
   n=n+1; f2nog[[n]] = list(family = "Ebenaceae"        , genus = "Ignotum.ebenus"        )
   n=n+1; f2nog[[n]] = list(family = "Elaeagnaceae"     , genus = "Ignotum.elaeagnus"     )
   n=n+1; f2nog[[n]] = list(family = "Elaeocarpaceae"   , genus = "Ignotum.elaeocarpus"   )
   n=n+1; f2nog[[n]] = list(family = "Emmotaceae"       , genus = "Ignotum.emmotum"       )
   n=n+1; f2nog[[n]] = list(family = "Ericaceae"        , genus = "Ignotum.erica"         )
   n=n+1; f2nog[[n]] = list(family = "Erythroxylaceae"  , genus = "Ignotum.erythroxylum"  )
   n=n+1; f2nog[[n]] = list(family = "Euphorbiaceae"    , genus = "Ignotum.euphorbia"     )
   n=n+1; f2nog[[n]] = list(family = "Fabaceae"         , genus = "Ignotum.faba"          )
   n=n+1; f2nog[[n]] = list(family = "Fagaceae"         , genus = "Ignotum.fagus"         )
   n=n+1; f2nog[[n]] = list(family = "Gentianaceae"     , genus = "Ignotum.gentiana"      )
   n=n+1; f2nog[[n]] = list(family = "Gleicheniaceae"   , genus = "Ignotum.gleichenia"    )
   n=n+1; f2nog[[n]] = list(family = "Goodeniaceae"     , genus = "Ignotum.goodenia"      )
   n=n+1; f2nog[[n]] = list(family = "Goupiaceae"       , genus = "Goupia"                )
   n=n+1; f2nog[[n]] = list(family = "Grossulariaceae"  , genus = "Ignotum.grossularia"   )
   n=n+1; f2nog[[n]] = list(family = "Haemodoraceae"    , genus = "Ignotum.haemodorum"    )
   n=n+1; f2nog[[n]] = list(family = "Hamamelidaceae"   , genus = "Ignotum.hamamelis"     )
   n=n+1; f2nog[[n]] = list(family = "Humiriaceae"      , genus = "Ignotum.humiria"       )
   n=n+1; f2nog[[n]] = list(family = "Hypericaceae"     , genus = "Ignotum.hypericum"     )
   n=n+1; f2nog[[n]] = list(family = "Icacinaceae"      , genus = "Ignotum.icacina"       )
   n=n+1; f2nog[[n]] = list(family = "Ignotaceae"       , genus = "Ignotum"               )
   n=n+1; f2nog[[n]] = list(family = "Iridaceae"        , genus = "Ignotum.iris"          )
   n=n+1; f2nog[[n]] = list(family = "Ixonanthaceae"    , genus = "Ignotum.ixonanthes"    )
   n=n+1; f2nog[[n]] = list(family = "Juglandaceae"     , genus = "Ignotum.juglans"       )
   n=n+1; f2nog[[n]] = list(family = "Lacistemataceae"  , genus = "Ignotum.lacistema"     )
   n=n+1; f2nog[[n]] = list(family = "Lamiaceae"        , genus = "Ignotum.lamium"        )
   n=n+1; f2nog[[n]] = list(family = "Lauraceae"        , genus = "Ignotum.laurus"        )
   n=n+1; f2nog[[n]] = list(family = "Lecythidaceae"    , genus = "Ignotum.lecythis"      )
   n=n+1; f2nog[[n]] = list(family = "Lianaceae"        , genus = "Liana"                 )
   n=n+1; f2nog[[n]] = list(family = "Linaceae"         , genus = "Ignotum.linum"         )
   n=n+1; f2nog[[n]] = list(family = "Loganiaceae"      , genus = "Ignotum.logania"       )
   n=n+1; f2nog[[n]] = list(family = "Lythraceae"       , genus = "Ignotum.lythrum"       )
   n=n+1; f2nog[[n]] = list(family = "Magnoliaceae"     , genus = "Ignotum.magnolia"      )
   n=n+1; f2nog[[n]] = list(family = "Malpighiaceae"    , genus = "Ignotum.malpighia"     )
   n=n+1; f2nog[[n]] = list(family = "Malvaceae"        , genus = "Ignotum.malva"         )
   n=n+1; f2nog[[n]] = list(family = "Marantaceae"      , genus = "Ignotum.maranta"       )
   n=n+1; f2nog[[n]] = list(family = "Marcgraviaceae"   , genus = "Ignotum.marcgravia"    )
   n=n+1; f2nog[[n]] = list(family = "Melastomataceae"  , genus = "Ignotum.melastoma"     )
   n=n+1; f2nog[[n]] = list(family = "Menispermaceae"   , genus = "Ignotum.menispermum"   )
   n=n+1; f2nog[[n]] = list(family = "Meliaceae"        , genus = "Ignotum.melia"         )
   n=n+1; f2nog[[n]] = list(family = "Moraceae"         , genus = "Ignotum.morus"         )
   n=n+1; f2nog[[n]] = list(family = "Myricaceae"       , genus = "Ignotum.myrica"        )
   n=n+1; f2nog[[n]] = list(family = "Myristicaceae"    , genus = "Ignotum.myristica"     )
   n=n+1; f2nog[[n]] = list(family = "Myrtaceae"        , genus = "Ignotum.myrtus"        )
   n=n+1; f2nog[[n]] = list(family = "Nothofagaceae"    , genus = "Nothofagus"            )
   n=n+1; f2nog[[n]] = list(family = "Nyctaginaceae"    , genus = "Ignotum.nyctaginia"    )
   n=n+1; f2nog[[n]] = list(family = "Ochnaceae"        , genus = "Ignotum.ochna"         )
   n=n+1; f2nog[[n]] = list(family = "Olacaceae"        , genus = "Ignotum.olax"          )
   n=n+1; f2nog[[n]] = list(family = "Oleaceae"         , genus = "Ignotum.olea"          )
   n=n+1; f2nog[[n]] = list(family = "Opiliaceae"       , genus = "Ignotum.opilia"        )
   n=n+1; f2nog[[n]] = list(family = "Paeoniaceae"      , genus = "Paeonia"               )
   n=n+1; f2nog[[n]] = list(family = "Passifloraceae"   , genus = "Ignotum.passiflora"    )
   n=n+1; f2nog[[n]] = list(family = "Penaeaceae"       , genus = "Ignotum.penaea"        )
   n=n+1; f2nog[[n]] = list(family = "Pentaphylacaceae" , genus = "Ignotum.pentaphylax"   )
   n=n+1; f2nog[[n]] = list(family = "Peridiscaceae"    , genus = "Ignotum.peridiscus"    )
   n=n+1; f2nog[[n]] = list(family = "Phrymaceae"       , genus = "Ignotum.phryma"        )
   n=n+1; f2nog[[n]] = list(family = "Phyllanthaceae"   , genus = "Ignotum.phyllanthus"   )
   n=n+1; f2nog[[n]] = list(family = "Phytolaccaceae"   , genus = "Ignotum.phytolacca"    )
   n=n+1; f2nog[[n]] = list(family = "Picrodendraceae"  , genus = "Ignotum.picrodendron"  )
   n=n+1; f2nog[[n]] = list(family = "Pinaceae"         , genus = "Ignotum.pinus"         )
   n=n+1; f2nog[[n]] = list(family = "Piperaceae"       , genus = "Ignotum.piper"         )
   n=n+1; f2nog[[n]] = list(family = "Plantaginaceae"   , genus = "Ignotum.plantago"      )
   n=n+1; f2nog[[n]] = list(family = "Plumbaginaceae"   , genus = "Ignotum.plumbago"      )
   n=n+1; f2nog[[n]] = list(family = "Poaceae"          , genus = "Ignotum.poa"           )
   n=n+1; f2nog[[n]] = list(family = "Podocarpaceae"    , genus = "Ignotum.podocarpus"    )
   n=n+1; f2nog[[n]] = list(family = "Polygalaceae"     , genus = "Ignotum.polygala"      )
   n=n+1; f2nog[[n]] = list(family = "Polygonaceae"     , genus = "Ignotum.polygonum"     )
   n=n+1; f2nog[[n]] = list(family = "Primulaceae"      , genus = "Ignotum.primula"       )
   n=n+1; f2nog[[n]] = list(family = "Proteaceae"       , genus = "Ignotum.protea"        )
   n=n+1; f2nog[[n]] = list(family = "Putranjivaceae"   , genus = "Ignotum.putranjiva"    )
   n=n+1; f2nog[[n]] = list(family = "Ranunculaceae"    , genus = "Ignotum.ranunculus"    )
   n=n+1; f2nog[[n]] = list(family = "Restionaceae"     , genus = "Ignotum.restio"        )
   n=n+1; f2nog[[n]] = list(family = "Rhabdodendraceae" , genus = "Ignotum.rhabdodendron" )
   n=n+1; f2nog[[n]] = list(family = "Rhabdodendraceae" , genus = "Ignotum.rhabdodendron" )
   n=n+1; f2nog[[n]] = list(family = "Rhamnaceae"       , genus = "Ignotum.rhamnus"       )
   n=n+1; f2nog[[n]] = list(family = "Rhizophoraceae"   , genus = "Ignotum.rhizophora"    )
   n=n+1; f2nog[[n]] = list(family = "Rosaceae"         , genus = "Ignotum.rosa"          )
   n=n+1; f2nog[[n]] = list(family = "Rubiaceae"        , genus = "Ignotum.rubia"         )
   n=n+1; f2nog[[n]] = list(family = "Rutaceae"         , genus = "Ignotum.ruta"          )
   n=n+1; f2nog[[n]] = list(family = "Salicaceae"       , genus = "Ignotum.salix"         )
   n=n+1; f2nog[[n]] = list(family = "Santalaceae"      , genus = "Ignotum.santalum"      )
   n=n+1; f2nog[[n]] = list(family = "Sapindaceae"      , genus = "Ignotum.sapindus"      )
   n=n+1; f2nog[[n]] = list(family = "Sapotaceae"       , genus = "Ignotum.sapota"        )
   n=n+1; f2nog[[n]] = list(family = "Scrophulariaceae" , genus = "Ignotum.scrophularia"  )
   n=n+1; f2nog[[n]] = list(family = "Simaroubaceae"    , genus = "Ignotum.simarouba"     )
   n=n+1; f2nog[[n]] = list(family = "Siparunaceae"     , genus = "Ignotum.siparuna"      )
   n=n+1; f2nog[[n]] = list(family = "Smilacaceae"      , genus = "Ignotum.smilax"        )
   n=n+1; f2nog[[n]] = list(family = "Solanaceae"       , genus = "Ignotum.solanum"       )
   n=n+1; f2nog[[n]] = list(family = "Staphyleaceae"    , genus = "Ignotum.staphylea"     )
   n=n+1; f2nog[[n]] = list(family = "Stemonuraceae"    , genus = "Ignotum.stemonurus"    )
   n=n+1; f2nog[[n]] = list(family = "Strelitziaceae"   , genus = "Ignotum.strelitzia"    )
   n=n+1; f2nog[[n]] = list(family = "Styracaceae"      , genus = "Ignotum.styrax"        )
   n=n+1; f2nog[[n]] = list(family = "Taxaceae"         , genus = "Ignotum.taxus"         )
   n=n+1; f2nog[[n]] = list(family = "Thymelaeaceae"    , genus = "Ignotum.thymelaea"     )
   n=n+1; f2nog[[n]] = list(family = "Ulmaceae"         , genus = "Ignotum.ulmus"         )
   n=n+1; f2nog[[n]] = list(family = "Urticaceae"       , genus = "Ignotum.urtica"        )
   n=n+1; f2nog[[n]] = list(family = "Verbenaceae"      , genus = "Ignotum.verbena"       )
   n=n+1; f2nog[[n]] = list(family = "Violaceae"        , genus = "Ignotum.viola"         )
   n=n+1; f2nog[[n]] = list(family = "Vitaceae"         , genus = "Ignotum.vitis"         )
   n=n+1; f2nog[[n]] = list(family = "Vochysiaceae"     , genus = "Ignotum.vochysia"      )
   n=n+1; f2nog[[n]] = list(family = "Winteraceae"      , genus = "Ignotum.wintera"       )
   n=n+1; f2nog[[n]] = list(family = "Xanthorrhoeaceae" , genus = "Ignotum.xanthorrhoea"  )
   n=n+1; f2nog[[n]] = list(family = "Zamiaceae"        , genus = "Ignotum.zamia"         )
   n=n+1; f2nog[[n]] = list(family = "Zygophyllaceae"   , genus = "Ignotum.zygophyllum"   )
   #---------------------------------------------------------------------------------------#





   #----- Convert the g2f list into a data frame. -----------------------------------------#
   f2nog = list.2.data.frame(f2nog)
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Check whether there are families missing in f2nog look up table.                  #
   #---------------------------------------------------------------------------------------#
   idx = which(! g2f$family %in% f2nog$family)
   if (length(idx) > 0){
     cat(" You must add families to f2nog...","\n")
     family.tofill = t(t(sort(unique(g2f$family[idx]))))
     browser()
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Fill in non-informative genus for individual with known family but unknown genus. #
   #---------------------------------------------------------------------------------------#
   yes.gen = ! is.na(datum$genus ) & regexpr(pattern="Ignotum"   ,text=datum$genus ) %<% 0
   yes.fam = ! is.na(datum$family) & regexpr(pattern="Ignotaceae",text=datum$family) %<% 0
   nothing = ! yes.gen & ! yes.fam
   no.gen  = ! yes.gen &   yes.fam
   #------ No information whatsoever.  Fill in with non-informative taxa. -----------------#
   datum$genus      [nothing] = "Ignotum"
   datum$scientific [nothing] = "Ignotum NA"
   datum$family     [nothing] = "Ignotaceae"
   #------ No genus, family only. ---------------------------------------------------------#
   idx                        = match(datum$family[no.gen],f2nog$family)
   datum$genus      [no.gen ] = f2nog$genus[idx]
   datum$scientific [no.gen ] = paste(datum$genus[no.gen],NA_character_,sep=" ")
   #------ Genus is known. ----------------------------------------------------------------#
   idx                        = match(datum$genus[yes.gen],g2f$genus)
   datum$family     [yes.gen] = g2f$family [idx]
   #---------------------------------------------------------------------------------------#

   return(datum)
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      This attributes scientific names based on common names for surveys carried out by   #
# the Sustainable Landscapes team.  It is not a good idea to use this function for any     #
# other data sets because common names may mean completely different things for other      #
# lads.                                                                                    #
#------------------------------------------------------------------------------------------#
scientific.lookup.SL <<- function(datum,lookup.path){
   #----- Append columns in case they don't exist. ----------------------------------------#
   if (! "scientific"    %in% names(datum)){
      datum$scientific    = rep(NA_character_, nrow(datum))
   }#end if (! "scientific"    %in% names(datum))
   if (! "genus"         %in% names(datum)){
      datum$genus         = rep(NA_character_, nrow(datum)) 
   }#end if (! "scientific"    %in% names(datum))
   if (! "gf.scientific" %in% names(datum)){
      datum$gf.scientific = rep(0, nrow(datum))
   }#end if (! "gf.scientific" %in% names(datum))
   #---------------------------------------------------------------------------------------#


   #----- Read in the look-up table. ------------------------------------------------------#
   lookup.file = paste(lookup.path,"SL_taxon_lookup.csv",sep="/")
   look.up = read.csv(file=lookup.file,stringsAsFactors=FALSE)
   look.up$common     = tolower(trim(look.up$common    ))
   look.up$scientific = trim(look.up$scientific)
   look.up$family     = trim(look.up$family    )
   #---------------------------------------------------------------------------------------#
        
        
        
   #----- Break into genus and species. ---------------------------------------------------#
   gs.list                  = sapply(X = tolower(look.up$scientific),FUN=strsplit,split=" ")
   gs.length                = sapply(X = gs.list, FUN = length)
   gs.mat                   = cbind( mapply(FUN="[",gs.list,MoreArgs=list(1))
                                     , mapply(FUN="[",gs.list,MoreArgs=list(2))
   )#end cbind
   g                        = capwords(gs.mat[,1],strict=TRUE)
   s                        = tolower(gs.mat[,2])
   g.s                      = paste(g,s,sep=" ")
   g.s[is.na(g) & is.na(s)] = NA
   look.up$scientific       = g.s
   look.up$genus            = g
   #---------------------------------------------------------------------------------------#
   
   
   
   #----- Trim the common names, and simplify/replace some names. -------------------------#
   datum$common                      = tolower(trim(datum$common))
   datum$common[is.na(datum$common)] = "mato"
   #---------------------------------------------------------------------------------------#
        
        
       
   #----- Find all unique common names. ---------------------------------------------------#
   unique.common = unique(datum$common)
   n.common      = length(unique.common)
   notfound      = NULL
   for (n in sequence(n.common)){
      #----- Find the trees that have the same common name as this one. -------------------#
      cat (" - ",n,"/",n.common," -- ",unique.common[n],"\n")
      w.dat  = which(datum$common %in% unique.common[n])
      n.dat  = length(w.dat)
      #------------------------------------------------------------------------------------#
      
      
      #----- Find the trees in the look-up table with the same common name. ---------------#
      w.look = which(look.up$common %in% unique.common[n])
      n.look = length(w.look)
      #------------------------------------------------------------------------------------#
      
      
      
      #------------------------------------------------------------------------------------#
      #      Check how many trees have the same common name in the look-up table.          #
      #------------------------------------------------------------------------------------#
      if (n.look >= 1){
         #---------------------------------------------------------------------------------#
         #   In case only one scientific name has been matched, this will allocate the     #
         # same scientific name for all selected individuals, otherwise this will randomly #
         # attribute the scientific names.                                                 #
         #---------------------------------------------------------------------------------#
         w.look                     = lit.sample(x=w.look,size=n.dat,replace=TRUE)
         #---------------------------------------------------------------------------------#


         #----- Only one.  Use it. --------------------------------------------------------#
         datum$scientific   [w.dat] = look.up$scientific[w.look]
         datum$genus        [w.dat] = look.up$genus     [w.look]
         datum$gf.scientific[w.dat] = 1
         #---------------------------------------------------------------------------------#
      }else{
         #----- Not found in the data base, warn the user. --------------------------------#
         notfound                   = c(notfound,unique.common[n])
         datum$scientific   [w.dat] = "Ignotum"
         datum$genus        [w.dat] = "Ignotum"
         datum$gf.scientific[w.dat] = 0
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#
   return(datum)
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      This function is a look-up table for scientific names given the common name.  The   #
# look up table came from the reference below and from the first survey in 1999, and it is #
# probably fine for Km 67 but not other places as they may be regional names and they may  #
# refer to a completely different plant elsewhere.  Only the trees that occurred in the    #
# survey were added, feel free to complete the list...                                     #
#                                                                                          #
# Silva, J. N. M.; Carvalho, J. O. P.; Lopes, J. C. A., 1985: Forest inventory of an       #
#     experimental area in the Tapajos National Forest.  Boletim de Pesquisa Florestal,    #
#     n. 10/11, p.38--110.                                                                 #
#------------------------------------------------------------------------------------------#
scientific.lookup.tnf <<- function(datum){
   n.datum = nrow(datum)

   #----- Look-up table, feel free to add more trees here. --------------------------------#
   tnf        = list()
   tnf[[  1]] = list( common     = "abiu"                      
                    , scientific = "Pouteria reticulata"             
                    )#end list
   tnf[[  2]] = list( common     = "abiu branco"               
                    , scientific = "Pradosia ptychandra"             
                    )#end list
   tnf[[  3]] = list( common     = "abiu cutite"               
                    , scientific = "Pouteria macrophylla"            
                    )#end list
   tnf[[  4]] = list( common     = "abiu cutite folha verde"   
                    , scientific = "Pouteria venosa"                 
                    )#end list
   tnf[[  5]] = list( common     = "abiu vermelho"             
                    , scientific = "Pouteria torta"                  
                    )#end list
   tnf[[  6]] = list( common     = "abiu-casca-grossa"         
                    , scientific = "Pouteria bilocularis"            
                    )#end list
   tnf[[  7]] = list( common     = "abiu-mangabinha"           
                    , scientific = "Micropholis egensis"             
                    )#end list
   tnf[[  8]] = list( common     = "abiurana"                  
                    , scientific = "Pouteria torta"                  
                    )#end list
   tnf[[  9]] = list( common     = "abiurana vermelha"         
                    , scientific = "Pouteria torta"                  
                    )#end list
   tnf[[ 10]] = list( common     = "acariquara"                
                    , scientific = "Minquartia guianensis"           
                    )#end list
   tnf[[ 11]] = list( common     = "acoita-cavalo"             
                    , scientific = "Lueheopsis duckeana"             
                    )#end list
   tnf[[ 12]] = list( common     = "amapa"                     
                    , scientific = "Brosimum"                        
                    )#end list
   tnf[[ 13]] = list( common     = "amapa amargoso"            
                    , scientific = "Brosimum guianense"              
                    )#end list
   tnf[[ 14]] = list( common     = "amapa doce"                
                    , scientific = "Brosimum parinarioides"          
                    )#end list
   tnf[[ 15]] = list( common     = "amapai"                    
                    , scientific = "Brosimum lactescens"             
                    )#end list
   tnf[[ 16]] = list( common     = "amaparana"                 
                    , scientific = "Thyrsodium spruceanum"           
                    )#end list
   tnf[[ 17]] = list( common     = "amarelao"                  
                    , scientific = "Apuleia leiocarpa"               
                    )#end list
   tnf[[ 18]] = list( common     = "anani"                     
                    , scientific = "Symphonia globulifera"           
                    )#end list
   tnf[[ 19]] = list( common     = "andiroba"                  
                    , scientific = "Carapa guianensis"               
                    )#end list
   tnf[[ 20]] = list( common     = "andirobarana"              
                    , scientific = "Guarea macrophylla"              
                    )#end list
   tnf[[ 21]] = list( common     = "angelim da mata"           
                    , scientific = "Hymenolobium excelsum"           
                    )#end list
   tnf[[ 22]] = list( common     = "angelim rajado"            
                    , scientific = "Zygia racemosa"                  
                    )#end list
   tnf[[ 23]] = list( common     = "angelim vermelho"          
                    , scientific = "Hymenolobium"                    
                    )#end list
   tnf[[ 24]] = list( common     = "apui"                      
                    , scientific = "Ficus broadwayi"                 
                    )#end list
   tnf[[ 25]] = list( common     = "araca da mata"             
                    , scientific = "Eugenia patrisii"                
                    )#end list
   tnf[[ 26]] = list( common     = "araracanga"                
                    , scientific = "Aspidosperma eteanum"            
                    )#end list
   tnf[[ 27]] = list( common     = "araticum"                  
                    , scientific = "Annona densicoma"                
                    )#end list
   tnf[[ 28]] = list( common     = "aroeira"                   
                    , scientific = "Astronium lecointei"             
                    )#end list
   tnf[[ 29]] = list( common     = "axixa"                     
                    , scientific = "Sterculia speciosa"              
                    )#end list
   tnf[[ 30]] = list( common     = "bacaba"                    
                    , scientific = "Oenocarpus"                      
                    )#end list
   tnf[[ 31]] = list( common     = "bacuri da mata"            
                    , scientific = "Garcinia madruno"                
                    )#end list
   tnf[[ 32]] = list( common     = "breu"                      
                    , scientific = "Protium apiculatum"              
                    )#end list
   tnf[[ 33]] = list( common     = "breu folha grande"         
                    , scientific = "Protium"                         
                    )#end list
   tnf[[ 34]] = list( common     = "breu manga"                
                    , scientific = "Tetragastris altissima"          
                    )#end list
   tnf[[ 35]] = list( common     = "breu sucuruba"             
                    , scientific = "Trattinnickia rhoifolia"         
                    )#end list
   tnf[[ 36]] = list( common     = "breu vermelho"             
                    , scientific = "Protium tenuifolium"             
                    )#end list
   tnf[[ 37]] = list( common     = "cabeca-de-urubu"           
                    , scientific = "Duroia macrophylla"              
                    )#end list
   tnf[[ 38]] = list( common     = "cacau da mata"             
                    , scientific = "Theobroma speciosum"             
                    )#end list
   tnf[[ 39]] = list( common     = "caferana"                  
                    , scientific = "Coussarea albescens"             
                    )#end list
   tnf[[ 40]] = list( common     = "caju da mata"              
                    , scientific = "Anacardium"                      
                    )#end list
   tnf[[ 41]] = list( common     = "canela"                    
                    , scientific = "Ocotea acutangula"               
                    )#end list
   tnf[[ 42]] = list( common     = "canela-de-jacamim"         
                    , scientific = "Rinorea flavescens"              
                    )#end list
   tnf[[ 43]] = list( common     = "canela-de-velho"           
                    , scientific = "Rinorea macrocarpa"              
                    )#end list
   tnf[[ 44]] = list( common     = "caneleira"                 
                    , scientific = "Casearia javitensis"             
                    )#end list
   tnf[[ 45]] = list( common     = "caneleira folha peluda"    
                    , scientific = "Casearia commersoniana"          
                    )#end list
   tnf[[ 46]] = list( common     = "caneleira vermelha"        
                    , scientific = "Matayba purgans"                 
                    )#end list
   tnf[[ 47]] = list( common     = "capitiu"                   
                    , scientific = "Siparuna cristata"               
                    )#end list
   tnf[[ 48]] = list( common     = "caqui"                     
                    , scientific = "Diospyros vestita"               
                    )#end list
   tnf[[ 49]] = list( common     = "carapanauba"               
                    , scientific = "Aspidosperma oblongum"           
                    )#end list
   tnf[[ 50]] = list( common     = "carapanauba amarela"       
                    , scientific = "Aspidosperma auriculatum"        
                    )#end list
   tnf[[ 51]] = list( common     = "castanha-do-para"          
                    , scientific = "Bertholletia excelsa"            
                    )#end list
   tnf[[ 52]] = list( common     = "castanha-sapucaia"         
                    , scientific = "Lecythis pisonis"                
                    )#end list
   tnf[[ 53]] = list( common     = "cedro"                     
                    , scientific = "Cedrela odorata"                 
                    )#end list
   tnf[[ 54]] = list( common     = "chichua"                   
                    , scientific = "Maytenus pruinosa"               
                    )#end list
   tnf[[ 55]] = list( common     = "cipo"                      
                    , scientific = "Liana"                           
                    )#end list
   tnf[[ 56]] = list( common     = "cocao"                     
                    , scientific = "Poecilanthe effusa"              
                    )#end list
   tnf[[ 57]] = list( common     = "copaiba"                   
                    , scientific = "Copaifera reticulata"            
                    )#end list
   tnf[[ 58]] = list( common     = "coracao-de-negro"          
                    , scientific = "Chamaecrista xinguensis"         
                    )#end list
   tnf[[ 59]] = list( common     = "cuiarana"                  
                    , scientific = "Buchenavia grandis"              
                    )#end list
   tnf[[ 60]] = list( common     = "cumaru"                    
                    , scientific = "Dipteryx odorata"                
                    )#end list
   tnf[[ 61]] = list( common     = "cumarui"                   
                    , scientific = "Prunus myrtifolia"               
                    )#end list
   tnf[[ 62]] = list( common     = "cumate preto"              
                    , scientific = "Calyptranthes lucida"            
                    )#end list
   tnf[[ 63]] = list( common     = "cunario"                   
                    , scientific = "Connarus perrottetii"            
                    )#end list
   tnf[[ 64]] = list( common     = "embauba"                   
                    , scientific = "Cecropia"                        
                    )#end list
   tnf[[ 65]] = list( common     = "embauba branca"            
                    , scientific = "Cecropia distachya"              
                    )#end list
   tnf[[ 66]] = list( common     = "embauba vermelha"          
                    , scientific = "Cecropia sciadophylla"           
                    )#end list
   tnf[[ 67]] = list( common     = "embaubarana"               
                    , scientific = "Pourouma guianensis"             
                    )#end list
   tnf[[ 68]] = list( common     = "envira"                    
                    , scientific = "Xylopia"                         
                    )#end list
   tnf[[ 69]] = list( common     = "envira branca"             
                    , scientific = "Guatteria amazonica"             
                    )#end list
   tnf[[ 70]] = list( common     = "envira-cana"               
                    , scientific = "Xylopia nitida"                  
                    )#end list
   tnf[[ 71]] = list( common     = "envira preta"              
                    , scientific = "Guatteria poeppigiana"           
                    )#end list
   tnf[[ 72]] = list( common     = "envira-surucucu"           
                    , scientific = "Duguetia echinophora"            
                    )#end list
   tnf[[ 73]] = list( common     = "envira vermelha"           
                    , scientific = "Xylopia ochrantha"               
                    )#end list
   tnf[[ 74]] = list( common     = "escorrega-macaco"          
                    , scientific = "Capirona decorticans"            
                    )#end list
   tnf[[ 75]] = list( common     = "farinha seca"              
                    , scientific = "Ampelocera edentula"             
                    )#end list
   tnf[[ 76]] = list( common     = "fava"                      
                    , scientific = "Abarema jupunba"                 
                    )#end list
   tnf[[ 77]] = list( common     = "fava amargosa"             
                    , scientific = "Vataireopsis speciosa"           
                    )#end list
   tnf[[ 78]] = list( common     = "fava-bolota"               
                    , scientific = "Parkia pendula"                  
                    )#end list
   tnf[[ 79]] = list( common     = "fava da rosca"             
                    , scientific = "Enterolobium schomburgkii"       
                    )#end list
   tnf[[ 80]] = list( common     = "fava folha fina"           
                    , scientific = "Pseudopiptadenia psilostachy"    
                    )#end list
   tnf[[ 81]] = list( common     = "fava-mapuxiqui"            
                    , scientific = "Albizia pedicellaris"            
                    )#end list
   tnf[[ 82]] = list( common     = "fava-saboeiro"             
                    , scientific = "Abarema"                         
                    )#end list
   tnf[[ 83]] = list( common     = "fava-timbauba"             
                    , scientific = "Enterolobium maximum"            
                    )#end list
   tnf[[ 84]] = list( common     = "fava-arara-tucupi"         
                    , scientific = "Parkia multijuga"                
                    )#end list
   tnf[[ 85]] = list( common     = "fava-barbatimao"           
                    , scientific = "Stryphnodendron pulcherrimum"    
                    )#end list
   tnf[[ 86]] = list( common     = "fava-saboeiro"             
                    , scientific = "Abarema"                         
                    )#end list
   tnf[[ 87]] = list( common     = "freijo"                    
                    , scientific = "Cordia exaltata"                 
                    )#end list
   tnf[[ 88]] = list( common     = "freijo branco"             
                    , scientific = "Cordia bicolor"                  
                    )#end list
   tnf[[ 89]] = list( common     = "freijorana"                
                    , scientific = "Cordia"                          
                    )#end list
   tnf[[ 90]] = list( common     = "geniparana"                
                    , scientific = "Gustavia poeppigiana"            
                    )#end list
   tnf[[ 91]] = list( common     = "ginja"                     
                    , scientific = "Eugenia omissa"                  
                    )#end list
   tnf[[ 92]] = list( common     = "goiabarana"                
                    , scientific = "Eugenia omissa"                  
                    )#end list
   tnf[[ 93]] = list( common     = "goiabinha"                 
                    , scientific = "Myrciaria tenella"               
                    )#end list
   tnf[[ 94]] = list( common     = "gombeira"                  
                    , scientific = "Swartzia laurifolia"             
                    )#end list
   tnf[[ 95]] = list( common     = "gombeira folha peluda"     
                    , scientific = "Swartzia laxiflora"              
                    )#end list
   tnf[[ 96]] = list( common     = "gombeira vermelha"         
                    , scientific = "Swartzia laurifolia"             
                    )#end list
   tnf[[ 97]] = list( common     = "guariuba"                  
                    , scientific = "Clarisia racemosa"               
                    )#end list
   tnf[[ 98]] = list( common     = "inajarana"                 
                    , scientific = "Quararibea guianensis"           
                    )#end list
   tnf[[ 99]] = list( common     = "inga"                      
                    , scientific = "Inga macrophylla"                
                    )#end list
   tnf[[100]] = list( common     = "inga branco"               
                    , scientific = "Inga alba"                       
                    )#end list
   tnf[[101]] = list( common     = "inga vermelho"             
                    , scientific = "Inga alba"                       
                    )#end list
   tnf[[102]] = list( common     = "inga xixica"               
                    , scientific = "Inga"                            
                    )#end list
   tnf[[103]] = list( common     = "itauba amarela"            
                    , scientific = "Mezilaurus itauba"               
                    )#end list
   tnf[[104]] = list( common     = "itaubarana"                
                    , scientific = "Casearia"                        
                    )#end list
   tnf[[105]] = list( common     = "amapa amargoso"                    
                    , scientific = "Brosimum guianense"              
                    )#end list
   tnf[[106]] = list( common     = "jarana"                    
                    , scientific = "Lecythis lurida"                 
                    )#end list
   tnf[[107]] = list( common     = "jatauba"                   
                    , scientific = "Matayba purgans"                 
                    )#end list
   tnf[[108]] = list( common     = "joao-mole"                 
                    , scientific = "Guapira venosa"                  
                    )#end list
   tnf[[109]] = list( common     = "joao-mole grande"          
                    , scientific = "Neea"                            
                    )#end list
   tnf[[110]] = list( common     = "jutai-pororoca"            
                    , scientific = "Dialium guianense"               
                    )#end list
   tnf[[111]] = list( common     = "jutai-acu"                 
                    , scientific = "Hymenaea courbaril"              
                    )#end list
   tnf[[112]] = list( common     = "jutai-mirim"               
                    , scientific = "Hymenaea parvifolia"             
                    )#end list
   tnf[[113]] = list( common     = "jutaiarana"                
                    , scientific = "Swartzia arborescens"            
                    )#end list
   tnf[[114]] = list( common     = "lacre-da-mata"             
                    , scientific = "Vismia"                          
                    )#end list
   tnf[[115]] = list( common     = "lacre vermelho"            
                    , scientific = "Vismia latifolia"                
                    )#end list
   tnf[[116]] = list( common     = "louro"                     
                    , scientific = "Nectandra pulverulenta"          
                    )#end list
   tnf[[117]] = list( common     = "louro-abacate"             
                    , scientific = "Ocotea glomerata"                
                    )#end list
   tnf[[118]] = list( common     = "louro amarelo"             
                    , scientific = "Licaria brasiliensis"            
                    )#end list
   tnf[[119]] = list( common     = "louro mole"                
                    , scientific = "Nectandra pulverulenta"          
                    )#end list
   tnf[[120]] = list( common     = "louro preto"               
                    , scientific = "Nectandra reticulata"            
                    )#end list
   tnf[[121]] = list( common     = "louro vermelho"            
                    , scientific = "Ocotea cuprea"                   
                    )#end list
   tnf[[122]] = list( common     = "macacauba"                 
                    , scientific = "Platymiscium trinitatis"         
                    )#end list
   tnf[[123]] = list( common     = "macaranduba"               
                    , scientific = "Manilkara huberi"                
                    )#end list
   tnf[[124]] = list( common     = "macucu"                    
                    , scientific = "Licania heteromorpha"            
                    )#end list
   tnf[[125]] = list( common     = "mamorana"                  
                    , scientific = "Eriotheca globosa"               
                    )#end list
   tnf[[126]] = list( common     = "mandioqueira ariana"       
                    , scientific = "Qualea grandiflora"              
                    )#end list
   tnf[[127]] = list( common     = "mandioqueira rosa"         
                    , scientific = "Qualea dinizii"                  
                    )#end list
   tnf[[128]] = list( common     = "maparana"                  
                    , scientific = "Drypetes variabilis"             
                    )#end list
   tnf[[129]] = list( common     = "marupa"                    
                    , scientific = "Simarouba amara"                 
                    )#end list
   tnf[[130]] = list( common     = "mata-calado"               
                    , scientific = "Lacistema aggregatum"            
                    )#end list
   tnf[[131]] = list( common     = "matamata"                  
                    , scientific = "Lecythis holcogyne"              
                    )#end list
   tnf[[132]] = list( common     = "matamata branco"           
                    , scientific = "Eschweilera pedicellata"         
                    )#end list
   tnf[[133]] = list( common     = "matamata jarani"           
                    , scientific = "Lecythis holcogyne"              
                    )#end list
   tnf[[134]] = list( common     = "matamata preto"            
                    , scientific = "Eschweilera coriacea"            
                    )#end list
   tnf[[135]] = list( common     = "matamata vermelho"         
                    , scientific = "Eschweilera bracteosa"           
                    )#end list
   tnf[[136]] = list( common     = "melancieira"               
                    , scientific = "Alexa grandiflora"               
                    )#end list
   tnf[[137]] = list( common     = "mirindiba doce"            
                    , scientific = "Glycydendron amazonicum"         
                    )#end list
   tnf[[138]] = list( common     = "morototo"                  
                    , scientific = "Schefflera morototoni"           
                    )#end list
   tnf[[139]] = list( common     = "muiauba"                   
                    , scientific = "Mouriri nigra"                   
                    )#end list
   tnf[[140]] = list( common     = "muiracatiara"              
                    , scientific = "Astronium lecointei"             
                    )#end list
   tnf[[141]] = list( common     = "muirapinima"               
                    , scientific = "Maquira calophylla"              
                    )#end list
   tnf[[142]] = list( common     = "muirarema"                 
                    , scientific = "Trichilia micrantha"             
                    )#end list
   tnf[[143]] = list( common     = "muiratinga"                
                    , scientific = "Maquira guianensis"              
                    )#end list
   tnf[[144]] = list( common     = "muiratinga folha larga"    
                    , scientific = "Maquira sclerophylla"            
                    )#end list
   tnf[[145]] = list( common     = "muiratinga folha lisa"     
                    , scientific = "Maquira sclerophylla"            
                    )#end list
   tnf[[146]] = list( common     = "muiratinga folha longa"    
                    , scientific = "Maquira sclerophylla"            
                    )#end list
   tnf[[147]] = list( common     = "muiratinga folha peluda"   
                    , scientific = "Perebea mollis"                  
                    )#end list
   tnf[[148]] = list( common     = "muirauba"                  
                    , scientific = "Mouriri brachyanthera"           
                    )#end list
   tnf[[149]] = list( common     = "murta"                     
                    , scientific = "Myrcia fallax"                   
                    )#end list
   tnf[[150]] = list( common     = "muruci-da-mata"            
                    , scientific = "Byrsonima arthropoda"            
                    )#end list
   tnf[[151]] = list( common     = "murure"                    
                    , scientific = "Brosimum acutifolium"            
                    )#end list
   tnf[[152]] = list( common     = "mututi"                    
                    , scientific = "Pterocarpus rohrii"              
                    )#end list
   tnf[[153]] = list( common     = "muuba"                     
                    , scientific = "Bellucia grossularioides"        
                    )#end list
   tnf[[154]] = list( common     = "pama"                      
                    , scientific = "Pseudolmedia macrophylla"        
                    )#end list
   tnf[[155]] = list( common     = "papa-terra"                 
                    , scientific = "Miconia ruficalyx"               
                    )#end list
   tnf[[156]] = list( common     = "papa-terra amarelo"         
                    , scientific = "Miconia lepidota"                
                    )#end list
   tnf[[157]] = list( common     = "papa-terra folha peluda"    
                    , scientific = "Miconia phanerostila"            
                    )#end list
   tnf[[158]] = list( common     = "parapara"                  
                    , scientific = "Jacaranda copaia"                
                    )#end list
   tnf[[159]] = list( common     = "passarinheira"             
                    , scientific = "Casearia ulmifolia"              
                    )#end list
   tnf[[160]] = list( common     = "pata-de-vaca"              
                    , scientific = "Bauhinia"                        
                    )#end list
   tnf[[161]] = list( common     = "papo-de-mutum"             
                    , scientific = "Lacunaria crenata"               
                    )#end list
   tnf[[162]] = list( common     = "pau-cobra"                 
                    , scientific = "Salacia impressifolia"           
                    )#end list
   tnf[[163]] = list( common     = "pau-de-arco amarelo"       
                    , scientific = "Handroanthus serratifolius"           
                    )#end list
   tnf[[164]] = list( common     = "pau-de-colher"             
                    , scientific = "Lacmellea aculeata"              
                    )#end list
   tnf[[165]] = list( common     = "pau-de-remo"               
                    , scientific = "Chimarrhis turbinata"            
                    )#end list
   tnf[[166]] = list( common     = "pau-jacare"                
                    , scientific = "Laetia procera"                  
                    )#end list
   tnf[[167]] = list( common     = "pau-para-tudo"             
                    , scientific = "Swartzia recurva"                
                    )#end list
   tnf[[168]] = list( common     = "pente-de-macaco"           
                    , scientific = "Apeiba glabra"                   
                    )#end list
   tnf[[169]] = list( common     = "molongo"            
                    , scientific = "Ambelania acida"                 
                    )#end list
   tnf[[170]] = list( common     = "piquia"                    
                    , scientific = "Caryocar villosum"               
                    )#end list
   tnf[[171]] = list( common     = "pitomba"                   
                    , scientific = "Talisia retusa"                  
                    )#end list
   tnf[[172]] = list( common     = "quariquarana"              
                    , scientific = "Rinorea"                         
                    )#end list
   tnf[[173]] = list( common     = "quaruba"                   
                    , scientific = "Vochysia"                        
                    )#end list
   tnf[[174]] = list( common     = "quaruba rosa"              
                    , scientific = "Vochysia surinamensis"           
                    )#end list
   tnf[[175]] = list( common     = "quaruba verdadeira"        
                    , scientific = "Vochysia maxima"                 
                    )#end list
   tnf[[176]] = list( common     = "quarubarana"               
                    , scientific = "Erisma uncinatum"                
                    )#end list
   tnf[[177]] = list( common     = "quinarana"                 
                    , scientific = "Geissospermum laeve"             
                    )#end list
   tnf[[178]] = list( common     = "saboeira"                  
                    , scientific = "Abarema jupunba"                 
                    )#end list
   tnf[[179]] = list( common     = "sucupira amarela"          
                    , scientific = "Bowdichia nitida"                
                    )#end list
   tnf[[180]] = list( common     = "sucupira preta"            
                    , scientific = "Diplotropis triloba"             
                    )#end list
   tnf[[181]] = list( common     = "tachi"                     
                    , scientific = "Tachigali"                       
                    )#end list
   tnf[[182]] = list( common     = "tachi branco"              
                    , scientific = "Tachigali alba"                  
                    )#end list
   tnf[[183]] = list( common     = "tachi preto"               
                    , scientific = "Tachigali myrmecophila"          
                    )#end list
   tnf[[184]] = list( common     = "tachi preto folha grauda"  
                    , scientific = "Tachigali myrmecophila"          
                    )#end list
   tnf[[185]] = list( common     = "tachi preto folha miuda"   
                    , scientific = "Tachigali paniculata"            
                    )#end list
   tnf[[186]] = list( common     = "tachi vermelho"            
                    , scientific = "Tachigali chrysophylla"          
                    )#end list
   tnf[[187]] = list( common     = "taperebarana"              
                    , scientific = "Heisteria laxiflora"             
                    )#end list
   tnf[[188]] = list( common     = "taruma"                    
                    , scientific = "Vitex triflora"                  
                    )#end list
   tnf[[189]] = list( common     = "tatajuba"                  
                    , scientific = "Bagassa guianensis"              
                    )#end list
   tnf[[190]] = list( common     = "tatapiririca"              
                    , scientific = "Tapirira guianensis"             
                    )#end list
   tnf[[191]] = list( common     = "tatapiririca folha peluda" 
                    , scientific = "Tapirira obtusa"                 
                    )#end list
   tnf[[192]] = list( common     = "tatapiririca vermelha"     
                    , scientific = "Tapirira guianensis"             
                    )#end list
   tnf[[193]] = list( common     = "tauari"                    
                    , scientific = "Couratari stellata"              
                    )#end list
   tnf[[194]] = list( common     = "tento"                     
                    , scientific = "Simaba polyphylla"               
                    )#end list
   tnf[[195]] = list( common     = "tento folha grauda"        
                    , scientific = "Ormosia paraensis"               
                    )#end list
   tnf[[196]] = list( common     = "tento folha miuda"         
                    , scientific = "Abarema mataybifolia"            
                    )#end list
   tnf[[197]] = list( common     = "tento preto"               
                    , scientific = "Abarema mataybifolia"            
                    )#end list
   tnf[[198]] = list( common     = "tres folhas"               
                    , scientific = "Allophylus punctatus"            
                    )#end list
   tnf[[199]] = list( common     = "ucuuba"                    
                    , scientific = "Virola"                          
                    )#end list
   tnf[[200]] = list( common     = "ucuuba folha peluda"       
                    , scientific = "Virola crebrinervia"             
                    )#end list
   tnf[[201]] = list( common     = "ucuuba vermelha"           
                    , scientific = "Virola elongata"                 
                    )#end list
   tnf[[202]] = list( common     = "ucuuba terra-firme"        
                    , scientific = "Virola michelii"                 
                    )#end list
   tnf[[203]] = list( common     = "ucuubarana"                
                    , scientific = "Iryanthera sagotiana"            
                    )#end list
   tnf[[204]] = list( common     = "uruazeiro"                 
                    , scientific = "Cordia exaltata"                 
                    )#end list
   tnf[[205]] = list( common     = "urucurana"                 
                    , scientific = "Aparisthmium cordatum"           
                    )#end list
   tnf[[206]] = list( common     = "uxi"                       
                    , scientific = "Endopleura"                      
                    )#end list
   tnf[[207]] = list( common     = "uxi liso"                  
                    , scientific = "Endopleura uchi"                 
                    )#end list
   tnf[[208]] = list( common     = "verdinho"                  
                    , scientific = "Ilex petiolaris"                 
                    )#end list
   tnf[[209]] = list( common     = "cipo"                  
                    , scientific = "Liana"                 
                    )#end list
   #---------------------------------------------------------------------------------------#



   #----- Convert the TNF list into a data frame. -----------------------------------------#
   tnf = data.frame( apply( X      = t(sapply(X = tnf , FUN = c))
                          , MARGIN = c(1,2)
                          , FUN    = unlist
                          )#end apply
                   , stringsAsFactors = FALSE
                   )#end data.frame
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Create an output data frame with the scientific names.  The default "no-name"     #
   # data does not count as gap-filling.                                                   #
   #---------------------------------------------------------------------------------------#
   out            = data.frame( common           = rep("mato"      ,times=n.datum)
                              , scientific       = rep("Ignotum"   ,times=n.datum)
                              , stringsAsFactors = FALSE
                              )#end data.frame
   idx            = match(datum$common,tnf$common)
   ok             = ! is.na(idx)
   out[ok,]       = tnf[idx[ok],]
   out$gfflg      = as.numeric(! out$common %in% c("cipo","mato"))
   #---------------------------------------------------------------------------------------#
   

   #---------------------------------------------------------------------------------------#
   #      Copy data to the main data frame.                                                #
   #---------------------------------------------------------------------------------------#
   sel = is.na(datum$scientific)
   datum$common       [sel] = out$common    [sel]
   datum$scientific   [sel] = out$scientific[sel]
   datum$gf.scientific[sel] = out$gfflg     [sel]
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Put NA to species if only genus is known.                                         #
   #---------------------------------------------------------------------------------------#
   gs.list          = sapply(X = tolower(datum$scientific),FUN=strsplit,split=" ")
   gs.length        = sapply(X = gs.list, FUN = length)
   gs.mat           = cbind( mapply(FUN="[",gs.list,MoreArgs=list(1))
                           , mapply(FUN="[",gs.list,MoreArgs=list(2))
                           )#end cbind
   datum$genus      = capwords(gs.mat[,1],strict=TRUE)
   datum$scientific = paste(datum$genus,tolower (gs.mat[,2]),sep=" ")
   #---------------------------------------------------------------------------------------#
   

   #---------------------------------------------------------------------------------------#
   #      Fill in the family.                                                              #
   #---------------------------------------------------------------------------------------#
   datum = standard.family.name(datum)
   #---------------------------------------------------------------------------------------#

   return(datum)
}#end if
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     Fill in the wood density for all individuals.  We do this in three stages, and       #
# save how the wood density was determined.  The numbers represent the flag given for      #
# each wood density.                                                                       #
#  0 -- Individuals have full identification and the species is listed in the              #
#       database; we used the reported density for that species.                           #
#  1 -- The species/genus is identified but it isn't listed in the database, we use an     #
#       average of all other species of that genus that exist in database.                 #
#  2 -- We could not identify the individual to the genus level, but we know the           #
#       family.  We use the average wood density of all individuals of this census that    #
#       belong to that family.                                                             #
#  3 -- No taxonomic information could be retrieved for this individual, so we filled      #
#       with random sampling.                                                              #
#                                                                                          #
#  We also add 10 when the scientific name was gap filled.                                 #
#                                                                                          #
#  INPUT variables:                                                                        #
#                                                                                          #
#  - datum       -- data frame with data.  This is going to be the output as well          #
#  - wood        -- wood density data base.                                                #
#  - region      -- which regions to use for genus average when species is not known or    #
#                   not available in the data set.                                         #
#  - fill.sample -- use random sampling to fill unidentified individuals, or individuals   #
#                   with genus that is not available at the wood density data base?        #
#                   If FALSE then it uses averages                                         #
#  - weight      -- A weighting factor to give either probability or to weight             #
#                   the average.  This could be a vector with weights, or a character with #
#                   the name of the variable in datum to use as the weight, or an integer  #
#                   with the column to be used as the weighting factor.                    #
#  - verbose     -- Flag to control the amount of information                              #
#------------------------------------------------------------------------------------------#
find.wood.density <<- function( datum
                              , wood
                              , region      = NULL
                              , fill.sample = TRUE
                              , weight      = NULL
                              , verbose     = FALSE
                              ){
   #---------------------------------------------------------------------------------------#
   #     Initialise gap filling flag and region in case they are not there.                #
   #---------------------------------------------------------------------------------------#
   if (! "gf.wood.dens" %in% names(datum)){
      datum$gf.wood.dens = rep(NA,times=nrow(datum))
   }#end (! "gf.wood.dens" %in% names(datum))
   if (! "region" %in% names(datum)){
      datum$region = rep(NA_character_,times=nrow(datum))
   }#end (! "gf.wood.dens" %in% names(datum))
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Check whether weight is a vector, or a character.  Make them a vector here.       #
   #---------------------------------------------------------------------------------------#
   if (is.null(weight)){
      #----- No weight provided, use equal weights. ---------------------------------------#
      wgtfac = rep(x=1/nrow(datum),times=nrow(datum))
      #------------------------------------------------------------------------------------#
   }else if (is.character(weight) && (length(weight) == 1)){
      #----- Character with column name was provided. -------------------------------------#
      if (weight %in% names(datum)){
         wgtfac = ifelse(datum[[weight]] %>% 0, datum[[weight]], 0)
      }else{
         stop(paste0(" Weight Variable name (",weight,") not found in datum!"))
      }#end if
      #------------------------------------------------------------------------------------#
   }else if (is.numeric(weight) && (length(weight) == 1)){
      #----- Column index provided. -------------------------------------------------------#
      if (! (weight %wr% c(1,ncol(datum)))){
         stop(paste0(" Weight column index (",weight,") doesn't make sense"))
      }else if (is.numeric(datum[,weight])){
         wgtfac = ifelse(datum[,weight] %>% 0, datum[,weight], 0)
      }else{
         stop(paste0(" Column ",weight," of data frame is not numeric!"))
      }#end if
      #------------------------------------------------------------------------------------#
   }else if (is.numeric(weight) && (length(weight) == nrow(datum))){
      wgtfac = ifelse(weight %>% 0, weight, 0)
   }else{
      stop("Variable weight is not properly set!")
   }#end if (is.null(weight))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Check that the weighting factor makes sense.                                      #
   #---------------------------------------------------------------------------------------#
   if (sum(wgtfac) > 0){
      wgtfac = wgtfac / sum(wgtfac)
   }else{
      stop(" Invalid weighting variable. Most numbers should be positive...")
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Find out which region to use.                                                    #
   #---------------------------------------------------------------------------------------#
   if (! "region" %in% names(wood)){
      stop(" You must be using an old wood density file.  You need region information")
   }else{
      all.regions = sort(unique(wood$region))
      if (is.null(region)){
         #----- No region has been given, use all regions. --------------------------------#
         region = all.regions
         #---------------------------------------------------------------------------------#
      }else{
         #----- Make sure all regions are valid. ------------------------------------------#
         if (! all(region %in% all.regions)){
            cat0("--------------------------------------------------------")
            cat0(" - Your region: ",region)
            cat0(" - Acceptable regions: ",paste(all.regions,collapse="; "))
            cat0("--------------------------------------------------------")
            stop(" Please only include regions listed in the data base.")
         }#end if
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     First loop, fill in information to all individuals for which the genus is         #
   # known.                                                                                #
   #---------------------------------------------------------------------------------------#
   #----- Separate all species. -----------------------------------------------------------#
   species          = unique(datum$scientific)
   nspecies         = length(species)
   
   sci.genus.mean   = matrix(nrow=0,ncol=3
                            ,dimnames=list(NULL,c("scientific","genus","family")))
   sci.loose.mean   = matrix(nrow=0,ncol=3
                            ,dimnames=list(NULL,c("scientific","genus","family")))
   for (s in sequence(nspecies)){
      if (! is.na(species[s]) && length(grep("Ignotum",species[s])) == 0){
         if (verbose) cat("   - ",s,"/",nspecies,"  -  ",species[s],"...","\n") 

         #----- Get the genus and family of this species, in case we need it. -------------#
         igen        = which (datum$scientific %in% species[s])
         this.genus  = unique(datum$genus[igen])
         this.family = unique(capwords(datum$family[igen],strict=TRUE))
         if (length(this.genus) != 1 || length(this.family) != 1){
            cat(" - Perhaps a bastard genus?","\n")
            browser()
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Check whether we have biomass for this species.                             #
         #---------------------------------------------------------------------------------#
         if (species[s] %in% wood$scientific){
            #------------------------------------------------------------------------------#
            #     Find the index of this species in the database, and assign the wood      #
            # density from there.                                                          #
            #------------------------------------------------------------------------------#
            iwood = intersect( which(wood$scientific %in% species[s] )
                             , intersect( which( wood$genus  %in% this.genus )
                                        , which( wood$family %in% this.family) ) )
            if (length(iwood) != 1){
               #----- Intersection was zero! Misidentification, maybe? --------------------#
               loose      = TRUE
               cat ("Weird, length(iwood)=",length(iwood),"...","\n",sep="")
               browser()
               #---------------------------------------------------------------------------#
            }else{
               if (any(c(FALSE,is.finite(wood$density[iwood])),na.rm=TRUE)){
                  sel       = datum$scientific %in% species[s]
                  datum$wood.dens   [sel] = wood$density[iwood]
                  datum$gf.wood.dens[sel] = 0
                  datum$region      [sel] = wood$region [iwood]
                  fill.genus              = FALSE
                  loose                   = FALSE
               }else{
                  loose                   = TRUE
                  fill.genus              = this.genus %in% wood$genus
               }#end if
               #---------------------------------------------------------------------------#
            }#end if
            #------------------------------------------------------------------------------#
         }else{
            loose      = TRUE
            fill.genus = this.genus %in% wood$genus
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     If this species is to be filled with genus average.                         #
         #---------------------------------------------------------------------------------#
         if (fill.genus){
            #------------------------------------------------------------------------------#
            #     Find all the plants from this genus, take the average, and attribute     #
            # to this species.  In this case, warn the user about the species, it may      #
            # be a typo in the scientific name that is easy to fix.                        #
            #------------------------------------------------------------------------------#
            iwood  = intersect( which( wood$genus  %in% this.genus )
                              , intersect( which( wood$family %in% this.family)
                                         , which( wood$region %in% region     )
                                         )#end intersect
                              )#end intersect
            if (length(iwood) == 0){
               #----- Intersection was zero! Misidentification, maybe? --------------------#
               loose      = TRUE
            }else{
               if (any(is.finite(wood$density[iwood]))){
                  sel                     = datum$scientific   %in% species[s]
                  overwrite               = sel & ! is.na(datum$wood.dens)
                  if (any(overwrite)) browser()
 
                  datum$wood.dens   [sel] = mean(wood$density[iwood],na.rm=TRUE)
                  datum$gf.wood.dens[sel] = 1
                  datum$region      [sel] = commonest(wood$region[iwood],na.rm=TRUE)
                  sci.genus.mean          = rbind(sci.genus.mean
                                                 ,c(species[s],this.genus,this.family))
                  if (verbose){
                     cat("     * Species",species[s]
                        ,"not found. Use genus average instead...","\n")
                  }#end if
                  loose = FALSE
               }else{
                  #---- The plant probably doesn't have any family. -----------------------#
                  loose = TRUE
               }#end if
               #---------------------------------------------------------------------------#
            }#end if
            #------------------------------------------------------------------------------#
         }#end if
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #    Append plants that have no family together.                                  #
         #---------------------------------------------------------------------------------#
         if (loose){
            #---- The plant probably doesn't have any family. -----------------------------#
            sci.loose.mean = rbind(sci.loose.mean
                                  ,c(species[s],this.genus,this.family))
            #------------------------------------------------------------------------------#
         }#end if
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     List all genera that were not filled with species information.                    #
   #---------------------------------------------------------------------------------------#
   sci.genus.mean = sci.genus.mean[order(sci.genus.mean[,"scientific"]),]
   genus.only     = grepl(pattern=" NA$",x=sci.genus.mean[,"scientific"],ignore.case=TRUE)
   if (verbose && (nrow(sci.genus.mean[!genus.only,,drop=FALSE]) > 0)){
      cat0("")
      cat0("-----------------------------------------------------------------------")
      cat0(" Found species that were not in Zanne's data base (check for synonyms)!")
      print(sci.genus.mean[! genus.only,,drop=FALSE],quote=FALSE)
      cat0("-----------------------------------------------------------------------")
      cat0("")
   }#end if (verbose && nrow(sci.genus.mean[!genus.only,,drop=FALSE]) > 0)
   if (verbose && (nrow(sci.genus.mean[genus.only,,drop=FALSE]) > 0)){
      cat0("")
      cat0("-----------------------------------------------------------------------")
      cat0("Only genus was provided: we used the genus mean wood density:")
      print (sci.genus.mean[genus.only,,drop=FALSE],quote=FALSE)
      cat0("-----------------------------------------------------------------------")
      cat0("")
   }#end if (verbose && nrow(sci.genus.mean[genus.only,,drop=FALSE]) > 0)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     List all genera that didn't belong to any known family.                           #
   #---------------------------------------------------------------------------------------#
   if (verbose && (nrow(sci.loose.mean) > 0)){
      cat0("")
      cat0("-----------------------------------------------------------------------")
      cat0(" Found genera that didn't have any data in Zanne's data base!")
      print(sci.loose.mean,quote=FALSE)
      cat0("-----------------------------------------------------------------------")
      cat0("")
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Second loop: we list all families, and look for individuals that have no genus    #
   # associated. We compute the mean wood density of the known individuals for that        #
   # family and use that as an estimate of wood density.                                   #
   #---------------------------------------------------------------------------------------#
   #----- Separate all families. ----------------------------------------------------------#
   families          = unique(datum$family)
   nfamilies         = length(families)
   sci.family.sample = matrix(nrow=0,ncol=3
                             ,dimnames=list(NULL,c("scientific","family","wood.dens")))
   #----- Loop over all families. ---------------------------------------------------------#
   for (f in sequence(nfamilies)){
      if (families[f] != "Ignotaceae"){
         if (verbose) cat("   - ",f,"/",nfamilies,"  -  ",families[f],"...","\n")

         #----- Get the individuals that belong to this family. ---------------------------#
         ifam  = which (datum$family %in% families[f] & is.finite(datum$wood.dens) )
         imiss = which (datum$family %in% families[f] & is.na    (datum$wood.dens) )
         if (length(imiss) > 0 && length(ifam) > 0){
            #----- Decide whether to use sample or average. -------------------------------#
            if (fill.sample){
               sample.wood.dens       = lit.sample( x       = datum$wood.dens[ifam]
                                                  , size    = length(imiss)
                                                  , replace = TRUE
                                                  , prob    = wgtfac[ifam]
                                                  )#end sample
               sample.region          = lit.sample( x       = datum$region[ifam]
                                                  , size    = length(imiss)
                                                  , replace = TRUE
                                                  , prob    = wgtfac[ifam]
                                                  )#end sample
            }else{
               sample.wood.dens       = weighted.mean( x     = datum$wood.dens[ifam]
                                                     , w     = wgtfac[ifam]
                                                     , na.rm = TRUE
                                                     )#end weighted.mean
               sample.region          = weighted.commonest( x     = datum$region[ifam]
                                                          , w     = wgtfac[ifam]
                                                          , na.rm = TRUE
                                                          )#end weighted.commonest
            }#end if
            #------------------------------------------------------------------------------#


            #------ Fill in with the sample/mean. -----------------------------------------#
            datum$wood.dens   [imiss] = sample.wood.dens
            datum$gf.wood.dens[imiss] = 2
            datum$region      [imiss] = sample.region
            gf2                       = cbind(datum$scientific[imiss]
                                             ,datum$family    [imiss]
                                             ,sprintf("%6.3f",sample.wood.dens)
                                             )#end cbind
            sci.family.sample         = rbind(sci.family.sample,gf2)
            #------------------------------------------------------------------------------#
         }#end if
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Stop if there was any genus that didn't belong to any known family.               #
   #---------------------------------------------------------------------------------------#
   if (nrow(sci.family.sample) > 0){
      if (verbose){
         cat (" Found families with unidentified genera!","\n")
         print(sci.family.sample,quote=FALSE)
      }#end if
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Final block.  We fill in the wood density for unknown individuals, by randomly     #
   # sampling from the individuals we know the density.                                    #
   #---------------------------------------------------------------------------------------#
   imiss                     = which(is.na(datum$wood.dens))
   nmiss                     = length(imiss)
   if (nmiss > 0){
      if (verbose){
         cat (" The following families are filled with global sampling: ","\n")
         fam.global.sampling = t(t(sort(unique(datum$family[imiss]))))
         print(fam.global.sampling)
      }#end if

      #----- Decide whether to use sample or averaged values. -----------------------------#
      if (fill.sample){
         sample.wood.dens = lit.sample( x       = datum$wood.dens[-imiss]
                                      , size    = nmiss
                                      , replace = TRUE
                                      , prob    = wgtfac[-imiss]
                                      )#end lit.sample
         sample.region    = lit.sample( x       = datum$region[-imiss]
                                      , size    = nmiss
                                      , replace = TRUE
                                      , prob    = wgtfac[-imiss]
                                      )#end lit.sample
      }else{
         sample.wood.dens = weighted.mean( x     = datum$wood.dens[-imiss]
                                         , w     = wgtfac         [-imiss]
                                         , na.rm = TRUE
                                         )#end weighted.mean
         sample.region    = weighted.commonest( x     = datum$region   [-imiss]
                                              , w     = wgtfac         [-imiss]
                                              , na.rm = TRUE
                                              )#end weighted.commonest
      }#end if
      #------------------------------------------------------------------------------------#


      #----- Fill the remaining gaps. -----------------------------------------------------#
      datum$wood.dens   [imiss] = sample.wood.dens
      datum$region      [imiss] = sample.region
      datum$gf.wood.dens[imiss] = 3
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Adjust the gap-filling flag for wood density by adding whether the scientific     #
   # name is gap-filled.                                                                   #
   #---------------------------------------------------------------------------------------#
   if ("gf.scientific" %in% names(datum)){
      datum$gf.wood.dens = datum$gf.wood.dens + 10 * datum$gf.scientific
   }#end if ("gf.scientific" %in% names(datum))
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Assign a plant functional type based on the wood density.                         #
   #---------------------------------------------------------------------------------------#
   pft.cut    = cut(datum$wood.dens,breaks=pft.breaks)
   pft.levels = levels(pft.cut)
   pft.idx    = match(pft.cut,pft.levels)
   datum$pft  = mypfts[pft.idx]
   #---------------------------------------------------------------------------------------#
   return(datum)
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      This attributes scientific names based on common names for surveys in Manaus.  It   #
# is not a good idea to use this anywhere else because common names may mean completely    #
# dif
#==========================================================================================#
#==========================================================================================#
scientific.lookup.mao <<- function(datum,lookup.path){

   #----- Read in the look-up table. ------------------------------------------------------#
   lookup.file = paste(lookup.path,"manaus_taxon_lookup.csv",sep="/")
   look.up = read.csv(file=lookup.file,stringsAsFactors=FALSE)
   look.up$common     = tolower(trim(look.up$common    ))
   look.up$scientific = trim(look.up$scientific)
   look.up$family     = trim(look.up$family    )
   #---------------------------------------------------------------------------------------#



   #----- Break into genus and species. ---------------------------------------------------#
   gs.list                  = sapply(X = tolower(look.up$scientific),FUN=strsplit,split=" ")
   gs.length                = sapply(X = gs.list, FUN = length)
   gs.mat                   = cbind( mapply(FUN="[",gs.list,MoreArgs=list(1))
                                   , mapply(FUN="[",gs.list,MoreArgs=list(2))
                                   )#end cbind
   g                        = capwords(gs.mat[,1],strict=TRUE)
   s                        = tolower(gs.mat[,2])
   g.s                      = paste(g,s,sep=" ")
   g.s[is.na(g) & is.na(s)] = NA
   look.up$scientific       = g.s
   look.up$genus            = g
   #---------------------------------------------------------------------------------------#



   #----- Trim the common names, and simplify/replace some names. -------------------------#
   datum$common                      = tolower(trim(datum$common))
   datum$common[is.na(datum$common)] = "mato"
   #---------------------------------------------------------------------------------------#



   #----- Find all unique common names. ---------------------------------------------------#
   unique.common = unique(datum$common)
   n.common      = length(unique.common)
   notfound      = NULL
   for (n in 1:n.common){
      #----- Find the trees that have the same common name as this one. -------------------#
      cat (" - ",n,"/",n.common," -- ",unique.common[n],"\n")
      w.dat  = which(datum$common %in% unique.common[n])
      n.dat  = length(w.dat)
      #------------------------------------------------------------------------------------#


      #----- Find the trees in the look-up table with the same common name. ---------------#
      w.look = which(look.up$common %in% unique.common[n])
      n.look = length(w.look)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Check how many trees have the same common name in the look-up table.          #
      #------------------------------------------------------------------------------------#
      if (n.look == 1){
         #----- Only one.  Use it. --------------------------------------------------------#
         datum$scientific   [w.dat] = look.up$scientific[w.look]
         datum$genus        [w.dat] = look.up$genus     [w.look]
         datum$gf.scientific[w.dat] = 1
      }else if (n.look > 1){
         datum$scientific   [w.dat] = sample( x       = look.up$scientific[w.look]
                                            , size    = n.dat
                                            , replace = TRUE
                                            )#end sample
         datum$genus        [w.dat] = look.up$genus     [w.look]
         datum$gf.scientific[w.dat] = 1
      }else{
         notfound = c(notfound,unique.common[n])
         datum$scientific   [w.dat] = "Ignotum"
         datum$genus        [w.dat] = "Ignotum"
         datum$gf.scientific[w.dat] = 0
      }#end if
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#
   return(datum)
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      This attributes scientific names based on common names for surveys in Rebio Jaru.   #
# It is NOT a good idea to use this anywhere else because common names may mean completely #
# diferent things...                                                                       #
#==========================================================================================#
#==========================================================================================#
scientific.lookup.rja <<- function(datum,lookup.path){

   #----- Read in the look-up table. ------------------------------------------------------#
   lookup.file = paste(lookup.path,"rondonia_taxon_lookup.csv",sep="/")
   look.up = read.csv(file=lookup.file,stringsAsFactors=FALSE)
   look.up$common     = tolower(trim(look.up$common    ))
   look.up$scientific = trim(look.up$scientific)
   look.up$family     = trim(look.up$family    )
   #---------------------------------------------------------------------------------------#



   #----- Break into genus and species. ---------------------------------------------------#
   gs.list                  = sapply(X = tolower(look.up$scientific),FUN=strsplit,split=" ")
   gs.length                = sapply(X = gs.list, FUN = length)
   gs.mat                   = cbind( mapply(FUN="[",gs.list,MoreArgs=list(1))
                                   , mapply(FUN="[",gs.list,MoreArgs=list(2))
                                   )#end cbind
   g                        = capwords(gs.mat[,1],strict=TRUE)
   s                        = tolower(gs.mat[,2])
   g.s                      = paste(g,s,sep=" ")
   g.s[is.na(g) & is.na(s)] = NA
   look.up$scientific       = g.s
   look.up$genus            = g
   #---------------------------------------------------------------------------------------#



   #----- Trim the common names, and simplify/replace some names. -------------------------#
   datum$common                      = tolower(trim(datum$common))
   datum$common[is.na(datum$common)] = "mato"
   #---------------------------------------------------------------------------------------#



   #----- Find all unique common names. ---------------------------------------------------#
   unique.common = unique(datum$common)
   n.common      = length(unique.common)
   notfound      = NULL
   for (n in 1:n.common){
      #----- Find the trees that have the same common name as this one. -------------------#
      cat (" - ",n,"/",n.common," -- ",unique.common[n],"\n")
      w.dat  = which(datum$common %in% unique.common[n])
      n.dat  = length(w.dat)
      #------------------------------------------------------------------------------------#


      #----- Find the trees in the look-up table with the same common name. ---------------#
      w.look = which(look.up$common %in% unique.common[n])
      n.look = length(w.look)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Check how many trees have the same common name in the look-up table.          #
      #------------------------------------------------------------------------------------#
      if (n.look == 1){
         #----- Only one.  Use it. --------------------------------------------------------#
         datum$scientific   [w.dat] = look.up$scientific[w.look]
         datum$genus        [w.dat] = look.up$genus     [w.look]
         datum$gf.scientific[w.dat] = 1
      }else if (n.look > 1){
         datum$scientific   [w.dat] = sample( x       = look.up$scientific[w.look]
                                            , size    = n.dat
                                            , replace = TRUE
                                            )#end sample
         datum$genus        [w.dat] = look.up$genus     [w.look]
         datum$gf.scientific[w.dat] = 1
      }else{
         notfound = c(notfound,unique.common[n])
         datum$scientific   [w.dat] = "Ignotum"
         datum$genus        [w.dat] = "Ignotum"
         datum$gf.scientific[w.dat] = 0
      }#end if
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#
   return(datum)
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function returns two characters, the genus and the scientific name.             #
#------------------------------------------------------------------------------------------#
keep.gen.spe.only <<- function(x,out=c("both","genus","species")){
   out = match.arg(out)

   if (is.matrix(x) | is.array(x)){
      ans = apply(X=x,MARGIN=dim(x),FUN=keep.gen.spe.only,out=out)
   }else if (is.list(x)){
      ans = lapply(X=x,FUN=keep.gen.spe.only,out=out)
   }else if (is.data.frame(x)){
      ans = sapply(X=x,FUN=keep.gen.spe.only,out=out)
   }else if (length(x) > 1){
      ans = sapply(X=x,FUN=keep.gen.spe.only,out=out)
   }else{
   
      if (is.na(x) | x %in% "") x="Ignotum"
   
      #----- Remove stuff that is not genus. ----------------------------------------------#
      x    = sub(pattern="cf\\."  ,replacement=""           ,ignore.case=TRUE,x=x)
      x    = sub(pattern="cf\\ "  ,replacement=""           ,ignore.case=TRUE,x=x)
      x    = sub(pattern="Ind\\ " ,replacement="deleteme\\ ",ignore.case=TRUE,x=x)
      #------------------------------------------------------------------------------------#

      #----- Split words. -----------------------------------------------------------------#
      ans  = unlist(strsplit(unlist(c(tolower(x))),split=" "))
      nans = length(ans)
      #------------------------------------------------------------------------------------#


      #----- Decide whether to keep the genus, species, or both. --------------------------#
      gen = NA_character_
      spe = NA_character_
      if (nans > 0){
         gen = capwords(ans[1],strict=TRUE)
         if (nans > 1 && (! gen %in% "Deleteme")){
            spe = tolower(ans[2])
            if (spe %in% c("sp.","sp.1","ni","ind","deleteme")) spe = NA_character_
         }else{
            gen = NA_character_
         }#end if
      }#end if
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Select the proper output.                                                      #
      #------------------------------------------------------------------------------------#
      if (out %in% "genus"){
         ans = gen
      }else if (out %in% "species"){
         ans = spe
      }else if (out %in% "both"){
         ans = ifelse(is.na(spe),gen,paste(gen,spe,sep=" "))
      }#end if
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#
   return(ans)
}#end function keep.gen.spe.only
#==========================================================================================#
#==========================================================================================#
