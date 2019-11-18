#----- List of flags for undetermined species. --------------------------------------------#
unknown.wildcard     <<- c("aff","cf","deleteme","ind","indet","na","ni","sp","spp"
                          ,"spnov","unknown","unidentified"
                          ,paste0("sp"  ,sequence(99))
                          ,paste0("spb" ,sequence(99))
                          ,paste0("spp" ,sequence(99))
                          ,paste0("sp"  ,sequence(99),"cay-atdn")
                          ,paste0("sp"  ,sequence(99),"guyafor")
                          ,paste0("spfg",sequence(99),"-holst")
                          )#end c
unknown.common       <<- "mato"
unk.liana.common     <<- "cipo"
unknown.phylum       <<- "Ignotophyta"
unknown.order        <<- "Ignotales"
unknown.family       <<- "Ignotaceae"
unknown.genus        <<- "Ignotum"
unknown.epithet      <<- "indet"
unk.liana.order      <<- "Lianales"
unk.liana.family     <<- "Lianaceae"
unk.liana.genus      <<- "Liana"
unknown.scientific   <<- paste(unknown.genus,unknown.epithet)
unk.liana.scientific <<- paste(unk.liana.genus,unknown.epithet)
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


   #----- Remove underscore marks. --------------------------------------------------------#
   x   = gsub(pattern="_",replacement=" ",x=x)
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
   sel = (x %in% "?"                              ); x[sel] = NA_character_
   sel = (x %in% "abacaba"                        ); x[sel] = "bacaba"
   sel = (x %in% "abicuiba"                       ); x[sel] = "ucuuba"
   sel = (x %in% "abirana rosadinha"              ); x[sel] = "abiu rosadinho"
   sel = (x %in% "abil fura fura"                 ); x[sel] = "abiu-fura-fura"
   sel = (x %in% "abiui"                          ); x[sel] = "abiu"
   sel = (x %in% "abiu acariquara"                ); x[sel] = "abiu-acariquara"
   sel = (x %in% "abiu cramuri"                   ); x[sel] = "abiu-cramuri"
   sel = (x %in% "abiu cutite"                    ); x[sel] = "abiu-cutite"
   sel = (x %in% "abiu cutite folha verde"        ); x[sel] = "abiu-cutite"
   sel = (x %in% "abiu cutiti"                    ); x[sel] = "abiu-cutite"
   sel = (x %in% "abiu-cutiti"                    ); x[sel] = "abiu-cutite"
   sel = (x %in% "abiu guajara"                   ); x[sel] = "abiu-guajara"
   sel = (x %in% "abiu goiabao"                   ); x[sel] = "abiu-goiabao"
   sel = (x %in% "abiu mangabinha"                ); x[sel] = "abiu-mangabinha"
   sel = (x %in% "abiu tauari"                    ); x[sel] = "tauari"
   sel = (x %in% "abiu vermelha"                  ); x[sel] = "abiu vermelho"
   sel = (x %in% "abiurana vermelho"              ); x[sel] = "abiurana vermelha"
   sel = (x %in% "abiurana vermlho"               ); x[sel] = "abiurana vermelha"
   sel = (x %in% "abiuarana"                      ); x[sel] = "abiurana"
   sel = (x %in% "abiurana rosadinha"             ); x[sel] = "abiu rosadinho"
   sel = (x %in% "acoita cavalo"                  ); x[sel] = "acoita-cavalo"
   sel = (x %in% "algodao brabo"                  ); x[sel] = "algodao-bravo"
   sel = (x %in% "algodao bravo"                  ); x[sel] = "algodao-bravo"
   sel = (x %in% "amalelinha"                     ); x[sel] = "amarelinho"
   sel = (x %in% "amapa tirana"                   ); x[sel] = "amapatirana"
   sel = (x %in% "amarelinha"                     ); x[sel] = "amarelinho"
   sel = (x %in% "ameixa"                         ); x[sel] = "ameixa-do-para"
   sel = (x %in% "amerelinho"                     ); x[sel] = "amarelinho"
   sel = (x %in% "amescla"                        ); x[sel] = "breu mescla"
   sel = (x %in% "amoninha"                       ); x[sel] = "mamoninha"
   sel = (x %in% "ananin"                         ); x[sel] = "anani"
   sel = (x %in% "angelim aroreira"               ); x[sel] = "angelim-aroeira"
   sel = (x %in% "angelim aroeira"                ); x[sel] = "angelim-aroeira"
   sel = (x %in% "angelim margoso"                ); x[sel] = "angelim amargoso"
   sel = (x %in% "angelim pedra"                  ); x[sel] = "angelim-pedra"
   sel = (x %in% "angelim pedro"                  ); x[sel] = "angelim-pedra"
   sel = (x %in% "angelim peroba"                 ); x[sel] = "angelim-pedra"
   sel = (x %in% "apuii"                          ); x[sel] = "apui"
   sel = (x %in% "aquariquara"                    ); x[sel] = "acariquara" 
   sel = (x %in% "aquariquarana"                  ); x[sel] = "acariquarana"
   sel = (x %in% "araca nego"                     ); x[sel] = "araca"
   sel = (x %in% "araca de anta"                  ); x[sel] = "araca-de-anta"
   sel = (x %in% "araracanga marfim"              ); x[sel] = "araracanga-marfim"
   sel = (x %in% "aratacio"                       ); x[sel] = "arataciu"
   sel = (x %in% "araticu"                        ); x[sel] = "araticum"
   sel = (x %in% "ata-menju"                      ); x[sel] = "atameju"
   sel = (x %in% "bacaba de leque"                ); x[sel] = "bacaba-de-leque"
   sel = (x %in% "bacabinha"                      ); x[sel] = "bacabi"
   sel = (x %in% "bacabinha quina"                ); x[sel] = "bacabinha-quina"
   sel = (x %in% "bacupari"                       ); x[sel] = "bacuripari"
   sel = (x %in% "bacuri-pari"                    ); x[sel] = "bacuripari"
   sel = (x %in% "bacuri pari"                    ); x[sel] = "bacuripari"
   sel = (x %in% "barbatimao"                     ); x[sel] = "fava-barbatimao"
   sel = (x %in% "bate puta"                      ); x[sel] = "batiputa"
   sel = (x %in% "baubarana"                      ); x[sel] = "embaubarana"
   sel = (x %in% "breu/louro preto?"              ); x[sel] = NA_character_
   sel = (x %in% "breu aroeira"                   ); x[sel] = "breu-aroeira"
   sel = (x %in% "breu sucuruba"                  ); x[sel] = "breu-sucuruba"
   sel = (x %in% "breu sucuruba branco"           ); x[sel] = "breu-sucuruba branco"
   sel = (x %in% "breu sucurubinha"               ); x[sel] = "breu-sucurubinha"
   sel = (x %in% "babao"                          ); x[sel] = "macauba"
   sel = (x %in% "bolao"                          ); x[sel] = "fava-bolota"
   sel = (x %in% "bombeira"                       ); x[sel] = "pau-pombo"
   sel = (x %in% "brau"                           ); x[sel] = "breu"
   sel = (x %in% "brejauba"                       ); x[sel] = "brejauva"              
   sel = (x %in% "breu sucuuba"                   ); x[sel] = "breu sucuruba"         
   sel = (x %in% "cabeca de urubu"                ); x[sel] = "cabeca-de-urubu"       
   sel = (x %in% "cabela"                         ); x[sel] = "louro canela"          
   sel = (x %in% "cabriuva"                       ); x[sel] = "cabreuva"              
   sel = (x %in% "cabriuna"                       ); x[sel] = "cabreuva"              
   sel = (x %in% "cabreu"                         ); x[sel] = "cabreuva"              
   sel = (x %in% "cabreuba"                       ); x[sel] = "cabreuva"              
   sel = (x %in% "caca piolho"                    ); x[sel] = "mata-piolho"           
   sel = (x %in% "cacau bravo"                    ); x[sel] = "cacaui"                
   sel = (x %in% "cacau da mata"                  ); x[sel] = "cacau"  
   sel = (x %in% "cacaurana"                      ); x[sel] = "cacau"                 
   sel = (x %in% "cachixa"                        ); x[sel] = "caxixa"
   sel = (x %in% "cachudinha"                     ); x[sel] = "cascudinha"            
   sel = (x %in% "cagaca"                         ); x[sel] = "abiurana-cagaca"
   sel = (x %in% "cajarana"                       ); x[sel] = "jarana"                
   sel = (x %in% "caju acu"                       ); x[sel] = "cajuacu"
   sel = (x %in% "cajuba"                         ); x[sel] = "caju"
   sel = (x %in% "calcho"                         ); x[sel] = "caucho"                
   sel = (x %in% "canafistula"                    ); x[sel] = "fava-marimari"
   sel = (x %in% "canela brava"                   ); x[sel] = "catuaba"
   sel = (x %in% "canela de anta"                 ); x[sel] = "canela-de-anta"        
   sel = (x %in% "canela de veado"                ); x[sel] = "canela-de-veado"       
   sel = (x %in% "canela de velho"                ); x[sel] = "canela-de-velho"       
   sel = (x %in% "canelha velha"                  ); x[sel] = "canela-de-velho"       
   sel = (x %in% "canela de jacamim"              ); x[sel] = "canela-de-jacamim"     
   sel = (x %in% "canella de jacami"              ); x[sel] = "canela-de-jacamim"     
   sel = (x %in% "canella ge jacami"              ); x[sel] = "canela-de-jacamim"     
   sel = (x %in% "canniela"                       ); x[sel] = "canela"                
   sel = (x %in% "capa bode"                      ); x[sel] = "capa-bode"
   sel = (x %in% "capoeiro"                       ); x[sel] = "capueiro"
   sel = (x %in% "capoeiro preto"                 ); x[sel] = "capueiro preto"
   sel = (x %in% "capoeiro branco"                ); x[sel] = "capueiro"
   sel = (x %in% "captiurana"                     ); x[sel] = "capitiurana"
   sel = (x %in% "capueiro branco"                ); x[sel] = "capueiro"
   sel = (x %in% "caqui branco"                   ); x[sel] = "caqui"
   sel = (x %in% "caqui folha grande"             ); x[sel] = "caqui"
   sel = (x %in% "carobia"                        ); x[sel] = "caroba"                
   sel = (x %in% "cascudinho"                     ); x[sel] = "cascudinha"            
   sel = (x %in% "cascudo"                        ); x[sel] = "cascudinha"            
   sel = (x %in% "castanha do brasil"             ); x[sel] = "castanha-do-para"
   sel = (x %in% "castanha do para"               ); x[sel] = "castanha-do-para"
   sel = (x %in% "castanha"                       ); x[sel] = "castanha-do-para"
   sel = (x %in% "castanheiro"                    ); x[sel] = "castanha-do-para"
   sel = (x %in% "castanha de galinha"            ); x[sel] = "castanha-de-galinha"   
   sel = (x %in% "castanha de morcego"            ); x[sel] = "castanha-de-morcego"
   sel = (x %in% "castanha de periquito"          ); x[sel] = "castanha-de-periquito" 
   sel = (x %in% "castanha de sapocaia"           ); x[sel] = "castanha-sapucaia"     
   sel = (x %in% "castanha de sapucaia"           ); x[sel] = "castanha-sapucaia"     
   sel = (x %in% "castanha sapocaia"              ); x[sel] = "castanha-sapucaia"     
   sel = (x %in% "castanheira"                    ); x[sel] = "castanha-do-para"      
   sel = (x %in% "cauba"                          ); x[sel] = "macacauba"
   sel = (x %in% "cauxo"                          ); x[sel] = "caucho"                
   sel = (x %in% "caxeta"                         ); x[sel] = "caixeta"               
   sel = (x %in% "caxicha"                        ); x[sel] = "caxixa"           
   sel = (x %in% "caximbeira"                     ); x[sel] = "cachimbeiro"           
   sel = (x %in% "caximbeiro"                     ); x[sel] = "cachimbeiro"           
   sel = (x %in% "caxudinha"                      ); x[sel] = "cascudinha"            
   sel = (x %in% "cedroarana"                     ); x[sel] = "cedrorana"
   sel = (x %in% "cega jumenta"                   ); x[sel] = "cega-jumento"
   sel = (x %in% "cega machado"                   ); x[sel] = "cega-machado"
   sel = (x %in% "chamaecrista"                   ); x[sel] = "coracao-de-negro"
   sel = (x %in% "chapeu de sol"                  ); x[sel] = "chapeu-de-sol"
   sel = (x %in% "chocolate"                      ); x[sel] = "cacau"                 
   sel = (x %in% "cip"                            ); x[sel] = unk.liana.common
   sel = (x %in% "cipo(dbh a 0.9m do chao)"       ); x[sel] = unk.liana.common
   sel = (x %in% "cipo+arapo"                     ); x[sel] = unk.liana.common
   sel = (x %in% "cipo abuta"                     ); x[sel] = "cipo-abuta"
   sel = (x %in% "cipo cocoloba"                  ); x[sel] = "cipo-cocoloba"
   sel = (x %in% "cipo de fogo"                   ); x[sel] = "cipo-de-fogo"
   sel = (x %in% "cipo de fogo folha peluda"      ); x[sel] = "cipo-de-fogo folha peluda"
   sel = (x %in% "cipo estopinha"                 ); x[sel] = "cipo-estopinha"
   sel = (x %in% "cipo falso curare"              ); x[sel] = "cipo falso-curare"
   sel = (x %in% "cipo folha de canela"           ); x[sel] = "cipo folha-de-canela"
   sel = (x %in% "cipo leiteiro"                  ); x[sel] = "cipo-leiteiro"
   sel = (x %in% "cipo timbo"                     ); x[sel] = "cipo-timbo"
   sel = (x %in% "cipo timborana"                 ); x[sel] = "cipo-timborana"
   sel = (x %in% "coca brava"                     ); x[sel] = "coca-brava"
   sel = (x %in% "coco pau"                       ); x[sel] = "coco-pau"
   sel = (x %in% "comida de jaboti"               ); x[sel] = "comida-de-jabuti"      
   sel = (x %in% "comida de jabuti"               ); x[sel] = "comida-de-jabuti"      
   sel = (x %in% "comida de pomba"                ); x[sel] = "comida-de-pombo"       
   sel = (x %in% "comida de pombo"                ); x[sel] = "comida-de-pombo"       
   sel = (x %in% "coracao  de negro"              ); x[sel] = "coracao-de-negro"      
   sel = (x %in% "coracao de nego"                ); x[sel] = "coracao-de-negro"      
   sel = (x %in% "coracao de negro"               ); x[sel] = "coracao-de-negro"      
   sel = (x %in% "corante de indio"               ); x[sel] = "urucum"                
   sel = (x %in% "coro preto"                     ); x[sel] = "louro preto"           
   sel = (x %in% "coussarea racemosa"             ); x[sel] = "caferana"              
   sel = (x %in% "crista de mutum"                ); x[sel] = "crista"  
   sel = (x %in% "cucutiteriba"                   ); x[sel] = "cucutitiriba"
   sel = (x %in% "cumari"                         ); x[sel] = "cumaru"                
   sel = (x %in% "cumaru/apui"                    ); x[sel] = "apui"                  
   sel = (x %in% "cumaru de ferro"                ); x[sel] = "cumaru-ferro"          
   sel = (x %in% "cumatezinho/goiabinha fm"       ); x[sel] = "cumatezinho"           
   sel = (x %in% "cupiu"                          ); x[sel] = "cupui"
   sel = (x %in% "cupuacu da mata"                ); x[sel] = "cupuacu-da-mata"
   sel = (x %in% "cutiti"                         ); x[sel] = "abiu cutite"
   sel = (x %in% "cutiti"                         ); x[sel] = "abiu cutite"
   sel = (x %in% "cutite"                         ); x[sel] = "abiu cutite"
   sel = (x %in% "cutiteriba"                     ); x[sel] = "cucutitiriba"
   sel = (x %in% "cutitiriba"                     ); x[sel] = "cucutitiriba"
   sel = (x %in% "cruzeiro"                       ); x[sel] = "quina-cruzeiro"
   sel = (x %in% "dulacia/cachaceira"             ); x[sel] = "cachaceiro"
   sel = (x %in% "dulacia/cachaceira"             ); x[sel] = "cachaceiro"
   sel = (x %in% "einvira preta"                  ); x[sel] = "envira preta" 
   sel = (x %in% "embauba branco"                 ); x[sel] = "embauba branca"        
   sel = (x %in% "embauba vick"                   ); x[sel] = "embauba"               
   sel = (x %in% "embauba torem"                  ); x[sel] = "embauba toren"  
   sel = (x %in% "embirata"                       ); x[sel] = "envira ata"            
   sel = (x %in% "embireira branca"               ); x[sel] = "envira"                
   sel = (x %in% "embireira rosa"                 ); x[sel] = "envira"                
   sel = (x %in% "envira biriba"                  ); x[sel] = "envira-biriba"         
   sel = (x %in% "envira cana"                    ); x[sel] = "envira-cana"           
   sel = (x %in% "envira cabo de rodo"            ); x[sel] = "envira cabo-de-rodo"   
   sel = (x %in% "envira caju"                    ); x[sel] = "envira-caju"           
   sel = (x %in% "envira conduru"                 ); x[sel] = "envira-conduru"        
   sel = (x %in% "envira cunauaru"                ); x[sel] = "envira-cunauaru"       
   sel = (x %in% "envira de caju"                 ); x[sel] = "envira-jacu"           
   sel = (x %in% "envira de jacu"                 ); x[sel] = "envira-jacu"           
   sel = (x %in% "envira-mao-de-onca"             ); x[sel] = "envira mao-de-onca"    
   sel = (x %in% "envira molia branco"            ); x[sel] = "envira-molia branco"
   sel = (x %in% "envira preto"                   ); x[sel] = "envira preta"          
   sel = (x %in% "envira quiabo"                  ); x[sel] = "axixa"
   sel = (x %in% "envira-quiabo"                  ); x[sel] = "axixa"
   sel = (x %in% "envira sombrera"                ); x[sel] = "envira-sombreiro"
   sel = (x %in% "envira sombreira"               ); x[sel] = "envira-sombreiro"
   sel = (x %in% "envira sombreiro"               ); x[sel] = "envira-sombreiro"
   sel = (x %in% "envira surucu"                  ); x[sel] = "envira-surucucu"
   sel = (x %in% "envira surucucu"                ); x[sel] = "envira-surucucu"
   sel = (x %in% "envira taia"                    ); x[sel] = "envira-taia"           
   sel = (x %in% "envira turi"                    ); x[sel] = "envira-turi"
   sel = (x %in% "envira turi duro"               ); x[sel] = "envira-turi"
   sel = (x %in% "envira vassourinha"             ); x[sel] = "envira-vassourinha"             
   sel = (x %in% "envira vassourinha folha peluda"); x[sel] = "envira-vassourinha folha peluda"
   sel = (x %in% "envira vermelho"                ); x[sel] = "envira vermelha"       
   sel = (x %in% "envireira"                      ); x[sel] = "envira"       
   sel = (x %in% "escorrega macaco"               ); x[sel] = "escorrega-macaco"      
   sel = (x %in% "escurrega macaco"               ); x[sel] = "escorrega-macaco"      
   sel = (x %in% "espeturana f. g."               ); x[sel] = "espeturana"            
   sel = (x %in% "espinheira santa"               ); x[sel] = "espinheira-santa"      
   sel = (x %in% "fava amarg"                     ); x[sel] = "fava amargosa"
   sel = (x %in% "fava arara"                     ); x[sel] = "fava-arara-tucupi"
   sel = (x %in% "fava arara tucupi"              ); x[sel] = "fava-arara-tucupi"     
   sel = (x %in% "fava atana"                     ); x[sel] = "fava-atana"
   sel = (x %in% "fiora preta"                    ); x[sel] = NA_character_
   sel = (x %in% "fava barbatimao"                ); x[sel] = "fava-barbatimao"       
   sel = (x %in% "fava bolacha"                   ); x[sel] = "fava-bolacha"          
   sel = (x %in% "fava bolota"                    ); x[sel] = "fava-bolota"           
   sel = (x %in% "fava cana"                      ); x[sel] = "fava-cana"           
   sel = (x %in% "fava core"                      ); x[sel] = "fava-core"           
   sel = (x %in% "fava de anta"                   ); x[sel] = "fava-de-anta"          
   sel = (x %in% "fava marimari"                  ); x[sel] = "fava-marimari"         
   sel = (x %in% "fava mapuxique"                 ); x[sel] = "fava-mapuxiqui"        
   sel = (x %in% "fava mapuxiqui"                 ); x[sel] = "fava-mapuxiqui"        
   sel = (x %in% "fava orelha"                    ); x[sel] = "fava orelha-de-macaco" 
   sel = (x %in% "fava orelha de macaco"          ); x[sel] = "fava orelha-de-macaco" 
   sel = (x %in% "fava paricarana"                ); x[sel] = "fava-paricana"         
   sel = (x %in% "fava saboeira"                  ); x[sel] = "fava-saboeira"         
   sel = (x %in% "fava tambori"                   ); x[sel] = "fava-tamboril"         
   sel = (x %in% "fava tamburi"                   ); x[sel] = "fava-tamboril"         
   sel = (x %in% "fava tamboril"                  ); x[sel] = "fava-tamboril" 
   sel = (x %in% "fava tamboriu"                  ); x[sel] = "fava-tamboril"  
   sel = (x %in% "favera amargosa"                ); x[sel] = "fava amargosa"         
   sel = (x %in% "faveira branca"                 ); x[sel] = "fava branca"           
   sel = (x %in% "feijo branco"                   ); x[sel] = "freijo branco"         
   sel = (x %in% "ferdinandusa elliptica"         ); x[sel] = "bacabinha quina"       
   sel = (x %in% "figado de preguisa"             ); x[sel] = "figado-de-preguica"
   sel = (x %in% "figueira brava"                 ); x[sel] = "figueira"              
   sel = (x %in% "gameleiro"                      ); x[sel] = "gameleira"             
   sel = (x %in% "gapeba"                         ); x[sel] = "guapeva"               
   sel = (x %in% "guapeba"                        ); x[sel] = "guapeva"               
   sel = (x %in% "gema de ovo"                    ); x[sel] = "gema-de-ovo"           
   sel = (x %in% "geniparana"                     ); x[sel] = "jeniparana"            
   sel = (x %in% "genipapo"                       ); x[sel] = "jenipapo"              
   sel = (x %in% "goiaba"                         ); x[sel] = "araca"
   sel = (x %in% "goiaba de anta"                 ); x[sel] = "goiaba-de-anta"
   sel = (x %in% "goibarana"                      ); x[sel] = "goiabarana"            
   sel = (x %in% "gombeira fg"                    ); x[sel] = "gombeira folha grande"
   sel = (x %in% "gombeira vermelho"              ); x[sel] = "gombeira vermelha"     
   sel = (x %in% "grao de galo"                   ); x[sel] = "grao-de-galo"          
   sel = (x %in% "grao de guariba"                ); x[sel] = "grao-de-guariba"       
   sel = (x %in% "grao de macaco"                 ); x[sel] = "grao-de-guariba" 
   sel = (x %in% "gravilola brava"                ); x[sel] = "graviola-brava"        
   sel = (x %in% "graviola brava"                 ); x[sel] = "graviola-brava"        
   sel = (x %in% "guaiba"                         ); x[sel] = "araca"                 
   sel = (x %in% "guaiaba"                        ); x[sel] = "araca"
   sel = (x %in% "guajara bolacha"                ); x[sel] = "guajara-bolacha"       
   sel = (x %in% "guajara ferro"                  ); x[sel] = "guajara-ferro"         
   sel = (x %in% "guajara pedra"                  ); x[sel] = "guajara-pedra"         
   sel = (x %in% "guajara mirim"                  ); x[sel] = "guajara-mirim"         
   sel = (x %in% "guariuva"                       ); x[sel] = "guariuba"              
   sel = (x %in% "guaruba"                        ); x[sel] = "quaruba"
   sel = (x %in% "ibirucu"                        ); x[sel] = "embirucu"              
   sel = (x %in% "imbauba torem"                  ); x[sel] = "embauba toren"         
   sel = (x %in% "imbirata"                       ); x[sel] = "envira ata"            
   sel = (x %in% "imbireira"                      ); x[sel] = "envira"                
   sel = (x %in% "imbireira rosa"                 ); x[sel] = "envira"                
   sel = (x %in% "imbiricu"                       ); x[sel] = "embirucu"              
   sel = (x %in% "imbirucu"                       ); x[sel] = "embirucu"              
   sel = (x %in% "inga a"                         ); x[sel] = "inga"
   sel = (x %in% "inga amarela"                   ); x[sel] = "inga amarelo"
   sel = (x %in% "inga branca"                    ); x[sel] = "inga branco"           
   sel = (x %in% "inga chichica"                  ); x[sel] = "inga-xixica"           
   sel = (x %in% "inga de orelha"                 ); x[sel] = "inga-de-orelha"        
   sel = (x %in% "inga de preguica"               ); x[sel] = "inga-de-preguica"      
   sel = (x %in% "inga de rodo"                   ); x[sel] = "inga-de-rodo"          
   sel = (x %in% "inga de rosca"                  ); x[sel] = "inga-de-rosca"         
   sel = (x %in% "inga f g"                       ); x[sel] = "inga folha grauda"
   sel = (x %in% "inga f.p."                      ); x[sel] = "inga folha peluda"     
   sel = (x %in% "inga folha grande"              ); x[sel] = "inga folha grauda"
   sel = (x %in% "inga folhao"                    ); x[sel] = "inga folha grauda"
   sel = (x %in% "inga peluda"                    ); x[sel] = "inga peludo"         
   sel = (x %in% "inga roceiro"                   ); x[sel] = "inga-roceiro"         
   sel = (x %in% "inga titica"                    ); x[sel] = "inga-xixica"           
   sel = (x %in% "inga vermelha"                  ); x[sel] = "inga vermelho"         
   sel = (x %in% "inga xixica"                    ); x[sel] = "inga-xixica"           
   sel = (x %in% "ipe amerelo"                    ); x[sel] = "ipe amarelo"           
   sel = (x %in% "jaboticaba"                     ); x[sel] = "jabuticaba"            
   sel = (x %in% "jacare"                         ); x[sel] = "pau-jacare"
   sel = (x %in% "jacariuba"                      ); x[sel] = "jacareuba"
   sel = (x %in% "jambo"                          ); x[sel] = "jambo-do-mato"
   sel = (x %in% "jara"                           ); x[sel] = "jarana"
   sel = (x %in% "jaruma"                         ); x[sel] = "taruma"
   sel = (x %in% "jauari"                         ); x[sel] = "tauari"
   sel = (x %in% "jenita"                         ); x[sel] = "janita"                
   sel = (x %in% "jito"                           ); x[sel] = "gito"
   sel = (x %in% "joao mole"                      ); x[sel] = "joao-mole"             
   sel = (x %in% "joao moleza"                    ); x[sel] = "joao-moleza"           
   sel = (x %in% "jotobazinho"                    ); x[sel] = "jatobazinho"           
   sel = (x %in% "jutai mirim"                    ); x[sel] = "jutai-mirim"           
   sel = (x %in% "jutai acu"                      ); x[sel] = "jutai-acu"             
   sel = (x %in% "jutai pororoca"                 ); x[sel] = "jutai-pororoca"        
   sel = (x %in% "kaqui"                          ); x[sel] = "caqui"        
   sel = (x %in% "lacre da mata"                  ); x[sel] = "lacre-da-mata"         
   sel = (x %in% "laranginha"                     ); x[sel] = "laranjinha"            
   sel = (x %in% "leiteiro"                       ); x[sel] = "leiteira"              
   sel = (x %in% "leitera"                        ); x[sel] = "leiteira"              
   sel = (x %in% "loro amarelo"                   ); x[sel] = "louro amarelo"         
   sel = (x %in% "louro?"                         ); x[sel] = "louro"
   sel = (x %in% "louro abacate"                  ); x[sel] = "louro-abacate"         
   sel = (x %in% "louro aritu"                    ); x[sel] = "louro-aritu"           
   sel = (x %in% "louro branco"                   ); x[sel] = "louro"                 
   sel = (x %in% "louro bosta"                    ); x[sel] = "louro-bosta"           
   sel = (x %in% "louro canela"                   ); x[sel] = "louro canelado" 
   sel = (x %in% "louro chumbo"                   ); x[sel] = "louro-chumbo"          
   sel = (x %in% "louro faia"                     ); x[sel] = "louro-faia"            
   sel = (x %in% "louro p"                        ); x[sel] = "louro preto"
   sel = (x %in% "louro pimenta"                  ); x[sel] = "louro-pimenta"         
   sel = (x %in% "louro tamanco"                  ); x[sel] = "louro-tamanco"         
   sel = (x %in% "louro seda"                     ); x[sel] = "louro-seda"            
   sel = (x %in% "macucu de sangue"               ); x[sel] = "macucu-de-sangue"
   sel = (x %in% "mafim"                          ); x[sel] = "marfim"
   sel = (x %in% "mamao jacatia"                  ); x[sel] = "jacaratia"             
   sel = (x %in% "mamica de porca"                ); x[sel] = "mamica-de-porca"       
   sel = (x %in% "mamonini"                       ); x[sel] = "mamoninha"             
   sel = (x %in% "mandioqueiro"                   ); x[sel] = "mandioqueira"
   sel = (x %in% "mandioqueiro escamoso"          ); x[sel] = "mandioqueira"
   sel = (x %in% "mangueira"                      ); x[sel] = "manguerana"
   sel = (x %in% "manguerano"                     ); x[sel] = "manguerana"
   sel = (x %in% "maparajuba"                     ); x[sel] = "parajuba"
   sel = (x %in% "maprounea"                      ); x[sel] = "caxixa"
   sel = (x %in% "mapuxique"                      ); x[sel] = "fava-mapuxiqui"
   sel = (x %in% "mapuxiqui"                      ); x[sel] = "fava-mapuxiqui"
   sel = (x %in% "maquira"                        ); x[sel] = "muiratinga"
   sel = (x %in% "maracata"                       ); x[sel] = "marassacaca"
   sel = (x %in% "maracatia"                      ); x[sel] = "muiracatiara"
   sel = (x %in% "maracacaca"                     ); x[sel] = "marassacaca"
   sel = (x %in% "marasacaca"                     ); x[sel] = "marassacaca"
   sel = (x %in% "marimari"                       ); x[sel] = "fava-marimari"
   sel = (x %in% "massaranduba"                   ); x[sel] = "macaranduba"
   sel = (x %in% "mata fome"                      ); x[sel] = "mata-fome"
   sel = (x %in% "matamata branca"                ); x[sel] = "matamata branco"
   sel = (x %in% "matamata ci"                    ); x[sel] = "matamata-ci"
   sel = (x %in% "matamata cinza"                 ); x[sel] = "matamata-ci"
   sel = (x %in% "matamata jiboia"                ); x[sel] = "matamata-jiboia"
   sel = (x %in% "matamata vermelha"              ); x[sel] = "matamata vermelho"
   sel = (x %in% "mata pau+jito"                  ); x[sel] = "gito"
   sel = (x %in% "mata caldo"                     ); x[sel] = "mata-calado"
   sel = (x %in% "melanciera"                     ); x[sel] = "melancieira"
   sel = (x %in% "morto"                          ); x[sel] = "morta"
   sel = (x %in% "morango de macaco"              ); x[sel] = "morango-de-macaco"
   sel = (x %in% "morango de morcego"             ); x[sel] = "morango-de-macaco"
   sel = (x %in% "muiratinga folha grande/amapa"  ); x[sel] = "muiratinga folha grande"
   sel = (x %in% "muiratinga fura fura"           ); x[sel] = "muiratinga fura-fura"
   sel = (x %in% "mulatero"                       ); x[sel] = "mulateiro"
   sel = (x %in% "mulugu"                         ); x[sel] = "mulungu"
   sel = (x %in% "murta da mata"                  ); x[sel] = "murta-da-mata"
   sel = (x %in% "murta vassourinha"              ); x[sel] = "murta-vassourinha"
   sel = (x %in% "murucidu mata"                  ); x[sel] = "muruci-da-mata"
   sel = (x %in% "muruci da mata"                 ); x[sel] = "muruci-da-mata"
   sel = (x %in% "muruci fp"                      ); x[sel] = "muruci folha peluda"
   sel = (x %in% "mutama"                         ); x[sel] = "mutambo"
   sel = (x %in% "mutamba"                        ); x[sel] = "mutambo"
   sel = (x %in% "mututiassu"                     ); x[sel] = "mututi-acu"
   sel = (x %in% "ni"                             ); x[sel] = NA_character_
   sel = (x %in% "orelha de burro"                ); x[sel] = "orelha-de-burro"
   sel = (x %in% "orelha de macaco"               ); x[sel] = "fava orelha-de-macaco"
   sel = (x %in% "olho de sapo"                   ); x[sel] = "olho-de-sapo"
   sel = (x %in% "olho de veado"                  ); x[sel] = "olho-de-veado"
   sel = (x %in% "olho de viado"                  ); x[sel] = "olho-de-veado"
   sel = (x %in% "ouro branco"                    ); x[sel] = "seringueira"
   sel = (x %in% "p bolacha"                      ); x[sel] = "guajara-bolacha"
   sel = (x %in% "paineira"                       ); x[sel] = "sumauma"
   sel = (x %in% "palmito"                        ); x[sel] = "acai"
   sel = (x %in% "palmito babosa"                 ); x[sel] = "acai"
   sel = (x %in% "pao de sangue"                  ); x[sel] = "pau-sangue"
   sel = (x %in% "papo de mutum"                  ); x[sel] = "papo-de-mutum"
   sel = (x %in% "papo de mutum ff"               ); x[sel] = "papo-de-mutum folha fina"
   sel = (x %in% "papo de mutum  ff"              ); x[sel] = "papo-de-mutum folha fina"
   sel = (x %in% "papo-de-mutum ff"               ); x[sel] = "papo-de-mutum folha fina"
   sel = (x %in% "papo-de-mutum  ff"              ); x[sel] = "papo-de-mutum folha fina"
   sel = (x %in% "paricarana"                     ); x[sel] = "fava-paricana"
   sel = (x %in% "passarinhiera"                  ); x[sel] = "passarinheira"
   sel = (x %in% "pata de vaca"                   ); x[sel] = "pata-de-vaca"
   sel = (x %in% "pata preta"                     ); x[sel] = "ata preta"
   sel = (x %in% "patua"                          ); x[sel] = "pataua"
   sel = (x %in% "pau de arco"                    ); x[sel] = "pau-de-arco"
   sel = (x %in% "pau d.arco"                     ); x[sel] = "pau-de-arco"
   sel = (x %in% "pau de bicho"                   ); x[sel] = "pau-de-bicho"
   sel = (x %in% "pau de cobra"                   ); x[sel] = "pau-cobra"
   sel = (x %in% "pau colher"                     ); x[sel] = "pau-de-colher"
   sel = (x %in% "pau de colher"                  ); x[sel] = "pau-de-colher"
   sel = (x %in% "pau de jacare"                  ); x[sel] = "pau-jacare"
   sel = (x %in% "pau de macaco"                  ); x[sel] = "pau-de-macaco"
   sel = (x %in% "pau de rego"                    ); x[sel] = "pau-de-remo"
   sel = (x %in% "pau de remo"                    ); x[sel] = "pau-de-remo"
   sel = (x %in% "pau de sangue"                  ); x[sel] = "pau-sangue"
   sel = (x %in% "pau doce"                       ); x[sel] = "pau-doce"
   sel = (x %in% "pau jacare"                     ); x[sel] = "pau-jacare"
   sel = (x %in% "pau marfim"                     ); x[sel] = "pau-marfim"
   sel = (x %in% "pau mulato"                     ); x[sel] = "pau-mulato"
   sel = (x %in% "pau para tudo"                  ); x[sel] = "pau-para-tudo"
   sel = (x %in% "pau pereira"                    ); x[sel] = "peroba mica"
   sel = (x %in% "pau-pra-tudo"                   ); x[sel] = "pau-para-tudo"
   sel = (x %in% "pau pra tudo"                   ); x[sel] = "pau-para-tudo"
   sel = (x %in% "paupratudo"                     ); x[sel] = "pau-para-tudo"
   sel = (x %in% "pau prego"                      ); x[sel] = "pau-de-remo"
   sel = (x %in% "pau purui"                      ); x[sel] = "purui"
   sel = (x %in% "pau rego"                       ); x[sel] = "pau-de-remo"
   sel = (x %in% "pau sangue"                     ); x[sel] = "pau-sangue"
   sel = (x %in% "pereauna"                       ); x[sel] = "perebuna"
   sel = (x %in% "pedra umi"                      ); x[sel] = "pedra ume-caa"
   sel = (x %in% "pelo de cutia"                  ); x[sel] = "pelo-de-cutia"
   sel = (x %in% "pente de macaco"                ); x[sel] = "pente-de-macaco"
   sel = (x %in% "pepino da mata"                 ); x[sel] = "pepino-do-mato"
   sel = (x %in% "pepino-da-mata"                 ); x[sel] = "pepino-do-mato"
   sel = (x %in% "pepino do mato"                 ); x[sel] = "pepino-do-mato"
   sel = (x %in% "pepino-do-mato"                 ); x[sel] = "pepino-do-mato"
   sel = (x %in% "perna de moca"                  ); x[sel] = "perna-de-moca"
   sel = (x %in% "piqui"                          ); x[sel] = "piquia"  
   sel = (x %in% "piqui rosa"                     ); x[sel] = "piquia"
   sel = (x %in% "piquiazeiro"                    ); x[sel] = "piquia"
   sel = (x %in% "pitomba da mata"                ); x[sel] = "pitomba-da-mata"
   sel = (x %in% "pororoca"                       ); x[sel] = "jutai-pororoca"
   sel = (x %in% "prapara"                        ); x[sel] = "parapara"
   sel = (x %in% "pratudo"                        ); x[sel] = "pau-para-tudo"
   sel = (x %in% "puruirana/purui branco"         ); x[sel] = "purui branco"
   sel = (x %in% "quaiquara"                      ); x[sel] = "acariquara"
   sel = (x %in% "quariquara"                     ); x[sel] = "acariquara"
   sel = (x %in% "quariquarana"                   ); x[sel] = "acariquara"
   sel = (x %in% "quariquari"                     ); x[sel] = "acariquara"            
   sel = (x %in% "quari quari"                    ); x[sel] = "acariquara"            
   sel = (x %in% "quaruba cedro"                  ); x[sel] = "quaruba-cedro"
   sel = (x %in% "quebrado"                       ); x[sel] = NA_character_
   sel = (x %in% "quina"                          ); x[sel] = "quinarana" 
   sel = (x %in% "quina cruzeiro"                 ); x[sel] = "quina-cruzeiro"        
   sel = (x %in% "rim de paca"                    ); x[sel] = "rim-de-paca"           
   sel = (x %in% "ripeiro"                        ); x[sel] = "ripeira"               
   sel = (x %in% "roxao"                          ); x[sel] = "roxinho"               
   sel = (x %in% "roxinao"                        ); x[sel] = "roxinho"               
   sel = (x %in% "sangra de agua"                 ); x[sel] = "sangra-de-agua"
   sel = (x %in% "sapucaia"                       ); x[sel] = "castanha-sapucaia"
   sel = (x %in% "saboeira"                       ); x[sel] = "fava-saboeiro"
   sel = (x %in% "saboeira amarela"               ); x[sel] = "fava-saboeiro amarela"
   sel = (x %in% "saboeiro"                       ); x[sel] = "fava-saboeiro"
   sel = (x %in% "saboiera"                       ); x[sel] = "fava-saboeiro"
   sel = (x %in% "sabueira"                       ); x[sel] = "fava-saboeiro"  
   sel = (x %in% "sabueiro"                       ); x[sel] = "fava-saboeiro"
   sel = (x %in% "sabuguero"                      ); x[sel] = "sabugueiro"
   sel = (x %in% "sajinera"                       ); x[sel] = NA_character_
   sel = (x %in% "samauma de terra firme"         ); x[sel] = "sumauma da terra firme" 
   sel = (x %in% "sangra d`agua"                  ); x[sel] = "sangra-de-agua" 
   sel = (x %in% "sardinheiro"                    ); x[sel] = "sardinheira"
   sel = (x %in% "segador"                        ); x[sel] = "cegador"               
   sel = (x %in% "seringa"                        ); x[sel] = "seringueira"
   sel = (x %in% "seringa branca"                 ); x[sel] = "seringueira"           
   sel = (x %in% "seringa branco"                 ); x[sel] = "seringueira"           
   sel = (x %in% "seringa verdadeira"             ); x[sel] = "seringueira"           
   sel = (x %in% "seringarana preta"              ); x[sel] = "seringarana"           
   sel = (x %in% "seritinga"                      ); x[sel] = "seringueira"           
   sel = (x %in% "sorveira"                       ); x[sel] = "sorva"
   sel = (x %in% "sorveira leite"                 ); x[sel] = "sorva"                 
   sel = (x %in% "sorvo"                          ); x[sel] = "sorva"
   sel = (x %in% "sova"                           ); x[sel] = "sorva"
   sel = (x %in% "sucuba"                         ); x[sel] = "sucuuba" 
   sel = (x %in% "sucupira p sapo"                ); x[sel] = "sucupira pele-de-sapo" 
   sel = (x %in% "sucupira pele de sapo"          ); x[sel] = "sucupira pele-de-sapo" 
   sel = (x %in% "sucuuba preta"                  ); x[sel] = "sucuuba" 
   sel = (x %in% "tachi branca"                   ); x[sel] = "tachi branco"          
   sel = (x %in% "tachi do csmpo"                 ); x[sel] = "tachi-do-campo"        
   sel = (x %in% "tachi preta"                    ); x[sel] = "tachi preto"           
   sel = (x %in% "tachi preto ???"                ); x[sel] = "tachi preto"
   sel = (x %in% "tachi preto folh"               ); x[sel] = "tachi preto"
   sel = (x %in% "tachi vermelha"                 ); x[sel] = "tachi vermelho"        
   sel = (x %in% "talquari"                       ); x[sel] = "tauari"                
   sel = (x %in% "tamaquarao"                     ); x[sel] = "tamaquare"             
   sel = (x %in% "tamarindu"                      ); x[sel] = "tamarindo"              
   sel = (x %in% "tamauma"                        ); x[sel] = "sumauma"             
   sel = (x %in% "tamboril"                       ); x[sel] = "fava-tamboril"
   sel = (x %in% "tamboriul"                      ); x[sel] = "fava-tamboril"
   sel = (x %in% "tanari roxo"                    ); x[sel] = "tauari"
   sel = (x %in% "tangarana"                      ); x[sel] = "tangirana"
   sel = (x %in% "tanimbuca"                      ); x[sel] = "tanibuca"
   sel = (x %in% "tapiririca"                     ); x[sel] = "tatapiririca"
   sel = (x %in% "tatapirirca"                    ); x[sel] = "tatapiririca"
   sel = (x %in% "tatapiririca verm."             ); x[sel] = "tatapiririca vermelha"
   sel = (x %in% "taturana"                       ); x[sel] = "taturuba"
   sel = (x %in% "tauri"                          ); x[sel] = "tauari"                
   sel = (x %in% "tento folha"                    ); x[sel] = "tento"
   sel = (x %in% "tento foha grauda"              ); x[sel] = "tento folha grauda"    
   sel = (x %in% "tintero"                        ); x[sel] = "tinteiro"              
   sel = (x %in% "titiriba"                       ); x[sel] = "cucutitiriba"              
   sel = (x %in% "tucuma acu"                     ); x[sel] = "tucuma-acu"
   sel = (x %in% "ucuarana"                       ); x[sel] = "urucurana"             
   sel = (x %in% "ucuuba da varzea"               ); x[sel] = "ucuuba-da-varzea"    
   sel = (x %in% "ucuuba da terra firme"          ); x[sel] = "ucuuba terra-firme"    
   sel = (x %in% "ucuuba terra firme"             ); x[sel] = "ucuuba terra-firme"    
   sel = (x %in% "ucuuba tf"                      ); x[sel] = "ucuuba terra-firme"    
   sel = (x %in% "ucuuba vermelho"                ); x[sel] = "ucuuba vermelha"       
   sel = (x %in% "umbia"                          ); x[sel] = "goiabarana"
   sel = (x %in% "unha de vaca"                   ); x[sel] = "pata-de-vaca"          
   sel = (x %in% "uruci"                          ); x[sel] = "muruci"
   sel = (x %in% "urucu"                          ); x[sel] = "urucum"
   sel = (x %in% "urucuri"                        ); x[sel] = "urucum"
   sel = (x %in% "uruucurana"                     ); x[sel] = "urucurana"             
   sel = (x %in% "virola"                         ); x[sel] = "ucuuba"
   sel = (x %in% "xaonoquito"                     ); x[sel] = "pau vermelho"
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
standard.scientific.name <<- function(dat,nafill=TRUE){
   #---------------------------------------------------------------------------------------#
   #     First we make sure all scientific names and families have the genus and family    #
   # capitalised (e.g.  scientific name: Araucaria angustifolia; family: Araucariaceae.    #
   #---------------------------------------------------------------------------------------#
   nplants = nrow(dat)
   #----- Make it case-insensitive. -------------------------------------------------------#
   dat$scientific = tolower(dat$scientific)
   #----- Break into genus and species. ---------------------------------------------------#
   g.s     = keep.gen.spe.only(dat$scientific)
   #---------------------------------------------------------------------------------------#

   #----- Remove the indetermined epithet for cleaning up typos. --------------------------#
   g.s     = gsub( pattern     = paste0("\\ ",unknown.epithet,"$")
                 , replacement = ""
                 , x           = g.s
                 )#end gsub
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #      Substitutions of scientific name.  We standardise these first so family becomes  #
   # easier.                                                                               #
   #---------------------------------------------------------------------------------------#
   synonym             = read.csv( file             = file.path(srcdir,"synonym_taxon.csv")
                                 , header           = TRUE
                                 , stringsAsFactors = FALSE
                                 )#end read.csv
   nsynonym            = nrow(synonym)
   idx                 = match(g.s,synonym$synonym)
   sel                 = is.finite(idx)
   g.s[sel]            = synonym$accepted[idx[sel]]
   #---------------------------------------------------------------------------------------#


   #----- Break again into genus and species. ---------------------------------------------#
   g.s = keep.gen.spe.only(x=g.s,is.debug=TRUE)
   #---------------------------------------------------------------------------------------#





   #----- Convert the g2f list into a data frame. -----------------------------------------#
   f2nog = read.csv( file             = file.path(srcdir,"family_list.csv")
                   , header           = TRUE
                   , stringsAsFactors = FALSE
                   )#end read.csv
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Replace family from genus with "unknown".  A few families have a single genus, so #
   # in this case use the genus.                                                           #
   #---------------------------------------------------------------------------------------#
   idx      = match(g.s,f2nog$family)
   fam      = ! is.na(idx)
   g.s[fam] = f2nog$genus[idx[fam]]
   #---------------------------------------------------------------------------------------#


   #----- Remove family names from scientific columns. ------------------------------------#
   bye                 = g.s %in% c("Ind","Unknown")
   g.s[bye]            = unknown.scientific
   liana               = g.s %in% c("Liane")
   g.s[liana]          = unk.liana.scientific
   #---------------------------------------------------------------------------------------#


   #------ Find genus. --------------------------------------------------------------------#
   g   = keep.gen.spe.only(x=g.s,out="genus")
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #      Copy the data back to the structure.                                             #
   #---------------------------------------------------------------------------------------#
   if (nafill){
      dat$scientific = g.s
      dat$genus      = g
   }else{
      dat$scientific = ifelse( test = g.s %in% unknown.scientific
                             , yes  = NA_character_,no=g.s
                             , no   = grepl( pattern     = paste0(" ",unknown.common)
                                           , replacement = ""
                                           , x           = g.s
                                           )#end grepl
                             )#end ifelse
      dat$genus      = ifelse(test=g.s %in% unknown.genus     ,yes=NA_character_,no=g  )
   }#end if (!nafill)
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
standard.family.name <<- function(datum,nafill=TRUE){
   #----- Make sure families are properly capitalised. ------------------------------------#
   if ("family" %in% datum){
      datum$family = capwords(datum$family,strict=TRUE)
   }else{
      datum$family = rep(unknown.family,times=nrow(datum))
   }#end if ("family" %in% datum)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Substitute obsolete/misspelt families that are straightforward.                   #
   #---------------------------------------------------------------------------------------#
   datum$family = sub("Bombacaceae"                 ,"Malvaceae"     ,x=datum$family)
   datum$family = sub("Caesalpinaceae"              ,"Fabaceae"      ,x=datum$family)
   datum$family = sub("Cecropiaceae"                ,"Urticaceae"    ,x=datum$family)
   datum$family = sub("Compositae"                  ,"Asteraceae"    ,x=datum$family)
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
   datum$family = sub("Myrsinaceae"                 ,"Primulaceae"   ,x=datum$family)
   datum$family = sub("Papilionaceae"               ,"Fabaceae"      ,x=datum$family)
   datum$family = sub("Tiliaceae"                   ,"Malvaceae"     ,x=datum$family)
   datum$family = sub("Sterculiaceae"               ,"Malvaceae"     ,x=datum$family)
   #---------------------------------------------------------------------------------------#



   #----- Read in the family look-up tables (including the dummy genera). -----------------#
   g2f = read.csv( file             = file.path(srcdir,"genus_family.csv")
                 , header           = TRUE
                 , stringsAsFactors = FALSE
                 )#end read.csv
   f2nog = read.csv( file             = file.path(srcdir,"family_list.csv")
                   , header           = TRUE
                   , stringsAsFactors = FALSE
                   )#end read.csv
   #---------------------------------------------------------------------------------------#


   #----- Append the dummy names to the family list. --------------------------------------#
   f2nog = f2nog[,c("genus","family")]
   g2f   = rbind(g2f,f2nog)
   uniq  = ! duplicated(g2f$genus)
   g2f   = g2f[uniq,,drop=FALSE]
   o     = order(g2f$genus)
   g2f   = g2f[o,,drop=FALSE]
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Check whether the family has all genera needed.                                   #
   #---------------------------------------------------------------------------------------#
   idx = which( ( ! (is.na(datum$genus) | (datum$genus %in% unknown.genus) ) ) 
              & ( ! datum$genus %in% g2f$genus ) 
              )#end which
   if (length(idx) > 0){
     tofill.gen = t(t(sort(unique(datum$genus[idx]))))
     vars.show  = c("trans","tag","x","y","scientific","genus","family")
     toshow     = which(names(datum) %in% vars.show)
     toshow     = datum[idx,toshow]
     toshow     = toshow[order(toshow[,"genus"]),]
     toshow     = toshow[! duplicated(toshow[,"genus"]),]
     cat0("-----------------------------------------------------------")
     cat0(" You must add genera to genus_family.csv:")
     cat0(" ")
     print(toshow,quote=FALSE)
     cat0("-----------------------------------------------------------")
     cat0(" ")
     browser()
   }#end if
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Check whether there are families missing in f2nog look up table.                  #
   #---------------------------------------------------------------------------------------#
   idx = which(! g2f$family %in% f2nog$family)
   if (length(idx) > 0){
     cat0(" You must add families to family_list.csv.")
     family.tofill = t(t(sort(unique(g2f$family[idx]))))
     browser()
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Fill in non-informative genus for individual with known family but unknown genus. #
   #---------------------------------------------------------------------------------------#
   yes.gen = ! ( is.na(datum$genus ) | (datum$genus %in% unknown.genus  ) )
   yes.fam = ! ( is.na(datum$family) | (datum$family %in% unknown.family) )
   nothing = ! yes.gen & ! yes.fam
   no.gen  = ! yes.gen &   yes.fam
   #------ No information whatsoever.  Fill in with non-informative taxa. -----------------#
   datum$genus      [nothing] = unknown.genus
   datum$scientific [nothing] = unknown.scientific
   datum$family     [nothing] = unknown.family
   #------ No genus, family only. ---------------------------------------------------------#
   idx                        = match(datum$family[no.gen],f2nog$family)
   datum$genus      [no.gen ] = f2nog$genus[idx]
   datum$scientific [no.gen ] = paste(datum$genus[no.gen],unknown.epithet,sep=" ")
   #------ Genus is known. ----------------------------------------------------------------#
   idx                        = match(datum$genus[yes.gen],g2f$genus)
   datum$family     [yes.gen] = g2f$family [idx]
   #---------------------------------------------------------------------------------------#

   #----- Replace placeholders with NA in case the user doesn't want dummy names. ---------#
   if (! nafill){
      datum$family = ifelse( test = datum$family %in% unknown.family
                           , yes  = NA_character_
                           , no   = datum$family
                           )#end ifelse
   }#end (! nafill)
   #---------------------------------------------------------------------------------------#

   return(datum)
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      This is the most up-to-date table to find orders and phyla given the family.  This  #
# contains.  In case the family is not found in the data base, we add non-informative      #
# values.                                                                                  #
#------------------------------------------------------------------------------------------#
standard.order.phylum.name <<- function(datum,nafill=TRUE){
   #----- Make sure orders and phyla are properly capitalised. ----------------------------#
   if ("order" %in% datum){
      datum$order = capwords(datum$order,strict=TRUE)
   }else{
      datum$order = rep(unknown.order,times=nrow(datum))
   }#end if ("order" %in% datum)
   if ("phylum" %in% datum){
      datum$phylum = capwords(datum$phylum,strict=TRUE)
   }else{
      datum$phylum = rep(unknown.phylum,times=nrow(datum))
   }#end if ("phylum" %in% datum)
   #---------------------------------------------------------------------------------------#



   #----- Convert the g2f list into a data frame. -----------------------------------------#
   f2o = read.csv( file             = file.path(srcdir,"family_order+phylum.csv")
                 , header           = TRUE
                 , stringsAsFactors = FALSE
                 )#end read.csv
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Check whether the family has all genera needed.                                   #
   #---------------------------------------------------------------------------------------#
   idx = which( ( ! (is.na(datum$family) | (datum$family %in% unknown.family) ) ) 
              & ( ! datum$family %in% f2o$family )
              & ( regexpr(pattern=unknown.family,text=datum$family,ignore.case=T) %==% -1)
              )#end which
   if (length(idx) > 0){
     tofill.fam = t(t(sort(unique(datum$family[idx]))))
     vars.show  = c("trans","tag","x","y","scientific","genus","family","order","phylum")
     toshow     = which(names(datum) %in% vars.show)
     toshow     = datum[idx,toshow]
     toshow     = toshow[order(toshow[,"genus"]),]
     toshow     = toshow[! duplicated(toshow[,"genus"]),]
     cat0("-----------------------------------------------------------")
     cat0(" You must add family to genus_order+phylum.csv:")
     cat0(" ")
     print(toshow,quote=FALSE)
     cat0("-----------------------------------------------------------")
     cat0(" ")
     browser()
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Fill in orders and phyla for individual with known families.                      #
   #---------------------------------------------------------------------------------------#
   tofill               = ( ( ! is.na(datum$family) )
                          & ( ! grepl(pattern=unknown.family,x=datum$family) )
                          )#end tofill
   idx                  = match(datum$family[tofill],f2o$family)
   datum$order [tofill] = f2o$order [idx]
   datum$phylum[tofill] = f2o$phylum[idx]
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
   lookup.file        = file.path(lookup.path,"SL_taxon_lookup.csv")
   look.up            = read.csv(file=lookup.file,stringsAsFactors=FALSE)
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
   g.s[is.na(g) & is.na(s)] = NA_character_
   look.up$scientific       = g.s
   look.up$genus            = g
   #---------------------------------------------------------------------------------------#
   
   
   
   #----- Trim the common names, and simplify/replace some names. -------------------------#
   datum$common                      = tolower(trim(datum$common))
   datum$common[is.na(datum$common)] = unknown.common
   #---------------------------------------------------------------------------------------#
        
        
       
   #----- Find all unique common names. ---------------------------------------------------#
   unique.common = unique(datum$common)
   n.common      = length(unique.common)
   notfound      = NULL
   for (n in sequence(n.common)){
      #----- Find the trees that have the same common name as this one. -------------------#
      cat0(" - ",n,"/",n.common," -- ",unique.common[n])
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
         datum$scientific   [w.dat] = unknown.scientific
         datum$genus        [w.dat] = unknown.genus
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
#      This attributes scientific names based on common names for surveys carried out at   #
# Barro Colorado Island.                                                                   #
#------------------------------------------------------------------------------------------#
scientific.lookup.BCI <<- function(datum,lookup.path){
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
   lookup.file        = file.path(lookup.path,"BCI_taxon_lookup.csv")
   look.up            = read.csv(file=lookup.file,stringsAsFactors=FALSE)
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
   g.s[is.na(g) & is.na(s)] = NA_character_
   look.up$scientific       = g.s
   look.up$genus            = g
   #---------------------------------------------------------------------------------------#
   
   
   
   #----- Trim the common names, and simplify/replace some names. -------------------------#
   datum$common                      = tolower(trim(datum$common))
   datum$common[is.na(datum$common)] = unknown.common
   #---------------------------------------------------------------------------------------#
        
        
       
   #----- Find all unique common names. ---------------------------------------------------#
   unique.common = unique(datum$common)
   n.common      = length(unique.common)
   notfound      = NULL
   for (n in sequence(n.common)){
      #----- Find the trees that have the same common name as this one. -------------------#
      cat0(" - ",n,"/",n.common," -- ",unique.common[n])
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
         datum$scientific   [w.dat] = unknown.scientific
         datum$genus        [w.dat] = unknown.genus
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
#      This attributes scientific names based on common names for surveys carried out at   #
# La Selva.                                                                                #
#------------------------------------------------------------------------------------------#
scientific.lookup.LSE <<- function(datum,lookup.path){
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
   lookup.file        = file.path(lookup.path,"LSE_taxon_lookup.csv")
   look.up            = read.csv(file=lookup.file,stringsAsFactors=FALSE)
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
   g.s[is.na(g) & is.na(s)] = NA_character_
   look.up$scientific       = g.s
   look.up$genus            = g
   #---------------------------------------------------------------------------------------#
   
   
   
   #----- Trim the common names, and simplify/replace some names. -------------------------#
   datum$common                      = tolower(trim(datum$common))
   datum$common[is.na(datum$common)] = unknown.common
   #---------------------------------------------------------------------------------------#
        
        
       
   #----- Find all unique common names. ---------------------------------------------------#
   unique.common = unique(datum$common)
   n.common      = length(unique.common)
   notfound      = NULL
   for (n in sequence(n.common)){
      #----- Find the trees that have the same common name as this one. -------------------#
      cat0(" - ",n,"/",n.common," -- ",unique.common[n])
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
         datum$scientific   [w.dat] = unknown.scientific
         datum$genus        [w.dat] = unknown.genus
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
#      This attributes scientific names based on common names for surveys carried out at   #
# Tanguro (not Sustainable Landscapes).                                                    #
#------------------------------------------------------------------------------------------#
scientific.lookup.TAN <<- function(datum,lookup.path){
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
   lookup.file        = file.path(lookup.path,"TAN_taxon_lookup.csv")
   look.up            = read.csv(file=lookup.file,stringsAsFactors=FALSE)
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
   g.s[is.na(g) & is.na(s)] = NA_character_
   look.up$scientific       = g.s
   look.up$genus            = g
   #---------------------------------------------------------------------------------------#
   
   
   
   #----- Trim the common names, and simplify/replace some names. -------------------------#
   datum$common                      = tolower(trim(datum$common))
   datum$common[is.na(datum$common)] = unknown.common
   #---------------------------------------------------------------------------------------#
        
        
       
   #----- Find all unique common names. ---------------------------------------------------#
   unique.common = unique(datum$common)
   n.common      = length(unique.common)
   notfound      = NULL
   for (n in sequence(n.common)){
      #----- Find the trees that have the same common name as this one. -------------------#
      cat0(" - ",n,"/",n.common," -- ",unique.common[n])
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
         datum$scientific   [w.dat] = unknown.scientific
         datum$genus        [w.dat] = unknown.genus
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
scientific.lookup.TNF <<- function(datum,lookup.path){
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
   lookup.file        = file.path(lookup.path,"TNF_taxon_lookup.csv")
   look.up            = read.csv(file=lookup.file,stringsAsFactors=FALSE)
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
   g.s[is.na(g) & is.na(s)] = NA_character_
   look.up$scientific       = g.s
   look.up$genus            = g
   #---------------------------------------------------------------------------------------#
   
   
   
   #----- Trim the common names, and simplify/replace some names. -------------------------#
   datum$common                      = tolower(trim(datum$common))
   datum$common[is.na(datum$common)] = unknown.common
   #---------------------------------------------------------------------------------------#
        
        
       
   #----- Find all unique common names. ---------------------------------------------------#
   unique.common = unique(datum$common)
   n.common      = length(unique.common)
   notfound      = NULL
   for (n in sequence(n.common)){
      #----- Find the trees that have the same common name as this one. -------------------#
      cat0(" - ",n,"/",n.common," -- ",unique.common[n])
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
         datum$scientific   [w.dat] = unknown.scientific
         datum$genus        [w.dat] = unknown.genus
         datum$gf.scientific[w.dat] = 0
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#
   return(datum)
}#end function scientific.lookup.TNF
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
      #----- Remove trailing spaces. ------------------------------------------------------#
      wood$region = trim(wood$region)
      #------------------------------------------------------------------------------------#

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
      if (! is.na(species[s]) && length(grep(unknown.genus,species[s])) == 0){
         if (verbose) cat0("   - ",s,"/",nspecies,"  -  ",species[s],".") 

         #----- Get the genus and family of this species, in case we need it. -------------#
         igen        = which (datum$scientific %in% species[s])
         this.genus  = unique(datum$genus[igen])
         this.family = unique(capwords(datum$family[igen],strict=TRUE))
         if (length(this.genus) != 1 || length(this.family) != 1){
            cat0(" - Perhaps a bastard genus?")
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
               cat0("Weird, length(iwood)=",length(iwood),"!")
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
                     cat0("     * Species",species[s]
                        ,"not found.  Use genus average instead.")
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
      if (! (families[f] %in% unknown.family)){
         if (verbose) cat0("   - ",f,"/",nfamilies,"  -  ",families[f],".")

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
         cat0(" Found families with unidentified genera!")
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
         cat0(" The following families are filled with global sampling: ")
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
#  - trait       -- trait to fill.                                                         #
#  - tdb.csv     -- csv file containing the trait data base.                               #
#  - country     -- restrict attribution to these countries (NULL uses all countries).     #
#  - continents  -- restrict attribution to these continents (NULL uses all continents).   #
#  - fsample     -- use random sampling to fill unidentified individuals, or individuals   #
#                   with genus that is not available at the trait data base?               #
#                   If FALSE then it uses averages                                         #
#  - weight      -- A weighting factor to give either probability or to weight             #
#                   the average.  This could be a vector with weights, or a character with #
#                   the name of the variable in datum to use as the weight, or an integer  #
#                   with the column to be used as the weighting factor.                    #
#  - verbose     -- Flag to control the amount of information                              #
#------------------------------------------------------------------------------------------#
find.trait <<- function( datum
                       , trait     = c("wood.dens","SLA","leaf.turnover"
                                      ,"c2n.leaf","c2p.leaf")
                       , tdb.csv   = file.path(srcdir,"TRY_species_table.csv")
                       , country   = NULL
                       , continent = NULL
                       , fsample   = TRUE
                       , weight    = NULL
                       , verbose   = FALSE
                       ){

   #----- Check that trait is valid. ------------------------------------------------------#
   trait    = match.arg(trait)
   n.trait  = paste0("n." ,trait)
   gf.trait = paste0("gf.",trait)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Initialise gap filling flag and region in case they are not there.                #
   #---------------------------------------------------------------------------------------#
   if (! trait %in% names(datum)){
      datum[[trait]]    = rep(NA_real_  ,times=nrow(datum))
   }#end (! trait %in% names(datum))
   if (! gf.trait %in% names(datum)){
      datum[[gf.trait]] = rep(NA_integer_,times=nrow(datum))
   }#end (! gf.trait %in% names(datum))
   if (! "country" %in% names(datum)){
      datum$country = rep(NA_character_,times=nrow(datum))
   }#end (! "country" %in% names(datum))
   if (! "continent" %in% names(datum)){
      datum$continent = rep(NA_character_,times=nrow(datum))
   }#end (! "continent" %in% names(datum))
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
   #     Read data base to fill traits.                                                    #
   #---------------------------------------------------------------------------------------#
   tdb  = read.csv(file=tdb.csv,header=TRUE,stringsAsFactors=FALSE)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Restrict data base to valid traits in selected regions, and make sure the data   #
   # base can be used.                                                                     #
   #---------------------------------------------------------------------------------------#
   if (! (trait %in% names(tdb))){
      stop(paste0(" Trait ",trait," is not available in trait file \""
                 ,basename(tdb.csv),"\"."))
   }else if (! all(c("country","continent") %in% names(tdb))){
      stop(paste0(" \"country\" and/or \"continent\" missing in trait file \""
                 ,basename(tdb.csv),"\"."))
   }else{
      #----- Restrict species to regions of interest. -------------------------------------#
      if (is.null(country)){
         sel.country = rep(x=TRUE,times=nrow(tdb))
      }else{
         sel.country = tdb$country %in% country
      }#end if (is.null(country))
      if (is.null(continent)){
         sel.continent = rep(x=TRUE,times=nrow(tdb))
      }else{
         sel.continent = tdb$continent %in% continent
      }#end if (is.null(continent))
      keep = is.finite(tdb[[trait]]) & sel.country & sel.continent
      tdb  = tdb[keep,,drop=FALSE]
      ntdb = nrow(tdb)
      if (ntdb == 0){
         cat0(" - Selected trait: ",trait)
         if (! is.null(country)){
            cat0(" - Selected countries: " ,paste(country  ,collapse="; "))
         }#end if (! is.null(country))
         if (! is.null(continent)){
            cat0(" - Selected continents: ",paste(continent,collapse="; "))
         }#end if (! is.null(continent))
         stop(paste0("No data available for trait ",trait," in file \""
                    ,basename(tdb.csv),"\"."))
      }#end if (ntdb == 0)
      #------------------------------------------------------------------------------------#
   }#end if (! (trait %in% names(tdb)))
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
      if ((! is.na(species[s])) && length(grep(unknown.genus,species[s])) == 0){
         if (verbose) cat0("   - ",s,"/",nspecies,"  -  ",species[s],".") 

         #----- Get the genus and family of this species, in case we need it. -------------#
         igen        = which (datum$scientific %in% species[s])
         this.genus  = unique(datum$genus[igen])
         this.family = unique(capwords(datum$family[igen],strict=TRUE))
         if (length(this.genus) != 1 || length(this.family) != 1){
            cat0(" - Perhaps a bastard genus?")
            browser()
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Check whether we have biomass for this species.                             #
         #---------------------------------------------------------------------------------#
         if (species[s] %in% tdb$scientific){
            #------------------------------------------------------------------------------#
            #     Find the index of this species in the database, and assign the trait     #
            # from there.                                                                  #
            #------------------------------------------------------------------------------#
            itrait = intersect( which(tdb$scientific %in% species[s] )
                              , intersect( which(tdb$genus  %in% this.genus )
                                         , which(tdb$family %in% this.family) ) )
            if (length(itrait) != 1){
               #----- Intersection was zero! Misidentification, maybe? --------------------#
               loose      = TRUE
               cat0("Weird, length(itrait)=",length(itrait),"!")
               browser()
               #---------------------------------------------------------------------------#
            }else{
               if (any(c(FALSE,is.finite(tdb[[trait]][itrait])),na.rm=TRUE)){
                  sel       = datum$scientific %in% species[s]
                  datum[[trait]]   [sel] = tdb[[trait]][itrait]
                  datum[[gf.trait]][sel] = 0L
                  datum$country    [sel] = tdb$country  [itrait]
                  datum$continent  [sel] = tdb$continent[itrait]
                  fill.genus              = FALSE
                  loose                   = FALSE
               }else{
                  loose                   = TRUE
                  fill.genus              = this.genus %in% tdb$genus
               }#end if
               #---------------------------------------------------------------------------#
            }#end if
            #------------------------------------------------------------------------------#
         }else{
            loose      = TRUE
            fill.genus = this.genus %in% tdb$genus
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
            itrait  = intersect( which( tdb$genus  %in% this.genus )
                               , which( tdb$family %in% this.family)
                              )#end intersect

            #----- Find individuals that belong to this species. --------------------------#
            sel  = datum$scientific %in% species[s]
            nsel = sum(sel)
            #------------------------------------------------------------------------------#

            if (length(itrait) == 0){
               #----- Intersection was zero! Misidentification, maybe? --------------------#
               loose      = TRUE
               #---------------------------------------------------------------------------#
            }else if (nsel > 0){
               #------ Select data that can be used for weighting. ------------------------#
               t.value     = tdb[[trait  ]][itrait]
               t.weight    = tdb[[n.trait]][itrait] / sum(tdb[[n.trait]][itrait])
               t.continent = tdb$continent [itrait]
               t.country   = tdb$country   [itrait]
               t.idx       = seq_along(t.value)
               #---------------------------------------------------------------------------#


               #----- Decide whether to use sample or average. ----------------------------#
               if (fsample){
                  smp.idx       = lit.sample(x=t.idx,size=nsel,replace=TRUE,prob=t.weight)
                  smp.trait     = t.value    [smp.idx]
                  smp.country   = t.country  [smp.idx]
                  smp.continent = t.continent[smp.idx]
               }else{
                  smp.trait     = weighted.mean     (x=t.value    ,w=t.weight,na.rm=TRUE)
                  smp.country   = weighted.commonest(x=t.country  ,w=t.weight,na.rm=TRUE)
                  smp.continent = weighted.commonest(x=t.continent,w=t.weight,na.rm=TRUE)
               }#end if (fsample)
               #---------------------------------------------------------------------------#

               datum[[trait]]   [sel] = smp.trait
               datum[[gf.trait]][sel] = 1L
               datum$country    [sel] = smp.country
               datum$continent  [sel] = smp.continent
               sci.genus.mean         = rbind(sci.genus.mean
                                             ,c(species[s],this.genus,this.family))
               if (verbose){
                  cat0("     * Species ",species[s]," not found."
                                        ,"  Use genus average instead.")
               }#end if
               loose = FALSE
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
   sci.genus.mean = sci.genus.mean[order(sci.genus.mean[,"scientific"]),,drop=FALSE]
   genus.only     = grepl( pattern     = " NA$"
                         , x           = sci.genus.mean[,"scientific",drop=FALSE]
                         , ignore.case = TRUE
                         )#end grepl
   if (verbose && (nrow(sci.genus.mean[!genus.only,,drop=FALSE]) > 0)){
      cat0("")
      cat0("-----------------------------------------------------------------------")
      cat0(" Found species that were not in trait data base (check for synonyms)!")
      print(sci.genus.mean[! genus.only,,drop=FALSE],quote=FALSE)
      cat0("-----------------------------------------------------------------------")
      cat0("")
   }#end if (verbose && nrow(sci.genus.mean[!genus.only,,drop=FALSE]) > 0)
   if (verbose && (nrow(sci.genus.mean[genus.only,,drop=FALSE]) > 0)){
      cat0("")
      cat0("-----------------------------------------------------------------------")
      cat0("Only genus was provided: filled with average genus value:")
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
      cat0(" Found genera from families with no valid data in the trait data base!")
      print(sci.loose.mean,quote=FALSE)
      cat0("-----------------------------------------------------------------------")
      cat0("")
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Second loop: we list all families, and look for individuals that have no genus    #
   # associated. We compute the mean trait of the known individuals for that family and    #
   # use that as an estimate of the trait.                                                 #
   #---------------------------------------------------------------------------------------#
   #----- Separate all families. ----------------------------------------------------------#
   families          = unique(datum$family)
   nfamilies         = length(families)
   sci.family.sample = matrix(nrow=0,ncol=3
                             ,dimnames=list(NULL,c("scientific","family",trait)))
   #----- Loop over all families. ---------------------------------------------------------#
   for (f in sequence(nfamilies)){
      if (! (families[f] %in% unknown.family)){
         if (verbose) cat0("   - ",f,"/",nfamilies,"  -  ",families[f],".")

         #----- Get the individuals that belong to this family. ---------------------------#
         ifam  = which (datum$family %in% families[f] & is.finite(datum[[trait]]) )
         imiss = which (datum$family %in% families[f] & is.na    (datum[[trait]]) )
         nmiss = length(imiss)
         #---------------------------------------------------------------------------------#

         if (length(imiss) > 0 && length(ifam) > 0){
            #------ Select data that can be used for weighting. ---------------------------#
            t.value     = datum[[trait]] [ifam]
            t.weight    = wgtfac         [ifam]
            t.continent = datum$continent[ifam]
            t.country   = datum$country  [ifam]
            t.idx       = seq_along(t.value)
            #------------------------------------------------------------------------------#



            #----- Decide whether to use sample or average. -------------------------------#
            if (fsample){
               smp.idx       = lit.sample(x=t.idx,size=nmiss,replace=TRUE,prob=t.weight)
               smp.trait     = t.value    [smp.idx]
               smp.country   = t.country  [smp.idx]
               smp.continent = t.continent[smp.idx]
            }else{
               smp.trait     = weighted.mean     (x=t.value    ,w=t.weight,na.rm=TRUE)
               smp.country   = weighted.commonest(x=t.country  ,w=t.weight,na.rm=TRUE)
               smp.continent = weighted.commonest(x=t.continent,w=t.weight,na.rm=TRUE)
            }#end if (fsample)
            #------------------------------------------------------------------------------#



            #------ Fill in with the sample/mean. -----------------------------------------#
            datum[[trait]]   [imiss] = smp.trait
            datum[[gf.trait]][imiss] = 2L
            datum$country    [imiss] = smp.country
            datum$continent  [imiss] = smp.continent
            gf2                      = cbind(datum$scientific[imiss]
                                            ,datum$family    [imiss]
                                            ,sprintf("%6.3f",smp.trait)
                                            )#end cbind
            sci.family.sample        = rbind(sci.family.sample,gf2)
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
         cat0(" Found families with unidentified genera!")
         print(sci.family.sample,quote=FALSE)
      }#end if
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Final block.  We fill in the trait for unknown individuals, by randomly sampling   #
   # from the individuals we know the density.                                             #
   #---------------------------------------------------------------------------------------#
   imiss  = which(is.na(datum[[trait]]))
   nmiss  = length(imiss)
   nvalid = nrow(datum) - nmiss

   if (nmiss > 0){
      if (verbose){
         cat0(" The following families are filled with global sampling: ")
         fam.global.sampling = t(t(sort(unique(datum$family[imiss]))))
         print(fam.global.sampling,quote=FALSE)
      }#end if
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Check whether this is going to be a complete guess or something slightly more #
      # elegant.                                                                           #
      #------------------------------------------------------------------------------------#
      if (nvalid == 0){
         #----- Use any data from the trait data base. ------------------------------------#
         warning(" None of the trees are known!  Trait-filling is going to be very crude!")
         itrait      = which(tdb$life.type %in% "T")
         t.value     = tdb[[trait  ]][itrait]
         t.weight    = tdb[[n.trait]][itrait] / sum(tdb[[n.trait]][itrait])
         t.continent = tdb$continent [itrait]
         t.country   = tdb$country   [itrait]
         t.idx       = seq_along(t.value)
         #---------------------------------------------------------------------------------#
      }else{

         #------ Select data that can be used for weighting. ------------------------------#
         t.value     = datum[[trait]] [-imiss]
         t.weight    = wgtfac         [-imiss]
         t.continent = datum$continent[-imiss]
         t.country   = datum$country  [-imiss]
         t.idx       = seq_along(t.value)
         #---------------------------------------------------------------------------------#
      }#end if (nvalid == 0)
      #------------------------------------------------------------------------------------#


      #----- Decide whether to use sample or averaged values. -----------------------------#
      if (fsample){
         smp.idx       = lit.sample(x=t.idx,size=nmiss,replace=TRUE,prob=t.weight)
         smp.trait     = t.value    [smp.idx]
         smp.country   = t.country  [smp.idx]
         smp.continent = t.continent[smp.idx]
      }else{
         smp.trait     = weighted.mean     (x=t.value    ,w=t.weight,na.rm=TRUE)
         smp.country   = weighted.commonest(x=t.country  ,w=t.weight,na.rm=TRUE)
         smp.continent = weighted.commonest(x=t.continent,w=t.weight,na.rm=TRUE)
      }#end if (fsample)
      #------------------------------------------------------------------------------------#


      #----- Fill the remaining gaps. -----------------------------------------------------#
      datum[[trait]]   [imiss] = smp.trait
      datum$country    [imiss] = smp.country
      datum$continent  [imiss] = smp.continent
      datum[[gf.trait]][imiss] = 3L
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Adjust the gap-filling flag for trait by adding whether the scientific name       #
   # itself was gap-filled.                                                                #
   #---------------------------------------------------------------------------------------#
   if ("gf.scientific" %in% names(datum)){
      datum[[gf.trait]] = datum[[gf.trait]] + as.integer(10 * datum$gf.scientific)
   }#end if ("gf.scientific" %in% names(datum))
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Assign a plant functional type based on the wood density.                         #
   #---------------------------------------------------------------------------------------#
   if (trait %in% "wood.dens"){
      pft.cut    = cut(datum[[trait]],breaks=pft.breaks)
      pft.levels = levels(pft.cut)
      pft.idx    = match(pft.cut,pft.levels)
      datum$pft  = mypfts[pft.idx]
   }#end if (trait %in% "wood.dens")
   #---------------------------------------------------------------------------------------#
   return(datum)
}#end function find.trait
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
   g.s[is.na(g) & is.na(s)] = NA_character_
   look.up$scientific       = g.s
   look.up$genus            = g
   #---------------------------------------------------------------------------------------#



   #----- Trim the common names, and simplify/replace some names. -------------------------#
   datum$common                      = tolower(trim(datum$common))
   datum$common[is.na(datum$common)] = unknown.common
   #---------------------------------------------------------------------------------------#



   #----- Find all unique common names. ---------------------------------------------------#
   unique.common = unique(datum$common)
   n.common      = length(unique.common)
   notfound      = NULL
   for (n in 1:n.common){
      #----- Find the trees that have the same common name as this one. -------------------#
      cat0(" - ",n,"/",n.common," -- ",unique.common[n],".")
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
         datum$scientific   [w.dat] = unknown.scientific
         datum$genus        [w.dat] = unknown.genus
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
   g.s[is.na(g) & is.na(s)] = NA_character_
   look.up$scientific       = g.s
   look.up$genus            = g
   #---------------------------------------------------------------------------------------#



   #----- Trim the common names, and simplify/replace some names. -------------------------#
   datum$common                      = tolower(trim(datum$common))
   datum$common[is.na(datum$common)] = unknown.common
   #---------------------------------------------------------------------------------------#



   #----- Find all unique common names. ---------------------------------------------------#
   unique.common = unique(datum$common)
   n.common      = length(unique.common)
   notfound      = NULL
   for (n in 1:n.common){
      #----- Find the trees that have the same common name as this one. -------------------#
      cat0(" - ",n,"/",n.common," -- ",unique.common[n],".")
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
         datum$scientific   [w.dat] = unknown.scientific
         datum$genus        [w.dat] = unknown.genus
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
keep.gen.spe.only <<- function(x,out=c("both","genus","species"),is.debug=FALSE){
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
      if (is.na(x) | x %in% c("",unknown.genus,unknown.scientific)) x=unknown.scientific


      #---- Temporarily make the name lower case, and remove any dots that may exist. -----#
      x = tolower(x)
      x = gsub(pattern="\\.",replacement="",ignore.case=TRUE,x=x)
      #------------------------------------------------------------------------------------#

      #----- In case there is an " x ", attach it to the epithet. -------------------------#
      x = gsub(pattern="\\ x\\ ",replacement=" x_",ignore.case=TRUE,x=x)
      #------------------------------------------------------------------------------------#

      #----- Separate the bogus _form2 that appears in some Paracou species. --------------#
      x = gsub(pattern="_form2$",replacement=" form2",x=x)
      #------------------------------------------------------------------------------------#


      #----- Split words, and remove stuff like sp, cf, etc. ------------------------------#
      ans  = unlist(strsplit(unlist(c(tolower(x))),split=" "))
      wild = ans %in% c(unknown.wildcard,unknown.epithet)
      ans  = ans[! wild]
      nans = length(ans)
      #------------------------------------------------------------------------------------#


      #----- Decide whether to keep the genus, species, or both. --------------------------#
      gen = unknown.genus
      spe = unknown.epithet
      if (nans > 1){
         gen = capwords(ans[1],strict=TRUE)
         spe = tolower(ans[2])
      }else if (nans == 1){
         gen = capwords(ans[1],strict=TRUE)
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

      #----- Separate the x from the epithet for cross-species. ---------------------------#
      ans = gsub(pattern="x_",replacement="x ",x=ans)
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#



   return(ans)
}#end function keep.gen.spe.only
#==========================================================================================#
#==========================================================================================#
