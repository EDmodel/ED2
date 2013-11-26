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

   #----- Full replacements. --------------------------------------------------------------#
   sel = (! is.na(x) & x == "araticu"                 ); x[sel] = "araticum"
   sel = (! is.na(x) & x == "baubarana"               ); x[sel] = "embaubarana"
   sel = (! is.na(x) & x == "breu/louro preto?"       ); x[sel] = NA
   sel = (! is.na(x) & x == "babao"                   ); x[sel] = "macauba"
   sel = (! is.na(x) & x == "bolao"                   ); x[sel] = "fava-bolota"
   sel = (! is.na(x) & x == "brau"                    ); x[sel] = "breu"
   sel = (! is.na(x) & x == "cauba"                   ); x[sel] = "macacauba"
   sel = (! is.na(x) & x == "cajuba"                  ); x[sel] = "caju"
   sel = (! is.na(x) & x == "castanha"                ); x[sel] = "castanha do para"
   sel = (! is.na(x) & x == "castanheiro"             ); x[sel] = "castanha do para"
   sel = (! is.na(x) & x == "cipo(dbh a 0.9m do chao)"); x[sel] = "cipo"
   sel = (! is.na(x) & x == "cipo+arapo"              ); x[sel] = "cipo"
   sel = (! is.na(x) & x == "cutiti"                  ); x[sel] = "abiu cutite"
   sel = (! is.na(x) & x == "cutite"                  ); x[sel] = "abiu cutite"
   sel = (! is.na(x) & x == "cruzeiro"                ); x[sel] = "quina-cruzeiro"
   sel = (! is.na(x) & x == "envira surucu"           ); x[sel] = "envira surucucu"
   sel = (! is.na(x) & x == "fiora preta"             ); x[sel] = NA
   sel = (! is.na(x) & x == "jambo"                   ); x[sel] = "jambo-do-mato"
   sel = (! is.na(x) & x == "jara"                    ); x[sel] = "jarana"
   sel = (! is.na(x) & x == "jito"                    ); x[sel] = "gito"
   sel = (! is.na(x) & x == "louro?"                  ); x[sel] = "louro"
   sel = (! is.na(x) & x == "maracatia"               ); x[sel] = "muiracatiara"
   sel = (! is.na(x) & x == "mata pau+jito"           ); x[sel] = "gito"
   sel = (! is.na(x) & x == "muruci"                  ); x[sel] = "muruci da mata"
   sel = (! is.na(x) & x == "palmito"                 ); x[sel] = "acai"
   sel = (! is.na(x) & x == "palmito babosa"          ); x[sel] = "acai"
   sel = (! is.na(x) & x == "patua"                   ); x[sel] = "pataua"
   sel = (! is.na(x) & x == "quariquara"              ); x[sel] = "acariquara"
   sel = (! is.na(x) & x == "tachi preto ???"         ); x[sel] = "tachi preto"
   sel = (! is.na(x) & x == "tachi preto folh"        ); x[sel] = "tachi preto"
   sel = (! is.na(x) & x == "tento folha"             ); x[sel] = "tento"
   sel = (! is.na(x) & x == "saboeiro"                ); x[sel] = "fava-saboeiro"
   sel = (! is.na(x) & x == "saboiera"                ); x[sel] = "fava-saboeiro"
   sel = (! is.na(x) & x == "seringa"                 ); x[sel] = "seringueira"
   sel = (! is.na(x) & x == "sajinera"                ); x[sel] = NA
   sel = (! is.na(x) & x == "sorveira"                ); x[sel] = "sorva"
   sel = (! is.na(x) & x == "sorvo"                   ); x[sel] = "sorva"
   sel = (! is.na(x) & x == "sova"                    ); x[sel] = "sorva"
   sel = (! is.na(x) & x == "tatapiririca verm."      ); x[sel] = "tatapiririca vermelha"
   sel = (! is.na(x) & x == "umbia"                   ); x[sel] = "goiabarana"
   #---------------------------------------------------------------------------------------#



   #----- Substitutions.  Only things that are not part of other words. -------------------#
   x = sub(pattern="abacaba"                ,replacement="bacaba"                ,x=x)
   x = sub(pattern="abicuiba"               ,replacement="ucuuba"                ,x=x)
   x = sub(pattern="abiui"                  ,replacement="abiu"                  ,x=x)
   x = sub(pattern="abiu casca grossa"      ,replacement="abiu-casca-grossa"     ,x=x)
   x = sub(pattern="abiu cutite folha verde",replacement="abiu cutite"           ,x=x)
   x = sub(pattern="abiu mangabinha"        ,replacement="abiu-mangabinha"       ,x=x)
   x = sub(pattern="abiu vermelha"          ,replacement="abiu vermelho"         ,x=x)
   x = sub(pattern="abiuarana vermelha"     ,replacement="abiurana vermelha"     ,x=x)
   x = sub(pattern="abiuarana vermelho"     ,replacement="abiurana vermelha"     ,x=x)
   x = sub(pattern="abiuarana vermlho"      ,replacement="abiurana vermelha"     ,x=x)
   x = sub(pattern="abiurana vermelho"      ,replacement="abiurana vermelha"     ,x=x)
   x = sub(pattern="abiuarana"              ,replacement="abiurana"              ,x=x)
   x = sub(pattern="acoita cavalo"          ,replacement="acoita-cavalo"         ,x=x)
   x = sub(pattern="algodoeira"             ,replacement="sumauma"               ,x=x)
   x = sub(pattern="amalelinha"             ,replacement="amarelinho"            ,x=x)
   x = sub(pattern="amarelinha"             ,replacement="amarelinho"            ,x=x)
   x = sub(pattern="amerelinho"             ,replacement="amarelinho"            ,x=x)
   x = sub(pattern="angelim margoso"        ,replacement="angelim amargoso"      ,x=x)
   x = sub(pattern="angelim pedro"          ,replacement="angelim pedra"         ,x=x)
   x = sub(pattern="apuii"                  ,replacement="apui"                  ,x=x)
   x = sub(pattern="araca nego"             ,replacement="araca"                 ,x=x)
   x = sub(pattern="barbatimao"             ,replacement="fava-barbatimao"       ,x=x)
   x = sub(pattern="bolao"                  ,replacement="fava bolota"           ,x=x)
   x = sub(pattern="brejauba"               ,replacement="brejauva"              ,x=x)
   x = sub(pattern="breu sucuuba"           ,replacement="breu sucuruba"         ,x=x)
   x = sub(pattern="cabeca de urubu"        ,replacement="cabeca-de-urubu"       ,x=x)
   x = sub(pattern="cabela"                 ,replacement="louro canela"          ,x=x)
   x = sub(pattern="cabriuna"               ,replacement="cabriuva"              ,x=x)
   x = sub(pattern="cacau bravo"            ,replacement="cacaui"                ,x=x)
   x = sub(pattern="cachudinha"             ,replacement="cascudinha"            ,x=x)
   x = sub(pattern="calcho"                 ,replacement="caucho"                ,x=x)
   x = sub(pattern="canela de velho"        ,replacement="canela-de-velho"       ,x=x)
   x = sub(pattern="canelha velha"          ,replacement="canela-de-velho"       ,x=x)
   x = sub(pattern="canella de jacami"      ,replacement="canela-de-jacamim"     ,x=x)
   x = sub(pattern="canella ge jacami"      ,replacement="canela-de-jacamim"     ,x=x)
   x = sub(pattern="canniela"               ,replacement="canela"                ,x=x)
   x = sub(pattern="carobia"                ,replacement="caroba"                ,x=x)
   x = sub(pattern="cascudinho"             ,replacement="cascudinha"            ,x=x)
   x = sub(pattern="cascudo"                ,replacement="cascudinha"            ,x=x)
   x = sub(pattern="castanha de sapocaia"   ,replacement="castanha sapucaia"     ,x=x)
   x = sub(pattern="castanha de sapucaia"   ,replacement="castanha sapucaia"     ,x=x)
   x = sub(pattern="castanha sapocaia"      ,replacement="castanha sapucaia"     ,x=x)
   x = sub(pattern="cauxo"                  ,replacement="caucho"                ,x=x)
   x = sub(pattern="caxeta"                 ,replacement="caixeta"               ,x=x)
   x = sub(pattern="caximbeira"             ,replacement="cachimbeiro"           ,x=x)
   x = sub(pattern="caximbeiro"             ,replacement="cachimbeiro"           ,x=x)
   x = sub(pattern="caxudinha"              ,replacement="cascudinha"            ,x=x)
   x = sub(pattern="chocolate"              ,replacement="cacau"                 ,x=x)
   x = sub(pattern="coracao  de negro"      ,replacement="coracao-de-negro"      ,x=x)
   x = sub(pattern="coracao de nego"        ,replacement="coracao-de-negro"      ,x=x)
   x = sub(pattern="coracao de negro"       ,replacement="coracao-de-negro"      ,x=x)
   x = sub(pattern="corante de indio"       ,replacement="urucum"                ,x=x)
   x = sub(pattern="coro preto"             ,replacement="louro preto"           ,x=x)
   x = sub(pattern="coussarea racemosa"     ,replacement="caferana"              ,x=x)
   x = sub(pattern="cutiti"                 ,replacement="cutite"                ,x=x)
   x = sub(pattern="cumari"                 ,replacement="cumaru"                ,x=x)
   x = sub(pattern="cumaru/apui"            ,replacement="apui"                  ,x=x)
   x = sub(pattern="embauba branco"         ,replacement="embauba branca"        ,x=x)
   x = sub(pattern="embauba vick"           ,replacement="embauba"               ,x=x)
   x = sub(pattern="embirata"               ,replacement="envira ata"            ,x=x)
   x = sub(pattern="embireira branca"       ,replacement="envira"                ,x=x)
   x = sub(pattern="embireira rosa"         ,replacement="envira"                ,x=x)
   x = sub(pattern="envira preto"           ,replacement="envira preta"          ,x=x)
   x = sub(pattern="envira vermelho"        ,replacement="envira vermelha"       ,x=x)
   x = sub(pattern="escorrega macaco"       ,replacement="escorrega-macaco"      ,x=x)
   x = sub(pattern="escurrega macaco"       ,replacement="escorrega-macaco"      ,x=x)
   x = sub(pattern="fava arara tucupi"      ,replacement="fava-arara-tucupi"     ,x=x)
   x = sub(pattern="fava saboeira"          ,replacement="fava-saboeira"         ,x=x)
   x = sub(pattern="feijo branco"           ,replacement="freijo branco"         ,x=x)
   x = sub(pattern="figueira brava"         ,replacement="figueira"              ,x=x)
   x = sub(pattern="gameleiro"              ,replacement="gameleira"             ,x=x)
   x = sub(pattern="gapeba"                 ,replacement="abiu"                  ,x=x)
   x = sub(pattern="guapeba"                ,replacement="abiu"                  ,x=x)
   x = sub(pattern="gema de ovo"            ,replacement="amarelao"              ,x=x)
   x = sub(pattern="genipapo"               ,replacement="jenipapo"              ,x=x)
   x = sub(pattern="goibarana"              ,replacement="goiabarana"            ,x=x)
   x = sub(pattern="gombeira vermelho"      ,replacement="gombeira vermelha"     ,x=x)
   x = sub(pattern="goiaba"                 ,replacement="araca"                 ,x=x)
   x = sub(pattern="guaiaba"                ,replacement="araca"                 ,x=x)
   x = sub(pattern="guaiba"                 ,replacement="araca"                 ,x=x)
   x = sub(pattern="guariuva"               ,replacement="guariuba"              ,x=x)
   x = sub(pattern="ibirucu"                ,replacement="embirucu"              ,x=x)
   x = sub(pattern="imbirata"               ,replacement="envira ata"            ,x=x)
   x = sub(pattern="imbireira"              ,replacement="envira"                ,x=x)
   x = sub(pattern="imbireira rosa"         ,replacement="envira"                ,x=x)
   x = sub(pattern="imbiricu"               ,replacement="embirucu"              ,x=x)
   x = sub(pattern="imbirucu"               ,replacement="embirucu"              ,x=x)
   x = sub(pattern="inga f.p."              ,replacement="inga"                  ,x=x)
   x = sub(pattern="inga vermelha"          ,replacement="inga vermelho"         ,x=x)
   x = sub(pattern="inga titica"            ,replacement="inga xixica"           ,x=x)
   x = sub(pattern="ipe amerelo"            ,replacement="ipe amarelo"           ,x=x)
   x = sub(pattern="jaboticaba"             ,replacement="jabuticaba"            ,x=x)
   x = sub(pattern="jaracatia"              ,replacement="jacaratia"             ,x=x)
   x = sub(pattern="jenita"                 ,replacement="janita"                ,x=x)
   x = sub(pattern="jotobazinho"            ,replacement="jatobazinho"           ,x=x)
   x = sub(pattern="jutai mirim"            ,replacement="jutai-mirim"           ,x=x)
   x = sub(pattern="jutai acu"              ,replacement="jutai-acu"             ,x=x)
   x = sub(pattern="laranginha"             ,replacement="laranjinha"            ,x=x)
   x = sub(pattern="leiteiro"               ,replacement="leiteira"              ,x=x)
   x = sub(pattern="leitera"                ,replacement="leiteira"              ,x=x)
   x = sub(pattern="louro branco"           ,replacement="louro"                 ,x=x)
   x = sub(pattern="loro amarelo"           ,replacement="louro amarelo"         ,x=x)
   x = sub(pattern="mamao jacatia"          ,replacement="jacaratia"             ,x=x)
   x = sub(pattern="maracatiara"            ,replacement="muiracatiara"          ,x=x)
   x = sub(pattern="mata-mata"              ,replacement="matamata"              ,x=x)
   x = sub(pattern="mata mata"              ,replacement="matamata"              ,x=x)
   x = sub(pattern="matamata branca"        ,replacement="matamata branco"       ,x=x)
   x = sub(pattern="matamata vermelha"      ,replacement="matamata vermelho"     ,x=x)
   x = sub(pattern="mata caldo"             ,replacement="mata-calado"           ,x=x)
   x = sub(pattern="melanciera"             ,replacement="melancieira"           ,x=x)
   x = sub(pattern="moratinga"              ,replacement="muiratinga"            ,x=x)
   x = sub(pattern="muratinga"              ,replacement="muiratinga"            ,x=x)
   x = sub(pattern="morta"                  ,replacement="defunta"               ,x=x)
   x = sub(pattern="murucidu mata"          ,replacement="muruci da mata"        ,x=x)
   x = sub(pattern="muiriatinga"            ,replacement="muiratinga"            ,x=x)
   x = sub(pattern="mutama"                 ,replacement="mutambo"               ,x=x)
   x = sub(pattern="mutamba"                ,replacement="mutambo"               ,x=x)
   x = sub(pattern="ocooba"                 ,replacement="ucuuba"                ,x=x)
   x = sub(pattern="ocuuba"                 ,replacement="ucuuba"                ,x=x)
   x = sub(pattern="ouro branco"            ,replacement="seringueira"           ,x=x)
   x = sub(pattern="papa terra"             ,replacement="papaterra"             ,x=x)
   x = sub(pattern="papa-terra"             ,replacement="papaterra"             ,x=x)
   x = sub(pattern="para para"              ,replacement="parapara"              ,x=x)
   x = sub(pattern="para-para"              ,replacement="parapara"              ,x=x)
   x = sub(pattern="papo de mutum"          ,replacement="pato-de-mutum"         ,x=x)
   x = sub(pattern="passarinhiera"          ,replacement="passarinheira"         ,x=x)
   x = sub(pattern="pata de vaca"           ,replacement="pata-de-vaca"          ,x=x)
   x = sub(pattern="paineira"               ,replacement="sumauma"               ,x=x)
   x = sub(pattern="pao de sangue"          ,replacement="pau-sangue"            ,x=x)
   x = sub(pattern="pau de arco"            ,replacement="pau-de-arco"           ,x=x)
   x = sub(pattern="pau d.arco"             ,replacement="pau-de-arco"           ,x=x)
   x = sub(pattern="pau de cobra"           ,replacement="pau-cobra"             ,x=x)
   x = sub(pattern="pau de colher"          ,replacement="pau-de-colher"         ,x=x)
   x = sub(pattern="pau de jacare"          ,replacement="pau-jacare"            ,x=x)
   x = sub(pattern="pau de remo"            ,replacement="pau-de-remo"           ,x=x)
   x = sub(pattern="pau de sangue"          ,replacement="pau-sangue"            ,x=x)
   x = sub(pattern="pau jacare"             ,replacement="pau-jacare"            ,x=x)
   x = sub(pattern="pau para tudo"          ,replacement="pau-para-tudo"         ,x=x)
   x = sub(pattern="pau pereira"            ,replacement="peroba mica"           ,x=x)
   x = sub(pattern="pau sangue"             ,replacement="pau-sangue"            ,x=x)
   x = sub(pattern="pente de macaco"        ,replacement="pente-de-macaco"       ,x=x)
   x = sub(pattern="perna de moca"          ,replacement="perna-de-moca"         ,x=x)
   x = sub(pattern="piquiazeiro"            ,replacement="piquia"                ,x=x)
   x = sub(pattern="piqui rosa"             ,replacement="piquia"                ,x=x)
   x = sub(pattern="prapara"                ,replacement="parapara"              ,x=x)
   x = sub(pattern="quaiquara"              ,replacement="acariquara"            ,x=x)
   x = sub(pattern="quariquari"             ,replacement="acariquara"            ,x=x)
   x = sub(pattern="quari quari"            ,replacement="acariquara"            ,x=x)
   x = sub(pattern="quina cruzeiro"         ,replacement="quina-cruzeiro"        ,x=x)
   x = sub(pattern="quina quina"            ,replacement="quinaquina"            ,x=x)
   x = sub(pattern="quina-quina"            ,replacement="quinaquina"            ,x=x)
   x = sub(pattern="muida"                  ,replacement="miuda"                 ,x=x)
   x = sub(pattern="roxao"                  ,replacement="roxinho"               ,x=x)
   x = sub(pattern="roxinao"                ,replacement="roxinho"               ,x=x)
   x = sub(pattern="segador"                ,replacement="cegador"               ,x=x)
   x = sub(pattern="seritinga"              ,replacement="seringueira"           ,x=x)
   x = sub(pattern="sorveira leite"         ,replacement="sorva"                 ,x=x)
   x = sub(pattern="taxi"                   ,replacement="tachi"                 ,x=x)
   x = sub(pattern="tachi branca"           ,replacement="tachi branco"          ,x=x)
   x = sub(pattern="tachi preta"            ,replacement="tachi preto"           ,x=x)
   x = sub(pattern="tachi vermelha"         ,replacement="tachi vermelho"        ,x=x)
   x = sub(pattern="tauri"                  ,replacement="tauari"                ,x=x)
   x = sub(pattern="talquari"               ,replacement="tauari"                ,x=x)
   x = sub(pattern="tamarindu"              ,replacement="azedinha"              ,x=x)
   x = sub(pattern="tamarino"               ,replacement="azedinha"              ,x=x)
   x = sub(pattern="tento foha grauda"      ,replacement="tento folha grauda"    ,x=x)
   x = sub(pattern="ucuarana"               ,replacement="urucurana"             ,x=x)
   x = sub(pattern="ucuuba terra firme"     ,replacement="ucuuba-terra-firme"    ,x=x)
   x = sub(pattern="ucuuba tf"              ,replacement="ucuuba-terra-firme"    ,x=x)
   x = sub(pattern="uruucurana"             ,replacement="urucurana"             ,x=x)
   x = sub(pattern="ucuuba vermelho"        ,replacement="ucuuba vermelha"       ,x=x)
   x = sub(pattern="uchi"                   ,replacement="uxi"                   ,x=x)
   x = sub(pattern="unha de vaca"           ,replacement="pata-de-vaca"          ,x=x)
   x = sub(pattern="verdadiera"             ,replacement="verdadeira"            ,x=x)
   x = sub(pattern="xixua"                  ,replacement="chichua"               ,x=x)
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
   dat$scientific = sub("cf.","",x=dat$scientific)
   #----- Break into genus and species. ---------------------------------------------------#
   gs.list                  = sapply(X = tolower(dat$scientific),FUN=strsplit,split=" ")
   gs.length                = sapply(X = gs.list, FUN = length)
   gs.mat                   = do.call(rbind,gs.list)
   g.only                   = gs.length < 2
   gs.mat[g.only,2]         = NA
   g                        = capwords(gs.mat[,1],strict=TRUE)
   s                        = tolower(gs.mat[,2])
   g.s                      = paste(g,s,sep=" ")
   g.s[is.na(g) & is.na(s)] = NA
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #      Full substitutions (only if the entire name matches).                            #
   #---------------------------------------------------------------------------------------#
   g.s[g.s == "Dialium sp" ] = "Dialium"
   g.s[g.s == "Eugenia sp."] = "Eugenia"
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #      Substitutions of scientific name.  We standardise these first so family becomes  #
   # easier.                                                                               #
   #---------------------------------------------------------------------------------------#
   g.s = sub("Agonadra sp."               ,"Agonandra"                    ,x=g.s)
   g.s = sub("Aiouea densiflora"          ,"Aiouea laevis"                ,x=g.s)
   g.s = sub("Allophyllus floribunda"     ,"Allophylus floribundus"       ,x=g.s)
   g.s = sub("Ampelocera endentula"       ,"Ampelocera edentula"          ,x=g.s)
   g.s = sub("Amphirrhox surinamensis"    ,"Amphirrhox longifolia"        ,x=g.s)
   g.s = sub("Anadenanthera falcata"      ,"Anadenanthera peregrina"      ,x=g.s)
   g.s = sub("Aniba roseodora"            ,"Aniba rosaeodora"             ,x=g.s)
   g.s = sub("Aniba sp."                  ,"Aniba"                        ,x=g.s)
   g.s = sub("Annona decicoma"            ,"Annona densicoma"             ,x=g.s)
   g.s = sub("Apeiba burchelii"           ,"Apeiba glabra"                ,x=g.s)
   g.s = sub("Aspidosperma aracanga"      ,"Aspidosperma araracanga"      ,x=g.s)
   g.s = sub("Aspidosperma auriculata"    ,"Aspidosperma auriculatum"     ,x=g.s)
   g.s = sub("Aspidosperma desmantum"     ,"Aspidosperma desmanthum"      ,x=g.s)
   g.s = sub("Aspidosperma desmathum"     ,"Aspidosperma desmanthum"      ,x=g.s)
   g.s = sub("Aspidosperma eteanun"       ,"Aspidosperma eteanum"         ,x=g.s)
   g.s = sub("Aspidosperma nitidum"       ,"Aspidosperma excelsum"        ,x=g.s)
   g.s = sub("Astronium le-cointei"       ,"Astronium lecointei"          ,x=g.s)
   g.s = sub("Austroplenckia populnea"    ,"Plenckia populnea"            ,x=g.s)
   g.s = sub("Balisia pedicelares"        ,"Albizia pedicellaris"         ,x=g.s)
   g.s = sub("Balizia pedicellaris"       ,"Albizia pedicellaris"         ,x=g.s)
   g.s = sub("Bauhinia jarensis"          ,"Bauhinia"                     ,x=g.s)
   g.s = sub("Bellucia grossulariodis"    ,"Bellucia grossularioides"     ,x=g.s)
   g.s = sub("Brosimum autifolium"        ,"Brosimum acutifolium"         ,x=g.s)
   g.s = sub("Brosimum sp."               ,"Brosimum"                     ,x=g.s)
   g.s = sub("Brosimum lactascens"        ,"Brosimum lactescens"          ,x=g.s)
   g.s = sub("Byrsonima schultesiana"     ,"Byrsonima arthropoda"         ,x=g.s)
   g.s = sub("Byrsonima estipulacea"      ,"Byrsonima stipulacea"         ,x=g.s)
   g.s = sub("Capirona ulei"              ,"Capirona decorticans"         ,x=g.s)
   g.s = sub("Capparis frondosa"          ,"Capparidastrum frondosum"     ,x=g.s)
   g.s = sub("Cecropia distachia"         ,"Cecropia distachya"           ,x=g.s)
   g.s = sub("Cedrela fistis"             ,"Cedrela fissilis"             ,x=g.s)
   g.s = sub("Cedrella odorata"           ,"Cedrela odorata"              ,x=g.s)
   g.s = sub("Chanouchiton kapleri"       ,"Chaunochiton kappleri"        ,x=g.s)
   g.s = sub("Chaunochiton kapleri"       ,"Chaunochiton kappleri"        ,x=g.s)
   g.s = sub("Clarisia elicifolia"        ,"Clarisia ilicifolia"          ,x=g.s)
   g.s = sub("Compomanesia"               ,"Campomanesia"                 ,x=g.s)
   g.s = sub("Connarus perrottetes"       ,"Connarus perrottetii"         ,x=g.s)
   g.s = sub("Connarus spermattetii"      ,"Connarus perrotettii"         ,x=g.s)
   g.s = sub("Cordia scrabida"            ,"Cordia exaltata"              ,x=g.s)
   g.s = sub("Coupeia sp."                ,"Couepia"                      ,x=g.s)
   g.s = sub("Coussarea racemosa"         ,"Coussarea albescens"          ,x=g.s)
   g.s = sub("Crepidospermum gondotiano"  ,"Crepidospermum goudotianum"   ,x=g.s)
   g.s = sub("Crepidospermum goudotiano"  ,"Crepidospermum goudotianum"   ,x=g.s)
   g.s = sub("Cupania hirta"              ,"Cupania hirsuta"              ,x=g.s)
   g.s = sub("Cybistax antisiphyllitica"  ,"Cybistax antisyphilitica"     ,x=g.s)
   g.s = sub("Dendrobangia sp."           ,"Dendrobangia"                 ,x=g.s)
   g.s = sub("Dendrobrangea boliviana"    ,"Dendrobangia boliviana"       ,x=g.s)
   g.s = sub("Dialium sp."                ,"Dialium"                      ,x=g.s)
   g.s = sub("Dialium guianensis"         ,"Dialium guianense"            ,x=g.s)
   g.s = sub("Diclinanonna matogrossensis","Diclinanona matogrossensis"   ,x=g.s)
   g.s = sub("Didymopanax vinosum"        ,"Schefflera vinosa"            ,x=g.s)
   g.s = sub("Diospyros praetermissa?"    ,"Diospyros vestita"            ,x=g.s)
   g.s = sub("Diospyros praetermissa"     ,"Diospyros vestita"            ,x=g.s)
   g.s = sub("Diploon venezuelano"        ,"Diploon cuspidatum"           ,x=g.s)
   g.s = sub("Dipteriyx odorata"          ,"Dipteryx odorata"             ,x=g.s)
   g.s = sub("Endlecheria"                ,"Endlicheria"                  ,x=g.s)
   g.s = sub("Enterolobium schobburgkii"  ,"Enterolobium schomburgkii"    ,x=g.s)
   g.s = sub("Enterolobium schombburgkii" ,"Enterolobium schomburgkii"    ,x=g.s)
   g.s = sub("Ephedrantus parviflorus"    ,"Ephedranthus parviflorus"     ,x=g.s)
   g.s = sub("Ephedrantus parvifolius"    ,"Ephedranthus parviflorus"     ,x=g.s)
   g.s = sub("Eriotheca ni"               ,"Eriotheca"                    ,x=g.s)
   g.s = sub("Eriotheca sp"               ,"Eriotheca"                    ,x=g.s)
   g.s = sub("Eriotreca globosa"          ,"Eriotheca globosa"            ,x=g.s)
   g.s = sub("Eschweilera amazomica"      ,"Eschweilera amazonica"        ,x=g.s)
   g.s = sub("Eschweilera apiculatum"     ,"Eschweilera apiculata"        ,x=g.s)
   g.s = sub("Eschweilera idatimom"       ,"Lecythis idatimon"            ,x=g.s)
   g.s = sub("Eschweilera observa"        ,"Eschweilera obversa"          ,x=g.s)
   g.s = sub("Eschweilera pedicelata"     ,"Eschweilera pedicellata"      ,x=g.s)
   g.s = sub("Eschweilera sp."            ,"Eschweilera"                  ,x=g.s)
   g.s = sub("Eugenia ni"                 ,"Eugenia"                      ,x=g.s)
   g.s = sub("Eugenia pacumensis"         ,"Austromyrtus ploumensis"      ,x=g.s)
   g.s = sub("Eugenia schumburgkii"       ,"Eugenia lambertiana"          ,x=g.s)
   g.s = sub("Geissospermum velosii"      ,"Geissospermum laeve"          ,x=g.s)
   g.s = sub("Geissospermum velozii"      ,"Geissospermum laeve"          ,x=g.s)
   g.s = sub("Glicoxilum"                 ,"Pouteria oppositifolia"       ,x=g.s)
   g.s = sub("Glycydendrom amazonicus"    ,"Glycydendron amazonicum"      ,x=g.s)
   g.s = sub("Glycydendron amazonicus"    ,"Glycydendron amazonicum"      ,x=g.s)
   g.s = sub("Guarea guianensis"          ,"Guarea"                       ,x=g.s)
   g.s = sub("Guarea subsessiliflora"     ,"Guarea macrophylla"           ,x=g.s)
   g.s = sub("Guatteria cardoniana"       ,"Guatteria recurvisepala"      ,x=g.s)
   g.s = sub("Hymenolobium flavium"       ,"Hymenolobium flavum"          ,x=g.s)
   g.s = sub("Ilex parviflora"            ,"Ilex petiolaris"              ,x=g.s)
   g.s = sub("Inga dibaldiana"            ,"Inga thibaudiana"             ,x=g.s)
   g.s = sub("Inga paraenses"             ,"Inga paraensis"               ,x=g.s)
   g.s = sub("Inga poliphylla"            ,"Inga"                         ,x=g.s)
   g.s = sub("Inga sp."                   ,"Inga"                         ,x=g.s)
   g.s = sub("Jacaraitia espinhosa"       ,"Jacaratia spinosa"            ,x=g.s)
   g.s = sub("Joannesia hevioides"        ,"Joannesia heveoides"          ,x=g.s)
   g.s = sub("Lacmelea aculeata"          ,"Lacmellea aculeata"           ,x=g.s)
   g.s = sub("Laetia procerra"            ,"Laetia procera"               ,x=g.s)
   g.s = sub("Licania densa"              ,"Licania densiflora"           ,x=g.s)
   g.s = sub("Licaria heteromorpha"       ,"Licania heteromorpha"         ,x=g.s)
   g.s = sub("Lucanarea"                  ,"Lacunaria"                    ,x=g.s)
   g.s = sub("Luheopsis duckeana"         ,"Lueheopsis duckeana"          ,x=g.s)
   g.s = sub("Lustema pubescei"           ,"Lacistema pubescens"          ,x=g.s)
   g.s = sub("Mabea sp."                  ,"Mabea"                        ,x=g.s)
   g.s = sub("Maluria"                    ,"Marlierea umbraticola"        ,x=g.s)
   g.s = sub("Maquira callophylla"        ,"Maquira calophylla"           ,x=g.s)
   g.s = sub("Marmaroxylon racemosum"     ,"Zygia racemosa"               ,x=g.s)
   g.s = sub("Maytenos guianensis"        ,"Maytenus guyanensis"          ,x=g.s)
   g.s = sub("Maytenus guianensis"        ,"Maytenus guyanensis"          ,x=g.s)
   g.s = sub("Meia maderensis"            ,"Neea"                         ,x=g.s)
   g.s = sub("Mezelaurus"                 ,"Mezilaurus"                   ,x=g.s)
   g.s = sub("Mezelaurus itauba"          ,"Mezilaurus itauba"            ,x=g.s)
   g.s = sub("Miconia chrysophyllum"      ,"Miconia chrysophylla"         ,x=g.s)
   g.s = sub("Michopholis sp."            ,"Micropholis"                  ,x=g.s)
   g.s = sub("Michopholis venulosa"       ,"Micropholis venulosa"         ,x=g.s)
   g.s = sub("Microphilis"                ,"Micropholis"                  ,x=g.s)
   g.s = sub("Micropholis guianensis"     ,"Micropholis guyanensis"       ,x=g.s)
   g.s = sub("Micropholis sp."            ,"Micropholis"                  ,x=g.s)
   g.s = sub("Microphylis acutangula"     ,"Micropholis acutangula"       ,x=g.s)
   g.s = sub("Microphylis sp."            ,"Micropholis"                  ,x=g.s)
   g.s = sub("Mouriri abnormis"           ,"Votomita guianensis"          ,x=g.s)
   g.s = sub("Myrciaria reticulata"       ,"Myrcia reticulata"            ,x=g.s)
   g.s = sub("Myrcia rutipula"            ,"Myrcia rufipila"              ,x=g.s)
   g.s = sub("Newtonia psilostachya"      ,"Pseudopiptadenia psilostachya",x=g.s)
   g.s = sub("Ocotea baturitensis"        ,"Ocotea"                       ,x=g.s)
   g.s = sub("Ocotea caudata"             ,"Ocotea cernua"                ,x=g.s)
   g.s = sub("Ocotea sp."                 ,"Ocotea"                       ,x=g.s)
   g.s = sub("Omedia perebea"             ,"Perebea mollis"               ,x=g.s)
   g.s = sub("Onichiopetalum amazonico"   ,"Onychopetalum amazonicum"     ,x=g.s)
   g.s = sub("Pausandra densiflora"       ,"Pausandra trianae"            ,x=g.s)
   g.s = sub("Pouroma guianensis"         ,"Pourouma guianensis"          ,x=g.s)
   g.s = sub("Pouruma guianensis"         ,"Pourouma guianensis"          ,x=g.s)
   g.s = sub("Pouteria biloculares"       ,"Pouteria bilocularis"         ,x=g.s)
   g.s = sub("Pouteria filipis"           ,"Pouteria filipes"             ,x=g.s)
   g.s = sub("Pouteria lasiocarpa"        ,"Pouteria caimito"             ,x=g.s)
   g.s = sub("Pouteria heterosepala"      ,"Pouteria polysepala"          ,x=g.s)
   g.s = sub("Pouteria paraensis"         ,"Pouteria macrocarpa"          ,x=g.s)
   g.s = sub("Pouteria sp."               ,"Pouteria"                     ,x=g.s)
   g.s = sub("Protium heptafilum"         ,"Protium heptaphyllum"         ,x=g.s)
   g.s = sub("Protium sp."                ,"Protium"                      ,x=g.s)
   g.s = sub("Prumes myrtifoliu"          ,"Prunus myrtifolia"            ,x=g.s)
   g.s = sub("Prunus myrtifolius"         ,"Prunus myrtifolia"            ,x=g.s)
   g.s = sub("Pseudolmedia murure"        ,"Pseudolmedia macrophylla"     ,x=g.s)
   g.s = sub("Ptecelobuim jucumba"        ,"Abarema jupunba"              ,x=g.s)
   g.s = sub("Pterocarpus amazonium"      ,"Pterocarpus santalinoides"    ,x=g.s)
   g.s = sub("Pterocarpus rhoire"         ,"Pterocarpus rohrii"           ,x=g.s)
   g.s = sub("Pterocarpus rhoiri"         ,"Pterocarpus rohrii"           ,x=g.s)
   g.s = sub("Ragala guianensis"          ,"Chrysophyllum sanguinolentum" ,x=g.s)
   g.s = sub("Rapanea ferruginea"         ,"Myrsine coriacea"             ,x=g.s)
   g.s = sub("Rapanea guianensis"         ,"Myrsine guianensis"           ,x=g.s)
   g.s = sub("Rheedia acuminata"          ,"Garcinia madruno"             ,x=g.s)
   g.s = sub("Rim de"                     ,"Crudia"                       ,x=g.s)
   g.s = sub("Rinorea pectino-squamata"   ,"Rinorea pectinosquamata"      ,x=g.s)
   g.s = sub("Rinoria guianensis"         ,"Rinorea guianensis"           ,x=g.s)
   g.s = sub("Rinoria racenosa"           ,"Rinorea racemosa"             ,x=g.s)
   g.s = sub("Rolinha esxuccar"           ,"Rollinia exsucca"             ,x=g.s)
   g.s = sub("Rollinia esxucca"           ,"Rollinia exsucca"             ,x=g.s)
   g.s = sub("Sacrogrotis guianensis"     ,"Sacoglottis guianensis"       ,x=g.s)
   g.s = sub("Sagotia sp."                ,"Sagotia"                      ,x=g.s)
   g.s = sub("Salacea"                    ,"Salacia"                      ,x=g.s)
   g.s = sub("Salacia imprissifolia"      ,"Salacia impressifolia"        ,x=g.s)
   g.s = sub("Sapium sp."                 ,"Sapium"                       ,x=g.s)
   g.s = sub("Sclerolobium chrysophylum"  ,"Tachigali chrysophylla"       ,x=g.s)
   g.s = sub("Sclerolobium chrysophyllum" ,"Tachigali chrysophylla"       ,x=g.s)
   g.s = sub("Sclerolobium guianensis"    ,"Tachigali guianensis"         ,x=g.s)
   g.s = sub("Sclerolobium guianense"     ,"Tachigali guianensis"         ,x=g.s)
   g.s = sub("Sclerolobium sp."           ,"Sclerolobium"                 ,x=g.s)
   g.s = sub("Simaba guianenses"          ,"Simaba guianensis"            ,x=g.s)
   g.s = sub("Simabacedron planch."       ,"Simaba cedron"                ,x=g.s)
   g.s = sub("Simarouba armara"           ,"Simarouba amara"              ,x=g.s)
   g.s = sub("Simaruba armara"            ,"Simarouba amara"              ,x=g.s)
   g.s = sub("Sipararuna cristata"        ,"Siparuna cristata"            ,x=g.s)
   g.s = sub("Sipararuna decipiens"       ,"Siparuna decipiens"           ,x=g.s)
   g.s = sub("Stryphnodendron pulchrrimum","Stryphnodendron pulcherrimum" ,x=g.s)
   g.s = sub("Stryphnodendron sp."        ,"Stryphnodendron"              ,x=g.s)
   g.s = sub("Styrax ferrugineum"         ,"Styrax ferrugineus"           ,x=g.s)
   g.s = sub("Swartzia arborensis"        ,"Swartzia arborescens"         ,x=g.s)
   g.s = sub("Swartzia flamingii"         ,"Swartzia flaemingii"          ,x=g.s)
   g.s = sub("Swartizia microcarpum"      ,"Swartzia microcarpa"          ,x=g.s)
   g.s = sub("Swartzia retusa"            ,"Swartzia recurva"             ,x=g.s)
   g.s = sub("Swartzia sp."               ,"Swartzia"                     ,x=g.s)
   g.s = sub("Tachigalia alba"            ,"Tachigali alba"               ,x=g.s)
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
   g.s = sub("Thysodium sp."              ,"Thyrsodium"                   ,x=g.s)
   g.s = sub("Trattinickia laurencei"     ,"Trattinnickia lawrancei"      ,x=g.s)
   g.s = sub("Trattinickia laurense"      ,"Trattinnickia lawrancei"      ,x=g.s)
   g.s = sub("Trattinickia rhoifolia"     ,"Trattinnickia rhoifolia"      ,x=g.s)
   g.s = sub("Trichilia lequente"         ,"Trichilia lecointei"          ,x=g.s)
   g.s = sub("Trichilia sp."              ,"Trichilia"                    ,x=g.s)
   g.s = sub("Triquilia lequente"         ,"Trichilia lecointei"          ,x=g.s)
   g.s = sub("Trymatococcus paraensis"    ,"Trymatococcus amazonicus"     ,x=g.s)
   g.s = sub("Vataicreopesis speciosa"    ,"Vataireopsis speciosa"        ,x=g.s)
   g.s = sub("Vataireopesis speciosa"     ,"Vataireopsis speciosa"        ,x=g.s)
   g.s = sub("Vernonia diffusa"           ,"Vernonia erigeroides"         ,x=g.s)
   g.s = sub("Vimia guianensis"           ,"Vismia guianensis"            ,x=g.s)
   g.s = sub("Virola crebinervia"         ,"Virola crebrinervia"          ,x=g.s)
   g.s = sub("Virola melinonii"           ,"Virola michelii"              ,x=g.s)
   g.s = sub("Virola melionii"            ,"Virola michelii"              ,x=g.s)
   g.s = sub("Virola michelli"            ,"Virola michelii"              ,x=g.s)
   g.s = sub("Vochisia surinamensis"      ,"Vochysia surinamensis"        ,x=g.s)
   #---------------------------------------------------------------------------------------#







   #----- Break again into genus and species. ---------------------------------------------#
   gs.mat      = t(as.data.frame(sapply(X=g.s,FUN=strsplit,split=" ")))
   del         = gs.mat == "na"
   gs.mat[del] = NA
   g           = capwords(gs.mat[,1],strict=TRUE)
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #      Copy the data back to the structure.                                             #
   #---------------------------------------------------------------------------------------#
   if ("scientific" %in% names(dat)) dat$scientific = g.s
   if ("genus"      %in% names(dat)) dat$genus      = g
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
   datum$family = capwords(datum$family,strict=TRUE)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Look-up table for all families.                                                   #
   #---------------------------------------------------------------------------------------#
   g2f        = list()

   g2f[[  1]] = list( genus = "Abarema"            , family = "Fabaceae"          )
   g2f[[  2]] = list( genus = "Acmanthera"         , family = "Malpighiaceae"     )
   g2f[[  3]] = list( genus = "Acrocomia"          , family = "Arecaceae"         )
   g2f[[  4]] = list( genus = "Acosmium"           , family = "Fabaceae"          )
   g2f[[  5]] = list( genus = "Adenophaedra"       , family = "Euphorbiaceae"     )
   g2f[[  6]] = list( genus = "Agonandra"          , family = "Opiliaceae"        )
   g2f[[  7]] = list( genus = "Aiouea"             , family = "Lauraceae"         )
   g2f[[  8]] = list( genus = "Albizia"            , family = "Fabaceae"          )
   g2f[[  9]] = list( genus = "Alchorneopsis"      , family = "Euphorbiaceae"     )
   g2f[[ 10]] = list( genus = "Aldina"             , family = "Fabaceae"          )
   g2f[[ 11]] = list( genus = "Alexa"              , family = "Fabaceae"          )
   g2f[[ 12]] = list( genus = "Allantoma"          , family = "Lecythidaceae"     )
   g2f[[ 13]] = list( genus = "Allophylus"         , family = "Sapindaceae"       )
   g2f[[ 14]] = list( genus = "Amaioua"            , family = "Rubiaceae"         )
   g2f[[ 15]] = list( genus = "Ambelania"          , family = "Apocynaceae"       )
   g2f[[ 16]] = list( genus = "Amburana"           , family = "Fabaceae"          )
   g2f[[ 17]] = list( genus = "Ampelocera"         , family = "Ulmaceae"          )
   g2f[[ 18]] = list( genus = "Amphirrhox"         , family = "Violaceae"         )
   g2f[[ 19]] = list( genus = "Anacardium"         , family = "Anacardiaceae"     )
   g2f[[ 20]] = list( genus = "Anadenanthera"      , family = "Fabaceae"          )
   g2f[[ 21]] = list( genus = "Anaxagorea"         , family = "Annonaceae"        )
   g2f[[ 22]] = list( genus = "Andira"             , family = "Fabaceae"          )
   g2f[[ 23]] = list( genus = "Aniba"              , family = "Lauraceae"         )
   g2f[[ 24]] = list( genus = "Anisophyllea"       , family = "Anisophylleaceae"  )
   g2f[[ 25]] = list( genus = "Annona"             , family = "Annonaceae"        )
   g2f[[ 26]] = list( genus = "Anomalocalyx"       , family = "Euphorbiaceae"     )
   g2f[[ 27]] = list( genus = "Antonia"            , family = "Loganiaceae"       )
   g2f[[ 28]] = list( genus = "Aparisthmium"       , family = "Euphorbiaceae"     )
   g2f[[ 29]] = list( genus = "Apeiba"             , family = "Malvaceae"         )
   g2f[[ 30]] = list( genus = "Aptandra"           , family = "Olacaceae"         )
   g2f[[ 31]] = list( genus = "Apuleia"            , family = "Fabaceae"          )
   g2f[[ 32]] = list( genus = "Aspidosperma"       , family = "Apocynaceae"       )
   g2f[[ 33]] = list( genus = "Astrocaryum"        , family = "Arecaceae"         )
   g2f[[ 34]] = list( genus = "Astronium"          , family = "Anacardiaceae"     )
   g2f[[ 35]] = list( genus = "Attalea"            , family = "Arecaceae"         )
   g2f[[ 36]] = list( genus = "Austromyrtus"       , family = "Myrtaceae"         )
   g2f[[ 37]] = list( genus = "Bactris"            , family = "Arecaceae"         )
   g2f[[ 38]] = list( genus = "Bagassa"            , family = "Moraceae"          )
   g2f[[ 39]] = list( genus = "Balizia"            , family = "Fabaceae"          )
   g2f[[ 40]] = list( genus = "Batesia"            , family = "Fabaceae"          )
   g2f[[ 41]] = list( genus = "Batocarpus"         , family = "Moraceae"          )
   g2f[[ 42]] = list( genus = "Bauhinia"           , family = "Fabaceae"          )
   g2f[[ 43]] = list( genus = "Bellucia"           , family = "Melastomataceae"   )
   g2f[[ 44]] = list( genus = "Bertholletia"       , family = "Lecythidaceae"     )
   g2f[[ 45]] = list( genus = "Blepharocalyx"      , family = "Myrtaceae"         )
   g2f[[ 46]] = list( genus = "Bixa"               , family = "Bixaceae"          )
   g2f[[ 47]] = list( genus = "Bocageopsis"        , family = "Annonaceae"        )
   g2f[[ 48]] = list( genus = "Bocoa"              , family = "Fabaceae"          )
   g2f[[ 49]] = list( genus = "Bombax"             , family = "Malvaceae"         )
   g2f[[ 50]] = list( genus = "Botryarrhena"       , family = "Rubiaceae"         )
   g2f[[ 51]] = list( genus = "Bowdichia"          , family = "Fabaceae"          )
   g2f[[ 52]] = list( genus = "Brachychiton"       , family = "Malvaceae"         )
   g2f[[ 53]] = list( genus = "Brosimum"           , family = "Moraceae"          )
   g2f[[ 54]] = list( genus = "Buchenavia"         , family = "Combretaceae"      )
   g2f[[ 55]] = list( genus = "Byrsonima"          , family = "Malpighiaceae"     )
   g2f[[ 56]] = list( genus = "Caesalpinia"        , family = "Fabaceae"          )
   g2f[[ 57]] = list( genus = "Calliandra"         , family = "Fabaceae"          )
   g2f[[ 58]] = list( genus = "Calophyllum"        , family = "Calophyllaceae"    )
   g2f[[ 59]] = list( genus = "Calyptranthes"      , family = "Myrtaceae"         )
   g2f[[ 60]] = list( genus = "Campomanesia"       , family = "Myrtaceae"         )
   g2f[[ 61]] = list( genus = "Capirona"           , family = "Rubiaceae"         )
   g2f[[ 62]] = list( genus = "Capparidastrum"     , family = "Capparaceae"       )
   g2f[[ 63]] = list( genus = "Capparis"           , family = "Capparaceae"       )
   g2f[[ 64]] = list( genus = "Caraipa"            , family = "Calophyllaceae"    )
   g2f[[ 65]] = list( genus = "Carapa"             , family = "Meliaceae"         )
   g2f[[ 66]] = list( genus = "Cariniana"          , family = "Lecythidaceae"     )
   g2f[[ 67]] = list( genus = "Caryocar"           , family = "Caryocaraceae"     )
   g2f[[ 68]] = list( genus = "Casearia"           , family = "Salicaceae"        )
   g2f[[ 69]] = list( genus = "Cassia"             , family = "Fabaceae"          )
   g2f[[ 70]] = list( genus = "Castilla"           , family = "Moraceae"          )
   g2f[[ 71]] = list( genus = "Catostemma"         , family = "Malvaceae"         )
   g2f[[ 72]] = list( genus = "Cecropia"           , family = "Urticaceae"        )
   g2f[[ 73]] = list( genus = "Cedrela"            , family = "Meliaceae"         )
   g2f[[ 74]] = list( genus = "Cedrelinga"         , family = "Fabaceae"          )
   g2f[[ 75]] = list( genus = "Ceiba"              , family = "Malvaceae"         )
   g2f[[ 76]] = list( genus = "Cereus"             , family = "Cactaceae"         )
   g2f[[ 77]] = list( genus = "Chaetocarpus"       , family = "Euphorbiaceae"     )
   g2f[[ 78]] = list( genus = "Chamaecrista"       , family = "Fabaceae"          )
   g2f[[ 79]] = list( genus = "Chaunochiton"       , family = "Olacaceae"         )
   g2f[[ 80]] = list( genus = "Cheiloclinium"      , family = "Celastraceae"      )
   g2f[[ 81]] = list( genus = "Chimarrhis"         , family = "Rubiaceae"         )
   g2f[[ 82]] = list( genus = "Chromolucuma"       , family = "Sapotaceae"        )
   g2f[[ 83]] = list( genus = "Chrysobalanus"      , family = "Chrysobalanaceae"  )
   g2f[[ 84]] = list( genus = "Chrysophyllum"      , family = "Sapotaceae"        )
   g2f[[ 85]] = list( genus = "Clarisia"           , family = "Moraceae"          )
   g2f[[ 86]] = list( genus = "Clusia"             , family = "Clusiaceae"        )
   g2f[[ 87]] = list( genus = "Cnidoscolus"        , family = "Euphorbiaceae"     )
   g2f[[ 88]] = list( genus = "Coccoloba"          , family = "Polygonaceae"      )
   g2f[[ 89]] = list( genus = "Conceveiba"         , family = "Euphorbiaceae"     )
   g2f[[ 90]] = list( genus = "Connarus"           , family = "Connaraceae"       )
   g2f[[ 91]] = list( genus = "Copaifera"          , family = "Fabaceae"          )
   g2f[[ 92]] = list( genus = "Cordia"             , family = "Boraginaceae"      )
   g2f[[ 93]] = list( genus = "Corythophora"       , family = "Lecythidaceae"     )
   g2f[[ 94]] = list( genus = "Couepia"            , family = "Chrysobalanaceae"  )
   g2f[[ 95]] = list( genus = "Couma"              , family = "Apocynaceae"       )
   g2f[[ 96]] = list( genus = "Couratari"          , family = "Lecythidaceae"     )
   g2f[[ 97]] = list( genus = "Coussapoa"          , family = "Urticaceae"        )
   g2f[[ 98]] = list( genus = "Coussarea"          , family = "Rubiaceae"         )
   g2f[[ 99]] = list( genus = "Coutarea"           , family = "Rubiaceae"         )
   g2f[[100]] = list( genus = "Crepidospermum"     , family = "Burseraceae"       )
   g2f[[101]] = list( genus = "Croton"             , family = "Euphorbiaceae"     )
   g2f[[102]] = list( genus = "Crudia"             , family = "Fabaceae"          )
   g2f[[103]] = list( genus = "Cupania"            , family = "Sapindaceae"       )
   g2f[[104]] = list( genus = "Cybianthus"         , family = "Primulaceae"       )
   g2f[[105]] = list( genus = "Cybistax"           , family = "Bignoniaceae"      )
   g2f[[106]] = list( genus = "Dacryodes"          , family = "Burseraceae"       )
   g2f[[107]] = list( genus = "Dalbergia"          , family = "Fabaceae"          )
   g2f[[108]] = list( genus = "Dendrobangia"       , family = "Cardiopteridaceae" )
   g2f[[109]] = list( genus = "Dialium"            , family = "Fabaceae"          )
   g2f[[110]] = list( genus = "Diclinanona"        , family = "Annonaceae"        )
   g2f[[111]] = list( genus = "Dicorynia"          , family = "Fabaceae"          )
   g2f[[112]] = list( genus = "Dicranostyles"      , family = "Convolvulaceae"    )
   g2f[[113]] = list( genus = "Dicypellium"        , family = "Lauraceae"         )
   g2f[[114]] = list( genus = "Dimorphandra"       , family = "Fabaceae"          )
   g2f[[115]] = list( genus = "Dinizia"            , family = "Fabaceae"          )
   g2f[[116]] = list( genus = "Diospyros"          , family = "Ebenaceae"         )
   g2f[[117]] = list( genus = "Diplotropis"        , family = "Fabaceae"          )
   g2f[[118]] = list( genus = "Dipteryx"           , family = "Fabaceae"          )
   g2f[[119]] = list( genus = "Discophora"         , family = "Stemonuraceae"     )
   g2f[[120]] = list( genus = "Doliocarpus"        , family = "Dilleniaceae"      )
   g2f[[121]] = list( genus = "Drypetes"           , family = "Putranjivaceae"    )
   g2f[[122]] = list( genus = "Duckeodendron"      , family = "Solanaceae"        )
   g2f[[123]] = list( genus = "Duckesia"           , family = "Humiriaceae"       )
   g2f[[124]] = list( genus = "Duguetia"           , family = "Annonaceae"        )
   g2f[[125]] = list( genus = "Duroia"             , family = "Rubiaceae"         )
   g2f[[126]] = list( genus = "Dystovomita"        , family = "Clusiaceae"        )
   g2f[[127]] = list( genus = "Ecclinusa"          , family = "Sapotaceae"        )
   g2f[[128]] = list( genus = "Elaeoluma"          , family = "Sapotaceae"        )
   g2f[[129]] = list( genus = "Elizabetha"         , family = "Fabaceae"          )
   g2f[[130]] = list( genus = "Emmotum"            , family = "Emmotaceae"        )
   g2f[[131]] = list( genus = "Endlicheria"        , family = "Lauraceae"         )
   g2f[[132]] = list( genus = "Endopleura"         , family = "Humiriaceae"       )
   g2f[[133]] = list( genus = "Enterolobium"       , family = "Fabaceae"          )
   g2f[[134]] = list( genus = "Eperua"             , family = "Fabaceae"          )
   g2f[[135]] = list( genus = "Ephedranthus"       , family = "Annonaceae"        )
   g2f[[136]] = list( genus = "Eriotheca"          , family = "Malvaceae"         )
   g2f[[137]] = list( genus = "Erisma"             , family = "Vochysiaceae"      )
   g2f[[138]] = list( genus = "Erythroxylum"       , family = "Erythroxylaceae"   )
   g2f[[139]] = list( genus = "Eschweilera"        , family = "Lecythidaceae"     )
   g2f[[140]] = list( genus = "Eugenia"            , family = "Myrtaceae"         )
   g2f[[141]] = list( genus = "Euterpe"            , family = "Arecaceae"         )
   g2f[[142]] = list( genus = "Euxylophora"        , family = "Rutaceae"          )
   g2f[[143]] = list( genus = "Faramea"            , family = "Rubiaceae"         )
   g2f[[144]] = list( genus = "Ferdinandusa"       , family = "Rubiaceae"         )
   g2f[[145]] = list( genus = "Ficus"              , family = "Moraceae"          )
   g2f[[146]] = list( genus = "Fusaea"             , family = "Annonaceae"        )
   g2f[[147]] = list( genus = "Garcinia"           , family = "Clusiaceae"        )
   g2f[[148]] = list( genus = "Geissospermum"      , family = "Apocynaceae"       )
   g2f[[149]] = list( genus = "Genipa"             , family = "Rubiaceae"         )
   g2f[[150]] = list( genus = "Gilbertiodendron"   , family = "Fabaceae"          )
   g2f[[151]] = list( genus = "Glycydendron"       , family = "Euphorbiaceae"     )
   g2f[[152]] = list( genus = "Goupia"             , family = "Goupiaceae"        )
   g2f[[153]] = list( genus = "Guapira"            , family = "Nyctaginaceae"     )
   g2f[[154]] = list( genus = "Guarea"             , family = "Meliaceae"         )
   g2f[[155]] = list( genus = "Guatteria"          , family = "Annonaceae"        )
   g2f[[156]] = list( genus = "Guazuma"            , family = "Malvaceae"         )
   g2f[[157]] = list( genus = "Gustavia"           , family = "Lecythidaceae"     )
   g2f[[158]] = list( genus = "Hebepetalum"        , family = "Linaceae"          )
   g2f[[159]] = list( genus = "Heisteria"          , family = "Olacaceae"         )
   g2f[[160]] = list( genus = "Helianthostylis"    , family = "Moraceae"          )
   g2f[[161]] = list( genus = "Helicostylis"       , family = "Moraceae"          )
   g2f[[162]] = list( genus = "Henriettea"         , family = "Melastomataceae"   )
   g2f[[163]] = list( genus = "Henriettella"       , family = "Melastomataceae"   )
   g2f[[164]] = list( genus = "Hevea"              , family = "Euphorbiaceae"     )
   g2f[[165]] = list( genus = "Himatanthus"        , family = "Apocynaceae"       )
   g2f[[166]] = list( genus = "Hippocratea"        , family = "Celastraceae"      )
   g2f[[167]] = list( genus = "Hirtella"           , family = "Chrysobalanaceae"  )
   g2f[[168]] = list( genus = "Humiria"            , family = "Humiriaceae"       )
   g2f[[169]] = list( genus = "Humiriastrum"       , family = "Humiriaceae"       )
   g2f[[170]] = list( genus = "Hura"               , family = "Euphorbiaceae"     )
   g2f[[171]] = list( genus = "Hymenaea"           , family = "Fabaceae"          )
   g2f[[172]] = list( genus = "Hymenolobium"       , family = "Fabaceae"          )
   g2f[[173]] = list( genus = "Ignotum"            , family = "Ignotaceae"        )
   g2f[[174]] = list( genus = "Ilex"               , family = "Aquifoliaceae"     )
   g2f[[175]] = list( genus = "Inga"               , family = "Fabaceae"          )
   g2f[[176]] = list( genus = "Iryanthera"         , family = "Myristicaceae"     )
   g2f[[177]] = list( genus = "Isertia"            , family = "Rubiaceae"         )
   g2f[[178]] = list( genus = "Jacaranda"          , family = "Bignoniaceae"      )
   g2f[[179]] = list( genus = "Jacaratia"          , family = "Caricaceae"        )
   g2f[[180]] = list( genus = "Jatropha"           , family = "Euphorbiaceae"     )
   g2f[[181]] = list( genus = "Joannesia"          , family = "Euphorbiaceae"     )
   g2f[[182]] = list( genus = "Justicia"           , family = "Acanthaceae"       )
   g2f[[183]] = list( genus = "Lacistema"          , family = "Lacistemataceae"   )
   g2f[[184]] = list( genus = "Lacmellea"          , family = "Apocynaceae"       )
   g2f[[185]] = list( genus = "Lacunaria"          , family = "Ochnaceae"         )
   g2f[[186]] = list( genus = "Ladenbergia"        , family = "Rubiaceae"         )
   g2f[[187]] = list( genus = "Laetia"             , family = "Salicaceae"        )
   g2f[[188]] = list( genus = "Lecythis"           , family = "Lecythidaceae"     )
   g2f[[189]] = list( genus = "Leonia"             , family = "Violaceae"         )
   g2f[[190]] = list( genus = "Liana"              , family = "Lianaceae"         )
   g2f[[191]] = list( genus = "Licania"            , family = "Chrysobalanaceae"  )
   g2f[[192]] = list( genus = "Licaria"            , family = "Lauraceae"         )
   g2f[[193]] = list( genus = "Lindackeria"        , family = "Achariaceae"       )
   g2f[[194]] = list( genus = "Lippia"             , family = "Verbenaceae"       )
   g2f[[195]] = list( genus = "Lonchocarpus"       , family = "Fabaceae"          )
   g2f[[196]] = list( genus = "Luehea"             , family = "Malvaceae"         )
   g2f[[197]] = list( genus = "Lueheopsis"         , family = "Malvaceae"         )
   g2f[[198]] = list( genus = "Mabea"              , family = "Euphorbiaceae"     )
   g2f[[199]] = list( genus = "Machaerium"         , family = "Fabaceae"          )
   g2f[[200]] = list( genus = "Macoubea"           , family = "Apocynaceae"       )
   g2f[[201]] = list( genus = "Macrolobium"        , family = "Fabaceae"          )
   g2f[[202]] = list( genus = "Malouetia"          , family = "Apocynaceae"       )
   g2f[[203]] = list( genus = "Manicaria"          , family = "Arecaceae"         )
   g2f[[204]] = list( genus = "Manihot"            , family = "Euphorbiaceae"     )
   g2f[[205]] = list( genus = "Manilkara"          , family = "Sapotaceae"        )
   g2f[[206]] = list( genus = "Maquira"            , family = "Moraceae"          )
   g2f[[207]] = list( genus = "Marcgravia"         , family = "Marcgraviaceae"    )
   g2f[[208]] = list( genus = "Marlierea"          , family = "Myrtaceae"         )
   g2f[[209]] = list( genus = "Matayba"            , family = "Sapindaceae"       )
   g2f[[210]] = list( genus = "Matisia"            , family = "Malvaceae"         )
   g2f[[211]] = list( genus = "Mauritia"           , family = "Arecaceae"         )
   g2f[[212]] = list( genus = "Mauritiella"        , family = "Arecaceae"         )
   g2f[[213]] = list( genus = "Maytenus"           , family = "Celastraceae"      )
   g2f[[214]] = list( genus = "Melicoccus"         , family = "Sapindaceae"       )
   g2f[[215]] = list( genus = "Mezilaurus"         , family = "Lauraceae"         )
   g2f[[216]] = list( genus = "Miconia"            , family = "Melastomataceae"   )
   g2f[[217]] = list( genus = "Micrandra"          , family = "Euphorbiaceae"     )
   g2f[[218]] = list( genus = "Micrandropsis"      , family = "Euphorbiaceae"     )
   g2f[[219]] = list( genus = "Micropholis"        , family = "Sapotaceae"        )
   g2f[[220]] = list( genus = "Mikania"            , family = "Asteraceae"        )
   g2f[[221]] = list( genus = "Mimosa"             , family = "Fabaceae"          )
   g2f[[222]] = list( genus = "Minquartia"         , family = "Olacaceae"         )
   g2f[[223]] = list( genus = "Misanteca"          , family = "Lauraceae"         )
   g2f[[224]] = list( genus = "Monopteryx"         , family = "Fabaceae"          )
   g2f[[225]] = list( genus = "Moronobea"          , family = "Clusiaceae"        )
   g2f[[226]] = list( genus = "Mouriri"            , family = "Melastomataceae"   )
   g2f[[227]] = list( genus = "Myrocarpus"         , family = "Fabaceae"          )
   g2f[[228]] = list( genus = "Myrcia"             , family = "Myrtaceae"         )
   g2f[[229]] = list( genus = "Myrciaria"          , family = "Myrtaceae"         )
   g2f[[230]] = list( genus = "Myrsine"            , family = "Primulaceae"       )
   g2f[[231]] = list( genus = "Naucleopsis"        , family = "Moraceae"          )
   g2f[[232]] = list( genus = "Nectandra"          , family = "Lauraceae"         )
   g2f[[233]] = list( genus = "Neea"               , family = "Nyctaginaceae"     )
   g2f[[234]] = list( genus = "Ocotea"             , family = "Lauraceae"         )
   g2f[[235]] = list( genus = "Oenocarpus"         , family = "Arecaceae"         )
   g2f[[236]] = list( genus = "Onychopetalum"      , family = "Annonaceae"        )
   g2f[[237]] = list( genus = "Ormosia"            , family = "Fabaceae"          )
   g2f[[238]] = list( genus = "Osteophloeum"       , family = "Myristicaceae"     )
   g2f[[239]] = list( genus = "Ouratea"            , family = "Ochnaceae"         )
   g2f[[240]] = list( genus = "Oxandra"            , family = "Annonaceae"        )
   g2f[[241]] = list( genus = "Pachira"            , family = "Malvaceae"         )
   g2f[[242]] = list( genus = "Palicourea"         , family = "Rubiaceae"         )
   g2f[[243]] = list( genus = "Parahancornia"      , family = "Apocynaceae"       )
   g2f[[244]] = list( genus = "Paraia"             , family = "Lauraceae"         )
   g2f[[245]] = list( genus = "Parinari"           , family = "Chrysobalanaceae"  )
   g2f[[246]] = list( genus = "Parkia"             , family = "Fabaceae"          )
   g2f[[247]] = list( genus = "Paullinia"          , family = "Sapindaceae"       )
   g2f[[248]] = list( genus = "Pausandra"          , family = "Euphorbiaceae"     )
   g2f[[249]] = list( genus = "Paypayrola"         , family = "Violaceae"         )
   g2f[[250]] = list( genus = "Peltogyne"          , family = "Fabaceae"          )
   g2f[[251]] = list( genus = "Pera"               , family = "Euphorbiaceae"     )
   g2f[[252]] = list( genus = "Perebea"            , family = "Moraceae"          )
   g2f[[253]] = list( genus = "Peridiscus"         , family = "Peridiscaceae"     )
   g2f[[254]] = list( genus = "Phenakospermum"     , family = "Strelitziaceae"    )
   g2f[[255]] = list( genus = "Piptadenia"         , family = "Fabaceae"          )
   g2f[[256]] = list( genus = "Platonia"           , family = "Clusiaceae"        )
   g2f[[257]] = list( genus = "Platymiscium"       , family = "Fabaceae"          )
   g2f[[258]] = list( genus = "Plenckia"           , family = "Celastraceae"      )
   g2f[[259]] = list( genus = "Poecilanthe"        , family = "Fabaceae"          )
   g2f[[260]] = list( genus = "Pogonophora"        , family = "Euphorbiaceae"     )
   g2f[[261]] = list( genus = "Poraqueiba"         , family = "Icacinaceae"       )
   g2f[[262]] = list( genus = "Posoqueria"         , family = "Rubiaceae"         )
   g2f[[263]] = list( genus = "Pourouma"           , family = "Urticaceae"        )
   g2f[[264]] = list( genus = "Pouteria"           , family = "Sapotaceae"        )
   g2f[[265]] = list( genus = "Pradosia"           , family = "Sapotaceae"        )
   g2f[[266]] = list( genus = "Protium"            , family = "Burseraceae"       )
   g2f[[267]] = list( genus = "Prunus"             , family = "Rosaceae"          )
   g2f[[268]] = list( genus = "Pseudolmedia"       , family = "Moraceae"          )
   g2f[[269]] = list( genus = "Pseudopiptadenia"   , family = "Fabaceae"          )
   g2f[[270]] = list( genus = "Pseudoxandra"       , family = "Annonaceae"        )
   g2f[[271]] = list( genus = "Psidium"            , family = "Myrtaceae"         )
   g2f[[272]] = list( genus = "Pterandra"          , family = "Malpighiaceae"     )
   g2f[[273]] = list( genus = "Pterocarpus"        , family = "Fabaceae"          )
   g2f[[274]] = list( genus = "Ptychopetalum"      , family = "Olacaceae"         )
   g2f[[275]] = list( genus = "Qualea"             , family = "Vochysiaceae"      )
   g2f[[276]] = list( genus = "Quararibea"         , family = "Malvaceae"         )
   g2f[[277]] = list( genus = "Quiina"             , family = "Ochnaceae"         )
   g2f[[278]] = list( genus = "Raputia"            , family = "Rutaceae"          )
   g2f[[279]] = list( genus = "Rauvolfia"          , family = "Apocynaceae"       )
   g2f[[280]] = list( genus = "Recordoxylon"       , family = "Fabaceae"          )
   g2f[[281]] = list( genus = "Rhodostemonodaphne" , family = "Lauraceae"         )
   g2f[[282]] = list( genus = "Richeria"           , family = "Phyllanthaceae"    )
   g2f[[283]] = list( genus = "Rinorea"            , family = "Violaceae"         )
   g2f[[284]] = list( genus = "Rollinia"           , family = "Annonaceae"        )
   g2f[[285]] = list( genus = "Roucheria"          , family = "Linaceae"          )
   g2f[[286]] = list( genus = "Ruizterania"        , family = "Vochysiaceae"      )
   g2f[[287]] = list( genus = "Ryania"             , family = "Salicaceae"        )
   g2f[[288]] = list( genus = "Sacoglottis"        , family = "Humiriaceae"       )
   g2f[[289]] = list( genus = "Sagotia"            , family = "Euphorbiaceae"     )
   g2f[[290]] = list( genus = "Salacia"            , family = "Celastraceae"      )
   g2f[[291]] = list( genus = "Sapium"             , family = "Euphorbiaceae"     )
   g2f[[292]] = list( genus = "Sarcaulus"          , family = "Sapotaceae"        )
   g2f[[293]] = list( genus = "Schefflera"         , family = "Araliaceae"        )
   g2f[[294]] = list( genus = "Schinopsis"         , family = "Anacardiaceae"     )
   g2f[[295]] = list( genus = "Schizolobium"       , family = "Fabaceae"          )
   g2f[[296]] = list( genus = "Sclerolobium"       , family = "Fabaceae"          )
   g2f[[297]] = list( genus = "Scleronema"         , family = "Malvaceae"         )
   g2f[[298]] = list( genus = "Senefeldera"        , family = "Euphorbiaceae"     )
   g2f[[299]] = list( genus = "Senna"              , family = "Fabaceae"          )
   g2f[[300]] = list( genus = "Sextonia"           , family = "Lauraceae"         )
   g2f[[301]] = list( genus = "Simaba"             , family = "Simaroubaceae"     )
   g2f[[302]] = list( genus = "Simarouba"          , family = "Simaroubaceae"     )
   g2f[[303]] = list( genus = "Siparuna"           , family = "Siparunaceae"      )
   g2f[[304]] = list( genus = "Sloanea"            , family = "Elaeocarpaceae"    )
   g2f[[305]] = list( genus = "Socratea"           , family = "Arecaceae"         )
   g2f[[306]] = list( genus = "Sorocea"            , family = "Moraceae"          )
   g2f[[307]] = list( genus = "Spondias"           , family = "Anacardiaceae"     )
   g2f[[308]] = list( genus = "Sterculia"          , family = "Malvaceae"         )
   g2f[[309]] = list( genus = "Strychnos"          , family = "Loganiaceae"       )
   g2f[[310]] = list( genus = "Stryphnodendron"    , family = "Fabaceae"          )
   g2f[[311]] = list( genus = "Styrax"             , family = "Styracaceae"       )
   g2f[[312]] = list( genus = "Swartzia"           , family = "Fabaceae"          )
   g2f[[313]] = list( genus = "Swietenia"          , family = "Meliaceae"         )
   g2f[[314]] = list( genus = "Symphonia"          , family = "Clusiaceae"        )
   g2f[[315]] = list( genus = "Syzygium"           , family = "Myrtaceae"         )
   g2f[[316]] = list( genus = "Tabebuia"           , family = "Bignoniaceae"      )
   g2f[[317]] = list( genus = "Tabernaemontana"    , family = "Apocynaceae"       )
   g2f[[318]] = list( genus = "Tachigali"          , family = "Fabaceae"          )
   g2f[[319]] = list( genus = "Talisia"            , family = "Sapindaceae"       )
   g2f[[320]] = list( genus = "Tapirira"           , family = "Anacardiaceae"     )
   g2f[[321]] = list( genus = "Tapura"             , family = "Dichapetalaceae"   )
   g2f[[322]] = list( genus = "Taralea"            , family = "Fabaceae"          )
   g2f[[323]] = list( genus = "Terminalia"         , family = "Combretaceae"      )
   g2f[[324]] = list( genus = "Tetragastris"       , family = "Burseraceae"       )
   g2f[[325]] = list( genus = "Theobroma"          , family = "Malvaceae"         )
   g2f[[326]] = list( genus = "Thyrsodium"         , family = "Anacardiaceae"     )
   g2f[[327]] = list( genus = "Tocoyena"           , family = "Rubiaceae"         )
   g2f[[328]] = list( genus = "Toulicia"           , family = "Sapindaceae"       )
   g2f[[329]] = list( genus = "Touroulia"          , family = "Ochnaceae"         )
   g2f[[330]] = list( genus = "Tovomita"           , family = "Clusiaceae"        )
   g2f[[331]] = list( genus = "Trattinnickia"      , family = "Burseraceae"       )
   g2f[[332]] = list( genus = "Trema"              , family = "Cannabaceae"       )
   g2f[[333]] = list( genus = "Trichilia"          , family = "Meliaceae"         )
   g2f[[334]] = list( genus = "Trymatococcus"      , family = "Moraceae"          )
   g2f[[335]] = list( genus = "Unonopsis"          , family = "Annonaceae"        )
   g2f[[336]] = list( genus = "Vantanea"           , family = "Humiriaceae"       )
   g2f[[337]] = list( genus = "Vatairea"           , family = "Fabaceae"          )
   g2f[[338]] = list( genus = "Vataireopsis"       , family = "Fabaceae"          )
   g2f[[339]] = list( genus = "Vernonia"           , family = "Asteraceae"        )
   g2f[[340]] = list( genus = "Virola"             , family = "Myristicaceae"     )
   g2f[[341]] = list( genus = "Vismia"             , family = "Hypericaceae"      )
   g2f[[342]] = list( genus = "Vitex"              , family = "Lamiaceae"         )
   g2f[[343]] = list( genus = "Vochysia"           , family = "Vochysiaceae"      )
   g2f[[344]] = list( genus = "Votomita"           , family = "Melastomataceae"   )
   g2f[[345]] = list( genus = "Vouacapoua"         , family = "Fabaceae"          )
   g2f[[346]] = list( genus = "Warszewiczia"       , family = "Rubiaceae"         )
   g2f[[347]] = list( genus = "Xylopia"            , family = "Annonaceae"        )
   g2f[[348]] = list( genus = "Zanthoxylum"        , family = "Rutaceae"          )
   g2f[[349]] = list( genus = "Zygia"              , family = "Fabaceae"          )


   #---------------------------------------------------------------------------------------#


   #----- Convert the g2f list into a data frame. -----------------------------------------#
   g2f = data.frame( apply( X      = t(sapply(X = g2f , FUN = c))
                          , MARGIN = c(1,2)
                          , FUN    = unlist
                          )#end apply
                   , stringsAsFactors = FALSE
                   )#end data.frame
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Check whether the family has all genera needed.                                   #
   #---------------------------------------------------------------------------------------#
   idx = which(! is.na(datum$genus) & ! datum$genus %in% g2f$genus)
   if (length(idx) > 0){
     cat(" You must add genera to g2f...","\n")
     tofill = t(t(sort(unique(datum$genus[idx]))))
     browser()
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     List the default name for all families if genera was unknown.                     #
   #---------------------------------------------------------------------------------------#
   f2nogenus = list()

   f2nogenus[[  1]] = list(family = "Acanthaceae"      , genus = "Ignotum.acanthus"     )
   f2nogenus[[  2]] = list(family = "Achariaceae"      , genus = "Ignotum.acharia"      )
   f2nogenus[[  3]] = list(family = "Anacardiaceae"    , genus = "Ignotum.anacardium"   )
   f2nogenus[[  4]] = list(family = "Anisophylleaceae" , genus = "Ignotum.anisophyllea" )
   f2nogenus[[  5]] = list(family = "Annonaceae"       , genus = "Ignotum.annona"       )
   f2nogenus[[  6]] = list(family = "Apocynaceae"      , genus = "Ignotum.apocynum"     )
   f2nogenus[[  7]] = list(family = "Aquifoliaceae"    , genus = "Ignotum.ilex"         )
   f2nogenus[[  8]] = list(family = "Araliaceae"       , genus = "Ignotum.aralia"       )
   f2nogenus[[  9]] = list(family = "Arecaceae"        , genus = "Ignotum.areca"        )
   f2nogenus[[ 10]] = list(family = "Asteraceae"       , genus = "Ignotum.aster"        )
   f2nogenus[[ 11]] = list(family = "Bignoniaceae"     , genus = "Ignotum.bignonia"     )
   f2nogenus[[ 12]] = list(family = "Bixaceae"         , genus = "Ignotum.bixa"         )
   f2nogenus[[ 13]] = list(family = "Boraginaceae"     , genus = "Ignotum.borago"       )
   f2nogenus[[ 14]] = list(family = "Burseraceae"      , genus = "Ignotum.bursera"      )
   f2nogenus[[ 15]] = list(family = "Cactaceae"        , genus = "Ignotum.cactus"       )
   f2nogenus[[ 16]] = list(family = "Calophyllaceae"   , genus = "Ignotum.calophyllum"  )
   f2nogenus[[ 17]] = list(family = "Cannabaceae"      , genus = "Ignotum.cannabis"     )
   f2nogenus[[ 18]] = list(family = "Capparaceae"      , genus = "Ignotum.capparis"     )
   f2nogenus[[ 19]] = list(family = "Cardiopteridaceae", genus = "Ignotum.cardiopteris" )
   f2nogenus[[ 20]] = list(family = "Caricaceae"       , genus = "Ignotum.carica"       )
   f2nogenus[[ 21]] = list(family = "Caryocaraceae"    , genus = "Ignotum.caryocar"     )
   f2nogenus[[ 22]] = list(family = "Celastraceae"     , genus = "Ignotum.celastrus"    )
   f2nogenus[[ 23]] = list(family = "Chrysobalanaceae" , genus = "Ignotum.chrysobalanus")
   f2nogenus[[ 24]] = list(family = "Clusiaceae"       , genus = "Ignotum.clusia"       )
   f2nogenus[[ 25]] = list(family = "Combretaceae"     , genus = "Ignotum.combretum"    )
   f2nogenus[[ 26]] = list(family = "Connaraceae"      , genus = "Ignotum.connarus"     )
   f2nogenus[[ 27]] = list(family = "Convolvulaceae"   , genus = "Ignotum.convolvulus"  )
   f2nogenus[[ 28]] = list(family = "Dichapetalaceae"  , genus = "Ignotum.dichapetalum" )
   f2nogenus[[ 29]] = list(family = "Dilleniaceae"     , genus = "Ignotum.dillenia"     )
   f2nogenus[[ 30]] = list(family = "Ebenaceae"        , genus = "Ignotum.ebenus"       )
   f2nogenus[[ 31]] = list(family = "Elaeocarpaceae"   , genus = "Ignotum.elaeocarpus"  )
   f2nogenus[[ 32]] = list(family = "Emmotaceae"       , genus = "Ignotum.emmotum"      )
   f2nogenus[[ 33]] = list(family = "Erythroxylaceae"  , genus = "Ignotum.erythroxylum" )
   f2nogenus[[ 34]] = list(family = "Euphorbiaceae"    , genus = "Ignotum.euphorbia"    )
   f2nogenus[[ 35]] = list(family = "Fabaceae"         , genus = "Ignotum.faba"         )
   f2nogenus[[ 36]] = list(family = "Goupiaceae"       , genus = "Goupia"               )
   f2nogenus[[ 37]] = list(family = "Humiriaceae"      , genus = "Ignotum.humiria"      )
   f2nogenus[[ 38]] = list(family = "Hypericaceae"     , genus = "Ignotum.hypericum"    )
   f2nogenus[[ 39]] = list(family = "Icacinaceae"      , genus = "Ignotum.icacina"      )
   f2nogenus[[ 40]] = list(family = "Ignotaceae"       , genus = "Ignotum"              )
   f2nogenus[[ 41]] = list(family = "Lacistemataceae"  , genus = "Ignotum.lacistema"    )
   f2nogenus[[ 42]] = list(family = "Lamiaceae"        , genus = "Ignotum.lamium"       )
   f2nogenus[[ 43]] = list(family = "Lauraceae"        , genus = "Ignotum.laurus"       )
   f2nogenus[[ 44]] = list(family = "Lecythidaceae"    , genus = "Ignotum.lecythis"     )
   f2nogenus[[ 45]] = list(family = "Lianaceae"        , genus = "Liana"                )
   f2nogenus[[ 46]] = list(family = "Linaceae"         , genus = "Ignotum.linum"        )
   f2nogenus[[ 47]] = list(family = "Loganiaceae"      , genus = "Ignotum.logania"      )
   f2nogenus[[ 48]] = list(family = "Malpighiaceae"    , genus = "Ignotum.malpighia"    )
   f2nogenus[[ 49]] = list(family = "Malvaceae"        , genus = "Ignotum.malva"        )
   f2nogenus[[ 50]] = list(family = "Marcgraviaceae"   , genus = "Ignotum.marcgravia"   )
   f2nogenus[[ 51]] = list(family = "Melastomataceae"  , genus = "Ignotum.melastoma"    )
   f2nogenus[[ 52]] = list(family = "Meliaceae"        , genus = "Ignotum.melia"        )
   f2nogenus[[ 53]] = list(family = "Moraceae"         , genus = "Ignotum.morus"        )
   f2nogenus[[ 54]] = list(family = "Myristicaceae"    , genus = "Ignotum.myristica"    )
   f2nogenus[[ 55]] = list(family = "Myrtaceae"        , genus = "Ignotum.myrtus"       )
   f2nogenus[[ 56]] = list(family = "Nyctaginaceae"    , genus = "Ignotum.nyctaginia"   )
   f2nogenus[[ 57]] = list(family = "Ochnaceae"        , genus = "Ignotum.ochna"        )
   f2nogenus[[ 58]] = list(family = "Olacaceae"        , genus = "Ignotum.olax"         )
   f2nogenus[[ 59]] = list(family = "Opiliaceae"       , genus = "Ignotum.opilia"       )
   f2nogenus[[ 60]] = list(family = "Peridiscaceae"    , genus = "Ignotum.peridiscus"   )
   f2nogenus[[ 61]] = list(family = "Phyllanthaceae"   , genus = "Ignotum.phyllanthus"  )
   f2nogenus[[ 62]] = list(family = "Primulaceae"      , genus = "Ignotum.primula"      )
   f2nogenus[[ 63]] = list(family = "Polygonaceae"     , genus = "Ignotum.polygonum"    )
   f2nogenus[[ 64]] = list(family = "Putranjivaceae"   , genus = "Ignotum.putranjiva"   )
   f2nogenus[[ 65]] = list(family = "Rosaceae"         , genus = "Ignotum.rosa"         )
   f2nogenus[[ 66]] = list(family = "Rubiaceae"        , genus = "Ignotum.rubia"        )
   f2nogenus[[ 67]] = list(family = "Rutaceae"         , genus = "Ignotum.ruta"         )
   f2nogenus[[ 68]] = list(family = "Salicaceae"       , genus = "Ignotum.salix"        )
   f2nogenus[[ 69]] = list(family = "Sapindaceae"      , genus = "Ignotum.sapindus"     )
   f2nogenus[[ 70]] = list(family = "Sapotaceae"       , genus = "Ignotum.sapota"       )
   f2nogenus[[ 71]] = list(family = "Simaroubaceae"    , genus = "Ignotum.simarouba"    )
   f2nogenus[[ 72]] = list(family = "Siparunaceae"     , genus = "Ignotum.siparuna"     )
   f2nogenus[[ 73]] = list(family = "Solanaceae"       , genus = "Ignotum.solanum"      )
   f2nogenus[[ 74]] = list(family = "Stemonuraceae"    , genus = "Ignotum.stemonurus"   )
   f2nogenus[[ 75]] = list(family = "Strelitziaceae"   , genus = "Ignotum.strelitzia"   )
   f2nogenus[[ 76]] = list(family = "Styracaceae"      , genus = "Ignotum.styrax"       )
   f2nogenus[[ 77]] = list(family = "Ulmaceae"         , genus = "Ignotum.ulmus"        )
   f2nogenus[[ 78]] = list(family = "Urticaceae"       , genus = "Ignotum.urtica"       )
   f2nogenus[[ 79]] = list(family = "Verbenaceae"      , genus = "Ignotum.verbena"      )
   f2nogenus[[ 80]] = list(family = "Violaceae"        , genus = "Ignotum.viola"        )
   f2nogenus[[ 81]] = list(family = "Vochysiaceae"     , genus = "Ignotum.vochysia"     )
   #---------------------------------------------------------------------------------------#





   #----- Convert the g2f list into a data frame. -----------------------------------------#
   f2nogenus = data.frame( apply( X      = t(sapply(X = f2nogenus , FUN = c))
                                , MARGIN = c(1,2)
                                , FUN    = unlist
                                )#end apply
                         , stringsAsFactors = FALSE
                         )#end data.frame
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Check whether there are families missing in f2nogenus look up table.              #
   #---------------------------------------------------------------------------------------#
   idx = which(! g2f$family %in% f2nogenus$family)
   if (length(idx) > 0){
     cat(" You must add families to f2nogenus...","\n")
     family.tofill = t(t(sort(unique(g2f$family[idx]))))
     browser()
   }#end if
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Fill in non-informative genus for individual with known family but unknown genus. #
   #---------------------------------------------------------------------------------------#
   nothing =   is.na(datum$genus) & is.na(datum$family)
   no.gen  =   is.na(datum$genus) & ! is.na(datum$family)
   yes.gen = ! is.na(datum$genus)
   #------ No information whatsoever.  Fill in with non-informative taxa. -----------------#
   datum$genus      [nothing] = "Ignotum"
   datum$scientific [nothing] = "Ignotum NA"
   datum$family     [nothing] = "Ignotaceae"
   #------ No genus, family only. ---------------------------------------------------------#
   idx                        = match(datum$family[no.gen],f2nogenus$family)
   datum$genus      [no.gen ] = f2nogenus$genus[idx]
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
   tnf[[ 51]] = list( common     = "castanha do para"          
                    , scientific = "Bertholletia excelsa"            
                    )#end list
   tnf[[ 52]] = list( common     = "castanha sapucaia"         
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
   tnf[[ 70]] = list( common     = "envira cana"               
                    , scientific = "Xylopia nitida"                  
                    )#end list
   tnf[[ 71]] = list( common     = "envira preta"              
                    , scientific = "Guatteria poeppigiana"           
                    )#end list
   tnf[[ 72]] = list( common     = "envira surucucu"           
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
   tnf[[ 78]] = list( common     = "fava bolota"               
                    , scientific = "Parkia pendula"                  
                    )#end list
   tnf[[ 79]] = list( common     = "fava da rosca"             
                    , scientific = "Enterolobium schomburgkii"       
                    )#end list
   tnf[[ 80]] = list( common     = "fava folha fina"           
                    , scientific = "Pseudopiptadenia psilostachy"    
                    )#end list
   tnf[[ 81]] = list( common     = "fava mapuxiqui"            
                    , scientific = "Albizia pedicellaris"            
                    )#end list
   tnf[[ 82]] = list( common     = "fava saboeiro"             
                    , scientific = "Abarema"                         
                    )#end list
   tnf[[ 83]] = list( common     = "fava timbauba"             
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
   tnf[[105]] = list( common     = "janita"                    
                    , scientific = "Brosimum guianense"              
                    )#end list
   tnf[[106]] = list( common     = "jarana"                    
                    , scientific = "Lecythis lurida"                 
                    )#end list
   tnf[[107]] = list( common     = "jatauba"                   
                    , scientific = "Matayba purgans"                 
                    )#end list
   tnf[[108]] = list( common     = "joao mole"                 
                    , scientific = "Guapira venosa"                  
                    )#end list
   tnf[[109]] = list( common     = "joao mole grande"          
                    , scientific = "Neea"                            
                    )#end list
   tnf[[110]] = list( common     = "jutai pororoca"            
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
   tnf[[114]] = list( common     = "lacre da mata"             
                    , scientific = "Vismia"                          
                    )#end list
   tnf[[115]] = list( common     = "lacre vermelho"            
                    , scientific = "Vismia latifolia"                
                    )#end list
   tnf[[116]] = list( common     = "louro"                     
                    , scientific = "Nectandra pulverulenta"          
                    )#end list
   tnf[[117]] = list( common     = "louro abacate"             
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
   tnf[[150]] = list( common     = "muruci da mata"            
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
   tnf[[155]] = list( common     = "papaterra"                 
                    , scientific = "Miconia ruficalyx"               
                    )#end list
   tnf[[156]] = list( common     = "papaterra amarelo"         
                    , scientific = "Miconia lepidota"                
                    )#end list
   tnf[[157]] = list( common     = "papaterra folha peluda"    
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
   tnf[[161]] = list( common     = "pato-de-mutum"             
                    , scientific = "Lacunaria crenata"               
                    )#end list
   tnf[[162]] = list( common     = "pau-cobra"                 
                    , scientific = "Salacia impressifolia"           
                    )#end list
   tnf[[163]] = list( common     = "pau-de-arco amarelo"       
                    , scientific = "Tabebuia serratifolia"           
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
   tnf[[169]] = list( common     = "pepino da mata"            
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
   tnf[[202]] = list( common     = "ucuuba-terra-firme"        
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
   gs.mat           = t(as.data.frame(gs.list))
   g.only           = gs.length < 2
   gs.mat[g.only,2] = NA
   datum$genus      = capwords(gs.mat[,1],strict=TRUE)
   datum$scientific = paste(datum$genus,tolower (gs.mat[,2]),sep=" ")
   #---------------------------------------------------------------------------------------#
   

   #---------------------------------------------------------------------------------------#
   #      Fill in the family.                                                              #
   #---------------------------------------------------------------------------------------#
   datum$family = standard.family.name(datum)
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
#------------------------------------------------------------------------------------------#
find.wood.density <<- function(datum,wood,verbose=FALSE){

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
   for (s in 1:nspecies){
      if (! is.na(species[s]) && length(grep("Ignotum",species[s])) == 0){
         if (verbose) cat("   - ",s,"/",nspecies,"  -  ",species[s],"...","\n") 

         #----- Get the genus and family of this species, in case we need it. -------------#
         igen        = which (datum$scientific == species[s])
         this.genus  = unique(datum$genus[igen])
         this.family = unique(datum$family[igen])
         if (length(this.genus) != length(this.family)){
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
            iwood = intersect( which(wood$scientific == species[s] )
                             , intersect( which( wood$genus     == this.genus )
                                        , which( wood$family    == this.family) ) )
            if (length(iwood) == 0){
               #----- Intersection was zero! Misidentification, maybe? --------------------#
               loose = TRUE
            }else{
               if (any(is.finite(wood$density[iwood]))){
                  sel                     = datum$scientific == species[s]
                  datum$wood.dens   [sel] = wood$density[iwood]
                  datum$gf.wood.dens[sel] = 0
                  fill.genus              = FALSE
                  loose                   = FALSE
               }else{
                  loose                   = TRUE
                  fill.genus              = this.genus %in% wood$genus
               }#end if
               #---------------------------------------------------------------------------#
            }#end if
         #---------------------------------------------------------------------------------#
         }else{
            loose      = TRUE
            fill.genus = this.genus %in% wood$genus
         }#end if
         #---------------------------------------------------------------------------------#

         if (fill.genus){
            #------------------------------------------------------------------------------#
            #     Find all the plants from this genus, take the average, and attribute     #
            # to this species.  In this case, warn the user about the species, it may      #
            # be a typo in the scientific name that is easy to fix.                        #
            #------------------------------------------------------------------------------#
            iwood              = intersect( which( wood$genus  == this.genus )
                                          , which( wood$family == this.family) )
            if (length(iwood) == 0){
               #----- Intersection was zero! Misidentification, maybe? --------------------#
               loose      = TRUE
            }else{
               if (any(is.finite(wood$density[iwood]))){
                  sel                     = datum$scientific == species[s]
                  datum$wood.dens   [sel] = mean(wood$density[iwood],na.rm=TRUE)
                  datum$gf.wood.dens[sel] = 1
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
   #     List all genera that didn't belong to any known family.                           #
   #---------------------------------------------------------------------------------------#
   if (verbose && nrow(sci.loose.mean) > 0){
      cat (" Found genera that didn't belong to any known family!","\n")
      print(sci.loose.mean)
      browser()
   }#end if


   genus.only = grep(" NA",sci.genus.mean[,"scientific"])

   if (verbose){
     cat ("Only genus was provided: we used the genus mean wood density:","\n")
     print (sci.genus.mean[genus.only,])

     cat   ("Species was provided but not found in the database:","\n")
     print (sci.genus.mean[-genus.only,])
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
   for (f in 1:nfamilies){
      if (families[f] != "Ignotaceae"){
         if (verbose) cat("   - ",f,"/",nfamilies,"  -  ",families[f],"...","\n")

         #----- Get the individuals that belong to this family. ---------------------------#
         ifam         = which (datum$family == families[f]  & is.finite(datum$wood.dens))
         imiss        = which (datum$family == families[f]  & is.na    (datum$wood.dens))
         if (length(imiss) > 0 && length(ifam) > 0){
            sample.wood.dens          = sample( x       = datum$wood.dens[ifam]
                                              , size    = length(imiss)
                                              , replace = TRUE
                                              )#end sample
            datum$wood.dens   [imiss] = sample.wood.dens
            datum$gf.wood.dens[imiss] = 2
            
            gf2               = cbind(datum$scientific[imiss]
                                     ,datum$family    [imiss]
                                     ,sprintf("%6.3f",sample.wood.dens)
                                     )#end cbind
            sci.family.sample = rbind(sci.family.sample,gf2)

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
         cat (" Found families that unidientified genera!","\n")
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

      datum$wood.dens   [imiss] = sample(x=datum$wood.dens[-imiss],size=nmiss,replace=TRUE)
      datum$gf.wood.dens[imiss] = 3
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Adjust the gap-filling flag for wood density by adding whether the scientific     #
   # name is gap-filled.                                                                   #
   #---------------------------------------------------------------------------------------#
   datum$gf.wood.dens = datum$gf.wood.dens + 10 * datum$gf.scientific
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
   gs.mat                   = t(as.data.frame(gs.list))
   g.only                   = gs.length < 2
   gs.mat[g.only,2]         = NA
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
   gs.mat                   = t(as.data.frame(gs.list))
   g.only                   = gs.length < 2
   gs.mat[g.only,2]         = NA
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
