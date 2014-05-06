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
   sel = (! is.na(x) & x == "?"                       ); x[sel] = NA
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
   sel = (! is.na(x) & x == "cip"                     ); x[sel] = "cipo"
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
   dat$scientific = tolower(dat$scientific)
   #----- Break into genus and species. ---------------------------------------------------#
   g.s     = keep.gen.spe.only(dat$scientific)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Substitutions of scientific name.  We standardise these first so family becomes  #
   # easier.                                                                               #
   #---------------------------------------------------------------------------------------#
   g.s = sub("Aiouea densiflora"          ,"Aiouea laevis"                ,x=g.s)
   g.s = sub("Allophyllus floribunda"     ,"Allophylus floribundus"       ,x=g.s)
   g.s = sub("Ampelocera endentula"       ,"Ampelocera edentula"          ,x=g.s)
   g.s = sub("Amphirrhox longiflora"      ,"Amphirrhox longifolia"        ,x=g.s)
   g.s = sub("Amphirrhox surinamensis"    ,"Amphirrhox longifolia"        ,x=g.s)
   g.s = sub("Anadenanthera falcata"      ,"Anadenanthera peregrina"      ,x=g.s)
   g.s = sub("Anartia"                    ,"Tabernaemontana"              ,x=g.s)
   g.s = sub("Aniba roseodora"            ,"Aniba rosaeodora"             ,x=g.s)
   g.s = sub("Annona decicoma"            ,"Annona densicoma"             ,x=g.s)
   g.s = sub("Apeiba burchelii"           ,"Apeiba glabra"                ,x=g.s)
   g.s = sub("Aspidosperma aracanga"      ,"Aspidosperma araracanga"      ,x=g.s)
   g.s = sub("Aspidosperma auriculata"    ,"Aspidosperma auriculatum"     ,x=g.s)
   g.s = sub("Aspidosperma cruentum"      ,"Aspidosperma desmanthum"      ,x=g.s)
   g.s = sub("Aspidosperma desmantum"     ,"Aspidosperma desmanthum"      ,x=g.s)
   g.s = sub("Aspidosperma desmathum"     ,"Aspidosperma desmanthum"      ,x=g.s)
   g.s = sub("Aspidosperma eteanun"       ,"Aspidosperma eteanum"         ,x=g.s)
   g.s = sub("Aspidosperma nitidum"       ,"Aspidosperma excelsum"        ,x=g.s)
   g.s = sub("Astronium le-cointei"       ,"Astronium lecointei"          ,x=g.s)
   g.s = sub("Austroplenckia populnea"    ,"Plenckia populnea"            ,x=g.s)
   g.s = sub("Balisia pedicelares"        ,"Albizia pedicellaris"         ,x=g.s)
   g.s = sub("Balizia pedicellaris"       ,"Albizia pedicellaris"         ,x=g.s)
   g.s = sub("Bauhinia jarensis"          ,"Bauhinia deleteme"            ,x=g.s)
   g.s = sub("Bellucia grossulariodis"    ,"Bellucia grossularioides"     ,x=g.s)
   g.s = sub("Brosimum autifolium"        ,"Brosimum acutifolium"         ,x=g.s)
   g.s = sub("Brosimum lactascens"        ,"Brosimum lactescens"          ,x=g.s)
   g.s = sub("Byrsonima schultesiana"     ,"Byrsonima arthropoda"         ,x=g.s)
   g.s = sub("Byrsonima estipulacea"      ,"Byrsonima stipulacea"         ,x=g.s)
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
   g.s = sub("Coussarea racemosa"         ,"Coussarea albescens"          ,x=g.s)
   g.s = sub("Crepidospermum gondotiano"  ,"Crepidospermum goudotianum"   ,x=g.s)
   g.s = sub("Crepidospermum goudotiano"  ,"Crepidospermum goudotianum"   ,x=g.s)
   g.s = sub("Cupania hirta"              ,"Cupania hirsuta"              ,x=g.s)
   g.s = sub("Cybistax antisiphyllitica"  ,"Cybistax antisyphilitica"     ,x=g.s)
   g.s = sub("Dendrobrangea boliviana"    ,"Dendrobangia boliviana"       ,x=g.s)
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
   g.s = sub("Eriotreca globosa"          ,"Eriotheca globosa"            ,x=g.s)
   g.s = sub("Eschweilera amazomica"      ,"Eschweilera amazonica"        ,x=g.s)
   g.s = sub("Eschweilera apiculatum"     ,"Eschweilera apiculata"        ,x=g.s)
   g.s = sub("Eschweilera idatimom"       ,"Lecythis idatimon"            ,x=g.s)
   g.s = sub("Eschweilera observa"        ,"Eschweilera obversa"          ,x=g.s)
   g.s = sub("Eschweilera pedicelata"     ,"Eschweilera pedicellata"      ,x=g.s)
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
   g.s = sub("Inga aff\\."                ,"Inga affinis"                 ,x=g.s)
   g.s = sub("Inga jenmanii"              ,"Inga sertulifera"             ,x=g.s)
   g.s = sub("Inga dibaldiana"            ,"Inga thibaudiana"             ,x=g.s)
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
   g.s = sub("Maluria"                    ,"Marlierea umbraticola"        ,x=g.s)
   g.s = sub("Maquira callophylla"        ,"Maquira calophylla"           ,x=g.s)
   g.s = sub("Marmaroxylon racemosum"     ,"Zygia racemosa"               ,x=g.s)
   g.s = sub("Maytenos guianensis"        ,"Maytenus guyanensis"          ,x=g.s)
   g.s = sub("Maytenus guianensis"        ,"Maytenus guyanensis"          ,x=g.s)
   g.s = sub("Meia maderensis"            ,"Neea"                         ,x=g.s)
   g.s = sub("Mezelaurus"                 ,"Mezilaurus"                   ,x=g.s)
   g.s = sub("Mezelaurus itauba"          ,"Mezilaurus itauba"            ,x=g.s)
   g.s = sub("Miconia chrysophyllum"      ,"Miconia chrysophylla"         ,x=g.s)
   g.s = sub("Michopholis venulosa"       ,"Micropholis venulosa"         ,x=g.s)
   g.s = sub("Microphilis"                ,"Micropholis"                  ,x=g.s)
   g.s = sub("Micropholis guianensis"     ,"Micropholis guyanensis"       ,x=g.s)
   g.s = sub("Microphylis acutangula"     ,"Micropholis acutangula"       ,x=g.s)
   g.s = sub("Mouriri abnormis"           ,"Votomita guianensis"          ,x=g.s)
   g.s = sub("Myrciaria reticulata"       ,"Myrcia reticulata"            ,x=g.s)
   g.s = sub("Myrcia rutipula"            ,"Myrcia rufipila"              ,x=g.s)
   g.s = sub("Newtonia psilostachya"      ,"Pseudopiptadenia psilostachya",x=g.s)
   g.s = sub("Ocotea baturitensis"        ,"Ocotea"                       ,x=g.s)
   g.s = sub("Ocotea caudata"             ,"Ocotea cernua"                ,x=g.s)
   g.s = sub("Omedia perebea"             ,"Perebea mollis"               ,x=g.s)
   g.s = sub("Onichiopetalum amazonico"   ,"Onychopetalum amazonicum"     ,x=g.s)
   g.s = sub("Pausandra densiflora"       ,"Pausandra trianae"            ,x=g.s)
   g.s = sub("Pouroma guianensis"         ,"Pourouma guianensis"          ,x=g.s)
   g.s = sub("Pouruma guianensis"         ,"Pourouma guianensis"          ,x=g.s)
   g.s = sub("Pourouma vilosa"            ,"Pourouma villosa"             ,x=g.s)
   g.s = sub("Pouteria ambelanifolia"     ,"Pouteria ambelaniifolia"      ,x=g.s)
   g.s = sub("Pouteria biloculares"       ,"Pouteria bilocularis"         ,x=g.s)
   g.s = sub("Pouteria filipis"           ,"Pouteria filipes"             ,x=g.s)
   g.s = sub("Pouteria lasiocarpa"        ,"Pouteria caimito"             ,x=g.s)
   g.s = sub("Pouteria gonggrijpii"       ,"Pouteria gongrijpii"          ,x=g.s)
   g.s = sub("Pouteria heterosepala"      ,"Pouteria polysepala"          ,x=g.s)
   g.s = sub("Pouteria paraensis"         ,"Pouteria macrocarpa"          ,x=g.s)
   g.s = sub("Protiumuceanum"             ,"Protium spruceanum"           ,x=g.s)
   g.s = sub("Protium heptafilum"         ,"Protium heptaphyllum"         ,x=g.s)
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
   g.s = sub("Salacea"                    ,"Salacia"                      ,x=g.s)
   g.s = sub("Salacia imprissifolia"      ,"Salacia impressifolia"        ,x=g.s)
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
   g.s = sub("Trattinickia laurencei"     ,"Trattinnickia lawrancei"      ,x=g.s)
   g.s = sub("Trattinickia laurense"      ,"Trattinnickia lawrancei"      ,x=g.s)
   g.s = sub("Trattinickia rhoifolia"     ,"Trattinnickia rhoifolia"      ,x=g.s)
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
   g.s = sub("Virola melinonii"           ,"Virola michelii"              ,x=g.s)
   g.s = sub("Virola melionii"            ,"Virola michelii"              ,x=g.s)
   g.s = sub("Virola michelli"            ,"Virola michelii"              ,x=g.s)
   g.s = sub("Vochisia surinamensis"      ,"Vochysia surinamensis"        ,x=g.s)
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
   datum$family = capwords(datum$family,strict=TRUE)
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
   datum$family = sub("Mimosaceae"                  ,"Fabaceae"      ,x=datum$family)
   datum$family = sub("Papilionaceae"               ,"Fabaceae"      ,x=datum$family)
   datum$family = sub("Tiliaceae"                   ,"Malvaceae"     ,x=datum$family)
   datum$family = sub("Sterculiaceae"               ,"Malvaceae"     ,x=datum$family)
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Look-up table for all families.                                                   #
   #---------------------------------------------------------------------------------------#
   n=0  ; g2f      = list()
   n=n+1; g2f[[n]] = list( genus = "Abarema"            , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Acmanthera"         , family = "Malpighiaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Acrocomia"          , family = "Arecaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Acosmium"           , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Adenophaedra"       , family = "Euphorbiaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Agonandra"          , family = "Opiliaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Aiouea"             , family = "Lauraceae"         )
   n=n+1; g2f[[n]] = list( genus = "Albizia"            , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Alchorneopsis"      , family = "Euphorbiaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Aldina"             , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Alexa"              , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Alibertia"          , family = "Rubiaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Allantoma"          , family = "Lecythidaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Allophylus"         , family = "Sapindaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Amaioua"            , family = "Rubiaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Ambelania"          , family = "Apocynaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Amburana"           , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Ampelocera"         , family = "Ulmaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Amphirrhox"         , family = "Violaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Anacardium"         , family = "Anacardiaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Anadenanthera"      , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Anaxagorea"         , family = "Annonaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Andira"             , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Aniba"              , family = "Lauraceae"         )
   n=n+1; g2f[[n]] = list( genus = "Anisophyllea"       , family = "Anisophylleaceae"  )
   n=n+1; g2f[[n]] = list( genus = "Annona"             , family = "Annonaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Anomalocalyx"       , family = "Euphorbiaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Antonia"            , family = "Loganiaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Aparisthmium"       , family = "Euphorbiaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Apeiba"             , family = "Malvaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Aptandra"           , family = "Olacaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Apuleia"            , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Aspidosperma"       , family = "Apocynaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Astrocaryum"        , family = "Arecaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Astronium"          , family = "Anacardiaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Attalea"            , family = "Arecaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Austromyrtus"       , family = "Myrtaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Bactris"            , family = "Arecaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Bagassa"            , family = "Moraceae"          )
   n=n+1; g2f[[n]] = list( genus = "Balizia"            , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Batesia"            , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Batocarpus"         , family = "Moraceae"          )
   n=n+1; g2f[[n]] = list( genus = "Bauhinia"           , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Bellucia"           , family = "Melastomataceae"   )
   n=n+1; g2f[[n]] = list( genus = "Bertholletia"       , family = "Lecythidaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Blepharocalyx"      , family = "Myrtaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Bixa"               , family = "Bixaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Bocageopsis"        , family = "Annonaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Bocoa"              , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Bombax"             , family = "Malvaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Botryarrhena"       , family = "Rubiaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Bowdichia"          , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Brachychiton"       , family = "Malvaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Brosimum"           , family = "Moraceae"          )
   n=n+1; g2f[[n]] = list( genus = "Buchenavia"         , family = "Combretaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Byrsonima"          , family = "Malpighiaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Caesalpinia"        , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Calliandra"         , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Calophyllum"        , family = "Calophyllaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Calyptranthes"      , family = "Myrtaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Campomanesia"       , family = "Myrtaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Capirona"           , family = "Rubiaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Capparidastrum"     , family = "Capparaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Capparis"           , family = "Capparaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Caraipa"            , family = "Calophyllaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Carapa"             , family = "Meliaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Cariniana"          , family = "Lecythidaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Caryocar"           , family = "Caryocaraceae"     )
   n=n+1; g2f[[n]] = list( genus = "Casearia"           , family = "Salicaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Cassia"             , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Castilla"           , family = "Moraceae"          )
   n=n+1; g2f[[n]] = list( genus = "Catostemma"         , family = "Malvaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Cecropia"           , family = "Urticaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Cedrela"            , family = "Meliaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Cedrelinga"         , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Ceiba"              , family = "Malvaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Cereus"             , family = "Cactaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Chaetocarpus"       , family = "Euphorbiaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Chaetochlamys"      , family = "Acanthaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Chamaecrista"       , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Chaunochiton"       , family = "Olacaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Cheiloclinium"      , family = "Celastraceae"      )
   n=n+1; g2f[[n]] = list( genus = "Chimarrhis"         , family = "Rubiaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Chloroleucon"       , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Chromolucuma"       , family = "Sapotaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Chrysobalanus"      , family = "Chrysobalanaceae"  )
   n=n+1; g2f[[n]] = list( genus = "Chrysophyllum"      , family = "Sapotaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Clarisia"           , family = "Moraceae"          )
   n=n+1; g2f[[n]] = list( genus = "Clusia"             , family = "Clusiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Cnidoscolus"        , family = "Euphorbiaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Coccoloba"          , family = "Polygonaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Commiphora"         , family = "Burseraceae"       )
   n=n+1; g2f[[n]] = list( genus = "Conceveiba"         , family = "Euphorbiaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Connarus"           , family = "Connaraceae"       )
   n=n+1; g2f[[n]] = list( genus = "Copaifera"          , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Cordia"             , family = "Boraginaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Corythophora"       , family = "Lecythidaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Couepia"            , family = "Chrysobalanaceae"  )
   n=n+1; g2f[[n]] = list( genus = "Couma"              , family = "Apocynaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Couratari"          , family = "Lecythidaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Coussapoa"          , family = "Urticaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Coussarea"          , family = "Rubiaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Coutarea"           , family = "Rubiaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Crepidospermum"     , family = "Burseraceae"       )
   n=n+1; g2f[[n]] = list( genus = "Croton"             , family = "Euphorbiaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Crudia"             , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Cupania"            , family = "Sapindaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Cybianthus"         , family = "Primulaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Cybistax"           , family = "Bignoniaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Dacryodes"          , family = "Burseraceae"       )
   n=n+1; g2f[[n]] = list( genus = "Dalbergia"          , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Dendrobangia"       , family = "Cardiopteridaceae" )
   n=n+1; g2f[[n]] = list( genus = "Dialium"            , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Diclinanona"        , family = "Annonaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Dicorynia"          , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Dicranostyles"      , family = "Convolvulaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Dicypellium"        , family = "Lauraceae"         )
   n=n+1; g2f[[n]] = list( genus = "Dimorphandra"       , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Dinizia"            , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Diospyros"          , family = "Ebenaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Diplotropis"        , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Dipteryx"           , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Discophora"         , family = "Stemonuraceae"     )
   n=n+1; g2f[[n]] = list( genus = "Doliocarpus"        , family = "Dilleniaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Drypetes"           , family = "Putranjivaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Duckeodendron"      , family = "Solanaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Duckesia"           , family = "Humiriaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Duguetia"           , family = "Annonaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Duroia"             , family = "Rubiaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Dystovomita"        , family = "Clusiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Ecclinusa"          , family = "Sapotaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Elaeoluma"          , family = "Sapotaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Elizabetha"         , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Emmotum"            , family = "Emmotaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Endlicheria"        , family = "Lauraceae"         )
   n=n+1; g2f[[n]] = list( genus = "Endopleura"         , family = "Humiriaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Enterolobium"       , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Eperua"             , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Ephedranthus"       , family = "Annonaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Eriotheca"          , family = "Malvaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Erisma"             , family = "Vochysiaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Erythroxylum"       , family = "Erythroxylaceae"   )
   n=n+1; g2f[[n]] = list( genus = "Eschweilera"        , family = "Lecythidaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Eugenia"            , family = "Myrtaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Euterpe"            , family = "Arecaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Euxylophora"        , family = "Rutaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Faramea"            , family = "Rubiaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Ferdinandusa"       , family = "Rubiaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Ficus"              , family = "Moraceae"          )
   n=n+1; g2f[[n]] = list( genus = "Fraunhofera"        , family = "Celastraceae"      )
   n=n+1; g2f[[n]] = list( genus = "Fusaea"             , family = "Annonaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Garcinia"           , family = "Clusiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Geissospermum"      , family = "Apocynaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Genipa"             , family = "Rubiaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Gilbertiodendron"   , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Glycydendron"       , family = "Euphorbiaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Goupia"             , family = "Goupiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Guapira"            , family = "Nyctaginaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Guarea"             , family = "Meliaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Guatteria"          , family = "Annonaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Guazuma"            , family = "Malvaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Gustavia"           , family = "Lecythidaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Hebepetalum"        , family = "Linaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Heisteria"          , family = "Olacaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Helianthostylis"    , family = "Moraceae"          )
   n=n+1; g2f[[n]] = list( genus = "Helicostylis"       , family = "Moraceae"          )
   n=n+1; g2f[[n]] = list( genus = "Henriettea"         , family = "Melastomataceae"   )
   n=n+1; g2f[[n]] = list( genus = "Henriettella"       , family = "Melastomataceae"   )
   n=n+1; g2f[[n]] = list( genus = "Hevea"              , family = "Euphorbiaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Himatanthus"        , family = "Apocynaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Hippocratea"        , family = "Celastraceae"      )
   n=n+1; g2f[[n]] = list( genus = "Hirtella"           , family = "Chrysobalanaceae"  )
   n=n+1; g2f[[n]] = list( genus = "Humiria"            , family = "Humiriaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Humiriastrum"       , family = "Humiriaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Hura"               , family = "Euphorbiaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Hymenaea"           , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Hymenolobium"       , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Ignotum"            , family = "Ignotaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Ilex"               , family = "Aquifoliaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Inga"               , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Iryanthera"         , family = "Myristicaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Isertia"            , family = "Rubiaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Jacaranda"          , family = "Bignoniaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Jacaratia"          , family = "Caricaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Jatropha"           , family = "Euphorbiaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Joannesia"          , family = "Euphorbiaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Justicia"           , family = "Acanthaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Lacistema"          , family = "Lacistemataceae"   )
   n=n+1; g2f[[n]] = list( genus = "Lacmellea"          , family = "Apocynaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Lacunaria"          , family = "Ochnaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Ladenbergia"        , family = "Rubiaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Laetia"             , family = "Salicaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Lafoensia"          , family = "Lythraceae"        )
   n=n+1; g2f[[n]] = list( genus = "Lecythis"           , family = "Lecythidaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Leonia"             , family = "Violaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Liana"              , family = "Lianaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Licania"            , family = "Chrysobalanaceae"  )
   n=n+1; g2f[[n]] = list( genus = "Licaria"            , family = "Lauraceae"         )
   n=n+1; g2f[[n]] = list( genus = "Lindackeria"        , family = "Achariaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Lippia"             , family = "Verbenaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Lonchocarpus"       , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Luehea"             , family = "Malvaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Lueheopsis"         , family = "Malvaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Mabea"              , family = "Euphorbiaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Machaerium"         , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Macoubea"           , family = "Apocynaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Macrolobium"        , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Malouetia"          , family = "Apocynaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Manicaria"          , family = "Arecaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Manihot"            , family = "Euphorbiaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Manilkara"          , family = "Sapotaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Maquira"            , family = "Moraceae"          )
   n=n+1; g2f[[n]] = list( genus = "Marcgravia"         , family = "Marcgraviaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Marlierea"          , family = "Myrtaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Matayba"            , family = "Sapindaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Matisia"            , family = "Malvaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Mauritia"           , family = "Arecaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Mauritiella"        , family = "Arecaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Maytenus"           , family = "Celastraceae"      )
   n=n+1; g2f[[n]] = list( genus = "Melicoccus"         , family = "Sapindaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Mezilaurus"         , family = "Lauraceae"         )
   n=n+1; g2f[[n]] = list( genus = "Miconia"            , family = "Melastomataceae"   )
   n=n+1; g2f[[n]] = list( genus = "Micrandra"          , family = "Euphorbiaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Micrandropsis"      , family = "Euphorbiaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Micropholis"        , family = "Sapotaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Mikania"            , family = "Asteraceae"        )
   n=n+1; g2f[[n]] = list( genus = "Mimosa"             , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Minquartia"         , family = "Olacaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Misanteca"          , family = "Lauraceae"         )
   n=n+1; g2f[[n]] = list( genus = "Monopteryx"         , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Moronobea"          , family = "Clusiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Mouriri"            , family = "Melastomataceae"   )
   n=n+1; g2f[[n]] = list( genus = "Myracrodruon"       , family = "Anacardiaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Myrocarpus"         , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Myrcia"             , family = "Myrtaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Myrciaria"          , family = "Myrtaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Myrsine"            , family = "Primulaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Naucleopsis"        , family = "Moraceae"          )
   n=n+1; g2f[[n]] = list( genus = "Nectandra"          , family = "Lauraceae"         )
   n=n+1; g2f[[n]] = list( genus = "Neea"               , family = "Nyctaginaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Ocotea"             , family = "Lauraceae"         )
   n=n+1; g2f[[n]] = list( genus = "Oenocarpus"         , family = "Arecaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Onychopetalum"      , family = "Annonaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Ormosia"            , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Osteophloeum"       , family = "Myristicaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Ouratea"            , family = "Ochnaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Oxandra"            , family = "Annonaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Pachira"            , family = "Malvaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Palicourea"         , family = "Rubiaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Parahancornia"      , family = "Apocynaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Paraia"             , family = "Lauraceae"         )
   n=n+1; g2f[[n]] = list( genus = "Parinari"           , family = "Chrysobalanaceae"  )
   n=n+1; g2f[[n]] = list( genus = "Parkia"             , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Paullinia"          , family = "Sapindaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Pausandra"          , family = "Euphorbiaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Paypayrola"         , family = "Violaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Peltogyne"          , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Pera"               , family = "Euphorbiaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Perebea"            , family = "Moraceae"          )
   n=n+1; g2f[[n]] = list( genus = "Peridiscus"         , family = "Peridiscaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Phenakospermum"     , family = "Strelitziaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Piptadenia"         , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Piptocarpha"        , family = "Asteraceae"        )
   n=n+1; g2f[[n]] = list( genus = "Pithecellobium"     , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Plathymenia"        , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Platonia"           , family = "Clusiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Platymiscium"       , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Platypodium"        , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Plenckia"           , family = "Celastraceae"      )
   n=n+1; g2f[[n]] = list( genus = "Poecilanthe"        , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Pogonophora"        , family = "Euphorbiaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Poraqueiba"         , family = "Icacinaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Posoqueria"         , family = "Rubiaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Pourouma"           , family = "Urticaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Pouteria"           , family = "Sapotaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Pradosia"           , family = "Sapotaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Protium"            , family = "Burseraceae"       )
   n=n+1; g2f[[n]] = list( genus = "Prunus"             , family = "Rosaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Pseudobombax"       , family = "Malvaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Pseudolmedia"       , family = "Moraceae"          )
   n=n+1; g2f[[n]] = list( genus = "Pseudopiptadenia"   , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Pseudoxandra"       , family = "Annonaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Psidium"            , family = "Myrtaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Pterandra"          , family = "Malpighiaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Pterocarpus"        , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Ptychopetalum"      , family = "Olacaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Qualea"             , family = "Vochysiaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Quararibea"         , family = "Malvaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Quiina"             , family = "Ochnaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Raputia"            , family = "Rutaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Rauvolfia"          , family = "Apocynaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Recordoxylon"       , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Rhodostemonodaphne" , family = "Lauraceae"         )
   n=n+1; g2f[[n]] = list( genus = "Richeria"           , family = "Phyllanthaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Rinorea"            , family = "Violaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Rollinia"           , family = "Annonaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Roucheria"          , family = "Linaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Roupala"            , family = "Proteaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Ruizterania"        , family = "Vochysiaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Ryania"             , family = "Salicaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Sacoglottis"        , family = "Humiriaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Sagotia"            , family = "Euphorbiaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Salacia"            , family = "Celastraceae"      )
   n=n+1; g2f[[n]] = list( genus = "Sapium"             , family = "Euphorbiaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Sarcaulus"          , family = "Sapotaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Schefflera"         , family = "Araliaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Schinopsis"         , family = "Anacardiaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Schizolobium"       , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Sclerolobium"       , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Scleronema"         , family = "Malvaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Senefeldera"        , family = "Euphorbiaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Senna"              , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Sextonia"           , family = "Lauraceae"         )
   n=n+1; g2f[[n]] = list( genus = "Simaba"             , family = "Simaroubaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Simarouba"          , family = "Simaroubaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Siparuna"           , family = "Siparunaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Sloanea"            , family = "Elaeocarpaceae"    )
   n=n+1; g2f[[n]] = list( genus = "Socratea"           , family = "Arecaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Sorocea"            , family = "Moraceae"          )
   n=n+1; g2f[[n]] = list( genus = "Spondias"           , family = "Anacardiaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Sterculia"          , family = "Malvaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Strychnos"          , family = "Loganiaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Stryphnodendron"    , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Styrax"             , family = "Styracaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Swartzia"           , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Swietenia"          , family = "Meliaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Symphonia"          , family = "Clusiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Syzygium"           , family = "Myrtaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Tabebuia"           , family = "Bignoniaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Tabernaemontana"    , family = "Apocynaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Tachigali"          , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Talisia"            , family = "Sapindaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Tapirira"           , family = "Anacardiaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Tapura"             , family = "Dichapetalaceae"   )
   n=n+1; g2f[[n]] = list( genus = "Taralea"            , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Terminalia"         , family = "Combretaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Tetragastris"       , family = "Burseraceae"       )
   n=n+1; g2f[[n]] = list( genus = "Theobroma"          , family = "Malvaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Thyrsodium"         , family = "Anacardiaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Tocoyena"           , family = "Rubiaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Toulicia"           , family = "Sapindaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Touroulia"          , family = "Ochnaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Tovomita"           , family = "Clusiaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Trattinnickia"      , family = "Burseraceae"       )
   n=n+1; g2f[[n]] = list( genus = "Trema"              , family = "Cannabaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Trichilia"          , family = "Meliaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Trymatococcus"      , family = "Moraceae"          )
   n=n+1; g2f[[n]] = list( genus = "Unonopsis"          , family = "Annonaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Vantanea"           , family = "Humiriaceae"       )
   n=n+1; g2f[[n]] = list( genus = "Vatairea"           , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Vataireopsis"       , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Vernonia"           , family = "Asteraceae"        )
   n=n+1; g2f[[n]] = list( genus = "Virola"             , family = "Myristicaceae"     )
   n=n+1; g2f[[n]] = list( genus = "Vismia"             , family = "Hypericaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Vitex"              , family = "Lamiaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Vochysia"           , family = "Vochysiaceae"      )
   n=n+1; g2f[[n]] = list( genus = "Votomita"           , family = "Melastomataceae"   )
   n=n+1; g2f[[n]] = list( genus = "Vouacapoua"         , family = "Fabaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Warszewiczia"       , family = "Rubiaceae"         )
   n=n+1; g2f[[n]] = list( genus = "Xylopia"            , family = "Annonaceae"        )
   n=n+1; g2f[[n]] = list( genus = "Zanthoxylum"        , family = "Rutaceae"          )
   n=n+1; g2f[[n]] = list( genus = "Zygia"              , family = "Fabaceae"          )
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
     cat(" You must add genera to g2f...","\n")
     tofill.gen = t(t(sort(unique(datum$genus[idx]))))
     toshow     = which(names(datum) %in% c("trans","tag","x","y","scientific","family"))
     toshow     = datum[idx,toshow]
     browser()
   }#end if
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     List the default name for all families if genera was unknown.                     #
   #---------------------------------------------------------------------------------------#
   n=0  ; f2nog      = list()
   n=n+1; f2nog[[n]] = list(family = "Acanthaceae"      , genus = "Ignotum.acanthus"     )
   n=n+1; f2nog[[n]] = list(family = "Achariaceae"      , genus = "Ignotum.acharia"      )
   n=n+1; f2nog[[n]] = list(family = "Anacardiaceae"    , genus = "Ignotum.anacardium"   )
   n=n+1; f2nog[[n]] = list(family = "Anisophylleaceae" , genus = "Ignotum.anisophyllea" )
   n=n+1; f2nog[[n]] = list(family = "Annonaceae"       , genus = "Ignotum.annona"       )
   n=n+1; f2nog[[n]] = list(family = "Apocynaceae"      , genus = "Ignotum.apocynum"     )
   n=n+1; f2nog[[n]] = list(family = "Aquifoliaceae"    , genus = "Ignotum.ilex"         )
   n=n+1; f2nog[[n]] = list(family = "Araliaceae"       , genus = "Ignotum.aralia"       )
   n=n+1; f2nog[[n]] = list(family = "Arecaceae"        , genus = "Ignotum.areca"        )
   n=n+1; f2nog[[n]] = list(family = "Asteraceae"       , genus = "Ignotum.aster"        )
   n=n+1; f2nog[[n]] = list(family = "Bignoniaceae"     , genus = "Ignotum.bignonia"     )
   n=n+1; f2nog[[n]] = list(family = "Bixaceae"         , genus = "Ignotum.bixa"         )
   n=n+1; f2nog[[n]] = list(family = "Boraginaceae"     , genus = "Ignotum.borago"       )
   n=n+1; f2nog[[n]] = list(family = "Burseraceae"      , genus = "Ignotum.bursera"      )
   n=n+1; f2nog[[n]] = list(family = "Cactaceae"        , genus = "Ignotum.cactus"       )
   n=n+1; f2nog[[n]] = list(family = "Calophyllaceae"   , genus = "Ignotum.calophyllum"  )
   n=n+1; f2nog[[n]] = list(family = "Cannabaceae"      , genus = "Ignotum.cannabis"     )
   n=n+1; f2nog[[n]] = list(family = "Capparaceae"      , genus = "Ignotum.capparis"     )
   n=n+1; f2nog[[n]] = list(family = "Cardiopteridaceae", genus = "Ignotum.cardiopteris" )
   n=n+1; f2nog[[n]] = list(family = "Caricaceae"       , genus = "Ignotum.carica"       )
   n=n+1; f2nog[[n]] = list(family = "Caryocaraceae"    , genus = "Ignotum.caryocar"     )
   n=n+1; f2nog[[n]] = list(family = "Celastraceae"     , genus = "Ignotum.celastrus"    )
   n=n+1; f2nog[[n]] = list(family = "Chrysobalanaceae" , genus = "Ignotum.chrysobalanus")
   n=n+1; f2nog[[n]] = list(family = "Clusiaceae"       , genus = "Ignotum.clusia"       )
   n=n+1; f2nog[[n]] = list(family = "Combretaceae"     , genus = "Ignotum.combretum"    )
   n=n+1; f2nog[[n]] = list(family = "Connaraceae"      , genus = "Ignotum.connarus"     )
   n=n+1; f2nog[[n]] = list(family = "Convolvulaceae"   , genus = "Ignotum.convolvulus"  )
   n=n+1; f2nog[[n]] = list(family = "Dichapetalaceae"  , genus = "Ignotum.dichapetalum" )
   n=n+1; f2nog[[n]] = list(family = "Dilleniaceae"     , genus = "Ignotum.dillenia"     )
   n=n+1; f2nog[[n]] = list(family = "Ebenaceae"        , genus = "Ignotum.ebenus"       )
   n=n+1; f2nog[[n]] = list(family = "Elaeocarpaceae"   , genus = "Ignotum.elaeocarpus"  )
   n=n+1; f2nog[[n]] = list(family = "Emmotaceae"       , genus = "Ignotum.emmotum"      )
   n=n+1; f2nog[[n]] = list(family = "Erythroxylaceae"  , genus = "Ignotum.erythroxylum" )
   n=n+1; f2nog[[n]] = list(family = "Euphorbiaceae"    , genus = "Ignotum.euphorbia"    )
   n=n+1; f2nog[[n]] = list(family = "Fabaceae"         , genus = "Ignotum.faba"         )
   n=n+1; f2nog[[n]] = list(family = "Goupiaceae"       , genus = "Goupia"               )
   n=n+1; f2nog[[n]] = list(family = "Humiriaceae"      , genus = "Ignotum.humiria"      )
   n=n+1; f2nog[[n]] = list(family = "Hypericaceae"     , genus = "Ignotum.hypericum"    )
   n=n+1; f2nog[[n]] = list(family = "Icacinaceae"      , genus = "Ignotum.icacina"      )
   n=n+1; f2nog[[n]] = list(family = "Ignotaceae"       , genus = "Ignotum"              )
   n=n+1; f2nog[[n]] = list(family = "Lacistemataceae"  , genus = "Ignotum.lacistema"    )
   n=n+1; f2nog[[n]] = list(family = "Lamiaceae"        , genus = "Ignotum.lamium"       )
   n=n+1; f2nog[[n]] = list(family = "Lauraceae"        , genus = "Ignotum.laurus"       )
   n=n+1; f2nog[[n]] = list(family = "Lecythidaceae"    , genus = "Ignotum.lecythis"     )
   n=n+1; f2nog[[n]] = list(family = "Lianaceae"        , genus = "Liana"                )
   n=n+1; f2nog[[n]] = list(family = "Linaceae"         , genus = "Ignotum.linum"        )
   n=n+1; f2nog[[n]] = list(family = "Loganiaceae"      , genus = "Ignotum.logania"      )
   n=n+1; f2nog[[n]] = list(family = "Lythraceae"       , genus = "Ignotum.lythrum"      )
   n=n+1; f2nog[[n]] = list(family = "Malpighiaceae"    , genus = "Ignotum.malpighia"    )
   n=n+1; f2nog[[n]] = list(family = "Malvaceae"        , genus = "Ignotum.malva"        )
   n=n+1; f2nog[[n]] = list(family = "Marcgraviaceae"   , genus = "Ignotum.marcgravia"   )
   n=n+1; f2nog[[n]] = list(family = "Melastomataceae"  , genus = "Ignotum.melastoma"    )
   n=n+1; f2nog[[n]] = list(family = "Meliaceae"        , genus = "Ignotum.melia"        )
   n=n+1; f2nog[[n]] = list(family = "Moraceae"         , genus = "Ignotum.morus"        )
   n=n+1; f2nog[[n]] = list(family = "Myristicaceae"    , genus = "Ignotum.myristica"    )
   n=n+1; f2nog[[n]] = list(family = "Myrtaceae"        , genus = "Ignotum.myrtus"       )
   n=n+1; f2nog[[n]] = list(family = "Nyctaginaceae"    , genus = "Ignotum.nyctaginia"   )
   n=n+1; f2nog[[n]] = list(family = "Ochnaceae"        , genus = "Ignotum.ochna"        )
   n=n+1; f2nog[[n]] = list(family = "Olacaceae"        , genus = "Ignotum.olax"         )
   n=n+1; f2nog[[n]] = list(family = "Opiliaceae"       , genus = "Ignotum.opilia"       )
   n=n+1; f2nog[[n]] = list(family = "Peridiscaceae"    , genus = "Ignotum.peridiscus"   )
   n=n+1; f2nog[[n]] = list(family = "Phyllanthaceae"   , genus = "Ignotum.phyllanthus"  )
   n=n+1; f2nog[[n]] = list(family = "Primulaceae"      , genus = "Ignotum.primula"      )
   n=n+1; f2nog[[n]] = list(family = "Proteaceae"       , genus = "Ignotum.protea"       )
   n=n+1; f2nog[[n]] = list(family = "Polygonaceae"     , genus = "Ignotum.polygonum"    )
   n=n+1; f2nog[[n]] = list(family = "Putranjivaceae"   , genus = "Ignotum.putranjiva"   )
   n=n+1; f2nog[[n]] = list(family = "Rosaceae"         , genus = "Ignotum.rosa"         )
   n=n+1; f2nog[[n]] = list(family = "Rubiaceae"        , genus = "Ignotum.rubia"        )
   n=n+1; f2nog[[n]] = list(family = "Rutaceae"         , genus = "Ignotum.ruta"         )
   n=n+1; f2nog[[n]] = list(family = "Salicaceae"       , genus = "Ignotum.salix"        )
   n=n+1; f2nog[[n]] = list(family = "Sapindaceae"      , genus = "Ignotum.sapindus"     )
   n=n+1; f2nog[[n]] = list(family = "Sapotaceae"       , genus = "Ignotum.sapota"       )
   n=n+1; f2nog[[n]] = list(family = "Simaroubaceae"    , genus = "Ignotum.simarouba"    )
   n=n+1; f2nog[[n]] = list(family = "Siparunaceae"     , genus = "Ignotum.siparuna"     )
   n=n+1; f2nog[[n]] = list(family = "Solanaceae"       , genus = "Ignotum.solanum"      )
   n=n+1; f2nog[[n]] = list(family = "Stemonuraceae"    , genus = "Ignotum.stemonurus"   )
   n=n+1; f2nog[[n]] = list(family = "Strelitziaceae"   , genus = "Ignotum.strelitzia"   )
   n=n+1; f2nog[[n]] = list(family = "Styracaceae"      , genus = "Ignotum.styrax"       )
   n=n+1; f2nog[[n]] = list(family = "Ulmaceae"         , genus = "Ignotum.ulmus"        )
   n=n+1; f2nog[[n]] = list(family = "Urticaceae"       , genus = "Ignotum.urtica"       )
   n=n+1; f2nog[[n]] = list(family = "Verbenaceae"      , genus = "Ignotum.verbena"      )
   n=n+1; f2nog[[n]] = list(family = "Violaceae"        , genus = "Ignotum.viola"        )
   n=n+1; f2nog[[n]] = list(family = "Vochysiaceae"     , genus = "Ignotum.vochysia"     )
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
            }else{
               if (any(c(FALSE,is.finite(wood$density[iwood])),na.rm=TRUE)){
                  sel       = datum$scientific %in% species[s]
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
            iwood              = intersect( which( wood$genus  %in% this.genus )
                                          , which( wood$family %in% this.family) )
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
   for (f in sequence(nfamilies)){
      if (families[f] != "Ignotaceae"){
         if (verbose) cat("   - ",f,"/",nfamilies,"  -  ",families[f],"...","\n")

         #----- Get the individuals that belong to this family. ---------------------------#
         ifam  = which (datum$family %in% families[f]  & is.finite(datum$wood.dens))
         imiss = which (datum$family %in% families[f]  & is.na    (datum$wood.dens))
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
keep.gen.spe.only <<- function(x,out="both"){
   if (is.matrix(x) | is.array(x)){
      ans = apply(X=x,MARGIN=dim(x),FUN=keep.gen.spe.only,out=out)
   }else if (is.list(x)){
      ans = lapply(X=x,FUN=keep.gen.spe.only,out=out)
   }else if (is.data.frame(x)){
      ans = sapply(X=x,FUN=keep.gen.spe.only,out=out)
   }else if (length(x) > 1){
      ans = sapply(X=x,FUN=keep.gen.spe.only,out=out)
   }else{
      #----- Remove stuff that is not genus. ----------------------------------------------#
      x    = sub(pattern="cf\\."  ,replacement=""           ,ignore.case=TRUE,x=x)
      x    = sub(pattern="cf\\ "  ,replacement=""           ,ignore.case=TRUE,x=x)
      x    = sub(pattern="Ind\\ " ,replacement="deleteme\\ ",ignore.case=TRUE,x=x)
      #------------------------------------------------------------------------------------#

      #----- Split words. -----------------------------------------------------------------#
      ans  = unlist(strsplit(unlist(c(tolower(x))),split=" "))
      nans = length(ans)
      #------------------------------------------------------------------------------------#

      #----- Simplify option so just the first letter matters. ----------------------------#
      out1 = substring(tolower(out),1,1)
      #------------------------------------------------------------------------------------#


      #----- Decide whether to keep the genus, species, or both. --------------------------#
      gen = NA_character_
      spe = NA_character_
      if (nans > 0){
         gen = capwords(ans[1],strict=TRUE)
         if (nans > 1 && (! gen %in% "Deleteme")){
            spe = tolower(ans[2])
            if (spe %in% c("sp.","sp.1","ni","ind")) spe = NA_character_
         }else{
            gen = NA_character_
         }#end if
      }#end if
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Select the proper output.                                                      #
      #------------------------------------------------------------------------------------#
      if (out1 == "g"){
         ans = gen
      }else if (out1 == "s"){
         ans = spe
      }else if (out1 == "b"){
         ans = ifelse(is.na(spe),gen,paste(gen,spe,sep=" "))
      }else{
         stop (paste(" Invalid out (",out,")",sep=""))
      }#end if
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#
   return(ans)
}#end function keep.gen.spe.only
#==========================================================================================#
#==========================================================================================#
