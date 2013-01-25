infoloc <- function(loctag,here=getwd()){
   loctag = tolower(loctag)

   if       (loctag == "rolim"){
      fullname = "Rolim de Moura, RO"
      iata     = "rmo"
      pathin   = paste(here,loctag,"input",sep="/")
      pathout  = paste(here,loctag,"output",sep="/")
      plotout  = paste(here,loctag,"plots",sep="/")
      prefix   = paste(pathin,iata,sep="/")
      suffix   = ".txt"
   }else if (loctag == "ourop"){
      fullname = "Ouro Preto d'Oeste, RO"
      iata     = "fns"
      pathin   = paste(here,loctag,"input",sep="/")
      pathout  = paste(here,loctag,"output",sep="/")
      plotout  = paste(here,loctag,"plots",sep="/")
      prefix   = paste(pathin,iata,sep="/")
      suffix   = ".txt"
   }else if (loctag == "rebio"){
      fullname = "Rebio Jaru, RO"
      iata     = "rbj"
      pathin   = paste(here,loctag,"input",sep="/")
      pathout  = paste(here,loctag,"output",sep="/")
      plotout  = paste(here,loctag,"plots",sep="/")
      prefix   = paste(pathin,iata,sep="/")
      suffix   = ".txt"
   }else if (loctag == "pvelho"){
      fullname = "Porto Velho, RO"
      iata     = "pvh"
      pathin   = paste(here,loctag,"input",sep="/")
      pathout  = paste(here,loctag,"output",sep="/")
      plotout  = paste(here,loctag,"plots",sep="/")
      prefix   = paste(pathin,iata,sep="/")
      suffix   = ".txt"
   }else if (loctag == "vilhena"){
      fullname = "Vilhena, RO"
      iata     = "bvh"
      pathin   = paste(here,loctag,"input",sep="/")
      pathout  = paste(here,loctag,"output",sep="/")
      plotout  = paste(here,loctag,"plots",sep="/")
      prefix   = paste(pathin,iata,sep="/")
      suffix   = ".txt"
   }else{
      fullname = NA
      iata     = NA
      pathin   = NA
      pathout  = NA
      plotout  = NA
      prefix   = NA
      suffix   = NA
      warning(paste("The place",loctag,"is not recognised by infoloc.r, NAs added..."))
   }#end if



   thisplace = list(fullname=fullname,iata=iata,pathin=pathin,pathout=pathout,
                    plotout=plotout,prefix=prefix,suffix=suffix)
   return(thisplace)
} #end function
