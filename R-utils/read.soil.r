#==========================================================================================#
#==========================================================================================#
#      This function reads in the soil texture maps.                                       #
#------------------------------------------------------------------------------------------#
read.soil <<- function(stext.type,mapdir=getwd(),fracexp=0.40){
   #------ Define some settings depending on the data set. --------------------------------#
   if (stext.type == "quesada"){
      na.strings = "-9999"
      byrow      = TRUE
      yrev       = TRUE
      reclass    = NULL
      desc       = "Quesada"
      sfile      = file.path(mapdir,"soilmap","quesada.txt")
   }else if (stext.type == "radam"){
      na.strings = "99"
      byrow      = FALSE
      yrev       = TRUE
      reclass    = c( 7, 9, 5, 2, 8, 9, 8, 9,16, 2
                    , 6, 1, 2, 2, 8, 7, 5, 3, 3
                    )#end c
      desc       = "RADAM (INPE)"
      sfile      = file.path(mapdir,"soilmap","radam.txt")
   }else if (stext.type == "igbp"){
      na.strings = "-9999"
      byrow      = TRUE
      yrev       = TRUE
      reclass    = NULL
      desc       = "IGBP"
      sfile      = file.path(mapdir,"soilmap","igbp.txt")
   }else if (stext.type == "fao"){
      na.strings = "NA"
      byrow      = FALSE
      yrev       = FALSE
      reclass    = c( 6, 8, 4, 7, 7, 8,16, 4, 4, 4
                    , 7, 4, 4, 4, 8, 4, 8, 4, 4, 8
                    , 4, 2, 4, 4, 4, 4, 6, 8, 8, 8
                    , 4, 8, 8, 2, 6, 4, 7, 4, 4, 3
                    , 4, 6, 7, 4, 4, 4, 4, 4, 4, 4
                    , 4, 4, 4, 4, 4, 4, 2, 4, 4, 2
                    , 4, 3, 4, 2, 7, 6, 4, 4, 6, 8
                    , 8, 7, 2, 5, 4, 5, 6, 6, 4, 2
                    , 2, 2, 4, 6, 2, 2, 2, 2, 2, 4
                    , 2, 2, 2, 4, 2, 4, 3, 6, 2, 7
                    , 4, 4, 4, 8, 8, 8, 3, 7, 4, 4
                    , 4, 3, 6, 4, 2, 4, 4, 4, 2, 2
                    , 2, 4, 6, 4, 4, 7, 7, 6, 3, 2
                    , 2, 6, 6, 6)
      desc       = "Zobler (FAO)"
      sfile      = file.path(mapdir,"soilmap","fao.txt")
   }else{
      stop(paste(" - Invalid class: ",which.class,"!!!"),sep="")
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Read the file. ------------------------------------------------------------------#
   myall = scan(file=sfile,na.strings=na.strings,quiet=TRUE)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Retrieve the coordinate information.                                             #
   #---------------------------------------------------------------------------------------#
   nlon = myall[1]
   lon0 = myall[2]
   dlon = myall[3]
   nlat = myall[4]
   lat0 = myall[5]
   dlat = myall[6]
   #---------------------------------------------------------------------------------------#



   #------ Save the rest of the data. -----------------------------------------------------#
   mydat = myall[7:length(myall)]
   sel = is.finite(mydat)
   if (! is.null(reclass)){
      mydat[sel] = reclass[mydat[sel]]
   }#end if
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #      Create a list that contains all the soil texture information.                    #
   #---------------------------------------------------------------------------------------#
   myres        = list()
   myres$name   = stext.type
   myres$desc   = desc
   myres$source = sfile
   myres$lon    = lon0 + seq(from=0,to=nlon-1,by=1) * dlon
   myres$lat    = lat0 + seq(from=0,to=nlat-1,by=1) * dlat
   myres$nlon   = nlon
   myres$nlat   = nlat
   myres$dlon   = dlon
   myres$dlat   = dlat
   myres$stext  = matrix(mydat,ncol=nlat,nrow=nlon,byrow=byrow)
   if (yrev) myres$stext[,1:nlat] = myres$stext[,rev(1:nlat)]

   limlon       = pretty.xylim(u=myres$lon,fracexp=0.0    ,is.log=FALSE)
   limlat       = pretty.xylim(u=myres$lat,fracexp=fracexp,is.log=FALSE)
   myres$size   = plotsize(proje=TRUE,limlon=limlon,limlat=limlat,paper="letter")
   #---------------------------------------------------------------------------------------#


   #----- Return the list. ----------------------------------------------------------------#
   return(myres)
   #---------------------------------------------------------------------------------------#
}#end function
#==========================================================================================#
#==========================================================================================#
