#------------------------------------------------------------------------------------------#
# Function readctl                                                                         #
# Developed by Marcos Longo - EPS/Harvard University                                       #
#                                                                                          #
#   This function reads the ctl file, and writes all needed information into a variable of #
# kind list (info).                                                                        #
#------------------------------------------------------------------------------------------#
readctl = function(ctlfile,size=4){
  library(chron)
  ctl0 = scan(file=ctlfile,what="character",sep=(c("\n")),multi.line=TRUE)
  ctl  = strsplit(x=ctl0,split=" ")
  lmax = length(ctl)
  for (l in 1:lmax) ctl[[l]] = ctl[[l]][ctl[[l]] != ""]

  template = FALSE
  endian = .Platform$endian
  for (l in 1:lmax){
     what = tolower(ctl[[l]][1])

     #----- Reading the binary file name --------------------------------------------------#
     if (what == "dset") binary = ctl[[l]][2]

     #----- Reading the undefined value ---------------------------------------------------#
     if (what == "undef") undef = as.numeric(ctl[[l]][2])

     #----- Reading the options line ------------------------------------------------------#
     if (what == "options"){
        if (any(tolower(ctl[[l]]) == "template")) template = TRUE
        if (any(tolower(ctl[[l]]) == "little_endian")) endian = "little"
        if (any(tolower(ctl[[l]]) == "big_endian"))    endian = "big"
     }#end if (what == "options")

     #----- Reading the xdef line ---------------------------------------------------------#
     if (what == "xdef") {
        glon = gridp(ctl[[l]],ctl[[l+1]])
        xmax = length(glon)
     } #end if (what == "xdef")

     #----- Reading the ydef line ---------------------------------------------------------#
     if (what == "ydef") {
        glat = gridp(ctl[[l]],ctl[[l+1]])
        ymax = length(glat)
     } #end if (what == "ydef")

     #----- Reading the zdef line ---------------------------------------------------------#
     if (what == "zdef") {
        glev = gridp(ctl[[l]],ctl[[l+1]])
        lmax = length(glev)
     } #end if (what == "zdef")

     #----- Reading the tdef line ---------------------------------------------------------#
     if (what == "tdef") {
        gtime = gridt(as.numeric(ctl[[l]][2]),ctl[[l]][4],ctl[[l]][5])
        tmax  = length(gtime)
     } #end if (what == "tdef")

     #----- Finding the variable names and their sizes ------------------------------------#
     if (what == "vars"){
        nvars = as.numeric(ctl[[l]][2])
        varname     = NULL
        zmax        = NULL
        description = NULL
        #----- Extracting the variable name and size --------------------------------------#
        for (v in 1:nvars){
           #----- Switching underscores by dots, so R won't mess up -----------------------#
           aux          = ctl[[l+v]][1]
           saux         = strsplit(x=aux,split="")[[1]]
           uscore       = which(saux == "_")
           saux[uscore] = "."
           aux          = paste(saux,collapse="")

           varname = c(varname,aux)
           zmax    = c(zmax,as.numeric(ctl[[l+v]][2]))
           #------ Retrieve the original line, since it has the description ---------------#
           aux = ctl0[[l+v]][1]
           saux         = strsplit(aux,split="")[[1]]
           uscore       = which(saux == "_")
           saux[uscore] = "."
           aux          = paste(saux,collapse="")
           description = c(description,aux)
        } #end for (v in 1:nvars)
        zmax[zmax == 0] = 1
        vmax = length(varname)
     } #end if (what == "vars")
  } #end for (l in 1:lmax)

  #----- Defining the binary name (or names if it is a template) --------------------------#
  path      = paste(dirname(ctlfile),"/",sep="")
  psplit    = strsplit(x=binary,split="")[[1]]
  ishat     = "^" %in% psplit
  hatp      = which(psplit == "^")[1]
  fullpath  = all(substring(binary,1,1) == "/")
  if (ishat){
     binary = paste(path,substring(binary,hatp+1),sep="/")
  }else if (! fullpath){
     binary = paste(path,binary,sep="/")
  }#end fi

  if (template) binary = bintemp(binary,gtime)

#----- Reading the variables --------------------------------------------------------------#
  pos = xmax*ymax*zmax
  if (template){
    posicz = cumsum(pos)
    posica = 1; for (p in 2:vmax) posica = c(posica,posicz[p-1]+1)
  }else{
    posicz = matrix(NA,ncol=vmax,nrow=tmax,dimnames=list(gtime,varname))
    posica = matrix(NA,ncol=vmax,nrow=tmax,dimnames=list(gtime,varname))
    for (tt in 1:tmax) {
      posicz[tt,] = (tt-1)*sum(pos)+cumsum(pos)
      posica[tt,1] = max(0,posicz[tt-1,vmax],na.rm=TRUE)+1
      if (vmax > 1) posica[tt,2:vmax] = posicz[tt,1:(vmax-1)]+1
    } #end for (tt in 1:tmax)
  } #end if (template)

  info = NULL
  info$ctl         = ctlfile
  info$binary      = binary
  info$size        = size
  info$template    = template
  info$undef       = undef
  info$endian      = endian
  info$glon        = glon
  info$glat        = glat
  info$glev        = glev
  info$gtime       = gtime
  info$xmax        = xmax
  info$ymax        = ymax
  info$lmax        = lmax
  info$tmax        = tmax
  info$varname     = varname
  info$zmax        = zmax
  info$description = description
  info$posica      = posica
  info$posicz      = posicz
  return(info)
} #end function readctl
#------------------------------------------------------------------------------------------#
