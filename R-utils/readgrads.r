#------------------------------------------------------------------------------------------#
# Function readgrads                                                                       #
# Developed by Marcos Longo - EPS/Harvard University                                       #
# This function extracts the variablem, within the asked interval, using the processed ctl #
# information.                                                                             #
#------------------------------------------------------------------------------------------#
readgrads <<- function(vari,info,coord="grid",xlim=NA,ylim=NA,zlim=NA,tlim=NA){

   #----- Assign the number of variables. -------------------------------------------------#
   nvars = length(vari)

   #----- Check whether all requested variables are present. ------------------------------#
   isthere = tolower(vari) %in% info$varname
   if (! all(isthere)){
     cat("Variable ",vari[! isthere]," was not found in the ctl file "
        ,info$ctl,"...",sep="","\n")
     stop("The variable list has invalid variables...")
   }else{
     varid = match(tolower(vari),info$varname)
   } #end if (! tolower(vari) in info$varname)


   #----- Find the boundaries of the subset to be retrieved. ------------------------------#
   xlim = sort(xlim)
   ylim = sort(ylim)
   zlim = sort(zlim)
   tlim = sort(tlim)
   if (tolower(coord) == "grid"){
     options(warn=-1)
     ta = max(tlim[1],1        ,na.rm=TRUE)
     tz = min(tlim[2],info$tmax,na.rm=TRUE)
     xa = max(xlim[1],1        ,na.rm=TRUE)
     xz = min(xlim[2],info$xmax,na.rm=TRUE)
     ya = max(ylim[1],1        ,na.rm=TRUE)
     yz = min(ylim[2],info$ymax,na.rm=TRUE)
     za = max(zlim[1],1        ,na.rm=TRUE)

     #-------------------------------------------------------------------------------------#
     #     The number of levels may vary for each variable.  Ensure that it never exceeds  #
     # the maximum number of levels for each variable.                                     #
     #-------------------------------------------------------------------------------------#
     zz = info$zmax[varid]
     zz[zz > zlim[2]] = zlim[2]
     options(warn=0)
   }else if (tolower(coord) == "coord"){
     options(warn=-1)
     ta = max(1,which.min(abs(info$gtime-tlim[1]),na.rm=TRUE))
     tz = min(tmax,which.min(abs(info$gtime-tlim[2]),na.rm=TRUE))

     xa = max(1,which.min(abs(info$glon-xlim[1]),na.rm=TRUE))
     xz = min(xmax,which.min(abs(info$glon-xlim[2]),na.rm=TRUE))
     ya = max(1,which.min(abs(info$glat-ylim[1]),na.rm=TRUE))
     yz = min(ymax,which.min(abs(info$glat-ylim[2]),na.rm=TRUE))
     za = max(1,which.min(abs(info$glev-zlim[1]),na.rm=TRUE))

     #-------------------------------------------------------------------------------------#
     #     The number of levels may vary for each variable.  Ensure that it never exceeds  #
     # the maximum number of levels for each variable.                                     #
     #-------------------------------------------------------------------------------------#
     ze          = which.min(abs(info$glev-zlim[2]),na.rm=TRUE)
     zz          = info$zmax[varid]
     zz[zz > ze] = ze

     options(warn=0)
   }else{
     stop(paste("Error in coord! The value",coord,"is invalid!",
                "Choose 'grid' for grid points or 'coord' for grid coordinates"))
   } #end if (tolower(coord) == "grid")

   #---------------------------------------------------------------------------------------#
   #     Initialise the output, which will contain the data and information about the      #
   # boundaries of the subset.                                                             #
   #---------------------------------------------------------------------------------------#
   outd = list()
   outd$tmax  = (tz-ta)+1
   outd$zmax  = (zz-za)+1
   outd$xmax  = (xz-xa)+1
   outd$ymax  = (yz-ya)+1
   outd$gtime = info$gtime[ta:tz]
   outd$glev  = info$glev[za:max(zz)]
   outd$glon  = info$glon[xa:xz]
   outd$glat  = info$glat[ya:yz]


   #---------------------------------------------------------------------------------------#
   #     Here we must consider the template and non-template cases slightly differently.   #
   #---------------------------------------------------------------------------------------#
   if (info$template){
      #------------------------------------------------------------------------------------#
      #     Template case.  Because we have a single variable per time, we can assign the  #
      # variables with the right size.  We first assign a template variable, and assign it #
      # to each variable.                                                                  #
      #------------------------------------------------------------------------------------#
      for (v in sequence(nvars)){
         data  = array(NA,c(length(ta:tz),length(za:zz[v]),length(xa:xz),length(ya:yz)))
         assign(x=vari[v],value=data)
         rm(data)
      } #end for (v in 1:nvars)


      #------------------------------------------------------------------------------------#
      #     Now we loop over the times, and for each time we retrieve the variable.  Some  #
      # files may be missing, and if this is the case, we skip the time, so the matrix     #
      # will remain with NA at that time.                                                  #
      #------------------------------------------------------------------------------------#
      tt = 0
      for (tabs in ta:tz){
         tt = tt + 1

         #----- Create the file name, and supported compressed files. ---------------------#
         this.file     = info$binary[tabs]
         this.file.gz  = paste(this.file,"gz" ,sep=".")
         this.file.bz2 = paste(this.file,"bz2",sep=".")
         #---------------------------------------------------------------------------------#

         #----- Set up a temporary file, and make sure it doesn't exist. ------------------#
         temp.file = file.path(tempdir(),basename(this.file))
         warn.orig = getOption("warn")
         options(warn=-1)
         file.remove(temp.file)
         options(warn=warn.orig)
         #---------------------------------------------------------------------------------#



         #----- Check whether the file (or the compressed file) exists. -------------------#
         if (file.exists(this.file)){
            dummy    = file.copy(from=this.file,to=temp.file)
            read.now = TRUE
         }else if (file.exists(this.file.gz)){
            dummy    = gunzip(filename=this.file.gz,destname=temp.file,remove=FALSE)
            read.now = TRUE
         }else if (file.exists(this.file.bz2)){
            dummy    = bunzip2(filename=this.file.gz,destname=temp.file,remove=FALSE)
            read.now = TRUE
         }else{
            read.now = FALSE
         }#end if
         #---------------------------------------------------------------------------------#


         if (read.now){
            cat("[+] Reading data from ",basename(info$binary[tabs]),"...",sep="","\n")

            mybin    = file(description=temp.file,open="rb")

            #----- Find the range of the full domain for this time and variable. ----------#
            xmx   = info$xmax
            ymx   = info$ymax
            xr    = 1:info$xmax
            yr    = 1:info$ymax
            nzall = sum(info$zmax)
            npts  = xmx*ymx*nzall
            aux   = readBin(mybin,n=npts,what="numeric",size=info$size,endian=info$endian)

            for (v in sequence(nvars)){
               zmx = info$zmax[varid[v]]
               zr  = 1:info$zmax[varid[v]]
               pr  = info$posica[varid[v]]:info$posicz[varid[v]]

               #---------------------------------------------------------------------------#
               #     Initialise the scratch data array that will temporarily hold the      #
               # dataset for this time and this variable.                                  #
               #---------------------------------------------------------------------------#
               thisdata = array(NA,c(xmx,ymx,zmx))
               thisdata[xr,yr,zr] = aux[pr]

               #----- Permute the dimensions so x and y go to the right. ------------------#
               thisdata = aperm(a=thisdata,perm=c(3,1,2))

               #----- Assign NA to missing data. ------------------------------------------#
               thisdata[abs((thisdata-(info$undef))/info$undef) < eps()] = NA

               dmin = min(thisdata,na.rm=TRUE)
               dmax = max(thisdata,na.rm=TRUE)
               cat("   [-] Reading variable ",vari[v],"...",sep="","\n")

               #---------------------------------------------------------------------------#
               #     Copy the variable to the scratch, copy the subset to this scratch,    #
               # and copy the scratch back to the variable.                                #
               #---------------------------------------------------------------------------#
               data        = get(vari[v])
               data[tt,,,] = thisdata[za:zz[v],xa:xz,ya:yz]
               assign(x=vari[v],value=data)

               #----- Free the memory used to store the scratch variables. ----------------#
               rm(thisdata,data)
            }#end for (v in 1:nvars)
            rm(aux)
            close(mybin)
            file.remove(temp.file)
         }else{
            warning(paste("Skipping the non-existent file",basename(info$binary[tabs]),
                          "and adding NA..."))
         } #end if
      }#end for (tt in ta:tz)
   }else{
      #------------------------------------------------------------------------------------#
      #     Non-template case.  We will now  assign scratch variables with the entire      #
      # domain first, and crop it to the subset afterwards.  Here because we only need to  #
      # open one file, we loop the variable outside the time loop.                         #
      #------------------------------------------------------------------------------------#

      #----- Create the file name, and supported compressed files. ------------------------#
      this.file     = info$binary[1]
      this.file.gz  = paste(this.file,"gz" ,sep=".")
      this.file.bz2 = paste(this.file,"bz2",sep=".")
      #------------------------------------------------------------------------------------#

      #----- Set up a temporary file, and make sure it doesn't exist. ---------------------#
      temp.file = file.path(tempdir(),basename(this.file))
      warn.orig = getOption("warn")
      options(warn=-1)
      file.remove(temp.file)
      options(warn=warn.orig)
      #------------------------------------------------------------------------------------#



      #----- Check whether the file (or the compressed file) exists. ----------------------#
      if (file.exists(this.file)){
         dummy    = file.copy(from=this.file,to=temp.file)
         read.now = TRUE
      }else if (file.exists(this.file.gz)){
         dummy    = gunzip(filename=this.file.gz,destname=temp.file,remove=FALSE)
         read.now = TRUE
      }else if (file.exists(this.file.bz2)){
         dummy    = bunzip2(filename=this.file.gz,destname=temp.file,remove=FALSE)
         read.now = TRUE
      }else{
         stop(paste("File",basename(info$binary[1]),"(or compressed form) doesn't exist!"))
      }#end if
      #------------------------------------------------------------------------------------#



      mybin = file(description=temp.file,open="rb")

      #----- Find the range of the full domain for this variable. -------------------------#
      xmx   = info$xmax
      ymx   = info$ymax
      nzall = sum(info$zmax)
      tmx   = info$tmax
      #------------------------------------------------------------------------------------#



      #----- Read all variables. ----------------------------------------------------------#
      aux   = readBin(mybin,n=xmx*ymx*nzall*tmx,what="numeric",size=info$size,
                      endian=info$endian)
      #------------------------------------------------------------------------------------#



      #----- Close connection only after all variables have been read. --------------------#
      close(mybin)
      file.remove(temp.file)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Loop over variables.                                                           #
      #------------------------------------------------------------------------------------#
      for (v in sequence(nvars)){
         cat("[+] Reading variable ",vari[v],"...",sep="","\n")
         zmx   = info$zmax[varid[v]]
         xr    = sequence(info$xmax)
         yr    = sequence(info$ymax)
         zr    = sequence(info$zmax[varid[v]])

         #----- Find the full dimensions of the output array of this variable. ------------#
         nxsub = xz    - xa + 1
         nysub = yz    - ya + 1
         nzsub = zz[v] - za + 1
         ntsub = tz    - ta + 1

         #----- Assign a temporary aray with the full dimension. --------------------------#
         datascr = array(NA,c(ntsub,xmx,ymx,zmx))

         #----- Now we loop over the time of interest and extract the data. ---------------#
         tt = 0
         for (tabs in ta:tz){
            pr          = info$posica[tabs,varid[v]]:info$posicz[tabs,varid[v]]
            tt          = tt + 1
            thisdata    = aux[pr]
            thisdata[abs((thisdata-info$undef)/info$undef) < eps()] = NA
            datascr[tt,,,] = thisdata
         }#end for (tt in ta:tz)

         #----- Define the output data, then copy the sub-set of the scratch there... -----#
         data = array(NA,c(ntsub,nxsub,nysub,nzsub))
         data[1:ntsub,1:nxsub,1:nysub,1:nzsub] = datascr[1:ntsub,xa:xz,ya:yz,za:zz[v]]
         rm(thisdata,datascr)
         #----- Rearrange the dimensions so t is to the left, and x and y to the right. ---#
         data = aperm(a=data,perm=c(1,4,2,3))
         #----- We now copy data to the array with the variable name. ---------------------#
         assign(x=vari[v],value=data)
      }#end for (v in 1:nvars)
      rm (aux)
   }#end if(info$template)
   #---------------------------------------------------------------------------------------#


   #----- Include the variables to the output list. ----------------------------------------#
   namesoutd = names(outd)
   for (v in 1:nvars){
      data = get(x=vari[v])
      outd$data   = data
      rm(data)
      namesoutd   = c(namesoutd,vari[v])
      names(outd) = namesoutd
   }#end for (v in 1:nvars)
   #---------------------------------------------------------------------------------------#
  return(outd)
}#end function readgrads
#------------------------------------------------------------------------------------------#
