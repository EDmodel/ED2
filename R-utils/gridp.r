#------------------------------------------------------------------------------------------#
# Function gridp                                                                           #
# Developed by Marcos Longo - EPS/Harvard University                                       #
# March 28, 2007
#                                                                                          #
#   This function determine the grid levels of each coordinate.                            #
#------------------------------------------------------------------------------------------#
gridp <- function(ctla,ctlb){
  ngrid <- as.numeric(ctla[2])
  how  <- tolower(ctla[3])
  nstr <- length(ctla)
  if(how == "linear"){
    grid0 <- as.numeric(ctla[4])
    dgrid <- as.numeric(ctla[5])
    gridp <- grid0+(1:ngrid-1)*dgrid         
  }else if(nstr == ngrid+3){
    gridp <- as.numeric(ctla[1:ngrid+3])
  }else{
    gridp <- as.numeric(c(ctla[4:nstr],ctlb[1:(ngrid-nstr+3)]))
  } #end if(how == "linear")
  return(gridp)
} #end function gridp
#------------------------------------------------------------------------------------------#
