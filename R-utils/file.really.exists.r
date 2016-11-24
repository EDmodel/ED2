############################################################################################
# Fuction file.really.exists                                                               #
#                                                                                          #
#   Marcos Longo - EPS/Harvard University.                                                 #
#   Cambridge, October 4, 2006.                                                            #
#                                                                                          #
#   This function tests whether a file exists or not (as file.exists does), but it also    #
# consider the possibility that the file can be large. The original R function usually     #
# fails in this case.                                                                      #
############################################################################################
file.really.exists <- function(filename){

###### Keeping it simple: if file.exists returns TRUE, then it's TRUE ######################
  foundit <- file.exists(filename)

###### In case foundit is FALSE, but I am in a Linux environment, I double check. ##########
  if (! foundit && Sys.info()["sysname"] == "Linux"){
    foundit <- as.numeric(system(command=paste("ls -1",filename,"2> /dev/null | wc -l"),
                                 intern=TRUE)) > 0
  } #end if (! foundit)
  return (foundit)
} #end function file.really.exists
############################################################################################
