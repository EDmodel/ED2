#------------------------------------------------------------------------------------------#
# Function eps                                                                             #
# Developed by Marcos Longo - EPS/Harvard University                                       #
#                                                                                          #
#   This function simply gives the machine epsilon for single precision                    #
#------------------------------------------------------------------------------------------#
 eps <<- function(){
   ans = 2^(-24)
   return(ans)
 }
#------------------------------------------------------------------------------------------#
