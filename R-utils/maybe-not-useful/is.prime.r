#------------------------------------------------------------------------------------------#
# Function is.prime                                                                        #
# Developed by Marcos Longo - EPS/Harvard                                                  #
#                                                                                          #
#   This function is among the most stupid I've ever done. It checks whether a number is   #
# prime or not. The way it does it? It just checks from a list I found using google. It    #
# checks only up to 104729, larger numbers will generate NA                                #
#------------------------------------------------------------------------------------------#
is.prime <- function(n){

if ((n %% 1) > 1.e-6){
  warning("The function is.prime works only for integer numbers.",immediate.=TRUE)
  numpri <- NA
}else if (abs(n) > 104729){
  warning("The function is.prime works only for numbers up to 104,729.",immediate.=TRUE)
  numpri <- NA
}else{ 
  numpri <- abs(n) %in% primes()
} #end if (!is.integer(n))
  return(numpri)
} #end function(is.prime)
#------------------------------------------------------------------------------------------#
