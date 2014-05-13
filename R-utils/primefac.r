#==========================================================================================#
#==========================================================================================#
# Function decompose                                                                       #
# Developed by Marcos Longo                                                                #
#                                                                                          #
# This function decomposes a number in prime factors, and check the power of each one...   #
#------------------------------------------------------------------------------------------#
primefac <<- function(n,tol=1.e-6){
  #----- Check whether the number is "integer" --------------------------------------------#
  if ((n %% 1) > tol){
     warning("Function decompose works only for integer numbers.",immediate.=TRUE)
  }else if (is.prime(n)){
    #----- If the number is prime, skip all steps and return the prime definition... ------#
    factors  = n
    powers   = 1
    allnumb  = c(1,n)
  }else{
    #----- Select prime number candidates. ------------------------------------------------#
    factors  = primes()[primes() <= n/2]
    #----- Determine how many factors are needed ------------------------------------------#
    cmax  = length(factors)
    fmax  = floor(logb(n,factors))
    powers  = rep(NA,times=cmax)
    #----- Find the maximum power which is still a divisor. -------------------------------#
    for (cc in 1:cmax) powers[cc]  = max(which(n %% factors[cc]^(seq(0,fmax[cc])) < tol))-1
    #----- Eliminate those which maximum power is 0 (non-divisors) ------------------------#
    factors  = factors[powers > 0]
    powers  = powers[powers > 0]
    allpowers  = c(1,n)
    for (cc in 1:length(factors)) allpowers  = c(allpowers,factors[cc]^seq(1:powers[cc]))
    allpowers  = sort(unique(allpowers))
    #----- Do the matrix of permutations --------------------------------------------------#
    aux  = matrixperm(allpowers[allpowers != n])
    candidates  = unique(apply(X=aux,MARGIN=1,FUN=prod))
    #----- Select the outcomes that indeed produce a divisor of n -------------------------#
    alldiv = sort(candidates[n %% candidates < tol])
  } #end if
  result  = list(factors=factors,powers=powers,alldiv=alldiv)
  return(result)
}#end function primefac
#------------------------------------------------------------------------------------------#
