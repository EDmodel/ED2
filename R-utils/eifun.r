#------------------------------------------------------------------------------------------#
#     EIFUN -- This function computes the exponential integral function, defined by        #
#                                                                                          #
#                           _                                                              #
#                      x   / \  exp(t)                                                     #
#              Ei(x) =     |   -------- dt                                                 #
#                      0 \_/      t                                                        #
#------------------------------------------------------------------------------------------#
eifun = function(x){

   euler    = 0.57721566490153286060651209008240243104215933593992
   tiny.off = 1.0e-30
   eps      = 1.0e-7

   #---------------------------------------------------------------------------------------#
   #    Initialise ei with the part that does not vary.  Notice that this will set the     #
   # result for any x value to infinity, but that's fine because the result should be      #
   # infinity anyway.                                                                      #
   #---------------------------------------------------------------------------------------#
   ei = euler + log(abs(x))
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Iterate until convergence.                                                         #
   #---------------------------------------------------------------------------------------#
   powermeth   = abs(x) <= 15
   #----- We only solve the numbers that are enough large. --------------------------------#
   large       = abs(x) > tiny.off & powermeth
   #----- Initialise iteration aux vars, keeping the tiny guys out of the process. --------#
   iter         = 0
   fact         = x
   term         = x
   fact[large]  = 1
   fact[!large] = 0
   while (any(large)){
      iter = iter + 1
      fact[large] = fact[large] * x[large] / iter
      large       = abs(fact) >= tiny.off
      term[large] = fact[large] / iter
      large       = large & abs(term) >= eps * abs(ei)
      ei[large]   = ei[large] + term[large]
   }#end while
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Use asymptotic expansion for large numbers.  This may diverge so we must check.    #
   #---------------------------------------------------------------------------------------#
   asymptot    = ! powermeth
   large       = asymptot
   #----- Initialise the sum for the terms that will play. --------------------------------#
   term        = x
   term[large] = 1
   ei[large]   = 1
   iter = 0
   while (any(large)){
      iter        = iter + 1
      prev        = term
      term[large] = term[large] * iter / x[large]
      
      #----- Check convergence. -----------------------------------------------------------#
      large    = large & abs(term) >= eps
      diverge  = large & abs(term) >  abs(prev)
      converge = large & abs(term) <= abs(prev)
      #------------------------------------------------------------------------------------#
      #   Remove the previous iteration for the diverging points, and add for the good     #
      # ones.                                                                              #
      #------------------------------------------------------------------------------------#
      ei[diverge]  = ei[diverge] - prev[diverge]
      ei[converge] = ei[converge] + term[converge]
      
      large   = converge
      
      nlarge  = sum(large)
      if (nlarge > 0){
         maxterm = max(abs(term[large]))
      }else{
         maxterm = 0.
      }#end if
   }#end while
   ei[asymptot] = ei[asymptot] * exp(x[asymptot]) / x[asymptot]

   return(ei)
}#end function
