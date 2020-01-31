#==========================================================================================#
#==========================================================================================#
#     This function aggregates error statistics when each entry is a collection of         #
# residuals, standard deviation, and root mean square error.                               #
#------------------------------------------------------------------------------------------#
error.aggr <<- function(n,bias,sigma,rmse,yobs,na.rm=TRUE){
   #----- Aggregated bias. ----------------------------------------------------------------#
   if (na.rm){
      keep   = is.finite(n) & is.finite(bias) & is.finite(sigma) & is.finite(yobs)
      n      = n     [keep]
      bias   = bias  [keep]
      sigma  = sigma [keep]
      rmse   = rmse  [keep]
      yobs   = yobs  [keep]
   }#end if(na.rm)
   #---------------------------------------------------------------------------------------#


   #----- Find standard deviation of observations. ----------------------------------------#
   o.sigma = sd(yobs)
   #---------------------------------------------------------------------------------------#


   #----- Find aggregated values. ---------------------------------------------------------#
   a.n     = sum(n)
   a.bias  = weighted.mean(x=bias,w=n)
   a.sigma = sqrt(sum((n-1)*sigma^2 + n*(a.bias-bias)^2)/(a.n-1))
   a.rmse  = sqrt(a.bias^2+a.sigma^2)
   a.rsqu  = 1. - a.rmse^2/o.sigma^2
   #---------------------------------------------------------------------------------------#


   #----- Return aggregated errors. -------------------------------------------------------#
   ans = c(n = a.n, bias = a.bias, sigma = a.sigma, rmse = a.rmse, rsqu = a.rsqu)
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end error.aggr
#==========================================================================================#
#==========================================================================================#
