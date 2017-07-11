#==========================================================================================#
#==========================================================================================#
#      This function predicts the mortality rate as a function of wood density using the   #
# linear fit proposed by                                                                   #
#                                                                                          #
# Kraft, N. J. B., M. R. Metz, R. S. Condit, and J. Chave. The relationship between wood   #
#     density and mortality in a global tropical forest data set. New Phytol.,             #
#     188(4):1124-1136, Dec 2010. doi:10.1111/j.1469-8137.2010.03444.x                     #
#                                                                                          #
# x is the set of parameters that are optimised for the dataset.                           #
#------------------------------------------------------------------------------------------#
mort.dens.kraft <<- function(x,datum){
   a    = x
   mort = a[1] * (datum$wd - mean.wood.dens.global) + a[2]
   return(mort)
}#end function
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#       This function finds the sum of the squares, which is the log likelihood if we      #
# asume the errors to be independent and log-normally distributed.                         #
#==========================================================================================#
#==========================================================================================#
mort.kraft.support <<- function(x,datum){

   mort.try   = mort.dens.kraft(x,datum)
   ndat       = nrow(datum)
   residual   = log(mort.try) - log(datum$mort)
   sigma      = sqrt(sum(residual)^2 / (ndat - 1))
   support    = - ndat * log(sigma) - sum(residual^2/(2*sigma^2),na.rm=TRUE)
   return(support)
}#end function mort.kraft.support
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#       This function optimises the mortality rate as a function of wood density.  The     #
# default first guess are the coefficients from Yasuni.                                    #
#==========================================================================================#
#==========================================================================================#
mort.kraft.optim <<- function(datum,x.1st=NULL){

   #----- Keep only the lines with information on both density and mortality. -------------#
   sel     = is.finite(datum$mort) & is.finite(datum$wd)
   data.in = datum[sel,]
   n.dat   = nrow(data.in)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    First guess: if none is given, assume non-informative, so 0 slope and the inter-   #
   # cept being the mean.                                                                  #
   #---------------------------------------------------------------------------------------#
   if(is.null(x.1st)) x.1st = c(0,mean(data.in$mort))
   #---------------------------------------------------------------------------------------#


   #----- Run the optimisation using log-normal distribution of residuals. ----------------#
   opt  = optim( par         = x.1st
               , fn          = mort.kraft.support
               , datum       = data.in
               , control     = list( trace   = verbose
                                   , fnscale = -1
                                   , maxit   = 20000
                                   , REPORT  = 1
                                   )#end list
               , hessian     = TRUE
               )#end optim
   if (opt$convergence != 0)  stop (" Solution of mort.kraft.optim didn't converge...")
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Save the output from the optimiser.                                               #
   #---------------------------------------------------------------------------------------#
   ans               = list()
   ans$hessian       = opt$hessian
   ans$df            = n.dat - 2
   ans$coefficients  = opt$par
   ans$std.err       = sqrt(diag(solve(-ans$hessian)))
   ans$t.value       = ans$coefficients / ans$std.err
   ans$p.value       = 2.0 * pt(-abs(ans$t.value),df=ans$df)
   ans$first.guess   = x.1st
   ans$support       = opt$value
   #----- Also save the fitted values and the residuals. ----------------------------------#
   ans$fitted.values = mort.dens.kraft(x=ans$coefficients,datum=datum)
   ans$residuals     = log(ans$fitted.values) - log(datum$mort)
   #---- Estimate the global standard error and R2. ---------------------------------------#
   ss.err            = sum(ans$residuals^2,na.rm=TRUE)
   df.err            = ans$df
   mean.y            = mean(log(datum$mort),na.rm=TRUE)
   ss.tot            = sum((log(datum$mort)-mean.y)^2,na.rm=TRUE)
   df.tot            = n.dat - 1
   ans$r.square      = 1.0 - ss.err * df.tot / ( ss.tot * df.err )
   ans$sigma         = sqrt(ss.err / ans$df)
   ans$ab.abs        = c(opt$par[1],opt$par[2]-opt$par[1]*mean.wood.dens.global)
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end mort.kraft
#==========================================================================================#
#==========================================================================================#
