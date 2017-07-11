#==========================================================================================#
#==========================================================================================#
#      This function predicts GPP for a given light profile.                               #
#------------------------------------------------------------------------------------------#
predict.gpp.from.par <<- function(x,par.in){
   gpp = x[1] + x[2] * par / (x[3] + par)
   return(gpp)
}#end function
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#       This function finds the sum of the squares, which is the log likelihood if we      #
# asume the errors to be independent and normally distributed (a big assumption).          #
#==========================================================================================#
#==========================================================================================#
rshort.in.wn.support <<- function(x,datum){

   rsbdown.try = predict.rshort.bdown( x        = x
                                     , rad.in   = rshort.in
                                     , atm.prss = atm.prss
                                     , cosz     = cosz
                                     , rad.type = "rshort"
                                     )#end predict.rshort.bdown

   residual.par.full = rsbdown.try$par.full - par.in
   residual.nir.full = rsbdown.try$nir.full - nir.in
   residual.par.diff = rsbdown.try$par.diff - par.diff

   residual          = c(residual.par.full,residual.nir.full,residual.par.diff)
   sigma             = c(sigma.par.full   ,sigma.nir.full   ,sigma.par.diff   )
   chi.square        = sum((residual/sigma)^2,na.rm=TRUE)

   support      = - chi.square
   return(support)
}#end function rlong.in.mmi.lnlike
#==========================================================================================#
#==========================================================================================#
