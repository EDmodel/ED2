#==========================================================================================#
#==========================================================================================#
#     This function fits a VPD response curve to given data.                               #
#------------------------------------------------------------------------------------------#
lm.vpd.response   <<- function( vpd
                              , gpp
                              , pred.vpd = seq( from = min(vpd,na.rm=TRUE)
                                              , to   = max(vpd,na.rm=TRUE)
                                              , length.out = 100
                                              )#end seq
                              , ...
                              ){#end function


   #----- Set the input dataset. ----------------------------------------------------------#
   use        = is.finite(vpd) & is.finite(gpp)
   data.in    = data.frame(vpd=vpd[use],gpp=gpp[use])
   data.pred  = data.frame(vpd=pred.vpd)
   n.use      = sum(use)
   n.coeff    = 3
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Main fit.                                                                         #
   #---------------------------------------------------------------------------------------#
   fit.gpp = try(lm(gpp ~ vpd + I(vpd^2), data = data.in,...),silent=TRUE)
   #---------------------------------------------------------------------------------------#


   if ("try-error" %in% is(fit.gpp)){
      #------------------------------------------------------------------------------------#
      #     Make all data sets NA.                                                         #
      #------------------------------------------------------------------------------------#
      coeff.gpp    = rep(NA,times=3)
      std.err.gpp  = rep(NA,times=3)
      tvalue.gpp   = rep(NA,times=3)
      pvalue.gpp   = rep(NA,times=3)
      df.gpp       = NA
      df.tot       = sum(is.finite(gpp)) - 1
      r2.gpp       = NA
      pred.gpp     = NA + pred.vpd
      expected.gpp = NA + pred.vpd
      q025.gpp     = NA + pred.vpd
      q975.gpp     = NA + pred.vpd
      #------------------------------------------------------------------------------------#
   }else{
      #------------------------------------------------------------------------------------#
      #     Save summary and predictions.                                                  #
      #------------------------------------------------------------------------------------#
      summ.gpp = summary(fit.gpp)
      pred.gpp = predict(fit.gpp,newdata=data.pred,se.fit=TRUE,interval="confidence")
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #     Copy data.                                                                     #
      #------------------------------------------------------------------------------------#
      coeff.gpp    = summ.gpp$coefficients[,1]
      std.err.gpp  = summ.gpp$coefficients[,2]
      tvalue.gpp   = summ.gpp$coefficients[,3]
      pvalue.gpp   = summ.gpp$coefficients[,4]
      df.gpp       = summ.gpp$df[2]
      df.tot       = sum(is.finite(gpp)) - 1
      #------------------------------------------------------------------------------------#



      #----- Find the adjusted R2. --------------------------------------------------------#
      residual.gpp   = summ.gpp$residuals
      df.gpp         = n.use - n.coeff
      ss.gpp         = sum(residual.gpp^2,na.rm=TRUE) / df.gpp
      sigma.gpp      = sqrt(ss.gpp)
      ss.tot         = var(data.in$gpp,na.rm=TRUE)
      r2.gpp         = 1.0 - ss.gpp / ss.tot
      #------------------------------------------------------------------------------------#


      #----- Find predicted values and their confidence band. -----------------------------#
      expected.gpp   = pred.gpp$fit[,1]
      q025.gpp       = pred.gpp$fit[,2]
      q975.gpp       = pred.gpp$fit[,3]
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Return answer. ------------------------------------------------------------------#
   ans = list( coefficients = coeff.gpp
             , std.err      = std.err.gpp
             , tvalue       = tvalue.gpp
             , pvalue       = pvalue.gpp
             , r2           = r2.gpp
             , vpd          = pred.vpd
             , gpp          = expected.gpp
             , q025         = q025.gpp
             , q975         = q975.gpp
             )#end list
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function
#==========================================================================================#
#==========================================================================================#
