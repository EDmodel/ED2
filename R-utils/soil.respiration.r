#==========================================================================================#
#==========================================================================================#
#     This function computes the root respiration based on ED2.                            #
#------------------------------------------------------------------------------------------#
find.root.respiration <<- function(broot,ipft,soil.temp,slz=-0.1){
   #----- Make sure that PFT and broot have the same size. --------------------------------#
   error = 0
   if ( (length(broot) > 1) && (length(ipft) == 1) ){
      #----- Turn ipft into a vector of the same length as broot. -------------------------#
      zpft = rep(ipft,times=length(broot))
      #------------------------------------------------------------------------------------#
   }else if ( (length(broot) > 1) && (length(ipft) != length(broot)) ){
      #----- Sizes do not match, stop. ----------------------------------------------------#
      cat0("==============================================================")
      cat0(" Size of \"broot\" = ",length(broot),"."                       )
      cat0(" Size of \"ipft\"  = ",length(ipft) ,"."                       )
      cat0(" "                                                             )
      cat0(" Variables \"broot\" and \"ipft\" must have the same size, or" )
      cat0("    \"ipft\" must be a scalar."                                )
      cat0("==============================================================")
      cat0(" "                                                             )
      error = error + 1
      #------------------------------------------------------------------------------------#
   }else{
      #----- Copy ipft to temporary vector. -----------------------------------------------#
      zpft = ipft
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Make sure that soil.temp and slz have the same size. ----------------------------#
   if ( length(soil.temp) != length(slz) ){
      #----- Sizes do not match, stop. ----------------------------------------------------#
      cat0("==============================================================")
      cat0(" Size of \"soil.temp\" = ",length(soil.temp),"."               )
      cat0(" Size of \"slz\"       = ",length(slz) ,"."                    )
      cat0(" "                                                             )
      cat0(" Variables \"soil.temp\" and \"slz\" must have the same size." )
      cat0("==============================================================")
      cat0(" "                                                             )
      error = error + 1
      #------------------------------------------------------------------------------------#
   }#end if
   if (error > 0) stop("Check argument dimensions in function \"root.respiration\".")
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #     Make sure the soil layers have negative depth, and go from deepest to shallowest. #
   #---------------------------------------------------------------------------------------#
   slz        = - abs(slz)
   oslz       = order(slz)
   slz        = slz[oslz]
   dslz       = diff(c(slz,0.))
   soil.depth = sum(dslz)
   soil.temp = soil.temp[oslz]
   #---------------------------------------------------------------------------------------#



   #----- Create matrix with rows being soil depths and columns being root biomass/PFT. ---#
   nslz      = length(slz)
   ncoh      = length(zpft)
   BROOT     = matrix(data=broot    ,nrow=nslz,ncol=ncoh,byrow=TRUE )
   ZPFT      = matrix(data=zpft     ,nrow=nslz,ncol=ncoh,byrow=TRUE )
   DSLZ      = matrix(data=dslz     ,nrow=nslz,ncol=ncoh,byrow=FALSE)
   SOIL.TEMP = matrix(data=soil.temp,nrow=nslz,ncol=ncoh,byrow=FALSE)
   #---------------------------------------------------------------------------------------#

   #----- Find the respiration with no correction. ----------------------------------------#
   RR.NOCORR = collatz(temp=SOIL.TEMP,refval=pft$rr0[ZPFT],q10=pft$rr.q10[ZPFT])
   #---------------------------------------------------------------------------------------#


   #----- Find the temperature function for each layer. -----------------------------------#
   TLOW.FUN  = pft$rr.decay.ecold[ZPFT] * (pft$rr.low.temp[ZPFT] + t00 - SOIL.TEMP)
   TLOW.FUN  = 0. * TLOW.FUN + pmax(lnexp.min,pmin(lnexp.max,TLOW.FUN))
   TLOW.FUN  = 1.0 + exp(TLOW.FUN)
   THIGH.FUN = pft$rr.decay.ecold[ZPFT] * (SOIL.TEMP - pft$rr.high.temp[ZPFT] - t00)
   THIGH.FUN = 0. * THIGH.FUN + pmax(lnexp.min,pmin(lnexp.max,THIGH.FUN))
   THIGH.FUN = 1.0 + exp(THIGH.FUN)
   #---------------------------------------------------------------------------------------#


   #---- Find root respiration from each layer. -------------------------------------------#
   RR.LYR = BROOT * RR.NOCORR / TLOW.FUN / THIGH.FUN * DSLZ / soil.depth
   ans    = colSums(RR.LYR)
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end find.root.respiration
