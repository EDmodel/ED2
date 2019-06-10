#==========================================================================================#
#==========================================================================================#
#      This function determines the plot aggregated error due to measurement and           #
# allometry.                                                                               #
#                                                                                          #
# Input:                                                                                   #
# ---------------------------------------------------------------------------------------- #
# datum       - a data frame with all trees in this plot (note, you can't run this for all #
#               plots at once, either use a for loop or mapply). The data frame must       #
#               contain the following variables:                                           #
#               * nplant   - 1/area sampled if tree is alive, zero if tree is dead [1/m2]  #
#               * ntotal   - 1/area sampled                                        [1/m2]  #
#                 For nplant/ntotal, if the sampled area was 2500 m2, the number should be #
#                 0.0004.  In case a differential sampling effort was used, then nplant    #
#                 varies.  For example, if trees with 10 <= DBH < 35 cm were measured in a #
#                 50x5 subplot and trees with DBH >= 35 cm were measured in the entire     #
#                 2500m2 plots, then nplant should be 0.004 for the trees with DBH < 35cm  #
#                 and 0.0004 for trees with DBH >= 35cm.                                   #
#               * AGC      - above-ground carbon                                   [ kgC]  #
#               * ME.AGC   - measurement uncertainty of above-ground carbon        [ kgC]  #
#               * LNAGC    - log of above-ground carbon                                    #
#               * SD.LNAGC - allometry uncertainty, using the log-scale                    #
#               * X        - x position in the plot (used only if epsilon.smp is NULL)     #
#               * Y        - y position in the plot (used only if epsilon.smp is NULL)     #
# xmax        - maximum size along the x axis (used only if epsilon.smp is NULL)           #
# ymax        - maximum size along the y axis (used only if epsilon.smp is NULL)           #
# subalong    - which axis has the subplot                                                 #
# epsilon.smp - in case epsilon.smp is null, the function will try to estimate within-plot #
#               sampling uncertainty.  This is unlikely to work unless you have large      #
#               plots or at the very least have plots without sub-sampling.                #
#               Alternatively, you may provide the number from previous studies.           #
# n.sub       - number of subplot samples                                                  #
# n.real      - number of replicates for estimating uncertainty.  Large numbers            #
#               (10000 or more) are needed for stable results.                             #
# ---------------------------------------------------------------------------------------- #
#                                                                                          #
#                                                                                          #
#                                                                                          #
# ---------------------------------------------------------------------------------------- #
# Output:                                                                                  #
# ---------------------------------------------------------------------------------------- #
# A vector with 8 numbers:                                                                 #
# se.agb.xxxxx - uncertainties for plot estimate of biomass (no standing dead) [kgC/m2]    #
# se.acd.xxxxx - uncertainties for plot estimate of biomass+necromass          [kgC/m2]    #
# se.xxx.measurement - contribution of measurement uncertainty to plot uncertainty         #
# se.xxx.allometry   - contribution of allometry uncertainty to plot uncertainty           #
# se.xxx.sampling    - contribution of sampling uncertainty to plot uncertainty            #
# se.xxx.census      - total uncertainty (combining the three terms above).                #
#------------------------------------------------------------------------------------------#
find.acd.error <<- function(datum,xmax,ymax,subalong=c("x","y"),epsilon.smp=NULL
                           ,n.sub=25,n.real=10000){

   #----- Standardise subalong. -----------------------------------------------------------#
   subalong = match.arg(subalong)
   #---------------------------------------------------------------------------------------#


   #----- Number of data points. ----------------------------------------------------------#
   ndatum = nrow(datum)
   #---------------------------------------------------------------------------------------#


   #----- First, find the combined error. -------------------------------------------------#
   se.use = data.frame( measurement = sqrt(log(1+(datum$ME.AGC / datum$AGC)^2))
                      , allometry   = datum$SD.LNAGC
                      )#end datum
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Initialise answer.                                                               #
   #---------------------------------------------------------------------------------------#
   se.agb        = rep(NA,length(se.use)+2)
   names(se.agb) = c(names(se.use),"sampling","census")
   se.acd = rep(NA,length(se.use)+2)
   names(se.acd) = c(names(se.use),"sampling","census")
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #      Create population matrices.                                                      #
   #---------------------------------------------------------------------------------------#
   NPLANT = matrix( data = rep(datum$nplant,times=n.real), nrow=ndatum,ncol=n.real)
   NTOTAL = matrix( data = rep(datum$ntotal,times=n.real), nrow=ndatum,ncol=n.real)
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #      Initialise answer.  We will create the vector before and loop through the        #
   # errors.                                                                               #
   #---------------------------------------------------------------------------------------#
   for (e in seq_along(se.use)){
      AGC       = matrix( data = rlnorm( n = n.real*ndatum
                                       , meanlog = rep(datum$LNAGC,times=n.real)
                                       , sdlog   = rep(se.use[[e]],times=n.real)
                                       )#end rlnorm
                        , nrow = ndatum
                        , ncol = n.real
                        )#end matrix
      AGB       = colSums(NPLANT * AGC)
      ACD       = colSums(NTOTAL * AGC)
      se.agb[e] = sd(AGB)
      se.acd[e] = sd(ACD)
   }#end for (e in seq_along(se.use))
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Find total biomass, to be used by the within plot sampling error.                #
   #---------------------------------------------------------------------------------------#
   agb.bar = sum(datum$nplant * datum$AGC)
   acd.bar = sum(datum$ntotal * datum$AGC)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Unless epsilon.smp is provided, we estimate the within-plot sampling error.      #
   #---------------------------------------------------------------------------------------#
   if (is.null(epsilon.smp)){
      #------------------------------------------------------------------------------------#
      #      Split domain into smaller subsplots.                                          # 
      #------------------------------------------------------------------------------------#
      if (subalong %in% "x"){
         xbreaks = seq(from=0,to=xmax,length.out=n.sub+1)
         xlwr    =       sqrt(.Machine$double.eps)  * xmax
         xupr    = (1. - sqrt(.Machine$double.eps)) * xmax
         xnow    = pmax(xlwr,pmin(xupr,datum$X))
         isub    = as.integer(cut(xnow,breaks=xbreaks))
      }else{
         ybreaks = seq(from=0,to=ymax,length.out=n.sub+1)
         ylwr    =       sqrt(.Machine$double.eps)  * ymax
         yupr    = (1. - sqrt(.Machine$double.eps)) * ymax
         ynow    = pmax(ylwr,pmin(yupr,datum$Y))
         isub    = as.integer(cut(ynow,breaks=ybreaks))
      }#end if
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Find the AGB/ACD for each subplot.                                            #
      #------------------------------------------------------------------------------------#
      agb.sub          = rep(0.,times=n.sub)
      acd.sub          = rep(0.,times=n.sub)
      agb.tmp          = tapply(X=n.sub*datum$nplant*datum$AGC,INDEX=isub,FUN=sum)
      acd.tmp          = tapply(X=n.sub*datum$ntotal*datum$AGC,INDEX=isub,FUN=sum)
      agb.idx          = as.integer(names(agb.tmp))
      acd.idx          = as.integer(names(acd.tmp))
      agb.sub[agb.idx] = agb.tmp
      acd.sub[acd.idx] = acd.tmp
      #------------------------------------------------------------------------------------#




      #------ Create replicates using bootstrap with replacement. -------------------------#
      AGB = colMeans( matrix( data = sample(agb.sub,size=n.real*n.sub,replace=TRUE)
                            , nrow = n.sub
                            , ncol = n.real
                            )#end matrix
                    )#end colMeans
      ACD = colMeans( matrix( data = sample(acd.sub,size=n.real*n.sub,replace=TRUE)
                            , nrow = n.sub
                            , ncol = n.real
                            )#end matrix
                    )#end colMeans
      se.agb["sampling"] = sd(AGB)
      se.acd["sampling"] = sd(ACD)
      #------------------------------------------------------------------------------------#
   }else{
      #----- Use pre-defined sampling error. ----------------------------------------------#
      se.agb["sampling"] = epsilon.smp * agb.bar
      se.acd["sampling"] = epsilon.smp * acd.bar
      #------------------------------------------------------------------------------------#
   }#end if (is.null(epsilon.smp))
   #---------------------------------------------------------------------------------------#



   #------ Find combined source of errors. ------------------------------------------------#
   se.agb["census"] = sqrt(sum(se.agb[c("measurement","allometry","sampling")]^2))
   se.acd["census"] = sqrt(sum(se.acd[c("measurement","allometry","sampling")]^2))
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Error is the standard deviation of all realisations.                             #
   #---------------------------------------------------------------------------------------#
   ans        = c(se.agb,se.acd)
   names(ans) = c(paste0("se.agb.",names(se.agb)),paste0("se.acd.",names(se.acd)))
   #---------------------------------------------------------------------------------------#


   #----- Return answer. ------------------------------------------------------------------#
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end find.acd.error
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#      This function determines the plot aggregated error due to measurement and           #
# allometry.                                                                               #
#                                                                                          #
# Input:                                                                                   #
# ---------------------------------------------------------------------------------------- #
# datum       - a data frame with all trees in this plot (note, you can't run this for all #
#               plots at once, either use a for loop or mapply). The data frame must       #
#               contain the following variables:                                           #
#               * nplant   - 1/area sampled if tree is alive, zero if tree is dead [1/m2]  #
#                 For nplant/ntotal, if the sampled area was 2500 m2, the number should be #
#                 0.0004.  In case a differential sampling effort was used, then nplant    #
#                 varies.  For example, if trees with 10 <= DBH < 35 cm were measured in a #
#                 50x5 subplot and trees with DBH >= 35 cm were measured in the entire     #
#                 2500m2 plots, then nplant should be 0.004 for the trees with DBH < 35cm  #
#                 and 0.0004 for trees with DBH >= 35cm.                                   #
#               * LA       - leaf area                                          [m2_leaf]  #
#               * AE.AGC   - allometry uncertainty of leaf area                 [m2_leaf]  #
#               * ME.AGC   - measurement uncertainty of leaf area               [m2_leaf]  #
#               * X        - x position in the plot (used only if epsilon.smp is NULL)     #
#               * Y        - y position in the plot (used only if epsilon.smp is NULL)     #
# xmax        - maximum size along the x axis (used only if epsilon.smp is NULL)           #
# ymax        - maximum size along the y axis (used only if epsilon.smp is NULL)           #
# subalong    - which axis has the subplot                                                 #
# epsilon.smp - in case epsilon.smp is null, the function will try to estimate within-plot #
#               sampling uncertainty.  This is unlikely to work unless you have large      #
#               plots or at the very least have plots without sub-sampling.                #
#               Alternatively, you may provide the number from previous studies.           #
# n.sub       - number of subplot samples                                                  #
# n.real      - number of replicates for estimating uncertainty.  Large numbers            #
#               (10000 or more) are needed for stable results.                             #
# ---------------------------------------------------------------------------------------- #
#                                                                                          #
#                                                                                          #
#                                                                                          #
# ---------------------------------------------------------------------------------------- #
# Output:                                                                                  #
# ---------------------------------------------------------------------------------------- #
# A vector with 4 numbers:                                                                 #
# se.lai.xxxxx - uncertainties for plot estimate of biomass (no standing dead) [kgC/m2]    #
# se.lai.measurement - contribution of measurement uncertainty to plot uncertainty         #
# se.lai.allometry   - contribution of allometry uncertainty to plot uncertainty           #
# se.lai.sampling    - contribution of sampling uncertainty to plot uncertainty            #
# se.lai.census      - total uncertainty (combining the three terms above).                #
#------------------------------------------------------------------------------------------#
find.lai.error <<- function(datum,xmax,ymax,subalong=c("x","y"),epsilon.smp=NULL
                           ,n.sub=25,n.real=10000){

   #----- Standardise subalong. -----------------------------------------------------------#
   subalong = match.arg(subalong)
   #---------------------------------------------------------------------------------------#


   #----- Number of data points. ----------------------------------------------------------#
   ndatum = nrow(datum)
   #---------------------------------------------------------------------------------------#


   #----- First, find the combined error. -------------------------------------------------#
   se.use             = data.frame( measurement = sqrt(log(1+(datum$ME.LA / datum$LA)^2))
                                  , allometry   = datum$AE.LA
                                  )#end datum
   se.use$measurement = ifelse(is.finite(se.use$measurement),se.use$measurement,NA_real_)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Initialise answer.                                                               #
   #---------------------------------------------------------------------------------------#
   se.lai        = rep(NA_real_,length(se.use)+2)
   names(se.lai) = c(names(se.use),"sampling","census")
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Create population matrices.                                                      #
   #---------------------------------------------------------------------------------------#
   NPLANT = matrix( data = rep(datum$nplant,times=n.real), nrow=ndatum,ncol=n.real)
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #      Initialise answer.  We will create the vector before and loop through the        #
   # errors.                                                                               #
   #---------------------------------------------------------------------------------------#
   for (e in seq_along(se.use)){
      if (e == 1){
         LA     = matrix( data = rlnorm( n       = n.real*ndatum
                                       , meanlog = rep(log(datum$LA),times=n.real)
                                       , sdlog   = rep(se.use[[e]]  ,times=n.real)
                                       )#end rlnorm
                        , nrow = ndatum
                        , ncol = n.real
                        )#end matrix
         LA     = ifelse(test=is.finite(LA),yes=LA,no=0)
      }else{
         LA     = matrix( data = rnorm( n    = n.real*ndatum
                                      , mean = rep(datum$LA   ,times=n.real)
                                      , sd   = rep(se.use[[e]],times=n.real)
                                      )#end rlnorm
                        , nrow = ndatum
                        , ncol = n.real
                        )#end matrix
      }#end if (e == 1)

      LAI        = colSums(NPLANT * LA)
      se.lai[e]  = sd(LAI)
   }#end for (e in seq_along(se.use))
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Find total biomass, to be used by the within plot sampling error.                #
   #---------------------------------------------------------------------------------------#
   lai.bar = sum(datum$nplant * datum$LA)
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #      Unless epsilon.smp is provided, we estimate the within-plot sampling error.      #
   #---------------------------------------------------------------------------------------#
   if (is.null(epsilon.smp)){
      #------------------------------------------------------------------------------------#
      #      Split domain into smaller subsplots.                                          # 
      #------------------------------------------------------------------------------------#
      if (subalong %in% "x"){
         xbreaks = seq(from=0,to=xmax,length.out=n.sub+1)
         xlwr    =       sqrt(.Machine$double.eps)  * xmax
         xupr    = (1. - sqrt(.Machine$double.eps)) * xmax
         xnow    = pmax(xlwr,pmin(xupr,datum$X))
         isub    = as.integer(cut(xnow,breaks=xbreaks))
      }else{
         ybreaks = seq(from=0,to=ymax,length.out=n.sub+1)
         ylwr    =       sqrt(.Machine$double.eps)  * ymax
         yupr    = (1. - sqrt(.Machine$double.eps)) * ymax
         ynow    = pmax(ylwr,pmin(yupr,datum$Y))
         isub    = as.integer(cut(ynow,breaks=ybreaks))
      }#end if
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Find the AGB/ACD for each subplot.                                            #
      #------------------------------------------------------------------------------------#
      lai.sub          = rep(0.,times=n.sub)
      lai.tmp          = tapply(X=n.sub*datum$nplant*datum$LA ,INDEX=isub,FUN=sum)
      lai.idx          = as.integer(names(lai.tmp))
      lai.sub[lai.idx] = lai.tmp
      #------------------------------------------------------------------------------------#




      #------ Create replicates using bootstrap with replacement. -------------------------#
      LAI = colMeans( matrix( data = sample(lai.sub,size=n.real*n.sub,replace=TRUE)
                            , nrow = n.sub
                            , ncol = n.real
                            )#end matrix
                    )#end colMeans
      se.lai["sampling"] = sd(LAI)
      #------------------------------------------------------------------------------------#
   }else{
      #----- Use pre-defined sampling error. ----------------------------------------------#
      se.lai["sampling"] = epsilon.smp * lai.bar
      #------------------------------------------------------------------------------------#
   }#end if (is.null(epsilon.smp))
   #---------------------------------------------------------------------------------------#



   #------ Find combined source of errors. ------------------------------------------------#
   se.lai["census"] = sqrt(sum(se.lai[c("measurement","allometry","sampling")]^2))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Error is the standard deviation of all realisations.                             #
   #---------------------------------------------------------------------------------------#
   ans        = se.lai
   names(ans) = paste0("se.lai.",names(se.lai))
   #---------------------------------------------------------------------------------------#


   #----- Return answer. ------------------------------------------------------------------#
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end find.lai.error
#==========================================================================================#
#==========================================================================================#
