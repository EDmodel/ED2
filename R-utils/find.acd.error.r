#==========================================================================================#
#==========================================================================================#
#      This function determines the plot aggregated error due to measurement and           #
# allometry.                                                                               #
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
