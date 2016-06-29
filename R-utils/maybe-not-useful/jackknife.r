#==========================================================================================#
#==========================================================================================#
#      Function jackknife                                                                  #
#      This function estimates the confidence interval for some statistic using the        #
# jackknifing resampling.  Inputs:                                                         #
#                                                                                          #
# * x        -- the data set.                                                              #
# * stat     -- the statistic.  Any function can be used, though only the statistic will   #
#               be applied to the first argument only.                                     #
# * method   -- which method to use.  Options are:                                         #
#               default:  this is the true jackknifing, or the "leave-m out" method. If    #
#                         m > 1, then data will be removed in blocks.                      #
#               randset:  a "bootstrappy jackknife", randomly sub-sample leaving m out.    #
#                         this is only useful if m is a significant fraction of the        #
#                         original length, otherwise the default is better.                #
#               combn:    choose m to create the closest number of replications as R.      #
#                         m will be ignored in this case.                                  #
# * m        -- Number of points to be removed.  Either an integer with the actual number  #
#               of points, or if 0 < m < 1, the fraction to be removed each time.          #
# * R        -- Number of realisations (ignored if method is "default")                    #
# * ci.level -- Confidence interval                                                        #
# * ...      -- additional arguments for function stat.                                    #
#------------------------------------------------------------------------------------------#
jackknife <<- function(x,stat,method="default",m=1,R=1000,ci.level=0.95,...){
   rcmax = 1000000


   #----- Save call for return (and for error evaluation). --------------------------------#
   mycall = match.call()
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Check data size.  Jackknife requires at least 3 points (not that it will be any   #
   # good with 2 points...).  If not, return NA with warning.
   #---------------------------------------------------------------------------------------#
   nx = length(x)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     check the method.                                                                 #
   #---------------------------------------------------------------------------------------#
   meth = substring(tolower(method),1,1)
   if (! meth %in% c("d","r","c")){
      stop(paste(" Method: ",method," is invalid.! Enter either default, randset, or combn!"
                ,sep=""))
   }else if (meth %in% "d"){
      R      = nx
      method = "default"
   }else if (meth %in% "c"){
      cc    = choose(n=nx,k=seq(from=1,to=floor(nx/2)))
      m     = which.min(ifelse(cc<R,Inf,cc-R))
      R.try = cc[m]
      if (R.try > rcmax){
         warning(paste(" Too many combinations (R=",R.try,"). Use random sampling..."
                      ,sep=""))
         meth   = "r"
         method = "randset"
      }else{
         R      = R.try
         method = "combn"
      }#end if
   }else{
      combmax = choose(nx,m)
      if (R > combmax && combmax > rcmax){
         stop(paste(" Too many realisations and they wouldn't be independent!","\n"
                   ," R =",R," combinations = ",combmax," rcmax = ",rcmax,sep=""))
      }else if (R > combmax){
         warning(paste("   Too many realisations for too few independent combinations!"
                      ,"\n   Using combination method...",sep=""))
         R      = combmax
         method = "combn"
         meth   = "c"
      }else{
         method = "randset"
      }#end if
   }#end if
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #      Check sub-sample size.                                                           #
   #---------------------------------------------------------------------------------------#
   if (m <= 0){
      stop(paste(" m =",m," is not valid! Choose a positive number!",sep=""))
   }else if (m < 1){
      m = max(1,min(nx-2,round(nx*m)))
   }else if (m >= nx-2){
      warning(paste(" m =",m," is too large.  Using ",nx-2," instead...",sep=""))
      m = nx-2
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Don't try if nx is too small.                                                     #
   #---------------------------------------------------------------------------------------#
   if (nx <= 2){
      warning(paste(" Vector x is too short (length=",nx,")",sep=""))
      realisation = rep(NA,times=R)
      expected    = NA
      std.err     = NA
      ci          = c(NA,NA)
   }else{
      #------------------------------------------------------------------------------------#
      #     Check size.  If the product of length of x and number of iterations is not too #
      # long, use apply, otherwise, for loop is actually faster.                           #
      #------------------------------------------------------------------------------------#
      S = nx - m
      use.apply = nx*(nx+1) < rcmax
      if (meth %in% "d"){
         if (use.apply){
            realisation = apply( X      = matrix(x,nrow=nx+1,ncol=R)[sequence(S),]
                               , MARGIN = 2
                               , FUN    = stat
                               , ...
                               )#end apply
         }else{
            realisation = rep(NA,times=R)
            out = seq(from=nx-m,to=nx-1)
            for (i in sequence(R)){
               out = (out%%nx)+1
               realisation[i] = stat(x[-out],...)
            }#end for
         }#end if
      }else if (meth %in% "c"){
         #---------------------------------------------------------------------------------#
         #     Run all possible combinations.                                              #
         #---------------------------------------------------------------------------------#
         realisation = combn(x,m=m,FUN=stat,...)
         #---------------------------------------------------------------------------------#
      }else{
         #---------------------------------------------------------------------------------#
         #     Random sampling.  Check whether apply family is usable.                     #
         #---------------------------------------------------------------------------------#
         if (use.apply){
            realisation = apply( X      = replicate( n    = R
                                                   , expr = sample(x,size=S,replace=FALSE)
                                                   )#end replicate
                               , MARGIN = 2
                               , FUN    = stat
                               , ...
                               )#end apply
         }else{
            realisation = rep(NA,times=R)
            out = seq(from=nx-m,to=nx-1)
            for (i in sequence(nx)){
               out            = sample(nx,size=m,replace=FALSE)
               realisation[i] = stat(x[-out],...)
            }#end for
         }#end if
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#


      #----- Find the statistics. ---------------------------------------------------------#
      expected = mean(realisation)
      std.err  = se(realisation)
      tci      = quantile(realisation,prob=(1+c(-1,1)*ci.level)/2)
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Build the list with the results. ------------------------------------------------#
   ans = list( call        = mycall
             , realisation = realisation
             , expected    = expected
             , std.err     = std.err
             , ci          = tci
             , R           = R
             , m           = m
             , method      = method
             )#end list
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end jackknife
#==========================================================================================#
#==========================================================================================#
