#==========================================================================================#
#==========================================================================================#
#     This is just a wrapper for boot.ci, in which we retrieve the lower bound of the      #
# confidence interval only.                                                                #
#------------------------------------------------------------------------------------------#
boot.expected <<- function(data,statistic,R,...){
    ans = boot(data,statistic,R)$t0
    return(ans)
}#end boot.ci.lower
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This is just a wrapper for boot.ci, in which we retrieve the lower bound of the      #
# confidence interval only.                                                                #
#------------------------------------------------------------------------------------------#
boot.ci.lower <<- function(boot.out,...){
    ci.now = boot.ci(boot.out,...)$percent
    if (length(ci.now) == 5){
       ans = ci.now[4]
    }else{
       ans = NA
       warning(" Failed using bootstrap...")
    }#end if
    return(ans)
}#end boot.ci.lower
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This is just a wrapper for boot.ci, in which we retrieve the upper bound of the      #
# confidence interval only.                                                                #
#------------------------------------------------------------------------------------------#
boot.ci.upper <<- function(boot.out,...){
    ci.now = boot.ci(boot.out,...)$percent
    if (length(ci.now) == 5){
       ans = ci.now[5]
    }else{
       ans = NA
       warning(" Failed using bootstrap...")
    }#end if
    return(ans)
}#end boot.ci.upper
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     Auxiliary functions that carry the index, useful for bootstrap.                      #
#------------------------------------------------------------------------------------------#
boot.mean   <<- function(x,idx) mean    ( x = x[idx]               , na.rm = TRUE)
boot.median <<- function(x,idx) median  ( x = x[idx]               , na.rm = TRUE)
boot.sd     <<- function(x,idx) sd      ( x = x[idx]               , na.rm = TRUE)
boot.var    <<- function(x,idx) var     ( x = x[idx]               , na.rm = TRUE)
boot.sum    <<- function(x,idx) sum     ( x = x[idx]               , na.rm = TRUE)
boot.q025   <<- function(x,idx) quantile( x = x[idx], probs = 0.025, na.rm = TRUE)
boot.q975   <<- function(x,idx) quantile( x = x[idx], probs = 0.975, na.rm = TRUE)
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     Demography-related functions.                                                        #
#------------------------------------------------------------------------------------------#
#----- Recruitment. -----------------------------------------------------------------------#
boot.recruit <<- function(dat,idx,dtime){
   n.idx         = length(idx)
   p.use         = rbinom(n=n.idx,size=1,prob=dat$p.use        [idx])
   p.established = rbinom(n=n.idx,size=1,prob=dat$p.established[idx])
   N             = sum( dat$property   [idx] * p.use                ,na.rm=TRUE)
   E             = sum( dat$property   [idx] * p.use * p.established,na.rm=TRUE)
   ans           = exp(log(N/E)/dtime) - 1.0

   return(ans)
}#end boot.recruit
#----- Mortality. -------------------------------------------------------------------------#
boot.mortality <<- function(dat,idx,dtime){
   n.idx         = length(idx)
   p.use         = rbinom(n=n.idx,size=1,prob=dat$p.use        [idx])
   p.survivor    = rbinom(n=n.idx,size=1,prob=dat$p.survivor   [idx])
   N             = sum( dat$property   [idx] * p.use                ,na.rm=TRUE)
   S             = sum( dat$property   [idx] * p.use * p.survivor   ,na.rm=TRUE)
   ans = 1.0 - exp(- log(N/S)/dtime)
   return(ans)
}#end boot.mortality
#----- Growth. ----------------------------------------------------------------------------#
boot.growth <<- function(dat,idx){
   n.idx         = length(idx)
   x             = dat$growth[idx]
   w             = dat$count [idx]
   ans           = -log(weighted.mean(x=exp(-x),w=w,na.rm=TRUE))
   return(ans)
}#end boot.growth
#----- Accumulated recruitment. -----------------------------------------------------------#
boot.acc.recruit <<- function(dat,idx,dtime){
   n.idx         = length(idx)
   p.use         = rbinom(n=n.idx,size=1,prob=dat$p.use        [idx])
   p.established = rbinom(n=n.idx,size=1,prob=dat$p.established[idx])

   p.recruit     = p.use * ( 1. - p.established)

   ans           = sum( dat$property   [idx] * p.recruit / dtime    ,na.rm=TRUE)
   return(ans)
}#end boot.recruit
#----- Accumulated mortality. -------------------------------------------------------------#
boot.acc.mortality <<- function(dat,idx,dtime){
   n.idx         = length(idx)
   p.use         = rbinom(n=n.idx,size=1,prob=dat$p.use        [idx])
   p.survivor    = rbinom(n=n.idx,size=1,prob=dat$p.survivor   [idx])

   p.dead        = p.use * (1. - p.survivor)
   ans           = sum( dat$property   [idx] * p.dead / dtime      ,na.rm=TRUE)
   return(ans)
}#end boot.mortality
#----- Above-ground net primary productivity. ---------------------------------------------#
boot.acc.growth <<- function(dat,idx){
   n.idx = length(idx)
   ans   = sum(dat$pop[idx] * ( dat$nok[idx] - dat$lok[idx] )  / dat$dtime[idx] )
   return(ans)
}#end boot.acc.growth
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     Binomial distribution function.  You can return either the expected value, or the    #
# lower or upper bounds of the confidence interval.                                        #
#------------------------------------------------------------------------------------------#
boot.binom <<- function(dat,idx,out="expected",conf=0.95){
   #---------------------------------------------------------------------------------------#
   #      Success is estimated using the probability.                                      #
   #---------------------------------------------------------------------------------------#
   N     = sum(dat$count[idx],na.rm=TRUE)
   S     = sum(rbinom(n=length(idx),size=1,prob=dat$probability[idx]),na.rm=TRUE)
   prob  = success / count
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Check which output to use.                                                       #
   #---------------------------------------------------------------------------------------#
   out.int = tolower(substring(out,1,1))
   if (out.int == "e"){
      #----- Expected value. --------------------------------------------------------------#
      ans = prob
      #------------------------------------------------------------------------------------#
   }else if (out.int %in% c("l","u","h")){
      #------------------------------------------------------------------------------------#
      #     Confidence interval.  Use the "exact" binomial confidence interval.            #
      #------------------------------------------------------------------------------------#
      pci95 = c( qbeta(p = (1. - conf)/2., shape1 = S  , shape2 = N-S+1 )
               , qbeta(p = (1. + conf)/2., shape1 = S+1, shape2 = N-S   ) )
      #------------------------------------------------------------------------------------#
      if(out.int %in% c("l")){
         #------ Lower bound. -------------------------------------------------------------#
         ans = min(pci95)
         #---------------------------------------------------------------------------------#
      }else{
         #------ Upper bound. -------------------------------------------------------------#
         ans = max(pci95)
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#
   }else{
      #----- Invalid option, quit. --------------------------------------------------------#
      cat(" - Requested output: ",out,"\n")
      cat(" - Accepted options: 'expected', 'lower', or 'upper' (or 'higher')","\n")
      stop(" Invalid output")
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end boot.binom
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      This function computes the fortnightly means using only days with full              #
# record, and using bootstrap to sample the days.                                          #
#------------------------------------------------------------------------------------------#
boot.fortnight.mean <<- function(data.in,index){

   data.use = data.in[index,]

   #----- Aggregate the data by hour and fortnightly period. ------------------------------#
   ta.fnmean   = tapply(X=data.use$x,INDEX=data.use$fortnight,FUN=mean,na.rm=TRUE)
   idx.fnmean  = as.numeric(names(ta.fnmean))
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Collapse fortnightly periods.  Make sure all periods are defined in the output,   #
   # if none of them were selected, make them NA.                                          #
   #---------------------------------------------------------------------------------------#
   fnmean                      = rep(NA,times=24)
   fnmean[idx.fnmean]          = ta.fnmean
   fnmean[! is.finite(fnmean)] = NA
   #---------------------------------------------------------------------------------------#

   return(fnmean)
}#end function
#==========================================================================================#
#==========================================================================================#

