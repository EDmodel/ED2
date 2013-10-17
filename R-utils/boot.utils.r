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
boot.mean   <<- function(x,idx) mean        ( x = x[idx]               , na.rm = TRUE)
boot.median <<- function(x,idx) median      ( x = x[idx]               , na.rm = TRUE)
boot.sd     <<- function(x,idx) sd          ( x = x[idx]               , na.rm = TRUE)
boot.var    <<- function(x,idx) var         ( x = x[idx]               , na.rm = TRUE)
boot.sum    <<- function(x,idx) sum         ( x = x[idx]               , na.rm = TRUE)
boot.q025   <<- function(x,idx) quantile    ( x = x[idx], probs = 0.025, na.rm = TRUE)
boot.q975   <<- function(x,idx) quantile    ( x = x[idx], probs = 0.975, na.rm = TRUE)
boot.moment <<- function(x,idx) four.moments( x = x[idx]               , na.rm = TRUE)
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
#      This function computes the fortnightly means sampling hours, but not the years and  #
# the fortnights.                                                                          #
#------------------------------------------------------------------------------------------#
boot.fortnight.mean <<- function(data.in,R,ci=0.95,...){
   call.now = match.call()


   #----- Save the number of hours. -------------------------------------------------------#
   lab.hour    = sort(unique(data.in$hour     ))
   lab.ftnight = sort(unique(data.in$fortnight))
   lab.year    = sort(unique(data.in$year     ))
   n.hour      = length(lab.hour)
   n.year      = length(lab.year)
   n.ftnight   = yr.ftnight
   #---------------------------------------------------------------------------------------#


   #------ Append bootstrap mean to the list of arguments to go to bootstrap. -------------#
   dotdotdot = modifyList( x   = list(...)
                         , val = list(R=R,stat=mean,realisation.only=TRUE,na.rm=TRUE)
                         )#end modifyList
   #---------------------------------------------------------------------------------------#


   #----- Split the data into lists. ------------------------------------------------------#
   list.use = split(x = data.in$x, f = list(data.in$hour,data.in$fortnight,data.in$year))
   n.list   = lapply(X=list.use,FUN=length)
   if (any(unlist(n.list) == 0)){
      stop("Empty elements in your list!")
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Run bootstrap for each class.                                                    #
   #---------------------------------------------------------------------------------------#
   boot.samples  = mapply( FUN = boot.sampling, x = list.use, MoreArgs  = dotdotdot)
   boot.arr4     = array(data=boot.samples,dim=c(R,n.hour,n.ftnight,n.year))
   boot.arr3     = apply(X=boot.arr4,MARGIN=c(1,3,4),FUN=mean)
   boot.mat      = apply(X=boot.arr3,MARGIN=c(1,2),FUN=mean,na.rm=TRUE)
   boot.expected = apply(X=boot.mat,MARGIN=2,FUN=mean,na.rm=TRUE)
   boot.qlow     = apply(X=boot.mat,MARGIN=2,FUN=quantile,prob=0.5*(1.-ci),na.rm=TRUE)
   boot.qhigh    = apply(X=boot.mat,MARGIN=2,FUN=quantile,prob=0.5*(1.+ci),na.rm=TRUE)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Collapse fortnightly periods.  Make sure all periods are defined in the output,   #
   # if none of them were selected, make them NA.                                          #
   #---------------------------------------------------------------------------------------#
   empty    = rep(x = NA, times = yr.ftnight)
   expected = empty
   qlow     = empty
   qhigh    = empty
   expected[lab.ftnight] = ifelse(is.finite(boot.expected),boot.expected,NA)
   qlow    [lab.ftnight] = ifelse(is.finite(boot.qlow    ),boot.qlow    ,NA)
   qhigh   [lab.ftnight] = ifelse(is.finite(boot.qhigh   ),boot.qhigh   ,NA)
   ans      = list(call=call.now,expected=expected,qlow=qlow,qhigh=qhigh)
   #---------------------------------------------------------------------------------------#


   #------ Return the statistics. ---------------------------------------------------------#
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function
#==========================================================================================#
#==========================================================================================#







#==========================================================================================#
#==========================================================================================#
#      This function computes the fortnightly means using hierarchical bootstrap on years  #
# and hours.                                                                               #
#------------------------------------------------------------------------------------------#
bhier.fortnight.mean <<- function(data.in,R,ci=0.95,...){
   call.now = match.call()


   #----- Save the number of hours. -------------------------------------------------------#
   lab.hour    = sort(unique(data.in$hour     ))
   lab.ftnight = sort(unique(data.in$fortnight))
   lab.year    = sort(unique(data.in$year     ))
   n.hour      = length(lab.hour)
   n.year      = length(lab.year)
   n.ftnight   = yr.ftnight
   n.data.in   = nrow(data.in)
   #---------------------------------------------------------------------------------------#


   #------ Append bootstrap mean to the list of arguments to go to bootstrap. -------------#
   dotdotdot = modifyList( x   = list(...)
                         , val = list(R=R,stat=mean,realisation.only=TRUE,na.rm=TRUE)
                         )#end modifyList
   #---------------------------------------------------------------------------------------#


   #----- Split the data into lists just to test that there aren't empty elements. --------#
   list.use = split(x = data.in$x, f = list(data.in$hour,data.in$fortnight,data.in$year))
   n.list   = lapply(X=list.use,FUN=length)
   if (any(unlist(n.list) == 0)){
      stop("Empty elements in your list!")
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Split data into three levels: fortnights, year, and hour.                        #
   #---------------------------------------------------------------------------------------#
   level.1.list  = split(x = data.in     , f = data.in$fortnight)
   level.2.list  = mapply( FUN      = function(dat) split(x=dat,f=dat$year)
                         , dat      = level.1.list
                         , SIMPLIFY = FALSE
                         )#End mapply
   level.3.list  = mapply( FUN      = function(dat){
                                         mapply( FUN = function(dat) split(x=dat,f=dat$hour)
                                               , dat = dat
                                               , SIMPLIFY = FALSE
                                               )#end mapply
                                      }#end function
                         , dat      = level.2.list
                         , SIMPLIFY = FALSE
                         )#End mapply
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Run bootstrap for each hour.                                                     #
   #---------------------------------------------------------------------------------------#
   level.2.sample = mapply( FUN = function(dat,...){
                                     mapply( FUN      = boot.list
                                           , dat      = dat
                                           , MoreArgs = list(...)
                                           , SIMPLIFY = FALSE
                                           )#end mapply
                                  }#end function(dat)
                          , dat = level.3.list
                          , MoreArgs = dotdotdot
                          , SIMPLIFY = FALSE
                          )#end level.3.sample
   #---------------------------------------------------------------------------------------#






   #---------------------------------------------------------------------------------------#
   #       Sample the years.                                                               #
   #---------------------------------------------------------------------------------------#
   level.1.sample = mapply( FUN      = high.sampler
                          , x        = level.2.sample
                          , MoreArgs = list(R=R)
                          , SIMPLIFY = FALSE
                          )#end mapply
   boot.samples   = mapply( FUN      = boot.collapse
                          , x        = level.1.sample
                          , SIMPLIFY = FALSE
                          )#end mapply
   #---------------------------------------------------------------------------------------#




   #------ Find expected value and confidence intervals. ----------------------------------#
   boot.expected = sapply(X=boot.samples,FUN=mean,na.rm=TRUE)
   boot.se       = sapply(X=boot.samples,FUN=sd  ,na.rm=TRUE)
   boot.qlow     = sapply(X=boot.samples,FUN=quantile,prob=0.5*(1.-ci),na.rm=TRUE)
   boot.qhigh    = sapply(X=boot.samples,FUN=quantile,prob=0.5*(1.+ci),na.rm=TRUE)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Collapse fortnightly periods.  Make sure all periods are defined in the output,   #
   # if none of them were selected, make them NA.                                          #
   #---------------------------------------------------------------------------------------#
   empty    = rep(x = NA, times = yr.ftnight)
   expected = empty
   std.err  = empty
   qlow     = empty
   qhigh    = empty
   expected[lab.ftnight] = ifelse(is.finite(boot.expected),boot.expected,NA)
   std.err [lab.ftnight] = ifelse(is.finite(boot.se      ),boot.se      ,NA)
   qlow    [lab.ftnight] = ifelse(is.finite(boot.qlow    ),boot.qlow    ,NA)
   qhigh   [lab.ftnight] = ifelse(is.finite(boot.qhigh   ),boot.qhigh   ,NA)
   ans      = list(call=call.now,expected=expected,se=std.err,qlow=qlow,qhigh=qhigh)
   #---------------------------------------------------------------------------------------#


   #------ Return the statistics. ---------------------------------------------------------#
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#    Boot.sampling is a simple sampler that doesn't require package boot.                  #
#------------------------------------------------------------------------------------------#
boot.sampling <<- function( x
                          , stat
                          , R                = 1000
                          , ci.level         = 0.95
                          , realisation.only = FALSE
                          , quiet            = TRUE
                          , not.finite.2.na  = TRUE
                          ,...
                          ){
   rcmax = 500000000


   #----- Save call for return (and for error evaluation). --------------------------------#
   mycall = match.call()
   #---------------------------------------------------------------------------------------#


   #---- Coerce x into a vector, and find its size. ---------------------------------------#
   x  = unlist(x)
   nx = length(x)
   #---------------------------------------------------------------------------------------#



   #---- Coerce x into a vector, and find its size. ---------------------------------------#
   if (nx <= 2){
      if (! quiet) warning(paste(" Vector x is too short (length=",nx,")",sep=""))
      realisation = rep(NA,times=R)


      #----- Find the statistics. ---------------------------------------------------------#
      if (! realisation.only){
         expected    = NA
         std.err     = NA
         tci         = c(NA,NA)
      }#end if
      #------------------------------------------------------------------------------------#
   }else{
      #------------------------------------------------------------------------------------#
      #     Check size.  If the product of length of x and number of iterations is not too #
      # long, use apply, otherwise, for loop is actually faster.                           #
      #------------------------------------------------------------------------------------#
      use.apply = nx*R < rcmax
      if (use.apply){
         idx.real = rep(x=sequence(R),each=nx)
         x.sample = sample(x=x,size=nx*R,replace=TRUE)
         realisation = tapply(X=x.sample,INDEX=idx.real,FUN=stat,...)
      }else{
         realisation = rep(NA,times=R)
         for (i in sequence(nx)){
            idx            = sample(x,size=nx,replace=TRUE)
            realisation[i] = stat(x[idx],...)
         }#end for
      }#end if
      if (not.finite.2.na) realisation[! is.finite(realisation)] = NA
      #------------------------------------------------------------------------------------#


      #----- Find the statistics. ---------------------------------------------------------#
      if (! realisation.only){
         expected = mean(realisation)
         std.err  = se(realisation)
         tci      = quantile(realisation,prob=(1+c(-1,1)*ci.level)/2)
      }#end if
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Build the list with the results. ------------------------------------------------#
   if (realisation.only){
      ans = realisation
   }else{
      ans = list( call        = mycall
                , realisation = realisation
                , expected    = expected
                , std.err     = std.err
                , ci          = tci
                , R           = R
                , m           = m
                , method      = method
                )#end list
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Return answer. ------------------------------------------------------------------#
   return(ans)
   #---------------------------------------------------------------------------------------#

}#end function boot.sampling
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      Additional auxiliary functions for hierarchical sampling.  These functions are not  #
# intended to be directly called, but if you find some use for it, feel free to use them.  #
#------------------------------------------------------------------------------------------#
#----- Bootstrap the statistics within the inner level of a list. -------------------------#
boot.list <<- function(dat,...){
   ans = mapply( FUN      = function(dat,...) boot.sampling(x=dat$x,...)
               , dat      = dat
               , MoreArgs = list(...)
               , SIMPLIFY = TRUE
               )#end mapply
   return(ans)
}#end function boot.list
#----- Select realisations for higher hierarchical level for when it has few points. ------#
high.sampler <<- function(x,R){
   ans = replicate( n    = R
                  , expr = replicate( length(x)
                                    , x[[sample(length(x),size=1)]][sample(R,size=1),]
                                    )#end replicate
                  )#end replicate
   return(ans)
}#end function high.sampler
#----- Collapse the realisations by hour and by year. -------------------------------------#
boot.collapse <<- function(x){
    tmp                   = apply(X=x  ,MARGIN=c(2,3),FUN=mean,na.rm=FALSE)
    tmp[! is.finite(tmp)] = NA
    ans                   = apply(X=tmp,MARGIN=2     ,FUN=mean,na.rm=TRUE )
    ans[! is.finite(ans)] = NA
    return(ans)
}#end function boot.collapse
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      This is a wrapper that computes multiple statistics using bootstrapping.            #
#------------------------------------------------------------------------------------------#
boot.six.summary <<- function(x,R,conf=0.95){

   #------ Initialise data with NA in case the function fails. ----------------------------#
   ans = rep(NA,times=n.six.summary)
   names(ans) = six.summary.names
   #---------------------------------------------------------------------------------------#


   #------ Delete non-finite values. ------------------------------------------------------#
   if (length(x) != 0) x  = x[is.finite(x)]
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Run bootstrap only if there are data.                                            #
   #---------------------------------------------------------------------------------------#
   if (length(x) > 0){
      #------ Run bootstrap and get the estimate of the four moments. ---------------------#
      bo       = boot(data=x,statistic=boot.moment,R=R)
      expected = bo$t0[1]
      variance = bo$t0[2]
      skewness = bo$t0[3]
      kurtosis = bo$t0[4]
      #------------------------------------------------------------------------------------#



      #----- Find confidence intervals (silent R because boot.ci has annoying messages). --#
      bci = shhh(fun=boot.ci,boot.out=bo,index=c(1,2),conf=conf,type=c("stud","perc"))
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Check whether confidence interval worked.                                      #
      #------------------------------------------------------------------------------------#
      if ("try-error" %in% is(bci)){
         ci.lower = quantile(x=bo$t[,1],prob=0.5*(1.0-conf),na.rm=TRUE)
         ci.upper = quantile(x=bo$t[,1],prob=0.5*(1.0+conf),na.rm=TRUE)
      }else if(length(bci$student) == 5 && all(is.finite(bci$student))){
         ci.lower = bci$student[4]
         ci.upper = bci$student[5]
      }else if(length(bci$percent) == 5 && all(is.finite(bci$percent))){
         ci.lower = bci$percent[4]
         ci.upper = bci$percent[5]
      }else{
         ci.lower = quantile(x=bo$t[,1],prob=0.5*(1.0-conf),na.rm=TRUE)
         ci.upper = quantile(x=bo$t[,1],prob=0.5*(1.0+conf),na.rm=TRUE)
      }#end if
      #------------------------------------------------------------------------------------#




      #----- Standardise non-finite values to NA. -----------------------------------------#
      expected = ifelse(is.finite(expected),expected,NA)
      variance = ifelse(is.finite(variance),variance,NA)
      skewness = ifelse(is.finite(skewness),skewness,NA)
      kurtosis = ifelse(is.finite(kurtosis),kurtosis,NA)
      ci.lower = ifelse(is.finite(ci.lower),ci.lower,NA)
      ci.upper = ifelse(is.finite(ci.upper),ci.upper,NA)
      #------------------------------------------------------------------------------------#



      #----- Make sure statistics go to the right place. ----------------------------------#
      ans[match("expected",six.summary.names)] = expected
      ans[match("variance",six.summary.names)] = variance
      ans[match("skewness",six.summary.names)] = skewness
      ans[match("kurtosis",six.summary.names)] = kurtosis
      ans[match("ci.lower",six.summary.names)] = ci.lower
      ans[match("ci.upper",six.summary.names)] = ci.upper
      #------------------------------------------------------------------------------------#
   }#end if (length(x) == 0)
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end function boot.six.summary
#==========================================================================================#
#==========================================================================================#
