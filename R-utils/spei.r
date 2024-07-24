vs10 <<- c( C0 = 2.515517, C1 = 0.802853, C2 = 0.010328
          , d0 = 1.000000, d1 = 1.432788, d2 = 0.189269, d3 = 0.001308
          )#end c

#==========================================================================================#
#==========================================================================================#
#     This function calculates the standardised precipitation evapotranspiration index     #
# (SPEI), based on VS10.  The original publication used the log-logistic distribution, but #
# this function allows the user to fit other distributions (normal, skew-normal, .                                                                      #
#                                                                                          #
# References:                                                                              #
#                                                                                          #
# Vicente-Serrano SM, Begueria S, Lopez-Moreno JI. 2010. A multiscalar drought index       #
#    sensitive to global warming: The standardized precipitation evapotranspiration index. #
#    J. Climate, 23: 1696-1718. doi:10.1175/2009JCLI2909.1 (VS10).                         #
#                                                                                          #
# Hosking JRM. 1990. L-moments: Analysis and estimation of distributions using linear      #
#    combinations of order statistics. J. R. Stat. Soc. B, 52: 105-124.                    #
#    doi:10.1111/j.2517-6161.1990.tb01775.x (H90).                                         #
#------------------------------------------------------------------------------------------#
spei <<- function( x
                 , toler  = 1.e-5
                 , itmax  = 50
                 , finite = TRUE
                 , distrib = c("llog3","sn")
                 ){

    #--------------------------------------------------------------------------------------#
    #     Select the distribution to use.                                                  #
    #--------------------------------------------------------------------------------------#
    distrib = match.arg(distrib)
    #--------------------------------------------------------------------------------------#



    #--------------------------------------------------------------------------------------#
    #    Decide what to do based on na.rm.                                                 #
    #--------------------------------------------------------------------------------------#
    if ( ((! finite) && any(! is.finite(x))) || (all(! is.finite(x))) ){
       #----- Invalid entries. Return nothing. --------------------------------------------#
       ans = NA_real_ * x
       #-----------------------------------------------------------------------------------#
    }else{
       #----- Exclude invalid entries. ----------------------------------------------------#
       fine  = is.finite(x)
       xfine = x[fine]
       nfine = length(xfine)
       cdf   = NA_real_ * x
       ans   = NA_real_ * x
       #-----------------------------------------------------------------------------------#


       #-----------------------------------------------------------------------------------#
       #     Find the first guesses for weights and parameter estimators.                  #
       #-----------------------------------------------------------------------------------#
       if (distrib %in% "llog3"){
          #----- Find first guess for parameters. -----------------------------------------#
          xlwr    = min(xfine) - 0.01 * diff(range(xfine))
          fit1st  = fitdistr( x       = x - xlwr
                            , densfun = dllogis
                            , start   = list(shape=1,scale=1)
                            )#end fit1st
          ln.scale = log(fit1st$estimate["scale"])
          ln.shape = log(fit1st$estimate["shape"])
          #--------------------------------------------------------------------------------#


          #----- Fit distribution. --------------------------------------------------------#
          fitnow   = try( fitdistr( x       = xfine
                                  , densfun = dologis
                                  , start   = list( thres    = xlwr
                                                  , ln.scale = ln.scale
                                                  , ln.shape = ln.shape
                                                  )#end list
                                )#end fitdistr
                      , silent = TRUE
                      )#end try
          if (! ("try-error" %in% is(fitnow))){
             #------ Find the statistics. -------------------------------------------------#
             ll.scale = exp(fitnow$estimate["ln.scale"])
             ll.shape = exp(fitnow$estimate["ln.shape"])
             ll.thres = fitnow$estimate["thres"]
             #-----------------------------------------------------------------------------#


             #----- Find the cumulative distribution function values. ---------------------#
             cdf[fine] = pllog3(q=xfine,thres=ll.thres,shape=ll.shape,scale=ll.scale)
             #-----------------------------------------------------------------------------#
          }#end if (! ("try-error" %in% is(fitnow)))
          #--------------------------------------------------------------------------------#
       }else if (distrib %in% "sn"){
          #----- Find first guess for parameters. -----------------------------------------#
          xmean   = mean(xfine)
          xsdev   = sd(xfine)
          #--------------------------------------------------------------------------------#


          #----- Fit distribution. --------------------------------------------------------#
          fitnow   = try( fitdistr( x       = xfine
                                  , densfun = dsn
                                  , start   = list( xi    = xmean
                                                  , omega = xsdev
                                                  , alpha = 0.
                                                  )#end list
                                  , tau     = 0.
                                  )#end fitdistr
                      , silent = TRUE
                      )#end try
          if (! ("try-error" %in% is(fitnow))){
             #----- Retrieve the statistics. ----------------------------------------------#
             sn.xi    = fitnow$estimate["xi"   ]
             sn.omega = fitnow$estimate["omega"]
             sn.alpha = fitnow$estimate["alpha"]
             #-----------------------------------------------------------------------------#


             #----- Find the cumulative distribution function values. ---------------------#
             cdf[fine] = psn(x=xfine,xi=sn.xi,omega=sn.omega,alpha=sn.alpha)
             #-----------------------------------------------------------------------------#
          }#end if ("try-error" %in% is(fitnow))
          #--------------------------------------------------------------------------------#
       }#end if (distrib %in% "llog3")
       #-----------------------------------------------------------------------------------#



       #-----------------------------------------------------------------------------------#
       #     Find the standardised values of p3cdf:                                        #
       #-----------------------------------------------------------------------------------#
       # P     = 1. - cdf
       # FS    = ifelse( test = P %le% 0.5, yes = 1.0, no = -1.0)
       # PUSE  = 0.5 - sqrt((P - 0.5)^2)
       # W     = sqrt(-2. * log(PUSE))
       #-----------------------------------------------------------------------------------#




       #-----------------------------------------------------------------------------------#
       #     Find the SPEI, following VS10:                                                #
       #-----------------------------------------------------------------------------------#
       # snom = vs10["C0"] + W * ( vs10["C1"] + W *   vs10["C2"] )
       #sden = vs10["d0"] + W * ( vs10["d1"] + W * ( vs10["d2"] + W * vs10["d3"] ) )
       # ans  = FS * ( W - snom / sden )
       #-----------------------------------------------------------------------------------#


       #-----------------------------------------------------------------------------------#
       #     To obtain the index, we convert the cumulative distribution function into     #
       # the equivalent quantile of a normal distribution of mean zero and standard        #
       # deviation one.                                                                    #
       #-----------------------------------------------------------------------------------#
       ans[fine] = qnorm(p=cdf[fine],mean=0,sd=1)
       #-----------------------------------------------------------------------------------#
    }#end if
    #--------------------------------------------------------------------------------------#


    return(ans)
}#end spei
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#      Wrapper for the three-parameter log-logistic distribution with log-transformed      #
# parameters (to help with optimisation).                                                  #
#------------------------------------------------------------------------------------------#
#------ Density. --------------------------------------------------------------------------#
dologis <<- function(x,ln.scale,ln.shape,thres,...){
   scale = exp(ln.scale)
   shape = exp(ln.shape)
   ans   = dllog3(x=pmax(x,thres),scale=scale,shape=shape,thres=thres,...)
   return(ans)
}#end dologis
#==========================================================================================#
#==========================================================================================#
