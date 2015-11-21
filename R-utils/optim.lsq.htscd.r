#==========================================================================================#
#==========================================================================================#
#     This function finds the best heteroscedastic model.  This method simultaneously fits #
# both the main formula (lsq.formula) and the error following the formula given by         #
# sig.formula.                                                                             #
#                                                                                          #
# INPUT                                                                                    #
#                                                                                          #
# - lsq.formula:   the model to be fitted.                                                 #
#                  e.g.  y ~ exp(a0*xvar0+a1*xvar1)                                        #
# - sig.formula:   the normalised model for the standard deviation of the residuals, i.e.  #
#                  sigma/sigma0, where sigma0 is the standard deviation of all residuals.  #
#                  e.g. ~ exp(s2*xvar2)                                                    #
#                                                                                          #
#                  Variables in sig.formula must also be passed through data.              #
#                  In addition,  you may use the following variables:                      #
#                     yhat - fitted values                                                 #
#                     yres - residuals                                                     #
#                  If you have any variable with these names in data, then the model will  #
#                  use the variables stored in data, not the fitted or residuals.          #
#                  If you want to run a homoscedastic fit, set sig.formula = NULL (which   #
#                  is the default)
#                                                                                          #
# - data:          Data frame with predictors for both lsq.formula and sig.formula         #
# - lsq.first:     First guess for all parameters of lsq.formula to be fitted              #
# - sig.first:     First guess for all parameters of sig.formula to be fitted              #
#                  Required unless sig.formula = NULL                                      #
# - err.method:    Which method to use to estimate error.  Acceptable values are:          #
#                  - "hessian": use the eigenvalues from the Hessian matrix                #
#                  - "bootstrap": use bootstrap                                            #
#                  Partial matching is fine.                                               #
# - maxit.optim    Maximum number of iteration within each optim call.                     #
# - tol.gain:      Tolerance gain for multiple calls of optim (full model).                #
# - tol.optim:     Aimed tolerance for optimisation (full model).                          #
# - maxit:         Maximum number of optim calls (full model).                             #
# - n.boot:        Number of bootstrap realisations.                                       #
# - ci.level:      Confidence interval for error estimate (bootstrap only)                 #
# - verbose:       Show information whilst it is optimising.                               #
#                                                                                          #
#  OUTPUT: an object of type lsq.htscd (a list really) with the following elements.        #
#                                                                                          #
# - call:          Call that generated the current object                                  #
# - err.method:    Method for error evaluation                                             #
# - lsq.formula:   Model to fit through least squares                                      #
# - sig.formula:   Formula to represent the error                                          #
# - df:            Degrees of freedom                                                      #
# - coefficients:  Coefficients and summary table with error estimate                      #
# - x.best:        Optimised values (same as first column of coefficients)                 #
# - hessian:       Hessian matrix                                                          #
# - coeff.boot:    Coefficients for each realisation (bootstrap only)                      #
# - support.boot:  Support function for each realisation (bootstrap only)                  #
# - support:       Support for optimised coefficients.                                     #
# - data:          Data frame with model predictors used for fitting                       #
# - actual.values: Vector with actual values used for model fitting                        #
# - fitted.values: Vector with fitted values                                               #
# - residuals:     Vector with residuals                                                   #
# - sigma:         Local scale for residuals                                               #
# - goodness:      Summary statistics for goodness of fit                                  #
# - AIC:           Akaike information criterion                                            #
# - BIC:           Bayesian information criterion                                          #
# - conf.int:      Confidence interval used for cross validation (bootstrap only)          #
# - cross.val:     Data frame with cross validation data.  Estimates were obtained for     #
#                  points that were not included in each bootstrap realisation             #
#                  (bootstrap only)                                                        #
#   * n:           Number of independent cross validation estimates.                       #
#   * expected:    Expected value.                                                         #
#   * qlow:        Lower quantile.                                                         #
#   * qhigh:       Upper quantile.                                                         #
#   * bias:        Mean bias.                                                              #
#   * sigma:       Standard error of the residuals.                                        #
#   * rmse:        root mean square error.                                                 #
#------------------------------------------------------------------------------------------#
optim.lsq.htscd <<- function( lsq.formula
                            , sig.formula = NULL
                            , data
                            , lsq.first
                            , sig.first
                            , err.method  = c("hessian","hess","bootstrap","boot")
                            , maxit.optim = 20000
                            , tol.gain    = 0.0001
                            , tol.optim   = 0.00001
                            , is.debug    = FALSE
                            , maxit       = 100
                            , n.boot      = 1000
                            , ci.level    = 0.95
                            , i1st.max    = 5
                            , verbose     = FALSE
                            , ...
                            ){


   #----- Choose the preferred method to obtain coefficients and errors. ------------------#
   err.method = match.arg(err.method)
   err.short  = tolower(substring(err.method,1,1))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #   Make sure mandatory variables aren't missing.                                       #
   #---------------------------------------------------------------------------------------#
   if (  missing(lsq.formula) || missing(data)       || missing(lsq.first  )
      || ( missing(sig.first) && (! is.null(sig.formula)) ) ){
      cat("---------------------------------------------------------------","\n",sep="")
      cat("   Some variables are missing from function call:"              ,"\n",sep="")
      cat(" - lsq.formula is missing:   ",missing(lsq.formula)             ,"\n",sep="")
      cat(" - data is missing:          ",missing(data       )             ,"\n",sep="")
      cat(" - lsq.first is missing:     ",missing(lsq.first  )             ,"\n",sep="")
      if (! is.null(sig.formula)){
         cat(" - sig.formula is NULL:      ",is.null(sig.formula)          ,"\n",sep="")
         cat(" - sig.first is missing:     ",missing(sig.first  )          ,"\n",sep="")
      }#end if (! is.null(sig.formula))
      cat("---------------------------------------------------------------","\n",sep="")
      stop(" Missing variables")
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Set variables for homoscedastic fit in case sig.formula is NULL.                 #
   #---------------------------------------------------------------------------------------#
   if (is.null(sig.formula)){
      sig.first  = NULL
      skip.sigma = TRUE
   }else{
      skip.sigma = FALSE
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Convert first guesses to lists, they will be eventually reverted to vectors.      #
   #---------------------------------------------------------------------------------------#
   lsq.first = as.list(lsq.first)
   sig.first = as.list(sig.first)
   #---------------------------------------------------------------------------------------#


   #----- Make sure lsq.formula is formula, and extract variable names from formula. ------#
   lsq.formula    = try(as.formula(lsq.formula),silent=TRUE)
   if ("try-error" %in% is(lsq.formula)){
      lsq.formula = as.formula(eval(lsq.formula))
   }#end if
   lsq.vars       = all.vars(lsq.formula)
   yname          = lsq.vars[1]
   names.lsq.1st  = names(lsq.first)
   n.lsq.1st      = length(lsq.first)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      In case we want to fit a heteroscedastic fit, make sure sig.formula is formula,  #
   # and extract variable names from formula.  Also, make sure sigma0 is one of the        #
   # variables.                                                                            #
   #---------------------------------------------------------------------------------------#
   if (! skip.sigma){
      sig.formula = try(as.formula(sig.formula),silent=TRUE)
      if ("try-error" %in% is(sig.formula)){
         sig.formula = as.formula(eval(sig.formula))
      }#end if
      #----- Append sigma0 to sig.formula. ------------------------------------------------#
      if (! "sigma0" %in% all.vars(sig.formula)){
         #----- Append sigma0 to formula and to first guess. ------------------------------#
         sig.formula = attr(x=terms(sig.formula),which="term.labels")
         if (grepl(pattern="^I\\(",x=sig.formula) && grepl(pattern="\\)$",x=sig.formula)){
            sig.formula = gsub(pattern="^I\\(",replacement="",x=sig.formula)
            sig.formula = gsub(pattern="\\)$" ,replacement="",x=sig.formula)
            sig.formula = as.formula(paste0("~ I(sigma0*",sig.formula,")"))
         }else{
            sig.formula = as.formula(paste0("~ I(sigma0*",sig.formula,")"))
         }#end if
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Evaluate the first guess then use residuals to estimate sigma0.             #
         #---------------------------------------------------------------------------------#
         zero      = optim.lsq.htscd( lsq.formula = lsq.formula
                                    , sig.formula = NULL
                                    , data        = data
                                    , lsq.first   = lsq.first
                                    , err.method  = "hess"
                                    )#end optim.lsq.htscd
         sig.first = c(list(sigma0=mean(zero$sigma)),sig.first)
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#

      #----- Get variable names. ----------------------------------------------------------#
      sig.vars      = all.vars(sig.formula)
      names.sig.1st = names(sig.first)
      n.sig.1st     = length(sig.first)
      #------------------------------------------------------------------------------------#
   }else{
      sig.vars      = NULL
      names.sig.1st = NULL
      n.sig.1st     = 0
   }#end if
   #---------------------------------------------------------------------------------------#


   

   #---------------------------------------------------------------------------------------#
   #      Sanity check: lsq.first and sig.first must be named lists or named vectors, and  #
   # all names must appear in their respective formulae.                                   #
   #---------------------------------------------------------------------------------------#
   if (is.null(names.lsq.1st) || ( (! skip.sigma) && is.null(names.sig.1st))){
      cat("------------------------------------------------------------","\n",sep="")
      cat(" Names missing from lsq.first: ",is.null(names.lsq.1st)      ,"\n",sep="")
      if (skip.sigma){
         cat("------------------------------------------------------------","\n",sep="")
         stop("Both lsq.first must be a named list or a named vector!")
      }else{
         cat(" Names missing from sig.first: ",is.null(names.sig.1st)      ,"\n",sep="")
         cat("------------------------------------------------------------","\n",sep="")
         stop("Both lsq.first and sig.first must be named lists or named vectors!")
      }#end if (skip.sigma)
   }else{
      #------------------------------------------------------------------------------------#
      #      Check that all names for first guess exist in the formula.                    #
      #------------------------------------------------------------------------------------#
      fine.lsq.1st  = names.lsq.1st %in% lsq.vars[-1]
      if (skip.sigma){
         fine.sig.1st  = NULL
      }else{
         fine.sig.1st  = names.sig.1st %in% sig.vars
      }#end if (skip.sigma)
      if (! all(c(fine.lsq.1st,fine.sig.1st))){
         #---- Print a message telling that first guesses are not properly set. -----------#
         cat("------------------------------------------------------------","\n",sep="")
         cat("    Some names in the 1st guess are missing from formulae"   ,"\n",sep="")
         cat("    - lsq.first: ",names.lsq.1st[! fine.lsq.1st]             ,"\n",sep="")
         if(! skip.sigma){
            cat("    - sig.first: ",names.sig.1st[! fine.sig.1st]          ,"\n",sep="")
         }#end if(! skip.sigma)
         cat("------------------------------------------------------------","\n",sep="")
         stop("All names in 1st guesses must appear in the respective formulae!")
      }#end if (! all(c(fine.lsq.1st,fine.sig.1st)))
      #------------------------------------------------------------------------------------#
   }#end if (is.null(names.lsq.1st) || is.null(names.sig.1st))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Keep only the variables in data that are needed for the fitting.                  #
   #---------------------------------------------------------------------------------------#
   names.data    = names(data)
   names.data    = names.data[names.data %in% c(lsq.vars,sig.vars)]
   opt.data      = data[,names.data,drop=FALSE]
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Coerce first guesses to a vector, and append "lsq." and "sig." to their names.     #
   #---------------------------------------------------------------------------------------#
   if (skip.sigma){
      x.1st        = unlist(lsq.first)
      names(x.1st) = paste0("lsq.",names.lsq.1st)
   }else{
      x.1st        = c(unlist(lsq.first),unlist(sig.first))
      names(x.1st) = c(paste0("lsq.",names.lsq.1st),paste0("sig.",names.sig.1st))
   }#end if (skip.sigma)
   n.par        = length(x.1st)
   names.par    = c(names.lsq.1st,names.sig.1st)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Data selection and total number of parameters.                                    #
   #---------------------------------------------------------------------------------------#
   use      = apply(X=opt.data,MARGIN=1,FUN=function(x) all(is.finite(x)))
   opt.data = opt.data[use,,drop=FALSE]
   n.use    = nrow(opt.data)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Check that there are at least four valid data points, otherwise, crash!           #
   #---------------------------------------------------------------------------------------#
   if (n.use <= n.par){
      cat (" - Number of valid points: ",n.use,"\n")
      cat (" - Minimum number of valid points: ",n.par+1,"\n")
      stop(" Too few valid data points!")
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Check how verbose to be.                                                          #
   #---------------------------------------------------------------------------------------#
   if (is.logical(verbose)){
      optim.verbose = FALSE
   }else{
      optim.verbose = verbose >= 2
      verbose       = verbose >= 1
   }#end if
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Set options to the optimiser.                                                     #
   #---------------------------------------------------------------------------------------#
   ctrl.optim = list( trace   = optim.verbose
                    , REPORT  = 1
                    , fnscale = -1
                    , maxit   = maxit.optim
                    , reltol  = tol.optim
                    , ndeps   = rep(sqrt(tol.optim),times=n.par)
                    )#end list
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #    Update the first guess with Hessian, so it is not too far from the answer, then    #
   # determine the scale for each coefficient.                                             #
   #---------------------------------------------------------------------------------------#
   success = FALSE
   it      = 0
   it.giveup = ifelse(test=skip.sigma,yes=1,no=i1st.max)
   #----- Loop until a stable first guess is found. ---------------------------------------#
   while (! success & (it < i1st.max)){
      it          = it + 1
      opt.1st     = try( optim( par         = x.1st
                              , fn          = support.lsq.htscd
                              , lsq.formula = lsq.formula
                              , sig.formula = sig.formula
                              , data        = opt.data
                              , control     = ctrl.optim
                              , hessian     = TRUE
                              , ...
                              )#end optim
                       , silent = TRUE
                       )#end try

      #----- Check whether it worked. --------------------------------------------------#
      if ("try-error" %in% is(opt.1st)){
         success             = FALSE
      }else{
         #---------------------------------------------------------------------------------#
         #     Accept step only if it converged.                                           #
         #---------------------------------------------------------------------------------#
         success             = opt.1st$convergence %==% 0
         #---------------------------------------------------------------------------------#
      }#end if ("try-error" %in% is(opt.1st))
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    In case the first guess of sigma0 is far off, make it smaller and try again.    #
      #------------------------------------------------------------------------------------#
      if ((! success) && (! skip.sigma)){
         x.1st["sig.sigma0"] = x.1st["sig.sigma0"] / 2
      }#end if
      #------------------------------------------------------------------------------------#
   }#end while(! success & (i1st < i1st.max))
   #---------------------------------------------------------------------------------------#


   #----- Update both the first guess and the scale in case it converged. -----------------#
   if (success){
      x.1st      = opt.1st$par
      ctrl.optim = modifyList( x   = ctrl.optim
                             , val = list(parscale=ifelse(x.1st==0,1,abs(x.1st)))
                             )#end modifyList
      #----- Copy first guess to opt in case this is the only optimisation that works. ----#
      opt = opt.1st
      #------------------------------------------------------------------------------------#
   }#end if (success)
   #---------------------------------------------------------------------------------------#



   #----- Show first guess. ---------------------------------------------------------------#
   if (success && verbose){
      cat("             > First guess:   "
         , paste(paste(names(x.1st),sprintf("%.3f",opt.1st$par),sep=" = "),collapse=";   ")
         ,"\n",sep="")
   }else if (! success){
      #------------------------------------------------------------------------------------#
      #     Stop in case the first guess didn't work.                                      #
      #------------------------------------------------------------------------------------#
      if (is.debug){
         browser()
      }else{
         stop (" Solution of change didn't converge...")
      }#end if
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Initialise the output list.                                                       #
   #---------------------------------------------------------------------------------------#
   names.qual       = c("Estimate","Std. Error","t value","Pr(>|t|)")
   ans              = list()
   ans$call         = match.call()
   ans$err.method   = ifelse(err.short %in% "h","Hessian","Bootstrap")
   ans$lsq.formula  = lsq.formula
   ans$sig.formula  = sig.formula
   ans$df           = n.use - n.par
   ans$coefficients = matrix(data=NA,nrow=n.par,ncol=4,dimnames=list(names.par,names.qual))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Optimise the parameters.  Call optim multiple times until results are stable.     #
   #---------------------------------------------------------------------------------------#
   bef.support = -Inf
   it          = 0
   iterate     = TRUE
   nsteps      = 0
   while (iterate){
      it      = it + 1
      #------------------------------------------------------------------------------------#
      #     Call optimisation function.                                                    #
      #------------------------------------------------------------------------------------#
      opt.hess = try( optim( par         = x.1st
                           , fn          = support.lsq.htscd
                           , lsq.formula = lsq.formula
                           , sig.formula = sig.formula
                           , data        = opt.data
                           , control     = ctrl.optim
                           , hessian     = TRUE
                           , ...
                           )#end optim
                    , silent = TRUE
                    )#end try
      #------------------------------------------------------------------------------------#
      

      #------------------------------------------------------------------------------------#
      #     Decide what to do based on success/failure.                                    #
      #------------------------------------------------------------------------------------#
      if ("try-error" %in% is(opt.hess)){
         #----- Optimisation step failed.  Nudge 1st guess and try again. -----------------#
         success = FALSE
         #---------------------------------------------------------------------------------#
      }else{
         #---------------------------------------------------------------------------------#
         #     Accept step only if it converged and if the log-likelihood is negative.     #
         # The negative requirement is to make sure the step did not converge to a bogus   #
         # solution, as the likelihood is a product of probabilities, which should be less #
         # than 1, hence the negative requirement.                                         #
         #---------------------------------------------------------------------------------#
         success = opt.hess$convergence %==% 0
         #---------------------------------------------------------------------------------#
      }#end if ("try-error" %in% is(opt.hess))
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #    Update optimisation step only if the optimisation step converged and the        #
      # support improved.                                                                  #
      #------------------------------------------------------------------------------------#
      if (success){
         now.support = opt.hess$value
         nsteps.now  = opt.hess$counts["function"]
         nsteps      = nsteps + nsteps.now

         #----- Check whether we are improving things. ------------------------------------#
         if (is.finite(bef.support)){
            gain = 200. * ( (now.support - bef.support) 
                          / (abs(now.support) + abs(bef.support) ) )
         }else{
            gain = 200.
         }#end if
         #---------------------------------------------------------------------------------#


         #----- Check whether we should iterate. ------------------------------------------#
         iterate = ( gain > (100. * tol.gain) ) && it < maxit
         if (verbose){
            cat("             > Iteration ",it
               ,";   Converged: ", ! iterate
               ,";   Success: "  , success
               ,";   Support: "  , sprintf("%.3f",now.support)
               ,";   Gain: "     , sprintf("%.3f",gain       )   
               ,";   Steps: "    , sprintf("%6i" ,nsteps.now )
               ,"\n")
         }#end if
         #---------------------------------------------------------------------------------#



         #----- Update support. -----------------------------------------------------------#
         if (iterate){
            bef.support = now.support
            x.1st       = opt.hess$par
            ctrl.optim  = modifyList( x   = ctrl.optim
                                    , val = list(parscale=ifelse(x.1st==0,1,abs(x.1st)))
                                    )#end modifyList
         }else if (gain > 0){
            #----- Update opt only if this is an improvement. -----------------------------#
            opt = opt.hess
            #------------------------------------------------------------------------------#
         }#end if
         #---------------------------------------------------------------------------------#
      }else{
         #----- Optimisation step failed, quit iteration. ---------------------------------#
         iterate = FALSE
         if (verbose){
            cat( "             > Iteration ",it
               , ";   Optim call failed.  Exit improvement loop."
               , "\n"
               , sep=""
               )
         }#end if (verbose)
         #---------------------------------------------------------------------------------#
      }#end if (success)
      #------------------------------------------------------------------------------------#
   }#end while
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     List best guesses that work directly on evaluate.lsq.htscd.                       #
   #---------------------------------------------------------------------------------------#
   ans$x.best           = opt$par
   names(ans$x.best)    = names(x.1st)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Copy the results to ans.                                                          #
   #---------------------------------------------------------------------------------------#
   ans$hessian          = opt$hessian
   ans$coefficients[,1] = opt$par
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find errors associated with coefficients; decide whether to use Hessian matrix    #
   # or bootstrap.                                                                         #
   #---------------------------------------------------------------------------------------#
   if (err.short %in% "h"){
      #----- Use diagonal elements of Hessian matrix to estimate standard error. ----------#
      se                   = try(sqrt(diag(solve(-ans$hessian))),silent=TRUE)
      if ("try-error" %in% is(se)){
         warning(" Hessian matrix cannot be inverted.  Optimisation probably failed.")
      }else{
         names(se)            = NULL
         ans$coefficients[,2] = se
         ans$coefficients[,3] = ans$coefficients[,1] / ans$coefficients[,2]
         ans$coefficients[,4] = 2.0 * pt(-abs(ans$coefficients[,3]),df=ans$df)
      }#end if
      #------------------------------------------------------------------------------------#
   }else{
      #------------------------------------------------------------------------------------#
      #      Use bootstrap to estimate parameters and errors.                              #
      #------------------------------------------------------------------------------------#



      #----- Initialise arrays that will store the data. ----------------------------------#
      coeff.boot   = matrix( data     = NA
                           , nrow     = n.boot
                           , ncol     = n.par
                           , dimnames = list(NULL,names(x.1st))
                           )#end matrix
      support.boot = rep(NA,times=n.boot)
      xval.boot    = matrix( data = NA
                           , ncol = n.boot
                           , nrow = nrow(opt.data)
                           , dimnames = list(rownames(opt.data),NULL)
                           )#end matrix
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Loop until we reach the sought number of bootstrap realisations.  In case some #
      # realisation doesn't work, skip it.                                                 #
      #------------------------------------------------------------------------------------#
      ib           = 0
      while (ib < n.boot){
         #----- Select samples for this realisation. --------------------------------------#
         idx         = sample.int(n=n.use,replace=TRUE)
         ixval       = which(! (sequence(n.use) %in% idx))
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Call the optimisation procedure.                                            #
         #---------------------------------------------------------------------------------#
         opt.boot    = try( optim( par         = x.1st
                                 , fn          = support.lsq.htscd
                                 , lsq.formula = lsq.formula
                                 , sig.formula = sig.formula
                                 , data        = opt.data[idx,,drop=FALSE]
                                 , control     = ctrl.optim
                                 , hessian     = FALSE
                                 , ...
                                 )#end optim
                          , silent = TRUE
                          )#end try
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Check whether this bootstrap realisation worked.                            #
         #---------------------------------------------------------------------------------#
         if ("try-error" %in% is(opt.boot)){
            success     = FALSE
         }else{
            #------------------------------------------------------------------------------#
            #     Accept step only if it converged and if the log-likelihood is negative.  #
            # The negative requirement is to make sure the step did not converge to a      #
            # bogus solution, as the likelihood is a product of probabilities, which       #
            # should be less than 1, hence the negative requirement.                       #
            #------------------------------------------------------------------------------#
            success     = opt.boot$convergence %==% 0
            nsteps.now  = opt.boot$counts["function"]
            #------------------------------------------------------------------------------#
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Check whether to append to the data set.                                   #
         #---------------------------------------------------------------------------------#
         if (success){
            ib                = ib + 1
            coeff.boot  [ib,] = opt.boot$par
            support.boot[ib ] = opt.boot$value
            #----- Run cross validation. --------------------------------------------------#
            if (length(ixval) > 0){
               ypred  = evaluate.lsq.htscd( x           = opt.boot$par
                                          , lsq.formula = lsq.formula
                                          , sig.formula = sig.formula
                                          , data        = opt.data[ixval,,drop=FALSE]
                                          , ...
                                          )#end evaluate.lsq.htscd

               xval.boot[ixval,ib] = ypred$yhat
            }#end if
            #------------------------------------------------------------------------------#

            if (verbose){
               cat( "             > Bootstrap  ",ib
                  , ";   Support: ",sprintf("%.3f",opt.boot$value )
                  , ";   Parameters:  "
                  ,  paste( paste(names(x.1st),sprintf("%.3f",opt.boot$par),sep=" = ")
                          , collapse=";   "
                          )#end paste
                  , "\n"
                  , sep=""
                  )#end cat
            }#end if (verbose)
            #------------------------------------------------------------------------------#
         }else if (verbose){
            cat("             > Bootstrap  realisation failed, skip it.","\n",sep="")
         }#end if (success)
         #---------------------------------------------------------------------------------#
      }#end while (ib < n.boot)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Copy the results to ans.                                                       #
      #------------------------------------------------------------------------------------#
      ans$coefficients[,2] = apply(X=coeff.boot,MARGIN=2,FUN=sd  )
      ans$coefficients[,3] = ans$coefficients[,1] / ans$coefficients[,2]
      ans$coefficients[,4] = 2.0 * pt(-abs(ans$coefficients[,3]),df=ans$df)
      ans$coeff.boot       = coeff.boot
      ans$support.boot     = support.boot
      #------------------------------------------------------------------------------------#
   }#end if (err.short %in% "h")
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Keep only the variables in data that are needed for the fitting.                  #
   #---------------------------------------------------------------------------------------#
   names.data   = names(data)
   names.data   = names.data[names.data %in% c(lsq.vars,sig.vars)]
   out.data     = data[,names.data,drop=FALSE]
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Find predicted values, sigma, and additional evaluations of goodness of fit.      #
   #---------------------------------------------------------------------------------------#
   ypred             = evaluate.lsq.htscd( x           = ans$x.best
                                         , lsq.formula = lsq.formula
                                         , sig.formula = sig.formula
                                         , data        = out.data
                                         , ...
                                         )#end evaluate.lsq.htscd
   ans$support       = support.lsq.htscd( x           = ans$x.best
                                        , lsq.formula = lsq.formula
                                        , sig.formula = sig.formula
                                        , data        = out.data
                                        , ...
                                        )#end support.lsq.htscd
   ans$data          = out.data
   ans$actual.values = ypred$yact
   ans$fitted.values = ypred$yhat
   ans$residuals     = ypred$yres
   ans$sigma         = ypred$sigma
   ans$goodness      = test.goodness( x.mod        = ans$fitted.values
                                    , x.obs        = ans$actual.values
                                    , n.parameters = length(lsq.first)
                                    )#end test.goodness
   ans$AIC           = 2*n.par - 2*ans$support + 2*n.par*(n.par+1)/(n.use-n.par-1)
   ans$BIC           = n.par*(log(n.use)-log(2*pi)) - 2*ans$support
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Cross-validation (Bootstrap only).                                                #
   #---------------------------------------------------------------------------------------#
   if (err.short %in% "b"){
      #------------------------------------------------------------------------------------#
      #     Find summary statistics for cross-validation.                                  #
      #------------------------------------------------------------------------------------#
      qlow          = 0.5 - 0.5 * ci.level
      qhigh         = 0.5 + 0.5 * ci.level
      xval.boot     = ifelse(is.finite(xval.boot),xval.boot,NA)
      xval.resid    = - apply(X=xval.boot ,MARGIN=2,FUN="-", opt.data[[yname]])
      n.xval        =   apply(X=xval.boot ,MARGIN=1,FUN=function(x) sum(! is.na(x)))
      expect.xval   =   apply(X=xval.boot ,MARGIN=1,FUN=mean                ,na.rm=TRUE)
      qlow.xval     =   apply(X=xval.boot ,MARGIN=1,FUN=quantile,probs=qlow ,na.rm=TRUE)
      qhigh.xval    =   apply(X=xval.boot ,MARGIN=1,FUN=quantile,probs=qhigh,na.rm=TRUE)
      bias.xval     = - apply(X=xval.resid,MARGIN=1,FUN=mean                ,na.rm=TRUE)
      sigma.xval    =   apply(X=xval.resid,MARGIN=1,FUN=sd                  ,na.rm=TRUE)
      rmse.xval     =   sqrt(bias.xval^2+sigma.xval^2)
      #------------------------------------------------------------------------------------#

      #------------------------------------------------------------------------------------#
      #     Attach to answer.                                                              #
      #------------------------------------------------------------------------------------#
      ans$conf.int  = ci.level
      ans$cross.val = data.frame( n        = n.xval
                                , expected = expect.xval
                                , qlow     = qlow.xval
                                , qhigh    = qhigh.xval
                                , bias     = bias.xval
                                , sigma    = sigma.xval
                                , rmse     = rmse.xval
                                )#end data.frame
      ans$xval.mat  = xval.boot
      #------------------------------------------------------------------------------------#
   }#end if (err.short %in% "b")
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Print results on standard output.                                                 #
   #---------------------------------------------------------------------------------------#
   if (verbose){
      cat("             > Optimised values:   "
         , paste(paste(names(x.1st),sprintf("%.3f",opt$par),sep=" = "),collapse=";   ")
         ,"\n",sep="")
   }#end if(verbose)
   #---------------------------------------------------------------------------------------#

   class(ans) = c("lsq.htscd")

   #----- Return answer -------------------------------------------------------------------#
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end optim.lsq.htscd
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#      Function to check whether the object came from heteroscedastic fit.                 #
#------------------------------------------------------------------------------------------#
is.lsq.htscd <<- function(x){
   ans = inherits(x,"lsq.htscd")
   return(ans)
}#end is.lsq.htscd
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#    predict.lsq.htscd -- This function predicts the model using the object.               #
#------------------------------------------------------------------------------------------#
predict.lsq.htscd <<- function( object
                              , newdata   = NULL
                              , confint   = FALSE
                              , level     = 0.95
                              , se.fit    = FALSE
                              , pred.boot = FALSE
                              , ...
                              ){
   #---------------------------------------------------------------------------------------#
   #    This must use an object created by optim.lsq.htscd.                                #
   #---------------------------------------------------------------------------------------#
   stopifnot(is.lsq.htscd(object))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     A few handy aliases.                                                              #
   #---------------------------------------------------------------------------------------#
   lsq.formula = object$lsq.formula
   sig.formula = object$sig.formula
   if (is.null(newdata)){
      data = object$data
   }else{
      data = newdata
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Prepare the vector to go to the model evaluation.                                 #
   #---------------------------------------------------------------------------------------#
   x     = object$x.best
   now   = evaluate.lsq.htscd( x           = x
                             , lsq.formula = lsq.formula
                             , sig.formula = sig.formula
                             , data        = data
                             , ...
                             )#end evaluate.lsq.htscd
   yhat  = now$yhat
   yact  = now$yact
   yres  = now$yres
   sigma = now$sigma
   #---------------------------------------------------------------------------------------#



   #----- Split the coefficients. ---------------------------------------------------------#
   x.lsq        = x[grepl(pattern="^lsq\\.",x=names(x))]
   x.sig        = x[grepl(pattern="^sig\\.",x=names(x))]
   names(x.lsq) = gsub(pattern="^lsq\\.",replacement="",x=names(x.lsq))
   names(x.sig) = gsub(pattern="^sig\\.",replacement="",x=names(x.sig))
   n.par        = length(x.lsq)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Calculate confidence intervals?                                                    #
   #---------------------------------------------------------------------------------------#
   if ( (confint || se.fit || pred.boot) && (object$err.method %in% "Hessian")){
      #------------------------------------------------------------------------------------#
      #    Haven't figured out how to estimate CI bands using Hessian, warn the user.      #
      #------------------------------------------------------------------------------------#
      if (confint){
         ylwr   = yhat * NA
         yupr   = yhat * NA
         yfit   = matrix( data = cbind(yhat,ylwr,yupr)
                        , ncol = 3
                        , nrow = length(yhat)
                        , dimnames = list(rownames(data),c("fit","lwr","upr"))
                        )#end matrix
      }#end if (confint)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    Haven't figured out how to estimate SE using Hessian, warn the user.            #
      #------------------------------------------------------------------------------------#
      if (se.fit){
         yste        = NA
         names(yste) = rownames(data)
      }#end if (sefit)
      #------------------------------------------------------------------------------------#

      warning("Confidence interval and SE are currently not available for Hessian solver.")

   }else if (confint || se.fit || pred.boot){
      coeff.boot = object$coeff.boot
      nboot      = nrow(coeff.boot)
      yhat.boot  = matrix(nrow=nrow(data),ncol=nboot,dimnames=list(rownames(data),NULL))
      for (ib in sequence(nboot)){
         xnow           = coeff.boot[ib,]
         now            = evaluate.lsq.htscd( x           = xnow
                                            , lsq.formula = lsq.formula
                                            , sig.formula = sig.formula
                                            , data        = data
                                            , ...
                                            )#end evaluate.lsq.htscd
         yhat.boot[,ib] = now$yhat
      }#end for (ib in sequence(nboot))
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Using quantiles to obtain estimates for lower and upper limits.                #
      #------------------------------------------------------------------------------------#
      if (confint){
         ci.lwr = 0.5 - 0.5 * level
         ci.upr = 0.5 + 0.5 * level
         ylwr   = apply(X=yhat.boot,MARGIN=1,FUN=quantile,probs=ci.lwr,na.rm=TRUE)
         yupr   = apply(X=yhat.boot,MARGIN=1,FUN=quantile,probs=ci.upr,na.rm=TRUE)
         yfit   = matrix( data = cbind(yhat,ylwr,yupr)
                        , ncol = 3
                        , nrow = length(yhat)
                        , dimnames = list(rownames(data),c("fit","lwr","upr"))
                        )#end matrix
      }#end if
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Standard error obtained by standard deviation of the realisation estimates.    #
      #------------------------------------------------------------------------------------#
      if (se.fit){
         yste        = apply(X=yhat.boot,MARGIN=1,FUN=sd,na.rm=TRUE)
         names(yste) = rownames(data)
      }#end if
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find number of degrees of freedom and residual scale.                             #
   #---------------------------------------------------------------------------------------#
   if (se.fit){
      y.df = nrow(data) - n.par
      yrsc = sqrt(sum((yres/sigma)^2)/y.df)
   }#end if (se.fit)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Define answer depending on the requests.                                          #
   #---------------------------------------------------------------------------------------#
   if (confint && se.fit && pred.boot){
      ans = list(fit=yfit,se.fit=yste,boot=yhat.boot,df=y.df,residual.scale=yrsc)
   }else if (confint && se.fit){
      ans = list(fit=yfit,se.fit=yste,df=y.df,residual.scale=yrsc)
   }else if (confint && pred.boot){
      ans = list(fit=yfit,boot=yhat.boot)
   }else if (se.fit && pred.boot){
      ans = list(fit=yhat,boot=yhat.boot,se.fit=yste,df=y.df,residual.scale=yrsc)
   }else if (confint){
      ans = yfit
   }else if (se.fit){
      ans = list(fit=yhat,se.fit=yste,df=y.df,residual.scale=yrsc)
   }else if (pred.boot){
      ans = list(fit=yhat,boot=yhat.boot)
   }else{
      ans = yhat
   }#end if (confint && se.fit)
   #---------------------------------------------------------------------------------------#


   return(ans)
}#end predict.lsq.htscd
#------------------------------------------------------------------------------------------#





#==========================================================================================#
#==========================================================================================#
#    fitted.lsq.htscd -- This function returns the fitted values using the object.         #
#------------------------------------------------------------------------------------------#
fitted.lsq.htscd <<- function(object,...){
   #---------------------------------------------------------------------------------------#
   #    This must use an object created by optim.lsq.htscd.                                #
   #---------------------------------------------------------------------------------------#
   stopifnot(is.lsq.htscd(object))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     A few handy aliases.                                                              #
   #---------------------------------------------------------------------------------------#
   lsq.formula = object$lsq.formula
   sig.formula = object$sig.formula
   data        = newdata
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Prepare the vector to go to the model evaluation.                                 #
   #---------------------------------------------------------------------------------------#
   x   = object$x.best
   now = evaluate.lsq.htscd( x           = x
                           , lsq.formula = lsq.formula
                           , sig.formula = sig.formula
                           , data        = data
                           , ...
                           )#end evaluate.lsq.htscd
   ans = now$yhat
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end fitted.lsq.htscd
#------------------------------------------------------------------------------------------#





#==========================================================================================#
#==========================================================================================#
#    residuals.lsq.htscd -- This function predicts the model using the object.             #
#------------------------------------------------------------------------------------------#
residuals.lsq.htscd <<- function(object,...){
   #---------------------------------------------------------------------------------------#
   #    This must use an object created by optim.lsq.htscd.                                #
   #---------------------------------------------------------------------------------------#
   stopifnot(is.lsq.htscd(object))
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     A few handy aliases.                                                              #
   #---------------------------------------------------------------------------------------#
   lsq.formula = object$lsq.formula
   sig.formula = object$sig.formula
   data        = object$data
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Prepare the vector to go to the model evaluation.                                 #
   #---------------------------------------------------------------------------------------#
   x   = object$x.best
   now = evaluate.lsq.htscd( x           = x
                           , lsq.formula = lsq.formula
                           , sig.formula = sig.formula
                           , data        = data
                           , ...
                           )#end evaluate.lsq.htscd
   ans = now$yres
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end residuals.lsq.htscd
#------------------------------------------------------------------------------------------#





#==========================================================================================#
#==========================================================================================#
#    coef.lsq.htscd -- This function retrieves the model coefficients using the object.    #
#------------------------------------------------------------------------------------------#
coef.lsq.htscd <<- function(object){
   #---------------------------------------------------------------------------------------#
   #    This must use an object created by optim.lsq.htscd.                                #
   #---------------------------------------------------------------------------------------#
   stopifnot(is.lsq.htscd(object))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Return the vector with coefficients.                                              #
   #---------------------------------------------------------------------------------------#
   ans = object$x.best
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end coef.lsq.htscd
#------------------------------------------------------------------------------------------#





#==========================================================================================#
#==========================================================================================#
#    formula.lsq.htscd -- This function retrieves the model formulae using the object.     #
#------------------------------------------------------------------------------------------#
formula.lsq.htscd <<- function(object){
   #---------------------------------------------------------------------------------------#
   #    This must use an object created by optim.lsq.htscd.                                #
   #---------------------------------------------------------------------------------------#
   stopifnot(is.lsq.htscd(object))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Return a list with the formulae.                                                  #
   #---------------------------------------------------------------------------------------#
   ans = list(lsq = object$lsq.formula, sig = object$sig.formula)
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end formula.lsq.htscd
#------------------------------------------------------------------------------------------#





#==========================================================================================#
#==========================================================================================#
#    logLik.lsq.htscd -- This function retrieves the support function.                     #
#------------------------------------------------------------------------------------------#
logLik.lsq.htscd <<- function(object){
   #---------------------------------------------------------------------------------------#
   #    This must use an object created by optim.lsq.htscd.                                #
   #---------------------------------------------------------------------------------------#
   stopifnot(is.lsq.htscd(object))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Return the log-likelihood (aka support).                                          #
   #---------------------------------------------------------------------------------------#
   ans = object$support
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end logLik.lsq.htscd
#------------------------------------------------------------------------------------------#





#==========================================================================================#
#==========================================================================================#
#    AIC.lsq.htscd -- This function retrieves the Akaike Information Criterion.            #
#------------------------------------------------------------------------------------------#
AIC.lsq.htscd <<- function(object){
   #---------------------------------------------------------------------------------------#
   #    This must use an object created by optim.lsq.htscd.                                #
   #---------------------------------------------------------------------------------------#
   stopifnot(is.lsq.htscd(object))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Return a list with the AIC value.                                                 #
   #---------------------------------------------------------------------------------------#
   ans = object$AIC
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end AIC.lsq.htscd
#------------------------------------------------------------------------------------------#





#==========================================================================================#
#==========================================================================================#
#    BIC.lsq.htscd -- This function retrieves the Bayes Information Criterion.             #
#------------------------------------------------------------------------------------------#
BIC.lsq.htscd <<- function(object){
   #---------------------------------------------------------------------------------------#
   #    This must use an object created by optim.lsq.htscd.                                #
   #---------------------------------------------------------------------------------------#
   stopifnot(is.lsq.htscd(object))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Return a list with the BIC value.                                                 #
   #---------------------------------------------------------------------------------------#
   ans = object$BIC
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end BIC.lsq.htscd
#------------------------------------------------------------------------------------------#





#==========================================================================================#
#==========================================================================================#
#    summary.lsq.htscd -- Dummy function, it returns the object itself.                    #
#------------------------------------------------------------------------------------------#
summary.lsq.htscd <<- function(object){
   #---------------------------------------------------------------------------------------#
   #    This must use an object created by optim.lsq.htscd.                                #
   #---------------------------------------------------------------------------------------#
   stopifnot(is.lsq.htscd(object))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Return the vector with coefficients.                                              #
   #---------------------------------------------------------------------------------------#
   return(object)
   #---------------------------------------------------------------------------------------#
}#end summary.lsq.htscd
#------------------------------------------------------------------------------------------#





#==========================================================================================#
#==========================================================================================#
#    print.lsq.htscd -- Prints information on screen.                                      #
#------------------------------------------------------------------------------------------#
print.lsq.htscd <<- function(object){
   #---------------------------------------------------------------------------------------#
   #    This must use an object created by optim.lsq.htscd.                                #
   #---------------------------------------------------------------------------------------#
   stopifnot(is.lsq.htscd(object))
   #---------------------------------------------------------------------------------------#


   p.value         = object$coefficients[,4]
   coeff.table     = as.data.frame(object$coefficients)
   coeff.table[,1] = sprintf("%g",signif(coeff.table[,1],5))
   coeff.table[,2] = sprintf("%g",signif(coeff.table[,2],5))
   coeff.table[,3] = sprintf("%g",signif(coeff.table[,3],4))
   coeff.table[,4] = ifelse( p.value %<% 1.e-16,"< 1e-16",sprintf("%g",signif(p.value,3)))

   #----- Append the significance test. ---------------------------------------------------#
   sig.brks           = c(-Inf,0.001,0.01,0.05,0.1,Inf)
   sig.symbs          = c("***","**","*",".","")
   coeff.table$signif = sig.symbs[as.numeric(cut(p.value,breaks=sig.brks))]
   kplot              = paste0("p-value: 0 \'***\' 0.001 \'**\' 0.01"
                              ," \'*\' 0.05 \'.\' 0.1 \' \' 1")
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Some metrics showing goodness of fit.                                             #
   #---------------------------------------------------------------------------------------#
   bias      = sprintf("%g",signif(object$goodness$bias     ,4))
   sigma     = sprintf("%g",signif(object$goodness$sigma    ,4))
   rmse      = sprintf("%g",signif(object$goodness$rmse     ,4))
   r.squared = sprintf("%g",signif(object$goodness$r.squared,3))
   n         = object$goodness$n
   dfres     = object$goodness$df.err
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      List what to show.                                                               #
   #---------------------------------------------------------------------------------------#
   cat0("")
   cat0("Call:")
   print(object$call,quote=FALSE)
   cat0("")
   cat0("Coefficients:")
   print(coeff.table,quote=FALSE)
   cat0(kplot)
   cat0("")
   cat0("Degrees of freedom: ",n," total; ",dfres," residual")
   cat0("Bias: ",bias,"; Sigma: ",sigma,"; RMSE: ",rmse,"; r.squared:",r.squared)
   cat0("")
   #---------------------------------------------------------------------------------------#




   #----- Set result to be invisible. -----------------------------------------------------#
   invisible()
   #---------------------------------------------------------------------------------------#
}#end summary.lsq.htscd
#------------------------------------------------------------------------------------------#




#==========================================================================================#
#==========================================================================================#
#    evaluate.lsq.htscd -- This function evaluates the model using the passed parameters.  #
#                          It returns a data frame with four variables:                    #
#                          yhat  -- The predicted values.                                  #
#                          yact  -- The actual values.                                     #
#                          yres  -- Residuals                                              #
#                          sigma -- Local estimate of sigma.                               #
#------------------------------------------------------------------------------------------#
evaluate.lsq.htscd <<- function(x,lsq.formula,sig.formula,data,...){

   #----- Split the coefficients. ---------------------------------------------------------#
   x.lsq        = x[grepl(pattern="^lsq\\.",x=names(x))]
   names(x.lsq) = gsub(pattern="^lsq\\.",replacement="",x=names(x.lsq))
   if (is.null(sig.formula)){
      x.sig        = NULL
   }else{
      x.sig        = x[grepl(pattern="^sig\\.",x=names(x))]
      names(x.sig) = gsub(pattern="^sig\\.",replacement="",x=names(x.sig))
   }#end if (! is.null(sig.formula))
   n.par        = length(x.lsq)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Coerce data to a data frame.                                                       #
   #---------------------------------------------------------------------------------------#
   data  = as.data.frame(data)
   n.use = sum(apply(X=data,MARGIN=1,FUN=function(x) all(is.finite(x))))
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Initialise a new environment where calculations will occur.                       #
   #---------------------------------------------------------------------------------------#
   modeval = new.env()
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Get additional arguments.                                                          #
   #---------------------------------------------------------------------------------------#
   dotdotdot = list(...)
   ndots     = length(dotdotdot)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Go through all variables in x.lsq, x.sig, data, and ... and assign them to the    #
   # new environment.                                                                      #
   #---------------------------------------------------------------------------------------#
   for (i in seq_along(x.lsq)){
      dummy = assign(x=names(x.lsq)[i],value=x.lsq[i]  ,envir=modeval)
   }#end for (i in seq_along(x.lsq))
   for (i in seq_along(x.sig)){
      dummy = assign(x=names(x.sig)[i],value=x.sig[i]  ,envir=modeval)
   }#end for (i in seq_along(x.sig))
   for (j in seq_along(data)){
      dummy = assign(x=names(data)[j],value=data[[j]]  ,envir=modeval)
   }#end for (j in seq_along(lsq.data))
   for (j in seq_along(dotdotdot)){
      dummy = assign(x=names(dotdotdot)[j],value=dotdotdot[[j]],envir=modeval)
   }#end for (j in seq_along(dotdotdot))
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Find the expressions we shall evaluate to determine yhat and sigma.               #
   #---------------------------------------------------------------------------------------#
   lsq.expr = attr(x=terms(lsq.formula),which="term.labels")
   if (! is.null(sig.formula)){
      sig.expr = attr(x=terms(sig.formula),which="term.labels")
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Find yhat (fitted) and yres (residuals). ----------------------------------------#
   yhat     = eval(expr=parse(text=lsq.expr),envir=modeval)
   resp     = all.vars(lsq.formula)[1]
   if (resp %in% names(data)){
      yact  = data[[resp]]
   }else{
      yact  = yhat * NA
   }#end if
   yres     = yact - yhat
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Add yhat and yres to modeval so the distribution of residuals can also use these  #
   # values.  The only exception is if variables with such names exist in the data frame,  #
   # in which case we don't overwrite them.                                                #
   #---------------------------------------------------------------------------------------#
   if (! "yhat" %in% c(names(data),names(dotdotdot))){
      dummy = assign("yhat",value=yhat,envir=modeval)
   }#end if (! "yhat" %in% c(names(data),names(dotdotdot)))
   if (! "yres" %in% c(names(data),names(dotdotdot))){
      dummy = assign("yres",value=yres,envir=modeval)
   }#end if (! "yres" %in% c(names(data),names(dotdotdot)))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the scale for sigma (aka sigma0), and add to the model evaluation            #
   # environment.                                                                          #
   #---------------------------------------------------------------------------------------#
   if (is.null(sig.formula)){
      sigma   = rep(sqrt( sum(yres^2,na.rm=TRUE) / (n.use - n.par) ),times=length(yres))
      dummy   = assign(x="sigma",value=sigma  ,envir=modeval)
   }else{
      sigma   = eval(expr=parse(text=sig.expr),envir=modeval)
   }#end if
   #---------------------------------------------------------------------------------------#

   #----- Append data to a data frame. ----------------------------------------------------#
   ans = data.frame( yhat = yhat, yact = yact, yres = yres, sigma = sigma)
   #---------------------------------------------------------------------------------------#


   return(ans)
}#end evaluate.lsq.htscd
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#    support.lsq.htscd -- The function that evaluates the support for a given set of       #
#                         parameters, using the least square estimator.  This is           #
#                         essentially the same as a "normal" least squares, except that it #
#                         uses different values of sigma for each data point.              #
#------------------------------------------------------------------------------------------#
support.lsq.htscd <<- function(x,lsq.formula,sig.formula,data,...){

   #---- Predict the values. --------------------------------------------------------------#
   ypred   = evaluate.lsq.htscd(x,lsq.formula,sig.formula,data,...)
   lnprob  = dnorm(x=ypred$yres,mean=0,sd=ypred$sigma,log=TRUE)
   support = sum(lnprob)
   return(support)
   #---------------------------------------------------------------------------------------#
}#end fit.nls.htscd
#==========================================================================================#
#==========================================================================================#
