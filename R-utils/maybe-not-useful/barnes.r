#==========================================================================================#
#==========================================================================================#
#     Define some global constants.                                                        #
#------------------------------------------------------------------------------------------#
gamma0     <<- 0.20
min.weight <<- 1.e-6
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     This is the main function that controls the Barnes' objective analysis based on the  #
# following paper:                                                                         #
#                                                                                          #
# Koch, S. E.; M. desJardins, P. J. Kocin, 1983: An interactive Barnes objective map       #
#    analysis scheme for use with satellite and conventional data. J. Climate and Appl.    #
#    Meteorol., 22, 1487-1503. (hereafter KJK83).                                          #
#                                                                                          #
# Inputs variables:                                                                        #
#    ref.lola - A vector containing the longitude and the latitude of the sought point.    #
#    rem.lola - A two-column matrix, each row containing the longitude and the latitude of #
#               every column site.                                                         #
#    rem.dat  - A multiple column matrix, each column containing the dataset to be         #
#               analysed.  The columns of rem.dat MUST correspond to the rows of rem.lola. #
#              (rows are times, columns are sites).                                        #
#    dist    - A vector containing the distances between the sites and the sought point.   #
#    gam     - The gamma parameter as in Koch et al. (1983).  If not given, we will use    #
#              the default value.                                                          #
#    verbose - A logical flag.  If TRUE then we will print more output to screen.          #
#                                                                                          #
# Output:                                                                                  #
#    gamma0        - The decaying term gamma as in KJK83.                                  #
#    delta.n       - Typical distance between sites.                                       #
#    dist.mat      - Matrix with distance between sites                                    #
#    first         - Result from the first step                                            #
#    second        - Result from the second step                                           #
#    nvalid        - Number of valid points for each time of the objective analysis.       #
#    fitted.values - Actual answer from the objective analysis.                            #
#------------------------------------------------------------------------------------------#
obj.analysis <<- function(ref.lola,rem.lola,rem.dat,gam=NULL,verbose=FALSE){

    #--------------------------------------------------------------------------------------#
    #     Fill in gamma and the typical dn with defaults in case nothing has been          #
    # provided.                                                                            #
    #--------------------------------------------------------------------------------------#
    if (is.null(gam)) gam = gamma0
    #--------------------------------------------------------------------------------------#



    #--------------------------------------------------------------------------------------#
    #     Check whether the dimensions of the dataset match with the vector of distances.  #
    #--------------------------------------------------------------------------------------#
    if (verbose) print(paste("    * Checking dimensions...",sep=""))
    nrdat  = nrow(rem.dat )
    ncdat  = ncol(rem.dat )
    nrlola = nrow(rem.lola)
    nclola = ncol(rem.lola)
    if (ncdat != nrlola | nclola != 2){
       print(paste("------------------------------------------"))
       print(paste(" In function obj.analysis:"                ))
       print(paste(" ---> # of rows of rem.dat:      ",nrdat ,sep=""))
       print(paste(" ---> # of columns of rem.dat:   ",ncdat ,sep=""))
       print(paste(" ---> # of rows of rem.lola:     ",nrlola,sep=""))
       print(paste(" ---> # of columns of rem.lola:  ",nrlola,sep=""))
       print(paste("------------------------------------------"))
       print(paste(" 1.  Number of columns of rem.dat must match the number"
                      ," of rows of rem.lola!",sep=""))
       print(paste(" 1.  Number of columns of rem.lola must be two!",sep=""))
       stop("Invalid input")
    }else if (length(ref.lola) != 2){
       stop(" Ref.lola must have two elements (longitude and latitude, in this order)!")
    }#end if
    #--------------------------------------------------------------------------------------#



    #----- Initialise the answer with NA. -------------------------------------------------#
    first.all   = rep(NA,times=nrdat)
    second.all  = rep(NA,times=nrdat)
    ans         = rep(NA,times=nrdat)
    #--------------------------------------------------------------------------------------#



    #--------------------------------------------------------------------------------------#
    #     Make sure that the reference coordinates are a matrix with two columns and one   #
    # row.                                                                                 #
    #--------------------------------------------------------------------------------------#
    ref.lola = matrix(ref.lola,ncol=2,nrow=1)
    #--------------------------------------------------------------------------------------#



    #--------------------------------------------------------------------------------------#
    #     Find the vector of distances between the remote sites and the reference.         #                        
    #--------------------------------------------------------------------------------------#
    if (verbose){
       print(paste("    * Finding distance between remote sites"
                      ,"and the reference site...",sep=""))
    }#end if
    rem.dist = rdist.earth(x1=rem.lola,x2=ref.lola,miles=FALSE) * 1000.
    rem.dist = c(rem.dist)
    #--------------------------------------------------------------------------------------#



    #--------------------------------------------------------------------------------------#
    #     Find the vector of distances amongst the remote sites, and find the mean minimum #
    # distance between the sites: that will be our proxy for delta.n.                      #
    #--------------------------------------------------------------------------------------#
    if (verbose)  print(paste("    * Finding typical distance between sites...",sep=""))
    mat.dist       = rdist.earth(x1=rem.lola,miles=FALSE) * 1000.
    #----- Make the diagonals NA because so the site itself cannot be its closest site. ---#
    diag(mat.dist) = NA
    delta.n  = mean(apply(X=mat.dist,MARGIN=2,FUN=min,na.rm=TRUE))
    #--------------------------------------------------------------------------------------#



    #----- Map the data that has at least one valid data point. ---------------------------#
    if (verbose){
       print(paste("    * Flag times when at least one point was valid...",sep=""))
    }#end if
    valid.dat = rowSums(is.finite(rem.dat))
    ok        = valid.dat != 0
    nok       = sum(ok)
    #--------------------------------------------------------------------------------------#



    #--------------------------------------------------------------------------------------#
    #     Delete the lines without enough data.                                            #
    #--------------------------------------------------------------------------------------#
    if (verbose) print(paste("    * Select lines with at least one valid datum...",sep=""))
    use.dat  = rem.dat [ok,]
    #--------------------------------------------------------------------------------------#



    #--------------------------------------------------------------------------------------#
    #     Convert the distance vector into a matrix of the same size as the data           #
    #--------------------------------------------------------------------------------------#
    if ( (ncdat*nok) %% length(rem.dist) != 0) browser(text="rem.dist")
    use.dist           = matrix(rem.dist,ncol=ncdat,nrow=nok,byrow=TRUE)
    dimnames(use.dist) = dimnames(use.dat)
    #--------------------------------------------------------------------------------------#



    #--------------------------------------------------------------------------------------#
    #     Find the default kappa0 and kappa1 based on equations 11 of KJK83.                           #
    #--------------------------------------------------------------------------------------#
    if (verbose) print(paste("    * Find the decay coefficients...",sep=""))
    kappa0 = kappa.barnes(gam=gam,dn=delta.n)
    kappa1 = gamma0 * kappa0
    #--------------------------------------------------------------------------------------#


    #--------------------------------------------------------------------------------------#
    #     Find the first guess, which is a simple weighted average, with weights defined   #
    # as functions of the distance to the point.                                           #
    #--------------------------------------------------------------------------------------#
    if (verbose){
       print(paste("    * Finding the first guess of the objective analysis...",sep=""))
    }#end if
    first = oa.weighted.mean(dat=use.dat,dist=use.dist,kap=kappa0)
    #--------------------------------------------------------------------------------------#



    #--------------------------------------------------------------------------------------#
    #     Now we must estimate the error that a first guess would cause to the other grid  #
    # points by estimating the first guess for each of the points that we have, and        #
    # replacing its measurement by our first guess.                                        #
    #--------------------------------------------------------------------------------------#
    if (verbose) print(paste("    * Finding the nudging term...",sep=""))
    #----- Initialise the nudge factor with NA. -------------------------------------------#
    use.nudge = NA * use.dat
    #----- Loop over the sites. -----------------------------------------------------------# 
    for (p in 1:nrlola){
       p.iata = dimnames(rem.lola)[[1]][p]
       if (verbose) print(paste("      # Replacing column ",p," (",p.iata,")...",sep=""))

       #-----------------------------------------------------------------------------------# 
       #      Pretend that the reference site is site [p] and site [p] is the reference    #
       # site.                                                                             #
       #-----------------------------------------------------------------------------------#
       p.ref.lola     = matrix(rem.lola[p,],ncol=2,nrow=1)
       p.rem.lola     = rem.lola
       p.rem.lola[p,] = ref.lola
       dimnames(p.rem.lola)[[1]][p] = "1st"
       #-----------------------------------------------------------------------------------#



       #-----------------------------------------------------------------------------------#
       #     Delete the lines without enough data.                                         #
       #-----------------------------------------------------------------------------------#
       p.use.dat     = use.dat
       p.use.dat[,p] = first
       dimnames(p.use.dat)[[2]][p] = "1st"
       #-----------------------------------------------------------------------------------#


       #-----------------------------------------------------------------------------------#
       #     Find the vector of distances between the remote sites and the reference, and  #
       # the corresponding matrix.                                                         #
       #-----------------------------------------------------------------------------------#
       p.rem.dist           = rdist.earth(x1=p.rem.lola,x2=p.ref.lola,miles=FALSE) * 1000.
       p.rem.dist           = c(p.rem.dist)
       if ( (ncdat*nok) %% length(p.rem.dist) != 0) browser(text="p.rem.dist")
       p.use.dist           = matrix(p.rem.dist,ncol=ncdat,nrow=nok,byrow=TRUE)
       dimnames(p.use.dist) = dimnames(p.use.dat)
       #-----------------------------------------------------------------------------------#
 

       #-----------------------------------------------------------------------------------#
       #     Find the first guess for the reference sites, which is a simple weighted      #
       # average, with weights defined as functions of the distance to the point.          #
       #-----------------------------------------------------------------------------------#
       p.first = oa.weighted.mean(dat=p.use.dat,dist=p.use.dist,kap=kappa0)
       #-----------------------------------------------------------------------------------#



       #-----------------------------------------------------------------------------------#
       #     Replace the column of the nudging factor with the difference between the      #
       # actual observation and the result from the first guess.                           #
       #-----------------------------------------------------------------------------------#
       use.nudge[,p] = use.dat[,p] - p.first
       #-----------------------------------------------------------------------------------#
    }#end for

    #--------------------------------------------------------------------------------------#
    #     The second step is just a weighted average of the nudging factor, with weights   #
    # proportional to the distance but with a lower kappa.                                 #
    #--------------------------------------------------------------------------------------#
    if (verbose){
       print(paste("    * Finding the second guess of the objective analysis...",sep=""))
    }#end if
    second = oa.weighted.mean(dat=use.nudge,dist=use.dist,kap=kappa1)
    #--------------------------------------------------------------------------------------#



    #--------------------------------------------------------------------------------------#
    #     The result from the objective analysis is the sum of the first and second        #
    # iterations.  Here we put the answers on the points that we could run, leaving the    #
    # other variables NA.                                                                  #
    #--------------------------------------------------------------------------------------#
    first.all  [ok] = first
    second.all [ok] = second
    ans        [ok] = first + second
    #--------------------------------------------------------------------------------------#




    #--------------------------------------------------------------------------------------#
    #      The result will be a list with some of the intermediate variables and the final #
    # answer for the objective analysis.                                                   #
    #--------------------------------------------------------------------------------------#
    res = list(gamma=gamma0,kappa=kappa0,delta.n=delta.n,dist=mat.dist
              ,first=first.all,second=second.all,nvalid=valid.dat,fitted.values=ans)
    return(res)
    #--------------------------------------------------------------------------------------#
}#end function
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     This function will define the best kappa0 based on D0(2 Delta_n) for the Barnes      #
# objective analysis as in Koch et al. (1983).  The value of D0 is chosen to give          #
# D1*=exp(-1) for a given gamma.                                                           #
#                                                                                          #
#        D1* = D0 * (1 + D0^(gamma-1) - D0^gamma) (Koch et al, 1983, equation 11.)         #
#------------------------------------------------------------------------------------------#
kappa.barnes = function(gam,dn){


   #----- Define the function for which we will seek the root. ----------------------------#
   zerofun = function(x,gam) x * (1 + x^(gam-1) - x^gam) - exp(-1.)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find a good combination of first guesses that have opposite signal.  D0 must be   #
   # between 0 and 1, so we define the interval as being between 1.e-10 and 1.             #
   #---------------------------------------------------------------------------------------#
   d0.barnes = uniroot(f=zerofun,interval=c(1.e-10,1.),gam=gam)$root
   lnd02dtn = log(d0.barnes)
   #---------------------------------------------------------------------------------------#


   #----- Define kappa_0. -----------------------------------------------------------------#
   kappa.zero = - (2. * dn / pi) * (2. * dn / pi) * lnd02dtn
   #---------------------------------------------------------------------------------------#


   #----- Define kappa_0. -----------------------------------------------------------------#
   return(kappa.zero)
   #---------------------------------------------------------------------------------------#
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#    This subroutine will perform the first pass of the objective analysis.                #
#------------------------------------------------------------------------------------------#
oa.weighted.mean = function(dat,dist,kap){

   #----- Find the weights. ---------------------------------------------------------------#
   weight                = exp(- dist * dist / kap)
   #---------------------------------------------------------------------------------------#



   #----- Set the weight to zero for those data that are too far or missing. --------------#
   toofar                  = weight < min.weight
   nodata                  = is.na(dat)
   if (length(dist) != length(toofar) | length(dist) != length(nodata)) browser()
   weight[toofar]          = 0.
   weight[nodata]          = 0.
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     The first step is simply a weighted average of the values.                        #
   #---------------------------------------------------------------------------------------#
   ans = rowSums(dat * weight,na.rm=TRUE) / rowSums(weight,na.rm=TRUE)
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end function oa.weighted.mean
#==========================================================================================#
#==========================================================================================#
