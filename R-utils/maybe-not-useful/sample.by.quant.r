#==========================================================================================#
#==========================================================================================#
#       This function creates a sample of the original data using the quantiles to make    #
# the data biased towards lower or higher quantiles.                                       #
#                                                                                          #
# Input variables:                                                                         #
#  ~ x      -- the vector to be sampled                                                    #
#  ~ v      -- values associated with the vector x (from which we draw the distribution)   #
#  ~ a      -- the factor to shift the probability of sampling any element.  Check method. #
#              For all methods, there are two special cases:                               #
#              a = NA - no sampling                                                        #
#              a = 0  - uniform sampling with replacement.                                 #
#  ~ method -- (case insensitive and partial match).  From v we determine the quantiles qq #
#              of each value, then we define the PDF according to the chosen method:       #
#                                                                                          #
#              -- "power": p (q) = k * abs(a) ^ (sign(a)*( 1 - 2 * q)); 0 <= q <= 1        #
#                          p (q) = 0  otherwise                                            #
#                                                                                          #
#              -- "logit": p (q) = k * inv.logit((a-sign(a))*(q-0.5-0.01*a)) 0 <= q <= 1   #
#                          p (q) = 0  otherwise                                            #
#                                                                                          #
#              -- "bias" : for each year, we sample one year giving more probability for   #
#                          quantiles that are lower or higher than the previous choice:    #
#                          a = 1 or -1 - equal chances.                                    #
#                          a < -1      - lower quantiles will be a^2 more likely to be     #
#                                        chosen.                                           #
#                          a >  1      - higher quantiles will be a^2 more likely to be    #
#                                        chosen.                                           #
#                                                                                          #
#              -- "skew" : the PDF is skewed towards lower or higher quantiles.            #
#                          using skew-normal probability distribution.                     #
#                          a  - the location is shifted by "a" units of scale              #
#                               parameter (sign tells the direction), unless a is NA or    #
#                               zero.                                                      #
#   k doesn't need to be determine because sample scales the PDF.                          #
#------------------------------------------------------------------------------------------#
sample.by.quant <<- function(x,v,method="skew",a=1){

   #----- First we check whether a is valid. ----------------------------------------------#
   if (length(a) != 1){
      stop (paste(" Factor a must be a scalar, and you entered one with length "
                 ,length(a),"...",sep=""))
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- First we check whether x or v are valid. ----------------------------------------#
   if (missing(x) || missing(v)){
      cat (" - x is missing: ",missing(x),"\n")
      cat (" - v is missing: ",missing(v),"\n")
      stop(" Either x or v (or both) is missing...")
   }else if (length(x) < length(v)){
      v.names   = as.numeric(names(v))
      ix        = match(x,v.names)
      va        = min(v.names)
      vz        = max(v.names)
      xa        = min(x)
      xz        = max(x)
      nv        = length(v)
      nx        = length(x)
      ncyc      = xz - xa + 1
      nfullcyc  = floor((xa-va)/ncyc)
      firstfull = xa - nfullcyc * ncyc
      ia        = match(xz - (firstfull - va),x) + 1
      idx       = c(seq(from=ia,to=nx,by=1),rep(1:nx,times=1+ceiling(nv/nx)))[1:nv]
      v         = v [v.names %in% x][idx]
      ix        = ix[idx]
      names(v)  = v.names
      x         = v.names
   }else if (length(x) > length(v)){
      keep     = intersect(x,as.numeric(names(v)))
      ix       = which(x %in% keep)
      x        = x[x                    %in% keep]
      v        = v[as.numeric(names(v)) %in% keep]
   }else{
      ix       = sequence(length(x))
   }#end if
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Simplify the name of the method.                                                  #
   #---------------------------------------------------------------------------------------#
   if (length(method) != 1){
      stop(" Variable method must be scalar!")
   }else{
      use = substring(tolower(method),1,1)
   }#end if
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #      Define some constants.                                                           #
   #---------------------------------------------------------------------------------------#
   nx = length(x)
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Find the quantiles of each value of v using the empirical distribution function.  #
   #---------------------------------------------------------------------------------------#
   v.ecdf = ecdf  (v)
   quant  = v.ecdf(v)
   quant  = (quant - min(quant)) / (max(quant)-min(quant))
   dquant = mean(diff(quant))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     If a is NA, we simply copy the values of X. If it is 0, then we just sample the   #
   # values with uniform distribution.                                                     #
   #---------------------------------------------------------------------------------------#
   if (is.na(a)){
      ans  = list ( orig   = x
                  , sample = x
                  , expect = mean(v)
                  , s.idx  = ix
                  , quant  = quant
                  , prob   = rep(1/nx,times=nx)
                  )#end list
   }else if (a == 0){
      s         = sample(x=x,size=nx,replace=TRUE)
      cnt       = table(s)
      prob      = rep(0,times=nx)
      idx       = match(as.numeric(names(cnt)),x)
      prob[idx] = cnt/nx
      ans       = list ( orig   = x
                       , sample = x
                       , expect = mean(v)
                       , s.idx  = ix[match(s,x)]
                       , quant  = quant
                       , prob   = prob
                       )#end list
   }else if (use == "p"){

      prob = abs(a) ^ (sign(a)*(1-2*quant))

      s    = sample(x=x,size=nx,replace=TRUE,prob=prob)
      ans  = list( orig   = x
                 , sample = s
                 , expect = mean(v)
                 , s.idx  = ix[match(s,x)]
                 , quant  = quant
                 , prob   = prob
                 )#end list
   }else if (use == "l"){
      prob = inv.logit((a-sign(a))*(quant-0.5-0.01*a))
      s    = sample(x=x,size=nx,replace=TRUE,prob=prob)
      ans  = list( orig   = x
                 , sample = s
                 , expect = mean(v)
                 , s.idx  = ix[match(s,x)]
                 , quant  = quant
                 , prob   = prob
                 )#end list
   }else if (use == "b"){
      s    = rep(NA,times=nx)
      prob = abs(a) ^ (sign(a)*quant > sign(a)*0.5)
      s[1] = sample(x=x,size=1,replace=TRUE,prob=prob)
      for (i in 2:nx){
         prob = abs(a) ^ (sign(a)*(2*(quant > quant[i-1])-1))
         s[i] = sample(x=x,size=1,replace=TRUE,prob=prob)
      }#end for

      cnt       = table(s)
      prob      = rep(0,times=nx)
      idx       = match(as.numeric(names(cnt)),x)
      prob[idx] = cnt/nx
      ans       = list( orig   = x
                      , sample = s
                      , expect = mean(v)
                      , s.idx  = ix[match(s,x)]
                      , quant  = quant
                      , prob   = prob
                      )#end list
   }else if (use == "s"){
      #----- Find the skew normal statistics. ---------------------------------------------#
      v.loc    = sn.location(x=v,na.rm=TRUE)
      v.scl    = sn.scale   (x=v,na.rm=TRUE)
      v.shp    = sn.shape   (x=v,na.rm=TRUE)
      #------------------------------------------------------------------------------------#



      #----- Modify the shape and location based on "a". ----------------------------------#
      s.loc   = v.loc + (round(a)-sign(a)) * 0.1 * v.scl
      s.scl   = v.scl
      s.shp   = v.shp
      d.shp   = round.log(abs(s.shp))
      q.100   = qsn(p=0.10,location=s.loc,scale=s.scl,shape=s.shp)
      q.900   = qsn(p=0.90,location=s.loc,scale=s.scl,shape=s.shp)
      v.min   = min(v)
      v.max   = max(v)
      #it = 0
      #while ( it < 50 && (q.100 < v.min || q.900 > v.max)){
      #   it = it+1
      #   s.scl = 0.99 * s.scl
      #   if (q.100 < v.min){
      #      s.shp = s.shp + 0.1 * d.shp
      #   }#end if
      #   if (q.900 > v.max){
      #      s.shp = s.shp - 0.1 * d.shp
      #   }#end if
      #   q.100 = qsn(p=0.10,location=s.loc,scale=s.scl,shape=s.shp)
      #   q.900 = qsn(p=0.90,location=s.loc,scale=s.scl,shape=s.shp)
      #}#end while
      #------------------------------------------------------------------------------------#


      #----- Use the new statistics to find the new distribution. -------------------------#
      v.goal    = rsn(n=nx,location=s.loc,scale=s.scl,shape=s.shp)
      s.idx     = mapply(FUN=which.closest,x=v.goal,MoreArgs=list(A=v))
      s         = x[s.idx]
      cnt       = table(s)
      prob      = rep(0,times=nx)
      idx       = match(as.numeric(names(cnt)),x)
      prob[idx] = cnt/nx
      ans       = list( orig   = x
                      , sample = s
                      , expect = mean(v)
                      , s.idx  = ix[match(s,x)]
                      , quant  = quant
                      , prob   = prob
                      )#end list
      #------------------------------------------------------------------------------------#

   }else{
       stop (paste(" Method ",method," is not recognised...",sep=""))
   }#end if
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end function
#------------------------------------------------------------------------------------------#
