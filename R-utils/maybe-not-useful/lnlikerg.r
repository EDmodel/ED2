#------------------------------------------------------------------------------------------#
#     This function will build the support, or log of the likelihood, based on rain gauge  #
# precipitation statistics.                                                                #
#                                                                                          #
#    Support is computed as:                                                               #
#                                                                                          #
#  S = k + n * p * {ln(PP) - ln(SS) - [s^2+(m-MM)^2]/[2*SS^2]} + n * (1-p) * ln(1-PP)      #
#                                                                                          #
# where:                                                                                   #
#                                                                                          #
# k  - arbitrary constant (it will be assigned later, here it will be 0.                   #
# n  - number of observations.                                                             #
# p  - fraction of observations with precipitation greater than 0.                         #
# m  - mean log of non-zero precipitation, based on observations.                          #
# s  - standard deviation of the log of non-zero precipitation, based on observations.     #
# NN - number of points used by the model statistics.                                      #
# PP - model probability of precipitation                                                  #
# MM - model mean of log of non-zero precipitation.                                        #
# SS - model standard deviation of log of non-zero precipitation.                          #
#------------------------------------------------------------------------------------------#
lnlikerg = function(brams,obs){

   #----- For the likelihood, the mean and standard deviation must be in log scale. -------#
   n  = obs$rg$dccnt
   p  = obs$rg$dprob
   m  = log(obs$rg$dmean)
   s  = log(obs$rg$dsdev)
   NN = brams$rg$dccnt
   PP = brams$rg$dprob
   MM = log(brams$rg$dmean)
   SS = log(brams$rg$dsdev)

   #---------------------------------------------------------------------------------------#
   #     Precipitation mean and standard deviation are both NA when it never rains at a    #
   # given point at a particular time.  In order to avoid not accounting for the obser-    #
   # vation, we impose both m and s to be zero when p is 0. The "rainfall"  term should be #
   # zero anyway because everything is multiplied by p.                                    #
   #---------------------------------------------------------------------------------------#
   sel = p == 0.
   m[sel]  =  0
   s[sel]  =  0

   #---------------------------------------------------------------------------------------#
   #     EXTREMELY CONTROVERSIAL POINT.  I DON'T THINK THIS IS REALLY RIGHT, BUT THAT WAS  #
   # THE BEST I COULD THINK...                                                             #
   #   There is nothing that prevents PP to be zero or one.  However, when PP is 0, the    #
   # support would be 0 if p is also 0, and -Infinity otherwise.  Likewise, when PP is 1,  #
   # support would be 0 if p is also 1, and -Infinity otherwise.  This is bad... My        #
   # argument is that if PP is 0, it doesn't necessarily mean that PP is really 0, because #
   # this was computed based on discrete data, and therefore it could be anything between  #
   # 0 and 0.5/NN.  Therefore I just make it 0.25/NN, and assume that the mean and         #
   # standard deviation are the average value for the entire domain.  Likewise, if PP is   #
   # 1, I assume that the PP is actually (NN-0.25)/NN, but keep the same distribution that #
   # was found.  The advantage of doing this is that it penalises the model by a lot if    #
   # NN is large and by a little when NN is small.                                         #
   #---------------------------------------------------------------------------------------#
   MMmean  = matrix(rep(apply(X=MM,MARGIN=2,FUN=mean,na.rm=TRUE),each=nrow(MM))
                   ,nrow=nrow(MM),ncol=ncol(MM))
   SSmean  = matrix(rep(apply(X=SS,MARGIN=2,FUN=mean,na.rm=TRUE),each=nrow(SS))
                   ,nrow=nrow(SS),ncol=ncol(SS))
   sel     = PP == 0
   MM[sel] = MMmean[sel]
   SS[sel] = SSmean[sel]
   PP[sel] = 0.25 / NN[sel]

   sel     = PP == 1 / NN
   SS[sel] = SSmean[sel]

   sel     = PP == 1
   PP[sel] = (NN[sel]-0.25)/NN[sel]

   #---------------------------------------------------------------------------------------#
   #   Now find the likelihood for each point independently.                               #
   #---------------------------------------------------------------------------------------#
   lnlike = n * (        p  * (log(PP) - log(SS) - (s^2 + (m - MM)^2)/(2*SS^2))
                + (1.0 - p) * log(1.0 - PP) )
   lnlike = lnlike - max(lnlike,na.rm=TRUE)

   return(lnlike)
}#end function
#------------------------------------------------------------------------------------------#

