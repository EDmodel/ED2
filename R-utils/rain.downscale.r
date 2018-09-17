#==========================================================================================#
#==========================================================================================#
#     This function downscales rainfall by taking into account the precipitation fraction  #
# (roughly speaking the average number of times when precipitation is recorded.            #
#------------------------------------------------------------------------------------------#
rain.downscale <<- function(lon,when.in,prate.in,mu.in.lut,nsub,rain.min=0.2/3600
                           ,mask.low=TRUE){

   #----- Hour blocks. --------------------------------------------------------------------#
   huniq  = sequence(24)-1
   nhidx  = ncol(mu.in.lut)
   bwidth = length(huniq) / nhidx
   hidx   = (sequence(nhidx) - 1) * day.hr / nhidx
   hblock = ( ceiling(huniq / bwidth) %% nhidx) + 1
   #---------------------------------------------------------------------------------------#


   #----- Initialise gridded look-up tables. ----------------------------------------------#
   g.rho.lut  = matrix(data=NA_real_,nrow=12,ncol=nhidx,dimnames=list(month.abb,hidx))
   g.rbar.lut = g.rho.lut
   g.mu.lut   = g.rho.lut
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     UTC correction: We use the geographic position as opposed to legal time zone as   #
   # legal time zones often do not match the solar one.  We keep it simple and adjust in   #
   # hourly shifts, as we rarely have sub-hourly data to make a more accurate correction.  #
   #---------------------------------------------------------------------------------------#
   utcoff.in = round(lon/15)
   #---------------------------------------------------------------------------------------#



   #------ Find precipitation density. ----------------------------------------------------#
   ndat.in         = length(when.in)
   month.in        = nummonths(when.in)
   hour.in         = hours(when.in)
   hidx.in         = hblock[match((hour.in+utcoff.in)%%24,huniq)]
   tmp             = aggregate( x    = prate.in
                              , by   = list(month.in,hidx.in)
                              , FUN  = mean.above
                              , xlwr = rain.min
                              , xnot = 0.
                              )#end aggregate
   ntmp            = nrow(tmp)
   idx             = cbind(tmp[[1]],tmp[[2]])
   g.rho.lut[idx]  = tmp[[3]]
   #---------------------------------------------------------------------------------------#



   #------ Find the expected rainfall (including zeroes). ---------------------------------#
   tmp              = aggregate( x     = prate.in
                               , by    = list(month.in,hidx.in)
                               , FUN   = mean
                               , na.rm = TRUE
                               )#end aggregate
   ntmp             = nrow(tmp)
   idx              = cbind(tmp[[1]],tmp[[2]])
   g.rbar.lut[idx]  = tmp[[3]]
   #---------------------------------------------------------------------------------------#


   #----- Find the ratio following Ryan's thesis. -----------------------------------------#
   g.mu.lut         = ifelse(test=g.rho.lut %>% 0,yes=g.rbar.lut/g.rho.lut,no=0.)
   #---------------------------------------------------------------------------------------#


   #----- Copy the gridded and the climate mu to vectors that match output. ---------------#
   idx    = cbind(month.in,hidx.in)
   mu.in  = mu.in.lut [idx]
   g.rbar = g.rbar.lut[idx]
   g.mu   = g.mu.lut  [idx]
   sg.mu  = ifelse(test = g.mu %>% 0., yes = mu.in / g.mu, no = 0.)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Build CDF matrices. The first column can range from 0 to 1.  Additional columns    #
   # are restricted to 1 - sum of each column.  After all columns are filled, we shuffle   #
   # the order so rainfall peaks can happen at any substep.                                #
   #---------------------------------------------------------------------------------------#
   eps.off   = sqrt(.Machine$double.eps)
   CDF.MAT   = matrix( data = NA_real_, nrow = ndat.in, ncol = nsub)
   for (n in sequence(nsub-1)){
      cdf.sum = pmin(1.-eps.off,pmax(eps.off,rowSums(CDF.MAT,na.rm=TRUE)))
      cdf.rnd = runif(n=ndat.in,min=eps.off,max=1-cdf.sum)
      if (! (all(is.finite(cdf.sum)) && all(is.finite(cdf.rnd)))){
         cat0(" NA values were spotted.")
         browser()
      }#End if 
      CDF.MAT[,n] = cdf.rnd
   }#end for
   CDF.MAT[,nsub] = 1 - rowSums(CDF.MAT,na.rm=TRUE)
   CDF.MAT        = t(apply(X=CDF.MAT,MARGIN=1,FUN=lit.sample,size=nsub,replace=FALSE))
   #---------------------------------------------------------------------------------------#


   #----- Make the downscaled precipitation matrix. ---------------------------------------#
   SG.MU     = matrix( data = sg.mu   , nrow = ndat.in, ncol = nsub)
   G.RBAR    = matrix( data = g.rbar  , nrow = ndat.in, ncol = nsub)
   PRATE.IN  = matrix( data = prate.in, nrow = ndat.in, ncol = nsub)
   WGT.OUT   = ifelse( test = SG.MU %>% 0. & PRATE.IN %>% 0.
                     , yes  = PRATE.IN / SG.MU * log(1.0/(1.0 - CDF.MAT))
                     , no   = 1 / nsub
                     )#end ifelse
   WGT.SUM   = matrix( data = rep(rowMeans(WGT.OUT),nsub)
                     , nrow = ndat.in
                     , ncol = nsub
                     )
   PRATE.OUT = PRATE.IN * WGT.OUT / WGT.SUM
   #---------------------------------------------------------------------------------------#



   #----- Expand the vector with rainfall. ------------------------------------------------#
   prate.out  = c(t(PRATE.OUT))
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Mask low precipitation events, but only when at least one rainfall event exceeds #
   # the minimum.                                                                          #
   #---------------------------------------------------------------------------------------#
   if (mask.low){
      ptotal.out = sum(prate.out)
      prate.mask  = ifelse(prate.out %>=% rain.min, prate.out, 0.)
      ptotal.mask = sum(prate.mask)
      if (ptotal.mask %>=% rain.min){
          prate.out = prate.mask * ptotal.out / ptotal.mask
      }#end if (ptotal.mask %>=% rain.min)
   }#end if (mask.low)
   #---------------------------------------------------------------------------------------#


   #------ Make sure the average output rainfall matches the average input rainfall. ------#
   prate.out.bar = mean(prate.out,na.rm=TRUE)
   if (prate.out.bar %>=% rain.min){
      prate.in.bar = mean(prate.in,na.rm=TRUE)
      prate.out    = prate.out * prate.in.bar / prate.out.bar
   }#end if (prate.out.bar %>=% rain.min)
   #---------------------------------------------------------------------------------------#

   return(prate.out)
   #---------------------------------------------------------------------------------------#
}#end rain.downscale
#==========================================================================================#
#==========================================================================================#
