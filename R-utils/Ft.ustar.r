#==========================================================================================#
#==========================================================================================#
#      Alternative method to determine the u* filter.  It is based on Bonal et al. (2008)  #
# figure, using a simple algorithm that attempts to meet the three idealisations.          #
#                                                                                          #
# 1. Ensure that remaining data are not locally or globally related to u*                  #
# 2. Ensure that most of the signal comes from eddy covariance, not storage.               #
# 3. Retain as much data as possible.                                                      #
#------------------------------------------------------------------------------------------#
Ft.ustar <<- function(ustar,cflxca,cflxst,nighttime,delta=0.01,nmin=10){

   #----- Check that all inputs are given. ------------------------------------------------#
   if (missing(ustar) || missing(cflxca) || missing(cflxst) || missing(nighttime)){
      cat(" Missing ustar:     ",missing(ustar    ),"\n",sep="")
      cat(" Missing cflxca:    ",missing(cflxca   ),"\n",sep="")
      cat(" Missing cflxst:    ",missing(cflxst   ),"\n",sep="")
      cat(" Missing nighttime: ",missing(nighttime),"\n",sep="")
   }#end if 
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    In case storage is always NA, it means storage was not measured. In this case, we  #
   # assume zero storage whenever the CO2 flux is not missing.                             #
   #---------------------------------------------------------------------------------------#
   if ( all(! is.finite(cflxst)) && any(is.finite(cflxca)) ){
      cflxst = ifelse( test = is.finite(cflxca), yes = 0., no = NA_real_ )
   }#end if ( all(! is.finite(cflxst)) && any(is.finite(cflxca)) )
   #---------------------------------------------------------------------------------------#



   #----- Delete missing data. ------------------------------------------------------------#
   keep   = is.finite(ustar) & is.finite(cflxca) & is.finite(cflxst) & nighttime
   keep   = ifelse(is.na(keep),FALSE,keep)
   ustar  = ustar [keep]
   cflxca = cflxca[keep]
   cflxst = cflxst[keep]
   nee    = cflxca + cflxst
   #---------------------------------------------------------------------------------------#


   #----- Split the ustar into quantiles. -------------------------------------------------#
   ustar.breaks    = seq(from=min(ustar),to=max(ustar+delta),by=delta)
   #---------------------------------------------------------------------------------------#



   #----- First guess for u* classes. -----------------------------------------------------#
   ustar.cut = as.numeric(cut(x=ustar,breaks=ustar.breaks,right=FALSE))
   #---------------------------------------------------------------------------------------#



   #------ Look for empty classes and merge them with classes with enough points. ---------#
   ustar.cnt = table(ustar.cut)
   ustar.use = ustar.cnt >= nmin
   ustar.num = as.numeric(names(ustar.cnt))
   ustar.idx = ustar.num[ mapply( FUN      = which.closest
                                , x        = ustar.num
                                , MoreArgs = list(A=ustar.num,mask=ustar.use)) ]
   #---------------------------------------------------------------------------------------#



   #----- Reduce the number of break points and re-categorise u*. -------------------------#
   ustar.breaks = ustar.breaks[c(sort(unique(ustar.idx)),length(ustar.breaks))]
   ustar.cut    = as.numeric(cut(x=ustar,breaks=ustar.breaks,right=FALSE))
   nbins        = length(unique(ustar.cut))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Split the data into the classes.                                                  #
   #---------------------------------------------------------------------------------------#
   sp.nee       = split(x=nee   ,f=ustar.cut)
   sp.cflxca    = split(x=cflxca,f=ustar.cut)
   sp.cflxst    = split(x=cflxst,f=ustar.cut)
   sp.ust       = split(x=ustar ,f=ustar.cut)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find out the first class that average cflxca is statistically significantly       #
   # higher than storage.                                                                  #
   #---------------------------------------------------------------------------------------#
   b      = 0
   iterate = TRUE
   while (iterate){
      b   = b + 1
      ttt     = t.test(x=sp.cflxca[[b]],y=sp.cflxst[[b]],alternative="greater")
      iterate = b < nbins && ttt$p.value %ge% 0.01
   }#end for
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #   Start from the class we just found, and continue until the end.  We look for the    #
   # first occurrence of 4 consecutive classes in which the mean NEE is the same. If it    #
   # gets near the end, we reduce to 3 and 2, and if it fails even at this point, u*       #
   # filter is rejected.                                                                   #
   #---------------------------------------------------------------------------------------#
   b       = b - 1
   iterate = TRUE
   success = FALSE
   while (iterate){
      b = b + 1
      if (b == nbins-1){
         nnn  = list(sp.nee[[b]],sp.nee[[b+1]])
         uuu  = list(sp.ust[[b]],sp.ust[[b+1]])
      }else if (b == nbins - 2){
         nnn  = list(sp.nee[[b]],sp.nee[[b+1]],sp.nee[[b+2]])
         uuu  = list(sp.ust[[b]],sp.ust[[b+1]],sp.ust[[b+2]])
      }else if (b == nbins - 3){
         nnn  = list(sp.nee[[b]],sp.nee[[b+1]],sp.nee[[b+2]],sp.nee[[b+3]])
         uuu  = list(sp.ust[[b]],sp.ust[[b+1]],sp.ust[[b+2]],sp.ust[[b+3]])
      }else{
         nnn  = list(sp.nee[[b]],sp.nee[[b+1]],sp.nee[[b+2]],sp.nee[[b+3]],sp.nee[[b+4]])
         uuu  = list(sp.ust[[b]],sp.ust[[b+1]],sp.ust[[b+2]],sp.ust[[b+3]],sp.ust[[b+4]])
      }#end if
      fff  = F.test(x=nnn)
      p.ft = fff$p.value
      p.lv = fff$levene$p.value
      #------------------------------------------------------------------------------------#



      #---- Check the p.value for a linear fit. -------------------------------------------#
      nnn       = unlist(nnn)
      uuu       = unlist(uuu)
      lm.now    = lm(nnn ~ uuu)
      summ.now  = summary(lm.now)
      p.lm      = summ.now$coefficients[2,4]
      #------------------------------------------------------------------------------------#



      #---- Check whether the p.value is sufficiently large. ------------------------------#
      success = (  ( p.ft %le% p.lm && p.lm %ge% 0.10 && p.ft %ge% 0.10 )
                || ( p.ft %ge% 0.50 && p.lm %ge% 0.50 ) )
      iterate = ( ! success ) && (b < (nbins - 2))
      cat(" u* = ",ustar.breaks[b],";   p.lm    = ",sprintf("%.2f",p.lm)
                                  ,";   p.ft    = ",sprintf("%.2f",p.ft)
                                  ,";   success = ",success
                                  ,"\n",sep="")
      #------------------------------------------------------------------------------------#
   }#end while
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #     Check whether this was successful.                                                #
   #---------------------------------------------------------------------------------------#
   if (success){
      ustar.out = ustar.breaks[b]
   }else{
      ustar.out = ustar.breaks[nbins]
   }#end if
   #---------------------------------------------------------------------------------------#
   return(ustar.out)
}#end function Ft.ustar
#==========================================================================================#
#==========================================================================================#
