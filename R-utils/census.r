#------------------------------------------------------------------------------------------#
#     This function finds the "census" from the population pop for each element of the     #
# list of possible values categ.                                                           #
#------------------------------------------------------------------------------------------#
census.survey <<- function(pop,categ){
   #---------------------------------------------------------------------------------------#
   #    First check whether the categ vector contains duplicates.  In case it does, crash. #
   #---------------------------------------------------------------------------------------#
   if (length(categ) != length(unique(categ))){
      stop("Vector categ shall not contain duplicates...")
   }#end if
   #---------------------------------------------------------------------------------------#

   #----- Find the number of categories. --------------------------------------------------#
   categ = as.vector(categ)
   ncat  = length(categ)

   #----- Make the population a vector and find its length. -------------------------------#
   popv = as.vector(pop)
   npop = length(popv)

   #----- Make a matrix with the population repeated ncat times. --------------------------#
   popmat = matrix(data=rep(x=popv,times=ncat),ncol=ncat)

   #----- Make a matrix of categories with each category repeated npop times. -------------#
   catmat = matrix(data=rep(x=categ,each=npop),nrow=npop)


   #----- The count is going to be the column sum. ----------------------------------------#
   mycount = colSums(popmat == catmat,na.rm=TRUE)

   #----- The answer is a vector with the count for each category. ------------------------#
   return(mycount)
}#end if
#------------------------------------------------------------------------------------------#
