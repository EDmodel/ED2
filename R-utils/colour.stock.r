#==========================================================================================#
#==========================================================================================#
#     This function creates a vector with colours for plots that are mostly CB friendly.   #
#  The quality decays as the number of lines increases.                                    #
#------------------------------------------------------------------------------------------#
colour.stock <<- function(n){
   #----- Build the colour stock.  Add pale and dark versions of the default colours. -----#
   cstock = c("#3B24B3","#A3CC52","#E65C17"
             ,"#990F0F","#306614","#2996CC"
             ,"#B49ED2","#F5C858","#00F3FB"
             )#end c
   pale   = round(0.5 * (col2rgb(cstock) + 255))
   pale   = RGB(R=pale[1,],G=pale[2,],B=pale[3,])
   dark   = round(0.5 * col2rgb(cstock))
   dark   = RGB(R=dark[1,],G=dark[2,],B=dark[3,])
   cstock = c(cstock,pale,dark)

   #---------------------------------------------------------------------------------------#
   #    Make sure does not exceed the maximum number of colours.                           #
   #---------------------------------------------------------------------------------------#
   if (n > length(cstock)){
      cat("-----------------------------------------","\n",sep="")
      cat(" Provided n      = ",n                    ,"\n",sep="")
      cat(" Maximum allowed = ",length(cstock)       ,"\n",sep="")
      cat("-----------------------------------------","\n",sep="")
      stop("Either reduce n or add more colours to the colour stock")
   }else{
      ans = cstock[sequence(n)]
      return(ans)
   }#end if
   #---------------------------------------------------------------------------------------#
}#end colour.stock
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     This function creates a vector with pch for plots to make it easier for CB folks.    #
#  The quality decays as the number of lines increases.                                    #
#------------------------------------------------------------------------------------------#
pch.stock <<- function(n){
   #----- Build the colour stock.  Add pale and dark versions of the default colours. -----#
   pstock = c(16, 4,13,17, 6, 7, 0, 5, 2
             , 3, 1,14,10,18, 9, 8,12,15
             ,11, 6,16, 5,13, 4,17, 3, 0
             )#end c

   #---------------------------------------------------------------------------------------#
   #    Make sure does not exceed the maximum number of colours.                           #
   #---------------------------------------------------------------------------------------#
   if (n > length(pstock)){
      cat("-----------------------------------------","\n",sep="")
      cat(" Provided n      = ",n                    ,"\n",sep="")
      cat(" Maximum allowed = ",length(cstock)       ,"\n",sep="")
      cat("-----------------------------------------","\n",sep="")
      stop("Either reduce n or add more pch types to the pch stock")
   }else{
      ans = pstock[sequence(n)]
      return(ans)
   }#end if
   #---------------------------------------------------------------------------------------#
}#end pch.stock
#==========================================================================================#
#==========================================================================================#
