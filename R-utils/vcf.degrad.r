#----- Polygon that defines the degradation classes. --------------------------------------#
vcf.lines     = list()
vcf.lines$VC0 = list( shrub = c(0.25,0.00,0.00,0.25)
                    , tree  = c(0.75,1.00,0.75,0.75)
                    , bare  = c(0.00,0.00,0.25,0.00)
                    )#end list
vcf.lines$VC1 = list( shrub = c(0.50,0.25,0.05,0.30,0.50)
                    , tree  = c(0.50,0.75,0.75,0.50,0.50)
                    , bare  = c(0.00,0.00,0.20,0.20,0.00)
                    )#end list
vcf.lines$VC2 = list( shrub = c(0.30,0.05,0.00,0.00,0.30)
                    , tree  = c(0.50,0.75,0.75,0.50,0.50)
                    , bare  = c(0.20,0.20,0.25,0.50,0.20)
                    )#end list
vcf.lines$VC3 = list( shrub = c(0.75,0.50,0.30,0.55,0.75)
                    , tree  = c(0.25,0.50,0.50,0.25,0.25)
                    , bare  = c(0.00,0.00,0.20,0.20,0.00)
                    )#end list
vcf.lines$VC4 = list( shrub = c(0.55,0.30,0.00,0.25,0.55)
                    , tree  = c(0.25,0.50,0.50,0.25,0.25)
                    , bare  = c(0.20,0.20,0.50,0.50,0.20)
                    )#end list
vcf.lines$VC5 = list( shrub = c(1.00,0.75,0.25,0.50,1.00)
                    , tree  = c(0.00,0.25,0.25,0.00,0.00)
                    , bare  = c(0.00,0.00,0.50,0.50,0.00)
                    )#end list
vcf.lines$VC6 = list( shrub = c(0.50,0.00,0.00,0.50)
                    , tree  = c(0.00,0.50,0.00,0.00)
                    , bare  = c(0.50,0.50,1.00,0.50)
                    )#end list
#----- List of centroids of each class. ---------------------------------------------------#
o253         = .25/3
o503         = .50/3
d292         = .375 - o253
a667         = .50+o503
a833         = .75+o253
vcf.centroid = data.frame( shrub = c(o253,0.275, o253,0.525,0.275,0.625,o503)
                         , tree  = c(a833,0.625,0.625,0.375,0.375,0.125,o503)
                         , bare  = c(o253,0.100, d292,0.100,0.350,0.250,a667)
                         )#end data.frame
rownames(vcf.centroid) = names(vcf.lines)
#----- Make variables global. -------------------------------------------------------------#
vcf.lines    <<- vcf.lines
vcf.centroid <<- vcf.centroid
n.vcf        <<- nrow(vcf.centroid)
vcf.number   <<- sequence(n.vcf)-1
vcf.names    <<- names(vcf.lines)
#------------------------------------------------------------------------------------------#




#==========================================================================================#
#==========================================================================================#
#     This function attributes a degradation class based on vegetation cover.              #
#------------------------------------------------------------------------------------------#
vcf.degrad <<- function(shrub,tree,bare){
   #----- It is fine to provide only two classes because the third is redundant. ----------#
   miss.shrub = missing(shrub)
   miss.tree  = missing(tree )
   miss.bare  = missing(bare )
   if ((miss.shrub + miss.tree + miss.bare) > 1){
      cat0("   --> Shrub is missing: ",miss.shrub)
      cat0("   --> Tree  is missing: ",miss.tree )
      cat0("   --> Bare  is missing: ",miss.bare )
      stop(" You must provide at least two elements to build vegetation cover!")
   }else if (miss.shrub){
      shrub = 1. - tree - bare
   }else if (miss.tree ){
      tree  = 1. - shrub - bare
   }else if (miss.bare ){
      bare  = 1. - shrub - tree
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Make sure that all non-NA elements add up to one and are bounded. ---------------#
   total    = shrub + tree + bare
   good     = is.na(total) | ( total %wr% (1 + c(-1,+1)*sqrt(.Machine$double.eps)) )
   ok.shrub = is.na(shrub) | ( shrub %wr% c(0,1) )
   ok.tree  = is.na(tree ) | ( tree  %wr% c(0,1) )
   ok.bare  = is.na(bare ) | ( bare  %wr% c(0,1) )
   if (any(! good)){
      stop(" shrub + tree + bare must add up to 1.")
   }else if (! all(c(ok.shrub,ok.tree,ok.bare))){
      stop(" shrub, tree, and bare must be all between 0 and 1.")
   }#end total
   #---------------------------------------------------------------------------------------#


   #----- Create vector by classes. -------------------------------------------------------#
   vc0 = tree %>=% 0.75
   vc1 = tree %>=% 0.50 & tree  %<%  0.75 & bare  %<%  0.20
   vc2 = tree %>=% 0.50 & tree  %<%  0.75 & bare  %>=% 0.20
   vc3 = tree %>=% 0.25 & tree  %<%  0.50 & bare  %<%  0.20
   vc4 = tree %>=% 0.25 & tree  %<%  0.50 & bare  %>=% 0.20 & bare %<% 0.50
   vc5 = tree %<%  0.25 & bare  %<%  0.50
   vc6 = bare %>=% 0.50
   #---------------------------------------------------------------------------------------#


   #----- Attribute classes. --------------------------------------------------------------#
   vcf      = NA + tree
   vcf[vc0] = 0
   vcf[vc1] = 1
   vcf[vc2] = 2
   vcf[vc3] = 3
   vcf[vc4] = 4
   vcf[vc5] = 5
   vcf[vc6] = 6
   #---------------------------------------------------------------------------------------#


   return(vcf)
   #---------------------------------------------------------------------------------------#
}#end vcover.degrad
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#       Convert the quality data into count.                                               #
#------------------------------------------------------------------------------------------#
vcf.qaqc.pass <<- function(x,nblock=10000){


   #---------------------------------------------------------------------------------------#
   #     Find blocks to process.                                                           #
   #---------------------------------------------------------------------------------------#
   nx = length(x)
   na = seq(from=1,to=nx,by=nblock)
   nz = pmin(na + nblock - 1,nx)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Initialise pass with 8 (everything failed), we will split large arrays into for   #
   # computation efficiency.                                                               #
   #---------------------------------------------------------------------------------------#
   pass = rep(8,times=nx)
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #      Go through each block.                                                           #
   #---------------------------------------------------------------------------------------#
   for (n in seq_along(na)){
      idx       = seq(from=na[n],to=nz[n],by=1)
      nnow      = length(idx)
      pass[idx] = colSums(matrix(data=as.numeric(intToBits(x[idx])),nrow=32,ncol=nnow))
   }#end for 
   #---------------------------------------------------------------------------------------#


   #----- Return answer. ------------------------------------------------------------------#
   return(pass)
   #---------------------------------------------------------------------------------------#
}#end function
#==========================================================================================#
#==========================================================================================#
