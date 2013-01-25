#------------------------------------------------------------------------------------------#
#    This function splits the number into number of rows and columns that will make the    #
# prettiest (prettiest here means the closest to a golden ratio rectangle).                #
#------------------------------------------------------------------------------------------#
pretty.box = function(n,horizontal=TRUE,angle.crit=atan2((1.+sqrt(5))/2.,1)*180./pi){


   #---------------------------------------------------------------------------------------#
   #      Check the size of n.  If it has a single value, then we must find the number of  #
   # rows and columns.  If it has lenght 2, then the first value becomes the number of     #
   # columns and the second becomes the number of rows.                                    #
   #---------------------------------------------------------------------------------------#
   n.n = length(n)
   if (n.n == 0){
      #------------------------------------------------------------------------------------#
      #      Empty! Make a single box.                                                     #
      #------------------------------------------------------------------------------------#
      nbrow  = 1
      nbcol  = 1
      nbox   = 1
      mat    = matrix(1,ncol=1,nrow=1)
      xangle = 45
      yangle = 45
      #------------------------------------------------------------------------------------#

   }else if (n.n == 1){
      #------------------------------------------------------------------------------------#
      #      Single site.  We must test whether we can split in nows and columns that      #
      # looks nice by the closest to square root.  If the number is prime, we must add one #
      # so it can be split.   We keep adding empty boxes until we get a number that has a  #
      # good ratio.                                                                        #
      #------------------------------------------------------------------------------------#
      iterate = TRUE
      n       = max(n,1)
      nbox    = n
      while (iterate){

         nbrow.pot = seq(from=1,to=floor(nbox),by = 1)
         nbcol.pot = nbox / nbrow.pot

         #----- Select only the cases that make a 2-D grid with the right orientation. ----#
         if (horizontal){
            ok  = nbcol.pot == as.integer(nbcol.pot) & nbcol.pot >= nbrow.pot
         }else{
            ok  = nbcol.pot == as.integer(nbcol.pot) & nbcol.pot <= nbrow.pot
         }#end if
         #---------------------------------------------------------------------------------#


         #----- Discard bad options. ------------------------------------------------------#
         nbrow.pot = nbrow.pot[ok]
         nbcol.pot = as.integer(nbcol.pot[ok])
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Choose the one that is the closest to the square.                          #
         #---------------------------------------------------------------------------------#
         nuse = which.min(abs(nbrow.pot-nbcol.pot))
         nbrow = nbrow.pot[nuse]
         nbcol = nbcol.pot[nuse]
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Make a matrix for the layout.                                              #
         #---------------------------------------------------------------------------------#
         mat          = matrix(sequence(nbrow*nbcol),ncol=nbcol,nrow=nbrow,byrow=TRUE)
         mat[mat > n] = 0
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Check whether the angle is acceptable...                                   #
         #---------------------------------------------------------------------------------#
         xangle = atan2(nbrow,nbcol) * 180. / pi
         yangle = atan2(nbcol,nbrow) * 180. / pi
         iterate = nbox > 2 && (xangle > angle.crit | yangle > angle.crit)
         if (iterate) nbox = nbox + 1
      }#end while
      #------------------------------------------------------------------------------------#

   }else if (n.n == 2){
      #------------------------------------------------------------------------------------#
      #      Two values were given.  No iteration needed.                                  #
      #------------------------------------------------------------------------------------#
      nbcol  = max(n[1],1)
      nbrow  = max(n[2],1)
      nbox   = nbcol*nbrow
      mat    = matrix(sequence(nbrow*nbcol),ncol=nbcol,nrow=nbrow,byrow=TRUE)
      xangle = atan2(nbrow,nbcol) * 180. / pi
      yangle = atan2(nbcol,nbrow) * 180. / pi
      #------------------------------------------------------------------------------------#

   }else{
      #------------------------------------------------------------------------------------#
      #      3-D boxes are not available...                                                #
      #------------------------------------------------------------------------------------#
      cat (" n           = (",paste(n,sep="; "),")","\n")
      cat (" length of n =  ",length(n)            ,"\n")
      stop(" Invalid n! It must have length between 0 and 2")
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Return a list with the information.                                               #
   #---------------------------------------------------------------------------------------#
   if (horizontal){
      ans = list(nrow=nbrow,ncol=nbcol,nbox=nbox,angle=yangle,mat=mat)
   }else{
      ans = list(nrow=nbrow,ncol=nbcol,nbox=nbox,angle=xangle,mat=mat)
   }#end if
   return(ans)
}#end function
#------------------------------------------------------------------------------------------#
