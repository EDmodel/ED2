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


   #----- Make a matrix with offset of one in case the user appends a legend. -------------#
   sel           = mat == 0
   mat.off       = mat + 1
   mat.off[sel]  = 0
   mat.off2      = mat + 2
   mat.off2[sel] = 0
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Create the margins for each box assuming that axes should be plotted only on the  #
   # bottom and left panels.                                                               #
   #---------------------------------------------------------------------------------------#
   n.seq = sequence(n)
   if (horizontal){
      left   = ( n.seq %% nbcol ) == 1 | nbcol == 1
      right  = ( n.seq %% nbcol ) == 0
      top    = n.seq <=   nbcol
      bottom = n.seq >  ( nbrow - 1 ) * nbcol
   }else{
      left   = n.seq <= nbrow
      right  = n.seq > ( nbcol - 1 ) * nbrow
      top    = ( n.seq %% nbrow ) == 1 | nbrow == 1
      bottom = ( n.seq %% nbrow ) == 0
   }#end if
   centre = ! ( left | right )
   middle = ! ( bottom | top )
   if (nbcol == 1){
      mar.left  = rep(5.1,times=n)
      mar.right = rep(1.1,times=n)
   }else if (nbcol == 2){
      mar.left  = 2.1 + 2 * left
      mar.right = 0.1 + 2 * right
   }else{
      mar.left  = 2.1 + 2 * left   - 2 * right
      mar.right = 0.1 + 2 * centre + 4 * right
   }#end if
   if (nbrow == 1){
      mar.bottom = rep(5.1,times=n)
      mar.top    = rep(4.1,times=n)
   }else if (nbrow == 2){
      mar.bottom = 2.1 + 2 * bottom
      mar.top    = 2.1 + 2 * top
   }else{
      mar.bottom = 2.1 + 2 * bottom - 2 * top
      mar.top    = 1.1 + 2 * middle + 4 * top
   }#end if
   mar                = cbind(mar.bottom,mar.left,mar.top,mar.right)
   dimnames(mar)[[2]] = c("bottom","left","top","right")

   #----- Margins for when all boxes get axis plots. --------------------------------------#
   if (nbcol == 1){
      mar0.left  = 5.1
      mar0.right = 1.1
   }else{
      mar0.left  = 4.1
      mar0.right = 0.1
   }#end if
   if (nbrow == 1){
      mar0.bottom = 5.1
      mar0.top    = 4.1
   }else{
      mar0.bottom = 4.1
      mar0.top    = 2.1
   }#end if
   mar0        = c(mar0.bottom,mar0.left,mar0.top,mar0.right)
   names(mar0) = c("bottom","left","top","right")
   #---------------------------------------------------------------------------------------#




   #----- Key margins for when all boxes get axis plots. ----------------------------------#
   key.left  = 0.6
   key.right = 4.1
   if (nbrow == 1){
      key.bottom = 5.1
      key.top    = 4.1
   }else{
      key.bottom = 4.1
      key.top    = 2.1
   }#end if
   mar.key        = c(key.bottom,key.left,key.top,key.right)
   names(mar.key) = c("bottom","left","top","right")
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Return a list with the information.                                               #
   #---------------------------------------------------------------------------------------#
   ans = list( nrow       = nbrow
             , ncol       = nbcol
             , nbox       = nbox
             , angle      = horizontal * yangle + (1 - horizontal) * xangle
             , mat        = mat
             , mat.off    = mat.off
             , mat.off2   = mat.off2
             , bottom     = bottom
             , left       = left
             , top        = top
             , right      = right
             , mar        = mar
             , mar0       = mar0
             , mar.key    = mar.key
             )#end list
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function
#------------------------------------------------------------------------------------------#
