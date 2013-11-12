#==========================================================================================#
#==========================================================================================#
#     This function plots global titles for multiple panel plots. Arguments are the same   #
# as function "title", plus optional argument telling which relative position should be    #used #
# for xlab, ylab, main, and sub.                                                           #
#                                                                                          #
# IMPORTANT:  This should be the last thing to be called before you close the device.  If  #
#             you use it before, it will probably mess up the plots.                       #
#------------------------------------------------------------------------------------------#
gtitle <<- function( main      = NULL
                   , sub       = NULL
                   , xlab      = NULL
                   , ylab      = NULL
                   , off.main  = 0       # Offset for main (relative, [0,1])
                   , off.sub   = 0       # Offset for sub  (relative, [0,1])
                   , off.xlab  = 0       # Offset for xlab (relative, [0,1])
                   , off.ylab  = 0       # Offset for ylab (relative, [0,1])
                   , off.right = 0       # Offset for the right margin (relative, [0,1])
                   , line.main = NA      # Line, same as title
                   , line.sub  = NA      # Line, same as title
                   , line.xlab = NA      # Line, same as title
                   , line.ylab = NA      # Line, same as title
                   , outer     = FALSE   # Added for compability, but it is ignored
                   ,...
                   ){
   #----- Save the original par settings. -------------------------------------------------#
   par.dims = par(no.readonly=FALSE)
   par.orig = par(no.readonly=TRUE )
   on.exit(par(par.orig))
   #---------------------------------------------------------------------------------------#


   #----- Grab the arguments from the ellipsis. -------------------------------------------#
   dotdotdot = list(...)
   #---------------------------------------------------------------------------------------#


   #----- Define outer margins based on device region. ------------------------------------#
   din = par.dims$din
   omi = c(max(off.xlab,off.sub)*din[2],off.ylab*din[1],off.main*din[2],off.right*din[1])
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Use the entire plotting area unless the title is to be placed somewhere inside.  #
   #---------------------------------------------------------------------------------------#
   if ("par.user" %in% ls()){
      if ("mar" %in% names(par.user)){
         mar = par.user$mar
      }else{
         mar = c(5.1,4.1,4.1,2.1)
      }#end if
      par.now     = modifyList(x=par.user,val=dotdotdot)
      par.now     = modifyList(x=par.now ,val=list(new=TRUE,bg=NULL,omi=omi,mar=mar))
   }else{
      mar     = c(5.1,4.1,4.1,2.1)
      par.now = modifyList(x=dotdotdot,list(new=TRUE,bg=NULL,omi=omi,mar=mar))
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Open a new window, then call title as usual, except that outer is forced to be     #
   # FALSE.                                                                                #
   #---------------------------------------------------------------------------------------#
   layout(mat=matrix(1,ncol=1,nrow=1))
   par(par.now)
   plot.new()
   plot.window(xlim=c(0,1),ylim=c(0,1))
   if (! is.null(main)) title(main = main,line=line.main       ,outer=FALSE,...)
   if (! is.null(sub )) title(sub  = sub ,line=line.sub ,font=2,outer=FALSE,...)
   if (! is.null(xlab)) title(xlab = xlab,line=line.xlab       ,outer=FALSE,...)
   if (! is.null(ylab)) title(ylab = ylab,line=line.ylab       ,outer=FALSE,...)
   #---------------------------------------------------------------------------------------#


   invisible()
}#end function gtitle
#==========================================================================================#
#==========================================================================================#
