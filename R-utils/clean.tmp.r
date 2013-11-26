#==========================================================================================#
#==========================================================================================#
#      This function is a simple workaround for R version 2.15.0.  It deletes the          #
# temporary PDF files from the temp directory.                                             #
#------------------------------------------------------------------------------------------#
clean.tmp <<- function(){
   unlink(list.files(path=tempdir(), pattern="^pdf.", full.names=TRUE))
}#end function clean.tmp
#------------------------------------------------------------------------------------------#
