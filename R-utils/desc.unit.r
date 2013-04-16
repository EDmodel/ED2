#==========================================================================================#
#==========================================================================================#
#     This function creates an expression object that has the description and the units    #
# with sub-scripts, superscripts and stuff.                                                #
#------------------------------------------------------------------------------------------#
desc.unit <<- function(desc,unit){
   if (missing(desc) | missing(unit)){
      cat(" - Missing \"desc\": ",missing(desc),"\n")
      cat(" - Missing \"unit\": ",missing(unit),"\n")
      stop(" Both desc and unit must be provided")
   }#end if

   if (is.null(desc) && is.null(unit)){
      answer = ""
   }else if (is.null(unit)){
      answer = desc
   }else if (is.null(desc)){
      answer = parse(text=paste("paste(\"[\",",unit,",\"]\")",sep=""))
   }else{
      answer = parse(text=paste("paste(\"",desc,"\",\" [\",",unit,",\"]\")",sep=""))
   }#end if

   return(answer)
}#end function
#==========================================================================================#
#==========================================================================================#
