#==========================================================================================#
#==========================================================================================#
#     This function diverts the standard output to nowhere.                                #
#------------------------------------------------------------------------------------------#
shhh <<- function(fun,skip.error=TRUE,...){
   dotdotdot = list(...)

   nowhere   = tempfile()
   dummy     = sink(file=nowhere)
   
   ans       = try(expr=do.call(what=fun,args=dotdotdot),silent=TRUE)
   dummy     = sink(file=NULL)
   dummy     = file.remove(nowhere)
   if ("try-error" %in% ans && ! skip.error ) do.call(what=fun,args=dotdotdot)
   return(ans)
}#end shhh
#==========================================================================================#
#==========================================================================================#
