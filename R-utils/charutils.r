#==========================================================================================#
#==========================================================================================#
#     This sub-routine makes the first letter of every entry capitalised, whilst leaving   #
# the other letters lower case.                                                            #
# Original subroutine comes from the help on toupper.  Here I modified it just to deal     #
# with NAs.                                                                                #
#------------------------------------------------------------------------------------------#
capwords <<- function(s, strict = FALSE) {

    #----- Function to be applied for each element of s. ----------------------------------#
    cap = function(x,strict=FALSE){

        #----- First letter is always upper case. -----------------------------------------#
        first  = toupper(substring(x,1,1))

        #----- Check whether to force the remainder to be lower case or not. --------------#
        if (strict){
           remainder = tolower(substring(x,2))
        }else{
           remainder = substring(x,2)
        }#end if
        #----------------------------------------------------------------------------------#
        ans      = paste(first,remainder,sep="",collapse=" ")
        return(ans)
    }#end function
    #--------------------------------------------------------------------------------------#


    #----- Remember which elements were NA, then we reset them to NA. ---------------------#
    sel    = is.na(s)
    #--------------------------------------------------------------------------------------#

    #----- Fix case for all dataset. ------------------------------------------------------#
    ans = sapply( X=strsplit(s, split = " "),FUN=cap,strict=strict
                , USE.NAMES = !is.null(names(s)))
    #--------------------------------------------------------------------------------------#
    
    #---- Force NAs to remain NAs. --------------------------------------------------------#
    ans[sel] = NA
    return(ans)
    #--------------------------------------------------------------------------------------#
}#end if
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#      This function deletes spaces from strings, trimming.  Default is to trim both       #
# sides, but you can also trim the left or the right only.                                 #
#------------------------------------------------------------------------------------------#
trim <<- function(x,side="both"){
   if (side %in% c("both","left") ) x = sub(pattern="^\\s+",replacement="",x=x)
   if (side %in% c("both","right")) x = sub(pattern="\\s+$",replacement="",x=x)
   return(x)
}#end function trim
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#       This function concatenates two strings, but skips the NAs.                         #
#------------------------------------------------------------------------------------------#
concatenate.message <<- function(x1,x2,sep="; "){
    if (length(x1) != length(x2)) stop("  Message vectors must have the same length!")

    only.x1 = ( ! is.na(x1) ) &     is.na(x2)
    only.x2 =     is.na(x1)   & ( ! is.na(x2) )
    both    = ( ! is.na(x1) ) & ( ! is.na(x2) )

    full          = rep(NA_character_,times=length(x1))
    full[only.x1] = x1[only.x1]
    full[only.x2] = x2[only.x2]
    full[both   ] = paste(x1[both],x2[both],sep=sep)

    return(full)
}#end function
#==========================================================================================#
#==========================================================================================#
