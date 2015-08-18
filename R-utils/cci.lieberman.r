#==========================================================================================#
#==========================================================================================#
#    Function cci.lieberman.                                                               #
#    This function computes the crown closure index for each individual, following the     #
# method presented by                                                                      #
#                                                                                          #
#     Lieberman, M., D. Lieberman, R. Peralta, G. S. Hartshorn, 1995: Canopy closure and   #
# the distribution of tropical forest tree species at La Selva, Costa Rica.  J. Trop.      #
# Ecol., 11 (2), 161--177.                                                                 #
#------------------------------------------------------------------------------------------#
cci.lieberman <<- function(xyz,radius=10,closure=TRUE){
   #---------------------------------------------------------------------------------------#
   #     Check whether the data are good to go.                                            #
   #---------------------------------------------------------------------------------------#
   if (! is.data.frame(xyz)){
      if (is.list(xyz)){
         xyz = try(list.2.data.frame(xyz),silent=TRUE)
         if ("try-error" %in% is(xyz)){
            cat0("1")
            browser()
            stop( paste0( "xyz must be a data frame, or an object that can "
                        , "be coerced to data frame"
                        )#end paste0
                )#end stop
         }#end if
      }else if (is.matrix(xyz)){
         xyz = data.frame(xyz)
      }else if (! is.data.frame(xyz)){
         cat0("3")
         browser()
         stop ("xyz must be a data frame, or an object that can be coerced to data frame")
      }#end if (is.list(xyz))
   }#end if (! is.data.frame(xyz))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Make up the names if they aren't provided.                                        #
   #---------------------------------------------------------------------------------------#
   if (! all(c("x","y","z") %in% names(xyz))){
      #------------------------------------------------------------------------------------#
      #    Check that we can make up x, y, and z                                           #
      #------------------------------------------------------------------------------------#
      if (length(xyz) < 3){
         stop("xyz must have at least three variables")
      }else{
         if (! "x" %in% names(xyz)) names(xyz)[which(! names(xyz) %in% "xyz")[1]] = "x"
         if (! "y" %in% names(xyz)) names(xyz)[which(! names(xyz) %in% "xyz")[1]] = "y"
         if (! "z" %in% names(xyz)) names(xyz)[which(! names(xyz) %in% "xyz")[1]] = "z"
      }#end if (length(xyz) < 3)
      #------------------------------------------------------------------------------------#
   }#end if (! all(c("x","y","z") %in% names(xyz)))
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Create a list with each element.                                                  #
   #---------------------------------------------------------------------------------------#
   i.xyz = split(x=xyz,f=sequence(nrow(xyz)))
   i.cci = mapply( FUN      = .Int.cci.lieberman
                 , i.xyz    = i.xyz
                 , MoreArgs = list(p.xyz=xyz,radius=radius,closure=closure)
                 , SIMPLIFY = TRUE
                 )#end mapply
   names(i.cci) = NULL
   return(i.cci)
   #---------------------------------------------------------------------------------------#
}#end cci.lieberman
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      This is the internal call of cci.lieberman.  It shouldn't be called directly by the #
# user, but from inside cci.lieberman (otherwise it wouldn't be called "internal").        #
#------------------------------------------------------------------------------------------#
.Int.cci.lieberman <<- function(i.xyz,p.xyz,radius,closure){


   #----- Make sure this function has been called by cloud.metrics or grid.metrics. -------#
   patt  = "^([A-Za-z0-9]+)(\\({1})(.*)(\\){1})$"
   repl  = "\\1"
   n     = 0
   mess  = TRUE
   top   = FALSE
   wcm   = list()
   while (! top){
      n = n + 1
      wcm[[n]] = try( gsub( pattern     = patt
                          , replacement = repl
                          , x           = deparse(sys.call(-n))
                          )#end gsub
                    , silent = TRUE
                    )#end try
      if ("try-error" %in% is(wcm[[n]])){
         wcm[[n]] = NA
         top      = TRUE
      }else{
         #----- Not an error.  Check whether this has been called by a friend function. ---#
         wcm[[n]] = paste(wcm[[n]],collapse="")
         top      = substring(wcm[[n]],1,4) %==% "NULL"
         mess     = mess && ! grepl("^cci.lieberman",wcm[[n]])
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#
   }#end while
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Do not allow this function to continue in case it was an externall call.          #
   #---------------------------------------------------------------------------------------#
   if (mess){
      wcm =  sapply(X=wcm,FUN=rbind)
      print(wcm)
      bye =  paste( " Function .Int.cci.lieberman is internal,"
                  , " and can only be called by cci.lieberman"
                  , sep = ""
                  )#end paste
      stop(bye)
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Calculate the index for this individual.                                          #
   #---------------------------------------------------------------------------------------#
   i.dr   = sqrt((i.xyz$x-p.xyz$x)^2 + (i.xyz$y-p.xyz$y)^2)
   i.dz   = p.xyz$z - i.xyz$z
   keep   = i.dr > 0 & i.dr <= radius
   p.xyz  = p.xyz[keep,]
   i.dr   = i.dr [keep ]
   i.dz   = i.dz [keep ]
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #    Find the index.                                                                    #
   #---------------------------------------------------------------------------------------#
   if (closure){
      sin.theta = pmax(0.,sin(i.dz/sqrt(i.dz^2+i.dr^2)))
   }else{
      sin.theta = pmax(0.,-sin(i.dz/sqrt(i.dz^2+i.dr^2)))
   }#end if (closure)
   cci = sum(sin.theta)
   #---------------------------------------------------------------------------------------#

   return(cci)
}#end function .Int.cci.lieberman
#==========================================================================================#
#==========================================================================================#
