#==========================================================================================#
#==========================================================================================#
#    cat0 is similar to cat, but it appends a "\n" and assumes sep="".                     #
#------------------------------------------------------------------------------------------#
cat0 <<- function(...,sep=""){
   cat(...,"\n",sep=sep)
   invisible()
}#end function cat0
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#   fill.bg fills in the background of pch that have background with some colour between   #
# col and the background.                                                                  #
# pch  - the pch                                                                           #
# col  - the main colour                                                                   #
# bg   - background colour                                                                 #
# dial - dial for background.  0 makes bg = col, 1 makes bg = bg, any numbers in between   #
#        will make something in between.  Values outside the range will extrapolate at     #
#        your own risk.                                                                    #
#                                                                                          #
#------------------------------------------------------------------------------------------#
fill.bg <<- function(pch,col,bg="white",dial=0.0){
   #----- Expand vectors so they all have the same size. ----------------------------------#
   nsz = max(length(pch),length(col),length(bg),length(dial))
   if (length(pch ) < nsz) pch  = rep(pch ,times=ceiling(nsz/length(pch )))[sequence(nsz)]
   if (length(col ) < nsz) col  = rep(col ,times=ceiling(nsz/length(col )))[sequence(nsz)]
   if (length(bg  ) < nsz) bg   = rep(bg  ,times=ceiling(nsz/length(bg  )))[sequence(nsz)]
   if (length(dial) < nsz) dial = rep(dial,times=ceiling(nsz/length(dial)))[sequence(nsz)]
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Make dial a data frame with three columns.                                       #
   #---------------------------------------------------------------------------------------#
   dial = data.frame(red=dial,green=dial,blue=dial)
   #---------------------------------------------------------------------------------------#


   #----- Turn col and bg into data frames. -----------------------------------------------#
   rgb.col = data.frame(t(col2rgb(col)))
   rgb.bg  = data.frame(t(col2rgb(bg )))
   #---------------------------------------------------------------------------------------#



   #----- Interpolate colours, make sure none of the intensities go outside range. --------#
   rgb.bg  = rgb.col + dial * (rgb.bg - rgb.col)
   rgb.bg  = data.frame(sapply(X=rgb.bg,FUN=function(x) pmax(0,pmin(255,round(x)))))
   #---------------------------------------------------------------------------------------#


   #----- Generate colour entry. ----------------------------------------------------------#
   rgb.bg  = mapply( FUN      = rgb
                   , red      = rgb.bg$red
                   , green    = rgb.bg$green
                   , blue     = rgb.bg$blue
                   , MoreArgs = list(maxColorValue=255)
                   )#end mapply
   #---------------------------------------------------------------------------------------#


   #----- Discard colours for those points without background support. --------------------#
   rgb.bg  = ifelse(test=pch %>% 20,yes=rgb.bg,no="transparent")
   #---------------------------------------------------------------------------------------#

   return(rgb.bg)
}#end fill.bg
#==========================================================================================#
#==========================================================================================#
