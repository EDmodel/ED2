#==========================================================================================#
#==========================================================================================#
#     Below are several colour palettes based on a single colour, useful for when one      #
# wants to plot net change.  You must say which colours to use for negative and positive.  #
#------------------------------------------------------------------------------------------#
two.palettes <<- function(x,n=20,white=1,low="blue",high="orangered"){

   #----- Make sure the name has an appropriate hue. --------------------------------------#
   bye = stopifnot(low  %in% c("orangered","green","brown","grey","blue","purple"))
   bye = stopifnot(high %in% c("orangered","green","brown","grey","blue","purple"))
   #---------------------------------------------------------------------------------------#



   #----- If x is a list, we use lapply to return the information for each element. -------#
   if (is.list(x) && ! is.data.frame(x)){
      ans = lapply(X=x,FUN=two.palettes,n=n,white=white,low=low,high=high)
   }else{

      #----- Break the values into bins, making 0 central. --------------------------------#
      if (sum(is.finite(x)) == 0){
         x.brks = c(0,0.5,1)
      }else{
         x.span = c(0,max(abs(x),na.rm=TRUE))
         if (all(x.span == 0)){
            x.brks = c(0,0.5,1)
         }else{
            x.brks = pretty(x.span,n=max(3,round(n/2)))
         }#end if
      }#end if
      x.brks = unique(c(-rev(x.brks),0,x.brks))
      n.brks = length(x.brks)
      x.cut  = cut(x,breaks = x.brks)
      x.lev  = levels(x.cut)

      #----- Find the indices that correspond to the colours. -----------------------------#
      x.idx  = match(x.cut,x.lev)
      if (is.data.frame(x)){
         x.idx = data.frame(x.idx,row.names=row.names(x))
         names(x.idx) = names(x)
      }else if (is.array(x)){
         x.idx = array(data=x.idx,dim=dim(x),dimnames=dimnames(x))
      }else{
         if (length(x) == length(x.idx)) names(x.idx) = names(x)
      }#end if
      #------------------------------------------------------------------------------------#



      #------ Find the size of each side. -------------------------------------------------#
      n.col   = n.brks-1
      n.each  = floor(n.col/2) - white
      n.white = n.col - 2*n.each
      #------------------------------------------------------------------------------------#


      #----- Get the correct colour palette. ----------------------------------------------#
      hue.low  = get(paste("hue",low ,sep="."))
      hue.high = get(paste("hue",high,sep="."))
      #------------------------------------------------------------------------------------#


      #----- Find the colours. ------------------------------------------------------------#
      col.low  = rev(hue.low(n=n.each))
      col.high = hue.high(n=n.each)
      col.out  = c(col.low,rep(background,n.white),col.high)
      #------------------------------------------------------------------------------------#



      #----- Build a list with the output. ------------------------------------------------#
      ans = list( idx       = x.idx
                , breaks    = x.brks
                , n.breaks  = n.brks
                , levels    = x.lev
                , colours   = col.out
                , n.colours = n.col
                )#end list
      #------------------------------------------------------------------------------------#
   }#end if

   return(ans)
}#end function two.palettes
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#      Below is the list of colour palettes by the hue.                                    #
#------------------------------------------------------------------------------------------#
#----- Blue. ------------------------------------------------------------------------------#
hue.blue <<- function(n,inv=FALSE){
   nodes     = c("skyblue","deepskyblue","dodgerblue","royalblue3","midnightblue")
   nodes     = data.frame(t(col2rgb(nodes)))
   pivot     = round(seq(from=1,to=n,length.out=nrow(nodes)),digits=0)
   rgb.out   = data.frame(t(mapply(FUN=spline,y=nodes,MoreArgs=list(x=pivot,n=n))))$y
   rgb.out   = lapply(X=rgb.out,FUN=as.integer)
   rgb.out   = lapply(X=rgb.out,FUN=pmax,  0)
   rgb.out   = lapply(X=rgb.out,FUN=pmin,255)
   rgb.out   = rgb(r=rgb.out$red,g=rgb.out$green,b=rgb.out$blue,maxColorValue=255)
   return(rgb.out)
}#end hue.blue
#----- Orange-Red. ------------------------------------------------------------------------#
hue.orangered <<- function(n,inv=FALSE){
   nodes     = c("gold","goldenrod","darkorange1","orangered","#7E0000")
   nodes     = data.frame(t(col2rgb(nodes)))
   pivot     = round(seq(from=1,to=n,length.out=nrow(nodes)),digits=0)
   rgb.out   = data.frame(t(mapply(FUN=spline,y=nodes,MoreArgs=list(x=pivot,n=n))))$y
   rgb.out   = lapply(X=rgb.out,FUN=as.integer)
   rgb.out   = lapply(X=rgb.out,FUN=pmax,  0)
   rgb.out   = lapply(X=rgb.out,FUN=pmin,255)
   rgb.out   = rgb(r=rgb.out$red,g=rgb.out$green,b=rgb.out$blue,maxColorValue=255)
   return(rgb.out)
}#end hue.orangered
#----- Green. -----------------------------------------------------------------------------#
hue.green <<- function(n,inv=FALSE){
   nodes     = c("darkolivegreen1","olivedrab3","chartreuse2","forestgreen","#004000")
   nodes     = data.frame(t(col2rgb(nodes)))
   pivot     = round(seq(from=1,to=n,length.out=nrow(nodes)),digits=0)
   rgb.out   = data.frame(t(mapply(FUN=spline,y=nodes,MoreArgs=list(x=pivot,n=n))))$y
   rgb.out   = lapply(X=rgb.out,FUN=as.integer)
   rgb.out   = lapply(X=rgb.out,FUN=pmax,  0)
   rgb.out   = lapply(X=rgb.out,FUN=pmin,255)
   rgb.out   = rgb(r=rgb.out$red,g=rgb.out$green,b=rgb.out$blue,maxColorValue=255)
   return(rgb.out)
}#end hue.orangered
#----- Brown. -----------------------------------------------------------------------------#
hue.brown <<- function(n){
   nodes     = c("lightgoldenrod2","sandybrown","sienna3","saddlebrown","#603000")
   nodes     = data.frame(t(col2rgb(nodes)))
   pivot     = round(seq(from=1,to=n,length.out=nrow(nodes)),digits=0)
   rgb.out   = data.frame(t(mapply(FUN=spline,y=nodes,MoreArgs=list(x=pivot,n=n))))$y
   rgb.out   = lapply(X=rgb.out,FUN=as.integer)
   rgb.out   = lapply(X=rgb.out,FUN=pmax,  0)
   rgb.out   = lapply(X=rgb.out,FUN=pmin,255)
   rgb.out   = rgb(r=rgb.out$red,g=rgb.out$green,b=rgb.out$blue,maxColorValue=255)
   return(rgb.out)
}#end hue.orangered
#----- Grey. ------------------------------------------------------------------------------#
hue.grey <<- function(n){
   nodes     = c("grey84","grey68","grey52","grey36","grey20")
   nodes     = data.frame(t(col2rgb(nodes)))
   pivot     = round(seq(from=1,to=n,length.out=nrow(nodes)),digits=0)
   rgb.out   = data.frame(t(mapply(FUN=spline,y=nodes,MoreArgs=list(x=pivot,n=n))))$y
   rgb.out   = lapply(X=rgb.out,FUN=as.integer)
   rgb.out   = lapply(X=rgb.out,FUN=pmax,  0)
   rgb.out   = lapply(X=rgb.out,FUN=pmin,255)
   rgb.out   = rgb(r=rgb.out$red,g=rgb.out$green,b=rgb.out$blue,maxColorValue=255)
   return(rgb.out)
}#end hue.orangered
#----- Purple. ----------------------------------------------------------------------------#
hue.purple <<- function(n){
   nodes     = c("#DEE1FF","#B0BCFF","#7E76FF","#6B4EFF"
                ,"#6337C9","#581BA2","#440084","#2A0053")
   nodes     = c("#DEE1FF","#B0BCFF","#7E76FF","#6B4EFF"
                ,"#581BA2","#2A0053")
   nodes     = data.frame(t(col2rgb(nodes)))
   pivot     = round(seq(from=1,to=n,length.out=nrow(nodes)),digits=0)
   rgb.out   = data.frame(t(mapply(FUN=spline,y=nodes,MoreArgs=list(x=pivot,n=n))))$y
   rgb.out   = lapply(X=rgb.out,FUN=as.integer)
   rgb.out   = lapply(X=rgb.out,FUN=pmax,  0)
   rgb.out   = lapply(X=rgb.out,FUN=pmin,255)
   rgb.out   = rgb(r=rgb.out$red,g=rgb.out$green,b=rgb.out$blue,maxColorValue=255)
   return(rgb.out)
}#end hue.orangered
#------ Inverted scale. -------------------------------------------------------------------#
ihue.blue      <<- function(n) rev(hue.blue     (n))
ihue.orangered <<- function(n) rev(hue.orangered(n))
ihue.green     <<- function(n) rev(hue.green    (n))
ihue.brown     <<- function(n) rev(hue.brown    (n))
ihue.grey      <<- function(n) rev(hue.grey     (n))
ihue.purple    <<- function(n) rev(hue.purple   (n))
#------------------------------------------------------------------------------------------#
