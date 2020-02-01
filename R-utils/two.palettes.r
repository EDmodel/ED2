#==========================================================================================#
#==========================================================================================#
#     Below are several colour palettes based on a single colour, useful for when one      #
# wants to plot net change.  You must say which colours to use for negative and positive.  #
#------------------------------------------------------------------------------------------#
two.palettes <<- function(x,n=20,white=1,low="blue",high="orangered",zero=NULL){

   #----- Make sure the name has an appropriate hue. --------------------------------------#
   plist = c("orangered","green","bluegreen","brown","grey","blue","purple")
   bye   = stopifnot(low  %in% plist)
   bye   = stopifnot(high %in% plist)
   #---------------------------------------------------------------------------------------#


   #----- Make default background colour in case none is given. ---------------------------#
   if (is.null(zero)) zero = background
   #---------------------------------------------------------------------------------------#



   #----- If x is a list, we use lapply to return the information for each element. -------#
   if (is.list(x) && ! is.data.frame(x)){
      ans = lapply(X=x,FUN=two.palettes,n=n,white=white,low=low,high=high,zero=zero)
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
      col.out  = c(col.low,rep(zero,n.white),col.high)
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
hue.blue <<- function(n){
   nodes     = c("#D1DCFF","#B4E6FF","#46B4FF","#0082C8","#003264")
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
hue.orangered <<- function(n){
   nodes     = c("#FEEEB8","#FDCE87","#FFB43C","#C85A0A","#960000")
   nodes     = data.frame(t(col2rgb(nodes)))
   pivot     = round(seq(from=1,to=n,length.out=nrow(nodes)),digits=0)
   rgb.out   = data.frame(t(mapply(FUN=spline,y=nodes,MoreArgs=list(x=pivot,n=n))))$y
   rgb.out   = lapply(X=rgb.out,FUN=as.integer)
   rgb.out   = lapply(X=rgb.out,FUN=pmax,  0)
   rgb.out   = lapply(X=rgb.out,FUN=pmin,255)
   rgb.out   = rgb(r=rgb.out$red,g=rgb.out$green,b=rgb.out$blue,maxColorValue=255)
   return(rgb.out)
}#end hue.orangered
#----- BlueGreen. -----------------------------------------------------------------------------#
hue.bluegreen <<- function(n){
   nodes     = c("#C7EAE5","#80CDC1","#53B3A8","#35978F","#01665E")
   nodes     = data.frame(t(col2rgb(nodes)))
   pivot     = round(seq(from=1,to=n,length.out=nrow(nodes)),digits=0)
   rgb.out   = data.frame(t(mapply(FUN=spline,y=nodes,MoreArgs=list(x=pivot,n=n))))$y
   rgb.out   = lapply(X=rgb.out,FUN=as.integer)
   rgb.out   = lapply(X=rgb.out,FUN=pmax,  0)
   rgb.out   = lapply(X=rgb.out,FUN=pmin,255)
   rgb.out   = rgb(r=rgb.out$red,g=rgb.out$green,b=rgb.out$blue,maxColorValue=255)
   return(rgb.out)
}#end hue.orangered
#----- BlueGreen. -----------------------------------------------------------------------------#
hue.green <<- function(n){
   nodes     = c("#D6ECD4","#B5E2B0","#83CC89","#3DAA5E","#008933")
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
   nodes     = c("#F6E8C3","#DFC27D","#D0A14D","#BF812D","#8C510A")
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
   nodes     = c("#E2D7E7","#CDB5D8","#B58DC4","#985EAE","#7C2D96")
   nodes     = data.frame(t(col2rgb(nodes)))
   pivot     = round(seq(from=1,to=n,length.out=nrow(nodes)),digits=0)
   rgb.out   = data.frame(t(mapply(FUN=spline,y=nodes,MoreArgs=list(x=pivot,n=n))))$y
   rgb.out   = lapply(X=rgb.out,FUN=as.integer)
   rgb.out   = lapply(X=rgb.out,FUN=pmax,  0)
   rgb.out   = lapply(X=rgb.out,FUN=pmin,255)
   rgb.out   = rgb(r=rgb.out$red,g=rgb.out$green,b=rgb.out$blue,maxColorValue=255)
   return(rgb.out)
}#end hue.orangered
#----- Green, colour-blind friendly. ------------------------------------------------------#
hue.green.cbf <<- function(n){
   nodes     = c("#D2EAE9","#B0F3F2","#50DFDC","#31A29F","#185957")
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
