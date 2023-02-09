##
##Custom functions for APA regulation in CB
##
##


PlotHeatmapRegulators <- function(mots=unique(unlist(global_positive)),mymats=list(ttest_APA_NGF_Ip[,-1],-ttest_APA_NGF_Id[,-1]),scaling=TRUE){
  
  if(scaling){
    mymats[[1]]<- mymats[[1]]/max(abs(mymats[[1]]),na.rm=TRUE)
    mymats[[2]]<- mymats[[2]]/max(abs(mymats[[2]]),na.rm=TRUE)
  }
  
  mymats                <- lapply(mymats,function(Z)return(Z[match(mots,rownames(Z)),]))
  GS                    <- unlist(lapply(mots,function(Z)return((unlist(strsplit(Z,split="_"))[1]))))
  mat                   <- cbind(mymats[[1]],mymats[[2]])[,]
  rownames(mymats[[1]]) <-rownames(mymats[[2]])<-rownames(mat)<-GS
  mat                   <- mat[!duplicated(GS),]
  dd                    <- as.matrix(dist(mymats[[1]][!duplicated(GS),],method="man"))
  diag(dd)              <- 0
  dd.row                <- as.dendrogram(hclust(as.dist(dd),method="ward.D"))
  row.ord               <- order.dendrogram(dd.row)
  require("RColorBrewer")
  require("gplots")
  mypalette             <-  colorRampPalette(colors= brewer.pal(n = 6, name = "PRGn"), bias = 1, space = c("rgb"), interpolate = c("spline"))
  mycols                <-  rev(mypalette(n=100))
  b.1                   <-  seq(from=-max(abs(mymats[[1]])),to=max(abs(mymats[[1]])),length.out=101)
  
  #b.1                   <-  seq(from=-max(abs(mymats[[2]])),to=max(abs(mymats[[2]])),length.out=101)
  
  
  #b.1                   <-  seq(from=-30,to=30,length.out=101)
  heatmap.2(mat[row.ord,], keysize=1,mar=c(6,5),col=mycols,breaks=b.1, scale="none",Rowv=FALSE,Colv=NA,key=TRUE,symkey=FALSE, density.info="none", trace="none",cexCol=0.7, cexRow=0.7, font.lab=1,dendrogram="row")
  
}


ExtractTopRegulatorsSimple <- function(sub2=t(standardize(t(cbind(as.numeric(ttest_transport_NGF[,colix1[2]]),as.numeric(Zscore_NGF_overtransported[,colix2[2]])))))){
  
  lims<- boxplot(sub2,plot=FALSE)$stats[5,]
  return(which(sub2[,1]>lims[1]&sub2[,2]>lims[2]))
}



ExtractTopRegulators <- function(n=10,sub2=t(standardize(t(cbind(as.numeric(ttest_transport_NGF[,colix1[2]]),as.numeric(Zscore_NGF_overtransported[,colix2[2]])))))){
  
  IX_sel <- which(apply(sub2>5.0,1,sum)==2)
  
  standardize <- function(z){
    rowmed <- apply(z, 1, mean)
    rowmad <- apply(z, 1, sd)
    rv <- sweep(z, 1, rowmed)
    rv <- sweep(rv, 1, rowmad, "/")
    return(rv)
  }
  sub2 <- t(standardize(t(sub2)))
  km5   <- kmeans(sub2, centers = n, nstart=25)
  Lev3  <- levels(as.factor(km5$cluster))
  Cols  <- c("grey","black","red","blue","green","cyan","magenta","yellow","darkgreen","lightpink")
  names(Cols)<- Lev3
  cols3 <- unlist(lapply(km5$cluster,function(x)return(Cols[match(x,names(Cols))])))
  
  plot(sub2[,1],sub2[,2],col=cols3,pch=19,cex=0.7,xlab="",ylab="",frame=FALSE,las=1)
  mtext(side=1,line=2,text="Ip")
  mtext(side=2,line=2,text="Id")
  grid()
  
  ix_max<- which(apply(sub2,1,sum)==max(apply(sub2,1,sum)))[1]
  ix_min<- which(apply(sub2,1,sum)==min(apply(sub2,1,sum)))[1]
  
  cl_max <- km5$cluster[ix_max]
  cl_min <- km5$cluster[ix_min]
  
  selmost <- intersect(which(km5$cluster==cl_max),IX_sel)
  selmin  <- intersect(which(km5$cluster==cl_min),IX_sel)
  return(list(selmost,selmin))
}


plotEnrichAlongUTR <- function(mots="CPSF160_iClip"){
  par(mfrow=c(3,2))
  plot(x=-unlist(lapply(myranges,mean))[-1],y=ttest_APA_NGF_Ip[match(mots,rownames(ttest_APA_NGF_Ip)),-1],main="Proximal",col=colsNGF[1],type="l",lwd=1.5,las=1,ylab="-log10(P-value)",frame=FALSE,xlab="distante to 3'end",ylim=c(min(0,min(ttest_APA_NGF_Ip[match(mots,rownames(ttest_APA_NGF_Ip)),-1],na.rm=TRUE)),max(0,max(ttest_APA_NGF_Ip[match(mots,rownames(ttest_APA_NGF_Ip)),-1],na.rm=TRUE))))
  points(x=-unlist(lapply(myranges,mean))[-1],y=ttest_APA_NGF_Ip[match(mots,rownames(ttest_APA_NGF_Ip)),-1],col=colsNGF[1],pch=19)
  mtext(side=2,line=2,text=mots,font=1,cex=0.5)
  abline(h=0,col="black")
  abline(v=0,col="darkgrey",lty=2)
  
  plot(x=-unlist(lapply(myranges,mean))[-1],y=-ttest_APA_NGF_Id[match(mots,rownames(ttest_APA_NGF_Id)),-1],main="Distal",col=colsNGF[1],type="l",lwd=1.5,las=1,ylab="-log10(P-value)",frame=FALSE,xlab="distante to 3'end",ylim=c(min(0,min(-ttest_APA_NGF_Id[match(mots,rownames(ttest_APA_NGF_Id)),-1],na.rm=TRUE)),max(0,max(-ttest_APA_NGF_Id[match(mots,rownames(ttest_APA_NGF_Id)),-1],na.rm=TRUE))))
  points(x=-unlist(lapply(myranges,mean))[-1],y=-ttest_APA_NGF_Id[match(mots,rownames(ttest_APA_NGF_Id)),-1],col=colsNGF[1],pch=19)
  mtext(side=2,line=2,text=mots,font=1,cex=0.5)
  abline(h=0,col="black")
  abline(v=0,col="darkgrey",lty=2)
  
  plot(x=-unlist(lapply(myranges,mean))[-1],y=ttest_APA_NT3_Ip[match(mots,rownames(ttest_APA_NT3_Ip)),-1],main="Proximal",col=colsNT3[1],type="l",lwd=1.5,las=1,ylab="-log10(P-value)",frame=FALSE,xlab="distante to 3'end",ylim=c(min(0,min(ttest_APA_NT3_Ip[match(mots,rownames(ttest_APA_NT3_Ip)),-1],na.rm=TRUE)),max(0,max(ttest_APA_NT3_Ip[match(mots,rownames(ttest_APA_NT3_Ip)),-1],na.rm=TRUE))))
  points(x=-unlist(lapply(myranges,mean))[-1],y=ttest_APA_NT3_Ip[match(mots,rownames(ttest_APA_NT3_Ip)),-1],col=colsNT3[1],pch=19)
  mtext(side=2,line=2,text=mots,font=1,cex=0.5)
  abline(h=0,col="black")
  abline(v=0,col="darkgrey",lty=2)
  
  plot(x=-unlist(lapply(myranges,mean))[-1],y=-ttest_APA_NT3_Id[match(mots,rownames(ttest_APA_NT3_Id)),-1],main="Distal",col=colsNT3[1],type="l",lwd=1.5,las=1,ylab="-log10(P-value)",frame=FALSE,xlab="distante to 3'end",ylim=c(min(0,min(-ttest_APA_NT3_Id[match(mots,rownames(ttest_APA_NT3_Id)),-1],na.rm=TRUE)),max(0,max(-ttest_APA_NT3_Id[match(mots,rownames(ttest_APA_NT3_Id)),-1],na.rm=TRUE))))
  points(x=-unlist(lapply(myranges,mean))[-1],y=-ttest_APA_NT3_Id[match(mots,rownames(ttest_APA_NT3_Id)),-1],col=colsNT3[1],pch=19)
  mtext(side=2,line=2,text=mots,font=1,cex=0.5)
  abline(h=0,col="black")
  abline(v=0,col="darkgrey",lty=2)
  
  
  plot(x=-unlist(lapply(myranges,mean))[-1],y=ttest_diffAPA_Ip[match(mots,rownames(ttest_diffAPA_Ip)),-1],main="Proximal",col="grey",type="l",lwd=1.5,las=1,ylab="-log10(P-value)",frame=FALSE,xlab="distante to 3'end",ylim=c(min(0,min(ttest_diffAPA_Ip[match(mots,rownames(ttest_diffAPA_Ip)),-1],na.rm=TRUE)),max(0,max(ttest_diffAPA_Ip[match(mots,rownames(ttest_diffAPA_Ip)),-1],na.rm=TRUE))))
  abline(h=0,col="black")
  abline(v=0,col="darkgrey",lty=2)
  points(x=-unlist(lapply(myranges,mean))[-1],y=ttest_diffAPA_Ip[match(mots,rownames(ttest_diffAPA_Ip)),-1],col="grey",pch=19)
  mtext(side=2,line=2,text=mots,font=1,cex=0.5)
  
  plot(x=-unlist(lapply(myranges,mean))[-1],y=ttest_diffAPA_Id[match(mots,rownames(ttest_diffAPA_Id)),-1],main="Distal",col="grey",type="l",lwd=1.5,las=1,ylab="-log10(P-value)",frame=FALSE,xlab="distante to 3'end",ylim=c(min(0,min(ttest_diffAPA_Id[match(mots,rownames(ttest_diffAPA_Id)),-1],na.rm=TRUE)),max(0,max(ttest_diffAPA_Id[match(mots,rownames(ttest_diffAPA_Id)),-1],na.rm=TRUE))))
  points(x=-unlist(lapply(myranges,mean))[-1],y=ttest_diffAPA_Id[match(mots,rownames(ttest_diffAPA_Id)),-1],col="grey",pch=19)
  mtext(side=2,line=2,text=mots,font=1,cex=0.5)
  abline(h=0,col="black")
  abline(v=0,col="darkgrey",lty=2)
  
}

plotEnrichAlongUTRfocus <- function(mots="CPSF160_iClip",sel=c(20:30)){
  par(mfrow=c(3,2))
  plot(x=-unlist(lapply(myranges,mean))[sel],y=ttest_APA_NGF_Ip[match(mots,rownames(ttest_APA_NGF_Ip)),sel],main="Proximal",col=colsNGF[1],type="l",lwd=1.5,las=1,ylab="-log10(P-value)",frame=FALSE,xlab="distante to 3'end",ylim=c(min(0,min(ttest_APA_NGF_Ip[match(mots,rownames(ttest_APA_NGF_Ip)),sel],na.rm=TRUE)),max(0,max(ttest_APA_NGF_Ip[match(mots,rownames(ttest_APA_NGF_Ip)),sel],na.rm=TRUE))))
  points(x=-unlist(lapply(myranges,mean))[sel],y=ttest_APA_NGF_Ip[match(mots,rownames(ttest_APA_NGF_Ip)),sel],col=colsNGF[1],pch=19)
  mtext(side=2,line=2,text=mots,font=1,cex=0.5)
  abline(v=0,lty=3)
  abline(h=0,col="black")
  plot(x=-unlist(lapply(myranges,mean))[sel],y=-ttest_APA_NGF_Id[match(mots,rownames(ttest_APA_NGF_Id)),sel],main="Distal",col=colsNGF[1],type="l",lwd=1.5,las=1,ylab="-log10(P-value)",frame=FALSE,xlab="distante to 3'end",ylim=c(min(0,min(-ttest_APA_NGF_Id[match(mots,rownames(ttest_APA_NGF_Id)),sel],na.rm=TRUE)),max(0,max(-ttest_APA_NGF_Id[match(mots,rownames(ttest_APA_NGF_Id)),sel],na.rm=TRUE))))
  points(x=-unlist(lapply(myranges,mean))[sel],y=-ttest_APA_NGF_Id[match(mots,rownames(ttest_APA_NGF_Id)),sel],col=colsNGF[1],pch=19)
  mtext(side=2,line=2,text=mots,font=1,cex=0.5)
  abline(v=0,lty=3)
  abline(h=0,col="black")
  
  plot(x=-unlist(lapply(myranges,mean))[sel],y=ttest_APA_NT3_Ip[match(mots,rownames(ttest_APA_NT3_Ip)),sel],main="Proximal",col=colsNT3[1],type="l",lwd=1.5,las=1,ylab="-log10(P-value)",frame=FALSE,xlab="distante to 3'end",ylim=c(min(0,min(ttest_APA_NT3_Ip[match(mots,rownames(ttest_APA_NT3_Ip)),sel],na.rm=TRUE)),max(0,max(ttest_APA_NT3_Ip[match(mots,rownames(ttest_APA_NT3_Ip)),sel],na.rm=TRUE))))
  points(x=-unlist(lapply(myranges,mean))[sel],y=ttest_APA_NT3_Ip[match(mots,rownames(ttest_APA_NT3_Ip)),sel],col=colsNT3[1],pch=19)
  mtext(side=2,line=2,text=mots,font=1,cex=0.5)
  abline(v=0,lty=3)
  abline(h=0,col="black")
  
  plot(x=-unlist(lapply(myranges,mean))[sel],y=-ttest_APA_NT3_Id[match(mots,rownames(ttest_APA_NT3_Id)),sel],main="Distal",col=colsNT3[1],type="l",lwd=1.5,las=1,ylab="-log10(P-value)",frame=FALSE,xlab="distante to 3'end",ylim=c(min(0,min(-ttest_APA_NT3_Id[match(mots,rownames(ttest_APA_NT3_Id)),sel],na.rm=TRUE)),max(0,max(-ttest_APA_NT3_Id[match(mots,rownames(ttest_APA_NT3_Id)),sel],na.rm=TRUE))))
  points(x=-unlist(lapply(myranges,mean))[sel],y=-ttest_APA_NT3_Id[match(mots,rownames(ttest_APA_NT3_Id)),sel],col=colsNT3[1],pch=19)
  mtext(side=2,line=2,text=mots,font=1,cex=0.5)
  abline(v=0,lty=3)
  abline(h=0,col="black")
  
  
  plot(x=-unlist(lapply(myranges,mean))[sel],y=ttest_diffAPA_Ip[match(mots,rownames(ttest_diffAPA_Ip)),sel],main="Proximal",col="grey",type="l",lwd=1.5,las=1,ylab="-log10(P-value)",frame=FALSE,xlab="distante to 3'end",ylim=c(min(0,min(ttest_diffAPA_Ip[match(mots,rownames(ttest_diffAPA_Ip)),sel],na.rm=TRUE)),max(0,max(ttest_diffAPA_Ip[match(mots,rownames(ttest_diffAPA_Ip)),sel],na.rm=TRUE))))
  points(x=-unlist(lapply(myranges,mean))[sel],y=ttest_diffAPA_Ip[match(mots,rownames(ttest_diffAPA_Ip)),sel],col="grey",pch=19)
  abline(v=0,lty=3)
  abline(h=0,col="black")
  
  mtext(side=2,line=2,text=mots,font=1,cex=0.5)
  plot(x=-unlist(lapply(myranges,mean))[sel],y=ttest_diffAPA_Id[match(mots,rownames(ttest_diffAPA_Id)),sel],main="Distal",col="grey",type="l",lwd=1.5,las=1,ylab="-log10(P-value)",frame=FALSE,xlab="distante to 3'end",ylim=c(min(0,min(ttest_diffAPA_Id[match(mots,rownames(ttest_diffAPA_Id)),sel],na.rm=TRUE)),max(0,max(ttest_diffAPA_Id[match(mots,rownames(ttest_diffAPA_Id)),sel],na.rm=TRUE))))
  points(x=-unlist(lapply(myranges,mean))[sel],y=ttest_diffAPA_Id[match(mots,rownames(ttest_diffAPA_Id)),sel],col="grey",pch=19)
  mtext(side=2,line=2,text=mots,font=1,cex=0.5)
  abline(v=0,lty=3)
  abline(h=0,col="black")
  
}



plotEnrichAlongUTRfocusMerged <- function(mots="CPSF160_iClip",sel=c(20:30),cond="NGF",scaling=TRUE){

  if(cond=="NGF"){
    myX  = -unlist(lapply(myranges,mean))[sel]
    myIp = ttest_APA_NGF_Ip[match(mots,rownames(ttest_APA_NGF_Ip)),sel]
    myId = ttest_APA_NGF_Id[match(mots,rownames(ttest_APA_NGF_Id)),sel]
    if(scaling){
      myIp <- myIp/max(abs(myIp))
      myId <- myId/max(abs(myId))
    }
    plot(myX,myIp,col="white",las=1,ylab="-log10(P-value)",frame=FALSE,xlab="distante to 3'end",ylim=c(min(c(myIp,myId),na.rm=TRUE),max(c(myIp,myId),na.rm=TRUE)))
    abline(h=0,col="black")
    abline(v=0,col="red",lty=2)
    lines(myX,myIp,col="black",lwd=1.5)
    lines(myX,myId,col="grey",lwd=1.5)
    points(myX,myIp,col="black",pch=19)
    points(myX,myId,col="grey",pch=19)
    mtext(side=2,line=2,text=mots,font=1,cex=0.5)
    mtext(side=3,line=0,text="NGF",font=1,cex=0.5)
  }
  
  
  if(cond=="NT3"){
    myX  = -unlist(lapply(myranges,mean))[sel]
    myIp = ttest_APA_NT3_Ip[match(mots,rownames(ttest_APA_NT3_Ip)),sel]
    myId = ttest_APA_NT3_Id[match(mots,rownames(ttest_APA_NT3_Id)),sel]
    plot(myX,myIp,col="white",las=1,ylab="-log10(P-value)",frame=FALSE,xlab="distante to 3'end",ylim=c(min(c(myIp,myId),na.rm=TRUE),max(c(myIp,myId),na.rm=TRUE)))
    abline(h=0,col="black")
    abline(v=0,col="red",lty=2)
    lines(myX,myIp,col="black",lwd=1.5)
    lines(myX,myId,col="grey",lwd=1.5)
    points(myX,myIp,col="black",pch=19)
    points(myX,myId,col="grey",pch=19)
    mtext(side=2,line=2,text=mots,font=1,cex=0.5)
    mtext(side=3,line=0,text="NT3",font=1,cex=0.5)
  }
  
  if(cond=="diff"){
    myX  = -unlist(lapply(myranges,mean))[sel]
    myIp = ttest_diffAPA_Ip[match(mots,rownames(ttest_diffAPA_Ip)),sel]
    myId = ttest_diffAPA_Ip[match(mots,rownames(ttest_diffAPA_Ip)),sel]
    plot(myX,myIp,col="white",las=1,ylab="-log10(P-value)",frame=FALSE,xlab="distante to 3'end",ylim=c(min(c(myIp,myId),na.rm=TRUE),max(c(myIp,myId),na.rm=TRUE)))
    abline(h=0,col="black")
    abline(v=0,col="red",lty=2)
    lines(myX,myIp,col="black",lwd=1.5)
    lines(myX,myId,col="grey",lwd=1.5)
    points(myX,myIp,col="black",pch=19)
    points(myX,myId,col="grey",pch=19)
    mtext(side=2,line=2,text=mots,font=1,cex=0.5)
    mtext(side=3,line=0,text="diff",font=1,cex=0.5)
  }
}



