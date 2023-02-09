
GetEnrichDiff<- function(goID="GO:0010574"){
  go.entrez       <- genesInTerm(mysampleGO, goID)
  if(length(go.entrez)>0){
    mySize          <- length(go.entrez[[1]])
    go.entrez       <- go.entrez[[1]]#To extract all genes related to this term
    go.gs           <- unique(as.character(t2g$external_gene_name)[which(t2g$entrezgene%in%go.entrez)])
    
    is.in.go        <- as.numeric(names(myVec[[2]])%in%go.gs)
    is.in.go        <- cbind(is.in.go,replicate(500,is.in.go[sample(c(1:length(is.in.go)))]))
    
    
    #Norm
    my.val.ngf.norm        <- apply(myVec[[2]]*is.in.go,2,sum)
    my.val.nt3.norm        <- apply(myVec[[6]]*is.in.go,2,sum)
    my.diff.norm           <- my.val.ngf.norm-my.val.nt3.norm
    zscore.norm            <- (my.diff.norm[1]-mean(my.diff.norm[-1]))/sd(my.diff.norm[-1])
    
    #Error
    my.val.ngf.error        <- apply(myVec[[4]]*is.in.go,2,sum)
    my.val.nt3.error        <- apply(myVec[[8]]*is.in.go,2,sum)
    my.diff.error           <- my.val.ngf.error-my.val.nt3.error
    zscore.error            <- (my.diff.error[1]-mean(my.diff.error[-1]))/sd(my.diff.error)
    
    
    return(c(my.diff.norm[1],zscore.norm,my.diff.error[1],zscore.error,length(go.gs)))
  }
  else{
    return(NA)
  }
}

GetEnrichTransportFull <- function(no.sites=mycount.all.penta.big[,1]){
  no.sites              <- cbind(no.sites,replicate(500,no.sites[sample(c(1:length(no.sites)))]))
  
  #Norm
  my.val.ngf.norm        <- apply(myVec[[2]]*no.sites,2,sum)
  zscore.ngf.norm        <- (my.val.ngf.norm[1]-mean(my.val.ngf.norm[-1]))/sd(my.val.ngf.norm[-1])    
  my.val.nt3.norm        <- apply(myVec[[6]]*no.sites,2,sum)
  zscore.nt3.norm        <- (my.val.nt3.norm[1]-mean(my.val.nt3.norm[-1]))/sd(my.val.nt3.norm[-1])    
  my.diff.norm           <- my.val.ngf.norm-my.val.nt3.norm
  zscore.norm            <- (my.diff.norm[1]-mean(my.diff.norm[-1]))/sd(my.diff.norm[-1])
  
  #Error
  my.val.ngf.error        <- apply(myVec[[4]]*no.sites,2,sum)
  zscore.ngf.error        <- (my.val.ngf.error[1]-mean(my.val.ngf.error[-1]))/sd(my.val.ngf.error[-1])
  my.val.nt3.error        <- apply(myVec[[8]]*no.sites,2,sum)
  zscore.nt3.error        <- (my.val.nt3.error[1]-mean(my.val.nt3.error[-1]))/sd(my.val.nt3.error[-1])
  my.diff.error           <- my.val.ngf.error-my.val.nt3.error
  zscore.error            <- (my.diff.error[1]-mean(my.diff.error[-1]))/sd(my.diff.error[-1])
  
  
  return(c(zscore.ngf.norm,zscore.nt3.norm,my.diff.norm[1],zscore.norm,
           zscore.ngf.error,zscore.nt3.error,my.diff.error[1],zscore.error))
}


PlotEnrichNoSites<- function(no.sites,colid="p.val.max",myName=MOT,myN=250){
  require(scales)
  par(mfrow=c(2,2),mar=c(3,4,4,3))
  myX <- myOut.ngf[no.sites==0,match(colid,colnames(myOut.ngf))]
  myY <- myOut.ngf[no.sites>0,match(colid,colnames(myOut.ngf))] 
  boxplot(list(myX,myY),outline=F,frame=F,las=1,col=c("white","grey"),xaxt="n")
  mtext(side=3,line=0,cex=0.6,text=paste("P(t-test)=",scientific_format(3)(t.test(x=myX,y=myY)$p.value)))
  mtext(side=3,line=1,cex=0.6,text=paste("P(MW)=",scientific_format(3)(wilcox.test(x=myX,y=myY)$p.value)))
  mtext(side=3,line=2,text="NGF",cex=0.7)
  mtext(side=1,line=3,text=colid,cex=0.7)
  legend("bottom",ncol=2,pch=15,col=c("white","grey"),leg=c(paste("without ",myName,sep=""),paste("with ",myName,sep="")),cex=0.7,bty="n")
  

  myX <- myOut.nt3[no.sites==0,match(colid,colnames(myOut.nt3))]
  myY <- myOut.nt3[no.sites>0,match(colid,colnames(myOut.nt3))] 
  boxplot(list(myX,myY),outline=F,frame=F,las=1,col=c("white","grey"),xaxt="n")
  mtext(side=3,line=0,cex=0.6,text=paste("P(t-test)=",scientific_format(3)(t.test(x=myX,y=myY)$p.value)))
  mtext(side=3,line=1,cex=0.6,text=paste("P(MW)=",scientific_format(3)(wilcox.test(x=myX,y=myY)$p.value)))
  mtext(side=3,line=2,text="NT3",cex=0.7)
  mtext(side=1,line=3,text=colid,cex=0.7)
  legend("bottom",ncol=2,pch=15,col=c("white","grey"),leg=c(paste("without ",myName,sep=""),paste("with ",myName,sep="")),cex=0.7,bty="n")
  

  
  myZ <-roll_sum(no.sites[sort(myOut.ngf[,match(colid,colnames(myOut.ngf))],decreasing=F,index.return=T)$ix],n=myN)
  plot(x=c(1:length(myZ)),y=myZ,col="blue",type="s",xaxt="n",las=1,frame=F,ylab="",ylim=c(0,max(myZ)))
  mtext(side=2,line=3,text="average no.motifs",cex=0.7)
  mtext(side=1,line=0,text=paste("sorted according to increasing ",colid,sep=""),cex=0.7)
  lines(apply(do.call(lapply(c(1:500),function(IX)return(roll_sum(no.sites[sort(myOut.ngf[sample(c(1:nrow(myOut.ngf))),match(colid,colnames(myOut.ngf))],decreasing=F,index.return=T)$ix],n=myN))),what=rbind),2,mean),col="red")
  
  myZ <-roll_sum(no.sites[sort(myOut.nt3[,match(colid,colnames(myOut.nt3))],decreasing=F,index.return=T)$ix],n=myN)
  plot(x=c(1:length(myZ)),y=myZ,col="blue",type="s",xaxt="n",las=1,frame=F,ylab="",ylim=c(0,max(myZ)))
  mtext(side=2,line=3,text="average no.motifs",cex=0.7)
  lines(apply(do.call(lapply(c(1:500),function(IX)return(roll_sum(no.sites[sort(myOut.nt3[sample(c(1:nrow(myOut.nt3))),match(colid,colnames(myOut.nt3))],decreasing=F,index.return=T)$ix],n=myN))),what=rbind),2,mean),col="red")
  mtext(side=1,line=0,text=paste("sorted according to increasing ",colid,sep=""),cex=0.7)
  

}


GetEnrichTransport <- function(no.sites=no.stau.sites,myvec=myVec,mymat.rdm=myMat.rdm){
  myout <- rep(NA,length(mymat.rdm))
  for(i in c(1:length(mymat.rdm))){
    my.val        <- sum(no.sites*myvec[[i]],na.rm=T)
    my.val.random <- apply(no.sites*mymat.rdm[[i]],2,sum)
    myout[i]      <- (my.val-mean(my.val.random))/sd(my.val.random)
  }
  return(myout)
}


#Create BINs where the distribution of axonal read counts is expected to be random i.e. does not depend on any covariate
CreateBins <- function(myout){
  bins.cb                        <- cut(myout$cb,breaks=seq(from=3,to=20,by=1),iclude.lowest=T)
  bins.txL                       <- cut(myout$txL,breaks=seq(from=2,to=4.5,by=0.25),iclude.lowest=T)
  bins.both                      <- paste(bins.cb,bins.txL,sep=".")
  temp                           <- data.frame(table(bins.both))
  tokeep                         <- temp$bins.both[temp$Freq>=5] 
  #tokeep                         <- temp$bins.both[temp$Freq>=1] 
  bins.both[!bins.both%in%tokeep]<-NA
  return(factor(as.character(bins.both)))
}



GetProbs <- function(mybins=bins.both.ngf,myfit=fit4.2.ngf,out=out.ngf){
  
  myOut <- do.call(lapply(c(1:length(levels(mybins))),function(IX){
    print(IX)
    myprediction         <- compute.prediction(IX,myfit=myfit,bins=mybins)
    subdat               <- out[which(mybins==levels(mybins)[IX]),]
    
    if(mean(myprediction)>0){
      fit.nb               <- fitdist(round(subdat$axons,digits=6)*10^6, "nbinom",method="mle",
                                      start=list(mu=max(c(round(mean(myprediction),digits=6)*10^6,0)),
                                                 size=round(var(myprediction),digits=6)*10^6))
    }
    else{
      fit.nb               <- fitdist(round(subdat$axons,digits=6)*10^6, "nbinom",method="mle")
    }
    fit.norm    <- fitdist(subdat$axons, "norm",start=list(mean=mean(myprediction),sd=sd(myprediction)))
    
    my.rdm.norm <- rnorm(n=10^4, mean=fit.norm$estimate["mean"], sd=fit.norm$estimate["sd"])
    my.rdm.nb   <- rnbinom(n=10^4, size=fit.nb$estimate["size"], mu=fit.nb$estimate["mu"])/10^6
    temp        <- do.call(what=rbind,lapply(subdat$axons,function(Z)return(c(sum(my.rdm.norm<=Z)/10^4,sum(my.rdm.norm>=Z)/10^4,
                                                                              sum(my.rdm.nb<=Z)/10^4,sum(my.rdm.nb>=Z)/10^4))))
    colnames(temp)<- c("prob(N).smaller","prob(N).bigger","prob(NB).smaller","prob(NB).bigger")
    temp          <-data.frame(temp,subdat)
    return(temp)
  }),what=rbind)
  
  
  p.value.norm <- apply(cbind(-log10(myOut$prob.N..smaller+min(myOut$prob.N..smaller[myOut$prob.N..smaller>0])),
                              -log10(myOut$prob.N..bigger+min(myOut$prob.N..bigger[myOut$prob.N..bigger>0]))),1,function(Z){
                                M<- max(Z)
                                if(which(Z==M)==1)return(-M)
                                else
                                  return(M)
                              })
  
  p.value.bn <- apply(cbind(-log10(myOut$prob.NB..smaller+min(myOut$prob.NB..smaller[myOut$prob.NB..smaller>0])),
                            -log10(myOut$prob.NB..bigger+min(myOut$prob.NB..bigger[myOut$prob.NB..bigger>0]))),1,function(Z){
                              M<- max(Z)
                              if(which(Z==M)==1)return(-M)
                              else
                                return(M)
                            })
  
  
  myOut <- data.frame(diff=p.value.norm-p.value.bn,p.value.norm,p.value.bn,myOut)
  return(myOut)
}



compute.prediction <- function(IX,myfit=fit4.2.ngf,bins=bins.both){
  temp        <- as.numeric(unlist(strsplit(gsub(gsub("\\].+?]","",levels(bins)[IX]),pattern="\\(",repl=""),split=",")))
  myrange.cb  <- seq(temp[1],temp[2],length=100)
  temp        <- as.numeric(unlist(strsplit(gsub(levels(bins)[IX],pattern="\\]",repl=""),split="[\\,\\(]"))[c(4,5)])
  myrange.txL <- seq(temp[1],temp[2],length=100)
  newData     <- data.frame(cb=rep(myrange.cb,100),txL=unlist(lapply(myrange.txL,function(Z)return(rep(Z,100)))))
  myprediction<- predict(myfit,newdata=newData)
  return(myprediction)
}

compute.posterior.mean <- function(IX,myfit){
  
  subdat               <- out[which(bins.both==levels(bins.both)[IX]),]
  myprediction         <- compute.prediction(IX,myfit)
  myrange.values.axons <- seq(from=0,to=22,length=1000)
  like                 <- dnorm(myrange.values.axons,mean=mean(subdat$axons),sd=sd(subdat$axons)) 
  #lines(myrange.values,like,col="blue")
  prior                <- dnorm(myrange.values, mean=mean(myprediction),sd=sd(myprediction))
  #lines(myrange.values,prior,col="green")
  post                 <- prior*like
  post                 <- post/sum(post)
  m                    <- sum(myrange.values.axons*post)
  s                    <- sqrt(sum(myrange.values.axons^2*post)-m^2)
  return(c(m,s))
}


require("plotly")
#Sys.setenv("plotly_username"="rluisier")
#Sys.setenv("plotly_api_key"="C2Xx7CqDQRgCVxeqLJSr")
#api_create(p, filename = "r-docs-midwest-boxplots")
getPlots <- function(subdat){
  p1 <- plot_ly(subdat, x = ~txL, y = ~cb, z = ~axons,marker = list(color = ~axons, colorscale = c('#FFE1A1', '#683531'), showscale = TRUE)) %>%
    add_markers() %>%
    layout(scene = list(xaxis = list(title = 'txL'),
                        yaxis = list(title = 'cb'),
                        zaxis = list(title = 'axons')),
           annotations = list(
             x = 1.13,
             y = 1.05,
             text = 'axons',
             xref = 'paper',
             yref = 'paper',
             showarrow = FALSE
           ))
  
  p2 <- plot_ly(subdat, x = ~cb, y = ~axons,marker = list(color = ~axons, colorscale = c('#FFE1A1', '#683531'), showscale = TRUE)) %>%
    add_markers() %>%
    layout(scene = list(xaxis = list(title = 'cb'),
                        yaxis = list(title = 'axons')),
           annotations = list(
             x = 1.13,
             y = 1.05,
             text = 'axons',
             xref = 'paper',
             yref = 'paper',
             showarrow = FALSE
           ))
  p3 <- plot_ly(subdat, x = ~txL, y = ~axons,marker = list(color = ~axons, colorscale = c('#FFE1A1', '#683531'), showscale = TRUE)) %>%
    add_markers() %>%
    layout(scene = list(xaxis = list(title = 'txL'),
                        yaxis = list(title = 'axons')),
           annotations = list(
             x = 1.13,
             y = 1.05,
             text = 'axons',
             xref = 'paper',
             yref = 'paper',
             showarrow = FALSE
           ))
  return(list(p1,p2,p3))
}



GetEnrichCLIPdata <- function(myclip=clip,tar.utr=all.utr.proxi.ngf,myMatRdm=list(myMat.ngf.dpud,myMat.ngf.drud),mycolid=c("diffPUD","dRUD"),myOutf=myOut.ngf){
  
  gOver                              <- findOverlaps(query=myclip,subject =tar.utr,ignore.strand=FALSE)
  no.stau.sites                      <- table(subjectHits(gOver))
  no.stau.sites                      <- no.stau.sites[match(as.character(c(1:length(tar.utr))),names(no.stau.sites))]
  no.stau.sites[is.na(no.stau.sites)]<-0
  names(no.stau.sites)               <- tar.utr$ID
  no.stau.sites                      <- as.vector(no.stau.sites)
  return(GetEnrichTransportFlex(no.sites=no.stau.sites,myOut=myOutf,myMat=myMatRdm,colid=mycolid))
}

#I may want to add the scores at each position instead of binary!
GetEnrichCLIPdataWithStrength <- function(myclip=clip,tar.utr=all.utr.proxi.ngf,myMatRdm=list(myMat.ngf.dpud,myMat.ngf.drud),mycolid=c("diffPUD","dRUD"),myOutf=myOut.ngf){
  
  gOver                              <- data.frame(findOverlaps(query=myclip,subject =tar.utr,ignore.strand=FALSE))
  no.stau.sites                      <- tapply(score(myclip)[gOver[,1]],INDEX=factor(paste(gOver[,2])),FUN=sum)
  no.stau.sites                      <- no.stau.sites[match(as.character(c(1:length(tar.utr))),names(no.stau.sites))]
  no.stau.sites[is.na(no.stau.sites)]<-0
  names(no.stau.sites)               <- tar.utr$ID
  no.stau.sites                      <- as.vector(no.stau.sites)
  
  return(GetEnrichTransportFlex(no.sites=no.stau.sites,myOut=myOutf,myMat=myMatRdm,colid=mycolid))
}



GetEnrichGOTransportAPA<-function(myOut=myOut.ngf,mysampleGO=mysampleGO.ngf,goID="GO:0051641",colid=c("diffPUD","dRUD"),myMat=list(myMat.ngf.dpud,myMat.ngf.drud)){
  go.entrez       <- genesInTerm(mysampleGO, goID)
  if(length(go.entrez)>0){
    mySize          <- length(go.entrez[[1]])
    go.entrez       <- go.entrez[[1]]#To extract all genes related to this term
    go.gs           <- unique(as.character(t2g$external_gene_name)[which(t2g$entrezgene%in%go.entrez)])
    
    is.in.go       <- as.numeric(myOut$geneSymbol%in%go.gs)
    
    if(sum(is.in.go)>20){
      my.val         <- sum(myOut[,match(colid[1],colnames(myOut))]*is.in.go)
      my.val.random  <- apply(is.in.go*myMat[[1]],2,sum)
      my.zval.dPUD   <- (my.val-mean(my.val.random))/sd(my.val.random)
      
      my.val         <- sum(myOut[,match(colid[2],colnames(myOut))]*is.in.go)
      my.val.random  <- apply(is.in.go*myMat[[2]],2,sum)
      my.zval.dRUD   <- (my.val-mean(my.val.random))/sd(my.val.random)
      
      return(c(my.zval.dPUD,my.zval.dRUD))
    }
    else{
      return(c(NA,NA))
    }
  }
  else{
    return(c(NA,NA))
  }
}


PlotEnrichNoSitesBINARY<- function(no.sites,colid=c("PUD.ngf.rel","PUD.nt3.rel","dPUD.rel"),myN=250,myName=colnames(mycount.all.tetra.big)[match(MOT,colnames(mycount.all.tetra.big))],YLAB="PUD.rel"){
  require(scales)
  par(mfrow=c(2,3),mar=c(3,4,4,3))
  myX <- myOut[no.sites==0,match(colid[1],colnames(myOut))]
  myY <- myOut[no.sites>0,match(colid[1],colnames(myOut))] 
  boxplot(list(myX,myY),outline=F,frame=F,las=1)
  mtext(side=3,line=0,cex=0.6,text=paste("P(t-test)=",scientific_format(3)(t.test(x=myX,y=myY)$p.value)))
  mtext(side=3,line=1,text="NGF",cex=0.7)
  mtext(side=2,line=3,text=YLAB)
  
  myX <- myOut[no.sites==0,match(colid[2],colnames(myOut))]
  myY <- myOut[no.sites>0,match(colid[2],colnames(myOut))] 
  boxplot(list(myX,myY),outline=F,frame=F,las=1)
  mtext(side=3,line=0,cex=0.6,text=paste("P(t-test)=",scientific_format(3)(t.test(x=myX,y=myY)$p.value)))
  mtext(side=3,line=1,text="NT3",cex=0.7)
  mtext(side=3,line=2,text=myName,cex=0.7)
  
  myX <- myOut[no.sites==0,match(colid[3],colnames(myOut))]
  myY <- myOut[no.sites>0,match(colid[3],colnames(myOut))] 
  boxplot(list(myX,myY),outline=F,frame=F,las=1)
  mtext(side=3,line=0,cex=0.6,text=paste("P(t-test)=",scientific_format(3)(t.test(x=myX,y=myY)$p.value)))
  mtext(side=3,line=1,text="NGF:NT3",cex=0.7)
  
  
  myZ <-roll_mean(no.sites[sort(myOut[,match(colid[1],colnames(myOut))],decreasing=T,index.return=T)$ix],n=myN)
  plot(x=c(1:length(myZ)),y=myZ,col="blue",type="s",xaxt="n",las=1,frame=F,ylab="",ylim=c(0,max(myZ)))
  mtext(side=2,line=3,text=YLAB,cex=0.7)
  lines(apply(do.call(lapply(c(1:500),function(IX)return(roll_mean(no.sites[sort(myOut[sample(c(1:nrow(myOut))),match(colid[1],colnames(myOut))],decreasing=T,index.return=T)$ix],n=myN))),what=rbind),2,mean),col="red")
  
  myZ <-roll_mean(no.sites[sort(myOut[,match(colid[2],colnames(myOut))],decreasing=T,index.return=T)$ix],n=myN)
  plot(x=c(1:length(myZ)),y=myZ,col="blue",type="s",xaxt="n",las=1,frame=F,ylab="",ylim=c(0,max(myZ)))
  mtext(side=2,line=3,text=YLAB,cex=0.7)
  lines(apply(do.call(lapply(c(1:500),function(IX)return(roll_mean(no.sites[sort(myOut[sample(c(1:nrow(myOut))),match(colid[2],colnames(myOut))],decreasing=T,index.return=T)$ix],n=myN))),what=rbind),2,mean),col="red")
  
  myZ <-roll_mean(no.sites[sort(myOut[,match(colid[3],colnames(myOut))],decreasing=T,index.return=T)$ix],n=myN)
  plot(x=c(1:length(myZ)),y=myZ,col="blue",type="s",xaxt="n",las=1,frame=F,ylab="",ylim=c(0,max(myZ)))
  mtext(side=2,line=3,text=YLAB,cex=0.7)
  lines(apply(do.call(lapply(c(1:500),function(IX)return(roll_mean(no.sites[sort(myOut[sample(c(1:nrow(myOut))),match(colid[3],colnames(myOut))],decreasing=T,index.return=T)$ix],n=myN))),what=rbind),2,mean),col="red")
  
}



GetEnrichGO <- function(goID="GO:0009966",myvec=myVec,myMat=myMat.rdm){
  go.entrez       <- genesInTerm(mysampleGO, goID)
  if(length(go.entrez)>0){
    mySize          <- length(go.entrez[[1]])
    go.entrez       <- go.entrez[[1]]#To extract all genes related to this term
    go.gs           <- unique(as.character(t2g$external_gene_name)[which(t2g$entrezgene%in%go.entrez)])
    
    is.in.go      <- as.numeric(names(myvec)%in%go.gs)
    my.val        <- sum(myvec*is.in.go)
    my.val.random <- apply(is.in.go*myMat,2,sum)
    return((my.val-mean(my.val.random))/sd(my.val.random))
  }
  else{
    return(NA)
  }
}

GetEnrichGOttest <- function(goID="GO:0009966",myvec=myVec){
  go.entrez       <- genesInTerm(mysampleGO, goID)
  if(length(go.entrez)>0){
    mySize          <- length(go.entrez[[1]])
    go.entrez       <- go.entrez[[1]]#To extract all genes related to this term
    go.gs           <- unique(as.character(t2g$external_gene_name)[which(t2g$entrezgene%in%go.entrez)])
    
    is.in.go      <- as.numeric(names(myvec)%in%go.gs)
    if(sum(is.in.go)>10){
    return(t.test(x=myvec[is.in.go==0],y=myvec[is.in.go==1])$p.value)
    }
    else{
      return(NA)
    }
  }
  else{
    return(NA)
  }
}

SelectPositive <- function(mat){
  mat[mat<0]<-NA
  return(mat)
}
SelectNegative <- function(mat){
  mat[mat>0]<-NA
  return(mat)
}

plotLineWithSD<-function(x=(-unlist(lapply(myranges,mean))[-1]),mymat=fisher_proxi_NGF_reg,col1="#81A4D6",col2=rgb(129/255,164/255,214/255,0.3),YLIM=c(-2,2),YLAB="Ip Fisher Enrichment",MAIN="Proximal"){
  mean_force=apply(mymat,2,function(W)return(median(W,na.rm=TRUE)))
  sd=apply(mymat,2,function(W)return(sd(W,na.rm=TRUE)))
  psd<-mean_force+sd
  nsd<-mean_force-sd
  plot(x, mean_force, ty="l", col=col1, ylab="", lty=1,lwd=3,las=1,frame=FALSE,ylim=YLIM,xlab="")
  lines(x, psd,col=col2)
  lines(x, nsd,col=col2)
  polygon(x=c(x, rev(x)), y=c(psd, rev(nsd)), col=col2,border=col2)
  lines(x, mean_force, col=col1,lwd=3)
  abline(h=0,col="black")
  mtext(side=3,line=0,text=MAIN,cex=0.7)
  mtext(side=1,line=2,text="distance from 3' end",cex=0.7)
  mtext(side=2,line=2,text=YLAB,cex=0.7)
}

#With PUD.rel
GetEnrichGOTransport <- function(goID="GO:0009966"){
  go.entrez       <- genesInTerm(mysampleGO, goID)
  if(length(go.entrez)>0){
    mySize          <- length(go.entrez[[1]])
    go.entrez       <- go.entrez[[1]]#To extract all genes related to this term
    go.gs           <- unique(as.character(t2g$external_gene_name)[which(t2g$entrezgene%in%go.entrez)])
    
    #Order this according to myOut$PUD.ngf.rel
    is.in.go      <- as.numeric(myOut$geneSymbol%in%go.gs)
    my.val        <- sum(myOut$PUD.ngf.rel*is.in.go)
    my.val.random <- apply(is.in.go*myMat.ngf.rel,2,sum)
    my.zval.ngf   <- (my.val-mean(my.val.random))/sd(my.val.random)
    
    my.val        <- sum(myOut$PUD.nt3.rel*is.in.go)
    my.val.random <- apply(is.in.go*myMat.nt3.rel,2,sum)
    my.zval.nt3   <- (my.val-mean(my.val.random))/sd(my.val.random)
    
    my.val        <- sum(myOut$dPUD.rel*is.in.go)
    my.val.random <- apply(is.in.go*myMat.diff.rel,2,sum)
    my.zval.diff   <- (my.val-mean(my.val.random))/sd(my.val.random)
    
    return(c(my.zval.ngf,my.zval.nt3,my.zval.diff))
  }
  else{
    return(c(NA,NA,NA))
  }
}

GetEnrichGOTransportFlex <- function(goID="GO:0009966",colid=c("RUD.ngf","RUD.nt3","dRUD"),myMat=list(myMat.ngf.rud,myMat.nt3,rud,myMat.dRUD)){
  go.entrez       <- genesInTerm(mysampleGO, goID)
  if(length(go.entrez)>0){
    mySize          <- length(go.entrez[[1]])
    go.entrez       <- go.entrez[[1]]#To extract all genes related to this term
    go.gs           <- unique(as.character(t2g$external_gene_name)[which(t2g$entrezgene%in%go.entrez)])
    
    #Order this according to myOut$PUD.ngf.rel
    is.in.go      <- as.numeric(myOut$geneSymbol%in%go.gs)
    my.val        <- sum(myOut[,match(colid[1],colnames(myOut))]*is.in.go)
    my.val.random <- apply(is.in.go*myMat[[1]],2,sum)
    my.zval.ngf   <- (my.val-mean(my.val.random))/sd(my.val.random)
    
    my.val        <- sum(myOut[,match(colid[2],colnames(myOut))]*is.in.go)
    my.val.random <- apply(is.in.go*myMat[[2]],2,sum)
    my.zval.nt3   <- (my.val-mean(my.val.random))/sd(my.val.random)
    
    my.val        <- sum(myOut[,match(colid[3],colnames(myOut))]*is.in.go)
    my.val.random <- apply(is.in.go*myMat[[3]],2,sum)
    my.zval.diff   <- (my.val-mean(my.val.random))/sd(my.val.random)
    
    return(c(my.zval.ngf,my.zval.nt3,my.zval.diff))
  }
  else{
    return(c(NA,NA,NA))
  }
}


#With PUD
GetEnrichGOTransport.init <- function(goID="GO:0009966"){
  go.entrez       <- genesInTerm(mysampleGO, goID)
  if(length(go.entrez)>0){
    mySize          <- length(go.entrez[[1]])
    go.entrez       <- go.entrez[[1]]#To extract all genes related to this term
    go.gs           <- unique(as.character(t2g$external_gene_name)[which(t2g$entrezgene%in%go.entrez)])
    
    #Order this according to myOut$PUD.ngf.rel
    is.in.go      <- as.numeric(myOut$geneSymbol%in%go.gs)
    my.val        <- sum(myOut$PUD.ngf*is.in.go)
    my.val.random <- apply(is.in.go*myMat.ngf,2,sum)
    my.zval.ngf   <- (my.val-mean(my.val.random))/sd(my.val.random)
    
    my.val        <- sum(myOut$PUD.nt3*is.in.go)
    my.val.random <- apply(is.in.go*myMat.nt3,2,sum)
    my.zval.nt3   <- (my.val-mean(my.val.random))/sd(my.val.random)
    
    my.val        <- sum(myOut$dPUD*is.in.go)
    my.val.random <- apply(is.in.go*myMat.diff,2,sum)
    my.zval.diff   <- (my.val-mean(my.val.random))/sd(my.val.random)
    
    return(c(my.zval.ngf,my.zval.nt3,my.zval.diff))
  }
  else{
    return(c(NA,NA,NA))
  }
}


PlotEnrichGO <- function(goID="GO:0008088",myvec=myVec[[1]],name="-log10(P-value) transport"){
  require(scales)
  go.entrez       <- genesInTerm(mysampleGO, goID)
  mySize          <- length(go.entrez[[1]])
  go.entrez       <- go.entrez[[1]]#To extract all genes related to this term
  go.gs           <- unique(as.character(t2g$external_gene_name)[which(t2g$entrezgene%in%go.entrez)])
  
  #Order this according to myOut$PUD.ngf.rel
  is.in.go      <- as.numeric(names(myvec)%in%go.gs)
  
  par(mfrow=c(1,3),mar=c(3,3,3,3))
  boxplot(list(myvec[!names(myvec)%in%go.gs],myvec[names(myvec)%in%go.gs]),outline=F,frame=F,las=1)
  mtext(side=3,line=0,cex=0.6,text=paste("P(t-test)=",scientific_format(3)(t.test(x=myvec[!names(myvec)%in%go.gs],y=myvec[names(myvec)%in%go.gs])$p.value)))
  mtext(side=2,line=2,text=name,cex=0.7)
  plot(is.in.go[sort(myvec,decreasing=T,index.return=T)$ix],col=rgb(0,0,1,0.8),type="h",lwd=0.5,xaxt="n",las=1,frame=F,xlab="",ylab="")
  mtext(side=3,line=1,text=paste(goterms[[goID]]," (",goID,")"),cex=0.7)
  l1 <- cumsum(is.in.go[sort(myvec,decreasing=T,index.return=T)$ix])
  l2 <- apply(replicate(500,cumsum(sample(is.in.go))),1,mean)
  plot(l1,type="l",lwd=0.5,xaxt="n",las=1,frame=F,xlab="")
  lines(l2,type="l",lwd=0.5,xaxt="n",las=1,frame=F,col=rgb(1,0,0))
}

PlotEnrichGOTransportBis <- function(goID="GO:0009966",YLAB="transport efficiency ",vec.ngf=myVec[[2]],vec.nt3=myVec[[6]],YLIM=c(-3,3)){
  require(scales)
  go.entrez       <- genesInTerm(mysampleGO, goID)
  mySize          <- length(go.entrez[[1]])
  go.entrez       <- go.entrez[[1]]#To extract all genes related to this term
  go.gs           <- unique(as.character(t2g$external_gene_name)[which(t2g$entrezgene%in%go.entrez)])
  is.in.go        <- as.numeric(names(vec.ngf)%in%go.gs)
  
  #par(mfrow=c(1,3),mar=c(3,3,3,3))
  myX <- vec.ngf[!names(vec.ngf)%in%go.gs]
  myY <- vec.ngf[names(vec.ngf)%in%go.gs]
  mylist<-list(myX,myY)
  names(mylist)<-c("not in GO","in GO")
  boxplot(mylist,outline=F,frame=F,las=1,ylim=YLIM,col="#81A4D6")
  mtext(side=3,line=0,cex=0.6,text=paste("P(t-test)=",scientific_format(3)(t.test(x=myX,y=myY)$p.value)))
  mtext(side=2,line=2,text=paste(YLAB,"NGF"),cex=0.7)

  myX <- vec.nt3[!names(vec.nt3)%in%go.gs]
  myY <- vec.nt3[names(vec.nt3)%in%go.gs]
  mylist<-list(myX,myY)
  names(mylist)<-c("not in GO","in GO")
  boxplot(mylist,outline=F,frame=F,las=1,ylim=YLIM,col="#AF72B0")
  mtext(side=3,line=0,cex=0.6,text=paste("P(t-test)=",scientific_format(3)(t.test(x=myX,y=myY)$p.value)))
  mtext(side=2,line=2,text=paste(YLAB,"NT3"),cex=0.7)
  mtext(side=3,line=1,text=paste(goterms[[goID]]," (",goID,")"),cex=0.7)
  

  diff <- vec.ngf-vec.nt3
  myX <- diff[!names(vec.nt3)%in%go.gs]
  myY <- diff[names(vec.nt3)%in%go.gs]
  mylist<-list(myX,myY)
  names(mylist)<-c("not in GO","in GO")
  boxplot(mylist,outline=F,frame=F,las=1,col="#C6C6C5")
  mtext(side=3,line=0,cex=0.6,text=paste("P(t-test)=",scientific_format(3)(t.test(x=myX,y=myY)$p.value)))  
  mtext(side=2,line=2,text=paste("diff",YLAB),cex=0.7)
  

}


PlotEnrichGOTransport <- function(goID="GO:0009966"){
  require(scales)
  go.entrez       <- genesInTerm(mysampleGO, goID)
  mySize          <- length(go.entrez[[1]])
  go.entrez       <- go.entrez[[1]]#To extract all genes related to this term
  go.gs           <- unique(as.character(t2g$external_gene_name)[which(t2g$entrezgene%in%go.entrez)])
  
  #Order this according to myOut$PUD.ngf.rel
  is.in.go      <- as.numeric(myOut$geneSymbol%in%go.gs)
  
  par(mfrow=c(3,3),mar=c(3,3,3,3))
  boxplot(list(myOut$PUD.ngf.rel[!myOut$geneSymbol%in%go.gs],myOut$PUD.ngf.rel[myOut$geneSymbol%in%go.gs]),outline=F,frame=F,las=1)
  mtext(side=3,line=0,cex=0.6,text=paste("P(t-test)=",scientific_format(3)(t.test(x=myOut$PUD.ngf.rel[myOut$geneSymbol%in%go.gs],y=myOut$PUD.ngf.rel[!myOut$geneSymbol%in%go.gs])$p.value)))
  plot(is.in.go[sort(myOut$PUD.ngf.rel,decreasing=T,index.return=T)$ix],col=rgb(0,0,1,0.01),type="h",lwd=0.5,xaxt="n",las=1,frame=F,xlab="",ylab="")
  mtext(side=3,line=1,text=paste(goterms[[goID]]," (",goID,")"),cex=0.7)
  mtext(side=3,line=0,text="NGF",cex=0.7)
  l1 <- cumsum(is.in.go[sort(myOut$PUD.ngf.rel,decreasing=T,index.return=T)$ix])
  l2 <- apply(replicate(500,cumsum(sample(is.in.go))),1,mean)
  plot(l1,type="l",lwd=0.5,xaxt="n",las=1,frame=F,xlab="")
  lines(l2,type="l",lwd=0.5,xaxt="n",las=1,frame=F,col=rgb(1,0,0))
  
  boxplot(list(myOut$PUD.nt3.rel[!myOut$geneSymbol%in%go.gs],myOut$PUD.nt3.rel[myOut$geneSymbol%in%go.gs]),outline=F,frame=F,las=1)
  mtext(side=3,line=0,cex=0.6,text=paste("P(t-test)=",scientific_format(3)(t.test(x=myOut$PUD.nt3.rel[myOut$geneSymbol%in%go.gs],y=myOut$PUD.nt3.rel[!myOut$geneSymbol%in%go.gs])$p.value)))
  plot(is.in.go[sort(myOut$PUD.nt3.rel,decreasing=T,index.return=T)$ix],col=rgb(0,0,1,0.01),type="h",lwd=0.5,xaxt="n",las=1,frame=F,xlab="",ylab="")
  mtext(side=3,line=0,text="NT3",cex=0.7)
  l1 <- cumsum(is.in.go[sort(myOut$PUD.nt3.rel,decreasing=T,index.return=T)$ix])
  l2 <- apply(replicate(500,cumsum(sample(is.in.go))),1,mean)
  plot(l1,type="l",lwd=0.5,xaxt="n",las=1,frame=F,xlab="")
  lines(l2,type="l",lwd=0.5,xaxt="n",las=1,frame=F,col=rgb(1,0,0))
  
  boxplot(list(myOut$dPUD.rel[!myOut$geneSymbol%in%go.gs],myOut$dPUD.rel[myOut$geneSymbol%in%go.gs]),outline=F,frame=F,las=1)
  mtext(side=3,line=0,cex=0.6,text=paste("P(t-test)=",scientific_format(3)(t.test(x=myOut$dPUD.rel[myOut$geneSymbol%in%go.gs],y=myOut$dPUD.rel[!myOut$geneSymbol%in%go.gs])$p.value)))
  plot(is.in.go[sort(myOut$dPUD.rel,decreasing=T,index.return=T)$ix],col=rgb(0,0,1,0.01),type="h",lwd=0.5,xaxt="n",las=1,frame=F,xlab="",ylab="")
  l1 <- cumsum(is.in.go[sort(myOut$dPUD.rel,decreasing=T,index.return=T)$ix])
  l2 <- apply(replicate(500,cumsum(sample(is.in.go))),1,mean)
  plot(l1,type="l",lwd=0.5,xaxt="n",las=1,frame=F,xlab="")
  lines(l2,type="l",lwd=0.5,xaxt="n",las=1,frame=F,col=rgb(1,0,0))
}

PlotEnrichGOTransportFlex <- function(myOut=myOut,goID="GO:0009966",mynames=c("NGF","NT3","NGF:NT3"),colid=c("RUD.ngf","RUD.nt3","dRUD"),mysampleGO=mysampleGO){
  require(scales)
  go.entrez       <- genesInTerm(mysampleGO, goID)
  go.entrez       <- go.entrez[[1]]#To extract all genes related to this term
  go.gs           <- unique(as.character(t2g$external_gene_name)[which(t2g$entrezgene%in%go.entrez)])
  mySize          <- length(go.gs)
  is.in.go        <- as.numeric(myOut$geneSymbol%in%go.gs)
  
  tempPlot <- function(myX,myY,myAll,Name="NGF"){
    boxplot(list(myX,myY),outline=F,frame=F,las=1)
    mtext(side=3,line=0,cex=0.6,text=paste("P(t-test)=",scientific_format(3)(t.test(x=myX,y=myY,var.equal = FALSE)$p.value)))
    plot(is.in.go[sort(myAll,decreasing=T,index.return=T)$ix],col=rgb(0,0,1,0.01),type="h",lwd=0.5,xaxt="n",las=1,frame=F,xlab="",ylab="")
    mtext(side=3,line=1,text=paste(goterms[[goID]]," (",goID,")"),cex=0.7)
    mtext(side=3,line=0,text=Name,cex=0.7)
    l1 <- cumsum(is.in.go[sort(myAll,decreasing=T,index.return=T)$ix])
    l2 <- apply(replicate(500,cumsum(sample(is.in.go))),1,mean)
    plot(l1,type="l",lwd=0.5,xaxt="n",las=1,frame=F,xlab="")
    lines(l2,type="l",lwd=0.5,xaxt="n",las=1,frame=F,col=rgb(1,0,0))
  }
  
  par(mfrow=c(length(mynames),3),mar=c(3,3,3,3))
  for(i in c(1:length(mynames))){
    tempPlot(myX=myOut[!myOut$geneSymbol%in%go.gs,match(colid[i],colnames(myOut))],
             myY=myOut[myOut$geneSymbol%in%go.gs,match(colid[i],colnames(myOut))],
             myAll=myOut[,match(colid[i],colnames(myOut))],
             Name=mynames[i])
  }
  mtext(side=1,line=1,text=paste("ntot genes=",length(unique(myOut$geneSymbol)),sep=""),cex=0.7)
  mtext(side=1,line=0,text=paste("# genes in GO=",sum(is.in.go),"( ntot.GO=",mySize,")",sep=""),cex=0.7)

}



GetGOI <- function(sampleGO=temp.sampleGO,mygoID=as.character(temp$GO.ID)[1]){
  entrez.oi       <- intersect(genesInTerm(sampleGO, mygoID)[[1]],sigGenes(sampleGO))
  gs.oi           <- unique(as.character(t2g$external_gene_name[match(entrez.oi,t2g$ensembl_transcript_id)]))
  return(paste(gs.oi, collapse = ', '))
}


GetEnrichTransportFlex <- function(no.sites=no.stau.sites,myOut=myOut.ngf,myMat=list(myMat.ngf,myMat.nt3,myMat.diff),colid=c("PUD.ngf","PUD.nt3","dPUD")){
  myout <- rep(NA,length(myMat))
  for(i in c(1:length(myMat))){
    my.val        <- sum(myOut[,match(colid[i],colnames(myOut))]*no.sites,na.rm=T)
    my.val.random <- apply(no.sites*myMat[[i]],2,sum)
    myout[i]      <- (my.val-mean(my.val.random))/sd(my.val.random)
  }
  return(myout)
}


myDrawVenn <- function(mylist=list(sel.gs.NT3[[1]],sel.gs.ngf[[1]]),mylabels=c("NT3","NGF"),myCols=c("green","blue")){
  myCounts<- c(length(setdiff(mylist[[1]],mylist[[2]])),length(intersect(mylist[[1]],mylist[[2]])),length(setdiff(mylist[[2]],mylist[[1]])))
  library(VennDiagram)
  grid.newpage()
  draw.pairwise.venn(area1 = length(mylist[[1]]), area2 = length(mylist[[2]]), cross.area = myCounts[2], category = mylabels,scaled=TRUE,fill=myCols,alpha=0.3,col=myCols)
  
}

GetGOI <- function(sampleGO=temp.sampleGO,mygoID=as.character(temp$GO.ID)[1]){
  entrez.oi       <- intersect(genesInTerm(sampleGO, mygoID)[[1]],sigGenes(sampleGO))
  gs.oi           <- unique(as.character(t2g$external_gene_name[match(entrez.oi,t2g$ensembl_transcript_id)]))
  return(paste(gs.oi, collapse = ', '))
}

ExtractRegion <- function(foi=nt3.utr,min=150,max=100){
  
  POS         <-  as.character(foi$strand)=="+"
  start       <-  as.numeric(foi$end)-min
  end         <-  as.numeric(foi$end)-max
  
  start[!POS] <-  as.numeric(foi$start[!POS])+max
  end[!POS]   <-  as.numeric(foi$start[!POS])+min
  foi$start   <- start
  foi$end     <- end
  print(sum(foi$end<=foi$start))
  return(foi)
}

ExtractRegionGR <- function(foi=all.utr,min=150,max=100){
  
  POS         <-  as.character(strand(foi))=="+"
  start       <-  as.numeric(end(foi))-min
  end         <-  as.numeric(end(foi))-max
  
  start[!POS] <-  as.numeric(start(foi)[!POS])+max
  end[!POS]   <-  as.numeric(start(foi)[!POS])+min

  ranges(foi)<-IRanges(start=start,end=end)
  
  return(foi)
}

CreateBedFiles <- function(foi,tarfile="~/Desktop/DataAnalsyisRiccio/Dec2016/axonal_remodelling/targets_for_catia.bed"){
  start        <- foi$start
  end          <- foi$end
  chr          <- as.character(foi$chr)
  strand       <- as.character(foi$strand)
  name         <- as.character(foi$ID)
  #start        <- as.numeric(format(start,scientific = FALSE))
  #end          <- as.numeric(format(end,scientific = FALSE))
  out          <- cbind(chr, start,end,name,score=rep(1,length(start)),strand)
  options(scipen=1000) 
  write.table(x=out,file=tarfile, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = FALSE,col.names = FALSE)
}




ExtractEnrichmentMotif <- function(set1=sel.uniqueID.ngf[[4]],set2=NA,mymatch=mycount.all.tetra.r1){
  
  myIX         <- match(set1,names(seq.all.r1))
  myIX         <- myIX[!is.na(myIX)]
  mymotifs     <- colnames(mymatch)
  mycount.seq1 <- apply(mymatch[myIX,]>0,2,sum)
  
  if(!is.na(set2)){
    myIX2        <- match(set2,names(seq.all.r1))
    myIX2        <- myIX2[!is.na(myIX2)]
    mycount.seq2 <- apply(mymatch[myIX2,]>0,2,sum)
    
    myp.value    <- unlist(lapply(c(1:length(mycount.seq1)),function(IX)return(fisher.test(rbind(c(mycount.seq1[IX],length(myIX)-mycount.seq1[IX]),
                                                                                                 c(mycount.seq2[IX],length(myIX2)-mycount.seq2[IX])))$p.value)))
    my.enrich    <- do.call(what=rbind,args=lapply(c(1:length(mycount.seq1)),function(IX)return(c(mycount.seq1[IX]/length(myIX)*length(myIX2)/mycount.seq2[IX],mycount.seq2[IX]/length(myIX2)*length(myIX)/mycount.seq1[IX]))))
    colnames(my.enrich)<- c("freq.set1/freq.set2","freq.set2/freq.set1")
    
    out        <- data.frame(count.seq1=mycount.seq1,freq.seq1=mycount.seq1/length(myIX),count.seq2=mycount.seq2,mymotifs,freq.sep2=mycount.seq2/length(myIX2),my.enrich,myp.value)
    return(out)
  }
  
  else{
    mycount.all  <- apply(mymatch>0,2,sum)
    IX           <- lapply(c(1:500),function(Z)return(sample(x=c(1:length(seq.all.r1)),size=length(myIX),replace=FALSE)))
    mycount.bg   <- do.call(what=rbind,args=lapply(IX,function(Z)return(apply(mymatch[Z,]>0,2,sum))))
    my.mean      <- apply(mycount.bg,2,mean)
    my.sd        <- apply(mycount.bg,2,sd)
    my.z         <- (mycount.seq1-my.mean)/my.sd
    names(my.z)  <- colnames(mymatch)
    outi         <- data.frame(Z=my.z,count.seq1=mycount.seq1,mean.bg=my.mean,sd.bg=my.sd)
    myp.value    <- unlist(lapply(c(1:length(mycount.seq1)),function(IX)return(fisher.test(rbind(c(mycount.seq1[IX],length(myIX)-mycount.seq1[IX]),
                                                                                                 c(mycount.all[IX],length(seq.all.r1)-mycount.all[IX])))$p.value)))
    my.enrich    <- do.call(what=rbind,args=lapply(c(1:length(mycount.seq1)),function(IX)return(c(mycount.seq1[IX]/length(myIX)*length(seq.all.r1)/mycount.all[IX],mycount.all[IX]/length(seq.all.r1)*length(myIX)/mycount.seq1[IX]))))
    colnames(my.enrich)<- c("freq.fg/freq.bg","freq.bg/freq.fg")
    IX          <-which(myp.value<=0.05)
    out        <- data.frame(count.seq1=mycount.seq1,freq.seq1=mycount.seq1/length(myIX),count.all=mycount.all,mymotifs,freq.all=mycount.all/length(seq.all.r1),my.enrich,myp.value,outi)
    return(out)
  }
}



PlotTrackUTR <- function(txOI="ENSRNOT00000036556",uniqueID=NA,genome="rn5",mytitle="test.pdf",GS=NA,LOG=FALSE,LOAD=FALSE){
  require(Gviz)
  require(GenomicRanges)
  # A. Load require files
  if(LOAD){load("/home/rluisier/data/Riccio/Exp_1/Dec2016/utrid/required_for_plotTrack_April_19.RData")}
  
  
  # B. Focus region based on txID or uniqueTXID (it may be important to plot overlapping features, especially if antisense)
  if(!is.na(txOI)){
    myGRTX      <- union(e75oi[myGeneModel$transcript==txOI],myGR[as.character(myGR$txID)==txOI])
    myRegion    <- e75oi[unique(subjectHits(findOverlaps(myGRTX,e75oi,ignore.strand=FALSE))),]#Select everything which overlaps with the txID
    uniqueID=NA
    print("txID known")
  }
  if(!is.na(uniqueID)){
    #uniqueID<- "ENSRNOT00000019820.3"
    txOI        <- unlist(strsplit(uniqueID,split="\\."))[1]
    myGRTX      <- union(e75oi[myGeneModel$transcript==txOI],myGR[as.character(myGR$txID)==uniqueID])
    myUTRoi     <- myGR[match(uniqueID,names(myGR)),]
    myRegion    <- e75oi[unique(subjectHits(findOverlaps(myGRTX,e75oi,ignore.strand=FALSE))),]#Select everything which overlaps with the txID
    #myRegion    <- e75oi[unique(subjectHits(findOverlaps(myUTRoi,e75oi,ignore.strand=FALSE))),]#Select everything which overlaps with the txID
    #start(myRegion) <- start(myRegion)-250
    #end(myRegion) <- end(myRegion)+250
    print("uniqueID known")
  }
  if(!is.na(GS)){
    print("GS known")
    txOI        <- unique(as.character(myGeneModel$transcript)[which(as.character(myGeneModel$symbol)==GS)])
    myGRTX      <- union(e75oi[myGeneModel$transcript==txOI],myGR[as.character(myGR$txID)==txOI])
    myRegion    <- e75oi[unique(subjectHits(findOverlaps(myGRTX,e75oi,ignore.strand=FALSE))),]
    uniqueID    <- NA
  }
  
  GS          <- unique(as.character(myGeneModel$symbol)[which(as.character(myGeneModel$transcript)==txOI)])
  
  # C. Create Tracks
  chr         <- as.character(unique(seqnames(myRegion)))
  # a) GeneModel
  tr1         <- as.data.frame(myRegion)
  tr1         <- tr1[,match(c("seqnames","start","end","width","strand","type","gene_id.","X.exon_number.","X.transcript_id.","X.gene_name."),colnames(tr1))]
  colnames(tr1)<- colnames(myGeneModel)
  tr1$exon    <- as.character(tr1$exon)
  tr1$exon[is.na(tr1$exon)]<-"utr"
  Tr1         <- GeneRegionTrack(tr1, genome=genome, chromosome=chr, name=GS,collapseTranscripts=FALSE,geneSymbols=TRUE,fill="black",transcriptAnnotation="symbol", background.title="lightgrey",cex.axis=0.8,cex.group=0.5,rot.title=0,cex.title=0.8,size=1.0,stackHeight=0.2,col.symbol="black",col.title="black",labelPos="above")
  print("GeneModel Track done")
  
  # b) De novo 3' UTR
  if(!is.na(uniqueID)){myUTRoi     <- myGR[match(uniqueID,names(myGR)),]}
  all.utr     <- myGR[as.character(myGR$txID)%in%txOI,]
  if(!is.na(uniqueID)){all.utr     <- all.utr[width(all.utr)<=width(myUTRoi),]}
  tr2        <- as.data.frame(all.utr)
  tr2        <- tr2[,c(1,2,3,4,5,11)]
  colnames(tr2)[1]<-"chromosome"
  colnames(tr2)[6]<-"transcript"
  if(tr2$strand[1]=="+"){
    tr2$start   <- rep(min(tr2$start),nrow(tr2))
    tr2         <- tr2[order(tr2$width),]
    tr2$transcript <- order(tr2$width)
  }
  
  if(tr2$strand[1]=="-"){
    tr2$end        <- rep(max(tr2$end),nrow(tr2))
    tr2            <- tr2[order(tr2$width),]
    tr2$transcript <- order(tr2$width)
  }
  
  
  Tr2         <- GeneRegionTrack(tr2, genome=genome, chromosome=chr, name="De novo 3'UTR",cex.axis=0.8,cex.group=0.5,rot.title=0,cex.title=0.8,size=0.1,stackHeight=0.5,transcriptAnnotation="transcript",col.symbol="black",col.title="black",fill="grey",showFeatureId=TRUE,labelPos="below",collapseTranscripts = FALSE)
  print("3'UTR Track done")
  
  # d) PolyA.db
  tr3             <- as.data.frame(intersect(myGRpaf,myRegion))
  colnames(tr3)[1]<-"chromosome"
  if(nrow(tr3)>0){
    Tr3             <- GeneRegionTrack(tr3, genome=genome, chromosome=chr, name="PA.db",fill="blue",background.title="white",cex.axis=0.8,cex.group=0.5,rot.title=0,cex.title=0.8,size=0.5,stackHeight=0.2,col.symbol="black",col.title="black",collapse=TRUE,collapseTranscripts = TRUE)
    IS.PA.DB=TRUE
  }
  else{
    IS.PA.DB=FALSE
  }
  print("PA.db Track done")
  
  
  # e) Axis track
  gtrack      <- GenomeAxisTrack(range=ranges(union(myRegion,myGR[as.character(myGR$txID)==txOI,])),name=chr)
  print("gtrack  done")
  
  # f) Data-tracks
  # Set the minimu and maximu ranges
  startC      <- min(c(tr1$start,tr2$start,tr3$start))
  endC        <- max(c(tr1$end,tr2$end,tr3$end))
  
  # Create data-tracks
  mydataTracks <- list()
  mynames      <- c("NGF.axon.1","NGF.axon.2","NGF.cb.2","NGF.cb.2","NT3.axon.1","NT3.axon.2","NT3.cb.2","Nt3.cb.2")
  mycols       <- c("chartreuse3","chartreuse3","chocolate1","chocolate1","blue","blue","lightblue","lightblue")
  if(as.character(tr2$strand)[1]=="+"){myCov<-myCovPos}
  else{myCov<-myCovNeg}
  
  for(IX in c(1:length(mynames))){
    focusCov            <-  which(start(myCov[[IX]])>=startC&end(myCov[[IX]])<=endC&as.character(seqnames(myCov[[IX]]))==chr)
    if(length(focusCov)>0){
      mydataTracks[[IX]]  <-  DataTrack(data=score(myCov[[IX]])[focusCov], start=start(myCov[[IX]])[focusCov], end=end(myCov[[IX]])[focusCov], chromosome=chr,
                                        genome=genome, name=mynames[IX],col.axis="black",col.symbol="black",col.title="black",cex.title=0.7,size=1.0,stackHeight=0.5,fill=mycols[IX],cex.axis=0.5)
    }
    else{
      lim    <- c(tr2$start[1], tr2$end[1])
      myL    <- tr2$end[1]-tr2$start[1]+1
      coords <- sort(c(lim[1], sample(seq(from = lim[1],to = lim[2]), myL-1), lim[2]))
      dat    <- rep(0,myL)
      mydataTracks[[IX]]  <-  DataTrack(data = dat, start = coords[-length(coords)], end = coords[-1], chromosome=chr,genome=genome, name=mynames[IX],
                                        col.axis="black",col.symbol="black",col.title="black",cex.title=0.7,size=1.0,stackHeight=0.5,fill=mycols[IX],cex.axis=0.5)
    }
    print(IX)
  }
  print("data tracks  done")
  
  # g) PAS track
  paoi       <- myPA[unique(subjectHits(findOverlaps(myGR[as.character(myGR$txID)%in%txOI,],myPA,ignore.strand=FALSE))),]#Select everything which overlaps with the 3' UTR
  if(length(paoi)>0){
    if(length(paoi)>10){
      paoi <- reduce(paoi[paoi$rep.names.myRes..IX...length.end..%in%"AATAAA",])
      Trpas             <- GeneRegionTrack(paoi, genome=genome, chromosome=chr, name="PAS reduced",fill="red",background.title="white",cex.axis=0.8,cex.group=0.5,rot.title=0,cex.title=0.8,size=0.5,stackHeight=0.2,col.symbol="black",col.title="black",collapse=TRUE)
    }
    else{
      Trpas             <- GeneRegionTrack(reduce(paoi), genome=genome, chromosome=chr, name="PAS",fill="red",background.title="white",cex.axis=0.8,cex.group=0.5,rot.title=0,cex.title=0.8,size=0.5,stackHeight=0.2,col.symbol="black",col.title="black",collapse = TRUE)
    }
    IS.PAS=TRUE
    print("there is an intersection with PAS")
  }
  else{
    IS.PAS=FALSE
  }
  
  
  
  # E. Set min and max for zooming on 3' UTR region
  if(as.character(tr2$strand)[1]=="+"){
    MIN  <- tr2$start[1]-500
    MAX  <- max(tr2$end)+500
  }
  if(as.character(tr2$strand)[1]=="-"){
    MIN  <- min(tr2$start)-500
    MAX  <- tr2$end[1]+500
  }
  
  WS1 <- max(1000,floor((endC-startC)/25))
  WS2 <- max(1000,floor((MAX-MIN)/10))
  
  
  # F. Plot
  if(!is.na(mytitle)){pdf(mytitle)}
  # Linear scale
  if((!IS.PA.DB)&(!IS.PAS)){
    myList<- list(mydataTracks[[3]],mydataTracks[[4]],mydataTracks[[1]],mydataTracks[[2]],mydataTracks[[7]],mydataTracks[[8]],mydataTracks[[5]],mydataTracks[[6]],Tr2,Tr1,gtrack)
  }
  if((!IS.PA.DB)&IS.PAS){
    myList<- list(mydataTracks[[3]],mydataTracks[[4]],mydataTracks[[1]],mydataTracks[[2]],mydataTracks[[7]],mydataTracks[[8]],mydataTracks[[5]],mydataTracks[[6]],Tr2,Tr1,gtrack,Trpas)
  }
  if(IS.PA.DB&(!IS.PAS)){
    myList<- list(mydataTracks[[3]],mydataTracks[[4]],mydataTracks[[1]],mydataTracks[[2]],mydataTracks[[7]],mydataTracks[[8]],mydataTracks[[5]],mydataTracks[[6]],Tr2,Tr1,gtrack,Tr3)
  }
  if(IS.PA.DB&IS.PAS){
    myList<- list(mydataTracks[[3]],mydataTracks[[4]],mydataTracks[[1]],mydataTracks[[2]],mydataTracks[[7]],mydataTracks[[8]],mydataTracks[[5]],mydataTracks[[6]],Tr2,Tr1,gtrack,Trpas,Tr3)
  }
  
  
  plotTracks(myList,type="hist",background.panel="white", background.title="white",add53=TRUE, add35=TRUE,window=WS1,fontcolor="black",thinBoxFeature="UTR")
  #Zooming
  plotTracks(myList,type="hist",background.panel="white", background.title="white",add53=TRUE, add35=TRUE,window=WS2,fontcolor="black",thinBoxFeature="UTR",from=MIN, to=MAX,littleTicks=TRUE,extend.left=1.2)
  
  if(LOG){
    #Log2 Scale
    plotTracks(myList,type="hist",background.panel="white", background.title="white",add53=TRUE, add35=TRUE,window=WS1,fontcolor="black",thinBoxFeature="UTR",transformation=function(x){x <- log2(x); x})
    #Zooming
    plotTracks(myList,type="hist",background.panel="white", background.title="white",add53=TRUE, add35=TRUE,window=WS2,fontcolor="black",thinBoxFeature="UTR",from=MIN, to=MAX,littleTicks=TRUE,transformation=function(x){x <- log2(x); x},extend.left=-1.2)
  }
  if(!is.na(mytitle)){dev.off()}
  
  rm(Tr1,Tr2,Tr3,txOI, myGRTX, myRegion,paoi,Trpas,mydataTracks,IS.PA.DB,IS.PAS,uniqueID,GS)
  
}


PlotEnrichGOCompareNGFNT3 <- function(goID="GO:0032784",alpha=1.0){
  require(scales)
  
  #NGF
  go.entrez           <- genesInTerm(mysampleGO.ngf, goID)
  go.entrez           <- go.entrez[[1]]#To extract all genes related to this term
  go.gs               <- unique(as.character(t2g$external_gene_name)[which(t2g$entrezgene%in%go.entrez)])
  mySize              <- length(go.gs)
  is.in.go.ngf        <- as.numeric(myOut.ngf$geneSymbol%in%go.gs)
  
  #NT3
  go.entrez           <- genesInTerm(mysampleGO.nt3, goID)
  go.entrez           <- go.entrez[[1]]#To extract all genes related to this term
  go.gs               <- unique(as.character(t2g$external_gene_name)[which(t2g$entrezgene%in%go.entrez)])
  mySize              <- length(go.gs)
  is.in.go.nt3        <- as.numeric(myOut.NT3$geneSymbol%in%go.gs)
  
  
  myVals <- list(myX.ngf= myOut.ngf$dPUD[is.in.go.ngf==0],
                 myY.ngf=myOut.ngf$dPUD[is.in.go.ngf==1],
                 myX.nt3=myOut.NT3$dPUD[is.in.go.nt3==0],
                 myY.nt3=myOut.NT3$dPUD[is.in.go.nt3==1],
                 myAll.ngf=myOut.ngf$dPUD,
                 myAll.nt3=myOut.NT3$dPUD)
  colsNGF <- c("#81A4D6","#2D598E","#083872")
  colsNT3 <- c("#AE73B1","#79387C","#57055B")
  
  layout(matrix(c(1,1,6,7,2,3,4,5),ncol=4,nrow=2,byrow = FALSE))
  par(mar=c(4,4,4,1))
  boxplot(myVals[-c(5,6)],col=c("white",colsNGF[1],"white",colsNT3[1]),las=1,frame=FALSE,outline=FALSE,xaxt="n")
  
  pvals <- scientific_format(3)(c(t.test(x=myVals[[1]],y=myVals[[2]],var.equal = FALSE)$p.value,
                                  t.test(x=myVals[[2]],y=myVals[[4]],var.equal = FALSE)$p.value,
                                  t.test(x=myVals[[3]],y=myVals[[4]],var.equal = FALSE)$p.value))
  
  mtext(side=3,line=0,cex=0.6,text=paste("P(t-test NGF)=",pvals[1]))
  mtext(side=3,line=1,cex=0.6,text=paste("P(t-test NT3)=",pvals[3]))
  mtext(side=3,line=2,cex=0.6,text=paste("P(t-test NGF:NT3)=",pvals[2]))
  mtext(side=1,line=0,at=c(1,2,3,4),cex=0.6,text=c("not in","in GO","not in","in GO"))
  mtext(side=1,line=1,at=c(1,2,3,4),cex=0.6,text=c(sum(is.in.go.ngf==0),sum(is.in.go.ngf),sum(is.in.go.nt3==0),sum(is.in.go.nt3)))
  mtext(side=2,line=3,cex=0.6,text="dPUD")
  
  
  plot(is.in.go.ngf[sort(myVals[[5]],decreasing=FALSE,index.return=T)$ix],col=colsNGF[1],type="h",lwd=0.5,xaxt="n",las=1,frame=F,xlab="",ylab="")
  mtext(side=3,line=3,text=goterms[[goID]],cex=0.7,font=2)
  mtext(side=3,line=2,text=goID,cex=0.7)
  mtext(side=2,line=2,text="is in GO",cex=0.7)
  
  plot(is.in.go.nt3[sort(myVals[[6]],decreasing=FALSE,index.return=T)$ix],col=colsNT3[1],type="h",lwd=0.5,xaxt="n",las=1,frame=F,xlab="",ylab="")
  mtext(side=2,line=2,text="is in GO",cex=0.7)
  mtext(side=1,line=2,text="increasing dPUD",cex=0.7)
  
  l1 <- cumsum(is.in.go.ngf[sort(myVals[[5]],decreasing=TRUE,index.return=T)$ix])
  l2 <- apply(replicate(500,cumsum(sample(is.in.go.ngf))),1,mean)
  plot(l1,type="l",lwd=0.5,xaxt="n",las=1,xlab="",col=colsNGF[1])
  lines(l2,type="l",lwd=0.5,xaxt="n",las=1,frame=F,col=rgb(1,0,0))
  mtext(side=2,line=2,text="cumulative sum of the number in GO",cex=0.7)
  mtext(side=1,line=2,text="decreasing dPUD",cex=0.7)
  
  
  l1 <- cumsum(is.in.go.nt3[sort(myVals[[6]],decreasing=TRUE,index.return=T)$ix])
  l2 <- apply(replicate(500,cumsum(sample(is.in.go.nt3))),1,mean)
  plot(l1,type="l",lwd=0.5,xaxt="n",las=1,xlab="",col=colsNT3[1])
  lines(l2,type="l",lwd=0.5,xaxt="n",las=1,frame=F,col=rgb(1,0,0))
  mtext(side=2,line=2,text="cumulative sum of the number in GO",cex=0.7)
  mtext(side=1,line=2,text="decreasing dPUD",cex=0.7)
  
  #Focus on the origin of the difference
  
  
  
  
  #par(mar=c(4,4,4,1),mfrow=c(2,2))
  
  myVals <- list(myX.ngf= myOut.ngf$mPUD.cb[is.in.go.ngf==0],
                 myY.ngf=myOut.ngf$mPUD.cb[is.in.go.ngf==1],
                 myX.nt3=myOut.ngf$mPUD.axon[is.in.go.ngf==0],
                 myY.nt3=myOut.ngf$mPUD.axon[is.in.go.ngf==1],
                 myAll.ngf=myOut.ngf$dPUD)
  boxplot(myVals[-c(5)],col=c("white",colsNGF[1],"white",colsNGF[2]),las=1,frame=FALSE,outline=FALSE,xaxt="n",ylim=c(0,1))
  
  
  pvals <- scientific_format(3)(c(t.test(x=myVals[[1]],y=myVals[[2]],var.equal = FALSE)$p.value,
                                  t.test(x=myVals[[3]],y=myVals[[4]],var.equal = FALSE)$p.value,
                                  t.test(x=myVals[[2]],y=myVals[[4]],var.equal = FALSE)$p.value))
  
  mtext(side=3,line=2,cex=0.6,text=paste("P(t-test CB)=",pvals[1]))
  mtext(side=3,line=1,cex=0.6,text=paste("P(t-test Axons)=",pvals[2]))
  mtext(side=3,line=0,cex=0.6,text=paste("P(t-test Axons:CB)=",pvals[3]))
  
  mtext(side=1,line=0,at=c(1,2,3,4),cex=0.6,text=c("not in","in GO","not in","in GO"))
  mtext(side=1,line=1,at=c(1,2,3,4),cex=0.6,text=c(sum(is.in.go.ngf==0),sum(is.in.go.ngf),sum(is.in.go.ngf==0),sum(is.in.go.ngf)))
  mtext(side=1,line=2,at=c(1.5,3.5),cex=0.6,text=c("cell body","axons"))
  mtext(side=2,line=3,cex=0.6,text="PUD")
  # mtext(side=3,line=3,text=paste(goterms[[goID]],goID),cex=0.7,font=2)
  
  myVals <- list(myX.ngf= myOut.NT3$mPUD.cb[is.in.go.nt3==0],
                 myY.ngf=myOut.NT3$mPUD.cb[is.in.go.nt3==1],
                 myX.nt3=myOut.NT3$mPUD.axon[is.in.go.nt3==0],
                 myY.nt3=myOut.NT3$mPUD.axon[is.in.go.nt3==1])
  boxplot(myVals[-c(5)],col=c("white",colsNT3[1],"white",colsNT3[2]),las=1,frame=FALSE,outline=FALSE,xaxt="n",ylim=c(0,1))
  
  
  pvals <- scientific_format(3)(c(t.test(x=myVals[[1]],y=myVals[[2]],var.equal = FALSE)$p.value,
                                  t.test(x=myVals[[3]],y=myVals[[4]],var.equal = FALSE)$p.value,
                                  t.test(x=myVals[[2]],y=myVals[[4]],var.equal = FALSE)$p.value))
  
  mtext(side=3,line=2,cex=0.6,text=paste("P(t-test CB)=",pvals[1]))
  mtext(side=3,line=1,cex=0.6,text=paste("P(t-test Axons)=",pvals[2]))
  mtext(side=3,line=0,cex=0.6,text=paste("P(t-test Axons:CB)=",pvals[3]))
  
  mtext(side=1,line=0,at=c(1,2,3,4),cex=0.6,text=c("not in","in GO","not in","in GO"))
  mtext(side=1,line=1,at=c(1,2,3,4),cex=0.6,text=c(sum(is.in.go.nt3==0),sum(is.in.go.nt3),sum(is.in.go.nt3==0),sum(is.in.go.nt3)))
  mtext(side=1,line=2,at=c(1.5,3.5),cex=0.6,text=c("cell body","axons"))
  mtext(side=2,line=3,cex=0.6,text="PUD")
  # mtext(side=3,line=3,text=paste(goterms[[goID]],goID),cex=0.7,font=2)
}



GetEnrichClipAxonalRemodelling <- function(clip=clipdat[,1],myres=candidate_axonal_remodeling[[1]],myout=myOut.ngf){
  
  
  
  #Proximal -- only use the proximal pool of 3' UTR (where we expect to find most sites)
  is_bound              <- factor(ifelse(clip[which(rownames(clipdat)%in%unique(myout$imIDp))]>0,"bound","unbound"))
  if(length(levels(is_bound))==1){
    fisher.proxi    <- NA
    welch.proxi.pud <- NA
    welch.proxi.rud <- NA
  }
  else{
    is_pool.proxi         <- factor(ifelse(rownames(clipdat[which(rownames(clipdat)%in%unique(myout$imIDp)),])%in%myres$imIDp,"positive","none"))#be careful because you here need clupdat
    p_bound_1             <- sum(is_bound=="bound",na.rm=TRUE)/sum(!is.na(is_bound))
    p_bound_proxi         <- sum(is_bound=="bound"&is_pool.proxi=="positive",na.rm=TRUE)/sum(is_pool.proxi=="positive")
    
    
    p_pool.proxi          <- sum(is_pool.proxi=="positive",na.rm=TRUE)/sum(!is.na(is_pool.proxi))
    p_expected.proxi      <- p_bound_1*p_pool.proxi
    p_observed.proxi      <- sum(is_bound=="bound"&is_pool.proxi=="positive",na.rm=TRUE)/sum(!is.na(is_bound))
    fisher.proxi          <- -log10(fisher.test(x=is_bound,y=is_pool.proxi)$p.value)*sign(p_observed.proxi-p_expected.proxi)
    
    if(length(which(myres$imIDp%in%rownames(clipdat)[clip==1]))>=3){
      welch.proxi.pud       <- -log10(t.test(x=myres$mPUD.axon[which(myres$imIDp%in%rownames(clipdat)[clip==0])],
                                             y=myres$mPUD.axon[which(myres$imIDp%in%rownames(clipdat)[clip>0])],var.equal = FALSE)$p.value)*sign(p_observed.proxi-p_expected.proxi)
      welch.proxi.rud       <- -log10(t.test(x=myres$mRUD.axon[which(myres$imIDp%in%rownames(clipdat)[clip==0])],
                                             y=myres$mRUD.axon[which(myres$imIDp%in%rownames(clipdat)[clip>0])],var.equal = FALSE)$p.value)*sign(p_observed.proxi-p_expected.proxi)
    }
    else{
      welch.proxi.pud<-0
      welch.proxi.rud<-0
    }
  }
  
  
  #Distal -- only use the proximal pool of 3' UTR (where we expect to find most sites)
  is_bound              <- factor(ifelse(clip[which(rownames(clipdat)%in%unique(myout$uniqueID))]>0,"bound","unbound"))
  if(length(levels(is_bound))==1){
    fisher.dist    <- NA
    welch.dist.pud <- NA
    welch.dist.rud <- NA
  }
  else{
    
    is_pool.dist          <- factor(ifelse(rownames(clipdat[which(rownames(clipdat)%in%unique(myout$uniqueID)),])%in%myres$uniqueID,"positive","none"))
    p_bound_2             <- sum(is_bound=="bound",na.rm=TRUE)/sum(!is.na(is_bound))
    p_bound_dist          <- sum(is_bound=="bound"&is_pool.dist=="positive",na.rm=TRUE)/sum(is_pool.dist=="positive")
    
    p_pool.dist           <- sum(is_pool.dist=="positive",na.rm=TRUE)/sum(!is.na(is_pool.dist))
    p_expected.dist       <- p_bound_2*p_pool.dist
    p_observed.proxi      <- sum(is_bound=="bound"&is_pool.dist=="positive",na.rm=TRUE)/sum(!is.na(is_bound))
    fisher.dist          <- -log10(fisher.test(x=is_bound,y=is_pool.dist)$p.value)*sign(p_observed.proxi-p_expected.dist)
    if(length(which(myres$uniqueID%in%rownames(clipdat)[clip==1]))>=3){
      
      welch.dist.pud       <- -log10(t.test(x=myres$mPUD.axon[which(myres$uniqueID%in%rownames(clipdat)[clip==0])],
                                            y=myres$mPUD.axon[which(myres$uniqueID%in%rownames(clipdat)[clip>0])],var.equal = FALSE)$p.value)*sign(p_observed.proxi-p_expected.dist)
      welch.dist.rud       <- -log10(t.test(x=myres$mRUD.axon[which(myres$uniqueID%in%rownames(clipdat)[clip==0])],
                                            y=myres$mRUD.axon[which(myres$uniqueID%in%rownames(clipdat)[clip>0])],var.equal = FALSE)$p.value)*sign(p_observed.proxi-p_expected.dist)
    }
    else{
      welch.dist.pud<-0
      welch.dist.rud<-0
      
    }
  }
  
  out<-c(p_bound_1,p_bound_proxi,p_bound_2,p_bound_dist,fisher.proxi,welch.proxi.pud,welch.proxi.rud,fisher.dist,welch.dist.pud,welch.dist.rud)  
  names(out)<- c("frac_bg_proxi","frac_proxi","frac_bg_dist","frac_dist","fisher.proxi","welch.proxi.pud","welch.proxi.rud","fisher.dist","welch.dist.pud","welch.dist.rud")
  return(out)
}

plotLineWithSD<-function(x=(-unlist(lapply(myranges,mean))[-1]),mymat=fisher_proxi_NGF_reg,col1="#81A4D6",col2=rgb(129/255,164/255,214/255,0.3),YLIM=c(-2,2),YLAB="Ip Fisher Enrichment",MAIN="Proximal"){
  mean_force=apply(mymat,2,median)
  sd=apply(mymat,2,sd)
  psd<-mean_force+sd
  nsd<-mean_force-sd
  plot(x, mean_force, ty="l", col=col1, ylab="", lty=1,lwd=3,las=1,frame=FALSE,ylim=YLIM,xlab="")
  lines(x, psd,col=col2)
  lines(x, nsd,col=col2)
  polygon(x=c(x, rev(x)), y=c(psd, rev(nsd)), col=col2,border=col2)
  lines(x, mean_force, col=col1,lwd=3)
  abline(h=0,col="black")
  mtext(side=3,line=0,text=MAIN,cex=0.7)
  mtext(side=1,line=2,text="distance from 3' end",cex=0.7)
  mtext(side=2,line=2,text=YLAB,cex=0.7)
}

plotEnrichAlongUTRGeneral <- function(mots="CPSF160_iClip",sel=c(20:30),cond="NGF",scaling=TRUE,mytests=list( ttest_APA_NGF_Ip,ttest_APA_NGF_Id)){
  
  
  myX  = -unlist(lapply(myranges,mean))[sel]
  myIp = mytests[[1]][match(mots,rownames(mytests[[1]])),sel]
  myId = -mytests[[2]][match(mots,rownames(mytests[[2]])),sel]
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
  mtext(side=3,line=0,text=cond,font=1,cex=0.5)
  mtext(side=1,line=2,text="distance from 3' end",font=1,cex=0.5)
}

plotEnrichAlongUTRSimple <- function(mots="TDP43_iClip",sel=c(15:30),scaling=TRUE,mytests=list(ttest_transport_NGF,ttest_transport_NT3),mytestname="Welsch test",YLAB="-log10(P-value)"){
  
  
  myX  = -unlist(lapply(myranges,mean))[sel]
  myIp = mytests[[1]][match(mots,rownames(mytests[[1]])),sel]
  myId = mytests[[2]][match(mots,rownames(mytests[[2]])),sel]
  if(scaling){
    myIp <- myIp/max(abs(myIp))
    myId <- myId/max(abs(myId))
  }
  plot(myX,myIp,col="white",las=1,ylab=YLAB,frame=FALSE,xlab="distante to 3'end",ylim=c(min(c(myIp,myId),na.rm=TRUE),max(c(myIp,myId),na.rm=TRUE)))
  abline(h=0,col="black")
  abline(v=0,col="red",lty=2)
  lines(myX,myIp,col=colsNGF[1],lwd=1.5)
  lines(myX,myId,col=colsNT3[1],lwd=1.5)
  points(myX,myIp,col=colsNGF[1],pch=19)
  points(myX,myId,col=colsNT3[1],pch=19)
  mtext(side=3,line=0,text=mots,font=1,cex=0.5)
  mtext(side=2,line=2,text=mytestname,font=1,cex=0.5)
  mtext(side=1,line=2,text="distance from  Ip 3' end",font=1,cex=0.5)
  
}

GetEnrichClipBinary <- function(clip=no_reads[,i],val_vec=myVal[[j]],val_mat=myMat[[j]]){
  is_bound <- factor(ifelse(clip>1,"bound","unbound"))
  is_pool  <- factor(ifelse(val_vec==1,"positive","none"))
  
  p_bound   <- sum(is_bound=="bound",na.rm=TRUE)/sum(!is.na(is_bound))
  p_pool    <- sum(is_pool=="positive",na.rm=TRUE)/sum(!is.na(is_pool))
  p_expected<- p_bound*p_pool
  p_observed<- sum(is_bound=="bound"&is_pool=="positive",na.rm=TRUE)/sum(!is.na(is_pool))
  
  my.pval       <- fisher.test(x=is_bound,y=is_pool)$p.value
  my.val        <- sum(clip*val_vec,na.rm=TRUE)
  my.val.random <- apply(clip*val_mat,2,function(Z)return(sum(Z,na.rm=TRUE)))
  my.zscore     <- (my.val-mean(my.val.random))/sd(my.val.random)
  #my.random     <- sum(apply(val_mat,2,function(W)return(fisher.test(x=factor(ifelse(clip>0,"bound","unbound")),y=factor(ifelse(W==1,"positive","none")))$p.value))<my.pval)/1000
  return(c(enr=p_observed/p_expected,p_e=p_expected,p_o=p_observed,myz=my.zscore,P_Fisher=my.pval)
         #return(c(enr=p_observed/p_expected,p_e=p_expected,p_o=p_observed,myz=my.zscore,P_Fisher=my.pval,P_perm=my.random)
  )
}

EnrichedInPool <- function(IX=1,pool=list(as.character(diff.transport.neurotrophins$uniqueID)[myVal[[1]]==1],as.character(diff.transport.neurotrophins$uniqueID)[myVal[[3]]==1])){
  
  sites   <- clipdat[match(diff.transport.neurotrophins$uniqueID,rownames(clipdat)),IX]
  
  frac_bg <- sum(sites>0,na.rm=TRUE)/sum(!is.na(sites))
  frac_NGF<- sum(sites[match(pool[[1]],diff.transport.neurotrophins$uniqueID)]>0,na.rm=TRUE)/length(pool[[1]])
  frac_NT3<- sum(sites[match(pool[[2]],diff.transport.neurotrophins$uniqueID)]>0,na.rm=TRUE)/length(pool[[2]])
  pFngf <- -log10(fisher.test(x=factor(ifelse(sites>0,"bound","unbound")),y=factor(ifelse(diff.transport.neurotrophins$uniqueID%in%pool[[1]],"inpool","not")))$p.value)*sign(frac_NGF-frac_bg)
  pFnt3 <- -log10(fisher.test(x=factor(ifelse(sites>0,"bound","unbound")),y=factor(ifelse(diff.transport.neurotrophins$uniqueID%in%pool[[2]],"inpool","not")))$p.value)*sign(frac_NT3-frac_bg)
  
  
  tE <- list(
    diff.transport.neurotrophins$transport.ngf[which(sites==0&diff.transport.neurotrophins$uniqueID%in%pool[[1]])],
    diff.transport.neurotrophins$transport.ngf[which(sites>0&diff.transport.neurotrophins$uniqueID%in%pool[[1]])],
    diff.transport.neurotrophins$transport.nt3[which(sites==0&diff.transport.neurotrophins$uniqueID%in%pool[[2]])],
    diff.transport.neurotrophins$transport.nt3[which(sites>0&diff.transport.neurotrophins$uniqueID%in%pool[[2]])]
  )
  
  if(min(unlist(lapply(tE[c(1,2)],length)))>10){
    pWngf <- -log10(t.test(tE[[1]],tE[[2]],var.equal=FALSE)$p.value)*sign(mean(tE[[2]],na.rm=TRUE)-mean(tE[[1]],na.rm=TRUE))
  }else{
    pWngf<-NA
  }
  
  if(min(unlist(lapply(tE[c(3,4)],length)))>10){
    pWnt3 <- -log10(t.test(tE[[3]],tE[[4]],var.equal=FALSE)$p.value)*sign(mean(tE[[4]],na.rm=TRUE)-mean(tE[[3]],na.rm=TRUE))
  }else{
    pWnt3<-NA
  }
  
  
  
  return(c(frac_bg,frac_NGF,frac_NT3,pFngf,pFnt3,pWngf,pWnt3))  
  
}

plotEnrichAlongUTRSimple <- function(mots="TDP43_iClip",sel=c(15:30),scaling=TRUE,mytests=list(ttest_transport_NGF,ttest_transport_NT3),mytestname="Welsch test",YLAB="-log10(P-value)",LWD=1.5){
  
  
  myX  = -unlist(lapply(myranges,mean))[sel]
  myIp = mytests[[1]][match(mots,rownames(mytests[[1]])),sel]
  myId = mytests[[2]][match(mots,rownames(mytests[[2]])),sel]
  if(scaling){
    myIp <- myIp/max(abs(myIp))
    myId <- myId/max(abs(myId))
  }
  plot(myX,myIp,col="white",las=1,ylab=YLAB,frame=FALSE,xlab="distante to 3'end",ylim=c(min(c(myIp,myId),na.rm=TRUE),max(c(myIp,myId),na.rm=TRUE)))
  abline(h=0,col="black")
  abline(v=0,col="red",lty=2)
  lines(myX,myIp,col=colsNGF[1],lwd=LWD,type="h")
  lines(myX+15,myId,col=colsNT3[1],lwd=LWD,type="h")
  mtext(side=3,line=0,text=mots,font=1,cex=0.5)
  mtext(side=2,line=2,text=mytestname,font=1,cex=0.5)
  mtext(side=1,line=2,text="distance from  Ip 3' end",font=1,cex=0.5)
  
  plot(myX,myIp,col="white",las=1,ylab=YLAB,frame=FALSE,xlab="distante to 3'end",ylim=c(min(c(myIp,myId),na.rm=TRUE),max(c(myIp,myId),na.rm=TRUE)))
  abline(h=0,col="black")
  abline(v=0,col="red",lty=2)
  lines(myX,myIp,col=colsNGF[1],lwd=LWD)
  lines(myX,myId,col=colsNT3[1],lwd=LWD)
  points(myX,myIp,col=colsNGF[1],pch=19,cex=1.0)
  points(myX,myId,col=colsNT3[1],pch=19,cex=1.0)
  mtext(side=3,line=0,text=mots,font=1,cex=0.5)
  mtext(side=2,line=2,text=mytestname,font=1,cex=0.5)
  mtext(side=1,line=2,text="distance from  Ip 3' end",font=1,cex=0.5)
  
}

plotEnrichAlongUTRGeneralNGF_NT3 <- function(mots="CPSF160_iClip",sel=c(20:30),cond="NGF",scaling=TRUE,mytests=list( list(ttest_APA_NGF_Ip,ttest_APA_NT3_Ip),list(ttest_APA_NGF_Id,ttest_APA_NT3_Id)),mytestname="Fisher P-value"){
  
  
  myX  = -unlist(lapply(myranges,mean))[sel]
  myIp = mytests[[1]][[1]][match(mots,rownames(mytests[[1]][[1]])),sel]
  myId = mytests[[1]][[2]][match(mots,rownames(mytests[[1]][[2]])),sel]
  if(scaling){
    myIp <- myIp/max(abs(myIp))
    myId <- myId/max(abs(myId))
  }
  plot(myX,myIp,col="white",las=1,ylab="-log10(P-value)",frame=FALSE,xlab="distante to 3'end",ylim=c(min(c(myIp,myId),na.rm=TRUE),max(c(myIp,myId),na.rm=TRUE)))
  abline(h=0,col="black")
  abline(v=0,col="red",lty=2)
  lines(myX,myIp,col=colsNGF[1],lwd=1.5)
  lines(myX,myId,col=colsNT3[1],lwd=1.5)
  points(myX,myIp,col=colsNGF[1],pch=19)
  points(myX,myId,col=colsNT3[1],pch=19)
  mtext(side=3,line=0,text=mots,font=1,cex=0.5)
  mtext(side=2,line=2,text=mytestname,font=1,cex=0.5)
  mtext(side=1,line=2,text="distance from  Ip 3' end",font=1,cex=0.5)
  
  myIp = -mytests[[2]][[1]][match(mots,rownames(mytests[[2]][[1]])),sel]
  myId = -mytests[[2]][[2]][match(mots,rownames(mytests[[2]][[2]])),sel]
  if(scaling){
    myIp <- myIp/max(abs(myIp))
    myId <- myId/max(abs(myId))
  }
  plot(myX,myIp,col="white",las=1,ylab="-log10(P-value)",frame=FALSE,xlab="distante to 3'end",ylim=c(min(c(myIp,myId),na.rm=TRUE),max(c(myIp,myId),na.rm=TRUE)))
  abline(h=0,col="black")
  abline(v=0,col="red",lty=2)
  lines(myX,myIp,col=colsNGF[1],lwd=1.5)
  lines(myX,myId,col=colsNT3[1],lwd=1.5)
  points(myX,myIp,col=colsNGF[1],pch=19)
  points(myX,myId,col=colsNT3[1],pch=19)
  mtext(side=3,line=0,text=mots,font=1,cex=0.5)
  mtext(side=2,line=2,text=mytestname,font=1,cex=0.5)
  mtext(side=1,line=2,text="distance from  Id 3' end",font=1,cex=0.5)
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


PlotFractionWithSitesImproved <- function(mymot="HuR_iClip",myregion="[0:-50]",myRes=list(NGF=res.sel.ngf[[3]],NT3=res.sel.NT3[[3]])){
  
  clipdat   <- myClips[[match(myregion,mynames_rages)]]

  #
  #NGF Ip
  #
  myres=myRes[[1]]
  myout=myOut.ngf[!duplicated(myOut.ngf$imIDp),]
  sel=match(myout$imIDp,rownames(clipdat))
  sites=clipdat[sel[!is.na(sel)],match(mymot,colnames(clipdat))]
  myout=myout[!is.na(sel),]
  is_remodelled               <- factor(ifelse(myout$imIDp%in%myres$imIDp,"remodelled","none"))
  is_bound                    <- factor(ifelse(sites>0,"bound","unbound"))
  fisher_ngf                  <- scientific_format(2)(fisher.test(x=is_bound,y=is_remodelled)$p.value)
  frac_bg_ngf                 <- sum(sites>0)/length(sites)
  frac_fg_ngf                 <- sum(sites>0&myout$imIDp%in%myres$imIDp)/nrow(myres)
  #
  #NT3 Ip
  #
  myres=myRes[[2]]
  myout=myOut.NT3[!duplicated(myOut.NT3$imIDp),]
  sel=match(myout$imIDp,rownames(clipdat))
  sites=clipdat[sel[!is.na(sel)],match(mymot,colnames(clipdat))]
  myout=myout[!is.na(sel),]
  
  is_remodelled               <- factor(ifelse(myout$imIDp%in%myres$imIDp,"remodelled","none"))
  is_bound                    <- factor(ifelse(sites>0,"bound","unbound"))
  fisher_nt3                  <- scientific_format(2)(fisher.test(x=is_bound,y=is_remodelled)$p.value)
  frac_bg_nt3                 <- sum(sites>0)/length(sites)
  frac_fg_nt3                 <- sum(sites>0&myout$imIDp%in%myres$imIDp)/nrow(myres)
  
  myvals                      <-matrix(c(frac_bg_ngf,frac_fg_ngf,frac_bg_nt3,frac_fg_nt3),nrow=2,ncol=2,byrow=TRUE)
  rownames(myvals)<- c("NGF","NT3")
  mp=barplot(t(myvals),beside=TRUE,col=c("grey",colsNGF[1],"grey",colsNT3[1]),las=1,frame=FALSE,ylim=c(0,max(myvals+0.2)))
  mtext(side=2,line=3,text="fraction of Ip with peaks",cex=0.7)
  mtext(side=3,line=1,text=paste(mymot,"--",myregion),cex=0.7,font=2)
  mtext(side=3,line=0,at=c(apply(mp,2,mean)),text=paste("P=",c(fisher_ngf,fisher_nt3)),cex=0.6,font=1)
  

  #
  #NGF Id
  #
  myres=myRes[[1]]
  myout=myOut.ngf
  sel=match(myout$uniqueID,rownames(clipdat))
  sites=clipdat[sel[!is.na(sel)],match(mymot,colnames(clipdat))]
  myout=myout[!is.na(sel),]
  is_remodelled               <- factor(ifelse(myout$uniqueID%in%myres$uniqueID,"remodelled","none"))
  is_bound                    <- factor(ifelse(sites>0,"bound","unbound"))
  fisher_ngf                  <- scientific_format(2)(fisher.test(x=is_bound,y=is_remodelled)$p.value)
  frac_bg_ngf                 <- sum(sites>0)/length(sites)
  frac_fg_ngf                 <- sum(sites>0&myout$uniqueID%in%myres$uniqueID)/nrow(myres)
  #
  #NT3 Ip
  #
  myres=myRes[[2]]
  myout=myOut.NT3
  sel=match(myout$uniqueID,rownames(clipdat))
  sites=clipdat[sel[!is.na(sel)],match(mymot,colnames(clipdat))]
  myout=myout[!is.na(sel),]
  
  is_remodelled               <- factor(ifelse(myout$uniqueID%in%myres$uniqueID,"remodelled","none"))
  is_bound                    <- factor(ifelse(sites>0,"bound","unbound"))
  fisher_nt3                  <- scientific_format(2)(fisher.test(x=is_bound,y=is_remodelled)$p.value)
  frac_bg_nt3                 <- sum(sites>0)/length(sites)
  frac_fg_nt3                 <- sum(sites>0&myout$uniqueID%in%myres$uniqueID)/nrow(myres)
  
  myvals                      <-matrix(c(frac_bg_ngf,frac_fg_ngf,frac_bg_nt3,frac_fg_nt3),nrow=2,ncol=2,byrow=TRUE)
  rownames(myvals)<- c("NGF","NT3")
  mp=barplot(t(myvals),beside=TRUE,col=c("grey",colsNGF[1],"grey",colsNT3[1]),las=1,frame=FALSE,ylim=c(0,max(myvals+0.2)))
  mtext(side=2,line=3,text="fraction of Id with peaks",cex=0.7)
  mtext(side=3,line=1,text=paste(mymot,"--",myregion),cex=0.7,font=2)
  mtext(side=3,line=0,at=c(apply(mp,2,mean)),text=paste("P=",c(fisher_ngf,fisher_nt3)),cex=0.6,font=1)
  
}



PlotFractionWithSites <- function(mymot="HuR_iClip",myregion="[0:-50]",myRes=list(NGF=res.sel.ngf[[3]],NT3=res.sel.NT3[[3]])){
  
  myOut     <- list(myOut.ngf,myOut.NT3)
  clipdat   <- myClips[[match(myregion,mynames_rages)]]
  clip      <- clipdat[,match(mymot,colnames(clipdat))]

  
  #Proximal -- only use the proximal pool of 3' UTR (where we expect to find most sites)
  is_bound              <- lapply(myOut,function(myout)return(factor(ifelse(clip[which(rownames(clipdat)%in%unique(myout$imIDp))]>0,"bound","unbound"))))#identify the ones in the pool of proximal 3' UTR that are bound by the RBPs
  is_pool.proxi         <- mapply(myout=myOut,myres=myRes,function(myout,myres)return(factor(ifelse(rownames(clipdat[which(rownames(clipdat)%in%unique(myout$imIDp)),])%in%myres$imIDp,"positive","none"))))#Identify in the pool of proximal those which are candidate of axonal remodeling
  p_bound_1             <- lapply(is_bound,function(test)return(sum(test=="bound",na.rm=TRUE)/sum(!is.na(test))))#fraction of proximal that are bound
  p_bound_proxi         <- mapply(test1=is_bound,test2=is_pool.proxi,function(test1,test2)return(sum(test1=="bound"&test2=="positive",na.rm=TRUE)/sum(test2=="positive")))#fraction of the candidate of axonal remodelling that are bound
  p_pool.proxi          <- lapply(is_pool.proxi,function(test)return(sum(test=="positive",na.rm=TRUE)/sum(!is.na(test))))#fraction which of proximal which are candidate of axonal remodelling
  p_expected.proxi      <- mapply(A=p_bound_1,B=p_pool.proxi,function(A,B)return(A*B))
  p_observed.proxi      <- mapply(A=is_bound,B=is_pool.proxi,function(A,B)return(sum(A=="bound"&B=="positive",na.rm=TRUE)/sum(!is.na(A))))
  fisher.proxi          <- mapply(A=is_bound,B=is_pool.proxi,function(A,B)return(scientific_format(3)(fisher.test(x=A,y=B)$p.value)))
   myvals          <-matrix(c(unlist(p_bound_1),unlist(p_bound_proxi)),nrow=2,ncol=2,byrow=FALSE)
   rownames(myvals)<- c("NGF","NT3")
   mp=barplot(t(myvals),beside=TRUE,col=c("grey",colsNGF[1],"grey",colsNT3[1]),las=1,frame=FALSE,ylim=c(0,max(myvals+0.2)))
  mtext(side=2,line=3,text="fraction of Ip with peaks",cex=0.7)
  mtext(side=3,line=1,text=paste(mymot,"--",myregion),cex=0.7,font=2)
  mtext(side=3,line=0,at=c(apply(mp,2,mean)),text=paste("P=",unlist(fisher.proxi)),cex=0.6,font=1)
  
  #Proximal -- only use the proximal pool of 3' UTR (where we expect to find most sites)
  is_bound              <- lapply(myOut,function(myout)return(factor(ifelse(clip[which(rownames(clipdat)%in%unique(myout$uniqueID))]>0,"bound","unbound"))))#identify the ones in the pool of proximal 3' UTR that are bound by the RBPs
  is_pool.proxi         <- mapply(myout=myOut,myres=myRes,function(myout,myres)return(factor(ifelse(rownames(clipdat[which(rownames(clipdat)%in%unique(myout$uniqueID)),])%in%myres$uniqueID,"positive","none"))))#Identify in the pool of proximal those which are candidate of axonal remodeling
  p_bound_1             <- lapply(is_bound,function(test)return(sum(test=="bound",na.rm=TRUE)/sum(!is.na(test))))#fraction of proximal that are bound
  p_bound_proxi         <- mapply(test1=is_bound,test2=is_pool.proxi,function(test1,test2)return(sum(test1=="bound"&test2=="positive",na.rm=TRUE)/sum(test2=="positive")))#fraction of the candidate of axonal remodelling that are bound
  p_pool.proxi          <- lapply(is_pool.proxi,function(test)return(sum(test=="positive",na.rm=TRUE)/sum(!is.na(test))))#fraction which of proximal which are candidate of axonal remodelling
  p_expected.proxi      <- mapply(A=p_bound_1,B=p_pool.proxi,function(A,B)return(A*B))
  p_observed.proxi      <- mapply(A=is_bound,B=is_pool.proxi,function(A,B)return(sum(A=="bound"&B=="positive",na.rm=TRUE)/sum(!is.na(A))))
  fisher.proxi          <- mapply(A=is_bound,B=is_pool.proxi,function(A,B)return(scientific_format(3)(fisher.test(x=A,y=B)$p.value)))
  myvals          <-matrix(c(unlist(p_bound_1),unlist(p_bound_proxi)),nrow=2,ncol=2,byrow=FALSE)
  rownames(myvals)<- c("NGF","NT3")
  mp=barplot(t(myvals),beside=TRUE,col=c("grey",colsNGF[1],"grey",colsNT3[1]),las=1,frame=FALSE,ylim=c(0,max(myvals+0.2)))
  mtext(side=2,line=3,text="fraction of Id with peaks",cex=0.7)
  mtext(side=3,line=1,text=paste(mymot,"--",myregion),cex=0.7,font=2)
  mtext(side=3,line=0,at=c(apply(mp,2,mean)),text=paste("P=",unlist(fisher.proxi)),cex=0.6,font=1)
  
}

PlotFractionWithSitesRegions <- function(mymot="HuR_iClip",myregions=c("[250:200]","[200:150]","[150:100]","[100:50]","[50:0]","[0:-50]","[-50:-100]","[-100:-150]"),myRes=list(NGF=res.sel.ngf[[3]],NT3=res.sel.NT3[[3]])){
  
  myOut     <- list(myOut.ngf,myOut.NT3)
  clipdat   <- myClips[match(myregions,mynames_rages)]
  clip      <- do.call(what=cbind,lapply(clipdat,function(Z)return(Z[,match(mymot,colnames(Z))])))
  rownames(clip)<- rownames(clipdat[[1]])
  
  pvals_Ip          <- matrix(ncol=ncol(clip),nrow=2,0)
  fracs_Ip          <- matrix(ncol=ncol(clip),nrow=4,0)
  colnames(fracs_Ip)<- myregions
  rownames(fracs_Ip) <-c("ct_NGF","NGF","ct_NT3","NT3") 
  colnames(pvals_Ip)<- myregions
  rownames(pvals_Ip) <-c("NGF","NT3") 
  
  for(IX in c(1:ncol(clip))){
  #Proximal -- only use the proximal pool of 3' UTR (where we expect to find most sites)
  is_bound              <- lapply(myOut,function(myout)return(factor(ifelse(clip[which(rownames(clip)%in%unique(myout$imIDp)),IX]>0,"bound","unbound"))))#identify the ones in the pool of proximal 3' UTR that are bound by the RBPs
  is_pool.proxi         <- mapply(myout=myOut,myres=myRes,function(myout,myres)return(factor(ifelse(rownames(clip[which(rownames(clip)%in%unique(myout$imIDp)),])%in%myres$imIDp,"positive","none"))))#Identify in the pool of proximal those which are candidate of axonal remodeling
  p_bound_1             <- lapply(is_bound,function(test)return(sum(test=="bound",na.rm=TRUE)/sum(!is.na(test))))#fraction of proximal that are bound
  p_bound_proxi         <- mapply(test1=is_bound,test2=is_pool.proxi,function(test1,test2)return(sum(test1=="bound"&test2=="positive",na.rm=TRUE)/sum(test2=="positive")))#fraction of the candidate of axonal remodelling that are bound
  p_pool.proxi          <- lapply(is_pool.proxi,function(test)return(sum(test=="positive",na.rm=TRUE)/sum(!is.na(test))))#fraction which of proximal which are candidate of axonal remodelling
  p_expected.proxi      <- mapply(A=p_bound_1,B=p_pool.proxi,function(A,B)return(A*B))
  p_observed.proxi      <- mapply(A=is_bound,B=is_pool.proxi,function(A,B)return(sum(A=="bound"&B=="positive",na.rm=TRUE)/sum(!is.na(A))))
  fracs_Ip[,IX]  <- as.vector(matrix(c(unlist(p_bound_1),unlist(p_bound_proxi)),nrow=2,ncol=2,byrow=TRUE))
  pvals_Ip[,IX]  <- mapply(A=is_bound,B=is_pool.proxi,function(A,B)return(scientific_format(3)(fisher.test(x=A,y=B)$p.value)))
  }
  YLIM=c(0,max(fracs_Ip+0.1))
  fracs_Ip_m <- rbind(apply(fracs_Ip[c(1,3),],2,mean),fracs_Ip[c(2,4),])
  mp=barplot(fracs_Ip_m,beside=TRUE,col=c("grey",colsNGF[1],colsNT3[1]),las=1,frame=FALSE,ylim=YLIM,xaxt="n")
  mtext(side=2,line=3,text="fraction of Ip with peaks",cex=0.7)
  mtext(side=3,line=2,text=paste(mymot,"--",colnames(clip)[IX]),cex=0.7,font=2)
  mtext(side=3,line=1,at=c(apply(mp,2,mean)),text=paste("P=",pvals_Ip[1,]),cex=0.5,font=1)
  mtext(side=3,line=0,at=c(apply(mp,2,mean)),text=paste("P=",pvals_Ip[2,]),cex=0.5,font=1)
  mtext(side=1,line=0,at=c(apply(mp,2,mean)),text=myregions,cex=0.5)
  
  
  
  pvals_Ip          <- matrix(ncol=ncol(clip),nrow=2,0)
  fracs_Ip          <- matrix(ncol=ncol(clip),nrow=4,0)
  colnames(fracs_Ip)<- myregions
  rownames(fracs_Ip) <-c("ct_NGF","NGF","ct_NT3","NT3") 
  colnames(pvals_Ip)<- myregions
  rownames(pvals_Ip) <-c("NGF","NT3") 
  
  for(IX in c(1:ncol(clip))){
    #Proximal -- only use the proximal pool of 3' UTR (where we expect to find most sites)
    is_bound              <- lapply(myOut,function(myout)return(factor(ifelse(clip[which(rownames(clip)%in%unique(myout$uniqueID)),IX]>0,"bound","unbound"))))#identify the ones in the pool of proximal 3' UTR that are bound by the RBPs
    is_pool.proxi         <- mapply(myout=myOut,myres=myRes,function(myout,myres)return(factor(ifelse(rownames(clip[which(rownames(clip)%in%unique(myout$uniqueID)),])%in%myres$uniqueID,"positive","none"))))#Identify in the pool of proximal those which are candidate of axonal remodeling
    p_bound_1             <- lapply(is_bound,function(test)return(sum(test=="bound",na.rm=TRUE)/sum(!is.na(test))))#fraction of proximal that are bound
    p_bound_proxi         <- mapply(test1=is_bound,test2=is_pool.proxi,function(test1,test2)return(sum(test1=="bound"&test2=="positive",na.rm=TRUE)/sum(test2=="positive")))#fraction of the candidate of axonal remodelling that are bound
    p_pool.proxi          <- lapply(is_pool.proxi,function(test)return(sum(test=="positive",na.rm=TRUE)/sum(!is.na(test))))#fraction which of proximal which are candidate of axonal remodelling
    p_expected.proxi      <- mapply(A=p_bound_1,B=p_pool.proxi,function(A,B)return(A*B))
    p_observed.proxi      <- mapply(A=is_bound,B=is_pool.proxi,function(A,B)return(sum(A=="bound"&B=="positive",na.rm=TRUE)/sum(!is.na(A))))
    fracs_Ip[,IX]  <- as.vector(matrix(c(unlist(p_bound_1),unlist(p_bound_proxi)),nrow=2,ncol=2,byrow=TRUE))
    pvals_Ip[,IX]  <- mapply(A=is_bound,B=is_pool.proxi,function(A,B)return(scientific_format(3)(fisher.test(x=A,y=B)$p.value)))
  }
  
  fracs_Ip_m <- rbind(apply(fracs_Ip[c(1,3),],2,mean),fracs_Ip[c(2,4),])
  mp=barplot(fracs_Ip_m,beside=TRUE,col=c("grey",colsNGF[1],colsNT3[1]),las=1,frame=FALSE,ylim=YLIM,xaxt="n")
  mtext(side=2,line=3,text="fraction of Id with peaks",cex=0.7)
  mtext(side=3,line=2,text=paste(mymot,"--",colnames(clip)[IX]),cex=0.7,font=2)
  mtext(side=3,line=1,at=c(apply(mp,2,mean)),text=paste("P=",pvals_Ip[1,]),cex=0.5,font=1)
  mtext(side=3,line=0,at=c(apply(mp,2,mean)),text=paste("P=",pvals_Ip[2,]),cex=0.5,font=1)
  mtext(side=1,line=0,at=c(apply(mp,2,mean)),text=myregions,cex=0.5)
  
}

AnalysisRegAxonalRemodelling <- function(tempot){
  par(mfrow=c(3,2))
  #Welch test on the PUD values -- need to plot both Welch in Axons and Welch in CB to give an understanding of the position-dependent positive association of the peaks with PUD
  plotEnrichAlongUTRGeneralNGF_NT3(mots=tempot,sel=c(2:30),cond="NGF",scaling=FALSE,mytests=list( list(ttest_APA_CB_NGF_Ip,ttest_APA_CB_NT3_Ip),list(ttest_APA_CB_NGF_Id,ttest_APA_CB_NT3_Id)),mytestname = "Welch test PUD CB")
  plotEnrichAlongUTRGeneralNGF_NT3(mots=tempot,sel=c(2:30),cond="NGF",scaling=FALSE,mytests=list( list(ttest_APA_NGF_Ip,ttest_APA_NT3_Ip),list(ttest_APA_NGF_Id,ttest_APA_NT3_Id)),mytestname = "Welch test PUD AX")
  plotEnrichAlongUTRGeneralNGF_NT3(mots=tempot,sel=c(2:30),cond="NGF",scaling=FALSE,mytests=list( list(fisher_proxi_NGF_reg,fisher_proxi_NT3_reg),list(fisher_dist_NGF_reg,fisher_dist_NT3_reg)),mytestname = "Fisher test")
  par(mfrow=c(2,1))
  PlotFractionWithSitesRegions(mymot=tempot,myregions=c("[250:200]","[200:150]","[150:100]","[100:50]","[50:0]","[0:-50]","[-50:-100]","[-100:-150]"),myRes=list(NGF=res.sel.ngf[[3]],NT3=res.sel.NT3[[3]]))
}



PlotHeatmapRegulators <- function(mots=unique(unlist(global_positive)),mymats=list(ttest_APA_NGF_Ip[,-1],-ttest_APA_NGF_Id[,-1]),scaling=TRUE){
  
  if(scaling){
    mymats[[1]]<- mymats[[1]]/max(abs(mymats[[1]]),na.rm=TRUE)
    mymats[[2]]<- mymats[[2]]/max(abs(mymats[[2]]),na.rm=TRUE)
  }
  
  mymats                <- lapply(mymats,function(Z)return(Z[match(mots,rownames(Z)),]))
  GS                    <- unlist(lapply(as.character(mots),function(Z)return((unlist(strsplit(Z,split="_"))[1]))))
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



PlotEffectMotifsAlongUTR <- function(mymot=pos_regs[[1]][1],cliplist=myClipOI,myN=500,YMAX=0.25){
  #Specific for the transport
  myrangesOI  <- c("[200:150]","[150:100]","[100:50]")
  require(scales)
  for(IX in c(1:length(myClipOI))){
    no.sites <- myClipOI[[IX]][,match(mymot,colnames(myClipOI[[IX]]))]
    par(mfrow=c(2,3),mar=c(4,4,4,2))
    myX <- diff.transport.neurotrophins$transport.ngf[no.sites==0]
    myY <- diff.transport.neurotrophins$transport.ngf[no.sites>0] 
    boxplot(list(unbound=myX,bound=myY),outline=F,frame=F,las=1,col="#81A4D6",ylim=c(-3,3))
    mtext(side=3,line=0,cex=0.6,text=paste("P(Welsh)=",scientific_format(3)(t.test(x=myX,y=myY,var.equal = FALSE)$p.value)))
    mtext(side=3,line=1,cex=0.6,text=paste("P(MW)=",scientific_format(3)(wilcox.test(x=myX,y=myY)$p.value)))
    mtext(side=3,line=2,text=paste("NGF",mymot),cex=0.7)
    mtext(side=2,line=3,text="transport efficiency",cex=0.7)
    abline(h=0)
    
    myX <- diff.transport.neurotrophins$transport.nt3[no.sites==0]
    myY <- diff.transport.neurotrophins$transport.nt3[no.sites>0] 
    boxplot(list(unbound=myX,bound=myY),outline=F,frame=F,las=1,col="#AE72B0",ylim=c(-3,3))
    mtext(side=3,line=0,cex=0.6,text=paste("P(Welsh)=",scientific_format(3)(t.test(x=myX,y=myY,var.equal = FALSE)$p.value)))
    mtext(side=3,line=1,cex=0.6,text=paste("P(MW)=",scientific_format(3)(wilcox.test(x=myX,y=myY)$p.value)))
    mtext(side=3,line=2,text=paste("NT3",mymot),cex=0.7)
    mtext(side=2,line=3,text="transport efficiency",cex=0.7)
    abline(h=0)
    
    myX <- diff.transport.neurotrophins$diff.transport[no.sites==0]
    myY <- diff.transport.neurotrophins$diff.transport[no.sites>0] 
    boxplot(list(unbound=myX,bound=myY),outline=F,frame=F,las=1,col="#C7C5C5",ylim=c(-3,3))
    mtext(side=3,line=0,cex=0.6,text=paste("P(Welsh)=",scientific_format(3)(t.test(x=myX,y=myY,var.equal = FALSE)$p.value)))
    mtext(side=3,line=1,cex=0.6,text=paste("P(MW)=",scientific_format(3)(wilcox.test(x=myX,y=myY)$p.value)))
    mtext(side=3,line=2,text=paste("NGF:NT3",mymot),cex=0.7)
    mtext(side=2,line=3,text="transport efficiency diff",cex=0.7)
    abline(h=0)
    
    require("roll")
    par(mar=c(4,4,1,2))
    trans_NGF<-  (no.sites>0)[sort(diff.transport.neurotrophins$transport.ngf,decreasing=T,index.return=T)$ix]
    myZ      <-  roll_mean(trans_NGF,width=myN)
    plot(x=c(1:length(myZ)),y=rev(myZ),type="s",las=1,frame=F,ylab="",ylim=c(0,YMAX),xlab="",col="#81A4D6")
    mtext(side=2,line=3,text="average fraction with sites",cex=0.7)
    mtext(side=1,line=3,text="transport efficiency",cex=0.7)
    abline(h=0)
    trans_NT3<-  (no.sites>0)[sort(diff.transport.neurotrophins$transport.nt3,decreasing=T,index.return=T)$ix]
    myZ      <-  roll_mean(trans_NT3,width=myN)
    plot(x=c(1:length(myZ)),y=rev(myZ),type="s",las=1,frame=F,ylab="",ylim=c(0,YMAX),xlab="",col="#AE72B0")
    mtext(side=2,line=3,text="average fraction with sites",cex=0.7)
    mtext(side=1,line=3,text="transport efficiency",cex=0.7)
    abline(h=0)
    mtext(side=3,line=0,text=myrangesOI[IX],cex=0.7)
    trans_diff<-  (no.sites>0)[sort(diff.transport.neurotrophins$diff.transport,decreasing=T,index.return=T)$ix]
    myZ      <-  roll_mean(trans_diff,width=myN)
    plot(x=c(1:length(myZ)),y=rev(myZ),type="s",las=1,frame=F,ylab="",xlab="",ylim=c(0,YMAX),col="#C7C5C5")
    mtext(side=2,line=3,text="average fraction with sites",cex=0.7)
    mtext(side=1,line=3,text="transport efficiency",cex=0.7)
    abline(h=0)
  }
  
}


PlotContinuousValuesAlongUTR <- function(motifs="SFPQ_2_HepG2"){
  myrangesL         <- list(c(100,-100),c(3000,2950),c(2750,2700),c(2500,2450),c(2250,2200),c(2000,1950),c(1750,1700),c(1500,1450),c(1400,1350),c(1300,1250),c(1200,1150),c(1100,1050),c(1000,950),c(900,850),c(800,750),c(700,650),c(600,550),c(500,450),c(450,400),c(400,350),c(350,300),c(300,250),c(250,200),c(200,150),c(150,100),c(100,50),c(50,0),c(0,-50),c(-50,-100),c(-100,-150))
  mynames_ragesL   <- unlist(lapply(myrangesL,function(Z)return(paste("[",Z[[1]],":",Z[[2]],"]",sep=""))))
  
  myranges         <- list(c(100,-100),c(3000,2950),c(2750,2700),c(2500,2450),c(1000,950),c(600,550),c(500,450),c(450,400),c(400,350),c(350,300),c(300,250),c(250,200),c(200,150),c(150,100),c(100,50),c(50,0),c(0,-50),c(-50,-100),c(-100,-150))
  mynames_rages    <- unlist(lapply(myranges,function(Z)return(paste("[",as.character(Z[[1]]),":",as.character(Z[[2]]),"]",sep=""))))
  
  myranges <-myrangesL 
  mynames_rages <- mynames_ragesL
  layout(matrix(c(1,1,2,3,4,5),ncol=2,nrow=3,byrow = TRUE))
  oival_welch <- list(ttest_transport_NGF[match(motifs,rownames(ttest_transport_NGF)),],
                      ttest_transport_NT3[match(motifs,rownames(ttest_transport_NT3)),])
  
  oival_ngf_over <- list(
    zscore=Zscore_NGF_overtransported[match(motifs,rownames(Zscore_NGF_overtransported)),],
    fisher=-log10(1E-30+PFisher_NGF_overtransported[match(motifs,rownames(PFisher_NGF_overtransported)),])
  )
  
  oival_nt3_over <- list(
    zscore=Zscore_NT3_overtransported[match(motifs,rownames(Zscore_NT3_overtransported)),],
    fisher=-log10(1E-30+PFisher_NT3_overtransported[match(motifs,rownames(PFisher_NT3_overtransported)),])
    
  )
  
  oival_ngf_under <- list(
    zscore=Zscore_NGF_undertransported[match(motifs,rownames(Zscore_NGF_undertransported)),],
    fisher=-log10(1E-30+PFisher_NGF_undertransported[match(motifs,rownames(PFisher_NGF_undertransported)),])
    
  )
  
  oival_nt3_under <- list(
    zscore=Zscore_NT3_undertransported[match(motifs,rownames(Zscore_NT3_undertransported)),],
    fisher=-log10(1E-30+PFisher_NT3_undertransported[match(motifs,rownames(PFisher_NT3_undertransported)),])
  )
  
  
  yLab=c("Z-score","-log10(Pvalue)-Fisher","Welch test")
  
  
  MIN=min(c(oival_welch[[1]][-1],oival_welch[[2]][-1]),na.rm=TRUE)
  MAX=max(oival_welch[[1]][-1],oival_welch[[2]][-1],na.rm=TRUE)
  
  plot(x=-unlist(lapply(myrangesL,mean))[-1],y=oival_welch[[1]][-1],las=1,ylab="Welch t-test [-log10(P)]",cex.axis=0.8,cex.lab=0.8,main=motifs,cex.main=0.7,col="white",pch=19,cex=0.6,ylim=c(MIN,MAX),frame=FALSE,xlab="distance from 3'end")
  points(x=-unlist(lapply(myrangesL,mean))[-1],y=oival_welch[[1]][-1],col="#81A4D6",pch=19,cex=0.6)
  lines(x=-unlist(lapply(myrangesL,mean))[-1],y=oival_welch[[1]][-1],col="#81A4D6",lwd=1.6)
  points(x=-unlist(lapply(myrangesL,mean))[-1],y=oival_welch[[2]][-1],col="#AE72B0",pch=19,cex=0.6)
  lines(x=-unlist(lapply(myrangesL,mean))[-1],y=oival_welch[[2]][-1],col="#AE72B0",lwd=1.6)
  abline(h=0,lty=1,col="black")
  abline(v=0,lty=2,col="red")
  
  
  
  
  for(i in c(1,2)){
    
    MIN=min(c(oival_ngf_over[[i]][-1],oival_nt3_over[[i]][-1],oival_ngf_under[[i]][-1],oival_nt3_under[[i]][-1]))
    MAX=max(c(oival_ngf_over[[i]][-1],oival_nt3_over[[i]][-1],oival_ngf_under[[i]][-1],oival_nt3_under[[i]][-1]))
    if(i==2&MIN>0){MIN<-0}
    if(i==2&MAX<0){MAX<-0}
    if(i==3){MIN<-0}
    plot(x=-unlist(lapply(myranges,mean))[-1],y=oival_ngf_over[[i]][-1],las=1,ylab=yLab[i],cex.axis=0.8,cex.lab=0.8,main=paste(motifs,"(overtransported)"),cex.main=0.7,col="white",pch=19,cex=0.6,ylim=c(MIN,MAX),frame=FALSE,xlab="distance from 3'end")
    points(x=-unlist(lapply(myranges,mean))[-1],y=oival_ngf_over[[i]][-1],col="#81A4D6",pch=19,cex=0.6)
    lines(x=-unlist(lapply(myranges,mean))[-1],y=oival_ngf_over[[i]][-1],col="#81A4D6",lwd=1.6)
    points(x=-unlist(lapply(myranges,mean))[-1],y=oival_nt3_over[[i]][-1],col="#AE72B0",pch=19,cex=0.6)
    lines(x=-unlist(lapply(myranges,mean))[-1],y=oival_nt3_over[[i]][-1],col="#AE72B0",lwd=1.6)
    if(i==1){legend("topleft",lty=1,pch=19,col=c("#81A4D6","#AE72B0"),leg=c("NGF","NT3"),cex=0.5,bty="n")}
    if(i==1){
      abline(h=0,lty=1,col="black")
      abline(h=c(2.32,-2.32),lty=2,col="grey")
    }
    if(i==2){
      abline(h=2,lty=2,col="grey")
    }
    abline(v=0,col="darkred",lty=3)
    
  }
  
  for(i in c(1,2)){
    #MIN=min(c(oival_ngf_under[[i]][-1],oival_nt3_under[[i]][-1]))
    #MAX=max(c(oival_ngf_under[[i]][-1],oival_nt3_under[[i]][-1]))
    MIN=min(c(oival_ngf_over[[i]][-1],oival_nt3_over[[i]][-1],oival_ngf_under[[i]][-1],oival_nt3_under[[i]][-1]))
    MAX=max(c(oival_ngf_over[[i]][-1],oival_nt3_over[[i]][-1],oival_ngf_under[[i]][-1],oival_nt3_under[[i]][-1]))
    if(i==2&MIN>0){MIN<-0}
    if(i==2&MAX<0){MAX<-0}
    if(i==3){MIN<-0}
    
    
    plot(x=-unlist(lapply(myranges,mean))[-1],y=oival_ngf_under[[i]][-1],las=1,ylab=yLab[i],cex.axis=0.8,cex.lab=0.8,main=paste(motifs,"(undertransported)"),cex.main=0.7,col="white",pch=19,cex=0.6,ylim=c(MIN,MAX),frame=FALSE,xlab="distance from 3'end")
    points(x=-unlist(lapply(myranges,mean))[-1],y=oival_ngf_under[[i]][-1],col="#81A4D6",pch=19,cex=0.6)
    lines(x=-unlist(lapply(myranges,mean))[-1],y=oival_ngf_under[[i]][-1],col="#81A4D6",lwd=1.6)
    points(x=-unlist(lapply(myranges,mean))[-1],y=oival_nt3_under[[i]][-1],col="#AE72B0",pch=19,cex=0.6)
    lines(x=-unlist(lapply(myranges,mean))[-1],y=oival_nt3_under[[i]][-1],col="#AE72B0",lwd=1.6)
    
    if(i==1){
      abline(h=0,lty=1,col="black")
      abline(h=c(2.32,-2.32),lty=2,col="grey")
    }
    if(i==2){
      abline(h=2,lty=2,col="grey")
    }
    abline(v=0,col="darkred",lty=3)
    
  }
  
}



PlotEnrichClips<- function(mymot="BUD13_1_HepG2",myN=250,YMAX=0.6,IX=26){
  require(scales)
  
  if(length(IX)==1){
    clipdat  <- myClips[[IX]]
    clipdat  <- clipdat[match(diff.transport.neurotrophins$uniqueIDs,rownames(clipdat)),]
    col_ix   <- match(mymot,colnames(clipdat))
    no.sites <- clipdat[,col_ix]
  }
  if(length(IX)>1){
    clipdat <- myClips[[IX[1]]]
    for(ix in c(2:length(IX))){
      clipdat <- clipdat+myClips[[IX[ix]]]
    }
    col_ix   <- match(mymot,colnames(clipdat))
    no.sites <- clipdat[,col_ix]
  }
  
  
  par(mfrow=c(2,3),mar=c(3,4,4,3))
  myX <- diff.transport.neurotrophins$transport.ngf[no.sites==0]
  myY <- diff.transport.neurotrophins$transport.ngf[no.sites>0] 
  boxplot(list(unbound=myX,bound=myY),outline=F,frame=F,las=1,col="#81A4D6",ylim=c(-3,3))
  mtext(side=3,line=0,cex=0.6,text=paste("P(Welsh)=",scientific_format(3)(t.test(x=myX,y=myY,var.equal = FALSE)$p.value)))
  mtext(side=3,line=1,cex=0.6,text=paste("P(MW)=",scientific_format(3)(wilcox.test(x=myX,y=myY)$p.value)))
  mtext(side=3,line=2,text=paste("NGF",mymot),cex=0.7)
  mtext(side=2,line=3,text="transport efficiency",cex=0.7)
  abline(h=0)
  
  myX <- diff.transport.neurotrophins$transport.nt3[no.sites==0]
  myY <- diff.transport.neurotrophins$transport.nt3[no.sites>0] 
  boxplot(list(unbound=myX,bound=myY),outline=F,frame=F,las=1,col="#AE72B0",ylim=c(-3,3))
  mtext(side=3,line=0,cex=0.6,text=paste("P(Welsh)=",scientific_format(3)(t.test(x=myX,y=myY,var.equal = FALSE)$p.value)))
  mtext(side=3,line=1,cex=0.6,text=paste("P(MW)=",scientific_format(3)(wilcox.test(x=myX,y=myY)$p.value)))
  mtext(side=3,line=2,text=paste("NT3",mymot),cex=0.7)
  mtext(side=2,line=3,text="transport efficiency",cex=0.7)
  abline(h=0)
  
  myX <- diff.transport.neurotrophins$diff.transport[no.sites==0]
  myY <- diff.transport.neurotrophins$diff.transport[no.sites>0] 
  boxplot(list(unbound=myX,bound=myY),outline=F,frame=F,las=1,col="#C7C5C5",ylim=c(-3,3))
  mtext(side=3,line=0,cex=0.6,text=paste("P(Welsh)=",scientific_format(3)(t.test(x=myX,y=myY,var.equal = FALSE)$p.value)))
  mtext(side=3,line=1,cex=0.6,text=paste("P(MW)=",scientific_format(3)(wilcox.test(x=myX,y=myY)$p.value)))
  mtext(side=3,line=2,text=paste("NGF:NT3",mymot),cex=0.7)
  mtext(side=2,line=3,text="transport efficiency diff",cex=0.7)
  abline(h=0)
  
  require("roll")
  trans_NGF<-  (no.sites>0)[sort(diff.transport.neurotrophins$transport.ngf,decreasing=T,index.return=T)$ix]
  myZ      <-  roll_mean(trans_NGF,width=myN)
  plot(x=c(1:length(myZ)),y=rev(myZ),type="s",xaxt="n",las=1,frame=F,ylab="",ylim=c(0,YMAX),col="#81A4D6")
  mtext(side=2,line=3,text="average fraction with sites",cex=0.7)
  mtext(side=1,line=0,text="transport efficiency",cex=0.7)
  
  trans_NT3<-  (no.sites>0)[sort(diff.transport.neurotrophins$transport.nt3,decreasing=T,index.return=T)$ix]
  myZ      <-  roll_mean(trans_NT3,width=myN)
  plot(x=c(1:length(myZ)),y=rev(myZ),type="s",xaxt="n",las=1,frame=F,ylab="",ylim=c(0,YMAX),col="#AE72B0")
  mtext(side=2,line=3,text="average fraction with sites",cex=0.7)
  mtext(side=1,line=0,text="transport efficiency",cex=0.7)
  
  trans_diff<-  (no.sites>0)[sort(diff.transport.neurotrophins$diff.transport,decreasing=T,index.return=T)$ix]
  myZ      <-  roll_mean(trans_diff,width=myN)
  plot(x=c(1:length(myZ)),y=rev(myZ),type="s",xaxt="n",las=1,frame=F,ylab="",ylim=c(0,YMAX),col="#C7C5C5")
  mtext(side=2,line=3,text="average fraction with sites",cex=0.7)
  mtext(side=1,line=0,text="transport efficiency diff",cex=0.7)
  
}



PlotValuesAlongUTR <- function(motifs="SFPQ_2_HepG2"){
  oival_ngf_over <- list(enr=enr_NGF_overtransported[match(motifs,rownames(enr_NGF_overtransported)),],
                         zscore=Zscore_NGF_overtransported[match(motifs,rownames(Zscore_NGF_overtransported)),],
                         fisher=-log10(1E-30+PFisher_NGF_overtransported[match(motifs,rownames(PFisher_NGF_overtransported)),]),
                         perm=-log10(1E-30+PPerm_NGF_overtransported[match(motifs,rownames(PPerm_NGF_overtransported)),]
                         ))
  
  oival_nt3_over <- list(enr=enr_NT3_overtransported[match(motifs,rownames(enr_NT3_overtransported)),],
                         zscore=Zscore_NT3_overtransported[match(motifs,rownames(Zscore_NT3_overtransported)),],
                         fisher=-log10(1E-30+PFisher_NGF_overtransported[match(motifs,rownames(PFisher_NT3_overtransported)),]),
                         perm=-log10(1E-30+PPerm_NT3_overtransported[match(motifs,rownames(PPerm_NT3_overtransported)),]
                         ))
  
  oival_ngf_under <- list(enr=enr_NGF_undertransported[match(motifs,rownames(enr_NGF_undertransported)),],
                          zscore=Zscore_NGF_undertransported[match(motifs,rownames(Zscore_NGF_undertransported)),],
                          fisher=-log10(1E-30+PFisher_NGF_undertransported[match(motifs,rownames(PFisher_NGF_undertransported)),]),
                          perm=-log10(1E-30+PPerm_NGF_undertransported[match(motifs,rownames(PPerm_NGF_undertransported)),]
                          ))
  
  oival_nt3_under <- list(enr=enr_NT3_undertransported[match(motifs,rownames(enr_NT3_undertransported)),],
                          zscore=Zscore_NT3_undertransported[match(motifs,rownames(Zscore_NT3_undertransported)),],
                          fisher=-log10(1E-30+PFisher_NT3_undertransported[match(motifs,rownames(PFisher_NT3_undertransported)),]),
                          perm=-log10(1E-30+PPerm_NT3_undertransported[match(motifs,rownames(PPerm_NT3_undertransported)),]
                          ))
  
  
  yLab=c("P(expected)/P(observed)","Z-score","-log10(Pvalue)-Fisher","-log10(P-value) --Permutation")
  for(i in c(2,3)){
    barplot(oival_ngf_over[[i]][-1],las=1,ylab=yLab[i],xlab="",cex.axis=0.8,cex.lab=0.8,main=paste(motifs,"(over NGF)"),cex.main=0.7,col="#81A4D6")
  }
  
  for(i in c(2,3)){
    barplot(oival_nt3_over[[i]][-1],las=1,ylab=yLab[i],xlab="",cex.axis=0.8,cex.lab=0.8,main=paste(motifs,"(over NT3)"),cex.main=0.7,col="#AE72B0")
  }
  
  for(i in c(2,3)){
    barplot(oival_ngf_under[[i]][-1],las=1,ylab=yLab[i],xlab="",cex.axis=0.8,cex.lab=0.8,main=paste(motifs,"(under  NGF)"),cex.main=0.7,col="#81A4D6")
  }
  for(i in c(2,3)){
    barplot(oival_nt3_under[[i]][-1],las=1,ylab=yLab[i],xlab="",cex.axis=0.8,cex.lab=0.8,main=paste(motifs,"(under  NT3)"),cex.main=0.7,col="#AE72B0")
  }
}

GetEnrichClip <- function(clip=no_reads[,IX],val_vec=no_reads[,IX],val_mat=myMat[[i]]){
  my.val        <- sum(clip*val_vec,na.rm=TRUE)
  my.val.random <- apply(clip*val_mat,2,function(Z)return(sum(Z,na.rm=TRUE)))
  return((my.val-mean(my.val.random))/sd(my.val.random))
}

standardize <- function(z){
  rowmed <- apply(z, 1, mean)
  rowmad <- apply(z, 1, sd)
  rv <- sweep(z, 1, rowmed)
  rv <- sweep(rv, 1, rowmad, "/")
  return(rv)
}

ExtractTopRegulators <- function(n=10,sub2=t(standardize(t(cbind(as.numeric(ttest_transport_NGF[,colix1[2]]),as.numeric(Zscore_NGF_overtransported[,colix2[2]]))))),XLAB="Welch",YLAB="Zscore over-transported"){
  

  
  km5   <- kmeans(sub2, centers = n, nstart=25)
  Lev3  <- levels(as.factor(km5$cluster))
  Cols  <- c("grey","black","red","blue","green","cyan","magenta","yellow","darkgreen","lightpink")
  names(Cols)<- Lev3
  cols3 <- unlist(lapply(km5$cluster,function(x)return(Cols[match(x,names(Cols))])))
  
  plot(sub2[,1],sub2[,2],col=cols3,pch=19,cex=0.7,xlab="",ylab="",frame=FALSE,las=1)
  mtext(side=1,line=2,text=XLAB)
  mtext(side=2,line=2,text=YLAB)
  grid()
  
  ix_max<- which(apply(sub2,1,sum)==max(apply(sub2,1,sum)))[1]
  ix_min<- which(apply(sub2,1,sum)==min(apply(sub2,1,sum)))[1]
  
  cl_max <- km5$cluster[ix_max]
  cl_min <- km5$cluster[ix_min]
  
  selmost <- which(km5$cluster==cl_max)
  selmin  <- which(km5$cluster==cl_min)
  return(list(selmost,selmin))
}



ExtractEnrichmentMotif <- function(set1=sel.uniqueID.ngf[[4]],set2=NA,mymatch=mycount.all.sfpq){
  
  myIX         <- match(set1,rownames(mymatch))
  myIX         <- myIX[!is.na(myIX)]
  mymotifs     <- colnames(mymatch)
  mycount.seq1 <- apply(mymatch[myIX,]>0,2,sum)
  
  if(!is.na(set2)){
    myIX2        <- match(set2,rownames(mymatch))
    myIX2        <- myIX2[!is.na(myIX2)]
    mycount.seq2 <- apply(mymatch[myIX2,]>0,2,sum)
    
    myp.value    <- unlist(lapply(c(1:length(mycount.seq1)),function(IX)return(fisher.test(rbind(c(mycount.seq1[IX],length(myIX)-mycount.seq1[IX]),
                                                                                                 c(mycount.seq2[IX],length(myIX2)-mycount.seq2[IX])))$p.value)))
    my.enrich    <- do.call(what=rbind,args=lapply(c(1:length(mycount.seq1)),function(IX)return(c(mycount.seq1[IX]/length(myIX)*length(myIX2)/mycount.seq2[IX],mycount.seq2[IX]/length(myIX2)*length(myIX)/mycount.seq1[IX]))))
    colnames(my.enrich)<- c("freq.set1/freq.set2","freq.set2/freq.set1")
    
    out        <- data.frame(count.seq1=mycount.seq1,freq.seq1=mycount.seq1/length(myIX),count.seq2=mycount.seq2,mymotifs,freq.sep2=mycount.seq2/length(myIX2),my.enrich,myp.value)
    return(out)
  }
  
  else{
    mycount.all  <- apply(mymatch>0,2,sum)
    IX           <- lapply(c(1:500),function(Z)return(sample(x=c(1:nrow(mymatch)),size=length(myIX),replace=FALSE)))
    mycount.bg   <- do.call(what=rbind,args=lapply(IX,function(Z)return(apply(mymatch[Z,]>0,2,sum))))
    my.mean      <- apply(mycount.bg,2,mean)
    my.sd        <- apply(mycount.bg,2,sd)
    my.z         <- (mycount.seq1-my.mean)/my.sd
    names(my.z)  <- colnames(mymatch)
    outi         <- data.frame(Z=my.z,count.seq1=mycount.seq1,mean.bg=my.mean,sd.bg=my.sd)
    myp.value    <- unlist(lapply(c(1:length(mycount.seq1)),function(IX)return(fisher.test(rbind(c(mycount.seq1[IX],length(myIX)-mycount.seq1[IX]),
                                                                                                 c(mycount.all[IX],nrow(mymatch)-mycount.all[IX])))$p.value)))
    my.enrich    <- do.call(what=rbind,args=lapply(c(1:length(mycount.seq1)),function(IX)return(c(mycount.seq1[IX]/length(myIX)*length(seq.all.r1)/mycount.all[IX],mycount.all[IX]/nrow(mymatch)*length(myIX)/mycount.seq1[IX]))))
    colnames(my.enrich)<- c("freq.fg/freq.bg","freq.bg/freq.fg")
    IX          <-which(myp.value<=0.05)
    out        <- data.frame(count.seq1=mycount.seq1,freq.seq1=mycount.seq1/length(myIX),count.all=mycount.all,mymotifs,freq.all=mycount.all/nrow(mymatch),my.enrich,myp.value,outi)
    return(out)
  }
}


CreateConsensusMatrix <- function(algn){
  temp           <- msaConvert(algn, "bios2mds::align")
  mymat          <- matrix(0,nrow=length(temp[[1]]),ncol=4)
  colnames(mymat)<-c("A","C","G","U")
  for(i in c(1:length(temp))){
    mymat[,1]<-mymat[,1]+as.integer(temp[[i]]=="A")
    mymat[,2]<-mymat[,2]+as.integer(temp[[i]]=="C")
    mymat[,3]<-mymat[,3]+as.integer(temp[[i]]=="G")
    mymat[,4]<-mymat[,4]+as.integer(temp[[i]]=="U")
  }
  
  mymat <- mymat/length(temp)
}


PlotEffectPerRangeCB<-function(motOI=MOT,BY=0.15){
  myX1      <- diff.transport.neurotrophins$cb.ngf
  myX2      <- diff.transport.neurotrophins$axons.ngf
  sites   <- apply(do.call(what=cbind,args=lapply(c("[250:200]","[200:150]","[150:100]","[100:50]"),function(Z)return(myClips[[match(Z,mynames_rages)]][,match(motOI,colnames(myClips[[1]]))]))),1,sum) 
  sites   <- sites[match(diff.transport.neurotrophins$uniqueIDs,rownames(myClips[[1]]))]>0
  myX1    <- myX1[!is.na(sites)]
  myX2    <- myX2[!is.na(sites)]
  myY     <- sites[!is.na(sites)]
  bins    <- cut(myX1,breaks=quantile(myX1,seq(from=0,by=BY,to=1.0)),include.lowest = FALSE)
  
  mycomps <- list()
  mywelch<- vector()
  q=1
  k=1
  for(LEV in levels(bins)){
    myaxons <- myX2[bins==LEV]
    mycomps[[k]] <- myaxons[myY[bins==LEV]==0]
    k=k+1
    mycomps[[k]] <- myaxons[myY[bins==LEV]>0]
    k=k+1
    mywelch[q]   <- -log10(t.test(mycomps[[k-1]],mycomps[[k-2]])$p.value)*sign(mean(mycomps[[k-1]],na.rm=TRUE)-mean(mycomps[[k-2]],na.rm=TRUE))
    q=q+1
  }
  #par(mfrow=c(2,2))
  boxplot(mycomps,outline=FALSE,las=1,frame=FALSE,col=c("white",colsNGF[1]))
  mtext(side=1,line=2.5,text="expression in CB",cex=0.7)
  mtext(side=2,line=2.5,text="expression in axons [log2]",cex=0.7)
  mtext(side=3,line=0,text=motOI,cex=0.7)
  barplot(mywelch,col=colsNGF[1],las=1,ylab="")
  mtext(side=2,line=2.5,text="P-value [-log10]",cex=0.7)
  grid()
  abline(h=2,col="red",lty=2)
  
  
  
  myX1      <- diff.transport.neurotrophins$cb.nt3
  myX2      <- diff.transport.neurotrophins$axons.nt3
  sites   <- apply(do.call(what=cbind,args=lapply(c("[250:200]","[200:150]","[150:100]","[100:50]"),function(Z)return(myClips[[match(Z,mynames_rages)]][,match(motOI,colnames(myClips[[1]]))]))),1,sum) 
  sites   <- sites[match(diff.transport.neurotrophins$uniqueIDs,rownames(myClips[[1]]))]>0
  myX1    <- myX1[!is.na(sites)]
  myX2    <- myX2[!is.na(sites)]
  myY     <- sites[!is.na(sites)]
  bins    <- cut(myX1,breaks=quantile(myX1,seq(from=0,by=BY,to=1.0)),include.lowest = FALSE)
  
  mycomps <- list()
  mywelch<- vector()
  q=1
  k=1
  for(LEV in levels(bins)){
    myaxons <- myX2[bins==LEV]
    mycomps[[k]] <- myaxons[myY[bins==LEV]==0]
    k=k+1
    mycomps[[k]] <- myaxons[myY[bins==LEV]>0]
    k=k+1
    mywelch[q]   <- -log10(t.test(mycomps[[k-1]],mycomps[[k-2]])$p.value)*sign(mean(mycomps[[k-1]],na.rm=TRUE)-mean(mycomps[[k-2]],na.rm=TRUE))
    q=q+1
  }
  #par(mfrow=c(2,2))
  boxplot(mycomps,outline=FALSE,las=1,frame=FALSE,col=c("white",colsNT3[1]))
  mtext(side=1,line=2.5,text="expression in CB",cex=0.7)
  mtext(side=2,line=2.5,text="expression in axons [log2]",cex=0.7)
  mtext(side=3,line=0,text=motOI,cex=0.7)
  barplot(mywelch,col=colsNT3[1],las=1,ylab="")
  mtext(side=2,line=2.5,text="P-value [-log10]",cex=0.7)
  grid()
  abline(h=2,col="red",lty=2)
}


PlotFractionBoundInFunctionTransport <- function(motOI=MOT,myLIM=1.0,BY=0.1){
  myX     <- diff.transport.neurotrophins$transport.ngf
  sites   <- apply(do.call(what=cbind,args=lapply(c("[250:200]","[200:150]","[150:100]","[100:50]"),function(Z)return(myClips[[match(Z,mynames_rages)]][,match(motOI,colnames(myClips[[1]]))]))),1,sum) 
  sites   <- sites[match(diff.transport.neurotrophins$uniqueIDs,rownames(myClips[[1]]))]>0
  myX     <- myX[!is.na(sites)]
  myY     <- sites[!is.na(sites)]
  bins    <- cut(myX,breaks=quantile(myX,seq(from=0,by=BY,to=1.0)),include.lowest = FALSE)
  palette <- colorRampPalette(colors=c("white", colsNGF[1]))
  barplot(tapply(X=myY,INDEX=bins,FUN=function(W)return(sum(W)/length(W))),col=palette(length(levels(bins))),las=1,frame=FALSE,ylim=c(0,myLIM),ylab="",xaxt="n")
  mtext(side=1,line=0.5,text="transport efficiency",cex=0.7)
  mtext(side=2,line=4,text="fraction of transcript with",cex=0.7)
  mtext(side=2,line=3,text=paste(motOI,"peaks"),cex=0.7)
  
  myX     <- diff.transport.neurotrophins$transport.nt3
  sites   <- apply(do.call(what=cbind,args=lapply(c("[250:200]","[200:150]","[150:100]","[100:50]"),function(Z)return(myClips[[match(Z,mynames_rages)]][,match(motOI,colnames(myClips[[1]]))]))),1,sum) 
  sites   <- sites[match(diff.transport.neurotrophins$uniqueIDs,rownames(myClips[[1]]))]>0
  myX     <- myX[!is.na(sites)]
  myY     <- sites[!is.na(sites)]
  bins    <- cut(myX,breaks=quantile(myX,seq(from=0,by=BY,to=1.0)),include.lowest = FALSE)
  palette <- colorRampPalette(colors=c("white", colsNT3[1]))
  barplot(tapply(X=myY,INDEX=bins,FUN=function(W)return(sum(W)/length(W))),col=palette(length(levels(bins))),las=1,frame=FALSE,ylim=c(0,myLIM),ylab="",xaxt="n")
  mtext(side=1,line=0.5,text="transport efficiency",cex=0.7)
  mtext(side=2,line=4,text="fraction of transcript with",cex=0.7)
  mtext(side=2,line=3,text=paste(motOI,"peaks"),cex=0.7)
  
}

PlotFractionWithSitesPoolsOI <- function(pool=list(
  as.character(diff.transport.neurotrophins$uniqueID)[myVal[[1]]==1],
  as.character(diff.transport.neurotrophins$uniqueID)[myVal[[3]]==1],
  as.character(diff.transport.neurotrophins$uniqueID)[myVal[[2]]==1],
  as.character(diff.transport.neurotrophins$uniqueID)[myVal[[4]]==1]),
  regionOI="[200:150]",motOI="CPSF160_iClip",myMAX=1.0,myLIM=c(-3,3)){
  
  
  clipdat <- myClips[[match(regionOI,mynames_rages)]]
  sites   <- clipdat[match(diff.transport.neurotrophins$uniqueID,rownames(clipdat)),match(motOI,colnames(clipdat))]
  
  frac_bg <- sum(sites>0,na.rm=TRUE)/sum(!is.na(sites))
  frac_NGF_o<- sum(sites[match(pool[[1]],diff.transport.neurotrophins$uniqueID)]>0,na.rm=TRUE)/length(pool[[1]])
  frac_NT3_o<- sum(sites[match(pool[[2]],diff.transport.neurotrophins$uniqueID)]>0,na.rm=TRUE)/length(pool[[2]])
  
  frac_NGF_u<- sum(sites[match(pool[[3]],diff.transport.neurotrophins$uniqueID)]>0,na.rm=TRUE)/length(pool[[3]])
  frac_NT3_u<- sum(sites[match(pool[[4]],diff.transport.neurotrophins$uniqueID)]>0,na.rm=TRUE)/length(pool[[4]])
  
  barplot(c(frac_bg,frac_NGF_o,frac_NT3_o,frac_NGF_u,frac_NT3_u),col=c("grey",colsNGF[1],colNT3[1],colsNGF[1],colNT3[1]),las=1,xlab="",ylab="",ylim=c(0,myMAX))
  mtext(side=2,line=3,text="fraction with sites",cex=0.7)
  mtext(side=3,line=1,text=motOI,cex=0.7)
  mtext(side=3,line=0,text=regionOI,cex=0.7)
  
  
  
  tE <- list(
    diff.transport.neurotrophins$transport.ngf[which(sites==0)],
    diff.transport.neurotrophins$transport.ngf[which(sites>0)],
    diff.transport.neurotrophins$transport.nt3[which(sites==0)],
    diff.transport.neurotrophins$transport.nt3[which(sites>0)]
  )
  boxplot(tE,col=c("grey",colsNGF[1],"grey",colsNT3[1]),ylab="",las=2,frame=FALSE,xlab="",xaxt="n",outline=FALSE,YLIM=myLIM)
  mtext(side=2,line=3,text="transport efficiency",cex=0.7)
  mtext(at=c(1,2,3,4),line=0,side=1,text=c("unbound","bound","unbound","bound"),cex=0.5)
  mtext(side=3,line=1,text=paste("P(NGF=",scientific(t.test(tE[[1]],tE[[2]],var.equal=FALSE)$p.value),sep=""),cex=0.7)
  mtext(side=3,line=0,text=paste("P(NT3=",scientific(t.test(tE[[3]],tE[[4]],var.equal=FALSE)$p.value),sep=""),cex=0.7)
}


PlotFractionWithSitesPoolOI <- function(pool=list(as.character(diff.transport.neurotrophins$uniqueID)[myVal[[1]]==1],as.character(diff.transport.neurotrophins$uniqueID)[myVal[[3]]==1]),
                                        regionOI="[200:150]",motOI="CPSF160_iClip",myMAX=1.0){
  
  
  clipdat <- myClips[[match(regionOI,mynames_rages)]]
  sites   <- clipdat[match(diff.transport.neurotrophins$uniqueID,rownames(clipdat)),match(motOI,colnames(clipdat))]
  
  frac_bg <- sum(sites>0,na.rm=TRUE)/sum(!is.na(sites))
  frac_NGF<- sum(sites[match(pool[[1]],diff.transport.neurotrophins$uniqueID)]>0,na.rm=TRUE)/length(pool[[1]])
  frac_NT3<- sum(sites[match(pool[[2]],diff.transport.neurotrophins$uniqueID)]>0,na.rm=TRUE)/length(pool[[2]])
  barplot(c(frac_bg,frac_NGF,frac_NT3),col=c("grey",colsNGF[1],colNT3[1]),las=1,xlab="",ylab="",ylim=c(0,myMAX))
  mtext(side=2,line=3,text="fraction with sites",cex=0.7)
  mtext(side=3,line=1,text=motOI,cex=0.7)
  mtext(side=3,line=0,text=regionOI,cex=0.7)
  
  
  
  tE <- list(
    diff.transport.neurotrophins$transport.ngf[which(diff.transport.neurotrophins$uniqueID%in%pool[[1]])],
    diff.transport.neurotrophins$transport.ngf[which(sites>0&diff.transport.neurotrophins$uniqueID%in%pool[[1]])],
    diff.transport.neurotrophins$transport.nt3[which(diff.transport.neurotrophins$uniqueID%in%pool[[2]])],
    diff.transport.neurotrophins$transport.nt3[which(sites>0&diff.transport.neurotrophins$uniqueID%in%pool[[2]])]
  )
  boxplot(tE,col=c("grey",colsNGF[1],"grey",colsNT3[1]),ylab="",las=2,frame=FALSE,xlab="",xaxt="n",outline=FALSE)
  mtext(side=2,line=3,text="transport efficiency")
  mtext(at=c(1,2,3,4),line=0,side=1,text=c("unbound","bound","unbound","bound"),cex=0.5)
}


PlotTransportAnalysisGlobal <- function(mot="KHSRP_2_K652",by_te=0.1,by_cb=0.1,limfrac=0.8,myMAX=0.6){
  
  
  PlotEnrichClips(mymot=mot,myN=250,YMAX=myMAX,IX=c(23:26))
  par(mfrow=c(2,2))
  PlotEffectPerRangeCB(motOI=mot,BY=by_cb)
  
  
  par(mar=c(3,4,4,1))
  layout(matrix(c(1,1,4,4,5,5,2,3,6,6,7,7),ncol=6,nrow=2,byrow=TRUE))
  mp=barplot(myFracsROI[match(mot,rownames(myFracsROI)),],col=c("grey",colsNGF[2],colsNGF[1],colsNT3[2],colsNT3[1],colsNGF[2],colsNT3[2]),las=1,xaxt="n")
  mtext(side=2,line=2.5,text="fraction with sites",cex=0.7)
  mtext(side=3,line=3,text=mot,font=2,cex=0.7)
  mtext(side=3,line=0,at=mp[-1],text=scientific(10^(-abs(myPvalsFisherROI[match(mot,rownames(myPvalsFisherROI)),]))[c(1,2,4,5,7,8)],digit=1),cex=0.5)
  mtext(side=3,line=1,at=apply(matrix(mp[-1],ncol=3,nrow=2),2,mean),text=scientific(10^(-abs(myPvalsFisherROI[match(mot,rownames(myPvalsFisherROI)),]))[c(3,6,9)],digit=1),cex=0.5)
  PlotFractionBoundInFunctionTransport(motOI=mot,myLIM=limfrac,BY=by_te)
  plotEnrichAlongUTRSimple(mots=mot,sel=c(2:30),scaling=FALSE,mytests=list(ttest_transport_NGF,ttest_transport_NT3),mytestname="Welsch test",YLAB="-log10(P-value)",LWD=1.5)
  plotEnrichAlongUTRSimple(mots=mot,sel=c(15:30),scaling=FALSE,mytests=list(ttest_transport_NGF,ttest_transport_NT3),mytestname="Welsch test",YLAB="-log10(P-value)",LWD=1.5)
}

ExtractRegulators <- function(thetest=fisher_proxi_NT3_reg,limdist=(-300),limsig=2.0){
  remodel_positive_NT3             <- lapply(c(1:ncol(thetest)),function(W){
    LIMS <- boxplot(thetest[,W],plot=FALSE)$stats[5,1]
    return(rownames(thetest)[(thetest[,W])>LIMS])
  })
  remodel_positive_NT3_names        <- lapply(remodel_positive_NT3,function(W)return(unique(unlist(lapply(W,function(Z)return((unlist(strsplit(Z,split="_"))[1])))))))
  names(remodel_positive_NT3_names) <- myranges
  
  myPositiveNT3                <- unique(unlist(remodel_positive_NT3_names))
  MOTS_POS                     <-unique(unlist(remodel_positive_NT3))
  START                        <- unlist(lapply(myranges,function(Z)return(-Z[1])))
  END                          <- unlist(lapply(myranges,function(Z)return(-Z[2])))
  
  ROIp                       <- mynames_rages[unlist(lapply(MOTS_POS,function(mymot)sort(thetest[match(mymot,rownames(thetest)),-1],decreasing=TRUE,index.return=TRUE)$ix[1]+1))]
  ValuesOI                   <- unlist(mapply(A=MOTS_POS,B=ROIp,function(A,B){
    thetest[match(A,rownames(thetest)),match(B,colnames(thetest))]
  }))
  
  
  
  Id_positive_NT3_regulators <- data.frame(start.Ip=START[match(ROIp,mynames_rages)],
                                           end.Ip=END[match(ROIp,mynames_rages)],
                                           Clip=MOTS_POS,
                                           vals_Ip=ValuesOI,
                                           GS=unlist(lapply(MOTS_POS,function(Z)return((unlist(strsplit(Z,split="_"))[1])))))
  
  Id_positive_NT3_regulators <- Id_positive_NT3_regulators[Id_positive_NT3_regulators$start.Ip>limdist,]
  Id_positive_NT3_regulators <- Id_positive_NT3_regulators[Id_positive_NT3_regulators$vals_Ip>limsig,]
  
  
  
  
  remodel_positive_NT3             <- lapply(c(1:ncol(thetest)),function(W){
    LIMS <- boxplot(-thetest[,W],plot=FALSE)$stats[5,1]
    return(rownames(-thetest)[(fisher_proxi_NT3_reg[,W])>LIMS])
  })
  remodel_positive_NT3_names        <- lapply(remodel_positive_NT3,function(W)return(unique(unlist(lapply(W,function(Z)return((unlist(strsplit(Z,split="_"))[1])))))))
  names(remodel_positive_NT3_names) <- myranges
  
  myPositiveNT3                <- unique(unlist(remodel_positive_NT3_names))
  MOTS_POS                     <-unique(unlist(remodel_positive_NT3))
  START                        <- unlist(lapply(myranges,function(Z)return(-Z[1])))
  END                          <- unlist(lapply(myranges,function(Z)return(-Z[2])))
  
  ROIp                       <- mynames_rages[unlist(lapply(MOTS_POS,function(mymot)sort(-thetest[match(mymot,rownames(thetest)),-1],decreasing=TRUE,index.return=TRUE)$ix[1]+1))]
  ValuesOI                   <- unlist(mapply(A=MOTS_POS,B=ROIp,function(A,B){
    thetest[match(A,rownames(thetest)),match(B,colnames(thetest))]
  }))
  
  
  
  Id_negative_NT3_regulators <- data.frame(start.Ip=START[match(ROIp,mynames_rages)],
                                           end.Ip=END[match(ROIp,mynames_rages)],
                                           Clip=MOTS_POS,
                                           vals_Ip=ValuesOI,
                                           GS=unlist(lapply(MOTS_POS,function(Z)return((unlist(strsplit(Z,split="_"))[1])))))
  
  Id_negative_NT3_regulators <- Id_negative_NT3_regulators[Id_negative_NT3_regulators$vals_Ip<(-limsig),]
  
  
  Id_negative_NT3_regulators <- Id_negative_NT3_regulators[Id_negative_NT3_regulators$start.Ip>limdist,]
  return(list(Id_positive_NT3_regulators,Id_negative_NT3_regulators))
  
}