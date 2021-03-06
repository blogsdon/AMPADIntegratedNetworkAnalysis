require(synapseClient)
synapseLogin()

res <- synQuery('select name, id from file where projectId=="syn2580853" and center=="UFL-Mayo-ISB" and disease=="Alzheimers Disease" and platform=="IlluminaHiSeq2000" and other=="geneCounts" and other=="Normalized"')
synObj <- synGet(res$file.id)

source('ensembl2symbol.R')
#Read in data
tcxAlzheimer <- read.table(synObj@filePath)

#log transformation + boundary condition to prevent log transformation from failing
tcxAlzheimer <- log(tcxAlzheimer+0.1)

#extract top 1000 varying genes
genesd <- apply(tcxAlzheimer,1,sd)
plot(sort(genesd))
top1000 <- order(genesd,decreasing=T)[1:1000]
tcxAlzheimer1000 <- data.matrix(t(as.matrix(tcxAlzheimer)[top1000,]))

library(metaNet)
sparrowNet <- matrix(0,1e3,1e3)
lassoNet <- matrix(0,1e3,1e3)
ssNet <- matrix(0,1e3,1e3)
ridgeNet <- matrix(0,1e3,1e3)
rfNet <- matrix(0,1e3,1e3)
xnorm <- scale(tcxAlzheimer1000)
for (i in 1:1000){
  y <- xnorm[,i]
  x <- xnorm[,-i]
  eigen <- svd(x)$d^2
  system.time(res <- metaReg(y = y,x = x,eigen = eigen))
  sparrowNet[i,-i] <- res[,1];
  lassoNet[i,-i] <- res[,2];
  ssNet[i,-i] <- res[,3];
  ridgeNet[i,-i] <- res[,4];
  rfNet[i,-i] <- res[,5];
  cat('finished ',i,'of 1000\n')
}

for (i in 1:1000){
  y <- xnorm[,i]
  x <- xnorm[,-i]
  res <- vbsrWrapperZ(y=y,x=x)
  sparrowNet[i,-i] <- res
  cat('finished ',i,'of 1000\n')
}

sparrowNetAvg <- sparrowNet/2+t(sparrowNet)/2
sparrowNetAvg <- sparrowNetAvg^2>qchisq(0.05/choose(1000,2),1,lower.tail=F)
lassoNetAvg <- lassoNet/2+t(lassoNet)/2
ssNetAvg <- ssNet/2+t(ssNet)/2
ridgeNetAvg <- ridgeNet/2+t(ridgeNet)/2
rfNetAvg <- rfNet/2+t(rfNet)/2

sparrowHub <- signif(scale(rowSums(sparrowNetAvg)),2)
lassoHub <- signif(scale(rowSums(lassoNetAvg)),2)
ssHub <- signif(scale(rowSums(ssNetAvg)),2)
ridgeHub <- signif(scale(rowSums(ridgeNetAvg)),2)
rfHub <- signif(scale(rowSums(rfNetAvg)),2)
aggHub <- signif(sparrowHub+lassoHub+ssHub+ridgeHub+rfHub,2)

geneId <- ensembl2symbol(colnames(tcxAlzheimer1000))
names(geneId) <- colnames(tcxAlzheimer1000)
dfHub <- data.frame(geneId,sparrowHub,lassoHub,ssHub,ridgeHub,rfHub,aggHub,row.names=colnames(tcxAlzheimer1000))



alzGenes <- scan('geneLists/geneCardsAlzheimer.txt',what='character')
alzGenetics <- unique(scan('geneLists//alzheimerGeneticsGenes.txt',what='character'))
alzheimerGene<-(dfHub$geneId%in%alzGenes)
alzheimerGWAS <- (dfHub$geneId%in%alzGenetics)
dfHub <- data.frame(dfHub,alzheimerGene)

png(file='UFL-Mayo-ISB/Figures/hubPairs.png',height=3.42,width=3.42,units='in',res=300,pointsize=10)
pairs(dfHub[,2:6],pch=16,cex=.8)
dev.off()

#enrichment(alzGenes,dfHub$geneId[order(dfHub$sparrowHub,decreasing=T)[1:500]],dfHub$geneId)
sparrowEnrich <- data.matrix(as.matrix(enrichmentPath(targetList = as.character(alzGenes), rankedList = as.character(dfHub$geneId[order(dfHub$sparrowHub,decreasing=T)]))))
lassoEnrich <- data.matrix(as.matrix(enrichmentPath(targetList = as.character(alzGenes), rankedList = as.character(dfHub$geneId[order(dfHub$lassoHub,decreasing=T)]))))
ssEnrich <- data.matrix(as.matrix(enrichmentPath(targetList = as.character(alzGenes), rankedList = as.character(dfHub$geneId[order(dfHub$ssHub,decreasing=T)]))))
ridgeEnrich <- data.matrix(as.matrix(enrichmentPath(targetList = as.character(alzGenes), rankedList = as.character(dfHub$geneId[order(dfHub$ridgeHub,decreasing=T)]))))
rfEnrich <- data.matrix(as.matrix(enrichmentPath(targetList = as.character(alzGenes), rankedList = as.character(dfHub$geneId[order(dfHub$rfHub,decreasing=T)]))))
aggEnrich <- data.matrix(as.matrix(enrichmentPath(targetList = as.character(alzGenes), rankedList = as.character(dfHub$geneId[order(dfHub$aggHub,decreasing=T)]))))

pvals <- apply(cbind(sparrowEnrich[,3],lassoEnrich[,3],ssEnrich[,3],ridgeEnrich[,3],rfEnrich[,3],aggEnrich[,3]),2,as.numeric)
enr<- apply(cbind(sparrowEnrich[,2],lassoEnrich[,2],ssEnrich[,2],ridgeEnrich[,2],rfEnrich[,2],aggEnrich[,2]),2,as.numeric)


png(file='UFL-Mayo-ISB//Figures/enrichmentPlot.png',height=3.42,width=3.42,units='in',res=300,pointsize=10)
matplot(matrix(rep(1:1000,6),1000,6),-log10(pvals),'l',lty=1:6,lwd=2,ylim=c(0,8),xlab='Ranked Expression Hubs',main='Enrichment for known Alz Genes\nin TCX RNAseq networks')
legend('topright',c('Sparrow','Lasso','Stability Selection','Ridge','Random Forest','Aggregate'),col=1:6,lty=1:6,lwd=3,cex=.6)
dev.off()
matplot(matrix(rep(1:1000,6),1000,6),enr,'l',lty=1:6,lwd=3)

tcresult<-as.tableColumns(dfHub)
cols<-tcresult$tableColumns
fileHandleId<-tcresult$fileHandleId
projectId<-"syn3243392"
schema<-TableSchema(name="TCX Alzheimer Preliminary Integrated Analysis", parent=projectId, columns=cols)
table<-Table(schema, fileHandleId)
table<-synStore(table, retrieveData=TRUE)
#onWeb(table)
