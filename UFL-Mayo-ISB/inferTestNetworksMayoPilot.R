require(synapseClient)
synapseLogin()

res <- synQuery('select name, id from file where projectId=="syn2580853" and center=="UFL-Mayo-ISB" and disease=="Alzheimers Disease" and platform=="IlluminaHiSeq2000" and other=="geneCounts" and other=="Normalized"')
synObj <- synGet(res$file.id)

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
x <- scale(tcxAlzheimer1000)
y <- tcxAlzheimer1000[,3]
x <- x[,-3]
eigen <- svd(x)$d^2
system.time(res <- metaReg(y = y,x = x,eigen = eigen))



