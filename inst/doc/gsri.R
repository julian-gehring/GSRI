###################################################
### chunk number 1: settings
###################################################
#line 53 "gsri.Rnw"
set.seed(1)
options(width=65, SweaveHooks=list(fig=function() par(mar=c(5.1, 5.1, 4.1, 1.1))))


###################################################
### chunk number 2: preparePdfCdfFig
###################################################
#line 81 "gsri.Rnw"
x <- seq(0, 1, len=50)
r <- 0.2
rate <- 10
d0 <- dunif(x)
d1 <- dexp(x, rate)
d <- (1-r)*d0 + r*d1
c0 <- punif(x)
c1 <- pexp(x, rate)
c <- (1-r)*c0 + r*c1
x <- seq(0, 1, len=50)
cex <- 1.5


###################################################
### chunk number 3: pdfFig
###################################################
#line 95 "gsri.Rnw"
plot(x, d, type="n", xaxt="n", yaxt="n", ylim=c(0, max(d)), xlab=expression(paste("p-value ", italic(p))), ylab=expression(paste("probability density ", rho(p))), cex.lab=cex)
lines(c(0, 1), rep(1-r, 2), lwd=2, col="darkgray")
lines(x, d, lwd=2)
axis(1, at=seq(0, 1, len=5), labels=c(0, NA, 0.5, NA, 1), cex.axis=cex)
axis(2, at=seq(0, max(d), by=0.5), labels=c(0, NA, 1, NA, 2, NA), cex.axis=cex)
axis(2, at=1-r, labels=expression(paste(1-italic(r))), cex.axis=cex)
text(0.5, (1-r)/2, labels=expression(1-italic(r)), cex=cex, offset=NULL, adj=c(0.5, 0.5))
text(0.09, 1.05, labels=expression(italic(r)), cex=cex, offset=NULL, adj=c(0, 0))


###################################################
### chunk number 4: cdfFig
###################################################
#line 106 "gsri.Rnw"
plot(x, c, type="n", xaxt="n", yaxt="n", xlim=c(0, 1), ylim=c(0, 1), xlab=expression(paste("p-value ", italic(p))), ylab=expression(paste("cumulative density ", F(p))), cex.lab=cex)
lines(c(0, 1), c(r, 1), lwd=2, col="darkgray")
lines(x, c, lwd=2)
axis(1, at=seq(0, 1, len=5), labels=c(0, NA, 0.5, NA, 1), cex.axis=cex)
axis(2, at=seq(0, 1, len=5), labels=c(0, NA, 0.5, NA, 1), cex.axis=cex)
axis(2, at=r, labels=expression(paste(italic(r))), cex.axis=cex)
text(0.75, 0.75, labels=expression(1-italic(r)), cex=cex, offset=NULL, adj=c(0, 0))


###################################################
### chunk number 5: extractData
###################################################
#line 154 "gsri.Rnw"
library(Biobase)
data(sample.ExpressionSet)
eset <- sample.ExpressionSet
eset
exprs <- exprs(eset)
phenotypes <- pData(phenoData(eset))
summary(phenotypes)


###################################################
### chunk number 6: gsriAllProbes
###################################################
#line 178 "gsri.Rnw"
library(GSRI)
gAllProbes <- gsri(eset, phenotypes$type)
gAllProbes


###################################################
### chunk number 7: gsriAllGenesSet
###################################################
#line 200 "gsri.Rnw"
library(GSEABase)
gs <- GeneSet(eset, setName="allGenes")
ind <- grep("^AFFX", geneIds(gs), invert=TRUE)
gsAllGenes <- gs[ind]
gsAllGenes


###################################################
### chunk number 8: gsriAllGenes
###################################################
#line 208 "gsri.Rnw"
gAllGenesType <- gsri(eset, phenotypes$type, gsAllGenes, name="allGenesType")
gAllGenesSex <- gsri(exprs, phenotypes$sex, gsAllGenes, name="allGenesSex")
gAllGenesType
gAllGenesSex


###################################################
### chunk number 9: gsriChr17
###################################################
#line 224 "gsri.Rnw"
gsChr17 <- getBroadSets(system.file("extdata", "c1chr17.xml", package="GSRI"))
gsChr17
gChr17 <- gsri(eset, phenotypes$type, gsChr17)
gChr17


###################################################
### chunk number 10: gsriCol5
###################################################
#line 246 "gsri.Rnw"
gmt <- getGmt(system.file("extdata", "c1c10.gmt", package="GSRI"))
gCol5 <- gsri(eset, phenotypes$type, gmt)
gCol5
gCol5Sort <- sortGsri(gCol5, c("nRegGenes", "pRegGenes"))
summary(gCol5Sort)
exportFile <- tempfile()
export(gCol5Sort, exportFile)


###################################################
### chunk number 11: gsriScoreFtest
###################################################
#line 269 "gsri.Rnw"
phenotypes$class <- cut(phenotypes$score, seq(0, 1, length.out=4), label=c("low", "medium", "high"))
summary(phenotypes$class)
g3 <- gsri(eset, phenotypes$class, gsChr17, test=rowF)
g3


###################################################
### chunk number 12: gsriAllGenesArgs
###################################################
#line 282 "gsri.Rnw"
g3arg2 <- gsri(eset, phenotypes$class, gsChr17, test=rowF, name="chr17_2", nBoot=50, alpha=0.1, grenander=FALSE)
g3arg2


###################################################
### chunk number 13: limmaTest
###################################################
#line 298 "gsri.Rnw"
library(limma)
limmaTest <- function(exprs, groups, id, index, testArgs) {
    design <- cbind(offset=1, diff=groups)
    fit <- lmFit(exprs[ ,index], design)
    fit <- eBayes(fit)
    pval <- fit$p.value[id,"diff"]
    return(pval)
}
g3Limma <- gsri(eset, phenotypes$type, gsChr17, test=limmaTest)
g3Limma


###################################################
### chunk number 14: plot1
###################################################
#line 319 "gsri.Rnw"
plot(g3)


###################################################
### chunk number 15: plot2
###################################################
#line 334 "gsri.Rnw"
plot(gCol5, 5, ecdf=list(type="o"), plot=list(xlab="p", main="GSRI results: chr9"), reg=list(col="lightblue", lty=1, lwd=1.5), gsri=list(col="darkblue"))


###################################################
### chunk number 16: weights
###################################################
#line 360 "gsri.Rnw"
library(hgu95av2.db)
gNames <- rownames(exprs(eset))
ind <- Lkeys(hgu95av2GO) %in% gNames
evidence <- factor(toTable(hgu95av2GO)[ind,"Evidence"])
summary(evidence)
l <- lapply(gNames, function(name, names, evidence) evidence[names %in% name], gNames, evidence)
expInd <- sapply(l, function(l) any(l %in% "EXP"))
goWeight <- rep(0.5, length.out=length(expInd))
goWeight[expInd] <- 1
gCol5go <- gsri(eset, phenotypes$type, weight=goWeight)
gCol5go
gCol5go2 <- gsri(eset, phenotypes$type, gmt, weight=goWeight)
gCol5go2


###################################################
### chunk number 17: sessionInfo
###################################################
#line 387 "gsri.Rnw"
toLatex(sessionInfo(), locale=FALSE)


