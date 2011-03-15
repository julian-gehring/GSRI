###################################################
### chunk number 1: settings
###################################################
#line 103 "gsri_new.Rnw"
set.seed(1)
options(width=65)


###################################################
### chunk number 2: attachPackages
###################################################
#line 119 "gsri_new.Rnw"
library(GSRI)
library(GSEABase)
library(Biobase)
data(sample.ExpressionSet)


###################################################
### chunk number 3: extractData
###################################################
#line 126 "gsri_new.Rnw"
eset <- sample.ExpressionSet
eset
exprs <- exprs(eset)
phenotypes <- pData(phenoData(eset))
summary(phenotypes)


###################################################
### chunk number 4: gsriAllProbes
###################################################
#line 142 "gsri_new.Rnw"
gAllProbes <- gsri(eset, phenotypes$type)
gAllProbes


###################################################
### chunk number 5: gsriAllGenesSet
###################################################
#line 161 "gsri_new.Rnw"
gs <- GeneSet(eset, setName="allGenes")
ind <- grep("^AFFX", geneIds(gs), invert=TRUE)
gsAllGenes <- gs[ind]
gsAllGenes


###################################################
### chunk number 6: gsriAllGenes
###################################################
#line 168 "gsri_new.Rnw"
gAllGenesType <- gsri(eset, phenotypes$type, gsAllGenes, name="allGenesType")
gAllGenesSex <- gsri(exprs, phenotypes$sex, gsAllGenes, name="allGenesSex")
gAllGenesType
gAllGenesSex


###################################################
### chunk number 7: gsriChr17
###################################################
#line 184 "gsri_new.Rnw"
gsChr17 <- getBroadSets(system.file("extdata", "c1chr17.xml", package="GSRI"))
gsChr17
gChr17 <- gsri(eset, phenotypes$type, gsChr17)
gChr17


###################################################
### chunk number 8: gsriCol5
###################################################
#line 206 "gsri_new.Rnw"
gmt <- getGmt(system.file("extdata", "c1c10.gmt", package="GSRI"))
gCol5 <- gsri(eset, phenotypes$type, gmt)
gCol5
gCol5Sort <- sortGsri(gCol5, c("nRegGenes", "pRegGenes"))
summary(gCol5Sort)
exportFile <- tempfile()
export(gCol5Sort, exportFile)


###################################################
### chunk number 9: gsriScoreFtest
###################################################
#line 229 "gsri_new.Rnw"
phenotypes$class <- cut(phenotypes$score, seq(0, 1, length.out=4), label=c("low", "medium", "high"))
summary(phenotypes$class)
g3 <- gsri(eset, phenotypes$class, gsChr17, test=rowF)
g3


###################################################
### chunk number 10: gsriAllGenesArgs
###################################################
#line 242 "gsri_new.Rnw"
g3arg2 <- gsri(eset, phenotypes$class, gsChr17, test=rowF, name="chr17_2", nBoot=50, alpha=0.1, grenander=FALSE)
g3arg2


###################################################
### chunk number 11: limmaTest
###################################################
#line 258 "gsri_new.Rnw"
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
### chunk number 12: plot1
###################################################
#line 279 "gsri_new.Rnw"
plot(g3)


###################################################
### chunk number 13: plot2
###################################################
#line 294 "gsri_new.Rnw"
plot(gCol5, 5, ecdf=list(type="o"), plot=list(xlab="p", main="GSRI results: chr9"), reg=list(col="lightblue", lty=1, lwd=1.5), gsri=list(col="darkblue"))


###################################################
### chunk number 14: weights
###################################################
#line 320 "gsri_new.Rnw"
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
### chunk number 15: sessionInfo
###################################################
#line 347 "gsri_new.Rnw"
toLatex(sessionInfo(), locale=FALSE)


