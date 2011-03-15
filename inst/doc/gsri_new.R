###################################################
### chunk number 1: settings
###################################################
#line 47 "gsri_new.Rnw"
set.seed(1)
options(width=65)


###################################################
### chunk number 2: attachPackages
###################################################
#line 63 "gsri_new.Rnw"
library(GSRI)
library(GSEABase)
library(Biobase)
data(sample.ExpressionSet)


###################################################
### chunk number 3: extractData
###################################################
#line 70 "gsri_new.Rnw"
eset <- sample.ExpressionSet
eset
exprs <- exprs(eset)
phenotypes <- pData(phenoData(eset))
summary(phenotypes)


###################################################
### chunk number 4: n1
###################################################
#line 84 "gsri_new.Rnw"
gAllProbes <- gsri(eset, phenotypes$type)
gAllProbes


###################################################
### chunk number 5: n2
###################################################
#line 89 "gsri_new.Rnw"
gs <- GeneSet(eset, setName="allGenes")
ind <- grep("^AFFX", geneIds(gs), invert=TRUE)
gsAllGenes <- gs[ind]
gsAllGenes
gAllGenesType <- gsri(eset, phenotypes$type, gsAllGenes, name="allGenesType")
gAllGenesSex <- gsri(exprs, phenotypes$sex, gsAllGenes, name="allGenesSex")
gAllGenesType
gAllGenesSex


###################################################
### chunk number 6: plot0
###################################################
#line 100 "gsri_new.Rnw"
plot(gAllGenesType)


###################################################
### chunk number 7: n1
###################################################
#line 113 "gsri_new.Rnw"
gmt <- getGmt(system.file("extdata", "c1c10.gmt", package="GSRI"))
gCol5 <- gsri(eset, phenotypes$type, gmt)
gCol5
gCol5Sort <- sortGsri(gCol5, c("nRegGenes", "pRegGenes"))
exportFile <- tempfile()
export(gCol5Sort, exportFile)


###################################################
### chunk number 8: n1
###################################################
#line 122 "gsri_new.Rnw"
gsChr17 <- getBroadSets(system.file("extdata", "c1chr17.xml", package="GSRI"))
gsChr17
gChr17 <- gsri(eset, phenotypes$type, gsChr17)


###################################################
### chunk number 9: n1
###################################################
#line 128 "gsri_new.Rnw"
phenotypes$class <- cut(phenotypes$score, seq(0, 1, length.out=4), label=c("low", "medium", "high"))
summary(phenotypes$class)
g3 <- gsri(eset, phenotypes$class, gsChr17, test=rowF)
g3


###################################################
### chunk number 10: plot1
###################################################
#line 135 "gsri_new.Rnw"
plot(g3)


###################################################
### chunk number 11: plot2
###################################################
#line 145 "gsri_new.Rnw"
plot(gCol5, 5, ecdf=list(type="s"))


###################################################
### chunk number 12: weights
###################################################
#line 155 "gsri_new.Rnw"
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
### chunk number 13: sessionInfo
###################################################
#line 190 "gsri_new.Rnw"
toLatex(sessionInfo(), locale=FALSE)


