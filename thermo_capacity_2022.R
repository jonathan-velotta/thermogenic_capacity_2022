# load packages
library(edgeR)
library(WGCNA)
library(reshape)
library(Hmisc)
library(lme4)
library(lmerTest)
library(multcomp)
library(gprofiler2)
library(igraph)
library(plyr)
library(vegan)
thermo0 <- read.csv("thermo_capacity_counts.csv") #raw counts table
colnames(thermo0)
thermo <- thermo0[, -c(2:6)] #omit chromosome, start/end, and length data
rownames(thermo) <- thermo$Geneid #replace rownames with GeneIDs
thermo$Geneid <- NULL #rownames MUST be GeneIDs. Columns must ONLY be individuals
colnames(thermo)
id <- read.csv("thermo_capacity_rnaseq_samples.csv") #ID names for individuals
id <- id[order(id$mouse_id),] #order sample ID's alphabetically
colnames(thermo)
dim(thermo) #should be 89 columns with individuals in alphabetical order
dim(id) #make sure 89 rows in alphabetical order
colnames(thermo) <- id$ID #WARNING. This replaces the individuals IDs, so make sure they are in order
head(thermo)
dim(thermo)
#######################################################################################################################################################
#gastroc analysis
gastroc.id <- subset(id, tissue=="gas") #subset ids by tissue
gastroc.id <- gastroc.id[order(gastroc.id$mouse_id),]
gas.id <- gastroc.id$ID #make a vector containing gastroc mouse ids
gastroc <- thermo[,colnames(thermo) %in% gas.id] #subset ontogeny data by gastroc mouse ids
head(gastroc)#only gastroc data
dim(gastroc)
match(gastroc.id$ID, colnames(gastroc)) #verify that colnames and sample id names are in the same order: 1-47
colnames(gastroc) <- gastroc.id$mouse_id ###WARNING: this step removes tissue identifier in ID names. Proceed with caution.
match(gastroc.id$mouse_id, colnames(gastroc)) #verify that colnames and sample id names are in the same order: 1-47
# write.csv(gastroc, "gastroc_raw_counts_thermocapacity.csv")
########
population <- gastroc.id$population #specify population
acclimation <- gastroc.id$acclimation #specify acclimation treatment: N=normoxia, H=hypoxia, C=cold, CH=cold, hypoxia
table <- data.frame(population, acclimation)
group <- factor(paste(table$population, table$acclimation, sep="_"))
cbind(table, group=group)
table$population = as.factor(table$population)
# table$acclimation <- relevel(table$acclimation, ref="N")
table$population <- relevel(table$population, ref="LN")
design <- model.matrix(~population*acclimation, data=table)
#filter data
gastroc$mean = rowMeans(gastroc) #rowMeans takes means of each row
keep_gastroc = subset(gastroc, mean >= 10) #filter by means
dim(keep_gastroc)
keep_gastroc$mean = NULL #clean up dataset
dim(keep_gastroc)
y0 <- DGEList(counts=keep_gastroc, group=group) # make a DGE list
y <- calcNormFactors(y0) # normalize
plotMDS(y) #exploration plots
y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
gastroc.norm <- cpm(y0, log=TRUE, prior.count=1, normalized.lib.sizes=TRUE) #cpm normalized and log transformed expression dataset
#write.csv(gastroc.norm, file="thermo_gastroc_norm_counts.csv")
plotMDS(gastroc.norm)

########## PCA separated by acclimation treatment
pca <- prcomp(t(gastroc.norm), scale=FALSE)
summary(pca)
pc = as.data.frame(pca$x)
pc$population <- population
pc$acclimation <- acclimation
pc$family <- gastroc.id$family
ln <- subset(pc, population == "LN")
lnN <- subset(ln, acclimation=="N")
lnH <- subset(ln, acclimation=="H")
lnC <- subset(ln, acclimation=="C")
lnCH <- subset(ln, acclimation=="CH")
me <- subset(pc, population == "ME")
meN <- subset(me, acclimation=="N")
meH <- subset(me, acclimation=="H")
meC <- subset(me, acclimation=="C")
meCH <- subset(me, acclimation=="CH")
# pdf("Figures/pca_gastroc.pdf", h=8, w=8)
par(mfrow=c(1,1),oma=c(4,4,1,1),mar=c(4,4,1,1))
plot(lnN$PC1,lnN$PC1, type='n', pch=19, xlim=c(-50,30), ylim = c(-50,60), cex=2, yaxt='n', xaxt='n',ylab="", xlab="",col="black", lwd=2, bg=as.numeric(pc$population))
axis(side=1, lwd=2, cex.axis=2, las=1)
mtext("PC1 (12.5%)", side=1, line=4, cex=2)
axis(side=2, lwd=2, cex.axis=2, las=1)
mtext("PC2 (7.6%)", side=2, cex=2, line=4)
ordiellipse(pca,pc$population,conf=0.95, draw="polygon", col="gray92", border = "white", lwd=3)
points(lnN$PC1, lnN$PC2, pch=21, col="darkred", bg="darkred",cex=2.5, lwd=2)
points(meN$PC1, meN$PC2, pch=21, col="darkred", bg="darkred", cex=2.5, lwd=2)
points(lnH$PC1, lnH$PC2, pch=21, col="red", bg="red", cex=2.5, lwd=2)
points(meH$PC1, meH$PC2, pch=21, col="red", bg="red", cex=2.5, lwd=2)
points(lnC$PC1, lnC$PC2, pch=21, col="darkblue", bg="darkblue", cex=2.5, lwd=2)
points(meC$PC1, meC$PC2, pch=21, col="darkblue", bg="darkblue", cex=2.5, lwd=2)
points(lnCH$PC1, lnCH$PC2, pch=21, col="blue", bg="blue", cex=2.5, lwd=2)
points(meCH$PC1, meCH$PC2, pch=21, col="blue", bg="blue", cex=2.5, lwd=3)
# text(pc$PC1, pc$PC2, pc$family, cex=1, pos=2)
box(which="plot", lty="solid", lwd=2)
# dev.off()

# population PCA
# pdf("Figures/pop_pca_gastroc.pdf", h=8, w=8)
par(mfrow=c(1,1),oma=c(3,4,1,1),mar=c(3,4,1,1), bg="white")
plot(ln$PC1,ln$PC2, type='n', pch=19, xlim=c(-50,40), ylim = c(-40,30), cex=2, yaxt='n', xaxt='n',ylab="", xlab="",col="black", lwd=2, bg=as.numeric(pc$population))
axis(side=1, lwd=2, cex.axis=2, las=1)
mtext("PC1 (12.8%)", side=1, line=4, cex=2)
axis(side=2, lwd=2, cex.axis=2, las=1)
mtext("PC2 (6.8%)", side=2, cex=2, line=4)
ordiellipse(pca,pc$population,conf=0.95, draw="polygon", col="gray92", border = "white", lwd=5)
points(ln$PC1, ln$PC2, pch=21, col="darkorange", bg="darkorange",cex=4, lwd=2)
points(me$PC1, me$PC2, pch=21, col="darkblue", bg="darkblue", cex=4, lwd=2)
box(which="plot", lty="solid", lwd=2)
# dev.off()

####################################################################################################################################################### GASTROC WGCNA #####################################################
#WGCNA for gastroc
head(gastroc.norm)#normalized read counts for gastroc only
dim(gastroc.norm)
Expr0 = as.data.frame(t(gastroc.norm)) #transpose expression data for further analysis
rownames(Expr0) = colnames(gastroc.norm)
gsg = goodSamplesGenes(Expr0, verbose = 3) #check for genes and samples with too many missing values
gsg$allOK
Expr = Expr0
#cluster the samples to check for outliers
sampleTree = hclust(dist(Expr), method = "average");
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

###Import Trait Data
traitData0 <- read.csv("thermo_capacity_traits_all.csv")#all trait data
dev.off()
boxplot(vo2_mass~population*acclimation, data=traitData0)
Samples <- rownames(Expr)
traitData <- traitData0[traitData0$mouse_id %in% Samples,] #We do not have trait data for all samples
remove <- subset(Samples, !(Samples %in% traitData$mouse_id))#list of samples for which there is no trait data associated
# re-plot vo2 data ###
plot(traitData0$mass, traitData0$vo2, log="xy") #log=xy makes this a log-log plot
ln <- subset(traitData0, population=="LN")
me <- subset(traitData0, population=="ME")
plot(ln$mass, ln$vo2, xlim=c(12,30), ylim=c(1,7), pch=21, bg="darkorange", cex=2, log="xy")
points(me$mass, me$vo2, pch=21, bg="darkblue", cex=2)
#
mod <- nls(vo2~a*(mass^b), data=traitData0, start=list(a=1, b=1)) #fit a nonlinear regression for the mass metabolic rate relationship
lines(sort(traitData0$mass), predict(mod, list(mass=sort(traitData0$mass))), lwd=3, col="black")
##### make trait dataframe for WGCNA
traitx <- traitData0[, c(1:10)] #omit all traits except mass and vo2
traitx <- na.omit(traitx)
########Create a new Expr dataframe with missing samples removed
Expr1 = Expr[!row.names(Expr)%in%remove,]
traitRows <- match(rownames(Expr1), traitData$mouse_id) # make sure IDs are in the same order as Expr dataset.
Traits0 <- traitData[traitRows, -1]
rownames(Traits0) <- traitData[traitRows, 1]
Traits <- Traits0
dim(Traits)
traity <- traitx[traitx$mouse_id %in% rownames(Traits),]
dim(traity)
match(rownames(Traits), rownames(Expr1))# should be in order if datasets are matched
collectGarbage()
#
mod <- nls(vo2~a*(mass^b), data=traity, start=list(a=1, b=1)) #fit a nonlinear regression for the mass metabolic rate relationship
traity$vo2.res <- resid(mod)
voplot <- traity[, c("population", "acclim1", "vo2.res")]
voplot <- voplot[, c(1, 2, 3)]
vo.melt <- melt.data.frame(voplot, c("population", "acclim1"))
vo.melt$variable <- NULL
colnames(vo.melt) <- c("population", "acclim1", "vo2.res")
vo.melt <- na.omit(vo.melt)
#
vo.plot <- melt((tapply(vo.melt$vo2.res, list(vo.melt$population, vo.melt$acclim1),mean)))
names(vo.plot) <- c('population','acclim1',"vo2.res")
vo.sd <- melt(tapply(vo.melt$vo2.res, list(vo.melt$population, vo.melt$acclim1),function(x) sd(x)/sqrt(length(x))))
vo.plot <- data.frame(vo.plot,vo.sd[,3])
names(vo.plot)[4]='sem'
head(vo.plot)
#vo2 plot
# pdf("Figures/vo2.pdf", h=8, w=8)
lnCol <- "#DB9E06"
meCol <- "#263D92"
par(mfrow=c(1,1),oma=c(4,7,1,1),mar=c(1,0.5,1,0.5))
x1=vo.plot[vo.plot$population=='LN','acclim1']
y1=vo.plot[vo.plot$population=='LN','vo2.res']
plot(x1,y1,type='p',pch=19,col=lnCol,lwd=5,ylim=c(-1,1.5),ylab=' ',xlab=' ',xaxt='n',xlim=c(-.5,3.5),yaxt='n', cex.main=1.5)
errbar(x1,y1,y1+vo.plot[vo.plot$population=='LN','sem'],y1-vo.plot[vo.plot$population=='LN','sem'],add=T, errbar.col=lnCol, lwd=3)
points(x1, y1, bg=lnCol, col="black", pch=21, cex=3, lwd=2.5)
axis(2 ,cex.axis=2, las=2)
axis(1,at=c(0,1,2,3), labels=c("N","H","C","CH"),cex.axis=2)
x2=vo.plot[vo.plot$population=='ME','acclim1']
y2=vo.plot[vo.plot$population=='ME','vo2.res']
# lines(x2,y2, col=meCol, lwd=5)
errbar(x2,y2,y2+vo.plot[vo.plot$population=='ME','sem'],y2-vo.plot[vo.plot$population=='ME','sem'],add=T, errbar.col=meCol, lwd=3)
points(x2,y2,pch=21,cex=3, col='black', bg=meCol, lwd=2.5)
box(which = "plot", lty = "solid", lwd=2)
dev.off()

## vo2 re-analysis
anova(lmer(vo2~population*po2*temp + mass + (1|family), data=traitData0)) #three-way lmer
#posthoc analysis one-way pop-effect
N <- subset(traitData0, acclimation=="N"); anova(lmer(vo2~population +mass + (1|family), data=N))
H <- subset(traitData0, acclimation=="H"); anova(lmer(vo2~population + mass + (1|family), data=H))
C <- subset(traitData0, acclimation=="C"); anova(lmer(vo2~population + mass + (1|family), data=C))
CH <- subset(traitData0, acclimation=="CH"); anova(lmer(vo2~population + mass + (1|family), data=CH))

#Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(Expr, powerVector = powers, verbose = 5)
#pdf('wgcna/rv_beta_plot.pdf', h=4, w=7)
# dev.off()
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.9,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#dev.off()

# # construct gene network
cor <- WGCNA::cor
Net <- blockwiseModules(Expr, power = 7, maxBlockSize = dim(Expr)[2],TOMType = "signed", networkType = "signed", minModuleSize = 30,reassignThreshold = 0, mergeCutHeight = 0.25,numericLabels = TRUE, pamRespectsDendro = FALSE,saveTOMs = TRUE,saveTOMFileBase = "ExprTOM",verbose = 3)
save(Net,file = "thermo_capacity_gastroc_Network.RData")
load gastroc network
load(file = "thermo_capacity_gastroc_Network.RData")#load saved Network file
table(Net$colors)
moduleLabels = Net$colors
moduleColors = labels2colors(Net$colors)
MEs = Net$MEs;
geneTree = Net$dendrograms[[1]];
table(moduleColors)
dim(table(moduleColors)) #should be 37 modules
write.csv(as.data.frame(table(moduleColors)),file="gastroc_module_table.csv")
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(Net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(Net$dendrograms[[1]], mergedColors[Net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

#######################################################################################################################
# Trait Module Associations
gastroc.traits <- Traits[, c(8:11)] # subset data to only vo2 and o2 sat
names(gastroc.traits)
# Define numbers of genes and samples
nGenes = ncol(Expr1);
nSamples = nrow(Expr1)
MEs0 = moduleEigengenes(Expr1, moduleColors)$eigengenes
MEs = orderMEs(MEs0) #the rownames of this dataset are equal to Expr
match(rownames(gastroc.traits), rownames(Expr1))
moduleTraitCor = round(cor(MEs, gastroc.traits, use = "p"), 3)
moduleTraitPvalue = round(corPvalueStudent(moduleTraitCor, nSamples), 3)
moduleTraitQvalue <- round(matrix(p.adjust(as.vector(as.matrix(moduleTraitPvalue)), method='fdr'),ncol=4), 3)
colnames(moduleTraitQvalue) <- colnames(gastroc.traits)
rownames(moduleTraitQvalue) <- colnames(MEs)
moduleTrait <- moduleTraitQvalue[apply(moduleTraitQvalue, MARGIN = 1, function(x) any(x <= 0.05)), ]#table of modules with significant correlations after FDR correction
# write.csv(moduleTrait, file="Tables/gastroc_module_trait_cor_FDR.csv")

### Module trait correlation plots ###############################
MEs0 = moduleEigengenes(Expr1, moduleColors)$eigengenes
MEs = orderMEs(MEs0) #the rownames of this dataset are equal to Expr
METable0 <- MEs
METable <- METable0[colnames(METable0) %in% "MEdarkgrey"]
match(rownames(METable), rownames(Traits0))
METable$vo2 <- Traits0$vo2
METable$population <- Traits0$population
METable$acclimation <- Traits0$acclimation
METable1 <- melt.data.frame(METable, id.vars=c("population", "acclimation","vo2"))
colnames(METable1) <- c("population", "acclimation", "vo2","module", "ME")
module <- c("MEdarkgrey")
# pdf("Figures/gastroc_darkgrey_correlation.pdf", h=8, w=8)
par(mfrow=c(1,1),oma=c(4,4,1,1),mar=c(4,4,1,1), bg="white")
for (i in module){
  x1=METable1[METable1$module==i & METable1$population=="LN", "ME"]
  y1=METable1[METable1$module==i & METable1$population=="LN", "vo2"]
  plot(x1,y1, pch=21, xlim=c(-0.5,0.4), ylim = c(1.5,5.5), cex=2,ylab="", xlab="",col="black", bg="darkorange",lwd=1, main=i)
  x=METable1[METable1$module==i, "ME"]
  y=METable1[METable1$module==i, "vo2"]
  abline(lm(y~x))
  points(x1,y1, col="black", bg="darkorange", cex=2, pch=21)
  x2=METable1[METable1$module==i & METable1$population=="ME", "ME"]
  y2=METable1[METable1$module==i & METable1$population=="ME", "vo2"]
  points(x2,y2, col="black", bg="darkblue", cex=2, pch=21)
}
# dev.off()

### Module trait linar mixed effects models  ###############################
MEs0 = moduleEigengenes(Expr1, moduleColors)$eigengenes
MEs = orderMEs(MEs0) #the rownames of this dataset are equal to Expr
cor.table0 <- MEs
match(rownames(cor.table0), rownames(Traits))
cor.table0$population <- Traits$population
cor.table0$family <- Traits$family
cor.table0$mass <- Traits$mass
cor.table0$vo2 <- Traits$vo2
cor.table0$vo2_mass <- Traits$vo2_mass
cor.table0$arterial_sat <- Traits$arterial_sat
cor.table <- melt.data.frame(cor.table0, id.vars=c("population", "family", "mass", "vo2", "vo2_mass", "arterial_sat"))
colnames(cor.table) <- c("population", "family", "mass", "vo2", "vo2_mass", "arterial_sat", "module", "ME")
### Linear mixed effects model 
cor.lmer <- dlply(cor.table, .(module), lmer, formula = vo2~ME + mass+ (1|family:population))
cor.out <- lapply(cor.lmer, summary)
cor.out.table <- data.frame(matrix(unlist(cor.out), nrow=length(cor.out), byrow=T))[, c(42:43)]
rownames(cor.out.table) <- names(cor.out)
colnames(cor.out.table) <- c("p_vo2", "p_mass")
sig.cor <- subset(cor.out.table, p_vo2<0.05) #there are no LMER effects of module on vo2

#############################################################################
# Create the starting data frame
#######################################################################################################################
# Define variable weight containing the weight column of datTrait
vo2_mass = as.data.frame(gastroc.traits$vo2_mass)
names(vo2_mass) <- "vo2_mass"
MEs = Net$MEs
MEs0 = moduleEigengenes(Expr, moduleColors)$eigengenes
MEs = orderMEs(MEs0) #the rownames of this dataset are equal to Expr
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(Expr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM.", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(Expr1, vo2_mass, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(vo2_mass), sep="");
names(GSPvalue) = paste("p.GS.", names(vo2_mass), sep="");
genes=names(Expr)
GSPvalue$qvalue <- p.adjust(GSPvalue$p.GS.vo2_mass, method='fdr')
geneInfoGastroc = data.frame(Gene = genes,
                              geneTraitSignificance,
                              GSPvalue,
                              moduleColor = moduleColors,
                              geneModuleMembership,
                              MMPvalue)
head(geneInfoGastroc)
dim(geneInfoGastroc)
# write.csv(geneInfoGastroc, file="geneInfoGastroc.csv")
# 


############################################## LMER on ME ############
head(gastroc.id)
match(rownames(Expr), gastroc.id$mouse_id)
setSamples = rownames(Expr)
MEs0<-moduleEigengenes(Expr, moduleColors)$eigengenes #redo for all samples
MEs<-orderMEs(MEs0) #the rownames of this dataset are equal to Expr
ME<-MEs
match(rownames(ME), gastroc.id$mouse_id)
gas.table0 <- MEs
gas.table0$population <- as.factor(gastroc.id$population)
gas.table0$acclimation <- as.factor(gastroc.id$acclimation)
gas.table0$po2 <- as.factor(gastroc.id$po2)
gas.table0$temp <- as.factor(gastroc.id$temp)
gas.table0$family <- as.factor(gastroc.id$family)
gas.table <- melt.data.frame(gas.table0, id.vars=c("population", "acclimation","po2", "temp", "family"))
colnames(gas.table) <- c("population", "acclimation","po2", "temp", "family","module", "ME")
# model
model <- dlply(gas.table, .(module), lmer, formula = ME~population*po2*temp + (1|family:population)) #run linear mixed model
out <- lapply(model, anova) # output anova table from model results
out.table <- round(data.frame(matrix(unlist(out), nrow=length(out), byrow=T))[, c(36:42)],3) #create dataframe from output table
rownames(out.table) <- names(out)
colnames(out.table) <- c("population", "po2", "temp", "pop:po2", "pop:temp", "po2:temp", "3-way")
qTable <- round(matrix(p.adjust(as.vector(as.matrix(out.table)), method='fdr'),ncol=7), 3)
rownames(qTable) <- rownames(out.table)
colnames(qTable) <- colnames(out.table)
qTable.sig <- qTable[apply(qTable[, ], MARGIN = 1, function(x) any(x <= 0.05)), ]#table of modules with significance after FDR correction
write.csv(qTable.sig, file="Tables/gastroc_lmer_sig_pvalues.csv")
# f.out <- round(data.frame(matrix(unlist(out), nrow=length(out), byrow=T))[, c(29:35)],3) #table of F values
qPop <- qTable[apply(qTable[, -c(2:3, 6)], MARGIN = 1, function(x) any(x <= 0.05)), ]; dim(qPop) # at least one pop effect
qTreat <-  qTable[apply(qTable[, c(2:3, 6)], MARGIN = 1, function(x) any(x <= 0.05)), ]; dim(qTreat) # at least a treat effect
inter <- intersect(rownames(qPop), rownames(qTreat))
qTreat.only <- qTreat[!(rownames(qTreat) %in% inter), ]; dim(qTreat.only)

# MEdarkgrey
darkgrey <- subset(gas.table, module=="MEdarkgrey")
anova(lmer(ME~population*po2 + population*temp + po2*temp + (1|family:population), data=darkgrey)) #three-way lmer
#posthoc analysis one-way pop-effect
N <- subset(darkgrey, acclimation=="N"); anova(lmer(ME~population + (1|family:population), data=N))
H <- subset(darkgrey, acclimation=="H"); anova(lmer(ME~population + (1|family:population), data=H))
C <- subset(darkgrey, acclimation=="C"); anova(lmer(ME~population + (1|family:population), data=C))
CH <- subset(darkgrey, acclimation=="CH"); anova(lmer(ME~population + (1|family:population), data=CH))
#posthoc analysis one-way treat effect
me <- subset(darkgrey, population == "ME")
summary(lmer(ME~acclimation + (1|family), data=me))
me.warm <- subset(me, temp=="warm"); anova(lmer(ME~po2 + (1|family:population), data=me.warm)) #r2 = -0.04 N v H
me.C <- subset(me, acclimation == "C")
me.CH <- subset(me, acclimation == "CH")
me.N <- subset(me, acclimation== "N")
me.NC <- rbind(me.N, me.C); anova(lmer(ME~acclimation + (1|family:population), data=me.NC)) # -0.32 N v C
me.NCH <- rbind(me.N, me.CH); anova(lmer(ME~acclimation + (1|family:population), data=me.NCH)) # -0.65 N v CH
#
ln <- subset(darkgrey, population == "LN")
ln.warm <- subset(ln, temp=="warm"); anova(lmer(ME~po2 + (1|family:population), data=ln.warm)) #r2 = 0.08 N v H
ln.C <- subset(ln, acclimation == "C"); summary(lm(ME~po2, data=ln.C)) #r2 = 0.08 N v H
ln.CH <- subset(ln, acclimation == "CH")
ln.N <- subset(ln, acclimation== "N")
NC <- rbind(ln.N, ln.C); anova(lmer(ME~acclimation + (1|family:population), data=NC)) # -0.09 N v C
NCH <- rbind(ln.N, ln.CH); anova(lmer(ME~acclimation + (1|family:population), data=NCH)) # -0.21 N v CH

# MEred
red <- subset(gas.table, module=="MEred")
anova(lmer(ME~population*po2*temp + (1|family:population), data=red)) #three-way lmer
#posthoc analysis one-way
N <- subset(red, acclimation=="N"); anova(lmer(ME~population + (1|family:population), data=N))
H <- subset(red, acclimation=="H"); anova(lmer(ME~population + (1|family:population), data=H))
C <- subset(red, acclimation=="C"); anova(lmer(ME~population + (1|family:population), data=C))
CH <- subset(red, acclimation=="CH"); anova(lmer(ME~population + (1|family:population), data=CH))
me <- subset(red, population == "ME")
summary(lmer(ME~acclimation + (1|family), data=me))
me.warm <- subset(me, temp=="warm"); anova(lmer(ME~po2 + (1|family:population), data=me.warm)) #r2 = -0.04 N v H
me.C <- subset(me, acclimation == "C")
me.CH <- subset(me, acclimation == "CH")
me.N <- subset(me, acclimation== "N")
me.NC <- rbind(me.N, me.C); anova(lmer(ME~acclimation + (1|family:population), data=me.NC)) # -0.32 N v C
me.NCH <- rbind(me.N, me.CH); anova(lmer(ME~acclimation + (1|family:population), data=me.NCH)) # -0.65 N v CH
#
ln <- subset(red, population == "LN")
ln.warm <- subset(ln, temp=="warm"); anova(lmer(ME~po2 + (1|family:population), data=ln.warm)) #r2 = 0.08 N v H
ln.C <- subset(ln, acclimation == "C")
ln.CH <- subset(ln, acclimation == "CH")
ln.N <- subset(ln, acclimation== "N")
NC <- rbind(ln.N, ln.C); anova(lmer(ME~acclimation + (1|family:population), data=NC)) # -0.09 N v C
NCH <- rbind(ln.N, ln.CH); anova(lmer(ME~acclimation + (1|family:population), data=NCH)) # -0.21 N v CH

# MEbrown
brown <- subset(gas.table, module=="MEbrown")
anova(lmer(ME~population*po2*temp + (1|family:population), data=brown)) #three-way lmer
N <- subset(brown, acclimation=="N"); anova(lmer(ME~population + (1|family:population), data=N)) #posthoc analysis one-way
H <- subset(brown, acclimation=="H"); anova(lmer(ME~population + (1|family:population), data=H))
C <- subset(brown, acclimation=="C"); anova(lmer(ME~population + (1|family:population), data=C))
CH <- subset(brown, acclimation=="CH"); anova(lmer(ME~population + (1|family:population), data=CH))
# MEdarkgreen
darkgreen <- subset(gas.table, module=="MEdarkgreen")
anova(lmer(ME~population*po2*temp + (1|family:population), data=darkgreen)) #three-way lmer
N <- subset(darkgreen, acclimation=="N"); anova(lmer(ME~population + (1|family:population), data=N)) #posthoc analysis one-way
H <- subset(darkgreen, acclimation=="H"); anova(lmer(ME~population + (1|family:population), data=H))
C <- subset(darkgreen, acclimation=="C"); anova(lmer(ME~population + (1|family:population), data=C))
CH <- subset(darkgreen, acclimation=="CH"); anova(lmer(ME~population + (1|family:population), data=CH))

# MElightgreen
lightgreen <- subset(gas.table, module=="MElightgreen")
anova(lmer(ME~population*po2*temp + (1|family:population), data=lightgreen)) #three-way lmer
N <- subset(lightgreen, acclimation=="N"); anova(lmer(ME~population + (1|family:population), data=N)) #posthoc analysis one-way
H <- subset(lightgreen, acclimation=="H"); anova(lmer(ME~population + (1|family:population), data=H))
C <- subset(lightgreen, acclimation=="C"); anova(lmer(ME~population + (1|family:population), data=C))
CH <- subset(lightgreen, acclimation=="CH"); anova(lmer(ME~population + (1|family:population), data=CH))
boxplot(ME~population*acclimation, data=lightgreen)
me <- subset(lightgreen, population == "ME")
me.warm <- subset(me, temp=="warm"); anova(lmer(ME~po2 + (1|family:population), data=me.warm)) #r2 = -0.04 N v H
me.C <- subset(me, acclimation == "C")
me.CH <- subset(me, acclimation == "CH")
me.N <- subset(me, acclimation== "N")
me.NC <- rbind(me.N, me.C); anova(lmer(ME~acclimation + (1|family:population), data=me.NC)) # -0.32 N v C
me.NCH <- rbind(me.N, me.CH); anova(lmer(ME~acclimation + (1|family:population), data=me.NCH)) # -0.65 N v CH
#
ln <- subset(lightgreen, population == "LN")
ln.warm <- subset(ln, temp=="warm"); anova(lmer(ME~po2 + (1|family:population), data=ln.warm)) #r2 = 0.08 N v H
ln.C <- subset(ln, acclimation == "C")
ln.CH <- subset(ln, acclimation == "CH")
ln.N <- subset(ln, acclimation== "N")
NC <- rbind(ln.N, ln.C); anova(lmer(ME~acclimation + (1|family:population), data=NC)) # -0.09 N v C
NCH <- rbind(ln.N, ln.CH); anova(lmer(ME~acclimation + (1|family:population), data=NCH)) # -0.21 N v CH

# MEpurple
purple <- subset(gas.table, module=="MEpurple")
anova(lmer(ME~population*po2*temp + (1|family:population), data=purple)) #three-way lmer
N <- subset(purple, acclimation=="N"); anova(lmer(ME~population + (1|family:population), data=N)) #posthoc analysis one-way
H <- subset(purple, acclimation=="H"); anova(lmer(ME~population + (1|family:population), data=H))
C <- subset(purple, acclimation=="C"); anova(lmer(ME~population + (1|family:population), data=C))
CH <- subset(purple, acclimation=="CH"); anova(lmer(ME~population + (1|family:population), data=CH))

me <- subset(purple, population == "ME")
me.warm <- subset(me, temp=="warm"); anova(lmer(ME~po2 + (1|family:population), data=me.warm)) #r2 = -0.04 N v H
me.C <- subset(me, acclimation == "C")
me.CH <- subset(me, acclimation == "CH")
me.N <- subset(me, acclimation== "N")
me.NC <- rbind(me.N, me.C); anova(lmer(ME~acclimation + (1|family:population), data=me.NC)) # -0.32 N v C
me.NCH <- rbind(me.N, me.CH); anova(lmer(ME~acclimation + (1|family:population), data=me.NCH)) # -0.65 N v CH
#
ln <- subset(purple, population == "LN")
ln.warm <- subset(ln, temp=="warm"); anova(lmer(ME~po2 + (1|family:population), data=ln.warm)) #r2 = 0.08 N v H
ln.C <- subset(ln, acclimation == "C")
ln.CH <- subset(ln, acclimation == "CH")
ln.N <- subset(ln, acclimation== "N")
NC <- rbind(ln.N, ln.C); anova(lmer(ME~acclimation + (1|family:population), data=NC)) # -0.09 N v C
NCH <- rbind(ln.N, ln.CH); anova(lmer(ME~acclimation + (1|family:population), data=NCH)) # -0.21 N v CH

################## functional enrichment ################## 
convert <- read.csv("//Users//jonathanvelotta1//Dropbox//RWork//convert_light.csv")
convert <- unique(convert)
geneInfoConvert <- merge(convert, geneInfoGastroc, by="Gene")
dim(geneInfoConvert)
dim(geneInfoGastroc) # row dimensions should match
geneInfoConvert$moduleColor <- paste("ME", geneInfoConvert$moduleColor, sep='')
write.table(geneInfoConvert$Mus, file="enrichment/background_gastroc.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)#background for gastroc dataset
geneInfoConvert1 <- geneInfoConvert[geneInfoConvert$moduleColor %in% enrichmentMods,]
### This loops writes MUS IDs to a text file
# modules <- unique(geneInfoConvert$moduleColor)
# for (i in modules){
#   setwd("//Users//jonathanvelotta1//Dropbox//RWork//thermo_capacity//enrichment/gastroc")
#   x <- subset(geneInfoConvert, geneInfoConvert$moduleColor==i)
#   write.table(x$Mus, paste(i, "_mus.txt"), row.names=FALSE, col.names=paste(">", i), quote=FALSE)#add FASTA-style colnames to run multi-query in gProfiler
#   setwd("//Users//jonathanvelotta1//Dropbox//RWork//thermo_capacity")
# }
# combined <- c(as.vector(rownames(output.adjust)), as.vector(rownames(moduleTrait)))
# combined[duplicated(combined)]
# enrichmentMods <- unique(combined)#cat modules with sig modules from ANOVA and correlation test

##### Module eigengene plots #######
ME.melt <- MEs
ME.melt$ID <- rownames(MEs)
ME.melt <- melt(ME.melt)
colnames(ME.melt) <- c("mouse_id", "module", "ME")
df0 <- merge(ME.melt, gastroc.id, by="mouse_id")
#
plot <- melt((tapply(df0$ME,list(df0$population, df0$acclim1, df0$module),mean)))
names(plot) <- c('population','acclim1','module', "ME")
sd <- melt(tapply((df0$ME),list(df0$population,df0$acclim1, df0$module),function(x) sd(x)/sqrt(length(x))))
plot <- data.frame(plot,sd[,4])
names(plot)[5]='sem'
head(plot)
##### candidates
mods <- c("MEdarkgrey","MElightgreen", "MEpurple")
sub.plot <- plot[plot$module %in% mods,]
pdf("Figures/gastroc_darkgrey-lightgreen-purple_modules.pdf", h=15, w=6)
par(mfrow=c(3,1),oma=c(4,7,1,1),mar=c(1,0.5,1,0.5))
for(i in mods) {
  x1=plot[plot$module==i & plot$population=='LN','acclim1']
  y1=plot[plot$module==i & plot$population=='LN','ME']
  plot(x1,y1,type='l',pch=19,col='darkorange',lwd=5,ylim=c(-.3,.4),ylab=' ',xlab=' ',xaxt='n',xlim=c(-0.5,3.5),yaxt='n', main=i, cex.main=1.5)
  errbar(x1,y1,y1+plot[plot$module==i & plot$population=='LN','sem'],y1-plot[plot$module==i & plot$population=='LN','sem'],add=T, errbar.col="darkorange", lwd=3)
  #if(i==plot[2]){axis(2,las=1,cex.axis=1.5)}
  points(x1,y1,pch=21,cex=4, col='black', bg="darkorange", lwd=2.5)
  if(i==mods[1]){axis(2,at=c(-0.3, -0.1, 0.1, 0.3, 0.5) ,cex.axis=2, las=2)}
  if(i==mods[2]){axis(2,at=c(-0.3, -0.1, 0.1, 0.3, 0.5) ,cex.axis=2, las=2)}
  if(i==mods[3]){axis(2,at=c(-0.3, -0.1, 0.1, 0.3, 0.5) ,cex.axis=2, las=2)}
  if(i==mods[3]){axis(1,at=c(0,1,2,3), labels=c("N","H","C","CH"),cex.axis=2)}
  # if(i==mods[2]){axis(1,at=c(0,1,2,3), labels=c("N","H","C","CH"),cex.axis=2)}
  x2=jitter(plot[plot$module==i & plot$population=='ME','acclim1'])
  y2=plot[plot$module==i & plot$population=='ME','ME']
  errbar(x2,y2,y2+plot[plot$module==i & plot$population=='ME','sem'],y2-plot[plot$module==i & plot$population=='ME','sem'],add=T, errbar.col="darkblue", lwd=3)
  points(x2,y2,pch=21,cex=4, col='black', bg="darkblue", lwd=2.5)
  lines(x2,y2, col='darkblue', lwd=5)
  box(which = "plot", lty = "solid", lwd=2)
}
dev.off()

# enrichmenet dotplot
darkgrey.go <- read.csv("enrichment/gastroc_darkgrey_enrichment.csv")
darkgrey.kegg <- subset(darkgrey.go, source=="KEGG")[,c(2:3,5,8)]
darkgrey.kegg <- darkgrey.kegg[order(darkgrey.kegg$negative_log10_of_adjusted_p_value),]
pdf("Figures/darkgrey.KEGG.pdf", h=6, w=5)
scale <- darkgrey.kegg$intersection_size/2
dotchart(darkgrey.kegg$negative_log10_of_adjusted_p_value, pt.cex=scale, labels=darkgrey.kegg$term_name, pch=21, col="black", bg="grey", xlim=c(1,7), lcolor="lightgrey", xlab="-log(p-value)")
dev.off()

# enrichmenet dotplot
lightgreen.go <- read.csv("enrichment/gastroc_lightgreen_enrichment.csv")
lightgreen.kegg <- subset(lightgreen.go, source=="KEGG")[,c(2:3,5,8)]
lightgreen.kegg <- lightgreen.kegg[order(lightgreen.kegg$negative_log10_of_adjusted_p_value),]
pdf("Figures/lightgreen.KEGG.pdf", h=6, w=6)
scale <- lightgreen.kegg$intersection_size/10
dotchart(lightgreen.kegg$negative_log10_of_adjusted_p_value, pt.cex=scale, labels=lightgreen.kegg$term_name, pch=21, col="black", bg="grey", xlim=c(1,50), lcolor="lightgrey", xlab="-log(p-value)")
dev.off()


# enrichmenet dotplot
purple.go <- read.csv("enrichment/gastroc_purple_enrichment.csv")
purple.kegg <- subset(purple.go, source=="KEGG")[,c(2:3,5,8)]
purple.kegg <- purple.kegg[order(purple.kegg$negative_log10_of_adjusted_p_value),]
pdf("Figures/purple.KEGG.pdf", h=6, w=3)
scale <- purple.kegg$intersection_size/2
dotchart(purple.kegg$negative_log10_of_adjusted_p_value, pt.cex=3, labels=purple.kegg$term_name, pch=21, col="black", bg="grey", xlim=c(1,15), lcolor="lightgrey", xlab="-log(p-value)")
dev.off()


########################### network figure
############################### Cytoscape Network
MEs0<-moduleEigengenes(Expr, moduleColors)$eigengenes #redo for all samples
MEs<-orderMEs(MEs0) #the rownames of this dataset are equal to Expr
ME<-MEs
me.matrix <- data.matrix(MEs) # convert ME data frame to a data matrix
ME_cor <- cor(me.matrix) # compute pearson correlation for ME matrix
ME_cor <- melt(ME_cor)
names(ME_cor) <- c("node1", "node2", "r")
pvals <- rcorr(me.matrix)$P
pvals <- melt(pvals)
names(pvals) <- c("node1", "node2", "pvalue")
ME_cor$pvalue <- pvals$pvalue #add pvalues to the correlations
ME_cor <- na.omit(ME_cor) #omit the 1:1 correlations and pvalues
cor <- subset(ME_cor,!duplicated(r))
cor$abs <- abs(cor$r)
cor$abs <- as.numeric(cor$abs)
labels <- data.frame(table(moduleColors))
labels <- labels[! (labels$moduleColors=="grey"),] #remove the grey module
labels$node1 <- as.factor(paste("ME", labels$moduleColors, sep=""))
cor1 <- merge(cor, labels, by="node1")
cor1$q <- p.adjust(cor1$pvalue, method="bonferroni")
cyto <- subset(cor1, q<=0.05)
cyto <- cyto[order(cyto$node1),]
dim(cyto)
links <- cyto
names(links) <- c("from", "to", "r", "pvalue", "abs", "moduleColors", "Freq", "q")
links <- links[order(links$abs, decreasing = TRUE),]
nodes <- labels[,c(3,1,2)]
names(nodes) <- c("id", "moduleColors", "Freq")

# igraph code
# net <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
# net
# plot(net)
# net <- simplify(net, remove.multiple = F, remove.loops = T)
# 
# plot(net, edge.arrow.size=.4,vertex.label=NA)
# V(net)$color <- as.character(nodes$moduleColors)
# # Compute node degrees (#links) and use that to set node size:
# deg <- degree(net, mode="all")
# V(net)$size <- deg*3
# # We could also use the size value:
# V(net)$size <- log(nodes$Freq)*5
# # The labels are currently node IDs.
# # Setting them to NA will render no labels:
# V(net)$label <- as.vector(nodes$moduleColors)
# # Set edge width based on weight:
# #E(net)$width <- ifelse(links$Pvalues > 0.05,0,abs(links$Correlation*5))
# E(net)$width <- abs(links$r*10)
# #change arrow size and edge color:
# E(net)$arrow.size <- 1
# E(net)$edge.color <- "black"
# # Set the network layout:
# par(mfrow=c(1,1))
# graph_attr(net, "layout") <- layout_with_lgl
# plot(net)
# ####Explore different Layouts
# layouts <- grep("^layout_", ls("package:igraph"), value=TRUE)[-1] 
# # Remove layouts that do not apply to our graph.
# layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama|tree", layouts)]
# V(net)$frame.color <- "black"
# par(mfrow=c(3,3), mar=c(1,1,1,1))
# for (layout in layouts) {
#   print(layout)
#   l <- do.call(layout, list(net))
#   plot(net, edge.arrow.mode=0, layout=l, main=layout) }

# Choose final layout
net <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
net <- simplify(net, remove.multiple = F, remove.loops = T)
# V(net)$size <- log(nodes$Freq)*5
V(net)$color <- as.character(nodes$moduleColors)
V(net)$label <- ''
E(net)$width <- abs(links$r*7)
V(net)$label <- as.vector(nodes$moduleColors)
V(net)$frame.color <- "black"
#change arrow size and edge color:
E(net)$arrow.size <- 0.5
E(net)$edge.color <- "black"
# save(net_groups, file="net_groups_gastroc_network.RData")
# load(file="net_groups_gastroc_network.RData")
# coords <- layout_in_circle(net, order=order(membership(net_groups)))
pdf("Figures/gastroc_igraph.pdf", h=10, w=10)
# plot(net, layout=coords)
plot(net, layout=layout_with_mds)
dev.off()

#barplot of modules sizes
pdf("Figures/gastroc_mod_size.pdf", h=10, w=10)
nodes1 <- nodes[order(nodes$Freq, decreasing = TRUE), ]
barplot(nodes1$Freq, col = as.character(nodes1$moduleColors), ylim = c(0,1000))
dev.off()


################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
#lung analysis
lung.id <- subset(id, tissue=="lung") #subset ids by tissue
lung.id <- lung.id[! (lung.id$mouse_id=="LN.F1.113.am"),]
lung.id <- lung.id[! (lung.id$mouse_id=="LN.F1.113.b"),]
lung.id <- lung.id[! (lung.id$mouse_id=="LN.F1.132.af"),]
lung.id1 <- lung.id$ID #make a vector containing lung mouse ids
lung <- thermo[,colnames(thermo) %in% lung.id1] #subset ontogeny data by lung mouse ids
head(lung)#only lung data
dim(lung)
match(lung.id$ID, colnames(lung)) #verify that colnames and sample id names are in the same order
colnames(lung) <- lung.id$mouse_id ###WARNING: this step removes tissue identifier in ID names. Proceed with caution.
match(lung.id$mouse_id, colnames(lung))
# write.csv(lung, "lung_raw_counts_thermocapacity.csv")
########
population <- lung.id$population
acclimation <- lung.id$acclimation
table <- data.frame(population, acclimation)
group <- factor(paste(table$population, table$acclimation, sep="_"))
cbind(table, group=group)
table$population = as.factor(table$population)
table$acclimation <- relevel(table$acclimation, ref="N")
table$population <- relevel(table$population, ref="LN")
design <- model.matrix(~population*acclimation, data=table)
#filter data
lung$mean = rowMeans(lung) #rowMeans takes means of each row
keep_lung = subset(lung, mean >= 10) #filter by means
dim(keep_lung)
keep_lung$mean = NULL #clean up dataset
dim(keep_lung)
y0 <- DGEList(counts=keep_lung, group=group) # make a DGE list
# dim(y0)
# y <- calcNormFactors(y0) # normalize
# plotMDS(y) #exploration plots
# y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
# y <- estimateGLMTrendedDisp(y, design)
# y <- estimateGLMTagwiseDisp(y, design)
lung.norm <- cpm(y0, log=TRUE, prior.count=1, normalized.lib.sizes=TRUE) #cpm normalized and log transformed expression dataset
#write.csv(lung.norm, file="thermo_lung_norm_counts.csv")
plotMDS(lung.norm)
dim(lung.norm)

pca <- prcomp(t(lung.norm), scale=FALSE)
summary(pca)
pc = as.data.frame(pca$x)
pc$population <- population
pc$acclimation <- acclimation
ln <- subset(pc, population == "LN")
lnN <- subset(ln, acclimation=="N")
lnH <- subset(ln, acclimation=="H")
lnC <- subset(ln, acclimation=="C")
lnCH <- subset(ln, acclimation=="CH")
me <- subset(pc, population == "ME")
meN <- subset(me, acclimation=="N")
meH <- subset(me, acclimation=="H")
meC <- subset(me, acclimation=="C")
meCH <- subset(me, acclimation=="CH")
pdf("Figures/pca_lung.pdf", h=8, w=8)
par(mfrow=c(1,1),oma=c(4,4,1,1),mar=c(4,4,1,1))
plot(lnN$PC1,lnN$PC2, type='n', pch=19, cex=2, xlim=c(-50, 50), ylim=c(-60, 70), yaxt='n', xaxt='n',ylab="", xlab="",col="black", lwd=2, bg=as.numeric(pc$population), bty='n')
axis(side=1, lwd=2, cex.axis=2, las=1)
mtext("PC1", side=1, line=4, cex=2)
axis(side=2, lwd=2, cex.axis=2, las=1)
mtext("PC2", side=4, cex=2, line=4)
ordiellipse(pca,pc$population,conf=0.95, draw="polygon", col="gray92", border = "white", lwd=5)
points(lnN$PC1, lnN$PC2, pch=21, col="darkred", bg="darkred",cex=2.5, lwd=2)
points(meN$PC1, meN$PC2, pch=21, col="darkred", bg="darkred", cex=2.5, lwd=2)
points(lnH$PC1, lnH$PC2, pch=21, col="red", bg="red", cex=2.5, lwd=2)
points(meH$PC1, meH$PC2, pch=21, col="red", bg="red", cex=2.5, lwd=2)
points(lnC$PC1, lnC$PC2, pch=21, col="darkblue", bg="darkblue", cex=2.5, lwd=2)
points(meC$PC1, meC$PC2, pch=21, col="darkblue", bg="darkblue", cex=2.5, lwd=2)
points(lnCH$PC1, lnCH$PC2, pch=21, col="blue", bg="blue", cex=2.5, lwd=2)
points(meCH$PC1, meCH$PC2, pch=21, col="blue", bg="blue", cex=2.5, lwd=3)
box(which="plot", lty="solid", lwd=2)
dev.off()

####################################### Lung WGCNA #####################################################
##############################################################################################################################################
head(lung.norm)#normalized read counts for lung only
dim(lung.norm)
Expr0 = as.data.frame(t(lung.norm)) #transpose expression data for further analysis
rownames(Expr0) = colnames(lung.norm)
gsg = goodSamplesGenes(Expr0, verbose = 3) #check for genes and samples with too many missing values
gsg$allOK
Expr = Expr0
#cluster the samples to check for outliers
sampleTree = hclust(dist(Expr), method = "average");
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 2.5,
     cex.axis = 1.5, cex.main = 2)

###Import Trait Data
traitData0 <- read.csv("thermo_capacity_traits_all.csv")#all trait daeta
dev.off()
# boxplot(vo2_mass~population*acclimation, data=traitData0)
plot(traitData0$mass, traitData0$vo2)
cor.test(traitData0$mass, traitData0$vo2)
ln <- subset(traitData0, population=="LN")
me <- subset(traitData0, population=="ME")
plot(ln$mass, ln$vo2, xlim=c(12,30), ylim=c(1,7), pch=21, bg="darkorange", cex=2)
points(me$mass, me$vo2, pch=21, bg="darkblue", cex=2)
traitx <- traitData0[, c(1:10)]
traitx <- na.omit(traitx)
vo2.res <- resid(lm(vo2~mass, data=traitx))
traitx$vo2.res <- vo2.res
boxplot(vo2.res~acclimation*population, data=traitx)
# 
Samples <- rownames(Expr)
traitData <- traitData0[traitData0$mouse_id %in% Samples,]#We do not have trait data for all samples
remove <- subset(Samples, !(Samples %in% traitData$mouse_id))#list of samples for which there is no trait data associatedrow
# ########Create a new Expr dataframe with missing samples removed
Expr1 = Expr[!row.names(Expr)%in%remove,]
# traitRows <- match(rownames(Expr1), traitData$mouse_id) # make sure IDs are in the same order as Expr dataset.
# Traits0 <- traitData[traitRows, -1]
# rownames(Traits0) <- traitData[traitRows, 1]
# Traits <- Traits0
# traity <- traitx[traitx$mouse_id %in% rownames(Traits),]
# match(traity$mouse_id, rownames(Traits))
# Traits$vo2.res <- traity$vo2.res
# collectGarbage();
# write.csv(Traits, file="lung_trait_data.csv")
# write.csv(traity, file="lung_vo2_resid_data.csv")
Traits <- read.csv("lung_trait_data.csv", row.names = 1) #this dataset contains all data associated with each lung RNA sample
match(rownames(Traits), rownames(Expr1))# should be in order if datasets are matched
collectGarbage()
#
# # Choose a set of soft-thresholding powers
# powers = c(c(1:10), seq(from = 12, to=20, by=2))
# sft = pickSoftThreshold(Expr, powerVector = powers, verbose = 5)
# # pdf('wgcna/rv_beta_plot.pdf', h=4, w=7)
# # dev.off()
# par(mfrow = c(1,2));
# cex1 = 0.9;
# plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
#      main = paste("Scale independence"));
# text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      labels=powers,cex=cex1,col="red");
# # this line corresponds to using an R^2 cut-off of h
# abline(h=0.9,col="red")
# plot(sft$fitIndices[,1], sft$fitIndices[,5],
#      xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
#      main = paste("Mean connectivity"))
# text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# #dev.off()
# 
# # construct gene network
# cor <- WGCNA::cor
# Net <- blockwiseModules(Expr, power = 7, maxBlockSize = dim(Expr)[2],TOMType = "signed", networkType = "signed", minModuleSize = 30,reassignThreshold = 0, mergeCutHeight = 0.25,numericLabels = TRUE, pamRespectsDendro = FALSE,saveTOMs = TRUE,saveTOMFileBase = "ExprTOM",verbose = 3)
# save(Net,file = "thermo_capacity_lung_Network.RData")
load(file = "thermo_capacity_lung_Network.RData")#load saved Network file
table(Net$colors)
moduleLabels = Net$colors
moduleColors = labels2colors(Net$colors)
length(moduleColors)
MEs = Net$MEs;
geneTree = Net$dendrograms[[1]];
table(moduleColors)
dim(table(moduleColors)) #63 modules
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(Net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(Net$dendrograms[[1]], mergedColors[Net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

#######################################################################################################################
# Trait Module Associations
lung.traits <- Traits[, c(9:11, 15, 17, 19, 21)] # subset data to only vo2 and o2 sat
names(lung.traits)
# Define numbers of genes and samples
nGenes = ncol(Expr1);
nSamples = nrow(Expr1)
MEs0 = moduleEigengenes(Expr1, moduleColors)$eigengenes
MEs = orderMEs(MEs0) #the rownames of this dataset are equal to Expr
match(rownames(lung.traits), rownames(Expr1))
moduleTraitCor = round(cor(MEs, lung.traits, use = "p"), 3)
moduleTraitPvalue = round(corPvalueStudent(moduleTraitCor, nSamples), 3)
moduleTraitQvalue <- round(matrix(p.adjust(as.vector(as.matrix(moduleTraitPvalue)), method='fdr'),ncol=dim(moduleTraitPvalue)[2]), 3)
colnames(moduleTraitQvalue) <- colnames(lung.traits)
rownames(moduleTraitQvalue) <- colnames(MEs)
moduleTraitP <- moduleTraitPvalue[apply(moduleTraitPvalue, MARGIN = 1, function(x) any(x <= 0.05)), ]#table of modules with significant correlations before FDR correction
moduleTraitQ <- moduleTraitQvalue[apply(moduleTraitQvalue, MARGIN = 1, function(x) any(x <= 0.05)), ]#table of modules with significant correlations after FDR correction

### Module trait correlation plots ###############################
MEs0 = moduleEigengenes(Expr1, moduleColors)$eigengenes
MEs = orderMEs(MEs0) #the rownames of this dataset are equal to Expr
match(rownames(lung.traits), rownames(Expr1))
METable <- MEs
METable$vo2 <- Traits$vo2
METable$population <- Traits$population
METable$acclimation <- Traits$acclimation
METable1 <- melt.data.frame(METable, id.vars=c("population", "acclimation","vo2"))
colnames(METable1) <- c("population", "acclimation", "vo2","module", "ME")
module <- "MEdarkred"
pdf("Figures/lung_darkred_vo2_ME_correlation.pdf", h=8, w=8)
par(mfrow=c(1,1),oma=c(4,4,1,1),mar=c(4,4,1,1), bg="white")
for (i in module){
  x1=METable1[METable1$module==i & METable1$population=="LN", "ME"]
  y1=METable1[METable1$module==i & METable1$population=="LN", "vo2"]
  plot(x1,y1, pch=21, xlim=c(-0.3,0.6), ylim = c(1.5, 7), cex=2,ylab="", xlab="",col="black", bg="darkorange",lwd=1, main=i)
  points(x1,y1, col="black", bg="darkorange", cex=2, pch=21)
  x=METable1[METable1$module==i, "ME"]
  y=METable1[METable1$module==i, "vo2"]
  abline(lm(y~x))
  x2=METable1[METable1$module==i & METable1$population=="ME", "ME"]
  y2=METable1[METable1$module==i & METable1$population=="ME", "vo2"]
  points(x2,y2, col="black", bg="darkblue", cex=2, pch=21)
}
dev.off()

darkred <- subset(METable1, module=="MEdarkred")
test <- MEs[, c(61, 62)]

### Module trait linar mixed effects models  ###############################
cor.table0 <- MEs
match(rownames(cor.table0), rownames(Traits))
cor.table0$population <- Traits$population
cor.table0$family <- Traits$family
cor.table0$mass <- Traits$mass
cor.table0$vo2 <- Traits$vo2
cor.table0$vo2_mass <- Traits$vo2_mass
cor.table0$arterial_sat <- Traits$arterial_sat
cor.table <- melt.data.frame(cor.table0, id.vars=c("population", "family", "mass", "vo2", "vo2_mass", "arterial_sat"))
colnames(cor.table) <- c("population", "family", "mass", "vo2", "vo2_mass", "arterial_sat", "module", "ME")
### Linear mixed effects model for ME on VO2
cor.lmer <- dlply(cor.table, .(module), lmer, formula = vo2~ME + mass + (1|family:population))
cor.out <- lapply(cor.lmer, anova)
cor.out.table <- data.frame(matrix(unlist(cor.out), nrow=length(cor.out), byrow=T))[, c(11:12)]
rownames(cor.out.table) <- names(cor.out)
colnames(cor.out.table) <- c("p_vo2", "p_mass")
cor.out.table$q_vo2 <- p.adjust(cor.out.table$p_vo2, method="fdr")
sig.cor.vo2 <- subset(cor.out.table, q_vo2<0.05) # a few modules are significantly associated with VO2 when taking family into account
### Linear mixed effects model for ME on arterial o2 sat
cor.lmer <- dlply(cor.table, .(module), lmer, formula = arterial_sat ~ ME + (1|family:population))
cor.out <- lapply(cor.lmer, anova)
cor.out.table <- data.frame(matrix(unlist(cor.out), nrow=length(cor.out), byrow=T))[, c(5:6)]
rownames(cor.out.table) <- names(cor.out)
colnames(cor.out.table) <- c("F", "p_sat")
cor.out.table$q_sat <- p.adjust(cor.out.table$p_sat, method="fdr")
sig.cor.sat <- subset(cor.out.table, q_sat<0.05) # a few modules are significantly associated with VO2 when taking family into account

#############################################################################
# Create the starting data frame
#######################################################################################################################
MEs = Net$MEs
MEs0 = moduleEigengenes(Expr, moduleColors)$eigengenes
MEs = orderMEs(MEs0) #the rownames of this dataset are equal to Expr
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(Expr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM.", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
#geneTraitSignificance = as.data.frame(cor(Expr, osmo, use = "p"));
#GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
#names(geneTraitSignificance) = paste("GS.", names(osmo), sep="");
#names(GSPvalue) = paste("p.GS.", names(osmo), sep="");
genes=names(Expr)
geneInfoLung = data.frame(Gene = genes,
                      moduleColor = moduleColors,
                      geneModuleMembership,
                      MMPvalue)
dim(geneInfoLung)
############################################## LMER on ME ############
head(lung.id)
match(rownames(Expr), lung.id$mouse_id)
setSamples = rownames(Expr)
MEs0<-moduleEigengenes(Expr, moduleColors)$eigengenes #redo for all samples
MEs<-orderMEs(MEs0) #the rownames of this dataset are equal to Expr
ME<-MEs
match(rownames(ME), lung.id$mouse_id)
lung.table0 <- MEs
lung.table0$population <- as.factor(lung.id$population)
lung.table0$acclimation <- as.factor(lung.id$acclimation)
lung.table0$po2 <- as.factor(lung.id$po2)
lung.table0$temp <- as.factor(lung.id$temp)
lung.table0$family <- as.factor(lung.id$family)
lung.table <- melt.data.frame(lung.table0, id.vars=c("population", "acclimation","po2", "temp", "family"))
colnames(lung.table) <- c("population", "acclimation","po2", "temp", "family","module", "ME")
# model
model <- dlply(lung.table, .(module), lmer, formula = ME~population*po2*temp + (1|family:population)) #run linear mixed model
out <- lapply(model, anova) # output anova table from model results
# write.csv(out, "lung_lmer_output.csv")
out.table <- round(data.frame(matrix(unlist(out), nrow=length(out), byrow=T))[, c(36:42)],3) #create dataframe from output table
# tail(out.table)
rownames(out.table) <- names(out)
colnames(out.table) <- c("population", "po2", "temp", "pop:po2", "pop:temp", "po2:temp", "3-way")
qTable <- round(matrix(p.adjust(as.vector(as.matrix(out.table)), method='fdr'),ncol=7), 3)
rownames(qTable) <- rownames(out.table)
colnames(qTable) <- colnames(out.table)
qTable.sig <- qTable[apply(qTable[, ], MARGIN = 1, function(x) any(x <= 0.05)), ]#table of modules with significance after FDR correction
qTable.sig <- qTable.sig[-18, ] #remove grey
write.csv(qTable.sig, file="Lung_lmer_pvals_table.csv")
#lightcyan, angiogensis, higher ME esp in Cold
#re check blue and brown
#lightsteelblue1 glycolosis, down in LN in H
#midnightblue, immune system, higher in LN in cold

# MEdarkred, weak population effect
darkred <- subset(lung.table, module=="MEdarkred")
anova(lmer(ME~population*po2*temp + (1|family:population), data=darkred)) #three-way lmer
N <- subset(darkred, acclimation=="N"); anova(lmer(ME~population + (1|family:population), data=N)) #posthoc analysis one-way
H <- subset(darkred, acclimation=="H"); anova(lmer(ME~population + (1|family:population), data=H))
C <- subset(darkred, acclimation=="C"); anova(lmer(ME~population + (1|family:population), data=C))
CH <- subset(darkred, acclimation=="CH"); anova(lmer(ME~population + (1|family:population), data=CH))
me <- subset(darkred, population == "ME")
me.warm <- subset(me, temp=="warm"); anova(lmer(ME~po2 + (1|family:population), data=me.warm)) #r2 = -0.04 N v H
me.C <- subset(me, acclimation == "C")
me.CH <- subset(me, acclimation == "CH")
me.N <- subset(me, acclimation== "N")
me.NC <- rbind(me.N, me.C); anova(lmer(ME~acclimation + (1|family:population), data=me.NC)) # -0.32 N v C
me.NCH <- rbind(me.N, me.CH); anova(lmer(ME~acclimation + (1|family:population), data=me.NCH)) # -0.65 N v CH

LN <- subset(darkred, population == "LN")
LN.warm <- subset(LN, temp=="warm"); anova(lmer(ME~po2 + (1|family:population), data=LN.warm)) #r2 = -0.04 N v H
LN.C <- subset(LN, acclimation == "C")
LN.CH <- subset(LN, acclimation == "CH")
LN.N <- subset(LN, acclimation== "N")
LN.NC <- rbind(LN.N, LN.C); anova(lmer(ME~acclimation + (1|family:population), data=LN.NC)) # -0.32 N v C
LN.NCH <- rbind(LN.N, LN.CH); anova(lmer(ME~acclimation + (1|family:population), data=LN.NCH)) # -0.65 N v CH


palevioletred3 <- subset(lung.table, module=="MEpalevioletred3")
anova(lmer(ME~population*po2*temp + (1|family:population), data=palevioletred3)) #three-way lmer
N <- subset(palevioletred3, acclimation=="N"); anova(lmer(ME~population + (1|family:population), data=N)) #posthoc analysis one-way
H <- subset(palevioletred3, acclimation=="H"); anova(lmer(ME~population + (1|family:population), data=H))
C <- subset(palevioletred3, acclimation=="C"); anova(lmer(ME~population + (1|family:population), data=C))
CH <- subset(palevioletred3, acclimation=="CH"); anova(lmer(ME~population + (1|family:population), data=CH))


##### Module eigengene plots #######
ME.melt <- MEs
ME.melt$ID <- rownames(MEs)
ME.melt <- melt(ME.melt)
colnames(ME.melt) <- c("mouse_id", "module", "ME")
df0 <- merge(ME.melt, lung.id, by="mouse_id")
plot <- melt((tapply(df0$ME,list(df0$population, df0$acclim1, df0$module),mean)))
names(plot) <- c('population','acclim1','module', "ME")
sd <- melt(tapply((df0$ME),list(df0$population,df0$acclim1, df0$module),function(x) sd(x)/sqrt(length(x))))
plot <- data.frame(plot,sd[,4])
names(plot)[5]='sem'
mods <- c("MEdarkred", "MEpalevioletred3")
sub.plot <- plot[plot$module %in% mods,]
module = unique(sub.plot$module)
pdf("Figures/lung_modules.pdf", h=8, w=8)
par(mfrow=c(3,1),oma=c(4,7,1,1),mar=c(2,2,3,2))
for(i in module) {
  x1=plot[plot$module==i & plot$population=='LN','acclim1']
  y1=plot[plot$module==i & plot$population=='LN','ME']
  plot(x1,y1,type='l',pch=19,col='darkorange',lwd=3,ylim=c(-.3,.4),ylab=' ',xlab=' ',xaxt='n',xlim=c(-0.5,3.5),yaxt='n', main=i, cex.main=1.5, las=1)
  errbar(x1,y1,y1+plot[plot$module==i & plot$population=='LN','sem'],y1-plot[plot$module==i & plot$population=='LN','sem'],add=T, errbar.col="darkorange", lwd=3)
  #if(i==plot[2]){axis(2,las=1,cex.axis=1.5)}
  points(x1,y1,pch=21,cex=3, col='black', bg="darkorange", lwd=1.5)
  axis(1,at=c(0,1,2,3), labels=c("N","H","C","CH"),cex.axis=1.5)
  axis(2,at=c(-0.3, -0.1, 0.1, 0.3) ,cex.axis=1.5, las=1)
  x2=jitter(plot[plot$module==i & plot$population=='ME','acclim1'])
  y2=plot[plot$module==i & plot$population=='ME','ME']
  errbar(x2,y2,y2+plot[plot$module==i & plot$population=='ME','sem'],y2-plot[plot$module==i & plot$population=='ME','sem'],add=T, errbar.col="darkblue", lwd=3)
  points(x2,y2,pch=21,cex=3, col='black', bg="darkblue", lwd=1.5)
  lines(x2,y2, col='darkblue', lwd=3)
}
dev.off()

# enrichmenet dotplot
darkred.go <- read.csv("enrichment/lung_darkred_enrichment.csv")
darkred.kegg <- subset(darkred.go, source=="REAC")[,c(2:3,5,8)]
darkred.kegg <- darkred.kegg[order(darkred.kegg$negative_log10_of_adjusted_p_value),]
pdf("Figures/darkred.REAC.pdf", h=6, w=7)
scale <- darkred.kegg$intersection_size/2
dotchart(darkred.kegg$negative_log10_of_adjusted_p_value, pt.cex=scale, labels=darkred.kegg$term_name, pch=21, col="black", bg="grey", xlim=c(1,5), lcolor="lightgrey", xlab="-log(p-value)")
dev.off()

################## functional enrichment ################## 
convert <- read.csv("//Users//jonathanvelotta1//Dropbox//RWork//convert_light.csv")
convert <- unique(convert)
geneInfoConvert <- merge(geneInfoLung, convert, by="Gene", all.y=FALSE)
dim(geneInfoConvert)
dim(geneInfoLung) # row dimensions should match 
geneInfoConvert$moduleColor <- paste("ME", geneInfoConvert$moduleColor, sep='')
write.table(geneInfoConvert$Mus, file="enrichment/background_lung.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)#background for lung dataset
# geneInfoConvert1 <- geneInfoConvert[geneInfoConvert$moduleColor %in% enrichmentMods,]
### This loops writes MUS IDs to a text file
modules <- unique(geneInfoConvert$moduleColor)
for (i in modules){
  x <- subset(geneInfoConvert, geneInfoConvert$moduleColor==i)
  write.table(x$Mus, paste(i,"_lung_mus.txt"), row.names=FALSE, col.names=paste(">", i), quote=FALSE)#add FASTA-style colnames to run multi-query in gProfiler
}

########################### network figure
############################### Cytoscape Network
head(MEs)
me.matrix <- data.matrix(MEs) # convert ME data frame to a data matrix
ME_cor <- cor(me.matrix) # compute pearson correlation for ME matrix
ME_cor <- melt(ME_cor)
names(ME_cor) <- c("node1", "node2", "r")
pvals <- rcorr(me.matrix)$P
pvals <- melt(pvals)
names(pvals) <- c("node1", "node2", "pvalue")
ME_cor$pvalue <- pvals$pvalue #add pvalues to the correlations
ME_cor <- na.omit(ME_cor) #omit the 1:1 correlations and pvalues
cor <- subset(ME_cor,!duplicated(r))
cor$abs <- abs(cor$r)
cor$abs <- as.numeric(cor$abs)
labels <- data.frame(table(moduleColors))
labels <- labels[! (labels$moduleColors=="grey"),] #remove the grey module
labels$node1 <- as.factor(paste("ME", labels$moduleColors, sep=""))
cor1 <- merge(cor, labels, by="node1")
cor1$q <- p.adjust(cor1$pvalue, method="bonferroni")
cyto <- subset(cor1, q<=0.05)
cyto <- cyto[order(cyto$node1),]
dim(cyto)

links <- cyto
names(links) <- c("from", "to", "r", "pvalue", "abs", "moduleColors", "Freq", "q")
links <- links[order(links$abs, decreasing = TRUE),]
nodes <- labels[,c(3,1,2)]
names(nodes) <- c("id", "moduleColors", "Freq")

########################### network figure
############################### Cytoscape Network
MEs0<-moduleEigengenes(Expr, moduleColors)$eigengenes #redo for all samples
MEs<-orderMEs(MEs0) #the rownames of this dataset are equal to Expr
ME<-MEs
me.matrix <- data.matrix(MEs) # convert ME data frame to a data matrix
ME_cor <- cor(me.matrix) # compute pearson correlation for ME matrix
ME_cor <- melt(ME_cor)
names(ME_cor) <- c("node1", "node2", "r")
pvals <- rcorr(me.matrix)$P
pvals <- melt(pvals)
names(pvals) <- c("node1", "node2", "pvalue")
ME_cor$pvalue <- pvals$pvalue #add pvalues to the correlations
ME_cor <- na.omit(ME_cor) #omit the 1:1 correlations and pvalues
cor <- subset(ME_cor,!duplicated(r))
cor$abs <- abs(cor$r)
cor$abs <- as.numeric(cor$abs)
labels <- data.frame(table(moduleColors))
labels <- labels[! (labels$moduleColors=="grey"),] #remove the grey module
labels$node1 <- as.factor(paste("ME", labels$moduleColors, sep=""))
cor1 <- merge(cor, labels, by="node1")
cor1$q <- p.adjust(cor1$pvalue, method="bonferroni")
cyto <- subset(cor1, q<=0.05)
cyto <- cyto[order(cyto$node1),]
dim(cyto)
links <- cyto
names(links) <- c("from", "to", "r", "pvalue", "abs", "moduleColors", "Freq", "q")
links <- links[order(links$abs, decreasing = TRUE),]
nodes <- labels[,c(3,1,2)]
names(nodes) <- c("id", "moduleColors", "Freq")

# igraph code
# net <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
# net
# plot(net)
# net <- simplify(net, remove.multiple = F, remove.loops = T)
# 
# plot(net, edge.arrow.size=.4,vertex.label=NA)
# V(net)$color <- as.character(nodes$moduleColors)
# # Compute node degrees (#links) and use that to set node size:
# deg <- degree(net, mode="all")
# V(net)$size <- deg*3
# # We could also use the size value:
# V(net)$size <- log(nodes$Freq)*5
# # The labels are currently node IDs.
# # Setting them to NA will render no labels:
# V(net)$label <- as.vector(nodes$moduleColors)
# # Set edge width based on weight:
# #E(net)$width <- ifelse(links$Pvalues > 0.05,0,abs(links$Correlation*5))
# E(net)$width <- abs(links$r*10)
# #change arrow size and edge color:
# E(net)$arrow.size <- 1
# E(net)$edge.color <- "black"
# # Set the network layout:
# par(mfrow=c(1,1))
# graph_attr(net, "layout") <- layout_with_lgl
# plot(net)
# ####Explore different Layouts
# layouts <- grep("^layout_", ls("package:igraph"), value=TRUE)[-1] 
# # Remove layouts that do not apply to our graph.
# layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama|tree", layouts)]
# V(net)$frame.color <- "black"
# par(mfrow=c(3,3), mar=c(1,1,1,1))
# for (layout in layouts) {
#   print(layout)
#   l <- do.call(layout, list(net))
#   plot(net, edge.arrow.mode=0, layout=l, main=layout) }

# Choose final layout
net <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
net <- simplify(net, remove.multiple = F, remove.loops = T)
# V(net)$size <- log(nodes$Freq)*5
V(net)$color <- as.character(nodes$moduleColors)
V(net)$label <- ''
E(net)$width <- abs(links$r*7)
V(net)$label <- as.vector(nodes$moduleColors)
V(net)$frame.color <- "black"
#change arrow size and edge color:
E(net)$arrow.size <- 0.5
E(net)$edge.color <- "black"
# save(net_groups, file="net_groups_gastroc_network.RData")
# load(file="net_groups_gastroc_network.RData")
# coords <- layout_in_circle(net, order=order(membership(net_groups)))
pdf("Figures/lung_igraph.pdf", h=10, w=10)
# plot(net, layout=coords)
plot(net, layout=layout_with_mds)
dev.off()

#barplot of modules sizes
pdf("Figures/lung_mod_size.pdf", h=10, w=10)
nodes1 <- nodes[order(nodes$Freq, decreasing = TRUE), ]
barplot(nodes1$Freq, col = as.character(nodes1$moduleColors), ylim = c(0,1400))
dev.off()







