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

options(stringsAsFactors = FALSE);
# note: use "allow" in RStudio, "enable" in R
allowWGCNAThreads()

#setwd("/Home/Documents/thermogenic_capacity_2022")
############ Read-in files and prepare expression datset #######################

# Read-in files:
thermo_counts <- read.table("thermo_capacity_counts.txt", header = TRUE, sep = "\t")

#thermo <- thermo0[, -c(2:6)] #omit chromosome, start/end, and length data
thermo <- thermo_counts[, -c(2:6)] 
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

##################
#### gastroc analysis #########


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
#clean up dataset
keep_gastroc$mean = NULL 
dim(keep_gastroc)
# make a DGE list
y0 <- DGEList(counts=keep_gastroc, group=group) 
# normalize
y <- calcNormFactors(y0) 
plotMDS(y) #exploration plots
y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
#cpm normalized and log transformed expression dataset
gastroc.norm <- cpm(y0, log=TRUE, prior.count=1, normalized.lib.sizes=TRUE) 
#write.csv(gastroc.norm, file="thermo_gastroc_norm_counts.csv")
plotMDS(gastroc.norm)
head(keep_gastroc)



####################### GASTROC WGCNA ###########################

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
#all trait data
traitData0 <- read.csv("thermo_capacity_traits_all.csv")
dev.off()
#boxplot(vo2_mass~population*acclimation, data=traitData0)
Samples <- rownames(Expr)
#We do not have trait data for all samples
traitData <- traitData0[traitData0$mouse_id %in% Samples,] 
# list of samples for which there is no trait data associated
remove <- subset(Samples, !(Samples %in% traitData$mouse_id))
# re-plot vo2 data ###
# log=xy makes this a log-log plot
plot(traitData0$mass, traitData0$vo2, log="xy") 
ln <- subset(traitData0, population=="LN")
me <- subset(traitData0, population=="ME")
plot(ln$mass, ln$vo2, xlim=c(12,30), ylim=c(1,7), pch=21, bg="darkorange", cex=2, log="xy")
points(me$mass, me$vo2, pch=21, bg="darkblue", cex=2)

#
# fit a nonlinear regression for the mass metabolic rate relationship
mod <- nls(vo2~a*(mass^b), data=traitData0, start=list(a=1, b=1)) 
lines(sort(traitData0$mass), predict(mod, list(mass=sort(traitData0$mass))), lwd=3, col="black")

##### make trait dataframe for WGCNA ########################################
# omit all traits except mass and vo2
traitx <- traitData0[, c(1:10)] 
traitx <- na.omit(traitx)

########Create a new Expr dataframe with missing samples removed #############
Expr1 = Expr[!row.names(Expr)%in%remove,]
# make sure IDs are in the same order as Expr dataset.
traitRows <- match(rownames(Expr1), traitData$mouse_id) 
Traits0 <- traitData[traitRows, -1]
rownames(Traits0) <- traitData[traitRows, 1]
Traits <- Traits0
dim(Traits)
traity <- traitx[traitx$mouse_id %in% rownames(Traits),]
dim(traity)
# should be in order if datasets are matched
match(rownames(Traits), rownames(Expr1))
collectGarbage()

#
#fit a nonlinear regression for the mass metabolic rate relationship
mod <- nls(vo2~a*(mass^b), data=traity, start=list(a=1, b=1)) 
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
# sft = 7
Net <- blockwiseModules(Expr, power = 7, maxBlockSize = dim(Expr)[2],TOMType = "signed", networkType = "signed", minModuleSize = 30,reassignThreshold = 0, mergeCutHeight = 0.25,numericLabels = TRUE, pamRespectsDendro = FALSE,saveTOMs = TRUE,saveTOMFileBase = "ExprTOM",verbose = 3)
save(Net,file = "thermo_capacity_gastroc_Network.RData")

#load gastroc network
load(file = "thermo_capacity_gastroc_Network.RData")#load saved Network file
table(Net$colors)
moduleLabels = Net$colors
moduleColors = labels2colors(Net$colors)
MEs = Net$MEs;
dim(MEs)
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


############## Trait Module Associations ###################
gastroc.traits <- Traits[, c(8:11)] # subset data to only vo2 and o2 sat
names(gastroc.traits)
# Define numbers of genes and samples
nGenes = ncol(Expr1);
nSamples = nrow(Expr1)
nGenes
nSamples
MEs0 = moduleEigengenes(Expr1, moduleColors)$eigengenes
MEs = orderMEs(MEs0) #the rownames of this dataset are equal to Expr
dim(MEs)
match(rownames(gastroc.traits), rownames(Expr1))
moduleTraitCor = round(cor(MEs, gastroc.traits, use = "p"), 3)
moduleTraitPvalue = round(corPvalueStudent(moduleTraitCor, nSamples), 3)
moduleTraitQvalue <- round(matrix(p.adjust(as.vector(as.matrix(moduleTraitPvalue)), method='fdr'),ncol=4), 3)
colnames(moduleTraitQvalue) <- colnames(gastroc.traits)
rownames(moduleTraitQvalue) <- colnames(MEs)
moduleTrait <- moduleTraitQvalue[apply(moduleTraitQvalue, MARGIN = 1, function(x) any(x <= 0.05)), ]#table of modules with significant correlations after FDR correction
write.csv(moduleTrait, file="Tables/gastroc_module_trait_cor_FDR.csv")



### Module trait linar mixed effects models  ###############################

MEs0 = moduleEigengenes(Expr1, moduleColors)$eigengenes
#the rownames of this dataset are equal to Expr
MEs = orderMEs(MEs0) 
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
cor.table


### Linear mixed effects model ###########################################
cor.lmer <- dlply(cor.table, .(module), lmer, formula = vo2~ME + mass+ (1|family:population))
cor.out <- lapply(cor.lmer, summary)
cor.out.table <- data.frame(matrix(unlist(cor.out), nrow=length(cor.out), byrow=T))[, c(42:43)]
rownames(cor.out.table) <- names(cor.out)
colnames(cor.out.table) <- c("p_vo2", "p_mass")
sig.cor <- subset(cor.out.table, p_vo2<0.05) 
cor.out.table
#there are no LMER effects of module on vo2


########## Create the starting data frame ####################################


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


#################################### LMER on ME #############################
head(gastroc.id)
match(rownames(Expr), gastroc.id$mouse_id)
setSamples = rownames(Expr)
#redo for all samples
MEs0<-moduleEigengenes(Expr, moduleColors)$eigengenes 
#the rownames of this dataset are equal to Expr
MEs<-orderMEs(MEs0) 
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
model <- dlply(gas.table, .(module), lmer, formula = ME~population*po2*temp + (1|family:population)) #run 



################# Linear Mixed Effects Model ##################################

# output anova table from model results
out <- lapply(model, anova) 
out

# create dataframe from output table
out.table <- round(data.frame(matrix(unlist(out), nrow=length(out), byrow=T))[, c(36:42)],3) 
rownames(out.table) <- names(out)
colnames(out.table) <- c("population", "po2", "temp", "pop:po2", "pop:temp", "po2:temp", "3-way")

qTable <- round(matrix(p.adjust(as.vector(as.matrix(out.table)), method='fdr'),ncol=7), 3)
rownames(qTable) <- rownames(out.table)
colnames(qTable) <- colnames(out.table)

# table of modules with significance after FDR correction
qTable.sig <- qTable[apply(qTable[, ], MARGIN = 1, function(x) any(x <= 0.05)), ]
write.csv(qTable.sig, file="Tables/gastroc_lmer_sig_pvalues.csv")
# table of F values
# f.out <- round(data.frame(matrix(unlist(out), nrow=length(out), byrow=T))[, c(29:35)],3) 

# at least one pop effect
qPop <- qTable[apply(qTable[, -c(2:3, 6)], MARGIN = 1, function(x) any(x <= 0.05)), ]; 
# at least a treat effect
qTreat <-  qTable[apply(qTable[, c(2:3, 6)], MARGIN = 1, function(x) any(x <= 0.05)), ]; 

dim(qPop) 
dim(qTreat) 


inter <- intersect(rownames(qPop), rownames(qTreat))
qTreat.only <- qTreat[!(rownames(qTreat) %in% inter), ]; dim(qTreat.only)


# a list of the colors used for module ID
modNames

############ TO DO ################
# rewrite code using for loop
#for (item in modNames)
#  {
#  item <- subset(gas.table, module = "ME"+item)
#}

#################### MEred ###############################################
red <- subset(gas.table, module=="MEred")
# three-way lmer
anova(lmer(ME~population*po2*temp + (1|family:population), data = red)) 
# results: p = 0.005 pop*po2, p = 0.047 pop*temp

# posthoc analysis one-way for RED

# posthoc all - RED
N <- subset(red, acclimation == "N")
H <- subset(red, acclimation == "H")
C <- subset(red, acclimation == "C")
CH <- subset(red, acclimation == "CH")

# posthoc Mt Evens - RED
me <- subset(red, population == "ME")
me.warm <- subset(me, temp=="warm")
me.C <- subset(me, acclimation == "C")
me.CH <- subset(me, acclimation == "CH")
me.N <- subset(me, acclimation== "N")

# posthoc Lincoln - RED
ln <- subset(red, population == "LN")
ln.warm <- subset(ln, temp=="warm")
ln.C <- subset(ln, acclimation == "C")
ln.CH <- subset(ln, acclimation == "CH")
ln.N <- subset(ln, acclimation== "N")

# Mt Events results - RED
me.NC <- rbind(me.N, me.C); 
me.NCH <- rbind(me.N, me.CH); 

# Lincoln results - RED
NC <- rbind(ln.N, ln.C); 
NCH <- rbind(ln.N, ln.CH); 

summary(lmer(ME~acclimation + (1|family), data=me))
# acclimationH p = 0.009, acclimationN = 0.04

anova(lmer(ME~population + (1|family:population), data=N))
anova(lmer(ME~population + (1|family:population), data=H))
anova(lmer(ME~population + (1|family:population), data=C))
anova(lmer(ME~population + (1|family:population), data=CH))
anova(lmer(ME~po2 + (1|family:population), data = me.warm))
anova(lmer(ME~acclimation + (1|family:population), data = me.NC))
anova(lmer(ME~acclimation + (1|family:population), data = me.NCH))
anova(lmer(ME~po2 + (1|family:population), data = ln.warm))
anova(lmer(ME~acclimation + (1|family:population), data = NC))
anova(lmer(ME~acclimation + (1|family:population), data = NCH))





modNames
#################### MEpink ###############################################
pink <- subset(gas.table, module=="MEpink")
# three-way lmer
anova(lmer(ME~population*po2*temp + (1|family:population), data = pink)) 
# results: p = 0.005 pop*po2, p = 0.047 pop*temp

# posthoc analysis one-way for PINK

# posthoc all - PINK
N <- subset(pink, acclimation == "N")
H <- subset(pink, acclimation == "H")
C <- subset(pink, acclimation == "C")
CH <- subset(pink, acclimation == "CH")

# posthoc Mt Evens - PINK
me <- subset(pink, population == "ME")
me.warm <- subset(me, temp=="warm")
me.C <- subset(me, acclimation == "C")
me.CH <- subset(me, acclimation == "CH")
me.N <- subset(me, acclimation== "N")

# posthoc Lincoln - PINK
ln <- subset(pink, population == "LN")
ln.warm <- subset(ln, temp=="warm")
ln.C <- subset(ln, acclimation == "C")
ln.CH <- subset(ln, acclimation == "CH")
ln.N <- subset(ln, acclimation== "N")

# Mt Events results - PINK
me.NC <- rbind(me.N, me.C); 
me.NCH <- rbind(me.N, me.CH); 

# Lincoln results - PINK
NC <- rbind(ln.N, ln.C); 
NCH <- rbind(ln.N, ln.CH); 

summary(lmer(ME~acclimation + (1|family), data=me))
# acclimationH p = 0.009, acclimationN = 0.04

anova(lmer(ME~population + (1|family:population), data=N))
anova(lmer(ME~population + (1|family:population), data=H))
anova(lmer(ME~population + (1|family:population), data=C))
anova(lmer(ME~population + (1|family:population), data=CH))
anova(lmer(ME~po2 + (1|family:population), data = me.warm))
anova(lmer(ME~acclimation + (1|family:population), data = me.NC))
anova(lmer(ME~acclimation + (1|family:population), data = me.NCH))
anova(lmer(ME~po2 + (1|family:population), data = ln.warm))
anova(lmer(ME~acclimation + (1|family:population), data = NC))
anova(lmer(ME~acclimation + (1|family:population), data = NCH))


modNames
#################### MEsalmon ###############################################
salmon <- subset(gas.table, module=="MEsalmon")
# three-way lmer
anova(lmer(ME~population*po2*temp + (1|family:population), data = salmon)) 
# results: p = 0.005 pop*po2, p = 0.047 pop*temp

# posthoc analysis one-way for SALMON

# posthoc all - SALMON
N <- subset(salmon, acclimation == "N")
H <- subset(salmon, acclimation == "H")
C <- subset(salmon, acclimation == "C")
CH <- subset(salmon, acclimation == "CH")

# posthoc Mt Evens - SALMON
me <- subset(salmon, population == "ME")
me.warm <- subset(me, temp=="warm")
me.C <- subset(me, acclimation == "C")
me.CH <- subset(me, acclimation == "CH")
me.N <- subset(me, acclimation== "N")

# posthoc Lincoln - SALMON
ln <- subset(salmon, population == "LN")
ln.warm <- subset(ln, temp=="warm")
ln.C <- subset(ln, acclimation == "C")
ln.CH <- subset(ln, acclimation == "CH")
ln.N <- subset(ln, acclimation== "N")

# Mt Events results - SALMON
me.NC <- rbind(me.N, me.C); 
me.NCH <- rbind(me.N, me.CH); 

# Lincoln results - SALMON
NC <- rbind(ln.N, ln.C); 
NCH <- rbind(ln.N, ln.CH); 

summary(lmer(ME~acclimation + (1|family), data=me))
# acclimationH p = 0.009, acclimationN = 0.04

### 
anova(lmer(ME~population + (1|family:population), data=N))
anova(lmer(ME~population + (1|family:population), data=H))
anova(lmer(ME~population + (1|family:population), data=C))
anova(lmer(ME~population + (1|family:population), data=CH))
anova(lmer(ME~po2 + (1|family:population), data = me.warm))
anova(lmer(ME~acclimation + (1|family:population), data = me.NC))
anova(lmer(ME~acclimation + (1|family:population), data = me.NCH))
anova(lmer(ME~po2 + (1|family:population), data = ln.warm))
anova(lmer(ME~acclimation + (1|family:population), data = NC))
anova(lmer(ME~acclimation + (1|family:population), data = NCH))



modNames
#################### MEbrown ###############################################
brown <- subset(gas.table, module=="MEbrown")
# three-way lmer
anova(lmer(ME~population*po2*temp + (1|family:population), data = brown)) 
# results: p = 0.005 pop*po2, p = 0.047 pop*temp

# posthoc analysis one-way for BROWN

# posthoc all - BROWN
N <- subset(brown, acclimation == "N")
H <- subset(brown, acclimation == "H")
C <- subset(brown, acclimation == "C")
CH <- subset(brown, acclimation == "CH")

# posthoc Mt Evens - BROWN
me <- subset(brown, population == "ME")
me.warm <- subset(me, temp=="warm")
me.C <- subset(me, acclimation == "C")
me.CH <- subset(me, acclimation == "CH")
me.N <- subset(me, acclimation== "N")

# posthoc Lincoln - BROWN
ln <- subset(brown, population == "LN")
ln.warm <- subset(ln, temp=="warm")
ln.C <- subset(ln, acclimation == "C")
ln.CH <- subset(ln, acclimation == "CH")
ln.N <- subset(ln, acclimation== "N")

# Mt Events results - BROWN
me.NC <- rbind(me.N, me.C); 
me.NCH <- rbind(me.N, me.CH); 

# Lincoln results - BROWN
NC <- rbind(ln.N, ln.C); 
NCH <- rbind(ln.N, ln.CH); 

summary(lmer(ME~acclimation + (1|family), data=me))
# acclimationH p = 0.009, acclimationN = 0.04

### 
anova(lmer(ME~population + (1|family:population), data=N))
anova(lmer(ME~population + (1|family:population), data=H))
anova(lmer(ME~population + (1|family:population), data=C))
anova(lmer(ME~population + (1|family:population), data=CH))
anova(lmer(ME~po2 + (1|family:population), data = me.warm))
anova(lmer(ME~acclimation + (1|family:population), data = me.NC))
anova(lmer(ME~acclimation + (1|family:population), data = me.NCH))
anova(lmer(ME~po2 + (1|family:population), data = ln.warm))
anova(lmer(ME~acclimation + (1|family:population), data = NC))
anova(lmer(ME~acclimation + (1|family:population), data = NCH))



modNames
#################### MEmagenta ###############################################
magenta <- subset(gas.table, module=="MEmagenta")
# three-way lmer
anova(lmer(ME~population*po2*temp + (1|family:population), data = magenta)) 
# results: p = 0.005 pop*po2, p = 0.047 pop*temp

# posthoc analysis one-way for MAGENTA

# posthoc all - MAGENTA
N <- subset(magenta, acclimation == "N")
H <- subset(magenta, acclimation == "H")
C <- subset(magenta, acclimation == "C")
CH <- subset(magenta, acclimation == "CH")

# posthoc Mt Evens - MAGENTA
me <- subset(magenta, population == "ME")
me.warm <- subset(me, temp=="warm")
me.C <- subset(me, acclimation == "C")
me.CH <- subset(me, acclimation == "CH")
me.N <- subset(me, acclimation== "N")

# posthoc Lincoln - MAGENTA
ln <- subset(magenta, population == "LN")
ln.warm <- subset(ln, temp=="warm")
ln.C <- subset(ln, acclimation == "C")
ln.CH <- subset(ln, acclimation == "CH")
ln.N <- subset(ln, acclimation== "N")

# Mt Events results - MAGENTA
me.NC <- rbind(me.N, me.C); 
me.NCH <- rbind(me.N, me.CH); 

# Lincoln results - MAGENTA
NC <- rbind(ln.N, ln.C); 
NCH <- rbind(ln.N, ln.CH); 

summary(lmer(ME~acclimation + (1|family), data=me))
# acclimationH p = 0.009, acclimationN = 0.04

### 
anova(lmer(ME~population + (1|family:population), data=N))
anova(lmer(ME~population + (1|family:population), data=H))
anova(lmer(ME~population + (1|family:population), data=C))
anova(lmer(ME~population + (1|family:population), data=CH))
anova(lmer(ME~po2 + (1|family:population), data = me.warm))
anova(lmer(ME~acclimation + (1|family:population), data = me.NC))
anova(lmer(ME~acclimation + (1|family:population), data = me.NCH))
anova(lmer(ME~po2 + (1|family:population), data = ln.warm))
anova(lmer(ME~acclimation + (1|family:population), data = NC))
anova(lmer(ME~acclimation + (1|family:population), data = NCH))


modNames
#################### MEgreen ###############################################
green <- subset(gas.table, module=="MEgreen")
# three-way lmer
anova(lmer(ME~population*po2*temp + (1|family:population), data = green)) 
# results: p = 0.005 pop*po2, p = 0.047 pop*temp

# posthoc analysis one-way for GREEN

# posthoc all - GREEN
N <- subset(green, acclimation == "N")
H <- subset(green, acclimation == "H")
C <- subset(green, acclimation == "C")
CH <- subset(green, acclimation == "CH")

# posthoc Mt Evens - GREEN
me <- subset(green, population == "ME")
me.warm <- subset(me, temp=="warm")
me.C <- subset(me, acclimation == "C")
me.CH <- subset(me, acclimation == "CH")
me.N <- subset(me, acclimation== "N")

# posthoc Lincoln - GREEN
ln <- subset(green, population == "LN")
ln.warm <- subset(ln, temp=="warm")
ln.C <- subset(ln, acclimation == "C")
ln.CH <- subset(ln, acclimation == "CH")
ln.N <- subset(ln, acclimation== "N")

# Mt Events results - GREEN
me.NC <- rbind(me.N, me.C); 
me.NCH <- rbind(me.N, me.CH); 

# Lincoln results - GREEN
NC <- rbind(ln.N, ln.C); 
NCH <- rbind(ln.N, ln.CH); 

summary(lmer(ME~acclimation + (1|family), data=me))
# acclimationH p = 0.009, acclimationN = 0.04

### 
anova(lmer(ME~population + (1|family:population), data=N))
anova(lmer(ME~population + (1|family:population), data=H))
anova(lmer(ME~population + (1|family:population), data=C))
anova(lmer(ME~population + (1|family:population), data=CH))
anova(lmer(ME~po2 + (1|family:population), data = me.warm))
anova(lmer(ME~acclimation + (1|family:population), data = me.NC))
anova(lmer(ME~acclimation + (1|family:population), data = me.NCH))
anova(lmer(ME~po2 + (1|family:population), data = ln.warm))
anova(lmer(ME~acclimation + (1|family:population), data = NC))
anova(lmer(ME~acclimation + (1|family:population), data = NCH))


modNames
#################### MEblack ###############################################
black <- subset(gas.table, module=="MEblack")
# three-way lmer
anova(lmer(ME~population*po2*temp + (1|family:population), data = black)) 
# results: p = 0.005 pop*po2, p = 0.047 pop*temp

# posthoc analysis one-way for BLACK

# posthoc all - BLACK
N <- subset(black, acclimation == "N")
H <- subset(black, acclimation == "H")
C <- subset(black, acclimation == "C")
CH <- subset(black, acclimation == "CH")

# posthoc Mt Evens - BLACK
me <- subset(black, population == "ME")
me.warm <- subset(me, temp=="warm")
me.C <- subset(me, acclimation == "C")
me.CH <- subset(me, acclimation == "CH")
me.N <- subset(me, acclimation== "N")

# posthoc Lincoln - BLACK
ln <- subset(black, population == "LN")
ln.warm <- subset(ln, temp=="warm")
ln.C <- subset(ln, acclimation == "C")
ln.CH <- subset(ln, acclimation == "CH")
ln.N <- subset(ln, acclimation== "N")

# Mt Events results - BLACK
me.NC <- rbind(me.N, me.C); 
me.NCH <- rbind(me.N, me.CH); 

# Lincoln results - BLACK
NC <- rbind(ln.N, ln.C); 
NCH <- rbind(ln.N, ln.CH); 

summary(lmer(ME~acclimation + (1|family), data=me))
# acclimationH p = 0.009, acclimationN = 0.04

### 
anova(lmer(ME~population + (1|family:population), data=N))
anova(lmer(ME~population + (1|family:population), data=H))
anova(lmer(ME~population + (1|family:population), data=C))
anova(lmer(ME~population + (1|family:population), data=CH))
anova(lmer(ME~po2 + (1|family:population), data = me.warm))
anova(lmer(ME~acclimation + (1|family:population), data = me.NC))
anova(lmer(ME~acclimation + (1|family:population), data = me.NCH))
anova(lmer(ME~po2 + (1|family:population), data = ln.warm))
anova(lmer(ME~acclimation + (1|family:population), data = NC))
anova(lmer(ME~acclimation + (1|family:population), data = NCH))





###########################################################################
###########################################################################

###################### MEbrown ##############################################
brown <- subset(gas.table, module=="MEbrown")
anova(lmer(ME~population*po2*temp + (1|family:population), data=brown)) #three-way lmer
N <- subset(brown, acclimation=="N"); anova(lmer(ME~population + (1|family:population), data=N)) #posthoc analysis one-way
H <- subset(brown, acclimation=="H"); anova(lmer(ME~population + (1|family:population), data=H))
C <- subset(brown, acclimation=="C"); anova(lmer(ME~population + (1|family:population), data=C))
CH <- subset(brown, acclimation=="CH"); anova(lmer(ME~population + (1|family:population), data=CH))

#################### MEdarkgreen ###########################################
darkgreen <- subset(gas.table, module=="MEdarkgreen")
anova(lmer(ME~population*po2*temp + (1|family:population), data=darkgreen)) #three-way lmer
N <- subset(darkgreen, acclimation=="N"); anova(lmer(ME~population + (1|family:population), data=N)) #posthoc analysis one-way
H <- subset(darkgreen, acclimation=="H"); anova(lmer(ME~population + (1|family:population), data=H))
C <- subset(darkgreen, acclimation=="C"); anova(lmer(ME~population + (1|family:population), data=C))
CH <- subset(darkgreen, acclimation=="CH"); anova(lmer(ME~population + (1|family:population), data=CH))


##################### MElightgreen ########################################
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

################ MEpurple ################################################################
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

