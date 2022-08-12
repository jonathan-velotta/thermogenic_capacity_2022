library(WGCNA)
library(reshape)
library(lmerTest)



options(stringsAsFactors = FALSE);

# note: use Allow in RStudio, enable in R
allowWGCNAThreads()

######## Read-in files and prepare expression datset ###########
# Read-in files:
thermo_counts <- read.table("thermo_capacity_counts.txt", header = TRUE, sep = "\t")
id <- read.csv("thermo_capacity_rnaseq_samples.csv") 
traitData0 <- read.csv("thermo_capacity_traits_all.csv")

# format expression data so col = mice, rows = genes
thermo <- thermo_counts[, -c(2:6)] 
rownames(thermo) <- thermo$Geneid 
thermo$Geneid <- NULL 
id <- id[order(id$mouse_id),] 
colnames(thermo) <- id$ID 

###### Prepare gastroc dataset ################

# subset gastroc specific data from expr dataset using sample id file ()
gastroc_id <- subset(id, tissue == "gas") 
gastroc_id <- gastroc_id[order(gastroc_id$mouse_id),]
gas_id <- gastroc_id$ID 
gastroc <- thermo[,colnames(thermo) %in% gas_id] 
match(gastroc_id$ID, colnames(gastroc)) 
colnames(gastroc) <- gastroc_id$mouse_id 
# only keep data with >=10 reads
gastroc$mean = rowMeans(gastroc) 
keep_gastroc = subset(gastroc, mean >= 10) 
keep_gastroc$mean = NULL 
# save gastroc data as csv file (pre and post removal of >= 10 reads data)
write.csv(gastroc, "gastroc_raw_counts_thermocapacity.csv")
write.csv(keep_gastroc, "gastroc_raw_counts_thermocapacity.csv")

######## Create dataframe of treatment conditions for each sample (=mouse) ##############

# vector of populations/acclimations - > data frame for treatment of each mouse
population <- gastroc_id$population 
acclimation <- gastroc_id$acclimation 
table <- data.frame(population, acclimation)
# Paste together the pop + acclimation data for each mouse,
# Separated by a "_" and set as the factor of the dataframe:
group <- factor(paste(table$population, table$acclimation, sep = "_"))
# Add the group vector to the table dataframe:
cbind(table, group = group)
table$population = as.factor(table$population)
table$population <- relevel(table$population, ref = "LN")
# Normalize gastroc data and save the csv file
y0 <- DGEList(counts = keep_gastroc, group = group) 
y <- calcNormFactors(y0) 
gastroc_norm <- cpm(y0, log = TRUE, prior.count=1, normalized.lib.sizes = TRUE) 
write.csv(gastroc_norm, file = "thermo_gastroc_norm_counts.csv")

############### Create dataframe for expression data ###########################
# Transpose and rename gastroc data frame to reflect data format needed for
# WGCNA:
Expr0 = as.data.frame(t(gastroc_norm)) 
rownames(Expr0) = colnames(gastroc_norm)
# Remove genes with too many missing values
gsg = goodSamplesGenes(Expr0, verbose = 3)  
gsg$allOK
Expr = Expr0
Samples <- rownames(Expr)
# Remove samples for which we do not have trait data
traitData <- traitData0[traitData0$mouse_id %in% Samples,] 
remove <- subset(Samples, !(Samples %in% traitData$mouse_id))
# Omit all traits except mass and vo2
traitx <- traitData0[, c(1:10)] 
traitx <- na.omit(traitx)
# Create a new Expr dataframe with missing samples removed
Expr1 <- Expr[!row.names(Expr)%in%remove,]
traitRows <- match(rownames(Expr1), traitData$mouse_id) 
Traits0 <- traitData[traitRows, -1]
rownames(Traits0) <- traitData[traitRows, 1]
Traits <- Traits0
traity <- traitx[traitx$mouse_id %in% rownames(Traits),]
match(rownames(Traits), rownames(Expr1))
collectGarbage()

############ Sample outliers ###################################################
# Detect sample outliers
sampleTree = hclust(dist(datExpr), method = "average");
sizeGrWindow(12,9)
pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

######### Softhresholding power ################################################
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
sft = pickSoftThreshold(Expr, powerVector = powers, verbose = 5, networkType = "signed")
#pdf('wgcna/rv_beta_plot.pdf', h=4, w=7)
#dev.off()

############### Gastroc blockwise modules network ##############################
Net <- blockwiseModules(Expr, power = 8, maxBlockSize = dim(Expr)[2],TOMType = "signed", networkType = "signed", minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25,numericLabels = TRUE, pamRespectsDendro = FALSE,saveTOMs = TRUE,saveTOMFileBase = "ExprTOM",verbose = 3)
# Turn gene module numbers into colors
moduleLabels = Net$colors
moduleColors = labels2colors(Net$colors)
MEs = Net$MEs;
geneTree = Net$dendrograms[[1]];
# Should be 15 modules)
dim(table(moduleColors)) 
# Create a csv of module colors and number of associated genes:
write.csv(as.data.frame(table(moduleColors)),file="gastroc_module_table.csv")


# Calculate number of genes (= columnes) and samples (= row/mouse)
nGenes <- ncol(Expr1);
nSamples <- nrow(Expr1)


########## Associating mouse traits with genes #################################


# Subset focal trais associate with gastroc data:
gastroc_traits <- Traits[, c(8:11)]
# Correlation of each gene with each ME:
MEs0 = moduleEigengenes(Expr1, moduleColors)$eigengenes
MEs = orderMEs(MEs0) 
# Make sure the rownames from gastroc_trais and Expr1 line up correctly
# Should be ascending numbers
match(rownames(gastroc_traits), rownames(Expr1))

###### Correlate trais with MEs ################################################
moduleTraitCor <- round(cor(MEs, gastroc_traits, use = "p"), 3)
moduleTraitPvalue <- round(corPvalueStudent(moduleTraitCor, nSamples), 3)
moduleTraitQvalue <- round(matrix(p.adjust(as.vector(as.matrix(moduleTraitPvalue)), method='fdr'),ncol=4), 3)
colnames(moduleTraitQvalue) <- colnames(gastroc_traits)
rownames(moduleTraitQvalue) <- colnames(MEs)
# Create table of modules with significant correlations (after FDR correction):
moduleTrait <- moduleTraitQvalue[apply(moduleTraitQvalue, MARGIN = 1, function(x) any(x <= 0.05)), ]
# Write csv of gastroc model trait correlation:
write.csv(moduleTrait, file="Tables/gastroc_module_trait_cor_FDR.csv")


####### Create correlation dataframe for linar mixed effects models  ###########
# Make sure the row names are equal in cor.table and Traits df
match(rownames(cor.table0), rownames(Traits))
cor.table0$population <- Traits$population
cor.table0$family <- Traits$family
cor.table0$mass <- Traits$mass
cor.table0$vo2 <- Traits$vo2
cor.table0$vo2_mass <- Traits$vo2_mass
cor.table0$arterial_sat <- Traits$arterial_sat
cor.table <- melt.data.frame(cor.table0, id.vars=c("population", "family", "mass", "vo2", "vo2_mass", "arterial_sat"))
colnames(cor.table) <- c("population", "family", "mass", "vo2", "vo2_mass", "arterial_sat", "module", "ME")


########### Linear mixed effects models ########################################
cor.lmer <- dlply(cor.table, .(module), lmer, formula = vo2~ME + mass+ (1|family:population))
cor.out <- lapply(cor.lmer, summary)
cor.out.table <- data.frame(matrix(unlist(cor.out), nrow=length(cor.out), byrow=T))[, c(42:43)]
rownames(cor.out.table) <- names(cor.out)
colnames(cor.out.table) <- c("p_vo2", "p_mass")
sig.cor <- subset(cor.out.table, p_vo2<0.05) #there are no LMER effects of module on vo2


############# Create the starting data frame #################################

# Define variable weight containing the weight column of datTrait
vo2_mass <- as.data.frame(gastroc_traits$vo2_mass)
names(vo2_mass) <- "vo2_mass"
modNames <- substring(names(MEs), 3)
geneModuleMembership <- as.data.frame(cor(Expr, MEs, use = "p"));
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) <- paste("MM.", modNames, sep="");
names(MMPvalue) <- paste("p.MM", modNames, sep="");
geneTraitSignificance <- as.data.frame(cor(Expr1, vo2_mass, use = "p"));
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) <- paste("GS.", names(vo2_mass), sep="");
names(GSPvalue) <- paste("p.GS.", names(vo2_mass), sep="");
genes <- names(Expr)
GSPvalue$qvalue <- p.adjust(GSPvalue$p.GS.vo2_mass, method='fdr')

geneInfoGastroc <- data.frame(Gene = genes, geneTraitSignificance, GSPvalue, moduleColor = moduleColors, geneModuleMembership, MMPvalue)
write.csv(geneInfoGastroc, file="geneInfoGastroc.csv")

