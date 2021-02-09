## Analysis for NKs's Drosphilia Cardiomyocytes
# contact: alwhiteh@ucsd.edu
# First we have taken the bam files that have been filtered for mitochondial reads, PCR dupicliates, and compensated for multimapped reads
# then MACS2 was used to call peaks for each bam file
# The resulting .xls files were used in this analysis

## Load Libaries
library(ChIPseeker)
#library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#BiocManager::install("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
library(clusterProfiler)
library(org.Dm.eg.db)

## We are going to take the extended summits from MACS2 + bedtools slop by 1 kb
library(DiffBind)
#samples <- read.csv("/Users/alexanderwhitehead/Documents/PhD/RNASeq_stuff/CF_ATAC/CF_ATAC_Sample_Sheet_Ext_Sum250.csv", stringsAsFactors = FALSE)
setwd("/Users/alexanderwhitehead/Library/Mobile Documents/com~apple~CloudDocs/Documents/PhD/RNASeq_stuff/kirkland_aging_ATAC/1_v_5wk")
samples_age<- read.csv("/Users/alexanderwhitehead/Library/Mobile Documents/com~apple~CloudDocs/Documents/PhD/RNASeq_stuff/kirkland_aging_ATAC/Kirkland_Master_Sample_Sheet_Age.csv", stringsAsFactors = FALSE)
samples_lamb<- read.csv("/Users/alexanderwhitehead/Library/Mobile Documents/com~apple~CloudDocs/Documents/PhD/RNASeq_stuff/kirkland_aging_ATAC/Kirkland_Master_Sample_Sheet_LamB.csv", stringsAsFactors = FALSE)
samples_lamc<- read.csv("/Users/alexanderwhitehead/Library/Mobile Documents/com~apple~CloudDocs/Documents/PhD/RNASeq_stuff/kirkland_aging_ATAC/Kirkland_Master_Sample_Sheet_LamC.csv", stringsAsFactors = FALSE)
samples_all <- read.csv("/Users/alexanderwhitehead/Library/Mobile Documents/com~apple~CloudDocs/Documents/PhD/RNASeq_stuff/kirkland_aging_ATAC/Kirkland_Master_Sample_Sheet.csv", stringsAsFactors = FALSE)

All_ATAC <- dba(sampleSheet=samples_all)
Age_ATAC <- dba(sampleSheet=samples_age)
LamB_ATAC <-dba(sampleSheet=samples_lamb)
LamC_ATAC <-dba(sampleSheet=samples_lamc)

#Generate a correlation heatmap using peak caller scores
pdf("Peak_Caller_scores_Correlations.pdf")
plot(All_ATAC)
dev.off()

#Generate counts for each sample
Age_ATAC <- dba.count(Age_ATAC, bParallel = T, filter = 1, summits = 100 ) #this takes a while
#Generate counts for each sample
LamB_ATAC <- dba.count(LamB_ATAC, bParallel = T, filter = 1, summits = 100 ) #this takes a while
#Generate counts for each sample
LamC_ATAC <- dba.count(LamC_ATAC, bParallel = T, filter = 1, summits = 100 ) #this takes a while
#Generate counts for each sample
All_ATAC <- dba.count(All_ATAC, bParallel = T, filter = 1, summits = 100 ) #this takes a while
# filter = 1 to remove consensus peaks where no peak has an RPKM value of at least 1 in any sample
# summits = 100 is recommended in vignette for ATAC-Seq, centers peaks around summits using 100 bp (rather than default 401 bp)
All_ATAC # should show Fractions of Reads in Peaks, this is from .21-.37 for the all samples analysis

# This time will plot correlation heatmap based on affinity scores (aka Read Density aka reads/length)
pdf("Read_densities_Correlations.pdf")
plot(All_ATAC)
dev.off()

#estblish a contrast 
#tell diffbind what criteria to contrast
All_ATAC <- dba.contrast(All_ATAC, group1 = c(17:19), group2 = c(20:22), minMembers = 2) #contrast 1 = Age of w1118 flies, grp1 = 1wk, grp2 = 5wk
All_ATAC <- dba.contrast(All_ATAC, group1 = c(1:3), group2 = c(12:16), minMembers = 2) #contrast 2 = LamB KO, grp1 = wt, grp2 = LambKO, we are excluding 
All_ATAC <- dba.contrast(All_ATAC, group1 = c(4:7), group2 = c(8:11), minMembers = 2) #contrast 3 = LamC KO grp1 = wt, grp2 = LamCKO
All_ATAC # View all of your contrasts

# Here's how to clear a contrast if you mess up: Age_ATAC$constrasts= NULL
#All_ATAC$contrasts = NULL

# here is test for only Age
Age_ATAC <- dba.contrast(Age_ATAC, group1 = c(1:3), group2 = c(4:6), minMembers = 2)
# here is test for only LamB
LamB_ATAC <- dba.contrast(LamB_ATAC, group1 = c(1:3), group2 = c(4:8), minMembers = 2)
# here is test for only LamC
LamC_ATAC <- dba.contrast(LamC_ATAC, group1 = c(1:4), group2 = c(5:8), minMembers = 2)

# Perfom the differential anlaysis
#All_ATAC_EdgeR <- dba.analyze(Age_ATAC, method = DBA_EDGER_GLM) 
#All_ATAC_DESEQ2 <- dba.analyze(Age_ATAC, method = DBA_DESEQ2)
All_ATAC <- dba.analyze(All_ATAC, method = DBA_ALL_METHODS)
Age_ATAC <- dba.analyze(Age_ATAC, method = DBA_ALL_METHODS) # Remove the method option when not trying to beat the no replicates rule
LamB_ATAC <- dba.analyze(LamB_ATAC, method = DBA_ALL_METHODS) # Remove the method option when not trying to beat the no replicates rule
LamC_ATAC <- dba.analyze(LamC_ATAC, method = DBA_ALL_METHODS) # Remove the method option when not trying to beat the no replicates rule

All_ATAC # the last column shows what # of binding sites are different based on FDR = 0.05

# Create dummy variable to for indexing the number of contrasts
num_cont = c(1:length(All_ATAC$contrasts)) #number of contrasts
name_cont = c("Age","LamBKO","LamCKO") # name of contrasts

#Generate MA, Volcano and Correlation Heatmaps for Each Contrast
pdf("MA_and_Volcano_p_0.05.pdf")
dba.plotMA(Age_ATAC, th = 0.05, contrast = 1, fold = 0.32, method = DBA_EDGER, bFlip = T) # generates MA plot, th = FDR cutoff, fold is LOG2, so 2^0.32 = 1.25 FC
dba.plotMA(Age_ATAC, th = 0.05, contrast =1, fold = 0.32,method = DBA_DESEQ2, bFlip = T) # generates MA plot, th = FDR cutoff
dba.plotMA(LamB_ATAC, th = 0.05, contrast = 1, fold = 0.32, method = DBA_EDGER, bFlip = T) # generates MA plot, th = FDR cutoff, fold is LOG2, so 2^0.32 = 1.25 FC
dba.plotMA(LamB_ATAC, th = 0.05, contrast =1, fold = 0.32,method = DBA_DESEQ2, bFlip = T)
dba.plotMA(LamC_ATAC, th = 0.05, contrast = 1, fold = 0.32, method = DBA_EDGER, bFlip = T) # generates MA plot, th = FDR cutoff, fold is LOG2, so 2^0.32 = 1.25 FC
dba.plotMA(LamC_ATAC, th = 0.05, contrast =1, fold = 0.32,method = DBA_DESEQ2, bFlip = T)
dba.plotVolcano(Age_ATAC, th = 0.05, contrast = 1,fold = 0.32, method = DBA_EDGER, bFlip = T) # generates volcano plot
dba.plotVolcano(Age_ATAC, th = 0.05, contrast = 1, fold = 0.32, method = DBA_DESEQ2, bFlip = T)
dba.plotVolcano(LamB_ATAC, th = 0.05, contrast = 1,fold = 0.32, method = DBA_EDGER, bFlip = T) # generates volcano plot
dba.plotVolcano(LamB_ATAC, th = 0.05, contrast = 1, fold = 0.32, method = DBA_DESEQ2, bFlip = T)
dba.plotVolcano(LamC_ATAC, th = 0.05, contrast = 1,fold = 0.32, method = DBA_EDGER, bFlip = T) # generates volcano plot
dba.plotVolcano(LamC_ATAC, th = 0.05, contrast = 1, fold = 0.32, method = DBA_DESEQ2, bFlip = T)
plot(Age_ATAC, contrast = 1, method = DBA_EDGER, th = 0.05)
plot(Age_ATAC, contrast = 1, method = DBA_DESEQ2, th = 0.05)
plot(LamB_ATAC, contrast = 1, method = DBA_EDGER, th = 0.05)
plot(LamB_ATAC, contrast = 1, method = DBA_DESEQ2, th = 0.05)
plot(LamC_ATAC, contrast = 1, method = DBA_EDGER, th = 0.05)
plot(LamC_ATAC, contrast = 1, method = DBA_DESEQ2, th = 0.05)
dev.off()

# Make PCA plots and check for outliers
pdf("ATAC_PCA.pdf")
dba.plotPCA(Age_ATAC,label = DBA_ID)
dba.plotPCA(LamB_ATAC,label = DBA_ID)
dba.plotPCA(LamC_ATAC,label = DBA_ID)
dba.plotPCA(All_ATAC, label = DBA_ID)
dev.off()

# Make box plot of all groups
pdf("ATAC_Boxplot_All_Counts.pdf")
dba.plotBox(Age_ATAC, bAll=T)
dba.plotBox(LamB_ATAC, bAll = T)
dba.plotBox(LamC_ATAC, bAll = T)
dev.off()

pdf("ATAC_Boxplot_Sig_Counts.pdf")
dba.plotBox(Age_ATAC, th=0.1)
dba.plotBox(LamB_ATAC, th = 0.1)
dba.plotBox(LamC_ATAC, th = 0.1)
dev.off()

#Create Heapmap based on affinities for differentially bound sites
pdf("ATAC_Heatmaps_No_Anno.pdf")
dba.plotHeatmap(Age_ATAC, correlations = FALSE, contrast = 1, th = 0.1, colScheme = myColor, 
                   score= DBA_SCORE_NORMALIZED, method = DBA_EDGER)
dba.plotHeatmap(Age_ATAC, correlations = FALSE, contrast = 1, th = 0.1, colScheme = myColor,
                score= DBA_SCORE_NORMALIZED, method = DBA_DESEQ2)
dba.plotHeatmap(LamB_ATAC, correlations = FALSE, contrast = 1, th = 0.1, colScheme = myColor,
                 score= DBA_SCORE_NORMALIZED, method = DBA_EDGER)
dba.plotHeatmap(LamB_ATAC, correlations = FALSE, contrast = 1, th = 0.1, colScheme = myColor,
                score= DBA_SCORE_NORMALIZED, method = DBA_DESEQ2)
dba.plotHeatmap(LamC_ATAC, correlations = FALSE, contrast = 1, th = 0.1, colScheme = myColor,
                 score= DBA_SCORE_NORMALIZED, method = DBA_EDGER)
dba.plotHeatmap(LamC_ATAC, correlations = FALSE, contrast = 1, th = 0.1, colScheme = myColor,
                score= DBA_SCORE_NORMALIZED, method = DBA_DESEQ2)
dev.off()

Age_ATAC_DESEQ2 <- dba.report(Age_ATAC, th = 0.01,contrast = 1,fold = 0.32, method = DBA_DESEQ2)
Age_ATAC_EdgeR <- dba.report(Age_ATAC, th = 0.01,contrast = 1,fold = 0.32, method = DBA_EDGER)
LamB_ATAC_DESEQ2 <- dba.report(LamB_ATAC, th = 0.01,contrast = 1,fold = 0.32, method = DBA_DESEQ2)
LamB_ATAC_EdgeR <- dba.report(LamB_ATAC, th = 0.01,contrast = 1,fold = 0.32, method = DBA_EDGER)
LamC_ATAC_DESEQ2 <- dba.report(LamC_ATAC, th = 0.01,contrast = 1,fold = 0.32, method = DBA_DESEQ2)
LamC_ATAC_EdgeR <- dba.report(LamC_ATAC, th = 0.01,contrast = 1,fold = 0.32, method = DBA_EDGER)


#quantitatively see how many DEGS up and down
sum(Age_ATAC_DESEQ2$Fold<0) # Down - here we have 29, # ~200 lower limit of motif finding input
sum(Age_ATAC_DESEQ2$Fold>0) # Up - here we have 440

#see how many peaks overlap between # of datasets
pdf("Histogram_Overlaps.pdf")
plot(dba.overlap(All_ATAC, mode=DBA_OLAP_RATE), type='b', ylab=" # of peaks", xlab = "Overlap at least this many peaksets")
dev.off()

# We will create separate reports that contain all the peaks (th =1) 
data.age.deseq2 <-dba.report(Age_ATAC, th = 1, contrast = 1, method = DBA_DESEQ2)
data.age.edger <-dba.report(Age_ATAC, th = 1, contrast = 1,method = DBA_EDGER)
data.lamb.deseq2 <-dba.report(LamB_ATAC, th = 1, contrast = 1,method = DBA_DESEQ2)
data.lamb.edger <-dba.report(LamB_ATAC, th = 1, contrast = 1,method = DBA_EDGER)
data.lamc.deseq2 <-dba.report(LamC_ATAC, th = 1, contrast = 1,method = DBA_DESEQ2)
data.lamc.edger <-dba.report(LamC_ATAC, th = 1, contrast = 1,method = DBA_EDGER)
  
# Need to convert to data frame first
data.age.deseq2_df <- as.data.frame(data.age.deseq2)
data.age.edger_df <- as.data.frame(data.age.edger)
data.lamb.deseq2_df <- as.data.frame(data.lamb.deseq2)
data.lamb.edger_df <- as.data.frame(data.lamb.edger)
data.lamc.deseq2_df <- as.data.frame(data.lamc.deseq2)
data.lamc.edger_df <- as.data.frame(data.lamc.edger)

#write the rownames as a seperate column, "peak"
data.age.deseq2_df$peaks <- rownames(data.age.deseq2_df)
data.age.edger_df$peaks <- rownames(data.age.edger_df)
data.lamb.deseq2_df$peaks <- rownames(data.lamb.deseq2_df)
data.lamb.edger_df$peaks <- rownames(data.lamb.edger_df)
data.lamc.deseq2_df$peaks <- rownames(data.lamc.deseq2_df)
data.lamc.edger_df$peaks <- rownames(data.lamc.edger_df)

# Then add column names to be able to pick which ones we want to write to .bed file (so that HOMER can read it)
bedfile_colnames <- c("seqnames","start","end","width","strand","Conc","Conc_group1","Conc_group2",
                    "Fold","p.value","FDR","peak")

colnames(data.age.deseq2_df) <- bedfile_colnames
colnames(data.age.edger_df) <- bedfile_colnames
colnames(data.lamb.deseq2_df) <- bedfile_colnames
colnames(data.lamb.edger_df) <- bedfile_colnames
colnames(data.lamc.deseq2_df) <- bedfile_colnames
colnames(data.lamc.edger_df) <- bedfile_colnames

# select only the columns (in order), that HOMER will use
homer.age.deseq2 <- data.age.deseq2_df[c("seqnames","start","end","peak","FDR","strand","Fold")] #need to be in this order
homer.age.edger <- data.age.edger_df[c("seqnames","start","end","peak","FDR","strand","Fold")] #need to be in this order
homer.lamb.deseq2 <- data.lamb.deseq2_df[c("seqnames","start","end","peak","FDR","strand","Fold")] #need to be in this order
homer.lamb.edger <- data.lamb.edger_df[c("seqnames","start","end","peak","FDR","strand","Fold")] #need to be in this order
homer.lamc.deseq2 <- data.lamc.deseq2_df[c("seqnames","start","end","peak","FDR","strand","Fold")] #need to be in this order
homer.lamc.edger <- data.lamc.edger_df[c("seqnames","start","end","peak","FDR","strand","Fold")] #need to be in this order

#finally export to .bed files for use in annotatePeaks.pl
write.table(homer.age.deseq2, file = paste("Kirkland_Age_DESEQ2.bed", sep="_"), sep = "\t", row.names=F, col.names = F, qmethod = "double")
write.table(homer.age.edger, file = paste("Kirkland_Age_EdgeR.bed", sep="_"), sep = "\t", row.names=F, col.names = F, qmethod = "double")
write.table(homer.lamb.deseq2, file = paste("Kirkland_LamB_DESEQ2.bed", sep="_"), sep = "\t", row.names=F, col.names = F, qmethod = "double")
write.table(homer.lamb.edger, file = paste("Kirkland_LamB_EdgeR.bed", sep="_"), sep = "\t", row.names=F, col.names = F, qmethod = "double")
write.table(homer.lamc.deseq2, file = paste("Kirkland_LamC_DESEQ2.bed", sep="_"), sep = "\t", row.names=F, col.names = F, qmethod = "double")
write.table(homer.lamc.edger, file = paste("Kirkland_LamC_EdgeR.bed", sep="_"), sep = "\t", row.names=F, col.names = F, qmethod = "double")

# you need to manually go in and change column 5 to "1" and then format as "number" in excel. hit save
# Then run annotatePeaks.pls

#Now we need to import the annotated text files and merge them with the _df files to have both the stats as well as annotations
library(dplyr)
homeranno.age.deseq2 <- read.delim2("/Users/alexanderwhitehead/Library/Mobile Documents/com~apple~CloudDocs/Documents/PhD/RNASeq_stuff/kirkland_aging_ATAC/for_HOMER_anno/Age_DESEQ2_HOMER_Anno.txt", sep ="\t", header = TRUE )
homeranno.age.edger <- read.delim2("/Users/alexanderwhitehead/Library/Mobile Documents/com~apple~CloudDocs/Documents/PhD/RNASeq_stuff/kirkland_aging_ATAC/for_HOMER_anno/Age_EDGER_HOMER_Anno.txt", sep ="\t", header = TRUE )
homeranno.lamb.deseq2 <- read.delim2("/Users/alexanderwhitehead/Library/Mobile Documents/com~apple~CloudDocs/Documents/PhD/RNASeq_stuff/kirkland_aging_ATAC/for_HOMER_anno/LamB_DESEQ2_HOMER_Anno.txt", sep ="\t", header = TRUE )
homeranno.lamb.edger <- read.delim2("/Users/alexanderwhitehead/Library/Mobile Documents/com~apple~CloudDocs/Documents/PhD/RNASeq_stuff/kirkland_aging_ATAC/for_HOMER_anno/LamB_EDGER_HOMER_Anno.txt", sep ="\t", header = TRUE )
homeranno.lamc.deseq2 <- read.delim2("/Users/alexanderwhitehead/Library/Mobile Documents/com~apple~CloudDocs/Documents/PhD/RNASeq_stuff/kirkland_aging_ATAC/for_HOMER_anno/LamC_DESEQ2_HOMER_Anno.txt", sep ="\t", header = TRUE )
homeranno.lamc.edger <- read.delim2("/Users/alexanderwhitehead/Library/Mobile Documents/com~apple~CloudDocs/Documents/PhD/RNASeq_stuff/kirkland_aging_ATAC/for_HOMER_anno/LamC_EDGER_HOMER_Anno.txt", sep ="\t", header = TRUE )

homer_colnames <- colnames(homeranno.age.deseq2)
homer_colnames[1] <- "peak" 

colnames(homeranno.age.deseq2) <- homer_colnames
colnames(homeranno.age.edger) <- homer_colnames
colnames(homeranno.lamb.deseq2) <- homer_colnames
colnames(homeranno.lamb.edger) <- homer_colnames
colnames(homeranno.lamc.deseq2) <- homer_colnames
colnames(homeranno.lamc.edger) <- homer_colnames

data.age.deseq2_df$peak <- as.integer(data.age.deseq2_df$peak)
data.age.edger_df$peak <- as.integer(data.age.edger_df$peak)
data.lamb.deseq2_df$peak <- as.integer(data.lamb.deseq2_df$peak)
data.lamb.edger_df$peak <- as.integer(data.lamb.edger_df$peak)
data.lamc.deseq2_df$peak <- as.integer(data.lamc.deseq2_df$peak)
data.lamc.edger_df$peak <- as.integer(data.lamc.edger_df$peak)

age_deseq2 <- inner_join(homeranno.age.deseq2, data.age.deseq2_df, c("peak"))
age_edger <- inner_join(homeranno.age.edger, data.age.edger_df, c("peak"))
lamb_deseq2 <- inner_join(homeranno.lamb.deseq2, data.lamb.deseq2_df, c("peak"))
lamb_edger <- inner_join(homeranno.lamb.edger, data.lamb.edger_df, c("peak"))
lamc_deseq2 <- inner_join(homeranno.lamc.deseq2, data.lamc.deseq2_df, c("peak"))
lamc_edger <- inner_join(homeranno.lamc.edger, data.lamc.edger_df, c("peak"))

# Remove redundant columns
# keep: 1:4, 8:14,16:19,25:30
age_deseq2 <- age_deseq2[c(1:5,8:14,16:19,25:30)]
age_edger <- age_edger[c(1:4,5:14,16:19,25:30)]
lamb_deseq2 <- lamb_deseq2[c(1:5,8:14,16:19,25:30)]
lamb_edger <- lamb_edger[c(1:5,8:14,16:19,25:30)]
lamc_deseq2 <- lamc_deseq2[c(1:5,8:14,16:19,25:30)]
lamc_edger <- lamc_edger[c(1:5,8:14,16:19,25:30)]

#output the final tables for easy viewing
write.table(age_deseq2, file = "Kirkland_Age_DESEQ2.bed", sep = "\t", row.names=F, col.names = T, qmethod = "double")
write.table(age_edger, file = "Kirkland_Age_EdgeR.bed", sep = "\t", row.names=F, col.names = T, qmethod = "double")
write.table(lamb_deseq2, file = "Kirkland_LamB_DESEQ2.bed", sep = "\t", row.names=F, col.names = T, qmethod = "double")
write.table(lamb_edger, file = "Kirkland_LamB_EdgeR.bed", sep = "\t", row.names=F, col.names = T, qmethod = "double")
write.table(lamc_deseq2, file = "Kirkland_LamC_DESEQ2.bed", sep = "\t", row.names=F, col.names = T, qmethod = "double")
write.table(lamc_edger, file = "Kirkland_LamC_EdgeR.bed", sep = "\t", row.names=F, col.names = T, qmethod = "double")


# Filtering and outputting for MOTIF analysis -----------------------------

Age_ATAC_DESEQ2_DARS_DOWN_W_AGE <-age_deseq2[which(age_deseq2$Fold>0 & age_deseq2$FDR<0.1),] 
Age_ATAC_DESEQ2_DARS_UP_W_AGE <- age_deseq2[which(age_deseq2$Fold<0 & age_deseq2$FDR<0.1),] 
Age_ATAC_EdgeR_DARS_DOWN_W_AGE <-age_edger[which(age_edger$Fold>0 & age_edger$FDR<0.1),] 
Age_ATAC_EdgeR_DARS_UP_W_AGE <- age_edger[which(age_edger$Fold<0 & age_edger$FDR<0.1),] 

LamB_ATAC_DESEQ2_DARS_DOWN_W_KO <-lamb_deseq2[which(lamb_deseq2$Fold>0 & lamb_deseq2$FDR<0.1),] 
LamB_ATAC_DESEQ2_DARS_UP_W_KO <- lamb_deseq2[which(lamb_deseq2$Fold<0 & lamb_deseq2$FDR<0.1),] 
LamB_ATAC_EdgeR_DARS_DOWN_W_KO <-lamb_edger[which(lamb_edger$Fold>0 & lamb_edger$FDR<0.1),] 
LamB_ATAC_EdgeR_DARS_UP_W_KO <- lamb_edger[which(lamb_edger$Fold<0 & lamb_edger$FDR<0.1),] 

LamC_ATAC_DESEQ2_DARS_DOWN_W_KO <-lamc_deseq2[which(lamc_deseq2$Fold>0 & lamc_deseq2$FDR<0.1),] 
LamC_ATAC_DESEQ2_DARS_UP_W_KO <- lamc_deseq2[which(lamc_deseq2$Fold<0 & lamc_deseq2$FDR<0.1),] 
LamC_ATAC_EdgeR_DARS_DOWN_W_KO <-lamc_edger[which(lamc_edger$Fold>0 & lamc_edger$FDR<0.1),] 
LamC_ATAC_EdgeR_DARS_UP_W_KO <- lamc_edger[which(lamc_edger$Fold<0 & lamc_edger$FDR<0.1),] 

# Create objects of DARS, agnostic of UP vs Down
Age_DARS_EdgeR <- age_edger[which(age_edger$FDR<0.1),] 
LamB_DARS_EdgeR <- lamb_edger[which(lamb_edger$FDR<0.1),] 
LamC_DARS_EdgeR <- lamc_edger[which(lamc_edger$FDR<0.1),] 

# Find common DARS between each condition by Nearest.Refseq
# Age and LamC
Age_and_LamC_EdgeR <- inner_join(Age_DARS_EdgeR, LamC_DARS_EdgeR, by = "Nearest.Refseq")
test1 <- unique(Age_and_LamC_EdgeR$Gene.Name.x)
write.table(test1, file = "Age_and_LamC_DARS_Name.csv",sep = "\t", row.names=F, col.names = F)
# Age and LamB
Age_and_LamB_EdgeR <- inner_join(Age_DARS_EdgeR, LamB_DARS_EdgeR, by = "Nearest.Refseq")
test2 <- unique(Age_and_LamB_EdgeR$Gene.Name.x)
write.table(test2, file = "Age_and_LamB_DARS_Name.csv",sep = "\t", row.names=F, col.names = F)
# LamC and Lamb
LamC_and_LamB_EdgeR <- inner_join(LamC_DARS_EdgeR, LamB_DARS_EdgeR, by = "Nearest.Refseq")
test3 <- unique(LamC_and_LamB_EdgeR$Gene.Name.x)
write.table(test3, file = "LamC_and_LamB_DARS_Name.csv",sep = "\t", row.names=F, col.names = F)

# Create blacklist of regions that are differentially accessible accross all comparisons
DARS_blacklist <- inner_join(Age_and_LamB_EdgeR,Age_and_LamC_EdgeR, by = "Nearest.Refseq")
DARS_blacklist <- inner_join(DARS_blacklist, LamC_and_LamB_EdgeR, by = "Nearest.Refseq")
test <- unique(DARS_blacklist$Gene.Name.x.x)

#test plotting FC age vs FC lamC
pdf("Correlation_Plots.pdf")
plot(-1*Age_and_LamC_EdgeR$Fold.x,-1*Age_and_LamC_EdgeR$Fold.y, main = "FC Age vs. LamCKO", xlab="Age FC (5wk/1wk)", ylab="LamCKO FC (KO/wt)")
plot(-1*Age_and_LamB_EdgeR$Fold.x,-1*Age_and_LamB_EdgeR$Fold.y, main = "FC Age vs. LamBKO", xlab="Age FC (5wk/1wk)",ylab="LamBKO FC (KO/wt)")
plot(-1*LamC_and_LamB_EdgeR$Fold.x,-1*LamC_and_LamB_EdgeR$Fold.y, main = "FC LamCKO vs. LamBKO",xlab="LamCKO FC (KO/wt)",ylab="LamBKO FC (KO/wt)")
ggplot(Age_and_LamC_EdgeR) +
         geom_point(aes(x =-1* Age_and_LamC_EdgeR$Fold.x, y =-1*Age_and_LamC_EdgeR$Fold.y))+
         ggtitle("FC Age vs. LamCKO, Zoomed-In") +
         xlab("Log2 Fold Change Age (5wk/1wk)") + 
         ylab("Log2 Fold Change LamC (KO/wt) ") +
         scale_y_continuous(limits = c(-2,2)) +
         scale_x_continuous(limits = c(-2,2)) +
         #geom_text_repel(aes(x = -1*Fold,
                             #y = -log10(FDR),
                             #label = ifelse(threshold == T, Gene.Name,"")),
                         #max.overlaps = 10,
                         #force = 2) +
         
         theme_few()+ 
         theme(legend.position = "none",
               plot.title = element_text(size = rel(1.5), hjust = 0.5),
               axis.title = element_text(size = rel(1.25))) 
dev.off()

#output the dataframes with FC and values
write.csv(Age_and_LamC_EdgeR, file = "Age_and_LamC_EdgeR.csv", row.names=F, col.names = T)
write.csv(Age_and_LamB_EdgeR, file = "Age_and_LamB_EdgeR.csv", row.names=F, col.names = T)
write.csv(LamC_and_LamB_EdgeR, file = "LamC_and_LamB_EdgeR.csv", row.names=F, col.names = T)

# we need to rearrange the columns to work with HOMER's findMotifsGenome.pl script
motif_colnames <- c("Chr","Start","End","peak","End","Strand")
Age_ATAC_DESEQ2_DARS_DOWN_W_AGE_1 <-Age_ATAC_DESEQ2_DARS_DOWN_W_AGE[,motif_colnames]
Age_ATAC_DESEQ2_DARS_UP_W_AGE_1 <-Age_ATAC_DESEQ2_DARS_UP_W_AGE[,motif_colnames]
Age_ATAC_EdgeR_DARS_DOWN_W_AGE_1 <-Age_ATAC_EdgeR_DARS_DOWN_W_AGE[,motif_colnames]
Age_ATAC_EdgeR_DARS_UP_W_AGE_1 <-Age_ATAC_EdgeR_DARS_UP_W_AGE[,motif_colnames]
LamB_ATAC_DESEQ2_DARS_DOWN_W_KO_1 <-LamB_ATAC_DESEQ2_DARS_DOWN_W_KO[,motif_colnames]
LamB_ATAC_DESEQ2_DARS_UP_W_KO_1 <-LamB_ATAC_DESEQ2_DARS_UP_W_KO[,motif_colnames]
LamB_ATAC_EdgeR_DARS_DOWN_W_KO_1 <-LamB_ATAC_EdgeR_DARS_DOWN_W_KO[,motif_colnames]
LamB_ATAC_EdgeR_DARS_UP_W_KO_1 <-LamB_ATAC_EdgeR_DARS_UP_W_KO[,motif_colnames]
LamC_ATAC_DESEQ2_DARS_DOWN_W_KO_1 <-LamC_ATAC_DESEQ2_DARS_DOWN_W_KO[,motif_colnames]
LamC_ATAC_DESEQ2_DARS_UP_W_KO_1 <-LamC_ATAC_DESEQ2_DARS_UP_W_KO[,motif_colnames]
LamC_ATAC_EdgeR_DARS_DOWN_W_KO_1 <-LamC_ATAC_EdgeR_DARS_DOWN_W_KO[,motif_colnames]
LamC_ATAC_EdgeR_DARS_UP_W_KO_1 <-LamC_ATAC_EdgeR_DARS_UP_W_KO[,motif_colnames]

# output
write.table(Age_ATAC_DESEQ2_DARS_DOWN_W_AGE_1, file = "Age_ATAC_DESEQ2_DARS_DOWN_W_AGE.bed",sep = "\t", row.names=F, col.names = F)
write.table(Age_ATAC_DESEQ2_DARS_UP_W_AGE_1, file = "Age_ATAC_DESEQ2_DARS_UP_W_AGE.bed",sep = "\t", row.names=F,col.names = F)
write.table(Age_ATAC_EdgeR_DARS_DOWN_W_AGE_1, file = "Age_ATAC_EdgeR_DARS_DOWN_W_AGE.bed",sep = "\t", row.names=F,col.names = F)
write.table(Age_ATAC_EdgeR_DARS_UP_W_AGE_1, file = "Age_ATAC_EdgeR_DARS_UP_W_AGE.bed",sep = "\t", row.names=F,col.names = F)

write.table(LamB_ATAC_DESEQ2_DARS_DOWN_W_KO_1, file = "LamB_ATAC_DESEQ2_DARS_DOWN_W_KO.bed",sep = "\t", row.names=F,col.names = F)
write.table(LamB_ATAC_DESEQ2_DARS_UP_W_KO_1, file = "LamB_ATAC_DESEQ2_DARS_UP_W_KO.bed",sep = "\t", row.names=F,col.names = F)
write.table(LamB_ATAC_EdgeR_DARS_DOWN_W_KO_1, file = "LamB_ATAC_EdgeR_DARS_DOWN_W_KO.bed",sep = "\t", row.names=F,col.names = F)
write.table(LamB_ATAC_EdgeR_DARS_UP_W_KO_1, file = "LamB_ATAC_EdgeR_DARS_UP_W_KO.bed",sep = "\t", row.names=F,col.names = F)


write.table(LamC_ATAC_DESEQ2_DARS_DOWN_W_KO_1, file = "LamC_ATAC_DESEQ2_DARS_DOWN_W_KO.bed",sep = "\t", row.names=F,col.names = F)
write.table(LamC_ATAC_DESEQ2_DARS_UP_W_KO_1, file = "LamC_ATAC_DESEQ2_DARS_UP_W_KO.bed",sep = "\t", row.names=F,col.names = F)
write.table(LamC_ATAC_EdgeR_DARS_DOWN_W_KO_1, file = "LamC_ATAC_EdgeR_DARS_DOWN_W_KO.bed",sep = "\t", row.names=F,col.names = F)
write.table(LamC_ATAC_EdgeR_DARS_UP_W_KO_1, file = "LamC_ATAC_EdgeR_DARS_UP_W_KO.bed",sep = "\t", row.names=F,col.names = F)

# For whatever reason, you need to open the files in excel and delete the last column manually,
# then change column 5 to all ones, and format as "number" in excel (should look like 1.00)
# After this, move the bed files into your ATAC seq directory for command line operations and run HOMER findMotifsGenome.pl

# Below we will do the meta analysis of our comparisons

# we are going to find the intersection of DARS that decrease/increase with age + lamCKO
DESEQ2_Meta_Dwn_w_age_LamCKO <- inner_join(LamC_ATAC_DESEQ2_DARS_DOWN_W_KO, Age_ATAC_DESEQ2_DARS_DOWN_W_AGE, by = "Nearest.Refseq")
DESEQ2_Meta_Up_w_age_LamCKO <- inner_join(LamC_ATAC_DESEQ2_DARS_UP_W_KO, Age_ATAC_DESEQ2_DARS_UP_W_AGE, by = "Nearest.Refseq")
EdgeR_Meta_Dwn_w_age_LamCKO <- inner_join(LamC_ATAC_EdgeR_DARS_DOWN_W_KO, Age_ATAC_EdgeR_DARS_DOWN_W_AGE,by = "Nearest.Refseq")
EdgeR_Meta_Up_w_age_LamCKO <- inner_join(LamC_ATAC_EdgeR_DARS_UP_W_KO, Age_ATAC_EdgeR_DARS_UP_W_AGE,by = "Nearest.Refseq")

#output
write.table(DESEQ2_Meta_Dwn_w_age_LamCKO, file = "LamC_and_Age_Down_w_Age_and_KO_DESEQ2.txt",sep = "\t", row.names=F)
write.table(DESEQ2_Meta_Up_w_age_LamCKO, file = "LamC_and_Age_Up_w_Age_and_KO_DESEQ2.txt",sep = "\t", row.names=F)
write.table(EdgeR_Meta_Dwn_w_age_LamCKO, file = "LamC_and_Age_Down_w_Age_and_KO_EdgeR.txt",sep = "\t", row.names=F)
write.table(EdgeR_Meta_Up_w_age_LamCKO, file = "LamC_and_Age_Up_w_Age_and_KO_EdgeR.txt",sep = "\t", row.names=F)

# now we will find the intersection of DARS that increase/decrease w age and lamBKO
DESEQ2_Meta_Dwn_w_age_LamBKO <- inner_join(LamB_ATAC_DESEQ2_DARS_DOWN_W_KO, Age_ATAC_DESEQ2_DARS_DOWN_W_AGE, by = "Nearest.Refseq")
DESEQ2_Meta_Up_w_age_LamBKO <- inner_join(LamB_ATAC_DESEQ2_DARS_UP_W_KO, Age_ATAC_DESEQ2_DARS_UP_W_AGE, by = "Nearest.Refseq")
EdgeR_Meta_Dwn_w_age_LamBKO <- inner_join(LamB_ATAC_EdgeR_DARS_DOWN_W_KO, Age_ATAC_EdgeR_DARS_DOWN_W_AGE,by = "Nearest.Refseq")
EdgeR_Meta_Up_w_age_LamBKO <- inner_join(LamB_ATAC_EdgeR_DARS_UP_W_KO, Age_ATAC_EdgeR_DARS_UP_W_AGE,by = "Nearest.Refseq")

write.table(DESEQ2_Meta_Dwn_w_age_LamBKO, file = "LamB_and_Age_Down_w_Age_and_KO_DESEQ2.txt",sep = "\t", row.names=F)
write.table(DESEQ2_Meta_Up_w_age_LamBKO, file = "LamB_and_Age_Up_w_Age_and_KO_DESEQ2.txt",sep = "\t", row.names=F)
write.table(EdgeR_Meta_Dwn_w_age_LamBKO, file = "LamB_and_Age_Down_w_Age_and_KO_EdgeR.txt",sep = "\t", row.names=F)
write.table(EdgeR_Meta_Up_w_age_LamBKO, file = "LamB_and_Age_Up_w_Age_and_KO_EdgeR.txt",sep = "\t", row.names=F)

# now we will find the intersection of DARS that increase/decrease w lamCKO and lamBKO
DESEQ2_Meta_Dwn_w_LamCKO_LamBKO <- inner_join(LamB_ATAC_DESEQ2_DARS_DOWN_W_KO, LamC_ATAC_DESEQ2_DARS_DOWN_W_KO, by = "Nearest.Refseq")
DESEQ2_Meta_Up_w_LamCKO_LamBKO <- inner_join(LamB_ATAC_DESEQ2_DARS_UP_W_KO, LamC_ATAC_DESEQ2_DARS_UP_W_KO, by = "Nearest.Refseq")
EdgeR_Meta_Dwn_w_LamCKO_LamBKO <- inner_join(LamB_ATAC_EdgeR_DARS_DOWN_W_KO, LamC_ATAC_EdgeR_DARS_DOWN_W_KO,by = "Nearest.Refseq")
EdgeR_Meta_Up_w_LamCKO_LamBKO <- inner_join(LamB_ATAC_EdgeR_DARS_UP_W_KO, LamC_ATAC_EdgeR_DARS_UP_W_KO,by = "Nearest.Refseq")

write.table(DESEQ2_Meta_Dwn_w_LamCKO_LamBKO, file = "LamB_and_LamC_Down_w_KO_DESEQ2.txt",sep = "\t", row.names=F)
write.table(DESEQ2_Meta_Up_w_LamCKO_LamBKO, file = "LamB_and_LamC_Up_w_KO_DESEQ2.txt",sep = "\t", row.names=F)
write.table(EdgeR_Meta_Dwn_w_LamCKO_LamBKO, file = "LamB_and_LamC_Down_w_KO_EdgeR.txt",sep = "\t", row.names=F)
write.table(EdgeR_Meta_Up_w_LamCKO_LamBKO, file = "LamB_and_LamC_Up_w_KO_EdgeR.txt",sep = "\t", row.names=F)


save.image(file = "Kirkland_Global_Objects.RData")
load("Kirkland_Global_Objects.RData")
#dba.save(Age_ATAC,file="Age_ATAC")
#writeLines(capture.output(sessionInfo()), "sessionInfo.txt")


# Make Annotated Heatmaps and Volcano Plots --------------------------------------------------

# Take these DARS and make a pheatmap
library(gplots)
library(pheatmap)
library(plyr)
library(tibble)
library(tidyr)


#Natalie wants a yellow and blue theme so we need to design the color palette
myColor <- colorRampPalette(c("#FBA919","white","#3C54A4"))(30)

# Convert the total data for each comparison into a tibble for convenience with ggplot2
age_deseq2_mat <- as_tibble(age_deseq2)

# sort by FDR from smallest to greatest
age_deseq2_mat <- arrange(age_deseq2_mat, FDR)

#Generate a heatmap of peak region densities for top 500 differentially accesible regions
pdf("Heatmap_Age_DESEQ2_scale.pdf", width = 10, height =15)
pheatmap((age_deseq2_mat[c(1:500),c(17:18)]),   #the matrix is already sorted by FDR from smallest to largest
         color = myColor,
         show_rownames = T,
         cluster_rows = T,
         cluster_cols = T,
         labels_row = age_deseq2_mat$Gene.Name,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         border_color = F,
         cex = 0.4,
         scale = "row",
         fontsize_row = 8,
         labels_col = c("1WK","5WK"),
         legend = T) 
# you can also use scale = row to normalize each row to median and plot std dev
dev.off()

pdf("Heatmap_Age_DESEQ2.pdf", width = 10, height =15)
pheatmap((age_deseq2_mat[c(1:500),c(17:18)]),   #the matrix is already sorted by FDR from smallest to largest
         color = myColor,
         show_rownames = T,
         cluster_rows = T,
         cluster_cols = T,
         labels_row = age_deseq2_mat$Gene.Name,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         border_color = F,
         cex = 0.4,
         fontsize_row = 8,
         labels_col = c("1WK","5WK"),
         legend = T) 
# you can also use scale = row to normalize each row to median and plot std dev
dev.off()

# Convert the total data for each comparison into a tibble for convenience with ggplot2
age_edger_mat <- as_tibble(age_edger)

# sort by FDR from smallest to greatest
age_edger_mat <- arrange(age_edger_mat, FDR)

#Generate a heatmap of peak region densities for top 500 differentially accesible regions
pdf("Heatmap_Age_EdgeR_scale.pdf", width = 10, height =15)
pheatmap((age_edger_mat[c(1:500),c(17:18)]),   #the matrix is already sorted by FDR from smallest to largest
         show_rownames = T,
         color = myColor,
         cluster_rows = T,
         cluster_cols = T,
         labels_row = age_edger_mat$Gene.Name,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         border_color = F,
         cex = 0.4,
         scale = "row",
         fontsize_row = 8,
         labels_col = c("1WK","5WK"),
         legend = T) 
# you can also use scale = row to normalize each row to median and plot std dev
dev.off()

pdf("Heatmap_Age_EdgeR.pdf", width = 10, height =15)
pheatmap((age_edger_mat[c(1:500),c(17:18)]),   #the matrix is already sorted by FDR from smallest to largest
         show_rownames = T,
         color = myColor,
         cluster_rows = T,
         cluster_cols = T,
         labels_row = age_edger_mat$Gene.Name,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         border_color = F,
         cex = 0.4,
         fontsize_row = 8,
         labels_col = c("1WK","5WK"),
         legend = T) 
# you can also use scale = row to normalize each row to median and plot std dev
dev.off()


# Now we will repeat this for the other comparisons
lamb_deseq2_mat <- as_tibble(lamb_deseq2)
# sort by FDR from smallest to greatest
lamb_deseq2_mat <- arrange(lamb_deseq2_mat, FDR)

#Generate a heatmap of peak region densities for top 500 differentially accesible regions
pdf("Heatmap_LamB_DESEQ2_scale.pdf", width = 10, height =15)
pheatmap((lamb_deseq2_mat[c(1:500),c(17:18)]),   #the matrix is already sorted by FDR from smallest to largest
         show_rownames = T,
         color = myColor,
         cluster_rows = T,
         cluster_cols = T,
         labels_row = lamb_deseq2_mat$Gene.Name,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         border_color = F,
         cex = 0.4,
         scale = "row",
         fontsize_row = 8,
         labels_col = c("WT","LamB iR"),
         legend = T) 
# you can also use scale = row to normalize each row to median and plot std dev
dev.off()

pdf("Heatmap_LamB_DESEQ2.pdf", width = 10, height =15)
pheatmap((lamb_deseq2_mat[c(1:500),c(17:18)]),   #the matrix is already sorted by FDR from smallest to largest
         show_rownames = T,
         color = myColor,
         cluster_rows = T,
         cluster_cols = T,
         labels_row = lamb_deseq2_mat$Gene.Name,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         border_color = F,
         cex = 0.4,
         fontsize_row = 8,
         labels_col = c("WT","LamB iR"),
         legend = T) 
# you can also use scale = row to normalize each row to median and plot std dev
dev.off()

# Next
lamb_edger_mat <- as_tibble(lamb_edger)
# sort by FDR from smallest to greatest
lamb_edger_mat <- arrange(lamb_edger_mat, FDR)

#Generate a heatmap of peak region densities for top 500 differentially accesible regions
pdf("Heatmap_LamB_EdgeR_scale.pdf", width = 10, height =15)
pheatmap((lamb_edger_mat[c(1:500),c(17:18)]),   #the matrix is already sorted by FDR from smallest to largest
         show_rownames = T,
         color = myColor,
         cluster_rows = T,
         cluster_cols = T,
         labels_row = lamb_edger_mat$Gene.Name,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         border_color = F,
         cex = 0.4,
         scale = "row",
         fontsize_row = 8,
         labels_col = c("WT","LamB iR"),
         legend = T) 
# you can also use scale = row to normalize each row to median and plot std dev
dev.off()

pdf("Heatmap_LamB_EdgeR.pdf", width = 10, height =15)
pheatmap((lamb_edger_mat[c(1:500),c(17:18)]),   #the matrix is already sorted by FDR from smallest to largest
         show_rownames = T,
         color = myColor,
         cluster_rows = T,
         cluster_cols = T,
         labels_row = lamb_edger_mat$Gene.Name,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         border_color = F,
         cex = 0.4,
         fontsize_row = 8,
         labels_col = c("WT","LamB iR"),
         legend = T) 
# you can also use scale = row to normalize each row to median and plot std dev
dev.off()

# Now we will repeat this for the other comparisons
lamc_deseq2_mat <- as_tibble(lamc_deseq2)
# sort by FDR from smallest to greatest
lamc_deseq2_mat <- arrange(lamc_deseq2_mat, FDR)

#Generate a heatmap of peak region densities for top 500 differentially accesible regions
pdf("Heatmap_LamC_DESEQ2_scale.pdf", width = 10, height =15)
pheatmap((lamc_deseq2_mat[c(1:500),c(17:18)]),   #the matrix is already sorted by FDR from smallest to largest
         show_rownames = T,
         color = myColor,
         cluster_rows = T,
         cluster_cols = T,
         labels_row = lamc_deseq2_mat$Gene.Name,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         border_color = F,
         cex = 0.4,
         scale = "row",
         fontsize_row = 8,
         labels_col = c("WT","LamC iR"),
         legend = T) 
# you can also use scale = row to normalize each row to median and plot std dev
dev.off()

pdf("Heatmap_LamC_DESEQ2.pdf", width = 10, height =15)
pheatmap((lamc_deseq2_mat[c(1:500),c(17:18)]),   #the matrix is already sorted by FDR from smallest to largest
         show_rownames = T,
         color = myColor,
         cluster_rows = T,
         cluster_cols = T,
         labels_row = lamc_deseq2_mat$Gene.Name,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         border_color = F,
         cex = 0.4,
         fontsize_row = 8,
         labels_col = c("WT","LamC iR"),
         legend = T) 
# you can also use scale = row to normalize each row to median and plot std dev
dev.off()

# Next
lamc_edger_mat <- as_tibble(lamc_edger)
# sort by FDR from smallest to greatest
lamc_edger_mat <- arrange(lamc_edger_mat, FDR)

#Generate a heatmap of peak region densities for top 500 differentially accesible regions
pdf("Heatmap_LamC_EdgeR_scale.pdf", width = 10, height =15)
pheatmap((lamc_edger_mat[c(1:500),c(17:18)]),   #the matrix is already sorted by FDR from smallest to largest
         show_rownames = T,
         color = myColor,
         cluster_rows = T,
         cluster_cols = T,
         labels_row = lamc_edger_mat$Gene.Name,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         border_color = F,
         cex = 0.4,
         scale = "row",
         fontsize_row = 8,
         labels_col = c("WT","LamC iR"),
         legend = T) 
# you can also use scale = row to normalize each row to median and plot std dev
dev.off()

pdf("Heatmap_LamC_EdgeR.pdf", width = 10, height =15)
pheatmap((lamc_edger_mat[c(1:500),c(17:18)]),   #the matrix is already sorted by FDR from smallest to largest
         show_rownames = T,
         color = myColor,
         cluster_rows = T,
         cluster_cols = T,
         labels_row = lamc_edger_mat$Gene.Name,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         border_color = F,
         cex = 0.4,
         fontsize_row = 8,
         labels_col = c("WT","LamC iR"),
         legend = T) 
# you can also use scale = row to normalize each row to median and plot std dev
dev.off()


# we DIY a volcano plot using ggplot2
library(ggplot2)
library(ggrepel)
library(ggthemes)

threshold_OE <- age_edger_mat$FDR < 0.1
length(which(threshold_OE))
age_edger_mat$threshold <- threshold_OE 

# Natalie preferences: standardize x (-3,3) and y (0,20) axes, color scale purple = down, green= up for volc
# and then yellow = low and blue = high for heatmaps

age_edger_mat <- age_edger_mat %>%
  mutate(threshold_OE = factor(case_when(Fold > 0.32 & FDR < 0.1 ~"Down",
                                         Fold < -0.32 & FDR < 0.1 ~ "Up",
                                         Threshold = TRUE ~ "NS")))

pdf("Age_Volc_EdgeR.pdf")
ggplot(age_edger_mat) +
  geom_point(aes(x = -1*Fold, y = -log10(FDR), colour = threshold_OE), size = 1) +
  scale_color_manual(name = "Threshold",
                     values=c("NS" = "black" , "Down" = "#C3A2CC" , "Up" = "#59CC8D")) + 
  ggtitle("Age, Right = Up with Aging") +
  xlab("Log2 Fold Change") + 
  ylab("-Log10 Adjusted P-Value") +
  scale_y_continuous(limits = c(0,20)) +
  scale_x_continuous(limits = c(-3,3)) +
  geom_text_repel(aes(x = -1*Fold,
                      y = -log10(FDR),
                      label = ifelse(threshold == T, Gene.Name,"")),
                  max.overlaps = 10,
                  force = 2) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) +
  geom_vline(xintercept = c(-0.32,0.32), linetype = "dotted") + 
  geom_hline(yintercept = 1, linetype = "dotted")+
  theme_few()
dev.off()



threshold_OE <- age_deseq2_mat$FDR < 0.1
length(which(threshold_OE))
age_deseq2_mat$threshold <- threshold_OE 

age_deseq2_mat <- age_deseq2_mat %>%
  mutate(threshold_OE = factor(case_when(Fold > 0.32 & FDR < 0.1 ~"Down",
                                         Fold < -0.32 & FDR < 0.1 ~ "Up",
                                         Threshold = TRUE ~ "NS")))

pdf("Age_Volc_DESEQ2.pdf")
ggplot(age_deseq2_mat) +
  geom_point(aes(x = -1*Fold, y = -log10(FDR), colour = threshold_OE), size = 1) +
  scale_color_manual(name = "Threshold",
                     values=c("NS" = "black" , "Down" = "#C3A2CC" , "Up" = "#59CC8D")) + 
  ggtitle("Age, Right = Up with Aging") +
  xlab("Log2 Fold Change") + 
  ylab("-Log10 Adjusted P-Value") +
  scale_y_continuous(limits = c(0,20)) +
  scale_x_continuous(limits = c(-3,3)) +
  geom_text_repel(aes(x = -1*Fold,
                      y = -log10(FDR),
                      label = ifelse(threshold == T, Gene.Name,"")),
                  max.overlaps = 10,
                  force = 2) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) +
  geom_vline(xintercept = c(-0.32,0.32), linetype = "dotted") + 
  geom_hline(yintercept = 1, linetype = "dotted")+
  theme_few()
dev.off()


# Repeat for LamB comparison
threshold_OE <- lamb_edger_mat$FDR < 0.1
length(which(threshold_OE))
lamb_edger_mat$threshold <- threshold_OE 

lamb_edger_mat <- lamb_edger_mat %>%
  mutate(threshold_OE = factor(case_when(Fold > 0.32 & FDR < 0.1 ~"Down",
                                         Fold < -0.32 & FDR < 0.1 ~ "Up",
                                         Threshold = TRUE ~ "NS")))

pdf("LamB_Volc_EdgeR.pdf")
ggplot(lamb_edger_mat) +
  geom_point(aes(x = -1*Fold, y = -log10(FDR), colour = threshold_OE), size = 1) +
  scale_color_manual(name = "Threshold",
                     values=c("NS" = "black" , "Down" = "#C3A2CC" , "Up" = "#59CC8D")) + 
  ggtitle("LamBKO, Right = Up with KO") +
  xlab("Log2 Fold Change") + 
  ylab("-Log10 Adjusted P-Value") +
  scale_y_continuous(limits = c(0,20)) +
  scale_x_continuous(limits = c(-3,3)) +
  geom_text_repel(aes(x = -1*Fold,
                      y = -log10(FDR),
                      label = ifelse(threshold == T, Gene.Name,"")),
                  max.overlaps = 10,
                  force = 2) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) +
  geom_vline(xintercept = c(-0.32,0.32), linetype = "dotted") + 
  geom_hline(yintercept = 1, linetype = "dotted")+
  theme_few()
dev.off()



threshold_OE <- lamb_deseq2_mat$FDR < 0.1 
length(which(threshold_OE))
lamb_deseq2_mat$threshold <- threshold_OE 

lamb_deseq2_mat <- lamb_deseq2_mat %>%
  mutate(threshold_OE = factor(case_when(Fold > 0.32 & FDR < 0.1 ~"Down",
                                         Fold < -0.32 & FDR < 0.1 ~ "Up",
                                         Threshold = TRUE ~ "NS")))

pdf("LamB_Volc_DESEQ2.pdf")
ggplot(lamb_deseq2_mat) +
  geom_point(aes(x = -1*Fold, y = -log10(FDR), colour = threshold_OE), size = 1) +
  scale_color_manual(name = "Threshold",
                     values=c("NS" = "black" , "Down" = "#C3A2CC" , "Up" = "#59CC8D")) + 
  ggtitle("LamBKO, Right = Up with KO") +
  xlab("Log2 Fold Change") + 
  ylab("-Log10 Adjusted P-Value") +
  scale_y_continuous(limits = c(0,20)) +
  scale_x_continuous(limits = c(-3,3)) +
  geom_text_repel(aes(x = -1*Fold,
                      y = -log10(FDR),
                      label = ifelse(threshold == T, Gene.Name,"")),
                  max.overlaps = 10,
                  force = 2) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) +
  geom_vline(xintercept = c(-0.32,0.32), linetype = "dotted") + 
  geom_hline(yintercept = 1, linetype = "dotted")+
  theme_few()
dev.off()


# Repeat for LamC comparison
threshold_OE <- lamc_edger_mat$FDR < 0.1
length(which(threshold_OE))
lamc_edger_mat$threshold <- threshold_OE 

lamc_edger_mat <- lamc_edger_mat %>%
  mutate(threshold_OE = factor(case_when(Fold > 0.32 & FDR < 0.1 ~"Down",
                                         Fold < -0.32 & FDR < 0.1 ~ "Up",
                                         Threshold = TRUE ~ "NS")))

pdf("LamC_Volc_EdgeR.pdf")
ggplot(lamc_edger_mat) +
  geom_point(aes(x = -1*Fold, y = -log10(FDR), colour = threshold_OE), size = 1) +
  scale_color_manual(name = "Threshold",
                     values=c("NS" = "black" , "Down" = "#C3A2CC" , "Up" = "#59CC8D")) + 
  ggtitle("LamCKO, Right = Up with KO") +
  xlab("Log2 Fold Change") + 
  ylab("-Log10 Adjusted P-Value") +
  scale_y_continuous(limits = c(0,20)) +
  scale_x_continuous(limits = c(-3,3)) +
  geom_text_repel(aes(x = -1*Fold,
                      y = -log10(FDR),
                      label = ifelse(threshold == T, Gene.Name,"")),
                  max.overlaps = 10,
                  force = 2) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) +
  geom_vline(xintercept = c(-0.32,0.32), linetype = "dotted") + 
  geom_hline(yintercept = 1, linetype = "dotted")+
  theme_few()
dev.off()


threshold_OE <- lamc_deseq2_mat$FDR < 0.1
length(which(threshold_OE))
lamc_deseq2_mat$threshold <- threshold_OE 

lamc_deseq2_mat <- lamc_deseq2_mat %>%
  mutate(threshold_OE = factor(case_when(Fold > 0.32 & FDR < 0.1 ~"Down",
                                         Fold < -0.32 & FDR < 0.1 ~ "Up",
                                         Threshold = TRUE ~ "NS")))

pdf("LamC_Volc_DESEQ2.pdf")
ggplot(lamc_deseq2_mat) +
  geom_point(aes(x = -1*Fold, y = -log10(FDR), colour = threshold_OE), size = 1) +
  scale_color_manual(name = "Threshold",
                     values=c("NS" = "black" , "Down" = "#C3A2CC" , "Up" = "#59CC8D")) + 
  ggtitle("LamCKO, Right = Up with KO") +
  xlab("Log2 Fold Change") + 
  ylab("-Log10 Adjusted P-Value") +
  scale_y_continuous(limits = c(0,20)) +
  scale_x_continuous(limits = c(-3,3)) +
  geom_text_repel(aes(x = -1*Fold,
                      y = -log10(FDR),
                      label = ifelse(threshold == T, Gene.Name,"")),
                  max.overlaps = 10,
                  force = 2) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) +
  geom_vline(xintercept = c(-0.32,0.32), linetype = "dotted") + 
  geom_hline(yintercept = 1, linetype = "dotted")+
  theme_few()
dev.off()




















# Functional Analysis Using CHIPSeeker ------------------------------------


#Now take these DEGS and send to KEGG + GO analysis
# We will want to use the GRanges structure for CHIPSeeker - here we called it data.peaksAnno
library(ChIPseeker)
library(ReactomePA)
library(clusterProfiler)
library(biomaRt)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
#get the TSS regions for your genome
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream =3000)

#see the structure of your DBA object from earlier - this will help you choose peaksets for CHIPSeeker
dba.show(Age_ATAC)

#Retrieve all of the peaksets from your previous list object and store in a new variable
allpeaks.consensus <- dba.peakset(Age_ATAC, bRetrieve = T) #bRetrieve=true just means that you are extracting peaksets instead of adding
allpeaks.1wk <-dba.peakset(Age_ATAC, 1:3, bRetrieve = T)
allpeaks.5wk <-dba.peakset(Age_ATAC, 4:6, bRetrieve = T)
allpeaks.lambwt <-dba.peakset(LamB_ATAC, 1:3, bRetrieve = T)
allpeaks.lambko <-dba.peakset(LamB_ATAC, 4:6, bRetrieve = T)
allpeaks.lamcwt <-dba.peakset(LamC_ATAC, 1:3, bRetrieve = T)
allpeaks.lamcko <-dba.peakset(LamC_ATAC, 4:6, bRetrieve = T)

# Calculate Tag Matrix
tagMatrix_all <- getTagMatrix(allpeaks.consensus, windows = promoter)
tagMat1_all <- getTagMatrix(allpeaks.1wk, windows = promoter)
tagMat5_all <- getTagMatrix(allpeaks.5wk, windows = promoter)
tagMatlambwt_all <- getTagMatrix(allpeaks.lambwt, windows = promoter)
tagMat5_lambko_all <- getTagMatrix(allpeaks.lambko, windows = promoter)
tagMatlamcwt_all <- getTagMatrix(allpeaks.lamcwt, windows = promoter)
tagMat5_lamcko_all <- getTagMatrix(allpeaks.lamcko, windows = promoter)


#Annotate Peaks
Age_ATAC_anno_all <- annotatePeak(allpeaks.consensus, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
Age_ATAC_anno1_all <- annotatePeak(allpeaks.1wk, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
Age_ATAC_anno5_all <- annotatePeak(allpeaks.5wk, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
LamB_ATAC_annowt_all <- annotatePeak(allpeaks.lambwt, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
LamB_ATAC_annoko_all <- annotatePeak(allpeaks.lambko, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
LamC_ATAC_annowt_all <- annotatePeak(allpeaks.lamcwt, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
LamC_ATAC_annoko_all <- annotatePeak(allpeaks.lamcko, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)

#Combine annotated files into a list to plot all at once
All_peaks_tagMatrixList <- list(Onewk=tagMat1_all,Fivewk=tagMat5_all, LamBWT =tagMatlambwt_all , 
                                LamBKO =tagMat5_lambko_all , LamCWT = tagMatlamcwt_all,LamCKO = tagMat5_lamcko_all)

#Generate ATAC-seq Heatmaps and Average Profiles (Histograms for Heatmap)
pdf("TagHeatmaps_All_Peaks.pdf")
tagHeatmap(All_peaks_tagMatrixList, xlim=c(-3000,3000), color=NULL)
plotAvgProf(All_peaks_tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=10, facet="row")
dev.off()

# I should try a covplot
# test_peaks_5 <-dba.peakset(Age_ATAC, 6, bRetrieve = T)
# test_peaks_1 <-dba.peakset(Age_ATAC, 1, bRetrieve = T)
# pdf("test_old.pdf")
# covplot(test_peaks_5, weightCol = "Score", title="ChIP Peaks over Chromosomes")
# dev.off()
# pdf("test_new.pdf")
# covplot(test_peaks_1, weightCol = "Score", title="ChIP Peaks over Chromosomes")
# dev.off()


# # create barplots, vennpie, upsetplots, and distance to TSS for the peaks of each group
pdf("CHIPSeeker_Anno_plots_by_Sample_All_Peaks.pdf")
plotAnnoBar(Age_ATAC_anno_all, title = "All Peaks Consensus")
plotAnnoBar(Age_ATAC_anno1_all, title = "1wk All Peaks Conensus")
plotAnnoBar(Age_ATAC_anno5_all, title = "5wk All Peaks Conensus")
vennpie(Age_ATAC_anno_all)
vennpie(Age_ATAC_anno1_all)
vennpie(Age_ATAC_anno5_all)
upsetplot(Age_ATAC_anno_all) #this shows overlapping annotations
upsetplot(Age_ATAC_anno1_all) 
upsetplot(Age_ATAC_anno5_all) 
dev.off()

# Now we will repeat this process with only the differentially expressed peaks
# Old will mean the peaks unique to 5wk, young will be unique to 1wk
age_diff_peaks <-dba.report(Age_ATAC, 
                            th = 0.1,
                            fold = 0.32,
                            bCounts= TRUE, 
                            bCalledDetail = TRUE,
                            contrast = 1,
                            method = DBA_EDGER)
# the filtering by bGain/bLoss doesn't work, so we will do it manually
age_diff_peaks_up <- age_diff_peaks[age_diff_peaks$Fold <0]
age_diff_peaks_down <- age_diff_peaks[age_diff_peaks$Fold >0]
# these are GRanges by the way :)

pdf("CovPlots_Age.pdf")
covplot(age_diff_peaks_up, weightCol = "Conc_Group2", title="ChIP Peaks over Chromosomes Enriched in 5wk")
covplot(age_diff_peaks_down, weightCol = "Conc_Group1", title="ChIP Peaks over Chromosomes Enriched in 1wk")
dev.off()



#Do annotation
Age_ATAC_anno1_diff <- annotatePeak(age_diff_peaks_down, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
Age_ATAC_anno5_diff <- annotatePeak(age_diff_peaks_up, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)

pdf("CHIPSeeker_Anno_Age_Sample_Diff_Peaks.pdf")
plotAnnoPie(Age_ATAC_anno1_diff)
plotAnnoPie(Age_ATAC_anno5_diff)
# plotAnnoBar(Age_ATAC_anno1_diff, title = "1wk Enriched Peaks")
# plotAnnoBar(Age_ATAC_anno5_diff, title = "5wk Enriched Peaks")
# vennpie(Age_ATAC_anno1_diff)
# vennpie(Age_ATAC_anno5_diff)
# upsetplot(Age_ATAC_anno1_diff) 
# upsetplot(Age_ATAC_anno5_diff) 
dev.off()

lamb_diff_peaks <-dba.report(LamB_ATAC, 
                            th = 0.1,
                            fold = 0.32,
                            bCounts= TRUE, 
                            bCalledDetail = TRUE,
                            contrast = 1,
                            method = DBA_EDGER)

lamb_diff_peaks_up <- lamb_diff_peaks[lamb_diff_peaks$Fold <0]
lamb_diff_peaks_down <- lamb_diff_peaks[lamb_diff_peaks$Fold >0]
tagMatLamBWT_diff <- getTagMatrix(lamb_diff_peaks_down, windows = promoter)
tagMatLamBKO_diff <- getTagMatrix(lamb_diff_peaks_up, windows = promoter)
LamB_ATAC_annoWT_diff <- annotatePeak(lamb_diff_peaks_down, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
LamB_ATAC_annoKO_diff <- annotatePeak(lamb_diff_peaks_up, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)


pdf("CHIPSeeker_Anno_LamB_Sample_Diff_Peaks.pdf")
plotAnnoPie(LamB_ATAC_annoWT_diff)
plotAnnoPie(LamB_ATAC_annoKO_diff)
dev.off()

lamc_diff_peaks <-dba.report(LamC_ATAC, 
                             th = 0.1,
                             fold = 0.32,
                             bCounts= TRUE, 
                             bCalledDetail = TRUE,
                             contrast = 1,
                             method = DBA_EDGER)

lamc_diff_peaks_up <- lamc_diff_peaks[lamc_diff_peaks$Fold <0]
lamc_diff_peaks_down <- lamc_diff_peaks[lamc_diff_peaks$Fold >0]
tagMatLamCWT_diff <- getTagMatrix(lamc_diff_peaks_down, windows = promoter)
tagMatLamCKO_diff <- getTagMatrix(lamc_diff_peaks_up, windows = promoter)
LamC_ATAC_annoWT_diff <- annotatePeak(lamc_diff_peaks_down, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
LamC_ATAC_annoKO_diff <- annotatePeak(lamc_diff_peaks_up, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)

pdf("CHIPSeeker_Anno_LamC_Sample_Diff_Peaks.pdf")
plotAnnoPie(LamC_ATAC_annoWT_diff)
plotAnnoPie(LamC_ATAC_annoKO_diff)
dev.off()

#Generate ATAC-seq Heatmaps and Average Profiles of differentially expressed peaks (Histograms for Heatmap)
#get tag matrix 
tagMat1_diff <- getTagMatrix(age_diff_peaks_down, windows = promoter)
tagMat5_diff <- getTagMatrix(age_diff_peaks_up, windows = promoter)
Diff_peaks_tagMatrixList <- list(One_Week=tagMat1_diff,Five_Week=tagMat5_diff, LamBWT =tagMatLamBWT_diff , 
                                 LamBKO =tagMatLamBKO_diff , LamCWT = tagMatLamCWT_diff,LamCKO = tagMatLamCKO_diff)




pdf("TagHeatmaps_Diff.pdf")
tagHeatmap(Diff_peaks_tagMatrixList, xlim=c(-3000,3000), color=NULL)
plotAvgProf(Diff_peaks_tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row")
dev.off()

# #repeat this for meta analysis genes
# GR_down_w_age_lamC <- makeGRangesFromDataFrame(EdgeR_Meta_Dwn_w_age_LamCKO,
#                          keep.extra.columns=FALSE,
#                          ignore.strand= TRUE,
#                          seqinfo=NULL,
#                          seqnames.field=c("seqnames", "seqname",
#                                           "chromosome", "chrom",
#                                           "Chr.x", "chromosome_name",
#                                           "seqid"),
#                          start.field="Start.x",
#                          end.field=c("End.x"),
#                          starts.in.df.are.0based=FALSE)
# GR_up_w_age_lamC <- makeGRangesFromDataFrame(EdgeR_Meta_Up_w_age_LamCKO,
#                                                keep.extra.columns=FALSE,
#                                                ignore.strand= TRUE,
#                                                seqinfo=NULL,
#                                                seqnames.field=c("seqnames", "seqname",
#                                                                 "chromosome", "chrom",
#                                                                 "Chr.x", "chromosome_name",
#                                                                 "seqid"),
#                                                start.field="Start.x",
#                                                end.field=c("End.x"),
#                                                starts.in.df.are.0based=FALSE)
# GR_down_w_age_lamC_anno <-annotatePeak(GR_down_w_age_lamC, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
# GR_up_w_age_lamC_anno <-annotatePeak(GR_up_w_age_lamC, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
# pdf("Anno_Down_w_Age_LamC_then_Up.pdf", width = 10, height = 15)
# plotAnnoPie(GR_down_w_age_lamC_anno)
# plotAnnoPie(GR_up_w_age_lamC_anno)
# plotAnnoBar(GR_down_w_age_lamC_anno, title = "Down with Age and LamC KO")
# plotAnnoBar(GR_up_w_age_lamC_anno, title = "Up with Age and LamC KO")
# vennpie(GR_down_w_age_lamC_anno)
# vennpie(GR_up_w_age_lamC_anno)
# upsetplot(GR_down_w_age_lamC_anno) 
# upsetplot(GR_up_w_age_lamC_anno) 
# dev.off()


# CHIPSeeker KEGG and GO - NOT WORKING YET, Just use Metascape Instead, better and easier --------------------------------------------------
#  0. make a list of matrices? GRanges?
# 1. Get filtered dataframe into GRanges object
# 2. Use dba.peakset(GRanges, bRetrieve=T) to 
#


# Need to find a way to make a DBA object from a dataframe
UPGranges <- makeGRangesFromDataFrame(DEGSUP, keep.extra.columns = T)
DOWNGranges <- makeGRangesFromDataFrame(DEGSDOWN, keep.extra.columns = T)

# Annotate the filtered sets individually
UPGranges_CSAnno <- annotatePeak(UPGranges, TxDb=txdb, tssRegion=c(-2000, 2000), verbose=FALSE)
DOWNGranges_CSAnno <- annotatePeak(DOWNGranges, TxDb=txdb, tssRegion=c(-2000, 2000), verbose=FALSE)

# Combine into a list 
list.annotations <- list(RR = UPGranges_CSAnno, RKO = DOWNGranges_CSAnno)


Diff_peaks_tagMatrixList
EdgeR_Meta_Dwn_w_age_LamCKO 
EdgeR_Meta_Up_w_age_LamCKO

## Functional Profiles Comparison 
genes = lapply(Diff_peaks_tagMatrixList, function(i) as.data.frame(i)$geneId)
names(genes) = sub("_", "\n", names(genes))

# The default for clustercompare will run a second set of statistical analysis
# - we will try another command, groupGO and see if we can avoid that
# - alternatively, retry clustercompare with p value of 1

testGO <- enrichGO(EdgeR_Meta_Dwn_w_age_LamCKO$Entrez.ID.x, OrgDb = org.Dm.eg.db, ont = "MF", readable =T)
pdf("test.pdf")
plotGOgraph(testGO)
#goplot(testGO)
upsetplot(testGO)
#barplot(testGO)
#cnetplot(testGO)
#dotplot(testGO)
#gseaplot(testGO)
#emapplot(testGO)
dev.off()

plotGOgraph(testGO, showCategory = 10, #How many top KEGG terms you want to display
            title = "KEGG Pathway Enrichment Analysis")


# Let's try reactome for enriched pathways - These are pretty generic, not helpful
# This approach will give you enrichment of genes defined by nearest feature to the peaks
pathway_Down_age <- enrichPathway(EdgeR_Meta_Dwn_w_age_LamBKO$Entrez.ID.x,organism = "fly")
#pathway_Primary <- enrichPathway(genes$Primary)
pdf("CHIPSeeker_Pathways.pdf", width = 10, height =15)
dotplot(pathway_Down_age)
#dotplot(pathway_Primary)
dev.off()

# This approach willl give enrichment defined as allowing many-to-many mappings - 
#  this function considers host gene (exon/intron), promoter and flanking gene from 
# intergentic region that may undergo control via cis-regulation aka Promoter/Enhancers/Silencers
# Also from CHIPseeker package

# First we will annotate many2many using seq2gene
up.m2m <- seq2gene(UPGranges, tssRegion = c(-2000,2000), flankDistance = 2000, TxDb = txdb)
down.m2m <- seq2gene(DOWNGranges, tssRegion = c(-2000,2000), flankDistance = 2000, TxDb = txdb)

# Then we will find enriched Reactome pathways given the many to many mapping
pathway.m2m.up <- enrichPathway(up.m2m)
pathway.m2m.down <- enrichPathway(down.m2m)

pdf("CHIPSeeker_m2m_Reactome.pdf")
dotplot(pathway.m2m.up, title = "RR Pathways")
dotplot(pathway.m2m.down, title = "RKO Pathways")
dev.off()




#download_KEGG(hsa) if you haven't already 
compKEGG <- compareCluster(geneCluster   = genes, #select which samples to use in the genes argument, format as [(1:2)]
                           fun           = "enrichKEGG",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH")

# Generate a dot plot for KEGG Terms 
pdf("CHIPSeeker_Dot_Plot_FunctionKEGG_All.pdf", width =10, height =15)
dotplot(compKEGG,
        showCategory = 10, #How many top KEGG terms you want to display
        title = "KEGG Pathway Enrichment Analysis")
dev.off()

#Plot counts in each KEGG term on X-axis for just one sample
deg_derived = genes$Derived
deg_primary = genes$Primary
KE_Derived_CF = enrichKEGG(deg_derived, pvalueCutoff = 0.05, minGSSize = 20, maxGSSize = 100, pAdjustMethod = "BH")
KE_Primary_CF = enrichKEGG(deg_primary, pvalueCutoff = 0.05, minGSSize = 20, maxGSSize = 100, pAdjustMethod = "BH")

#Generate Count X KEGG Plots
# Here count is the # of annotated genes, ratio is the count divided by all genes user provided in dataset (not in KEGG term)
pdf("CHIPSeeker_Dot_Plot_FunctionKEGG_byCounts.pdf", width =10, height =15)
dotplot(KE_Derived_CF, x="Count",showCategory=20, 
        title = "Derived CF Specific KEGG Terms")
dotplot(KE_Primary_CF, x="Count",showCategory=20, 
        title = "Primary Specific KEGG Terms")
dev.off()

# You can also do this for GO terms, GE stands for GENE ENRICHMENT
GE_Derived_CF = enrichGO(deg_derived,OrgDb = "org.Hs.eg.db", pvalueCutoff = 0.05, pAdjustMethod = "BH", readable =TRUE)
GE_Primary_CF = enrichGO(deg_primary,OrgDb = "org.Hs.eg.db", pvalueCutoff = 0.05, pAdjustMethod = "BH", readable =TRUE)
compGO <- compareCluster(geneCluster   = genes[(1:2)], #select which samples to use in the genes argument 
                         fun           = "enrichGO",
                         pvalueCutoff  = 0.05,
                         OrgDb = "org.Hs.eg.db",
                         pAdjustMethod = "BH")

# Create Dot Plots for GO Terms
pdf("CHIPSeeker_Dot_Plot_GO_15_Only_Counts.pdf", width =10, height =15)
dotplot(compGO,showCategory=20, 
        title = "Derived CF Specific GO Terms")
dotplot(GE_Derived_CF, x="Count",showCategory=20, 
        title = "Derived CF Specific GO Terms")
dotplot(GE_Primary_CF, x="Count",showCategory=20, 
        title = "Primary CF Specific GO Terms")
dev.off()

#You can extract the list of genes from a specific GO term:
# Data structure goes: enrichResult@geneSets$GO_ID
GO_results_Derived_CF <- GE_Derived_CF@result
GO_results_Derived_CF <- GO_results_Derived_CF[order(GO_results_Derived_CF$p.adjust),]
Cadherin_Derived_CF <- as.array(GO_results_Derived_CF[1,8])  # Look in the previous object to see what term to get based on p values
Cadherin_Derived_CF <- as.data.frame(unlist(strsplit(Cadherin_Derived_CF,"/")))
names(Cadherin_Derived_CF) <- c("Cadherin_Derived_GO")
head(Cadherin_Derived_CF)


#### Come back to this !!!
# We will try to visualize the GSEA result
library(ReactomePA)
deg15 = genes$`Unique
to_15_CF`
x <- enrichPathway(gene=deg15,pvalueCutoff=0.05, readable=T)
barplot(x, showCategory=8)
dotplot(x, showCategory=15)
emapplot(x)

deg19 = genes$`Unique
to_19_CF`
y <- enrichPathway(gene=deg19,pvalueCutoff=0.05, readable=T)
barplot(y, showCategory=8)
dotplot(y, showCategory=15)
emapplot(y)
cnetplot(y, categorySize = "pvalue", foldChange = genes$`19
         CF`)
# This makes me think that the genome math does not really help to pick DEGS/DARs - no FC shown


