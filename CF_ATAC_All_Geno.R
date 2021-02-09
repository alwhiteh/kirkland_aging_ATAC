## Analysis for Alex Whitehead's iPSC-Derived Cardiac Fibroblasts
# contact: alwhiteh@ucsd.edu
# CHIPSeeker Tutorial adapted from: https://www.bioconductor.org/packages/release/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html
setwd("/Users/alexanderwhitehead/Documents/PhD/RNASeq_stuff/CF_ATAC/All_Geno")
#renv::init()# initialize renv package for version control
#renv::restore() #use this to load version-specific packages
renv::snapshot()
# First we have taken the bed files that have been filtered for mitochondial reads, PCR dupicliates, and compensated for multimapped reads
# then MACS2 was used to call peaks for each bam file
# Next, bedtools slop was used to extend each feature by 1kb
# The resulting bed files were used in this analysis


## Load Libaries
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(clusterProfiler)

## We are going to take the extended summits from MACS2 + bedtools slop by 250 base pairs
library(DiffBind)
#samples <- read.csv("/Users/alexanderwhitehead/Documents/PhD/RNASeq_stuff/CF_ATAC/CF_ATAC_Sample_Sheet_Ext_Sum250.csv", stringsAsFactors = FALSE)
samples <- read.csv("/Users/alexanderwhitehead/Documents/PhD/RNASeq_stuff/CF_ATAC/EpiC_Geno/CF_ATAC_Sample_Sheet_Geno.csv", stringsAsFactors = FALSE)
CF_ATAC <- dba(sampleSheet=samples)

#Generate a correlation heatmap using peak caller scores
pdf("Peak_Caller_scores_Correlations.pdf")
plot(CF_ATAC)
dev.off()

#Generate counts for each sample
CF_ATAC <- dba.count(CF_ATAC)
CF_ATAC # should show Fractions of Reads in Peaks, this is from .32-.54 for the all samples analaysis

# Here is using windows generated from the default .xls output from MACS2
#samples2 <- read.csv("/Users/alexanderwhitehead/Documents/PhD/RNASeq_stuff/CF_ATAC/CF_ATAC_Sample_Sheet.csv", stringsAsFactors = FALSE)
#CF_ATAC2 <- dba(sampleSheet = samples2)
#CF_ATAC2 <- dba.count(CF_ATAC2)
#CF_ATAC2 # Here FRiP is from .19 to .46


# This time will plot correlation heatmap based on affinity scores (aka Read Density aka reads/length)
pdf("Read_densities_Correlations.pdf")
plot(CF_ATAC)
dev.off()

#estblish a contrast 
#CF_ATAC <- dba.contrast(CF_ATAC, categories = DBA_CONDITION,block = c(1:4, 9:16), minMembers = 2) #tell diffbind what criteria to contrast
CF_ATAC <- dba.contrast(CF_ATAC, group1 = c(1,2), group2 = c(3,4), minMembers = 2) #contrast 1 = CFs
CF_ATAC <- dba.contrast(CF_ATAC, group1 = c(5,6), group2 = c(7,8), minMembers = 2) #contrast 2 = EpiC
CF_ATAC <- dba.contrast(CF_ATAC, group1 = c(9,10), group2 = c(11,12), minMembers = 2) #contrast 3 = CPC
CF_ATAC <- dba.contrast(CF_ATAC, group1 = c(13,14), group2 = c(15,16), minMembers = 2) #contrast 4 = iPSC
CF_ATAC # View all of your contrasts

# Here's how to clear a contrast if you mess up: CF_ATAC$constrasts= NULL

# Perfom the differential anlaysis
CF_ATAC <- dba.analyze(CF_ATAC, method = DBA_EDGER) # Remove the method option when not trying to beat the no replicates rule
CF_ATAC # the last column shows what # of binding sites are different based on FDR = 0.05

# Create dummy variable to for indexing the number of contrasts
num_cont = c(1:length(CF_ATAC$contrasts))
name_cont = c("CFs","EpiCs","CPCs","iPSCs")

#Make mean average and volcano plots
for(i in num_cont) {
        pdf(paste("MA_and_Volcano_",name_cont[i],".pdf",sep=""))
        dba.plotMA(CF_ATAC, th = 0.01, contrast = i, method = DBA_EDGER) # generates MA plot, th = FDR cutoff
        dba.plotVolcano(CF_ATAC, th = 0.01, method = DBA_EDGER) # generates volcano plot
        plot(CF_ATAC, contrast =1, method = DBA_EDGER) # this generates a heatmap based on the differential genes from criterion defined in contrast 1, kind of recursive imho
        dev.off()
}

#Make some more plots
for(i in num_cont) {
        pdf(paste("ATAC_Venn_PCA_",name_cont[i],".pdf",sep=""))
        dba.plotVenn(CF_ATAC, mask =1:4)
        dba.plotPCA(CF_ATAC, label = DBA_ID)
        dba.plotBox(CF_ATAC, method = DBA_EDGER)
        dev.off()
}

#Create Heapmap based on affinities for differentially bound sites
for(i in num_cont) {
        pdf(paste("ATAC_Heatmap_No_Anno_",name_cont[i],".pdf",sep=""))
        dba.plotHeatmap(CF_ATAC, correlations = FALSE, 
                score= DBA_SCORE_TMM_MINUS_FULL_CPM)
        dev.off()
}
        
#Here is how to make a 3D plot that you can view through XQuartz (must be installed in command line)
dba.plotPCA(CF_ATAC, label=DBA_CONDITION, attributes=DBA_ID, b3D=TRUE, component=1:3)

# Retrieving the differentially bound sites aka sort from the entire dataset to select peaks based on statistical results and save to GRanges object
for(i in num_cont) {
        db_nam <- paste("ATAC_Report",name_cont[i],sep="_")
        assign(db_nam, dba.report(CF_ATAC, 
                         th = 0.01,
                         contrast = i,
                         method = DBA_EDGER)) #th is the FDR value that you would like to use to filter, you can also do fold change with "fold"
}

#ATAC_Report_CFs #view your GRanges object

# quantitatively see how many DEGS up and down
#sum(CF_ATAC.DB$Fold<0) # Down - here we have 232, # this is lower limit of motif finding input
#sum(CF_ATAC.DB$Fold>0) # Up - here we have 2412

#see how many peaks overlap between # of datasets
#pdf("Histogram_Overlaps.pdf")
#plot(dba.overlap(CF_ATAC, mode=DBA_OLAP_RATE), type='b', ylab=" # of peaks", xlab = "Overlap at least this many peaksets")
#dev.off()

# Annotate CHIP Peaks, this package is very selective with annotations, later we will use the CHIPSeeker package to get better annotations
library(ChIPpeakAnno)
data("TSS.human.GRCh38") #load TSS annotation for Hg 38 from BiomaRt
for(i in num_cont) {
        assign(paste("all.data.peaks",name_cont[i],sep="_"), dba.report(CF_ATAC, th = 1, contrast = i, method = DBA_EDGER)) # change th into 1 if you want everything, not just differential genes
}

#all.data.peaks 

#Annotate Peaks with information on closest TSS 
anno_manager <- as.character(rep(NA,length(num_cont))) # create vector of variable names to set up lappy for annotation
anno_manager <- mget(ls(pattern = "all.data.peaks"))
anno_manager <- sapply(anno_manager, annotatePeakInBatch, AnnotationData = TSS.human.GRCh38)

# Now we need to add GeneIDs
library(org.Hs.eg.db)

anno_manager <- sapply(anno_manager,addGeneIDs, orgAnn = "org.Hs.eg.db", IDs2Add = c("symbol","genename"))
for(i in num_cont) {
        write.table(anno_manager[i], file = paste("Genotype", name_cont[i],"ATAC.txt", sep ="_"), sep = "\t", row.names=F)
}

#Change from a GRanges object into a dataframe and create BED files for motif finding
deseq_colnames <- c("seqnames","start","end","width","strand","Conc","Conc_group1","Conc_group2",
                    "Fold","p.value","FDR","peak","feature","start_position","end_position",
                    "feature_strand","insideFeature","distancetoFeature","shortestDistance","fromOverlappingOrNearest","symbol","genename")

for(i in num_cont) {
        out_deseq <- as.data.frame(anno_manager[i])
        colnames(out_deseq) <- deseq_colnames
        out_deseq <- out_deseq[c("seqnames","start","end","peak","FDR","feature_strand","Fold")] #need to be in this order
        DEGSUP <- out_deseq[which(out_deseq$Fold>0 & out_deseq$FDR<0.01),] 
        DEGSDOWN <- out_deseq[which(out_deseq$Fold<0 & out_deseq$FDR<0.01),] 
        DEGS <- out_deseq[which(out_deseq$FDR<0.01),] 
        #DEGSUP <- DEGSUP[1:6]
        #DEGSDOWN <- DEGSDOWN[1:6]
        #DEGS <- DEGS[1:6]
        #DEGSUP[,5] <- trunc(DEGSUP[,5]+1) 
        #DEGSDOWN[,5] <- trunc(DEGSDOWN[,5]+1)
        #DEGS[,5] <- trunc(DEGS[,5]+1)
        write.table(DEGSUP, file = paste("RR",name_cont[i],"DARS.bed", sep="_"), sep = "\t", row.names=F, col.names = F, qmethod = "double")
        write.table(DEGSDOWN, file = paste("RKO",name_cont[i],"DARS.bed", sep="_"), sep = "\t", row.names=F, col.names =F,qmethod = "double")
        write.table(DEGS, file = paste("Geno",name_cont[i],"DARS.bed",sep="_"), sep = "\t", row.names=F, col.names = F,qmethod = "double")
}

# For whatever reason, you need to open the files in excel and delete the last column manually,
        # then change column 5 to all ones, and format as "number" in excel (should look like 1.00)



# After this, move the bed files into your ATAC seq directory for command line operations and run HOMER findMotifsGenome.pl

# Now we will combine all of the motifs for each sample into a matrix and plot to compare
# Here we are using the known results
library(tidyverse) # using tidyverse helps us read in data with repeat row names and better format column names - formate as tibble
RR_CF_motif <- read_tsv('/Users/alexanderwhitehead/Documents/PhD/RNASeq_stuff/CF_ATAC/All_Geno/RR_CFs_DARS.bedmotifs/knownResults.txt')
RR_CF_top20 <-RR_CF_motif[c(1:20),c(1,4)]
RR_EpiC_motif <- read_tsv('/Users/alexanderwhitehead/Documents/PhD/RNASeq_stuff/CF_ATAC/All_Geno/RR_EpiCs_DARS.bedmotifs/knownResults.txt')


# Make Annotated Heatmap --------------------------------------------------


# Take these DARS and make a pheatmap
library(gplots)
library(pheatmap)
library(plyr)
library(tibble)
library(tidyr)
#mat <- as_tibble(read.table("Derived_vs_Primary_CFs_ATAC.txt", sep ="\t", header = TRUE ))
mat <- as_tibble(read.table("Genotype_CFs_ATAC.txt", sep ="\t", header = TRUE ))


# you cannot move factors easily between objects, so we will convert them into characters
mat$symbol <- as.character(mat$symbol)
mat$feature <- as.character(mat$feature)
mat$symbol[is.na(mat$symbol)] <- mat$feature[is.na(mat$symbol)] #replace any NAs with ENSEMBL ID

# need to wrangle the rows with the same annotations and combine counts
#mat <- as.data.frame(mat)

#Generate a heatmap of peak region densities for top 500 differentially accesible regions
pdf("Heatmap_annotate_top500_scale.pdf", width = 10, height =15)
pheatmap((mat[c(1:500),c(7:8)]),   #the matrix is already sorted by FDR from smallest to largest
         show_rownames = T,
         cluster_rows = T,
         cluster_cols = T,
         labels_row = mat$symbol,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         border_color = F,
         cex = 0.4,
         scale = "row",
         fontsize_row = 8,
         legend = T) 
# you can also use scale = row to normalize each row to median and plot std dev
dev.off()

pdf("Heatmap_annotate_top500.pdf", width = 10, height =15)
pheatmap((mat[c(1:500),c(7:8)]),   #the matrix is already sorted by FDR from smallest to largest
         show_rownames = T,
         cluster_rows = T,
         cluster_cols = T,
         labels_row = mat$symbol,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         border_color = F,
         cex = 0.4,
         fontsize_row = 8,
         legend = T) 
dev.off()

# Functional Analysis Using CHIPSeeker ------------------------------------


#Now take these DEGS and send to KEGG + GO analysis
# We will want to use the GRanges structure for CHIPSeeker - here we called it data.peaksAnno
library(ChIPseeker)
library(ReactomePA)
library(clusterProfiler)
library(biomaRt)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#get the TSS regions for your genome
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream =3000)

#see the structure of your DBA object from earlier - this will help you choose peaksets for CHIPSeeker
dba.show(CF_ATAC)

#Retrieve all of the peaksets from your previous list object and store in a new variable
peaks.consensus <- dba.peakset(CF_ATAC, bRetrieve = T) #bRetrieve=true just means that you are extracting peaksets instead of adding

# break each condition out into individual variables - you can also do this with just one sample, just remove the & {sample 2}
peaks.derived <- peaks.consensus[CF_ATAC$called[,1]==1 & CF_ATAC$called[,2]==1]
peaks.primary <- peaks.consensus[CF_ATAC$called[,3]==1 & CF_ATAC$called[,4]==1]

# calculate tagMatrix
tagMatrix.der <- getTagMatrix(peaks.derived, windows = promoter)
tagMatrix.pri <- getTagMatrix(peaks.primary, windows = promoter)

# Combine the tagMatrices into a list for downstream analysis
tagMatrixList <- list(RR=tagMatrix.der,RKO=tagMatrix.pri)

#Plot tagMatrix heatmaps for each condition, not really helpful but here it is
pdf("CHIPSeeker_ATAC_Heatmap.pdf")
tagHeatmap(tagMatrixList, xlim=c(-3000,3000), color=NULL)
dev.off()

# Generic CHIPSeeker QC Stats ---------------------------------------------


#Plot the Average Profiles of ATAC Peaks on top of each other
pdf("CHIPSeeker_Average_Profiles_Stacked.pdf")
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))
dev.off()

#Plot the Average Profiles of ATAC Peaks Seperately, the resample value chooses how many peaks to sample
pdf("CHIPSeeker_Average_Peak_Profiles_by_Sample.pdf")
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row")
dev.off()

# Create Peak Heatmaps for each sample
pdf("Average_Peak_Heatmaps_by_Sample.pdf", width = 20, height = 20)
tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL)
dev.off()

# Annotate the peaks, here I added "Anno2" as a suffix to represent that these were annotated by CHIPSeeker, not Diffbind/CHIPPeakAnno
peaks.derived_CSAnno <- annotatePeak(peaks.derived, TxDb=txdb, tssRegion=c(-2000, 2000), verbose=FALSE)
peaks.primary_CSAnno <- annotatePeak(peaks.primary, TxDb=txdb, tssRegion=c(-2000, 2000), verbose=FALSE)

# create a list of annotations
list.annotations <- list(RR = peaks.derived_CSAnno, RKO = peaks.primary_CSAnno)

# This can be used to visualize samples side by side - doesn't work for vennpie or upsetplot
pdf("CHIPSeeker_Anno_Bar_TSS.pdf")
plotAnnoBar(list.annotations)
plotDistToTSS(list.annotations)
dev.off()

# create barplots, vennpie, upsetplots, and distance to TSS for the peaks of each group
pdf("CHIPSeeker_Anno_plots_by_Sample.pdf")
plotAnnoBar(peaks.derived_CSAnno, title = "RR Cells Feature Dist")
vennpie(peaks.derived_CSAnno)
upsetplot(peaks.derived_CSAnno) #this shows overlapping annotations
plotAnnoBar(peaks.primary_CSAnno, title = "RKO Cells Feature Dist")
vennpie(peaks.primary_CSAnno)
upsetplot(peaks.primary_CSAnno)
dev.off()

# CHIPSeeker KEGG and GO --------------------------------------------------
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




## Functional Profiles Comparison 
genes = lapply(list.annotations, function(i) as.data.frame(i)$geneId)
names(genes) = sub("_", "\n", names(genes))

# The default for clustercompare will run a second set of statistical analysis
        # - we will try another command, groupGO and see if we can avoid that
        # - alternatively, retry clustercompare with p value of 1

testGO <- enrichGO(genes$Derived, OrgDb = org.Hs.eg.db, ont = "cc", readable =T)
plotGOgraph(testGO, showCategory = 10, #How many top KEGG terms you want to display
        title = "KEGG Pathway Enrichment Analysis")





# Let's try reactome for enriched pathways - These are pretty generic, not helpful
# This approach will give you enrichment of genes defined by nearest feature to the peaks
pathway_Derived <- enrichPathway(genes$Derived)
pathway_Primary <- enrichPathway(genes$Primary)
pdf("CHIPSeeker_Pathways.pdf", width = 10, height =15)
dotplot(pathway_Derived)
dotplot(pathway_Primary)
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


