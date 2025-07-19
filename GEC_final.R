## Create the options object, initial starting code and other code throughout taken from pipeline
opt <- list() 
opt$project    <- "TCGA-LUSC" #changed project to TGCA-LUSC

## Setup package and GDC directory options
opt$wd <- path.expand(file.path("~","bio321g","DiffExp")) #file path changed to mine
opt$pckDir  <- c(
  file.path(opt$wd,"R_gdc"),              
  "/stor/scratch/Bio321G_NB_Spring2024/R" 
)
opt$gdcPath <- "/stor/scratch/Bio321G_NB_Spring2024/GDCdata" 
opt$biocPackages <- c(
  "ggplot2","ggrepel","TCGAbiolinks","DESeq2"
) 

##### Setup location to store packages #####
## Setup the working directory
if(!dir.exists(opt$wd)){dir.create(opt$wd,recursive = T,showWarnings = T)}
setwd(opt$wd)

opt$backupPath    <- file.path(opt$wd,"R_gdc")
opt$backupGdcPath <- file.path(opt$wd,"GDCdata")

##### Install packages #####

if(!dir.exists(opt$pckDir[1])){
  message("Creating home package directory:\n",opt$backupPath)
  dir.create(opt$backupPath,recursive = T,showWarnings = T)
}
if(!all(dir.exists(opt$pckDir))){
  message("Changing package directory to:\n",opt$backupPath)
  opt$pckDir <- opt$backupPath
}

# Add paths to the list of locations
.libPaths(opt$pckDir)
message(paste0(.libPaths(),sep = "\n"))

## checking for success
if(!all(opt$biocPackages%in%installed.packages())){
  if(!"BiocManager"%in%installed.packages()){
    install.packages("BiocManager", lib = opt$pckDir)
  }
  
  ## Update the installed packages
  update.packages(instlib = opt$pckDir,ask = F)
  
  BiocManager::install(
    opt$biocPackages, lib = opt$pckDir,
    ask = F,update = T
  )
}

# Check and load packages 
if(!all(opt$pckDir%in%.libPaths())){
  print("Setting lib paths...")
  .libPaths(opt$pckDir)
}
library(ggplot2)
library(ggrepel)
library(TCGAbiolinks) 
library(DESeq2) 


# Search for TCGA Files 
query1 <- GDCquery(
  project = opt$project,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  access = "open"
)

# Extract sample IDs from initial query 
samDf <- query1$results[[1]]
table(samDf$sample_type) 


#Subsetting and filtering data

groupHigh <- read.csv("~/bio321g/DiffExp/R_gdc/high mutated lrpb1/lusc LRP1B HIGH mpact mutation - gdc_sample_sheet.2024-03-26_HIGH_LRP1B.csv",fill=T,header = 1) #edit to change to my gdc files and file locations
groupNon <- read.csv("R_gdc/non mutated lrpb1/gdc_sample_sheet.2024-03-26_NONMUTATED_LRPB1 - gdc_sample_sheet.2024-03-26_NONMUTATED.csv", fill=T,header = 1)

groupHigh <- subset(groupHigh, Sample.Type != "Solid Tissue Normal") #removing solid tissue normal, only wanted primary tumor
groupNon <- subset(groupNon, Sample.Type != "Solid Tissue Normal")

groupHigh_sampleIDs <- groupHigh$Sample.ID #assigning sample.id to variable
groupNon_sampleIDs <- groupNon$Sample.ID


subset_Df <- samDf[samDf$sample.submitter_id %in% groupHigh_sampleIDs | samDf$sample.submitter_id %in% groupNon_sampleIDs, ] #subsetting to contain only the submitter ids that match high and non mutation

groupHigh_df <- samDf[samDf$sample.submitter_id %in% groupHigh_sampleIDs, ] #splitting the subset up into high and non df
groupNon_df <- samDf[samDf$sample.submitter_id %in% groupNon_sampleIDs, ]

set.seed(6481) #given seed
if(nrow(groupNon_df)>50){ #had to randomly select 50 rows from the Non df, this was the max allowed value
  groupNon_df <- groupNon_df[sample(nrow(groupNon_df), size = 50),]
}

samDf2 <- rbind(groupHigh_df, groupNon_df) #combining the two data sets together based on rows

samDf2$group <- ifelse(1:nrow(samDf2) <= 6, "sampled_groupNon", "groupHigh_df") #ensuring the first 6 rows of Sam Df2 are the more wild type sample (the Non mutated group), the rest/bottom are the high impact mutated
samDf2$analysis_workflow_type <- "STAR - Counts" #new column called analysis workflow type, filled using star counts
samDf2$sample_type <- "Primary Tumor" #new column named sample type, adding primary tumor 
write.csv(samDf2, file = "samDf2.csv", row.names = FALSE) #creates csv file



## checking if groups overlap
sum(groupHigh_df$Sample.ID%in%groupNon_df$Sample.ID) #Should be 0
sum(groupNon_df$Sample.ID%in%groupHigh_df$Sample.ID) #Should be 0

## checking if data has values
nrow(samDf2)             #Should be >0 
length(groupHigh$Sample.ID) #Should be >0
length(groupNon$Sample.ID) #Should be >0

## checking if subsetting has worked
sum(samDf2$sample.submitter_id%in%groupHigh$Sample.ID) #Should be >=6 and not >>>50
sum(samDf2$sample.submitter_id%in%groupNon$Sample.ID) #Should be >=6 and not >>>50
length(unique(samDf2$cases)) #This should be length of grp1 + grp2 and is used in next step


#### Download TCGA Files ####
query2 <- GDCquery(
  project = opt$project,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  access = "open",
  barcode = samDf2$cases
)


if(!file.exists(opt$gdcPath)){
  warning(
    "GDC download directory (",opt$gdcPath,") not found...\n",
    "Creating directory to store GDC data in the wd:\n",
    opt$backupGdcPath
  )
  opt$gdcPath <- opt$backupGdcPath
  dir.create(opt$gdcPath,recursive = T,showWarnings = T)
}


# Downloading data
GDCdownload(
  query = query2,
  method = "api",
  files.per.chunk = 10,
  directory = opt$gdcPath
)


# organizing multi dimensional data
dds <- GDCprepare(query = query2,directory = opt$gdcPath)

## Check if dds is the right class of data
class(dds) 

## Checking clinical data
dim(colData(dds))   
dim(as.data.frame(rowRanges(dds))) 

## Check sample data
assays(dds) 
dim(assays(dds)[[1]]) #Should be the number of human genes by the number of samples

## Checking for multiple observations / undesired tissues, etc
if(sum(duplicated(dds$sample_submitter_id))>0){
  warning("Some IDs present more than once. Consider subsetting if this was unintentional.")
}
if(length(unique(dds$sample_type))>1){
  warning("More than one type of sample tissue type present. Consider subsetting if this was unintentional.")
}


# Filtering the samples and loci analyzed
# Identify the groupings of samples 
# Create new "columns" in the sample data based on grouping
dds$groupHigh <- dds$sample_submitter_id%in%groupHigh$Sample.ID
dds$groupNon <- dds$sample_submitter_id%in%groupNon$Sample.ID


if(!all(dds$groupHigh==!dds$groupNon)){
  stop("Your groupings are not mutually exclusive")
}
  
dds$comp <- factor(dds$groupNon,levels = c("FALSE","TRUE")) #creates new column in dds, factor levels to determine order

# DESeq2 to create the filter data 
# Convert to DESeq object
dds <- DESeqDataSet(dds, design = ~ comp)

## Normalizeing the read depth values
dds <- estimateSizeFactors(dds)


# More Filtering 

rawCounts <- as.data.frame(counts(dds, normalized = F)) #turning dds into dataframe

counts_High <- rawCounts[,dds$groupHigh] #extracting grouphigh raw counts into counts high object, we had to split these raw counts into 2 different groups
counts_Non <- rawCounts[,dds$groupNon]

zcounts_High <- counts_High > 0 #seeing which raw counts are greater than 0 in high group
zcounts_Non <- counts_Non > 0

meanHigh <- rowMeans(zcounts_High) #finding mean of rows, will use these values to compare to 0.1 aka 10% we want to remove --> we want to keep the means that are 90% above zero
meanNon <- rowMeans(zcounts_Non)

mHighfilter <- meanHigh >= 0.1 #will provide true/false values for the mean row values greater or equal to 1 
mNonfilter <- meanNon >= 0.1
mBothfilter <- mHighfilter & mNonfilter #combining both to have a dds that has filters for both high and non

dd_t <- dds[which(mBothfilter),] #filtering out only the true values
dt <- rowMeans(as.data.frame(counts(dd_t, normalized = F))==0) #need to turn the dds back into numeric to make a histogram
hist(dt) #visualizing if the group above 0.9 is removed (should no longer show peak at 1.0)


normCounts <- as.data.frame(counts(dds, normalized = TRUE))
## Determining which loci to include in analysis

normCnt_var <- apply(X = normCounts,MARGIN = 1,FUN = var)
#Calculating  variance per row in a loop
normCnt_varCutoff <-max(c(0.01,median(normCnt_var)))      
normCnt_logic <- normCnt_var > normCnt_varCutoff          
normCounts_sub<- normCounts[which(normCnt_logic),]        


# Transpose to place samples in rows and loci in columns
transposeNorm <- t(normCounts_sub) #so swapping rows and columns?

## Calculating scaled PCA with prcomp
scaled_PCA <- prcomp(transposeNorm, scale. = TRUE)

pve <- scaled_PCA$sdev^2/sum(scaled_PCA$sdev^2) #proportion variance calculation using scaled_pca
xLab <- paste0("PC1 (",round(pve[1]*100,2),"%)") #variance pc1
yLab <- paste0("PC2 (",round(pve[2]*100,2),"%)") #variance pc2

visual <- data.frame(obsnames=rownames(scaled_PCA$x), scaled_PCA$x) #creating dataframe from scaled pca


#Creating PCA plot

install.packages("viridis") #installing package for color in plots


PCA_plot <- ggplot(
  data = scaled_PCA$x, #data used for ggplot is scaled pca
  mapping = aes(x = scaled_PCA$x[,1], y = scaled_PCA$x[,2], color= dds$comp, shape = dds$comp) #color and shape based off of dds$comp true false values, dds$comp label does appear on legend
)+
  geom_point()+
  theme_bw()+
  coord_fixed()+
  labs(x=xLab,y=yLab) #xlab and ylab used from variance calculations earlier, will show up on labels

print(PCA_plot) #prints pca plot into plots pane, quicker to visualize and see what needs to be fixed


#Differential expression analysis

#Convert to DESeq object
dds_seq <- DESeqDataSet(dd_t, design = ~ comp)

#recalculating normalization

ddsr <- estimateSizeFactors(dds) #used new name instead of dds just in case of errors

## Run the Deseq2 package analysis
ddseq <- DESeq(ddsr) #reassigned to new name
res <- results(ddseq) # Organizes the results


#Accumulating data into a data.frame 
#Add in gene data from the rowRanges section of the SummarizedExperiment object

colsAdded <- !colnames(as.data.frame(rowRanges(ddsr)))%in%colnames(as.data.frame(res))
resOutput <- cbind(
  as.data.frame(res),
  as.data.frame(rowRanges(ddsr))[,colsAdded]
)


#Volcano Plot

Significance <- resOutput$padj < 0.01 & resOutput$log2FoldChange >= 1.5 #subsetting where fdr adjusted p value is <0.01 and log2 fold change is >= 1.5

sigplot <- ggplot(
  data = resOutput, 
  mapping = aes(x= log2FoldChange, y = -log10(padj), color = baseMean)
)+
  geom_point(aes(shape = Significance))+ #put the shape specifically in geom point because of layer errors
  scale_color_viridis_c(trans = "log10")+
  theme_bw()+
  labs(x = "log2(foldchange)", y = "-log10(FDR-adjusted p-values)")+
  ggtitle("Differential Gene Expression Analysis of TGCA-LUSC Tumors
          with High Impact Mutations in LRP1B vs No Mutations in LRP1B, 
          p-Value = 0.05")+
  geom_hline(yintercept = -log10(0.05), color = "black", linetype = "dashed")+ #used p value of 0.05 as cutoff
  geom_vline(xintercept = c(-1.5), color = "blue", linetype = "dashed")+ #using 1.5 and -1.5 as cutoffs for FDR adjusted p value
  geom_vline(xintercept = c(1.5), color = "red", linetype = "dashed")+
  geom_text_repel(data = head(resOutput[order(resOutput$pvalue),],10), #layer that labels the significant genes taken from the top 10 rows of the ordered p values, different genes labeled if 10 changed to 20
                  aes(label = gene_name), #labeling genes by name
                  size = 3, box.padding = 3) #considering adding max.overlaps, but volcano plot looks messier, especially with the lines being used to label

print(sigplot) #prints volcano plot


#####  Code to visualize and summarize results #####


nmsig <- sum(complete.cases(resOutput$padj))
less0.01 <- sum(resOutput$padj < 0.01, na.rm = T) #error if na.rm isnt included, this code shows how many padj values have fdr adjusted less than 0.01, removed NA values
proportion <- (less0.01)/nmsig * 100 #calculating proportion, multiplied by 100 to turn into percentage
completel2fc <- complete.cases(resOutput$log2FoldChange) #had to remove missing or na values in l2fc 
l2fc <- sum(abs(resOutput$log2FoldChange[completel2fc]) > 1.5) #outputs the value of significant and absolute l2fc



