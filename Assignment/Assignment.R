
package_needed <- c("tidyverse", "factoextra", "ggpubr","edgeR","DESeq2")
install.packages(package_needed)
for (pkg in package_needed) {
  if(!require(pkg, character.only = T)) library(pkg, character.only = T)
} 


citation("edgeR")
citation("clusterProfiler")
citation("DESeq2")

count_matrix <- read_tsv("B-ALL.txt")# Please place in directory!

#1. Description of the dataset (number of samples, groups, other metadata)
dim(count_matrix)
print(paste("Number of genes: ", nrow(count_matrix)))
sample.names <- colnames(count_matrix[-1]) # Gets sample names from all samples(excluding gene)
split.column.names <- str_split_fixed(sample.names, "_", 2) # Splits Patient by sample type
tissue<-split.column.names[,2]
patient<-split.column.names[,1]
print(paste0("What cell types are these samples from? ", paste(unique(split.column.names[, 2]), collapse = " ")))# Outputs cell types 
print(paste0("How many patients are there samples from? ", length(unique(split.column.names[,1]))))

sample.meta.data<-data.frame(sample.names,tissue,patient)#Contains patient Number and tissue type 
gene.meta.data <- count_matrix[1]  # Contains the unfiltered Gene names


#Removal of genes with an occurance of 0 across all samples to allow for downstream signalling
summary(rowSums(count_matrix[-1]==0)==6)

filtered_data <- count_matrix[rowSums(count_matrix[-1])>0,]
#Ordering shows that 3 genes have much higher values than the rest of the dataset, removal of these genes aids in removal of highly varaible data
sort(rowSums(filtered_data[-1]),decreasing = TRUE)
filtered_data <- filtered_data[rowSums(filtered_data[-1])<400000,]





filtered_gene_names <- filtered_data[1]#Contains the filtered gene names

filtered_data <- filtered_data[-1]#Removes gene names from the matrix

row.names(filtered_data) <- filtered_gene_names$Gene#Assigns gene names to the columns in martix

#Creates a DGEList S4 class which takes the argument of a count martix, and allows division with sample meta data

dgelist <- DGEList(filtered_data,samples = patient,group = tissue,genes = filtered_gene_names)
#Why do I not assign cpm in to the DGEList? https://www.biostars.org/p/9475236/, https://www.biostars.org/p/491767/ -  Contains statements by one of the developers on why
#Most edgeR DE pipelines never modify the original counts in any way.
#Normalization for library size is instead implicit as part of the TMM model-fitting. 
#edgeR does not use cpm or rpkm values internally in its DE pipelines, rather they are only for export or for graphical purposes.
dgelist$samples$lib.size
#~~~~~~~~~~~~~~~~~~~~~~~Filtering~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Counts per million calculations 
cpm<-cpm(dgelist)
lcpm <- cpm(dgelist, log=TRUE)
#This function transforms a matrix into a long format containing gene, logCPM, and patients/tissue numbers 
#I use this function multipul times for checking filtering techniques 
long_transormation <- function(x,value){    #used this code alot so decided to make a function, allows df to be plotted easier 
  long<-as.data.frame(x) %>% 
    rownames_to_column("Gene") %>%
    pivot_longer(-Gene, values_to= value, names_to="Patients")
  return(long)
}

long_transormation(lcpm,"LogCPM") %>% ggplot()+geom_boxplot(aes(x = Patients, y = LogCPM, fill = Patients))#Variance!

index.filter<-filterByExpr(dgelist, group=tissue, min.count = 15, min.total.count = 20, min.prop = 0.7, design_matrix)#Modifed parameters to remove lowely expressed genes

dgelist.filtered<-dgelist[index.filter, , keep.lib.sizes=FALSE]
dgelist.filtered$samples$lib.size <- colSums(dgelist.filtered$counts)
lcpm.filtered<-cpm(dgelist.filtered, log=TRUE)

#Piping long transformed data into a density plot 
long_transormation(lcpm,"LogCPM") %>% ggplot() + geom_density(aes(LogCPM, colour=Patients))+
  ylab("Count Density")+ggtitle("Unfiltered LogCPM count density across patient samples")+
  theme_pubr()+
  theme(plot.title = element_text(hjust = 0.5))

long_transormation(lcpm.filtered, "LogCPM") %>% ggplot() + geom_density(aes(LogCPM, colour=Patients))+
  theme_pubr()+ylab("Count Density")+ggtitle("Customised Filtered LogCPM count density across patient samples")+
  theme(plot.title = element_text(hjust = 0.5))

plotMeanVar(dgelist.filtered)#views the mean varaince at the gene-level 

plotMeanVar(dgelist.filtered.norm)


#~~~~~~~~~~~~~~~~~~~~~~~Normalisation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Normalisation 
dgelist.filtered.norm<-calcNormFactors(dgelist.filtered, method="TMM")
dgelist.filtered.norm$samples$norm.factors
dgelist.filtered.norm$samples$lib.size <- (dgelist.filtered.norm$samples$lib.size) * (dgelist.filtered.norm$samples$norm.factors) # Re-calculate Lib size

#Normalisation doesn't make a massive difference since the data has already been normalised(but not for library size)

#Gene Names are removed from column after normalisation!
row.names(dgelist.filtered.norm$counts) <- dgelist.filtered.norm$genes[,1]
#check TMM effects on lib sizes 
unnorm_lib <- data.frame(
  Sample = colnames(dgelist.filtered),
  LibSize = dgelist.filtered$samples$lib.size,
  tissue = sample.meta.data$tissue)


norm_lib <- data.frame(
  Sample = colnames(dgelist.filtered.norm),
  LibSize = dgelist.filtered.norm$samples$lib.size,
  tissue = sample.meta.data$tissue
)

ggplot(unnorm_lib, aes(x = Sample, y = LibSize, color = tissue)) +
  geom_bar(stat = "identity", fill = "white",alpha = 0.7) +
  theme_pubr()+
  labs(title = "Barplot of Un-normalised Library Sizes",
       x = "Sample",
       y = "Library Size") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  
  
  ggplot(norm_lib, aes(x = Sample, y = LibSize, color = tissue)) +
  geom_bar(stat = "identity", fill = "white",alpha = 0.7) +
  theme_pubr()+
  labs(title = "Barplot of TMM normalised Library Sizes",
       x = "Sample",
       y = "Library Size") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#---------------------------------------------------------------------------------#

cpm.filtered.norm<-cpm(dgelist.filtered.norm) # Format the filtered.norm dgeobject with CPM and LogCPM for plotting
lcpm.filtered.norm<-cpm(dgelist.filtered.norm, log=TRUE)

long.cpm.filtered.norm <- long_transormation(cpm.filtered.norm,"CPM")

long.lcpm.filtered <- long_transormation(lcpm.filtered,"LogCPM")

long.lcpm.filtered.norm <- long_transormation(lcpm.filtered.norm,"LogCPM")

df.plotting<-full_join(long.cpm.filtered.norm, long.lcpm.filtered.norm)
df.plotting<-left_join(df.plotting, sample.meta.data, by=c("Patients"="sample.names"))
df.plotting<-left_join(df.plotting, gene.meta.data)


long.lcpm.filtered %>% left_join(sample.meta.data, by=c("Patients"="sample.names")) %>% left_join(gene.meta.data) %>% 
  ggplot()+geom_density(aes(LogCPM, colour=Patients))+
  ggtitle("LogCPM density across all patient samples before TMM normalisation")+
  theme_pubr()+
  theme(plot.title = element_text(hjust = 0.5))+

ggplot(df.plotting) + geom_density(aes(LogCPM, colour=Patients))+
  ggtitle("LogCPM density across all patient samples after TMM normalisation")+
  theme_pubr()+
  theme(plot.title = element_text(hjust = 0.5))

long.lcpm.filtered %>% left_join(sample.meta.data, by=c("Patients"="sample.names")) %>% left_join(gene.meta.data) %>% 
  ggplot()+geom_boxplot(aes(x=Patients, y = LogCPM, color = tissue))+
  ggtitle("Varaince before TMM normalisation")+
  theme_pubr()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
ggplot(df.plotting)+geom_boxplot(aes(x=Patients, y = LogCPM, color = tissue))+
  ggtitle("Varaince after TMM normalisation")+
  theme_pubr()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


ggplot(df.plotting %>% filter(Gene=="VEGFA")) + 
  geom_point(aes(x=tissue, y=CPM, color = patient)) +
  ggtitle("VEGFA expression in BM and CNS")+
  theme_pubr()+
  theme(plot.title = element_text(hjust = 0.5))




#3. PCA analysis (plot for dimensionality reduction and scree plot)

svd_PCA <-
  lcpm.filtered.norm %>% #Takes filtered and normalised logCPM
  t() %>% # Transposes martix, since PCA takes features on the column axis 
  prcomp(scale = T) # Genenerates a prcomp object, I scale the data using scale = T, This normalises the s.d to 1 and mean to 0  

summary(svd_PCA) # Shows S.D, variance and proportion of each Principle component

fviz_eig(svd_PCA, addlabels = T)+ # Genereates a scree plot showing how much variance is captured by each PC, shows that 99.9% of the variance is observed in PC1 & PC2
  theme_pubr()+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("Principle Components")+
  ylab("Percentage of variance captured")
  
#fviz_pca_ind plots a biplot of individals and creates a graph 
fviz_pca_ind(svd_PCA, 
             habillage = sample.meta.data$tissue,
             label = "none",
             invisible = "quali") + 
             ggtitle("Tissue's grouping")+
             theme_pubr()+
             theme(plot.title = element_text(hjust = 0.5))

fviz_pca_ind(svd_PCA, 
             habillage = sample.meta.data$patient,
             label = "none" ,
             invisible = "quali") + 
             ggtitle("Patient Grouping")+
             theme_pubr()+
             theme(plot.title = element_text(hjust = 0.5))
#PCA analysis shows that PC1 accoutns for 68.5% of the variance in the datasets and seperates the tissue type(+Ve for CNS and -ve for BM)
#PC2 accounts for 31.3% of the variance in the dataset and seperates varaince between patient 6 


rownames(svd_PCA$rotation) <- rownames(lcpm.filtered.norm)

fviz_pca_var(svd_PCA, select.var = list(cos2 = 15), repel = TRUE)+ # Shows eignen vector loadings for the top 20 gene in biplot, the majority of variation is in PC1 dimension
  theme_pubr()+
  ggtitle("Top 15 Eigenvector Gene Loading")+
  theme(plot.title = element_text(hjust = 0.5))

as.data.frame(svd_PCA$rotation[,1:2]) %>% # Top +VE contributing genes to PC1
  slice_max(n=10, order_by = PC1)

as.data.frame(svd_PCA$rotation[,1:2]) %>% #Top +ve contributing genes to PC2
  slice_max(n = 10, order_by = PC2)


as.data.frame(svd_PCA$rotation["VEGFA",])#Loading for VEGFA

#4. Identification of DEGs between the groups
library("ComplexHeatmap")

Patients_fac<-factor(sample.meta.data$patient) # Generates a factor object for patient
tissue<-factor(sample.meta.data$tissue) # Generates a factor object for tissue type(this is what i use for DGE analysis)

design_matrix <- model.matrix(~tissue) # Generates a martix which captures information about the experiemnt which is used for model fitting

dgelist.filtered.norm<-estimateDisp(dgelist.filtered.norm, design=design_matrix, robust = TRUE)# estimate the dispersion of read counts for each gene in RNA-Seq data by
#maximising the negative binomial likelyhood
plotBCV(dgelist.filtered.norm)
title('Genewise biological coefficient of variation (BCV) against gene abundance')

linear_model <- glmQLFit(dgelist.filtered.norm, design_matrix, robust = TRUE)#fitting quasi-likelihood (QL) negative binomial generalized linear models (GLM) to RNA-Seq count data
#The negative binomial distribution is commonly used to model count data, and it allows for overdispersion. The quasi-likelihood approach in edgeR is based on fitting a negative binomial liner model to the count data.


# estimated coefficients for each gene in the model. The coefficients correspond to the log-fold changes in gene expression associated with change in the corresponding variable.
head(linear_model$coefficients)
linear_model$coefficients["VEGFA",]

stat_test <- glmQLFTest(linear_model, coef = 2) # This function is used to perform likelihood ratio tests (LRTs) for differential gene expression between different conditions. 
#In our case the coef variable disinusies which groups are being compared. coef = 1 is BMvsCNS and coef = 2 is CNS vs BM


DE<-topTags(stat_test, n=Inf, sort.by = "logFC", adjust.method = "fdr")# A function from the limmia package and is used to test edgeR functions glmLRT, glmTreat, glmQLFTest
library(ggrepel)

ggplot() + #volcano plot showing DE gene before FDR and P-value filtering, shows the top 10 +ve and -ve DE genes 
  geom_point(data=filter(DE$table, FDR<0.05), aes(x=logFC, y=-log10(PValue)), colour = "black") +
  geom_point(data=filter(DE$table, FDR>=0.05), aes(x=logFC, y=-log10(PValue)), colour="grey") +
  geom_text_repel(data=slice_head(DE$table[order(DE$table$logFC, decreasing = FALSE), ], n=10), aes(x=logFC, y=-log10(PValue), label=Gene))+
  geom_text_repel(data=slice_head(DE$table[order(DE$table$logFC, decreasing = TRUE), ], n=10), aes(x=logFC, y=-log10(PValue), label=Gene))+
  ggtitle("Volcano plot of Differentally Expressed Genes before filtering")+
  theme_pubr()+
  theme(plot.title = element_text(hjust = 0.5))

DE$table <- DE$table[DE$table$PValue < 0.05,] #filters for pvalue

DE$table <- filter(DE$table, FDR<0.05) # filteres for FDR
DE$table$diffexpressed <- "NO" # Creates a new column with all "no" values
DE$table$diffexpressed[DE$table$logFC > 1] <- "UP" # All +Ve DE genes over this threshold are marked as Up regulated
DE$table$diffexpressed[DE$table$logFC < -1] <- "DOWN" #All -ve DE genes under this threshold are marked as Down regulated 

ggplot(data=DE$table, aes(x=logFC, y=-log10(PValue), col=diffexpressed)) + geom_point(size = 1) + theme_pubr()+ # Volcano plot that looks cool :)
  geom_text_repel(data=slice_head(DE$table[order(DE$table$logFC, decreasing = FALSE), ], n=10), aes(x=logFC, y=-log10(PValue), label=Gene))+
  geom_text_repel(data=slice_head(DE$table[order(DE$table$logFC, decreasing = TRUE), ], n=10), aes(x=logFC, y=-log10(PValue), label=Gene))+
  ggtitle("Volcano plot of positively and negatively Differentally Expressed Genes")+
  theme(plot.title = element_text(hjust = 0.5))

slice_head(DE$table[order(DE$table$logFC, decreasing = TRUE), ], n=20)


filtered_DEList <- DE$table[DE$table$PValue < 0.05,]
filtered_DEList <- filtered_DEList[DE$table$FDR < 0.05,]

edgeR_top <- slice_head(filtered_DEList[order(filtered_DEList$logFC, decreasing = TRUE),], n = 20)
view(slice_head(filtered_DEList[order(filtered_DEList$logFC, decreasing = TRUE),], n = 50)$Gene)
view(slice_head(filtered_DEList[order(filtered_DEList$logFC, decreasing = FALSE),], n = 50)$Gene)
edgeR_top[order(edgeR_top$PValue, decreasing = FALSE),]
#----------------------------------------------------------------------------------------------------------------

#COMPARE EdgeR with DESeq2!
  
library("DESeq2")

filtered_data[index.filter,]

dds <- DESeqDataSetFromMatrix(filtered_data, colData = dgelist.filtered.norm$samples, design = design_matrix)
res <- DESeq(dds)

keep <- rowSums(counts(dds)) >= 10
res <- res[keep, ]


res <- results(res, alpha = 0.05, lfcThreshold = 1) # Returns results from DESeq2 analysis
summary(res) # 637 Up regulated, 317 down regulated
res <- na.omit(res) # Removes rows with NA values

filtered_DESeq2 <- res[res$pvalue < 0.05,]
filtered_DESeq2 <- res[res$padj < 0.05,]
filtered_DESeq2 <- head(filtered_DESeq2[order(filtered_DESeq2$log2FoldChange, decreasing = TRUE),], n = 20)

top_DESeq2 <- filtered_DESeq2[order(filtered_DESeq2$pvalue, decreasing = FALSE),]

view(cbind(edgeR_top$Gene, row.names(top_DESeq2)))

#5. A heatmap of top 15 DEGs in each group

up_regulated <- slice_head(filtered_DEList[order(filtered_DEList$logFC, decreasing = TRUE), ], n = 15) # top 15 upregulated
down_regulated <- slice_head(DE$table[order(DE$table$logFC, decreasing = FALSE), ], n = 15) # Top 15 downregulated 

joined_genes<- c(up_regulated$Gene, down_regulated$Gene) # Joins the most up and down reuglated genes 

z.score.scaled.genes<-t(cpm.filtered.norm) %>% scale() %>% t() # Creates a z.score martix of the top 15 up and down regulated genes with patients and tissue
head(z.score.scaled.genes)

#Heatmap of Z-scores
Heatmap(matrix=z.score.scaled.genes[joined_genes,], 
        cluster_rows=TRUE, cluster_columns = TRUE, 
        show_row_names = TRUE,show_row_dend = TRUE,
        show_column_dend = TRUE,column_dend_side  = "top", 
        name = "Z-score", column_title = "Top 15 Postively and Negatively Differnetially Expressed Genes in CNS condition",
        row_split = 2, row_title = c("Up", "Down"))


