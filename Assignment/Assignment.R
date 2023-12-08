library("tidyverse")
library("edgeR")
library("clusterProfiler")
library("ggpubr")
library("Seurat")


count_matrix <- read_tsv("B-ALL.txt")

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

#~~~~~~~~~~~~~~~~~~~~~~~Filtering~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Counts per million calculations 
cpm<-cpm(dgelist)
lcpm <- cpm(dgelist, log=TRUE)
#This function transforms a matrix into a long format containing gene, logCPM, and patients/tissue numbers 
#I use this function multipul times for checking filtering techniques 
long_transormation <- function(x,value){
  long<-as.data.frame(x) %>% 
    rownames_to_column("Gene") %>%
    pivot_longer(-Gene, values_to= value, names_to="Patients")
  return(long)
}

#Piping long transformed data into a density plot 
long_transormation(lcpm,"LogCPM") %>% ggplot() + geom_density(aes(LogCPM, colour=Patients))


index.filter<-filterByExpr(dgelist, group=tissue, min.count = 10, min.total.count = 15, min.prop = 0.7)

dgelist.filtered<-dgelist[index.filter, , keep.lib.sizes=FALSE]
dim(dgelist.filtered)

lcpm.filtered<-cpm(dgelist.filtered, log=TRUE)
long_transormation(lcpm.filtered, "LogCPM") %>% ggplot() + geom_density(aes(LogCPM, colour=Patients))

plotMeanVar(dgelist.filtered)#views the mean varaince at the gene-level 


#~~~~~~~~~~~~~~~~~~~~~~~Normalisation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Normalisation 
dgelist.filtered.norm<-calcNormFactors(dgelist.filtered, method="TMM")
dgelist.filtered.norm$samples$norm.factors
plotMeanVar(dgelist.filtered.norm) #Normalisation doesn't make a massive difference since the data has already been normalised(but not for library size)

#Gene Names are removed from column after normalisation!
row.names(dgelist.filtered.norm$counts) <- dgelist.filtered.norm$genes[,1]



cpm.filtered.norm<-cpm(dgelist.filtered.norm)
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
  theme(plot.title = element_text(hjust = 0.5))+
  theme_pubr()

ggplot(df.plotting) + geom_density(aes(LogCPM, colour=Patients))+
  ggtitle("LogCPM density across all patient samples after TMM normalisation")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_pubr()


long.lcpm.filtered %>% left_join(sample.meta.data, by=c("Patients"="sample.names")) %>% left_join(gene.meta.data) %>% 
  ggplot()+geom_boxplot(aes(x=Patients, y = LogCPM, color = tissue))+
  ggtitle("Varaince before TMM filtering")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_pubr()


ggplot(df.plotting)+geom_boxplot(aes(x=Patients, y = LogCPM, color = tissue))+
  ggtitle("Varaince after TMM filtering")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_pubr()



ggplot(df.plotting %>% filter(Gene=="VEGFA")) + 
  geom_point(aes(x=tissue, y=CPM, color = patient)) +
  ggtitle("VEGFA expression in BM and CNS")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_pubr()



#3. PCA analysis (plot for dimensionality reduction and scree plot)
library("factoextra")

svd_PCA <-
  lcpm.filtered.norm %>%
  t() %>%
  prcomp(scale = F)

summary(svd_PCA)

fviz_eig(svd_PCA, addlabels = T)
fviz_pca_ind(svd_PCA, repel = F) 
#fviz_pca_ind plots a graph of individals and creates a graph 
fviz_pca_ind(svd_PCA, 
             habillage = sample.meta.data$patient,
             label = "none" ,
             invisible = "quali") + ggtitle("Sample")

fviz_pca_ind(svd_PCA, 
             habillage = sample.meta.data$tissue,
             label = "none",
             invisible = "quali") + ggtitle("Tissue.type")

rownames(svd_PCA$rotation) <- rownames(lcpm.filtered.norm)

selected_var <- svd_PCA$rotation["VEGFA", , drop = FALSE]


fviz_pca_var(svd_PCA, select.var = list(cos2 = 20), repel = TRUE)

as.data.frame(svd_PCA$rotation[,1:2]) %>%
  slice_max(n=10, order_by = PC1)


sort(svd_PCA$rotation[,1:2], decreasing = TRUE)
head(svd_PCA)


#4. Identification of DEGs between the groups
library("ComplexHeatmap")

Patients_fac<-factor(sample.meta.data$patient)
tissue<-factor(sample.meta.data$tissue)

design_matrix <- model.matrix(~tissue)


dgelist.filtered.norm<-estimateDisp(dgelist.filtered.norm, design=design_matrix)


linear_model <- glmQLFit(dgelist.filtered.norm, design_matrix)#Uses as quasi-likelyhood test 


head(linear_model$coefficients)
sort(linear_model$coefficients)
linear_model$coefficients["VEGFA",]

stat_test <- glmQLFTest(linear_model,coef = 2) # Creates a DGELRT object storing 

head(stat_test)

DE<-topTags(stat_test, n=Inf, sort.by = "logFC")# A function from the limmia package and is used to test edgeR functions glmLRT, glmTreat, glmQLFTest
library(ggrepel)

ggplot() + 
  geom_point(data=filter(DE$table, FDR<0.05), aes(x=logFC, y=-log10(PValue)), colour = "black") +
  geom_point(data=filter(DE$table, FDR>=0.05), aes(x=logFC, y=-log10(PValue)), colour="grey") +
  geom_text_repel(data=slice_head(DE$table[order(DE$table$logFC, decreasing = FALSE), ], n=10), aes(x=logFC, y=-log10(PValue), label=Gene))+
  geom_text_repel(data=slice_head(DE$table[order(DE$table$logFC, decreasing = TRUE), ], n=10), aes(x=logFC, y=-log10(PValue), label=Gene))+
  ggtitle("Volcano plot of Differentally Expressed Genes before filtering")+
  theme(plot.title = element_text(hjust = 0.5))

DE$table <- DE$table[DE$table$PValue < 0.05,]
DE$table <- filter(DE$table, FDR<0.05)
DE$table$diffexpressed <- "NO" # Creates a new column with all "no" values
DE$table$diffexpressed[DE$table$logFC > 0.6] <- "UP"
DE$table$diffexpressed[DE$table$logFC < -0.6] <- "DOWN"

ggplot(data=DE$table, aes(x=logFC, y=-log10(PValue), col=diffexpressed)) + geom_point() + theme_pubr()+
  geom_text_repel(data=slice_head(DE$table[order(DE$table$logFC, decreasing = FALSE), ], n=10), aes(x=logFC, y=-log10(PValue), label=Gene))+
  geom_text_repel(data=slice_head(DE$table[order(DE$table$logFC, decreasing = TRUE), ], n=10), aes(x=logFC, y=-log10(PValue), label=Gene))+
  ggtitle("Volcano plot of positively and negatively Differentally Expressed Genes")+
  theme(plot.title = element_text(hjust = 0.5))







dgelist.filtered.norm$counts["VEGFA",]

#5. A heatmap of top 15 DEGs in each group

top_15 <- slice_head(DE$table, n=15)

DE$table["TNFRSF4",]

DE$table["VEGFA",]

filtered_DEList <- DE$table[(DE$table$logFC >= 1) & (DE$table$PValue < 0.05),]

dim(filtered_DEList)

top_15 <- slice_head(filtered_DEList[order(filtered_DEList$logFC, decreasing = TRUE), ], n = 15)

p_value_filtered_top_15 <- slice_head(top_15[order(top_15$PValue, decreasing = FALSE), ], n=15)



head(DE$table)

joined_genes<- c(CNSvBM$Gene, BMvCNS$Gene)

DE$table[DE$table$Gene %in% joined_genes, ] %>% left_join(df.plotting, join_by("Gene" == "Gene"))

subset_DE <- data.frame(DE$table[joined_genes,])


Heatmap(matrix=subset_DE[2], cluster_rows=FALSE, cluster_columns = FALSE, show_row_names = TRUE)

CNSvBM <- slice_head(filtered_DEList[order(filtered_DEList$logFC, decreasing = TRUE), ], n = 15)
BMvCNS <- slice_head(DE$table[order(DE$table$logFC, decreasing = FALSE), ], n = 15)

merged <- data.frame(CNSvBM, BMvCNS)


Heatmap()




#DESEQ 2 analysis
