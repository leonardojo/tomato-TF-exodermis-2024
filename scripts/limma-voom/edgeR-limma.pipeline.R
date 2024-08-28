## --------------------------- DEG Analysis ---------------------------
#This script performs differentially expression analysis using the limma-voom package
#Author: Leonardo Jo
#Adapted from: https://github.com/plant-plasticity/tomato-root-atlas-2020/tree/master/Scripts/Limma-voom

library(edgeR)
library(tidyverse)
library(Glimma)
library(cowplot)
options(stringsAsFactors = F)
##install.packages("BiocManager", repos = "https://cloud.r-project.org") ## if you need to install edgeR
##BiocManager::install("edgeR") ## if you need to install edgeR

## Setting up the working directory to the folder where this R script is saved
path <- rstudioapi::getActiveDocumentContext()$path
Encoding(path) <- "UTF-8"
setwd(dirname(path))

## Then we load in your metadata
meta <- read.table("meta-tomato.txt",sep="\t",header=T)
sample.names <- meta %>% pull(Name)
head(meta)

## And load in your count data
counts <- read.table("mastercounts.tsv",sep="\t",header=T) %>%
  select(target_id,all_of(sample.names))
head(counts)

## Preparing Count matrix
count.matrix <- counts %>%
  column_to_rownames(var = "target_id")

## Removing non-expressed genes
zero.expressed.genes <- count.matrix %>%
  rownames_to_column(var = "Geneid") %>%
  pivot_longer(-Geneid) %>%
  group_by(Geneid) %>%
  summarise(sum = sum(value)) %>%
  filter(sum == 0)
  
count.matrix.filtered <- count.matrix %>%
  rownames_to_column(var = "Geneid") %>%
  filter(!Geneid %in% zero.expressed.genes$Geneid) %>%
  column_to_rownames(var = "Geneid")
  
cat("Initial number of genes:",nrow(count.matrix),"\n",
    "Genes not expressed in any tissues:",nrow(zero.expressed.genes),"\n",
    "Remaining number of genes:",nrow(count.matrix.filtered),"\n")

# Create a folder where your output is going. 
# You can change this name if you want to re-do analyses and not overwrite.
outDir = "Plots/"
outDir2 = "DEGs/"
outDir3 = "Tables/"
dir.create(outDir, showWarnings=T)
dir.create(outDir2, showWarnings=T)
dir.create(outDir3, showWarnings=T)
## Ignore warnings, if the table already exists, it won't change the outcome of the script


###### Design matrix
## Convert experimental metadata to factors for the design
experimentFactors <- lapply(apply(meta,2,split,""),unlist)
experimentFactors <- as.data.frame(lapply(experimentFactors,as.factor))
colnames(experimentFactors)

####Simplest design taking into account all possible interactions
Groups <- as.factor(str_sub(names(count.matrix.filtered),1,-3))
design <- model.matrix(~0+Groups) 
colnames(design) <- levels(Groups)
levels(Groups)

## Running edgeR
y <- DGEList(counts=count.matrix.filtered, group=Groups)

## one sample has 3 biological replicates, therefore, to be considered,
## a gene needs to be expressed (0.5 cpm) in at least 3 samples
keep <- rowSums(cpm(y)>0.5) >= 3
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]

## recalculate norm factor for each library after filtering
y <- calcNormFactors(y,method="TMM")
normalizedExpression <- cpm(y) ## get CPMs with the TMM normalization method

## saving TMM normalized expression table
normalizedExpression %>%
  data.frame() %>%
  rownames_to_column("Geneid") %>%
  write.table(paste0(outDir3,"TMMnormalizedcpm.txt"),sep="\t",row.names=F)


## saving TMM normalized expression table, averaged between bioreps
normalizedExpression %>%
  data.frame() %>%
  rownames_to_column("Geneid") %>%
  pivot_longer(-Geneid) %>%
  separate(name,into=c("sample","biorep")) %>%
  group_by(Geneid,sample) %>%
  summarise(average.cpm=mean(value,na.rm=T)) %>%
  pivot_wider(names_from=sample,values_from=average.cpm) %>%
  write.table(paste0(outDir3,"TMMnormalizedcpm.AVG.txt"),sep="\t",row.names=F)


## MDS plot
MDS.plot <- plotMDS(y)
tibble(sample = rownames(y$samples),
       dim1=MDS.plot$x,
       dim2=MDS.plot$y) %>%
  separate(sample,into=c("sample","biorep"),sep="[.]",remove = F) %>%
  ggplot() +
  geom_point(aes(dim1,dim2,colour=sample),size=7.5) +
  coord_fixed() +
  theme_bw(base_size = 18) +
  theme(panel.grid=element_blank()) +
  xlab(paste0("dim1 (",round(MDS.plot$var.explained[1]*100,2),"%)")) +
  ylab(paste0("dim2 (",round(MDS.plot$var.explained[2]*100,2),"%)"))
ggsave(paste0(outDir,"MDSplot.png"),dpi=600,width=8,height=8)

## Dendogram
scaled <- as.matrix(y$counts) %>%
  t() %>%
  scale() %>%
  t()

hc <- hclust(as.dist(1-cor(scaled, method="spearman")), method="complete") # Clusters columns by Spearman correlation.
sampleTree = as.dendrogram(hc, method="average")

tmpSave <- paste(outDir,"Dendogram",".png",sep="")
png(tmpSave,width=1600,height=800)
par(cex=1, font=2)
plot(sampleTree, xlab="", ylab="", main="", sub="", axes=FALSE)
par(cex=1)
title(ylab="Height", main="Sample Clustering")
axis(2)
dev.off()



# Variance stabilization with Voom - draw a graph of mean-variance trend
tmpSave <- paste(outDir,"Mean-variance trend",".png",sep="")
png(tmpSave)
v <- voom(y,design,plot = T)
dev.off()

# boxplot with stabilized mean-variance across samples
a <- y$counts %>%
  data.frame() %>%
  rownames_to_column("Geneid") %>%
  pivot_longer(-Geneid) %>%
  mutate(value = log(value,2)) %>%
  ggplot() +
  geom_boxplot(aes(name,value,fill=name),show.legend = F) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1),
        panel.grid = element_blank()) +
  ylab("log2[counts]") +
  xlab("Library") +
  ggtitle("Non normalized counts")

b <- v$E %>% 
  data.frame() %>%
  rownames_to_column(var="Geneid") %>%
  pivot_longer(-Geneid) %>%
  ggplot() +
  geom_boxplot(aes(name,value,fill=name),show.legend = F) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1),
        panel.grid = element_blank()) +
  ylab("log2[counts]") +
  xlab("Library") +
  ggtitle("TMM normalized counts with voom variance stabilization")

plot_grid(a,b,nrow=1)
ggsave(paste0(outDir,"NormalizedCounts-boxplot.png"),dpi=600,width=12,height=8)


## Let's start doing the differential analysis
### Set up contrasts
levels(Groups)  #This command lists the groups that you can use in your comparisons

# Change the contrasts to be ones that answer YOUR RESEARCH QUESTION
cont.matrix= makeContrasts(
  "myb2KO" = myb2KO-wt,
  "myb2OEX" = myb2OEX-wt,
  "myb5KO" = myb5KO-wt,
  "myb5OEX" = myb5OEX-wt,
  "wrky4KO" = wrky4KO-wt,
  "wrky5KO" = wrky5KO-wt,
  "wrkyOEX" = wrkyOEX-wt,
  levels=design)

# Fit your data and design to a linear model
fit <- lmFit(v, design)
fit

# Estimate contrast for each gene
fit2 <- contrasts.fit(fit, cont.matrix)
fit2

# Empirical Bayes smoothing of standard errors (shrinks standard errors that are much larger 
# or smaller than those from other genes towards the average standard error)
fit2 <- eBayes(fit2)

# Multiple Testing Across Genes and Contrasts
results <- decideTests(fit2,p.value = 0.05,adjust.method="BH") 
## adjust.method for p.value correction, BH: Barnes-Hutch
## p.value threshold for corrected pvalues
## lfc minimum absolute log2 fold change (1 means at least 2 fold higher or lower)

# What is in the "results" object? Quick check to look at first six rows:
head(results)
## 1 is upregulated, -1 is downregulated

# Check the summary of the results
summary(results)

# Save the summary as a .txt file
summary(results) %>% data.frame() %>%
  pivot_wider(names_from="Var2",values_from="Freq") %>%
  write.table(paste0(outDir3,"DESummary.txt"),sep="\t",quote = F,row.names=F)

# Draw a bar graph of the summary
summary(results) %>% data.frame() %>%
  filter(Var1 != "NotSig") %>%
  ggplot() +
  geom_bar(aes(Var2,Freq,fill=Var1),stat="identity",colour="black",position=position_dodge2(1)) +
  geom_text(aes(Var2,Freq+15,group=Var1,label=Freq),position=position_dodge2(0.9)) +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle=90,vjust=0.5,hjust=1)) +
  scale_fill_manual(values=c("navy","red3")) +
  scale_y_continuous(expand = expansion(mult = c(0,.1))) +
  xlab("Timepoint") +
  ylab("Number of DEGs")
  ggsave(paste0(outDir,"summarysDEGs.png"),dpi=300,width=8,height=5)

  
### Next, we will generate output files with all the DE statistics for each gene and each comparison
contrasts.names <-  colnames(cont.matrix)
for (k in contrasts.names) {
  tmp <- topTable(fit2, coef=k,number = Inf,sort.by = "none") %>%
    rownames_to_column("Geneid")
  ##allDEGs
    all.degs <- tmp %>%
      filter(adj.P.Val <= 0.05)
    ##UP.regulated
    UP.reg <- tmp %>%
      filter(adj.P.Val <= 0.05 & logFC >= 0)
    ##DOWN.regulated
    DOWN.reg <- tmp %>%
      filter(adj.P.Val <= 0.05 & logFC <= 0)
    write.table(all.degs,paste0(outDir2,k,"_allDEGs.txt"),sep="\t",row.names=F)
    write.table(UP.reg,paste0(outDir2,k,"_UPreg.txt"),sep="\t",row.names=F)
    write.table(DOWN.reg,paste0(outDir2,k,"_DOWNreg.txt"),sep="\t",row.names=F)
    write.table(tmp,paste0(outDir2,k,"_allGenes.txt"),sep="\t",row.names=F)
  }

### Create a single file with the logFC and p.value of all DEGs into one file
contrasts.names <-  colnames(cont.matrix)
final.tmp <- NULL
for (k in contrasts.names) {
  tmp <- topTable(fit2, coef=k,number = Inf,sort.by = "none") %>%
    rownames_to_column("Geneid") %>%
    select(Geneid,adj.P.Val,logFC) %>%
    pivot_longer(-Geneid) %>%
    mutate(dataset = k) 
    final.tmp <- rbind(final.tmp,tmp)
}

## Putting all in one table
final.tmp2 <- final.tmp %>% 
  unite(name,c(dataset,name)) %>%
  pivot_wider(names_from=name,values_from=value)

## filtering DEGs 
all.degs <- NULL
for (k in contrasts.names) {
  tmp <- topTable(fit2, coef=k,number = Inf,sort.by = "none") %>%
    rownames_to_column("Geneid") %>%
    filter(adj.P.Val <= 0.05) %>%
    select(Geneid)
  all.degs <- rbind(all.degs,tmp)
}

all.degs <- all.degs %>%
  select(Geneid) %>%
  unique

final.tmp2 %>%
  filter(Geneid %in% all.degs$Geneid)

write.table(final.tmp2,paste0(outDir3,"summaryDEGs.txt"),sep="\t",row.names=F)

