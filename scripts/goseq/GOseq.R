## --------------------------- GENE ONTOLOGY ENRICHMENT ANALYSIS ---------------------------
#This scripts perform GO terms enrichment analysis in a given list of genes.
#three type of GO terms test independent (BP,MF and CC).
#Do hypergeometric and FDR adjust.
#Author: Leonardo Jo

## Library
##if (!requireNamespace("BiocManager", quietly = TRUE))
##  install.packages("BiocManager")
##BiocManager::install("goseq")
library(goseq)
library(tidyverse)

#read all the tables
all.genes<-read.table("ITAG4.1_Gene_Model.lengths.txt", sep="\t",stringsAsFactors=FALSE, header=T)
GO_BP<-read.table("GOterms.phytozome.splitBP.summary.txt", sep="\t",stringsAsFactors=FALSE)
GO_MF<-read.table("GOterms.phytozome.splitMF.summary.txt", sep="\t",stringsAsFactors=FALSE)
GO_CC<-read.table("GOterms.phytozome.splitCC.summary.txt", sep="\t",stringsAsFactors=FALSE)


#this creates a vector with all the genes to assay
assayed.genes<-all.genes$Sol.ID
gene.length<-as.integer(all.genes$size)
names(gene.length)<-assayed.genes

#List all the gene files in the director
filelist=list.files(path="GENESETS/",pattern=".txt")


#Start analysis.
for (i in filelist){
  de.genes=read.table(paste0("GENESETS/",i),header=T,sep="\t") %>%
    pull(Geneid)

  gene.vector<-as.integer(assayed.genes%in%de.genes)
  names(gene.vector)=assayed.genes

  pwf<-nullp(gene.vector,bias.data=gene.length,plot.fit=F)
  
  #GO terms enrichment of Biological Process
  go.wall.BP <- goseq(pwf,gene2cat=GO_BP,method="Hypergeometric", use_genes_without_cat=T)
  go.wall.BP$qval <- p.adjust(go.wall.BP$over_represented_pvalue,method="BH")
  enriched.go.BP <-subset(go.wall.BP, go.wall.BP$qval <= 0.05)
  
  #if have Go terms qval lower than 0.05, draw a tree.
  if(nrow(enriched.go.BP)>0){
    enriched.go.wall.annot.BP=enriched.go.BP

    #for table
    Data1=NULL
    for(j in c(1:nrow(enriched.go.wall.annot.BP))){
      referencesGenes=GO_BP[GO_BP[,2]==enriched.go.wall.annot.BP[j,1],1]
      NumberOfReferencesGenes=length(referencesGenes)
      Genes=intersect(referencesGenes,de.genes)
      entries=NULL
      for(k in c(1:length(Genes))){
        entries=paste(entries,Genes[k],sep=",")
      }
      entries=substr(entries,2,nchar(entries))
      NumberOfGenes=length(Genes)
      Data2=cbind(NumberOfGenes,NumberOfReferencesGenes,entries)
      Data1=rbind(Data1,Data2)
    }
    enriched.go.wall.annot.BP=cbind(enriched.go.wall.annot.BP,Data1)
    enriched.go.wall.annot.BP=enriched.go.wall.annot.BP[order(enriched.go.wall.annot.BP$qval),]
  }
  if(nrow(enriched.go.BP)==0){
    enriched.go.wall.annot.BP=NULL
  }
  
  #GO terms enrichment of molecular function.
  go.wall.MF <- goseq(pwf,gene2cat=GO_MF,method="Hypergeometric", use_genes_without_cat=T)
  go.wall.MF$qval <- p.adjust(go.wall.MF$over_represented_pvalue,method="BH")
  enriched.go.MF <-subset(go.wall.MF, go.wall.MF$qval <= 0.05)
  
  #if have Go terms qval lower than 0.05, draw a tree.
  if(nrow(enriched.go.MF)>0){
    enriched.go.wall.annot.MF=enriched.go.MF

    
    #for table
    Data1=NULL
    for(j in c(1:nrow(enriched.go.wall.annot.MF))){
      referencesGenes=GO_MF[GO_MF[,2]==enriched.go.wall.annot.MF[j,1],1]
      NumberOfReferencesGenes=length(referencesGenes)
      Genes=intersect(referencesGenes,de.genes)
      entries=NULL
      for(k in c(1:length(Genes))){
        entries=paste(entries,Genes[k],sep=",")
      }
      entries=substr(entries,2,nchar(entries))
      NumberOfGenes=length(Genes)
      Data2=cbind(NumberOfGenes,NumberOfReferencesGenes,entries)
      Data1=rbind(Data1,Data2)
    }
    enriched.go.wall.annot.MF=cbind(enriched.go.wall.annot.MF,Data1)
    enriched.go.wall.annot.MF=enriched.go.wall.annot.MF[order(enriched.go.wall.annot.MF$qval),]
  }
  if(nrow(enriched.go.MF)==0){
    enriched.go.wall.annot.MF=NULL
  }
  
  #GO terms enrichment of Cellular component.
  go.wall.CC <- goseq(pwf,gene2cat=GO_CC,method="Hypergeometric", use_genes_without_cat=T)
  go.wall.CC$qval <- p.adjust(go.wall.CC$over_represented_pvalue,method="BH")
  enriched.go.CC <-subset(go.wall.CC, go.wall.CC$qval <= 0.05)
  
  #   #if have Go terms qval lower than 0.05, draw a tree.
  if(nrow(enriched.go.CC)>0){
    enriched.go.wall.annot.CC=enriched.go.CC

    
    #for table
    Data1=NULL
    for(j in c(1:nrow(enriched.go.wall.annot.CC))){
      referencesGenes=GO_CC[GO_CC[,2]==enriched.go.wall.annot.CC[j,1],1]
      NumberOfReferencesGenes=length(referencesGenes)
      Genes=intersect(referencesGenes,de.genes)
      entries=NULL
      for(k in c(1:length(Genes))){
        entries=paste(entries,Genes[k],sep=",")
      }
      entries=substr(entries,2,nchar(entries))
      NumberOfGenes=length(Genes)
      Data2=cbind(NumberOfGenes,NumberOfReferencesGenes,entries)
      Data1=rbind(Data1,Data2)
    }
    enriched.go.wall.annot.CC=cbind(enriched.go.wall.annot.CC,Data1)
    enriched.go.wall.annot.CC=enriched.go.wall.annot.CC[order(enriched.go.wall.annot.CC$qval),]
  }
  if(nrow(enriched.go.CC)==0){
    enriched.go.wall.annot.CC=NULL
  }
  
  if(nrow(enriched.go.BP)+nrow(enriched.go.CC)+nrow(enriched.go.MF)>0){
    enriched.go.wall.annot=NULL
    enriched.go.wall.annot=rbind(enriched.go.wall.annot.BP,enriched.go.wall.annot.MF,enriched.go.wall.annot.CC)
    enriched.go.wall.annot<- enriched.go.wall.annot[,-c(2:5)]
    name=paste(substr(i,1,nchar(i)-4),"GOseq.enrichment.SOL.txt",sep=".")
    write.table(enriched.go.wall.annot,name,sep="\t",quote=FALSE,row.names=F)
  }
}
