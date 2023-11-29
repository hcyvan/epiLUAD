library(ggplot2)
library(ggpubr)
library(ape)
library(DESeq2)
library(ChIPseeker)
library(org.Hs.eg.db)
library(DOSE)
library(KEGGREST)
library(clusterProfiler)
library(biomaRt)

# hclust------------------------------------------------------------------------
# data must have colname, and sampleID file
hclustplot <- function(data, sampleID){
  color_map <- data.frame(
    groups=c('LNM'   ,'AIS'    ,'MIA'    ,'IAC', 'LSQ'),
    colors=c('#00FF00','#00BFFF','#FFB90F','#FF0000','#800080'))
  #col <- color_map$colors ; names(col) <- color_map$groups
  pcolor <- color_map$colors[match(sampleID$cancer, color_map$groups)]
  
  dist <- dist(t(data),method = "euclidean")
  hclu <- hclust(dist,method = "ward.D")
  phyl <- as.phylo(hclu)
  
  # filename <- paste0(deparse(substitute(data)),'.pdf')
  # 
  # pdf(filename,width = 6,height = 6)
  plot(phyl, direction="downward",tip.color=pcolor, no.margin=TRUE,
       label.offset=2, cex=0.8,font=2,plot=TRUE)
  
  legend("bottom", legend =color_map$groups,col = color_map$colors,
         pch=15, pt.cex=1, cex=1, bty='n',horiz=TRUE)
  # dev.off()
}

# PCA---------------------------------------------------------------------------
# prepare data as hclustplot
pcaplot <- function(data,sampleID){
  color_map <- data.frame(
    groups=c('LNM'   ,'AIS'    ,'MIA'    ,'IAC', 'LSQ'),
    colors=c('#00FF00','#00BFFF','#FFB90F','#FF0000','#800080'))
  # col <- color_map$colors ; names(col) <- color_map$groups
  pcolor <- color_map$colors[match(sampleID$cancer, color_map$groups)]
  
  pca <- prcomp(t(data))
  df <- pca$x
  summ <- summary(pca)
  
  xlab_wgbs <- paste("PC1(",round(summ$importance[2,1]*100,2),"%)",sep="")
  ylab_wgbs <- paste("PC2(",round(summ$importance[2,2]*100,2),"%)",sep="")
  p <- ggplot(data=as.data.frame(df),aes(x=PC1,y=PC3))+
    geom_point(size = 3.0, alpha = 0.4 ,color=pcolor)+
    ggrepel::geom_text_repel(aes(label=colnames(data)),size=4.0,
                             segment.color=pcolor,color=pcolor,max.overlaps=10)+
    #scale_color_manual(values=color_map$colors,name=color_map$groups)+
    labs(x=xlab_wgbs,y=ylab_wgbs,color="",cex=2)+
    geom_hline(yintercept = 0,lty=2,col="black",lwd=1)+
    geom_vline(xintercept = 0,lty=2,col="black",lwd=1)+
    theme_bw()
  # filename <- paste0(deparse(substitute(data)),'.pdf')
  # ggsave(filename,plot = p,width = 6,height = 6)
}

# DESeq2------------------------------------------------------------------------
# data的列名为样品，行名为基因或者区间名，meta的列名为groups，行名为基因或者区间名
DEG_seq <- function(data,meta){
  dds <- DESeqDataSetFromMatrix(countData = round(data),
                                colData = meta,
                                design = ~ groups)
  
  keep <- rowSums(counts(dds)) >= 25000
  dds <- dds[keep,]
  dim(dds)
  dep <- DESeq(dds)
  resultsNames(dep)
  return(dep)
  # res <- results(dep,contrast = c("groups","LNM","LSQ"),
  #                alpha=0.001,tidy=T,pAdjustMethod = 'BH') 
  # # pAdjustMethod = c("holm", "hochberg", "hommel", "bonferroni", "BH"-default, "BY", "fdr", "none")
  # diff_gene <-subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
  # gene_up <- diff_gene_deseq2[diff_gene_deseq2$log2FoldChange < -1,]
  # gene_down <- diff_gene_deseq2[diff_gene_deseq2$log2FoldChange > 1,]
}

# Volcano plot------------------------------------------------------------------
volcanplot <- function(input){
  ggplot(input,aes(logFC,-log10(padj),color=change)) + 
    geom_point() + theme_bw() + 
    labs(title="ATAC-seq ATAC-seq differential accessibility",
         x=expression(log[2](LSQ/LNM)), y=expression(-log[10](q-value))) +
    scale_color_manual(values = c("red3","blue3","black")) +
    # scale_x_continuous(limits = c(-3, 3))+
    # xlim(-3,3)+
    geom_vline(xintercept=-1, linetype=2, colour="#B8423A",lwd=1) +
    geom_vline(xintercept=1, linetype=2, colour="#B8423A",lwd=1) +
    geom_hline(yintercept=-log10(0.01), linetype=2, colour="#B8423A",lwd=1) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
}


# Gene annotation---------------------------------------------------------------
annoGene <- function(input){
  require(TxDb.Hsapiens.UCSC.hg38.knownGene)
  txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
  
  input <- makeGRangesFromDataFrame(input,keep.extra.columns=T)
  anno <- annotatePeak(input,tssRegion = c(-5000,5000),TxDb = txdb,
                       annoDb = "org.Hs.eg.db")
  anno_input <- as.data.frame(anno)
  plotAnnoPie(anno)
  return(anno_input)
  # return(anno_read)
  # save(anno,anno_read,file = "annopeak.rdata")
}

# write .bed file---------------------------------------------------------------
wrtbed <- function(input){
  filename <- paste0(deparse(substitute(input)),'.bed')
  write.table(input,file = filename, 
              row.names = F, col.names = F,quote = F,sep = "\t")
}

# enrich------------------------------------------------------------------------
enrich <- function(gene_symbol){
  # change SYMBOL to ENTREZID
  gene <- bitr(gene_symbol,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
  
  # 富集分析
  GO <<- enrichGO(gene$ENTREZID,
                  OrgDb = 'org.Hs.eg.db',
                  keyType = "ENTREZID",
                  ont = "ALL",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = T)
  
  KEGG <<- enrichKEGG(gene$ENTREZID,
                      organism = 'hsa',
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05)
}

# boxplot-----------------------------------------------------------------------

boxp <- function(input,col){
  color_map <- data.frame(
    groups=c('LNM'   ,'AIS'    ,'MIA'    ,'IAC'     ,'LSQ'     ),
    colors=c('#00FF00','#00BFFF','#FFB90F','#FF0000','#800080'))
  col <- color_map$colors ; names(col) <- color_map$groups
  
  ggplot(input, aes(x=group, y=value,fill=group)) +
    scale_fill_manual(values=col)+
    geom_boxplot() +
    geom_jitter(shape=17, position=position_jitter(0.2), colour='black', size=1.5)+
    #stat_summary(fun.data=mean_sdl, mult=1, 
    #             geom="pointrange", color="black")+
    stat_summary(fun.data=mean_sdl, geom="pointrange", color="black")+
    theme_classic()+
    stat_compare_means( comparisons = list(c('LNM','AIS'),c('LNM','MIA'),c('LNM','IAC'),c('LNM','LSQ')),
                        label = 'p.signif', method = "t.test")+
    xlab("")+
    ylab("Mean Methylation Level")+
    theme(
      legend.position="right",
      axis.title.x = element_text(size=0),
      axis.title.y = element_text(size=12),
      axis.text = element_text(size = 12,colour="black"),
      legend.title = element_blank(),
      legend.text = element_text(size=12)
    )+
    guides(colour = guide_legend(override.aes = list(shape = 12,size=10)))
}

# scatter plot -----------------------------------------------------------------
scatterp <- function(df,col){
  
  filename <- paste0(deparse(substitute(df)),'.pdf')
  rho <- cor.test(df$wgbs, df$atac) #Calculate correlation
  
  ggplot(df,aes(x=wgbs,y=atac,color =group))+
    geom_smooth(method = "lm", se=FALSE, 
                color="black", formula = y~x) +
    geom_point(size=3,
               alpha=0.7)+
    scale_color_manual(values = col)+
    geom_text(aes(x = mean(wgbs), y = max(atac), 
                  label = paste("Rho =",round(rho$estimate,digits = 3),
                                "\n","p-value=",round(rho$p.value,digits = 3))),
              color = "black", size = 4)+
    ggrepel::geom_text_repel(aes(label=rownames(df)),size=3.0,
                             segment.color=pcolor,color=pcolor,max.overlaps=20)
  #theme(panel.grid = element_blank())
  
  ggsave(filename,width = 7, height = 6)
}

# heatmap plot------------------------------------------------------------------


# ensambleID to SYMBOL----------------------------------------------------------
esmbTosymb <- function(input){
  load('ensembl.rds')
  attr <- c("chromosome_name", "start_position", "end_position", "strand",
            "ensembl_gene_id","hgnc_symbol")
  input_info <- getBM(attributes = attr,
                          filters = "ensembl_gene_id",
                          values = input,
                          mart = ensembl)
  input_info$chromosome_name <- paste0("chr",input_info$chromosome_name)
  
  # chromosome_order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
  #                       "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
  #                       "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
  #                       "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
  # input_info <- input_info[order(match(input_info$chromosome_name,chromosome_order), 
  #                                        input_info$start_position), ]
  # input_info <- input_info[input_info$hgnc_symbol != '',]
  return(input_info)
}