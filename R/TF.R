source('./R/base.R')
source('./R/local/base.R')
library(dplyr)
library(reshape2)
library(ComplexHeatmap)
library(VennDiagram)


#----------------------------------------------------------------------------------------------------------------------
# Figure 4A. Common DARs and SC-DMRs Homer TFs
#----------------------------------------------------------------------------------------------------------------------
motifDAR<-read.csv(file.path(CONFIG$dataResult, 'atac.dar.homer.motifs.csv'))
motifSRDMR<-read.csv(file.path(CONFIG$dataResult, 'srdmr.homer.motifs.csv'))
darHomerTFs<-unique(motifDAR$TFs)
srdmrHomerTFs<-unique(motifSRDMR$TFs)
tf124<-intersect(darHomerTFs,srdmrHomerTFs)
cat(sprintf("Common TFs: %d", length(tf124)))
write.table(tf124, file.path(CONFIG$dataIntermediate,'tf', 'tf124.txt'), col.names = FALSE,row.names = FALSE,quote = FALSE)


venn.plot<-venn.diagram(x = list(DAR = darHomerTFs, `SR-DMR` = srdmrHomerTFs),
                        fill = c("red","blue"),
                        cat.cex = 2,
                        cex = 2,
                        filename = NULL)
saveImage2("tf124.venn.pdf",width = 4,height = 4)
grid.draw(venn.plot)
dev.off()

#----------------------------------------------------------------------------------------------------------------------
# Common DARs and SC-DMRs Homer TFs
#----------------------------------------------------------------------------------------------------------------------
darDeseq2<-saveRDS(file.path(CONFIG$dataIntermediate,'atac', 'darDeseq2.rds'))

diffPeak<-lapply(darDeseq2, function(x){
  rownames(x$diff)
})%>%unlist()%>%unique()

tf124<-loadData2(file.path(CONFIG$dataIntermediate,'tf', 'tf124.txt'),file.format = 'bed',header = FALSE)$V1
allMotif<-loadData2(file.path(CONFIG$dataIntermediate, 'atac','homer.mask','all.motif.region.txt'),file.format = 'bed')
allMotifFilter<-allMotif[allMotif$MotifScore>=5,]
allMotifFilter<-allMotifFilter[allMotifFilter$PositionID%in%diffPeak,]
tf124Motif<-allMotifFilter[sapply(strsplit(allMotifFilter$`Motif Name`, '\\('),function(x){x[1]})%in%tf124,]
saveRDS(tf124Motif,file.path(CONFIG$dataIntermediate, 'atac','homer.mask','all.motif.124.dar.rds'))
#----------------------------------------------------------------------------------------------------------------------
# Figure 4B. 124 TFs cluster
#----------------------------------------------------------------------------------------------------------------------
allMotif<-readRDS(file.path(CONFIG$dataIntermediate, 'atac','homer.mask','all.motif.124.dar.rds'))
allMotifList<-split(allMotif, allMotif$`Motif Name`)
motifs<-unique(allMotif$`Motif Name`)
motifOverlapList<-lapply(motifs, function(x){
  print(x)
  lapply(motifs, function(y){
    motifX<-allMotifList[[x]]
    motifY<-allMotifList[[y]]
    xy<-intersect(unique(motifX$PositionID),unique(motifY$PositionID))
    length(xy)
  })
})
motifOverlap<-sapply(motifOverlapList, function(x){
  unlist(x)
})
colnames(motifOverlap)<-motifs
rownames(motifOverlap)<-motifs

addHyper<-function(motifOverlapDf){
  motifOvelapMap<-diag(motifOverlap)
  names(motifOvelapMap)<-motifs
  motifOverlapDf$tf1Num<-motifOvelapMap[match(motifOverlapDf$motif1,names(motifOvelapMap))]
  motifOverlapDf$tf2Num<-motifOvelapMap[match(motifOverlapDf$motif2,names(motifOvelapMap))]
  motifOverlapDf$Total<-length(unique(allMotif$PositionID))
  p<-apply(motifOverlapDf, 1, function(x){
    M <- as.numeric(x[5])
    N <- as.numeric(x[6]) - M
    n <- as.numeric(x[4])
    k <- as.numeric(x[3])
    phyper(k, M, N, n, lower.tail = FALSE)
  })
  motifOverlapDf$hyperP<-p
  motifOverlapDf$hyperFDR<-p.adjust(p, method = "BH")
  motifOverlapDf<-arrange(motifOverlapDf, hyperFDR)
  motifOverlapDf
}

motifOverlapAll<-melt(motifOverlap)
colnames(motifOverlapAll)<-c('motif1','motif2', 'overlapNum')
motifOverlapAll<-addHyper(motifOverlapAll)
motifOverlapAll$tag <- -log10(motifOverlapAll$hyperFDR+1e-322)

m<-dcast(motifOverlapAll,motif1~motif2,value.var = 'tag')
m<-as.matrix(m[,-1])
colnames(m)<-sapply(strsplit(motifs,'\\('),function(x){x[1]})
rownames(m)<-sapply(strsplit(motifs,'\\('),function(x){x[1]})

myDist <- function(x) dist(x, method = "euclidean")
myHclust <- function(x) hclust(x, method = "ward.D2")
saveImage2("tf.124.heatmap.pdf",width = 9,height = 8)
Heatmap(m,
        name='-log10(FDR)',
        column_names_gp = grid::gpar(fontsize = 4),
        row_names_gp = grid::gpar(fontsize = 4),
        clustering_method_rows = 'ward.D2',
        clustering_method_columns = 'ward.D2',
        clustering_distance_rows = myDist,
        clustering_distance_columns = myDist
)
dev.off()
#----------------------------------------------------------------------------------------------------------------------
# Figure S3A-D. 124 TFs cluster
#----------------------------------------------------------------------------------------------------------------------
hc <- myHclust(myDist(m))
clusters <- cutree(hc, k = 4)
tfCluster<-lapply(1:max(clusters), function(i){
  a<-clusters[clusters==i]
  srdmrMotifTable<-loadData2(file.path(CONFIG$dataResult, 'srdmr.homer.motifs.csv'),force.refresh = TRUE)
  darMotifTable<-loadData2(file.path(CONFIG$dataResult, 'atac.dar.homer.motifs.csv'),force.refresh = TRUE)
  clusterTf<-unique(names(a))
  srdmr<-srdmrMotifTable[match(clusterTf,srdmrMotifTable$TFs),2:6]
  dar<-darMotifTable[match(clusterTf,darMotifTable$TFs),3:8]
  cbind(srdmr,dar)
})
names(tfCluster) <- c('cluster4','cluster3','cluster2','cluster1')
plotTFCluster<-function(cluster){
  m<-as.matrix(cluster[,2:ncol(cluster)])
  rownames(m)<-cluster$TFs
  colnames(m)<-c('EarlyHyperDMR', 'EarlyHypoDMR', 'LateHyperDMR','LateHypoDMR',
                 'AISHypoDARs','AISHyperDARs', 'MIAHypoDARs','MIAHyperDARs', 'IACHypoDARs', 'IACHyperDARs')
  Heatmap(m,
          cluster_rows=TRUE,
          cluster_columns = FALSE,
          show_row_names=TRUE,
          show_column_names=TRUE,
          border = TRUE,
          show_heatmap_legend = FALSE,
          col = colorRamp2(breaks = c(0,1.3), colors = c('white', 'black')),
          # top_annotation = column_annotation,
          column_names_gp = grid::gpar(fontsize = 10),
          row_names_gp = grid::gpar(fontsize = 10)
  )
}


saveImage2("tf124.cluster1.heatmap.pdf",width = 4,height = 3.5)
plotTFCluster(tfCluster$cluster3)
dev.off()
saveImage2("tf124.cluster2.heatmap.pdf",width = 3,height = 2.5)
plotTFCluster(tfCluster$cluster2)
dev.off()
saveImage2("tf124.cluster3.heatmap.pdf",width = 3,height = 3)
plotTFCluster(tfCluster$cluster1)
dev.off()
saveImage2("tf124.cluster4.heatmap.pdf",width = 3,height = 12)
plotTFCluster(tfCluster$cluster4)
dev.off()

