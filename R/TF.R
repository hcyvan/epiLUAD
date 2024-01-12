source('./R/base.R')
source('./R/local/base.R')
library(dplyr)
library(reshape2)
library(ComplexHeatmap)
library(VennDiagram)
library(circlize)


#----------------------------------------------------------------------------------------------------------------------
# Figure 3I. Common DARs and SCDMRs Homer TFs: EpiTFs
#----------------------------------------------------------------------------------------------------------------------
srdmrTFs<-readRDS(file.path(CONFIG$dataIntermediate,'wgbs', 'srdmr.tfs.rds'))
srdarTFs<-readRDS(file.path(CONFIG$dataIntermediate,'atac', 'srdar.tfs.rds'))
venn.plot<-venn.diagram(x = list(DAR = srdarTFs$tf, `SRDMR` = srdmrTFs$tf),
                        fill = c("red","blue"),
                        cat.cex = 2,
                        cex = 2,
                        filename = NULL)
saveImage2("tf.epiTFs.venn.pdf",width = 4,height = 4)
grid.draw(venn.plot)
dev.off()

tfs<-intersect(srdmrTFs$tf, srdarTFs$tf)
epiTFs<-list(
  tf=tfs,
  stageDAR=srdarTFs$stage[match(tfs, rownames(srdarTFs$stage)),],
  stageDMR=srdmrTFs$stage[match(tfs, rownames(srdmrTFs$stage)),],
  map=rbind(srdarTFs$map,srdmrTFs$map)%>%distinct(gene,motif,.keep_all = TRUE)%>%filter(gene%in%tfs)
)
saveRDS(epiTFs, file.path(CONFIG$dataIntermediate,'tf', 'epiTFs.rds'))
#----------------------------------------------------------------------------------------------------------------------
# Common DARs and SCDMRs Homer TFs: epiTFs 
#----------------------------------------------------------------------------------------------------------------------
darDeseq2<-readRDS(file.path(CONFIG$dataIntermediate,'atac', 'darDeseq2.rds'))
diffPeak<-lapply(darDeseq2, function(x){
  rownames(x$diff)
})%>%unlist()%>%unique()
epiTFs<-readRDS(file.path(CONFIG$dataIntermediate,'tf', 'epiTFs.rds'))
allMotif<-loadData2(file.path(CONFIG$dataIntermediate, 'atac','homer.mask','all.motif.region.txt'),file.format = 'bed')
allMotifFilter0<-allMotif[allMotif$MotifScore>=5,]
allMotifFilter<-allMotifFilter0[allMotifFilter0$PositionID%in%diffPeak,]

epiTFsDAR<-left_join(allMotifFilter, epiTFs$map,by=c(`Motif Name`='motif'))%>%
    filter(!is.na(gene))%>%
    dplyr::select(feature=PositionID,tf=gene)%>%
    distinct(feature, tf)
epiTFsAllRegion<-left_join(allMotifFilter0, epiTFs$map,by=c(`Motif Name`='motif'))%>%
  filter(!is.na(gene))%>%
  dplyr::select(feature=PositionID,tf=gene)%>%
  distinct(feature, tf)

saveRDS(epiTFsDAR,file.path(CONFIG$dataIntermediate, 'atac','homer.mask','tf.epiTFs.DAR.rds'))
saveRDS(epiTFsAllRegion,file.path(CONFIG$dataIntermediate, 'atac','homer.mask','tf.epiTFs.AllR.rds'))
#----------------------------------------------------------------------------------------------------------------------
# Figure 3H. Figure S3B. EpiTFs heatmap
#----------------------------------------------------------------------------------------------------------------------
epiTFsDAR<-readRDS(file.path(CONFIG$dataIntermediate, 'atac','homer.mask','tf.epiTFs.DAR.rds'))
epiTFsAllRegion<-readRDS(file.path(CONFIG$dataIntermediate, 'atac','homer.mask','tf.epiTFs.AllR.rds'))

getOverlapMatrix<-function(tfRegion){
  tfRegion<-split(tfRegion, tfRegion$tf)
  # total<-length(unique(allMotifFilter0$PositionID))
  total<-144452
  tfOverlap<-sapply(tfRegion, function(x){
    sapply(tfRegion, function(y){
      xy<-intersect(unique(x$feature),unique(y$feature))
      k<-length(xy)
      n<-nrow(x)
      M<-nrow(y)
      N<-total-M
      # print(sprintf("%d,%d,%d,%d", k, M, N,n))
      phyper(k, M, N, n, lower.tail = FALSE)
    })
  })
  m <- -log10(tfOverlap+1e-322)
}
plotTFClusterHeatmap<-function(m){
  myHclust <- function(x) hclust(x, method = "ward.D2")
  myDist <- function(x) dist(x, method = "euclidean")
  Heatmap(m,
          name='-log10(PValue)',
          column_names_gp = grid::gpar(fontsize = 6),
          row_names_gp = grid::gpar(fontsize = 6),
          clustering_method_rows = 'ward.D2',
          clustering_method_columns = 'ward.D2',
          clustering_distance_rows = myDist,
          clustering_distance_columns = myDist
  )
}

overlapMatrix<-getOverlapMatrix(epiTFsDAR)
saveImage2("tf.epiTFs.dar.heatmap.cluster.pdf",width = 7,height = 6)
plotTFClusterHeatmap(overlapMatrix)
dev.off()

overlapMatrix2<-getOverlapMatrix(epiTFsAllRegion)
saveImage2("tf.epiTFs.allRegion.heatmap.cluster.pdf",width = 7,height = 6)
plotTFClusterHeatmap(overlapMatrix2)
dev.off()

#----------------------------------------------------------------------------------------------------------------------
# Figure S3A. EpiTFs cluster
#----------------------------------------------------------------------------------------------------------------------
epiTFsDAR<-readRDS(file.path(CONFIG$dataIntermediate, 'atac','homer.mask','tf.epiTFs.DAR.rds'))
epiTFs<-readRDS(file.path(CONFIG$dataIntermediate,'tf', 'epiTFs.rds'))
overlapMatrix<-getOverlapMatrix(epiTFsDAR)

myDist <- function(x) dist(x, method = "euclidean")
myHclust <- function(x) hclust(x, method = "ward.D2")
hc <- myHclust(myDist(overlapMatrix))
clusters <- cutree(hc, k = 6)
sort(clusters)

tfcluster<-lapply(unique(clusters), function(i){
  tfs<-names(clusters[clusters==i])
  dar<-epiTFs$stageDAR[match(tfs, rownames(epiTFs$stageDAR)),]
  dmr<-epiTFs$stageDMR[match(tfs, rownames(epiTFs$stageDMR)),]
  cbind(dmr,dar)
})

plotCluster<-function(i){
  m<-as.matrix(tfcluster[[i]])
  print(dim(m))
  column_split <- factor(c(rep("SRDMRs", 6), rep("SRDARs", 6)))
  column_annotation <- HeatmapAnnotation(df = data.frame(
    ColAnno = c(rep("SRDMRs", 6), rep("SRDARs", 6))),
    col = list(ColAnno = c("SRDMRs" = "blue", "SRDARs" = "red")),
    which = "col")
  Heatmap(m,
             cluster_rows=TRUE,
             cluster_columns = FALSE,
             show_row_names=TRUE,
             show_column_names=TRUE,
             border = TRUE,
             show_heatmap_legend = FALSE,
             col = colorRamp2(breaks = c(0,1.3), colors = c('white', 'gray1')),
             column_names_gp = grid::gpar(fontsize = 10),
             row_names_gp = grid::gpar(fontsize = 10),
             cell_fun = function(j, i, x, y, width, height, fill) {
               grid::grid.rect(x, y, width, height, 
                               gp = grid::gpar(fill = fill, col = 'black'))
             },
          column_split = column_split
  )
}

c1<-plotCluster(1)
c2<-plotCluster(2)
c3<-plotCluster(3)
c4<-plotCluster(4)
c5<-plotCluster(5)
c6<-plotCluster(6)
saveImage2("tf.epiTFs.dar.cluster.pdf",width = 3,height = 11)
c1%v%c2%v%c3%v%c4%v%c5%v%c6
dev.off()
