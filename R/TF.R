source('./R/base.R')
source('./R/local/base.R')
library(dplyr)
library(reshape2)
library(ComplexHeatmap)
library(VennDiagram)
library(circlize)
library(MASS)


#----------------------------------------------------------------------------------------------------------------------
# Figure 3I. Common SRDARs and SRDMRs Homer TFs: EpiTFs
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
SRDAR<-loadSRDAR()
epiTFs<-readRDS(file.path(CONFIG$dataIntermediate,'tf', 'epiTFs.rds'))
allMotif<-loadData2(file.path(CONFIG$dataIntermediate, 'atac','homer.mask','all.motif.region.mask.given.txt'),file.format = 'bed')
allMotifFilter0<-allMotif[allMotif$MotifScore>=5,]
allMotifFilter<-allMotifFilter0[allMotifFilter0$PositionID%in%bed2Feature(SRDAR),]

epiTFsDAR<-left_join(allMotifFilter, epiTFs$map,by=c(`Motif Name`='motif'))%>%
    filter(!is.na(gene))%>%
    dplyr::select(feature=PositionID,tf=gene)%>%
    distinct(feature, tf)
epiTFsAllRegion<-left_join(allMotifFilter0, epiTFs$map,by=c(`Motif Name`='motif'))%>%
  filter(!is.na(gene))%>%
  dplyr::select(feature=PositionID,tf=gene)%>%
  distinct(feature, tf)
SRDARClass<-data.frame(feature=bed2Feature(SRDAR), class=SRDAR$class)
epiTFsDARClass<-left_join(epiTFsDAR,SRDARClass,by = 'feature')
saveRDS(epiTFsDARClass,file.path(CONFIG$dataIntermediate, 'atac','homer.mask','tf.epiTFs.DAR.rds'))
saveRDS(epiTFsAllRegion,file.path(CONFIG$dataIntermediate, 'atac','homer.mask','tf.epiTFs.AllR.rds'))

out<-data.frame(feature2Bed(epiTFsAllRegion$feature), tf=epiTFsAllRegion$tf)
saveTsv(out, file.path(CONFIG$dataIntermediate, 'atac','homer.mask','all.motif.region.mask.given.bed'))
#----------------------------------------------------------------------------------------------------------------------
# Figure 3J. Figure S3B. EpiTFs heatmap
#----------------------------------------------------------------------------------------------------------------------
epiTFsDAR<-readRDS(file.path(CONFIG$dataIntermediate, 'atac','homer.mask','tf.epiTFs.DAR.rds'))
epiTFsAllRegion<-readRDS(file.path(CONFIG$dataIntermediate, 'atac','homer.mask','tf.epiTFs.AllR.rds'))
SRDAR<-loadSRDAR()

getOverlapMatrixLogP<-function(epiTFsDARTarget){
  epiTFsAllRegionList<-split(epiTFsAllRegion,epiTFsAllRegion$tf)
  epiTFsDARList<-split(epiTFsDARTarget,epiTFsDARTarget$tf)
  tfs<-names(epiTFsAllRegionList)
  N0<-length(unique(epiTFsAllRegion$feature))
  N1<-length(unique(epiTFsDARTarget$feature))
  sapply(tfs, function(x){
    sapply(tfs, function(y){
      C0x<-epiTFsAllRegionList[[x]]$feature%>%unique
      C0y<-epiTFsAllRegionList[[y]]$feature%>%unique
      C1x<-epiTFsDARList[[x]]$feature%>%unique
      C1y<-epiTFsDARList[[y]]$feature%>%unique
      NC0xy<-length(intersect(C0x, C0y))
      NC1xy<-length(intersect(C1x, C1y))
      # NC0x<-length(union(C0x,C0y))
      # NC1x<-length(union(C1x,C1y))
      # print(sprintf("%s-%s: %d,%d,%d,%d", x, y,NC0x, NC0xy,  NC1x, NC1xy))
      # test_result <- prop.test(c(NC0xy, NC1xy), c(NC0x, NC1x))
      test_result <- prop.test(c(NC0xy, NC1xy), c(N0, N1))
      test_result$p.value
    })
  })->m
  mm <- -log10(m)
  mm
}
getOverlapMatrixLogFC<-function(epiTFsDARTarget,calMean=FALSE){
  epiTFsAllRegionList<-split(epiTFsAllRegion,epiTFsAllRegion$tf)
  epiTFsDARList<-split(epiTFsDARTarget,epiTFsDARTarget$tf)
  tfs<-names(epiTFsAllRegionList)
  N0<-length(unique(epiTFsAllRegion$feature))
  N1<-length(unique(epiTFsDARTarget$feature))
  sapply(tfs, function(x){
    sapply(tfs, function(y){
      C0x<-epiTFsAllRegionList[[x]]$feature%>%unique
      C0y<-epiTFsAllRegionList[[y]]$feature%>%unique
      C1x<-epiTFsDARList[[x]]$feature%>%unique
      C1y<-epiTFsDARList[[y]]$feature%>%unique
      NC0xy<-length(intersect(C0x, C0y))
      NC1xy<-length(intersect(C1x, C1y))
      (NC1xy/N1)/(NC0xy/N0)
    })
  })->m
  print(min(m))
  print(min(log2(m)))
  if (calMean){
    rowMeans(log2(m))
  }else{
    log2(m)
  }
}

getOverlapMatrixRatio<-function(epiTFsDARTarget){
  epiTFsDARList<-split(epiTFsDARTarget,epiTFsDARTarget$tf)
  tfs<-names(epiTFsDARList)
  N1<-length(unique(epiTFsDARTarget$feature))
  sapply(tfs, function(x){
    sapply(tfs, function(y){
      C1x<-epiTFsDARList[[x]]$feature%>%unique
      C1y<-epiTFsDARList[[y]]$feature%>%unique
      NC1xy<-length(intersect(C1x, C1y))
      NC1xy/N1
    })
  })->m
  m
}

plotTFClusterHeatmapLogP<-function(mm){
  myHclust <- function(x) hclust(x, method = "ward.D2")
  myDist <- function(x) dist(x, method = "euclidean")
  Heatmap(mm,
          name='-log10(PValue)',
          column_names_gp = grid::gpar(fontsize = 6),
          row_names_gp = grid::gpar(fontsize = 6),
          clustering_method_rows = 'ward.D2',
          clustering_method_columns = 'ward.D2',
          clustering_distance_rows = myDist,
          clustering_distance_columns = myDist
  )
}

plotTFClusterHeatmapLogFC<-function(mm){
  myHclust <- function(x) hclust(x, method = "ward.D2")
  myDist <- function(x) dist(x, method = "euclidean")
  Heatmap(mm,
          name='log2FC',
          column_names_gp = grid::gpar(fontsize = 6),
          row_names_gp = grid::gpar(fontsize = 6),
          # col=colorRamp2(c(min(mm), 0, max(mm)), c("green", "white", "red")),
          clustering_method_rows = 'ward.D2',
          clustering_method_columns = 'ward.D2',
          clustering_distance_rows = myDist,
          clustering_distance_columns = myDist
  )
}

plotTFClusterHeatmapRatio<-function(mm){
  myHclust <- function(x) hclust(x, method = "ward.D2")
  myDist <- function(x) dist(x, method = "euclidean")
  row_ha = rowAnnotation(foo2 = runif(10), bar2 = anno_barplot(runif(10)))
  Heatmap(mm,
          name='Ratio',
          column_names_gp = grid::gpar(fontsize = 6),
          row_names_gp = grid::gpar(fontsize = 6),
          right_annotation = row_ha,
          clustering_method_rows = 'ward.D2',
          clustering_method_columns = 'ward.D2',
          clustering_distance_rows = myDist,
          clustering_distance_columns = myDist
  )
}

doCluster<-function(mm,k=6){
  Heatmap(mm,
          name='Ratio',
          column_names_gp = grid::gpar(fontsize = 6),
          row_names_gp = grid::gpar(fontsize = 6),
          clustering_method_rows = 'ward.D2',
          clustering_method_columns = 'ward.D2',
          clustering_distance_rows = myDist,
          clustering_distance_columns = myDist
  )%>%draw()%>%row_dend()->tree
  colors<-c('1'="#E41A1CA0",
            '2'="#377EB8A0",
            '3'="#4DAF4AA0",
            '4'="#984EA3A0",
            '5'="#FF7F00A0",
            '6'="#FFFF33A0")
  list(
    tree=tree,
    clusters= cutree(as.hclust(tree), k=k),
    colors= colors[1:k]
  )
}

plotTFClusterHeatmapRatio2<-function(mm,fc,cluster){
  row_ha = rowAnnotation(cluster = cluster$clusters, 
                         log2FC = anno_barplot(fc,gp=gpar(border =NA,fill="gray6",lty="blank")),
                         col=list(cluster=cluster$colors))
  a<-Heatmap(mm,
             name='Ratio',
             column_names_gp = grid::gpar(fontsize = 6),
             row_names_gp = grid::gpar(fontsize = 6),
             right_annotation = row_ha,
             clustering_method_rows = 'ward.D2',
             clustering_method_columns = 'ward.D2',
             clustering_distance_rows = myDist,
             clustering_distance_columns = myDist
  )%>%draw()
  decorate_row_dend("Ratio", {
    ind = cluster$clusters[order.dendrogram(cluster$tree)]
    ind<-rev(ind)
    x1<-c(0)
    x2<-c(0)
    print(ind)
    for (i in unique(ind)){
      x1<-c(x1,x2[length(x2)])
      x2<-c(x2,x2[length(x2)]+sum(ind==i))
    }
    x1<-x1[-1]
    x2<-x2[-1]
    grid.rect(y = x1/length(ind), height = (x2 - x1)/length(ind),just = "bottom",
              default.units = "npc", gp = gpar(fill = cluster$colors[unique(ind)],col=NA))
  })
}
  

mmLogP<-getOverlapMatrixLogP(epiTFsDAR)
mmLogFC<-getOverlapMatrixLogFC(epiTFsDAR)
mmLogFCMean<-getOverlapMatrixLogFC(epiTFsDAR, calMean = TRUE)
mmRatio<-getOverlapMatrixRatio(epiTFsDAR)

cluster<-doCluster(mmRatio,k=6)

saveImage2("tf.epiTFs.srdar.cluster.heatmap.pvalue.pdf",width = 7,height = 6)
plotTFClusterHeatmapLogP(mmLogP)
dev.off()
saveImage2("tf.epiTFs.srdar.cluster.heatmap.log2fc.pdf",width = 7,height = 6)
plotTFClusterHeatmapLogFC(mmLogFC)#Figure S3B
dev.off()
saveImage2("tf.epiTFs.srdar.cluster.heatmap.ratio.pdf",width = 7,height = 6)
plotTFClusterHeatmapRatio2(mmRatio, mmLogFCMean, cluster)#Figure 3J
dev.off()
#----------------------------------------------------------------------------------------------------------------------
# Figure S3A. EpiTFs cluster
#----------------------------------------------------------------------------------------------------------------------
epiTFsDAR<-readRDS(file.path(CONFIG$dataIntermediate, 'atac','homer.mask','tf.epiTFs.DAR.rds'))
epiTFs<-readRDS(file.path(CONFIG$dataIntermediate,'tf', 'epiTFs.rds'))
clusters<-cluster$clusters
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
