source('./R/base.R')
library(ComplexHeatmap)
library(circlize)
library(ggplot2)


#----------------------------------------------------------------------------------------------------------------------
# The Average DNA methylation level of some Genomic Regions
#----------------------------------------------------------------------------------------------------------------------
genomicRegionMethyLevel<-readRDS(file.path(CONFIG$dataExternal, 'genomicRegionMethyLevel.rds'))
m<-do.call(rbind,lapply(genomicRegionMethyLevel, function(x){
  x<-x[,4:ncol(x)]
  colMeans(x,na.rm=TRUE)
}))
rownames(m)<-c('CpG Islands', 'CpG Shelves', 'CpG Shores', 'Exons', 'Intergenic', 'Intron', 'Promoter 1K', 'Promoter 5k', 'TSS', "3'UTR", "5'UTR")
samples<-groups$WGBS$selectBySample(colnames(m))

color.map<-groups$WGBS$getColorMapVec()
column_annotation <-HeatmapAnnotation(
  df=data.frame(Group=samples$Group),
  col = list(Group =color.map),
  show_annotation_name =FALSE,
  annotation_name_side='left'
)
saveImage2("genomeRegion.heatmap.pdf",width = 16,height = 2.5)
Heatmap(m,
        cluster_rows=TRUE,
        cluster_columns = FALSE,
        show_column_names=FALSE,
        bottom_annotation = column_annotation,
        # column_names_gp = gpar(col = samples$colors)
        heatmap_legend_param = list(
          title = "DNA Methylation Levels",
          legend_height = unit(4, "cm"),
          at = c(0,0.5,1),
          labels = c('0','0.5','1'),
          title_position = "lefttop-rot"
        ),
)
dev.off()
