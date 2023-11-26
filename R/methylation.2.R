source('./R/base.R')
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(GenomicRanges)
library(rtracklayer)

#----------------------------------------------------------------------------------------------------------------------
# The Average DNA methylation level of some Genomic Regions
#----------------------------------------------------------------------------------------------------------------------
genomicRegionMethyLevel<-readRDS(file.path(CONFIG$dataExternal, 'genomicRegionMethyLevel.rds'))
m<-do.call(rbind,lapply(genomicRegionMethyLevel, function(x){
  x<-x[,4:ncol(x)]
  colMeans(x,na.rm=TRUE)
}))
rownames(m)<-c('CpG Islands', 'CpG Sea','CpG Shelves', 'CpG Shores', 'Exons', 'Intergenic', 'Intron', 'Promoter 1K', 'Promoter 5k', 'TSS', "3'UTR", "5'UTR")
samples<-groups$WGBS$selectBySample(colnames(m))

color.map<-groups$WGBS$getColorMapVec()
column_annotation <-HeatmapAnnotation(
  df=data.frame(Group=samples$Group),
  col = list(Group =color.map),
  show_annotation_name =FALSE,
  annotation_name_side='left'
)
saveImage2("genomicRegion.heatmap.pdf",width = 16,height = 2.5)
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
saveImage2("genomicRegion.scale.heatmap.pdf",width = 16,height = 2.5)
Heatmap(t(scale(t(m))),
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
#----------------------------------------------------------------------------------------------------------------------
genomicRegionStatistic<-do.call(rbind,lapply(1:nrow(m),function(i){
  x<-m[i,]
  g.ctl<-match(groups$WGBS$select('CTL')$SampleName,colnames(m))
  g.ais<-match(groups$WGBS$select('AIS')$SampleName,colnames(m))
  g.mia<-match(groups$WGBS$select('MIA')$SampleName,colnames(m))
  g.iac<-match(groups$WGBS$select('IAC')$SampleName,colnames(m))
  
  level.mean<-c(mean(x[g.ctl]),mean(x[g.ais]),mean(x[g.mia]),mean(x[g.iac]))
  level.sd<-c(sd(x[g.ctl]),sd(x[g.ais]),sd(x[g.mia]),sd(x[g.iac]))
  #Coefficient of Variation
  level.cv<-level.sd/level.mean
  #Interquartile Range, IQR
  level.iqr <- c(IQR(x[g.ctl]),IQR(x[g.ais]),IQR(x[g.mia]),IQR(x[g.iac]))
  data.frame(GenomicRange=rownames(m)[i],Group=names(color.map), Mean=level.mean,SD=level.sd, CV=level.cv,IQR=level.iqr)
}))
write.csv(genomicRegionStatistic, file.path(CONFIG$dataResult, 'genomicRegionStatistic.csv'),row.names  = FALSE,quote = FALSE)


#----------------------------------------------------------------------------------------------------------------------
# Table S4. Hyper-DMCs and Hypo-DMCs of different group in LUAD				
#----------------------------------------------------------------------------------------------------------------------
dmcCTLvsAIS<-loadData2(file.path(CONFIG$dataExternal,'dmc', 'dmc.CTL.vs.AIS.txt'),file.format='bed')
dmcCTLvsMIA<-loadData2(file.path(CONFIG$dataExternal,'dmc', 'dmc.CTL.vs.MIA.txt'),file.format='bed')
dmcCTLvsIAC<-loadData2(file.path(CONFIG$dataExternal,'dmc', 'dmc.CTL.vs.IAC.txt'),file.format='bed')
sum(table(dmcCTLvsAIS$class))
sum(table(dmcCTLvsMIA$class))
sum(table(dmcCTLvsIAC$class))
CTLvsAIS<-data.frame(group="CTL vs. AIS", class=dmcCTLvsAIS$class)
CTLvsMIA<-data.frame(group="CTL vs. MIA", class=dmcCTLvsMIA$class)
CTLvsIAC<-data.frame(group="CTL vs. IAC", class=dmcCTLvsIAC$class)
data<-Reduce(rbind, list(CTLvsAIS,CTLvsMIA,CTLvsIAC))
target<-table(data)
test <- chisq.test(target)
print(target)
print(test$residuals)
print(test)

#----------------------------------------------------------------------------------------------------------------------
# Table S5. Hyper-DMCs and Hypo-DMCs in different genomic regions
#----------------------------------------------------------------------------------------------------------------------
genomicRegionMethyLevel<-readRDS(file.path(CONFIG$dataExternal, 'genomicRegionMethyLevel.rds'))
dmcCTLvsAIS<-loadData2(file.path(CONFIG$dataExternal,'dmc', 'dmc.CTL.vs.AIS.txt'),file.format='bed')
dmcCTLvsMIA<-loadData2(file.path(CONFIG$dataExternal,'dmc', 'dmc.CTL.vs.MIA.txt'),file.format='bed')
dmcCTLvsIAC<-loadData2(file.path(CONFIG$dataExternal,'dmc', 'dmc.CTL.vs.IAC.txt'),file.format='bed')
countDmcInGonomicRegion<-function(dmc) {
  hyper<-bed2GRanges(filter(dmc, class=='strongHyper'))
  hypo<-bed2GRanges(filter(dmc, class=='strongHypo'))
  count<-sapply(names(genomicRegionMethyLevel),function(x){
    bed<-genomicRegionMethyLevel[[x]]
    gr1<-bed2GRanges(bed)
    # rbind(sum(countOverlaps(gr1, hyper)), sum(countOverlaps(gr1, hypo)), sum(width(gr1))/(1000*1000))
    rbind(sum(countOverlaps(gr1, hyper)), sum(countOverlaps(gr1, hypo)))
  })
  # rownames(count) <- c('hyper', 'hypo', 'regionSizeMb')
  rownames(count) <- c('hyper', 'hypo')
  data.frame(count)
}
statistic<-list(
  AIS=countDmcInGonomicRegion(dmcCTLvsAIS),
  MIA=countDmcInGonomicRegion(dmcCTLvsMIA),
  AIC=countDmcInGonomicRegion(dmcCTLvsIAC)
)

genomicRegionDmc<-do.call(rbind,lapply(names(genomicRegionMethyLevel), function(region){
  data<-t(do.call(cbind,lapply(names(statistic), function(group){
    item<-statistic[[group]]
    item[region]
  })))
  rownames(data)<-paste(c('CTL.vs.AIS', 'CTL.vs.MIA', 'CTL.vs.IAC'),region)
  test <- chisq.test(data)
  out<-data.frame(cbind(data, test$residuals))
  out$region <- region
  out
  print(test)
}))
genomicRegionDmc
write.csv(genomicRegionDmc, file.path(CONFIG$dataResult, 'genomicRegionDmc.csv'),quote = FALSE)

