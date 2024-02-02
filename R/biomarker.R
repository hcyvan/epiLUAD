source('./R/base.R')
source('./R/local/base.R')
library(data.table)
library(ComplexHeatmap)
library(VennDiagram)
library(survminer)
library(survival)
library(glmnet)


#----------------------------------------------------------------------------------------------------------------------
# Figure 6A. The SRDARs overlap with SRDMRs
#----------------------------------------------------------------------------------------------------------------------
SRDMR<-loadSRDMR()
SRDAR<-loadSRDAR()
srdmr<-bed2GRanges(SRDMR)
srdar<-bed2GRanges(SRDAR)
dist<-distanceToNearest(srdar, srdmr)
d<-dist@elementMetadata@listData$distance
dist2<-dist[which(d<=0)]
dar<-srdar[dist2@from]
dmr<-srdmr[dist2@to]
darDmrOverlap<-data.frame(dar=GRanges2Feature(dar), dmr=GRanges2Feature(dmr))
saveTsv(darDmrOverlap, file.path(CONFIG$dataIntermediate,'panel', 'darDmrOverlap.bed'))
overlapNum<-length(unique(darDmrOverlap$dar))
darNum<-length(srdar)
dmrNum<-length(srdmr)
cat(sprintf("SRDMR: %s; SRDAR: %s; Overlap: %s", dmrNum, darNum, overlapNum))
#----------------------------------------------------------------------------------------------------------------------
# Figure 6B. The SRDARs overlap with SRDMRs
#----------------------------------------------------------------------------------------------------------------------
atacPeak<-AtacPeak('WGBS.ATAC')
methyLevel<-MethyLevel('WGBS.ATAC')
darDmrOverlap<-loadData2(file.path(CONFIG$dataIntermediate,'panel', 'darDmrOverlap.bed'))
mmDmr<-methyLevel$getMethy(darDmrOverlap$dmr, rmIfcontainNA = TRUE)
mmDarMethy<-atacPeak$getMatchData(darDmrOverlap$dar)$methy
pDmr<-Heatmap(mmDmr,
              name='SRDMR',
              top_annotation = getHeatmapAnnotatio(atacPeak$sample),
              cluster_rows=TRUE,
              cluster_columns = TRUE,
              show_row_names=FALSE,
              show_column_names=FALSE)


pDar<-Heatmap(mmDarMethy,
              name='SRDAR',
              cluster_rows=TRUE,
              cluster_columns = TRUE,
              show_row_names=FALSE,
              show_column_names=FALSE)
saveImage2("marker.dar.dmr.heatmap.pdf",width = 5,height = 3)
pDmr%v%pDar
dev.off()
#----------------------------------------------------------------------------------------------------------------------
# Figure 6C. Cluster of patients with DAR methylation level
#----------------------------------------------------------------------------------------------------------------------
darDmrOverlap<-loadData2(file.path(CONFIG$dataIntermediate,'panel', 'darDmrOverlap.bed'))
atacPeak<-AtacPeak('WGBS')
darMethy<-atacPeak$getAllMethy(unique(darDmrOverlap$dar))
saveImage2("marker.darMethy.cluster.pdf",width = 4,height = 4)
plotCircleCluster(darMethy, groups$WGBS$selectBySample(colnames(darMethy))$color)
dev.off()
#----------------------------------------------------------------------------------------------------------------------
# Figure 6D,E DAR methylation level and surivival
#----------------------------------------------------------------------------------------------------------------------
clinicalData<-readRDS(file.path(CONFIG$dataIntermediate,'marker', 'clinicalData.rds'))
atacPeak<-AtacPeak('WGBS')
darDmrOverlap<-loadData2(file.path(CONFIG$dataIntermediate,'marker', 'darDmrOverlap.bed'))
methy<-atacPeak$getAllMethy(darDmrOverlap$dar, sampleName = data$SampleName)
methy<-t(methy)
coxP<-sapply(1:ncol(methy), function(i){
  df<-data.frame(data, methy=unlist(methy[,i]))
  fit <- coxph(Surv(DFS, progression)~methy,data=df)
  p<-summary(fit)$coefficients[5]
  p
})

peakDmrSur<-data.frame(darDmrOverlap, coxP=coxP)
peakDmrSur<-filter(peakDmrSur, coxP<0.01)

mm<-atacPeak$getAllMethy(peakDmrSur$dar, sampleName = clinicalData$SampleName)
column_annotation <-HeatmapAnnotation(
  df=data.frame(Stage=groups$WGBS$selectBySample(clinicalData$SampleName)$Group,
                Event=clinicalData$event,
                Death=clinicalData$dead,
                Recurrence=clinicalData$recurrence,
                Age=clinicalData$Age,
                Sex=clinicalData$Sex),
  col = list(Stage =colorMapStage,
             Event=c('1'='red3','0'='green3'),
             Death=c('1'='red3','0'='green3'),
             Recurrence=c('1'='red3','0'='green3'),
             Sex=c('female'='red','male'='blue')),
  show_annotation_name =TRUE,
  annotation_name_side='right'
)
saveImage2("marker.dardmr.clincal.heatmap.pdf",width = 15,height = 3)
Heatmap(t(scale(t(mm))),
        cluster_rows=TRUE,
        name="Scaled Methylation Level",
        top_annotation = column_annotation,
        clustering_method_columns = 'ward.D2',
        clustering_method_rows = 'ward.D2',
        cluster_columns = TRUE,
        show_row_names=FALSE,
        show_column_names=FALSE)
dev.off()

df<-data.frame(clinicalData, methy=colMeans(mm,na.rm = TRUE))
df$group<-ifelse(df$methy>median(df$methy), 'High', 'Low')
fit <- survfit(Surv(DFS, event) ~group, data = df)
saveImage2("marker.dardmr.clincal.survplot.pdf",width = 3.5,height = 3.5)
ggsurvplot(fit,pval=T)
dev.off()
#----------------------------------------------------------------------------------------------------------------------
# Figure 6F,G. 
#----------------------------------------------------------------------------------------------------------------------
atacPeak<-AtacPeak('WGBS')
darDmrOverlap<-loadData2(file.path(CONFIG$dataIntermediate,'marker', 'darDmrOverlap.bed'))
wgbs<-atacPeak$getAllMethy(unique(darDmrOverlap$dar))
saveTsv(wgbs, file.path(CONFIG$dataIntermediate,'marker', 'wgbs.txt'),col.names = FALSE)
group<-groups$WGBS$selectBySample(colnames(wgbs))$Group
saveTsv(as.numeric(group)-1, file.path(CONFIG$dataIntermediate,'marker', 'wgbs.group.txt'),col.names = FALSE)
saveTsv(rownames(wgbs), file.path(CONFIG$dataIntermediate,'marker', 'wgbs.feature.txt'),col.names = FALSE)
# Random forest modeling using python sklearn
# ./doc/marker.ipynb
#----------------------------------------------------------------------------------------------------------------------
# Figure S4A,B
#----------------------------------------------------------------------------------------------------------------------
dataClinical<-loadData2(file.path(CONFIG$dataIntermediate,'marker', 'tcgaLUAD450kClinical.tsv'),force.refresh = TRUE)
atacPeak<-AtacPeak('WGBS')
tcgaLUAD450k<-readRDS(file.path(CONFIG$dataIntermediate,'marker', 'tcgaLUAD450k.rds'))
darDmrHM450Overlap<-loadData2(file.path(CONFIG$dataIntermediate,'marker', 'darDmrHM450Overlap.bed'),force.refresh = TRUE)
dataTCGA<-tcgaLUAD450k[match(intersect(darDmrHM450Overlap$hm450, rownames(tcgaLUAD450k)),rownames(tcgaLUAD450k)),]
label<-dataClinical$label[match(colnames(dataTCGA), dataClinical$id)]
saveRDS(list(data=dataTCGA,label=label), file.path(CONFIG$dataIntermediate,'marker', 'HM450kTCGAforModel.rds'))
saveTsv(dataTCGA, file.path(CONFIG$dataIntermediate,'marker', 'HM450k.txt'),col.names = FALSE)
saveTsv(label, file.path(CONFIG$dataIntermediate,'marker', 'HM450k.group.2class.txt'),col.names = FALSE)
saveTsv(rownames(dataTCGA), file.path(CONFIG$dataIntermediate,'marker', 'HM450k.feature.txt'),col.names = FALSE)

saveImage2("marker.tcgaLUAD450kClinical.tcga.cluster.pdf",width = 8,height = 8)
plotCircleCluster(dataTCGA, ifelse(label==1, 'red', 'green'))
dev.off()
# Random forest modeling using python sklearn
# ./doc/marker.ipynb
#----------------------------------------------------------------------------------------------------------------------
# Figure 6H.
#----------------------------------------------------------------------------------------------------------------------
HM450kTCGAforModel<-readRDS(file.path(CONFIG$dataIntermediate,'marker', 'HM450kTCGAforModel.rds'))
hm450kClass2<-loadData2(file.path(CONFIG$dataIntermediate,'marker', 'HM450k.featureImportance.class2.csv'),force.refresh = TRUE)
hm450kClass2<-hm450kClass2[hm450kClass2$inportance>0,]
match(hm450kClass2$feature, dataTCGAFeature)
match(hm450kClass2$feature, rownames(HM450kTCGAforModel$data))
mm=HM450kTCGAforModel$data[match(hm450kClass2$feature, rownames(HM450kTCGAforModel$data)),]
column_annotation <-HeatmapAnnotation(
  df=data.frame(Stage=ifelse(HM450kTCGAforModel$label==1, 'Tumor', 'Normal')),
  col = list(Stage =c('Tumor'='red3', 'Normal'='green3')),
  show_annotation_name =FALSE,
  annotation_name_side='left'
)
saveImage2("marker.hm450k.rf.feature.heatmap.pdf",width = 10,height = 4)
Heatmap(t(scale(t(mm))),
        name='Scaled Mehtylation Level',
        top_annotation = column_annotation,
        clustering_method_columns = 'ward.D',
        # clustering_method_rows = 'ward.D2',
        cluster_rows=TRUE,
        cluster_columns = TRUE,
        show_row_names=FALSE,
        show_column_names=FALSE)
dev.off()
#----------------------------------------------------------------------------------------------------------------------
# Figure 6I.
#----------------------------------------------------------------------------------------------------------------------
darDmrHM450Overlap<-loadData2(file.path(CONFIG$dataIntermediate,'marker', 'darDmrHM450Overlap.bed'))
myDMP<-readRDS(file.path(CONFIG$dataIntermediate,'open', 'GSE122126_epic.dmp.rds'))
GSE122126<-readRDS(file.path(CONFIG$dataIntermediate,'open', 'GSE122126_epic.beta.nrom.rds'))
label<-c(rep(0,8),rep(1,4))
darDmrHM450OverlapcfDNADMP<-filter(darDmrHM450Overlap, hm450Site%in%rownames(myDMP$N_to_T))
mm<-GSE122126[match(darDmrHM450OverlapcfDNADMP$hm450Site, rownames(GSE122126)),]
column_annotation <-HeatmapAnnotation(
  df=data.frame(Stage=ifelse(label==1, 'Tumor', 'Normal')),
  col = list(Stage =c('Tumor'='red3', 'Normal'='green3')),
  show_annotation_name =FALSE,
  annotation_name_side='left'
)
saveImage2("marker.cfDNA.heatmap.pdf",width = 3,height = 4)
Heatmap(t(scale(t(mm))),
        cluster_rows=TRUE,
        top_annotation = column_annotation,
        clustering_method_columns = 'ward.D2',
        clustering_method_rows = 'ward.D2',
        cluster_columns = TRUE,
        show_row_names=FALSE,
        show_column_names=FALSE)
dev.off()



