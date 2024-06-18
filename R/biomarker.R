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
darDmrOverlap<-loadData2(file.path(CONFIG$dataIntermediate,'marker', 'darDmrOverlap.bed'))
mmDmr<-methyLevel$getMethy(unique(darDmrOverlap$dmr), rmIfcontainNA = TRUE)
mmDarMethy<-atacPeak$getMatchData(darDmrOverlap$dar)$methy
mmDarAccess<-atacPeak$getMatchData(darDmrOverlap$dar)$access
overlapNew<-darDmrOverlap[match(rownames(mmDmr), darDmrOverlap$dmr),]
mmDarMethy<-atacPeak$getMatchData(overlapNew$dar)$methy
mmDarAccess<-atacPeak$getMatchData(overlapNew$dar)$access

calCor<-function(d1,d2,group){
  mmCor<-sapply(1:20, function(i){
    t<-cor.test(d1[,i], d2[,i])
    t$estimate
  })
  mmP<-sapply(1:20, function(i){
    t<-cor.test(d1[,i], d2[,i])
    t$p.value
  })
  data<-data.frame(sample=colnames(d1),group=group,value=mmCor,p=mmP,nlogp=-log10(mmP),q=p.adjust(mmP))
  data
}

df1<-calCor(mmDmr,mmDarMethy,'DmrMethy.vs.DarMethy')
df2<-calCor(mmDarMethy,mmDarAccess,'DarMethy.vs.DarAccess')
df3<-calCor(mmDmr,mmDarAccess,'DmrMethy.vs.DarAccess')
data<-rbind(df1,df2,df3)

colorMap<-c('#bb9727','#f27970','#32b897')
names(colorMap)<-c('DmrMethy.vs.DarMethy','DarMethy.vs.DarAccess','DmrMethy.vs.DarAccess')

yMax<-max(data$nlogp)
b1<-df1$nlogp
b2<-df2$nlogp
b3<-df3$nlogp
names(b1)<-df1$sample
names(b2)<-df2$sample
names(b3)<-df3$sample
saveImage2("marker.dar.dmr.heatmap.3.pdf",width =4,height = 5)
par(mfrow = c(3, 1))
barplot(b1,ylab="-Log10(P-Value)",las=2,cex.names =0.8, col='#bb9727',ylim = c(0, yMax))
barplot(b2,ylab="-Log10(P-Value)",las=2,cex.names =0.8, col='#f27970',ylim = c(0, yMax))
barplot(b3,ylab="-Log10(P-Value)",las=2,cex.names =0.8, col='#32b897',ylim = c(0, yMax))
dev.off()

p1<-ggplot(data=filter(data, group=='DmrMethy.vs.DarMethy'),aes(x=group,y=nlogp,fill=group))+
  scale_fill_manual(values=colorMap) +
  geom_violin( colour ='NA',trim=FALSE)+
  theme_bw()
p2<-ggplot(data=filter(data, group=='DarMethy.vs.DarAccess'),aes(x=group,y=nlogp,fill=group))+
  scale_fill_manual(values=colorMap) +
  geom_violin( colour ='NA',trim=FALSE)+
  theme_bw()
p3<-ggplot(data=filter(data, group=='DmrMethy.vs.DarAccess'),aes(x=group,y=nlogp,fill=group))+
  scale_fill_manual(values=colorMap) +
  geom_violin( colour ='NA',trim=FALSE)+
  theme_bw()
p4<-ggplot(data=filter(data, group=='DarMethy.vs.DarAccess'|group=='DmrMethy.vs.DarAccess'),aes(x=group,y=nlogp,fill=group))+
  scale_fill_manual(values=colorMap) +
  geom_violin( colour = 'NA',trim=FALSE)+
  geom_line(aes(group=sample))+
  theme_bw()
saveImage2("marker.dar.dmr.heatmap.2.pdf",width =4,height = 6.6)
p1/p2/p3/p4
dev.off()

pDmr<-Heatmap(mmDmr,
              col=colorRamp2(c(0, 0.5,1), c("#4574b6", "#fdfec2", "#d83127")),
              name='SRDMR Methylation Level',
              top_annotation = getHeatmapAnnotatio(atacPeak$sample),
              cluster_rows=TRUE,
              cluster_columns = TRUE,
              show_row_names=FALSE,
              show_column_names=FALSE)

pDarMethy<-Heatmap(mmDarMethy,
              col=colorRamp2(c(0, 0.5,1), c("#4574b6", "#fdfec2", "#d83127")),
              name='SRDAR Methylation Level',
              cluster_rows=TRUE,
              cluster_columns = TRUE,
              show_row_names=FALSE,
              show_column_names=FALSE)
pDarAccess<-Heatmap(mmDarAccess,
              name='SRDAR Accessibility',
              col=colorRamp2(c(min(m1.2),quantile(unlist(m1.2), 0.25),quantile(unlist(m1.2), 0.5),quantile(unlist(m1.2), 0.75) ,max(m1.2)), c("#333cac","#1d9cbb","#c5ba59", "#f5cc2f","#f6f803")),
              cluster_rows=TRUE,
              cluster_columns = TRUE,
              show_row_names=FALSE,
              show_column_names=FALSE)
saveImage2("marker.dar.dmr.heatmap.pdf",width =3.5,height = 5)
pDmr%v%pDarMethy%v%pDarAccess
dev.off()
#----------------------------------------------------------------------------------------------------------------------
# Figure 6C. Cluster of patients with DAR methylation level
#----------------------------------------------------------------------------------------------------------------------
darDmrOverlap<-loadData2(file.path(CONFIG$dataIntermediate,'marker', 'darDmrOverlap.bed'))
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
methy<-atacPeak$getAllMethy(darDmrOverlap$dar, sampleName = clinicalData$SampleName)
methy<-t(methy)
coxP<-sapply(1:ncol(methy), function(i){
  df<-data.frame(clinicalData, methy=unlist(methy[,i]))
  fit <- coxph(Surv(DFS, event)~methy,data=df)
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
  col = list(Stage =colorMapStage2,
             Event=c('1'='#f0988c','0'='#1f77b4'),
             Death=c('1'='#f0988c','0'='#1f77b4'),
             Recurrence=c('1'='#f0988c','0'='#1f77b4'),
             Sex=c('female'='#f0988c','male'='#1f77b4')),
  show_annotation_name =TRUE,
  annotation_name_side='right'
)
saveImage2("marker.dardmr.clincal.heatmap.pdf",width = 8,height = 2.5)
# mm<-t(scale(t(mm)))
mm<-na.omit(as.matrix(mm))
Heatmap(mm,
        col=colorRamp2(c(min(mm),median(mm), max(mm)), c("#4574b6", "#fdfec2", "#d83127")),
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
# Figure 6F,G. ROC curve of Random Forest Modeling
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
# Figure S4A,B Random Forest Modeling of TCGA LUAD
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
# Figure 6H. methyltion level of important site of TCGA LUAD select by random forest
#----------------------------------------------------------------------------------------------------------------------
HM450kTCGAforModel<-readRDS(file.path(CONFIG$dataIntermediate,'marker', 'HM450kTCGAforModel.rds'))
hm450kClass2<-loadData2(file.path(CONFIG$dataIntermediate,'marker', 'HM450k.featureImportance.class2.csv'),force.refresh = TRUE)
hm450kClass2<-hm450kClass2[hm450kClass2$inportance>0,]
mm=HM450kTCGAforModel$data[match(hm450kClass2$feature, rownames(HM450kTCGAforModel$data)),]
column_annotation <-HeatmapAnnotation(
  df=data.frame(Stage=ifelse(HM450kTCGAforModel$label==1, 'Tumor', 'Normal')),
  col = list(Stage =c('Tumor'='#c82423', 'Normal'='#2878b5')),
  show_annotation_name =FALSE,
  annotation_name_side='left'
)
mm<-t(scale(t(mm)))
saveImage2("marker.hm450k.rf.feature.heatmap.pdf",width = 3,height = 4)
Heatmap(mm,
        col=colorRamp2(c(-2,0, 2), c("#4574b6", "#fdfec2", "#d83127")),
        top_annotation = column_annotation,
        clustering_method_columns = 'ward.D',
        cluster_rows=TRUE,
        cluster_columns = TRUE,
        show_row_names=FALSE,
        show_column_names=FALSE)
dev.off()
#----------------------------------------------------------------------------------------------------------------------
# Figure 6I. Methylation level in 1416 regions of healthy and lung cancer cfDNA (GSE122126)
#----------------------------------------------------------------------------------------------------------------------
darDmrHM450Overlap<-loadData2(file.path(CONFIG$dataIntermediate,'marker', 'darDmrHM450Overlap.bed'))
myDMP<-readRDS(file.path(CONFIG$dataIntermediate,'open', 'GSE122126_epic.dmp.rds'))
GSE122126<-readRDS(file.path(CONFIG$dataIntermediate,'open', 'GSE122126_epic.beta.nrom.rds'))
label<-c(rep(0,8),rep(1,4))
darDmrHM450OverlapcfDNADMP<-filter(darDmrHM450Overlap, hm450Site%in%rownames(myDMP$N_to_T))
mm<-GSE122126[match(darDmrHM450OverlapcfDNADMP$hm450Site, rownames(GSE122126)),]
column_annotation <-HeatmapAnnotation(
  df=data.frame(Stage=ifelse(label==1, 'Tumor', 'Normal')),
  col = list(Stage =c('Tumor'='#c82423', 'Normal'='#2878b5')),
  show_annotation_name =FALSE,
  annotation_name_side='left'
)
saveImage2("marker.cfDNA.heatmap.pdf",width = 3,height = 4)
mm<-t(scale(t(mm)))
Heatmap(mm,
        col=colorRamp2(c(-2,0,2), c("#4574b6", "#fdfec2", "#d83127")),
        cluster_rows=TRUE,
        top_annotation = column_annotation,
        clustering_method_columns = 'ward.D2',
        clustering_method_rows = 'ward.D2',
        cluster_columns = TRUE,
        show_row_names=FALSE,
        show_column_names=FALSE)
dev.off()



