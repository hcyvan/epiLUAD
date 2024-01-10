source('./R/base.R')
source('./R/local/base.R')
library(dplyr)
library(reshape2)
library(ComplexHeatmap)
library(ggExtra)
library(ggplot2)
library(gridExtra)

#----------------------------------------------------------------------------------------------------------------------
# DARs and SC-DMRs Correlation
#----------------------------------------------------------------------------------------------------------------------
atacPeakTPM<-readRDS(file.path(CONFIG$dataIntermediate,'atac', 'atacPeakTPM.rds'))
atacPeakMethyLevel<-readRDS(file.path(CONFIG$dataIntermediate,'atac', 'atacPeakMethyLevel.rds'))
samples<-groups$WGBS.ATAC$select(c('CTL','AIS','MIA','IAC'))

methy<-atacPeakMethyLevel[,samples$SampleName]
keep.methy<-which(rowSums(is.na(methy))==0)
methy<-methy[keep.methy,]
access<-atacPeakTPM[keep.methy,samples$SampleName]

corDiag<-sapply(1:nrow(methy), function(i){
  t<-cor.test(unlist(methy[i,]), unlist(access[i,]), method = 'kendall')
  c(t$estimate, t$p.value)
})
corDiag<-data.frame(t(corDiag))
colnames(corDiag)<-c('tau', 'pvalue')
corDiag<-data.frame(atacPeakTPM[keep.methy,1:3],corDiag)
saveRDS(corDiag, file.path(CONFIG$dataIntermediate, 'atac.wgbs','atac.wgbs.kendall.rds'))
#----------------------------------------------------------------------------------------------------------------------
# Figure 5A. Methylation Levels vs Accessibility in positive, negative and no correlation regions
#----------------------------------------------------------------------------------------------------------------------
corDiag<-readRDS(file.path(CONFIG$dataIntermediate, 'atac.wgbs','atac.wgbs.kendall.rds'))
keep.negative<-which(corDiag$pvalue<0.05&corDiag$tau<0)
keep.positive<-which(corDiag$pvalue<0.05&corDiag$tau>0)
keep.no<-order(abs(corDiag$tau))[1:10000]
cat(sprintf("Negative: %d\n",length(keep.negative)))
cat(sprintf("Positive: %d\n",length(keep.positive)))
plot.methy.access.scatter<-function(keep, title){
  data<-data.frame(x=rowMeans(methy[keep,]), y=rowMeans(access[keep,]))
  p_scatter <- ggplot(data, aes(x = x, y = y)) +
    geom_point(size=0.1,alpha = 0.5) +
    ylab("Average Accessibility")+
    xlab("Average Methylation Level")+
    ggtitle(title)+
    theme_bw()
  p<-ggMarginal(p_scatter, type="density")
  p
}

p1<-plot.methy.access.scatter(keep.negative, 'Negative Correlation')
p2<-plot.methy.access.scatter(keep.positive, 'Positive Correlation')
p3<-plot.methy.access.scatter(keep.no, 'No Correlation')
saveImage2("atac.wgbs.cor.scatter.pdf",width = 5,height = 5)
grid.arrange(p1, p2,p3, nrow = 2)
dev.off()
#----------------------------------------------------------------------------------------------------------------------
# Figure 5B. Correlation across samples
#----------------------------------------------------------------------------------------------------------------------
cal.cor<-function(s){
  x<-methy[,s]
  y<-access[,s]
  set.seed(123)
  a<-sapply(1:5000, function(i){
    keep<-sample(1:length(x), 100)
    x2<-x[keep]
    y2<-y[keep]
    cor(x2,y2, method = 'kendall')
  })
  hist(a,breaks = 100,main=s)
  cat(sprintf('%s: %f\n',s,mean(a)))
  a
}
samples<-groups$WGBS.ATAC$selectBySample(data$variable)
match(data$variable, names(colorMapStage))
corEstList<-lapply(colnames(methy), function(s){
  cal.cor(s)
})
df<-data.frame(do.call(cbind,corEstList))
colnames(df)<-colnames(methy)
data<-melt(df)
data$stage<-samples$Group
data$colors<-samples$colors

saveImage2("atac.wgbs.cor.samples.pdf",width = 5,height = 2.5)
ggplot(data,aes(variable,value))+
  scale_fill_manual(values = colorMapStage)+
  stat_summary(mapping=aes(fill = stage),fun=mean,geom = "bar",fun.args = list(mult=1),width=0.7)+
  stat_summary(fun.data=mean_sdl,fun.args = list(mult=1),geom="errorbar",width=0.2) +
  labs(x = "",y = "Correlation Coefficient") +
  theme_classic()+
  theme(axis.text.x = element_blank())
dev.off()
#----------------------------------------------------------------------------------------------------------------------
# Table S14. DARs and Correlation Regions Overlap
#----------------------------------------------------------------------------------------------------------------------
corDiag<-readRDS(file.path(CONFIG$dataIntermediate, 'atac.wgbs','atac.wgbs.kendall.rds'))
darDeseq2<-readRDS(file.path(CONFIG$dataIntermediate,'atac', 'darDeseq2.rds'))
intersect.cor.dar<-function(diff){
  dar<-feature2Bed(rownames(diff))
  dar<-data.frame(dar, diff)
  corRegion<-filter(corDiag, pvalue<=0.05)
  corDiff<-inner_join(corRegion, dar, by=c('chrom','start','end'))
  a<-data.frame(cor=ifelse(corDiff$tau<0,'negative','positive'), dar=ifelse(corDiff$log2FoldChange>0,'hyper','hypo'))
  Nhyper<-sum(a$cor=='negative'&a$dar=='hyper')
  Nhypo<-sum(a$cor=='negative'&a$dar=='hypo')
  Phyper<-sum(a$cor=='positive'&a$dar=='hyper')
  Phypo<-sum(a$cor=='positive'&a$dar=='hypo')
  NAll<-sum(ifelse(corRegion$tau<0,1,0))
  PAll<-sum(ifelse(corRegion$tau>0,1,0))
  hyperAll<-sum(ifelse(dar$log2FoldChange>0,1,0))
  hypoAll<-sum(ifelse(dar$log2FoldChange<0,1,0))
  out<-t(data.frame(
    negative=c(Nhyper, Nhypo, NAll),
    positive=c(Phyper, Phypo, PAll),
    all=c(hyperAll, hypoAll, 0)
  ))
  colnames(out)<-c('hyper','hypo','all')
  out
}
taAIS<-intersect.cor.dar(darDeseq2$darCTLvsAIS$diff)
taMIA<-intersect.cor.dar(darDeseq2$darCTLvsMIA$diff)
taIAC<-intersect.cor.dar(darDeseq2$darCTLvsIAC$diff)
taAIS
taMIA
taIAC
#----------------------------------------------------------------------------------------------------------------------
# Figure 5C. Methylation Levels and Accessibility Heatmap in DARs
#----------------------------------------------------------------------------------------------------------------------
darDeseq2<-readRDS(file.path(CONFIG$dataIntermediate,'atac', 'darDeseq2.rds'))
atacPeakCount<-readRDS(file.path(CONFIG$dataIntermediate,'atac', 'atacPeakCount.rds'))
atacPeakTPM<-readRDS(file.path(CONFIG$dataIntermediate,'atac', 'atacPeakTPM.rds'))
atacPeakMethyLevel<-readRDS(file.path(CONFIG$dataIntermediate,'atac', 'atacPeakMethyLevel.rds'))
a<-do.call(rbind,lapply(names(darDeseq2), function(x){
  hypo<-subset(darDeseq2[[x]]$diff,log2FoldChange<0)
  hyper<-subset(darDeseq2[[x]]$diff,log2FoldChange>0)
  rbind(
    data.frame(feature=rownames(hypo),group=x, class='Hypo'),
    data.frame(feature=rownames(hyper),group=x, class='Hyper')
  )
}))
darAll<-dcast(a,feature ~ group, value.var='class')
darAll[,2:4][is.na(darAll[,2:4])] <- "NC"
colnames(darAll)[2:4]<-c('AIS', 'IAC', 'MIA')
darAll<-darAll[,c(1,2,4,3)]
dplyr::count(darAll,AIS,MIA,IAC)%>%arrange(desc(n))
bed<-feature2Bed(darAll$feature)
bedPeak<-dplyr::left_join(bed, atacPeakTPM, by = c('chrom'='chrom', 'start'='start', 'end'='end'))
bedMethy<-dplyr::left_join(bed, atacPeakMethyLevel, by = c('chrom'='chrom', 'start'='start', 'end'='end'))
m1<-bedPeak[,4:ncol(bedPeak)]
m2<-bedMethy[,4:ncol(bedMethy)]
m2.keep<-rowSums(is.na(m2))==0

v<-rep(FALSE, length(m2.keep))
set.seed(123)
v[sample(length(m2.keep), 2000)]<-TRUE
m2.keep.2<-m2.keep&v

darAll.2<-darAll[m2.keep.2,]
color.map<-c('red','white','green')
names(color.map)<-c('Hyper','NC','Hypo')
row_annotation <-rowAnnotation(
  show_annotation_name = TRUE,
  annotation_name_side='top',
  df=data.frame(AIS=darAll.2$AIS,MIA=darAll.2$MIA,IAC=darAll.2$IAC),
  col = list(AIS=color.map,MIA=color.map, IAC=color.map)
)
samples<-groups$WGBS.ATAC$selectBySample(colnames(m1))
color.map.col<-groups$WGBS.ATAC$getColorMapVec()
column_annotation <-HeatmapAnnotation(
  show_annotation_name = FALSE,
  df=data.frame(Stage=samples$Group),
  col = list(Stage =color.map.col),
  annotation_name_side='left'
)
m2.2<-m2[m2.keep.2,]
m1.2<-m1[m2.keep.2,]
h1<-Heatmap(m1.2,
            cluster_rows=TRUE,
            cluster_columns = FALSE,
            show_row_names=FALSE,
            show_column_names=FALSE,
            name = "Accessibility",
            column_title = "Accessibility",
            right_annotation  = row_annotation,
            top_annotation = column_annotation
)
h2<-Heatmap(m2.2,
            cluster_rows=TRUE,
            cluster_columns = FALSE,
            show_row_names=FALSE,
            show_column_names=FALSE,
            name = "Methylation Level",
            column_title = "Methylation Level",
            top_annotation = column_annotation
)
saveImage2("atac.wgbs.heatmap.dar.pdf",width = 6,height = 6)
h2+h1
dev.off()
#----------------------------------------------------------------------------------------------------------------------
# Figure 5D. Methylation Levels and Accessibility correlation in epiTFs' binding DARs
#----------------------------------------------------------------------------------------------------------------------
atacPeakTPM<-readRDS(file.path(CONFIG$dataIntermediate,'atac', 'atacPeakTPM.rds'))
atacPeakMethyLevel<-readRDS(file.path(CONFIG$dataIntermediate,'atac', 'atacPeakMethyLevel.rds'))
rnaTPM<-readRDS(file.path(CONFIG$dataIntermediate,'rna', 'rnaTPM.rds'))
rnaTPMSymbol<-mapEnsemble2Symbol(rnaTPM$ensemble)
darDeseq2<-readRDS(file.path(CONFIG$dataIntermediate,'atac', 'darDeseq2.rds'))
data.atacPeakTPM<-groups$WGBS.RNA.ATAC$pickColumnsByGroup(names(colorMapStage), atacPeakTPM)
data.atacPeakMethyLevel<-groups$WGBS.RNA.ATAC$pickColumnsByGroup(names(colorMapStage), atacPeakMethyLevel)
data.rnaTPM<-groups$WGBS.RNA.ATAC$pickColumnsByGroup(names(colorMapStage), rnaTPM)

allMotif<-readRDS(file.path(CONFIG$dataIntermediate, 'atac','homer.mask','tf.epiTFs.DAR.rds'))
epiTfs <- loadData2(file.path(CONFIG$dataIntermediate, 'tf','tf.epTFs.csv'))
motifs<-unique(allMotif$`Motif Name`)
epiTFsMotif<-data.frame(
  motif=motifs,
  tf=sapply(strsplit(motifs, '\\('),function(x){x[1]})
)
epiTFsMotif<-left_join(epiTFsMotif, epiTfs, by = 'tf')
epiTFsMotif<-filter(epiTFsMotif, !is.na(gene))

getCorAcessMethyTpm<-function(tf, gene, motif){
  selectTPM<-unlist(data.rnaTPM[match(gene, rnaTPMSymbol),])
  if(sum(is.na(selectTPM))>0){
    return(NA)
  }
  selectMotifRegion<-intersect(rownames(darDeseq2$darCTLvsIAC$hyper),unique(allMotif[allMotif$`Motif Name`%in%motif,]$PositionID))
  selectMotifIdx<-match(selectMotifRegion,bed2Feature(atacPeakTPM))
  access<-data.atacPeakTPM[selectMotifIdx,]
  methy<-data.atacPeakMethyLevel[selectMotifIdx,]
  rownames(access)<-selectMotifRegion
  rownames(methy)<-selectMotifRegion
  keep.methy<-which(rowSums(is.na(methy))==0)
  methy<-methy[keep.methy,]
  access<-access[keep.methy,]
  
  corAcessMethyTpm<-data.frame(t(sapply(1:nrow(access), function(i){
    x<-unlist(access[i,])
    y<-unlist(methy[i,])
    z<-unlist(selectTPM)
    t.xy<-cor.test(x, y,method='kendall')
    t.xz<-cor.test(x, z,method='kendall')
    t.yz<-cor.test(y, z,method='kendall')
    c(t.xy$estimate,t.xz$estimate,t.yz$estimate)
  })))
  colnames(corAcessMethyTpm)<-c('AccessMethy','AccessTPM','MethyTPM')
  rownames(corAcessMethyTpm)<-rownames(access)
  corAcessMethyTpm$order<-1:nrow(corAcessMethyTpm)
  corAcessMethyTpm<-arrange(corAcessMethyTpm, desc(abs(AccessMethy*MethyTPM*AccessTPM)))
  return(list(
    corAcessMethyTpm=corAcessMethyTpm,
    access=access,
    methy=methy,
    selectTPM=selectTPM
  ))
}

epiTFsCorAMTinHyperDarIAC<-lapply(split(epiTFsMotif, epiTFsMotif$gene), function(x){
  tf<-x[1,2]
  gene<-x[1,3]
  motif<-unlist(x[,1])
  print(tf)
  getCorAcessMethyTpm(tf,gene,motif)
})
epiTFsCorAMTinHyperDarIAC$MEF2B<-NULL
saveRDS(epiTFsCorAMTinHyperDarIAC, file.path(CONFIG$dataIntermediate,'atac', 'epiTFsCorAMTinHyperDarIAC.rds'))

plotCor<-function(key,color, label){
  accessMethCor<-do.call(rbind,lapply(names(epiTFsCorAMTinHyperDarIAC), function(x){
    tfAMY<-epiTFsCorAMTinHyperDarIAC[[x]]
    data.frame(
      cor=tfAMY$corAcessMethyTpm[[key]],
      tf=x
    )
  }))
  data<-accessMethCor
  ggplot(data=data,aes(x=tf,y=cor))+
    # geom_bar()+
    stat_summary(fun=mean,geom = "bar",fun.args = list(mult=1),fill = color)+
    stat_summary(fun.data=mean_sdl,fun.args = list(mult=1),geom="errorbar",width=0.2) +
    labs(x="",y=label)+
    scale_y_continuous(labels = function(x) format(x, nsmall = 2))+
    theme_classic()+
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
      axis.title.x = element_text(size=15,angle = 90, vjust = 0.5, hjust=1),
      axis.title.y = element_text(size=15),
      axis.text = element_text(size = 12,colour="black"))+
    guides(colour = guide_legend(override.aes = list(shape = 12,size=10)))
}

saveImage2("atac.wgbs.epiTFs.CorAcessMethyTpm.pdf",width = 16,height = 6)
p1<-plotCor('AccessMethy','green3', 'Correlation')
p2<-plotCor('AccessTPM','red3', 'Correlation ')
p3<-plotCor('MethyTPM','orange3', 'Correlation ')
grid.arrange(p1, p2,p3, nrow = 3)
dev.off()

colors <- c('green3','red3','orange3')
type <- c('Accessibility vs. Methylation','Accessibility vs. TPM','Methylation vs. TPM')
saveImage2("atac.wgbs.epiTFs.CorAcessMethyTpm.legend.pdf",width = 16,height = 7)
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend = type, pch=15, pt.cex=3, cex=1.5, bty='n',col = alpha(colors),horiz=FALSE)
dev.off()

# 
# TFinfo<-amt$ASCL1
# par(mfrow = c(9, 3))
# par(mar = c(2, 4, 1, 1))
# for(i in 1:9){
#   corAcessMethyTpm<-TFinfo$corAcessMethyTpm
#   access<-TFinfo$access
#   methy<-TFinfo$methy
#   selectTPM<-TFinfo$selectTPM
#   
#   idx<-corAcessMethyTpm[i,]$order
#   x<-unlist(access[idx,])
#   y<-unlist(methy[idx,])
#   z<-unlist(selectTPM)
#   plot(x,y, xlab="Accessibility", ylab='Methylation Level',pch=19,col=(groups$WGBS.RNA.ATAC$select(names(colorMapStage)))$colors)
#   abline(lm(y ~ x))
#   plot(z,x, ylab="Accessibility", xlab='Motif TPM',pch=19,col=(groups$WGBS.RNA.ATAC$select(names(colorMapStage)))$colors)
#   abline(lm(x ~ z))
#   plot(z,y, ylab="Methylation Level", xlab='Motif TPM',pch=19,col=(groups$WGBS.RNA.ATAC$select(names(colorMapStage)))$colors)
#   abline(lm(y ~ z))
# }
#----------------------------------------------------------------------------------------------------------------------
# DARs and SC-DMRs Correlation
#----------------------------------------------------------------------------------------------------------------------

       