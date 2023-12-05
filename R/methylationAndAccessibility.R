source('./R/base.R')
source('./R/local/base.R')
library(dplyr)
library(reshape2)
library(ComplexHeatmap)
library(ggExtra)
library(ggplot2)
library(gridExtra)


#----------------------------------------------------------------------------------------------------------------------
# Common DARs and SC-DMRs Homer TFs
#----------------------------------------------------------------------------------------------------------------------
motifDAR<-read.csv(file.path(CONFIG$dataResult, 'atac.dar.homer.motifs.csv'))
motifSRDMR<-read.csv(file.path(CONFIG$dataResult, 'srdmr.homer.motifs.csv'))
darHomerTFs<-unique(motifDAR$TSs)
srdmrHomerTFs<-unique(motifSRDMR$TSs)
cat(sprintf("Common TFs: %d", length((intersect(darHomerTFs,srdmrHomerTFs)))))


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


corDiag<-readRDS(file.path(CONFIG$dataIntermediate, 'atac.wgbs','atac.wgbs.kendall.rds'))
keep.negative<-which(corDiag$pvalue<0.05&corDiag$tau<0)
keep.positive<-which(corDiag$pvalue<0.05&corDiag$tau>0)
keep.no<-order(abs(corDiag$tau))[1:10000]

cat(sprintf("Negative: %d\n",length(keep.negative)))
cat(sprintf("Positive: %d\n",length(keep.positive)))
#----------------------------------------------------------------------------------------------------------------------
# Methylation Levels vs Accessibility in positive, negative and no correlation regions
#----------------------------------------------------------------------------------------------------------------------
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
# Correlation across samples
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
# DARs and Correlation Regions Overlap
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
# DARs and SC-DMRs Correlation
#----------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------
# DARs and SC-DMRs Correlation
#----------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------
# DARs and SC-DMRs Correlation
#----------------------------------------------------------------------------------------------------------------------

       