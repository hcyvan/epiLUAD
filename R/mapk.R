source('./R/base.R')
source('./R/local/base.R')
library(dplyr)
library(reshape2)
library(ComplexHeatmap)
library(stringr)
library(gridExtra)
library(RCircos)


plotFeatureRCircos<-function(features){
  atacPeak<-AtacPeak()
  SRDEG<-loadSRDEG()
  srdeg<-Srdeg()
  anno<-annoPeak(bed2GRanges(feature2Bed(features)))
  diff<-anno[anno$symbol%in%srdeg$all()]
  genes<-diff$symbol%>%sort%>%unique()
  names(genes)<-NULL
  print(genes)
  diff<-diff[diff$symbol!='KPNA7',] # Blocks genes of interest 
  dfDiff<-data.frame(
    Chromosome=GRanges2bed(diff)$chrom,
    chromStart=GRanges2bed(diff)$start,
    chromEnd=GRanges2bed(diff)$end,
    Gene=diff$symbol
  )
  dfDiff<-distinct(dfDiff,Gene,.keep_all = TRUE)
  data1<-atacPeak$getMatchData(features, 'CTL')
  data2<-atacPeak$getMatchData(features, 'AIS')
  data3<-atacPeak$getMatchData(features, 'MIA')
  data4<-atacPeak$getMatchData(features, 'IAC')
  df1<-data.frame(
    feature2Bed(rownames(data1$methy)),
    methy=rowMeans(data1$methy),
    access=rowMeans(data1$access)
  )
  df2<-data.frame(
    feature2Bed(rownames(data2$methy)),
    methy=rowMeans(data2$methy),
    access=rowMeans(data2$access)
  )
  df3<-data.frame(
    feature2Bed(rownames(data3$methy)),
    methy=rowMeans(data3$methy),
    access=rowMeans(data3$access)
  )
  df4<-data.frame(
    feature2Bed(rownames(data4$methy)),
    methy=rowMeans(data4$methy),
    access=rowMeans(data4$access)
  )
  df<-list(df1,df2,df3,df4)
  
  data(UCSC.HG38.Human.CytoBandIdeogram)
  cyto.info <- UCSC.HG38.Human.CytoBandIdeogram
  RCircos.Set.Core.Components(cyto.info, chr.exclude=NULL,tracks.inside=5, tracks.outside=6 )
  RCircos.Set.Plot.Area()
  RCircos.Chromosome.Ideogram.Plot()
  RC.param <- RCircos.Get.Plot.Parameters()
  RC.param['track.background'] <- "white"
  RC.param['grid.line.color'] <- NULL
  RCircos.Reset.Plot.Parameters(RC.param)
  
  for(i in 1:4){
    RC.param <- RCircos.Get.Plot.Parameters()
    RC.param['hist.color'] <- colorMapStage2[i]
    RCircos.Reset.Plot.Parameters(RC.param)
    RCircos.Histogram.Plot(df[[i]], 4,i, 'in', 1);
  }
  
  for(i in 1:4){
    RC.param <- RCircos.Get.Plot.Parameters()
    RC.param['hist.color'] <- colorMapStage2[i]
    RCircos.Reset.Plot.Parameters(RC.param)
    RCircos.Histogram.Plot(df[[i]], 5,i, 'out', 1)
  }
  RCircos.Gene.Connector.Plot(dfDiff, 5, 'out')
  RC.param <- RCircos.Get.Plot.Parameters()
  # RC.param['text.size'] <- 0.5
  RCircos.Reset.Plot.Parameters(RC.param)
  RCircos.Gene.Name.Plot(dfDiff, 4,6, 'out')
}
#----------------------------------------------------------------------------------------------------------------------
# Figure 7B. mapk deg
#----------------------------------------------------------------------------------------------------------------------
rnaTPM<-RnaTPM('RNA')
saveImage2("rna.mapk.up.heatmap.pdf",width = 8,height = 4)
mapkGenes<-list(
  rtk=c('EGFR','ERBB2','ERBB3','INSR','MET'),
  gf=c('EGF','FGF1','FGFBP1','TGFA'),
  vgcc=c('CACNA1F','CACNB1','CACNB3'),
  kinase=c('RAC3','MAP3K13','MAP3K9','RRAS2','BRAF')
)

p1<-rnaTPM$plotStageHeatmap(mapkGenes$rtk, colAnno = TRUE, colNames = FALSE)
p2<-rnaTPM$plotStageHeatmap(mapkGenes$gf, colAnno = FALSE,colNames = FALSE)
p3<-rnaTPM$plotStageHeatmap(mapkGenes$vgcc, colAnno = FALSE,colNames = FALSE)
p4<-rnaTPM$plotStageHeatmap(mapkGenes$kinase, colAnno = FALSE,colNames = FALSE)
p1%v%p2%v%p3%v%p4
dev.off()
#----------------------------------------------------------------------------------------------------------------------
# Figure 7C. MAPK activity
#----------------------------------------------------------------------------------------------------------------------
rnaTPM<-RnaTPM('RNA')
p1<-rnaTPM$plotStageBar('FOS')
p2<-rnaTPM$plotStageBar('JUN')
p3<-rnaTPM$plotStageBar('MMP3')
p4<-rnaTPM$plotStageBar('MMP13')
saveImage2("rna.mapk.target.barplot.pdf",width = 5,height = 3.5)
grid.arrange(p1,p2,p3,p4,nrow=2)
dev.off()
#----------------------------------------------------------------------------------------------------------------------
# Figure 7D. AP1 target
#----------------------------------------------------------------------------------------------------------------------
atacPeak<-AtacPeak()
rnaTPM<-RnaTPM()

getFeatures<-function(TF,chipseq){
  region<-readRDS(file.path(CONFIG$dataIntermediate, 'atac','homer.mask','tf.epiTFs.DAR.rds'))
  gr2<-bed2GRanges(loadData2(chipseq,header = FALSE))
  chipH3k27ac=file.path(CONFIG$dataIntermediate, 'mapk','H3K27AC_ENCFF614QOP.bed')
  h3k27ac<-bed2GRanges(loadData2(chipH3k27ac,header = FALSE))
  atacTFAnno<-filter(region,tf==TF)
  atacTFBed<-feature2Bed(unique(atacTFAnno$feature))
  featureTF<-bed2Feature(atacTFBed)
  atacTF<-bed2GRanges(atacTFBed)
  peak<-atacTF[findOverlaps(atacTF, gr2)@from]%>%unique()
  tfPeakNum<-length(peak)
  peak<-peak[findOverlaps(peak, h3k27ac)@from]%>%unique()
  tfPeakInterH3k27acNum<-length(peak)
  print(table(atacTFBed$chrom)) ##########
  print(gr2@seqnames%>%table()) #########
  print(h3k27ac@seqnames%>%table()) #########
  print(peak@seqnames%>%table()) #########
  cat()
  cat(sprintf("tfPeakNum: %d; tfPeakInterH3k27acNum: %d", tfPeakNum, tfPeakInterH3k27acNum))
  features<-bed2Feature(GRanges2bed(peak))
  features
}
chipseqFOS=file.path(CONFIG$dataIntermediate, 'mapk','FOS_ENCFF459DPT.bed')
chipseqJUN=file.path(CONFIG$dataIntermediate, 'mapk','JUN_ENCFF280RQS.bed')
featureFOS<-getFeatures('FOS', chipseqFOS)
featureJUN<-getFeatures('JUN', chipseqJUN)
featureAP1<-unique(union(featureFOS, featureJUN)) # 819 AP-1 binding region

saveImage2("rna.ap1.target.circle.pdf",width = 5,height = 5)
plotFeatureRCircos(featureAP1)
dev.off()
#----------------------------------------------------------------------------------------------------------------------
# Figure 7F. Peak Accessibility/Methylation vs. gene expression TPM
#----------------------------------------------------------------------------------------------------------------------
atacPeak<-AtacPeak()
saveImage2("gene.access.methy.tpm.met.egfr.fgfbp1.scatter.pdf",width = 3.5,height = 6)
par(mfrow=c(3,2))
atacPeak$plotCorAccessVsExpression('chr7:116671720-116673040', 'MET',withTest = TRUE,featureTitle = TRUE)#*
atacPeak$plotCorMethyVsExpression('chr7:116671720-116673040', 'MET',withTest = TRUE)
atacPeak$plotCorAccessVsExpression('chr7:55018238-55020917', 'EGFR',withTest = TRUE,featureTitle = TRUE)
atacPeak$plotCorMethyVsExpression('chr7:55018238-55020917', 'EGFR',withTest = TRUE)#*
atacPeak$plotCorAccessVsExpression('chr4:15938055-15938795', 'FGFBP1',withTest = TRUE,featureTitle = TRUE)
atacPeak$plotCorMethyVsExpression('chr4:15938055-15938795', 'FGFBP1',withTest = TRUE)
dev.off()
#----------------------------------------------------------------------------------------------------------------------
# Figure S5C. AP-1 target gene which is related to survival
#----------------------------------------------------------------------------------------------------------------------
rnaTPM<-RnaTPM('RNA')
p1<-rnaTPM$plotStageBar('EFNA5')
p2<-rnaTPM$plotStageBar('PERP')
saveImage2("rna.mapk.ap1.target.efna5,perp.barplot.pdf",width = 5,height = 3.5)
grid.arrange(p1,p2,nrow=2, ncol=2)
dev.off()
