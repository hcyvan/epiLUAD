source('./R/base.R')
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(GenomicRanges)
library(rtracklayer)
library(reshape2)
library(ggExtra)
library(gridExtra)
library(ggpubr)
library(Hmisc)
library(readxl)


drawDensity<-function(ratio1,ratio2,s1,s2,legend.position='topleft'){
  d1 <- density(ratio1)
  d2 <- density(ratio2)
  dens <- list(a=d1,b=d2)
  plot(NA, xlim=range(sapply(dens, "[", "x")),
       ylim=range(sapply(dens, "[", "y")),xlab="",ylab="Density",cex.lab=2, cex.axis=1.5)
  mapply(lines, dens, col=1:length(dens))
  polygon(d1, col=rgb(0, 1, 0,0.5), border=NA)
  polygon(d2, col=rgb(1, 0, 0,0.5), border=NA)
  legend(legend.position, legend=c(s1,s2), fill=c("green","red"), bty = "n")
}

#----------------------------------------------------------------------------------------------------------------------
# Figure 1A: The Average DNA methylation level of samples in each group
#----------------------------------------------------------------------------------------------------------------------
tableS3 <- read_excel(file.path(CONFIG$dataExternal, 'SupplementaryData.xlsx'),sheet = 'Table S3')

samples<-groups$WGBS$select(groupFactorLevel)
samplesMatch<-tableS3[match(samples$SampleName,tableS3$SampleName),]
data<-data.frame(
  group=samples$Group,
  ratio=samplesMatch$MCALL_MeanRatioCG_3X
)
saveImage2("methylation.level.mean.pdf",width = 5,height = 3)
ggplot(data=data,aes(x=group,y=ratio,fill=group))+
  scale_fill_manual(values=colorMapStage) +
  # geom_violin(trim=FALSE) +
  geom_violin(draw_quantiles = NULL, colour = NA,trim=FALSE)+
  geom_jitter(shape=17, position=position_jitter(0.2), colour='black', size=0.5)+
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="black")+
  theme_classic()+
  stat_compare_means( comparisons = list(c('CTL','AIS'),c('CTL','MIA'),c('CTL','IAC')),
                      label = 'p.signif', method = "t.test")+
  labs(x='',y='Mean Methylation Level')+
  theme(legend.position="right",
        axis.title.x = element_text(size=0),
        axis.title.y = element_text(size=15),
        axis.text = element_text(size = 12,colour="black"),
        legend.title = element_blank(),
        legend.text = element_text(size=12))+
  guides(colour = guide_legend(override.aes = list(shape = 12,size=10)))
dev.off()
#----------------------------------------------------------------------------------------------------------------------
# Figure 1B: The methylation level density distribution of CTL vs. AIS, CTL vs. MIA and CTL vs. IAC
#----------------------------------------------------------------------------------------------------------------------
ratio <- loadData2(file.path(CONFIG$dataExternal, 'LAD.group.ratio.bed'))
ratio<-removeNegativeOne(ratio)
saveImage2("methylation.level.density.pdf",width = 4,height = 3.5)
layout(1:3)
par(mar = c(3,5,1,1))
drawDensity(ratio$CTL, ratio$AIS, "CTL", "AIS")
drawDensity(ratio$CTL, ratio$MIA, "CTL", "MIA")
drawDensity(ratio$CTL, ratio$IAC, "CTL", "IAC")
dev.off()
#----------------------------------------------------------------------------------------------------------------------
# Figure 1C: The Average DNA methylation level of some Genomic Regions
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
  df=data.frame(Stage=samples$Group),
  col = list(Stage =color.map),
  show_annotation_name =FALSE,
  annotation_name_side='left'
)
saveImage2("genomicRegion.heatmap.pdf",width = 12,height = 2.5)
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
# Figure 1D. Methylation level of CpGs within 5,000 bp upstream and downstream relative to TSS
#----------------------------------------------------------------------------------------------------------------------
smooth2 <- function(hw=51) {
  function(arr){
    for(i in 1:length(arr)){
      left<-i-hw
      if(left < 1){
        left<-1
      }
      right<-i+hw
      if(right > length(arr)){
        right<-length(arr)
      }
      
      arr[i]<-mean(arr[left:right])
    }
    arr
  }
}


plot.group.methy.profile<- function(UP,DWON,signal.file,ylab="Methylation Level",xlab="Distance to TSS (bp)",cex=1.5,horiz=FALSE, draw.legend=TRUE){
  x<-seq(UP, DWON)
  signalMatrixGroup<-loadData2(signal.file)
  mmData<-signalMatrixGroup[,-1]
  mm<-mmData[(15001+UP):(15001+DWON),]
  # mm<-t(scale(t(mm)))
  mm<-as.data.frame(apply(mm,2,smooth2(51)))
  mm<-as.data.frame(mm)
  par(mar = c(5,5,1,1))
  if (draw.legend){
    plot(NA, xlim=c(UP, DWON), xlab=xlab,ylab=ylab,cex.lab=1.5, cex.axis=1.5,ylim=c(min(mm),max(mm)), bty='n')
  }else{
    plot(NA, xlim=c(UP, DWON),cex.lab=1.5, cex.axis=1.5,ylim=c(min(mm),max(mm)), xaxt='n',yaxt='n',xlab = "",ylab = "")
  }
  lines(x,mm[["CTL"]],col=colorMapStage[1])
  lines(x,mm[["AIS"]],col=colorMapStage[2])
  lines(x,mm[["MIA"]],col=colorMapStage[3])
  lines(x,mm[["IAC"]],col=colorMapStage[4])
  if (draw.legend){
    legend("bottomleft",legend=names(colorMapStage),fill=colorMapStage,bty = "n",cex=cex,horiz = horiz)
  }
}
saveImage2("methylation.level.profile.cgi.pdf",width = 5,height = 5)
plot.group.methy.profile(-5000, 5000,file.path(CONFIG$dataIntermediate, "signalRoundCgiMatrixGroup.bed"),xlab="Distance to CGI Center (bp)")
dev.off()
saveImage2("methylation.level.profile.cgi.2.pdf",width = 4,height = 4)
plot.group.methy.profile(-500, 500,file.path(CONFIG$dataIntermediate, "signalRoundCgiMatrixGroup.bed"),xlab="Distance to CGI Center (bp)",draw.legend=FALSE)
dev.off()
saveImage2("methylation.level.profile.cgi.3.pdf",width = 4,height = 4)
plot.group.methy.profile(2000, 3000,file.path(CONFIG$dataIntermediate, "signalRoundCgiMatrixGroup.bed"),xlab="Distance to CGI Center (bp)",draw.legend=FALSE)
dev.off()

saveImage2("methylation.level.profile.tss.pdf",width = 5,height = 5)
plot.group.methy.profile(-5000, 5000,file.path(CONFIG$dataIntermediate, "signalRoundTssMatrixGroup.bed"),xlab="Distance to TSS (bp)")
dev.off()
saveImage2("methylation.level.profile.tss.2.pdf",width = 4,height = 4)
plot.group.methy.profile(-500, 500,file.path(CONFIG$dataIntermediate, "signalRoundTssMatrixGroup.bed"),xlab="Distance to TSS (bp)",draw.legend=FALSE)
dev.off()
saveImage2("methylation.level.profile.tss.3.pdf",width = 4,height = 4)
plot.group.methy.profile(2000, 3000,file.path(CONFIG$dataIntermediate, "signalRoundTssMatrixGroup.bed"),xlab="Distance to TSS (bp)",draw.legend=FALSE)
dev.off()

#----------------------------------------------------------------------------------------------------------------------
# Table S5. Hyper-DMCs and Hypo-DMCs of different group in LUAD				
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
# Table S6. Hyper-DMCs and Hypo-DMCs in different genomic regions
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
    rbind(sum(countOverlaps(gr1, hyper)), sum(countOverlaps(gr1, hypo)))
  })
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
#----------------------------------------------------------------------------------------------------------------------
# Table S7. Stage-Related DMCs
#----------------------------------------------------------------------------------------------------------------------
dmcCTLvsAIS<-loadData2(file.path(CONFIG$dataExternal,'dmc', 'dmc.CTL.vs.AIS.txt'),file.format='bed')
dmcCTLvsMIA<-loadData2(file.path(CONFIG$dataExternal,'dmc', 'dmc.CTL.vs.MIA.txt'),file.format='bed')
dmcCTLvsIAC<-loadData2(file.path(CONFIG$dataExternal,'dmc', 'dmc.CTL.vs.IAC.txt'),file.format='bed')

extractDmc <- function(dmc,group) {
  dmc2<-dmc[,match(c('#chrom','start','end', 'class'),colnames(dmc))]
  colnames(dmc2)<-c('chrom','start','end', 'class')
  dmc2$class<-'Hypo'
  dmc2$class[dmc$class=='strongHyper'] <- 'Hyper'
  dmc2$group<-group
  dmc2
}
dmc<-list(
  AIS=extractDmc(dmcCTLvsAIS,'AIS'),
  MIA=extractDmc(dmcCTLvsMIA,'MIA'),
  IAC=extractDmc(dmcCTLvsIAC,'IAC')
)
lData<-Reduce(rbind, dmc)
srdmc<-dcast(lData,chrom + start+ end ~ group,value.var = 'class')
srdmc[,4:6][is.na(srdmc[,4:6])] <- "NC"
HyperDmcInAIS<-Reduce(rbind,list(
  filter(srdmc,AIS=='Hyper',MIA=="Hyper", IAC=="Hyper"),
  filter(srdmc,AIS=='Hyper',MIA=="Hyper", IAC=="NC"),
  filter(srdmc,AIS=='Hyper',MIA=="NC", IAC=="NC")
))%>%mutate(class='HyperInAIS')%>%dplyr::select(chrom, start,end,class)
HyperDmcInMIA<-Reduce(rbind,list(
  filter(srdmc,AIS=='NC',MIA=="Hyper", IAC=="NC"),
  filter(srdmc,AIS=='NC',MIA=="Hyper", IAC=="Hyper")
))%>%mutate(class='HyperInMIA')%>%dplyr::select(chrom, start,end,class)
HyperDmcInIAC<-Reduce(rbind,list(
  filter(srdmc,AIS=='NC',MIA=="NC", IAC=="Hyper")
))%>%mutate(class='HyperInIAC')%>%dplyr::select(chrom, start,end,class)
HypoDmcInAIS<-Reduce(rbind,list(
  filter(srdmc,AIS=='Hypo',MIA=="Hypo", IAC=="Hypo"),
  filter(srdmc,AIS=='Hypo',MIA=="Hypo", IAC=="NC"),
  filter(srdmc,AIS=='Hypo',MIA=="NC", IAC=="NC")
))%>%mutate(class='HypoInAIS')%>%dplyr::select(chrom, start,end,class)
HypoDmcInMIA<-Reduce(rbind,list(
  filter(srdmc,AIS=='NC',MIA=="Hypo", IAC=="NC"),
  filter(srdmc,AIS=='NC',MIA=="Hypo", IAC=="Hypo")
))%>%mutate(class='HypoInMIA')%>%dplyr::select(chrom, start,end,class)
HypoDmcInIAC<-Reduce(rbind,list(
  filter(srdmc,AIS=='NC',MIA=="NC", IAC=="Hypo")
))%>%mutate(class='HypoInIAC')%>%dplyr::select(chrom, start,end,class)
SRDMC<-Reduce(rbind,list(
  HyperDmcInAIS=HyperDmcInAIS,
  HyperDmcInMIA=HyperDmcInMIA,
  HyperDmcInIAC=HyperDmcInIAC,
  HypoDmcInAIS=HypoDmcInAIS,
  HypoDmcInMIA=HypoDmcInMIA,
  HypoDmcInIAC=HypoDmcInIAC
))
SRDMC$chrom <- factor(SRDMC$chrom, levels=chromFactorLevel)
SRDMC<-arrange(SRDMC, chrom, start)
srdmc.count<-count(srdmc,AIS,IAC,MIA)%>%dplyr::select(AIS, MIA, IAC, n)%>%arrange(AIS, MIA, IAC)
SRDMC.count<-count(SRDMC, class)%>%arrange(desc(n))

saveTsv(SRDMC, file.path(CONFIG$dataIntermediate, 'wgbs','srdmc.bed'),col.names=FALSE)
saveCsv(srdmc.count, file.path(CONFIG$dataIntermediate, 'wgbs','srdmc.count.detail.csv'))
saveCsv(SRDMC.count, file.path(CONFIG$dataIntermediate, 'wgbs','srdmc.count.csv'))
#----------------------------------------------------------------------------------------------------------------------
# Figure 2B-D. SRDMC Mehtylation Levels Heatmap
#----------------------------------------------------------------------------------------------------------------------
srdmcMethyLevel<-readRDS(file.path(CONFIG$dataIntermediate, 'wgbs','dmc.methyLevel.rds'))
# dataj<-readRDS(file.path(CONFIG$dataIntermediate, 'dmc.methyLevel.rds'))
srdmcMethyLevel<-split(srdmcMethyLevel,srdmcMethyLevel$class)
plot.srdmc<-function(data){
  m<-as.matrix(data[,5:ncol(data)])
  m1<-m[rowSums(is.na(m))==0,]
  samples<-groups$WGBS$selectBySample(colnames(m))
  color.map<-groups$WGBS$getColorMapVec()
  column_annotation <-HeatmapAnnotation(
    df=data.frame(Stage=samples$Group),
    col = list(Stage =color.map),
    show_annotation_name =FALSE,
    annotation_name_side='left'
  )
  nrow1<-nrow(m1)
  if (nrow(m1) >=1000){
    set.seed(123)
    m1<-m1[sample(1:nrow(m1),1000),]
  }
  nrow2<-nrow(m1)
  print(sprintf("%d => %d", nrow1, nrow2))
  Heatmap(m1,
          top_annotation = column_annotation,
          cluster_rows=TRUE,
          cluster_columns = FALSE,
          show_row_names=FALSE,
          show_column_names=FALSE,
          heatmap_legend_param = list(
            title = "DNA Methylation Levels",
            legend_height = unit(4, "cm"),
            at = c(0,0.5,1),
            labels = c('0','0.5','1'),
            title_position = "lefttop-rot"
          ),
  )
}

saveImage2("srdmr.heatmap.HyperInAIS.pdf",width = 3.5,height = 2)
plot.srdmc(srdmcMethyLevel$HyperInAIS)
dev.off()
saveImage2("srdmr.heatmap.HyperInMIA.pdf",width = 3.5,height = 2)
plot.srdmc(srdmcMethyLevel$HyperInMIA)
dev.off()
saveImage2("srdmr.heatmap.HyperInIAC.pdf",width = 3.5,height = 2)
plot.srdmc(srdmcMethyLevel$HyperInIAC)
dev.off()
saveImage2("srdmr.heatmap.HypoInAIS.pdf",width = 3.5,height = 2)
plot.srdmc(srdmcMethyLevel$HypoInAIS)
dev.off()
saveImage2("srdmr.heatmap.HypoInMIA.pdf",width = 3.5,height = 2)
plot.srdmc(srdmcMethyLevel$HypoInMIA)
dev.off()
saveImage2("srdmr.heatmap.HypoInIAC.pdf",width = 3.5,height = 2)
plot.srdmc(srdmcMethyLevel$HypoInIAC)
dev.off()
#----------------------------------------------------------------------------------------------------------------------
# SR-DMR
# eg: python ./script/dmc2dmr.py -i ./data/intermediate/wgbs/srdmc.bed -o ./data/intermediate/wgbs/srdmr.bed
#----------------------------------------------------------------------------------------------------------------------
SRDMR<-loadSRDMR()
SRDMR.count<-count(SRDMR, class)%>%arrange(desc(n))
saveCsv(SRDMR.count, file.path(CONFIG$dataIntermediate, 'wgbs','srdmr.count.csv'))
saveTsv(filter(SRDMR,class=='HyperInAIS')[,1:3], file.path(CONFIG$dataIntermediate, 'wgbs','srdmr.HyperInAIS.bed'),col.names = FALSE)
saveTsv(filter(SRDMR,class=='HypoInAIS')[,1:3], file.path(CONFIG$dataIntermediate, 'wgbs','srdmr.HypoInAIS.bed'),col.names = FALSE)
saveTsv(filter(SRDMR,class=='HyperInMIA')[,1:3], file.path(CONFIG$dataIntermediate, 'wgbs','srdmr.HyperInMIA.bed'),col.names = FALSE)
saveTsv(filter(SRDMR,class=='HypoInMIA')[,1:3], file.path(CONFIG$dataIntermediate, 'wgbs','srdmr.HypoInMIA.bed'),col.names = FALSE)
saveTsv(filter(SRDMR,class=='HyperInIAC')[,1:3], file.path(CONFIG$dataIntermediate, 'wgbs','srdmr.HyperInIAC.bed'),col.names = FALSE)
saveTsv(filter(SRDMR,class=='HypoInIAC')[,1:3], file.path(CONFIG$dataIntermediate, 'wgbs','srdmr.HypoInIAC.bed'),col.names = FALSE)
#----------------------------------------------------------------------------------------------------------------------
# Table S8. SR-DMR in genomic regions
#----------------------------------------------------------------------------------------------------------------------
SRDMR<-loadSRDMR()
genomicRegion<-readRDS(file.path(CONFIG$dataIntermediate, 'genomicRegion.rds'))
SRDMR.genomicRegion<-sapply(split(SRDMR,SRDMR$class),function(x){
  gr2<-bed2GRanges(x)
  sapply(genomicRegion, function(gr1){
    sum(countOverlaps(gr2, gr1)!=0)
  })
})
SRDMR.genomicRegion
plot.srdmr.barplot.genomicRegion <- function(regions){
  SRDMRgenomicRegion<-data.frame(SRDMR.genomicRegion,check.names = FALSE)
  SRDMRgenomicRegion$region<-rownames(SRDMRgenomicRegion)
  SRDMRgenomicRegion$region<-factor(regions[match(SRDMRgenomicRegion$region,names(regions))],levels = regions)
  SRDMRgenomicRegion<-SRDMRgenomicRegion[!is.na(SRDMRgenomicRegion$region),]
  df<-melt(SRDMRgenomicRegion,variable.name='Stage')
  ggplot(data=df, aes(y=value, x=region, fill=Stage))+
    scale_fill_manual(values=colorMapSRDMR)+
    geom_bar(stat="identity", position=position_dodge())+
    ylab("SRDMRs Number")+
    xlab("Genomic Regions")+
    theme_classic()
}
saveImage2("srdmr.barplot.genomicRegion.cpg.pdf",width = 6,height = 2.5)
regions<-c('cgIslands'='CpG Islands', 'cgShores'='CpG Shores', 'cgShelves'='CpG Shelves', 'cgSea'='CpG Sea')
plot.srdmr.barplot.genomicRegion(regions)
dev.off()
saveImage2("srdmr.barplot.genomicRegion.promoter.pdf",width = 8,height = 2.5)
regions<-c('promoter.1k'='Promoter.1k', 'promoter.5k'='Promoter.5k', 'utr5'="5'UTR", 'utr3'="3'UTR",'exons'='Exons','intron'='Intron','intergenic'='Intergenic')
plot.srdmr.barplot.genomicRegion(regions)
dev.off()

out<-data.frame(
  SRDMR=colnames(SRDMR.genomicRegion),
  sum=colSums(SRDMR.genomicRegion),
  t(SRDMR.genomicRegion)
)
saveCsv(out, file.path(CONFIG$dataIntermediate, 'wgbs','srdmr.count.genomicRegion.csv'))
#----------------------------------------------------------------------------------------------------------------------
# Figure 2F. SRDMRs Density near TSS and CGI
#----------------------------------------------------------------------------------------------------------------------
genomicRegion<-readRDS(file.path(CONFIG$dataIntermediate, 'genomicRegion.rds'))
cgIslands.gr<-genomicRegion$cgIslands
tss.gr<-genomicRegion$tss
SRDMR<-bed2GRanges(loadSRDMR())

annoDistQueryToSubject<-function(query, subject){
  midpoint<-function(gr){
    mid<-floor((start(gr) + end(gr)) / 2)
    start(gr)<-mid
    end(gr)<-mid
    gr
  }
  subject<-midpoint(subject)
  hits<-distanceToNearest(query, subject)
  query.hits<-query[queryHits(hits)]
  subject.hits<-subject[subjectHits(hits)]
  dist.sign<-ifelse((start(query.hits) > end(subject.hits))&(mcols(hits)$distance!=0),-1,1)
  query.hits$distanceToSubject<-mcols(hits)$distance * dist.sign
  query.hits
}
distToCGI<-annoDistQueryToSubject(SRDMR,cgIslands.gr)
distToTSS<-annoDistQueryToSubject(SRDMR,tss.gr)

set.seed(123)
lim<-50000
dist<-distToCGI
dist<-dist[abs(dist$distanceToSubject)<lim,]
dens<-lapply(split(dist, dist$class), function(x){
  density(x$distanceToSubject,bw = 'nrd0')
})
xlim<-range(sapply(dens, "[", "x"))
ylim<-range(sapply(dens, "[", "y"))

plot.scdmr.density.tss.vs.cgi<-function(g) {
  dist1<-distToTSS
  dist2<-distToCGI
  dist1<-dist1[dist1$class==g,]
  dist1<-dist1[abs(dist1$distanceToSubject)<lim,]
  dist2<-dist2[dist2$class==g,]
  dist2<-dist2[abs(dist2$distanceToSubject)<lim,]
  dens<-list(
    a=density(dist1$distanceToSubject),
    b=density(dist2$distanceToSubject)
  )
  # xlim<-range(sapply(dens, "[", "x"))
  # ylim<-range(sapply(dens, "[", "y"))
  colorMap<-c("#0000FF80","#FF000080")
  names(colorMap)<-c("TSS","CGI")
  plot(NA, main = g,
       xlab = "Distance to TSS/CGI",
       ylab = "SCDMR Density",
       xlim=xlim,
       ylim=ylim,
       col = "blue")
  lines(dens$a, lwd = 2,col=colorMap[1])
  lines(dens$b, lwd = 2,col=colorMap[2])
  legend("topleft", legend=names(colorMap), fill=colorMap, bty = "n")
}

saveImage2("srdmr.density.genomicRegion.pdf",width = 3.5,height = 6)
par(mfrow = c(3, 2))
plot.scdmr.density.tss.vs.cgi('HyperInAIS')
plot.scdmr.density.tss.vs.cgi('HypoInAIS')
plot.scdmr.density.tss.vs.cgi('HyperInMIA')
plot.scdmr.density.tss.vs.cgi('HypoInMIA')
plot.scdmr.density.tss.vs.cgi('HyperInIAC')
plot.scdmr.density.tss.vs.cgi('HypoInIAC')
dev.off()
#----------------------------------------------------------------------------------------------------------------------
# Figure 2F,G. GREAT analysis of SRDMRs
#----------------------------------------------------------------------------------------------------------------------
saveImage2("srdmr.great.HyperInAIS.BP.pdf",width = 4,height = 3.5)
plot.great(file.path(CONFIG$dataIntermediate,'wgbs', "srdmr.HyperInAIS.GOBiologicalProcess.tsv"),title="")
dev.off()
saveImage2("srdmr.great.HypoInIAC.BP.pdf",width = 4,height = 3.5)
plot.great(file.path(CONFIG$dataIntermediate,'wgbs', "srdmr.HypoInIAC.GOBiologicalProcess.tsv"),title="")
dev.off()
#----------------------------------------------------------------------------------------------------------------------
# Figure S1A. SRDMRs Length and CpG Number
#----------------------------------------------------------------------------------------------------------------------
SRDMR<-bed2GRanges(loadSRDMR())
SRDMR.list<-split(SRDMR,SRDMR$class)
plot.scdmr.status<-function(srdmr.type){
  SRDMR.list<-split(SRDMR,SRDMR$class)
  gr<-SRDMR.list[[srdmr.type]]
  data<-data.frame(x=gr$length, y=gr$cpg)
  p_scatter <- ggplot(data, aes(x = x, y = y)) +
    geom_point(size=0.3,alpha = 0.6) +
    xlab("SC-DMR size (bp)")+
    ylab("SC-DMR Number")+
    ggtitle(srdmr.type)+
    theme_bw()
  p<-ggMarginal(p_scatter, type="histogram",fill='white',bins = 100,size = 8)
  p
}
p1<-plot.scdmr.status('HyperInAIS')
p2<-plot.scdmr.status('HypoInAIS')
p3<-plot.scdmr.status('HyperInMIA')
p4<-plot.scdmr.status('HypoInMIA')
p5<-plot.scdmr.status('HyperInIAC')
p6<-plot.scdmr.status('HypoInIAC')
saveImage2("srdmr.status.scatter.cpg.length.pdf",width = 4,height = 6.5)
grid.arrange(p1,p2,p3,p4,p5,p6, nrow = 3)
dev.off()
srdmr.status<-sapply(SRDMR.list, function(x){
  out<-c(mean(x$cpg),mean(x$length))
  names(out)<-c('Average CpG Number', 'Average Size')
  out
})
srdmr.status<-data.frame(srdmr.status)
write.csv(srdmr.status, file.path(CONFIG$dataResult, 'srdmr.status.cpg.length.csv'),row.names  = TRUE,quote = FALSE)
#----------------------------------------------------------------------------------------------------------------------
# Table S9. SRDMRs Homer
#----------------------------------------------------------------------------------------------------------------------
HyperInAIS=homerKnownTFs(file.path(CONFIG$dataIntermediate, 'wgbs','homer.mask','HyperInAIS'), 'HyperInAIS')
HyperInMIA=homerKnownTFs(file.path(CONFIG$dataIntermediate, 'wgbs','homer.mask','HyperInMIA'), 'HyperInMIA')
HyperInIAC=homerKnownTFs(file.path(CONFIG$dataIntermediate, 'wgbs','homer.mask','HyperInIAC'), 'HyperInIAC')
HypoInIAC=homerKnownTFs(file.path(CONFIG$dataIntermediate, 'wgbs','homer.mask','HypoInIAC'), 'HypoInIAC')

tfs<-do.call(rbind,list(HyperInAIS,HyperInMIA,HyperInIAC,HypoInIAC))
srdmrTFs<-getSRTFS(tfs)
saveRDS(srdmrTFs,file.path(CONFIG$dataIntermediate,'wgbs', 'srdmr.tfs.rds'))
#----------------------------------------------------------------------------------------------------------------------
# Figure S1B. SR-DMR Homer Figure. heatmap
#----------------------------------------------------------------------------------------------------------------------
srdmrTFs<-readRDS(file.path(CONFIG$dataIntermediate,'wgbs', 'srdmr.tfs.rds'))
TFS<-dplyr::arrange(srdmrTFs$stage, desc(HyperInAIS),desc(HypoInAIS),desc(HyperInMIA),desc(HypoInMIA),desc(HyperInIAC),desc(HypoInIAC))
m<-as.matrix(TFS)
column_annotation <-HeatmapAnnotation(
  df=data.frame(SRDMR=names(colorMapGroup)),
  col = list(SRDMR =colorMapGroup),
  show_annotation_name =FALSE,
  annotation_name_side='left'
)
saveImage2("srdmr.homer.tfs.heatmap.pdf",width = 3,height = 8)
Heatmap(m,
        cluster_rows=FALSE,
        cluster_columns = FALSE,
        show_row_names=TRUE,
        show_column_names=FALSE,
        border = TRUE,
        show_heatmap_legend = FALSE,
        col = colorRamp2(breaks = c(0,1.3), colors = c('white', 'black')),
        top_annotation = column_annotation,
        column_names_gp = grid::gpar(fontsize = 7),
        row_names_gp = grid::gpar(fontsize = 6)
)
dev.off()
#----------------------------------------------------------------------------------------------------------------------
# Figure S1C. SR-DMR Homer Figure. heatmap count
#----------------------------------------------------------------------------------------------------------------------
srdmrTFs<-readRDS(file.path(CONFIG$dataIntermediate,'wgbs', 'srdmr.tfs.rds'))
data<-srdmrTFs$stage
mat<-sapply(names(data), function(x){
  sapply(names(data), function(y){
    sum(data[[x]]&data[[y]])
  })
})
mat[upper.tri(mat)] <- NA

saveImage2("srdmr.homer.tfs.heatmap.count.pdf",width = 4,height = 4)
Heatmap(mat,
        cluster_rows=FALSE,
        cluster_columns = FALSE,
        na_col = "white",
        show_heatmap_legend = FALSE,
        col = colorRamp2(breaks = c(min(mat,na.rm=TRUE),max(mat,na.rm=TRUE)), colors = c('#D3D3D3', '#A9A9A9')),
        cell_fun = function(j, i, x, y, width, height, fill) {
          if(!is.na(mat[i, j])){
            grid.text(sprintf("%d", mat[i, j]), x, y, gp = gpar(fontsize = 13, col='black'))
          }
        })
dev.off()
