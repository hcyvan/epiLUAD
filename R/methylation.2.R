source('./R/base.R')
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(GenomicRanges)
library(rtracklayer)
library(reshape2)
library(ComplexHeatmap)

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

#----------------------------------------------------------------------------------------------------------------------
# SR-DMC
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
EarlyHyperDmc<-Reduce(rbind,list(
  filter(srdmc,AIS=='Hyper',MIA=="Hyper", IAC=="Hyper"),
  filter(srdmc,AIS=='Hyper',MIA=="Hyper", IAC=="NC"),
  filter(srdmc,AIS=='Hyper',MIA=="NC", IAC=="NC"),
  filter(srdmc,AIS=='NC',MIA=="Hyper", IAC=="NC"),
  filter(srdmc,AIS=='NC',MIA=="Hyper", IAC=="Hyper")
))%>%mutate(class='Early-Hyper-DMC')%>%dplyr::select(chrom, start,end,class)
EarlyHypoDmc<-Reduce(rbind,list(
  filter(srdmc,AIS=='Hypo',MIA=="Hypo", IAC=="Hypo"),
  filter(srdmc,AIS=='Hypo',MIA=="Hypo", IAC=="NC"),
  filter(srdmc,AIS=='Hypo',MIA=="NC", IAC=="NC"),
  filter(srdmc,AIS=='NC',MIA=="Hypo", IAC=="NC"),
  filter(srdmc,AIS=='NC',MIA=="Hypo", IAC=="Hypo")
))%>%mutate(class='Early-Hypo-DMC')%>%dplyr::select(chrom, start,end,class)


LateHyperDmc<-filter(srdmc,AIS=='NC',MIA=="NC", IAC=="Hyper")%>%mutate(class='Late-Hyper-DMC')%>%dplyr::select(chrom, start,end,class)
LateHypoDmc<-filter(srdmc,AIS=='NC',MIA=="NC", IAC=="Hypo")%>%mutate(class='Late-Hypo-DMC')%>%dplyr::select(chrom, start,end,class)

SRDMC<-Reduce(rbind,list(EarlyHyperDmc,EarlyHypoDmc,LateHyperDmc,LateHypoDmc))
SRDMC$chrom <- factor(SRDMC$chrom, levels=chromFactorLevel)
SRDMC<-arrange(SRDMC, chrom, start)

srdmc.count<-count(srdmc,AIS,IAC,MIA)%>%arrange(desc(n))
SRDMC.count<-count(SRDMC, class)%>%arrange(desc(n))

write.table(SRDMC, file.path(CONFIG$dataIntermediate, 'srdmc.s2.bed'), quote = FALSE, sep = "\t", row.names = FALSE,col.names=FALSE)
write.table(srdmc.count, file.path(CONFIG$dataIntermediate, 'srdmc.s2.count.detail.bed'), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(SRDMC.count, file.path(CONFIG$dataIntermediate, 'srdmc.s2.count.bed'), quote = FALSE, sep = "\t", row.names = FALSE)
#----------------------------------------------------------------------------------------------------------------------
# SR-DMC Mehtylation Levels Heatmap
#----------------------------------------------------------------------------------------------------------------------
select.and.rename <- function(x) {
  df.wgbs<-groups$WGBS$select(groupFactorLevel)
  select.data<-x[,match(idNew2Old(df.wgbs$SampleName),colnames(x))]
  colnames(select.data)<-df.wgbs$SampleName
  data.fix<-cbind(x[,1:4],select.data)
  data.fix
}
SRDMC<-loadData2(file.path(CONFIG$dataIntermediate, 'srdmc.s2.bed'),force.refresh = TRUE,header = FALSE)
SRDMC.ratio<-loadData2(file.path(CONFIG$dataIntermediate, 'srdmc.s2.ratio.bed'))
dataj<-left_join(SRDMC, SRDMC.ratio, by=c('V1'='#chrom','V2'='start'))
dataj<-dataj[,c(1:4, 6:ncol(dataj))]
dataj<-select.and.rename(dataj)
dataj<-split(dataj,dataj$V4)

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
  print(dim(m1))
  if (nrow(m1) >=10000){
    set.seed(123)
    m1<-m1[sample(1:nrow(m1),10000),]
  }
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
saveImage2("srdmr.heatmap.Early-Hyper-DMC.pdf",width = 5,height = 4)
plot.srdmc(dataj$`Early-Hyper-DMC`)
dev.off()
saveImage2("srdmr.heatmap.Early-Hypo-DMC.pdf",width = 5,height = 4)
plot.srdmc(dataj$`Early-Hypo-DMC`)
dev.off()
saveImage2("srdmr.heatmap.Late-Hyper-DMC.pdf",width = 5,height = 4)
plot.srdmc(dataj$`Late-Hyper-DMC`)
dev.off()
saveImage2("srdmr.heatmap.Late-Hypo-DMC.pdf",width = 5,height = 4)
plot.srdmc(dataj$`Late-Hypo-DMC`)
dev.off()
#----------------------------------------------------------------------------------------------------------------------
# SR-DMR
# eg: python ./script/dmc2dmr.py -i ./data/intermediate/srdmc.s2.bed -o ./data/intermediate/srdmr.s2.bed
#----------------------------------------------------------------------------------------------------------------------
loadSRDMR<-function(){
  SRDMR<-loadData2(file.path(CONFIG$dataIntermediate, 'srdmr.s2.bed'),header = FALSE)
  colnames(SRDMR)<-c('chrom','start', 'end', 'class','cpg','length')
  SRDMR$class<-sub('DMC','DMR',SRDMR$class)
  SRDMR
}

SRDMR<-loadSRDMR()
SRDMR.count<-count(SRDMR, class)%>%arrange(desc(n))
write.table(SRDMR.count, file.path(CONFIG$dataIntermediate, 'srdmr.s2.count.bed'), quote = FALSE, sep = "\t", row.names = FALSE)


SRDMR<-loadSRDMR()
genomicRegion<-readRDS(file.path(CONFIG$dataIntermediate, 'genomicRegion.rds'))
SRDMR.genomicRegion<-sapply(split(SRDMR,SRDMR$class),function(x){
  gr2<-bed2GRanges(x)
  sapply(genomicRegion, function(gr1){
    sum(countOverlaps(gr1, gr2))
  })
})
#----------------------------------------------------------------------------------------------------------------------
# SR-DMR in genomic regions
#----------------------------------------------------------------------------------------------------------------------
plot.srdmr.barplot.genomicRegion <- function(regions){
  SRDMRgenomicRegion<-data.frame(SRDMR.genomicRegion,check.names = FALSE)
  SRDMRgenomicRegion$region<-rownames(SRDMRgenomicRegion)
  SRDMRgenomicRegion$region<-factor(regions[match(SRDMRgenomicRegion$region,names(regions))],levels = regions)
  SRDMRgenomicRegion<-SRDMRgenomicRegion[!is.na(SRDMRgenomicRegion$region),]
  df<-melt(SRDMRgenomicRegion,variable.name='Stage')
  stage.colormap<-c("#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd")
  names(stage.colormap)<-levels(df$variable)
  df$value<-log2(df$value)
  ggplot(data=df, aes(y=value, x=region, fill=Stage))+
    scale_fill_manual(values=stage.colormap)+
    geom_bar(stat="identity", position=position_dodge())+
    ylab("Log2(SR-DMR Number)")+
    xlab("Genomic Regions")+
    theme_classic()
}
saveImage2("srdmr.barplot.genomicRegion.cpg.pdf",width = 5,height = 3)
regions<-c('cgIslands'='CpG Islands', 'cgShores'='CpG Shores', 'cgShelves'='CpG Shelves', 'cgSea'='CpG Sea')
plot.srdmr.barplot.genomicRegion(regions)
dev.off()
saveImage2("srdmr.barplot.genomicRegion.promoter.pdf",width = 8,height = 3)
regions<-c('promoter.1k'='Promoter.1k', 'promoter.5k'='Promoter.5k', 'utr5'="5'UTR", 'utr3'="3'UTR",'exons'='Exons','intron'='Intron','intergenic'='Intergenic')
plot.srdmr.barplot.genomicRegion(regions)
dev.off()






