source('./R/base.R')
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(GenomicRanges)
library(rtracklayer)
library(reshape2)
library(ggExtra)
library(gridExtra)

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
write.csv(srdmc.count, file.path(CONFIG$dataResult, 'srdmc.s2.count.detail.csv'), quote = FALSE, row.names = FALSE)
write.csv(SRDMC.count, file.path(CONFIG$dataResult, 'srdmc.s2.count.csv'), quote = FALSE, row.names = FALSE)
#----------------------------------------------------------------------------------------------------------------------
# SR-DMC Mehtylation Levels Heatmap
#----------------------------------------------------------------------------------------------------------------------
dataj<-readRDS(file.path(CONFIG$dataIntermediate, 'dmc.methyLevel.rds'))
dataj<-split(dataj,dataj$class)
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
  print(dim(m1))
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
SRDMR<-loadSRDMR()
SRDMR.count<-count(SRDMR, class)%>%arrange(desc(n))
write.csv(SRDMR.count, file.path(CONFIG$dataResult, 'srdmr.s2.count.csv'), quote = FALSE, row.names = FALSE)
#----------------------------------------------------------------------------------------------------------------------
# SR-DMR in genomic regions
#----------------------------------------------------------------------------------------------------------------------
SRDMR<-loadSRDMR()
genomicRegion<-readRDS(file.path(CONFIG$dataIntermediate, 'genomicRegion.rds'))
SRDMR.genomicRegion<-sapply(split(SRDMR,SRDMR$class),function(x){
  gr2<-bed2GRanges(x)
  sapply(genomicRegion, function(gr1){
    sum(countOverlaps(gr1, gr2))
  })
})
plot.srdmr.barplot.genomicRegion <- function(regions){
  SRDMRgenomicRegion<-data.frame(SRDMR.genomicRegion,check.names = FALSE)
  SRDMRgenomicRegion$region<-rownames(SRDMRgenomicRegion)
  SRDMRgenomicRegion$region<-factor(regions[match(SRDMRgenomicRegion$region,names(regions))],levels = regions)
  SRDMRgenomicRegion<-SRDMRgenomicRegion[!is.na(SRDMRgenomicRegion$region),]
  # row.sum<-rowSums(SRDMRgenomicRegion[,1:4])
  # SRDMRgenomicRegion[,1:4]<-apply(SRDMRgenomicRegion[,1:4], 2, function(x) {x/row.sum} )
  # row.sum<-colSums(SRDMRgenomicRegion[,1:4])
  # SRDMRgenomicRegion[,1:4]<-t(apply(SRDMRgenomicRegion[,1:4], 1, function(x) {x/row.sum} ))
  SRDMRgenomicRegion[,1:4]<-log2(SRDMRgenomicRegion[,1:4])
  print(SRDMRgenomicRegion)
  df<-melt(SRDMRgenomicRegion,variable.name='Stage')
  ggplot(data=df, aes(y=value, x=region, fill=Stage))+
    scale_fill_manual(values=colorMapSRDMR)+
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

write.csv(data.frame(t(SRDMR.genomicRegion)), file.path(CONFIG$dataResult, 'srdmr.s2.count.genomicRegion.csv'), quote = FALSE)
#----------------------------------------------------------------------------------------------------------------------
# SR-DMR Density near TSS and CGI
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
       ylab = "SC-DMR Density",
       xlim=xlim,
       ylim=ylim,
       col = "blue")
  lines(dens$a, lwd = 2,col=colorMap[1])
  lines(dens$b, lwd = 2,col=colorMap[2])
  legend("topleft", legend=names(colorMap), fill=colorMap, bty = "n")
}

saveImage2("srdmr.density.genomicRegion.pdf",width = 6,height = 5)
par(mfrow = c(2, 2))
plot.scdmr.density.tss.vs.cgi('Early-Hyper-DMR')
plot.scdmr.density.tss.vs.cgi('Early-Hypo-DMR')
plot.scdmr.density.tss.vs.cgi('Late-Hyper-DMR')
plot.scdmr.density.tss.vs.cgi('Late-Hypo-DMR')
dev.off()

#----------------------------------------------------------------------------------------------------------------------
# SR-DMR Length and CpG Number
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
p1<-plot.scdmr.status('Early-Hyper-DMR')
p2<-plot.scdmr.status('Early-Hypo-DMR')
p3<-plot.scdmr.status('Late-Hyper-DMR')
p4<-plot.scdmr.status('Late-Hypo-DMR')
saveImage2("srdmr.status.scatter.cpg.length.pdf",width = 6,height = 6)
grid.arrange(p1, p2,p3,p4, nrow = 2)
dev.off()
srdmr.status<-sapply(SRDMR.list, function(x){
  out<-c(mean(x$cpg),mean(x$length))
  names(out)<-c('Average CpG Number', 'Average Size')
  out
})
srdmr.status<-data.frame(srdmr.status)
write.csv(srdmr.status, file.path(CONFIG$dataResult, 'srdmr.status.cpg.length.csv'),row.names  = TRUE,quote = FALSE)

#----------------------------------------------------------------------------------------------------------------------
# SR-DMR Homer
#----------------------------------------------------------------------------------------------------------------------
homerEarlyHyper=homerKnownTFs(file.path(CONFIG$dataIntermediate, 'homer.s2.mask','Early-Hyper-DMC'))
homerEarlyHypo=homerKnownTFs(file.path(CONFIG$dataIntermediate, 'homer.s2.mask','Early-Hypo-DMC'))
homerLateHyper=homerKnownTFs(file.path(CONFIG$dataIntermediate, 'homer.s2.mask','Late-Hyper-DMC'))
homerLateHypo=homerKnownTFs(file.path(CONFIG$dataIntermediate, 'homer.s2.mask','Late-Hypo-DMC'))

tfs<-do.call(rbind,list(
  data.frame(tf=homerEarlyHyper, motif=names(homerEarlyHyper),class='Early-Hyper-DMR'),
  data.frame(tf=homerEarlyHypo, motif=names(homerEarlyHypo),class='Early-Hypo-DMR'),
  data.frame(tf=homerLateHyper, motif=names(homerLateHyper),class='Late-Hyper-DMR'),
  data.frame(tf=homerLateHypo, motif=names(homerLateHypo),class='Late-Hypo-DMR')
))
tfWidthData<-dcast(tfs, tf~class,fun.aggregate = length)
motifWidthData<-dcast(tfs, motif~class,fun.aggregate = length)

motifTable<-data.frame(Motifs=motifWidthData$motif,
                       TSs=sapply(strsplit(motifWidthData$motif,'\\('),function(x){x[1]}),
                       motifWidthData[,2:5])
write.csv(motifTable, file.path(CONFIG$dataResult, 'srdmr.homer.motifs.csv'),row.names  = FALSE)
#----------------------------------------------------------------------------------------------------------------------
# SR-DMR Homer Figure. heatmap
#----------------------------------------------------------------------------------------------------------------------
TFS<-dplyr::arrange(tfWidthData, desc(`Early-Hyper-DMR`),desc(`Early-Hypo-DMR`),desc(`Late-Hyper-DMR`),desc(`Late-Hypo-DMR`))
m<-as.matrix(TFS[,2:ncol(TFS)])
rownames(m)<-TFS[,1]
m<-ifelse(m>0,1,0)

column_annotation <-HeatmapAnnotation(
  df=data.frame(SRDMR=names(colorMapSRDMR)),
  col = list(SRDMR =colorMapSRDMR),
  show_annotation_name =FALSE,
  annotation_name_side='left'
)
saveImage2("srdmr.homer.tfs.heatmap.pdf",width = 3,height = 12)
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
        row_names_gp = grid::gpar(fontsize = 5)
)
dev.off()
#----------------------------------------------------------------------------------------------------------------------
# SR-DMR Homer Figure. heatmap count
#----------------------------------------------------------------------------------------------------------------------
data<-tfWidthData[2:5]
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
