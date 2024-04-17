source('./R/base.R')
library(ChIPseeker)
library(DESeq2)
library(reshape2)
library(ggplot2)
library(ggExtra)
library(gridExtra)
library(circlize)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

#----------------------------------------------------------------------------------------------------------------------
# Figure 3A. ATAC-seq Peak
#----------------------------------------------------------------------------------------------------------------------
atacPeak<-loadData2(file.path(CONFIG$dataIntermediate, 'atac','atacPeak.bed'),header = FALSE, force.refresh = TRUE)
set.seed(123)
atacPeakCountToHeatmap<-atacPeak[sample(1:nrow(atacPeak), 6000),]
peak<-bed2GRanges(atacPeakCountToHeatmap)
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(peak, windows=promoter)
saveImage2("atac.peak.heatmap.pdf",width = 3,height = 6)
tagHeatmap(tagMatrix)
dev.off()

#----------------------------------------------------------------------------------------------------------------------
# ATAC-seq: Identify DARs
#----------------------------------------------------------------------------------------------------------------------
atacPeakCount<-readRDS(file.path(CONFIG$dataIntermediate,'atac', 'atacPeakCount.rds'))
peaks<-paste(paste(atacPeakCount[,1],atacPeakCount[,2],sep = ':'),atacPeakCount[,3], sep = '-')
findDAR<-function(s0,s1,minCount=10) {
  samples0<-groups$ATAC$select(s0)
  samples1<-groups$ATAC$select(s1)
  count<-atacPeakCount[,4:ncol(atacPeakCount)]
  g0<-match(samples0$SampleName,colnames(count))
  g1<-match(samples1$SampleName,colnames(count))
  countData<-cbind(count[,g0], count[,g1])
  rownames(countData)<-peaks
  keep<-rowSums(countData>=minCount)>=2
  countData<-countData[keep,]
  colData <- data.frame(condition = groups$ATAC$select(c(s0,s1))$Group)
  dds <- DESeqDataSetFromMatrix(countData = countData,
                                colData = colData,
                                design = ~ condition)
  dds <- DESeq(dds)
  res <- results(dds)
  res <- res[order(res$pvalue), ]
  diff <- subset(res, pvalue <= 0.001, abs(log2FoldChange)>=log2(1.8))
  hyper<-subset(diff, log2FoldChange>0)
  hypo<-subset(diff, log2FoldChange<0)
  hyperNum<-sum(diff$log2FoldChange > 0)
  hypoNum<-sum(diff$log2FoldChange < 0)
  print(paste('keep', sum(keep),'features'))
  print(paste("Hyper:",hyperNum, ";","Hypo:",hypoNum))
  list(
    diff=diff,
    hyper=hyper,
    hypo=hypo,
    res=res
  )
}

darCTLvsAIS<-findDAR('CTL','AIS')
darCTLvsMIA<-findDAR('CTL','MIA')
darCTLvsIAC<-findDAR('CTL','IAC')

darDeseq2<-list(
  darCTLvsAIS=darCTLvsAIS,
  darCTLvsMIA=darCTLvsMIA,
  darCTLvsIAC=darCTLvsIAC
)
saveRDS(darDeseq2, file.path(CONFIG$dataIntermediate,'atac', 'darDeseq2.rds'))
#----------------------------------------------------------------------------------------------------------------------
# Figure 3B. DAR numbers and Upset plot
#----------------------------------------------------------------------------------------------------------------------
darDeseq2<-readRDS(file.path(CONFIG$dataIntermediate,'atac', 'darDeseq2.rds'))
dar<-list(
  AISHyper=rownames(darDeseq2$darCTLvsAIS$hyper),
  AISHypo=rownames(darDeseq2$darCTLvsAIS$hypo),
  MIAHyper=rownames(darDeseq2$darCTLvsMIA$hyper),
  MIAHypo=rownames(darDeseq2$darCTLvsMIA$hypo),
  IACHyper=rownames(darDeseq2$darCTLvsIAC$hyper),
  IACHypo=rownames(darDeseq2$darCTLvsIAC$hypo)
)
sapply(names(dar), function(x){
  sapply(names(dar), function(y){
    length(intersect(dar[[x]],dar[[y]]))
  })
})
darDeseq2<-readRDS(file.path(CONFIG$dataIntermediate,'atac', 'darDeseq2.rds'))
m = make_comb_mat(dar)
setOrder=c("AISHyper", "AISHypo", "MIAHyper", "MIAHypo","IACHyper", "IACHypo")
saveImage2("atac.dar.upset.pdf",width = 4,height = 3)
UpSet(m, set_order = setOrder, comb_order = order(-comb_size(m)))
dev.off()
#----------------------------------------------------------------------------------------------------------------------
# SRDAR
#----------------------------------------------------------------------------------------------------------------------
darDeseq2<-readRDS(file.path(CONFIG$dataIntermediate,'atac', 'darDeseq2.rds'))
extractDar<-function(diff,group){
  data.frame(
    feature2Bed(rownames(diff)),
    class=ifelse(diff$log2FoldChange>0,'Hyper','Hypo'),
    group=group
  )
}
lData<-Reduce(rbind,list(
  AIS=extractDar(darDeseq2$darCTLvsAIS$diff,'AIS'),
  MIA=extractDar(darDeseq2$darCTLvsMIA$diff,'MIA'),
  IAC=extractDar(darDeseq2$darCTLvsIAC$diff,'IAC')
))
srdar<-dcast(lData,chrom + start+ end ~ group,value.var = 'class')
srdar[,4:6][is.na(srdar[,4:6])] <- "NC"

HyperInAIS<-Reduce(rbind,list(
  filter(srdar,AIS=='Hyper',MIA=="Hyper", IAC=="Hyper"),
  filter(srdar,AIS=='Hyper',MIA=="Hyper", IAC=="NC"),
  filter(srdar,AIS=='Hyper',MIA=="NC", IAC=="NC")
))%>%mutate(class='HyperInAIS')%>%dplyr::select(chrom, start,end,class)
HyperInMIA<-Reduce(rbind,list(
  filter(srdar,AIS=='NC',MIA=="Hyper", IAC=="NC"),
  filter(srdar,AIS=='NC',MIA=="Hyper", IAC=="Hyper")
))%>%mutate(class='HyperInMIA')%>%dplyr::select(chrom, start,end,class)
HyperInIAC<-Reduce(rbind,list(
  filter(srdar,AIS=='NC',MIA=="NC", IAC=="Hyper")
))%>%mutate(class='HyperInIAC')%>%dplyr::select(chrom, start,end,class)
HypoInAIS<-Reduce(rbind,list(
  filter(srdar,AIS=='Hypo',MIA=="Hypo", IAC=="Hypo"),
  filter(srdar,AIS=='Hypo',MIA=="Hypo", IAC=="NC"),
  filter(srdar,AIS=='Hypo',MIA=="NC", IAC=="NC")
))%>%mutate(class='HypoInAIS')%>%dplyr::select(chrom, start,end,class)
HypoInMIA<-Reduce(rbind,list(
  filter(srdar,AIS=='NC',MIA=="Hypo", IAC=="NC"),
  filter(srdar,AIS=='NC',MIA=="Hypo", IAC=="Hypo")
))%>%mutate(class='HypoInMIA')%>%dplyr::select(chrom, start,end,class)
HypoInIAC<-Reduce(rbind,list(
  filter(srdar,AIS=='NC',MIA=="NC", IAC=="Hypo")
))%>%mutate(class='HypoInIAC')%>%dplyr::select(chrom, start,end,class)
SRDAR<-Reduce(rbind,list(
  HyperInAIS=HyperInAIS,
  HyperInMIA=HyperInMIA,
  HyperInIAC=HyperInIAC,
  HypoInAIS=HypoInAIS,
  HypoInMIA=HypoInMIA,
  HypoInIAC=HypoInIAC
))
SRDAR$chrom <- factor(SRDAR$chrom, levels=chromFactorLevel)
SRDAR<-arrange(SRDAR, chrom, start)
srdar.count<-dplyr::count(srdar,AIS,IAC,MIA)%>%dplyr::select(AIS, MIA, IAC, n)%>%arrange(AIS, MIA, IAC)
SRDAR.count<-dplyr::count(SRDAR, class)%>%arrange(desc(n))
saveTsv(SRDAR, file.path(CONFIG$dataIntermediate, 'atac','srdar.bed'),col.names=FALSE)
saveTsv(filter(SRDAR,class=='HyperInAIS')[,1:3], file.path(CONFIG$dataIntermediate, 'atac','srdar.HyperInAIS.bed'),col.names = FALSE)
saveTsv(filter(SRDAR,class=='HypoInAIS')[,1:3], file.path(CONFIG$dataIntermediate, 'atac','srdar.HypoInAIS.bed'),col.names = FALSE)
saveTsv(filter(SRDAR,class=='HyperInMIA')[,1:3], file.path(CONFIG$dataIntermediate, 'atac','srdar.HyperInMIA.bed'),col.names = FALSE)
saveTsv(filter(SRDAR,class=='HypoInMIA')[,1:3], file.path(CONFIG$dataIntermediate, 'atac','srdar.HypoInMIA.bed'),col.names = FALSE)
saveTsv(filter(SRDAR,class=='HyperInIAC')[,1:3], file.path(CONFIG$dataIntermediate, 'atac','srdar.HyperInIAC.bed'),col.names = FALSE)
saveTsv(filter(SRDAR,class=='HypoInIAC')[,1:3], file.path(CONFIG$dataIntermediate, 'atac','srdar.HypoInIAC.bed'),col.names = FALSE)
#----------------------------------------------------------------------------------------------------------------------
# Figure 3C-E. SRDAR heatmap
#----------------------------------------------------------------------------------------------------------------------
SRDAR<-loadSRDAR()
atacPeakTPM<-readRDS(file.path(CONFIG$dataIntermediate,'atac', 'atacPeakTPM.rds'))
dar<-split(SRDAR,SRDAR$class)
plot.srdar<-function(region){
  data<-atacPeakTPM[match(bed2Feature(region),bed2Feature(atacPeakTPM)),-(1:3)]
  samples<-groups$ATAC$selectBySample(colnames(data))
  column_annotation <-HeatmapAnnotation(
    df=data.frame(Stage=samples$Group),
    col = list(Stage =colorMapStage),
    show_annotation_name =FALSE,
    annotation_name_side='left'
  )
  m1<-data
  nrow1<-nrow(m1)
  if (nrow(m1) >=1000){
    set.seed(123)
    m1<-m1[sample(1:nrow(m1),1000),]
  }
  nrow2<-nrow(m1)
  print(sprintf("%d => %d", nrow1, nrow2))
  Heatmap(t(scale(t(m1))),
          top_annotation = column_annotation,
          cluster_rows=TRUE,
          cluster_columns = FALSE,
          show_row_names=FALSE,
          show_column_names=FALSE,
          heatmap_legend_param = list(
            title = "Scaled Accessibility",
            legend_height = unit(4, "cm"),
            at = c(0,0.5,1),
            labels = c('0','0.5','1'),
            title_position = "lefttop-rot"
          ),
  )
}
saveImage2("atac.srdar.heatmap.HyperInAIS.pdf",width = 3.5,height = 2)
plot.srdar(dar$HyperInAIS)
dev.off()
saveImage2("atac.srdar.heatmap.HypoInAIS.pdf",width = 3.5,height = 2)
plot.srdar(dar$HypoInAIS)
dev.off()
saveImage2("atac.srdar.heatmap.HyperInMIA.pdf",width = 3.5,height = 2)
plot.srdar(dar$HyperInMIA)
dev.off()
saveImage2("atac.srdar.heatmap.HypoInMIA.pdf",width = 3.5,height = 2)
plot.srdar(dar$HypoInMIA)
dev.off()
saveImage2("atac.srdar.heatmap.HyperInIAC.pdf",width = 3.5,height = 2)
plot.srdar(dar$HyperInIAC)
dev.off()
saveImage2("atac.srdar.heatmap.HypoInIAC.pdf",width = 3.5,height = 2)
plot.srdar(dar$HypoInIAC)
dev.off()

#----------------------------------------------------------------------------------------------------------------------
# Figure 3F. SRDAR Density near TSS and CGI
#----------------------------------------------------------------------------------------------------------------------
genomicRegion<-readRDS(file.path(CONFIG$dataIntermediate, 'genomicRegion.rds'))
cgIslands.gr<-genomicRegion$cgIslands
tss.gr<-genomicRegion$tss
SRDAR<-loadSRDAR()
DAR<-split(SRDAR,SRDAR$class)
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
distToCGI<-annoDistQueryToSubject(bed2GRanges(SRDAR),cgIslands.gr)
distToTSS<-annoDistQueryToSubject(bed2GRanges(SRDAR),tss.gr)

set.seed(123)
lim<-50000
dist<-distToTSS
dist<-dist[abs(dist$distanceToSubject)<lim,]
dens<-lapply(split(dist, dist$class), function(x){
  density(x$distanceToSubject,bw = 'nrd0')
})
xlim<-range(sapply(dens, "[", "x"))
ylim<-range(sapply(dens, "[", "y"))

plot.dar.density.tss.vs.cgi<-function(g) {
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
  colorMap<-c("#0000FF80","#FF000080")
  names(colorMap)<-c("TSS","CGI")
  plot(NA, main = g,
       xlab = "Distance to TSS/CGI",
       ylab = "DARs Density",
       xlim=xlim,
       ylim=ylim,
       col = "blue")
  lines(dens$a, lwd = 2,col=colorMap[1])
  lines(dens$b, lwd = 2,col=colorMap[2])
  legend("topleft", legend=names(colorMap), fill=colorMap, bty = "n")
}

saveImage2("atac.dar.density.genomicRegion.pdf",width = 4,height = 6.5)
par(mfrow = c(3,2))
plot.dar.density.tss.vs.cgi('HyperInAIS')
plot.dar.density.tss.vs.cgi('HypoInAIS')
plot.dar.density.tss.vs.cgi('HyperInMIA')
plot.dar.density.tss.vs.cgi('HypoInMIA')
plot.dar.density.tss.vs.cgi('HyperInIAC')
plot.dar.density.tss.vs.cgi('HypoInIAC')
dev.off()
#----------------------------------------------------------------------------------------------------------------------
# SRDAR in genomic regions
#----------------------------------------------------------------------------------------------------------------------
SRDAR<-loadSRDAR()
DAR<-split(SRDAR, SRDAR$class)
genomicRegion<-readRDS(file.path(CONFIG$dataIntermediate, 'genomicRegion.rds'))
DAR.genomicRegion<-sapply(names(DAR),function(x){
  gr2<-bed2GRanges(DAR[[x]])
  sapply(genomicRegion, function(gr1){
    sum(countOverlaps(gr2, gr1)!=0)
  })
})
plot.dar.barplot.genomicRegion <- function(regions){
  allRegion<-data.frame(DAR.genomicRegion,check.names = FALSE)
  allRegion$region<-rownames(allRegion)
  allRegion$region<-factor(regions[match(allRegion$region,names(regions))],levels = regions)
  allRegion<-allRegion[!is.na(allRegion$region),]
  df<-melt(allRegion,variable.name='Stage')
  df$Stage<-factor(as.vector(df$Stage), levels=names(colorMapGroup))
  ggplot(data=df, aes(y=value, x=region, fill=Stage))+
    scale_fill_manual(values=colorMapGroup)+
    geom_bar(stat="identity", position=position_dodge())+
    ylab("DARs Number")+
    xlab("Genomic Regions")+
    theme_classic()
}
saveImage2("atac.dar.barplot.genomicRegion.cpg.pdf",width = 6,height = 2.5)
regions<-c('cgIslands'='CpG Islands', 'cgShores'='CpG Shores', 'cgShelves'='CpG Shelves', 'cgSea'='CpG Sea')
plot.dar.barplot.genomicRegion(regions)
dev.off()
saveImage2("atac.dar.barplot.genomicRegion.promoter.pdf",width = 8,height = 2.5)
regions<-c('promoter.1k'='Promoter.1k', 'promoter.5k'='Promoter.5k', 'utr5'="5'UTR", 'utr3'="3'UTR",'exons'='Exons','intron'='Intron','intergenic'='Intergenic')
plot.dar.barplot.genomicRegion(regions)
dev.off()
#----------------------------------------------------------------------------------------------------------------------
# Figure 3G,H. GREAT analysis of SRDARs
#----------------------------------------------------------------------------------------------------------------------
saveImage2("dar.great.AIS.hypo.BP.pdf",width = 4,height = 3)
plot.great(file.path(CONFIG$dataIntermediate, 'atac',"srdar.HypoInAIS.GOBiologicalProcess.tsv"),title="")
dev.off()
saveImage2("dar.great.IAC.hyper.HPO.pdf",width = 4,height = 6)
plot.great(file.path(CONFIG$dataIntermediate, 'atac',"srdar.HyperInIAC.GOBiologicalProcess.tsv"),title="")
dev.off()
#----------------------------------------------------------------------------------------------------------------------
# Figure S2A. ATAC-seq: DAR Homer Figure. heatmap
#----------------------------------------------------------------------------------------------------------------------
atacPeakTPM<-readRDS(file.path(CONFIG$dataIntermediate,'atac', 'atacPeakTPM.rds'))
atacPeakTPMFeatures<-bed2Feature(atacPeakTPM)
darDeseq2<-readRDS(file.path(CONFIG$dataIntermediate,'atac', 'darDeseq2.rds'))
plot.dar.status<-function(dar,dar.type){
  features<-rownames(dar)
  bed<-feature2Bed(features)
  tpmSelect<-atacPeakTPM[match(features, atacPeakTPMFeatures),]
  tpm<-groups$WGBS.ATAC$pickColumnsByGroup(c('IAC'), tpmSelect)
  data<-data.frame(
    x=bed$end - bed$start,
    y=dar$log2FoldChange
  )
  p_scatter <- ggplot(data, aes(x = x, y = y)) +
    geom_point(size=0.3,alpha = 0.6) +
    xlab("DARs size (bp)")+
    ylab("Log2FoldChange of DARs")+
    ggtitle(dar.type)+
    theme_bw()
  p<-ggMarginal(p_scatter, type="histogram",fill='white',bins = 100,size = 8)
  p
}

p1<-plot.dar.status(darDeseq2$darCTLvsAIS$hyper, 'HyperDARs In AIS')
p2<-plot.dar.status(darDeseq2$darCTLvsAIS$hypo, 'HypoDARs In AIS')
p3<-plot.dar.status(darDeseq2$darCTLvsMIA$hyper, 'HyperDARs In MIA')
p4<-plot.dar.status(darDeseq2$darCTLvsMIA$hypo, 'HypoDARs In MIA')
p5<-plot.dar.status(darDeseq2$darCTLvsIAC$hyper, 'HyperDARs In IAC')
p6<-plot.dar.status(darDeseq2$darCTLvsIAC$hypo, 'HypoDARs In IAC')
saveImage2("dar.status.scatter.logfc.length.pdf",width = 6,height = 9)
grid.arrange(p1, p2,p3,p4, p5,p6,nrow = 3)
dev.off()
#----------------------------------------------------------------------------------------------------------------------
# Table S11. ATAC-seq: DAR Homer
#----------------------------------------------------------------------------------------------------------------------
HypoInAIS=homerKnownTFs(file.path(CONFIG$dataIntermediate, 'atac','homer.mask','AIS.hypo'),'HypoInAIS')
HyperInMIA=homerKnownTFs(file.path(CONFIG$dataIntermediate, 'atac','homer.mask','MIA.hyper'),'HyperInMIA')
HyperInIAC=homerKnownTFs(file.path(CONFIG$dataIntermediate, 'atac','homer.mask','IAC.hyper'),'HyperInIAC')
HypoInIAC=homerKnownTFs(file.path(CONFIG$dataIntermediate, 'atac','homer.mask','IAC.hypo'),'HypoInIAC')
tfs<-do.call(rbind,list(HypoInAIS,HyperInMIA,HyperInIAC, HypoInIAC))
srdarTFs<-getSRTFS(tfs)
saveRDS(srdarTFs,file.path(CONFIG$dataIntermediate,'atac', 'srdar.tfs.rds'))
#----------------------------------------------------------------------------------------------------------------------
# Figure S2B. SRDAR Homer Figure. heatmap
#----------------------------------------------------------------------------------------------------------------------
srdarTFs<-readRDS(file.path(CONFIG$dataIntermediate,'atac', 'srdar.tfs.rds'))
TFS<-dplyr::arrange(srdarTFs$stage, desc(HyperInAIS),desc(HypoInAIS),desc(HyperInMIA),desc(HypoInMIA),desc(HyperInIAC),desc(HypoInIAC))
m<-as.matrix(TFS)
column_annotation <-HeatmapAnnotation(
  df=data.frame(SRDAR=factor(colnames(m), levels = colnames(m))),
  col = list(SRDAR =colorMapGroup),
  show_annotation_name =FALSE,
  annotation_name_side='left'
)
saveImage2("atac.dar.homer.tfs.heatmap.pdf",width = 3,height = 11)
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
# Figure S2C. ATAC-seq: DAR Homer Figure. heatmap count
#----------------------------------------------------------------------------------------------------------------------
srdarTFs<-readRDS(file.path(CONFIG$dataIntermediate,'atac', 'srdar.tfs.rds'))

data<-srdarTFs$stage
mat<-sapply(names(data), function(x){
  sapply(names(data), function(y){
    sum(data[[x]]&data[[y]])
  })
})
mat[upper.tri(mat)] <- NA

saveImage2("atac.dar.homer.tfs.heatmap.count.pdf",width = 4,height = 4)
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




