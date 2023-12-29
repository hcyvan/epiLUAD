source('./R/base.R')
library(ChIPseeker)
library(DESeq2)
library(reshape2)
library(ggplot2)
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
saveDAR<-function(dar,hypo.file,hyper.file) {
  hypo<-dar$hypo
  hyper<-dar$hyper
  cat(sprintf("hypo: %d; hyper: %d\n", dim(hypo)[1], dim(hyper)[1]))
  write.table(data.frame(feature2Bed(rownames(hypo)),hypo), file.path(CONFIG$dataIntermediate,'atac', hypo.file), sep='\t', row.names = FALSE, col.names = FALSE,quote = FALSE)
  write.table(data.frame(feature2Bed(rownames(hyper)),hyper), file.path(CONFIG$dataIntermediate,'atac', hyper.file), sep='\t', row.names = FALSE, col.names = FALSE,quote = FALSE)
}
saveDAR(darDeseq2$darCTLvsAIS,'dar.p400.AIS.hypo.bed','dar.p400.AIS.hyper.bed')
saveDAR(darDeseq2$darCTLvsMIA,'dar.p400.MIA.hypo.bed','dar.p400.MIA.hyper.bed')
saveDAR(darDeseq2$darCTLvsIAC,'dar.p400.IAC.hypo.bed','dar.p400.IAC.hyper.bed')
#----------------------------------------------------------------------------------------------------------------------
# Figure 3B. ATAC-seq: DAR numbers and Upset plot
#----------------------------------------------------------------------------------------------------------------------
darDeseq2<-readRDS(file.path(CONFIG$dataIntermediate,'atac', 'darDeseq2.rds'))
dar<-list(
  AISHyperDARs=rownames(darDeseq2$darCTLvsAIS$hyper),
  AISHypoDARs=rownames(darDeseq2$darCTLvsAIS$hypo),
  MIAHyperDARs=rownames(darDeseq2$darCTLvsMIA$hyper),
  MIAHypoDARs=rownames(darDeseq2$darCTLvsMIA$hypo),
  IACHyperDARs=rownames(darDeseq2$darCTLvsIAC$hyper),
  IACHypoDARs=rownames(darDeseq2$darCTLvsIAC$hypo)
)
sapply(names(dar), function(x){
  sapply(names(dar), function(y){
    length(intersect(dar[[x]],dar[[y]]))
  })
})
darDeseq2<-readRDS(file.path(CONFIG$dataIntermediate,'atac', 'darDeseq2.rds'))
m = make_comb_mat(dar)
setOrder=c("AISHyperDARs", "AISHypoDARs", "MIAHyperDARs", "MIAHypoDARs","IACHyperDARs", "IACHypoDARs")
saveImage2("atac.dar.upset.pdf",width = 6,height = 4)
UpSet(m, set_order = setOrder, comb_order = order(-comb_size(m)))
dev.off()
#----------------------------------------------------------------------------------------------------------------------
# Figure 3C. ATAC-seq: DAR Density near TSS and CGI
#----------------------------------------------------------------------------------------------------------------------
genomicRegion<-readRDS(file.path(CONFIG$dataIntermediate, 'genomicRegion.rds'))
cgIslands.gr<-genomicRegion$cgIslands
tss.gr<-genomicRegion$tss
darDeseq2<-readRDS(file.path(CONFIG$dataIntermediate,'atac', 'darDeseq2.rds'))
dar<-list(
  AISHyperDARs=rownames(darDeseq2$darCTLvsAIS$hyper),
  AISHypoDARs=rownames(darDeseq2$darCTLvsAIS$hypo),
  MIAHyperDARs=rownames(darDeseq2$darCTLvsMIA$hyper),
  MIAHypoDARs=rownames(darDeseq2$darCTLvsMIA$hypo),
  IACHyperDARs=rownames(darDeseq2$darCTLvsIAC$hyper),
  IACHypoDARs=rownames(darDeseq2$darCTLvsIAC$hypo)
)

DAR<-do.call(rbind,lapply(names(dar), function(x){
  data.frame(feature2Bed(dar[[x]]),class=x)
}))

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
distToCGI<-annoDistQueryToSubject(bed2GRanges(DAR),cgIslands.gr)
distToTSS<-annoDistQueryToSubject(bed2GRanges(DAR),tss.gr)

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

saveImage2("atac.dar.density.genomicRegion.pdf",width = 5,height = 6)
par(mfrow = c(3,2))
plot.dar.density.tss.vs.cgi('AISHyperDARs')
plot.dar.density.tss.vs.cgi('AISHypoDARs')
plot.dar.density.tss.vs.cgi('MIAHyperDARs')
plot.dar.density.tss.vs.cgi('MIAHypoDARs')
plot.dar.density.tss.vs.cgi('IACHyperDARs')
plot.dar.density.tss.vs.cgi('IACHypoDARs')
dev.off()

#----------------------------------------------------------------------------------------------------------------------
# Figure 3D. ATAC-seq: DAR in genomic regions
#----------------------------------------------------------------------------------------------------------------------
darDeseq2<-readRDS(file.path(CONFIG$dataIntermediate,'atac', 'darDeseq2.rds'))
DAR<-list(
  AISHyperDARs=rownames(darDeseq2$darCTLvsAIS$hyper),
  AISHypoDARs=rownames(darDeseq2$darCTLvsAIS$hypo),
  MIAHyperDARs=rownames(darDeseq2$darCTLvsMIA$hyper),
  MIAHypoDARs=rownames(darDeseq2$darCTLvsMIA$hypo),
  IACHyperDARs=rownames(darDeseq2$darCTLvsIAC$hyper),
  IACHypoDARs=rownames(darDeseq2$darCTLvsIAC$hypo)
)
genomicRegion<-readRDS(file.path(CONFIG$dataIntermediate, 'genomicRegion.rds'))
DAR.genomicRegion<-sapply(names(DAR),function(x){
  gr2<-bed2GRanges(feature2Bed(DAR[[x]]))
  sapply(genomicRegion, function(gr1){
    sum(countOverlaps(gr2, gr1)!=0)
  })
})
plot.dar.barplot.genomicRegion <- function(regions){
  allRegion<-data.frame(DAR.genomicRegion,check.names = FALSE)
  allRegion$region<-rownames(allRegion)
  allRegion$region<-factor(regions[match(allRegion$region,names(regions))],levels = regions)
  allRegion<-allRegion[!is.na(allRegion$region),]
  # print(allRegion)
  # allRegion[,1:6]<-log2(allRegion[,1:6]+1)
  # print(allRegion)
  df<-melt(allRegion,variable.name='Stage')
  print(df)
  ggplot(data=df, aes(y=value, x=region, fill=Stage))+
    scale_fill_manual(values=colorMapDAR)+
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
# write.csv(data.frame(t(DAR.genomicRegion)), file.path(CONFIG$dataResult, 'atac.dar.count.genomicRegion.csv'), quote = FALSE)
#----------------------------------------------------------------------------------------------------------------------
# Figure 3E,F. GREAT analysis of DMRs
#----------------------------------------------------------------------------------------------------------------------
plot.great<-function(tsv,title="") {
  enrich<-loadData2(tsv, comment.char = "#",header=FALSE)
  enrich.col<-c("TermName","BinomRank","BinomRawPValue","BinomFDRQVal","BinomFoldEnrichment","BinomObservedRegionHits","BinomRegionSetCoverage","HyperRank","HyperFDRQVal","HyperFoldEnrichment","HyperObservedGeneHits","HyperTotalGenes","HyperGeneSetCoverage")
  colnames(enrich)<-enrich.col
  data<-data.frame(x=enrich$TermName,y=-log10(enrich$BinomFDRQVal))
  data<-data[order(data$y),]
  ggplot(data, aes(x=x, y=y))+
    geom_bar(stat="identity", width=0.3, fill='black')+
    coord_flip()+
    scale_x_discrete(limits=data$x,labels = NULL )+
    theme_classic()+
    theme(legend.position="none")+
    scale_y_continuous(expand = c(0,0))+
    annotate("text", x=seq(1, nrow(data))+0.4, y=0,
             hjust = 0, cex=3,
             label= data$x)+
    labs(x="", y="-Log10(q-value)",title=title)+
    theme( axis.ticks.y = element_blank(),
           axis.line.y = element_blank(),
           axis.text.y = element_blank())
}
saveImage2("dar.great.AIS.hypo.BP.pdf",width = 4,height = 4)
plot.great(file.path(CONFIG$dataIntermediate, 'atac',"dar.p400.AIS.hypo.great.GOBiologicalProcess.tsv"),title="")
dev.off()
saveImage2("dar.great.AIS.hypo.CC.pdf",width = 4,height = 3)
plot.great(file.path(CONFIG$dataIntermediate, 'atac',"dar.p400.AIS.hypo.great.GOCellularComponent.tsv"),title="")
dev.off()
saveImage2("dar.great.MIA.hyper.HPO.pdf",width = 4,height = 1.2)
plot.great(file.path(CONFIG$dataIntermediate, 'atac',"dar.p400.MIA.hyper.great.HumanPhenotypeOntology.tsv"),title="")
dev.off()
saveImage2("dar.great.IAC.hyper.BP.pdf",width = 4,height = 6)
plot.great(file.path(CONFIG$dataIntermediate, 'atac',"dar.p400.IAC.hyper.great.GOBiologicalProcess.tsv"),title="")
dev.off()
saveImage2("dar.great.IAC.hyper.CC.pdf",width = 4,height = 1.2)
plot.great(file.path(CONFIG$dataIntermediate, 'atac',"dar.p400.IAC.hyper.great.GOCellularComponent.tsv"),title="")
dev.off()
saveImage2("dar.great.IAC.hyper.HPO.pdf",width = 4,height = 1.5)
plot.great(file.path(CONFIG$dataIntermediate, 'atac',"dar.p400.IAC.hyper.great.HumanPhenotypeOntology.tsv"),title="")
dev.off()
saveImage2("dar.great.IAC.hypo.BP.pdf",width = 4,height = 1.2)
plot.great(file.path(CONFIG$dataIntermediate, 'atac',"dar.p400.IAC.hypo.great.GOBiologicalProcess.tsv"),title="")
dev.off()
#----------------------------------------------------------------------------------------------------------------------
# ATAC-seq: DAR Homer
#----------------------------------------------------------------------------------------------------------------------
homerAISHyper=homerKnownTFs(file.path(CONFIG$dataIntermediate, 'atac','homer.mask','AIS.hyper'))
homerAISHypo=homerKnownTFs(file.path(CONFIG$dataIntermediate, 'atac','homer.mask','AIS.hypo'))
homerMIAHyper=homerKnownTFs(file.path(CONFIG$dataIntermediate, 'atac','homer.mask','MIA.hyper'))
homerMIAHypo=homerKnownTFs(file.path(CONFIG$dataIntermediate, 'atac','homer.mask','MIA.hypo'))
homerIACHyper=homerKnownTFs(file.path(CONFIG$dataIntermediate, 'atac','homer.mask','IAC.hyper'))
homerIACHypo=homerKnownTFs(file.path(CONFIG$dataIntermediate, 'atac','homer.mask','IAC.hypo'))

tfs<-do.call(rbind,list(
  # data.frame(tf=homerAISHyper, motif=names(homerAISHyper),class='AISHyperDARs'), # dim()[1]==0
  data.frame(tf=homerAISHypo, motif=names(homerAISHypo),class='AISHypoDARs'),
  data.frame(tf=homerMIAHyper, motif=names(homerMIAHyper),class='MIAHyperDARs'),
  data.frame(tf=homerMIAHypo, motif=names(homerMIAHypo),class='MIAHypoDARs'),
  data.frame(tf=homerIACHyper, motif=names(homerIACHyper),class='IACHyperDARs'),
  data.frame(tf=homerIACHypo, motif=names(homerIACHypo),class='IACHypoDARs')
))
tfs$class<-factor(tfs$class, levels = names(colorMapDAR))
motifWidthData<-dcast(tfs, motif~class,fun.aggregate = length)

motifTable<-data.frame(Motifs=motifWidthData$motif,
                       TFs=sapply(strsplit(motifWidthData$motif,'\\('),function(x){x[1]}),
                       AISHypoDARs=motifWidthData[,2],
                       AISHyperDARs=0,
                       motifWidthData[,3:ncol(motifWidthData)])
write.csv(motifTable, file.path(CONFIG$dataResult, 'atac.dar.homer.motifs.csv'),row.names  = FALSE)
#----------------------------------------------------------------------------------------------------------------------
# ATAC-seq: DAR Homer Figure. heatmap
#----------------------------------------------------------------------------------------------------------------------
tfWidthData.2<-data.frame(tfWidthData[,1:2],AISHyperDARs=0 ,tfWidthData[3:ncol(tfWidthData)])

TFS<-dplyr::arrange(tfWidthData.2, desc(AISHypoDARs),desc(AISHyperDARs),desc(MIAHypoDARs),desc(MIAHyperDARs),desc(IACHypoDARs),desc(IACHyperDARs))
m<-as.matrix(TFS[,2:ncol(TFS)])
rownames(m)<-TFS[,1]
m<-ifelse(m>0,1,0)

column_annotation <-HeatmapAnnotation(
  df=data.frame(DAR=names(colorMapDAR)),
  col = list(DAR =colorMapDAR),
  show_annotation_name =FALSE,
  annotation_name_side='left'
)
saveImage2("atac.dar.homer.tfs.heatmap.pdf",width = 3,height = 12)
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
        row_names_gp = grid::gpar(fontsize = 4)
)
dev.off()
#----------------------------------------------------------------------------------------------------------------------
# ATAC-seq: DAR Homer Figure. heatmap count
#----------------------------------------------------------------------------------------------------------------------
data<-tfWidthData.2[2:ncol(tfWidthData.2)]
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




