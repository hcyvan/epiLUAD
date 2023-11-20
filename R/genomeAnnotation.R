source('./R/base.R')
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(dplyr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

fixGranges<- function(gr,out.format='GRanges') {
  gr<-gr[seqnames(gr) %in% chromFactorLevel]
  seqlevels(gr) <- chromFactorLevel
  gr<-sort(gr, by = ~seqnames + start)
  if (out.format=='bed3') {
    data.frame(
      chrom = seqnames(gr),
      start = start(gr)-1,
      end = end(gr)
    )
  }else if(out.format=='bed4') {
    data.frame(
      chrom = seqnames(gr),
      start = start(gr)-1,
      end = end(gr),
      strand = strand(gr)
    )
  }else {
    gr
  }
}
#----------------------------------------------------------------------------------------------------------------------
# Get the region of CpG Islands, CpG Shores, CpG Shelves and CpG Sea.
#----------------------------------------------------------------------------------------------------------------------
cgiPath<-file.path(CONFIG$dataDir,'cpgIslandExt.hg38.ucsc.bed')

cpgIslandExt<-loadData('cpgIslandExt.hg38.ucsc', ext='bed', header=FALSE)
cpgIslandExt<-cpgIslandExt[cpgIslandExt$V1%in%chromFactorLevel,]

chromInfo<-loadData('chromInfo.hg38.ucsc',ext='tsv')
chromInfo<-chromInfo[chromInfo$`#chrom`%in%chromFactorLevel,]

cgIslands <- GRanges(
  seqnames = cpgIslandExt$V1,
  ranges = IRanges(start = cpgIslandExt$V2+1, end = cpgIslandExt$V3),
)
cgShoresPre<-resize(cgIslands, width = width(cgIslands) + 2000*2, fix = "center")
cgShelvesPre<-resize(cgIslands, width = width(cgIslands) + 4000*2, fix = "center")
cgSeaPre<-GRanges(
  seqnames = chromInfo$`#chrom`,
  ranges = IRanges(start = 1, end = chromInfo$size),
)
cgShores<-setdiff(cgShoresPre, cgIslands)
cgShelves<-setdiff(cgShelvesPre, cgShoresPre)
cgSea<-setdiff(cgSeaPre, cgShelvesPre)




saveRDS(fixGranges(cgIslands), file.path(CONFIG$dataIntermediate, 'cgIslands.rds'))
saveRDS(fixGranges(cgShores), file.path(CONFIG$dataIntermediate, 'cgShores.rds'))
saveRDS(fixGranges(cgShelves), file.path(CONFIG$dataIntermediate, 'cgShelves.rds'))
saveRDS(fixGranges(cgSea), file.path(CONFIG$dataIntermediate, 'cgSea.rds'))

write.table(fixGranges(cgIslands,out.format='bed3'), file.path(CONFIG$dataAnnotation, 'cgIslands.bed'), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(fixGranges(cgShores,out.format='bed3'), file.path(CONFIG$dataAnnotation, 'cgShores.bed'), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(fixGranges(cgShelves,out.format='bed3'), file.path(CONFIG$dataAnnotation, 'cgShelves.bed'), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(fixGranges(cgSea,out.format='bed3'), file.path(CONFIG$dataAnnotation, 'cgSea.bed'), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#----------------------------------------------------------------------------------------------------------------------
# Get the region of TSS, Promoters, 5'UTR, 3'UTR, Introns and Exons
#----------------------------------------------------------------------------------------------------------------------
region.5utr<-unlist(fiveUTRsByTranscript(txdb))
region.3utr<-unlist(threeUTRsByTranscript(txdb))
region.intron<-unlist(intronsByTranscript(txdb))
region.exons<-unlist(exonsBy(txdb))
region.promoter.1k<-promoters(txdb, upstream=1000, downstream=1000)
region.promoter.5k<-promoters(txdb, upstream=5000, downstream=5000)
region.transcripts<-transcripts(txdb)
region.tss<-resize(region.transcripts, width=1, fix="start")



write.table(fixGranges(region.tss,out.format = 'bed4'), file.path(CONFIG$dataAnnotation, 'tss.TxDb.Hsapiens.UCSC.hg38.knownGene.bed'),row.names = FALSE,sep = '\t',quote = FALSE, col.names = FALSE)
write.table(fixGranges(region.5utr,out.format = 'bed4'), file.path(CONFIG$dataAnnotation, 'utr5.TxDb.Hsapiens.UCSC.hg38.knownGene.bed'),row.names = FALSE,sep = '\t',quote = FALSE, col.names = FALSE)
write.table(fixGranges(region.3utr,out.format = 'bed4'), file.path(CONFIG$dataAnnotation, 'utr3.TxDb.Hsapiens.UCSC.hg38.knownGene.bed'),row.names = FALSE,sep = '\t',quote = FALSE, col.names = FALSE)
write.table(fixGranges(region.intron,out.format = 'bed4'), file.path(CONFIG$dataAnnotation, 'intron.TxDb.Hsapiens.UCSC.hg38.knownGene.bed'),row.names = FALSE,sep = '\t',quote = FALSE, col.names = FALSE)
write.table(fixGranges(region.exons,out.format = 'bed4'), file.path(CONFIG$dataAnnotation, 'exons.TxDb.Hsapiens.UCSC.hg38.knownGene.bed'),row.names = FALSE,sep = '\t',quote = FALSE, col.names = FALSE)
write.table(fixGranges(region.promoter.1k,out.format = 'bed4'), file.path(CONFIG$dataAnnotation, 'promoter.1k.TxDb.Hsapiens.UCSC.hg38.knownGene.bed'),row.names = FALSE,sep = '\t',quote = FALSE, col.names = FALSE)
write.table(fixGranges(region.promoter.5k,out.format = 'bed4'), file.path(CONFIG$dataAnnotation, 'promoter.5k.TxDb.Hsapiens.UCSC.hg38.knownGene.bed'),row.names = FALSE,sep = '\t',quote = FALSE, col.names = FALSE)


