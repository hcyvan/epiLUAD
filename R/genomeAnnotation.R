source('./R/base.R')
library(rtracklayer)
library(GenomicRanges)

cgiPath<-file.path(CONFIG$dataDir,'cpgIslandExt.hg38.ucsc.bed')


cpgIslandExt<-loadData('cpgIslandExt.hg38.ucsc', ext='bed', header=FALSE)
cpgIslandExt<-cpgIslandExt[cpgIslandExt$V1%in%chromFactorLevel,]

chromInfo<-loadData('chromInfo.hg38.ucsc',ext='tsv')
chromInfo<-chromInfo[chromInfo$`#chrom`%in%chromFactorLevel,]


cgIslands <- GRanges(
  seqnames = cpgIslandExt$V1,
  ranges = IRanges(start = cpgIslandExt$V2, end = cpgIslandExt$V3),
)
cgShoresPre<-resize(cgIslands, width = width(cgIslands) + 2000*2, fix = "center")
cgShelvesPre<-resize(cgIslands, width = width(cgIslands) + 4000*2, fix = "center")
cgSeaPre<-GRanges(
  seqnames = chromInfo$`#chrom`,
  ranges = IRanges(start = 0, end = chromInfo$size),
)
cgShores<-setdiff(cgShoresPre, cgIslands)
cgShelves<-setdiff(cgShelvesPre, cgShoresPre)
cgSea<-setdiff(cgSeaPre, cgShelvesPre)

sortGranges<- function(gr) {
  seqlevels(gr) <- chromFactorLevel
  sort(gr, by = ~seqnames + start)

}
toBed3<-function(gr){
  data.frame(
    seqnames = seqnames(gr),
    start = start(gr) - 1,
    end = end(gr)
  )
}

cgIslands<-sortGranges(cgIslands)
cgShores<-sortGranges(cgShores)
cgShelves<-sortGranges(cgShelves)
cgSea<-sortGranges(cgSea)

saveRDS(cgIslands, file.path(CONFIG$dataIntermediate, 'cgIslands.rds'))
saveRDS(cgShores, file.path(CONFIG$dataIntermediate, 'cgShores.rds'))
saveRDS(cgShelves, file.path(CONFIG$dataIntermediate, 'cgShelves.rds'))
saveRDS(cgSea, file.path(CONFIG$dataIntermediate, 'cgSea.rds'))

write.table(toBed3(cgIslands), file.path(CONFIG$dataAnnotation, 'cgIslands.bed'), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(toBed3(cgShores), file.path(CONFIG$dataAnnotation, 'cgShores.bed'), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(toBed3(cgShelves), file.path(CONFIG$dataAnnotation, 'cgShelves.bed'), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(toBed3(cgSea), file.path(CONFIG$dataAnnotation, 'cgSea.bed'), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

