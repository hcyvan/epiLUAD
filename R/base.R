library(dplyr)
library(yaml)
library(tools)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(circlize)

configDefaultPath <- './config.default.yaml'
configPath <- './config.yaml'
if (file.exists(configPath)) {
  
} else if(file.exists(configDefaultPath)) {
  configPath <- configDefaultPath
} else {
  stop('please add configure file [config.yaml] to this project root directory')
}

CONFIG <- yaml::read_yaml(configPath)
CONFIG['dataAnnotation']<-'./data/annotation'
CONFIG['dataIntermediate']<-'./data/intermediate'
CONFIG['dataExternal']<-'./data/external'
CONFIG['dataResult']<-'./data/result'

DATA_DIR <- CONFIG$dataDir
IMAGE_DIR <- CONFIG$imageDir
########################################### helper function ###########################################
annoPeak<-function(peakGr) {
  peakAnno <- annotatePeak(peakGr, tssRegion = c(-5000, 5.000), TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene)
  peakAnno@anno$symbol<-mapIds(org.Hs.eg.db, keys = peakAnno@anno$geneId, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
  peakAnno@anno$ensembl<-mapIds(org.Hs.eg.db, keys = peakAnno@anno$geneId, column = "ENSEMBL", keytype = "ENTREZID", multiVals = "first")
  anno<-peakAnno@anno
  anno[abs(anno$distanceToTSS)<100000,]
}

writeBed<-function(bed,file.name){
  write.table(bed, file.name,row.names = FALSE,sep = '\t',quote = FALSE, col.names = FALSE)
}

mapEnsemble2Symbol<-function(ensemble){
  mapIds(org.Hs.eg.db, keys = ensemble, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
}

ensemble2Symbol2<-function(ensemble){
  ensemblMart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  # saveRDS(ensemblMart, file.path(CONFIG$dataIntermediate, 'biomaRtEnsemblMart.rds'))
  genes_info <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), 
                      filters = 'ensembl_gene_id', 
                      values = ensemble, 
                      mart = ensemblMart)
  symbols<-genes_info$external_gene_name
  names(symbols)<-genes_info$ensembl_gene_id
  symbols
}
ensemble2Symbol3<-function(ensemble,keepIfNotMatch=FALSE){
  ensembl2symbolMap<-readRDS(file.path(CONFIG$dataIntermediate, 'ensembl2symbolMap.rds'))
  df<-data.frame(
    a=ensemble,
    b=ensembl2symbolMap[match(ensemble,names(ensembl2symbolMap))]
  )
  sapply(1:nrow(df),function(i){
    if (is.na(df[i,2])){
      ifelse(keepIfNotMatch, df[i,1], "")
    }else if (nchar(df[i,2])==0) {
      ifelse(keepIfNotMatch, df[i,1], "")
    }else {
      df[i,2]
    }
  })
}

homerKnownTFs<-function(dir,stage='LUAD',cutoff=0.001){
  result.txt<-file.path(dir,'knownResults.txt')
  result<-read.csv(result.txt,sep='\t')
  result<-result[result$q.value..Benjamini.< 0.02,]
  result<-filter(result, q.value..Benjamini.<=cutoff)
  homer<-sapply(result$Motif.Name, function(x){
    c(strsplit(x,"\\(")[[1]][1])
  })
  out<-data.frame(tf=homer, motif=names(homer),class=stage)
  rownames(out)<-NULL
  out
}


fixHomerTFs<-function(tf0){
  tf<-toupper(tf0)
  for (s in c('-FUSION','-HALFSITE', '-DISTAL',':EBOX', ':E-BOX','-AML')){
    tf<-gsub(s,"",tf)
  }
  tfDimerMap<-c(
    "JUN-AP1"='JUN',
    'C-JUN-CRE'='JUN',
    'EWS:ERG'='ERG',
    'EWS:FLI1'='FLI1',
    'NF1:FOXA1'='NF1',
    'PU.1-IRF'='PU.1',
    'PU.1:IRF8'='IRF8',
    'ARNT:AHR'='AHR',
    'ETS:RUNX'='ETS'
  )
  tfMap<-c(
    'PU.1'='SPI1',
    'NKX2.1'='NKX2-1',
    'NKX2.2'='NKX2-2',
    'NKX2.5'='NKX2-5',
    'NKX3.1'='NKX3-1',
    'NKX6.1'='NKX6-1',
    'BAPX1'='NKX3-2',
    'FRA1'='FOSL1',
    'FRA2'='FOSL2',
    'TLX?'='NR2E1',
    'CHOP'='DDIT3',
    'E2A'='TCF3',
    'HNF1'='HNF1A',
    'NF-E2'='NFE2',
    'NRF2'='NFE2L2',
    'P53'='TP53',
    'P63'='TP63',
    'P73'='TP73',
    'SNAIL1'='SNAI1',
    'SLUG'='SNAI2',
    'EWS'='EWSR1',
    'EKLF'='KLF1',
    'AP-2GAMMA'='TFAP2C',
    'AP-2ALPHA'='TFAP2A',
    'HIF-1B'='ARNT'
  )
  tfRemove<-c('AP-1','CRE', 'ETS', 'RUNX', 'FOX', 'HEB', 'TEAD', 'ZFP809','ETS:RUNX','FOX:EBOX','RUNX-AML')
  tf<-sapply(tf, function(x){
    ifelse(x%in%names(tfDimerMap), tfDimerMap[match(x, names(tfDimerMap))],x)
  })
  tf<-sapply(tf, function(x){
    ifelse(x%in%names(tfMap), tfMap[match(x, names(tfMap))],x)
  })
  tf<-sapply(tf, function(x){
    ifelse(x%in%tfRemove, NA,x)
  })
  tf
}

getSRTFS<-function(tfs){
  tfGeneMotif<-data.frame(gene=fixHomerTFs(tfs$tf),tf=tfs$tf,motif=tfs$motif)
  tfGeneMotif<-filter(tfGeneMotif,!is.na(gene))%>%distinct(gene,tf,motif)
  statusStage<-c("HyperInAIS","HypoInAIS", "HyperInMIA","HypoInMIA", "HyperInIAC", "HypoInIAC")
  tfStage<-sapply(split(tfGeneMotif, tfGeneMotif$gene), function(x){
    data<-left_join(x, tfs, by='motif')
    ifelse(statusStage%in%data$class,1,0)
  })%>%t%>%data.frame()
  colnames(tfStage)<-statusStage
  srTFs<-list(
    tf=rownames(tfStage),
    map=tfGeneMotif,
    stage=tfStage
  )
  srTFs
}

feature2Bed<-function(feature) {
  bed<-do.call(rbind,lapply(feature, function(x){
    tmp<-strsplit(x,':')[[1]]
    se<-strsplit(tmp[2],'-')[[1]]
    c(tmp[1], se[1], se[2])
  }))
  bed<-data.frame(bed)
  bed[,2]<-as.numeric(bed[,2])
  bed[,3]<-as.numeric(bed[,3])
  colnames(bed)<-c('chrom','start','end')
  bed
}

feature2GRanges<-function(feature){
  bed2GRanges(feature2Bed(feature))
}

bed2Feature<-function(bed){
  paste(paste(bed[,1], bed[,2], sep=':'),bed[,3],sep='-')
}

bed2GRanges<-function(bed){
  gr<-GRanges(seqnames = bed[,1], ranges = IRanges(start = bed[,2]+1, end =  bed[,3]))
  if (ncol(bed)>=4){
    .<-lapply(colnames(bed)[4:ncol(bed)],function(x){
      mcols(gr)[[x]]<<-bed[[x]]
    })
  }
  gr
}

GRanges2bed<-function(gr){
  data.frame(chrom=seqnames(gr), start=start(gr)-1, end=end(gr))
}
GRanges2Feature<-function(gr){
  bed<-data.frame(chrom=seqnames(gr), start=start(gr)-1, end=end(gr))
  bed2Feature(bed)
}
percent2Numeric <- function(x){
  as.numeric(substr(x,0,nchar(x)-1))/100
}

printf <- function(format, ...){
  cat(sprintf(format, ...))
}

removeNegativeOne <- function(m){
  m[rowSums(m[,4:ncol(m)]==-1)==0,]
}

loadData <- function(name, ext="csv", header=FALSE, force.refresh=FALSE){
  ext.path <- paste(file.path(DATA_DIR, name),ext,sep = '.')
  rds.path <- paste(file.path(DATA_DIR, name),'rds',sep = '.')
  rds.external.path <- paste(file.path(DATA_DIR, 'external', name),'rds',sep = '.')
  if(file.exists(rds.path) && !force.refresh){
    readRDS(rds.path)
  }else if(file.exists(ext.path)) {
    if(ext=='csv'){
      data<-read.csv(ext.path)
      saveRDS(data, rds.path)
      data
    }else if(ext=='tsv'){
      data<-read.csv(ext.path,sep = '\t',check.names = FALSE)
      saveRDS(data, rds.path)
      data
    }else if(ext=='bed'){
      data<-read.csv(ext.path,sep = '\t',check.names = FALSE, header=header)
      saveRDS(data, rds.path)
      data
    } else {
      stop(sprintf("The file format <%s> is not supported", ext)); 
    }
  } else if(file.exists(rds.external.path)){
    readRDS(rds.external.path)
  } else{
    stop(sprintf("File not exist: %s and %s", rds.path, ext.path))
  }
}

loadDataBed<-function(name){
  loadData(name, 'bed')
}

#' load data function
#' 
#' This function while load data and generate rds cache.
#' 
#' @param filename The file path of input data
#' @param file.format The file format. If not set, it will be automatically set based on the extension.
#' @param force.refresh Set TRUE will over write the rds cache
#' @param header see read.csv
#'
loadData2<-function(filename, file.format=NULL, force.refresh=FALSE, header=TRUE,comment.char = "") {
  filename.rds<-paste0(file_path_sans_ext(filename),'.rds')
  if (is.null(file.format)){
    file.format=file_ext(filename)
  }
  
  if (file.exists(filename.rds) && !force.refresh) {
    readRDS(filename.rds)
  } else {
    if(file.format=='bed'|file.format=='tsv'){
      data<-read.csv(filename,sep = '\t',check.names = FALSE, header=header, comment.char=comment.char)
      saveRDS(data, filename.rds)
      data
    } else if(file.format=='csv'){
      data<-read.csv(filename,check.names = FALSE, header=header)
      saveRDS(data, filename.rds)
      data
    }
  }
}
########################################### plot function ###########################################
plotEnrich<-function(term, qval,title=""){
  data<-data.frame(x=term,y=-log10(qval))
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
plot.great<-function(tsv,cutoff=0.01,title="") {
  enrich<-loadData2(tsv, comment.char = "#", header=FALSE)
  enrich.col<-c("TermName","BinomRank","BinomRawPValue","BinomFDRQVal","BinomFoldEnrichment","BinomObservedRegionHits","BinomRegionSetCoverage","HyperRank","HyperFDRQVal","HyperFoldEnrichment","HyperObservedGeneHits","HyperTotalGenes","HyperGeneSetCoverage")
  colnames(enrich)<-enrich.col
  enrich<-enrich[,-ncol(enrich)]
  enrich<-filter(enrich, BinomFDRQVal<=cutoff)
  plotEnrich(enrich$TermName, enrich$BinomFDRQVal, title)
}
plotBarError<-function(valueStage, xlabel="",ylabel="TPM"){
  ggplot(data=valueStage,aes(x=stage,y=value,fill=stage))+
    scale_fill_manual(values=colorMapStage) +
    stat_summary(mapping=aes(fill = stage),fun=mean,geom = "bar",fun.args = list(mult=1),width=0.7)+
    stat_summary(fun.data=mean_sdl,fun.args = list(mult=0.3),geom="errorbar",width=0.2) +
    labs(x=xlabel,y=ylabel)+
    theme_classic()+
    theme(legend.position="none",
          axis.title.x = element_text(size=15),
          axis.title.y = element_text(size=15),
          axis.text = element_text(size = 12,colour="black"))+
    guides(colour = guide_legend(override.aes = list(shape = 12,size=10)))
}
plotDot<-function(x,y, col,xlab,ylab,text.title=NA, withTest=FALSE){
  model <- lm(y ~ x)
  plot(x, y,col=col,pch=20,xlab=xlab,ylab=ylab)
  abline(model, col="gray3")
  if(!is.na(text.title)){
    title(text.title)
  }
  x0<-max(x)*0.66
  y0<-max(y)*0.66
  if(withTest){
    test<-cor.test(x,y, method = 'spearman')
    # text(x0, y0, labels=sprintf('P-value=%.3f',test$p.value), pos=3)
    text(x0, y0, labels=sprintf('R=%.2f',test$estimate))
  }
}

getHeatmapAnnotatio<-function(sample){
  column_annotation <-HeatmapAnnotation(
    df=data.frame(Stage=sample$Group),
    col = list(Stage =colorMapStage),
    show_annotation_name =FALSE,
    annotation_name_side='left'
  )
  column_annotation
}
############################################################################################
chromFactorLevel<-c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14',
     'chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY')
groupFactorLevel<-c('CTL', 'AIS', 'MIA', 'IAC')
gonomicRegionFactorLevel<-c('cgIslands', 'cgShores', 'cgShelves', 'cgSea',
                            'tss','promoter.1k', 'promoter.5k', 'utr5', 'utr3','exons','intron','intergenic')
# stageColor<-c("#31a354", "#addd8e", "#fd8d3c", "#e31a1c")
colorMapDAR <- c("#1f77b4", "#ff7f0e", "#17becf", "#e377c2", "#2ca02c", "#9467bd")
names(colorMapDAR)<-c('AISHypoDARs','AISHyperDARs','MIAHypoDARs','MIAHyperDARs','IACHypoDARs','IACHyperDARs')
colorMapGroup<-c("#1f77b4", "#ff7f0e", "#17becf", "#e377c2", "#2ca02c", "#9467bd")
names(colorMapGroup)<-c('HypoInAIS','HyperInAIS','HypoInIMA','HyperInMIA','HypoInIAC','HyperInIAC')
colorMapStage<-c('#00FF00','#00BFFF','#FFB90F','#FF0000')
names(colorMapStage)<-groupFactorLevel
colorMapStage2<-c('green','cyan','orange','red')
names(colorMapStage2)<-groupFactorLevel


loadSRDMR<-function(){
  SRDMR<-loadData2(file.path(CONFIG$dataIntermediate,'wgbs', 'srdmr.bed'),header = FALSE)
  colnames(SRDMR)<-c('chrom','start', 'end', 'class','cpg','length')
  SRDMR
}
loadSRDAR<-function(){
  SRDAR<-loadData2(file.path(CONFIG$dataIntermediate,'atac', 'srdar.bed'),header = FALSE)
  colnames(SRDAR)<-c('chrom','start', 'end', 'class')
  SRDAR 
}

loadSRDEG<-function(){
  SRDEG<-readRDS(file.path(CONFIG$dataIntermediate,'rna', 'srdeg.rds'))
  SRDEG
}
#############################################################################################
Group <- setRefClass(
  "group",
  fields = list(table = "data.frame", list='list',color.map='data.frame'),
  methods = list(
    initialize = function(data) {
      data$Group<-factor(data$Group, levels = groupFactorLevel)
      table<<-data
      list <<-split(table, table$Group)
      color.map<<-data.frame(
        group=names(colorMapStage),
        colors=colorMapStage
      )
    },
    select = function(groups,colors=NULL) {
      out<-data.frame()
      for(i in 1:length(groups)){
        tmp<-list[[groups[i]]]
        if(!is.null(colors)){
        }else{
          colors<-color.map$colors[match(groups,color.map$group)]
        }
        tmp$colors<-colors[i]
        out <- rbind(out,tmp)
      }
      out$Group<-factor(out$Group, levels = groupFactorLevel)
      out
    },
    selectBySample=function(samples) {
      tmp<-table[match(samples, table$SampleName),]
      tmp<-left_join(tmp, color.map, by=c('Group'='group'))
      tmp$Group<-factor(tmp$Group, levels = groupFactorLevel)
      tmp
    },
    pickColumnsByGroup=function(groups, data, na.rm=FALSE) {
      tmp<-select(groups)
      keep<-match(tmp$SampleName, colnames(data))
      if (na.rm){
        keep=na.omit(keep)
      }
      data[,keep]
    },
    getColorMapVec=function(){
      colorMapStage
    },
    show = function() {
      print(table(table$Group))
    }
  )
)

getGroups <- function() {
  data<-loadData2(file.path(CONFIG$dataExternal, 'samples.csv'))
  ret <- list()
  ret$WGBS <- Group$new(filter(data,WGBS=='Yes'))
  ret$RNA <-  Group$new(filter(data,RNA.seq=='Yes'))
  ret$ATAC <-  Group$new(filter(data,ATAC.seq=='Yes'))
  ret$WGBS.RNA <-  Group$new(filter(data,WGBS=='Yes', RNA.seq=='Yes'))
  ret$WGBS.ATAC <-  Group$new(filter(data,WGBS=='Yes', ATAC.seq=='Yes'))
  ret$WGBS.RNA.ATAC <- Group$new(filter(data,WGBS=='Yes', RNA.seq=='Yes', ATAC.seq=='Yes'))
  ret
}

saveImage <- function(file,...){
  file.path=file.path(IMAGE_DIR, file)
  if (endsWith(file, '.pdf')){
    pdf(file=file.path, ...)
  }
}

saveImage2 <- function(file,...){
  file.path=file.path(CONFIG$dataResult, file)
  if (endsWith(file, '.pdf')){
    pdf(file=file.path, ...)
  }
}

saveTsv<-function(data, filename,col.names = TRUE) {
  write.table(data, filename, quote = FALSE, sep = "\t", row.names = FALSE,col.names=col.names)
}
saveCsv<-function(data, filename) {
  write.csv(data, filename, quote = FALSE, row.names = FALSE)
}
################################################################################################
groups <- getGroups()
################################################################################################
RnaTPM <- setRefClass(
  "RnaTPM",
  fields = list(data = "data.frame", symbol='vector', ensemble='vector',sample='data.frame'),
  methods = list(
    initialize = function(batch='WGBS.RNA.ATAC') {
      rnaTPM<-readRDS(file.path(CONFIG$dataIntermediate,'rna', 'rnaTPM.wo.mt.rds'))
      symbol<<-ensemble2Symbol3(rnaTPM$ensemble)
      ensemble<<-rnaTPM$ensemble
      data<<-groups[[batch]]$pickColumnsByGroup(names(colorMapStage), rnaTPM)
      sample<<-groups[[batch]]$selectBySample(colnames(data))
    },
    getTPM = function(geneSymbol, reshapeIfOneRow=TRUE, rmIfcontainNA=FALSE) {
      out<-data[match(geneSymbol,symbol),]
      if(nrow(out)==1){
        if (reshapeIfOneRow){
          unlist(out)
        }else{
          rownames(out)<-geneSymbol
          out
        }
      }else{
        if (any(duplicated(geneSymbol))){
          rownames(out)<-NULL
        }else{
          rownames(out)<-geneSymbol
        }
        if(rmIfcontainNA){
          keep<-apply(out, 1, function(x){
            !(sum(is.na(x))>0)
          })
          out<-out[keep,]
        }
        out
      }
    },
    plotStageBar=function(geneSymbol){
      tpm<-getTPM(geneSymbol)
      dataInput<-data.frame(
        stage=sample$Group,
        value=tpm
      )
      plotBarError(dataInput, xlabel=geneSymbol,ylabel="TPM")
    },
    plotStageHeatmap=function(geneSymbols, colAnno=TRUE, colNames=TRUE){
      mm<-rnaTPM$getTPM(geneSymbols,reshapeIfOneRow=FALSE)
      column_annotation=NULL
      if(colAnno){
        column_annotation <-HeatmapAnnotation(
          df=data.frame(Stage=sample$Group),
          col = list(Stage =colorMapStage),
          show_annotation_name =FALSE,
          annotation_name_side='left'
        )
      }
      mm<-mm[!is.na(mm[,1]),]
      Heatmap(t(scale(t(mm))),
              cluster_rows=TRUE,
              cluster_columns = FALSE,
              show_row_names=TRUE,
              show_column_names=colNames,
              col=colorRamp2(c(-2, 0, 2), c("green3", "black", "red3")),
              top_annotation = column_annotation
      )
    },
    show = function() {
      print(data)
    }
  )
)


AtacPeak <- setRefClass(
  "AtacPeak",
  fields = list(dataMethy = "data.frame",dataAllMethy = "data.frame",dataAccess='data.frame', feature='vector',sample='data.frame',batch='character'),
  methods = list(
    initialize = function(batchSelect='WGBS.RNA.ATAC') {
      batch<<-batchSelect
      atacPeakMethyLevel<-readRDS(file.path(CONFIG$dataIntermediate,'atac', 'atacPeakMethyLevel.rds'))
      atacPeakAllMethyLevel<-readRDS(file.path(CONFIG$dataIntermediate,'atac', 'atacPeakAllMethyLevel.rds'))
      atacPeakTPM<-readRDS(file.path(CONFIG$dataIntermediate,'atac', 'atacPeakTPM.rds'))
      
      dataAllMethy<<-groups[[batch]]$pickColumnsByGroup(names(colorMapStage), atacPeakAllMethyLevel,na.rm=TRUE)
      dataMethy<<-groups[[batch]]$pickColumnsByGroup(names(colorMapStage), atacPeakMethyLevel,na.rm=TRUE)
      dataAccess<<-groups[[batch]]$pickColumnsByGroup(names(colorMapStage), atacPeakTPM, na.rm = TRUE)
      feature<<-bed2Feature(atacPeakTPM)
      sample<<-groups[[batch]]$selectBySample(colnames(dataMethy))
    },
    getMatchData = function(peakFeature,group=NULL,doScale=FALSE) {
      access<-getAccess(peakFeature,group,doScale)
      methy<-getMethy(peakFeature,group,doScale)
      if (length(peakFeature)>=2){
        keep<-apply(methy, 1, function(x){
          sum(x,na.rm = TRUE)!=0
        })
        access<-access[keep,]
        methy<-methy[keep,]
      }
      list(
        access=access,
        methy=methy
      )
    },
    getAccess = function(peakFeature, group=NULL,doScale=FALSE) {
      getData(peakFeature, dataAccess,group=group,doScale=doScale)
    },
    getMethy = function(peakFeature, group=NULL,doScale=FALSE) {
      getData(peakFeature, dataMethy,group=group,doScale=doScale)
    },
    getAllMethy = function(peakFeature, group=NULL,doScale=FALSE) {
      getData(peakFeature, dataAllMethy,group=group,doScale=doScale)
    },
    getData = function(peakFeature,data, group=NULL,doScale=FALSE) {
      out<-data[match(peakFeature,feature),]
      if (doScale){
        out<-t(scale(t(out)))
      }
      if(!is.null(group)) {
        out<-out[,sample$Group%in%group]
      }
      if(nrow(out)==1){
        unlist(out)
      }else{
        if (any(duplicated(peakFeature))){
          rownames(out)<-NULL
        }else{
          rownames(out)<-peakFeature
        }
        out
      }
    },
    plotCorAccessVsExpression=function(feature, gene, withTest=FALSE, featureTitle=FALSE){
      epi<-getMatchData(c(feature))
      rnaTPM<-RnaTPM()
      rna<-rnaTPM$getTPM(gene)
      if(featureTitle){
        plotDot(epi$access, rna, sample$colors, 'Accessibility', 'TPM',feature,withTest=withTest)
      }else{
        plotDot(epi$access, rna, sample$colors, 'Accessibility', 'TPM',gene,withTest=withTest)
      }
    },
    plotCorMethyVsExpression=function(feature, gene, withTest=FALSE, featureTitle=FALSE){
      epi<-getMatchData(c(feature))
      rnaTPM<-RnaTPM()
      rna<-rnaTPM$getTPM(gene)
      if(featureTitle){
        plotDot(epi$methy, rna, sample$colors, 'Methylation', 'TPM',feature,withTest=withTest)
      }else{
        plotDot(epi$methy, rna, sample$colors, 'Methylation', 'TPM',gene,withTest=withTest)
      }
    },
    show = function() {
      print('AtacPeak')
    }
  )
)
MethyLevel <- setRefClass(
  "MethyLevel",
  fields = list(data = "data.frame",feature='vector',sample='data.frame'),
  methods = list(
    initialize = function(batch='WGBS') {
      ratio<-readRDS(file.path(CONFIG$dataIntermediate,'wgbs', 'srdmr.ratio.rds'))
      feature<<-bed2Feature(ratio)
      data<<-groups[[batch]]$pickColumnsByGroup(names(colorMapStage), ratio)
      sample<<-groups[[batch]]$selectBySample(colnames(data))
      rownames(data)<<-feature
    },
    getMethy=function(selectFeature, rmIfcontainNA=FALSE){
      out<-data[match(selectFeature, feature),]
      if(rmIfcontainNA){
        keep<-apply(out, 1, function(x){
          !(sum(is.na(x))>0)
        })
      }else{
        keep<-apply(out, 1, function(x){
          sum(x,na.rm = TRUE)!=0
        })
      }
      
      out[keep,]
    },
    show = function() {
      print(dim(data))
    }
  )
)

Srdeg <- setRefClass(
  "Srdeg",
  fields = list(deg = "list",dfs='data.frame'),
  methods = list(
    initialize = function() {
      deg<<-loadSRDEG()
    },
    all = function(){
      unique(unlist(deg$deg))
    }
  )
)

Survival <- setRefClass(
  "Survival",
  fields = list(os = "data.frame",dfs='data.frame'),
  methods = list(
    initialize = function() {
      os<<-read.csv(file.path(CONFIG$dataIntermediate,'gepia2', 'table_survival.txt'),sep = '\t')
      dfs<<-read.csv(file.path(CONFIG$dataIntermediate,'gepia2', 'table_df_survival.txt'),sep = '\t')
    },
    pickSurvival = function(genes) {
      list(
        os=intersect(genes, os$Gene.Symbol),
        dfs=intersect(genes,dfs$Gene.Symbol)
      )
    },
    pickSurvivalDeg = function(genes) {
      srdeg<-loadSRDEG()
      a<-lapply(srdeg$deg, function(x){
        intersect(intersect(x, os$Gene.Symbol),genes)
      })
      b<-lapply(srdeg$deg, function(x){
        intersect(intersect(x, dfs$Gene.Symbol),genes)
      })
      ret<-list()
      list(
        os=unlist(a)%>%unique(),
        dfs=unlist(b)%>%unique()
      )
    }
  )
)
