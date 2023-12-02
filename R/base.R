library(dplyr)
library(yaml)
library(tools)


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
homerKnownTFs<-function(dir){
  result.txt<-file.path(dir,'knownResults.txt')
  result<-read.csv(result.txt,sep='\t')
  result<-result[result$q.value..Benjamini.< 0.05,]
  result
  sapply(result$Motif.Name, function(x){
    c(strsplit(x,"\\(")[[1]][1])
  })
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
loadData2<-function(filename, file.format=NULL, force.refresh=FALSE, header=TRUE) {
  filename.rds<-paste0(file_path_sans_ext(filename),'.rds')
  if (is.null(file.format)){
    file.format=file_ext(filename)
  }
  
  if (file.exists(filename.rds) && !force.refresh) {
    readRDS(filename.rds)
  } else {
    if(file.format=='bed'){
      data<-read.csv(filename,sep = '\t',check.names = FALSE, header=header)
      saveRDS(data, filename.rds)
      data
    } else if(file.format=='csv'){
      data<-read.csv(filename)
      saveRDS(data, filename.rds)
      data
    }
  }
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
colorMapSRDMR<-c("#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd")
names(colorMapSRDMR)<-c('Early-Hyper-DMR','Early-Hypo-DMR','Late-Hyper-DMR','Late-Hypo-DMR')
colorMapStage<-c('#00FF00','#00BFFF','#FFB90F','#FF0000')
names(colorMapStage)<-groupFactorLevel

loadSRDMR<-function(){
  SRDMR<-loadData2(file.path(CONFIG$dataIntermediate, 'srdmr.s2.bed'),header = FALSE)
  colnames(SRDMR)<-c('chrom','start', 'end', 'class','cpg','length')
  SRDMR$class<-sub('DMC','DMR',SRDMR$class)
  SRDMR
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

################################################################################################
groups <- getGroups()

