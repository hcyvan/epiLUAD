library(dplyr)

DATA_DIR='./tmp_data'
IMAGE_DIR='./tmp_image'

########################################### helper function ###########################################
percent2Numeric <- function(x){
  as.numeric(substr(x,0,nchar(x)-1))/100
}

printf <- function(format, ...){
  cat(sprintf(format, ...))
}

loadData <- function(name, ext='csv'){
  ext.path <- paste(file.path(IMAGE_DIR, name),ext,sep = '.')
  rds.path <- paste(file.path(IMAGE_DIR, name),'rds',sep = '.')
  rds.external.path <- paste(file.path(IMAGE_DIR, 'external', name),'rds',sep = '.')
  if(file.exists(rds.path)){
    readRDS(rds.path)
  }else if(file.exists(ext.path)) {
    if(ext=='csv'){
      data<-read.csv(ext.path)
      saveRDS(data, rds.path)
      data
    }else{
      stop(sprintf("The file format <%s> is not supported", ext)); 
    }
  }else if(file.exists(rds.external.path)){
    readRDS(rds.external.path)
  }else{
    stop(sprintf("File not exist: %s and %s", rds.path, ext.path))
  }
}

#############################################################################################


Group <- setRefClass(
  "group",
  fields = list(table = "data.frame", list='list'),
  methods = list(
    initialize = function(data) {
      table<<-data
      list <<-split(table, table$Group)
    },
    select = function(groups,colors=NULL) {
      color_map<-data.frame(
        group=c('L0','L1','L2','L3'),
        colors=c('green3','cyan','orange','red')
      )
      out<-data.frame()
      for(i in 1:length(groups)){
        tmp<-list[[groups[i]]]
        if(!is.null(colors)){
        }else{
          colors<-color_map$colors[match(groups,color_map$group)]
        }
        tmp$colors<-colors[i]
        out <- rbind(out,tmp)
      }
      out
    },
    show = function() {
      print(table(table$Group))
    }
  )
)

getGroups <- function() {
  data<-loadData('sampleInfo')
  ret <- list()
  ret$WGBS <- Group$new(filter(data,WGBS=='Yes'))
  ret$RNA <-  Group$new(filter(data,RNA_seq=='Yes'))
  ret$ATAC <-  Group$new(filter(data,ATAC_seq=='Yes'))
  ret$WGBS.RNA <-  Group$new(filter(data,WGBS=='Yes', RNA_seq=='Yes'))
  ret$WGBS.ATAC <-  Group$new(filter(data,WGBS=='Yes', ATAC_seq=='Yes'))
  ret$WGBS.RNA.ATAC <- Group$new(filter(data,WGBS=='Yes', RNA_seq=='Yes', ATAC_seq=='Yes'))
  ret
}

saveImage <- function(file,...){
  file.path=file.path(IMAGE_DIR, file)
  if (endsWith(file, '.pdf')){
    pdf(file=file.path, ...)
  }
}


################################################################################################
groups <- getGroups()





