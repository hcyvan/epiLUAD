library("dplyr")
library("stringr")
library("ggplot2")
source('./base.R')

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
    chose <- 'new' #
    thre <- 5000 #0.25
} else {
    chose <- args[1]
    thre <- args[2]
}
num <- 9
figdir <- "/PUBLIC/gomics/renqiaoling/analyze/fig/"
dmcdir <- paste0('/PUBLIC/gomics/renqiaoling/2023lungcancer/wgbs/4.mcomp',num,'/result/dmc.LNM.vs.')

select_cpg <- function(methy,chose){
    #groups=c('LSQ','LCC','SCLC')
    groups=c('AIS','MIA','IAC')
    if (chose %in% groups) {
        filname=paste(dmcdir,chose,".txt",sep='')
        dmc <- read.csv(filname,sep="\t")
        upset_list <- paste(dmc[,1],dmc[,2],dmc[,3])
    } else {
        filname=paste(dmcdir,groups[1],".txt",sep='')
        dmc <- read.csv(filname,sep="\t")
        if (thre>1) {
            dmc <- dmc[order(abs(dmc$credibleDif_1-0),decreasing=TRUE),][1:thre,]
        } else {
            dmc <- dmc[abs(dmc$credibleDif_1-0)>=thre,]
        }
        print(summary(dmc$credibleDif_1-0))
        upset_list <- paste(dmc[,1],dmc[,2],dmc[,3])
        for (samp in groups[2:length(groups)]) {
          filname=paste(dmcdir,samp,".txt",sep='')
          dmc <- read.csv(filname,sep="\t")
          if (thre>1) {
              dmc <- dmc[order(abs(dmc$credibleDif_1-0),decreasing=TRUE),][1:thre,]
          } else {
              dmc <- dmc[abs(dmc$credibleDif_1-0)>=thre,]
          }
          print(summary(dmc$credibleDif_1-0))
          if (chose=='intersect') {
            upset_list <- intersect(upset_list, paste(dmc[,1],dmc[,2],dmc[,3]))
          } else {
            upset_list <- union(upset_list, paste(dmc[,1],dmc[,2],dmc[,3]))
          }
          rm(dmc);gc()
        }
    }
    print(length(upset_list))
    cpg <- paste(methy[,1],methy[,2],methy[,3])
    return(methy[cpg %in% upset_list,]) #-1:-3
}

na_stat <- function(d){
    stat <- apply(d,1,function(x){sum(is.na(x))})
    print(c(sum(is.na(d)),sum(stat==0),sum(stat>2),sum(stat>100)))
}

rmna <- function(dat){
    dat_clean <- na.omit(dat)
    print(paste('raw rows:',nrow(dat),'; clean rows:',nrow(dat_clean)))
    return(dat_clean)
}

dheatmap <- function(dat,fign){
    if(chose=='union') {return()}
    library("pheatmap")
    dat <- rmna(dat)
    anno_col <- data.frame(group=salia, row.names=colnames(dat))
    a <- color_map$colors ; names(a)<-color_map$groups; anno_color <- list(group=a)
    
    ph <- pheatmap(dat, annotation_col=anno_col,annotation_names_col=FALSE, 
        show_rownames=FALSE, annotation_colors = anno_color,
        fontsize_col=5,width=16,height=9,
        filename = paste0(figdir,fign,".png"))
}

dtsne <- function(dat,fign,perp,miter=2000) {
    library("Rtsne")
    tsne_out = Rtsne(t(dat),perplexity = perp, verbose = TRUE, max_iter=miter)
    pdat = data.frame(tsne_out$Y)
    colnames(pdat) = c("D1","D2")
    
    ip2 <- ggplot(data=pdat,aes(x=D1,y=D2,color=salia))+
        geom_point(size = 3.0, alpha = 0.4 ,color=pcolor)+
        ggrepel::geom_text_repel(aes(label=colnames(dat)),size=4.0,segment.color=pcolor,color=pcolor,max.overlaps=10)+
        geom_hline(yintercept = 0,lty=2,col="black",lwd=1)+
        geom_vline(xintercept = 0,lty=2,col="black",lwd=1)+
        labs(title=sprintf("depth 10x perplexity=%2d iter=%4d",perp,miter))
        theme_bw()
    ggsave(paste0(figdir,fign,perp,".png"),width=20,height=10,dpi=300)
}

dviolin <- function(dat,fign){
    library("ggpubr")
    library("Hmisc")
    #dat <- rmna(dat)
    mdat <- data.frame(ratio=colMeans(dat,na.rm=TRUE),gtype=salia)
    print(head(mdat[order(mdat$ratio),],10))
    print(tail(mdat[order(mdat$ratio),],30))
    mdat$gtype <- factor(mdat$gtype,levels=color_map$groups)
    print(str(mdat))
    a <- color_map$colors ; names(a) <- color_map$groups
    print(summary(a))
    p2<-ggplot(data=mdat,aes(x=gtype,y=ratio,fill=gtype))+
        scale_fill_manual(values=a) +
        geom_violin(trim=FALSE) +
        geom_jitter(shape=17, position=position_jitter(0.2), colour='black', size=0.5)+
        stat_summary(fun.data=mean_sdl, geom="pointrange", color="black")+
        theme_classic()+
        stat_compare_means( comparisons = list(c('LNM','AIS'),c('LNM','MIA'),c('LNM','IAC')),
                           label = 'p.signif', method = "t.test")+
        labs(x='',y='Mean Methylation Level')+
        theme(legend.position="right",
            axis.title.x = element_text(size=0),
            axis.title.y = element_text(size=12),
            axis.text = element_text(size = 12,colour="black"),
            legend.title = element_blank(),
            legend.text = element_text(size=12))+
        guides(colour = guide_legend(override.aes = list(shape = 12,size=10)))
    ggsave(paste0(figdir,fign,".pdf"),width=6,height=3.5)
    #ggsave(paste0(figdir,fign,".png"),width=6,height=3.5,dpi=300)
}

dpca <- function(dat,fign){
    library("ggrepel")
    dat <- rmna(dat)
    pca1 <- prcomp(t(dat),scale.=FALSE,rank.=3)
    print(summary(pca1))
    df1 <- pca1$x
    # ps1 <- predict(pca1) |> as.data.frame()
    summ<-summary(pca1)
    xlab<-paste("pc1(",round(summ$importance[2,1]*100,2),"%)",sep="")
    ylab<-paste("pc2(",round(summ$importance[2,2]*100,2),"%)",sep="")
    p2<-ggplot(data=as.data.frame(df1),aes(x=PC1,y=PC2,color=salia))+
        geom_point(size = 3.0, alpha = 0.4 ,color=pcolor)+
        ggrepel::geom_text_repel(aes(label=colnames(dat)),size=4.0,segment.color=pcolor,color=pcolor,max.overlaps=10)+
        #scale_color_manual(values=color_map$colors,labels=color_map$groups)+
        labs(x=xlab,y=ylab,color="",cex=2)+
        geom_hline(yintercept = 0,lty=2,col="black",lwd=1)+
        geom_vline(xintercept = 0,lty=2,col="black",lwd=1)+
        theme_bw()
    ggsave(paste0(figdir,fign,".png"),width=20,height=10,dpi=300)
}

dclust <- function(dat,fign,clmd){
#dclust <- function(dist_mm,fign,clmd='ward.D'){
    library("ape")
    dist_mm<-dist(t(dat))
    print(paste('na number in dist_mm=',sum(is.na(dist_mm))))
    hclust_avg <- hclust(dist_mm,method = clmd)
    print(paste('na number in hclust_avg=',sum(is.na(hclust_avg$height))))
    
    figname <- paste0(figdir,fign,clmd,".png")
    png(figname,width=20,height=5,unit='in',res=300)
    
    par(mar = c(0,0,0,0))
    phyl<-as.phylo(hclust_avg)
    plot(phyl, direction="downward",tip.color=pcolor, no.margin=TRUE,
         label.offset=1, cex=0.4,font=2,plot=TRUE,main=clmd) 
    tiplabels(pch=21, col=pcolor, bg=pcolor, cex=0.8)
    #tiplabels(pch=21, col=hcolor, bg=hcolor, cex=0.8,adj=c(0.5,8.0))
    #tiplabels(pch=21, col=agecol, bg=agecol, cex=0.8,adj=c(0.5,4.0))
    #tiplabels(pch=21, col=sexcol, bg=sexcol, cex=0.8,adj=c(0.5,2.0))
    legend("bottom", legend =color_map$groups,col = color_map$colors,
           pch=15, pt.cex=1, cex=1, bty='n',horiz=TRUE)
    
    dev.off()
    system(paste('convert',figname,'-trim',figname),intern = FALSE)
}

#groups=c('LNM' ,'LAD' ,'LSQ'   ,'SCLC','TRANS','AAH'),
#deltype=c("#N/A","",'LSQ','SCLC','TRANS','LAD','AAH') #
#color_map <- data.frame(
#  group1=c(3,4,5,6,7,8,1),
#  group2=c('female','male','','','','',''),
#  group3=c('N'     ,'H'   ,'M'      ,'1'   ,'1'    ,'L'      ,'1'),
#  groups=c('LNM'   ,'AIS' ,'MIA'    ,'IAC' ,'AAH'  ,'LPA'    ,'LAD'  ),
#  colors=c('green3','blue','#9B59B6','red' ,'black','#FFA500','black'))
color_map <- data.frame(
  #groups=c('LNM'   ,'LSQ'    ,'LCC'    ,'SCLC'),
  groups=c('LNM'   ,'AIS'   ,'MIA'    ,'IAC'),
  colors=c('green3','#00BFFF','#FFB90F','red'))
  #groups=c('LNM'    ,'AIS'   ,'MIA'    ,'IAC' ,'LCC'    ,'LSQ'  ,'SCLC' ),
  #colors=c('green3','#00BFFF','#FFB90F','red' ,'salmon2','magenta1','cyan3' ))
#selctype=color_map$groups[-1]
selctype=c('FIB','HAM','INF','TB',color_map$groups)
print(selctype)

read_rna <- function(){
    group <- read.csv("/PUBLIC/gomics/renqiaoling/2023lungcancer/sampleinfo_rql.csv")
    group$type_rql2[group$type_rql2 %in% c('FIB','HAM','INF','TB')] <- 'LNM'
    #group$note_rql[group$ID %in% c(261,184,127,98)] <- 1
    group <- filter(group, type_rql2 %in% groups, note_rql==0)[,c('ID','RNA','type_rql2')]
    
    rnacount <- read.csv("/PUBLIC/gomics/chengyihang/yixing/rna_seq/merge_feature_count/4.feature_count.fpkm",check.names = FALSE)
    colnames(rnacount) <- str_split_fixed(colnames(rnacount),'_',n=2)[,1]
    print(dim(rnacount))

    dat <- rnacount[,colnames(rnacount) %in% group$RNA]
    idx <- match(colnames(dat),group$RNA)
    colnames(dat) <- paste(group$type_rql2[idx],group$ID[idx],sep="_")
    rm(rnacount,group);gc()
    return(dat)
}

read_methy <- function(){
    group <- read.csv("/PUBLIC/gomics/renqiaoling/2023lungcancer/sampleinfo_rql.csv",fileEncoding = "GB2312")
    #info <- read.csv("/PUBLIC/gomics/renqiaoling/2023lungcancer/clinical_information.csv",sep=",")
    #group$age <- info$age[match(group$ID,info$ID)]
    #group$sex <- info$sex[match(group$ID,info$ID)]

    methy <- readRDS('/PUBLIC/gomics/renqiaoling/2023lungcancer/wgbs/3.moabs.sort.merge.bed/merge_ratio_bed_d3_all.rds')
    #methy <- readRDS('/PUBLIC/gomics/renqiaoling/2023lungcancer/wgbs/3.moabs.sort.merge.bed/merge_ratio_bed_d10_one.rds')
    colnames(methy) <- gsub("\\.","-",colnames(methy))

    #dat  <- methy[,colnames(methy) %in% group$WGBS[
    #  ! str_split_fixed(group$type_rql2,"_",n=2)[,1] %in% deltype]] # also delete the first three columns
    methy <- methy[,colnames(methy) %in% c("X-chrom","start","end",group$WGBS[group$type_rql2 %in% selctype & group$note_rql==0])] # also delete the first three columns
    print(dim(methy)) # how many dmc are selected to draw
    if ('LSQ' %in% selctype | 'IAC' %in% selctype){
    methy2 <- readRDS('/PUBLIC/gomics/renqiaoling/2023lungcancer/wgbs/3.moabs.sort.merge.bed/merge_ratio_bed_d3_all_shanghai.rds')
    colnames(methy2) <- gsub("\\.","-",colnames(methy2))
    methy2 <- methy2[,colnames(methy2) %in% c("X-chrom","start","end",group$WGBS[group$type_rql2 %in% selctype & group$note_rql==0])] # also delete the first three columns
    methy <- full_join(methy,methy2,by=c("X-chrom","start","end" ))
    print(dim(methy)) # how many dmc are selected to draw
    }
    if(chose=='new'){
        dat <- methy#[,-1:-3]
    }else{
        dat <- select_cpg(methy,chose)
    }
    rm(methy);gc()

    idx <- match(colnames(dat)[-1:-3],group$WGBS)
    colnames(dat)[-1:-3] <- paste(group$type_rql2[idx],group$ID[idx],sep="_")
    #write.csv(dat,file="/PUBLIC/gomics/renqiaoling/analyze/mdata/AIS_union_dmc_methy_3x.csv",quote=FALSE,row.names = FALSE,col.names = TRUE)
    return(dat[,-1:-3])
}

dat <- read_methy()
salia <- str_split_fixed(colnames(dat),"_",n=2)[,1]
#salia[salia %in% c('AAH','AIS','MIA','LPA','IAC')] <- 'LAD'
salia[salia %in% c('FIB','HAM','INF','TB')] <- 'LNM'
print(table(salia))
pcolor <- color_map$colors[match(salia,color_map$groups)]
#hcolor <- color_map$colors[match(group$degree[idx],color_map$group3)]
#sexcol <- color_map$colors[match(group$sex[idx],color_map$group2)]
#agecol <- color_map$colors[match(group$age[idx]%/%10,color_map$group1)]

print(dim(dat))
dat <- replace(dat,dat==-1.000,NA)
dviolin(dat,paste0(chose,'violin'))
#dheatmap(dat,paste0(chose,'heatmapNA'))

#dat <- rmna(dat)
#dclust(dat,paste0(num,chose,thre,color_map$groups[2],'clust'),'ward.D')
#dist_mm<-dist(t(dat))
#for (clmd in c('average','mcquitty','median','centroid','ward.D')) { 
#    dclust(dist_mm,paste0(chose,'clust'),clmd)
#}
#dpca(dat,paste0(num,chose,thre,color_map$groups[2],'pca'))

#dtsne(dat,paste0(num,chose,thre,color_map$groups[2],'d10_tsne'),3,3000)
#for (perp in seq(2,52,10)) {
#    dtsne(dat,'d10_tsne_5k',perp)
#}
