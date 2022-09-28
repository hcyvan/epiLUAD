library(reshape2)
library(dplyr)
library(ggplot2)
library(ggalluvial)
library(ComplexHeatmap)
library(circlize)
library(ChIPseeker)
# library(patchwork)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


source('./R/base.R')



#----------------------------------------------------------------------------------------------------------------------
# Figure 2 a(left). Obtain SR-DMCs from DMCs of L1-L3 groups
#----------------------------------------------------------------------------------------------------------------------
extractDmc <- function(dmcTag,group) {
  dmc<-loadDataBed(dmcTag)
  dmc2<-dplyr::select(dmc,chrom=X.chrom,start,end, class)
  dmc2$class<-'Hypo'
  dmc2$class[dmc$class=='strongHyper'] <- 'Hyper'
  dmc2$group<-group
  dmc2
}
dmcL0vsL1<-extractDmc('dmcL0vsL1','L1')
dmcL0vsL2<-extractDmc('dmcL0vsL2','L2')
dmcL0vsL3<-extractDmc('dmcL0vsL3','L3')

lData<-Reduce(rbind,list(dmcL0vsL1,dmcL0vsL2,dmcL0vsL3))
srdmc<-dcast(lData,chrom + start+ end ~ group,value.var = 'class')
srdmc[,4:6][is.na(srdmc[,4:6])] <- "NC"
df<-srdmc%>%count(L1,L2,L3)
df<-(df%>%arrange(desc(n)))
df$L1 = factor(df$L1,levels = c('Hyper','NC','Hypo'))
df$L2 = factor(df$L2,levels = c('Hyper','NC','Hypo'))
df$L3 = factor(df$L3,levels = c('Hyper','NC','Hypo'))
df$fill<-as.factor(1:23)
df$col<-rep('gray', 23)
df$col[1]<-'green3'
df$col[2]<-'red'
df$col[3] <- 'darkorchid1'
df$col[8] <- 'darkorchid1'
df$col[10] <- 'darkorchid1'
df$col[5] <- 'cyan1'
df$col[13] <- 'cyan1'
df$col[14] <- 'cyan1'
df=df[1:14,]

saveImage("SRDMC.alluvium.pdf",width = 4,height = 9)
ggplot(df, aes(y = n, axis1 = L1, axis2 = L2, axis3=L3)) +
  geom_alluvium(aes(fill = fill)) +
  geom_stratum(width = 1/6, fill = "white") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),angle=90,size=5)+
  coord_cartesian(clip = 'off')+
  scale_fill_manual(values=df$col)+
  theme_bw()+
  theme(
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.x = element_text(size=20,vjust = 5,colour = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  )+
  guides(fill = F)+
  scale_x_continuous(breaks = 1:3, labels = c("L1", "L2", "L3"),position = 'top')
dev.off()

colors <- c('cyan1','darkorchid1','green3','red')
type <- c('Early-Hypo-DMC','Early-Hyper-DMC','Late-Hypo-DMC','Late-Hyper-DMC')
saveImage("SRDMC.legned.pdf",width = 5,height = 4)
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend = type, pch=15, pt.cex=3, cex=1.5, bty='n',col = alpha(colors, 0.6),horiz=FALSE)
dev.off()
#------------------------- Prepare SRDMC (stage related DMCs) -----------------------------------------------<<<
earlyHyperDmc<-Reduce(rbind,list(
  filter(srdmc,L1=='Hyper',L2=="Hyper", L3=="Hyper"),
  filter(srdmc,L1=='Hyper',L2=="Hyper", L3=="NC"),
  filter(srdmc,L1=='Hyper',L2=="NC", L3=="NC")
))%>%mutate(class='Early-Hyper-DMC')%>%select(chrom, start,end,class)
earlyHypoDmc<-Reduce(rbind,list(
  filter(srdmc,L1=='Hypo',L2=="Hypo", L3=="Hypo"),
  filter(srdmc,L1=='Hypo',L2=="Hypo", L3=="NC"),
  filter(srdmc,L1=='Hypo',L2=="NC", L3=="NC")
))%>%mutate(class='Early-Hypo-DMC')%>%select(chrom, start,end,class)
lateHyperDmc<-filter(srdmc,L1=='NC',L2=="NC", L3=="Hyper")%>%mutate(class='Late-Hyper-DMC')%>%select(chrom, start,end,class)
lateHypoDmc<-filter(srdmc,L1=='NC',L2=="NC", L3=="Hypo")%>%mutate(class='Late-Hypo-DMC')%>%select(chrom, start,end,class)
SRDMC<-Reduce(rbind,list(earlyHyperDmc,earlyHypoDmc,lateHyperDmc,lateHypoDmc))
SRDMC$chrom <- factor(SRDMC$chrom, levels=chromFactorLevel)
SRDMC<-arrange(SRDMC, chrom, start)
saveRDS(SRDMC, file.path(DATA_DIR, 'SRDMC.rds'))
#------------------------- Prepare SRDMC (stage related DMCs) ----------------------------------------------->>>

#----------------------------------------------------------------------------------------------------------------------
# Figure 2 a(right). The heatmap of SR-DMCs			
#----------------------------------------------------------------------------------------------------------------------
sampleMethyLevelDepth3xSRDMC<-loadData('sampleMethyLevelDepth3xSRDMC')
sampleMethyLevelDepth3xSRDMC<-removeNegativeOne(sampleMethyLevelDepth3xSRDMC)
SRDMC<-loadData('SRDMC')
type<-groups$WGBS$select(c('L0','L1','L2','L3'))


srdmcFrac<-sample_frac(SRDMC, 0.1)%>%arrange(chrom, start)
srdmcFrac<-inner_join(srdmcFrac, sampleMethyLevelDepth3xSRDMC,by=c('chrom', 'start'))%>%select(class, starts_with('L'))
srdmcFracList<-split(srdmcFrac,srdmcFrac$class)
srdmcFrac<-Reduce(rbind,list(
  srdmcFracList$`Early-Hyper-DMC`,
  srdmcFracList$`Early-Hypo-DMC`,
  srdmcFracList$`Late-Hyper-DMC`,
  srdmcFracList$`Late-Hypo-DMC`
))
sc<-table(srdmcFrac$class)
mm <- srdmcFrac[,-1]
row_type <- factor(
  Reduce(c,lapply(1:length(sc), function(x){
    rep(names(sc[x]), sc[x])
  })),
  levels = names(sc)
)
col_type<-factor(type[match(colnames(mm),type$SampleName),'Group'], levels = c('L3','L2','L1','L0'))
color_type<-distinct(type, Group,colors)
colColorMap<-color_type$colors
names(colColorMap)<-color_type$Group
top_ha<-HeatmapAnnotation(
  type=type$Group,
  col = list(type = colColorMap),
  annotation_legend_param = list(
    type=list(
      title = "Sample Type",
      title_position = "lefttop-rot"
    )
  ),
  show_legend = FALSE,
  show_annotation_name =FALSE
)
rowColorMap <- alpha(c('darkorchid1','cyan1','red','green3'), 0.6)
names(rowColorMap)<-names(sc)
left_ha<-rowAnnotation(
  stage=row_type,
  col = list(stage = rowColorMap),
  show_legend = FALSE,
  show_annotation_name =FALSE
)
col_fun = colorRamp2(c( -2,0,2),c("green", "black", "red"))
mm<-t(scale(t(mm)))
saveImage("SRDMC.heatmap.pdf",width = 5,height = 8)
ht<-Heatmap(mm,
            cluster_column_slices = FALSE,
            cluster_columns = T,
            cluster_rows  = F,
            column_names_gp = gpar(fontsize = 10),
            column_names_rot = 45,
            show_row_names = F,
            show_column_names = F,
            col = col_fun,
            column_split = type$Group,
            row_split = row_type,
            top_annotation = top_ha,
            left_annotation = left_ha,
            heatmap_legend_param = list(
              title = "Scaled DNA Methylation Levels",
              legend_height = unit(6, "cm"),
              at = c(-2,0,2),
              labels = c("-2",'0', "2"),
              title_position = "lefttop-rot"),
            
)
draw(ht,heatmap_legend_side='right')
dev.off()


#----------------------------------------------------------------------------------------------------------------------
# Figure S3 a. The distribution of SR-DMCs in the genome
#----------------------------------------------------------------------------------------------------------------------
SRDMC<-loadData('SRDMC')
SRDMC <- GRanges(SRDMC)
SRDMCList<-split(SRDMC, SRDMC$class)
options(ChIPseeker.ignore_1st_exon = T)
options(ChIPseeker.ignore_1st_intron = T)
options(ChIPseeker.ignore_downstream = T)
options(ChIPseeker.ignore_promoter_subcategory = T)
peakAnnoList<- lapply(SRDMCList, annotatePeak, tssRegion=c(-5000, 3000),TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene)
saveImage("SRDMC.distribution.Early-Hypo-DMC.pdf",width = 5,height = 3)
plotAnnoPie(peakAnnoList$`Early-Hypo-DMC`)
dev.off()
saveImage("SRDMC.distribution.Early-Hyper-DMC.pdf",width = 5,height = 3)
plotAnnoPie(peakAnnoList$`Early-Hyper-DMC`)
dev.off()
saveImage("SRDMC.distribution.Late-Hypo-DMC.pdf",width = 5,height = 3)
plotAnnoPie(peakAnnoList$`Late-Hypo-DMC`)
dev.off()
saveImage("SRDMC.distribution.Late-Hyper-DMC.pdf",width = 5,height = 3)
plotAnnoPie(peakAnnoList$`Late-Hyper-DMC`)
dev.off()

saveImage("SRDMC.distribution.distToTss.pdf",width = 7,height = 2)
plotDistToTSS(peakAnnoList,title ="", ylab = "Percentage of SR-DMC (%)")
dev.off()

#----------------------------------------------------------------------------------------------------------------------
# Figure S3 b. The distribution of SR-DMRs in the genome
#----------------------------------------------------------------------------------------------------------------------
SRDMR<-loadData('SRDMR')
SRDMR <- GRanges(SRDMR)
SRDMRList<-split(SRDMR, SRDMR$class)
options(ChIPseeker.ignore_1st_exon = T)
options(ChIPseeker.ignore_1st_intron = T)
options(ChIPseeker.ignore_downstream = T)
options(ChIPseeker.ignore_promoter_subcategory = T)
peakAnnoList<- lapply(SRDMRList, annotatePeak, tssRegion=c(-5000, 3000),TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene)
saveImage("SRDMR.distribution.Early-Hypo-DMR.pdf",width = 5,height = 3)
plotAnnoPie(peakAnnoList$`Early-Hypo-DMR`)
dev.off()
saveImage("SRDMR.distribution.Early-Hyper-DMR.pdf",width = 5,height = 3)
plotAnnoPie(peakAnnoList$`Early-Hyper-DMR`)
dev.off()
saveImage("SRDMR.distribution.Late-Hypo-DMR.pdf",width = 5,height = 3)
plotAnnoPie(peakAnnoList$`Late-Hypo-DMR`)
dev.off()
saveImage("SRDMR.distribution.Late-Hyper-DMR.pdf",width = 5,height = 3)
plotAnnoPie(peakAnnoList$`Late-Hyper-DMR`)
dev.off()

saveImage("SRDMR.distribution.distToTss.pdf",width = 7,height = 2)
plotDistToTSS(peakAnnoList,title ="", ylab = "Percentage of SR-DMR (%)")
dev.off()



##############################################################








###############################################################
colors <- c('cyan1','darkorchid1','green3','red')
type <- c('Early-Hypo-DMR','Early-Hyper-DMR','Late-Hypo-DMR','Late-Hyper-DMR')
saveImage("legend.SRDMR.pdf",width = 5,height = 4)
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend = type, pch=15, pt.cex=3, cex=1.5, bty='n',col = alpha(colors, 0.6),horiz=FALSE)
dev.off()
















