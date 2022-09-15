library(ggplot2)
library(ggpubr)
library(ape)
library(FactoMineR)
library(factoextra)



source('./R/base.R')


#----------------------------------------------------------------------------------------------------------------------
# The mean bisulfite conversion and  CpG depth
#----------------------------------------------------------------------------------------------------------------------
wgbsInfo <- loadData('wgbsInfo')
conversion <- percent2Numeric(wgbsInfo$MCALL_bisulfite_conversion.ratio)
cpgDepth <- wgbsInfo$MCALL_CG_depth
printf("The mean bisulfite conversion is: %s", mean(conversion))
printf("The mean CpG depth is: %s", mean(cpgDepth))

#----------------------------------------------------------------------------------------------------------------------
# Figure S1 a. The mean DNA methylation of all samples				
#----------------------------------------------------------------------------------------------------------------------
type<-groups$WGBS$select(c('L0','L1','L2','L3'))
wgbsInfo <- loadData('wgbsInfo')
wgbsInfo$MCALL_MeanRatioCG_3X<-percent2Numeric(wgbsInfo$MCALL_MeanRatioCG_3X)
methy<-wgbsInfo[match(type$SampleName ,wgbsInfo$SampleName),]$MCALL_MeanRatioCG_3X
saveImage("mean.methy.level.all.sample.pdf",width = 10,height = 4)
barplot(methy,col=type$colors,names=type$SampleName,las=2,cex.names =0.5,ylim=c(0,1.0),ylab="Methylation Level")
legend("top", 
       legend = c('L0','L1','L2', 'L3'), 
       fill = unique(type$colors),bty = "n",ncol = 4)
dev.off()

#----------------------------------------------------------------------------------------------------------------------
# Figure 1 a. The mean DNA methylation of different groups	
#----------------------------------------------------------------------------------------------------------------------
wgbsInfo <- loadData('wgbsInfo')
wgbsInfo$MCALL_MeanRatioCG_3X<-percent2Numeric(wgbsInfo$MCALL_MeanRatioCG_3X)
methy<-wgbsInfo[match(type$SampleName ,wgbsInfo$SampleName),]$MCALL_MeanRatioCG_3X

data<-data.frame(ratio=methy, group=type$Group, color=type$Group)
col<-unique(type$colors)
names(col)<-unique(type$type)
saveImage("mean.methy.level.all.sample.violin.plot.pdf",width = 6,height = 3.5)

ggplot(data, aes(y=ratio, x=group,fill=color)) +
  scale_fill_manual(values=col)+
  geom_violin(trim=FALSE) +
  geom_jitter(shape=17, position=position_jitter(0.2), colour='black', size=0.5)+
  stat_summary(fun.data=mean_sdl, mult=1, 
               geom="pointrange", color="black")+
  theme_classic()+
  stat_compare_means( comparisons = list(c('L0','L1'),c('L0','L3')),label = 'p.signif', method = "wilcox.test")+
  xlab("")+
  ylab("Mean Methylation Level")+
  theme(
    legend.position="right",
    axis.title.x = element_text(size=0),
    axis.title.y = element_text(size=12),
    axis.text = element_text(size = 12,colour="black"),
    legend.title = element_blank(),
    legend.text = element_text(size=12)
  )+
  guides(colour = guide_legend(override.aes = list(shape = 12,size=10)))
dev.off()

#----------------------------------------------------------------------------------------------------------------------
# Figure S1 b. Unsupervised hierarchical clustering
#----------------------------------------------------------------------------------------------------------------------

bed <- loadData('sampleMethyLevelDepth10x')
data <- bed[,-1:-3]
type<-groups$WGBS$select(c('L0','L1','L2','L3'))
type <- type[match(colnames(data),type$SampleName),]
dist_mm<-dist(t(data))
hclust_avg <- hclust(dist_mm)
saveImage("methylationLevel.cluster.pdf",width = 10,height = 10)
par(mar = c(0,0,0,0))
phyl<-as.phylo(hclust_avg)
plot(phyl, type = "fan",tip.color=type$colors,label.offset=3.5)
tiplabels(pch=21, col="black", adj=0, bg=type$colors, cex=2.5)
dev.off()

saveImage("methylationLevel.cluster.legend.pdf",width = 5,height = 3)
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =distinct(type,Group,.keep_all = TRUE)$Group, pch=15, pt.cex=3, cex=1.5, bty='n',col = distinct(type,Group,.keep_all = TRUE)$colors,horiz=TRUE)
dev.off()


#----------------------------------------------------------------------------------------------------------------------
# Figure S1 c. principal component analysis
#----------------------------------------------------------------------------------------------------------------------
bed <- loadData('sampleMethyLevelDepth10x')
data <- bed[,-1:-3]
type<-groups$WGBS$select(c('L0','L1','L2','L3'))
type <- type[match(colnames(data),type$SampleName),]
res.pca <- PCA(t(((as.matrix(data)))), graph = FALSE)
saveImage("methylationLevel.pca.pdf",width = 10,height = 10)
fviz_pca_ind(res.pca,
             repel = TRUE,
             label="none",
             col.ind=type$Group,
             palette=distinct(type,Group,.keep_all = TRUE)$colors,
             addEllipses=TRUE,
             # ellipse.type='t',
             # ellipse.level=0.75,
             ellipse.alpha = 0.01,
             pointsize =1.5,
             col.var="black",
             axes=c(1,2),
             pointshape = 19)+
  theme_classic()+
  theme(plot.title = element_blank(),axis.text=element_text(size=20,colour="black"),axis.title=element_text(size=20),legend.position = "none")
dev.off()




