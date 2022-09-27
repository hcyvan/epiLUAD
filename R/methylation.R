library(ggplot2)
library(ggpubr)
library(ape)
library(FactoMineR)
library(factoextra)



source('./R/base.R')

drawDensity<-function(ratio1,ratio2,s1,s2,legend.position='topleft'){
  d1 <- density(ratio1)
  d2 <- density(ratio2)
  dens <- list(a=d1,b=d2)
  plot(NA, xlim=range(sapply(dens, "[", "x")),
       ylim=range(sapply(dens, "[", "y")),xlab="",ylab="Density",cex.lab=2, cex.axis=1.5)
  mapply(lines, dens, col=1:length(dens))
  polygon(d1, col=rgb(0, 1, 0,0.5), border=NA)
  polygon(d2, col=rgb(1, 0, 0,0.5), border=NA)
  legend(legend.position, legend=c(s1,s2), fill=c("green","red"), bty = "n")
}

drawDensityAll<-function(r0,r1,r2,r3,legend.position='topleft'){
  par(mar = c(3,5,1,1))
  d0<-density(r0)
  d1<-density(r1)
  d2<-density(r2)
  d3<-density(r3)
  dens <- list(a=d0,b=d1,c=d2,d=d3)
  plot(NA, xlim=range(sapply(dens, "[", "x")),
       ylim=range(sapply(dens, "[", "y")),xlab="",ylab="Density",cex.lab=2, cex.axis=1.5)
  polygon(d0, border='green3',lwd=1)
  polygon(d1, border='cyan',lwd=1)
  polygon(d2, border='orange',lwd=1)
  polygon(d3, border='red',lwd=1)
  legend(legend.position, legend=c('L0','L1','L2','L3'), fill=c('green3','cyan','orange','red'), bty = "n",horiz = T)
}


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
legend("topleft", legend =distinct(type,Group,.keep_all = TRUE)$Group,
       pch=15, pt.cex=3, cex=1.5, bty='n',col = distinct(type,Group,.keep_all = TRUE)$colors,horiz=TRUE)
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
  theme(plot.title = element_blank(),
        axis.text=element_text(size=20,colour="black"),axis.title=element_text(size=20),legend.position = "none")
dev.off()

#----------------------------------------------------------------------------------------------------------------------
# Figure 1 b. principal component analysis
#----------------------------------------------------------------------------------------------------------------------
groupMethyLevelDepth3xAllLine1<-loadData('groupMethyLevelDepth3xAllLine1')
line<-removeNegativeOne(groupMethyLevelDepth3xAllLine1)
saveImage("methyLevel.density.line1.pdf",width = 4,height = 3.5)
layout(1:3)
par(mar = c(3,5,1,1))
drawDensity(line$L0, line$L1, "L0", "L1")
drawDensity(line$L0, line$L2, "L0", "L2")
drawDensity(line$L0, line$L3, "L0", "L3")
dev.off()

#----------------------------------------------------------------------------------------------------------------------
# Figure S2 a. principal component analysis
#----------------------------------------------------------------------------------------------------------------------
groupMethyLevelDepth3xAllLine1<-loadData('groupMethyLevelDepth3xAllLine1')
line<-removeNegativeOne(groupMethyLevelDepth3xAllLine1)
saveImage("methyLevel.density.all.line1.pdf",width = 8,height = 3)
drawDensityAll(line$L0,line$L1,line$L2,line$L3,'topleft')
dev.off()


#----------------------------------------------------------------------------------------------------------------------
# Figure S2 a. principal component analysis
#----------------------------------------------------------------------------------------------------------------------
saveImage("methyLevel.density.cgi.pdf",width = 4,height = 3.5)
layout(1:3)
par(mar = c(3,5,1,1))
drawDensity(cgi$L0, cgi$L1, "L0", "L1",'topright')
drawDensity(cgi$L0, cgi$L2, "L0", "L2",'topright')
drawDensity(cgi$L0, cgi$L3, "L0", "L3",'topright')
dev.off()

#----------------------------------------------------------------------------------------------------------------------
# Figure S2 a. principal component analysis
#----------------------------------------------------------------------------------------------------------------------
groupMethyLevelDepth3xAllCgi<-loadData('groupMethyLevelDepth3xAllCgi')
cgi<-removeNegativeOne(groupMethyLevelDepth3xAllCgi)
saveImage("methyLevel.density.all.cgi.pdf",width = 8,height = 3)
drawDensityAll(cgi$L0,cgi$L1,cgi$L2,cgi$L3,'topright')
dev.off()

#----------------------------------------------------------------------------------------------------------------------
# Figure 1 e. The DMC density distribute
#----------------------------------------------------------------------------------------------------------------------
draw.scatter <- function(dmc,xlab,ylab){
  percent <- function(x, digits = 2, format = "f", ...) {
    paste0(formatC(x * 1000, format = format, digits = digits, ...), "‰")
  }
  data<-data.frame(x=dmc$nominalRatio_0, y=dmc$nominalRatio_1)
  palette <- colorRampPalette(c("blue", "yellow", "red"))
  smoothScatter(data,colramp = palette,xlab=xlab,ylab=ylab,cex.lab=2, cex.axis=1.5)
  total<-29401360
  hyper<-table(dmc$class)[1]
  hypo<-table(dmc$class)[2]
  text(x = 0.15, y = 0.85, labels = percent(hyper/total),cex=1.5)
  text(x = 0.85, y = 0.15, labels = percent(hypo/total),cex=1.5)
}


dmcL0vsL1<-loadDataBed('dmcL0vsL1')
dmcL0vsL2<-loadDataBed('dmcL0vsL2')
dmcL0vsL3<-loadDataBed('dmcL0vsL3')

saveImage("dmc.L0.vs.L1.pdf",width = 4.5,height = 4)
draw.scatter(dmcL0vsL1, 'L0','L1')
dev.off()
saveImage("dmc.L0.vs.L2.pdf",width = 4.5,height = 4)
draw.scatter(dmcL0vsL2, 'L0','L2')
dev.off()
saveImage("dmc.L0.vs.L3.pdf",width = 4.5,height = 4)
draw.scatter(dmcL0vsL3, 'L0','L3')
dev.off()

#----------------------------------------------------------------------------------------------------------------------
# Figure 1 d. Methylation level of CpGs within 5,000 bp upstream and downstream relative to TSS
#----------------------------------------------------------------------------------------------------------------------
smooth2 <- function(hw=51) {
  function(arr){
    for(i in 1:length(arr)){
      left<-i-hw
      if(left < 1){
        left<-1
      }
      right<-i+hw
      if(right > length(arr)){
        right<-length(arr)
      }
      
      arr[i]<-mean(arr[left:right])
    }
    arr
  }
}

plot.group.methy.profile<- function(UP,DWON,ylab="Methylation Level",cex=1.5,horiz=FALSE){
  x<-seq(UP, DWON)
  signalMatrixGroup<-loadDataBed('signalMatrixGroup')
  mmData<-signalMatrixGroup[,-1]
  mm<-mmData[(15001+UP):(15001+DWON),]
  mm<-as.data.frame(apply(mm,2,smooth2(51)))
  mm<-as.data.frame(mm)
  par(mar = c(5,5,1,1))
  plot(NA, xlim=c(UP, DWON), xlab="bp to TSS",ylab=ylab,cex.lab=2, cex.axis=1.5,ylim=c(min(mm),max(mm)))
  lines(x,mm[["L0"]],col=alpha("green3", 0.5))
  lines(x,mm[["L1"]],col=alpha("cyan", 0.5))
  lines(x,mm[["L2"]],col=alpha("orange", 0.5))
  lines(x,mm[["L3"]],col=alpha("red", 0.5))
  legend("bottomleft",legend=c("L0","L1","L2","L3"),fill=c("green3","cyan","orange","red"),bty = "n",cex=cex,horiz = horiz)
}
saveImage("methyLevel.profile.pdf",width = 6,height = 4)
plot.group.methy.profile(-5000, 5000)
dev.off()

#----------------------------------------------------------------------------------------------------------------------
# Figure S1 c. The variation of methylation level of CpGs within 5,000 bp upstream and downstream relative to TSS
#----------------------------------------------------------------------------------------------------------------------
plot.sample.methy.profile <- function(UP,DWON){
  signalMatrixSample<-loadDataBed('signalMatrixSample')
  mmData<-signalMatrixSample[,-1]
  type<-groups$WGBS$select(c('L0','L1','L2','L3'))
  samples<-colnames(mmData)
  colors<-type$colors[match(samples,type$SampleName)]
  ylab="Methylation Level"
  x<-seq(UP, DWON)
  mm<-mmData[(15001+UP):(15001+DWON),]
  par(mar = c(5,5,1,1))
  plot(NA, xlim=c(UP, DWON), ylim=c(min(mm),max(mm)), xlab="bp to TSS",ylab=ylab,cex.lab=2, cex.axis=1.5)
  
  for(i in 1:length(samples)){
    lines(x,mm[[samples[i]]],col=alpha(colors[i],0.3))
  }
  legend("bottomleft",legend=c("L0","L1","L2","L3"),fill=c("green3","cyan","orange","red"),bty = "n",cex=1,horiz = T)
}
saveImage("methyLevel.profile.1.pdf",width = 5,height = 5)
plot.sample.methy.profile(-5000,-4900)
dev.off()
saveImage("methyLevel.profile.2.pdf",width = 5,height = 5)
plot.sample.methy.profile(-1000,-0)
dev.off()
saveImage("methyLevel.profile.3.pdf",width = 5,height = 5)
plot.sample.methy.profile(4900,5000)
dev.off()

UP<--5000
DWON<-5000
x<-seq(UP, DWON)
mm<-mmData[(15001+UP):(15001+DWON),]
SD<-apply(mm,1, function(x){
  sd(x)
})
SD<--log10(SD)
saveImage("methyLevel.profile.sd.pdf",width = 7,height = 2.5)
par(mar = c(5,5,1,1))
plot(UP:DWON,SD,pch=16,cex=0.1,xlab="bp to TSS",ylab="-Log10(SD)",cex.lab=1, cex.axis=1)
dev.off()



