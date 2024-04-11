source('./R/base.R')
library(edgeR)
library(stats)
library(VennDiagram)

plotBarErrorOfA549T5224<-function(gene){
  tpm<-readRDS(file.path(CONFIG$dataIntermediate,'rna_a549', 'rnaTPM.rds'))
  genes<-ensemble2Symbol3(tpm$ensemble,TRUE)
  line<-tpm[match(gene, genes),]
  valueStage <- data.frame(
    value=unlist(line[,2:10]),
    stage=factor(c(rep('0uM',3),rep('100uM',3),rep('200uM',3)))
  )
  xlabel=""
  ylabel="TPM"
  colorMapStageTmp<-c('#3587c8','#e73217','#66b62e')
  names(colorMapStageTmp)<-c('0uM','100uM','200uM')
  ggplot(data=valueStage,aes(x=stage,y=value,fill=stage))+
    scale_fill_manual(values=colorMapStageTmp) +
    stat_summary(mapping=aes(fill = stage),fun=mean,geom = "bar",fun.args = list(mult=1),width=0.7)+
    stat_summary(fun.data=mean_sdl,fun.args = list(mult=0.5),geom="errorbar",width=0.2) +
    labs(x=xlabel,y=ylabel)+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    theme(legend.position="none",
          axis.title.x = element_text(size=15),
          axis.title.y = element_text(size=15),
          axis.text = element_text(size = 12,colour="black"))+
    guides(colour = guide_legend(override.aes = list(shape = 12,size=10)))
}
#----------------------------------------------------------------------------------------------------------------------
# RNA-seq: Differential Expression Analysis of RNA-seq data of a549 cell line treated with different concentrations of T5224
#----------------------------------------------------------------------------------------------------------------------
count<-readRDS(file.path(CONFIG$dataIntermediate,'rna_a549', 'rnaCount.rds'))
genes<-ensemble2Symbol3(count$ensemble)
count<-count[which(genes!=""),]
dat<-count[,-1]
genes<-ensemble2Symbol3(count$ensemble)
salia<-c('s0','s0','s0','s100','s100','s100','s200','s200','s200')
y <- DGEList(counts = dat,genes = genes,group = salia)
keep <- rowSums(cpm(dat) > 0.5) >= 1
y <- y[keep, keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
plotMDS(y, col = as.numeric(salia))
design <- model.matrix( ~ 0 + salia)
system.time(y <- estimateDisp(y, design, robust = TRUE))
system.time(fit <- glmQLFit(y, design, robust = TRUE))
system.time(qlf.0.vs.100 <- glmQLFTest(fit, contrast =  c(-1,1,0)))
system.time(qlf.0.vs.200 <- glmQLFTest(fit, contrast =  c(-1,0,1)))
cutoff<-1.15
diff.0.vs.100<-topTags(qlf.0.vs.100, 10000)$table%>%arrange(abs(logFC))%>%filter(abs(logFC) > log2(cutoff) &PValue < 0.05)
diff.0.vs.200<-topTags(qlf.0.vs.200, 10000)$table%>%arrange(abs(logFC))%>%filter(abs(logFC) > log2(cutoff) &PValue< 0.05)

d100<-table(diff.0.vs.100$logFC>0)
d200<-table(diff.0.vs.200$logFC>0)

deg<-list(
  c100=diff.0.vs.100,
  c200=diff.0.vs.200
)
#----------------------------------------------------------------------------------------------------------------------
# Figure S5D. DEGs in a549 cell line treated with different concentrations of T5224
#----------------------------------------------------------------------------------------------------------------------
df<-data.frame(
  T5224=c(rep('100uM',2),rep('200uM',2)),
  DEG=factor(c('Up','Down','Up','Down'),levels = c('Up','Down')),
  Number=c(d100[1],d100[2],d200[1],d200[2])
)
saveImage2("rna.deg.a549.t5224.pdf",width = 2.5,height = 2.5)
ggplot(data=df, aes(x=T5224, y=Number, fill=DEG)) +
  geom_bar(stat="identity")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
#----------------------------------------------------------------------------------------------------------------------
# Figure S5E. Venn plot of AP-1 target genes and DEGs in a549 
#----------------------------------------------------------------------------------------------------------------------
ap1Target<-readRDS(file.path(CONFIG$dataIntermediate,'mapk', 'ap1.targets.rds'))
a<-ap1Target$Gene
b<-diff.0.vs.100$genes
c<-diff.0.vs.200$genes
cmp<-c('#ffb90f','#e73217','#66b62e')
names(cmp)<-c('AP-1', '100uM', '200uM')
venn.plot<-venn.diagram(x = list(`AP1` = a, `100uM` = b,`200uM` =c),
                        fill = cmp,
                        cat.cex = 2,
                        cex = 2,
                        filename = NULL)
saveImage2("rna.deg.venn.ap1.target.a549.t5224.pdf",width = 5,height = 5)
grid.draw(venn.plot)
dev.off()
#----------------------------------------------------------------------------------------------------------------------
# Figure 7E. FGFBP1 and EGFR in a549 cell line treated with different concentrations of T5224
#----------------------------------------------------------------------------------------------------------------------
saveImage2("rna.ap1.target.a549.FGFBP1.pdf",width = 2,height = 2.5)
plotBarErrorOfA549T5224('FGFBP1')
dev.off()
saveImage2("rna.ap1.target.a549.EGFR.pdf",width = 2,height = 2.5)
plotBarErrorOfA549T5224('EGFR')
dev.off()
## Statistical Test
gene<-'FGFBP1'
tpm<-readRDS(file.path(CONFIG$dataIntermediate,'rna_a549', 'rnaTPM.rds'))
genes<-ensemble2Symbol3(tpm$ensemble,TRUE)
line<-tpm[match(gene, genes),]
wilcox.test(unlist(line[,2:4]), unlist(line[,8:10]),paired = FALSE,alternative = 'greater')
wilcox.test(unlist(line[,2:4]), unlist(line[,5:7]),paired = FALSE,alternative = 'greater')
data <- data.frame(
  y=unlist(line[,2:10]),
  group=factor(c(1,1,1,2,2,2,3,3,3))
)
shapiro.test(data$y)
bartlett.test(y ~ group, data = data)
summary(aov(y~group,data=data))

gene<-'EGFR'
tpm<-readRDS(file.path(CONFIG$dataIntermediate,'rna_a549', 'rnaTPM.rds'))
genes<-ensemble2Symbol3(tpm$ensemble,TRUE)
line<-tpm[match(gene, genes),]
wilcox.test(unlist(line[,2:4]), unlist(line[,8:10]),paired = FALSE,alternative = 'greater')
wilcox.test(unlist(line[,2:4]), unlist(line[,5:7]),paired = FALSE,alternative = 'greater')

data <- data.frame(
  y=unlist(line[,2:10]),
  group=factor(c(1,1,1,2,2,2,3,3,3))
)
shapiro.test(data$y)
bartlett.test(y ~ group, data = data)
summary(aov(y~group,data=data))
#----------------------------------------------------------------------------------------------------------------------
# Figure 7F. DEGs in IAC and A549 with T-5224
#----------------------------------------------------------------------------------------------------------------------
deg<-loadSRDEG()
diff<-inner_join(diff.0.vs.100, deg$detail$IAC, by='genes')
uu<-filter(diff, logFC.x>0,logFC.y>0)%>%nrow
dd<-filter(diff, logFC.x<0,logFC.y<0)%>%nrow
ud<-filter(diff, logFC.x>0,logFC.y<0)%>%nrow
du<-filter(diff, logFC.x<0,logFC.y>0)%>%nrow
cat(paste('uu', uu))
cat(paste('dd', dd))
cat(paste('ud', ud))
cat(paste('du', du))

diff<-inner_join(diff.0.vs.200, deg$detail$IAC, by='genes')
uu<-filter(diff, logFC.x>0,logFC.y>0)%>%nrow
dd<-filter(diff, logFC.x<0,logFC.y<0)%>%nrow
ud<-filter(diff, logFC.x>0,logFC.y<0)%>%nrow
du<-filter(diff, logFC.x<0,logFC.y>0)%>%nrow
cat(paste('uu', uu))
cat(paste('dd', dd))
cat(paste('ud', ud))
cat(paste('du', du))
#----------------------------------------------------------------------------------------------------------------------
# Figure S5G. MMP3 and MMP7 in a549 cell line treated with different concentrations of T5224
#----------------------------------------------------------------------------------------------------------------------
saveImage2("rna.ap1.target.a549.MMP3.pdf",width = 2,height = 2.5)
plotBarErrorOfA549T5224('MMP3')
dev.off()
saveImage2("rna.ap1.target.a549.MMP7.pdf",width = 2,height = 2.5)
plotBarErrorOfA549T5224('MMP7')
dev.off()
## Statistical Test
gene<-'MMP3'
tpm<-readRDS(file.path(CONFIG$dataIntermediate,'rna_a549', 'rnaTPM.rds'))
genes<-ensemble2Symbol3(tpm$ensemble,TRUE)
line<-tpm[match(gene, genes),]
wilcox.test(unlist(line[,2:4]), unlist(line[,8:10]),paired = FALSE,alternative = 'greater')
wilcox.test(unlist(line[,2:4]), unlist(line[,5:7]),paired = FALSE,alternative = 'greater')
data <- data.frame(
  y=unlist(line[,2:10]),
  group=factor(c(1,1,1,2,2,2,3,3,3))
)
shapiro.test(data$y)
bartlett.test(y ~ group, data = data)
summary(aov(y~group,data=data))

gene<-'MMP7'
tpm<-readRDS(file.path(CONFIG$dataIntermediate,'rna_a549', 'rnaTPM.rds'))
genes<-ensemble2Symbol3(tpm$ensemble,TRUE)
line<-tpm[match(gene, genes),]
wilcox.test(unlist(line[,2:4]), unlist(line[,8:10]),paired = FALSE,alternative = 'greater')
wilcox.test(unlist(line[,2:4]), unlist(line[,5:7]),paired = FALSE,alternative = 'greater')

data <- data.frame(
  y=unlist(line[,2:10]),
  group=factor(c(1,1,1,2,2,2,3,3,3))
)
shapiro.test(data$y)
bartlett.test(y ~ group, data = data)
summary(aov(y~group,data=data))




