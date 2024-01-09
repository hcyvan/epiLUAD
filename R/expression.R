source('./R/base.R')
library(dplyr)
library(reshape2)
library(ComplexHeatmap)
library(edgeR)
library(VennDiagram)
library(circlize)
library(ggplot2)
library(GenomicRanges)
library(rtracklayer)
library(reshape2)
library(ggExtra)
library(gridExtra)
library(ggpubr)
library(clusterProfiler)
library(org.Hs.eg.db)

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

plotBarTPM<-function(geneSymbol){
  rnaTPM<-readRDS(file.path(CONFIG$dataIntermediate,'rna', 'rnaTPM.rds'))
  symbol<-ensemble2Symbol3(rnaTPM$ensemble, keepIfNotMatch = TRUE)
  gene<-unlist(rnaTPM[match(geneSymbol, symbol),][-1])
  samples<-groups$RNA$selectBySample(names(exp))
  data<-data.frame(
    stage=samples$Group,
    value=gene
  )
  plotBarError(data, xlabel=geneSymbol,ylabel="TPM")
}

#----------------------------------------------------------------------------------------------------------------------
# RNA-seq: Differential Expression Analysis
#----------------------------------------------------------------------------------------------------------------------
dat<-readRDS(file.path(CONFIG$dataIntermediate,'rna', 'rnaCountForDeg.rds'))
samples<-groups$RNA$selectBySample(colnames(dat))
salia<-samples$Group
y <- DGEList(counts = dat,genes = rownames(dat),group = salia)
keep <- rowSums(cpm(dat) > 0.5) >= 1
y <- y[keep, keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
plotMDS(y, col = as.numeric(salia))
design <- model.matrix( ~ 0 + salia)
system.time(y <- estimateDisp(y, design, robust = TRUE))
system.time(fit <- glmQLFit(y, design, robust = TRUE))
system.time(qlf.CTL.vs.AIS <- glmQLFTest(fit, contrast =  c(-1,1,0,0)))
system.time(qlf.CTL.vs.MIA <- glmQLFTest(fit, contrast =  c(-1,0,1,0)))
system.time(qlf.CTL.vs.IAC <- glmQLFTest(fit, contrast =  c(-1,0,0,1)))
cutoff<-1.5
CTL.vs.AIS<-topTags(qlf.CTL.vs.AIS, 10000)$table%>%arrange(abs(logFC))%>%filter(abs(logFC) > log2(cutoff) &FDR < 0.05)
CTL.vs.MIA<-topTags(qlf.CTL.vs.MIA, 10000)$table%>%arrange(abs(logFC))%>%filter(abs(logFC) > log2(cutoff) &FDR < 0.05)
CTL.vs.IAC<-topTags(qlf.CTL.vs.IAC, 10000)$table%>%arrange(abs(logFC))%>%filter(abs(logFC) > log2(cutoff) &FDR < 0.05)


table(CTL.vs.AIS$logFC>0)
table(CTL.vs.MIA$logFC>0)
table(CTL.vs.IAC$logFC>0)

deg<-list(
  AIS=CTL.vs.AIS,
  MIA=CTL.vs.MIA,
  IAC=CTL.vs.IAC
)
sapply(deg, nrow)
degAll<-list(
  AIS=qlf.CTL.vs.AIS,
  MIA=qlf.CTL.vs.MIA,
  IAC=qlf.CTL.vs.IAC
)
saveRDS(deg, file.path(CONFIG$dataIntermediate,'rna', 'deg.rds'))
saveRDS(degAll, file.path(CONFIG$dataIntermediate,'rna', 'deg.all.rds'))

write.csv(CTL.vs.AIS,file.path(CONFIG$dataIntermediate,'rna', 'rna.deg.CTL.vs.AIS.csv'),row.names =FALSE)
write.csv(CTL.vs.MIA,file.path(CONFIG$dataIntermediate,'rna', 'rna.deg.CTL.vs.MIA.csv'),row.names =FALSE)
write.csv(CTL.vs.IAC,file.path(CONFIG$dataIntermediate,'rna', 'rna.deg.CTL.vs.IAC.csv'),row.names =FALSE)

degDF<-do.call(rbind,list(data.frame(deg$AIS, group='AISvsCTL'),data.frame(deg$MIA, group='MIAvsCTL'),data.frame(deg$IAC, group='IACvsCTL')))
write.csv(degDF,file.path(CONFIG$dataResult, 'rna.deg.detail.csv'),row.names =FALSE)

#----------------------------------------------------------------------------------------------------------------------
# Figure 4A. DEG volcano plot
#----------------------------------------------------------------------------------------------------------------------
degAll<-readRDS(file.path(CONFIG$dataIntermediate,'rna', 'deg.all.rds'))
plot.Volcano <- function(qlf,main){
  dd <-data.frame(qlf$table,FDR = p.adjust(qlf$table$PValue, method = 'BH'),genes = qlf$genes)
  colors<-sapply(split(dd,seq(nrow(dd))),function(x){
    if(x$logFC > log2(1.5) && x$FDR < 0.05){
      'red'
    }else if(x$logFC < -log2(1.5) && x$FDR < 0.05) {
      'green3'
    }else{
      'black'
    }
  })
  str(colors)
  plot(dd$logFC,-log2(dd$FDR),pch=20,cex=0.5,col=colors,cex.axis=1.3,cex.lab=1.3,xlab = "log2FC",ylab="-log2FDR",main=main)
  text(x=6,y=0.5,labels=sum(colors=="red"),cex=1.5,col="red")
  text(x=-6,y=0.5,labels=sum(colors=="green3"),cex=1.5,col="green3")
  legend("topright", legend=c("UP","DOWN"),col = c("red","green3"),bty = "n",pch = 16,cex=0.8)
}

saveImage2("rna.deg.volcano.pdf",width = 8,height = 2.5)
par(mfrow=c(1,3),mar = c(4,5,2,1))
plot.Volcano(degAll$AIS, main="AIS vs CTL")
plot.Volcano(degAll$MIA, main="MIA vs CTL")
plot.Volcano(degAll$IAC, main="IAC vs CTL")
dev.off()
#----------------------------------------------------------------------------------------------------------------------
# Figure 4B. DEG Venn
#----------------------------------------------------------------------------------------------------------------------
deg<-readRDS(file.path(CONFIG$dataIntermediate,'rna', 'deg.rds'))
# deg<-lapply(deg, function(x){
#   x[x$logFC>0,]
# })
a<-deg$AIS$genes
b<-deg$MIA$genes
c<-deg$IAC$genes
venn.plot<-venn.diagram(x = list(AIS = a, MIA = b,IAC =c),
                        fill = colorMapStage[2:4],
                        cat.cex = 2,
                        cex = 2,
                        filename = NULL)
saveImage2("rna.deg.venn.pdf",width = 5,height = 5)
grid.draw(venn.plot)
dev.off()
#----------------------------------------------------------------------------------------------------------------------
# Figure 4C
#----------------------------------------------------------------------------------------------------------------------
deg<-readRDS(file.path(CONFIG$dataIntermediate,'rna', 'deg.rds'))
degStage<-list(
  AIS=deg$AIS$genes,
  MIA=setdiff(deg$MIA$genes, deg$AIS$genes),
  IAC=setdiff(deg$IAC$genes, union(deg$AIS$genes,deg$MIA$genes))
)
sapply(degStage, function(x){length(x)})
saveRDS(degStage, file.path(CONFIG$dataIntermediate,'rna', 'deg.stage.rds'))
plot.degStage<-function(gene) {
  rnaTPM<-readRDS(file.path(CONFIG$dataIntermediate,'rna', 'rnaTPM.rds'))
  symbol<-ensemble2Symbol3(rnaTPM$ensemble, keepIfNotMatch = TRUE)
  keep<-match(gene, symbol)[!is.na(match(gene, symbol))]
  m<-as.matrix(rnaTPM[keep,-1])
  rownames(m)<-symbol[keep]
  mm<-t(scale(t(m)))
  mm<-mm[!is.na(mm[,1]),]
  samples<-groups$RNA$selectBySample(colnames(mm))
  column_annotation <-HeatmapAnnotation(
    df=data.frame(Stage=samples$Group),
    col = list(Stage =colorMapStage),
    show_annotation_name =FALSE,
    annotation_name_side='left'
  )
  Heatmap(mm,
          top_annotation = column_annotation,
          cluster_rows=TRUE,
          cluster_columns = FALSE,
          show_row_names=FALSE,
          show_column_names=FALSE,
          heatmap_legend_param = list(
            title = "Scaled TPM",
            legend_height = unit(4, "cm"),
            at = c(0,0.5,1),
            title_position = "lefttop-rot"
          ),
  )
}
saveImage2("rna.deg.SR.AIS.heatmap.pdf",width = 4.5,height = 3)
plot.degStage(degStage$AIS)
dev.off()
saveImage2("rna.deg.SR.MIA.heatmap.pdf",width = 4.5,height = 3)
plot.degStage(degStage$MIA)
dev.off()
saveImage2("rna.deg.SR.IAC.heatmap.pdf",width = 4.5,height = 3)
plot.degStage(degStage$IAC)
dev.off()


#----------------------------------------------------------------------------------------------------------------------
# Figure S3E Expression of 24 significant differential expression TFs
#----------------------------------------------------------------------------------------------------------------------
tf124gene<-loadData2(file.path(CONFIG$dataResult,'tf124.symbol.csv'))
deg<-readRDS(file.path(CONFIG$dataIntermediate,'rna', 'deg.rds'))
tf124<-loadData2(file.path(CONFIG$dataIntermediate,'tf', 'tf124.txt'),file.format = 'bed',header = FALSE)$V1
geneTf124<-na.omit(unique(tf124gene$gene))
tf124DEG<-lapply(deg,function(x){
  intersect(x$genes, geneTf124)
})
tf124DEGList<-unique(unlist(tf124DEG))
p1<-plotBarTPM(tf124DEGList[1])
p2<-plotBarTPM(tf124DEGList[2])
p3<-plotBarTPM(tf124DEGList[3])
p4<-plotBarTPM(tf124DEGList[4])
p5<-plotBarTPM(tf124DEGList[5])
p6<-plotBarTPM(tf124DEGList[6])
p7<-plotBarTPM(tf124DEGList[7])
p8<-plotBarTPM(tf124DEGList[8])
p9<-plotBarTPM(tf124DEGList[9])
p10<-plotBarTPM(tf124DEGList[10])
p11<-plotBarTPM(tf124DEGList[11])
p12<-plotBarTPM(tf124DEGList[12])
p13<-plotBarTPM(tf124DEGList[13])
p14<-plotBarTPM(tf124DEGList[14])
p15<-plotBarTPM(tf124DEGList[15])
p16<-plotBarTPM(tf124DEGList[16])
p17<-plotBarTPM(tf124DEGList[17])
p18<-plotBarTPM(tf124DEGList[18])
p19<-plotBarTPM(tf124DEGList[19])
p20<-plotBarTPM(tf124DEGList[20])
p21<-plotBarTPM(tf124DEGList[21])
p22<-plotBarTPM(tf124DEGList[22])
p23<-plotBarTPM(tf124DEGList[23])
p24<-plotBarTPM(tf124DEGList[24])

saveImage2("tf124.deg.tf24.pdf",width = 14,height = 8)
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,nrow = 4)
dev.off()

#----------------------------------------------------------------------------------------------------------------------
# Figure 4D,E,F DEGs Enrichment Analysis
#----------------------------------------------------------------------------------------------------------------------
degStage<-readRDS(file.path(CONFIG$dataIntermediate,'rna', 'deg.stage.rds'))
doErich<-function(genes){
  pcg = bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  ego.CC <- enrichGO(gene          = pcg$ENTREZID,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "CC",
                     pvalueCutoff  = 0.005,
                     qvalueCutoff  = 0.005,
                     readable      = TRUE)
  ego.MF <- enrichGO(gene          = pcg$ENTREZID,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "MF",
                     pvalueCutoff  = 0.005,
                     qvalueCutoff  = 0.005,
                     readable      = TRUE)
  ego.BP <- enrichGO(gene          = pcg$ENTREZID,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "BP",
                     pvalueCutoff  = 0.005,
                     qvalueCutoff  = 0.005,
                     readable      = TRUE)
  
  kegg <- enrichKEGG(gene         = pcg$ENTREZID,
                           organism     = 'hsa',
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.005,
                           qvalueCutoff  = 0.005)
  
  ego<-list(
    cc=simplify(ego.CC,cutoff=0.7, by="p.adjust", select_fun=min),
    mf=simplify(ego.MF,cutoff=0.7, by="p.adjust", select_fun=min),
    bp=simplify(ego.BP,cutoff=0.7, by="p.adjust", select_fun=min),
    kegg=kegg
  )
  ego
}
goAIS<-doErich(degStage$AIS)
goMIA<-doErich(degStage$MIA)
goIAC<-doErich(degStage$IAC)


outputGo<-function(go,out){
  enrich<-head(go,n=1000)
  write.csv(enrich, file.path(CONFIG$dataIntermediate, "rna", out), row.names = FALSE,quote = FALSE)
}
outputGo(goAIS$bp, 'deg.AIS.go.bp.csv')
outputGo(goMIA$bp, 'deg.MIA.go.bp.csv')
outputGo(goIAC$bp, 'deg.IAC.go.bp.csv')
outputGo(goAIS$kegg, 'deg.AIS.kegg.csv')
outputGo(goMIA$kegg, 'deg.MIA.kegg.csv')
outputGo(goIAC$kegg, 'deg.IAC.kegg.csv')

plot.clusterprofiler<-function(go,title=""){
  enrich<-head(go,n=15)
  plotEnrich(enrich$Description, enrich$qvalue,title=title)

}
plot.clusterprofiler(goAIS$cc)
saveImage2("rna.deg.SR.go.AIS.BP.pdf",width = 4,height = 5)
plot.clusterprofiler(goAIS$bp)
dev.off()
plot.clusterprofiler(goAIS$mf)

plot.clusterprofiler(goMIA$cc)
saveImage2("rna.deg.SR.go.MIA.BP.pdf",width = 4,height = 4)
plot.clusterprofiler(goMIA$bp)
dev.off()
plot.clusterprofiler(goMIA$mf)

plot.clusterprofiler(goIAC$cc)
saveImage2("rna.deg.SR.go.IAC.BP.pdf",width = 4,height = 5)
plot.clusterprofiler(goIAC$bp)
dev.off()
plot.clusterprofiler(goIAC$mf)
#----------------------------------------------------------------------------------------------------------------------
# Figure 4G TF-DEG pair correlation
#----------------------------------------------------------------------------------------------------------------------
deg<-readRDS(file.path(CONFIG$dataIntermediate,'rna', 'deg.rds'))
rnaTPM<-readRDS(file.path(CONFIG$dataIntermediate,'rna', 'rnaTPM.rds'))
symbol<-ensemble2Symbol3(rnaTPM$ensemble, keepIfNotMatch = TRUE)
tf124gene<-loadData2(file.path(CONFIG$dataResult,'tf124.symbol.csv'))
geneTf124<-na.omit(unique(tf124gene$gene))
genes<-c()
genes<-unique(c(genes, unlist(sapply(deg, function(x){x$genes}))))
genes<-unique(c(genes, geneTf124))
keep<-match(genes, symbol)[!is.na(match(genes, symbol))]
m<-as.matrix(rnaTPM[keep,-1])
rownames(m)<-symbol[keep]
mm<-m
mm<-mm[!is.na(mm[,1]),]
cor(t(mm))->a
a[lower.tri(a)]<-1
df0<-melt(a)
df0<-filter(df0, value<1)
df<-filter(df0,abs(value)>=0.7)
df<-df[(df$Var1%in%geneTf124|df$Var2%in%geneTf124),]
df$Var1<-as.vector(df$Var1)
df$Var2<-as.vector(df$Var2)
df$tf<-ifelse(df$Var1%in%geneTf124, df$Var1, df$Var2)
df$gene<-ifelse(df$Var1%in%geneTf124, df$Var2, df$Var1)
tfCorNumber<-sapply(split(df, df$tf),function(x){
  nrow(x)
})
tfCorNumber<-sort(tfCorNumber, decreasing = TRUE)
saveImage2("rna.deg.tf124.coexp.pdf",width = 14,height = 4)
barplot(tfCorNumber,ylab="Co-expression gene number",las=2,cex.names =0.8, col='red')
dev.off()
















############
rnaTPM<-readRDS(file.path(CONFIG$dataIntermediate,'rna', 'rnaTPM.rds'))
symbol<-ensemble2Symbol3(rnaTPM$ensemble, keepIfNotMatch = TRUE)
keep<-match(geneTf124, symbol)[!is.na(match(geneTf124, symbol))]
tfs<-rnaTPM[keep,]
m<-as.matrix(tfs[,-1])
rownames(m)<-symbol[keep]
mm<-t(scale(t(m)))
mm<-mm[!is.na(mm[,1]),]
saveImage2("tf124.exp.heatmap.pdf",width = 8,height = 16)
Heatmap(mm,
        heatmap_legend_param = list(title = "Scaled TPM"),
        cluster_rows=TRUE,
        cluster_columns = FALSE,
        show_row_names=TRUE,
        show_column_names=TRUE)
dev.off()
#########
degStageMAPK<-readRDS(file.path(CONFIG$dataIntermediate,'rna', 'degStageMAPK.rds'))



geneList<-degStageMAPK$MIA$gene
geneList<-unlist(lapply(degStageMAPK, function(x){x$gene}))

rnaTPM<-readRDS(file.path(CONFIG$dataIntermediate,'rna', 'rnaTPM.rds'))
symbol<-ensemble2Symbol3(rnaTPM$ensemble, keepIfNotMatch = TRUE)
keep<-match(geneList, symbol)[!is.na(match(geneList, symbol))]
tfs<-rnaTPM[keep,]
m<-as.matrix(tfs[,-1])
rownames(m)<-symbol[keep]
mm<-t(scale(t(m)))
# mm<-m
mm<-mm[!is.na(mm[,1]),]
Heatmap(mm,
        heatmap_legend_param = list(title = "Scaled TPM"),
        cluster_rows=TRUE,
        cluster_columns = FALSE,
        show_row_names=TRUE,
        show_column_names=TRUE)


geneList<-c('RAC3','RRAS','RRAS2','MRAS','PDGFRB','CACNA1C')
plotBarTPM('RAC3')
plotBarTPM('TGFA')
plotBarTPM('ANGPT4')
plotBarTPM('RASGRP2')
plotBarTPM('RASGRF2')
plotBarTPM('PRKCB')

# RTK
plotBarTPM('MET') # up;;target
plotBarTPM('NTRK3')# down;NTRK2 sub; target
plotBarTPM('ERBB2')# up;HER2 alias; target
plotBarTPM('ERBB3')# up;MIA
plotBarTPM('RET')# no;; target
plotBarTPM('ALK')# no;; target;ALK -> RAS/ERK
plotBarTPM('EGFR')# up;MIA;target
# GF
plotBarTPM('FGF2')# down
plotBarTPM('FGF1')# up
plotBarTPM('TGFA')# up
plotBarTPM('ANGPT4')# down
plotBarTPM('IGF1')# down;MIA
plotBarTPM('EGF')# up;MIA
plotBarTPM('PGF')# down;MIA

# 

plotBarTPM('BRAF')# up;;target
plotBarTPM('DUSP6')# up;;target
#
plotBarTPM('PDCD1')# no;;target
plotBarTPM('CD274')# no;;target

############


geneList<-c('EGFR','IL1R2','BRAF','PDCD1','CD274','JUN','FOSL1','FOSL2','JUND','BATF','FOS','MMP9')

rnaTPM<-readRDS(file.path(CONFIG$dataIntermediate,'rna', 'rnaTPM.rds'))
symbol<-ensemble2Symbol3(rnaTPM$ensemble, keepIfNotMatch = TRUE)
keep<-match(geneList, symbol)[!is.na(match(geneList, symbol))]
tfs<-rnaTPM[keep,]
m<-as.matrix(tfs[,-1])
rownames(m)<-symbol[keep]
mm<-t(scale(t(m)))
# mm<-m
mm<-mm[!is.na(mm[,1]),]
Heatmap(mm,
        heatmap_legend_param = list(title = "Scaled TPM"),
        cluster_rows=TRUE,
        cluster_columns = FALSE,
        show_row_names=TRUE,
        show_column_names=TRUE)
heatmap(abs(cor(t(mm),method = 'spearman')))

cor(t(mm))

# a<-mm[match('JUN',rownames(mm)),]
# b<-mm[match('FOS',rownames(mm)),]
# plot(a,b)
#######################



geneList<-c('EGFR','IL1R2','BRAF','PDCD1','CD274')
rnaTPM<-readRDS(file.path(CONFIG$dataIntermediate,'rna', 'rnaTPM.rds'))
symbol<-ensemble2Symbol3(rnaTPM$ensemble, keepIfNotMatch = TRUE)
keep<-match(geneList, symbol)[!is.na(match(geneList, symbol))]
tfs<-rnaTPM[keep,]
m<-as.matrix(tfs[,-1])
rownames(m)<-symbol[keep]
mm<-m
samples<-groups$RNA$selectBySample(colnames(mm))
mmSingle<-lapply(1:nrow(mm), function(i){
  data.frame(
    value=mm[i,],
    stage=samples$Group
  )
})
names(mmSingle)<-rownames(mm)

plotBarTPM(mmSingle$EGFR, 'EGFR')
plotBarTPM(mmSingle$IL1R2, 'IL1R2')
plotBarTPM(mmSingle$BRAF, 'BRAF')
plotBarTPM(mmSingle$PDCD1, 'PDCD1')
plotBarTPM(mmSingle$CD274, 'CD274')


