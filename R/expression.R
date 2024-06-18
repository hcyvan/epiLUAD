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

srdegDetail<-lapply(names(degStage), function(x){
  genes<-degStage[[x]]
  DEG<-deg[[x]]
  DEG[match(genes, DEG$genes),]
})
names(srdegDetail)<-names(degStage)
srdegDetail$AIS$class<-ifelse(srdegDetail$AIS$logFC>0, 'UpInAIS','DownInAIS')
srdegDetail$MIA$class<-ifelse(srdegDetail$MIA$logFC>0, 'UpInMIA','DownInMIA')
srdegDetail$IAC$class<-ifelse(srdegDetail$IAC$logFC>0, 'UpInIAC','DownInIAC')

SRDEG<-list(
  deg=degStage,
  detail=srdegDetail
)
saveRDS(SRDEG, file.path(CONFIG$dataIntermediate,'rna', 'srdeg.rds'))
saveRDS(degStage, file.path(CONFIG$dataIntermediate,'rna', 'deg.stage.rds'))
plot.degStage<-function(gene, do_top_annotation=FALSE) {
  rnaTPM<-readRDS(file.path(CONFIG$dataIntermediate,'rna', 'rnaTPM.rds'))
  symbol<-ensemble2Symbol3(rnaTPM$ensemble, keepIfNotMatch = TRUE)
  keep<-match(gene, symbol)[!is.na(match(gene, symbol))]
  m<-as.matrix(rnaTPM[keep,-1])
  rownames(m)<-symbol[keep]
  mm<-t(scale(t(m)))
  mm<-mm[!is.na(mm[,1]),]
  samples<-groups$RNA$selectBySample(colnames(mm))
  if (do_top_annotation) {
    column_annotation <-HeatmapAnnotation(
      df=data.frame(Stage=samples$Group),
      col = list(Stage =colorMapStage2),
      show_annotation_name =FALSE,
      annotation_name_side='left'
    )
  }else{
    column_annotation<-NULL
  }
  mm<-2 * (mm - min(mm)) / (max(mm) - min(mm)) - 1
  print(sprintf("max: %f, min: %f", max(mm), min(mm)))
  Heatmap(mm,
          # col=colorRamp2(c(-1,-0.5,0,0.5 ,1), c("#333cac","#1d9cbb","#c5ba59", "#f5cc2f","#f6f803")),
          col=colorRamp2(c(-1, -0.6, 1), c("#0808ff", "#ffffff", "#ff0808")),
          top_annotation = column_annotation,
          cluster_rows=TRUE,
          cluster_columns = FALSE,
          show_row_names=FALSE,
          show_column_names=FALSE,
          heatmap_legend_param = list(
            title = "Scaled TPM",
            legend_height = unit(4, "cm"),
            title_position = "lefttop-rot"
          ),
  )
}
saveImage2("rna.deg.SR.AIS.heatmap.pdf",width = 3,height = 3)
up<-plot.degStage(filter(SRDEG$detail$AIS,logFC>0)$genes, do_top_annotation=TRUE)
down<-plot.degStage(filter(SRDEG$detail$AIS,logFC<0)$genes)
up%v%down
dev.off()
saveImage2("rna.deg.SR.MIA.heatmap.pdf",width = 3,height = 3)
up<-plot.degStage(filter(SRDEG$detail$MIA,logFC>0)$genes, do_top_annotation=TRUE)
down<-plot.degStage(filter(SRDEG$detail$MIA,logFC<0)$genes)
up%v%down
dev.off()
saveImage2("rna.deg.SR.IAC.heatmap.pdf",width = 3,height = 3)
up<-plot.degStage(filter(SRDEG$detail$IAC,logFC>0)$genes, do_top_annotation=TRUE)
down<-plot.degStage(filter(SRDEG$detail$IAC,logFC<0)$genes)
up%v%down
dev.off()

#----------------------------------------------------------------------------------------------------------------------
# Figure S3B Expression of significant differential expression epiTFs
#----------------------------------------------------------------------------------------------------------------------
rnaTPM<-RnaTPM('RNA')
epiTFs<-readRDS(file.path(CONFIG$dataIntermediate,'tf', 'epiTFs.rds'))
deg<-readRDS(file.path(CONFIG$dataIntermediate,'rna', 'deg.rds'))
genes<-lapply(deg,function(x){
  intersect(x$genes, epiTFs$tf)
})%>%unlist()%>%unique()

p1<-rnaTPM$plotStageBar(genes[1])
p2<-rnaTPM$plotStageBar(genes[2])
p3<-rnaTPM$plotStageBar(genes[3])
p4<-rnaTPM$plotStageBar(genes[4])
p5<-rnaTPM$plotStageBar(genes[5])
p6<-rnaTPM$plotStageBar(genes[6])
p7<-rnaTPM$plotStageBar(genes[7])
p8<-rnaTPM$plotStageBar(genes[8])
p9<-rnaTPM$plotStageBar(genes[9])
p10<-rnaTPM$plotStageBar(genes[10])
p11<-rnaTPM$plotStageBar(genes[11])
p12<-rnaTPM$plotStageBar(genes[12])
p13<-rnaTPM$plotStageBar(genes[13])
p14<-rnaTPM$plotStageBar(genes[14])
p15<-rnaTPM$plotStageBar(genes[15])
p16<-rnaTPM$plotStageBar(genes[16])
p17<-rnaTPM$plotStageBar(genes[17])
p18<-rnaTPM$plotStageBar(genes[18])
p19<-rnaTPM$plotStageBar(genes[19])
p20<-rnaTPM$plotStageBar(genes[20])

saveImage2("rna.deg.epiTFs.barplot.pdf",width = 8,height = 8)
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,nrow = 5)
dev.off()
#----------------------------------------------------------------------------------------------------------------------
# Figure 4D,E,F DEGs Enrichment Analysis
#----------------------------------------------------------------------------------------------------------------------
degStage<-readRDS(file.path(CONFIG$dataIntermediate,'rna', 'deg.stage.rds'))
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
  enrich<-head(go,n=8)
  plotEnrich2(enrich$Description, enrich$qvalue,title=title)

}
# plot.clusterprofiler(goAIS$cc)
saveImage2("rna.deg.SR.go.AIS.BP.pdf",width = 4,height = 4)
plot.clusterprofiler(goAIS$bp)
dev.off()
# plot.clusterprofiler(goAIS$mf)

# plot.clusterprofiler(goMIA$cc)
saveImage2("rna.deg.SR.go.MIA.BP.pdf",width = 4,height = 4)
plot.clusterprofiler(goMIA$bp)
dev.off()
# plot.clusterprofiler(goMIA$mf)

# plot.clusterprofiler(goIAC$cc)
saveImage2("rna.deg.SR.go.IAC.BP.pdf",width = 4.5,height = 4)
plot.clusterprofiler(goIAC$bp)
dev.off()
# plot.clusterprofiler(goIAC$mf)
#----------------------------------------------------------------------------------------------------------------------
# Figure 4G TF-DEG pair correlation
#----------------------------------------------------------------------------------------------------------------------
deg<-readRDS(file.path(CONFIG$dataIntermediate,'rna', 'deg.rds'))
srdeg<-loadSRDEG()
rnaTPM<-readRDS(file.path(CONFIG$dataIntermediate,'rna', 'rnaTPM.rds'))
symbol<-ensemble2Symbol3(rnaTPM$ensemble, keepIfNotMatch = TRUE)
epiTFs<-readRDS(file.path(CONFIG$dataIntermediate,'tf', 'epiTFs.rds'))
#epiTFsDEG<-lapply(deg,function(x){intersect(x$genes, epiTFs$tf)})%>%unlist()%>%unique()
epiTFsDEG<-do.call(rbind,lapply(srdeg$detail,function(x){
  a<-filter(x,genes%in%epiTFs$tf)
  data.frame(genes=a$genes, deg=ifelse(startsWith(a$class, 'Up'), 'up', 'down'))
}))



genesTF<-unique(epiTFs$tf)
genesDEG<-unique(unlist(sapply(deg,function(x){x$genes})))


rnaTPM<-RnaTPM()
mTF<-rnaTPM$getTPM(genesTF,rmIfcontainNA=TRUE)
mDEG<-rnaTPM$getTPM(genesDEG,rmIfcontainNA=TRUE)

mCorP<-sapply(data.frame(t(mTF),check.names = FALSE), function(x){
  sapply(data.frame(t(mDEG)), function(y){
    t<-cor.test(x,y)
    t$p.value
  })
})

mCorCor<-sapply(data.frame(t(mTF), check.names = FALSE), function(x){
  sapply(data.frame(t(mDEG)), function(y){
    t<-cor.test(x,y)
    t$estimate
  })
})

dfCorP<-melt(mCorP)
dfCorCor<-melt(mCorCor)

dfCor<-data.frame(TF=dfCorP$Var2, DEG=dfCorP$Var1, cor=dfCorCor$value, p=dfCorP$value)
dfCorTotal<-filter(dfCor,p<=0.05)
tfCorNumber<-sapply(split(dfCorTotal, dfCorTotal$TF),function(x){nrow(x)})
dfCorPositive<-filter(dfCor,p<0.05, cor>0)
dfCorNegative<-filter(dfCor,p<0.05, cor<0)
tfCorNumberPositive<-sapply(split(dfCorPositive, dfCorPositive$TF),function(x){nrow(x)})
tfCorNumberNegative<-sapply(split(dfCorNegative, dfCorNegative$TF),function(x){nrow(x)})


df<-sapply(split(dfCorTotal, dfCorTotal$TF),function(x){
  t<-nrow(x)
  p<-sum(x$cor>0)
  n<-sum(x$cor<0)
  c(t,p,n)
})%>%t
dfPlot<-data.frame(tf=rownames(df), total=df[,1], positive=df[,2], negative=df[,3])
dfPlot<-arrange(dfPlot, desc(total))


tfCorNumberNegative<-dfPlot$negative
tfCorNumberPositive<-dfPlot$positive
names(tfCorNumberNegative)<-rownames(dfPlot)
names(tfCorNumberPositive)<-rownames(dfPlot)
yMax<-max(c(tfCorNumberNegative, tfCorNumberPositive))

cm<-c('#c82423','gray','#2878b5')
names(cm)<-c('up','nc','down')
tfDegDir<-epiTFsDEG[match(names(tfCorNumber),epiTFsDEG$genes),]
tfDegDir$deg[is.na(tfDegDir$deg)]<-'nc'
tfColors<-cm[match(tfDegDir$deg, names(cm))]
saveImage2("rna.deg.epiTFs.coexp.anno.pdf",width = 10,height = 2.5)
dev.off()



saveImage2("rna.deg.epiTFs.coexp.pdf",width = 10,height = 6)
par(mfrow = c(3, 1))
barplot(tfCorNumberNegative,ylab="Co-Expression Gene",las=2,cex.names =0.8, col='#9ac9db',ylim = c(0, yMax))
barplot(tfCorNumberPositive,ylab="Co-Expression Gene",las=2,cex.names =0.8, col='#f8ac8c',ylim = c(0, yMax))
barplot(rep(1,length(tfColors)),ylab="Co-Expression Gene",las=2,cex.names =0.8, col=tfColors)
dev.off()
