library(dplyr)

source('./R/base.R')


#----------------------------------------------------------------------------------------------------------------------
# Table S2. Summary of clinical characteristics associated with different subtypes of LUAD				
#----------------------------------------------------------------------------------------------------------------------

# Age
## Kruskal-Wallis H test
data <- loadData('sampleInfo')
kruskal.test(Age ~ Group, data = data)
## Statistical parameters: median and IQR
group_by(data, Group) %>% summarise(median = median(Age),IQR = IQR(Age))

# Smoking
## two-way Chi-squared test
data <- loadData('sampleInfo')
data <- filter(data, Smoking!='unknown')
car.data<-data.frame(data$Group, data$Smoking)
target <- table(car.data)
chisq.test(target)
print(target)

# Sex
## two-way Chi-squared test
data <- loadData('sampleInfo')
car.data<-data.frame(data$Group, data$Sex)
target <- table(car.data)
chisq.test(target)
print(target)

data <- loadData('sampleInfo')
data <- filter(data, Group!='L0')
car.data<-data.frame(data$Group, data$Sex)
target <- table(car.data)
chisq.test(target)
print(target)

#----------------------------------------------------------------------------------------------------------------------
# Table S5. The distribution of SR-DMCs in genome
#----------------------------------------------------------------------------------------------------------------------
SRDMC<-loadData('SRDMC')
SRDMC <- GRanges(SRDMC)
SRDMCList<-split(SRDMC, SRDMC$class)
options(ChIPseeker.ignore_1st_exon = T)
options(ChIPseeker.ignore_1st_intron = T)
options(ChIPseeker.ignore_downstream = T)
options(ChIPseeker.ignore_promoter_subcategory = T)
peakAnnoList<- lapply(SRDMCList, annotatePeak, tssRegion=c(-5000, 3000),TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene)
parseAnno<- function(gr){
  stat<-gr@annoStat
  stat$num<-gr@annoStat$Frequency*gr@peakNum/100
  stat
}
eHyperDmc<-parseAnno(peakAnnoList$`Early-Hyper-DMC`)
eHypoDmc<-parseAnno(peakAnnoList$`Early-Hypo-DMC`)
lHyperDmc<-parseAnno(peakAnnoList$`Late-Hyper-DMC`)
lHypoDmc<-parseAnno(peakAnnoList$`Late-Hypo-DMC`)
tab<-data.frame(earlyHyper=eHyperDmc$num, earlyHypo=eHypoDmc$num, lateHyper=lHyperDmc$num, lateHypo=lHypoDmc$num)
rownames(tab)<-eHyperDmc$Feature
chisq.test(tab)

#----------------------------------------------------------------------------------------------------------------------
# Table S6. The distribution of SR-DMRs in genome
#----------------------------------------------------------------------------------------------------------------------
SRDMR<-loadData('SRDMR')
SRDMR <- GRanges(SRDMR)
SRDMRList<-split(SRDMR, SRDMR$class)
options(ChIPseeker.ignore_1st_exon = T)
options(ChIPseeker.ignore_1st_intron = T)
options(ChIPseeker.ignore_downstream = T)
options(ChIPseeker.ignore_promoter_subcategory = T)
peakAnnoList<- lapply(SRDMRList, annotatePeak, tssRegion=c(-5000, 3000),TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene)
parseAnno<- function(gr){
  stat<-gr@annoStat
  stat$num<-gr@annoStat$Frequency*gr@peakNum/100
  stat
}
eHyperDmr<-parseAnno(peakAnnoList$`Early-Hyper-DMR`)
eHypoDmr<-parseAnno(peakAnnoList$`Early-Hypo-DMR`)
lHyperDmr<-parseAnno(peakAnnoList$`Late-Hyper-DMR`)
lHypoDmr<-parseAnno(peakAnnoList$`Late-Hypo-DMR`)
tab<-data.frame(earlyHyper=eHyperDmr$num, earlyHypo=eHypoDmr$num, lateHyper=lHyperDmr$num, lateHypo=lHypoDmr$num)
rownames(tab)<-eHyperDmc$Feature
chisq.test(tab)
