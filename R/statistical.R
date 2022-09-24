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
# Table S4. Hyper-DMCs and Hypo-DMCs of different group in LUAD				
#----------------------------------------------------------------------------------------------------------------------
dmcL0vsL1<-loadDataBed('dmcL0vsL1')
dmcL0vsL2<-loadDataBed('dmcL0vsL2')
dmcL0vsL3<-loadDataBed('dmcL0vsL3')
sum(table(dmcL0vsL1$class))
sum(table(dmcL0vsL2$class))
sum(table(dmcL0vsL3$class))
L1.vs.L0<-data.frame(group="L1.vs.L0", class=dmcL0vsL1$class)
L2.vs.L0<-data.frame(group="L2.vs.L0", class=dmcL0vsL2$class)
L3.vs.L0<-data.frame(group="L3.vs.L0", class=dmcL0vsL3$class)
data<-Reduce(rbind, list(L1.vs.L0,L2.vs.L0,L3.vs.L0))
target<-table(data)
chisq.test(target)
print(target)

