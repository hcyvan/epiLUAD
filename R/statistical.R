library(dplyr)

source('./R/base.R')


#----------------------------------------------------------------------------------------------------------------------
# Table S2. Summary of clinical characteristics associated with different subtypes of LUAD				
#----------------------------------------------------------------------------------------------------------------------

# Age
## Kruskal-Wallis H test
data <- loadData('sampleInfo', './tmp_data')
kruskal.test(Age ~ Group, data = data)
## Statistical parameters: median and IQR
group_by(data, Group) %>% summarise(median = median(Age),IQR = IQR(Age))

# Smoking
## two-way Chi-squared test
data <- loadData('sampleInfo', './tmp_data')
data <- filter(data, Smoking!='unknown')
car.data<-data.frame(data$Group, data$Smoking)
target <- table(car.data)
chisq.test(target)
print(target)

# Sex
## two-way Chi-squared test
data <- loadData('sampleInfo', './tmp_data')
car.data<-data.frame(data$Group, data$Sex)
target <- table(car.data)
chisq.test(target)
print(target)

data <- loadData('sampleInfo', './tmp_data')
data <- filter(data, Group!='L0')
car.data<-data.frame(data$Group, data$Sex)
target <- table(car.data)
chisq.test(target)
print(target)
