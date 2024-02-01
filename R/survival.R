

os<-read.csv(file.path(CONFIG$dataIntermediate,'gepia2', 'table_survival.txt'),sep = '\t')
dfs<-read.csv(file.path(CONFIG$dataIntermediate,'gepia2', 'table_df_survival.txt'),sep = '\t')
# survival<-list(
#   os=os$Gene.Symbol,
#   dfs=dfs$Gene.Symbol
# )
# saveRDS(survival,file.path(CONFIG$dataIntermediate,'gepia2', 'survival.rds'))



survival<-readRDS(file.path(CONFIG$dataIntermediate,'gepia2', 'survival.rds'))







survival<-Survival()

survival$pickSurvivalDeg(genes)



