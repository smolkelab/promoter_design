# Test whether sequence designs are actually sequence-diverse
# used script facs-seq/mean_extraction/design_diversity_testing.sh on instance 3
# to generate table with fields Read_1,Read_2,Score,Filename
# For each filename, start by getting a boxplot of the 'Score' (score is higher if alignment better)

library(data.table)
library(ggplot2)
setwd('D:/Promoter Design Data/Diversity Test/')

fn = 'diversity_test.csv'
fn.key = 'designs_filename_map.csv'
x = fread(fn)
x$Key = ''
x$Promoter = ''
x$Design = ''
x$Target = ''
x$GC_Filter = ''
x$Merge = ''
#x.key = fread(fn.key, colClasses = list(character=2))
x.key = fread(fn.key, colClasses = 'character')
for(i in 1:nrow(x.key)) {
  rows.to.change = which(x$Filename ==  x.key$Filename[i])
  #x[rows.to.change,]$Key = x.key$ID_merged[i]
  x[rows.to.change,]$Key = x.key$ID_final[i]
  x[rows.to.change,]$Promoter = x.key$Promoter[i]
  x[rows.to.change,]$Design = x.key$Design[i]
  x[rows.to.change,]$Target = x.key$Target[i]
  x[rows.to.change,]$GC_Filter = x.key$GC_Filter[i]
  x[rows.to.change,]$Merge = x.key$Merge[i]
}

# With padding, all these sequences are 363 bases long - subtract parts of the alignments
# that correspond to this padding. For GPD, subtract 50; for ZEV, subtract 116
x$Score.adj = x$Score - 50 - 66*(x$Promoter == 'ZEV')

png(filename = 'D:/Promoter Design Data/Figures/PNGs/S19_GPD.png',
    units = 'cm', width = 10, height = 5, res = 600)

y = x[x$Promoter == 'GPD']
p = ggplot(data = y, aes(x = Key, y = Score.adj, fill = Design)) + theme_bw() + geom_boxplot(outlier.size = 0.5) + 
  theme(plot.title = element_text(hjust = 0.5, size = 8, face='bold')) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold')) +
  theme(legend.text = element_text(size=6), legend.title = element_text(size=8, face='bold')) +
  labs(x = 'Experiment ID', y = 'Alignment Score\nGPD')
print(p)
dev.off()

png(filename = 'D:/Promoter Design Data/Figures/PNGs/S19_ZEV-I.png',
    units = 'cm', width = 10, height = 5, res = 600)

y = x[x$Promoter == 'ZEV' & x$Target %in% c('Control', 'Induced')]
p = ggplot(data = y, aes(x = Key, y = Score.adj, fill = Design)) + theme_bw() + geom_boxplot(outlier.size = 0.5) + 
  theme(plot.title = element_text(hjust = 0.5, size = 8, face='bold')) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold')) +
  theme(legend.text = element_text(size=6), legend.title = element_text(size=8, face='bold')) +
  labs(x = 'Experiment ID', y = 'Alignment Score\nZEV-Induced')
print(p)
dev.off()

png(filename = 'D:/Promoter Design Data/Figures/PNGs/S19_ZEV-AR.png',
    units = 'cm', width = 10, height = 5, res = 600)

y = x[x$Promoter == 'ZEV' & x$Target %in% c('Control', 'AR')]
p = ggplot(data = y, aes(x = Key, y = Score.adj, fill = Design)) + theme_bw() + geom_boxplot(outlier.size = 0.5) + 
  theme(plot.title = element_text(hjust = 0.5, size = 8, face='bold')) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold')) +
  theme(legend.text = element_text(size=6), legend.title = element_text(size=8, face='bold')) +
  labs(x = 'Experiment ID', y = 'Alignment Score\nZEV-Activation Ratio')
print(p)
dev.off()

# tests for alignment distance between groups
# Depending on the p-value, just about any pairing is "significant"

gps = sort(unique(x$Key))

# pval.mat format: each cell is pval for test "row < column"
pval.mat = matrix(nrow = length(gps), ncol = length(gps))
rownames(pval.mat) = gps; colnames(pval.mat) = gps
for(i in 1:nrow(pval.mat)) { for(j in 1:ncol(pval.mat)) {
  #p = (t.test(x$Score.adj[x$Key == gps[i]], x$Score.adj[x$Key == gps[j]], alternative='less'))$p.value
  p = (wilcox.test(x$Score.adj[x$Key == gps[i]], x$Score.adj[x$Key == gps[j]], alternative='two.sided'))$p.value
  pval.mat[i,j] = p #min(-log10(p), 10)
}}

# pGPD: control vs. screen
max(apply(pval.mat[19:20,1:6],1,max)) # 0.01629639
# control vs. evolution
max(apply(pval.mat[21:24,1:6],1,max)) # 5.047073e-06
# gradient vs. control
max(apply(pval.mat[1:6,25:29],1,max)) # 4.083107e-106

#pZEV: control vs. screen
max(apply(pval.mat[30:31,7:18],1,max)) # 0.2912263
# control vs. evolution
max(apply(pval.mat[32:35,7:18],1,max)) # 0.006133204
# gradient vs. control
max(apply(pval.mat[7:18, 36:40],1,max)) # 1.151119e-73

#pZEV-AR: control vs. screen
max(apply(pval.mat[41:42,7:18],1,max)) # 0.06252029
# control vs. evolution
max(apply(pval.mat[43:46,7:18],1,max)) # 0.0006066611
# gradient vs. control
max(apply(pval.mat[7:18, 47:51],1,max)) # 6.860136e-50
