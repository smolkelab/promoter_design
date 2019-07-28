library(data.table)
library(ggplot2)
setwd('D:/Promoter Design Data/Read Group Fates/')

max.pre = 20

process.fates = function(fn, max.pre) {
x = fread(fn)
# Possible fates: Singleton, N, Ambiguous, Aligned - N never appears though
fates = vector(mode = 'character', length = nrow(x))
for(i in 1:nrow(x)) {
  pre = x$Pre[i]; ambig = x$Ambig[i]; n = x$N[i]; aln = x$Aln[i]
  if(pre == 1) { tmp = 'Singleton' }
  if(pre == ambig) { tmp = 'Ambiguous' }
  if(pre == n) { tmp = 'N' }
  if(pre == aln) { tmp = 'Aligned' }
  fates[i] = tmp
}
x$Fate = fates
x$Pre.adj = sapply(x$Pre, function(x) min(x,max.pre))

fate.mat = matrix(nrow = max.pre, ncol = 4)
colnames(fate.mat) = c('Singleton', 'N', 'Ambiguous', 'Aligned')
for(i in 1:nrow(fate.mat)) {
  for(j in 1:ncol(fate.mat)) {
    fate.mat[i,j] = sum(x$Pre.adj == i & x$Fate == colnames(fate.mat)[j])
}}
return(t(fate.mat))
}

fates.mat = lapply(dir(), process.fates, max.pre = max.pre)
names(fates.mat) = c('GPD', 'ZEV')
fates.mat[[1]] = melt(fates.mat[[1]]); fates.mat[[1]]$Promoter = 'GPD'
fates.mat[[2]] = melt(fates.mat[[2]]); fates.mat[[2]]$Promoter = 'ZEV'
fates = rbind(fates.mat[[1]], fates.mat[[2]])

png(filename = 'D:/Promoter Design Data/Figures/Final PNGs/S26.png',
    units = 'cm', res = 600, width = 16, height = 8)

p = ggplot(data = fates, aes(x = Var2, fill = Var1, weight= value)) + geom_bar() + facet_grid(.~Promoter) + 
  theme_bw() + 
  labs(x = 'Cluster Size', y = 'Number of Clusters', fill = 'Read Fate') + 
  theme(plot.title = element_text(hjust = 0.4, size = 12, face='bold')) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=12, face='bold'),
        legend.text = element_text(size=10), legend.title = element_text(size=12, face='bold'),
        legend.key.size = unit(0.5,'cm')) 
print(p)
dev.off()








