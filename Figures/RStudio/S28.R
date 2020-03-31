library(data.table)
library(ggplot2)

setwd('D:/Promoter Design Data/Clustering Test/')
fn = 'intergroup_alignment_test_GPD.csv'


table.from.fns = function(fns) {
  tables = list()
  for(i in 1:length(fns)) {
    f = fns[i]
    scores = sapply(readLines(f), function(x) strsplit(x,',')[[1]][2] )
    scores = as.numeric(scores)
    t = data.table(scores = scores, fn = f)
    #print(t)
    tables[[i]] = t
  }
  return(rbindlist(tables))
}

gpd = table.from.fns(c('internal_alignment_test_GPD.csv', 'intergroup_alignment_test_GPD.csv'))
gpd$Library = 'pGPD'
zev = table.from.fns(c('internal_alignment_test_ZEV.csv', 'intergroup_alignment_test_ZEV.csv'))
zev$Library = 'pZEV'

seqs = rbind(gpd, zev)
seqs$fn = sapply(seqs$fn, function(x) strsplit(x, '_')[[1]][1])
setnames(seqs, c('scores','fn'), c('Scores','Comparison'))

png(filename = 'D:/Promoter Design Data/Figures/PNGs/S28.png',
    units = 'in', res = 144, width = 7, height = 5.25)
ggplot(seqs, aes(x = Scores, col = Comparison, fill = Comparison)) + geom_bar() + facet_grid(Library ~ .) + 
  theme_bw() + labs(x = 'Alignment Score', y = 'Number of Read Clusters')
dev.off()

