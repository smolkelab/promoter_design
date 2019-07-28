library(data.table)
library(ggplot2)

load.csvs = function(dir.in) {
  fns = dir(dir.in)
  fns.load = paste0(dir.in, '/', fns)
  tags = sapply(fns, function(x) strsplit(x, '_')[[1]][1])
  loaded.tables = list()
  for(i in 1:length(fns.load)) {
    x = fread(fns.load[i])
    x$Table.ID = tags[i]
    loaded.tables[[i]] = x
  }
  return(rbindlist(loaded.tables))
}
dirs.in = c('AR_1','AR_2','AR_3')

ys = list()
setwd('D:/Promoter Design Data/Motif Search/')
for(i in 1:length(dirs.in)) {
  y = load.csvs(dirs.in[i])
  y$Table.label = 'Screening'
  y$Table.label[y$Table.ID == 27] = 'Evolution-GC'
  y$Motif.simple = substr(y$Motif,1,4) == 'GCTA'
  y$Motif.ID = paste0('ZEV ', i, collapse = '')
  ys[[i]] = y
}

motifs.eg = data.table(ZEV.1 = ys[[1]][ys[[1]]$Table.ID == 27,]$Motif.simple, 
                       ZEV.2 = ys[[2]][ys[[2]]$Table.ID == 27,]$Motif.simple, 
                       ZEV.3 = ys[[3]][ys[[3]]$Table.ID == 27,]$Motif.simple,
                       score.1 = ys[[1]][ys[[1]]$Table.ID == 27,]$Mut_score, 
                       score.2 = ys[[2]][ys[[2]]$Table.ID == 27,]$Mut_score, 
                       score.3 = ys[[3]][ys[[3]]$Table.ID == 27,]$Mut_score)


# Effect of motif redundancy for motifs in EG designs
motifs.eg$is.unique = apply(cbind(motifs.eg$ZEV.1, motifs.eg$ZEV.2, motifs.eg$ZEV.3), 1, sum) == 1

# Effect of motif redundancy for motifs in EG designs
motifs.eg$is.unique = apply(cbind(motifs.eg$ZEV.1, motifs.eg$ZEV.2, motifs.eg$ZEV.3), 1, sum) == 1

u.1 = data.table(score = motifs.eg$score.1[motifs.eg$ZEV.1], unique = motifs.eg$is.unique[motifs.eg$ZEV.1], pos = 1)
u.2 = data.table(score = motifs.eg$score.2[motifs.eg$ZEV.2], unique = motifs.eg$is.unique[motifs.eg$ZEV.2], pos = 2)
u.3 = data.table(score = motifs.eg$score.3[motifs.eg$ZEV.3], unique = motifs.eg$is.unique[motifs.eg$ZEV.3], pos = 3)
uniques = rbind(u.1, u.2, u.3)
uniques$pos = paste0('ZEV ', uniques$pos)

wilcox.test(uniques$score[uniques$unique], uniques$score[!(uniques$unique)], 'less') # p = 1.05e-07

png(filename = 'D:/Promoter Design Data/Figures/Final PNGs/S22.png',
    units = 'cm', width = 13.5, height = 10.5, res = 600)
#boxplot(score ~ unique + pos, data = uniques)
p = ggplot(data = uniques, aes(x = unique, y = score)) + geom_boxplot() + facet_grid(. ~ pos)  +
  theme_bw() + labs(x = 'Site Is Unique', y = 'Median Score Differential')
print(p)
dev.off()



