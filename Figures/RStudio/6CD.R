library(data.table)
library(ggplot2)
library(UpSetR)

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

setwd('D:/Promoter Design Data/Motif Search/')
dirs.in = c('AR_1','AR_2','AR_3')
fns = paste0('6C_scores_', dirs.in, '.png')
ys = list()
for(i in 1:length(dirs.in)) {
  y = load.csvs(dirs.in[i])
  y$Table.label = 'Screening'
  y$Table.label[y$Table.ID == 27] = 'Evolution-GC'
  y$Motif.simple = substr(y$Motif,1,4) == 'GCTA'
  y$Motif.ID = paste0('ZEV ', i, collapse = '')
  ys[[i]] = y
}
y = rbindlist(ys)
y$Table.label = factor(y$Table.label)
y$Table.label = factor(y$Table.label, levels = rev(levels(y$Table.label)))

setwd('D:/Promoter Design Data/Figures/PNGs')
png(filename = '6C_scores.png', units = 'cm', width = 9, height = 7, res = 600)

p = ggplot(data = y, aes(x = Mut_score, color = Motif.simple, fill = Motif.simple, pch = Table.label)) + 
  stat_bin(binwidth = 0.01, alpha = 0.7) + theme_bw() +
  facet_grid(Table.label ~ Motif.ID) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold'),
        legend.text = element_text(size=6), legend.title = element_text(size=8, face='bold'),
        legend.key.size = unit(0.2,'cm')) +
  labs(x = 'Median Score Differential', y = 'Motif Counts', color = 'GCTA', fill = 'GCTA')
print(p)
dev.off()

# Couldn't get these plots to size nicely automatically - use the RStudio dialog instead :(

#png(filename = '6D_upset_screening.png',
#        units = 'cm', width = 6, height = 12, res = 600)
motifs = list(ZEV.1 = which(ys[[1]][ys[[1]]$Table.ID == 23,]$Motif.simple), 
              ZEV.2 = which(ys[[2]][ys[[2]]$Table.ID == 23,]$Motif.simple), 
              ZEV.3 = which(ys[[3]][ys[[3]]$Table.ID == 23,]$Motif.simple))

upset(fromList(motifs), order.by = "degree", nintersects = NA, empty.intersections = TRUE,
      mainbar.y.label = '', 
      sets = c('ZEV.3', 'ZEV.2', 'ZEV.1'), keep.order = TRUE, text.scale = 1.5)
#dev.off()

#png(filename = '6D_upset_evolution-GC.png',
#    units = 'cm', width = 7, height = 7, res = 600)
motifs = list(ZEV.1 = which(ys[[1]][ys[[1]]$Table.ID == 27,]$Motif.simple), 
              ZEV.2 = which(ys[[2]][ys[[2]]$Table.ID == 27,]$Motif.simple), 
              ZEV.3 = which(ys[[3]][ys[[3]]$Table.ID == 27,]$Motif.simple))


upset(fromList(motifs), order.by = "degree", nintersects = NA, empty.intersections = TRUE,
      mainbar.y.label = '', 
      sets = c('ZEV.3', 'ZEV.2', 'ZEV.1'), keep.order = TRUE, text.scale = 1.5)
#dev.off()

### For text ###
# Median differentials
# non-GCTA
# Screening
# ZEV 1
median(y[y$Motif.ID == 'ZEV 1' & y$Table.label == 'Screening' & !(y$Motif.simple),]$Mut_score)
# ZEV 2
median(y[y$Motif.ID == 'ZEV 2' & y$Table.label == 'Screening' & !(y$Motif.simple),]$Mut_score)
# ZEV 3
median(y[y$Motif.ID == 'ZEV 3' & y$Table.label == 'Screening' & !(y$Motif.simple),]$Mut_score)
# Evolution-GC
# ZEV 1
median(y[y$Motif.ID == 'ZEV 1' & y$Table.label == 'Evolution-GC' & !(y$Motif.simple),]$Mut_score)
# ZEV 2
median(y[y$Motif.ID == 'ZEV 2' & y$Table.label == 'Evolution-GC' & !(y$Motif.simple),]$Mut_score)
# ZEV 3
median(y[y$Motif.ID == 'ZEV 3' & y$Table.label == 'Evolution-GC' & !(y$Motif.simple),]$Mut_score)

# GCTA
# Screening
# ZEV 1
median(y[y$Motif.ID == 'ZEV 1' & y$Table.label == 'Screening' & (y$Motif.simple),]$Mut_score)
# ZEV 2
median(y[y$Motif.ID == 'ZEV 2' & y$Table.label == 'Screening' & (y$Motif.simple),]$Mut_score)
# ZEV 3
median(y[y$Motif.ID == 'ZEV 3' & y$Table.label == 'Screening' & (y$Motif.simple),]$Mut_score)
# Evolution-GC
# ZEV 1
median(y[y$Motif.ID == 'ZEV 1' & y$Table.label == 'Evolution-GC' & (y$Motif.simple),]$Mut_score)
# ZEV 2
median(y[y$Motif.ID == 'ZEV 2' & y$Table.label == 'Evolution-GC' & (y$Motif.simple),]$Mut_score)
# ZEV 3
median(y[y$Motif.ID == 'ZEV 3' & y$Table.label == 'Evolution-GC' & (y$Motif.simple),]$Mut_score)

# Wilcox tests: GCTA vs. not
wilcox.test(y[y$Motif.ID == 'ZEV 1' & y$Table.label == 'Screening' & !(y$Motif.simple),]$Mut_score,
            y[y$Motif.ID == 'ZEV 1' & y$Table.label == 'Screening' & (y$Motif.simple),]$Mut_score, 'greater')
wilcox.test(y[y$Motif.ID == 'ZEV 2' & y$Table.label == 'Screening' & !(y$Motif.simple),]$Mut_score,
            y[y$Motif.ID == 'ZEV 2' & y$Table.label == 'Screening' & (y$Motif.simple),]$Mut_score, 'greater')
wilcox.test(y[y$Motif.ID == 'ZEV 3' & y$Table.label == 'Screening' & !(y$Motif.simple),]$Mut_score,
            y[y$Motif.ID == 'ZEV 3' & y$Table.label == 'Screening' & (y$Motif.simple),]$Mut_score, 'greater')

# Evolution-GC
wilcox.test(y[y$Motif.ID == 'ZEV 1' & y$Table.label == 'Evolution-GC' & !(y$Motif.simple),]$Mut_score,
            y[y$Motif.ID == 'ZEV 1' & y$Table.label == 'Evolution-GC' & (y$Motif.simple),]$Mut_score, 'greater')
wilcox.test(y[y$Motif.ID == 'ZEV 2' & y$Table.label == 'Evolution-GC' & !(y$Motif.simple),]$Mut_score,
            y[y$Motif.ID == 'ZEV 2' & y$Table.label == 'Evolution-GC' & (y$Motif.simple),]$Mut_score, 'greater')
wilcox.test(y[y$Motif.ID == 'ZEV 3' & y$Table.label == 'Evolution-GC' & !(y$Motif.simple),]$Mut_score,
            y[y$Motif.ID == 'ZEV 3' & y$Table.label == 'Evolution-GC' & (y$Motif.simple),]$Mut_score, 'greater')

# 6D

motifs.sc = data.table(ZEV.1 = ys[[1]][ys[[1]]$Table.ID == 23,]$Motif.simple, 
              ZEV.2 = ys[[2]][ys[[2]]$Table.ID == 23,]$Motif.simple, 
              ZEV.3 = ys[[3]][ys[[3]]$Table.ID == 23,]$Motif.simple,
              score.1 = ys[[1]][ys[[1]]$Table.ID == 23,]$Mut_score, 
              score.2 = ys[[2]][ys[[2]]$Table.ID == 23,]$Mut_score, 
              score.3 = ys[[3]][ys[[3]]$Table.ID == 23,]$Mut_score)

motifs.eg = data.table(ZEV.1 = ys[[1]][ys[[1]]$Table.ID == 27,]$Motif.simple, 
              ZEV.2 = ys[[2]][ys[[2]]$Table.ID == 27,]$Motif.simple, 
              ZEV.3 = ys[[3]][ys[[3]]$Table.ID == 27,]$Motif.simple,
              score.1 = ys[[1]][ys[[1]]$Table.ID == 27,]$Mut_score, 
              score.2 = ys[[2]][ys[[2]]$Table.ID == 27,]$Mut_score, 
              score.3 = ys[[3]][ys[[3]]$Table.ID == 27,]$Mut_score)


table(apply(motifs.sc,1,sum))
nrow(motifs.sc)
table(apply(motifs.eg,1,sum))
nrow(motifs.eg)

