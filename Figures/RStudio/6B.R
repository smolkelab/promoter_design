library(data.table)
library(ggplot2)

setwd('D:/Promoter Design Data/Motif Search/')

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

dir.in = 'TATATA'
x = load.csvs(dir.in)
y = x[x$Table.ID %in% c(12,16),]
y$Table.label = 'Screening'
y$Table.label[y$Table.ID == 16] = 'Evolution-GC'
y$Motif.simple = y$Motif
y$Motif.simple[!(y$Motif %in% c('TATATA','TTATAT','TTTATA'))] = 'Other'
y$Table.label = factor(y$Table.label)
y$Table.label = factor(y$Table.label, levels = rev(levels(y$Table.label)))

setwd('D:/Promoter Design Data/Figures/PNGs')

png(filename = '6B_scores.png',
    units = 'cm', width = 7.7, height = 5, res = 600)

p = ggplot(data = y, aes(x = Mut_score, color = Motif.simple, fill = Motif.simple, pch = Table.label)) + 
  stat_bin(binwidth = 0.01, alpha = 0.7) + theme_bw() + facet_wrap( ~ Table.label) +

  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold'),
        legend.text = element_text(size=6), legend.title = element_text(size=8, face='bold'),
        legend.key.size = unit(0.2,'cm')) +
  labs(x = 'Median Score Differential', y = 'Motif Counts', color = 'Motif', fill = 'Motif')
print(p)
dev.off()


# True sequence starts after 58 bp of padding for pZEV designs
png(filename = '6B_pos.png',
    units = 'cm', width = 7.7, height = 3, res = 600)
p = ggplot(data = y, aes(x = Pos - 58, fill = Table.label, color = Table.label)) + theme_bw() + 
  geom_density(alpha = 0.5) + #stat_bin(binwidth = 1, position = 'dodge') +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold'),
        legend.text = element_text(size=6), legend.title = element_text(size=8, face='bold'),
        legend.key.size = unit(0.2,'cm')) +
  labs(x = 'Sequence Position', y = 'Density', color = 'Strategy', fill = 'Strategy')
print(p)
dev.off()

### For text ###
# Median drop due to mutation
a = median(y$Mut_score[y$Table.label == 'Screening'])
b = median(y$Mut_score[y$Table.label == 'Evolution-GC'])
print(1 - 10^a)
print(1 - 10^b)




