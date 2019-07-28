setwd('D:/Promoter Design Data/FACS-Seq')
require(data.table)
#x = fread('FS9_means_just_designs.csv')
x = fread('final_validation_FACS-Seq_means_just_designs_corrected.csv')

get.min.gc = function(seq, window) {
  gcs = scan.for.gc(seq, window)
  return(min(gcs))
}

get.max.gc = function(seq, window) {
  gcs = scan.for.gc(seq, window)
  return(max(gcs))
}

scan.for.gc = function(seq, window) {
  gcs = vector(mode='numeric', length = nchar(seq) - window + 1)
  for(i in 1:(nchar(seq) - window + 1)) {
    subseq = substr(seq, i, i+window-1)
    gcs[i] = sum(strsplit(subseq, '')[[1]] %in% c('G','C'))/nchar(subseq)
  }
  return(gcs)
}

x$GC.min = sapply(x$Seqs, get.min.gc, window = 20)
x$GC.max = sapply(x$Seqs, get.max.gc, window = 20)
gc.exps = c(22,23,26,27,33,34,37,38,44,45,48,49)
x$isfilt = x$Experiment %in% gc.exps
xmin = x; xmin$val = x$GC.min; xmin$Cond = 'min'
xmax = x; xmax$val = x$GC.max; xmax$Cond = 'max'

xm = rbind(xmin, xmax)

xf = data.table(Cond = as.factor(xm$Cond), val = xm$val, Experiment = as.factor(xm$Experiment), isfilt = xm$isfilt)

xf$Cond.fac = factor(xf$Cond, levels = levels(factor(xf$Cond))[c(2,1)])

png(filename = 'D:/Promoter Design Data/Figures/Final PNGs/S18.png',
    units = 'cm', width = 18, height = 5, res = 600)

p = ggplot(data = xf, aes(x = as.factor(Experiment), y = val, color = isfilt)) + geom_boxplot() + theme_bw() + 
  facet_grid(Cond.fac ~ .) + 
  theme(plot.title = element_text(hjust = 0.5, size = 8, face='bold')) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold'),
        legend.text = element_text(size=6), legend.title = element_text(size=8, face='bold')) +
  labs(x = 'Experiment ID', y = 'GC Content', color = 'GC Constraint')
print(p)
dev.off()






