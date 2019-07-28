library(data.table)
library(ggplot2)
library(gplots)
library(beanplot)

setwd('D:/Promoter Design Data')
source('Code/FCS file analysis.R')

fns = 'pGPD Library Panel/Data/BK2018-03-22CSY3-12_0nM_1.0001.fcs'

files.dat = lapply(fns, load.fromfilename)
files.dat = lapply(files.dat, extract.data, vals = c('mCherry','GFP'))

# stack by replicate
mC = list()
gfp = list()
replicate = list()
for(i in 1:length(files.dat)) { 
  replicate[[i]] = rep(i, nrow(files.dat[[i]]))
  mC[[i]] = log10(files.dat[[i]][,1])
  gfp[[i]] = log10(files.dat[[i]][,2])
}

replicate = as.factor(unlist(replicate))
files.dat = data.frame(mCherry = unlist(mC), GFP = unlist(gfp))
files.dat$Replicate = replicate
files.dat$Activity = files.dat$GFP - files.dat$mCherry

png(filename = 'D:/Promoter Design Data/Figures/PNGs/S3 scatter - raw.png',
    units = 'cm', width = 6, height = 4, res = 600)
p = ggplot(files.dat, aes(x=mCherry, y=GFP, color=Activity)) + geom_point(alpha = 0.75, pch = '.', stroke = NA) + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.4, size = 8, face='bold')) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold'),
        legend.text = element_text(size=4), legend.title = element_text(size=6, face='bold'),
        legend.key.size = unit(0.2,'cm')) +
  scale_colour_gradient2(low = 'darkred', mid = 'gray', high = 'darkgreen') + 
  labs(x = 'mCherry (log10)', y = 'GFP (log10)', title = 'Two-Color Fluorescence:\nExample Library') +
  geom_abline(slope = 1, intercept = seq(from = -3, to = 3, by = 1), lty = 2, lwd = 0.5, alpha = 0.5)
print(p)
dev.off()

png(filename = 'D:/Promoter Design Data/Figures/PNGs/S3 bean - raw.png',
    units = 'cm', width = 2, height = 4, res = 600)

p = ggplot(files.dat, aes(x=Replicate, y=Activity)) + geom_violin(color='darkgray', fill = 'darkgray') + theme_bw() + 
  theme(axis.text = element_text(size=8), axis.title = element_text(size=10, face='bold'), legend.position="none",
        axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  geom_hline(yintercept = seq(from = -3, to = 3, by = 1), lty = 2, lwd = 0.5, alpha = 0.5)
print(p)
dev.off()




