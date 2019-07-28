# Visualize pGPD dataset
library(data.table)
library(ggplot2)
library(viridis)

setwd('D:/Promoter Design Data/FACS-Seq/')
dat.tra = fread('means_trainable_GPD.csv')

print(range(dat.tra$Strength))
print(mean(dat.tra$Strength))
print(sd(dat.tra$Strength))

print(median(dat.tra$Strength))
print(IQR(dat.tra$Strength))

nbins = 40
cols = colorRampPalette(c('darkred','gray','darkgreen'))(nbins)

png(filename = 'D:/Promoter Design Data/Figures/PNGs/1C.png',
    units = 'cm', width = 8, height = 6.4, res = 600)

p = ggplot(dat.tra, aes(x = Strength)) + 
  stat_bin(bins = nbins, fill = cols) +
  theme_bw() +
  theme(axis.title.x=element_blank()) + 
  theme(axis.title.y=element_blank()) +
  theme(axis.text.x = element_text(size=8), axis.text.y = element_text(size=8))
print(p)
dev.off()


