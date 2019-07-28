# Plots to show how FACS-seq works.
setwd('D:/Promoter Design Data/')
source('Code/FCS file analysis.R')
library(data.table)
library(ggplot2)
library(MASS)

fns = 'pGPD Library Panel/Data/BK2018-03-22CSY3-12_0nM_1.0001.fcs'

files.dat = lapply(lapply(fns, load.fromfilename), extract.data, vals = c('mCherry','GFP'))

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
files.dat$Strength = files.dat$GFP - files.dat$mCherry
# add a fake point to have even scales with pGPD
files.fake = data.frame(mCherry = NA, GFP = NA, Replicate = 1, Strength = 2.623037)
files.dat = rbind(files.dat, files.fake)

# add in points from an example FCS (for illustrative purposes)
fn.example = 'Final Flow Validation FCS/BK2018-12-16P2C07.0001.fcs'
dat.example = extract.data(load.fromfilename(fn.example), vals = c('mCherry','GFP'))
dat.example = data.frame(mCherry = log10(dat.example[,1]), GFP = log10(dat.example[,2]), Replicate = 1)
dat.example$Strength = dat.example$GFP - dat.example$mCherry
# limit the number of dat.example points shown
dat.example = dat.example[dat.example$Strength > -1,]
dat.example = dat.example[1:50,]

bin.edges = c(-1.39794001,-1.21043388,-1.02292776,-0.83542163,-0.6479155,-0.46040938,-0.27290325,-0.08539712,0.102109,0.28961513,0.47712125)
bin.edges = bin.edges + 0.17 

# Scatterplot of library and single-sequence cells (is P2C07 from 12-16-18 plate validation)

png(filename = 'D:/Promoter Design Data/Figures/PNGs/1B - scatter.png',
    units = 'cm', width = 4, height = 4, res = 600)
p = ggplot(files.dat, aes(x=mCherry, y=GFP, color=Strength)) + geom_point(alpha = 0.75, pch = '.', stroke = NA) + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.4, size = 8, face='bold')) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold')) + theme(legend.position = 'none') +
  scale_colour_gradient2(low = 'darkred', mid = 'gray', high = 'darkgreen') + 
  #labs(x = 'mCherry (log10)', y = 'GFP (log10)', title = 'Two-Color Fluorescence:\nExample Library') +
  xlim(0,4) + ylim(-1,5) + 
  geom_abline(slope = 1, intercept = bin.edges, lty = 2, lwd = 0.5, alpha = 0.5) + 
  geom_point(data = dat.example, size = 1, color = 'blue', fill = 'white', stroke = 0.05) + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  theme(axis.text.x=element_blank(), axis.text.y=element_blank(),
        axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) +
  theme(axis.text = element_text(size=10))
print(p)
dev.off()

scale.fac = 100
dat.example.scaled = data.frame(mCherry = rep(dat.example$mCherry, scale.fac), GFP = rep(dat.example$GFP, scale.fac),
                                Replicate = rep(dat.example$Replicate, scale.fac), Strength = rep(dat.example$Strength, scale.fac))

scale.fac.gauss = 1000
points.test = seq(from = -2, to = 2, by = 0.001)
fits = fitdistr(dat.example$Strength, 'normal')$estimate
points.test = data.frame(x = points.test, y = scale.fac.gauss*dnorm(points.test, mean = fits[1], sd = fits[2]))

breaks = c(-1.5,bin.edges,0.9)
cols = colorRampPalette(c('darkred','gray','darkgreen'))(12)


png(filename = 'D:/Promoter Design Data/Figures/PNGs/1B - histogram.png',
    units = 'cm', width = 4, height = 4, res = 600)

p = ggplot(files.dat, aes(x = Strength)) + theme_bw() +
  stat_bin(breaks = breaks, fill = cols, alpha = 0.5) +
  stat_bin(data = dat.example.scaled, breaks = breaks, fill = 'blue') + #cols) + 
  geom_line(data = points.test, aes(x = x, y = y)) + 
  geom_vline(xintercept = fits[1], lwd = 1) + 
  geom_vline(xintercept = bin.edges, lty = 2) +
  xlim(-1.5, 1) +
  theme(axis.title.x=element_blank()) + 
  theme(axis.title.y=element_blank()) +
  theme(axis.text.x=element_blank(), axis.text.y=element_blank(),
        axis.ticks.x=element_blank(), axis.ticks.y=element_blank())
print(p)
dev.off()










