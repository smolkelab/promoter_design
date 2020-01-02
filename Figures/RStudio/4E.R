setwd('D:/Datasets/20190110 VYB (pBK71 prelim)/')
setwd('D:/Promoter Design Data/Final Flow Validation FCS/Baseline')
fns = dir()
meds = sapply(fns, get.one.median)
SEM.scale = 1/sqrt(3)
meds.min = 10^(mean(meds)-sd(meds)*SEM.scale)
meds.max = 10^(mean(meds)+sd(meds)*SEM.scale)
allmeds = 10^meds
meds = 10^mean(meds)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

dat.zev.ar = dat[(is.na(dat$Seqs) & dat$Name %in% c('ZEV_Pr3', 'ZEV_Pr4', 'ZEV_Pr8')) | 
                   dat$Experiment.designed == 44,]# | dat$ZEV_AR_offscale_cand | dat$ZEV_I_AR_cand,]
dat.zev.ar$Name = factor(dat.zev.ar$Name, levels = dat.zev.ar$Name[order(dat.zev.ar$Experiment.designed, dat.zev.ar$VYB_AR)])
dat.zev.ar$`Design Type` = dat.zev.ar$Experiment.designed
dat.zev.ar$`Design Type`[is.na(dat.zev.ar$`Design Type`)] = 'Control'
dat.zev.ar$`Design Type`[dat.zev.ar$Experiment.designed == 44] = 'Evolution-GC'
dat.zev.ar$`Design Type` = factor(dat.zev.ar$`Design Type`)

dat.zev.bars = data.frame(Name=dat.zev.ar$Name, VYB_Mean_A=dat.zev.ar$VYB_Mean_A,
                          `Design Type`=dat.zev.ar$`Design Type`, bar.min.a=dat.zev.ar$bar.min.a,
                          bar.max.a=dat.zev.ar$bar.max.a)

dat.zev.bars.m = melt.for.jitter(dat.zev.ar, 'Uninduced')
neg = data.frame(Name = 'Background', VYB_Mean_A=mean(meds), Design.Type='Background',
                 bar.min.a=meds.min, bar.max.a=meds.max)

neg.m = data.frame(Name = 'Background', Design.Type='Background', variable = c('R1', 'R2', 'R3'),
                   value = allmeds)
#dat.zev.bars = rbind(neg, dat.zev.bars)
dat.zev.bars = rbind(dat.zev.bars[dat.zev.bars$Design.Type == 'Evolution-GC',], neg, dat.zev.bars[dat.zev.bars$Design.Type == 'Control',])

dat.zev.bars$Design.Type = factor(dat.zev.bars$Design.Type, 
                                  levels = levels(factor(dat.zev.bars$Design.Type))[c(2,3,1)] )

setnames(dat.zev.bars.m, 'Design Type', 'Design.Type')
#dat.zev.bars.m = rbind(neg.m, dat.zev.bars.m)
dat.zev.bars.m = rbind(dat.zev.bars.m[dat.zev.bars.m$Design.Type == 'Evolution-GC',], neg.m, dat.zev.bars.m[dat.zev.bars.m$Design.Type == 'Control',])

dat.zev.bars$Name = factor(dat.zev.bars$Name, levels = levels(dat.zev.bars$Name)[c(1:11,14,12:13)])
#dat.zev.bars.m$Name = factor(dat.zev.bars.m$Name, levels = dat.zev.bars.m$Name[c(1:11,14,12:13)])

dat.zev.bars = dat.zev.bars[dat.zev.bars$Name != 'ZEV_Pr8',]
dat.zev.bars.m = dat.zev.bars.m[dat.zev.bars.m$Name != 'ZEV_Pr8',]
xstr = c(paste0('pZEV-AR-', 1:10), 'Background', 'P4', 'P3')

png(filename = 'D:/Promoter Design Data/Figures/PNGs/4E.png',
    units = 'cm', width = w, height = h, res = 600)
p = ggplot(dat.zev.bars, aes(x = Name, y = VYB_Mean_A, color = Design.Type, fill = Design.Type,
                           ymin = bar.min.a, ymax = bar.max.a)) + 
  theme_bw() +  geom_col(alpha = 0.5, width = 0.6) + geom_errorbar(width = 0.6) +
  theme(plot.title = element_text(hjust = 0.4, size = 8, face='bold')) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold'),
        legend.text = element_text(size=6), legend.title = element_text(size=8, face='bold'),
        legend.key.size = unit(0.2,'cm')) +
  theme(axis.text.x=element_text(angle=270, hjust = 0),
        axis.ticks.x=element_blank()) +
  labs(x = '', y='Uninduced promoter activity \npZEV-AR designs', color = 'Design type', fill = 'Design type') + 
  scale_x_discrete(labels=xstr) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  scale_fill_manual(values=cols) + 
  scale_color_manual(values=cols)

q = p + geom_jitter(data=dat.zev.bars.m, aes(x=Name, y=value, color = Design.Type,
                                           ymin=NULL, ymax=NULL), width = 0.1, height = 0, size = 0.5)

print(q)

dev.off()

##### Is it expected to have values "below background"?
nx = names(table(as.character(dat.zev.bars.m$Name)))
for(i in 2:length(nx)) {
  v = dat.zev.bars.m$value
  n = as.character(dat.zev.bars.m$Name)
  vals = v[n == nx[i]]
  print(nx[i])
  print(vals)
  print(t.test(v[n == 'Background'], vals, 'greater')$p.value)
}





