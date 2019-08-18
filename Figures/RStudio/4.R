# Do setup - select control sequences,
# run analysis for single-sequence validation,
# use those results to align validation FACS-Seq data to single-sequence strengths.

library(data.table)
library(ggplot2)
setwd('D:/Promoter Design Data/FACS-Seq')
fn = 'final_means_ids_added.csv'
x = fread(fn)
x$AR = x$Means_B - x$Means_A
x$Means_Avg = apply(cbind(x$Means_A, x$Means_B),1,mean)
x = x[x$Experiment != 'None',]
write.csv(x, file = 'final_validation_FACS-Seq_means_just_designs.csv')

setwd('D:/Promoter Design Data/Code')
# identify control sequences
source('Select_control_seqs.R')

# select sequences for single-sequence validation
setwd('D:/Promoter Design Data/Code')
source('Merging FS9 data and design predictions.R')
setwd('D:/Promoter Design Data/Code')
source('Sequence selection for FS9 validation.R')

# Get single-sequence validation results
setwd('D:/Promoter Design Data/Code')
source('FS9 validation VYB analysis - get medians, attach metadata.R')

# Use single-sequence validation results to correct validation FACS-Seq means
setwd('D:/Promoter Design Data/Code')
source('Correct FS9 data with single sequence data.R')

# Validation FACS-seq results: selected; best sets only
setwd('D:/Promoter Design Data/FACS-Seq/')
dat.key = fread('Fig 4 key.csv')

gpd.strong.ex = c(1,5,19,21,23,25,27) 
zev.induced.ex = c(9, 16, 30, 32, 34, 36, 38, 39)
zev.ar.ex = c(10, 17, 41,43,45,47,49,50)
zev.2D.ex = zev.ar.ex #c(9,10,16,17,30, 32, 34, 36, 38, 39, 41,43,45,47,49,50)

dat.mer.use = fread('D:/Promoter Design Data/FACS-Seq/final_validation_FACS-Seq_means_just_designs_corrected.csv')

offscale.val = 25
dat.mer.use = dat.mer.use[dat.mer.use$Experiment %in% gpd.strong.ex,]
dat.mer.use$Offscale = dat.mer.use$Means_Avg == max(dat.mer.use$Means_Avg)
dat.mer.use$Means_Avg[dat.mer.use$Offscale] = log10(offscale.val)

dat.mer.use$Name = sapply(dat.mer.use$Experiment, function(x,y) y$Name[y$ID == x], y = dat.key )
dat.mer.use$F6_group = as.factor(sapply(dat.mer.use$Experiment, function(x,y) y$F6_group[y$ID == x], y = dat.key ))
group_text_key = c('Training Data', 'Screening', 'Evolution', 'Evolution-GC', 'Gradient', 'Gradient-GC')
dat.mer.use$F6_group_text = group_text_key[dat.mer.use$F6_group]
dat.mer.use$F6_group_text = factor(dat.mer.use$F6_group_text, 
                                   levels = levels(factor(dat.mer.use$F6_group_text))[c(6,5,1:4)] )

png(filename = 'D:/Promoter Design Data/Figures/PNGs/4A.png',
    units = 'cm', width = 9, height = 6, res = 600)

p = ggplot(dat.mer.use, aes(F6_group_text, 10^Means_Avg, color = F6_group_text)) + #, shape = Offscale)) + 
  
  geom_boxplot(data = dat.mer.use, outlier.size = 0, coef = 0, outlier.shape = NA) + 
  geom_jitter(data = dat.mer.use, width = 0.1, height = 0, size = 0.5) +
  geom_hline(yintercept = offscale.val, lty = 2) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 8, face='bold')) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold'),
        legend.title = element_blank(), axis.title.x=element_blank(), legend.position = 'none') +
  labs(y='Promoter Activity \npGPD Designs') + ylim(c(0,27)) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0))
print(p)
dev.off()

fix.m = function(x) {10^x[!is.na(x)]}
print('Training Data vs. Screening (GPD)')
td = fix.m(dat.mer.use$Means_Avg[dat.mer.use$F6_group_text == 'Training Data'])
sc = fix.m(dat.mer.use$Means_Avg[dat.mer.use$F6_group_text == 'Screening'])
wilcox.test(td, sc, 'less')
median(td)
median(sc)
max(td)
max(sc)
print('Evolution and Gradient medians (GPD)')
ev = fix.m(dat.mer.use$Means_Avg[dat.mer.use$F6_group_text == 'Evolution'])
eg = fix.m(dat.mer.use$Means_Avg[dat.mer.use$F6_group_text == 'Evolution-GC'])
gr = fix.m(dat.mer.use$Means_Avg[dat.mer.use$F6_group_text == 'Gradient'])
gg = fix.m(dat.mer.use$Means_Avg[dat.mer.use$F6_group_text == 'Gradient-GC'])
median(ev)
median(eg)
median(gr)
median(gg)
wilcox.test(ev, gr, 'less')
wilcox.test(eg, gg, 'less')
wilcox.test(ev, eg, 'less')
wilcox.test(gr, gg, 'less')

# 4B: pZEV-I
dat.mer.use = fread('D:/Promoter Design Data/FACS-Seq/final_validation_FACS-Seq_means_just_designs_corrected.csv')
dat.mer.use = dat.mer.use[dat.mer.use$Experiment %in% zev.induced.ex,]
offscale.val = 130
dat.mer.use$Offscale = dat.mer.use$Means_B == max(dat.mer.use$Means_B)
dat.mer.use$Means_B[dat.mer.use$Offscale] = log10(offscale.val)

dat.mer.use$Name = sapply(dat.mer.use$Experiment, function(x,y) y$Name[y$ID == x], y = dat.key )
dat.mer.use$F6_group = as.factor(sapply(dat.mer.use$Experiment, function(x,y) y$F6_group[y$ID == x], y = dat.key ))
group_text_key = c('Training Data', 'Screening', 'Evolution', 'Evolution-GC', 
                   'Gradient', 'Gradient-GC', 'Gradient*')
dat.mer.use$F6_group_text = group_text_key[dat.mer.use$F6_group]
dat.mer.use$F6_group_text = factor(dat.mer.use$F6_group_text, 
                                   levels = levels(factor(dat.mer.use$F6_group_text))[c(7,6,1:5)] )

png(filename = 'D:/Promoter Design Data/Figures/PNGs/4B.png',
    units = 'cm', width = 9, height = 6, res = 600)
p = ggplot(dat.mer.use, aes(F6_group_text, 10^Means_B, color = F6_group_text)) + 
  geom_boxplot(outlier.size = 0, coef = 0, outlier.shape = NA) + 
  geom_jitter(width = 0.1, height = 0, size = 0.5) +
  geom_hline(yintercept = offscale.val, lty = 2) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 8, face='bold')) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold'),
        legend.title = element_blank(), axis.title.x=element_blank(), legend.position = 'none') + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  labs(y='Induced Promoter Activity \npZEV-Induced Designs') + ylim(c(0,140))
print(p)
dev.off()

print('Training Data vs. Screening (ZEV-I)')
td = fix.m(dat.mer.use$Means_B[dat.mer.use$F6_group_text == 'Training Data'])
sc = fix.m(dat.mer.use$Means_B[dat.mer.use$F6_group_text == 'Screening'])
wilcox.test(td, sc, 'greater')
median(td)
median(sc)
max(td)
max(sc)
print('Evolution and Gradient medians (ZEV-I)')
ev = fix.m(dat.mer.use$Means_B[dat.mer.use$F6_group_text == 'Evolution'])
eg = fix.m(dat.mer.use$Means_B[dat.mer.use$F6_group_text == 'Evolution-GC'])
gr = fix.m(dat.mer.use$Means_B[dat.mer.use$F6_group_text == 'Gradient'])
gg = fix.m(dat.mer.use$Means_B[dat.mer.use$F6_group_text == 'Gradient-GC'])
g_star = fix.m(dat.mer.use$Means_B[dat.mer.use$F6_group_text == 'Gradient*'])
median(ev)
median(eg)
median(gr)
median(gg)
median(g_star)
wilcox.test(ev, gr, 'less')
wilcox.test(eg, gg, 'less')
wilcox.test(ev, eg, 'less')
wilcox.test(gr, gg, 'greater')
wilcox.test(gr, g_star, 'less')


# 4C: pZEV-AR
dat.mer.use = fread('D:/Promoter Design Data/FACS-Seq/final_validation_FACS-Seq_means_just_designs_corrected.csv')
dat.mer.use = dat.mer.use[dat.mer.use$Experiment %in% zev.ar.ex,]
dat.mer.use$Name = sapply(dat.mer.use$Experiment, function(x,y) y$Name[y$ID == x], y = dat.key )
dat.mer.use$F6_group_AR = as.factor(sapply(dat.mer.use$Experiment, function(x,y) y$F6_group_AR[y$ID == x], y = dat.key ))
group_text_key = c('Training Data', 'Screening', 'Evolution', 'Evolution-GC', 
                   'Gradient', 'Gradient-GC', 'Gradient*')
dat.mer.use$F6_group_text = group_text_key[dat.mer.use$F6_group_AR]
dat.mer.use$F6_group_text = factor(dat.mer.use$F6_group_text, 
                                   levels = levels(factor(dat.mer.use$F6_group_text))[c(7,6,1:5)] )

png(filename = 'D:/Promoter Design Data/Figures/PNGs/4C.png',
    units = 'cm', width = 9, height = 6, res = 600)
p = ggplot(dat.mer.use, aes(F6_group_text, 10^AR, color = F6_group_text)) + 
  geom_boxplot(outlier.size = 0, coef = 0, outlier.shape = NA) + 
  geom_jitter(width = 0.1, height = 0, size = 0.5) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 8, face='bold')) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold'),
        legend.title = element_blank(), axis.title.x=element_blank(), legend.position = 'none') +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  labs(y='Activation Ratio \npZEV-Activation Ratio Designs')
print(p)
dev.off()

over.thresh = function(x,thresh) { sum(x > thresh)/length(x)}

print('Training Data vs. Screening (ZEV-AR)')
td = fix.m(dat.mer.use$AR[dat.mer.use$F6_group_text == 'Training Data'])
sc = fix.m(dat.mer.use$AR[dat.mer.use$F6_group_text == 'Screening'])
wilcox.test(td, sc, 'less')
median(td)
median(sc)
max(td)
max(sc)
print('Evolution and Gradient medians (ZEV-AR)')
ev = fix.m(dat.mer.use$AR[dat.mer.use$F6_group_text == 'Evolution'])
eg = fix.m(dat.mer.use$AR[dat.mer.use$F6_group_text == 'Evolution-GC'])
gr = fix.m(dat.mer.use$AR[dat.mer.use$F6_group_text == 'Gradient'])
gg = fix.m(dat.mer.use$AR[dat.mer.use$F6_group_text == 'Gradient-GC'])
g_star = fix.m(dat.mer.use$AR[dat.mer.use$F6_group_text == 'Gradient*'])
median(ev)
median(eg)
median(gr)
median(gg)
median(g_star)
wilcox.test(ev, gr, 'less')
wilcox.test(eg, gg, 'less')
wilcox.test(ev, eg, 'less')
wilcox.test(gr, gg, 'less')
wilcox.test(gr, g_star, 'less')
over.thresh(td, 100)
over.thresh(sc, 100)
over.thresh(ev, 100)
over.thresh(eg, 100)
over.thresh(gr, 100)
over.thresh(gg, 100)
over.thresh(g_star, 100)

# 4D: off-scale issue for pZEV designs

dat.mer.use$ismin = dat.mer.use$Means_A == min(dat.mer.use$Means_A)
dat.mer.use$`Uninduced Strength` = 'Measured'
dat.mer.use$`Uninduced Strength`[dat.mer.use$ismin] = 'Lower Bound'

png(filename = 'D:/Promoter Design Data/Figures/PNGs/4D.png',
    units = 'cm', width = 9.5, height = 6, res = 600)
p = ggplot(dat.mer.use, aes(x = F6_group_text, fill = `Uninduced Strength`)) + geom_bar() + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face='bold')) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold'),
        legend.text = element_text(size=6), legend.title = element_text(size=8, face='bold'),
        legend.key.size = unit(0.1,'cm')) +
  labs(y='Number of Sequences', fill = 'Uninduced\nActivity in\npZEV Constructs') +
  theme(axis.text.x = element_text(angle = -45, hjust = 0), axis.title.x=element_blank())
print(p)
dev.off()

# Counts
# how many training data?
at = length(dat.mer.use$`Uninduced Strength`[dat.mer.use$F6_group_text == 'Training Data'])
# how many training data not measured?
bt = sum(dat.mer.use$`Uninduced Strength`[dat.mer.use$F6_group_text == 'Training Data'] == 'Lower Bound')

print(at); print(bt); print(bt/at)

a = length(dat.mer.use$`Uninduced Strength`[dat.mer.use$F6_group_text == 'Screening'])
b = sum(dat.mer.use$`Uninduced Strength`[dat.mer.use$F6_group_text == 'Screening'] == 'Lower Bound')
print(b/a) # 54.9%
fisher.test(matrix(nrow = 2, ncol = 2, c(at - bt, bt, a - b, b)), 'less') # 3e-5

a = length(dat.mer.use$`Uninduced Strength`[dat.mer.use$F6_group_text == 'Evolution'])
b = sum(dat.mer.use$`Uninduced Strength`[dat.mer.use$F6_group_text == 'Evolution'] == 'Lower Bound')
print(b/a) # 46.7%
fisher.test(matrix(nrow = 2, ncol = 2, c(at - bt, bt, a - b, b)), 'less') # 0.01

a = length(dat.mer.use$`Uninduced Strength`[dat.mer.use$F6_group_text == 'Evolution-GC'])
b = sum(dat.mer.use$`Uninduced Strength`[dat.mer.use$F6_group_text == 'Evolution-GC'] == 'Lower Bound')
print(b/a) # 59.8%
fisher.test(matrix(nrow = 2, ncol = 2, c(at - bt, bt, a - b, b)), 'less') # 3e-6

a = length(dat.mer.use$`Uninduced Strength`[dat.mer.use$F6_group_text == 'Gradient'])
b = sum(dat.mer.use$`Uninduced Strength`[dat.mer.use$F6_group_text == 'Gradient'] == 'Lower Bound')
print(b/a) # 71.7%
fisher.test(matrix(nrow = 2, ncol = 2, c(at - bt, bt, a - b, b)), 'less') # 2e-7

a = length(dat.mer.use$`Uninduced Strength`[dat.mer.use$F6_group_text == 'Gradient-GC'])
b = sum(dat.mer.use$`Uninduced Strength`[dat.mer.use$F6_group_text == 'Gradient-GC'] == 'Lower Bound')
print(b/a) # 60.9%
fisher.test(matrix(nrow = 2, ncol = 2, c(at - bt, bt, a - b, b)), 'less') # 2e-6

a = length(dat.mer.use$`Uninduced Strength`[dat.mer.use$F6_group_text == 'Gradient*'])
b = sum(dat.mer.use$`Uninduced Strength`[dat.mer.use$F6_group_text == 'Gradient*'] == 'Lower Bound')
print(b/a) # 68.2%
fisher.test(matrix(nrow = 2, ncol = 2, c(at - bt, bt, a - b, b)), 'less') # 2e-8









