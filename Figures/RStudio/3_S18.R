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
dat.key = fread('Fig 3 key.csv')

gpd.strong.ex = c(0, 1,5,19,21,23,25,27) # added 2019-11-08
zev.induced.ex = c(9, 16, 30, 32, 34, 36, 38, 39)
zev.ar.ex = c(10, 17, 41,43,45,47,49,50)
zev.2D.ex = zev.ar.ex #c(9,10,16,17,30, 32, 34, 36, 38, 39, 41,43,45,47,49,50)

dat.mer.use = fread('D:/Promoter Design Data/FACS-Seq/final_validation_FACS-Seq_means_just_designs_corrected.csv')

#offscale.val = 25
dat.mer.use = dat.mer.use[dat.mer.use$Experiment %in% gpd.strong.ex,]
dat.mer.use$Offscale = dat.mer.use$Means_Avg == max(dat.mer.use$Means_Avg)
#dat.mer.use$Means_Avg[dat.mer.use$Offscale] = log10(offscale.val)
print('3A Offscale')
print(sum(dat.mer.use$Offscale))
print(table(dat.mer.use$Experiment[dat.mer.use$Offscale]))
dat.mer.use = dat.mer.use[!dat.mer.use$Offscale,]

dat.mer.use$Name = sapply(dat.mer.use$Experiment, function(x,y) y$Name[y$ID == x], y = dat.key )
dat.mer.use$F6_group = as.factor(sapply(dat.mer.use$Experiment, function(x,y) y$F6_group[y$ID == x], y = dat.key ))
group_text_key = c('Training data', 'Screening', 'Evolution', 'Evolution-GC', 'Gradient', 'Gradient-GC')
dat.mer.use$F6_group_text = group_text_key[dat.mer.use$F6_group]
dat.mer.use$F6_group_text = factor(dat.mer.use$F6_group_text, 
                                   levels = levels(factor(dat.mer.use$F6_group_text))[c(6,5,1:4)] )

stat_box_data = function(y) { return(data.frame(y = max(y) + 3, label = length(y))) }

png(filename = 'D:/Promoter Design Data/Figures/PNGs/3A.png',
    units = 'cm', width = 8.8, height = 5, res = 600)

p = ggplot(dat.mer.use, aes(F6_group_text, 10^Means_Avg)) + #, color = F6_group_text)) + #, shape = Offscale)) + 
  
  geom_boxplot(data = dat.mer.use, outlier.size = 0, outlier.shape = NA) + # coef = 0, 
  geom_jitter(data = dat.mer.use, width = 0.1, height = 0, size = 0.5) +
  #geom_hline(yintercept = offscale.val, lty = 2) +
  # leave numbers off in final version so they can be added in Illustrator
  #stat_summary(fun.data = stat_box_data, geom = "text", hjust = 0.5, vjust = 0.9, size=3) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 8, face='bold')) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold'),
        legend.title = element_blank(), axis.title.x=element_blank(), legend.position = 'none') +
  labs(y='Promoter activity \npGPD designs') + ylim(c(0,20)) + # ylim(c(0,27)) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0))
print(p)
dev.off()

fix.m = function(x) {10^x[!is.na(x)]}
print('Training Data vs. Screening (GPD)')
td = fix.m(dat.mer.use$Means_Avg[dat.mer.use$F6_group_text == 'Training data'])
sc = fix.m(dat.mer.use$Means_Avg[dat.mer.use$F6_group_text == 'Screening'])
wilcox.test(td, sc, 'two.sided')$p.value # 1.178708e-06
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
#wilcox.test(ev, gr, 'less')
#wilcox.test(eg, gg, 'less')
#wilcox.test(ev, eg, 'less')
#wilcox.test(gr, gg, 'less')

# 3B: pZEV-I
dat.mer.use = fread('D:/Promoter Design Data/FACS-Seq/final_validation_FACS-Seq_means_just_designs_corrected.csv')
dat.mer.use = dat.mer.use[dat.mer.use$Experiment %in% zev.induced.ex,]
#offscale.val = 130
dat.mer.use$Offscale = dat.mer.use$Means_B == max(dat.mer.use$Means_B)
#dat.mer.use$Means_B[dat.mer.use$Offscale] = log10(offscale.val)
print('3B Offscale')
print(sum(dat.mer.use$Offscale))
print(table(dat.mer.use$Experiment[dat.mer.use$Offscale]))
dat.mer.use = dat.mer.use[!dat.mer.use$Offscale,]

dat.mer.use$Name = sapply(dat.mer.use$Experiment, function(x,y) y$Name[y$ID == x], y = dat.key )
dat.mer.use$F6_group = as.factor(sapply(dat.mer.use$Experiment, function(x,y) y$F6_group[y$ID == x], y = dat.key ))
group_text_key = c('Training data', 'Screening', 'Evolution', 'Evolution-GC', 
                   'Gradient', 'Gradient-GC', 'Gradient*')
dat.mer.use$F6_group_text = group_text_key[dat.mer.use$F6_group]
dat.mer.use$F6_group_text = factor(dat.mer.use$F6_group_text, 
                                   levels = levels(factor(dat.mer.use$F6_group_text))[c(7,6,1:5)] )

stat_box_data = function(y) { return(data.frame(y = max(y) + 15, label = length(y))) }

png(filename = 'D:/Promoter Design Data/Figures/PNGs/3B.png',
    units = 'cm', width = 8.8, height = 5, res = 600)
p = ggplot(dat.mer.use, aes(F6_group_text, 10^Means_B)) + #, color = F6_group_text)) + 
  geom_boxplot(outlier.size = 0, outlier.shape = NA) + # coef = 0, 
  geom_jitter(width = 0.1, height = 0, size = 0.5) +
  #geom_hline(yintercept = offscale.val, lty = 2) +
  # leave numbers off in final version so they can be added in Illustrator
  #stat_summary(fun.data = stat_box_data, geom = "text", hjust = 0.5, vjust = 0.9, size=3) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 8, face='bold')) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold'),
        legend.title = element_blank(), axis.title.x=element_blank(), legend.position = 'none') + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  labs(y='Induced promoter activity \npZEV-induced designs') + ylim(c(0,100)) # + ylim(c(0,140))
print(p)
dev.off()

print('Training Data vs. Screening (ZEV-I)')
td = fix.m(dat.mer.use$Means_B[dat.mer.use$F6_group_text == 'Training data'])
sc = fix.m(dat.mer.use$Means_B[dat.mer.use$F6_group_text == 'Screening'])
wilcox.test(td, sc, 'two.sided')$p.value
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
wilcox.test(ev, gr, 'two.sided')$p.value # 0.0004085065
wilcox.test(eg, gg, 'two.sided')$p.value # 0.511687
wilcox.test(ev, eg, 'two.sided')$p.value # 0.04484599
wilcox.test(gr, gg, 'two.sided')$p.value # 0.01216047
wilcox.test(gr, g_star, 'two.sided')$p.value # 1.697595e-05

# 3C: pZEV-AR
dat.mer.use = fread('D:/Promoter Design Data/FACS-Seq/final_validation_FACS-Seq_means_just_designs_corrected.csv')
dat.mer.use = dat.mer.use[dat.mer.use$Experiment %in% zev.ar.ex,]
dat.mer.use$Name = sapply(dat.mer.use$Experiment, function(x,y) y$Name[y$ID == x], y = dat.key )
dat.mer.use$F6_group_AR = as.factor(sapply(dat.mer.use$Experiment, function(x,y) y$F6_group_AR[y$ID == x], y = dat.key ))
group_text_key = c('Training data', 'Screening', 'Evolution', 'Evolution-GC', 
                   'Gradient', 'Gradient-GC', 'Gradient*')
dat.mer.use$F6_group_text = group_text_key[dat.mer.use$F6_group_AR]
dat.mer.use$F6_group_text = factor(dat.mer.use$F6_group_text, 
                                   levels = levels(factor(dat.mer.use$F6_group_text))[c(7,6,1:5)] )

stat_box_data = function(y) { return(data.frame(y = min(max(y) + 40, 200), label = length(y))) }

png(filename = 'D:/Promoter Design Data/Figures/PNGs/3C.png',
    units = 'cm', width = 8.8, height = 5, res = 600)
p = ggplot(dat.mer.use, aes(F6_group_text, 10^AR)) + #, color = F6_group_text)) + 
  geom_boxplot(outlier.size = 0, outlier.shape = NA) + # coef = 0, 
  geom_jitter(width = 0.1, height = 0, size = 0.5) +
  theme_bw() + 
  # leave numbers off in final version so they can be added in Illustrator
  #stat_summary(fun.data = stat_box_data, geom = "text", hjust = 0.5, vjust = 0.9, size=3) +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face='bold')) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold'),
        legend.title = element_blank(), axis.title.x=element_blank(), legend.position = 'none') +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  labs(y='Activation ratio \npZEV-AR designs')
print(p)
dev.off()

over.thresh = function(x,thresh) { sum(x > thresh)/length(x)}

print('Training Data vs. Screening (ZEV-AR)')
td = fix.m(dat.mer.use$AR[dat.mer.use$F6_group_text == 'Training data'])
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
wilcox.test(ev, gr, 'less')$p.value
wilcox.test(eg, gg, 'less')$p.value
wilcox.test(ev, eg, 'less')$p.value
wilcox.test(gr, gg, 'less')$p.value
wilcox.test(gr, g_star, 'less')$p.value
over.thresh(td, 100)
over.thresh(sc, 100)
over.thresh(ev, 100)
over.thresh(eg, 100)
over.thresh(gr, 100)
over.thresh(gg, 100)
over.thresh(g_star, 100)

# S18: off-scale issue for pZEV designs

dat.mer.use$ismin = dat.mer.use$Means_A == min(dat.mer.use$Means_A)
dat.mer.use$`Uninduced Strength` = 'Measured'
dat.mer.use$`Uninduced Strength`[dat.mer.use$ismin] = 'Lower Bound'
stat_box_data = function(y) { return(data.frame(y = max(y), label = length(y))) }

# df for geom_text to add bar counts
df.small = data.frame(F6_group_text = dat.mer.use$F6_group_text,
                      us=dat.mer.use$`Uninduced Strength`, val = 1)
df.lab = melt(dcast(df.small, us ~ F6_group_text, value.var = 'val'))
setnames(df.lab, c('us', 'variable'), c('Uninduced Strength', 'F6_group_text'))
df.lab$pos = 0
is.m = df.lab$`Uninduced Strength` == 'Measured'
df.lab$pos[is.m] = df.lab$value[is.m]/2
isnt.m = which(!is.m)
for(i in 1:length(isnt.m)) {
  idx = which(df.lab$F6_group_text == df.lab$F6_group_text[isnt.m[i]])
  lb = df.lab$value[idx[1]]
  me = df.lab$value[idx[2]]
  ans = me + lb/2
  df.lab$pos[isnt.m[i]] = ans
}

# Source Data
dat.source = dat.mer.use[order(dat.mer.use$F6_group_text, dat.mer.use$`Uninduced Strength`),]
dat.source = data.frame(Seqs=dat.source$Seqs, `Design Type`=dat.source$F6_group_text, Measured=dat.source$`Uninduced Strength`)
write.csv(dat.source, 'D:/Promoter Design Data/Source Data/S18.csv', quote = FALSE, row.names = FALSE)

png(filename = 'D:/Promoter Design Data/Figures/PNGs/S18.png',
    units = 'cm', width = 9.5, height = 6, res = 600)
p = ggplot(dat.mer.use, aes(x = F6_group_text, fill = `Uninduced Strength`)) + geom_bar() + theme_bw() +
  geom_text(data=df.lab,aes(x=F6_group_text,y=pos,label=value),vjust=0.5) +
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
at = length(dat.mer.use$`Uninduced Strength`[dat.mer.use$F6_group_text == 'Training data'])
# how many training data not measured?
bt = sum(dat.mer.use$`Uninduced Strength`[dat.mer.use$F6_group_text == 'Training data'] == 'Lower Bound')

print(at); print(bt); print(bt/at)

a = length(dat.mer.use$`Uninduced Strength`[dat.mer.use$F6_group_text == 'Screening'])
b = sum(dat.mer.use$`Uninduced Strength`[dat.mer.use$F6_group_text == 'Screening'] == 'Lower Bound')
print(b/a) # 54.9%
fisher.test(matrix(nrow = 2, ncol = 2, c(at - bt, bt, a - b, b)), 'two.sided')$p.value # 3.567337e-05

a = length(dat.mer.use$`Uninduced Strength`[dat.mer.use$F6_group_text == 'Evolution'])
b = sum(dat.mer.use$`Uninduced Strength`[dat.mer.use$F6_group_text == 'Evolution'] == 'Lower Bound')
print(b/a) # 46.7%
fisher.test(matrix(nrow = 2, ncol = 2, c(at - bt, bt, a - b, b)), 'two.sided')$p.value # 0.002200129

a = length(dat.mer.use$`Uninduced Strength`[dat.mer.use$F6_group_text == 'Evolution-GC'])
b = sum(dat.mer.use$`Uninduced Strength`[dat.mer.use$F6_group_text == 'Evolution-GC'] == 'Lower Bound')
print(b/a) # 59.8%
fisher.test(matrix(nrow = 2, ncol = 2, c(at - bt, bt, a - b, b)), 'two.sided')$p.value # 3.174218e-06

a = length(dat.mer.use$`Uninduced Strength`[dat.mer.use$F6_group_text == 'Gradient'])
b = sum(dat.mer.use$`Uninduced Strength`[dat.mer.use$F6_group_text == 'Gradient'] == 'Lower Bound')
print(b/a) # 71.7%
fisher.test(matrix(nrow = 2, ncol = 2, c(at - bt, bt, a - b, b)), 'two.sided')$p.value # 2.233967e-07

a = length(dat.mer.use$`Uninduced Strength`[dat.mer.use$F6_group_text == 'Gradient-GC'])
b = sum(dat.mer.use$`Uninduced Strength`[dat.mer.use$F6_group_text == 'Gradient-GC'] == 'Lower Bound')
print(b/a) # 60.9%
fisher.test(matrix(nrow = 2, ncol = 2, c(at - bt, bt, a - b, b)), 'two.sided')$p.value # 2.542183e-06

a = length(dat.mer.use$`Uninduced Strength`[dat.mer.use$F6_group_text == 'Gradient*'])
b = sum(dat.mer.use$`Uninduced Strength`[dat.mer.use$F6_group_text == 'Gradient*'] == 'Lower Bound')
print(b/a) # 68.2%
fisher.test(matrix(nrow = 2, ncol = 2, c(at - bt, bt, a - b, b)), 'two.sided')$p.value # 2.566722e-08









