# Single-sequence testing results:
# Show concordance w/ validation FACS-seq; show performance of designed promoters
# vs. benchmarks

setwd('D:/Promoter Design Data/')
library(data.table)
library(ggplot2)
source('Code/FCS file analysis.R')

get.one.median = function(fn) {
  x = extract.data(load.fromfilename(fn))
  rats = log10(x[,10]/x[,8])
  return(median(rats))
}
w = 9.5
h = 6

dat = fread('Validation/validation_output_fs_means_corrected.csv'); dat$V1 = NULL
dat$VYB_Mean_A = apply(cbind(dat$R1_FALSE, dat$R2_FALSE, dat$R3_FALSE), 1, mean)
dat$VYB_SD_A = apply(cbind(dat$R1_FALSE, dat$R2_FALSE, dat$R3_FALSE), 1, sd)
dat$VYB_Mean_B = apply(cbind(dat$R1_TRUE, dat$R2_TRUE, dat$R3_TRUE), 1, mean)
dat$VYB_SD_B = apply(cbind(dat$R1_TRUE, dat$R2_TRUE, dat$R3_TRUE), 1, sd)
dat$VYB_AR = dat$VYB_Mean_B - dat$VYB_Mean_A
dat$VYB_AR_SD = sqrt(dat$VYB_SD_A^2 + dat$VYB_SD_B^2)
dat$bar.min.a = dat$VYB_Mean_A-dat$VYB_SD_A
dat$bar.max.a = dat$VYB_Mean_A+dat$VYB_SD_A
dat$bar.min.b = dat$VYB_Mean_B-dat$VYB_SD_B
dat$bar.max.b = dat$VYB_Mean_B+dat$VYB_SD_B
dat$bar.min.ar = dat$VYB_AR-dat$VYB_AR_SD
dat$bar.max.ar = dat$VYB_AR+dat$VYB_AR_SD

# Did we choose this sequence because we expect it to be unreliable?
dat$trusted_a = dat$Scale_A == 0
dat$trusted_b = dat$Scale_B == 0

# for building final sequence-strength table: controls
#write.csv(dat[1:9,], 'C:/Users/Ben/Dropbox/Lab/Lab Deliverables/Promoter Manuscript/single_seqs_ctrls.csv')

# 5A
dat.use.a = data.frame(FS = dat$Means_A, VYB = dat$VYB_Mean_A, ymin = dat$bar.min.a, ymax = dat$bar.max.a, 
                       Promoter = dat$Promoter, rep = 'A', Offscale = !(dat$trusted_a))

# don't double-count pGPD promoters
dat.use.b = data.frame(FS = dat$Means_B, VYB = dat$VYB_Mean_B, ymin = dat$bar.min.b, ymax = dat$bar.max.b,
                       Promoter = dat$Promoter, rep = 'B', Offscale = !(dat$trusted_b))
dat.use.b = dat.use.b[dat.use.b$Promoter == 'ZEV',]

dat.use = rbind(dat.use.a, dat.use.b)
dat.use = dat.use[apply(dat.use, 1, function(x) all(!is.na(x))),]

dat.use$Type = 'pGPD'
dat.use$Type[dat.use$Promoter == 'ZEV' & dat.use$rep == 'A'] = 'pZEV - Uninduced'
dat.use$Type[dat.use$Promoter == 'ZEV' & dat.use$rep == 'B'] = 'pZEV - Induced'

mod = lm(dat.use$VYB[!dat.use$Offscale] ~ dat.use$FS[!dat.use$Offscale])
summary(mod) # R2 = 0.92
coefs = mod$coef


val.A.low = min(dat.use.a$FS[!is.na(dat.use.a$FS)])
val.B.low = min(dat.use.b$FS[!is.na(dat.use.b$FS)])
val.B.hi = max(dat.use.b$FS[!is.na(dat.use.b$FS)])
vals.offscale = c(val.A.low, val.B.low, val.B.hi)

png(filename = 'D:/Promoter Design Data/Figures/PNGs/5A.png',
    units = 'cm', width = w, height = h, res = 600)

p = ggplot(dat.use, aes(x = FS, y = VYB, color = Type, shape = Offscale, ymin = ymin, ymax = ymax)) + 
  geom_point(size = 1, stroke = 0) + geom_errorbar(size = 0.3) + theme_bw() + 
  geom_abline(slope = coefs[2], intercept = coefs[1], lty = 2) +
  geom_vline(xintercept = vals.offscale, lty = 2, lwd = 0.5, alpha = 0.5) + 
  theme(plot.title = element_text(hjust = 0.4, size = 8, face='bold')) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold'),
        legend.text = element_text(size=6), legend.title = element_text(size=8, face='bold'),
        legend.key.size = unit(0.2,'cm')) +
  labs(x = 'FACS-seq\nActivity Measurements (log10)', y = 'Individual Testing\nActivity Measurements (log10)')
print(p)
dev.off()

# 5B
# break log scaling
dat = fread('Validation/validation_output_fs_means_corrected.csv'); dat$V1 = NULL

dat$VYB_Mean_A = apply(cbind(dat$R1_FALSE, dat$R2_FALSE, dat$R3_FALSE), 1, mean)
dat$VYB_SD_A = apply(cbind(dat$R1_FALSE, dat$R2_FALSE, dat$R3_FALSE), 1, sd)
dat$VYB_Mean_B = apply(cbind(dat$R1_TRUE, dat$R2_TRUE, dat$R3_TRUE), 1, mean)
dat$VYB_SD_B = apply(cbind(dat$R1_TRUE, dat$R2_TRUE, dat$R3_TRUE), 1, sd)
dat$VYB_AR = dat$VYB_Mean_B - dat$VYB_Mean_A
dat$VYB_AR_SD = sqrt(dat$VYB_SD_A^2 + dat$VYB_SD_B^2)
dat$bar.min.a = 10^(dat$VYB_Mean_A-dat$VYB_SD_A)
dat$bar.max.a = 10^(dat$VYB_Mean_A+dat$VYB_SD_A)
dat$bar.min.b = 10^(dat$VYB_Mean_B-dat$VYB_SD_B)
dat$bar.max.b = 10^(dat$VYB_Mean_B+dat$VYB_SD_B)
dat$bar.min.ar = 10^(dat$VYB_AR-dat$VYB_AR_SD)
dat$bar.max.ar = 10^(dat$VYB_AR+dat$VYB_AR_SD)

dat$VYB_Mean_A = 10^(dat$VYB_Mean_A)
dat$VYB_Mean_B = 10^(dat$VYB_Mean_B)
dat$VYB_AR = 10^(dat$VYB_AR)

# Did we choose this sequence because we expect it to be unreliable?
dat$trusted_a = dat$Scale_A == 0
dat$trusted_b = dat$Scale_B == 0

dat.gpd = dat[dat$Name %in% c('GPD', 'TEF', 'CYC1') | 
                    dat$Experiment.designed %in% c(23,28),]
dat.gpd = dat.gpd[(dat.gpd$trusted_a & dat.gpd$trusted_b) | is.na(dat.gpd$Experiment.designed),]

dat.gpd$`Design Type` = dat.gpd$Experiment.designed
dat.gpd$`Design Type`[is.na(dat.gpd$`Design Type`)] = 'Control'
#dat.gpd$`Design Type`[dat.gpd$`Design Type` == 22] = 'Evolution; GC Filter; Mean'
dat.gpd$`Design Type`[dat.gpd$`Design Type` == 23] = 'Evolution-GC'
#dat.gpd$`Design Type`[dat.gpd$`Design Type` == 24] = 'Gradient; No Filter; Mean'
dat.gpd$`Design Type`[dat.gpd$`Design Type` == 28] = 'Gradient*'
dat.gpd$`Design Type` = factor(dat.gpd$`Design Type`)

dat.gpd$Name = factor(dat.gpd$Name, levels = dat.gpd$Name[order(dat.gpd$Experiment.designed, dat.gpd$VYB_Mean_A)])

# for building final sequence-strength table: pGPD
# Use the full set from S21 instead!
#write.csv(dat.gpd, 'C:/Users/Ben/Dropbox/Lab/Lab Deliverables/Promoter Manuscript/single_seqs_pGPD.csv')

png(filename = 'D:/Promoter Design Data/Figures/PNGs/5B.png',
    units = 'cm', width = w, height = h, res = 600)
xstr = c(paste0('pGPD-', 1:15), 'pCYC1', 'pTEF1', 'pGPD')
p = ggplot(dat.gpd, aes(Name, VYB_Mean_A, color = `Design Type`, ymin = bar.min.a, ymax = bar.max.a)) + 
  geom_point(size = 1, stroke = 0) + geom_errorbar(size = 0.3) + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.4, size = 8, face='bold')) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold'),
        legend.text = element_text(size=6), legend.title = element_text(size=8, face='bold'),
        legend.key.size = unit(0.2,'cm')) +
  theme(axis.text.x=element_text(angle=315),
        axis.ticks.x=element_blank()) +
  scale_x_discrete(labels=xstr) +
  labs(x = '', y='Promoter Activity\npGPD Designs', color = 'Design Type')
print(p)
dev.off()

# 5C
png(filename = 'D:/Promoter Design Data/Figures/PNGs/5C.png',
    units = 'cm', width = w, height = h, res = 600)
xstr = c(paste0('pZEV-I-', 1:11), 'P8', 'P4', 'P3')
dat.zev.i = dat[dat$Name %in% c('ZEV_Pr3', 'ZEV_Pr4', 'ZEV_Pr8') | 
                  dat$Experiment.designed == 34,] #c(34,39)

dat.zev.i = dat.zev.i[(dat.zev.i$trusted_a & dat.zev.i$trusted_b) | is.na(dat.zev.i$Experiment.designed),]

dat.zev.i$Name = factor(dat.zev.i$Name, levels = dat.zev.i$Name[order(dat.zev.i$Experiment.designed, dat.zev.i$VYB_Mean_B)])
dat.zev.i$`Design Type` = dat.zev.i$Experiment.designed
dat.zev.i$`Design Type`[is.na(dat.zev.i$`Design Type`)] = 'Control'
dat.zev.i$`Design Type`[dat.zev.i$`Design Type` == 34] = 'Evolution-GC'
#dat.zev.i$`Design Type`[dat.zev.i$`Design Type` == 39] = 'Gradient; No Filter'
dat.zev.i$`Design Type` = factor(dat.zev.i$`Design Type`)

# for building final sequence-strength table: pZEV-I
dat.zev.i$Seqs.final = sapply(dat.zev.i$Seqs, function(x) substr(x,59,nchar(x)-59))
#write.csv(dat.zev.i, 'C:/Users/Ben/Dropbox/Lab/Lab Deliverables/Promoter Manuscript/single_seqs_pZEV_I.csv')

p = ggplot(dat.zev.i, aes(Name, VYB_Mean_B, color = `Design Type`, 
                      ymin = bar.min.b, ymax = bar.max.b)) + 
  geom_point(size = 1, stroke = 0) + geom_errorbar() + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.4, size = 8, face='bold')) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold'),
        legend.text = element_text(size=6), legend.title = element_text(size=8, face='bold'),
        legend.key.size = unit(0.2,'cm')) +
  theme(axis.text.x=element_text(angle=315),
        axis.ticks.x=element_blank()) +
  scale_x_discrete(labels=xstr) +
  labs(x = '', y='Induced Promoter Activity \npZEV-Induced Designs', color = 'Design Type')
print(p)
dev.off()

# 5D
dat.zev.ar = dat[(is.na(dat$Seqs) & dat$Name %in% c('ZEV_Pr3', 'ZEV_Pr4', 'ZEV_Pr8')) | 
                   dat$Experiment.designed == 44,]# | dat$ZEV_AR_offscale_cand | dat$ZEV_I_AR_cand,]
xstr = c(paste0('pZEV-AR-', 1:10), 'P8', 'P4', 'P3')
dat.zev.ar$Name = factor(dat.zev.ar$Name, levels = dat.zev.ar$Name[order(dat.zev.ar$Experiment.designed, dat.zev.ar$VYB_AR)])
dat.zev.ar$`Design Type` = dat.zev.ar$Experiment.designed
dat.zev.ar$`Design Type`[is.na(dat.zev.ar$`Design Type`)] = 'Control'
#dat.zev.ar$`Design Type`[dat.zev.ar$ZEV_AR_offscale_cand] = 'Off-scale AR in FACS-seq'
#dat.zev.ar$`Design Type`[dat.zev.ar$ZEV_I_AR_cand] = 'High AR in FACS-seq'
dat.zev.ar$`Design Type`[dat.zev.ar$Experiment.designed == 44] = 'Evolution-GC'
dat.zev.ar$`Design Type` = factor(dat.zev.ar$`Design Type`)

# for building final sequence-strength table: pZEV-AR
dat.zev.ar$Seqs.final = sapply(dat.zev.ar$Seqs, function(x) substr(x,59,nchar(x)-59))
#write.csv(dat.zev.ar, 'C:/Users/Ben/Dropbox/Lab/Lab Deliverables/Promoter Manuscript/single_seqs_pZEV_AR.csv')

png(filename = 'D:/Promoter Design Data/Figures/PNGs/5D.png',
    units = 'cm', width = w, height = h, res = 600)

p = ggplot(dat.zev.ar, aes(Name, VYB_AR, color = `Design Type`, ymin = bar.min.ar, ymax = bar.max.ar)) + 
  geom_point(size = 1, stroke = 0) + geom_errorbar() + theme_bw() +
  theme(plot.title = element_text(hjust = 0.4, size = 8, face='bold')) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold'),
        legend.text = element_text(size=6), legend.title = element_text(size=8, face='bold'),
        legend.key.size = unit(0.2,'cm')) +
  theme(axis.text.x=element_text(angle=315),
        axis.ticks.x=element_blank()) +
  scale_x_discrete(labels=xstr) +
  labs(x = '', y='Activation Ratio \npZEV-Activation Ratio Designs', color = 'Design Type')
print(p)
dev.off()

# 5E
# get the pBK71 baseline
setwd('D:/Datasets/20190110 VYB (pBK71 prelim)/')
setwd('D:/Promoter Design Data/Final Flow Validation FCS/Baseline')
fns = dir()
meds = sapply(fns, get.one.median)
meds = 10^meds
cutoff = mean(meds)
cutoff.sd = sd(meds)

png(filename = 'D:/Promoter Design Data/Figures/PNGs/5E.png',
    units = 'cm', width = w, height = h, res = 600)

dat.zev.ar = dat.zev.ar[dat.zev.ar$Name != 'ZEV_Pr8',]

p = ggplot(dat.zev.ar, aes(x = VYB_Mean_A, y = VYB_Mean_B, color = `Design Type`, #Target, 
                       xmin = bar.min.a, xmax = bar.max.a, 
                       ymin = bar.min.b, ymax = bar.max.b)) + 
  geom_point(size = 1, stroke = 0, alpha = 0.8) + geom_errorbar(size = 0.3) + geom_errorbarh(size = 0.3) + 
  geom_abline(slope = 1, intercept = 0, lty = 2, size = 0.3) + theme_bw() + 
  geom_vline(xintercept = cutoff - cutoff.sd, lty = 2, size = 0.3) +
  geom_vline(xintercept = cutoff + cutoff.sd, lty = 2, size = 0.3) +
  theme(plot.title = element_text(hjust = 0.4, size = 8, face='bold')) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold'),
        legend.text = element_text(size=6), legend.title = element_text(size=8, face='bold'),
        legend.key.size = unit(0.2,'cm')) +
  labs(x = 'High-AR pZEV Design\n Uninduced Activity', 
       y = 'High-AR pZEV Design\n Induced Activity', color = 'Design Type')
print(p)
dev.off()

### Information for text ###
# 95% confidence intervals, with the correct log-normality
get.cids = function(df) {
  df$CI_A_low = 10^(log10(df$VYB_Mean_A) - qnorm(0.975)*df$VYB_SD_A/sqrt(3)) # qnorm(0.975) is about 1.96
  df$CI_A_hi = 10^(log10(df$VYB_Mean_A) + qnorm(0.975)*df$VYB_SD_A/sqrt(3))
  df$CI_B_low = 10^(log10(df$VYB_Mean_B) - qnorm(0.975)*df$VYB_SD_B/sqrt(3))
  df$CI_B_hi = 10^(log10(df$VYB_Mean_B) + qnorm(0.975)*df$VYB_SD_B/sqrt(3))
  return(df)
}

# for constant controls, we have 6 measurements: use all of them
get.cids.ctrl = function(df) {
  tmp = cbind(df$R1_FALSE, df$R1_TRUE, df$R2_FALSE, df$R2_TRUE, df$R3_FALSE, df$R3_TRUE)
  m = apply(tmp,1,mean)
  s = apply(tmp,1,sd)
  df$CI_low = 10^(m - qnorm(0.975)*s/sqrt(6)) # qnorm(0.975) is about 1.96
  df$CI_hi = 10^(m + qnorm(0.975)*s/sqrt(6))
  return(df)
}

# Designed pGPD
dat.tmp = dat.gpd[dat.gpd$Experiment.designed == 23,]
dat.tmp = dat.tmp[dat.tmp$VYB_Mean_A %in% range(dat.tmp$VYB_Mean_A),]
dat.tmp = get.cids(dat.tmp)
print(dat.tmp$VYB_Mean_A)
print(dat.tmp$CI_A_low)
print(dat.tmp$CI_A_hi)

# pTEF and pGPD
dat.tmp = dat[dat$Name %in% c('TEF','GPD'),];dat.tmp = get.cids.ctrl(dat.tmp)
#print(dat.tmp$VYB_Mean_A)
print(10^apply(cbind(dat.tmp$R1_FALSE, dat.tmp$R1_TRUE,dat.tmp$R2_FALSE, dat.tmp$R2_TRUE,dat.tmp$R3_FALSE, dat.tmp$R3_TRUE),
            1,mean))
print(dat.tmp$CI_low)
print(dat.tmp$CI_hi)

# Gradient
dat.tmp = dat.gpd[dat.gpd$Experiment.designed == 28,]; dat.tmp = get.cids(dat.tmp)
print(dat.tmp$VYB_Mean_A)
print(dat.tmp$CI_A_low)
print(dat.tmp$CI_A_hi)
# Gradient t-test vs. GPD
df = dat.gpd[dat.gpd$Experiment.designed == 28,]
g = dat.gpd[dat.gpd$Name == 'GPD',]
tmp = cbind(df$R1_FALSE, df$R2_FALSE, df$R3_FALSE)
g = cbind(g$R1_FALSE, g$R1_TRUE, g$R2_FALSE, g$R2_TRUE, g$R3_FALSE, g$R3_TRUE)
for(i in 1:nrow(tmp)) {
  print(10^mean(tmp[i,]))
  print(t.test(tmp[i,], g, 'greater')$p.value)
  print(t.test(tmp[i,], g, 'less')$p.value)
}



# Designed pZEV-I
dat.tmp = dat.zev.i[dat.zev.i$Experiment.designed == 34,]
dat.tmp = dat.tmp[dat.tmp$VYB_Mean_B %in% range(dat.tmp$VYB_Mean_B),]
dat.tmp = get.cids(dat.tmp)
print(dat.tmp$VYB_Mean_B)
print(dat.tmp$CI_B_low)
print(dat.tmp$CI_B_hi)

# p-values vs. P4 and P8
dat.tmp = dat.zev.i[dat.zev.i$Experiment.designed == 34,]
dat.ctrl = dat[dat$Name %in% c('ZEV_Pr3', 'ZEV_Pr4', 'ZEV_Pr8'),]
z = cbind(dat.ctrl$R1_TRUE, dat.ctrl$R2_TRUE, dat.ctrl$R3_TRUE)
tmp = cbind(dat.tmp$R1_TRUE, dat.tmp$R2_TRUE, dat.tmp$R3_TRUE)
pvals = matrix(nrow = nrow(z), ncol = nrow(tmp))
for(i in 1:nrow(z)) { for(j in 1:nrow(tmp)) {
  pvals[i,j] = t.test(z[i,], tmp[j,], 'less')$p.value
}}
colnames(pvals) = dat.tmp$VYB_Mean_B; rownames(pvals) = dat.ctrl$Name



# CI for Pr3
dat.ctrl = get.cids(dat.ctrl)
print(dat.ctrl$Name)
print(dat.ctrl$VYB_Mean_B)
print(dat.ctrl$CI_B_low)
print(dat.ctrl$CI_B_hi)

# Designed pZEV-AR
get.cids.ar = function(df) {
  tmp = cbind(df$R1_TRUE - df$R1_FALSE, df$R2_TRUE - df$R2_FALSE, df$R3_TRUE - df$R3_FALSE)
  m = apply(tmp,1,mean)
  s = apply(tmp,1,sd)
  df$CI_AR_low = 10^(m - qnorm(0.975)*s/sqrt(3)) # qnorm(0.975) is about 1.96
  df$CI_AR_hi = 10^(m + qnorm(0.975)*s/sqrt(3))
  return(df)
}


dat.zev.ar = dat[(is.na(dat$Seqs) & dat$Name %in% c('ZEV_Pr3', 'ZEV_Pr4', 'ZEV_Pr8')) | 
                   dat$Experiment.designed == 44,]# | dat$ZEV_AR_offscale_cand | dat$ZEV_I_AR_cand,]

dat.zev.ar$Name = factor(dat.zev.ar$Name, levels = dat.zev.ar$Name[order(dat.zev.ar$Experiment.designed, dat.zev.ar$VYB_AR)])
dat.zev.ar$`Design Type` = dat.zev.ar$Experiment.designed
dat.zev.ar$`Design Type`[is.na(dat.zev.ar$`Design Type`)] = 'Control'
dat.zev.ar$`Design Type`[dat.zev.ar$Experiment.designed == 44] = 'Evolution-GC'
dat.zev.ar$`Design Type` = factor(dat.zev.ar$`Design Type`)

dat.tmp = dat.zev.ar[dat.zev.ar$Experiment.designed == 44,]
dat.tmp = dat.tmp[dat.tmp$VYB_AR %in% range(dat.tmp$VYB_AR),]
dat.tmp = get.cids.ar(dat.tmp)
print(dat.tmp$VYB_AR)
print(dat.tmp$CI_AR_low)
print(dat.tmp$CI_AR_hi)



# AR for controls
dat.ctrl = dat[dat$Name %in% c('ZEV_Pr3', 'ZEV_Pr4', 'ZEV_Pr8'),]
dat.ctrl = get.cids.ar(dat.ctrl)
print(dat.ctrl$VYB_AR)
print(dat.ctrl$CI_AR_low)
print(dat.ctrl$CI_AR_hi)

# t-tests for uninduced and induced
# p-values vs. P4 and P8

# UNINDUCED
dat.tmp = dat.zev.ar[dat.zev.ar$Experiment.designed == 44,]
dat.ctrl = dat[dat$Name %in% c('ZEV_Pr3', 'ZEV_Pr4'),]
z = cbind(dat.ctrl$R1_FALSE, dat.ctrl$R2_FALSE, dat.ctrl$R3_FALSE)
tmp = cbind(dat.tmp$R1_FALSE, dat.tmp$R2_FALSE, dat.tmp$R3_FALSE)
pvals = matrix(nrow = nrow(z), ncol = nrow(tmp))
for(i in 1:nrow(z)) { for(j in 1:nrow(tmp)) {
  pvals[i,j] = t.test(z[i,], tmp[j,], 'greater')$p.value
}}
colnames(pvals) = dat.tmp$VYB_Mean_A; rownames(pvals) = dat.ctrl$Name



# INDUCED
dat.tmp = dat.zev.ar[dat.zev.ar$Experiment.designed == 44,]
dat.ctrl = dat[dat$Name %in% c('ZEV_Pr3', 'ZEV_Pr4'),]
z = cbind(dat.ctrl$R1_TRUE, dat.ctrl$R2_TRUE, dat.ctrl$R3_TRUE)
tmp = cbind(dat.tmp$R1_TRUE, dat.tmp$R2_TRUE, dat.tmp$R3_TRUE)
pvals = matrix(nrow = nrow(z), ncol = nrow(tmp))
for(i in 1:nrow(z)) { for(j in 1:nrow(tmp)) {
  pvals[i,j] = t.test(z[i,], tmp[j,], 'greater')$p.value
}}
colnames(pvals) = dat.tmp$VYB_Mean_B; rownames(pvals) = dat.ctrl$Name



