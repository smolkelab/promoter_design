# Single-sequence testing results:
# Show concordance w/ validation FACS-seq; show performance of designed promoters
# vs. benchmarks

setwd('D:/Promoter Design Data/')
library(data.table)
library(ggplot2)
library(reshape2)
source('Code/FCS file analysis.R')

get.one.median = function(fn) {
  x = extract.data(load.fromfilename(fn))
  rats = log10(x[,10]/x[,8])
  return(median(rats))
}

melt.for.jitter = function(df, mode) {
  if (mode == 'Uninduced') {df = data.frame(Name=df$Name, `Design Type`=df$`Design Type`, R1=df$R1_FALSE, R2=df$R2_FALSE, R3=df$R3_FALSE)}
  if (mode == 'Induced') {df = data.frame(Name=df$Name, `Design Type`=df$`Design Type`, R1=df$R1_TRUE, R2=df$R2_TRUE, R3=df$R3_TRUE)}
  if (mode == 'AR') { df = data.frame(Name=df$Name, `Design Type`=df$`Design Type`, R1=df$R1_TRUE - df$R1_FALSE, 
                                      R2=df$R2_TRUE - df$R2_FALSE, R3=df$R3_TRUE - df$R3_FALSE) }
  if (mode == 'All') {df = data.frame(Name=df$Name, `Design Type`=df$`Design Type`, R1=df$R1_FALSE, R2=df$R2_FALSE, R3=df$R3_FALSE,
                                      R4=df$R1_TRUE, R5=df$R2_TRUE, R6=df$R3_TRUE)}
  if (mode == 'All') {df = melt(df, id.vars = c('Name', 'Design.Type'), measure.vars = c('R1', 'R2', 'R3', 'R4', 'R5', 'R6'))}
  else { df = melt(df, id.vars = c('Name', 'Design.Type'), measure.vars = c('R1', 'R2', 'R3')) }
  df = data.table(df)
  setnames(df, 'Design.Type', 'Design Type')
  df$value = 10^df$value
  return(df)
}


w = 9.5
h = 6

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(2)
cols = c(cols, '#aaaaaa')


dat = fread('Validation/validation_output_fs_means_corrected.csv'); dat$V1 = NULL
dat$SEM.scale = 1/sqrt(3)
dat$SEM.scale[dat$Name %in% c('GPD', 'CYC1', 'TEF', 'ADH1', 'PGK1', 'TPI1')] = 1/sqrt(6)
dat$VYB_Mean_A = apply(cbind(dat$R1_FALSE, dat$R2_FALSE, dat$R3_FALSE), 1, mean)
dat$VYB_SD_A = apply(cbind(dat$R1_FALSE, dat$R2_FALSE, dat$R3_FALSE), 1, sd)
dat$VYB_Mean_B = apply(cbind(dat$R1_TRUE, dat$R2_TRUE, dat$R3_TRUE), 1, mean)
dat$VYB_SD_B = apply(cbind(dat$R1_TRUE, dat$R2_TRUE, dat$R3_TRUE), 1, sd)
dat$VYB_AR = dat$VYB_Mean_B - dat$VYB_Mean_A
dat$VYB_AR_SD = sqrt(dat$VYB_SD_A^2 + dat$VYB_SD_B^2)
dat$bar.min.a = dat$VYB_Mean_A-(dat$VYB_SD_A*dat$SEM.scale)
dat$bar.max.a = dat$VYB_Mean_A+(dat$VYB_SD_A*dat$SEM.scale)
dat$bar.min.b = dat$VYB_Mean_B-(dat$VYB_SD_B*dat$SEM.scale)
dat$bar.max.b = dat$VYB_Mean_B+(dat$VYB_SD_B*dat$SEM.scale)
dat$bar.min.ar = dat$VYB_AR-(dat$VYB_AR_SD*dat$SEM.scale)
dat$bar.max.ar = dat$VYB_AR+(dat$VYB_AR_SD*dat$SEM.scale)

# Did we choose this sequence because we expect it to be unreliable?
dat$trusted_a = dat$Scale_A == 0
dat$trusted_b = dat$Scale_B == 0

# for building final sequence-strength table: controls
#write.csv(dat[1:9,], 'C:/Users/Ben/Dropbox/Lab/Lab Deliverables/Promoter Manuscript/single_seqs_ctrls.csv')

# 4A
dat.use.a = data.frame(FS = dat$Means_A, VYB = dat$VYB_Mean_A, ymin = dat$bar.min.a, ymax = dat$bar.max.a, 
                       Promoter = dat$Promoter, rep = 'A', Offscale = !(dat$trusted_a))

# don't double-count pGPD promoters
dat.use.b = data.frame(FS = dat$Means_B, VYB = dat$VYB_Mean_B, ymin = dat$bar.min.b, ymax = dat$bar.max.b,
                       Promoter = dat$Promoter, rep = 'B', Offscale = !(dat$trusted_b))
dat.use.b = dat.use.b[dat.use.b$Promoter == 'ZEV',]

dat.use = rbind(dat.use.a, dat.use.b)
dat.use = dat.use[apply(dat.use, 1, function(x) all(!is.na(x))),]

dat.use$Type = 'pGPD'
dat.use$Type[dat.use$Promoter == 'ZEV' & dat.use$rep == 'A'] = 'pZEV - uninduced'
dat.use$Type[dat.use$Promoter == 'ZEV' & dat.use$rep == 'B'] = 'pZEV - induced'

dat.use = dat.use[!dat.use$Offscale,]
#mod = lm(dat.use$VYB[!dat.use$Offscale] ~ dat.use$FS[!dat.use$Offscale])
mod = lm(dat.use$VYB ~ dat.use$FS)
summary(mod) # R2 = 0.92
coefs = mod$coef

val.A.low = min(dat.use.a$FS[!is.na(dat.use.a$FS)])
val.B.low = min(dat.use.b$FS[!is.na(dat.use.b$FS)])
val.B.hi = max(dat.use.b$FS[!is.na(dat.use.b$FS)])
vals.offscale = c(val.A.low, val.B.low, val.B.hi)

png(filename = 'D:/Promoter Design Data/Figures/PNGs/4A.png',
    units = 'cm', width = w, height = h, res = 600)

#p = ggplot(dat.use, aes(x = FS, y = VYB, color = Type, shape = Offscale, ymin = ymin, ymax = ymax)) + 
p = ggplot(dat.use, aes(x = FS, y = VYB, color = Type, ymin = ymin, ymax = ymax)) + 
  geom_point(size = 1, stroke = 0) + geom_errorbar(size = 0.3) + theme_bw() + 
  geom_abline(slope = coefs[2], intercept = coefs[1], lty = 2) +
  theme(plot.title = element_text(hjust = 0.4, size = 8, face='bold')) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold'),
        legend.text = element_text(size=6), legend.title = element_text(size=8, face='bold'),
        legend.key.size = unit(0.2,'cm')) +
  scale_x_continuous(breaks=c(-1,0,1,2), limits = c(-1.01, 2.01)) + 
  labs(x = 'FACS-seq\nactivity measurements (log10)', y = 'Individual testing\nactivity measurements (log10)')
print(p)
dev.off()

# 4B
# break log scaling
dat = fread('Validation/validation_output_fs_means_corrected.csv'); dat$V1 = NULL
dat$SEM.scale = 1/sqrt(3)
dat$SEM.scale[dat$Name %in% c('GPD', 'CYC1', 'TEF', 'ADH1', 'PGK1', 'TPI1')] = 1/sqrt(6)
dat$VYB_Mean_A = apply(cbind(dat$R1_FALSE, dat$R2_FALSE, dat$R3_FALSE), 1, mean)
dat$VYB_SD_A = apply(cbind(dat$R1_FALSE, dat$R2_FALSE, dat$R3_FALSE), 1, sd)
dat$VYB_Mean_B = apply(cbind(dat$R1_TRUE, dat$R2_TRUE, dat$R3_TRUE), 1, mean)
dat$VYB_SD_B = apply(cbind(dat$R1_TRUE, dat$R2_TRUE, dat$R3_TRUE), 1, sd)
dat$VYB_AR = dat$VYB_Mean_B - dat$VYB_Mean_A
dat$VYB_AR_SD = sqrt(dat$VYB_SD_A^2 + dat$VYB_SD_B^2)
dat$bar.min.a = 10^(dat$VYB_Mean_A-dat$VYB_SD_A*dat$SEM.scale)
dat$bar.max.a = 10^(dat$VYB_Mean_A+dat$VYB_SD_A*dat$SEM.scale)
dat$bar.min.b = 10^(dat$VYB_Mean_B-dat$VYB_SD_B*dat$SEM.scale)
dat$bar.max.b = 10^(dat$VYB_Mean_B+dat$VYB_SD_B*dat$SEM.scale)
dat$bar.min.ar = 10^(dat$VYB_AR-dat$VYB_AR_SD*dat$SEM.scale)
dat$bar.max.ar = 10^(dat$VYB_AR+dat$VYB_AR_SD*dat$SEM.scale)

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
dat.gpd = dat.gpd[order(dat.gpd$VYB_Mean_A),]
dat.gpd = dat.gpd[order(dat.gpd$`Design Type`),]
# drop superfluous rows
dat.gpd = dat.gpd[c(1:3, 4,6,8,11,13,14:18),]

dat.gpd$`Design Type` = factor(dat.gpd$`Design Type`, 
                               levels = levels(factor(dat.gpd$`Design Type`))[c(2,3,1)] )

dat.gpd.me = melt.for.jitter(dat.gpd[dat.gpd$`Design Type` != 'Control',], 'Uninduced')
dat.gpd.mc = melt.for.jitter(dat.gpd[dat.gpd$`Design Type` == 'Control',], 'All')
dat.gpd.m = rbind(dat.gpd.me, dat.gpd.mc)

#dat.gpd$xstr = c(paste0('pGPD-', 1:15), 'pCYC1', 'pTEF1', 'pGPD')
xstr = c(paste0('pGPD-', c(1,3,5,8,10,11:15)), 'pCYC1', 'pTEF1', 'pGPD')

png(filename = 'D:/Promoter Design Data/Figures/PNGs/4B.png',
    units = 'cm', width = w, height = h + 0.3, res = 600)

p = ggplot(dat.gpd, aes(Name, VYB_Mean_A, color = `Design Type`, fill = `Design Type`,
                        ymin = bar.min.a, ymax = bar.max.a)) + 
  geom_errorbar(width = 0.6) + 
  geom_col(alpha = 0.5, width = 0.6) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.4, size = 8, face='bold')) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold'),
        legend.text = element_text(size=6), legend.title = element_text(size=8, face='bold'),
        legend.key.size = unit(0.2,'cm')) +
  theme(axis.text.x=element_text(angle=270, hjust = 0),
        axis.ticks.x=element_blank()) +
  scale_x_discrete(labels=xstr) +
  ylim(0, 25) +
  labs(x = '', y='Promoter activity\npGPD designs', color = 'Design type', fill = 'Design type') + 
  guides(colour = guide_legend(override.aes = list(alpha = 1)))  +
  scale_fill_manual(values=cols) + 
  scale_color_manual(values=cols)
q = p + geom_jitter(data=dat.gpd.m, aes(x=Name, y=value, color = `Design Type`,
                                        ymin=NULL, ymax=NULL), width = 0.1, height = 0, size = 0.5)
print(q)
dev.off()

# 4C
png(filename = 'D:/Promoter Design Data/Figures/PNGs/4C.png',
    units = 'cm', width = w, height = h + 0.3, res = 600)
dat.zev.i = dat[dat$Name %in% c('ZEV_Pr3', 'ZEV_Pr4', 'ZEV_Pr8') | 
                  dat$Experiment.designed == 34,] #c(34,39)

dat.zev.i = dat.zev.i[(dat.zev.i$trusted_a & dat.zev.i$trusted_b) | is.na(dat.zev.i$Experiment.designed),]
dat.zev.i$`Design Type` = dat.zev.i$Experiment.designed
dat.zev.i$`Design Type`[is.na(dat.zev.i$`Design Type`)] = 'Control'
dat.zev.i$`Design Type`[dat.zev.i$`Design Type` == 34] = 'Evolution-GC'
#dat.zev.i$`Design Type`[dat.zev.i$`Design Type` == 39] = 'Gradient*'
dat.zev.i$`Design Type` = factor(dat.zev.i$`Design Type`)

dat.zev.i = dat.zev.i[order(dat.zev.i$VYB_Mean_B),]
dat.zev.i = dat.zev.i[order(dat.zev.i$`Design Type`),]
dat.zev.i.all = dat.zev.i

# drop superfluous rows
dat.zev.i = dat.zev.i[c(1:3, 4,6,8,10,12,14),]
xstr = c(paste0('pZEV-I-', c(1,3,5,7,9,11)), 'P8', 'P4', 'P3')


dat.zev.i$Name = factor(dat.zev.i$Name, levels = dat.zev.i$Name[order(dat.zev.i$Experiment.designed, dat.zev.i$VYB_Mean_B)])

# for building final sequence-strength table: pZEV-I
dat.zev.i$Seqs.final = sapply(dat.zev.i$Seqs, function(x) substr(x,59,nchar(x)-59))
#write.csv(dat.zev.i, 'C:/Users/Ben/Dropbox/Lab/Lab Deliverables/Promoter Manuscript/single_seqs_pZEV_I.csv')
dat.zev.i$`Design Type` = factor(dat.zev.i$`Design Type`, 
                               levels = levels(factor(dat.zev.i$`Design Type`))[c(2,1)] )

dat.zev.i.m = melt.for.jitter(dat.zev.i, 'Induced')
p = ggplot(dat.zev.i, aes(Name, VYB_Mean_B, color = `Design Type`, fill = `Design Type`, 
                      ymin = bar.min.b, ymax = bar.max.b)) + 
  geom_errorbar(width = 0.6) + 
  geom_col(alpha = 0.5, width = 0.6) + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.4, size = 8, face='bold')) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold'),
        legend.text = element_text(size=6), legend.title = element_text(size=8, face='bold'),
        legend.key.size = unit(0.2,'cm')) +
  theme(axis.text.x=element_text(angle=270, hjust = 0),
        axis.ticks.x=element_blank()) +
  scale_x_discrete(labels=xstr) +
  labs(x = '', y='Induced promoter activity \npZEV-induced designs', color = 'Design type', fill = 'Design type') + 
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  scale_fill_manual(values=cols[c(1,3)]) + 
  scale_color_manual(values=cols[c(1,3)])
q = p + geom_jitter(data=dat.zev.i.m, aes(x=Name, y=value, color = `Design Type`,
                                        ymin=NULL, ymax=NULL), width = 0.1, height = 0, size = 0.5)
print(q)
dev.off()

# 4D: ARs
dat.zev.ar = dat[(is.na(dat$Seqs) & dat$Name %in% c('ZEV_Pr3', 'ZEV_Pr4', 'ZEV_Pr8')) | 
                   dat$Experiment.designed == 44,]# | dat$ZEV_AR_offscale_cand | dat$ZEV_I_AR_cand,]
xstr = c(paste0('pZEV-AR-', 1:10), 'P8', 'P4', 'P3')
dat.zev.ar$Name = factor(dat.zev.ar$Name, levels = dat.zev.ar$Name[order(dat.zev.ar$Experiment.designed, dat.zev.ar$VYB_AR)])
dat.zev.ar$`Design Type` = dat.zev.ar$Experiment.designed
dat.zev.ar$`Design Type`[is.na(dat.zev.ar$`Design Type`)] = 'Control'
dat.zev.ar$`Design Type`[dat.zev.ar$Experiment.designed == 44] = 'Evolution-GC'
dat.zev.ar$`Design Type` = factor(dat.zev.ar$`Design Type`)

# for building final sequence-strength table: pZEV-AR
dat.zev.ar$Seqs.final = sapply(dat.zev.ar$Seqs, function(x) substr(x,59,nchar(x)-59))

dat.zev.ar$`Design Type` = factor(dat.zev.ar$`Design Type`, 
                                 levels = levels(factor(dat.zev.ar$`Design Type`))[c(2,1)] )

dat.zev.ar.m = melt.for.jitter(dat.zev.ar, 'AR')
#write.csv(dat.zev.ar, 'C:/Users/Ben/Dropbox/Lab/Lab Deliverables/Promoter Manuscript/single_seqs_pZEV_AR.csv')

png(filename = 'D:/Promoter Design Data/Figures/PNGs/4D.png',
    units = 'cm', width = w, height = h, res = 600)

p = ggplot(dat.zev.ar, aes(Name, VYB_AR, fill = `Design Type`, color = `Design Type`,
                           ymin = bar.min.ar, ymax = bar.max.ar)) + 
  theme_bw() +
  geom_errorbar(width = 0.6) + 
  geom_col(alpha = 0.5, width = 0.6) +
  theme(plot.title = element_text(hjust = 0.4, size = 8, face='bold')) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold'),
        legend.text = element_text(size=6), legend.title = element_text(size=8, face='bold'),
        legend.key.size = unit(0.2,'cm')) +
  theme(axis.text.x=element_text(angle=270, hjust = 0),
        axis.ticks.x=element_blank()) +
  scale_x_discrete(labels=xstr) +
  labs(x = '', y='Activation ratio \npZEV-Activation ratio designs', color = 'Design type', fill = 'Design type') +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  scale_fill_manual(values=cols[c(1,3)]) + 
  scale_color_manual(values=cols[c(1,3)])
q = p + geom_jitter(data=dat.zev.ar.m, aes(x=Name, y=value, color = `Design Type`,
                                          ymin=NULL, ymax=NULL), width = 0.1, height = 0, size = 0.5)
print(q)
dev.off()

# NEW 4E: "4D_new_redesign.R"

### Information for text ###
# 95% confidence intervals, with the correct log-normality
get.cids = function(df, n = 3) {
  df$CI_A_low = 10^(log10(df$VYB_Mean_A) - qt(0.975,df=n-1)*df$VYB_SD_A/sqrt(n))
  df$CI_A_hi = 10^(log10(df$VYB_Mean_A) + qt(0.975,df=n-1)*df$VYB_SD_A/sqrt(n))
  df$CI_B_low = 10^(log10(df$VYB_Mean_B) - qt(0.975,df=n-1)*df$VYB_SD_B/sqrt(n))
  df$CI_B_hi = 10^(log10(df$VYB_Mean_B) + qt(0.975,df=n-1)*df$VYB_SD_B/sqrt(n))
  return(df)
}

# for constant controls, we have 6 measurements: use all of them
get.cids.ctrl = function(df, n = 6) {
  tmp = cbind(df$R1_FALSE, df$R1_TRUE, df$R2_FALSE, df$R2_TRUE, df$R3_FALSE, df$R3_TRUE)
  m = apply(tmp,1,mean)
  s = apply(tmp,1,sd)
  df$CI_low = 10^(m - qt(0.975,df=n-1)*s/sqrt(n))
  df$CI_hi = 10^(m + qt(0.975,df=n-1)*s/sqrt(n))
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
  print(i)
  print(10^mean(tmp[i,]))
  print(t.test(tmp[i,], g, 'greater')$p.value)
  print(t.test(tmp[i,], g, 'less')$p.value)
  print(t.test(tmp[i,], g, 'two.sided'))
}

# Designed pZEV-I
dat.tmp = dat.zev.i[dat.zev.i$Experiment.designed == 34,]
dat.tmp = dat.tmp[dat.tmp$VYB_Mean_B %in% range(dat.tmp$VYB_Mean_B),]
dat.tmp = get.cids(dat.tmp)
print(dat.tmp$VYB_Mean_B)
print(dat.tmp$CI_B_low)
print(dat.tmp$CI_B_hi)

# p-values vs. P4 and P8
dat.tmp = dat.zev.i.all[dat.zev.i.all$Experiment.designed == 34,]
dat.ctrl = dat[dat$Name %in% c('ZEV_Pr3', 'ZEV_Pr4', 'ZEV_Pr8'),]
z = cbind(dat.ctrl$R1_TRUE, dat.ctrl$R2_TRUE, dat.ctrl$R3_TRUE)
tmp = cbind(dat.tmp$R1_TRUE, dat.tmp$R2_TRUE, dat.tmp$R3_TRUE)
pvals = matrix(nrow = nrow(z), ncol = nrow(tmp))
for(i in 1:nrow(z)) { for(j in 1:nrow(tmp)) {
  pvals[i,j] = t.test(z[i,], tmp[j,], 'two.sided')$p.value
}}
colnames(pvals) = dat.tmp$VYB_Mean_B; rownames(pvals) = dat.ctrl$Name

apply(pvals, 1, max)
# ZEV_Pr3      ZEV_Pr4      ZEV_Pr8 
# 9.932181e-01 2.863542e-02 9.995833e-05 

# CI for Pr3
dat.ctrl = get.cids(dat.ctrl)
print(dat.ctrl$Name)
print(dat.ctrl$VYB_Mean_B)
print(dat.ctrl$CI_B_low)
print(dat.ctrl$CI_B_hi)

# Designed pZEV-AR
get.cids.ar = function(df, n = 3) {
  tmp = cbind(df$R1_TRUE - df$R1_FALSE, df$R2_TRUE - df$R2_FALSE, df$R3_TRUE - df$R3_FALSE)
  m = apply(tmp,1,mean)
  s = apply(tmp,1,sd)
  df$CI_AR_low = 10^(m - qt(0.975, n-1)*s/sqrt(n)) # qnorm(0.975) is about 1.96
  df$CI_AR_hi = 10^(m + qnorm(0.975, n-1)*s/sqrt(n))
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
  pvals[i,j] = t.test(z[i,], tmp[j,], 'two.sided')$p.value
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

# lastly, get the CI table that has the raw numbers for Supplementary Table 3
dat = fread('Validation/validation_output_fs_means_corrected.csv'); dat$V1 = NULL
dat$VYB_Mean_A = apply(cbind(dat$R1_FALSE, dat$R2_FALSE, dat$R3_FALSE), 1, mean)
dat$VYB_SD_A = apply(cbind(dat$R1_FALSE, dat$R2_FALSE, dat$R3_FALSE), 1, sd)
dat$VYB_Mean_B = apply(cbind(dat$R1_TRUE, dat$R2_TRUE, dat$R3_TRUE), 1, mean)
dat$VYB_SD_B = apply(cbind(dat$R1_TRUE, dat$R2_TRUE, dat$R3_TRUE), 1, sd)
dat$VYB_AR = dat$VYB_Mean_B - dat$VYB_Mean_A
dat$VYB_Mean_A = 10^(dat$VYB_Mean_A)
dat$VYB_Mean_B = 10^(dat$VYB_Mean_B)
dat$VYB_AR = 10^(dat$VYB_AR)
# create columns CI_A_low, CI_A_hi, CI_B_low, CI_B_hi
dat = get.cids(dat)
# create columns CI_AR_low, CI_AR_hi
dat = get.cids.ar(dat)
# special handling for constant controls: take advantage of the extra replicates
# create columns CI_low, CI_hi
dat = get.cids.ctrl(dat)
write.csv(dat, 'D:/Promoter Design Data/Validation/validation_data_with_CIs.csv')


