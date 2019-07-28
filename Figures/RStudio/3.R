library(data.table)
library(ggplot2)

setwd('D:/Promoter Design Data/Preds')

dat.orig = fread('preds_GPD.csv')
setkey(dat.orig, 'Seqs')

get.r2s = function(str.A, str.B, pred.str.A, pred.str.B) {
  ans = rep(0, 4)
  means.avg = apply(cbind(str.A, str.B), 1, mean)
  pred.means.avg = apply(cbind(pred.str.A, pred.str.B), 1, mean)
  ar = str.B - str.A
  pred.ar = pred.str.B - pred.str.A
  dat = data.table(str.A = str.A, str.B = str.B, pred.str.A = pred.str.A, pred.str.B = pred.str.B, is.gpd = str.A == str.B,
                   means.avg = means.avg, pred.means.avg = pred.means.avg, pred.ar = pred.ar, ar = ar)
  
  gpd = dat[dat$is.gpd,]
  zev = dat[!(dat$is.gpd),]
  ans[1] = summary(lm(pred.means.avg ~ means.avg, data = gpd))$r.squared
  ans[2] = summary(lm(pred.str.A ~ str.A, data = zev))$r.squared
  ans[3] = summary(lm(pred.str.B ~ str.B, data = zev))$r.squared
  ans[4] = summary(lm(pred.ar ~ ar, data = zev))$r.squared
  return(ans)
}

# test sets are different between original pGPD and merged datasets
setwd('D:/Promoter Design Data/Preds/Ensemble Preds/')
fns = dir()
r2s = matrix(nrow = 4, ncol = length(fns))
nr = dim(fread(fns[1]))[1]
ans.a = matrix(nrow = nr, ncol = length(fns))
ans.b = matrix(nrow = nr, ncol = length(fns))
for(i in 1:length(fns)) {
  tmp = fread(fns[i])
  if(i == 1) {
    sa = tmp$Strength_A; sb = tmp$Strength_B; seqs = tmp$Seqs
  }
  ans.a[,i] = tmp$Pred_Strength_A
  ans.b[,i] = tmp$Pred_Strength_B
  r2s[,i] = get.r2s(sa, sb, tmp$Pred_Strength_A, tmp$Pred_Strength_B)
}
merged.a = apply(ans.a,1,mean)
merged.b = apply(ans.b,1,mean)
dat.mer = data.table(Pred_Strength_A = merged.a, Pred_Strength_B = merged.b, 
                     Strength_A = sa, Strength_B = sb, Seqs = seqs)
dat.mer$Library = 'ZEV'
dat.mer$Library[dat.mer$Strength_A == dat.mer$Strength_B] = 'GPD'

setwd('D:/Promoter Design Data/Preds/')
write.csv(dat.mer, 'preds_joined_all.csv')

# for text

# medians
meds = apply(r2s, 1, median)
names(meds) = c('med.GPD', 'med.ZEV-U', 'med.ZEV-I', 'med.ZEV-AR')
print(meds)

# GPD
dat.mer$Pred_Strength_Mer = apply(cbind(dat.mer$Pred_Strength_A, dat.mer$Pred_Strength_B), 1, mean)
dat.mer$Strength_Mer = apply(cbind(dat.mer$Strength_A, dat.mer$Strength_B), 1, mean)
summary(lm(Pred_Strength_Mer ~ Strength_Mer, data = dat.mer[dat.mer$Library == 'GPD',])) # 0.80
summary(lm(Pred_Strength_A ~ Strength_A, data = dat.mer[dat.mer$Library == 'ZEV',])) # 0.84
summary(lm(Pred_Strength_B ~ Strength_B, data = dat.mer[dat.mer$Library == 'ZEV',])) # 0.79
dat.mer$AR = dat.mer$Strength_B - dat.mer$Strength_A
dat.mer$Pred_AR = dat.mer$Pred_Strength_B - dat.mer$Pred_Strength_A
summary(lm(Pred_AR ~ AR, data = dat.mer[dat.mer$Library == 'ZEV',])) # 0.82

# shared parameters for figure properties
fig.w = 8
fig.h = 6
ax.text = 6
ax.titl = 8
lg.text = 6
lg.titl = 8
lg.size = 0.2

# 3A
dat.mer =  fread('preds_joined_all.csv'); dat.mer = dat.mer[dat.mer$Library == 'GPD',]
dat.mer$Pred_Strength_Mer = apply(cbind(dat.mer$Pred_Strength_A, dat.mer$Pred_Strength_B), 1, mean)
setkey(dat.mer, 'Seqs')

# 12123 merged sequences
dat = merge(dat.orig, dat.mer)
dat$Orig_Err = (dat$Pred_Strength - dat$Strength)^2
dat$Merged_Err = (dat$Pred_Strength_Mer - dat$Strength)^2
dat$merged_better = abs(dat$Pred_Strength_Mer - dat$Strength) < abs(dat$Pred_Strength - dat$Strength)
dat$merged_better_txt = 'Original'
dat$merged_better_txt[dat$merged_better] = 'Merged'

png(filename = 'D:/Promoter Design Data/Figures/PNGs/3A.png',
    units = 'cm', width = fig.w, height = fig.h, res = 600)
p = ggplot(dat, aes(x = abs(dat$Pred_Strength - dat$Strength), y = abs(dat$Pred_Strength_Mer - dat$Strength), color = merged_better_txt)) +
  geom_point(alpha = 0.2) + #scale_color_viridis() + 
  theme_bw() + #xlim(c(0, 0.05)) + ylim(c(0, 0.05))
  geom_abline(intercept = 0, slope = 1, lty = 2) + 
  theme(axis.text = element_text(size=ax.text), axis.title = element_text(size=ax.titl, face='bold'),
        legend.text = element_text(size=lg.text), legend.title = element_text(size=lg.titl, face='bold'),
        legend.key.size = unit(lg.size,'cm')) +
  labs(x = 'Original Residuals (pGPD)', y='Final Residuals (pGPD)', color = 'Best Model')
print(p)
dev.off()

# 3B
setwd('D:/Promoter Design Data/Preds')
dat.orig = fread('preds_ZEV.csv')
setkey(dat.orig, 'Seqs')
dat.orig$AR = dat.orig$Induced - dat.orig$Uninduced
dat.orig$Pred_AR = dat.orig$Pred_Induced - dat.orig$Pred_Uninduced

setwd('D:/Promoter Design Data/Preds')
dat.mer =  fread('preds_joined_all.csv'); dat.mer = dat.mer[dat.mer$Library == 'ZEV',]
setnames(dat.mer, c('Pred_Strength_A', 'Pred_Strength_B'), c('Pred_Uninduced_Mer', 'Pred_Induced_Mer'))
dat.mer$Pred_AR_Mer = dat.mer$Pred_Induced_Mer - dat.mer$Pred_Uninduced_Mer
dat.mer$Seqs = sapply(dat.mer$Seqs, function(x) substr(x,34, nchar(x) - 33))
setkey(dat.mer, 'Seqs')
dat = merge(dat.orig, dat.mer)

dat$merged_better_uninduced = abs(dat$Pred_Uninduced_Mer - dat$Uninduced) < abs(dat$Pred_Uninduced - dat$Uninduced)
dat$merged_better_induced = abs(dat$Pred_Induced_Mer - dat$Induced) < abs(dat$Pred_Induced - dat$Induced)
dat$merged_better_AR = abs(dat$Pred_AR_Mer - dat$AR) < abs(dat$Pred_AR - dat$AR)

dat$merged_better_txt = 'Original'
dat$merged_better_txt[dat$merged_better_uninduced] = 'Merged'
png(filename = 'D:/Promoter Design Data/Figures/PNGs/3B.png',
    units = 'cm', width = fig.w, height = fig.h, res = 600)
p = ggplot(dat, aes(x = abs(dat$Pred_Uninduced - dat$Uninduced), y = abs(dat$Pred_Uninduced_Mer - dat$Uninduced), color = merged_better_txt)) +
  geom_point(alpha = 0.2) + #scale_color_viridis() + 
  theme_bw() + #xlim(c(0, 0.05)) + ylim(c(0, 0.05))
  geom_abline(intercept = 0, slope = 1, lty = 2) + 
  theme(axis.text = element_text(size=ax.text), axis.title = element_text(size=ax.titl, face='bold'),
        legend.text = element_text(size=lg.text), legend.title = element_text(size=lg.titl, face='bold'),
        legend.key.size = unit(lg.size,'cm')) +
  labs(x = 'Original Residuals\n(pZEV-Uninduced)', y='Final Residuals\n(pZEV-Uninduced)', color = 'Best Model')
print(p)
dev.off()

# 3C
dat$merged_better_txt = 'Original'
dat$merged_better_txt[dat$merged_better_induced] = 'Merged'
png(filename = 'D:/Promoter Design Data/Figures/PNGs/3C.png',
    units = 'cm', width = fig.w, height = fig.h, res = 600)
p = ggplot(dat, aes(x = abs(dat$Pred_Induced - dat$Induced), y = abs(dat$Pred_Induced_Mer - dat$Induced), color = merged_better_txt)) +
  geom_point(alpha = 0.2) + #scale_color_viridis() + 
  theme_bw() + #xlim(c(0, 0.05)) + ylim(c(0, 0.05))
  geom_abline(intercept = 0, slope = 1, lty = 2) + 
  theme(axis.text = element_text(size=ax.text), axis.title = element_text(size=ax.titl, face='bold'),
        legend.text = element_text(size=lg.text), legend.title = element_text(size=lg.titl, face='bold'),
        legend.key.size = unit(lg.size,'cm')) +
  labs(x = 'Original Residuals\n(pZEV-Induced)', y='Final Residuals\n(pZEV-Induced)', color = 'Best Model')
print(p)
dev.off()

# 3D
dat$merged_better_txt = 'Original'
dat$merged_better_txt[dat$merged_better_AR] = 'Merged'
png(filename = 'D:/Promoter Design Data/Figures/PNGs/3D.png',
    units = 'cm', width = fig.w, height = fig.h, res = 600)
p = ggplot(dat, aes(x = abs(dat$Pred_AR - dat$AR), y = abs(dat$Pred_AR_Mer - dat$AR), color = merged_better_txt)) +
  geom_point(alpha = 0.2) + #scale_color_viridis() + 
  theme_bw() + #xlim(c(0, 0.05)) + ylim(c(0, 0.05))
  geom_abline(intercept = 0, slope = 1, lty = 2) + 
  theme(axis.text = element_text(size=ax.text), axis.title = element_text(size=ax.titl, face='bold'),
        legend.text = element_text(size=lg.text), legend.title = element_text(size=lg.titl, face='bold'),
        legend.key.size = unit(lg.size,'cm')) +
  labs(x = 'Original Residuals\n(pZEV-Activation Ratio)', y='Final Residuals\n(pZEV-Activation Ratio)', color = 'Best Model')
print(p)
dev.off()




