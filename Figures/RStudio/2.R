# Fig. 2: plots for GPD and ZEV separate modeling
library(data.table)
library(ggplot2)
library(MASS)
library(viridis)
setwd('D:/Promoter Design Data')

#get_density = function(x, y, n = 200) {
get_density = function(x, y, n = 100) {
  dens = MASS::kde2d(x = x, y = y, n = n)
  ix = findInterval(x, dens$x)
  iy = findInterval(y, dens$y)
  ii = cbind(ix, iy)
  return(dens$z[ii])
}

# shared parameters for figure properties
fig.w = 8
fig.h = 5.5
ax.text = 6
ax.titl = 8
lg.text = 6
lg.titl = 8
lg.size = 0.2

# 2A
setwd('D:/Promoter Design Data/Preds/')
dat.mer =  fread('preds_joined_all.csv'); dat.mer.use = dat.mer[dat.mer$Library == 'GPD',]
dat.mer.use$Pred_Strength_Mer = apply(cbind(dat.mer.use$Pred_Strength_A, dat.mer.use$Pred_Strength_B), 1, mean)
setkey(dat.mer.use, 'Seqs')

setnames(dat.mer.use, c('Strength_A', 'Pred_Strength_Mer'), c('Strength.A', 'Strength.B'))
dat.mer.use$Density = get_density(dat.mer.use$Strength.A, dat.mer.use$Strength.B)
dat.mer.use$`Cell Count` = dat.mer.use$Density

png(filename = 'D:/Promoter Design Data/Figures/PNGs/2A.png',
    units = 'cm', width = fig.w, height = fig.h, res = 600)

# R2 = 0.80
print(summary(lm(Strength.B ~ Strength.A, data = dat.mer.use)))
p = ggplot(dat.mer.use, aes(x = Strength.A, y = Strength.B, color = `Cell Count`)) + 
  geom_point(alpha = 0.1, size = 0.2, stroke = 0, shape = 16) + theme_bw() + 
  theme(axis.text = element_text(size=ax.text), axis.title = element_text(size=ax.titl, face='bold'),
        legend.text = element_text(size=lg.text), legend.title = element_text(size=lg.titl, face='bold'),
        legend.key.size = unit(lg.size,'cm')) + 
  labs(x = 'Activity\n(pGPD)', y='Predicted activity\n(pGPD)', color = 'Density') +
  xlim(-0.55, 0.55) + ylim(-0.55, 0.55) + 
  scale_color_viridis()
print(p)
dev.off()

# 2B
setwd('D:/Promoter Design Data/Preds/')
dat.mer =  fread('preds_joined_all.csv'); dat.mer.use = dat.mer[dat.mer$Library == 'ZEV',]
setkey(dat.mer.use, 'Seqs')

setnames(dat.mer.use, c('Strength_A', 'Pred_Strength_A'), c('Strength.A', 'Strength.B'))
dat.mer.use$Density = get_density(dat.mer.use$Strength.A, dat.mer.use$Strength.B)
dat.mer.use$`Cell Count` = dat.mer.use$Density

png(filename = 'D:/Promoter Design Data/Figures/PNGs/2B.png',
    units = 'cm', width = fig.w, height = fig.h, res = 600)

# R2 = 0.84
print(summary(lm(Strength.B ~ Strength.A, data = dat.mer.use)))
p = ggplot(dat.mer.use, aes(x = Strength.A, y = Strength.B, color = `Cell Count`)) + 
  geom_point(alpha = 0.1, size = 0.2, stroke = 0, shape = 16) + theme_bw() + 
  theme(axis.text = element_text(size=ax.text), axis.title = element_text(size=ax.titl, face='bold'),
        legend.text = element_text(size=lg.text), legend.title = element_text(size=lg.titl, face='bold'),
        legend.key.size = unit(lg.size,'cm')) + 
  labs(x = 'Activity\n(pZEV-uninduced)', y='Predicted activity\n(pZEV-uninduced)', color = 'Density') +
  xlim(-0.60, 0.60) + ylim(-0.60, 0.60) + 
  scale_color_viridis()
print(p)
dev.off()

# 2C
setwd('D:/Promoter Design Data/Preds/')
dat.mer =  fread('preds_joined_all.csv'); dat.mer.use = dat.mer[dat.mer$Library == 'ZEV',]
setkey(dat.mer.use, 'Seqs')

setnames(dat.mer.use, c('Strength_B', 'Pred_Strength_B'), c('Strength.A', 'Strength.B'))
dat.mer.use$Density = get_density(dat.mer.use$Strength.A, dat.mer.use$Strength.B, 200)
dat.mer.use$`Cell Count` = dat.mer.use$Density

png(filename = 'D:/Promoter Design Data/Figures/PNGs/2C.png',
    units = 'cm', width = fig.w, height = fig.h, res = 600)

# R2 = 0.79
print(summary(lm(Strength.B ~ Strength.A, data = dat.mer.use)))
p = ggplot(dat.mer.use, aes(x = Strength.A, y = Strength.B, color = `Cell Count`)) + 
  geom_point(alpha = 0.1, size = 0.2, stroke = 0, shape = 16) + theme_bw() + 
  theme(axis.text = element_text(size=ax.text), axis.title = element_text(size=ax.titl, face='bold'),
        legend.text = element_text(size=lg.text), legend.title = element_text(size=lg.titl, face='bold'),
        legend.key.size = unit(lg.size,'cm')) + 
  labs(x = 'Activity\n(pZEV-induced)', y='Predicted activity\n(pZEV-induced)', color = 'Density') +
  xlim(1.0, 1.5) + ylim(1.0, 1.5) + 
  scale_color_viridis()
print(p)
dev.off()

# 2D
setwd('D:/Promoter Design Data/Preds/')
dat.mer =  fread('preds_joined_all.csv'); dat.mer.use = dat.mer[dat.mer$Library == 'ZEV',]
setkey(dat.mer.use, 'Seqs')
dat.mer.use$Strength.A = dat.mer.use$Strength_B - dat.mer.use$Strength_A
dat.mer.use$Strength.B = dat.mer.use$Pred_Strength_B - dat.mer.use$Pred_Strength_A

dat.mer.use$Density = get_density(dat.mer.use$Strength.A, dat.mer.use$Strength.B, 200)
dat.mer.use$`Cell Count` = dat.mer.use$Density

png(filename = 'D:/Promoter Design Data/Figures/PNGs/2D.png',
    units = 'cm', width = fig.w, height = fig.h, res = 600)

# R2 = 0.79
print(summary(lm(Strength.B ~ Strength.A, data = dat.mer.use)))
p = ggplot(dat.mer.use, aes(x = Strength.A, y = Strength.B, color = `Cell Count`)) + 
  geom_point(alpha = 0.1, size = 0.2, stroke = 0, shape = 16) + theme_bw() + 
  theme(axis.text = element_text(size=ax.text), axis.title = element_text(size=ax.titl, face='bold'),
        legend.text = element_text(size=lg.text), legend.title = element_text(size=lg.titl, face='bold'),
        legend.key.size = unit(lg.size,'cm')) + 
  labs(x = 'Activation ratio\n(pZEV)', y='Predicted activation\nratio (pZEV)', color = 'Density') +
  xlim(0.8, 2.0) + ylim(0.8, 2.0) + 
  scale_color_viridis()
print(p)
dev.off()

