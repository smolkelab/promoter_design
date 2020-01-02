# Fig. 2: plots for GPD and ZEV separate modeling
library(data.table)
library(ggplot2)
library(MASS)
library(viridis)
setwd('D:/Promoter Design Data')

get_density = function(x, y, n = 200) {
  dens = MASS::kde2d(x = x, y = y, n = n)
  ix = findInterval(x, dens$x)
  iy = findInterval(y, dens$y)
  ii = cbind(ix, iy)
  return(dens$z[ii])
}

process.raw.log = function(fn) {
  lines = suppressWarnings(readLines(fn))
  lines = lines[(which(lines == 'None')+1):length(lines)]
  lines = lines[seq(from = 2, to = length(lines), by = 2)]
  lines = lines[lines != 'Predictions written']
  lines = sapply(lines, function(x) strsplit(x, 'loss:')[[1]])
  lines = t(lines[2:3,])
  lines[,1] = sapply(lines[,1], function(x) strsplit(x, '- val_')[[1]][1])
  ans = matrix(nrow = nrow(lines), ncol = ncol(lines))
  for(i in 1:ncol(ans)) { ans[,i] = as.numeric(lines[,i]) }
  colnames(ans) = c('Training', 'Validation'); return(ans)
}

# shared parameters for figure properties
fig.w = 6.5
fig.h = 4.5
ax.text = 6
ax.titl = 8
lg.text = 6
lg.titl = 8
lg.size = 0.2

# S8A
setwd('D:/Promoter Design Data/Preds')
dat.mer.use = fread('preds_GPD.csv')
setnames(dat.mer.use, c('Strength', 'Pred_Strength'), c('Strength.A', 'Strength.B'))
dat.mer.use$Density = get_density(dat.mer.use$Strength.A, dat.mer.use$Strength.B)
dat.mer.use$`Cell Count` = dat.mer.use$Density

png(filename = 'D:/Promoter Design Data/Figures/PNGs_new/S8A.png',
    units = 'cm', width = fig.w, height = fig.h, res = 600)

# R2 = 0.65
print(summary(lm(Strength.B ~ Strength.A, data = dat.mer.use)))
p = ggplot(dat.mer.use, aes(x = Strength.A, y = Strength.B, color = `Cell Count`)) + 
  geom_point(alpha = 0.1, size = 0.2, stroke = 0, shape = 16) + theme_bw() + 
  theme(axis.text = element_text(size=ax.text), axis.title = element_text(size=ax.titl, face='bold'),
        legend.text = element_text(size=lg.text), legend.title = element_text(size=lg.titl, face='bold'),
        legend.key.size = unit(lg.size,'cm')) + 
  labs(x = 'Activity (pGPD)', y='Predicted Activity (pGPD)', color = 'Density') +
  scale_color_viridis()
print(p)
dev.off()

# S8B
setwd('D:/Promoter Design Data/Preds')
dat.mer.zev = fread('preds_ZEV.csv')
dat.mer.zev$AR = dat.mer.zev$Induced - dat.mer.zev$Uninduced
dat.mer.zev$Pred_AR = dat.mer.zev$Pred_Induced - dat.mer.zev$Pred_Uninduced
dat.mer.use = dat.mer.zev
setnames(dat.mer.use, c('Uninduced', 'Pred_Uninduced'), c('Strength.A', 'Strength.B'))
dat.mer.use$Density = get_density(dat.mer.use$Strength.A, dat.mer.use$Strength.B)
dat.mer.use$`Cell Count` = dat.mer.use$Density

png(filename = 'D:/Promoter Design Data/Figures/PNGs_new/S8B.png',
    units = 'cm', width = fig.w, height = fig.h, res = 600)

# R2 = 0.75
print(summary(lm(Strength.B ~ Strength.A, data = dat.mer.use)))
p = ggplot(dat.mer.use, aes(x = Strength.A, y = Strength.B, color = `Cell Count`)) + 
  geom_point(alpha = 0.1, size = 0.2, stroke = 0, shape = 16) + theme_bw() + 
  theme(axis.text = element_text(size=ax.text), axis.title = element_text(size=ax.titl, face='bold'),
        legend.text = element_text(size=lg.text), legend.title = element_text(size=lg.titl, face='bold'),
        legend.key.size = unit(lg.size,'cm')) + 
  labs(x = 'Activity\n(pZEV-Uninduced)', y='Predicted Activity\n(pZEV-Uninduced)', color = 'Density') +
  scale_color_viridis()
print(p)
dev.off()

# S8C
dat.mer.zev = fread('preds_ZEV.csv')
dat.mer.zev$AR = dat.mer.zev$Induced - dat.mer.zev$Uninduced
dat.mer.zev$Pred_AR = dat.mer.zev$Pred_Induced - dat.mer.zev$Pred_Uninduced
dat.mer.use = dat.mer.zev
setnames(dat.mer.use, c('Induced', 'Pred_Induced'), c('Strength.A', 'Strength.B'))
dat.mer.use$Density = get_density(dat.mer.use$Strength.A, dat.mer.use$Strength.B)
dat.mer.use$`Cell Count` = dat.mer.use$Density
png(filename = 'D:/Promoter Design Data/Figures/PNGs_new/S8C.png',
    units = 'cm', width = fig.w, height = fig.h, res = 600)

# R2 = 0.73
print(summary(lm(Strength.B ~ Strength.A, data = dat.mer.use)))
p = ggplot(dat.mer.use, aes(x = Strength.A, y = Strength.B, color = `Cell Count`)) + 
  geom_point(alpha = 0.1, size = 0.2, stroke = 0, shape = 16) + theme_bw() + 
  theme(axis.text = element_text(size=ax.text), axis.title = element_text(size=ax.titl, face='bold'),
        legend.text = element_text(size=lg.text), legend.title = element_text(size=lg.titl, face='bold'),
        legend.key.size = unit(lg.size,'cm')) + 
  labs(x = 'Activity\n(pZEV-Induced)', y='Predicted Activity\n(pZEV-Induced)', color = 'Density') +
  scale_color_viridis()
print(p)
dev.off()

# S8D
dat.mer.zev = fread('preds_ZEV.csv')
dat.mer.zev$AR = dat.mer.zev$Induced - dat.mer.zev$Uninduced
dat.mer.zev$Pred_AR = dat.mer.zev$Pred_Induced - dat.mer.zev$Pred_Uninduced
dat.mer.use = dat.mer.zev
setnames(dat.mer.use, c('AR', 'Pred_AR'), c('Strength.A', 'Strength.B'))
dat.mer.use$Density = get_density(dat.mer.use$Strength.A, dat.mer.use$Strength.B)
dat.mer.use$`Cell Count` = dat.mer.use$Density
png(filename = 'D:/Promoter Design Data/Figures/PNGs_new/S8D.png',
    units = 'cm', width = fig.w, height = fig.h, res = 600)

# R2 = 0.72
print(summary(lm(Strength.B ~ Strength.A, data = dat.mer.use)))
p = ggplot(dat.mer.use, aes(x = Strength.A, y = Strength.B, color = `Cell Count`)) + 
  geom_point(alpha = 0.1, size = 0.2, stroke = 0, shape = 16) + theme_bw() + 
  theme(axis.text = element_text(size=ax.text), axis.title = element_text(size=ax.titl, face='bold'),
        legend.text = element_text(size=lg.text), legend.title = element_text(size=lg.titl, face='bold'),
        legend.key.size = unit(lg.size,'cm')) + 
  labs(x = 'Activation Ratio\n(pZEV)', y='Predicted Activation\nRatio (pZEV)', color = 'Density') +
  scale_color_viridis()
print(p)
dev.off()


