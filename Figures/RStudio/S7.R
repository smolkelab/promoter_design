library(data.table)
library(ggplot2)
library(MASS)
library(viridis)
setwd('D:/Promoter Design Data/FACS-Seq/')
dat.mer = fread('means_merged_GPD.csv'); dat.mer$Fate = 'Accepted'

setnames(dat.mer, c('V1', 'V2', 'V3'), c('Seq', 'Strength.A', 'Strength.B'))
dat.mer$Reps.Consistent = abs(dat.mer$Strength.A - dat.mer$Strength.B) < 0.2
dat.mer$Extreme = 0
dat.mer$Strength = apply(cbind(dat.mer$Strength.A, dat.mer$Strength.B), 1, mean)
dat.mer$Extreme[dat.mer$Strength <= -0.8] = -1
dat.mer$Extreme[dat.mer$Strength >= 0.7] = 1

dat.mer$Accepted = dat.mer$Reps.Consistent & dat.mer$Extreme == 0
dat.mer$Fate = 1
dat.mer$Fate[dat.mer$Extreme == -1] = 2
dat.mer$Fate[dat.mer$Extreme == 1] = 3
dat.mer$Fate[!dat.mer$Reps.Consistent] = 4

dat.mer.use = dat.mer[dat.mer$Fate == 1,]

# cf. https://slowkow.com/notes/ggplot2-color-by-density/

get_density = function(x, y, n = 200) {
  dens = MASS::kde2d(x = x, y = y, n = n)
  ix = findInterval(x, dens$x)
  iy = findInterval(y, dens$y)
  ii = cbind(ix, iy)
  return(dens$z[ii])
}

dat.mer.use$Density = get_density(dat.mer.use$Strength.A, dat.mer.use$Strength.B)

png(filename = 'D:/Promoter Design Data/Figures/PNGs/S7.png',
    units = 'cm', width = 10, height = 8, res = 600)

p = ggplot(dat.mer.use, aes(x = Strength.A, y = Strength.B, color = Density)) + 
  geom_point(alpha = 0.02, size = 0.2, stroke = 0, shape = 16) + theme_bw() + 
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold'),
        legend.text = element_text(size=6), legend.title = element_text(size=8, face='bold'),
        legend.key.size = unit(0.2,'cm')) +
  labs(x = 'pGPD Promoter Activity (Replicate A) (log10)', y='pGPD Promoter Activity (Replicate B) (log10)',
       color = 'Density') +
  geom_abline(slope = 1, intercept = c(-0.2, 0.2), lwd = 0.2, lty = 2) +
  scale_color_viridis()
print(p)
dev.off()

# for text
summary(lm(Strength.B ~ Strength.A, data = dat.mer.use)) # R2 = 0.94

means.avg = apply(cbind(dat.mer.use$Strength.A, dat.mer.use$Strength.B), 1, mean)
range(means.avg)
mean(means.avg)
sd(means.avg)



