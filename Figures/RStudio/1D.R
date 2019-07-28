# Visualize pZEV dataset
require(data.table)
require(ggplot2)
require(MASS)
require(viridis)

# from old 'Fig 3 - ALL except A.R'
get_density = function(x, y, n = 200) {
  dens = MASS::kde2d(x = x, y = y, n = n)
  ix = findInterval(x, dens$x)
  iy = findInterval(y, dens$y)
  ii = cbind(ix, iy)
  return(dens$z[ii])
}

setwd('D:/Promoter Design Data/FACS-Seq/')
dat.mer.use = fread('means_trainable_ZEV.csv')
print(range(dat.mer.use$Uninduced))
print(mean(dat.mer.use$Uninduced))
print(sd(dat.mer.use$Uninduced))
print(range(dat.mer.use$Induced))
print(mean(dat.mer.use$Induced))
print(sd(dat.mer.use$Induced))

print(median(dat.mer.use$Uninduced))
print(IQR(dat.mer.use$Uninduced))
print(median(dat.mer.use$Induced))
print(IQR(dat.mer.use$Induced))


setnames(dat.mer.use, c('Seq', 'Uninduced', 'Induced'), c('Seq', 'Strength.A', 'Strength.B'))

# cf. https://slowkow.com/notes/ggplot2-color-by-density/


dat.mer.use$Density = get_density(dat.mer.use$Strength.A, dat.mer.use$Strength.B)
dat.mer.use$`Cell Count` = dat.mer.use$Density

png(filename = 'D:/Promoter Design Data/Figures/PNGs/1D.png',
    units = 'cm', width = 8, height = 6.4, res = 600)

p = ggplot(dat.mer.use, aes(x = Strength.A, y = Strength.B, color = `Cell Count`)) + 
  geom_point(alpha = 0.1, size = 0.2, stroke = 0, shape = 16) + theme_bw() + 
  theme(axis.title.x=element_blank()) + 
  theme(axis.title.y=element_blank()) +
  theme(axis.text.x = element_text(size=8), axis.text.y = element_text(size=8)) +
  #xlim(c(-2.0, 0.8)) + ylim(c(-0.4, 1.8)) + 
  theme(legend.title = element_blank()) + 
  theme(legend.key.size = unit(0.5, 'cm')) + 
  scale_color_viridis()
print(p)
dev.off()

# R2 of ZEV conditions
summary(lm(Strength.B ~ Strength.A, data = dat.mer.use)) # R2 = 0.21


