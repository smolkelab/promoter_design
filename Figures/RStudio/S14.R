# Adjust ZEV means with validation data
setwd('D:/Promoter Design Data/')
library(data.table)
library(ggplot2)

dat.ref = fread('D:/Promoter Design Data/FACS-Seq/final_validation_FACS-Seq_means_just_designs_corrected.csv'); dat.ref$V1 = NULL
reliable = c(9,10,11,15,16,17)
dat.ref = dat.ref[dat.ref$Experiment %in% reliable]
dat.zev = fread('D:/Promoter Design Data/FACS-Seq/means_trainable_ZEV.csv')

setnames(dat.ref, 'Seqs', 'Seq')
setkey(dat.ref, 'Seq')
setkey(dat.zev, 'Seq')

dat.mer = merge(dat.ref, dat.zev)

coefs.a = summary(lm(Means_A ~ Uninduced, data = dat.mer[dat.mer$Uninduced > 0,]))$coef[,1] # R2 = 0.85
dat.mer$Uninduced_Corr = dat.mer$Uninduced*coefs.a[2] + coefs.a[1]
coefs.b = summary(lm(Means_B ~ Induced, data = dat.mer[dat.mer$Induced > 0.9 & dat.mer$Induced < 1.8,]))$coef[,1] # R2 = 0.79
dat.mer$Induced_Corr = dat.mer$Induced*coefs.b[2] + coefs.b[1]

#print(ggplot(data = dat.mer, aes(x = Uninduced, y = Means_A)) + theme_bw() + geom_point() + 
#  geom_abline(slope = 1, intercept = 0,lty = 2))
#print(ggplot(data = dat.mer, aes(x = Uninduced_Corr, y = Means_A)) + theme_bw() + geom_point() + 
#  geom_abline(slope = 1, intercept = 0,lty = 2))

#print(ggplot(data = dat.mer, aes(x = Induced, y = Means_B)) + theme_bw() + geom_point() + 
#  geom_abline(slope = 1, intercept = 0,lty = 2))
#print(ggplot(data = dat.mer, aes(x = Induced_Corr, y = Means_B)) + theme_bw() + geom_point() + 
#  geom_abline(slope = 1, intercept = 0,lty = 2))



grid.y = c(dat.mer$Uninduced, dat.mer$Uninduced_Corr, dat.mer$Induced, dat.mer$Induced_Corr)
grid.x = c(dat.mer$Means_A, dat.mer$Means_A, dat.mer$Means_B, dat.mer$Means_B)
grid.c = c(rep('Uninduced', nrow(dat.mer)), rep('Uninduced', nrow(dat.mer)), rep('Induced', nrow(dat.mer)), rep('Induced', nrow(dat.mer)))
grid.r = c(rep('Original', nrow(dat.mer)), rep('Rescaled', nrow(dat.mer)), rep('Original', nrow(dat.mer)), rep('Rescaled', nrow(dat.mer)))

dat.plot = data.table(x = grid.x, y = grid.y, ro = grid.r, co = grid.c)

png(filename = 'D:/Promoter Design Data/Figures/Final PNGs/S14.png',
    units = 'cm', width = 12, height = 8, res = 600)

p = ggplot(data = dat.plot, aes(x = x, y = y, color = co)) + theme_bw() + geom_point() + 
  geom_abline(slope = 1, intercept = 0, lty = 2) + 
  facet_wrap(ro ~ .) + 
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold')) +
  labs(x = 'Promoter Activity - Validation FACS-seq (log10)', y = 'Promoter Activity - ZEV FACS-seq (log10)', color = '')
print(p)
dev.off()

dat.zev.orig = fread('D:/Promoter Design Data/FACS-Seq/means_trainable_ZEV.csv')
a_corr = dat.zev.orig$Uninduced*coefs.a[2] + coefs.a[1]
b_corr = dat.zev.orig$Induced*coefs.b[2] + coefs.b[1]
dat.zev.out = data.table(Seq = dat.zev.orig$Seq, Uninduced_Corr = a_corr, Induced_Corr = b_corr)
write.csv(dat.zev.out, 'D:/Promoter Design Data/FACS-Seq/means_trainable_ZEV_rescaled_to_validation.csv')



