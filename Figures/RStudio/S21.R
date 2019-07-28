library(data.table)
library(ggplot2)
setwd('D:/Promoter Design Data')

dat = fread('Validation/validation_output_fs_means_corrected.csv'); dat$V1 = NULL

w = 16
h = 7.5

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

png(filename = 'D:/Promoter Design Data/Figures/PNGs/S21.png',
    units = 'cm', width = w, height = h, res = 600)

dat.gpd = dat[dat$Name %in% c('GPD', 'TEF', 'CYC1') | 
                dat$Experiment.designed %in% c(22,23,24,28),]
dat.gpd = dat.gpd[(dat.gpd$trusted_a & dat.gpd$trusted_b) | is.na(dat.gpd$Experiment.designed),]

dat.gpd$`Design Type` = dat.gpd$Experiment.designed
dat.gpd$`Design Type`[is.na(dat.gpd$`Design Type`)] = 'Control'
dat.gpd$`Design Type`[dat.gpd$`Design Type` == 22] = 'Evolution-GC, No Exp. Penalty'
dat.gpd$`Design Type`[dat.gpd$`Design Type` == 23] = 'Evolution-GC'
dat.gpd$`Design Type`[dat.gpd$`Design Type` == 24] = 'Gradient, No Exp. Penalty'
dat.gpd$`Design Type`[dat.gpd$`Design Type` == 28] = 'Gradient*'
dat.gpd$`Design Type` = factor(dat.gpd$`Design Type`)

# move the non-tested designs to the end, hackily
dat.gpd$Experiment.designed[dat.gpd$Experiment.designed == 22] = 29
dat.gpd$Experiment.designed[dat.gpd$Experiment.designed == 24] = 30

dat.gpd$Name = factor(dat.gpd$Name, levels = dat.gpd$Name[order(dat.gpd$Experiment.designed, dat.gpd$VYB_Mean_A)])

# for building final sequence-strength table: pGPD
dat.gpd$Seqs.final = sapply(dat.gpd$Seqs, function(x) substr(x,26,nchar(x)-26))
#write.csv(dat.gpd, 'C:/Users/Ben/Dropbox/Lab/Lab Deliverables/Promoter Manuscript/single_seqs_pGPD.csv')
xstr = c(paste0('pGPD-', 1:25), 'pCYC1', 'pTEF1', 'pGPD')
p = ggplot(dat.gpd, aes(Name, VYB_Mean_A, color = `Design Type`, ymin = bar.min.a, ymax = bar.max.a)) + 
  geom_point(size = 1, stroke = 0) + geom_errorbar(size = 0.3) + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.4, size = 10, face='bold')) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold'),
        legend.text = element_text(size=6), legend.title = element_text(size=8, face='bold'),
        legend.key.size = unit(0.5,'cm')) +
  theme(axis.text.x=element_text(angle=315),
        axis.ticks.x=element_blank()) +
  scale_x_discrete(labels=xstr) +
  labs(x = '', y='Promoter Activity\nAll Tested pGPD Designs', color = 'Design Type')
print(p)
dev.off()


