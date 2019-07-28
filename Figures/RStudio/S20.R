library(data.table)
library(ggplot2)
setwd('D:/Promoter Design Data/Validation')

dat = fread('validation_output_fs_means_corrected.csv'); dat$V1 = NULL

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

dat = dat[is.na(dat$Seqs),]
dat$Source = 'Native'
dat$Source[grepl('ZEV', dat$Name)] = 'ZEV'

png(filename = 'D:/Promoter Design Data/Figures/PNGs/S20_raw.png',
    units = 'cm', width = 16, height = 12, res = 600)

p = ggplot(dat, aes(x = VYB_Mean_A, y = VYB_Mean_B, color = Source, 
                xmin = bar.min.a, xmax = bar.max.a,
                ymin = bar.min.b, ymax = bar.max.b)) + 
  geom_point(size = 1, stroke = 0) + geom_errorbarh(size = 0.3) + geom_errorbar(size = 0.3) + 
  theme_bw() + geom_abline(slope = 1, intercept = 0, lty = 2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face='bold')) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10, face='bold'),
        legend.text = element_text(size=10), legend.title = element_text(size=10, face='bold'),
        legend.key.size = unit(0.5,'cm')) +
  labs(x = 'Uninduced Promoter Activity\nComparison Sequences', 
       y = 'Induced Promoter Activity\nComparison Sequences')
print(p)
dev.off()


