setwd('D:/Promoter Design Data/FACS-Seq/')
library(data.table)
library(ggplot2)

dat.mer.use = fread('final_validation_FACS-Seq_means_just_designs_corrected.csv')
setnames(dat.mer.use, c('Means_A', 'Means_B'), c('Strength.A', 'Strength.B'))

png(filename = 'D:/Promoter Design Data/Figures/Final PNGs/S12.png',
    units = 'cm', width = 16, height = 32/3, res = 600)

p = ggplot(dat.mer.use, aes(x = Strength.A, y = Strength.B)) +#, color = Density)) + 
  geom_point(alpha = 1, size = 0.5, stroke = 0, shape = 16) + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 8, face='bold')) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold')) +
  theme(legend.position='none') +
  labs(x = 'Validation FACS-seq\nUninduced Promoter Activity (log10)', y='Validation FACS-seq\nInduced Promoter Activity (log10)') +
  scale_x_continuous(breaks = seq(from = -1, to = 2, by = 0.5))
print(p)
dev.off()


