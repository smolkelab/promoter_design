setwd('D:/Promoter Design Data/FACS-Seq')
require(data.table)
require(ggplot2)

gpd.strong.ex = c(1,5,18:28)

dat.mer.use = fread('final_validation_FACS-Seq_means_just_designs_corrected.csv')
dat.mer.use = dat.mer.use[dat.mer.use$Experiment %in% gpd.strong.ex,]

label.key = fread('supplementary_design_key.csv')
label.key$Final = sapply(label.key$Final, function(x) gsub('\\n', '\n', x, fixed = TRUE))

for(i in 1:nrow(label.key)) {
  dat.mer.use$Experiment[dat.mer.use$Experiment == label.key$Orig[i]] = label.key$Final[i]
}
dat.mer.use$Experiment = as.factor(dat.mer.use$Experiment)
offscale.val = max(dat.mer.use$Means_Avg)
#offscale.val = 25
#dat.mer.use$Means_Avg[dat.mer.use$Means_Avg > log10(offscale.val)] = log10(offscale.val)
dat.mer.use$Experiment = factor(dat.mer.use$Experiment, 
                                   levels = levels(factor(dat.mer.use$Experiment))[c(10, 6, 9, 5, 7, 1, 3, 8, 2, 4)] )

png(filename = 'D:/Promoter Design Data/Figures/Final PNGs/S15.png',
    units = 'cm', width = 12, height = 8, res = 600)

p = ggplot(dat.mer.use, aes(Experiment, 10^Means_Avg, color = Experiment)) + #, shape = Offscale)) + 
  
  geom_boxplot(data = dat.mer.use, outlier.size = 0, coef = 0, outlier.shape = NA) + 
  geom_jitter(data = dat.mer.use, width = 0.1, height = 0, size = 0.5) +
  geom_hline(yintercept = 10^offscale.val, lty = 2) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 8, face='bold')) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold'),
        legend.title = element_blank(), axis.title.x=element_blank(), legend.position = 'none') +
  labs(y='Promoter Activity\nAll Tested pGPD Designs')
print(p)
dev.off()

# For text: Tests for penalty/no-penalty sets
fix.m = function(x) {10^x[!is.na(x)]}
#NC/S
wilcox.test(fix.m(dat.mer.use$Means_Avg[dat.mer.use$Experiment == 'NC\nNP\nS']), 
            fix.m(dat.mer.use$Means_Avg[dat.mer.use$Experiment == 'NC\nP\nS']), 'less') # 5.9e-6
#NC/E
wilcox.test(fix.m(dat.mer.use$Means_Avg[dat.mer.use$Experiment == 'NC\nNP\nE']), 
            fix.m(dat.mer.use$Means_Avg[dat.mer.use$Experiment == 'NC\nP\nE']), 'less') # 0.003
#C/G
wilcox.test(fix.m(dat.mer.use$Means_Avg[dat.mer.use$Experiment == 'C\nNP\nG']), 
            fix.m(dat.mer.use$Means_Avg[dat.mer.use$Experiment == 'C\nP\nG']), 'less') # 7e-14









