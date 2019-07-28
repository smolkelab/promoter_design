setwd('D:/Promoter Design Data/FACS-Seq')
require(data.table)
require(ggplot2)

zev.i.ex = c(9,16,29:39)

dat.mer.use = fread('final_validation_FACS-Seq_means_just_designs_corrected.csv')
dat.mer.use = dat.mer.use[dat.mer.use$Experiment %in% zev.i.ex,]

label.key = fread('supplementary_design_key.csv')
label.key$Final = sapply(label.key$Final, function(x) gsub('\\n', '\n', x, fixed = TRUE))

for(i in 1:nrow(label.key)) {
  dat.mer.use$Experiment[dat.mer.use$Experiment == label.key$Orig[i]] = label.key$Final[i]
}

dat.mer.use$Experiment = as.factor(dat.mer.use$Experiment)
#dat.mer.use$Offscale = dat.mer.use$Means_B == max(dat.mer.use$Means_B)
offscale.val = max(dat.mer.use$Means_B)
#offscale.val = 130
dat.mer.use$Means_B[dat.mer.use$Means_B > offscale.val] = log10(offscale.val)
dat.mer.use$Experiment = factor(dat.mer.use$Experiment, 
                                   levels = levels(factor(dat.mer.use$Experiment))[c(12, 7, 11, 5, 8, 1, 3, 6, 9, 2, 4, 10)] )

png(filename = 'D:/Promoter Design Data/Figures/Final PNGs/S16.png',
    units = 'cm', width = 12, height = 8, res = 600)

p = ggplot(dat.mer.use, aes(Experiment, 10^Means_B, color = Experiment)) + #, shape = Offscale)) + 
  
  geom_boxplot(data = dat.mer.use, outlier.size = 0, coef = 0, outlier.shape = NA) + 
  geom_jitter(data = dat.mer.use, width = 0.1, height = 0, size = 0.5) +
  geom_hline(yintercept = 10^offscale.val, lty = 2) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 8, face='bold')) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold'),
        legend.title = element_blank(), axis.title.x=element_blank(), legend.position = 'none') +
  labs(y='Induced Promoter Activity\nAll Tested pZEV-Induced Designs')
print(p)
dev.off()

# For text: Tests for penalty/no-penalty sets
fix.m = function(x) {10^x[!is.na(x)]}
#NC/S
wilcox.test(fix.m(dat.mer.use$Means_B[dat.mer.use$Experiment == 'NC\nNP\nS']), 
            fix.m(dat.mer.use$Means_B[dat.mer.use$Experiment == 'NC\nP\nS']), 'less') # 2.6e-14
#NC/E
wilcox.test(fix.m(dat.mer.use$Means_B[dat.mer.use$Experiment == 'NC\nNP\nE']), 
            fix.m(dat.mer.use$Means_B[dat.mer.use$Experiment == 'NC\nP\nE']), 'less') # 0.26

median(fix.m(dat.mer.use$Means_B[dat.mer.use$Experiment == 'NC\nP\nE'])) # with extrapolation penalty
median(fix.m(dat.mer.use$Means_B[dat.mer.use$Experiment == 'NC\nNP\nE'])) # without

#C/E
wilcox.test(fix.m(dat.mer.use$Means_B[dat.mer.use$Experiment == 'C\nNP\nE']), 
            fix.m(dat.mer.use$Means_B[dat.mer.use$Experiment == 'C\nP\nE']), 'less') # 0.0009

#NC/G
wilcox.test(fix.m(dat.mer.use$Means_B[dat.mer.use$Experiment == 'NC\nNP\nG']), 
            fix.m(dat.mer.use$Means_B[dat.mer.use$Experiment == 'NC\nP\nG']), 'less') # 1.08e-6

#C/G
wilcox.test(fix.m(dat.mer.use$Means_B[dat.mer.use$Experiment == 'C\nNP\nG']), 
            fix.m(dat.mer.use$Means_B[dat.mer.use$Experiment == 'C\nP\nG']), 'less') # 0.0012


