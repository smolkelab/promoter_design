setwd('D:/Promoter Design Data/FACS-Seq')
require(data.table)
require(ggplot2)

zev.ar.ex = c(10, 17, 40:50)

dat.mer.use = fread('final_validation_FACS-Seq_means_just_designs_corrected.csv')
dat.mer.use = dat.mer.use[dat.mer.use$Experiment %in% zev.ar.ex,]

label.key = fread('supplementary_design_key.csv')
label.key$Final = sapply(label.key$Final, function(x) gsub('\\n', '\n', x, fixed = TRUE))

for(i in 1:nrow(label.key)) {
  dat.mer.use$Experiment[dat.mer.use$Experiment == label.key$Orig[i]] = label.key$Final[i]
}

dat.mer.use$Experiment = as.factor(dat.mer.use$Experiment)
dat.mer.use$Experiment = factor(dat.mer.use$Experiment, 
                                   levels = levels(factor(dat.mer.use$Experiment))[c(12, 7, 11, 5, 8, 1, 3, 6, 9, 2, 4, 10)] )

# for source data
dat.source = data.table(Experiment=dat.mer.use$Experiment, Value=dat.mer.use$AR)
dat.source$Experiment = sapply(dat.source$Experiment, function(x) gsub('\n', '-', x, fixed = TRUE))
dat.source = dat.source[order(dat.source$Experiment, dat.source$Value),]
write.csv(dat.source, 'D:/Promoter Design Data/Source Data/3C_S17.csv', quote = FALSE, row.names = FALSE)

stat_box_data = function(y) { return(data.frame(y = min(max(y) + 40, 200), label = length(y))) }

png(filename = 'D:/Promoter Design Data/Figures/Final PNGs/S17.png',
    units = 'cm', width = 12, height = 8, res = 600)

p = ggplot(dat.mer.use, aes(Experiment, 10^AR)) + #, color = Experiment)) + #, shape = Offscale)) + 
  
  geom_boxplot(data = dat.mer.use, outlier.size = 0, outlier.shape = NA) + # coef = 0, 
  geom_jitter(data = dat.mer.use, width = 0.1, height = 0, size = 0.5) +
  theme_bw() + 
  stat_summary(fun.data = stat_box_data, geom = "text", hjust = 0.5, vjust = 0.9, size=3) +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face='bold')) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold'),
        legend.title = element_blank(), axis.title.x=element_blank(), legend.position = 'none') +
  labs(y='Activation Ratio\nAll Tested pZEV-Activation Ratio Designs')
print(p)
dev.off()

