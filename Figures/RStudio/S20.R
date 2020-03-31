library(data.table)
library(ggplot2)
setwd('D:/Promoter Design Data/Validation')

dat = fread('validation_output_fs_means_corrected.csv'); dat$V1 = NULL
dat$SEM.scale = 1/sqrt(3)
dat$VYB_Mean_A = apply(cbind(dat$R1_FALSE, dat$R2_FALSE, dat$R3_FALSE), 1, mean)
dat$VYB_SD_A = apply(cbind(dat$R1_FALSE, dat$R2_FALSE, dat$R3_FALSE), 1, sd)
dat$VYB_Mean_B = apply(cbind(dat$R1_TRUE, dat$R2_TRUE, dat$R3_TRUE), 1, mean)
dat$VYB_SD_B = apply(cbind(dat$R1_TRUE, dat$R2_TRUE, dat$R3_TRUE), 1, sd)
dat$VYB_AR = dat$VYB_Mean_B - dat$VYB_Mean_A
dat$VYB_AR_SD = sqrt(dat$VYB_SD_A^2 + dat$VYB_SD_B^2)
dat$bar.min.a = dat$VYB_Mean_A-(dat$VYB_SD_A*dat$SEM.scale)
dat$bar.max.a = dat$VYB_Mean_A+(dat$VYB_SD_A*dat$SEM.scale)
dat$bar.min.b = dat$VYB_Mean_B-(dat$VYB_SD_B*dat$SEM.scale)
dat$bar.max.b = dat$VYB_Mean_B+(dat$VYB_SD_B*dat$SEM.scale)
dat$bar.min.ar = dat$VYB_AR-(dat$VYB_AR_SD*dat$SEM.scale)
dat$bar.max.ar = dat$VYB_AR+(dat$VYB_AR_SD*dat$SEM.scale)

dat = dat[is.na(dat$Seqs),]
dat$Source = 'Native'
dat$Source[grepl('ZEV', dat$Name)] = 'ZEV'

melt.for.jitter.s20 = function(df, is.induced) {
  if (is.induced) {df = data.frame(Name=df$Name, Source = df$Source, R1=df$R1_TRUE, R2=df$R2_TRUE, R3=df$R3_TRUE)}
  else {df = data.frame(Name=df$Name, Source = df$Source, R1=df$R1_FALSE, R2=df$R2_FALSE, R3=df$R3_FALSE)}
  df = melt(df, id.vars = c('Name', 'Source'), measure.vars = c('R1', 'R2', 'R3'))
  df = data.table(df)
  return(df)
}

pts.u = melt.for.jitter.s20(dat, FALSE)
setnames(pts.u, 'value', 'Uninduced')
pts.i = melt.for.jitter.s20(dat, TRUE)
setnames(pts.i, 'value', 'Induced')
pts = cbind(pts.u, pts.i)
pts = pts[,c(1,2,4,8)]

# Source Data
dat.source = data.frame(Name=dat$Name, R1_FALSE = dat$R1_FALSE, R2_FALSE = dat$R2_FALSE, R3_FALSE = dat$R3_FALSE,
                        R1_TRUE = dat$R1_TRUE, R2_TRUE = dat$R2_TRUE, R3_TRUE = dat$R3_TRUE)
write.csv(dat.source, 'D:/Promoter Design Data/Source Data/S20.csv', quote = FALSE, row.names = FALSE)

png(filename = 'D:/Promoter Design Data/Figures/PNGs/S20_raw.png',
    units = 'cm', width = 16, height = 12, res = 600)

p = ggplot(dat, aes(x = VYB_Mean_A, y = VYB_Mean_B, color = Source, 
                xmin = bar.min.a, xmax = bar.max.a,
                ymin = bar.min.b, ymax = bar.max.b)) + 
  geom_point(size = 1, stroke = 0) +
  geom_errorbarh(size = 0.3, height = 0.05) + geom_errorbar(size = 0.3, width = 0.05) + 
  theme_bw() + geom_abline(slope = 1, intercept = 0, lty = 2) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face='bold')) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10, face='bold'),
        legend.text = element_text(size=10), legend.title = element_text(size=10, face='bold'),
        legend.key.size = unit(0.5,'cm')) +
  labs(x = 'Uninduced Promoter Activity\nComparison Sequences', 
       y = 'Induced Promoter Activity\nComparison Sequences')
q = p + geom_point(data=pts, aes(x=Uninduced, y=Induced, color = Source), 
                   alpha = 0.5, size = 0.5, inherit.aes = FALSE)
print(q)
dev.off()


