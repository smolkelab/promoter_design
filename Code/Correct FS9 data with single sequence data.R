# Single-sequence testing results:
# Show concordance w/ validation FACS-seq; show performance of designed promoters
# vs. benchmarks

setwd('D:/Promoter Design Data/')
library(data.table)
library(ggplot2)
source('Code/FCS file analysis.R')

get.one.median = function(fn) {
  x = extract.data(load.fromfilename(fn))
  rats = log10(x[,10]/x[,8])
  return(median(rats))
}

dat = fread('Validation/validation_output.csv'); dat$V1 = NULL

# copy VYB results for 'a' to 'b' for GPD designs
dat[dat$Promoter == 'GPD']$R1_TRUE = dat[dat$Promoter == 'GPD']$R1_FALSE
dat[dat$Promoter == 'GPD']$R2_TRUE = dat[dat$Promoter == 'GPD']$R2_FALSE
dat[dat$Promoter == 'GPD']$R3_TRUE = dat[dat$Promoter == 'GPD']$R3_FALSE

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

# Did we choose this sequence because we expect it to be unreliable?
dat$trusted_a = dat$Scale_A == 0
dat$trusted_b = dat$Scale_B == 0

dat.use.a = data.frame(FS = dat$Means_A, VYB = dat$VYB_Mean_A, ymin = dat$bar.min.a, ymax = dat$bar.max.a, 
                       Promoter = dat$Promoter, rep = 'A', Offscale = !(dat$trusted_a))
dat.use.a = dat.use.a[apply(dat.use.a, 1, function(x) all(!is.na(x))),]

dat.use.b = data.frame(FS = dat$Means_B, VYB = dat$VYB_Mean_B, ymin = dat$bar.min.b, ymax = dat$bar.max.b,
                       Promoter = dat$Promoter, rep = 'B', Offscale = !(dat$trusted_b))

dat.use.b = dat.use.b[apply(dat.use.b, 1, function(x) all(!is.na(x))),]


# get linear models for reliable points
mod.a = lm(dat.use.a$VYB[!dat.use.a$Offscale & dat.use.a$VYB > -1] ~ dat.use.a$FS[!dat.use.a$Offscale & dat.use.a$VYB > -1])

dat.use.a$FS_Corr = dat.use.a$FS*coef(mod.a)[2] + coef(mod.a)[1]

mod.b = lm(dat.use.b$VYB[!dat.use.b$Offscale] ~ dat.use.b$FS[!dat.use.b$Offscale])

dat.use.b$FS_Corr = dat.use.b$FS*coef(mod.b)[2] + coef(mod.b)[1]

mod.a.corr = lm(dat.use.a$VYB[!dat.use.a$Offscale & dat.use.a$VYB > -1] ~ dat.use.a$FS_Corr[!dat.use.a$Offscale & dat.use.a$VYB > -1])
mod.b.corr = lm(dat.use.b$VYB[!dat.use.b$Offscale] ~ dat.use.b$FS_Corr[!dat.use.b$Offscale])



dat.use = rbind(dat.use.a, dat.use.b)
ggplot(dat.use, aes(x = FS, y = VYB, color = rep, shape = Promoter)) + theme_bw() + geom_point() +
  geom_abline(slope = coef(mod.a)[2], intercept = coef(mod.a)[1], col = 'red') + 
  geom_abline(slope = coef(mod.b)[2], intercept = coef(mod.b)[1], col = 'blue') +
  geom_abline(slope = 1, intercept = 0, lty = 2)

ggplot(dat.use, aes(x = FS_Corr, y = VYB, color = rep, shape = Promoter)) + theme_bw() + geom_point() +
  geom_abline(slope = coef(mod.a.corr)[2], intercept = coef(mod.a.corr)[1], col = 'red') + 
  geom_abline(slope = coef(mod.b.corr)[2], intercept = coef(mod.b.corr)[1], col = 'blue') +
  geom_abline(slope = 1, intercept = 0, lty = 2)

# also correct the validation results means
dat = fread('Validation/validation_output.csv'); dat$V1 = NULL
dat$Means_A = dat$Means_A*coef(mod.a)[2] + coef(mod.a)[1]
dat$Means_B = dat$Means_B*coef(mod.b)[2] + coef(mod.b)[1]
dat$Means_Avg = apply(cbind(dat$Means_A, dat$Means_B), 1, mean)
dat$AR = dat$Means_B - dat$Means_A
write.csv(dat, 'Validation/validation_output_fs_means_corrected.csv')

# use these models to correct the original validation FACS-Seq data
means.orig = fread('D:/Promoter Design Data/FACS-Seq/final_validation_FACS-Seq_means_just_designs.csv')
means.orig$V1 = NULL; means.orig$AR = NULL; means.orig$Means_Avg = NULL
means.orig$Means_A = means.orig$Means_A*coef(mod.a)[2] + coef(mod.a)[1]
means.orig$Means_B = means.orig$Means_B*coef(mod.b)[2] + coef(mod.b)[1]
means.orig$Means_Avg = apply(cbind(means.orig$Means_A, means.orig$Means_B), 1, mean)
means.orig$AR = means.orig$Means_B - means.orig$Means_A

write.csv(means.orig, 'D:/Promoter Design Data/FACS-Seq/final_validation_FACS-Seq_means_just_designs_corrected.csv')
  



