setwd('D:/Promoter Design Data/')
library(data.table)
library(reshape2)
source('Code/FCS file analysis.R')
dat = fread('Validation/validation_output_fs_means_corrected.csv'); dat$V1 = NULL
dat.u = melt(dat, id.vars = 'Name', measure.vars = c('R1_FALSE', 'R2_FALSE', 'R3_FALSE'))
res.aov = aov(value ~ Name, data = dat.u)

png(filename = 'D:/Promoter Design Data/Figures/PNGs/S31_A.png',
    units = 'cm', width = 12, height = 12, res = 600)
plot(res.aov, 2)
dev.off()
dat.i = melt(dat, id.vars = 'Name', measure.vars = c('R1_TRUE', 'R2_TRUE', 'R3_TRUE'))
dat.i = dat.i[!is.na(dat.i$value),]
res.aov = aov(value ~ Name, data = dat.i)
png(filename = 'D:/Promoter Design Data/Figures/PNGs/S31_B.png',
    units = 'cm', width = 12, height = 12, res = 600)
plot(res.aov, 2)
dev.off()
