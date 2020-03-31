library(data.table)
library(ggplot2)
library(ggseqlogo)

png.head = 'D:/Promoter Design Data/Figures/PNGs/'
df = fread('D:/Promoter Design Data/Designs/3_GPD_strong_nofilter_1sd_evolve-thresh_0.6_selected.txt')
df$Seqs = substr(df$Seqs, 26, 340)
# 312 bp = 6.4084 cm
# 15 bp = 0.31 cm
# 20 bp = 0.41 cm
# 13 bp = 0.27 cm

# whole sequence
#png(filename = paste0(png.head, 'S22A.png', collapse = ''),
#    units = 'cm', width = 21, height = 3.5, res = 600)
#ggplot() + geom_logo( df$Seqs ) + theme_logo() + ylim(0,2) + theme(axis.title.x=element_blank(),
#                                                                   axis.text.x=element_blank(),
#                                                                   axis.ticks.x=element_blank())
#dev.off()


# first GCR1: 15 bp
png(filename = paste0(png.head, '5B_pGPD_LOGO_GCR1-1.png', collapse = ''),
    units = 'cm', width = 8, height = 3.5, res = 600)
ggplot() + geom_logo( substr(df$Seqs, 56, 70) ) + theme_logo() + ylim(0,2) + theme(axis.title.x=element_blank(),
                                                                                   axis.text.x=element_blank(),
                                                                                   axis.ticks.x=element_blank())
dev.off()

# second GCR1: 15 bp
png(filename = paste0(png.head, '5B_pGPD_LOGO_GCR1-2.png', collapse = ''),
    units = 'cm', width = 8, height = 3.5, res = 600)
ggplot() + geom_logo( substr(df$Seqs, 89, 103) ) + theme_logo() + ylim(0,2) + theme(axis.title.x=element_blank(),
                                                                                    axis.text.x=element_blank(),
                                                                                    axis.ticks.x=element_blank())
dev.off()

# core upstream: 20 bp
png(filename = paste0(png.head, '5B_pGPD_LOGO_core.png', collapse = ''),
    units = 'cm', width = 8, height = 3.5, res = 600)
ggplot() + geom_logo( substr(df$Seqs, 200, 219) ) + theme_logo() + ylim(0,2) + theme(axis.title.x=element_blank(),
                                                                                     axis.text.x=element_blank(),
                                                                                     axis.ticks.x=element_blank())
dev.off()

# 5' UTR downstream: 13 bp
png(filename = paste0(png.head, '5B_pGPD_LOGO_5UTR.png', collapse = ''),
    units = 'cm', width = 8, height = 3.5, res = 600)
ggplot() + geom_logo( substr(df$Seqs, 303, 315) ) + theme_logo() + ylim(0,2) + theme(axis.title.x=element_blank(),
                                                                                     axis.text.x=element_blank(),
                                                                                     axis.ticks.x=element_blank())
dev.off()


#### ZEV-I

# 3' context
png.head = 'D:/Promoter Design Data/Figures/PNGs/'
df = fread('D:/Promoter Design Data/Designs/14_ZEV_induced_nofilter_1sd_evolve-thresh_1.6_selected.txt')
df$Seqs = substr(df$Seqs, 59, 307)
png(filename = paste0(png.head, 'S22A.png', collapse = ''),
    units = 'cm', width = 10, height = 3.5, res = 600)
ggplot() + geom_logo( substr(df$Seqs, 237, 249) ) + theme_logo() + ylim(0,2) + theme(axis.title.x=element_blank(),
                                                                   axis.text.x=element_blank(),
                                                                   axis.ticks.x=element_blank())
dev.off()

### ZEV-I upstream TATA
df = fread('D:/Promoter Design Data/Designs/14_ZEV_induced_nofilter_1sd_evolve-thresh_1.6_selected.txt')
df$Seqs = substr(df$Seqs, 59, 98)
png(filename = paste0(png.head, 'S22C.png', collapse = ''),
    units = 'cm', width = 21, height = 3.5, res = 600)
ggplot() + geom_logo( df$Seqs ) + theme_logo() + ylim(0,2) + theme(axis.title.x=element_blank(),
                                                                   axis.text.x=element_blank(),
                                                                   axis.ticks.x=element_blank())
dev.off()

### ZEV-AR

# whole sequence
df = fread('D:/Promoter Design Data/Designs/25_ZEV_AR_nofilter_1sd_evolve-thresh_1.85_selected.txt')
df$Seqs = substr(df$Seqs, 59, 307)
png(filename = paste0(png.head, 'S22B.png', collapse = ''),
    units = 'cm', width = 10, height = 3.5, res = 600)
ggplot() + geom_logo( substr(df$Seqs, 237, 249) ) + theme_logo() + ylim(0,2) + theme(axis.title.x=element_blank(),
                                                                   axis.text.x=element_blank(),
                                                                   axis.ticks.x=element_blank())
dev.off()

# GCTA extensions
df = fread('D:/Promoter Design Data/Designs/25_ZEV_AR_nofilter_1sd_evolve-thresh_1.85_selected.txt')
df$Seqs = substr(df$Seqs, 100, 150)
png(filename = paste0(png.head, 'S22D.png', collapse = ''),
    units = 'cm', width = 21, height = 3.5, res = 600)
ggplot() + geom_logo( df$Seqs ) + theme_logo() + ylim(0,2) + theme(axis.title.x=element_blank(),
                                                                   axis.text.x=element_blank(),
                                                                   axis.ticks.x=element_blank())
dev.off()

