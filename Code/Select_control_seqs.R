# Find desired sequences for synthesis as controls in FS9.
# Return them padded to 363 bp for compatibility with designed sequences.
setwd('D:/Promoter Design Data/Code/')
library(data.table)
source('Functions_select_control_seqs.R')

rs = 2017
nc = 40 # number of controls per set - careful if reducing, some sets are fixed size (all seqs beyond a threshold)
range.reps = 5 # when sampling a range of strengths, how many at each strength?
zev.grid.a = 10 # number of points to sample in the ZEV grid, uninduced dimension (even sampling of ZEV data)
zev.grid.b = 10 # number of points to sample in the ZEV grid, induced dimension (even sampling of ZEV data)
zev.grid.reps = 3 # number of sequences to sample from each point in the ZEV grid
outlier.thresh.gpd = 0.3
outlier.thresh.zev = 0.2
if( nc %% range.reps != 0) { warning('range.reps does not divide nc evenly')}
# pads, for rejected means: means.all has the pads already
pads.gpd = c('TACGTAAATAATTAATAGTAGTGAC','TGTCTAAAGGTGAAGAATTATTCAC')
pads.zev = c('TTTATCATTATCAATACTCGCCATTTCAAAGAATACGTAAATAATTAATAGTAGTGAC',
             'TGTCTAAAGGTGAAGAATTATTCACTGGTGTTGTCCCAATTTTGGTTGAATTAGATGG')

setwd('D:/Promoter Design Data/FACS-Seq/')
means.all = fread('means_trainable_joined.csv')
means.gpd.r = fread('means_rejected_GPD.csv')
means.zev.r = fread('means_rejected_ZEV.csv')
preds.test = fread('D:/Promoter Design Data/Preds/preds_joined_all.csv')
preds.test$GPD = preds.test$Strength_A == preds.test$Strength_B
preds.gpd = preds.test[preds.test$GPD,]
preds.zev = preds.test[!preds.test$GPD,]
preds.zev$AR = preds.zev$Strength_B - preds.zev$Strength_A
preds.zev$Pred_AR = preds.zev$Pred_Strength_B - preds.zev$Pred_Strength_A
setnames(means.gpd.r, c('V1','V2','V3'), names(means.all))
setnames(means.zev.r, c('V1','V2','V3'), names(means.all))
means.gpd.r$Seq = paste0(pads.gpd[1], means.gpd.r$Seq, pads.gpd[2])
means.zev.r$Seq = paste0(pads.zev[1], means.zev.r$Seq, pads.zev[2])
means.all$GPD = means.all$Strength_A == means.all$Strength_B
means.gpd = means.all[means.all$GPD,]
means.zev = means.all[!means.all$GPD,]
means.zev$AR = means.zev$Strength_B - means.zev$Strength_A

setwd('D:/Promoter Design Data/Validation/')
# GPD
#  High outlier
get.outlier(means.gpd.r, 'controls/00_GPD_top_outlier.csv', floor(nc/2), rs, TRUE)
# Top - no outlier
means.use = copy(means.gpd)
means.use$Strength = apply(cbind(means.use$Strength_A, means.use$Strength_B),1,mean)
sample.top(means.use, 'controls/01_GPD_top_no-outlier.csv', nc, rs)
# Range
means.use = copy(means.gpd)
means.use$Strength = apply(cbind(means.use$Strength_A, means.use$Strength_B),1,mean)
sample.range(means.use, 'controls/02_GPD_range.csv', nc*2, rs, range.reps)
#Low outlier
get.outlier(means.gpd.r, 'controls/03_GPD_bottom_outlier.csv',floor(nc/2),rs,FALSE)

#Model extreme error (test data only)
preds.use = copy(preds.gpd)
setnames(preds.use, c('Seqs','Strength_A'), c('Seq','Strength'))
preds.use$Pred_Strength = apply(cbind(preds.use$Pred_Strength_A, preds.use$Pred_Strength_B),1,mean)
get.outliers.worst.and.range(preds.use, 'controls/04_GPD_test-sample-outliers.csv',
                             nc, rs, outlier.thresh.gpd)

#Reliable strong predictions: in test data, get top 'nc' by sum of actual and predicted strength
preds.use$Strength.plus.Pred = preds.use$Strength + preds.use$Pred_Strength
set.seed(rs); preds.use = preds.use[sample(nrow(preds.use)),]
preds.use = preds.use[order(preds.use$Strength.plus.Pred)]
tail(preds.use,nc)
write.csv(data.frame(Seqs = tail(preds.use,nc)$Seq, Scores = tail(preds.use,nc)$Strength.plus.Pred), 
          file = 'controls/05_GPD_test-good-strong-preds.csv', row.names = FALSE, quote = FALSE)

# ZEV
# Outlier distribution is different for ZEV: e.g. only 1 'max' on the induced-positive end
# Same questions though: What's going on with outliers? What are the best sequences in the original data like?
# When the model and data disagree, who's right?
# Don't worry about "useful" too much - these are a lot stronger than the GPD anyway
# (important for design most likely, though)
# -outliers-
# Populations to sample:
# means.zev.r: Strength_A == min, Strength_B == min, Strength_A ~ max (incl. the one weird one with high Strength_B also)
# Strength_A == min
# hack for compatibility with 'get.outlier'
means.zev.r.use = copy(means.zev.r)
means.zev.r.use$Strength_B = means.zev.r.use$Strength_A
get.outlier(means.zev.r.use, 'controls/06_ZEV-uninduced_bottom_outlier.csv', floor(nc/2), rs, FALSE)

# Strength_B == min
means.zev.r.use = copy(means.zev.r)
means.zev.r.use$Strength_A = means.zev.r.use$Strength_B
get.outlier(means.zev.r.use, 'controls/07_ZEV-induced_bottom_outlier.csv', floor(nc/2), rs, FALSE)

# Strength_A ~ max
means.zev.r.use = means.zev.r[means.zev.r$Strength_A > 0.5]
set.seed(rs); means.zev.r.use = means.zev.r.use[sample(nrow(means.zev.r.use)),]
if(nrow(means.zev.r.use) > floor(nc/2)) { means.zev.r.use = means.zev.r.use[1:floor(nc/2),] }
write.csv(data.frame(Seqs = means.zev.r.use$Seq, Scores = means.zev.r.use$Strength_A), 
          file = 'controls/08_ZEV-uninduced_top_outlier.csv', row.names = FALSE, quote = FALSE)

# -best-
# means.zev highest Strength_B (all)
means.zev.use = copy(means.zev)
set.seed(rs); means.zev.use = means.zev.use[sample(nrow(means.zev.use)),]
means.zev.use = means.zev.use[order(means.zev.use$Strength_B),]
means.zev.use.top = tail(means.zev.use, floor(nc))
means.zev.use = means.zev.use[1:(nrow(means.zev.use) - floor(nc)),]
means.zev.use = means.zev.use[means.zev.use$Strength_B > 1.5,]
set.seed(rs); means.zev.use = means.zev.use[sample(nrow(means.zev.use), floor(nc)),]
means.zev.use = rbind(means.zev.use.top, means.zev.use)
write.csv(data.frame(Seqs = means.zev.use$Seq, Scores = means.zev.use$Strength_B), 
          file = 'controls/09_ZEV-Induced_top_no-outlier.csv', row.names = FALSE, quote = FALSE)

#  highest AR (NB: check for overlap): 40 sequences greater than 1.9; 60 sampled from the ones between 1.8 and 1.9
means.zev.use = copy(means.zev)
set.seed(rs); means.zev.use = means.zev.use[sample(nrow(means.zev.use)),]
means.zev.use = means.zev.use[order(means.zev.use$AR),]
means.zev.use.top = tail(means.zev.use, floor(nc/2))
means.zev.use = means.zev.use[1:(nrow(means.zev.use) - floor(nc/2)),]
means.zev.use = means.zev.use[means.zev.use$AR > 1.8,]
set.seed(rs); means.zev.use = means.zev.use[sample(nrow(means.zev.use), floor(nc/2)),]
means.zev.use = rbind(means.zev.use.top, means.zev.use)
write.csv(data.frame(Seqs = means.zev.use$Seq, Scores = means.zev.use$AR), 
          file = 'controls/10_ZEV-AR_high-selected-and-highest_no-outlier.csv', row.names = FALSE, quote = FALSE)

# -general-
# range of data (grid sample in main strength A/strength B space)
a_range = seq(from = -0.55, to = 0.6, length.out = zev.grid.a)
b_range = seq(from = 0.7, to = 1.55, length.out = zev.grid.b)
pts = expand.grid(a_range, b_range); names(pts) = c('A','B')
pts = pts[pts$B - pts$A > 0.8,]

set.seed(rs); means.zev.use = means.zev[sample(nrow(means.zev)),]
for(i in 1:nrow(pts)) {
  pt = c(pts$A[i], pts$B[i])
  means.zev.use$dist = (means.zev.use$Strength_A - pt[1])^2 + (means.zev.use$Strength_B - pt[2])^2
  means.zev.use = means.zev.use[order(means.zev.use$dist),]
  tmp = means.zev.use[1:zev.grid.reps,]
  if(i == 1) { ans = tmp }
  else { ans = rbind(ans, tmp) }
}
write.csv(data.frame(Seqs = ans$Seq, Scores_A = ans$Strength_A, Scores_B = ans$Strength_B), 
          file = 'controls/11_ZEV_grid.csv', row.names = FALSE, quote = FALSE)

# -modeling-
# Strength A outliers: worst 'nc' and 'nc' outside ZEV threshold 
preds.zev.use = copy(preds.zev)
setnames(preds.zev.use, c('Seqs','Strength_A','Pred_Strength_A'), c('Seqs','Strength','Pred_Strength'))
get.outliers.worst.and.range(preds.zev.use, 'controls/12_ZEV-A_test-sample-outliers.csv',
                             nc, rs, outlier.thresh.zev)

# Strength B outliers; worst 'nc' and 'nc' outside ZEV threshold
preds.zev.use = copy(preds.zev)
setnames(preds.zev.use, c('Seqs','Strength_B','Pred_Strength_B'), c('Seqs','Strength','Pred_Strength'))
get.outliers.worst.and.range(preds.zev.use, 'controls/13_ZEV-B_test-sample-outliers.csv',
                             nc, rs, outlier.thresh.zev)

#AR outliers; worst 'nc' and 'nc' outside ZEV threshold
preds.zev.use = copy(preds.zev)
setnames(preds.zev.use, c('Seqs','AR','Pred_AR'), c('Seq','Strength','Pred_Strength'))
get.outliers.worst.and.range(preds.zev.use, 'controls/14_ZEV-AR_test-sample-outliers.csv',
                             nc, rs, outlier.thresh.zev)

# Good predictions of strength for A, B, and AR
#A
preds.zev.use = copy(preds.zev)
set.seed(rs); preds.zev.use = preds.zev.use[sample(nrow(preds.zev.use)),]
setnames(preds.zev.use, c('Seqs','Strength_A','Pred_Strength_A'), c('Seq','Strength','Pred_Strength'))
preds.zev.use$Strength.plus.Pred = preds.zev.use$Strength + preds.zev.use$Pred_Strength
preds.zev.use = preds.zev.use[order(preds.zev.use$Strength.plus.Pred),]
preds.zev.use = tail(preds.zev.use, nc)
write.csv(data.frame(Seqs = preds.zev.use$Seq, Scores = preds.zev.use$Strength.plus.Pred), 
          file = 'controls/15_ZEV-A_test-good-strong-preds.csv', row.names = FALSE, quote = FALSE)

#B
preds.zev.use = copy(preds.zev)
set.seed(rs); preds.zev.use = preds.zev.use[sample(nrow(preds.zev.use)),]
setnames(preds.zev.use, c('Seqs','Strength_B','Pred_Strength_B'), c('Seq','Strength','Pred_Strength'))
preds.zev.use$Strength.plus.Pred = preds.zev.use$Strength + preds.zev.use$Pred_Strength
preds.zev.use = preds.zev.use[order(preds.zev.use$Strength.plus.Pred),]
preds.zev.use = tail(preds.zev.use, nc)
write.csv(data.frame(Seqs = preds.zev.use$Seq, Scores = preds.zev.use$Strength.plus.Pred), 
          file = 'controls/16_ZEV-B_test-good-strong-preds.csv', row.names = FALSE, quote = FALSE)

#AR
preds.zev.use = copy(preds.zev)
set.seed(rs); preds.zev.use = preds.zev.use[sample(nrow(preds.zev.use)),]
setnames(preds.zev.use, c('Seqs','AR','Pred_AR'), c('Seq','Strength','Pred_Strength'))
preds.zev.use$Strength.plus.Pred = preds.zev.use$Strength + preds.zev.use$Pred_Strength
preds.zev.use = preds.zev.use[order(preds.zev.use$Strength.plus.Pred),]
preds.zev.use = tail(preds.zev.use, nc)
write.csv(data.frame(Seqs = preds.zev.use$Seq, Scores = preds.zev.use$Strength.plus.Pred), 
          file = 'controls/17_ZEV-AR_test-good-strong-preds.csv', row.names = FALSE, quote = FALSE)

