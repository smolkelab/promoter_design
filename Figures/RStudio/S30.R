# given CSVs containing MLE scores and RMSDs vs. final means,
# generate graphs showing how sensitive mean predictions are to hyperparameter choice.
library(data.table)
setwd('D:/Promoter Design Data/Mean Fitting Scores/')

# manually pulled out of the fit logs
# On instance #3, ~/facs-seq_test/[GPD or ZEV]/nextseq/fit_means.log and ~/facs-seq_test/design_testing/fit_means.log
# GPD: global_ests: [{'fuzz': 0.06999999999999999, 'sigma': 0.038}, {'fuzz': 0.06999999999999999, 'sigma': 0.042}]
# ZEV: global_ests: [{'fuzz': 0.062000000000000006, 'sigma': 0.046000000000000006}, {'fuzz': 0.054, 'sigma': 0.02600000000
# 0000002}]
# FS9: global_ests: [{'fuzz': 0.046000000000000006, 'sigma': 0.082}, {'fuzz': 0.054, 'sigma': 0.058}]

contour.from.3col = function(dat, prefix, xpos, ypos, ...) {
  z = reshape(dat, timevar='sigma', idvar='fuzz',direction='wide')
  z = as.matrix(z)
  rownames(z) = z[,1]
  z = z[,2:ncol(z)]
  colvals = as.numeric(sapply(strsplit(colnames(z),prefix), '[[', 2))
  contour(z, x = as.numeric(rownames(z)), y = colvals, xlab = 'Fuzz', ylab = 'Sigma', ...)
  points(xpos,ypos,pch=16, col = 'red')
}

wrapper.2plot = function(fn, maintxt = '', xpos,ypos) {
  x = fread(fn)
  y = x; y$score = NULL
  y2 = x; y2$rmsd = NULL
  contour.from.3col(dat=y, prefix = 'rmsd.', xpos=xpos,ypos=ypos,
                    main = paste0(maintxt, 'RMSD', collapse = ''))
  contour.from.3col(dat=y2, prefix = 'score.', xpos=xpos,ypos=ypos,
                    main = paste0(maintxt, 'Score', collapse = ''))
}


png(filename = 'D:/Promoter Design Data/Figures/PNGs/S30.png',
  units = 'in', res = 144, width = 7.5, height = 5)
par(mfrow = c(3,4), mar = c(4.1, 4, 3.1, 1.1) )
#wrapper.2plot('score_A_GPD.txt', maintxt='pGPD-A: ',xpos = 0.07, ypos = 0.038)
#wrapper.2plot('score_B_GPD.txt', maintxt='pGPD-B: ',xpos = 0.07, ypos = 0.042)
#wrapper.2plot('score_A_ZEV.txt', maintxt='pZEV-Uninduced: ', xpos = 0.062, ypos = 0.046)
#wrapper.2plot('score_B_ZEV.txt', maintxt='pZEV-Induced: ' ,xpos = 0.054, ypos = 0.026)
#wrapper.2plot('score_A_validation.txt', maintxt='Validation-Uninduced: ', xpos = 0.046, ypos = 0.082)
#wrapper.2plot('score_B_validation.txt', maintxt='Validation-Induced: ' ,xpos = 0.054, ypos = 0.058)
wrapper.2plot('score_A_GPD.txt', maintxt='A: ',xpos = 0.07, ypos = 0.038)
wrapper.2plot('score_B_GPD.txt', maintxt='B: ',xpos = 0.07, ypos = 0.042)
wrapper.2plot('score_A_ZEV.txt', maintxt='Uninduced: ', xpos = 0.062, ypos = 0.046)
wrapper.2plot('score_B_ZEV.txt', maintxt='Induced: ' ,xpos = 0.054, ypos = 0.026)
wrapper.2plot('score_A_validation.txt', maintxt='Uninduced: ', xpos = 0.046, ypos = 0.082)
wrapper.2plot('score_B_validation.txt', maintxt='Induced: ' ,xpos = 0.054, ypos = 0.058)
dev.off()


