setwd('D:/Promoter Design Data/Model Ensemble Logs/')

process.raw.log = function(fn) {
  lines = suppressWarnings(readLines(fn))
  lines = lines[(which(lines == 'None')+1):length(lines)]
  lines = lines[seq(from = 2, to = length(lines), by = 2)]
  lines = lines[lines != 'Predictions written']
  lines = sapply(lines, function(x) strsplit(x, 'loss:')[[1]])
  lines = t(lines[2:3,])
  lines[,1] = sapply(lines[,1], function(x) strsplit(x, '- val_')[[1]][1])
  ans = matrix(nrow = nrow(lines), ncol = ncol(lines))
  for(i in 1:ncol(ans)) { ans[,i] = as.numeric(lines[,i]) }
  colnames(ans) = c('Training', 'Validation'); return(ans)
}

png(filename = 'D:/Promoter Design Data/Figures/Final PNGs/S8.png',
    units = 'in', res = 144, width = 6, height = 6)
par(mfrow = c(3,3))
for(i in 1:length(dir())) {
  xlab.txt = paste0('Model ', i-1, ' Epoch', collapse = '')
  ylab.txt = paste0('Model ', i-1, ' Loss', collapse = '')
  
  lines = process.raw.log(dir()[i])
  plot(lines[,1], type = 'l', lwd = 2,
       xlim = c(0, 50), ylim = c(0, 0.02), xlab = xlab.txt, ylab = ylab.txt)
  legend('topright', fill = c('black', 'red'), legend = c('Training', 'Validation'))
  points(lines[,2], type = 'l', col = 'red', lwd = 2)
  abline(v = nrow(lines) - 5, lty = 2)
}
dev.off()

