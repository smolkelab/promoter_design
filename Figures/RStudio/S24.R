source('D:/Promoter Design Data/Code/Bin edge calculation.R')

png('D:/Promoter Design Data/Figures/PNGs/S24.png',
    units = 'in', width = 6, height = 4, res = 600)
### Generate a final plot ###
par(mfrow = c(1,3), cex.lab = 0.75, cex.main = 0.75)

### Plot pZEV ###
plot(from_rats, to_rats, pch = '.', xlab = 'Promoter Activity (log10) - As Measured', 
     ylab = 'Promoter Activity (log10) - Reference Equivalent',
     xlim = c(-2, 1), ylim = c(-0.5, 2.0))
abline(lm(to_rats ~ from_rats), lty = 2)
abline(v = pzev.orig, lty = 2)
abline(h = pzev.cal, lty = 2)

### Plot uninduced ###

plot(output_0[[1]], rats_0, xlab = 'Promoter Activity (log10) - As Measured', 
     ylab = 'Promoter Activity (log10) - Reference Equivalent', pch = '.',
     xlim = c(-4, 1), ylim = c(-4, 3))
for(i in 2:length(output_0)) { points(output_0[[i]], rats_0, col = i, pch = '.') }
abline(lm(Uninduced_Mean ~ orig, data = final.frame), lty = 2)
for(i in 1:nrow(final.frame)) { abline(v = final.frame$orig[i], lty = 2) }
for(i in 1:nrow(final.frame)) { abline(h = final.frame$Uninduced_Mean[i], lty = 2) }

### Plot induced ###
plot(output_1[[1]], rats_1, xlab = 'Promoter Activity (log10) - As Measured', 
     ylab = 'Promoter Activity (log10) - Reference Equivalent', pch = '.',
     xlim = c(-4, 1), ylim = c(-4, 3))
for(i in 2:length(output_1)) { points(output_1[[i]], rats_1, col = i, pch = '.') }
abline(lm(Induced_Mean ~ orig, data = final.frame), lty = 2)
for(i in 1:nrow(final.frame)) { abline(v = final.frame$orig[i], lty = 2) }
for(i in 1:nrow(final.frame)) { abline(h = final.frame$Induced_Mean[i], lty = 2) }

dev.off()
# Final output: ZEV
print(pzev.cal)
# Final output: Validation Uninduced
print(final.frame$Uninduced_Mean)
# Final output: Validation Induced
print(final.frame$Induced_Mean)
