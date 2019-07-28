setwd('D:/Promoter Design Data/FACS-Seq/')
library(data.table)
library(ggplot2)

#tables = fread('FS9_results_validate_means.csv')
tables = fread('final_validation_FACS-Seq_means_just_designs_corrected.csv')
tables$Experiment.designed = tables$Experiment
tables = tables[!is.na(tables$Means_A) & !is.na(tables$Means_B) & tables$Experiment.designed < 18,]

exp.names = c('GPD top outlier','GPD top inlier','GPD range','GPD bottom outlier','GPD test-data outlier','GPD test-strong, good preds','ZEV uninduced bottom outlier','ZEV induced bottom outlier','ZEV uninduced top outlier','ZEV induced top inlier','ZEV AR high inlier','ZEV grid','ZEV uninduced test outlier','ZEV induced test outlier','ZEV AR test outlier','ZEV uninduced test-strong, good preds','ZEV induced test-strong, good preds','ZEV AR test-strong, good preds')
# NB: some sequences appear in tables more than once: this is OK

pads.gpd = c('TACGTAAATAATTAATAGTAGTGAC','TGTCTAAAGGTGAAGAATTATTCAC')
pads.zev = c('TTTATCATTATCAATACTCGCCATTTCAAAGAATACGTAAATAATTAATAGTAGTGAC',
             'TGTCTAAAGGTGAAGAATTATTCACTGGTGTTGTCCCAATTTTGGTTGAATTAGATGG')
means.all = fread('means_trainable_joined.csv')
means.gpd.r = fread('means_rejected_GPD.csv')
means.zev.r = fread('means_rejected_ZEV.csv')
setnames(means.gpd.r, c('V1','V2','V3'), names(means.all))
setnames(means.zev.r, c('V1','V2','V3'), names(means.all))
means.gpd.r$Seq = paste0(pads.gpd[1], means.gpd.r$Seq, pads.gpd[2])
means.zev.r$Seq = paste0(pads.zev[1], means.zev.r$Seq, pads.zev[2])
means.all = rbind(means.all, means.gpd.r, means.zev.r)

tables$Seqs[sapply(tables$Seqs, nchar) == 313] = paste0(pads.gpd[1], tables$Seqs[sapply(tables$Seqs, nchar) == 313], pads.gpd[2])
tables$Seqs[sapply(tables$Seqs, nchar) == 247] = paste0(pads.zev[1], tables$Seqs[sapply(tables$Seqs, nchar) == 247], pads.zev[2])

tables.use = data.table(Seq = tables$Seqs, Means_A_FS9 = tables$Means_A, Means_B_FS9 = tables$Means_B, Experiment.designed = tables$Experiment.designed)
means.use = data.table(Seq = means.all$Seq, Means_A_Orig = means.all$Strength_A, Means_B_Orig = means.all$Strength_B)
setkey(tables.use, Seq)
setkey(means.use, Seq)
table.final = merge(tables.use, means.use)

# change experiment names to match terminology used in final manuscript,
# and fit on figures
new.exp.names = c("GPD top outlier",
                  "GPD top inlier",
                  "GPD range",                            
                  "GPD bottom outlier",
                  "GPD test-data outlier",
                  "GPD test-active, good preds",        
                  "ZEV-U bottom outlier",
                  "ZEV-I bottom outlier",
                  "ZEV-U top outlier",            
                  "ZEV-I top inlier",
                  "ZEV AR high inlier",
                  "ZEV grid",                             
                  "ZEV-U test outlier",
                  "ZEV-I test outlier",
                  "ZEV AR test outlier",                
                  "ZEV-U test-active, good preds",
                  "ZEV-I test-active, good preds",
                  "ZEV AR test-active, good preds")

trim.name = function(name, refname = 'ZEV AR test outlier') {
  if(nchar(name) < nchar(refname)) {return(name)}
  name = substr(name,1,(nchar(refname)-3))
  name = paste0(name, '...', collapse = '')
  return(name)
}

new.exp.names = sapply(new.exp.names, trim.name)

setwd('D:/Promoter Design Data/Figures/PNGs')
exp.names.final = paste0(0:17, ': ', new.exp.names)
for(i in 0:17) {
  rs = which(table.final$Experiment.designed == i)
  table.tmp = copy(table.final)
  table.tmp = table.tmp[rs,]
  Original = c(table.tmp$Means_A_Orig, table.tmp$Means_B_Orig)
  Validation = c(table.tmp$Means_A_FS9, table.tmp$Means_B_FS9)
  Condition = rep('B', length(Original))
  Condition[1:(length(Condition)/2)] = rep('A', length(Condition)/2)
  table.tmp = data.frame(Original=Original, Validation=Validation, Condition=Condition)
  if(nrow(table.tmp) > 0) {
    png(filename = paste0('S13_', i, '.png', collapse=''),
        units = 'cm', width = 4, height = 4, res = 600)
    
      p = ggplot(data = table.tmp, aes(x = Original, y = Validation, color = Condition)) + geom_point(alpha = 1, size = 0.5, stroke = 0) + 
        theme_bw() + 
      xlim(c(-2,2)) + ylim(c(-2,2)) + 
      geom_hline(yintercept = -0.5, lty = 2, lwd = 0.5) + 
        theme(plot.title = element_text(hjust = 0.4, size = 8, face='bold')) +
        theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold'),
              legend.position = 'none') +
      labs(title = exp.names.final[i+1],
           x = 'Original (log10)',
           y = 'Validation (log10)')
    print(p)
  dev.off()
  }
}

