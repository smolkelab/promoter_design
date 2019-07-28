# Merge the predicted scores with actual scores (from FS9)
# for all designed sequences (designs + controls both).
# From the controls, figure out how FS9 results compare to original data
# (reproducible? Rescaling strengths?)
# From the experiments, figure out how accurate the model's predictions are
# in this extrapolated regime.

setwd('D:/Promoter Design Data/FACS-Seq/')
# for converting between padded and unpadded sequences later on
pads.gpd = c('TACGTAAATAATTAATAGTAGTGAC', 'TGTCTAAAGGTGAAGAATTATTCAC')
pads.zev = c('TTTATCATTATCAATACTCGCCATTTCAAAGAATACGTAAATAATTAATAGTAGTGAC', 'TGTCTAAAGGTGAAGAATTATTCACTGGTGTTGTCCCAATTTTGGTTGAATTAGATGG')

require(data.table)
fn = 'final_validation_FACS-Seq_means_just_designs.csv'
fs9.means = fread(fn)
fs9.means$V1 = NULL

# Generate a table like 'fs9_means' but with all designed sequences represented,
# whether measured or not; additionally, have column for predicted score
# and for actual score (e.g. for GPD this will be equal to mean strength, for ZEV AR, equal to AR)

get.table.by.idx = function(idx, dir.name) {
  fns = dir(dir.name)
  fn.tags = as.numeric(tstrsplit(fns,'_')[[1]])
  ans = fns[fn.tags == idx]
  # make sure exactly one output
  stopifnot(length(ans) == 1)
  return(fread(paste0(c(dir.name, '/', ans),collapse='')))
}

get.tables.in.dir = function(dir.name) {
  idxes = 1:length(dir(dir.name)) - 1 # 0-indexed
  tables = list()
  for(i in 1:length(idxes)) { tables[[i]] = get.table.by.idx(idxes[i], dir.name) }
  return(tables)
}

# Get designed sequences
setwd('D:/Promoter Design Data/Validation//')
tables.controls = get.tables.in.dir('controls')
tables.designs = get.tables.in.dir('new_promoters')

# tables.controls has a special case to consider - the ZEV grid has 'Scores_A' and 'Scores_B',
# and we want both.
# Merged table will have fields, 'Seqs', 'Experiment', 'Scores,'Scores_B' - 
# ZEV grid will have its 'Scores_A' in 'Scores', and will be the only sequence with actual values
# in 'Scores_B'. Rest will have NA.
id.grid = which(sapply(tables.controls, function(x) length(names(x))) > 2 )
setnames(tables.controls[[id.grid]], 'Scores_A', 'Scores')
for(i in 1:length(tables.controls)) {
  if(i != id.grid) { tables.controls[[i]]$Scores_B = NA}
  tables.controls[[i]]$Experiment = i-1
}
for(i in 1:length(tables.designs)) { 
  tables.designs[[i]]$Scores_B = NA
  tables.designs[[i]]$Experiment = i-1
}
tables.controls = rbindlist(tables.controls)
tables.designs = rbindlist(tables.designs)
# prepare to merge the Experiment fields
tables.designs$Experiment = tables.designs$Experiment + max(tables.controls$Experiment) + 1
tables = rbindlist(list(tables.controls, tables.designs))

# Now we need to merge in experimental results. Problem: 'tables$Seqs' is padded to 363;
# 'fs9.means$Seqs' is not (247 or 313).
# Solution: copy 'fs9.means$Seqs' to 'fs9.means$Seqs.orig'; create new 'fs9.means$Seqs'
# with padding.
setnames(fs9.means, 'Seqs', 'Seqs.orig')
new.seqs = vector(mode = 'character', length = nrow(fs9.means))
for(i in 1:nrow(fs9.means)) {
  seq = fs9.means$Seqs.orig[i]
  if(nchar(seq) == 313) {seq = paste0(c(pads.gpd[1], seq, pads.gpd[2]), collapse = '')}
  if(nchar(seq) == 247) {seq = paste0(c(pads.zev[1], seq, pads.zev[2]), collapse = '')}
  new.seqs[i] = seq
}
fs9.means$Seqs = new.seqs

setkey(fs9.means, Seqs)
setkey(tables, Seqs)

tables.merged = merge(fs9.means, tables, all=TRUE)
setnames(tables.merged, c('Experiment.x', 'Experiment.y'), c('Experiment.fs9means', 'Experiment.designed'))
# there are duplicate sequences in 'tables', but the merge appears to work as desired - 
# all possible matches appear to be covered
fn.key = 'designs_filename_map.csv'
x.key = fread(fn.key)
x.key$ID = x.key$ID + 18*(x.key$Type == 'D')
setnames(x.key, c('ID', 'Design'), c('Experiment.designed', 'Strategy'))
setkey(x.key, 'Experiment.designed')
setkey(tables.merged, 'Experiment.designed')
tables = merge(tables.merged, x.key, all = TRUE)

# Set measured score - again, ZEV grid requires special handling
setnames(tables, c('Scores', 'Scores_B'), c('Pred_Scores', 'Pred_Scores_B'))
set.measured.score = function(frame.in) {
  measured.score = vector(mode = 'numeric', length = nrow(frame.in))
  measured.score.B = vector(mode = 'numeric', length = nrow(frame.in))
  for(i in 1:nrow(frame.in)) {
  measured.score[i] = NA; measured.score.B[i] = NA
  handler = frame.in$Score_Handling[i]
    if(handler == 'AR') { measured.score[i] = frame.in$AR[i] }
    if(handler == 'Grid') { measured.score[i] = frame.in$Means_A[i]; measured.score.B[i] = frame.in$Means_B[i] }
    if(handler == 'Induced') { measured.score[i] = frame.in$Means_B[i] }
    if(handler == 'Strength') { measured.score[i] =  frame.in$Means_Avg[i]  }
    if(handler == 'Uninduced') { measured.score[i] = frame.in$Means_A[i] }
  }
  frame.in$Scores = measured.score
  frame.in$Scores_B = measured.score.B
  return(frame.in)
}

tables = set.measured.score(tables)

# Plot predicted score vs. measured score for each category
tables.notna = tables[!is.na(tables$Scores),]
u = unique(tables.notna$Experiment.designed) # length(u) = 48
par(mfrow = c(4,3))
for(i in 1:24) {
  tmp = tables.notna[tables.notna$Experiment.designed == u[i],]
  plot(tmp$Scores, tmp$Pred_Scores, pch = 20, main = u[i])
}

for(i in 25:48) {
  tmp = tables.notna[tables.notna$Experiment.designed == u[i],]
  plot(tmp$Scores, tmp$Pred_Scores, pch = 20, main = u[i])
}

# compare predicted vs. actual for ZEV-grid
par(mfrow = c(1,2))
tables.grid = tables[tables$Score_Handling == 'Grid',]
tables.grid = tables.grid[!is.na(tables.grid$Scores) & !is.na(tables.grid$Scores_B),]
plot(tables.grid$Scores, tables.grid$Pred_Scores)
abline(a = 0, b = 1, lty = 2, col = 'red')
plot(tables.grid$Scores_B, tables.grid$Pred_Scores_B)
abline(a = 0, b = 1, lty = 2, col = 'red')

# save the results
write.csv(tables, 'merged_FS9_results.csv')

