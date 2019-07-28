get.closest = function(target, means.use, range.reps) {
  diffs = abs(means.use$Strength - target)
  ord = order(diffs)
  return( means.use[ord,][1:range.reps,] )
}

# expect keys 'Seq','Strength_A','Strength_B'
get.outlier = function(dat.in, fn.out, nc, rs, get.top) {
  if(get.top) {
    dat.in$outlier = (dat.in$Strength_A == max(dat.in$Strength_A) & dat.in$Strength_B == max(dat.in$Strength_B))
  }
  else {
    dat.in$outlier = (dat.in$Strength_A == min(dat.in$Strength_A) & dat.in$Strength_B == min(dat.in$Strength_B))
  }
  set.seed(rs); rows.use = sample(which(dat.in$outlier), nc)
  seqs.use = dat.in$Seq[rows.use]
  scores.use = apply(cbind(dat.in$Strength_A, dat.in$Strength_B),1,mean)[rows.use]
  write.csv(data.frame(Seqs = seqs.use, Scores = scores.use), file = fn.out, row.names = FALSE, quote = FALSE)
}

# expect keys 'Seq', 'Strength'
sample.top = function(dat.in, fn.out, nc, rs) {
  # randomize dat.in to break alphabetization
  set.seed(rs); dat.in = dat.in[sample(nrow(dat.in)),]
  dat.in = tail(dat.in[order(dat.in$Strength),],nc)
  seqs.use = dat.in$Seq
  scores.use = dat.in$Strength
  write.csv(data.frame(Seqs = seqs.use, Scores = scores.use), file = fn.out, row.names = FALSE, quote = FALSE)
}

# expect keys 'Seq', 'Strength'
sample.range = function(dat.in, fn.out, nc, rs, range.reps) {
  set.seed(rs); dat.in = dat.in[sample(nrow(dat.in)),] # randomize to break alphabetization
  strength.targets = seq(from = min(dat.in$Strength), to = max(dat.in$Strength), length.out = nc/range.reps)
  closest.tables = sapply(strength.targets, get.closest, means.use = dat.in, range.reps = range.reps)
  seqs.use = unlist(closest.tables[which(names(dat.in) == 'Seq'),])
  scores.use = unlist(closest.tables[which(names(dat.in) == 'Strength'),])
  write.csv(data.frame(Seqs = seqs.use, Scores = scores.use), file = fn.out, row.names = FALSE, quote = FALSE)
}

# expect keys 'Pred_Strength','Strength', 'Seq'
# Get the 'nc' worst mis-predictions (no threshold)
get.prediction.outliers = function(dat.in, fn.out, nc, rs, outlier.type) {
  set.seed(rs); dat.in = dat.in[sample(nrow(dat.in)),]
  dat.in$Err = dat.in$Pred_Strength - dat.in$Strength
  dat.in = dat.in[order(dat.in$Err),] # smallest at top
  if(outlier.type == 'positive') { dat.in.use = tail(dat.in, nc) }
  else if(outlier.type == 'negative') { dat.in.use = head(dat.in, nc)}
  else { if(nc %% 2 != 0) { warning('Number of controls not divisible by 2') }
    dat.in.pos = tail(dat.in, floor(nc/2))
    dat.in.neg = head(dat.in, floor(nc/2))
    dat.in.use = rbind(dat.in.pos, dat.in.neg)
  }
  write.csv(data.frame(Seqs = dat.in.use$Seq, Scores = dat.in.use$Err), file = fn.out, row.names = FALSE, quote = FALSE)
}

# expect keys 'Pred_Strength','Strength', 'Seq'
# Get 'nc' mis-predictions worse than 'thresh', randomly chosen.
sample.prediction.outliers = function(dat.in, fn.out, nc, rs, thresh, outlier.type) {
  if(thresh < 0) {stop('Negative outlier threshold')}
  set.seed(rs); dat.in = dat.in[sample(nrow(dat.in)),]
  dat.in$Err = dat.in$Pred_Strength - dat.in$Strength
  if(outlier.type == 'positive') { 
    dat.in.use = dat.in[dat.in$Err > thresh,][1:nc,]
    }
  else if(outlier.type == 'negative') { 
    dat.in.use = dat.in[dat.in$Err < thresh,][1:nc,]
    }
  else { 
    dat.in.use = dat.in[abs(dat.in$Err) > thresh,][1:nc,]
  } 
  write.csv(data.frame(Seqs = dat.in.use$Seq, Scores = dat.in.use$Err), file = fn.out, row.names = FALSE, quote = FALSE)
}

# expect keys 'Pred_Strength','Strength', 'Seq'
# Get the 'nc' worst outliers, and then sample 'nc' other outliers worse than 'thresh'
# Always sample both positive and negative outliers
get.outliers.worst.and.range = function(dat.in, fn.out, nc, rs, thresh) {
  set.seed(rs); dat.in = dat.in[sample(nrow(dat.in)),]
  dat.in$Err = dat.in$Pred_Strength - dat.in$Strength
  dat.in$Abs.err = abs(dat.in$Err)
  dat.in = dat.in[rev(order(dat.in$Abs.err)),]
  dat.in.worst = dat.in[1:nc,]
  dat.in = dat.in[(nc+1):nrow(dat.in),]
  dat.in = dat.in[dat.in$Abs.err > thresh,]
  set.seed(rs); dat.in = dat.in[sample(nrow(dat.in),nc),]
  dat.in = rbind(dat.in.worst, dat.in)
  write.csv(data.frame(Seqs = dat.in$Seq, Scores = dat.in$Err), file = fn.out, row.names = FALSE, quote = FALSE)
}

vis = function(dat, xname, yname, ctrl.names, do.legend = TRUE, ...) {
  setnames(dat, c(xname, yname), c('x','y'))
  plot(dat$x, dat$y, pch = '.', xlab = xname, ylab = yname, ...)
  dat = dat[!is.na(dat$Last.Control.Use),]
  points(dat$x, dat$y, pch = 20, col = dat$Last.Control.Use + 1)
  if(do.legend) {
  ctrls.used = unique(dat$Last.Control.Use)
  ctrls.used = ctrls.used[!is.na(ctrls.used)]
  l.text = ctrl.names[ctrls.used + 1]
  legend('bottomright', legend = l.text, fill = ctrls.used + 1)
}}

vis.hist = function(dat, xname, ctrl.names, do.legend = TRUE, ...) {
  setnames(dat, xname, 'x')
  h = hist(dat$x, ...)
  dat = dat[!is.na(dat$Last.Control.Use),]
  points(dat$x, (0.8 + runif(length(dat$x), -0.02, 0.02))*max(h$counts), pch = 20, col = dat$Last.Control.Use + 1)
  if(do.legend) {
  ctrls.used = unique(dat$Last.Control.Use)
  ctrls.used = ctrls.used[!is.na(ctrls.used)]
  l.text = ctrl.names[ctrls.used + 1]
  legend('bottomright', legend = l.text, fill = ctrls.used + 1)
}}

sanitize.ctrl.names = function(dir.in) {
  n = dir(dir.in)
  san.one = function(x) {
    x = strsplit(x, '_')[[1]]
    x = x[2:length(x)]
    x = paste0(x, collapse = '_')
    x = strsplit(x, '.', fixed = TRUE)[[1]][1]
    return(x)
  }
  n = sapply(n, san.one)
  return(n)
}

# Individually plot each control set against the appropriate subset of data, on the correct axes.
summarize.ctrls = function(dat, ctrl.names, x.axes, y.axes, thresholds, panels, ...) {
  par(mfrow = panels)
  for(i in 1:length(ctrl.names)) {
    xname = x.axes[i]; yname = y.axes[i]
    thresh = thresholds[i]
    tmp = copy(dat)
    if(xname == yname) { setnames(tmp, c(toString(i-1), xname), c('this.control','x'))  }
    else {
      setnames(tmp, c(toString(i-1), xname, yname), c('this.control','x','y'))
    }
    tmp.libs = unique(tmp$Library[tmp$this.control]); tmp.libs = tmp.libs[!is.na(tmp.libs)]
    tmp.mods = unique(tmp$Model[tmp$this.control]); tmp.mods = tmp.mods[!is.na(tmp.mods)]
    tmp = tmp[(tmp$Library %in% tmp.libs & tmp$Model %in% tmp.mods),]
    if(xname == yname) {
      h = hist(tmp$x, xlab = xname, main = ctrl.names[i], ...)
      tmp = tmp[tmp$this.control,]
      points(tmp$x, (0.8 + runif(length(tmp$x), -0.02, 0.02))*max(h$counts), pch = 20, col = 'red')
    }
    else {
      plot(tmp$x, tmp$y, pch = '.', xlab = xname, ylab = yname, main = ctrl.names[i])
      tmp = tmp[tmp$this.control,]
      points(tmp$x, tmp$y, pch = 20, col = 'red')
      if(thresh != 0) {
      abline(a = thresh, b = 1, lty = 2, col = 'red')
      abline(a = -1*thresh, b = 1, lty = 2, col = 'red')
    }}
}}

  
  
  
