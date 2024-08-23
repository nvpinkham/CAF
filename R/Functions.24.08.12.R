medain4agg <- function(x){
  if(class(x) == "numeric" | class(x) == "integer"){
    return(median(x))
  }else{
    return(unique(x))
  }
}

plot.fam <- function(surv.map,
                     families, 
                     xlim = c(0, 90), 
                     lty.user = c(1, 1, 1, 1, 1, 2, 3, 1, 1, 1)){
  
  surv.map  <- surv.map[surv.map$family %in% families, ]

  if(class( surv.map$family) == "character"){
    surv.map$family <- factor(surv.map$family, levels = families)
  }
  
  surv <- Surv(surv.map$mouse.residency, surv.map$mouse.residency < 90)
  
  sur.tax <- survfit(surv ~ surv.map$family)
  mod.tax <- coxph(surv ~ surv.map$family)
  mod.aic <- AIC(mod.tax)
  
  sf <- fortify(sur.tax)
  leg <- aggregate(sf$n.risk, list(sf$strata), max)

  leg$col <- surv.map$fam.col[match(leg$Group.1, surv.map$family)]
  
  res <- summary(mod.tax)
  
  coef <- res$coefficients
  row.names(coef) <- sapply(strsplit(row.names(coef), "family"), "[[", 2)

  leg$hz <-  leg$hz.pretty <- NA  
  leg$p <-  leg$p.pretty <- NA  
  
  for(i in 1 : nrow(leg)){
    
    here <- which(row.names(coef) == leg$Group.1[i] )
    
    if(length(here) > 0){
      leg$hz[i] <- round(coef[here,2], 3)
      leg$p[i] <- round(coef[here,5], 3)
    }
  }
  

  leg$hz.pretty <- paste(", hazard ratio=", leg$hz)
  leg$hz.pretty[1] <- ", Baseline"
  leg$p.pretty <- paste(", p=", leg$p)
  leg$p.pretty[1] <- ""
  
  plot(sur.tax, col = leg$col, 
       xlim = xlim,
       lwd = 7, lty = lty.user, 
       font = 4,
       font.lab = 4,
       ylab = "Survival Probability", 
       xlab = paste("Days outside of Host\nModel AIC =", round(mod.aic)),
       main = paste( "\nCox model logtest p=", 
                    round(res$logtest[3], 3), 
                    ", Concordance=", round(res$concordance[1], 3)))
       
  
  legend("bottomright", 
         legend = paste0(leg$Group.1, " (n=", leg$x, leg$hz.pretty, leg$p.pretty, ")"), 
         fill =  leg$col,
         bty = "n", 
         text.font = 4,
         cex = 0.75)
  
  return(mod.tax)
}

plot.spore <- function(surv.map, var, colors, conf.int = T){
  
  var.name <- deparse(substitute(var))
  
  surv <- Surv(surv.map$mouse.residency, surv.map$mouse.residency < 90)
  
  cox.res <- coxph(surv ~ var)
  sur.fit <- survfit(surv ~ var)
  
  res <- summary(cox.res)
  
  mod.aic <- AIC(cox.res)
  
  sf <- fortify(sur.fit)
  leg <- aggregate(sf$n.risk, list(sf$strata), max)
  
  leg$col <- colors
  
  coef <- res$coefficients

  row.names(coef) <- gsub("var", "", row.names(coef))
  
  leg$hz <-  leg$hz.pretty <- NA  
  leg$p <-  leg$p.pretty <- NA  
  
  for(i in 1 : nrow(leg)){
    
    here <- which(row.names(coef) == leg$Group.1[i] )
    
    if(length(here) > 0){
      leg$hz[i] <- round(coef[here,2], 3)
      leg$p[i] <- round(coef[here,5], 3)
    }
  }
  
  leg$hz.pretty <- paste(", hazard ratio=", leg$hz)
  leg$hz.pretty[1] <- ", Baseline"
  leg$p.pretty <- paste(", p=", leg$p)
  leg$p.pretty[1] <- ""
  
  mod.aic <- AIC(cox.res)
  
  plot(sur.fit, 
       conf.int = conf.int,  
       main = paste0(# unique(surv.map$family),
         "\nCox model logtest p=", 
         round(res$logtest[3], 3), 
         paste(", Concordance=", round(res$concordance[1], 3))),
       xlab = paste("Days outside of Host\nModel AIC =", round(mod.aic)),
       ylab = "Survival Probability",
       font = 4,
       font.lab = 4,
       col = colors, lwd = 3)
  
  legend("bottomleft", 
         legend = paste0(leg$Group.1, " (n=", leg$x, leg$hz.pretty, leg$p.pretty, ")"), 
         fill =  colors,
         bty = "n", 
         text.font = 4,
         cex = 0.75)
  
  return(cox.res)
}

taxplot <- function (otu, tax, var, 
                     cutoff = 0.05, log = F) {
  par(mar = c(8, 4, 4, 2))
  if (log) {
    otu <- log(otu + 1)
    label.y <- "log transformed proportion of community"
  }
  else {
    label.y <- "Proportion of Community"
  }
  otu <- otu/rowSums(otu)
  
  tax.agg <- aggregate(t(otu), list(tax), sum)
  rownames(tax.agg) <- tax.agg$Group.1
  tax.agg$Group.1 <- NULL
  
  tax.agg <- aggregate(t(tax.agg), list(var), mean)
  rownames(tax.agg) <- tax.agg$Group.1
  tax.agg$Group.1 <- NULL
  tax.agg <- t(tax.agg)
  tax.agg <- tax.agg[order(rowSums(tax.agg), decreasing = T), 
  ]
  tax.agg.pick <- tax.agg[rowSums(tax.agg)/sum(tax.agg) >= 
                            cutoff, ]
  tax.agg.other <- tax.agg[rowSums(tax.agg)/sum(tax.agg) < 
                             cutoff, ]
  other <- colSums(tax.agg.other)
  tax.agg.pick <- rbind(tax.agg.pick, other)
  color = grDevices::colors()[grep("gr(a|e)y", grDevices::colors(), 
                                   invert = T)]
  set.seed(145941)
  
  color <-c("pink1", "green", "mediumpurple1", "slateblue1", "gold", "orchid",
            "turquoise2", "skyblue", "steelblue","tan2", "navyblue",
            "orange", "orangered", "coral2", "palevioletred", "violetred", "darkred",
            "springgreen2", "yellowgreen", "palegreen4",
            "wheat2", "tan", "magenta", "tan3", "brown",
            "grey70", "grey50", "blue")
  
  cols <- sample(color, nrow(tax.agg.pick))
  col.key <- cbind(rownames(tax.agg.pick), cols)
  col.key[nrow(col.key), 2] <- "grey50"
  barplot(tax.agg.pick,
          width = 1,
          beside = F, 
          xlim = c(0, ncol(tax.agg.pick) * 1.75),
          space = 0.1,
          col = col.key[, 2],
          las = 2, 
          ylab = label.y, 
          font = 4,
          font.lab = 4,
          cex.axis=1,
          cex.names=1)
  
  legend("right",
         bty = "n",
         text.font = 4,
         fill = rev(col.key[, 2]), 
         legend = rev(col.key[, 1]))
  
  return(col.key)
}

rowMedians <- function(x){
  apply(x, 1, median)
}

tax.shared.table <- function (otu, var, group1, group2, tax.level = tax$family, 
                              PresAbs = F, log = T, taxa2include, round.to = 3, 
                              add.n = T) {
  if (PresAbs & log) {
    warning("Log transformation not informative when using PresAbs")
  }
  var <- as.factor(var)
  group1.otus <- otu[var == group1, ]
  group2.otus <- otu[var == group2, ]
  if (PresAbs) {
    gg <- "number of discrete OTUs"
    group1.otus[group1.otus > 0] <- 1
    group2.otus[group2.otus > 0] <- 1
  }else{
    gg <- ""
  }
  group1.otus.shared <- group1.otus.unique <- group1.otus
  group2.otus.shared <- group2.otus.unique <- group2.otus
  group1.otus.unique[, colSums(group2.otus) > 0] <- 0
  group1.otus.shared[, colSums(group1.otus.unique) > 0] <- 0
  group2.otus.unique[, colSums(group1.otus) > 0] <- 0
  group2.otus.shared[, colSums(group2.otus.unique) > 0] <- 0
  
  group1.fam.unique <- aggregate(t(group1.otus.unique), list(tax.level),  sum)
  group1.fam.shared <- aggregate(t(group1.otus.shared), list(tax.level),  sum)
  group2.fam.unique <- aggregate(t(group2.otus.unique), list(tax.level), sum)
  group2.fam.shared <- aggregate(t(group2.otus.shared), list(tax.level), sum)
  
  g1.u <- rowMedians(group1.fam.unique[, -1])
  names(g1.u) <- group1.fam.unique$Group.1
  g1.s <- rowMedians(group1.fam.shared[, -1])
  names(g1.s) <- group1.fam.shared$Group.1
  g2.u <- rowMedians(group2.fam.unique[, -1])
  names(g2.u) <- group2.fam.unique$Group.1
  g2.s <- rowMedians(group2.fam.shared[, -1])
  names(g2.s) <- group2.fam.shared$Group.1
  
  ps.unique.l <- ps.unique.g <- ps.shared <- NULL
  
  for (i in 1:nrow(group1.fam.unique)) {
    ps.unique.l[i] <- wilcox.test(t(group1.fam.unique[i, -1]), 
                                 t(group2.fam.unique[i, -1]), alternative = "less", exact = F)$p.value
    ps.shared[i] <- wilcox.test(t(group1.fam.shared[i, -1]), 
                                t(group2.fam.shared[i, -1]), exact = F)$p.value
    ps.unique.g[i] <- wilcox.test(t(group1.fam.unique[i, -1]), 
                                 t(group2.fam.unique[i, -1]), alternative = "greater", exact = F)$p.value
  }
  ps.unique.l <- p.adjust(ps.unique.l, method = "fdr")
  ps.shared <- p.adjust(ps.shared, method = "fdr")
  ps.unique.g <- p.adjust(ps.unique.g, method = "fdr")
  
  ps.unique.l <- round(ps.unique.l, round.to)
  ps.shared <- round(ps.shared, round.to)
  ps.unique.g <- round(ps.unique.g, round.to)

  hh <- which(ps.unique.g > 0.05 | is.na(ps.unique.g))
  ps.unique.g[hh] <- " "
  ps.unique.g[-hh] <- paste0("(unique greater p=", ps.unique.g[-hh], ")")
  
  kk <- which(ps.shared > 0.05 | is.na(ps.shared))
  ps.shared[kk] <- " "
  ps.shared[-kk] <- paste0("(shared p=", ps.shared[-kk], ")")
  
  jj <- which(ps.unique.l > 0.05 | is.na(ps.unique.l))
  ps.unique.l[jj] <- " "
  ps.unique.l[-jj] <- paste0("(unique less p=", ps.unique.l[-jj], ")")
  

  tax.res <- rbind(g1.u, g1.s, g2.s, g2.u)
  
  if(add.n){
    otu.pick <- otu[var %in% c(group1, group2), ]
    ns <- table(tax.level[ colSums(otu.pick) > 0]) # number of taxa in both groups
    colnames(tax.res)[match(names(ns), colnames(tax.res))] <- paste0(names(ns), "(nOTUs=", ns, ")")
  }
  
  colnames(tax.res) <- paste0(colnames(tax.res), ps.unique.l, 
                              ps.shared, ps.unique.g)
  rownames(tax.res) <- c(paste0("unique to", group1), 
                         paste0("shared - ", gg, "in", group1), 
                         paste0("shared - ", gg, "in", group2), 
                         paste0("unique to ", group2))
  
  tax.res <- tax.res[, colSums(tax.res) != 0]
  
  
  if (log) {
    tax.res <- log10(as.matrix(tax.res) + 1)
  }
  
  tax.res <- tax.res[, order(colSums(tax.res), decreasing = F)]
  
  if (!missing(taxa2include)) {
    print(taxa2include)
    tax.res2 <- as.data.frame(matrix(nrow = 4, ncol = length(taxa2include)))
    row.names(tax.res2) <- row.names(tax.res)
    for (i in 1:length(taxa2include)) {
      aa <- grep(taxa2include[i], colnames(tax.res))
      if (length(aa) > 0) {
        tax.res2[, i] <- tax.res[, aa]
        colnames(tax.res2)[i] <- colnames(tax.res)[aa]
      }
      else {
        colnames(tax.res2)[i] <- taxa2include[i]
        tax.res2[, i] <- 0
      }
    }
    tax.res <- as.matrix(tax.res2)
  }
  
  colnames(tax.res) <- gsub("(", " (", colnames(tax.res), fixed = T)
  
  return(tax.res)
}

tax.shared <- function (otu, var, group1, group2, tax.level = tax$family, PresAbs = F, 
                        log = T, taxa2include) {
  tax.res <- tax.shared.table(otu, var, group1, group2, tax.level, 
                              PresAbs, log, taxa2include)
  par(mar = c(3, 20, 4, 4))
  par(xpd = TRUE)
  tax.res <- tax.res[, order(colSums(tax.res), decreasing = F)]
  here <- ceiling(max(colSums(tax.res))) - 0.5
  if (!missing(taxa2include)) {
    print(taxa2include)
    tax.res2 <- as.data.frame(matrix(nrow = 4, ncol = length(taxa2include)))
    row.names(tax.res2) <- row.names(tax.res)
    for (i in 1:length(taxa2include)) {
      aa <- grep(taxa2include[i], colnames(tax.res))
      if (length(aa) > 0) {
        tax.res2[, i] <- tax.res[, aa]
        colnames(tax.res2)[i] <- colnames(tax.res)[aa]
      }
      else {
        colnames(tax.res2)[i] <- taxa2include[i]
        tax.res2[, i] <- 0
      }
    }
    tax.res <- as.matrix(tax.res2)
  }
  
  a <- colnames(tax.res)
  b <- paste(colnames(tax.res), "dummy")
  
  while(all(a != b)){
    a <- colnames(tax.res) <- gsub(") ", ")", colnames(tax.res))
    b <- gsub(") ", ")", a)
  }
  labs <- gsub(")(", ") (", b, fixed = T)
  
  barplot(tax.res,
          cex.axis = 0.75,
          cex.names = 0.75,
          font.lab = 4,
          xaxt = "n", 
          yaxt = "n", 
          horiz = T, 
          main = "", 
          xlab = "median population size (log 10)", 
          las = 1, 
          col = c("orange", "#b38a61", "#887382", "purple"))
  
  axis(2,
       at= (1:length(labs) * 1.2) - .5,
       labels=labs, 
       las = 2, 
       cex.axis=.5,
       font = 4)
  
  
  if (log ) {
    
    segments(x0 = here - log10(10), x1 = here, y0 = 14, y1 = 14,  lwd = 2)
    text(here - 0.5, 14, "Counts\n10\n\n", cex = 0.5, font = 4)
    
    segments(x0 = here - log10(100), x1 = here, y0 = 13,  y1 = 13, lwd = 2)
    text(here - 0.5, 13, "100\n", cex = 0.5, font = 4)
    
    segments(x0 = here - 3, x1 = here, y0 = 12, y1 = 12, lwd = 2)
    text(here - 0.5, 12, "1000\n", cex = 0.5, font = 4)
  }
  legend("bottomright",
         fill = c("orange", "#b38a61", "#887382", "purple"),
         legend = row.names(tax.res), 
         bty = "n", 
         text.font = 4,
         cex = 0.75)
  return(tax.res)
}



get.blast.alignment <- function(protein, genome){
  
  protein.path <- paste0("Spore_Proteins_representatives/", protein, ".fasta")
  
  blast.1 <- blast.FUNction(query.fasta = protein.path, 
                            subject.fasta = genome,
                            algorithm = "tblastn",
                            outfmt = 6, 
                            max_tar = 1,
                            num_threads = 16)
  
  pick <- which(blast.1$bitscore == max(blast.1$bitscore))[1]
  blast.1 <- blast.1[pick , ]
  aa.pick <- blast.1$qseqid
  
  protien.aa <- microseq::readFasta(protein.path)
  protien.aa <- protien.aa[protien.aa$Header ==  aa.pick, ]
  
  blast.1$q.length <- nchar(protien.aa$Sequence)
  
  microseq::writeFasta(protien.aa, out.file = "temp.fasta")
  
  blast.2 <- blast.FUNction(query.fasta = "temp.fasta", 
                            subject.fasta = genome,
                            algorithm = "tblastn",
                            out = paste0(protein, "_blast_align.txt"),
                            outfmt = 1, 
                            max_tar = 1,
                            num_threads = 16)
  return(blast.1)
}

get.group.summary <- function(list = diffs, p_ident = p_ident, genome.paths = genome.paths){
  
  result <- NULL
  for(i in 1 : length(list)){
    print(i)
    
    for(j in list[[i]]){
      
      protein.i.j <- colnames(p_ident)[j]
      
      print(protein.i.j)
      
      genome.i <- genome.paths[grep(row.names(p_ident)[i], genome.paths)]
      
      res <- get.blast.alignment(protein = protein.i.j, genome = genome.i)
      
      result <- rbind(result, res)
    }
  }
  return(result)
}

plot.otus1 <- function(human.info = otu.info$p3.human ,  
                       mice.info = otu.info$otu.m3.p3 , 
                       groups = c("resident", "migrant")){
  # groups = c("resident", "migrant", "trans")){
  
  
  tab <- table(human.info[human.info %in% groups], 
               mice.info[human.info %in% groups] > 0)
  
  res <- fisher.test(tab)
  
  try(
    res <- paste("\n\nFisher Test p value =", round(res$p.value, 5), 
                 ", odds ratio = ", round(res$estimate, 5))
  )
  
  plot(tab, col = c(8, 3), ylab = "Recolonized", cex = 1, main = "")
  title(res)
  
  return(res)
}

rowMaxs <- function(x){
  return(apply(x, 1, max))
}


remove_rare <- function( table , cutoff_pro = .2) {
  # cutoff is the proportion of samples that a feature must be found in. 
  table <- t(table)# transpose to keep "tidy" ; easier than rewriting function...
  row2keep <- c()
  cutoff <- ceiling( cutoff_pro * ncol(table) )
  for ( i in 1:nrow(table) ) {
    row_nonzero <- length( which( table[ i , ]  > 0 ) )
    if ( row_nonzero > cutoff ) {
      row2keep <- c( row2keep , i)
    }
  }
  res <- t(table [ row2keep , , drop=F ])
  res <- as.data.frame(res)
  
  return(res )
}

pairwiseRF <- function(mat, 
                       map, 
                       dis_1 = "T 0", 
                       dis_2 = "B 7",
                       num_of_col = 50,
                       Accuray.cutoff = 0,
                       Factor = 1, 
                       plot = T){
  # returns features of a matrix that differentiate between two treatments based on RF analysis 
  # dis_1 and dis_2 are the treatment levels being compared
  # map needs a column named  "Description"  
  set.seed(42)
  mat.pick <- mat[map$Description %in% c(dis_1, dis_2 ), ]
  map.pick <- map[map$Description %in% c(dis_1, dis_2 ), ]
  
  #mat.pick$Group <- NULL
  
  mat.pick.rare_removed <- remove_rare(table = mat.pick, cutoff_pro = 0.2)
  
  spliter <- as.factor(as.character(map.pick$Description))
  
  fit <- randomForest::randomForest(x = mat.pick.rare_removed,
                                    y = spliter,
                                    importance=TRUE,
                                    proximities=TRUE,
                                    ntree = 5001)
  RF.sig <- rfUtilities::rf.significance(x = fit,
                                         xdata = mat.pick.rare_removed,
                                         nperm = 101,
                                         ntree = 501)
  print(fit)
  print(RF.sig)
  
  if(RF.sig$pValue < 0.05){
    
    fit_imp <- as.data.frame( fit$importance )
    fit_imp$features <- rownames(fit_imp )
    fit_imp_sorted <- dplyr::arrange( fit_imp  , desc(fit_imp$MeanDecreaseGini)  )
    colors <- vector(length = ncol(mat.pick.rare_removed))
    
    for(j in 1:ncol(mat.pick.rare_removed)){
      i <- fit_imp_sorted$features[j]
      t1.mean <- mean(mat.pick[map.pick$Description == dis_1, which(colnames(mat.pick) == i)])
      t2.mean <- mean(mat.pick[map.pick$Description == dis_2, which(colnames(mat.pick) == i)])
      if( t1.mean >  t2.mean){
        colors[j] <- "cadetblue3"
      }else{
        colors[j] <- "orange3"
      }
    }
    
    if(plot){
      barplot(fit_imp_sorted$MeanDecreaseGini[1:num_of_col],
              names.arg = fit_imp_sorted[1:num_of_col,"features"],
              ylab= "Mean Decrease in Accuracy (Variable Importance)",
              las = 2,
              col = colors,
              main= paste(dis_1,"or", dis_2, "Classification RF"))
      
      legend("topright",
             legend = paste("elivated in", c(dis_1, dis_2)),
             fill = c("cadetblue3", "orange3"))
      
    }
    res <- fit_imp_sorted$MeanDecreaseGini
    names(res) <- fit_imp_sorted$features
    
    if(plot){
      abline(h = median(res) + sd(res) * Factor,
             col = 2, lty = 2)
      text(50, median(res) + sd(res) * Factor + .01, paste(" median Feature impotance +", Factor, "sd"),
           col = 2, font = 2)
    }
    res.pick <- res[res > median(res) + (sd(res) * Factor)]
    colors.pick <- colors[res > median(res) + (sd(res) * Factor)]
    
    mat.bloom <- res.pick[colors.pick == "orange3"]
    mat.bust <- res.pick[colors.pick == "cadetblue3"]
    
    Result <- list(  mat.bust, mat.bloom)
    names(Result) <- c(dis_1, dis_2)
    
    return(Result)
  }else{
    return("VARIABLE p val > 0.05")
  }
}

add.line <- function(y = line, x = seq(2, 21, 3)){
  for(i in 1 : (length(y) - 1)){
    segments(x0 = x[i], y0 = y[i], x1 = x[i + 1], y1 = y[i+1], col = 2)
  }
}  

TemporalPerm <- function(otu = otu.mice, 
                         time = map.mice$days.between.hosts, 
                         groups = map.mice$sample.type){
  
  p <- !is.na(time) & !is.na(groups)
  
  time <- time[p]
  groups <- groups[p]
  otu <- otu[p,]
  
  times <- sort(unique(time))
  
  res <- matrix(nrow =  length(times), ncol = 2)
  
  for(i in 1 : length(times)){
    
    p.i <- time == times[i]
    
    res.i <- adonis2(otu[p.i,] ~ groups[p.i])
    res[i,] <-c(res.i$`Pr(>F)`[1], res.i$F[1])
  }
  colnames(res) <- c("p.value", "f.stat")
  
  return(res)
}

TemporalT <- function(var = map.mice$invsimp, 
                      time = map.mice$days.between.hosts, 
                      groups = map.mice$sample.type){
  
  p <- !is.na(time) & !is.na(groups) & !is.na(var)
  
  time <- time[p]
  groups <- groups[p]
  var <- var[p]
  
  times <- sort(unique(time))
  
  res <- matrix(nrow =  length(times), ncol = 2)
  for(i in 1 : length(times)){
    
    p.i <- time == times[i]
    
    res.i <- t.test(var[p.i] ~ groups[p.i])
    res[i,] <- c(res.i$p.value, diff(res.i$estimate))
  }
  colnames(res) <- c("p.value", "effect.size")
  return(res)
}

resident <- function(otu = otu.pick, dates, min.res = 14, max.abs = 30){
  # determines if an OTU is a resident in a temporal dataset 
  
  pres <- as.Date("1991-01-10")
  last.resident <- as.Date("1991-01-10")
  resident <- NULL
  
  otu.bi <- otu > 0
  
  for(i in 1 : length(otu.bi)){
    
    if(otu.bi[i]){
      
      # if the OTU is present 🎁 begin counting for how long 
      pres <- c(pres, dates[i])
      
      if(diff(range(pres)) >= min.res){
        
        last.resident <- tail(pres, 1) # keep track of when it was called resident
        resident[dates >= min(pres) & dates <= max(pres)] <- T # call true ever since the pres streak began
        
      }else{
        # not around long enough
        # will be changed if pres streak goes longer than min.res
        resident[i] <- F
      }
      
    }else{
      # OTU not seen
      time.since.last <- diff(c(last.resident, dates[i]))# how long since it was called a resident
      
      if(time.since.last <= max.abs){
        # within window still called a res, for now...
        resident[i] <- F
        
      }else{
        # not a resident or no longer a redisent 
        pres <- as.Date(character(0))# reset Present 🎁 counter
        resident[i] <- F
        resident[dates > last.resident] <- F # use last date seen to stop the residency counter
      }
    }
  }
  
  return(resident)
}

classify.residency <- function(otu.pick, 
                               dates = map.i$Date.Collected, 
                               cutoffs = c(60, 300), 
                               min.res = 14){
  # map.p must be for just one participant 
  
  resident.res <- resident(otu = otu.pick, 
                           dates = dates,
                           min.res = min.res,
                           max.abs = 30)
  
  time.elapsed <- sapply(dates, function(x) x - dates[1])
  
  resident.rle <- rle(resident.res)
  resident.rle2 <- rle(resident.res[time.elapsed >= 14])
  
  if(all(resident.rle2$values)){
    
    invasions <- 0
  }else{
    
    aa <- (which(!resident.rle2$values) + 1) %in% which(resident.rle2$values)
    invasions <- length(aa[aa])
  }
  
  beginnings <- c(1,1+which(diff(resident.res)!=0))
  ends <- c(which(diff(resident.res)!=0), length(resident.res))
  
  days <- NULL
  
  for(i in 1 : length(resident.rle$lengths)){
    
    days[i] <- dates[ends[i]] - dates[beginnings[i]]
  }
  
  if(!all(resident.res == F)){
    
    heres <- days[resident.rle$values == T]
    
    max.residency <- max(heres, na.rm = T)
    where.max <- which(max.residency == days & resident.rle$values)
    
    censor <- "no"
    if(max(where.max) == 1){
      censor <- "left"
    }
    if(min(where.max) == length(days)){
      censor <- "right"
    }
    if(min(where.max) == length(days) & max(where.max) == 1){
      censor <- "both"
    }
    
  }else{
    max.residency <- 0
    censor <- "no"
  }
  
  ######## BIN OTU RESIDENCY ########
  if(any(resident.res)){
    migration <- "Resident_1"
  }
  if(!any(resident.res) & max(otu.pick > 0)){
    migration <- "Transient"
  }
  #  if(!any(resident.res) & otu.pick[day.pick] == 0){
  #    migration <- "Transient_at_humanization"
  # }
  for(i in 1 : length(cutoffs)){
    # assign resdient to the highest level of cutoff
    if(max.residency >= cutoffs[i]){
      # print(paste("resident for", cutoffs[i]))
      migration <- paste0("Resident_", i + 1)
    }
  }
  if(all( otu.pick == 0)){
    migration <- "Missing in human"
  }
  
  
  res <- data.frame(migration, max.residency, censor)
  
  colnames(res) <- c("migration", "residency", "censor")
  
  return(res)
}

classify.residency.by.sample <- function(otu.pick, 
                                         dates = map.i$Date.Collected, 
                                         date.pick = "2017-03-07",
                                         cutoffs = c(50 ,100, 300), 
                                         min.res = 14){
  # map.p must be for just one participant 
  
  resident.res <- resident(otu = otu.pick, 
                           dates = dates,
                           min.res = min.res,
                           max.abs = 30)
  
 # cbind(resident.res, otu.pick)
  
  day.pick <- which(dates == as.Date(date.pick))
  
  time.elapsed <- sapply(dates, function(x) x - dates[1])
  
  resident.rle <- rle(resident.res)
  
  beginnings <- c(1, 1+which(diff(resident.res)!=0))
  beginnings <- beginnings[resident.rle$values]
  
  ends <- c(which(diff(resident.res)!=0), length(resident.res))
  ends <- ends[resident.rle$values]
  
  beginnings.pick <- rev(beginnings[which(beginnings <= day.pick)])[1]
  ends.pick <- ends[which(ends > day.pick)][1]
  
  otu.dummy <- otu.pick
  
  if(!is.na(beginnings.pick) & !is.na(ends.pick)){
    # Run Again with other colonization streaks blocked out
    block <- beginnings.pick:ends.pick
    otu.dummy[ ! 1:length(otu.dummy) %in% block ] <- 0
    
  }else{
    # set every day but the day in the sample to 0, this is to distinguish between trans and missing
    otu.dummy[ -day.pick] <- 0
  }
  
  resident.res <- resident(otu = otu.dummy, 
                           dates = dates,
                           min.res = min.res,
                           max.abs = 30)
  
  resident.rle <- rle(resident.res)
  
  beginnings <- c(1, 1+which(diff(resident.res)!=0))
  ends <- c(which(diff(resident.res)!=0), length(resident.res))
  
  days <- NULL
  
  for(i in 1 : length(resident.rle$lengths)){
    days[i] <- dates[ends[i]] - dates[beginnings[i]]
  }
  
  if(!all(resident.res == F)){
    
    heres <- days[resident.rle$values == T]
    
    residency <- max(heres, na.rm = T)
    where.max <- which(residency == days & resident.rle$values)
    
    censor <- "no"
    if(max(where.max) == 1){
      censor <- "left"
    }
    if(min(where.max) == length(days)){
      censor <- "right"
    }
    if(min(where.max) == length(days) & max(where.max) == 1){
      censor <- "both"
    }
    
  }else{
    residency <- 0
    censor <- "no"
  }
  
  ######## BIN OTU RESIDENCY ########
  if(any(resident.res)){
    migration <- "Resident_1"
  }
  if(!any(resident.res) & max(otu.pick > 0)){
    migration <- "Transient"
  }
#  if(!any(resident.res) & otu.pick[day.pick] == 0){
#    migration <- "Transient_at_humanization"
 # }
  for(i in 1 : length(cutoffs)){
    # assign resdient to the highest level of cutoff
    if(residency >= cutoffs[i]){
      # print(paste("resident for", cutoffs[i]))
      migration <- paste0("Resident_", i + 1)
    }
  }
  if(all( otu.pick == 0)){
    migration <- "Missing in human"
  }
  
  
  human.day <- otu.pick[day.pick]
  
  res <- data.frame(migration, residency, censor, human.day)
  
  colnames(res) <- c("migration", "residency", "censor", "count.humanize")
  
  return(res)
}

max.life <- function(time = otu.surv$days.between.hosts, event = otu.surv$recip.event){
  
  time.alive <- time[event == 0]
  a <- max(c(time.alive, 0))
  print(a)
  return(which(time == a))
}

plot.kap <- function(mod){
  survfit2(mod) %>% 
    ggsurvfit() +
    labs(
      x = "Days",
      y = "Overall survival probability by residency"
    ) +  
    add_confidence_interval()+
    add_risktable(times = c(0, 1, 7, 30, 90))+
    add_quantile(y_value = .5,
                 color = "gray50", 
                 linewidth = 0.75)
}

otu.summarize.samples <- function(otu, 
                                  map, 
                                  cutoffs = 100,
                                  by.site = T, 
                                  alsways.there = T, 
                                  use.median = T, 
                                  samples.used = samples.used){
  # residency calculate based on OTUs present in the samples used to humanize donor mice. 
  #   if an OTU recolonized a human donor multiple times. use the residency time containing the date used for humanization
  
  # by.site allows you to specify different cutoff by site/donor/individual - lust use list
  
  # alsways.there as in alwas a resident in the human to define resident_2 
  
  # use.median if median is in GF recipiant than the event occured, else max is used (only one opservation needed to say event has NOT occured)
  
  all(row.names(otu) == map$sample.id)
  
  otu.human <- otu[map$group == "human" , ]
  map.human <- map[map$group == "human" , ]
  otu.human <- otu.human[order(map.human$Date.Collected) , ]
  map.human <- map.human[order(map.human$Date.Collected) , ]
  
  map.mice <- map[map$group != "human" , ]
  otu.mice <- otu[map$group != "human" , ]
  
  participants <- sort(unique(map$participant))
  otus <- colnames(otu.mice)
  otu.sum <- NULL
  
  for(i in 1 : length(participants)){
    
    print(paste("Calculating summary for participant", i))
    
    otu.i <- otu.human[map.human$participant == participants[i] , ]
    map.i <- map.human[map.human$participant == participants[i] , ]
    otu.mice.i <- otu.mice[map.mice$participant == participants[i] ,]
    map.mice.i <- map.mice[map.mice$participant == participants[i] , ]
    
    sample_humanization <- samples.used[names(samples.used) == participants[i]]
    date4humanization <- map.i$Date.Collected[map.i$sample.id == sample_humanization]
    
    if(length(date4humanization) == 0){
      date4humanization <- map.i$Date.Collected[1]
    }
    
    cutoffs.using <- NULL
    
    if(by.site & alsways.there){
      cutoffs.using <- diff(range(map.i$Date.Collected))
    }
    
    else if(by.site & !alsways.there){
      cutoffs.using <- cutoffs[[i]]
    }
    
    else if(!by.site & !alsways.there){
      cutoffs.using <- cutoffs
      
    }else{
      
      stop("alsways.there only availbe if by site = T")
    }
    
    
    for(j in 1 : ncol(otu.i)){
      
      ###################### OTUs in the HUMANs ########################################
      otu.sum.i <- classify.residency.by.sample(otu.pick = otu.i[,j], 
                                                dates = map.i$Date.Collected, 
                                                date.pick = date4humanization,
                                                cutoffs = cutoffs.using)
      otu.sum.i$participant <- participants[i]
      otu.sum.i$otu <- colnames(otu.i)[j]
      otu.sum.i$otu.id <- paste0(otu.sum.i$participant,".", otu.sum.i$otu)
     # otu.sum.i$family <- tax.df$family[j]
      
      if(nrow(otu.mice.i) == 0){
        # print("no mouse info")
      }else{
        ###################### OTUs in the MICE ########################################
        
        donor.pop <- otu.mice[which(map.mice$participant == participants[i] &
                                      map.mice$group == "donor" ) , j]
        
        #otu.sum.i$in.donor.mouse <- max(donor.pop) > 0
        otu.sum.i$in.donor.mouse <- median(donor.pop) > 0      
        
        otu.sum.i$mouse.max <- max(otu.mice.i[map.mice.i$group == "donor" , j])
        otu.sum.i$mouse.mean <- mean(otu.mice.i[map.mice.i$group == "donor" , j])
        
        map.mice.i.recip <- map.mice.i[map.mice.i$group == "recipiant" & map.mice.i$sample.type == "scruff" , ]
        otu.mice.i.recip <- otu.mice.i[map.mice.i$group == "recipiant" & map.mice.i$sample.type == "scruff" , ]
        # in recip mouse scruff
        if(use.median){
          otu.agg <- aggregate(otu.mice.i.recip[, j], 
                               by = list(map.mice.i.recip$days.between.hosts), 
                               median)
        }else{
          otu.agg <- aggregate(otu.mice.i.recip[, j], 
                               by = list(map.mice.i.recip$days.between.hosts), 
                               max)
        }
        
        otu.agg$x <- otu.agg$x > 0
        otu.agg[1,2] <- T # call 0 true since the event (loss of viability) occurred at day 0
        
        max.viable <- otu.agg$Group.1[max(which(otu.agg$x))]
        otu.sum.i$mouse.residency <- max.viable
      }
      otu.sum <-merge(otu.sum, otu.sum.i, all = T)
    }
  }
  ################################################################################
  otu.sum$status <- otu.sum$mouse.residency == 90
  
  otu.sum$human <- "missing from humanization sample"
  otu.sum$mouse <- "missing from donor mouse"
  
  otu.sum$human[which(otu.sum$count.humanize > 0)] <- "in humanization sample"
  otu.sum$mouse[otu.sum$in.donor.mouse] <- "in donor mouse"
  
  otu.sum$censor2 <- otu.sum$censor
  otu.sum$censor2[otu.sum$censor == "left" | otu.sum$censor == "right"] <- "one side"
  otu.sum$censor2 <- factor(otu.sum$censor2, levels = c("both", "one side", "no"))
  
  otu.sum$migration.censor <- paste(otu.sum$migration, otu.sum$censor2)# don't use "+" because not all levels are meaningful
  
  return(otu.sum)
}

add.bins <- function(map, otu, otu.sum, otu.groups = "migration.censor"){
  
  var <- otu.sum[ , colnames(otu.sum) == otu.groups ]
  
#  map <- map[map$participant %in% otu.sum$participant , ]
  
  bins <- factor(var)
  bins <- levels(bins)
  
  for(i in 1 : length(bins)){
    map <- cbind(NA, map)
  }
  
  colnames(map)[1: length(bins)] <- bins
  
  for(i in 1 : nrow(map)){
    
    var.i <- var[otu.sum$participant == map$participant[i]]
    
    sample <- as.numeric(t(otu[i,]))
    
    if(length(var.i) > 0){
      
      bins.i <- aggregate(sample, list(var.i), sum)
      bins.i$x <-   bins.i$x / sum(bins.i$x)
      
      map[i , match(bins.i$Group.1, colnames(map)) ] <- bins.i$x
      
    }else{
      
      bins.i <- vector(length = length(bins))
      map[i , 1 : length(bins) ] <- NA
    }
  }
  
  return(map)
}

add.bins.bi <- function(map, otu, otu.sum, otu.groups = "migration.censor"){
  
  var <- otu.sum[ , colnames(otu.sum) == otu.groups ]
  
  bins <- factor(var)
  bins <- levels(bins)
  otu <- otu > 0
  
  for(i in 1 : length(bins)){
    map <- cbind(NA, map)
  }
  
  colnames(map)[1: length(bins)] <- bins
  
  for(i in 1 : nrow(map)){
    
    var.i <- var[otu.sum$participant == map$participant[i]]
    
    sample <- as.numeric(t(otu[i,]))
    bins.i <- aggregate(sample, list(var.i), sum)
    
    bins.i$x <-   bins.i$x / sum(bins.i$x)
    
    map[i , match(bins.i$Group.1, colnames(map)) ] <- bins.i$x
  }
  
  return(map)
}

color_palette <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#9370DB",
                   "#F0E442", "#56B4E9", "#FFAB30", "#6ACC64", "#C44E52",
                   "#B47CC7", "#FFFF99", "#00BFC4", "#9A82B8", "#8DB255",
                   "#FF6B6B", "#FFD700", "#1A9850", "#290333", "#909009",
                   "#FDAE61", "#F46D43", "#D73027", "#A50026")

# cols = c("#990090", "orange", "#909009", "#970900", "#109330", "blue","#290333", "#296330", "#995330")

color.groups <- function(var =  otu.sum$migration, 
                         cols = color_palette){
  
  groups <- unique(sort(var))
  cols.assigned <- NULL
  for(i in 1 : length(groups)){
    cols.assigned[var == groups[i]] <- cols[i]
  }
  return(cols.assigned)
}

color.key <- function(var = map$participant, 
                      key){
  # key <- aggregate(map.human$color, list(map.human$donor_ids), unique)

  groups <- unique(sort(var))
  cols.assigned <- NULL
  for(i in 1 : length(groups)){
    cols.assigned[var == key$Group.1[i]] <- key$x[i]
  }
  return(cols.assigned)
}


mod.get.groups <- function(mod){
  
  aa <- names(mod$strata)
  
  bb <- strsplit(aa, ",")
  
  groups <- NULL
  for(i in 1 : length(bb)){
    
    group.i <- NULL
    for(j in 1 : length(bb[[i]])){
      
      group.i[j] <- sapply(strsplit(bb[[i]][j], "="), "[[", 2)
      
    }
    groups[i] <- paste(group.i, collapse = " ")
  }
  
  groups <- gsub("\\s+", " ", stringr::str_trim(groups))
  
  return(groups)
}

plot.reslengths <- function(who = "p2", cutoffs = cutoffs.used){
  
  otu.sum.p <- otu.sum[otu.sum$participant == who & otu.sum$migration != "missing" , ]
  o <- order(otu.sum.p$max.residency)
  
  qqnorm(otu.sum.p$max.residency[o], 
         main = who,
         cex = 3, 
         pch = 21, 
         lwd = .3,
         col = 8, 
         ylab = "Days observed in human donor", 
         bg = otu.sum.p$col[o])
  abline(h = cutoffs, col = 2, lty = 2)
  
  leg <- aggregate(otu.sum.p$migration[o], list(otu.sum.p$col[o]), unique)
  
  tab <- table(otu.sum.p$migration)
  
  leg$x <- paste0(leg$x, " (n=", tab[match(leg$x, names(tab))], ")")
  
  legend("topleft", 
         fill = leg$Group.1, 
         legend = leg$x)
}

plot2d <- function(a){
  
  row.names(a) <- 1 : nrow(a)
  colnames(a) <- 1 : ncol(a)
  
  a <- a[rowSums(a, na.rm = T) > 0 , ]
  a <- a[ , colSums(a, na.rm = T) > 0]
  b <- t(a)
  b[b > 0.05] <- NA
  #b[b < 0.005] <- -1
  
  gplots::heatmap.2(b, 
                    #na.color = "grey70",
                    dendrogram='none', 
                    trace='none',
                    Rowv=FALSE, Colv=FALSE)
}

biome.proportion <- function(var = map.mice$migrant,
                             map.mice = map.mice,
                             col = "green", 
                             by.donor = F){
  
  if(by.donor){
    
    cols <- c(rgb(  34/255, 139/255,  34/255, alpha = 0.35),
              rgb( 255/255, 127/255,  80/255, alpha = 0.35),
              rgb( 255/255,  20/255, 147/255, alpha = 0.35))
    
    cols2 <- c(rgb(  34/255, 139/255,  34/255, alpha = 1),
               rgb( 255/255, 127/255,  80/255, alpha = 1),
               rgb( 255/255,  20/255, 147/255, alpha = 1))
    
    
    participants <- sort(unique(map.mice$participant))
    aa <- NULL
    aa[1 : length(participants)] <- T
    aa[1] <- F
    
    for(i in 1 : length(participants)){
      
      p <- !is.na(map.mice$days.between.hosts) & map.mice$participant == participants[i]
      
      vioplot::vioplot(var[p] ~ map.mice$days.between.host[p], 
                       xlab = "days outside host",
                       # ylab = "proportion of microbiome", 
                       at = log(sort(unique(map.mice$days.between.hosts[p]) + 1)), 
                       xlim = c(-1, 6), 
                       plotCentre = "none",
                       ylim = range(var), 
                       drawRect = F,
                       col = cols[i],
                       add = aa[i], 
      )
      
      # points(log((map.mice$days.between.hosts[p]) + 1), var[p])
      
      stripchart(var[p] ~ log((map.mice$days.between.hosts[p]) + 1), at = log(sort(unique(map.mice$days.between.hosts[p]) + 1)), 
                 method = "jitter", vertical = T, 
                 pch = 21, add = TRUE, col = cols2[i], 
                 bg = cols[i])
      
      
      anova(lm(var[p] ~ log(map.mice$days.between.host[p] + 1)))
      print(cor.test(var[p], log(map.mice$days.between.host[p] + 1), method = "spearman"))
      
      abline(lm(var[p] ~ log(map.mice$days.between.host[p] + 1)),  
             col = cols2[i], 
             lwd = 1.5,
             lty = 2)
    }
    
    legend("topleft", fill = cols2, legend = participants, )
    
  }else{
    p <- !is.na(map.mice$days.between.hosts) 
    
    vioplot::vioplot(var[p] ~ map.mice$days.between.host[p], col = col, 
                     xlab = "days outside host", ylab = "proportion of microbiome", 
                     at = log(sort(unique(map.mice$days.between.hosts[p]) + 1)), type = "n")
    
  }
  
  p <- !is.na(map.mice$days.between.hosts) 
  
  mod <- lmer(var[p] ~ log(map.mice$days.between.host[p] + 1) + map.mice$participant[p] + (1|map.mice$mouse.id[p]))
  mod.null <- lmer(var[p] ~  map.mice$participant[p] + (1|map.mice$mouse.id[p]))
  res <- anova(mod, mod.null)
  
  print(summary(mod))
  
  return(res)
}

biome.proportion.cat <- function(var = map.mice$migrant, map.mice = map.mice, col = "green"){
  
  p <- !is.na(map.mice$days.between.hosts) 
  #plot(var[p] ~ map.mice$days.between.host[p], pch = 21, bg = as.numeric(as.factor(map.mice$participant)[p]) + 1)
  
  vioplot::vioplot(var[p] ~ map.mice$days.between.host[p], col = col, 
                   xlab = "days outside host", ylab = "proportion of microbiome")
  
  print(anova(lm(var[p] ~ map.mice$days.between.host[p])))
  
  mod <- lmer(var[p] ~ as.character(map.mice$days.between.host[p]) + map.mice$participant[p] + (1|map.mice$mouse.id[p]))
  mod.null <- lmer(var[p] ~ map.mice$participant[p] + (1|map.mice$mouse.id[p]))
  res <- anova(mod, mod.null)
  
  print(summary(mod))
  
  return(res)
}

taxplot.migration <- function(group = "Resident_2 both"){
  
  par(mar = c(8, 4, 4, 2))
  
  otu.dummy <- otu  
  ps <- sort(unique(map$participant))
  
  otu.sum
  
  for(i in 1 : length(ps)){
    
    otu.sum.i <- otu.sum[otu.sum$participant == ps[i] , ]
    
    pick.1 <- which(map$participant == ps[i])
    pick.2 <- otu.sum.i$otu[otu.sum.i$migration.censor != group]
    pick.2 <- which(colnames(otu.dummy) %in% pick.2)
    
    otu.dummy[pick.1 , pick.2] <- 0
  }
  
  otu.dummy[,3]
  otu[,3]
  
  otu.sum$migration.censor == "Resident_2 both"
  
  all(row.names(otu.dummy) == map$sample.id)
  
  fam <- aggregate(t(otu.dummy), list(tax$family), sum)
  fam$Group.1
  row.names(fam) <- fam$Group.1
  fam$Group.1 <- NULL
  
  all(colnames(fam) == map$sample.id)
  
  tab <- table(tax$family)
  
  row.names(fam)  <- paste0(row.names(fam),  " (n=", tab,")")
  fam <- t(fam)
  
  bar <- aggregate(fam, 
                   list(map$time.between.hosts, 
                        map$participant), 
                   mean)
  
  bar.num <- bar[,-c(1:2)]
  
  #bar.num <- bar.num + 1
  #bar.num <- log(bar.num)
  
  bar.num <- bar.num[ , order(colSums(bar.num), decreasing = T)]
  bar.pick <- bar.num[ , colSums(bar.num) > sum(bar.num) / 100]
  bar.pick$other <-  rowSums(bar.num[,colSums(bar.num) <= sum(bar.num) / 100])
  
  row.names(bar.pick) <- paste(bar$Group.1, bar$Group.2)
  
  mat <- t(as.matrix(bar.pick))
  
  color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  set.seed(140404)
  cols <- sample(color, ncol(bar))
  col.key <- cbind(colnames(bar), cols, sapply(strsplit(colnames(bar), " "), "[[", 1))
  col.key <- rbind(col.key, c("other", "grey", "other"))
  
  col.key.pick <- col.key[match(row.names(mat), col.key[,1]) , ]
  
  group.sums <- aggregate(rowSums(otu), 
                          list(map$time.between.hosts, 
                               map$participant),
                          mean)
  
  non.migs <- group.sums$x - colSums(mat)  
  barplot(rbind(mat, non.migs), width = 1, beside = F, 
          space = c(1, .2, .2, 0, 0, 0, 0),
          xlim = c(0, ncol(bar.num) * 1.5),
          col =  c(col.key.pick[,2], "white"), las = 2, 
          main = paste("Taxonomy of", group, "bin\nbased on human longitudinal study"), 
          ylab = "mean count")
  
  legend("right", 
         bty = "n",
         fill = rev(col.key.pick[,2]), 
         legend = rev(col.key.pick[,1]))
  
  legend("topright", 
         bty = "n",
         fill = "white", 
         legend = paste("Non", group))
}

lmer4 <- function(response = 1:10, 
                  predictor =1:10,
                  mice = c(rep("A", 5), 
                           rep("B", 5)),
                  donors = c(rep("A", 5), 
                             rep("B", 5)), 
                  plot = T,
                  col.points = c("forestgreen", "coral", "deeppink"), 
                  alpha = .35, 
                  log = T){
  
  response.name <- deparse(substitute(response))
  response.name <- gsub(".*[$]","", response.name)
  
  predictor.name <- deparse(substitute(predictor))
  predictor.name <- gsub(".*[$]","", predictor.name)
  
  r_colors <- data.frame(color = colors())
  r_colors <- cbind(r_colors, t(col2rgb(colors())))
  r_colors.pick <- r_colors[ match(col.points, r_colors$color) , ]
  
  r_colors.pick <- rgb(r_colors.pick[2:4]/255, alpha = alpha)
  
  p <- !is.na(response) & !is.na(predictor) & !is.na(donors)
  
  response <- response[p]
  mice <- mice[p]
  donors <- donors[p]
  
  predictor <- predictor[p]
  try(predictor <- droplevels(predictor), silent = T)
  

  donors.unique <- sort(unique(donors))
  
  add <- NULL
  add[1 : length(donors.unique)] <- T
  add[1] <- F
  
  mod <- lme4::lmer(response ~ predictor + (1|mice) + (1|donors))
  mod.null <- lme4::lmer(response ~ (1|mice) + (1|donors))
  anova.res <- anova(mod, mod.null)
  mod.sum <- summary(mod)
  
  tail(anova.res$`Pr(>Chisq)`, 1)
  
  res <- paste("Chisq =", round(tail(anova.res$Chisq, 1), 3), 
               "p =", round(tail(anova.res$`Pr(>Chisq)`, 1), 5))
  
  if(plot){
    
    for(i in 1 : length(donors.unique)){
      
      p.i <- donors == donors.unique[i]
      
      if(log){
        
        names <- sort(unique(predictor))
        predictor.i <- log(predictor + 1)
        here <- sort(unique(predictor.i))
        
        predictor.name <- paste(predictor.name, "\n log transformed")
        
      }else{
        
        predictor.i <- predictor
        names <- sort(unique(predictor))# CHECK THIS!!
        here <- sort(unique(predictor.i))
        
      }
      
      vioplot::vioplot(response[p.i] ~ predictor.i[p.i], 
                       at = as.numeric(here),
                       names=names,
                       xlab = predictor.name,
                       ylab = response.name,
                       ylim = range(response, na.rm = T), 
                       plotCentre = "none",
                       drawRect = F, 
                       col = r_colors.pick[i],
                       add = add[i], 
                       main = res)
      
      stripchart(response[p.i] ~ predictor[p.i], 
                 at = here,
                 method = "jitter",
                 vertical = T, 
                 pch = 21, 
                 add = TRUE,
                 col = 1, 
                 bg = col.points[i],)
      
      abline(lm(response[p.i] ~ predictor.i[p.i]),  
             lwd = 1.5,
             lty = 2, 
             col = col.points[i])
    
    }
  }
  
  print(mod.sum)
  return(anova.res)
}


bc.pairs <- function(otu, 
                     variable = map$Description2, 
                     group1 = "p1 donor d0",
                     group2 = "p1 donor d1"){
  # variable is categorical and must be same length as nrow otu
  
  if(nrow(otu) != length(variable)){
    
    stop("ERROR #1")
    
  }else{
    
    variable.pick <- variable[variable %in% c(group1, group2)]
    otu.pick <- otu[variable %in% c(group1, group2) , ]

    otu.pick <- otu.pick[ , colSums(otu.pick) > 0 ]
    
    otu.dist <- vegdist(otu.pick)

    res <- NULL
    
    adonis.res <- adonis2(otu.dist ~ variable.pick)
    res$p.value <- adonis.res$`Pr(>F)`[1]
    res$f.stat <- adonis.res$F[1]
    res$r2 <- adonis.res$R2[1]
    
    centroids <- aggregate(otu.pick, list(variable.pick), mean)
    centroids$Group.1 <- NULL
   
    res$bc <- vegdist(centroids)

    res <- as.data.frame(res)
    
    return(res)
  }
}


bc.pairs.nested <- function(otu, 
                            variable = map$Description, 
                            variable.constant = map$participant, 
                            group1 = "donor d0",
                            group2 = "donor d1"){
  # variable is categorical and must be same length as nrow otu
  
  if(nrow(otu) != length(variable)){
    
    stop("ERROR #1")
    
  }else{
    
    p <- variable %in% c(group1, group2)
    
    variable.pick <- variable[p]
    constant.pick <- variable.constant[p]
    otu.pick <- otu[ p , ]
    
    otu.pick <- otu.pick[ , colSums(otu.pick) > 0 ]
    
    otu.dist <- vegdist(otu.pick)
    
    res <- NULL
    
    adonis.res <- adonis2(otu.dist ~ variable.pick + constant.pick)
    
    print(adonis.res)
    
    res$p.value <- adonis.res$`Pr(>F)`[1]
    res$f.stat <- adonis.res$F[1]
    res$r2 <- adonis.res$R2[1]
    
    centroids <- aggregate(otu.pick, list(variable.pick), mean)
    centroids$Group.1 <- NULL
    
    res$bc <- vegdist(centroids)
    
    res <- as.data.frame(res)
    
    return(res)
  }
}

otu.summarize.donors <- function(otu, 
                                 map, 
                                 cutoffs = 100,
                                 by.site = T, 
                                 alsways.there = T, 
                                 use.median = T){
  # residency calculate based on OTUs present in the samples used to humanize donor mice. 
  #   if an OTU recolonized a human donor multiple times. use the residency time containing the date used for humanization
  
  # by.site allows you to specify different cutoff by site/donor/individual - lust use list
  
  # alsways.there as in alwas a resident in the human to define resident_2 
  
  # use.median if median is in GF recipiant than the event occured, else max is used (only one opservation needed to say event has NOT occured)
  
  all(row.names(otu) == map$sample.id)
  
  otu.human <- otu[map$group == "human" , ]
  map.human <- map[map$group == "human" , ]
  otu.human <- otu.human[order(map.human$Date.Collected) , ]
  map.human <- map.human[order(map.human$Date.Collected) , ]
  
  participants <- sort(unique(map$participant))
  otus <- colnames(otu)
  otu.sum <- NULL
  
  for(i in 1 : length(participants)){
    
    print(paste("Calculating summary for participant", i))
    
    otu.i <- otu.human[map.human$participant == participants[i] , ]
    map.i <- map.human[map.human$participant == participants[i] , ]

    cutoffs.using <- NULL
    
    if(by.site & alsways.there){
      cutoffs.using <- diff(range(map.i$Date.Collected))
    }
    
    else if(by.site & !alsways.there){
      cutoffs.using <- cutoffs[[i]]
    }
    
    else if(!by.site & !alsways.there){
      cutoffs.using <- cutoffs
      
    }else{
      
      stop("alsways.there only availbe if by site = T")
    }
    
    
    for(j in 1 : ncol(otu.i)){
      
      ###################### OTUs in the HUMANs ########################################
      otu.sum.i <- classify.residency(otu.pick = otu.i[,j], 
                                      dates = map.i$Date.Collected, 
                                      cutoffs = cutoffs.using)
      
      otu.sum.i$participant <- participants[i]
      otu.sum.i$otu <- colnames(otu.i)[j]
      otu.sum.i$otu.id <- paste0(otu.sum.i$participant,".", otu.sum.i$otu)

      otu.sum <- rbind(otu.sum, otu.sum.i)
    }
  }
  ################################################################################
  
  otu.sum$human <- "missing from humanization sample"
  otu.sum$human[otu.sum$count.humanize > 0] <- "in humanization sample"

  otu.sum$censor2 <- otu.sum$censor
  otu.sum$censor2[otu.sum$censor == "left" | otu.sum$censor == "right"] <- "one side"
  otu.sum$censor2 <- factor(otu.sum$censor2, levels = c("both", "one side", "no"))
  
  otu.sum$migration.censor <- paste(otu.sum$migration, otu.sum$censor2)# don't use "+" because not all levels are meaningful
  
  return(otu.sum)
}


blastn <- function(query.fasta, genome.fasta){
  
  makedb <- paste("makeblastdb", "-dbtype nucl", "-in", genome.fasta)
  system(makedb)
  
  blast <- paste("blastn", "-query", query.fasta, 
                 "-db", genome.fasta, 
                 "-out", "temp.txt", "-outfmt 6")
  system(blast)
  
  if(file.size("temp.txt") > 0){
    blast.res <- read.table("temp.txt")
    colnames(blast.res) <- c("qseqid", 
                             "sseqid", 
                             "pident", 
                             "length", 
                             "mismatch", 
                             "gapopen", 
                             "qstart", 
                             "qend", 
                             "sstart", 
                             "send", 
                             "evalue", 
                             "bitscore")
    return(blast.res)
  }else{
    return(NULL)
  }
}

complement <- function(dna){
  
  dna.inter <- toupper(dna)
  
  dna.inter <- gsub("T", "a", dna.inter)
  dna.inter <- gsub("G", "c", dna.inter)
  dna.inter <- gsub("A", "t", dna.inter)
  dna.inter <- gsub("C", "g", dna.inter)
  
  res <- toupper(dna.inter)
  
  return(res)
}

filter.outliers <- function(var){
  
  Q <- quantile(var, probs=c(.25, .75), na.rm = FALSE)
  
  iqr <- IQR(var)
  
  up <-  Q[2]+1.5*iqr # Upper Range  
  low<- Q[1]-1.5*iqr # Lower Range﻿
  
  return(var >= low & var <= up)
}

plot.1cat <- function (var = Invsimp, groups = map$antibiotic.group[p], col = map$col.group[p], 
                       exclude.outliers = T, p.adjust = T, symbols = T, plot = T, ranked.sum = F, ylim){
  if (missing(col)) {
    col <- color.groups(groups)
  }
  if (missing(ylim)) {
    ylim <- c(min(var, na.rm = T), max(var, na.rm = T) * 
                1.1)
  }
  var.name <- deparse(substitute(var))
  groups.name <- deparse(substitute(groups))
  if (typeof(groups) == "character") {
    groups <- factor(groups)
  }
  try({
    col.agg <- aggregate(col, list(groups), unique)
  }, silent = T)
  if (max(lengths(col.agg$x)) > 1) {
    col.agg$x <- "grey"
    message("colors did not correspond to groups")
  }
  var1 <- var[!is.na(groups) & !is.na(var)]
  groups1 <- groups[!is.na(groups) & !is.na(var)]
  var <- var1
  groups <- droplevels(groups1)
  n <- table(groups)
  n <- paste0(names(n), " (n=", n, ")")
  
  groups.unique <- levels(groups)
  groups.unique <- groups.unique[!is.na(groups.unique)]
  
  if(plot){
    vioplot::vioplot(var ~ groups, 
                     col = col.agg$x,
                     ylab = "", 
                     xlab = "",
                     drawRect = F,
                     las = 2, 
                     font.axis = 4, 
                     font.lab=4,
                     names =  n,
                     ylim = ylim)
    title(ylab = var.name, font.lab = 4)
    
    for(i in 1 : length(groups.unique)){
      
      p.i <- which(groups == groups.unique[i])
      n.i <- length(p.i)
      
      points(i + sample(-100 : 100,  n.i, replace = T)/2000, 
             var[p.i], pch = 21, bg = 8)
    }
  }
  outs <- aggregate(var, list(groups), filter.outliers)
  
  aa = 0
  groups.unique <- levels(groups)
  groups.unique <- groups.unique[!is.na(groups.unique)]
  n.comps <- sum(1:(length(groups.unique) - 1))
  t.res.all <- as.data.frame(matrix(nrow = n.comps, ncol = 8))
  colnames(t.res.all) <- c("group1", "group2", "p", "effect_size", 
                           "t_stat", "df", "95_lo", "95_hi")
  counts = 0
  for (i in 1:length(groups.unique)) {
    for (j in 1:i) {
      if (j != i) {
        counts <- counts + 1
        var.i <- var[which(groups == groups.unique[i])]
        var.j <- var[which(groups == groups.unique[j])]
        if (exclude.outliers) {
          t.res <- t.test(var.i[outs$x[[i]]], var.j[outs$x[[j]]])
        }
        if (exclude.outliers) {
          t.res <- t.test(var.i[outs$x[[i]]], var.j[outs$x[[j]]])
        }
        
        if (ranked.sum) {
          # dummy
          t.res <- wilcox.test(var.i, var.j)
          t.res$estimate <- c(NA, NA)
          t.res$statistic <- NA
          t.res$parameter <- NA
          t.res$conf.int[1:2] <- NA
        }
        
        else {
          t.res <- t.test(var.i, var.j)
        }
        t.res.all[counts, ] <- c(as.character(groups.unique[i]), 
                                 as.character(groups.unique[j]), t.res$p.value, 
                                 diff(t.res$estimate), t.res$statistic, t.res$parameter, 
                                 t.res$conf.int[1:2])
      }
    }
  }
  if (p.adjust) {
    t.res.all$p.fdr <- p.adjust(t.res.all$p)
  }
  else {
    t.res.all$p.fdr <- NA
  }
  t.res.all$symbol <- NA
  t.res.all$symbol[t.res.all$p.fdr < 0.05] <- "*"
  t.res.all$symbol[t.res.all$p.fdr < 0.01] <- "**"
  t.res.all$symbol[t.res.all$p.fdr < 0.001] <- "***"
  for (i in 1:nrow(t.res.all)) {
    if (t.res.all$p.fdr[i] < 0.05) {
      g1 <- which(groups.unique == t.res.all$group1[i])
      g2 <- which(groups.unique == t.res.all$group2[i])
      h <- max(var) + i * diff(range(var, na.rm = T))/50
      segments(x0 = g1, y0 = h, x1 = g2, y1 = h)
      if (symbols) {
        text(x = mean(c(g1, g2)), y = h, t.res.all$symbol[i])
      }
      else {
        text(x = mean(c(g1, g2)), y = h, paste("p = ", 
                                               round(t.res.all$p.fdr[i], 5), "\n"), col = 1, 
             cex = 0.5)
      }
    }
  }
  if (!p.adjust) {
    t.res.all$p.fdr <- NULL
  }
  return(t.res.all)
}

num_gene.surv <- function(surv.map, col){
  
  if (missing(col)) {
    col <- "red"
  }
  
  surv <- Surv(surv.map$mouse.residency, surv.map$mouse.residency < 90)
  
  efs <- ps <- NULL
  
  for(i in 1 : 10){
    
    p <- !is.na(surv.map$gene.mean.50)
    
    if(length(table(surv.map$gene.mean.50[p] > i)) > 1){
      
      res.i <- t.test(surv.map$residency[p] ~ surv.map$gene.mean.50[p] > i)
      efs[i] <- diff(res.i$estimate)
      ps[i] <- res.i$p.value
    }
  }
  
  efs <- as.numeric(efs)
  ps <- as.numeric(ps)
  
  ps.col <-  NULL
  ps.col[ps <= 0.05] <- col 
  ps.col[ps > 0.05] <- "grey"
  
  plot(efs, 
       pch = 21, 
       bg = ps.col, 
       cex = -log10(ps),
       xlab = "Number of Spore Genes in Annotation",
       ylab = "Increase in median human\nresidancy time (days)")
  
  coef <- ps <- NULL
  
  for(i in 1 : 10){
    
    sur.fit <- survival::survfit(surv ~ surv.map$gene.mean.50 < i)
    cox.res <- survival::coxph(surv ~ surv.map$gene.mean.50 < i)
    cox.res <- summary(cox.res)
    coef[i] <- cox.res$coefficients[2]
    ps[i]   <- cox.res$logtest[3]
  }
  
  ps.col <-  NULL
  ps.col[ps <= 0.05] <- col
  ps.col[ps > 0.05] <- "grey"
  
  plot(coef, 
       pch = 21,
       bg = ps.col, 
       cex = -log10(p.adjust(ps)), 
       #  ylim = c(min(coef, na.rm = T) * .9, max(coef, na.rm = T) * 1.1),
       ylab = "Cox coefficient",
       xlab = "Number of spore genes in annotations")
  #legend("topleft", fill = c("red", "grey"), legend = c("p <= 0.05", "p > 0.05"), bty = "n")
  legend("bottomleft", fill = c(col, "grey"), legend = c("p <= 0.05", "p > 0.05"), bty = "n")
}


tax.shared.multicomp <- function (otu, 
                                  var1, 
                                  var2, 
                                  var1.groups,
                                  var2.groups, 
                                  tax.level = tax$family, 
                                  group1.color = "orange", 
                                  group2.color = "purple",
                                  PresAbs = F, 
                                  log = T,
                                  taxa2include) {
  
  # var1.groups must have 2 groups 
  # var1.groups can have more any number of groups
  if (length(group1.color) == 1) {
    group1.color <- rep(group1.color, nrow(otu))
  }
  if (length(group2.color) == 1) {
    group2.color <- rep(group2.color, nrow(otu))
  }
  
  n <- length(var2.groups)
  n.taxa <- length(taxa2include)
  add <- vector(length = n)
  add[-1] <- T
  
  res.all <- list()
  n <- length(var2.groups)
  
  for(i in 1 : n){
    
    var1.p <- var1[var2 == var2.groups[i]]
    otu.p <- otu[var2 == var2.groups[i],]
    
    res.i <- tax.shared.table(otu.p, 
                              var1.p,
                              group1 = var1.groups[1], 
                              group2 = var1.groups[2], 
                              tax.level, 
                              PresAbs,
                              log,
                              taxa2include)
    
    
    res.all[[i]] <- res.i
  }
  
  names(res.all) <- var2.groups
  here = max(sapply(res.all, function(x) colSums(x)))
  unique.1.col <- NULL
  shared.1.col <- NULL
  shared.2.col <- NULL
  unique.2.col <- NULL
  
  par(mar = c(3, 20, 4, 20))
  par(xpd = TRUE)
  
  for (i in 1:n) {
    p <- var2 == var2.groups[i]
    unique.1.col[i] <- group1.color[p][1]
    unique.2.col[i] <- group2.color[p][1]
    if (i != 1) {
      new.colnames <- gsub(paste0(taxa2include, collapse = "|"), 
                           "", colnames(res.all[[i]]))
      new.colnames <- gsub("  ", "", new.colnames)
    }
    else {
      new.colnames <- colnames(res.all[[i]])
    }
    shared.1 <- colorspace::mixcolor(0.2, colorspace::sRGB(t(col2rgb(group1.color[p][1]))), 
                                     colorspace::sRGB(t(col2rgb(group2.color[p][1]))))
    shared.2 <- colorspace::mixcolor(0.5, colorspace::sRGB(t(col2rgb(group1.color[p][1]))), 
                                     colorspace::sRGB(t(col2rgb(group2.color[p][1]))))
    shared.1.col[i] <- rgb(shared.1@coords, maxColorValue = 255)
    shared.2.col[i] <- rgb(shared.2@coords, maxColorValue = 255)
    colnames(res.all[[i]]) <- new.colnames
    barplot(res.all[[i]], ylim = c(0, ncol(res.all[[i]]) * 
                                     (n + 1)), xlim = c(0, here), cex.axis = 0.5, cex.names = 0.5, 
            xaxt = "n", horiz = T, main = "", xlab = "median population size (log 10)", 
            las = 1, add = add[i], space = c(n - i, rep(n, ncol(res.all[[i]]))[-1]), 
            col = c(unique.1.col[i], shared.1.col[i], shared.2.col[i], 
                    unique.2.col[i]))
  }
  if (log) {
    segments(x0 = here - log10(11), x1 = here, y0 = 14, y1 = 14, 
             lwd = 2)
    text(here - 0.5, 15, "10", cex = 0.5)
    segments(x0 = here - log10(101), x1 = here, y0 = 12, 
             y1 = 12, lwd = 2)
    text(here - 0.5, 13, "100", cex = 0.5)
    segments(x0 = here - 3, x1 = here, y0 = 10, y1 = 10, 
             lwd = 2)
    text(here - 0.5, 11, "1000", cex = 0.5)
    segments(x0 = here - 3, x1 = here, y0 = 10, y1 = 10, 
             lwd = 2)
    text(here - 0.5, 11, "1000", cex = 0.5)
  }
  text(here + 0.3, n.taxa * n/2, srt = 0, "u1", pos = 4)
  text(here + 0.6, n.taxa * n/2, srt = 0, "s1", pos = 4)
  text(here + 0.9, n.taxa * n/2, srt = 0, "s2", pos = 4)
  text(here + 1.2, n.taxa * n/2, srt = 0, "u2", pos = 4)
  legend(here + 0.3, n.taxa * n/2, fill = unique.1.col, legend = rep("", 
                                                                     n), cex = 1, bty = "n")
  legend(here + 0.6, n.taxa * n/2, fill = shared.1.col, legend = rep("", 
                                                                     n), cex = 1, bty = "n")
  legend(here + 0.9, n.taxa * n/2, fill = shared.2.col, legend = rep("", 
                                                                     n), cex = 1, bty = "n")
  legend(here + 1.2, n.taxa * n/2, fill = unique.2.col, legend = var2.groups, 
         cex = 1, bty = "n")
}
                    
read.mothur.taxonomy <- function (cons.taxonomy) {
  tax <- read.table(cons.taxonomy, header = T, row.names = 1)
  tax.list <- strsplit(tax$Taxonomy, ";|[(]|[)]")
  tax.df <- matrix(nrow = nrow(tax), ncol = 18)
  for (i in 1:18) {
    tax.df[, i] <- sapply(tax.list, "[[", i)
  }
  tax.df <- as.data.frame(tax.df)
  tax.df <- tax.df[, tax.df[1, ] != ""]
  colnames(tax.df) <- c("Kingdom", "King.conf", "Phylum", "Phy.conf", 
                        "Class", "Class.conf", "Order", "Order.conf", "family", 
                        "fam.conf", "genus", "gen.conf")
  return(tax.df)
}
                    


