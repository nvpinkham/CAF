#' Blast wrapper
#'
#' Run local blast and import the results as text file
#' @param query.fasta fasta file of sequence to search for
#' @param genome.fasta fasta file where the query.fasta sequences are searched for
#' @return hit table
#' @examples
#' @export
blast.FUNction <- function(query.fasta = "WADDL_DiagnosticPCR_Primers_dummy.fas",
                           genome.fasta = "test_temp.fasta",
                           word_size = 11, max_tar = 1){

  makedb <- paste("makeblastdb", "-dbtype nucl", "-in", genome.fasta)
  system(makedb)

  blast <- paste("blastn", "-query", query.fasta,
                 "-db", genome.fasta,
                 "-out", "temp.txt",
                 "-outfmt 6",
                 "-subject_besthit",
                 "-max_target_seqs", max_tar,
                 "-evalue 1000",
                 "-word_size", word_size)
  system(blast)

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
}

#' Remove outliers
#'
#' Remove outliers from vector
#' @param var numeric vector
#' @return numeric vector without outliers
#' @export
filter.outliers <- function(var){

  Q <- quantile(var, probs=c(.25, .75), na.rm = T)

  iqr <- IQR(var, na.rm = T)

  up <-  Q[2]+1.5*iqr # Upper Range
  low<- Q[1]-1.5*iqr # Lower Range﻿

  return(var > low & var < up)
}

#' Plots 1 numeric vector against 1 categorical vector
#'
#'
#' @param var numeric vector
#' @param groups categories vector
#' @return t test results
#' @export
plot.1var <- function(var = map$dist2a1,
                      groups = map$discription,
                      col = map$col,
                      exclude.outliers = TRUE){
  
  if(missing(col)) {
    col <- color.groups(groups)
  }
  
  var.name <- deparse(substitute(var))
  groups.name <- deparse(substitute(groups))
  
  if(typeof(groups) == "character"){
    groups <- factor(groups)
  }

  try({col.agg <- aggregate(col, list(groups), unique)}, silent = T)
          
  if(max(lengths(col.agg$x)) > 1){
    col.agg$x <- "grey"
    message("colors did not correspond to groups")
  }
  
  var1 <- var[!is.na(groups) & !is.na(var)]
  groups1 <- groups[!is.na(groups) & !is.na(var)]
  
  var <- var1
  groups <- droplevels(groups1)
  
  boxplot(var ~ groups,
          col = col.agg$x,
          xlab = groups.name,
          ylab = var.name,
          ylim = range(var, na.rm = T) * c(.9, 1.1))
  
  outs <- aggregate(var, list(groups), filter.outliers)
  aa = 0
  
  groups.unique <- levels(groups)
  groups.unique <- groups.unique[!is.na(groups.unique)]
  
  n.comps <- sum(1 : (length(groups.unique) - 1))
  t.res.all <- as.data.frame(matrix(nrow = n.comps, ncol = 8))
  
  colnames(t.res.all) <- c("group1", "group2",
                           "p", "effect_size", "t_stat", "df", "95_lo", "95_hi")
  counts = 0
  
  for(i in 1 : length(groups.unique)){
    for(j in 1 : i){
      if(j != i){
        
        counts <- counts + 1
        
        var.i <- var[which(groups == groups.unique[i])]
        var.j <- var[which(groups == groups.unique[j])]
        
        if(exclude.outliers){
          
          t.res <- t.test(var.i[outs$x[[i]]], var.j[outs$x[[j]]])
        }else{
          t.res <- t.test(var.i, var.j)
        }
        
        t.res.all[counts,] <- c(as.character(groups.unique[i]),
                                as.character(groups.unique[j]),
                                t.res$p.value,
                                diff(t.res$estimate),
                                t.res$statistic,
                                t.res$parameter,
                                t.res$conf.int[1:2])
        
        
        if(t.res$p.value < 0.05){
          aa <- aa + (diff(range(var, na.rm = T)) / 50)
          
          print(i)
          print(j)
          
          #boxplot(var.i, var.j, main = t.res$p.value)
          
          h <- max(c(var.i, var.j), na.rm = T) + (aa * 2)
          
          
          segments(x0 = i, y0 = h, x1 = j, y1 = h)
          text(mean(c(i,j)), h,
               paste("p = ", round(t.res$p.value, 5), "\n"),
               col = 1, cex = .5)
        }
      }
    }
  }
  return(t.res.all)
}

#' Plots NMDS from OTU table
#'
#'
#' @param otu site abundance table
#' @param groups categories vector
#' @return t test results
#' @export
group.nmds <- function(otu = otu.nasal,
                       predictor = map.nasal$`MOVI PCR`,
                       groups = c("positive", "negative"),
                       size = 1,
                       color = map.nasal$col1,
                       return.table = F,
                       spider = T,
                       legend.location = "bottomleft"){
  # predictor is a vector of the catagorical predicting variable
  # groups are the paticular catagories you want to explore
  # col is a vertor coresponding to the groups in the predictor

  all.groups <-  unique(predictor)
  others <- all.groups [!all.groups %in% groups]
  others <- paste(  others , collapse = "\nand ")

  otu.dist <- vegan::vegdist(otu)

  nmds <- metaMDS(otu.dist)
  plot(nmds$points, pch = 19, col = 8)#, cex  = map.fecal$invsimp / 10)

  nmds.pick <- nmds$points[ predictor %in% groups , ]
  col.pick <- color[ predictor %in% groups]
  otu.pick <- otu[ predictor %in% groups , ]

  predictor.pick <- predictor[ which(predictor %in% groups)]

  try(
    {
      col.key <- aggregate(col.pick, list(predictor.pick), unique)
      ns <- table(predictor.pick)
      m <- match(col.key$Group.1, names(ns))
      col.key$ns <- ns[m]
      legend(legend.location, fill = col.key$x,
             legend = paste0(col.key$Group.1, " (n=", col.key$ns, ")"),
             bty = "n")

      if(spider){
        ordispider(nmds.pick, predictor.pick,
                   label = F, lwd = 1, col = "darkgrey")
      }else{
        ordiellipse(nmds.pick, predictor.pick,
                    label = F, lwd = 1, col = col.key$x)
      }
    }
  )

  points(nmds.pick, pch = 21, bg = col.pick, cex = size)

  p <- !is.na(predictor.pick)
  res <- vegan::adonis2(vegdist(otu.pick[  p , ]) ~ predictor.pick[p], permutations = 9999)

  bres <- betadisper(vegdist(otu.pick[p , ]), predictor.pick[p])
  bres$group.distances

  predictor.unique <- unique(predictor.pick)

  otu.agg <- aggregate(otu.pick, by = list(predictor.pick), FUN = mean)
  row.names(otu.agg) <- otu.agg$Group.1
  otu.agg$Group.1 <-NULL
  centroid.dist <- vegdist(otu.agg)

  print(centroid.dist)

  res1 <- c( round(res$`Pr(>F)`[1], 4),
             round(res$R2[1], 3),
             round(res$F[1], 3),
             round(max(centroid.dist), 3),
             nrow(otu.pick))
  names(res1) <- c("p.val", "R2", "f.stat", "distance.between.centroids", "n")

  title(paste0("PERMANOVA\np val =", res1[1],
               ", R2 =", res1[2],
               ", f stat =", res1[3],
               ", n=", res1[5],
               "\n max BC between centroids =", res1[4]))

  if(return.table){

    return(res1)

  }

  return(nmds)
}


group.t <- function(response = map$invsimp,
                    predictor = map.nasal$`MOVI PCR`,
                    groups = c("positive", "negative"),
                    color = map.nasal$col1, var.id = "invsimp"){
  # predictor is a vector of the catagorical predicting variable
  # groups are the paticular catagories you want to explore
  # col is a vertor coresponding to the groups in the predictor

  all.groups <-  unique(predictor)
  others <- all.groups [!all.groups %in% groups]
  others <- paste(  others , collapse = "\nand ")

  col.pick <- color[ which(predictor %in% groups) ]
  response.pick <- response[  which(predictor %in% groups) ]
  predictor.pick <- predictor[  which(predictor %in% groups) ]

  boxplot(response.pick ~ predictor.pick, col = unique(col.pick), xlab = "", ylab = var.id)

  res <- t.test(response.pick ~ predictor.pick)

  res <- paste("t test\np value =", round(res$p.value, 3),
               "\n diff between means", round(diff(res$estimate), 3))

  title(res)
}


#' makes figure illustrating taxa that are shared between two sites and unique to either site
#'
#'
#' @param otu site abundance table
#' @param map environmental data coresponding to each site in the the site abundance table
#' @return
#' @export
tax.shared <- function(otu, map,  group1, group2,
                       tax.level = tax$family, binary = F, use.max = F, log = T){

  tax.familes <- sort(unique(tax.level))
  fam.count <- aggregate(t(otu), list(tax.level), sum)
  row.names(fam.count) <- fam.count$Group.1
  fam.count$Group.1 <- NULL
  fam.count <- as.data.frame(t(fam.count))

  fam.count <- fam.count[,colSums(fam.count) > 0]
  #### Run the thing #####

  map$discription <- as.factor(map$discription)

  p <- otu[map$discription == group1, ]
  r <- otu[map$discription == group2, ]

  if(binary){
    p[p > 0] <- 1
    r[r > 0] <- 1
  }

  p.dummy2 <- p.dummy <- p
  p.dummy[,colSums(r) > 0] <- 0# shared OTUs removed
  p.dummy2[,colSums(  p.dummy) > 0] <- 0# unique OTUs removed

  r.dummy2 <-  r.dummy <- r
  r.dummy[,colSums(p) > 0] <- 0# shared OTUs removed
  r.dummy2[,colSums(  r.dummy) > 0] <- 0# unique OTUs removed

  if(use.max){

    p.unique <- aggregate(t(p.dummy), list(tax.level), max)
    r.unique <- aggregate(t(r.dummy), list(tax.level), max)

    p.shared <- aggregate(t(p.dummy2), list(tax.level), max)
    r.shared <- aggregate(t(r.dummy2), list(tax.level), max)

  }else{

    p.unique <- aggregate(t(p.dummy), list(tax.level), sum)
    p.shared <- aggregate(t(p.dummy2), list(tax.level), sum)

    r.unique <- aggregate(t(r.dummy), list(tax.level), sum)
    r.shared <- aggregate(t(r.dummy2), list(tax.level), sum)
  }

  p.u <-  rowMeans(p.unique[,-1])
  names(p.u) <- p.unique$Group.1

  r.u <-  rowMeans(r.unique[,-1])
  names(r.u) <- r.unique$Group.1

  p.s <-  rowMeans(p.shared[,-1])
  names(p.s) <- p.shared$Group.1

  r.s <-  rowMeans(r.shared[,-1])
  names(r.s) <- r.shared$Group.1

  tax.res <- rbind(  p.u,
                     p.s,
                     r.s,
                     r.u)

  rownames(tax.res) <- c(paste("unique to", group1),
                         paste("shared - pop size in", group1),
                         paste("shared - pop size in", group2),
                         paste("unique to ", group2))

  tax.res <- tax.res[ , colSums(tax.res) != 0 ]

  par(mar = c(3, 15, 4, 4))

  par(xpd = TRUE) #Draw outside plot area

  if(log){
    tax.res <- log10(as.matrix(tax.res) + 1)#+ .000001)
  }
  tax.res <- tax.res[ , order(colSums(tax.res), decreasing = F) ]

  here <- ceiling(max(colSums(tax.res))) - .5

  barplot(tax.res,
          # xlim = c(-3.3, 16),
          cex.axis=.75,
          cex.names=.75,
          xaxt='n',
          horiz=T ,
          main = "",
          xlab = "median population size (log 10)",
          las = 1,
          col =c( "orange","#b38a61", "#887382", "purple"))


  #segments(x0 = here - log10(2), x1 = here, y0 = 16, y1 = 16)
  #text(here - .1, 17, "5", cex = 0.5)

  segments(x0 = here - log10(11), x1 = here, y0 = 14, y1 = 14, lwd = 2)
  text(here -.5, 15, "10", cex = 0.5)

  segments(x0 = here - log10(101), x1 = here, y0 = 12, y1 = 12,  lwd = 2)
  text(here - .5, 13, "100", cex = 0.5)

  segments(x0 = here - 3, x1 = here, y0 = 10, y1 = 10, lwd = 2)
  text(here - .5, 11, "1000", cex = 0.5)


  legend("bottomright",
         fill =c( "orange","#b38a61", "#887382", "purple"),
         bty = "n",
         legend = row.names(tax.res))


  return(tax.res)
}

#' makes figure illustrating taxa that are shared between two sites and unique to either site
#'
#'
#' @param table site abundance table
#' @param cutoff_pro proportion of sites that an organism must be observed in ranges from 0 to 1
#' @return filtered site abundance table
#' @export
remove_rare <- function( table , cutoff_pro ) {

  table <- t(table)# transpose to keep "tidy" ; easier that rewriting function...
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

#map$discription <- map$`MOVI PCR`

#' makes figure illustrating taxa that are shared between two sites and unique to either site
#'
#'
#' @param mat matrix or OTU table
#' @param map met data, need column named "discription"
#' @return RF results
#' @export
pairwiseRF <- function(mat = otu.fecal,
                       map = map.fecal,
                       dis_1 = "positive",
                       dis_2 = "negative",
                       num_of_col = 50, Accuray.cutoff = 0,
                       Factor = 1, plot = T){
  # report mats that are different between two treatments based on RF analysis
  set.seed(42)
  mat.pick <- mat[map$discription %in% c(dis_1, dis_2 ), ]
  map.pick <- map[map$discription %in% c(dis_1, dis_2 ), ]

  #mat.pick$Group <- NULL

  mat.pick.rare_removed <- remove_rare(table = mat.pick, cutoff_pro = 0.2)

  spliter <- as.factor(as.character(map.pick$discription))

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

  if(RF.sig$pValue < 0.1){

    fit_imp <- as.data.frame( fit$importance )
    fit_imp$features <- rownames(fit_imp )
    fit_imp_sorted <- dplyr::arrange( fit_imp  , desc(fit_imp$MeanDecreaseGini)  )
    colors <- vector(length = ncol(mat.pick.rare_removed))

    for(j in 1:ncol(mat.pick.rare_removed)){
      i <- fit_imp_sorted$features[j]
      t1.mean <- mean(mat.pick[map.pick$discription == dis_1, which(colnames(mat.pick) == i)])
      t2.mean <- mean(mat.pick[map.pick$discription == dis_2, which(colnames(mat.pick) == i)])
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
    return("VARIABLE p val > 0.1")
  }
}



color.groups <- function(vec, cols = c("darkorange", #(Red)
                                       "lightgreen", #(Lime)
                                       "#800080", #(Purple)
                                       "#6A5ACD", #(Slate Blue)
                                       "#0080FF", #(Blue)
                                       "#FFFF00", #(Yellow)
                                       "#00FFFF", #(Cyan)
                                       "#FF00FF", #(Magenta)
                                       "#800000", #(Maroon)
                                       "#008000", #(Green)
                                       "#000080", #(Navy)
                                       "#808000", #(Olive)
                                       "#008080", #(Teal)
                                       "red", #(Orange)
                                       "#A52A2A", #(Brown)
                                       "#FFC0CB", #(Pink)
                                       "#FFD700", #(Gold)
                                       "#008B8B", #(Dark Cyan)
                                       "#2E8B57", #(Sea Green)
                                       "#CD5C5C")) #(Indian Red)))
{
  levels <- unique(vec)
  vec.col <- vector(length = length(vec))

  for(i in 1 : length(levels)){
    vec.col[vec == levels[i]] <- cols[i]
  }

  vec.col[is.na(vec)] <- "grey40"

  return(vec.col)
}

taxplot <- function(otu = otu.nasal,
                    tax = tax.nasal$family,
                    var = map.nasal$discription.MOV2,
                    cutoff = 0.05){

  par(mar = c(8, 4, 4, 2))

  otu <- otu / rowSums(otu)
  tax.agg <- aggregate(t(otu), list(tax), sum)
  rownames(tax.agg) <- tax.agg$Group.1
  tax.agg$Group.1 <- NULL

  tax.agg <- aggregate(t(tax.agg), list(var), mean)
  rownames(tax.agg) <- tax.agg$Group.1
  tax.agg$Group.1 <- NULL

  tax.agg <- t(tax.agg)

  tax.agg <- tax.agg[order(rowSums(tax.agg), decreasing = T) ,  ]

  tax.agg.pick <- tax.agg[rowSums(tax.agg) / sum(tax.agg) >= cutoff , ]
  tax.agg.other <- tax.agg[rowSums(tax.agg) / sum(tax.agg) < cutoff , ]

  other <- colSums(tax.agg.other)
  tax.agg.pick <- rbind(tax.agg.pick, other)

  color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  set.seed(140404)
  cols <- sample(color, nrow(tax.agg.pick))
  col.key <- cbind(rownames(tax.agg.pick), cols)
  # col.key <- rbind(col.key, c("other", "grey", "other"))
  col.key[nrow(col.key),2] <- "grey50"

  barplot(tax.agg.pick *100, width = 1, beside = F,
          xlim = c(0, ncol(tax.agg.pick) * 1.5),
          space = .1,
          col =  col.key[,2], las = 2,
          ylab = "mean % of community within group")

  legend("right",
         bty = "n",
         fill = rev(col.key[,2]),
         legend = rev(col.key[,1]))
}
