tab_chrom_mod<-function (phy, table, missing = NA, filename = "ChromEvol_analysis",
                         folder = NA, percentage = TRUE, mainType = "All_Models",
                         fix.root = FALSE, maxChrNum = "-10", minChrNum = "1", branchMul = "1",
                         rootAt = "N1", numberSim = NA, ladderize = TRUE)
{


  ##This function not authored by me. I made some changes to Natalie Cusimano's script.

  if (class(phy) != "phylo")
    phy <- read.tree(phy)
  if (!is.null(phy$posterior) && is.null(phy$node.label))
    phy$node.label <- round(phy$posterior, digits = 2)
  if (ladderize)
    phy <- ladderize(phy)
  phy <- fixNodes(phy)
  if (!is.matrix(table)) {
    tableS <- scan(table, what = "list")
    if (length(grep(">", tableS)) > 1) {
      table <- read.cE.tab(table)
      print(table)
    }
    else table <- read.csv("sub_table.csv", stringsAsFactors = FALSE)[,c(2:3)]
  }
  table2 <- match.table2tip(sub, table)
  table2[, 1] <- phy$tip.label
  if (length(which(is.na(table2[, 2]))) != 0) {
    if (is.na(missing)) {
      phy2 <- drop.tip(phy, which(is.na(table2[, 2])),
                       1)
      table3 <- table2[-which(is.na(table2[, 2])), ]
    }
    else {
      table2[which(is.na(table2[, 2])), 2] <- "X"
      table3 <- table2
      phy2 <- phy
    }
  }
  else {
    table3 <- table2
    phy2 <- phy
  }
  if (percentage) {
    multi <- grep("_", table3[, 2])
    if (length(multi) > 0) {
      perc <- grep("=", table3[, 2])
      if (length(perc) == 0) {
        multi2 <- strsplit(table3[multi, 2], "_")
        for (i in 1:length(multi2)) {
          p <- round(1/length(multi2[[i]]), digits = 2)
          if (length(multi2[[i]]) * p != 1)
            p1 <- (1 - length(multi2[[i]]) * p) + p
          else p1 <- p
          ll <- paste(multi2[[i]][1], "=", p1, sep = "")
          for (j in 2:(length(multi2[[i]]))) {
            ll <- paste(ll, "_", multi2[[i]][j], "=",
                        p, sep = "")
          }
          multi2[[i]] <- ll
          table3[multi[i], 2] <- multi2[[i]]
        }
      }
    }
  }
  if (is.na(numberSim)) {
    if (mainType == "All_Models")
      numberSim <- 0
    else numberSim <- 10000
  }
  if (!is.na(folder))
    setwd(folder)
  write.table(paste0(">", table3[1, 1], "\n", table3[1, 2]),
              file = paste(filename, "_input_table.txt", sep = ""),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  for (i in 2:length(table3[, 1])) write.table(paste0(">",
                                                      table3[i, 1], "\n", table3[i, 2]), file = paste(filename,
                                                                                                      "_input_table.txt", sep = ""), quote = FALSE, row.names = FALSE,
                                               col.names = FALSE, append = TRUE)
  fnp <- paste(filename, "_params.txt", sep = "")
  write.tree(phy2, file = paste(filename, "_with_bootstraps.tree",
                                sep = ""))
  phy2$node.label <- NULL
  write.tree(phy2, file = paste(filename, "_input.tree", sep = ""))
  ChromEvol::write.param(filename, mainType, maxChrNum, minChrNum, branchMul,
                         rootAt, numberSim, fix.root)
  if (mainType != "All_Models")
    specify.model(mainType, filename)
  cat("\tStructure of the input table\n")
  # system(paste("open", fnp))
  list(table3, phy2)
  cat("number of tips of the original tree:", length(phy$tip.label),
      "\n", "number of tips of the ChromEvol input tree:",
      length(phy2$tip.label), "\n", length(table[, 1]) - length(table3[,
                                                                       1]), "species have been dropped from the table\n",
      length(phy$tip.label) - length(phy2$tip.label), "tips have been dropped from tree:\n")
  cat(phy$tip.label[which(is.na(table2[, 2]))], "\n", sep = "\n")
  if (!is.logical(fix.root))
    write.fix.root(fix.root)
  # setwd("~/R")
  invisible(list(table3, phy2))
}
