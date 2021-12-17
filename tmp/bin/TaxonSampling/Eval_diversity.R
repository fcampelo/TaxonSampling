# For parallel processing. for a serial run, do "cores <- 1"
suppressMessages(library("foreach"))
suppressMessages(library("doParallel"))
library("ggplot2")
library("ggpubr")
source("bin/TaxonSampling/TaxonSampling.R")

#Number of bootstraps
n <- 10

#x axis values when plotting results
x <- c(50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 750, 1000)

#where should we start looking?
root_taxon <- 1

#taxon ids to sequence names
idsFile <- "data/validation/metadata/TaxID2SeqID.txt"

#fasta files
multifasta <- "data/validation/fasta/mit_vertebrata.fasta"

#NCBI taxonomy files
taxondir <- "data/validation/taxdump/"

#loading node structure from NCBI Taxonomy
nodes <- suppressMessages(getnodes(taxondir))

#number of taxonomic IDs per node
countIDs <- TS_TaxonomyData(idsFile, nodes)
nodes <- Simplify_Nodes(nodes, countIDs)

cores <- 7
if (cores > 1) {
  cl <- makeCluster(cores)
  registerDoParallel(cl)
}
#
comb <- function(x, ...) {
  lapply(seq_along(x),
    function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

oper <- foreach(i=1:35, .combine='comb', .multicombine=TRUE,
                .init=list(list(), list())) %dopar% {
  list(i+2, i+3)
}

foreach (i = 1:n, .export = c("TS_Algorithm", "RandomSampling",
                              "Evaluate_TS", "TS_TaxonomyData",
                              "TS_Algorithm")) %dopar% {
                              }

# Reduce the node information to the necessary only, reduces search time.
nodes <- nodes[is.element(nodes$id, names(countIDs)), 1:2]


method <- "diversity"
randomize <- "no"

output_diversity_no_random <- list("1" = numeric(0), "2" = numeric(0), "3" = numeric(0),
                                   "4" = numeric(0), "5" = numeric(0), "6" = numeric(0),
                                   "7" = numeric(0), "8" = numeric(0), "9" = numeric(0), 
                                   "10" = numeric(0), "11" = numeric(0), "12" = numeric(0),
                                   "13" = numeric(0), "14" = numeric(0), "15" = numeric(0),
                                   "16" = numeric(0), "17" = numeric(0), "18" = numeric(0),
                                   "19" = numeric(0), "20" = numeric(0), "21" = numeric(0),
                                   "22" = numeric(0), "23" = numeric(0), "24" = numeric(0),
                                   "25" = numeric(0), "26" = numeric(0), "27" = numeric(0),
                                   "28" = numeric(0), "29" = numeric(0), "30" = numeric(0),
                                   "31" = numeric(0), "32" = numeric(0), "33" = numeric(0),
                                   "34" = numeric(0), "35" = numeric(0))


for (i in 1:n) {
  print(paste0("Bootstrap ", i, " of ", n))
  listOutput <- list()
  listOutput$"50" <- TS_Algorithm(root_taxon, 50, nodes, countIDs, method, randomize)
  listOutput$"100" <- TS_Algorithm(root_taxon, 100, nodes, countIDs, method, randomize)
  listOutput$"150" <- TS_Algorithm(root_taxon, 150, nodes, countIDs, method, randomize)
  listOutput$"200" <- TS_Algorithm(root_taxon, 200, nodes, countIDs, method, randomize)
  listOutput$"250" <- TS_Algorithm(root_taxon, 250, nodes, countIDs, method, randomize)
  listOutput$"300" <- TS_Algorithm(root_taxon, 300, nodes, countIDs, method, randomize)
  listOutput$"350" <- TS_Algorithm(root_taxon, 350, nodes, countIDs, method, randomize)
  listOutput$"400" <- TS_Algorithm(root_taxon, 400, nodes, countIDs, method, randomize)
  listOutput$"450" <- TS_Algorithm(root_taxon, 450, nodes, countIDs, method, randomize)
  listOutput$"500" <- TS_Algorithm(root_taxon, 500, nodes, countIDs, method, randomize)
  listOutput$"750" <- TS_Algorithm(root_taxon, 750, nodes, countIDs, method, randomize)
  listOutput$"1000" <- TS_Algorithm(root_taxon, 1000, nodes, countIDs, method, randomize)
  
  evalOutput <- list()
  evalOutput$"50" <- Evaluate_TS(listOutput$"50", nodes, countIDs)
  evalOutput$"100" <- Evaluate_TS(listOutput$"100", nodes, countIDs)
  evalOutput$"150" <- Evaluate_TS(listOutput$"150", nodes, countIDs)
  evalOutput$"200" <- Evaluate_TS(listOutput$"200", nodes, countIDs)
  evalOutput$"250" <- Evaluate_TS(listOutput$"250", nodes, countIDs)
  evalOutput$"300" <- Evaluate_TS(listOutput$"300", nodes, countIDs)
  evalOutput$"350" <- Evaluate_TS(listOutput$"350", nodes, countIDs)
  evalOutput$"400" <- Evaluate_TS(listOutput$"400", nodes, countIDs)
  evalOutput$"450" <- Evaluate_TS(listOutput$"450", nodes, countIDs)
  evalOutput$"500" <- Evaluate_TS(listOutput$"500", nodes, countIDs)
  evalOutput$"750" <- Evaluate_TS(listOutput$"750", nodes, countIDs)
  evalOutput$"1000" <- Evaluate_TS(listOutput$"1000", nodes, countIDs)
  
  
  for (level in 1:35) {
    no_random <- numeric(0)
    for (number in names(evalOutput)) {
      no_random <- c(no_random, length(evalOutput[[number]][[level]]))
    }
    output_diversity_no_random[[level]] <- rbind(output_diversity_no_random[[level]], no_random)
  }
}

method <- "diversity"
randomize <- "yes"

output_diversity_random <- list("1" = numeric(0), "2" = numeric(0), "3" = numeric(0),
                                "4" = numeric(0), "5" = numeric(0), "6" = numeric(0),
                                "7" = numeric(0), "8" = numeric(0), "9" = numeric(0), 
                                "10" = numeric(0), "11" = numeric(0), "12" = numeric(0),
                                "13" = numeric(0), "14" = numeric(0), "15" = numeric(0),
                                "16" = numeric(0), "17" = numeric(0), "18" = numeric(0),
                                "19" = numeric(0), "20" = numeric(0), "21" = numeric(0),
                                "22" = numeric(0), "23" = numeric(0), "24" = numeric(0),
                                "25" = numeric(0), "26" = numeric(0), "27" = numeric(0),
                                "28" = numeric(0), "29" = numeric(0), "30" = numeric(0),
                                "31" = numeric(0), "32" = numeric(0), "33" = numeric(0),
                                "34" = numeric(0), "35" = numeric(0))


for (i in 1:n) {
  print(paste0("Bootstrap ", i, " of ", n))
  listOutput <- list()
  listOutput$"50" <- TS_Algorithm(root_taxon, 50, nodes, countIDs, method, randomize)
  listOutput$"100" <- TS_Algorithm(root_taxon, 100, nodes, countIDs, method, randomize)
  listOutput$"150" <- TS_Algorithm(root_taxon, 150, nodes, countIDs, method, randomize)
  listOutput$"200" <- TS_Algorithm(root_taxon, 200, nodes, countIDs, method, randomize)
  listOutput$"250" <- TS_Algorithm(root_taxon, 250, nodes, countIDs, method, randomize)
  listOutput$"300" <- TS_Algorithm(root_taxon, 300, nodes, countIDs, method, randomize)
  listOutput$"350" <- TS_Algorithm(root_taxon, 350, nodes, countIDs, method, randomize)
  listOutput$"400" <- TS_Algorithm(root_taxon, 400, nodes, countIDs, method, randomize)
  listOutput$"450" <- TS_Algorithm(root_taxon, 450, nodes, countIDs, method, randomize)
  listOutput$"500" <- TS_Algorithm(root_taxon, 500, nodes, countIDs, method, randomize)
  listOutput$"750" <- TS_Algorithm(root_taxon, 750, nodes, countIDs, method, randomize)
  listOutput$"1000" <- TS_Algorithm(root_taxon, 1000, nodes, countIDs, method, randomize)
  
  evalOutput <- list()
  evalOutput$"50" <- Evaluate_TS(listOutput$"50", nodes, countIDs)
  evalOutput$"100" <- Evaluate_TS(listOutput$"100", nodes, countIDs)
  evalOutput$"150" <- Evaluate_TS(listOutput$"150", nodes, countIDs)
  evalOutput$"200" <- Evaluate_TS(listOutput$"200", nodes, countIDs)
  evalOutput$"250" <- Evaluate_TS(listOutput$"250", nodes, countIDs)
  evalOutput$"300" <- Evaluate_TS(listOutput$"300", nodes, countIDs)
  evalOutput$"350" <- Evaluate_TS(listOutput$"350", nodes, countIDs)
  evalOutput$"400" <- Evaluate_TS(listOutput$"400", nodes, countIDs)
  evalOutput$"450" <- Evaluate_TS(listOutput$"450", nodes, countIDs)
  evalOutput$"500" <- Evaluate_TS(listOutput$"500", nodes, countIDs)
  evalOutput$"750" <- Evaluate_TS(listOutput$"750", nodes, countIDs)
  evalOutput$"1000" <- Evaluate_TS(listOutput$"1000", nodes, countIDs)
  
  
  for (level in 1:35) {
    random <- numeric(0)
    for (number in names(evalOutput)) {
      random <- c(random, length(evalOutput[[number]][[level]]))
    }
    output_diversity_random[[level]] <- rbind(output_diversity_random[[level]], random)
  }
}


outputRandom <- list("1" = numeric(0), "2" = numeric(0), "3" = numeric(0),
                     "4" = numeric(0), "5" = numeric(0), "6" = numeric(0),
                     "7" = numeric(0), "8" = numeric(0), "9" = numeric(0), 
                     "10" = numeric(0), "11" = numeric(0), "12" = numeric(0),
                     "13" = numeric(0), "14" = numeric(0), "15" = numeric(0),
                     "16" = numeric(0), "17" = numeric(0), "18" = numeric(0),
                     "19" = numeric(0), "20" = numeric(0), "21" = numeric(0),
                     "22" = numeric(0), "23" = numeric(0), "24" = numeric(0),
                     "25" = numeric(0), "26" = numeric(0), "27" = numeric(0),
                     "28" = numeric(0), "29" = numeric(0), "30" = numeric(0),
                     "31" = numeric(0), "32" = numeric(0), "33" = numeric(0),
                     "34" = numeric(0), "35" = numeric(0))

for (i in 1:n) {
  print(paste0("Bootstrap ", i, " of ", n))
  listRandom <- list()
  listRandom$"50" <- RandomSampling(idsFile, nodes, 50)
  listRandom$"100" <- RandomSampling(idsFile, nodes, 100)
  listRandom$"150" <- RandomSampling(idsFile, nodes, 150)
  listRandom$"200" <- RandomSampling(idsFile, nodes, 200)
  listRandom$"250" <- RandomSampling(idsFile, nodes, 250)
  listRandom$"300" <- RandomSampling(idsFile, nodes, 300)
  listRandom$"350" <- RandomSampling(idsFile, nodes, 350)
  listRandom$"400" <- RandomSampling(idsFile, nodes, 400)
  listRandom$"450" <- RandomSampling(idsFile, nodes, 450)
  listRandom$"500" <- RandomSampling(idsFile, nodes, 500)
  listRandom$"750" <- RandomSampling(idsFile, nodes, 750)
  listRandom$"1000" <- RandomSampling(idsFile, nodes, 1000)
  
  evalRandom <- list()
  evalRandom$"50" <- Evaluate_TS(listRandom$"50", nodes, countIDs)
  evalRandom$"100" <- Evaluate_TS(listRandom$"100", nodes, countIDs)
  evalRandom$"150" <- Evaluate_TS(listRandom$"150", nodes, countIDs)
  evalRandom$"200" <- Evaluate_TS(listRandom$"200", nodes, countIDs)
  evalRandom$"250" <- Evaluate_TS(listRandom$"250", nodes, countIDs)
  evalRandom$"300" <- Evaluate_TS(listRandom$"300", nodes, countIDs)
  evalRandom$"350" <- Evaluate_TS(listRandom$"350", nodes, countIDs)
  evalRandom$"400" <- Evaluate_TS(listRandom$"400", nodes, countIDs)
  evalRandom$"450" <- Evaluate_TS(listRandom$"450", nodes, countIDs)
  evalRandom$"500" <- Evaluate_TS(listRandom$"500", nodes, countIDs)
  evalRandom$"750" <- Evaluate_TS(listRandom$"750", nodes, countIDs)
  evalRandom$"1000" <- Evaluate_TS(listRandom$"1000", nodes, countIDs)
  
  for (level in 1:35) {
    diversity <- numeric(0)
    for (number in names(evalRandom)) {
      diversity <- c(diversity, length(evalRandom[[number]][[level]]))
    }
    outputRandom[[level]] <- rbind(outputRandom[[level]], diversity)
  }
}


totalTaxa <- numeric(0)
children <- 1
for (i in 1:35) {
  taxon <- as.integer(children)
  children <- nodes$id[is.element(nodes$parent, taxon) & 
                       !is.element(nodes$id, taxon)]
  children <- intersect(children, names(countIDs))
  totalTaxa <- c(totalTaxa, length(children))
}


save.image("Eval_mito.RData")

confidence <- .995  # 99% = .995, 95% = .975

diversity_random_Means <- list()
diversity_random_CI <- list()

diversity_no_random_Means <- list()
diversity_no_random_CI <- list()

random_Means <- list()
random_CI <- list()

for (level in 1:35) {
  diversity_random_Means[[level]] <- colMeans(output_diversity_random[[level]])
  diversity_no_random_Means[[level]] <- colMeans(output_diversity_no_random[[level]])
  random_Means[[level]] <- colMeans(outputRandom[[level]])

  diversity_random_CI[[level]] <- apply(output_diversity_random[[level]], 2, sd)
  diversity_no_random_CI[[level]] <- apply(output_diversity_no_random[[level]], 2, sd)
  random_CI[[level]] <- apply(outputRandom[[level]], 2, sd)

  diversity_random_CI[[level]] <- diversity_random_CI[[level]]/sqrt(n)
  diversity_no_random_CI[[level]] <- diversity_no_random_CI[[level]]/sqrt(n)
  random_CI[[level]] <- random_CI[[level]]/sqrt(n)

  diversity_random_CI[[level]] <- qt(confidence, df = n-1) * diversity_random_CI[[level]]
  diversity_no_random_CI[[level]] <- qt(confidence, df = n-1) * diversity_no_random_CI[[level]]
  random_CI[[level]] <- qt(confidence, df = n-1) * random_CI[[level]]
}


level_15 <- ggplot()
level_17 <- ggplot()
level_19 <- ggplot()
level_21 <- ggplot()

level <- 15

  df <- data.frame(x, diversity_random_Means = diversity_random_Means[[level]], 
                   diversity_no_random_Means = diversity_no_random_Means[[level]],
                   random_Means = random_Means[[level]],
                   totalTaxa = rep(totalTaxa[level], 12))
  
  level_15 <- (ggplot(df, aes(x)) + 
          geom_point(aes(y=diversity_random_Means, colour="TS_div_rand")) +
          geom_line(aes(y=diversity_random_Means, colour="TS_div_rand")) + 
          geom_errorbar(aes(ymin=diversity_random_Means - diversity_random_CI[[level]],
                            ymax=diversity_random_Means + diversity_random_CI[[level]],
                            colour = "TS_div_rand"), width=1) +
          geom_point(aes(y=diversity_no_random_Means, colour="TS_div_non_rand")) +
          geom_line(aes(y=diversity_no_random_Means, colour="TS_div_non_rand")) + 
          geom_errorbar(aes(ymin=diversity_no_random_Means - diversity_no_random_CI[[level]],
                            ymax=diversity_no_random_Means + diversity_no_random_CI[[level]],
                            colour = "TS_div_non_rand"), width=1) +
          geom_point(aes(y=random_Means, colour="Random_sampling")) +
          geom_line(aes(y=random_Means, colour="Random_sampling")) +
          geom_errorbar(aes(ymin=random_Means - random_CI[[level]],
                            ymax=random_Means + random_CI[[level]],
                            colour = "Random_sampling"), width=1) +
          geom_point(aes(y=totalTaxa, colour="max")) +
          geom_line(aes(y=totalTaxa, colour="max")) +
          xlab("m") + ylab(paste0("# taxa (level = ", level, ")")) +
          scale_color_manual("Method",
                             values = c("TS_div_rand" = "orange",
                                        "TS_div_non_rand" = "red",
                                        "Random_sampling" = "blue",
                                        "max" = "black"))) 

  
  level <- 17
  
  df <- data.frame(x, diversity_random_Means = diversity_random_Means[[level]], 
                   diversity_no_random_Means = diversity_no_random_Means[[level]],
                   random_Means = random_Means[[level]],
                   totalTaxa = rep(totalTaxa[level], 12))
  
  level_17 <- (ggplot(df, aes(x)) + 
                 geom_point(aes(y=diversity_random_Means, colour="TS_div_rand")) +
                 geom_line(aes(y=diversity_random_Means, colour="TS_div_rand")) + 
                 geom_errorbar(aes(ymin=diversity_random_Means - diversity_random_CI[[level]],
                                   ymax=diversity_random_Means + diversity_random_CI[[level]],
                                   colour = "TS_div_rand"), width=1) +
                 geom_point(aes(y=diversity_no_random_Means, colour="TS_div_non_rand")) +
                 geom_line(aes(y=diversity_no_random_Means, colour="TS_div_non_rand")) + 
                 geom_errorbar(aes(ymin=diversity_no_random_Means - diversity_no_random_CI[[level]],
                                   ymax=diversity_no_random_Means + diversity_no_random_CI[[level]],
                                   colour = "TS_div_non_rand"), width=1) +
                 geom_point(aes(y=random_Means, colour="Random_sampling")) +
                 geom_line(aes(y=random_Means, colour="Random_sampling")) +
                 geom_errorbar(aes(ymin=random_Means - random_CI[[level]],
                                   ymax=random_Means + random_CI[[level]],
                                   colour = "Random_sampling"), width=1) +
                 geom_point(aes(y=totalTaxa, colour="max")) +
                 geom_line(aes(y=totalTaxa, colour="max")) +
                 xlab("m") + ylab(paste0("# taxa (level = ", level, ")")) +
                 scale_color_manual("Method",
                                    values = c("TS_div_rand" = "orange",
                                               "TS_div_non_rand" = "red",
                                               "Random_sampling" = "blue",
                                               "max" = "black"))) 
  
  
  level <- 19
  
  df <- data.frame(x, diversity_random_Means = diversity_random_Means[[level]], 
                   diversity_no_random_Means = diversity_no_random_Means[[level]],
                   random_Means = random_Means[[level]],
                   totalTaxa = rep(totalTaxa[level], 12))

  level_19 <- (ggplot(df, aes(x)) + 
                 geom_point(aes(y=diversity_random_Means, colour="TS_div_rand")) +
                 geom_line(aes(y=diversity_random_Means, colour="TS_div_rand")) + 
                 geom_errorbar(aes(ymin=diversity_random_Means - diversity_random_CI[[level]],
                                   ymax=diversity_random_Means + diversity_random_CI[[level]],
                                   colour = "TS_div_rand"), width=1) +
                 geom_point(aes(y=diversity_no_random_Means, colour="TS_div_non_rand")) +
                 geom_line(aes(y=diversity_no_random_Means, colour="TS_div_non_rand")) + 
                 geom_errorbar(aes(ymin=diversity_no_random_Means - diversity_no_random_CI[[level]],
                                   ymax=diversity_no_random_Means + diversity_no_random_CI[[level]],
                                   colour = "TS_div_non_rand"), width=1) +
                 geom_point(aes(y=random_Means, colour="Random_sampling")) +
                 geom_line(aes(y=random_Means, colour="Random_sampling")) +
                 geom_errorbar(aes(ymin=random_Means - random_CI[[level]],
                                   ymax=random_Means + random_CI[[level]],
                                   colour = "Random_sampling"), width=1) +
                 geom_point(aes(y=totalTaxa, colour="max")) +
                 geom_line(aes(y=totalTaxa, colour="max")) +
                 xlab("m") + ylab(paste0("# taxa (level = ", level, ")")) +
                 scale_color_manual("Method",
                                    values = c("TS_div_rand" = "orange",
                                               "TS_div_non_rand" = "red",
                                               "Random_sampling" = "blue",
                                               "max" = "black")))
  
  
  level <- 21
  
  df <- data.frame(x, diversity_random_Means = diversity_random_Means[[level]], 
                   diversity_no_random_Means = diversity_no_random_Means[[level]],
                   random_Means = random_Means[[level]],
                   totalTaxa = rep(totalTaxa[level], 12))
  
  level_21 <- (ggplot(df, aes(x)) + 
                 geom_point(aes(y=diversity_random_Means, colour="TS_div_rand")) +
                 geom_line(aes(y=diversity_random_Means, colour="TS_div_rand")) + 
                 geom_errorbar(aes(ymin=diversity_random_Means - diversity_random_CI[[level]],
                                   ymax=diversity_random_Means + diversity_random_CI[[level]],
                                   colour = "TS_div_rand"), width=1) +
                 geom_point(aes(y=diversity_no_random_Means, colour="TS_div_non_rand")) +
                 geom_line(aes(y=diversity_no_random_Means, colour="TS_div_non_rand")) + 
                 geom_errorbar(aes(ymin=diversity_no_random_Means - diversity_no_random_CI[[level]],
                                   ymax=diversity_no_random_Means + diversity_no_random_CI[[level]],
                                   colour = "TS_div_non_rand"), width=1) +
                 geom_point(aes(y=random_Means, colour="Random_sampling")) +
                 geom_line(aes(y=random_Means, colour="Random_sampling")) +
                 geom_errorbar(aes(ymin=random_Means - random_CI[[level]],
                                   ymax=random_Means + random_CI[[level]],
                                   colour = "Random_sampling"), width=1) +
                 geom_point(aes(y=totalTaxa, colour="max")) +
                 geom_line(aes(y=totalTaxa, colour="max")) +
                 xlab("m") + ylab(paste0("# taxa (level = ", level, ")")) +
                 scale_color_manual("Method",
                                    values = c("TS_div_rand" = "orange",
                                               "TS_div_non_rand" = "red",
                                               "Random_sampling" = "blue",
                                               "max" = "black")))


figure <- ggarrange(level_15, level_17, level_19, level_21, ncol=2, nrow=2, common.legend = TRUE)  

pdf("teste.pdf", width=8, height=8)
print(figure)
dev.off()

for (level in 1:35) {
  df <- data.frame(x, diversity_random_Means = diversity_random_Means[[level]], 
                      diversity_no_random_Means = diversity_no_random_Means[[level]],
                      random_Means = random_Means[[level]],
                      totalTaxa = rep(totalTaxa[level], 12))
  
  imageName <- paste0("level", level, ".pdf")
  pdf(imageName, width=5, height=5)
  print(ggplot(df, aes(x)) + 
          geom_point(aes(y=diversity_random_Means, colour="TS_div_rand")) +
          geom_line(aes(y=diversity_random_Means, colour="TS_div_rand")) + 
          geom_errorbar(aes(ymin=diversity_random_Means - diversity_random_CI[[level]],
                            ymax=diversity_random_Means + diversity_random_CI[[level]],
                            colour = "TS_div_rand"), width=1) +
          geom_point(aes(y=diversity_no_random_Means, colour="TS_div_non_rand")) +
          geom_line(aes(y=diversity_no_random_Means, colour="TS_div_non_rand")) + 
          geom_errorbar(aes(ymin=diversity_no_random_Means - diversity_no_random_CI[[level]],
                            ymax=diversity_no_random_Means + diversity_no_random_CI[[level]],
                            colour = "TS_div_non_rand"), width=1) +
          geom_point(aes(y=random_Means, colour="Random_sampling")) +
          geom_line(aes(y=random_Means, colour="Random_sampling")) +
          geom_errorbar(aes(ymin=random_Means - random_CI[[level]],
                            ymax=random_Means + random_CI[[level]],
                            colour = "Random_sampling"), width=1) +
          geom_point(aes(y=totalTaxa, colour="max")) +
          geom_line(aes(y=totalTaxa, colour="max")) +
          xlab("m") + ylab(paste0("# taxa (level = ", level, ")")) +
          scale_color_manual("Method",
                             values = c("TS_div_rand" = "orange",
                                        "TS_div_non_rand" = "red",
                                        "Random_sampling" = "blue",
                                        "max" = "black")))
  dev.off()
}
