# For parallel processing. for a serial run, do "cores <- 1"
suppressMessages(library("foreach"))
suppressMessages(library("doParallel"))
library("ggplot2")
source("bin/TaxonSampling/TaxonSampling.R")

#Number of bootstraps
n <- 10

#x axis values when plotting results
x <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)

#where should we start looking?
root_taxon <- 40674

#taxon ids to sequence names
idsFile <- "data/validation/taxonIDs_2_sequenceIDs.txt"

#fasta files
multifasta <- "data/validation/sequences_with_taxonIDs.fasta"

#NCBI taxonomy files
taxondir <- "taxdump/"

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

oper <- foreach(i=1:15, .combine='comb', .multicombine=TRUE,
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
                        "4" = numeric(0))

for (i in 1:n) {
  listOutput <- list()
  listOutput$"10" <- TS_Algorithm(root_taxon, 10, nodes, countIDs, method, randomize)
  listOutput$"20" <- TS_Algorithm(root_taxon, 20, nodes, countIDs, method, randomize)
  listOutput$"30" <- TS_Algorithm(root_taxon, 30, nodes, countIDs, method, randomize)
  listOutput$"40" <- TS_Algorithm(root_taxon, 40, nodes, countIDs, method, randomize)
  listOutput$"50" <- TS_Algorithm(root_taxon, 50, nodes, countIDs, method, randomize)
  listOutput$"60" <- TS_Algorithm(root_taxon, 60, nodes, countIDs, method, randomize)
  listOutput$"70" <- TS_Algorithm(root_taxon, 70, nodes, countIDs, method, randomize)
  listOutput$"80" <- TS_Algorithm(root_taxon, 80, nodes, countIDs, method, randomize)
  listOutput$"90" <- TS_Algorithm(root_taxon, 90, nodes, countIDs, method, randomize)
  listOutput$"100" <- TS_Algorithm(root_taxon, 100, nodes, countIDs, method, randomize)

  evalOutput <- list()
  evalOutput$"10" <- Evaluate_TS(listOutput$"10", nodes, countIDs)
  evalOutput$"20" <- Evaluate_TS(listOutput$"20", nodes, countIDs)
  evalOutput$"30" <- Evaluate_TS(listOutput$"30", nodes, countIDs)
  evalOutput$"40" <- Evaluate_TS(listOutput$"40", nodes, countIDs)
  evalOutput$"50" <- Evaluate_TS(listOutput$"50", nodes, countIDs)
  evalOutput$"60" <- Evaluate_TS(listOutput$"60", nodes, countIDs)
  evalOutput$"70" <- Evaluate_TS(listOutput$"70", nodes, countIDs)
  evalOutput$"80" <- Evaluate_TS(listOutput$"80", nodes, countIDs)
  evalOutput$"90" <- Evaluate_TS(listOutput$"90", nodes, countIDs)
  evalOutput$"100" <- Evaluate_TS(listOutput$"100", nodes, countIDs)
  
  
  for (level in 1:4) {
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
                                   "4" = numeric(0))

for (i in 1:n) {
  listOutput <- list()
  listOutput$"10" <- TS_Algorithm(root_taxon, 10, nodes, countIDs, method, randomize)
  listOutput$"20" <- TS_Algorithm(root_taxon, 20, nodes, countIDs, method, randomize)
  listOutput$"30" <- TS_Algorithm(root_taxon, 30, nodes, countIDs, method, randomize)
  listOutput$"40" <- TS_Algorithm(root_taxon, 40, nodes, countIDs, method, randomize)
  listOutput$"50" <- TS_Algorithm(root_taxon, 50, nodes, countIDs, method, randomize)
  listOutput$"60" <- TS_Algorithm(root_taxon, 60, nodes, countIDs, method, randomize)
  listOutput$"70" <- TS_Algorithm(root_taxon, 70, nodes, countIDs, method, randomize)
  listOutput$"80" <- TS_Algorithm(root_taxon, 80, nodes, countIDs, method, randomize)
  listOutput$"90" <- TS_Algorithm(root_taxon, 90, nodes, countIDs, method, randomize)
  listOutput$"100" <- TS_Algorithm(root_taxon, 100, nodes, countIDs, method, randomize)
  
  evalOutput <- list()
  evalOutput$"10" <- Evaluate_TS(listOutput$"10", nodes, countIDs)
  evalOutput$"20" <- Evaluate_TS(listOutput$"20", nodes, countIDs)
  evalOutput$"30" <- Evaluate_TS(listOutput$"30", nodes, countIDs)
  evalOutput$"40" <- Evaluate_TS(listOutput$"40", nodes, countIDs)
  evalOutput$"50" <- Evaluate_TS(listOutput$"50", nodes, countIDs)
  evalOutput$"60" <- Evaluate_TS(listOutput$"60", nodes, countIDs)
  evalOutput$"70" <- Evaluate_TS(listOutput$"70", nodes, countIDs)
  evalOutput$"80" <- Evaluate_TS(listOutput$"80", nodes, countIDs)
  evalOutput$"90" <- Evaluate_TS(listOutput$"90", nodes, countIDs)
  evalOutput$"100" <- Evaluate_TS(listOutput$"100", nodes, countIDs)
  
  
  for (level in 1:4) {
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
                     "13" = numeric(0), "14" = numeric(0), "15" = numeric(0))
for (i in 1:n) {
  listRandom <- list()
  listRandom$"10" <- RandomSampling(idsFile, nodes, 50)
  listRandom$"20" <- RandomSampling(idsFile, nodes, 100)
  listRandom$"30" <- RandomSampling(idsFile, nodes, 150)
  listRandom$"40" <- RandomSampling(idsFile, nodes, 200)
  listRandom$"50" <- RandomSampling(idsFile, nodes, 250)
  listRandom$"60" <- RandomSampling(idsFile, nodes, 300)
  listRandom$"70" <- RandomSampling(idsFile, nodes, 350)
  listRandom$"80" <- RandomSampling(idsFile, nodes, 400)
  listRandom$"90" <- RandomSampling(idsFile, nodes, 450)
  listRandom$"100" <- RandomSampling(idsFile, nodes, 500)
  
  evalRandom <- list()
  evalRandom$"10" <- Evaluate_TS(listRandom$"10", nodes, countIDs)
  evalRandom$"20" <- Evaluate_TS(listRandom$"20", nodes, countIDs)
  evalRandom$"30" <- Evaluate_TS(listRandom$"30", nodes, countIDs)
  evalRandom$"40" <- Evaluate_TS(listRandom$"40", nodes, countIDs)
  evalRandom$"50" <- Evaluate_TS(listRandom$"50", nodes, countIDs)
  evalRandom$"60" <- Evaluate_TS(listRandom$"60", nodes, countIDs)
  evalRandom$"70" <- Evaluate_TS(listRandom$"70", nodes, countIDs)
  evalRandom$"80" <- Evaluate_TS(listRandom$"80", nodes, countIDs)
  evalRandom$"90" <- Evaluate_TS(listRandom$"90", nodes, countIDs)
  evalRandom$"100" <- Evaluate_TS(listRandom$"100", nodes, countIDs)
  
  for (level in 3:15) {
    diversity <- numeric(0)
    for (number in names(evalRandom)) {
      diversity <- c(diversity, length(evalRandom[[number]][[level]]))
    }
    outputRandom[[level]] <- rbind(outputRandom[[level]], diversity)
  }
}


totalTaxa <- numeric(0)
children <- 1
for (i in 1:15) {
  taxon <- as.integer(children)
  children <- nodes$id[is.element(nodes$parent, taxon) & 
                       !is.element(nodes$id, taxon)]
  children <- intersect(children, names(countIDs))
  totalTaxa <- c(totalTaxa, length(children))
}


save.image("Eval_final.RData")

confidence <- .995  # 99% = .995, 95% = .975

diversity_random_Means <- list()
diversity_random_CI <- list()

diversity_no_random_Means <- list()
diversity_no_random_CI <- list()

random_Means <- list()
random_CI <- list()

for (level in 3:15) {
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


for (level in 3:15) {
  df <- data.frame(x, diversity_random_Means = diversity_random_Means[[level]], 
                      diversity_no_random_Means = diversity_no_random_Means[[level]],
                      random_Means = random_Means[[level]],
                      totalTaxa = rep(totalTaxa[level], 10))
  
  imageName <- paste0("level", level, ".pdf")
  pdf(imageName)
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
