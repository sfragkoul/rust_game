library(data.table)


consensus_strain <- function(strains) {
  # strains: list of character vectors, each of length 6
  dt <- as.data.table(do.call(rbind, strains))
  
  result <- character(ncol(dt))
  
  for (j in seq_len(ncol(dt))) {
    col_values <- dt[[j]]
    counts <- as.data.table(table(col_values))
    setnames(counts, c("gene", "N"))
    
    max_count <- max(counts$N)
    top_genes <- counts[N == max_count, gene]
    
    if (max_count >= 3) {
      # Clear majority (3 or 4 times)
      result[j] <- top_genes[1]
      
    } else if (max_count == 2) {
      if (length(top_genes) == 1) {
        # One gene with 2, others all different (XXYW → X)
        result[j] <- top_genes
      } else {
        # True tie (WWYY → W/Y)
        result[j] <- paste(top_genes, collapse = "/")
      }
      
    } else {
      # All unique or evenly split → ambiguous
      result[j] <- paste(top_genes, collapse = "/")
    }
  }
  
  return(result)
}

make_strain <- function(s) {
  
  return(strsplit(s, "")[[1]])
}

combine_strains <- function(string1, string2, string3, string4){
  
  str1 <- make_strain(string1)
  str2 <- make_strain(string2)
  str3 <- make_strain(string3)
  str4 <- make_strain(string4)
  
  strains <- list(str1, str2, str3, str4)
  cons <- consensus_strain(strains)
  
  return(cons)
  
}


resutls = combine_strains("YYYWGX", 
                          "YYYWGX", 
                          "XYYGGW", 
                          "XYYGGW")



#------------------------------------------------------------------------------

score_strain <- function(strain) {
  g_count <- sum(strain == "G")
  y_count <- sum(strain == "Y")
  dist <- abs(g_count - 2) + abs(y_count - 4)  # Manhattan distance
  return(list(G = g_count, Y = y_count, score = dist))
}



# ----------------------------
# Generate all repeated 4-strain combinations
strain_group <- c("YYYWGX", "XYYGGW", "GWGYWY")
all_combos <- expand.grid(strain_group, strain_group, strain_group, strain_group,
                          stringsAsFactors = FALSE)

# ----------------------------
# Run combine_strains() and score each
results <- apply(all_combos, 1, function(row) {
  cons <- combine_strains(row[1], row[2], row[3], row[4])
  sc <- score_strain(cons)
  list(
    combination = row,
    consensus = cons,
    G = sc$G,
    Y = sc$Y,
    score = sc$score
  )
})

# ----------------------------
# Find best results (lowest distance to 2G + 4Y)
scores <- sapply(results, function(x) x$score)
best_results <- results[scores == min(scores)]

# ----------------------------
# Print best results
for (res in best_results) {
  cat("Combination:", paste(res$combination, collapse = ", "), "\n")
  cat("Consensus:  ", paste(res$consensus, collapse = " "), "\n")
  cat("G count:", res$G, "Y count:", res$Y, "Score:", res$score, "\n\n")
}