# step111_cellchat_GSE167363_netP_recalc.R
# Purpose: give the CD14+ monocyte cohort (GSE167363) the same quantitative CellChat
# support that the PBMC cohort (GSE151263) already has (step110).
# Objects: 08_singlecell/GSE167363/cellchat/CellChat_Sepsis.rds (sepsis),
#          CellChat_HC.rds (healthy control), CellChat_merged.rds
#
# Object slots actually available (verified 2026-10-19):
#   netP$pathways, netP$prob (7 x 7 x 35 array: celltype x celltype x pathway)
#   net$prob, net$pval, net$LRs, net$count, net$weight   (LR x cell-pair matrices)
#   DB$interaction                                        (LR -> pathway mapping)
#   idents: ClassicalMono, InflammatoryMono, IntermediateMono, ISGMono,
#           MonoDC, NonclassicalMono, pDC
#
# netP$strength / netP$weight are NULL, so:
#   strength[p] = sum over all cell-type pairs of netP$prob[,,p]   (information flow)
#   weight[p]   = number of significant (pval <= 0.05) L-R pairs assigned to pathway p
# Both are written back into the rds objects and exported as CSV.
# Date: 2026-10-19

suppressMessages(library(CellChat))

cell_dir <- "08_singlecell/GSE167363/cellchat"

recalc_167363 <- function(rds_file, label) {
  obj <- readRDS(file.path(cell_dir, rds_file))
  prob_arr <- obj@netP$prob
  if (is.null(prob_arr)) stop("netP$prob missing in ", rds_file)
  pathways <- dimnames(prob_arr)[[3]]
  dms <- dim(prob_arr)

  # --- strength: total information flow per pathway (sum over all cell-type pairs) ---
  strength <- vapply(seq_len(dms[3]), function(k) sum(prob_arr[, , k]), numeric(1))
  names(strength) <- pathways

  # --- weight: number of significant L-R pairs per pathway ---
  # net$prob / net$pval are 7 x 7 x nLR arrays (celltype x celltype x L-R pair)
  lr_names  <- dimnames(obj@net$prob)[[3]]
  lr_pvals  <- obj@net$pval
  sig_lr <- character(0)
  if (!is.null(lr_pvals) && !is.null(lr_names)) {
    min_p <- apply(lr_pvals, 3, min, na.rm = TRUE)      # LR significant in any cell pair
    sig_lr <- lr_names[min_p <= 0.05]
  }
  path_map <- obj@DB$interaction
  idx <- match(toupper(sig_lr), toupper(path_map$interaction_name))
  pw_of_lr <- path_map$pathway_name[idx]
  pw_of_lr <- pw_of_lr[!is.na(pw_of_lr)]
  weight <- as.numeric(table(factor(pw_of_lr, levels = pathways)))
  names(weight) <- pathways

  # --- write slots back ---
  obj@netP$strength <- strength
  obj@netP$weight   <- weight
  names(obj@netP$strength) <- pathways
  names(obj@netP$weight)   <- pathways
  saveRDS(obj, file.path(cell_dir, rds_file))

  tab <- data.frame(pathway = pathways,
                    strength = unname(round(strength, 10)),
                    weight_LR_pairs = weight,
                    stringsAsFactors = FALSE)
  tab <- tab[order(-tab$strength), ]
  out_csv <- file.path(cell_dir, paste0("netP_strength_weight_", label, ".csv"))
  write.csv(tab, out_csv, row.names = FALSE)
  cat("===== GSE167363 /", label, "=====\n")
  print(tab)
  cat("written:", out_csv, "\n")

  # --- cross-check: sum of LR-level probabilities per pathway ---
  lr_prob <- obj@net$prob
  lr_sum  <- apply(lr_prob, 3, sum, na.rm = TRUE)       # total probability per L-R pair
  pw_all  <- path_map$pathway_name[match(toupper(lr_names),
                                         toupper(path_map$interaction_name))]
  df <- data.frame(pathway = pw_all, lr_sum = unname(lr_sum), stringsAsFactors = FALSE)
  df <- df[!is.na(df$pathway), ]
  if (nrow(df) > 0) {
    by_p <- aggregate(lr_sum ~ pathway, data = df, FUN = sum)
    cat("cross-check (sum of L-R pair probabilities per pathway), top 8:\n")
    print(head(by_p[order(-by_p$lr_sum), ], 8))
  } else {
    cat("cross-check: no L-R pairs could be mapped to pathways\n")
  }
  cat("\n")
  return(tab)
}

tabS <- recalc_167363("CellChat_Sepsis.rds", "Sepsis")
tabH <- recalc_167363("CellChat_HC.rds",     "HC")

# --- differential network: sepsis vs healthy control ---
allp <- union(tabS$pathway, tabH$pathway)
mS <- setNames(tabS$strength, tabS$pathway); mH <- setNames(tabH$strength, tabH$pathway)
sS <- mS[allp]; sS[is.na(sS)] <- 0
sH <- mH[allp]; sH[is.na(sH)] <- 0
wS <- setNames(tabS$weight_LR_pairs, tabS$pathway); wH <- setNames(tabH$weight_LR_pairs, tabH$pathway)

diff_net <- data.frame(
  pathway = allp,
  strength_Sepsis = unname(sS),
  strength_HC     = unname(sH),
  diff_Sepsis_minus_HC = unname(sS - sH),
  ratio_Sepsis_over_HC = ifelse(sH > 0, unname(sS / sH), NA),
  dominant_side = ifelse(sS > sH, "Sepsis", ifelse(sH > sS, "HC", "tie")),
  weight_LR_pairs_Sepsis = unname(ifelse(is.na(wS[allp]), 0, wS[allp])),
  weight_LR_pairs_HC     = unname(ifelse(is.na(wH[allp]), 0, wH[allp])),
  stringsAsFactors = FALSE
)
diff_net <- diff_net[order(-abs(diff_net$diff_Sepsis_minus_HC)), ]
diff_net$rank_by_abs_diff <- seq_len(nrow(diff_net))
out_diff <- file.path(cell_dir, "pathway_diff_network_GSE167363_recalc_2026-10-19.csv")
write.csv(diff_net, out_diff, row.names = FALSE)
cat("\n===== GSE167363 differential network (sepsis vs HC) =====\n")
print(diff_net)
cat("written:", out_diff, "\n")

cat("\n--- top-5 pathways by strength, each condition ---\n")
cat("Sepsis:\n"); print(head(tabS, 5))
cat("HC:\n");     print(head(tabH, 5))

cat("\n--- MIF / MHC-I / MHC-II focus ---\n")
for (pw in c("MIF", "MHC-I", "MHC-II", "GALECTIN", "CypA")) {
  r <- diff_net[diff_net$pathway == pw, ]
  if (nrow(r) == 1) {
    cat(sprintf("%-10s Sepsis=%.4e  HC=%.4e  ratio=%.2f  dominant=%s\n",
                pw, r$strength_Sepsis, r$strength_HC,
                r$ratio_Sepsis_over_HC, r$dominant_side))
  }
}
