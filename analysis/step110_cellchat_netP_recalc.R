# step110_cellchat_netP_recalc.R
# Purpose: netP$strength / netP$weight are NULL in CellChat_ARDS.rds / CellChat_Sepsis.rds.
# This script recomputes both slots from the objects themselves:
#   netP$strength[p] = sum of netP$prob[,,p]  (total communication probability / information flow)
#   netP$weight[p]  = number of significant L-R pairs assigned to pathway p (from interactions CSV mapping)
# Results are written back into the rds objects AND exported as CSV tables.
# Differential network between the two conditions is also rebuilt from the Desktop-verified
# interaction tables (ARDS_interactions.csv / Sepsis_interactions.csv).
# Date: 2026-10-19

suppressMessages(library(CellChat))

cell_dir <- "08_singlecell/cellchat"
out_dir  <- cell_dir

recalc <- function(rds_file, csv_file, label) {
  obj <- readRDS(file.path(cell_dir, rds_file))
  prob_arr <- obj@netP$prob          # cell x cell x pathway
  pathways <- dimnames(prob_arr)[[3]]
  stopifnot(!is.null(pathways))
  dms <- dim(prob_arr)

  # --- strength: total information flow per pathway (sum over all cell-type pairs) ---
  strength <- vapply(seq_len(dms[3]), function(k) sum(prob_arr[, , k]), numeric(1))
  names(strength) <- pathways

  # --- weight: number of significant L-R pairs per pathway (mapping from interactions CSV) ---
  inter <- read.csv(file.path(cell_dir, csv_file), stringsAsFactors = FALSE)
  weight <- table(inter$pathway_name)
  weight <- weight[pathways]
  weight[is.na(weight)] <- 0
  names(weight) <- pathways

  # --- write slots back into object ---
  obj@netP$strength <- strength
  obj@netP$weight   <- as.numeric(weight)
  names(obj@netP$weight) <- pathways
  saveRDS(obj, file.path(cell_dir, rds_file))

  # --- export table ---
  tab <- data.frame(pathway = pathways,
                    strength = unname(round(strength, 10)),
                    weight_LR_pairs = as.numeric(weight),
                    stringsAsFactors = FALSE)
  tab <- tab[order(-tab$strength), ]
  out_csv <- file.path(out_dir, paste0("netP_strength_weight_", label, ".csv"))
  write.csv(tab, out_csv, row.names = FALSE)
  cat("=====", label, "=====\n")
  print(tab)
  cat("written:", out_csv, "\n")
  return(tab)
}

tabA <- recalc("CellChat_ARDS.rds",  "ARDS_interactions.csv",  "ARDS")
tabS <- recalc("CellChat_Sepsis.rds", "Sepsis_interactions.csv", "Sepsis")

# --- sanity check: strength from array vs sum of L-R probs in interactions CSV ---
chk <- function(csv_file, label) {
  inter <- read.csv(file.path(cell_dir, csv_file), stringsAsFactors = FALSE)
  agg <- aggregate(prob ~ pathway_name, data = inter, FUN = sum)
  cat("cross-check (LR-level sum) for", label, ":\n")
  print(agg[order(-agg$prob), ])
}
chk("ARDS_interactions.csv", "ARDS")
chk("Sepsis_interactions.csv", "Sepsis")

# --- differential network: per-pathway strength comparison ---
allp <- union(tabA$pathway, tabS$pathway)
mA <- setNames(tabA$strength, tabA$pathway); mS <- setNames(tabS$strength, tabS$pathway)
strength_ards  <- mA[allp]; strength_ards[is.na(strength_ards)] <- 0
strength_sepsis <- mS[allp]; strength_sepsis[is.na(strength_sepsis)] <- 0
wA <- setNames(tabA$weight_LR_pairs, tabA$pathway); wS <- setNames(tabS$weight_LR_pairs, tabS$pathway)

diff_net <- data.frame(
  pathway = allp,
  strength_ARDS  = unname(strength_ards),
  strength_Sepsis = unname(strength_sepsis),
  diff_Sepsis_minus_ARDS = unname(strength_sepsis - strength_ards),
  ratio_Sepsis_over_ARDS = ifelse(strength_ards > 0, unname(strength_sepsis / strength_ards), NA),
  dominant_side = ifelse(strength_sepsis > strength_ards, "Sepsis",
                  ifelse(strength_ards > strength_sepsis, "ARDS", "tie")),
  weight_LR_pairs_ARDS = unname(ifelse(is.na(wA[allp]), 0, wA[allp])),
  weight_LR_pairs_Sepsis = unname(ifelse(is.na(wS[allp]), 0, wS[allp])),
  stringsAsFactors = FALSE
)
diff_net <- diff_net[order(-abs(diff_net$diff_Sepsis_minus_ARDS)), ]
diff_net$rank_by_abs_diff <- seq_len(nrow(diff_net))
out_diff <- file.path(out_dir, "pathway_diff_network_recalc_2026-10-19.csv")
write.csv(diff_net, out_diff, row.names = FALSE)
cat("\n===== differential network (recalc) =====\n")
print(diff_net)
cat("written:", out_diff, "\n")

# --- MIF focus ---
mif <- diff_net[diff_net$pathway == "MIF", ]
cat("\nMIF: ARDS strength =", sprintf("%.6e", mif$strength_ARDS),
    "| Sepsis strength =", sprintf("%.6e", mif$strength_Sepsis),
    "| ratio =", sprintf("%.2f", mif$ratio_Sepsis_over_ARDS),
    "| dominant side =", mif$dominant_side, "\n")
