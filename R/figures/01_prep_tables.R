# ---------------------------------------------------------------------------
# 01_prep_tables.R -- build the small derived tables the figure scripts read
#
# The analysis pipeline writes one directory per module; the figure scripts
# should not depend on the multi-gigabyte raw GEO downloads. This script
# extracts the few matrices the figures need and stores them as tidy CSV under
# data/derived. It is the only place where raw expression matrices are touched.
#
# Inputs  (ARDS_ROOT, default ../../ARDS_Sepsis_reanalysis):
#   01_DEG/{GSE76293,GSE69528}_{gene_expr,grp}.rds
#   00_rawdata/validation/GSE100159/{GSE100159_series_matrix.txt,GPL6884.annot}
#   00_rawdata/validation/GSE243066/*_raw_count.txt.gz
# Output  data/derived/hub6_expression_long.csv
# ---------------------------------------------------------------------------

source("R/00_setup.R")

ARDS_ROOT <- Sys.getenv("ARDS_ROOT", unset = "../../ARDS_Sepsis_reanalysis")
DERIVED   <- file.path("data", "derived", "hub6_expression_long.csv")
if (!dir.exists(ARDS_ROOT)) {
  # 原始矩阵不在本仓库内（体积大且可重新下载）。data/derived 里已附带由它
  # 生成的成品表，图件层可以直接跑；需要重建时才用 ARDS_ROOT 指到原始项目。
  if (file.exists(DERIVED)) {
    message("ARDS_ROOT not found; keeping the shipped ", DERIVED)
    quit(save = "no", status = 0)
  }
  stop("ARDS_ROOT is not set and ", DERIVED, " is missing")
}
dir.create(file.path("data", "derived"), recursive = TRUE, showWarnings = FALSE)

# --- discovery cohorts: normalised matrices are already stored ---------------

load_discovery <- function(ds) {
  expr <- readRDS(file.path(ARDS_ROOT, "01_DEG", paste0(ds, "_gene_expr.rds")))
  grp  <- readRDS(file.path(ARDS_ROOT, "01_DEG", paste0(ds, "_grp.rds")))$group
  list(expr = expr, grp = grp)
}

tidy <- function(x, cohort) {
  m <- x$expr[HUB6, , drop = FALSE]
  data.frame(
    cohort = cohort,
    gene   = rep(rownames(m), times = ncol(m)),
    sample = rep(colnames(m), each = nrow(m)),
    group  = rep(x$grp, each = nrow(m)),
    expr   = as.vector(m),
    stringsAsFactors = FALSE
  )
}

d76293 <- load_discovery("GSE76293")
d69528 <- load_discovery("GSE69528")

# --- GSE100159: series matrix + GPL6884 probe annotation ---------------------

read_gse100159 <- function() {
  dd <- file.path(ARDS_ROOT, "00_rawdata/validation/GSE100159")

  lines <- readLines(file.path(dd, "GPL6884.annot"), warn = FALSE)
  b <- which(lines == "!platform_table_begin")
  hdr <- strsplit(lines[b + 1], "\t")[[1]]
  ann <- read.delim(file.path(dd, "GPL6884.annot"), header = FALSE,
                    stringsAsFactors = FALSE, skip = b, quote = "",
                    comment.char = "", check.names = FALSE)
  colnames(ann) <- hdr
  ann <- ann[ann[["Gene symbol"]] != "", c("ID", "Gene symbol")]
  colnames(ann) <- c("probe", "gene")

  sm <- readLines(file.path(dd, "GSE100159_series_matrix.txt"), warn = FALSE)
  tb <- which(sm == "!series_matrix_table_begin")
  te <- which(sm == "!series_matrix_table_end")
  hdr2 <- gsub('"', "", strsplit(sm[tb + 1], "\t")[[1]])
  m4 <- read.delim(textConnection(sm[(tb + 2):(te - 1)]), header = FALSE,
                   sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  colnames(m4) <- hdr2
  rownames(m4) <- m4[[1]]
  e4 <- as.matrix(m4[, -1, drop = FALSE])
  storage.mode(e4) <- "double"

  i <- which(grepl("sample class:", sm, fixed = TRUE))[1]
  cls <- strsplit(gsub("!Sample_characteristics_ch1\t|\"", "", sm[i]), "\t")[[1]]
  cls <- cls[cls != ""][seq_len(ncol(e4))]
  grp <- ifelse(cls %in% c("sample class: melioidosis", "sample class: other infection"),
                "Sepsis", ifelse(cls == "sample class: Control", "Control", "Exclude"))

  keep <- grp %in% c("Sepsis", "Control")
  e4 <- e4[, keep, drop = FALSE]; grp <- grp[keep]

  # collapse probes to genes, keeping the probe with the highest mean signal
  common <- intersect(rownames(e4), ann$probe)
  e4 <- e4[common, , drop = FALSE]
  gmap <- ann$gene[match(common, ann$probe)]
  ord <- order(rowMeans(e4), decreasing = TRUE)
  dfp <- data.frame(probe = common[ord], gene = gmap[ord], stringsAsFactors = FALSE)
  dfp <- dfp[!duplicated(dfp$gene), ]
  e4 <- e4[dfp$probe, , drop = FALSE]
  rownames(e4) <- dfp$gene
  # the series matrix carries linear intensities; log2 puts this cohort on the
  # same scale as the other three
  if (max(e4, na.rm = TRUE) > 50) e4 <- log2(pmax(e4, 0) + 1)
  list(expr = e4, grp = grp)
}

# --- GSE243066: per-sample raw counts, CPM then log2 -----------------------

read_gse243066 <- function() {
  dd <- file.path(ARDS_ROOT, "00_rawdata/validation/GSE243066")
  files <- list.files(dd, pattern = "raw_count.txt.gz$", full.names = TRUE)
  mat <- NULL
  for (f in files) {
    x <- read.delim(gzfile(f), header = TRUE, check.names = FALSE,
                    stringsAsFactors = FALSE)
    colnames(x)[2] <- basename(f)
    mat <- if (is.null(mat)) x else merge(mat, x, by = "Gene_symbol", all = TRUE)
  }
  rownames(mat) <- mat$Gene_symbol
  mat <- as.matrix(mat[, -1, drop = FALSE])
  mat[is.na(mat)] <- 0
  storage.mode(mat) <- "double"
  smp <- colnames(mat)
  grp <- ifelse(grepl("_A[0-9]", smp), "ARDS",
                ifelse(grepl("_H[0-9]", smp), "Control", "Exclude"))
  cpm <- t(t(mat) / colSums(mat) * 1e6)
  lc <- log2(cpm + 1)
  keep <- grp %in% c("ARDS", "Control")
  list(expr = lc[, keep, drop = FALSE], grp = grp[keep])
}

long <- rbind(
  tidy(d76293, "GSE76293"),
  tidy(d69528, "GSE69528"),
  tidy(read_gse100159(), "GSE100159"),
  tidy(read_gse243066(), "GSE243066")
)
long$group <- ifelse(long$group %in% c("Control", "Healthy", "HC"), "Control", "Disease")

write.csv(long, file.path("data", "derived", "hub6_expression_long.csv"),
          row.names = FALSE)
msg("hub6_expression_long.csv: %d rows", nrow(long))
print(with(long, table(cohort, group)))
