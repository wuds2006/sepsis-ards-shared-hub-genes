# ---------------------------------------------------------------------------
# 99_check_figures.R -- pre-submission checks on the rendered figure files
#
# For every figure in Figures/main and Figures/supplementary:
#   1. all four formats present
#   2. SVG declares one font family only (Helvetica)
#   3. SVG font sizes are limited to the two allowed sizes (7 pt, 9 pt)
#   4. SVG carries live <text> nodes rather than outlined glyphs
#   5. PNG is 300 dpi at the intended physical width
#   6. PDF page box matches the PNG geometry
# Usage: Rscript R/99_check_figures.R [Figure1 Figure2 ...]
# ---------------------------------------------------------------------------

source("R/00_setup.R")

args <- commandArgs(trailingOnly = TRUE)
fig_dirs <- c(file.path("..", "Figures", "main"), file.path("..", "Figures", "supplementary"))

files <- unlist(lapply(fig_dirs, function(d) {
  if (!dir.exists(d)) return(character(0))
  sub("\\.pdf$", "", list.files(d, pattern = "\\.pdf$"))
}))
if (length(args)) files <- files[files %in% args]

if (!length(files)) { cat("no figures to check\n"); quit(status = 0) }

prob <- character(0)
note <- function(...) prob <<- c(prob, sprintf(...))

for (f in files) {
  dir <- fig_dirs[file.exists(file.path(fig_dirs, paste0(f, ".svg")))]
  base <- file.path(dir, f)

  for (ext in c(".pdf", ".svg", ".png", ".tiff")) {
    if (!file.exists(paste0(base, ext))) note("%s: missing %s", f, ext)
  }

  # PNG geometry: walk the chunk list for IHDR and pHYs
  png_file <- paste0(base, ".png")
  if (file.exists(png_file)) {
    raw <- readBin(png_file, "raw", 8192)
    be32 <- function(i) sum(as.integer(raw[i:(i + 3)]) * 256^(3:0))
    i <- 9L; px <- NA; ppm <- NA
    while (i + 11 < length(raw)) {
      len <- be32(i); typ <- rawToChar(raw[(i + 4):(i + 7)])
      if (typ == "IHDR") px <- be32(i + 8)
      if (typ == "pHYs") ppm <- be32(i + 8)
      if (typ == "IDAT") break
      i <- i + 12 + len
    }
    if (!is.na(px) && !is.na(ppm)) {
      w_mm <- px / ppm * 1000
      if (abs(w_mm - FIG_W_MM) > 0.6) note("%s: PNG width %.1f mm", f, w_mm)
      if (abs(ppm / 1000 * 25.4 - 300) > 3)
        note("%s: PNG resolution %.0f dpi", f, ppm / 1000 * 25.4)
    } else {
      note("%s: PNG has no pHYs chunk", f)
    }
  }

  # PDF page box, read with Ghostscript (the file carries binary streams)
  pdf_file <- paste0(base, ".pdf")
  if (file.exists(pdf_file)) {
    out <- suppressWarnings(system2("gs", c("-q", "-dNODISPLAY", "-dBATCH", "-c",
             shQuote(sprintf("(%s) (r) file runpdfbegin 1 pdfgetpage /MediaBox get == quit", pdf_file))),
             stdout = TRUE, stderr = FALSE))
    num <- suppressWarnings(as.numeric(gsub("[^0-9. ]", "", out)))
    num <- num[!is.na(num) & num > 0]
    if (length(num) >= 2) {
      w_mm <- max(num[1:2]) / 72 * 25.4
      if (abs(w_mm - FIG_W_MM) > 0.6) note("%s: PDF width %.1f mm", f, w_mm)
    }
  }

  # SVG font family, sizes, live text
  svg_file <- paste0(base, ".svg")
  if (file.exists(svg_file)) {
    s <- paste(readLines(svg_file, warn = FALSE), collapse = "\n")
    fam <- unique(gsub('"', "", sub(".*font-family: *", "", 
                  regmatches(s, gregexpr("font-family: *[^;]*", s))[[1]])))
    if (any(fam != "" & fam != FONT_BODY))
      note("%s: SVG fonts %s", f, paste(setdiff(fam, c("", FONT_BODY)), collapse = ","))
    sz <- as.numeric(gsub("[^0-9.]", "", 
          regmatches(s, gregexpr("font-size: *[0-9.]+px", s))[[1]]))
    bad <- setdiff(unique(sz), c(PT_BODY, PT_TAG))
    if (length(bad)) note("%s: SVG sizes %s", f, paste(bad, collapse = ","))
    n_text <- length(gregexpr("<text", s)[[1]])
    if (n_text < 10) note("%s: only %d text nodes (outlined glyphs?)", f, n_text)
  }
}

if (length(prob)) {
  cat("CHECKS WITH ISSUES\n"); cat(paste0("  - ", prob, collapse = "\n"), "\n")
} else {
  cat(sprintf("ALL CHECKS PASSED (%d figures)\n", length(files)))
}
