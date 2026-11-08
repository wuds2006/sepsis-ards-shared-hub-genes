# ---------------------------------------------------------------------------
# 00_graphical_abstract.R -- graphical abstract
#
# Style follows the flow-chart reference supplied by the author: a solid deep
# green title bar, dashed group frames carrying a solid pill header that sits
# on the frame's top edge, white rounded cards inside each frame, a colour icon
# at the left of every card, and heavy solid arrows down the columns and across
# between them.
#
# Drawn in grid on a millimetre canvas rather than in ggplot: round corners,
# dashed frames, pill headers straddling a frame edge and 1 mm arrow shafts
# cannot be expressed faithfully in a ggplot scale. Text still goes through
# Helvetica at 6.0-8.6 pt so that the file matches the other ten figures.
#
# Content is carried over unchanged from the previous version. Every number is
# lifted from Results and keeps its qualifier in the same card: 0.977 is
# labelled as post-ComBat, and MR is labelled as null with a reverse-only
# powered direction. Nothing is rounded up and nothing is added.
#
# Output : Figures/graphical_abstract/Graphical_Abstract.{pdf,svg,png,tiff}
# ---------------------------------------------------------------------------

source("R/00_setup.R")

GA_DIR <- file.path("..", "Figures", "graphical_abstract")
dir.create(GA_DIR, recursive = TRUE, showWarnings = FALSE)

W <- 180; H <- 110

# --- style -----------------------------------------------------------------

# Every colour below comes from the manuscript palette defined in 00_setup.R,
# or is a tint (mixw) / shade (dk) of one. Nothing is invented, so the
# graphical abstract cannot drift away from the ten main figures. The three
# column accents are the three palette colours that already carry meaning in
# the manuscript: brick for disease, steel for comparison, plum for the
# fourth category.

TH <- pal

mixw <- function(col, p) {          # p parts colour, the rest white
  r <- grDevices::col2rgb(col) / 255
  grDevices::rgb(r[1] * p + (1 - p), r[2] * p + (1 - p), r[3] * p + (1 - p))
}
dk <- function(col, p) {            # p < 1 darkens
  r <- grDevices::col2rgb(col) / 255
  grDevices::rgb(r[1] * p, r[2] * p, r[3] * p)
}

COL <- list(
  bar   = dk(TH[["steel"]], 0.76),  # title bar: steel, one step deeper
  arrow = TH[["steel"]],            # flow arrows, same hue as the pipeline
  card  = "#FFFFFF",
  ink   = TH[["ink"]],
  sub   = mixw(TH[["ink"]], 0.78),  # secondary text
  white = "#FFFFFF"
)

# --- grid primitives, all in millimetres -----------------------------------

seg <- function(x0, y0, x1, y1, col, lwd = 1, lend = "round", lty = "solid") {
  grid.segments(unit(x0, "mm"), unit(y0, "mm"), unit(x1, "mm"), unit(y1, "mm"),
                gp = gpar(col = col, lwd = lwd, lineend = lend, lty = lty))
}

rr <- function(x0, y0, x1, y1, r, fill = NA, col = NA, lwd = 1, lty = "solid") {
  grid.roundrect(unit((x0 + x1) / 2, "mm"), unit((y0 + y1) / 2, "mm"),
                 width = unit(x1 - x0, "mm"), height = unit(y1 - y0, "mm"),
                 r = unit(r, "mm"),
                 gp = gpar(fill = fill, col = col, lwd = lwd, lty = lty))
}

ell <- function(cx, cy, a, b, fill = NA, col = NA, lwd = 1) {
  th <- seq(0, 2 * pi, length.out = 60)
  grid.polygon(unit(cx + a * cos(th), "mm"), unit(cy + b * sin(th), "mm"),
               gp = gpar(fill = fill, col = col, lwd = lwd))
}

blob <- function(x, y, r, fill = NA, col = NA, lwd = 1) {
  grid.circle(unit(x, "mm"), unit(y, "mm"), unit(r, "mm"),
              gp = gpar(fill = fill, col = col, lwd = lwd))
}

# --- rich text: gene symbols italic, per AMA -------------------------------

GA_GENES <- c("CDK5RAP2", "CD177", "GRB10", "UPP1", "EXOSC4", "VNN1",
              "CEBPB", "MIF")

tokenize <- function(s) {
  pat <- paste0("\\b(?:", paste(GA_GENES, collapse = "|"), ")\\b")
  m <- gregexpr(pat, s, perl = TRUE)[[1]]
  if (length(m) == 1L && m[1] == -1L) return(list(list(t = s, it = FALSE)))
  ml <- attr(m, "match.length")
  out <- list(); pos <- 1L
  for (k in seq_along(m)) {
    if (m[k] > pos)
      out[[length(out) + 1]] <- list(t = substr(s, pos, m[k] - 1L), it = FALSE)
    out[[length(out) + 1]] <- list(t = substr(s, m[k], m[k] + ml[k] - 1L), it = TRUE)
    pos <- m[k] + ml[k]
  }
  if (pos <= nchar(s))
    out[[length(out) + 1]] <- list(t = substr(s, pos, nchar(s)), it = FALSE)
  out
}

tw <- function(s, size, face = "plain") {
  convertWidth(grobWidth(textGrob(s, gp = gpar(fontsize = size,
                                               fontfamily = FONT_BODY,
                                               fontface = face))),
               "mm", valueOnly = TRUE)
}

.miss <- character(0)

rich <- function(x, y, s, size, face = "plain", col = COL$ink,
                 just = "left", maxw = NULL, tag = "") {
  tk <- tokenize(s)
  ws <- vapply(tk, function(z) tw(z$t, size, if (z$it) "italic" else face),
               numeric(1))
  tot <- sum(ws)
  if (!is.null(maxw) && tot > maxw)
    .miss <<- c(.miss, sprintf("%-22s %5.1f > %5.1f mm  | %s", tag, tot, maxw, s))
  x0 <- switch(just, centre = x - tot / 2, right = x - tot, x)
  cx <- x0
  for (k in seq_along(tk)) {
    grid.text(tk[[k]]$t,
              x = unit(cx, "mm"), y = unit(y, "mm"),
              just = c("left", "centre"),
              gp = gpar(fontsize = size, fontfamily = FONT_BODY,
                        fontface = if (tk[[k]]$it) "italic" else face,
                        col = col))
    cx <- cx + ws[k]
  }
  invisible(tot)
}

# --- icons -----------------------------------------------------------------
# Each icon lives in a 6 mm box centred on (cx, cy). c1 carries the mass, c2
# the accent, so that a card is readable at a glance without a legend.

icon <- function(type, cx, cy, s, c1, c2) {
  h <- s / 2
  switch(type,
    db = {
      a <- h * 0.74; b <- h * 0.24
      rr(cx - a, cy - h * 0.60, cx + a, cy + h * 0.62, 0.05,
         fill = adjustcolor(c1, 0.62), col = NA)
      for (yy in c(cy - 1.15, cy + 1.25)) {
        th <- seq(pi, 2 * pi, length.out = 40)
        grid.lines(unit(cx + a * cos(th), "mm"),
                   unit(yy + b * sin(th), "mm"),
                   gp = gpar(col = adjustcolor(c1, 0.95), lwd = 0.8))
      }
      th <- seq(pi, 2 * pi, length.out = 40)
      grid.lines(unit(cx + a * cos(th), "mm"),
                 unit(cy - h * 0.60 + b * sin(th), "mm"),
                 gp = gpar(col = adjustcolor(c1, 0.95), lwd = 0.8))
      ell(cx, cy + h * 0.62, a, b, fill = adjustcolor(c1, 0.95), col = c1, lwd = 0.9)
    },
    funnel = {
      grid.polygon(unit(cx + c(-h * 0.82, h * 0.82, h * 0.26, -h * 0.26), "mm"),
                   unit(cy + c(h * 0.8, h * 0.8, -h * 0.2, -h * 0.2), "mm"),
                   gp = gpar(fill = adjustcolor(c1, 0.45), col = c1, lwd = 1))
      seg(cx, cy - h * 0.2, cx, cy - h * 0.9, c1, 1.3)
    },
    bars = {
      hs <- c(2.0, 3.3, 4.6)
      for (i in 1:3) {
        x <- cx - h * 0.72 + (i - 1) * h * 0.56
        rr(x, cy - h * 0.85, x + 0.95, cy - h * 0.85 + hs[i], 0.22,
           fill = if (i == 3) c1 else adjustcolor(c1, 0.66), col = NA)
      }
    },
    venn = {
      blob(cx - 1.05, cy, h * 0.72, fill = adjustcolor(c1, 0.55), col = c1, lwd = 0.9)
      blob(cx + 1.05, cy, h * 0.72, fill = adjustcolor(c2, 0.45), col = c2, lwd = 0.9)
    },
    gear = {
      for (a in seq(0, 2 * pi, length.out = 9)[-9]) {
        x <- cx + 2.15 * cos(a); y <- cy + 2.15 * sin(a)
        rr(x - 0.5, y - 0.5, x + 0.5, y + 0.5, 0.28, fill = c1, col = NA)
      }
      blob(cx, cy, 1.85, fill = NA, col = c1, lwd = 1.7)
      blob(cx, cy, 0.72, fill = c2, col = NA)
    },
    roc = {
      seg(cx - 2.2, cy - 2.1, cx - 2.2, cy + 2.3, c2, 1)
      seg(cx - 2.2, cy - 2.1, cx + 2.4, cy - 2.1, c2, 1)
      grid.lines(unit(cx + c(-2.2, -1.0, 0.4, 2.35), "mm"),
                 unit(cy + c(-2.1, -0.5, 0.8, 2.2), "mm"),
                 gp = gpar(col = c1, lwd = 1.7, lineend = "round", linejoin = "round"))
    },
    cell = {
      blob(cx, cy, h * 0.86, fill = adjustcolor(c1, 0.18), col = c1, lwd = 1.1)
      blob(cx + 0.35, cy + 0.15, h * 0.36, fill = c2, col = NA)
      blob(cx - 1.55, cy - 1.15, 0.3, fill = adjustcolor(c1, 0.9), col = NA)
      blob(cx + 1.45, cy + 1.25, 0.26, fill = adjustcolor(c1, 0.9), col = NA)
    },
    net = {
      nx <- cx + c(-1.85, 0.45, 2.05, -0.55)
      ny <- cy + c(1.65, 2.15, 0.85, -1.95)
      e <- rbind(c(1, 2), c(2, 3), c(1, 4), c(2, 4), c(3, 4))
      for (k in seq_len(nrow(e)))
        seg(nx[e[k, 1]], ny[e[k, 1]], nx[e[k, 2]], ny[e[k, 2]],
            adjustcolor(c1, 0.75), 0.9)
      for (k in 1:4) blob(nx[k], ny[k], 0.58, fill = c2, col = c1, lwd = 0.8)
    },
    target = {
      blob(cx, cy, h * 0.92, fill = NA, col = c1, lwd = 1.2)
      blob(cx, cy, h * 0.56, fill = NA, col = c1, lwd = 1.2)
      blob(cx, cy, h * 0.22, fill = c2, col = NA)
    },
    flask = {
      grid.polygon(unit(cx + c(-0.78, 0.78, h * 0.76, -h * 0.76), "mm"),
                   unit(cy + c(h * 0.5, h * 0.5, -h * 0.86, -h * 0.86), "mm"),
                   gp = gpar(fill = adjustcolor(c1, 0.35), col = c1, lwd = 1))
      rr(cx - 0.78, cy + h * 0.46, cx + 0.78, cy + h * 0.94, 0.2, fill = c1, col = NA)
      seg(cx - 1.55, cy - h * 0.32, cx + 1.55, cy - h * 0.32, c2, 1.1)
    },
    ban = {
      blob(cx, cy, h * 0.86, fill = adjustcolor(c1, 0.14), col = c1, lwd = 1.4)
      seg(cx - 1.65, cy - 1.65, cx + 1.65, cy + 1.65, c1, 1.4)
    },
    bound = {
      seg(cx, cy - h * 0.86, cx, cy + h * 0.86, c1, 1.5)
      # grid.lines needs two points; a single point draws nothing and the
      # arrow head silently disappears
      for (sgn in c(-1, 1)) {
        grid.lines(unit(c(cx + sgn * h * 0.88, cx + sgn * 0.34), "mm"),
                   unit(c(cy, cy), "mm"),
                   arrow = arrow(length = unit(1.4, "mm"), type = "closed",
                                 angle = 26),
                   gp = gpar(col = c2, lwd = 1.6, lineend = "round"))
      }
    }
  )
}

# --- layout ----------------------------------------------------------------

# three columns, four cards each; the column gap has to be wide enough for a
# heavy arrow to read as an arrow and not as a tick
COL_X <- list(c(5.50, 55.17), c(65.17, 114.83), c(124.83, 174.50))
BOX_Y <- c(4.5, 93.0)          # frame bottom, frame top (= pill centre)
CARD_H <- 15.0
CARD_GAP <- 5.0
PAD <- 4.0                     # card inset from the frame side

COLS <- list(
  list(tag = "Discovery", pill = TH[["brick"]], ic2 = TH[["ochre"]], cards = list(
    list(ic = "db",     head = "Two discovery cohorts",
         sub = c("ARDS 12 vs 12", "sepsis 83 vs 28")),
    list(ic = "funnel", head = "Re-curation",
         sub = c("off-target samples removed")),
    list(ic = "bars",   head = "Differential expression",
         sub = c("557 ARDS \u00b7 1,458 sepsis", "13 + 11 modules")),
    list(ic = "venn",   head = "Intersection",
         sub = c("179 shared candidates"))
  )),
  list(tag = "Evidence", pill = TH[["steel"]], ic2 = TH[["sage"]], cards = list(
    list(ic = "gear",   head = "11 classifiers",
         sub = c("LASSO \u00b7 RF \u00b7 SHAP", "6 core genes retained")),
    list(ic = "roc",    head = "External validation",
         sub = c("ARDS AUC 0.984", "sepsis 0.708, then 0.977",
                 "after ComBat adjustment")),
    list(ic = "cell",   head = "Blood and single cell",
         sub = c("monocytes up, CD8+ T down", "MIF top; complement up",
                 "CEBPB upstream of UPP1")),
    list(ic = "net",    head = "Causal tests",
         sub = c("MR null; reverse only", "no colocalization",
                 "LDSC rg = -0.05"))
  )),
  list(tag = "Translation", pill = TH[["plum"]], ic2 = TH[["slate"]], cards = list(
    list(ic = "target", head = "Two tractable targets",
         sub = c("VNN1 and UPP1 only")),
    list(ic = "flask",  head = "Docking",
         sub = c("luteolin on VNN1", "-10.01 kcal/mol")),
    list(ic = "ban",    head = "No approved drug",
         sub = c("none of the six")),
    list(ic = "bound",  head = "What this supports",
         sub = c("shared monocyte state", "not a causal pathway"), dark = TRUE)
  ))
)

# --- draw ------------------------------------------------------------------

ga_draw <- function() {
  grid.newpage()
  pushViewport(viewport(width = unit(W, "mm"), height = unit(H, "mm")))
  on.exit(popViewport(), add = TRUE)

  # title bar ---------------------------------------------------------------
  rr(5.5, 99.75, 174.5, 107.25, 3.75, fill = COL$bar, col = NA)
  rich(90, 103.5, "ARDS and sepsis share a monocyte-centred host-response state",
       8.6, "bold", COL$white, "centre", maxw = 163, tag = "title")

  # columns -----------------------------------------------------------------
  for (ci in seq_along(COLS)) {
    co <- COLS[[ci]]
    x0 <- COL_X[[ci]][1]; x1 <- COL_X[[ci]][2]
    accent     <- co$pill
    frame_fill <- mixw(accent, 0.09)
    frame_line <- mixw(accent, 0.44)
    card_edge  <- mixw(accent, 0.26)
    ic1        <- accent
    ic2        <- co$ic2
    dark_fill  <- dk(accent, 0.72)

    # dashed frame with the pill header straddling its top edge
    rr(x0, BOX_Y[1], x1, BOX_Y[2], 4.2,
       fill = frame_fill, col = frame_line, lwd = 1.0, lty = "33")

    pw <- tw(co$tag, 7.2, "bold") + 9.0
    rr(x0 + 5.0, BOX_Y[2] - 3.75, x0 + 5.0 + pw, BOX_Y[2] + 3.75, 3.75,
       fill = accent, col = NA)
    rich(x0 + 5.0 + pw / 2, BOX_Y[2], co$tag, 7.2, "bold", COL$white, "centre",
         maxw = pw - 3, tag = paste0("pill", ci))

    # cards
    cw <- (x1 - PAD) - (x0 + PAD)
    top <- BOX_Y[2] - 8.5
    for (ki in seq_along(co$cards)) {
      cd <- co$cards[[ki]]
      ytop <- top - (ki - 1) * (CARD_H + CARD_GAP)
      ybot <- ytop - CARD_H
      dark <- isTRUE(cd$dark)

      rr(x0 + PAD, ybot, x1 - PAD, ytop, 2.9,
         fill = if (dark) dark_fill else COL$card, col = if (dark) NA else card_edge,
         lwd = 0.9)

      cy <- (ytop + ybot) / 2
      icon(cd$ic, x0 + PAD + 5.5, cy, 6.5,
           if (dark) COL$white else ic1,
           if (dark) mixw(accent, 0.5) else ic2)

      tx  <- x0 + PAD + 10.8
      mw  <- (x1 - PAD) - 2.0 - tx
      tcol <- if (dark) COL$white else COL$ink
      scol <- if (dark) adjustcolor(COL$white, 0.86) else COL$sub

      k <- length(cd$sub)
      bh <- 2.42 + 2.50 + (k - 1) * 2.30
      yh <- cy + bh / 2 - 1.21
      rich(tx, yh, cd$head, 6.8, "bold", tcol, "left", maxw = mw,
           tag = paste0("c", ci, ki, "h"))
      for (si in seq_len(k))
        rich(tx, yh - 2.50 - (si - 1) * 2.30, cd$sub[si], 6.0, "plain", scol,
             "left", maxw = mw, tag = paste0("c", ci, ki, "s", si))

      # arrow down to the next card in the same column
      if (ki < length(co$cards)) {
        grid.lines(unit(x0 + PAD + 5.5, "mm"),
                   unit(c(ybot - 0.6, ybot - CARD_GAP + 0.8), "mm"),
                   arrow = arrow(length = unit(1.9, "mm"), type = "closed",
                                 angle = 28),
                   gp = gpar(col = COL$arrow, lwd = 2.6, lineend = "round"))
      }
    }

    # arrow across to the next column
    if (ci < length(COLS)) {
      gx0 <- x1 + 0.5; gx1 <- COL_X[[ci + 1]][1] - 1.1
      grid.lines(unit(c(gx0, gx1), "mm"), unit(c(47.0, 47.0), "mm"),
                 arrow = arrow(length = unit(2.4, "mm"), type = "closed",
                               angle = 30),
                 gp = gpar(col = COL$arrow, lwd = 3.4, lineend = "round"))
    }
  }
}

# --- devices ---------------------------------------------------------------
# save_fig() is not reused here: it draws a ggplot, and grid needs par(mar = 0)
# so that the 180 x 110 mm viewport coincides with the page.

ga_save <- function(name, dir, width_mm = 180, height_mm = 110) {
  w <- mm(width_mm); h <- mm(height_mm)
  prep <- function() par(mar = rep(0, 4), oma = rep(0, 4), xpd = NA)

  grDevices::pdf(file.path(dir, paste0(name, ".pdf")), width = w, height = h,
                 family = FONT_BODY, useDingbats = FALSE, onefile = FALSE,
                 bg = "white", colormodel = "srgb")
  prep(); ga_draw(); grDevices::dev.off()

  svglite::svglite(file.path(dir, paste0(name, ".svg")), width = w, height = h,
                   bg = "white")
  prep(); ga_draw(); grDevices::dev.off()

  png_file <- file.path(dir, paste0(name, ".png"))
  ragg::agg_png(png_file, width = w, height = h, units = "in", res = 300,
                background = "white")
  prep(); ga_draw(); grDevices::dev.off()
  img <- magick::image_read(png_file)
  magick::image_write(img, png_file, format = "png", density = "300x300")

  magick::image_write(img, file.path(dir, paste0(name, ".tiff")),
                      format = "tiff", compression = "lzw", density = "300x300")
  invisible(png_file)
}

ga_save("Graphical_Abstract", GA_DIR)

if (length(.miss)) {
  cat("\n--- text wider than its box ---\n")
  cat(paste(.miss, collapse = "\n"), "\n")
  stop(sprintf("%d text run(s) overflow", length(.miss)))
} else {
  msg("Graphical_Abstract written: %s (no overflow)", GA_DIR)
}
