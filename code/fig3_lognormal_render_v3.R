#!/usr/bin/env Rscript
# fig3_lognormal_render_v3.R — Fig.3 margin/export + panel-title fix (v3).
# v3 CHANGES vs v2 (reviewer-finaled 2026-10-16; values/curve/band/vlines/points/
# labels/axis ranges/size 9.2x8.2in @300dpi ALL UNCHANGED):
#   (1) panel (c,d) titles: "population median" -> "organism/population-level
#       median" (reviewer: match the actual organism+population scope of the
#       median variant); Fig.3 caption updated to the same labeling (main text).
#   (2) plot.margin(26,10,4,4) -> margin(30,12,8,14): wider left/top/bottom
#       margins so the y-axis title, tick labels, panel tags, and the 2-line
#       x-axis title cannot be clipped (reviewer: left axis/ticks/titles cut).
#   (3) ggsave: units = "in" explicit (export-setting documentation only).
# Prior header (v2): Fig.3 label-fix re-render.
# v2 CHANGES vs v1 (user-approved 2026-10-16; values/curve/band/vlines/axis UNCHANGED):
#   (1) species spelling: Oplegnathus fasciatu -> fasciatus, corrected in the fig3
#       render inputs (02_data/py_out/fig3_sizeclass_v2.csv + fig23_delivery/10_fig3
#       ssdtools_input/srcclass/merged CSVs, 2026-10-16; numbers diff-identical);
#   (2) Skeletonema costatum (panel b) label nudged 0.9 decades LEFT (ggrepel
#       nudge_x) so it stays fully inside the panel (v1 clipped it at the right
#       panel edge); all other species labels identical to v1.
# Prior header (v1): Fig.3 2x2 SSD re-render after user round-3 approval.
# PURPOSE (user-approved, 2026-10-16): Fig.3's curve and CI band must VISUALLY be
# the SAME fit that produced the reported TR1-TR4 thresholds. The previous
# ssdtools model-average curve (v1-v4 lineage) differs from the reported
# single-lognormal TR values by ~0.65-0.70x, which misleads readers into
# thinking curve and TR vlines share one fit.
# THIS RENDER (v1, new lineage):
#   - curve(p) = exp(mu + z_p * sigma), mu = mean(log Conc), sigma = sd(log Conc)
#     ddof = 0  ("engine recipe"; ssd_tomex.py / fig3_recompute_check.py lock)
#   - CI band = parametric bootstrap 95% (nboot = 1000, seed = 20261011) of the
#     SAME single lognormal fit; per-draw curve re-fit, 2.5/97.5 quantiles
#   - band+curve precomputed by py/fig3_lognormal_band_v1.py into
#     03_figures/output_archive/fig3_lognormal_build/ (curve/band/meta CSVs),
#     with assertions curve(0.05) = TR2/TR3, curve(0.10) = TR4 (1e-9 log tol)
#     and band lcl(0.05) = TR1 (q1 panels)
#   - R re-verifies here against the locked TR csv AND against the fit input
#     (mu/sigma from ds$Conc must equal the python-side meta values)
# UNCHANGED vs ssdtools_render_v4.R (user rule: HC5/TR locked, points/theme only
# for the size-class legend v2; species points, pp positions, vlines, tags,
# axis scheme, ggrepel labels, 9.2x8.2in @300dpi all identical):
#   - species points + v2 size-class colours (fig3_sizeclass_v2.csv hard gate)
#   - TR1-TR4 dashed vlines (row-direct from the locked csv, independent check)
#   - native (r-0.5)/n pp positions (convention unchanged; package name dropped
#     from the label since ssdtools no longer produces the curve)
# Output: 03_figures/output/fig3_lognormal_v2.png (final-only rule);
#         audit/log -> 03_figures/output_archive/fig3_lognormal_build/.
# New files only: nothing from the v1-v4 ssdtools lineage is overwritten; the
# retired MA curve remains reproducible via fig3_ssdtools_render_v4.R.

suppressMessages({library(ggplot2); library(patchwork); library(ggrepel)})

INPD <- "data"
BANDD <- "data"
OUTPNG_DIR <- "output"
OUTD <- "output"
dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
LOGC <- character(0)
PP_SRC <- "unseed"
PANELS_META <- list()
P <- function(s = "") { cat(s, "\n"); LOGC <<- c(LOGC, s) }
set.seed(20261011)

sd_pop <- function(x) sqrt(mean((x - mean(x))^2))

# ---- size-data composition legend (INFO-E) ----------------------------------
COL_MULT <- "#009E73"   # Multiple-size data - teal (unified redesign 2026-10-15)
COL_SING <- "#4C78A8"   # Single-size - steel blue
COL_BOTH <- "#CC79A7"   # Both size-data types (Okabe-Ito reddish purple)
COL_UNCL <- "#999999"   # Unclassified         (grey; shown only when present)
CLS_LV <- c("Multiple-size data", "Single-size data",
            "Both size-data types", "Unclassified")
cols_cls <- c("Multiple-size data" = COL_MULT, "Single-size data" = COL_SING,
              "Both size-data types" = COL_BOTH, "Unclassified" = COL_UNCL)

CLS_CSV <- "data/fig3_sizeclass_v2.csv"
CLS_SUM <- "data/fig3_sizeclass_v2_summary.csv"
sc_all <- read.csv(CLS_CSV, fileEncoding = "UTF-8-BOM")
stopifnot(all(c("panel", "Species", "class_v2", "run_type_ours") %in% names(sc_all)))
stopifnot(all(sc_all$class_v2 %in% CLS_LV))
names(sc_all)[names(sc_all) == "class_v2"] <- "class"
sc_sum <- read.csv(CLS_SUM, fileEncoding = "UTF-8-BOM")
stopifnot(ncol(sc_sum) >= 3)
for (pnl in unique(sc_all$panel)) {
  tab <- table(factor(sc_all$class[sc_all$panel == pnl], levels = CLS_LV))
  nsp <- length(unique(sc_all$Species[sc_all$panel == pnl]))
  for (lv in CLS_LV) {
    s <- sc_sum[[3]][sc_sum$panel == pnl & sc_sum$level_1 == lv]
    stopifnot(length(s) <= 1, as.integer(tab[[lv]]) == as.integer(if (length(s)) s else 0))
  }
  sn <- sc_sum[[3]][sc_sum$panel == pnl & sc_sum$level_1 == "n_species"]
  stopifnot(length(sn) == 1, nsp == as.integer(sn))
  rt <- table(sc_all$run_type_ours[sc_all$panel == pnl])
  P(sprintf("size-class v2 gate %s: %s ; n_spp=%d ; run_type_ours: %s", pnl,
            paste(sprintf("%s=%d", CLS_LV, as.integer(tab[CLS_LV])), collapse = " "),
            nsp, paste(sprintf("%s=%d", names(rt), as.integer(rt)), collapse = " ")))
}

# ---- TR values from the locked CSV (row-direct; independent ALL_MATCH check) --
tr <- read.csv(file.path(INPD, "TR_abcd_merged_T12_20_300.csv"),
               fileEncoding = "UTF-8-BOM")
pick <- function(mech, ssd, trn) {
  v <- tr$particles_L[tr$Mechanism == mech & tr$SSD_Collapse == ssd &
                        tr$TR == trn]
  stopifnot(length(v) == 1, is.finite(v))
  v
}
V <- list(
  fd_q1_v1 = pick("food_dilution", "q1_all", "TR1"),
  fd_q1_v2 = pick("food_dilution", "q1_all", "TR2"),
  fd_md_v1 = pick("food_dilution", "median_organism_population", "TR3"),
  fd_md_v2 = pick("food_dilution", "median_organism_population", "TR4"),
  tt_q1_v1 = pick("translocation", "q1_all", "TR1"),
  tt_q1_v2 = pick("translocation", "q1_all", "TR2"),
  tt_md_v1 = pick("translocation", "median_organism_population", "TR3"),
  tt_md_v2 = pick("translocation", "median_organism_population", "TR4"))
chk <- c(0.025481, 1.069965, 11.526777, 69.439538,
         0.691026, 39.142795, 144.866194, 1400.640967)
stopifnot(isTRUE(all.equal(unname(unlist(V)), chk, tolerance = 1e-4)))
P(sprintf("vlines locked: %s", paste(names(V), sprintf("%.5g", unlist(V)),
                                     sep = "=", collapse = " ; ")))

fmt <- function(v) sprintf("%.5g", v)

# power-of-3 log breaks; label "1" for 10^0, expression 10^k otherwise
make_breaks <- function(lo, hi, extra_dec = NA) {
  p3 <- seq(ceiling(lo / 3) * 3, floor(hi / 3) * 3, by = 3)
  if (!length(p3)) p3 <- 0
  if (!is.na(extra_dec)) p3 <- sort(unique(c(p3, extra_dec)))  # v3 axis-extension: user end tick
  brk <- 10^p3
  lab <- lapply(p3, function(k) if (k == 0) "1" else bquote(10^.(k)))
  dec <- 10^(min(floor(lo)):max(ceiling(hi)))
  minor <- unique(as.vector(sapply(dec, function(k) (2:9) * k)))
  minor <- minor[minor >= 10^lo & minor <= 10^hi]
  list(breaks = brk, labels = do.call(expression, lab), minor = minor)
}

panel_fn <- function(mech, coll, panel_id, tagtext, v1, v2, lab1, lab2,
                     n_exp, eff_exp, hc_target, key) {
  # v4 (user spec 2026-10-16): panel ends on 3-multiples, grid-native end ticks
  XEND_DEC <- switch(panel_id, "FD-Q1" = 12, "TT-Q1" = 18,
                     "FD-Median" = 12, "TT-Median" = 15)
  stopifnot(!is.na(XEND_DEC))
  f_in <- file.path(INPD, sprintf("ssdtools_input_marine_%s_%s.csv", mech, coll))
  f_sc <- file.path(INPD, sprintf("srcclass_marine_%s_%s.csv", mech, coll))
  inp <- read.csv(f_in, fileEncoding = "UTF-8-BOM")
  sc  <- read.csv(f_sc, fileEncoding = "UTF-8-BOM")
  stopifnot(nrow(inp) == n_exp, nrow(sc) == n_exp,
            sum(inp$n_effects) == eff_exp,
            isTRUE(all.equal(inp$Conc, sc$Conc, tolerance = 1e-6)),
            setequal(inp$Species, sc$Species))
  sc_p <- sc_all[sc_all$panel == panel_id, c("Species", "class")]
  stopifnot(setequal(sc_p$Species, inp$Species))
  dv <- merge(inp[, c("Species", "Conc")], sc_p, by = "Species")
  stopifnot(!anyNA(dv$class))
  dv$class <- factor(dv$class, levels = CLS_LV)
  n <- nrow(dv)
  ds <- dv[order(dv$Conc), ]
  tb <- table(ds$class)
  P(sprintf("-- %s/%s: n=%d spp, %d effects; conc=%.4g..%.4g; size-class: %s",
            mech, coll, n, sum(ds$n_effects_ <- inp$n_effects[
              match(ds$Species, inp$Species)]),
            min(ds$Conc), max(ds$Conc),
            paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse = " ")))
  # ---- curve/band = the engine single lognormal (precomputed) ----------------
  band_raw <- read.csv(file.path(BANDD,
      sprintf("fig3_lognormal_band_%s_%s.csv", mech, coll)))
  bmeta <- read.csv(file.path(BANDD,
      sprintf("fig3_lognormal_meta_%s_%s.csv", mech, coll)))
  stopifnot(nrow(band_raw) == 99,
            all(is.finite(band_raw$curve)), all(is.finite(band_raw$lcl)),
            all(is.finite(band_raw$ucl)), all(band_raw$lcl < band_raw$ucl))
  stopifnot(all(is.finite(bmeta$mu)), all(is.finite(bmeta$sigma_ddof0)),
            bmeta$n_species == n)
  stopifnot(abs(bmeta$mu - mean(log(ds$Conc))) < 1e-9,
            abs(bmeta$sigma_ddof0 - sd_pop(log(ds$Conc))) < 1e-9)  # same fit input
  ci <- data.frame(proportion = band_raw$pr, est = band_raw$curve,
                   lcl = band_raw$lcl, ucl = band_raw$ucl)
  i5 <- which(abs(ci$proportion - 0.05) < 1e-12); i10 <- which(abs(ci$proportion - 0.10) < 1e-12)
  stopifnot(length(i5) == 1, length(i10) == 1)
  hc_r <- ci$est[i5]
  rat <- hc_r / hc_target
  stopifnot(abs(log(rat)) < 1e-9)   # curve(0.05) IS the reported TR (the fix)
  # cross-panel checks: median panels' curve(0.10) = TR4 = v2
  if (isTRUE(hc_target == v1)) {
    stopifnot(abs(log(ci$est[i10] / v2)) < 1e-9)
    P(sprintf("   engine lognormal: HC5=%.4g curve(0.10)=%.4g=TR4 | band lcl(0.05)=%.4g ucl(0.05)=%.4g (95%% boot of the SAME fit)",
              hc_r, ci$est[i10], ci$lcl[i5], ci$ucl[i5]))
  } else {
    stopifnot(abs(log(ci$lcl[i5] / v1)) < 1e-9)  # q1 panels: band lcl(0.05) = TR1
    P(sprintf("   engine lognormal: HC5=%.4g=TR-locked | band lcl(0.05)=%.4g=TR1 | ucl(0.05)=%.4g | lcl(0.10)=%.4g ucl(0.10)=%.4g (95%% boot of the SAME fit)",
              hc_r, ci$lcl[i5], ci$ucl[i5], ci$lcl[i10], ci$ucl[i10]))
  }
  PP_SRC <<- "native (r-0.5)/n"
  pp  <- (seq_len(n) - 0.5) / n
  obs <- data.frame(Conc = ds$Conc, pp = pp, size_cls = as.character(ds$class),
                    Species = ds$Species)
  lo  <- floor(min(log10(c(obs$Conc, ci$est, v1, v2, ci$lcl))))
  hi  <- ceiling(max(log10(c(obs$Conc, ci$est, v1, v2))))
  hi  <- hi + 1.35    # room for right-side species labels
  hi  <- max(hi, XEND_DEC + 0.75)  # v4 (user spec): range; end tick comes from standard 3-decade grid
  bk  <- make_breaks(lo, hi)
  band <- data.frame(x = c(ci$lcl, rev(ci$ucl)),
                     y = c(ci$proportion, rev(ci$proportion)))
  gp <- ggplot()
  gp <- gp + geom_polygon(data = band, aes(x = x, y = y, fill = "95% CI"),
                          alpha = 0.45, color = NA)
  gp <- gp + geom_line(data = ci, aes(x = est, y = proportion),
                       color = "grey50", linewidth = 0.9)
  gp <- gp + geom_point(data = obs,
                        aes(x = Conc, y = pp, color = size_cls, fill = size_cls),
                        size = 2.7, stroke = 0.7)
  gp <- gp + geom_text_repel(
    data = obs, aes(x = Conc, y = pp, label = Species),
    size = 2.1, fontface = "italic", color = "grey15",
    min.segment.length = 0.25, max.overlaps = Inf,
    box.padding = 0.28, point.padding = 0.22, force = 1.2,
    nudge_x = ifelse(obs$Species == "Skeletonema costatum", -0.9,
              ifelse(obs$Species == "Parvocalanus crassirostris", -0.5,
              ifelse(obs$Species == "Karenia mikimotoi", -1.7, 0))),  # v3 label-clip fix C,D (decades)
    xlim = c(NA, 10^hi), segment.color = "grey62", segment.size = 0.3)
  gp <- gp + geom_vline(xintercept = v1, linetype = "dashed",
                        color = "grey25", linewidth = 0.65)
  gp <- gp + geom_vline(xintercept = v2, linetype = "dashed",
                        color = "grey25", linewidth = 0.65)
  gp <- gp + scale_x_log10(breaks = bk$breaks, labels = bk$labels,
                           minor_breaks = bk$minor,
                           limits = c(10^lo, 10^hi),
                           expand = expansion(mult = 0.015))
  gp <- gp + scale_y_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.2),
                                minor_breaks = seq(0.1, 0.9, 0.1),
                                labels = function(v) paste0(v * 100, "%"),
                                expand = expansion(mult = 0.01))
  pres <- CLS_LV[CLS_LV %in% as.character(unique(ds$class))]
  stopifnot(length(pres) >= 1)
  brk <- c("95% CI", pres)
  P(sprintf("   legend keys: %s", paste(brk, collapse = ", ")))
  gp <- gp + scale_fill_manual(values = c("95% CI" = "grey82", cols_cls[pres]),
                               breaks = brk, labels = brk, name = NULL)
  gp <- gp + scale_colour_manual(values = cols_cls[pres], name = NULL)
  gp <- gp + guides(fill = guide_legend(order = 1,
                     override.aes = list(colour = c("grey85",
                                                    unname(cols_cls[pres])))),
                    colour = "none")
  gp <- gp + labs(x = paste0("Concentration (particles/L)\n",
                             lab1, " = ", fmt(v1), " ; ", lab2, " = ", fmt(v2)),
                  y = "Percent of species affected", tag = tagtext)
  gp <- gp + theme_classic(base_size = 11) +
    theme(panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
          panel.grid.minor.y = element_line(colour = "grey94", linewidth = 0.2),
          panel.grid.minor.x = element_line(colour = "grey95", linewidth = 0.15))
  gp <- gp + theme(
    plot.tag = element_text(size = 11.5, face = "plain", hjust = 0, vjust = 0),
    plot.tag.position = c(0, 1.02),
    legend.position = c(0.02, 0.97), legend.justification = c(0, 1),
    legend.title = element_blank(),
    legend.background = element_rect(fill = "white", colour = "grey65",
                                     linewidth = 0.3),
    legend.key.size = unit(0.3, "cm"), legend.text = element_text(size = 8.5),
    axis.text.x = element_text(size = 9.5), axis.text.y = element_text(size = 9.5),
    plot.margin = margin(30, 12, 8, 14))
  P(sprintf("   pp=%s ; x-range=10^%.1f..10^%.1f ; breaks=%s", PP_SRC, lo, hi,
            paste(sprintf("%.0e", bk$breaks), collapse = ",")))
  PANELS_META[[key]] <<- sprintf(
    '{"n_species":%d,"n_effects":%d,"curve_hc5":%s,"tr_locked":%s,"ratio":%.6f,"v1":%s,"v2":%s,"mu":%.6f,"sigma_ddof0":%.6f,"size_class_counts":{%s}}',
    n, sum(inp$n_effects), format(hc_r, digits = 6, trim = TRUE),
    format(hc_target, digits = 6, trim = TRUE), rat,
    format(v1, digits = 6, trim = TRUE), format(v2, digits = 6, trim = TRUE),
    mean(log(ds$Conc)), sd_pop(log(ds$Conc)),
    paste(sprintf('"%s":%d', names(tb), as.integer(tb)), collapse = ","))
  list(plot = gp, hc5_curve = hc_r)
}

pa <- panel_fn("food_dilution", "q1_all", "FD-Q1",
               "(a) Food dilution - 1st quartile", V$fd_q1_v1, V$fd_q1_v2,
               "TR1 (HC5 95% CI lower)", "TR2 (HC5)", 20, 367, V$fd_q1_v2, "fd_q1")
pb <- panel_fn("translocation", "q1_all", "TT-Q1",
               "(b) Tissue translocation - 1st quartile", V$tt_q1_v1, V$tt_q1_v2,
               "TR1 (HC5 95% CI lower)", "TR2 (HC5)", 26, 288, V$tt_q1_v2, "tt_q1")
pc <- panel_fn("food_dilution", "median_organism_population", "FD-Median",
               "(c) Food dilution - organism/population-level median", V$fd_md_v1, V$fd_md_v2,
               "TR3 (HC5)", "TR4 (HC10)", 19, 184, V$fd_md_v1, "fd_md")
pd <- panel_fn("translocation", "median_organism_population", "TT-Median",
               "(d) Tissue translocation - organism/population-level median", V$tt_md_v1, V$tt_md_v2,
               "TR3 (HC5)", "TR4 (HC10)", 23, 141, V$tt_md_v1, "tt_md")
fig <- (pa$plot | pb$plot) / (pc$plot | pd$plot)  # true 2x2 (see fig2 v9 note)

png_path <- file.path(OUTPNG_DIR, "fig3_lognormal_v3.png")
set.seed(20261011)   # ggrepel label layout reproducibility
ggsave(png_path, fig, width = 9.2, height = 8.2, units = "in", dpi = 300, bg = "white")

vjson <- paste0('"', names(V), '":', sprintf("%.6g", unlist(V)),
                collapse = ",")
audit <- paste0(
  '{"figure":"fig3_lognormal_v3",',
  '"engine":"parametric lognormal (engine recipe locked 2026-10-13: mu=mean(log), sigma=sd ddof=0; curve(p)=exp(mu+z_p*sigma))",',
  '"ci":"parametric bootstrap 95% of the SAME single lognormal fit, nboot=1000, seed=20261011, 2.5/97.5 quantiles (default linear interp)",',
  '"data":"02_data/fig23_delivery/10_fig3 confirmed delivery (9 input files md5-SAME vs output_archive/fig3_merge/ssdtools_inputs, 2026-10-14): 4 ssdtools_input CSV; TR_abcd_merged_T12_20_300.csv row-direct; ',
  'band+curve precomputed by 03_figures/code/py/fig3_lognormal_band_v1.py into 03_figures/output_archive/fig3_lognormal_build/ (band/meta CSVs per panel); ',
  'size classes from 02_data/py_out/fig3_sizeclass_v2.csv (fig3_panel_sizeclass_v2.py study-series re-adjudication v2, gate: per-panel counts locked to fig3_sizeclass_v2_summary.csv)",',
  '"spec":"2x2 FD/TT x Q1/median; window 20-300 um; particles/L; alpha=2.07 engine descriptors",',
  '"pp_convention":"', PP_SRC, '",',
  '"vline_values":{', vjson, '},',
  '"panels":{',
    '"a_fd_q1":', PANELS_META[["fd_q1"]], ',',
    '"b_tt_q1":', PANELS_META[["tt_q1"]], ',',
    '"c_fd_md":', PANELS_META[["fd_md"]], ',',
    '"d_tt_md":', PANELS_META[["tt_md"]],
  '},',
  '"legend_scheme":"particle-size data composition - approved study-series v2 adjudication (2026-10-14): same-species comparable-test >=2 independent sizes -> Multiple-size (dates irrelevant); Within/Across-run recorded in run_type_ours column of fig3_sizeclass_v2.csv; Unclassified key omitted (0 in all panels)",',
  '"legend_colors":{"Multiple-size data":"', COL_MULT, '","Single-size data":"', COL_SING,
  '","Both size-data types":"', COL_BOTH, '","Unclassified":"', COL_UNCL, '"},',
  '"labels":"all species, ggrepel italic, seed=20261011; v2: Skeletonema costatum (panel b) nudge_x=-0.9 decades (label-clip fix, user-approved 2026-10-16); Oplegnathus fasciatus spelling fixed in the 10 render-input CSVs (numbers diff-identical)",',
  '"png":"fig3_lognormal_v3.png",',
  '"note":"user round-3 approved re-render (2026-10-16): Fig.3 curve+band now = the SINGLE parametric lognormal fit underlying the reported TR1-TR4 (engine recipe, seed 20261011); R asserts curve(0.05)=TR2/TR3, curve(0.10)=TR4, band lcl(0.05)=TR1 (q1) at <1e-9 log tol, and mu/sigma re-derived from the panel input CSV equal the python-side meta values; the retired ssdtools model-average curve lineage (fig3_ssdtools_render_v4.R, v4 audit: MA HC5 = 0.654-0.701x of TR-locked values) is kept reproducible in output_archive; points, vlines, sizes, tags, theme unchanged; v2 label-only fix (user-approved 2026-10-16): (i) Oplegnathus fasciatus spelling corrected in fig3 render inputs (10 CSVs, numbers diff-identical), (ii) Skeletonema costatum panel-(b) label nudged 0.9 decades left via ggrepel nudge_x; all other labels, points, curve, band, vlines, axis ranges identical to v1; v1 png kept in output_archive/fig3_lognormal_build/"}')
writeLines(audit, file.path(OUTD, "fig3_lognormal_v3_audit.json"))
writeLines(LOGC, file.path(OUTD, "fig3_lognormal_v3_log.txt"))
P(paste0("[SAVED] ", png_path))