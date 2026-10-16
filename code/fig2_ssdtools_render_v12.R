#!/usr/bin/env Rscript
# fig2_ssdtools_render_v12.R — fig2 (ssdtools 2.6.0) 4-panel SSD: GM & MIN tracks
# = v11 CLASS-REBUILD (2026-10-14, INFO v2 승인):
#   분류 = fig2_sizeclass_v2.csv (v2: 판정 단위=비교 가능한 연구계열(study series),
#   날짜=run 기록 보조변수 — Within-run/Across-run은 보충자료 기록).
#     급성: 적격 68행 (Sheet,No)->Domestic_Extended_Long 조인, miss=0
#           => Multiple 8 / Single 5 (보정② near-dup 112/112.1 적격행 미출현)
#     만성: ours 8종(candidates n_size 8/8 대조) M6/S2 + ToMEx_fix1000 3종 S3
#           + Jung2021-marine 4종 U4(원문연구 source 특정, 크기=범위 표기)
#           => Multiple 6 / Single 5 / Both 0 / Unclassified 4
#   범례 4명칭은 동일 표기(의미만 v2 개정), Both=0이어서 패널 범례에서 자동 생략.
#   유일 변경 = 클래스 매핑/범례; 값·트랙(GM/MIN)·HC5 락 검사·부트스트랩·식 종
#   구성 전부 v9/v10 그대로 불변 (집계 규칙 미변경 -> HC5·TR 유지).
# Data : $BASE/ssd_pooled_reps_v9.csv (28 rows: Acute 13 + Chronic 15) ;
#        HC5  = $BASE/ssd_pooled_hc5_v9.csv (lock, rel<=1e-6 검증 완료)
# Curve = lnorm+llogis model average; ribbon = parametric bootstrap 95% CI
# (nboot=1000, seed 20260929); dashed vline = python-locked HC5.
# Style = approved fig2 v9/v10; PNG untouched (new file only: fig2_ssdtools_v12.png).

suppressMessages({library(ssdtools); library(ggplot2); library(patchwork)
library(ggrepel)})

BASE <- "data"   # repo-relative: run Rscript from the repo root (see README)
CSVP <- "data/fig2_sizeclass_v2.csv"
OUTD <- "output"
LOGC <- character(0)
PP_SRC <- "unset"
P <- function(s = "") { cat(s, "\n"); LOGC <<- c(LOGC, s) }
set.seed(20260929)

# 4-class palette (INFO-E 확정색) — order fixed: Multiple/Single/Both/Unclassified
# (이 표기는 범례 라벨용; 분류 의미는 v2 연구계열 기준으로 개정됨)
ORDER_ALL <- c("Multiple-size data", "Single-size data",
               "Both size-data types", "Unclassified")
COL_MULTI   <- "#009E73"   # Multiple-size data - teal (통일 리디자인 2026-10-15: warm 제거)
COL_SINGLE  <- "#4C78A8"   # Single-size data - steel blue (TT #0072B2와 구분)
COL_BOTH    <- "#CC79A7"
COL_UNCLASS <- "#999999"
CLASS_COLS  <- c("Multiple-size data" = COL_MULTI, "Single-size data" = COL_SINGLE,
                 "Both size-data types" = COL_BOTH, "Unclassified" = COL_UNCLASS)
# v2 CSV 짧은 클래스명 -> 범례 라벨 매핑
LABEL_MAP <- c("Multiple" = "Multiple-size data", "Single" = "Single-size data",
               "Both" = "Both size-data types", "Unclassified" = "Unclassified")

d <- read.csv(file.path(BASE, "ssd_pooled_reps_v9.csv"))
h <- read.csv(file.path(BASE, "ssd_pooled_hc5_v9.csv"))
P(sprintf("ssdtools %s | R %s", as.character(packageVersion("ssdtools")), getRversion()))
P(sprintf("data rows=%d ; hc5 rows=%d", nrow(d), nrow(h)))
stopifnot(sum(d$Tier == "Acute") == 13, sum(d$Tier == "Chronic") == 15, nrow(d) == 28)

# ---- 4-class particle-size data composition (v2 study-series, legend/colour ONLY) ----
cls_raw <- read.csv(CSVP, fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
P(sprintf("sizeclass csv rows=%d ; species uniq=%d ; cols=%s", nrow(cls_raw),
          length(unique(paste(cls_raw$Tier, cls_raw$Species))),
          paste(names(cls_raw), collapse = ",")))
stopifnot(all(c("Tier", "Species", "size_class_v2") %in% names(cls_raw)))
mrg <- merge(d, cls_raw[, c("Tier", "Species", "size_class_v2")],
             by = c("Tier", "Species"), all.x = TRUE)
stopifnot(nrow(mrg) == nrow(d), !anyNA(mrg$size_class_v2))
stopifnot(all(mrg$size_class_v2 %in% names(LABEL_MAP)))
d <- mrg   # merged row order (panel_fn re-sorts by Conc, HC5 join by name)
d$size_cls <- factor(LABEL_MAP[d$size_class_v2], levels = ORDER_ALL)
for (t in c("Acute", "Chronic")) {
  tb <- table(d$size_cls[d$Tier == t])
  P(sprintf("size-class split %s: %s", t,
            paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse = " ")))
}

fmat <- function(v) formatC(v, format = "f", digits = 1, big.mark = ",")

# power-of-3 log breaks; label "1" for 10^0, expression 10^k otherwise
make_breaks <- function(lo, hi) {
  p3 <- seq(ceiling(lo / 3) * 3, floor(hi / 3) * 3, by = 3)
  if (!length(p3)) p3 <- 0
  brk <- 10^p3
  lab <- lapply(p3, function(k) if (k == 0) "1" else bquote(10^.(k)))
  # minor ticks: 2-9 within each decade inside the plotted range
  dec <- 10^(min(floor(lo)):max(ceiling(hi)))
  minor <- unique(as.vector(sapply(dec, function(k) (2:9) * k)))
  minor <- minor[minor >= 10^lo & minor <= 10^hi]
  list(breaks = brk, labels = do.call(expression, lab), minor = minor)
}

# observed-point plotting positions: ssdtools native convention (r-0.5)/n
get_pp <- function(ds, fit, n) {
  PP_SRC <<- "ssdtools native (r-0.5)/n (probe6-verified)"
  (seq_len(n) - 0.5) / n
}

panel_fn <- function(tier, track, tagtext, hc5_locked, n_exp) {
  col_name <- if (track == "GM") "gm" else "mn"
  dv <- d[d$Tier == tier, c("Species", col_name, "size_cls")]
  names(dv) <- c("Species", "Conc", "size_cls")
  n <- nrow(dv)
  stopifnot(n == n_exp)
  ds <- dv[order(dv$Conc), ]
  fit <- ssd_fit_dists(ds, dists = c("lnorm", "llogis"))
  gl <- as.data.frame(glance(fit, wt = TRUE))
  P(sprintf("-- %s-%s fit weights: %s", tier, track,
            paste(sprintf("%s=%.3f", gl$dist, gl$wt), collapse = " ")))
  pr  <- seq(0.01, 0.99, by = 0.01)
  pd0 <- as.data.frame(predict(fit, ci = FALSE, proportion = pr))
  ci <- tryCatch({
    set.seed(20260929)
    as.data.frame(predict(fit, ci = TRUE, nboot = 1000, proportion = pr))
  }, error = function(e) { P(paste("ci ERR:", conditionMessage(e))); NULL })
  stopifnot(!is.null(ci))
  hc_r <- pd0$est[pd0$proportion == 0.05]
  P(sprintf("%s-%s: ssdtools HC5=%.4g | python-locked HC5=%.4g | ratio=%.3f",
            tier, track, hc_r, hc5_locked, hc_r / hc5_locked))
  stopifnot(abs(log(hc_r) - log(hc5_locked)) < 0.03)
  stopifnot(all(is.finite(ci$lcl)), all(ci$lcl < ci$ucl))
  pp  <- get_pp(ds, fit, n)
  obs <- data.frame(Conc = ds$Conc, pp = pp, size_cls = ds$size_cls,
                    Species = ds$Species)
  lo  <- floor(min(log10(c(obs$Conc, pd0$est, hc5_locked))))
  hi  <- ceiling(max(log10(c(obs$Conc, pd0$est, hc5_locked * 1.5))))
  hi  <- hi + 1.35    # room for right-side species labels
  bk  <- make_breaks(lo, hi)
  band <- data.frame(x = c(ci$lcl, rev(ci$ucl)),
                     y = c(ci$proportion, rev(ci$proportion)))
  gp <- ggplot()
  gp <- gp + geom_polygon(data = band, aes(x = x, y = y, fill = "95% CI"),
                          alpha = 0.45, color = NA)
  gp <- gp + geom_line(data = ci, aes(x = est, y = proportion,
                 color = "Model-averaged SSD"), linewidth = 0.9)
  gp <- gp + geom_point(data = obs,
                        aes(x = Conc, y = pp, color = size_cls, fill = size_cls),
                        size = 2.7, stroke = 0.7)
  gp <- gp + geom_text_repel(
    data = obs, aes(x = Conc, y = pp, label = Species),
    size = 2.1, fontface = "italic", color = "grey15",
    min.segment.length = 0.25, max.overlaps = Inf,
    box.padding = 0.28, point.padding = 0.22, force = 1.2,
    xlim = c(NA, 10^hi), segment.color = "grey62", segment.size = 0.3)
  gp <- gp + geom_vline(xintercept = hc5_locked, linetype = "dashed",
                        color = "grey25", linewidth = 0.65)
  gp <- gp + scale_x_log10(breaks = bk$breaks, labels = bk$labels,
                           minor_breaks = bk$minor,
                           limits = c(10^lo, 10^hi),
                           expand = expansion(mult = 0.015))
  gp <- gp + scale_y_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.2),
                                minor_breaks = seq(0.1, 0.9, 0.1),
                                labels = function(v) paste0(v * 100, "%"),
                                expand = expansion(mult = 0.01))
  pres <- ORDER_ALL[ORDER_ALL %in% as.character(unique(ds$size_cls))]   # fixed order, presence-based
  brk2 <- c("95% CI", pres)
  P(sprintf("%s-%s legend keys: %s", tier, track, paste(brk2, collapse = ", ")))
  gp <- gp + scale_fill_manual(values = c("95% CI" = "grey82", CLASS_COLS),
                               breaks = brk2, labels = brk2, name = NULL, drop = FALSE)
  gp <- gp + scale_colour_manual(values = CLASS_COLS, name = NULL, drop = FALSE)
  gp <- gp + guides(fill = guide_legend(order = 1,
                     override.aes = list(colour = c("grey85", unname(CLASS_COLS[pres])))),
                    colour = "none")
  gp <- gp + labs(x = paste0("Concentration (particles/L)\nHC5 = ",
                             fmat(hc5_locked)),
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
    plot.margin = margin(26, 10, 4, 4))
  P(sprintf("%s-%s pp=%s ; x-range=10^%.1f..10^%.1f ; breaks=%s", tier, track,
            PP_SRC, lo, hi,
            paste(sprintf("%.0e", bk$breaks), collapse = ",")))
  list(plot = gp, hc5_ssdc = hc_r)
}

pa <- panel_fn("Acute",    "GM",  "(a) Acute - species GM",
               h$hc5[h$Tier == "Acute"    & h$Track == "GM"],  13)
pb <- panel_fn("Acute",    "MIN", "(b) Acute - species MIN",
               h$hc5[h$Tier == "Acute"    & h$Track == "MIN"], 13)
pc <- panel_fn("Chronic",  "GM",  "(c) Chronic - species GM",
               h$hc5[h$Tier == "Chronic"  & h$Track == "GM"],  15)
pd <- panel_fn("Chronic",  "MIN", "(d) Chronic - species MIN",
               h$hc5[h$Tier == "Chronic"  & h$Track == "MIN"], 15)
fig <- (pa$plot | pb$plot) / (pc$plot | pd$plot)  # true 2x2: / = vertical stack of
# two horizontal pairs. Do NOT add outer plot_layout(ncol=2): on a /-nested
# patchwork it arranges the two subassemblies SIDE BY SIDE -> 1x4 single row.
png_path <- file.path(OUTD, "fig2_ssdtools_v12.png")
set.seed(20260929)   # ggrepel label layout reproducibility
ggsave(png_path, fig, width = 9.2, height = 8.2, dpi = 300, bg = "white")
audit <- sprintf(paste0('{"figure":"fig2_ssdtools_v12","engine":"ssdtools %s (R %s)",',
                        '"data":"ssd_pooled_reps_v9.csv (28 rows verbatim; HC5 lock ssd_pooled_hc5_v9.csv) + fig2_sizeclass_v2.csv (v2 comparable-study-series re-judgment, unit=species x SSD panel; dates=run qualifiers)",',
                        '"panels":["Acute GM n=13","Acute MIN n=13",',
                        '"Chronic GM n=15","Chronic MIN n=15"],',
                        '"ci":"parametric bootstrap 95%% nboot=1000 seed=20260929",',
                        '"hc5_displayed":{"acute_gm_particles_per_L":%.1f,',
                        '"acute_min_particles_per_L":%.1f,',
                        '"chronic_gm_particles_per_L":%.1f,',
                        '"chronic_min_particles_per_L":%.1f},',
                        '"hc5_ssdc_recomputed":{"acute_gm":%.4g,"acute_min":%.4g,',
                        '"chronic_gm":%.4g,"chronic_min":%.4g},',
                        '"size_class_colors":{"Multiple-size data":"%s","Single-size data":"%s","Both size-data types":"%s","Unclassified":"%s"},',
                        '"labels":"all species, ggrepel italic, seed=20260929",',
                        '"pp_convention":"%s","png":"fig2_ssdtools_v12.png",',
                        '"note":"v11 class rebuild per INFO v2 2026-10-14: classification unit = comparable study series (dates = run qualifiers, Within/Across-run to supplement). Acute M8/S5/B0/U0, Chronic M6/S5/B0/U4. Values/tracks/HC5/CI identical to v9/v10"}'),
                 as.character(packageVersion("ssdtools")), getRversion(),
                 h$hc5[h$Tier == "Acute"   & h$Track == "GM"],
                 h$hc5[h$Tier == "Acute"   & h$Track == "MIN"],
                 h$hc5[h$Tier == "Chronic" & h$Track == "GM"],
                 h$hc5[h$Tier == "Chronic" & h$Track == "MIN"],
                 pa$hc5_ssdc, pb$hc5_ssdc, pc$hc5_ssdc, pd$hc5_ssdc,
                 COL_MULTI, COL_SINGLE, COL_BOTH, COL_UNCLASS, PP_SRC)
writeLines(audit, file.path(OUTD, "fig2_ssdtools_v12_audit.json"))
writeLines(LOGC, file.path(OUTD, "fig2_ssdtools_v12_log.txt"))
P(paste0("[SAVED] ", png_path))