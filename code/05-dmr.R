#!/usr/bin/env Rscript
# Differential methylation between vaccinated (Vac, 600 ATU) and
# pre-vaccination (T0, 0 ATU) Atlantic salmon parr with DSS.
#
# Usage:
#   Rscript code/05-dmr.R                       # defaults below
#   MIN_COV=10 DELTA=0.2 Rscript code/05-dmr.R  # override settings
#
# Input: the strand-merged CpG coverage files written by
#   code/04-dedup-methylation.sh extract
#   output/04-methylation/<lib>.CpG.CpG_report.merged_CpG_evidence.cov.gz
#   (chr, start, end, % methylation, methylated count, unmethylated count)
#
# Steps:
#   1. Read each library. Keep CpGs with >= MIN_COV reads in all four
#      libraries and drop CpGs above the HI_PERC coverage percentile in any
#      library (likely repeats or collapsed regions).
#   2. Sample QC: per-library summary, correlation, PCA. SsParrT0_1 was
#      sequenced on a different run; if it separates from the other three,
#      treat the results with care.
#   3. DSS DMLtest with smoothing (Feng et al. 2014). Smoothing borrows
#      information from neighboring CpGs, which is what makes n = 2 per
#      group workable.
#   4. Differentially methylated loci (callDML): posterior probability that
#      |Vac - T0| > DELTA is >= 1 - P_DML.
#      Differentially methylated regions (callDMR): >= DMR_MIN_CPG CpGs,
#      >= DMR_MIN_LEN bp, P_DMR on the per-CpG test.
#   5. Every DML/DMR gets a 'consistent' flag: TRUE when both Vac replicates
#      lie on the same side of both T0 replicates (raw, unsmoothed values).
#      With two replicates per group this is the more trustworthy subset.
#
# The batch effect cannot be modeled: the only library from the other run is
# in the T0 group, so a batch covariate would be confounded with condition.
# Differences are Vac minus T0 (positive = higher methylation after
# vaccination), on a 0-1 scale.
#
# Output: output/05-dmr/
#   sample-summary.tsv, correlation.tsv, pca.tsv, pca-variance.tsv, qc-plots.pdf
#   dmltest.tsv.gz                    per-CpG test results (all tested CpGs)
#   dml.tsv, dml.bed                  differentially methylated loci
#   dmr.tsv, dmr.bed                  differentially methylated regions
#   summary.tsv                       counts for the dashboard
#   dss-objects.rds                   BSseq object and DMLtest result

suppressPackageStartupMessages({
  library(DSS)
  library(bsseq)
  library(data.table)
  library(GenomicRanges)
})

env_num <- function(name, default) as.numeric(Sys.getenv(name, default))

REPO    <- Sys.getenv("REPO", "/home/shared/8TB_HDD_03/sr320/github/project-salmon-vaccine-wgbs")
IN_DIR  <- Sys.getenv("IN_DIR", file.path(REPO, "output/04-methylation"))
OUT_DIR <- Sys.getenv("OUT_DIR", file.path(REPO, "output/05-dmr"))
MIN_COV <- env_num("MIN_COV", 5)        # reads per CpG in every library
HI_PERC <- env_num("HI_PERC", 99.9)     # drop CpGs above this coverage percentile
DELTA   <- env_num("DELTA", 0.25)       # methylation difference (0-1) for DMLs
P_DML   <- env_num("P_DML", 0.001)      # DML: 1 - posterior probability
P_DMR   <- env_num("P_DMR", 0.01)       # DMR: per-CpG p-value threshold
DMR_DELTA   <- env_num("DMR_DELTA", 0.1)  # minimum mean difference for a DMR
DMR_MIN_LEN <- env_num("DMR_MIN_LEN", 50)
DMR_MIN_CPG <- env_num("DMR_MIN_CPG", 3)
DMR_MERGE   <- env_num("DMR_MERGE", 100)  # merge DMRs closer than this (bp)
CORES   <- env_num("CORES", 8)

samples <- data.table(
  id        = c("SsParrT0_1", "SsParrT0_2", "SsParrVac_1", "SsParrVac_2"),
  condition = c("T0", "T0", "Vac", "Vac")
)
samples[, file := file.path(IN_DIR, paste0(id, ".CpG.CpG_report.merged_CpG_evidence.cov.gz"))]
missing <- samples[!file.exists(file), file]
if (length(missing)) stop("missing input:\n  ", paste(missing, collapse = "\n  "))
vac <- samples[condition == "Vac", id]
t0  <- samples[condition == "T0", id]

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) cat(format(Sys.time(), "%F %T"), "|", ..., "\n")
log_msg("settings: MIN_COV", MIN_COV, "HI_PERC", HI_PERC, "DELTA", DELTA,
        "P_DML", P_DML, "P_DMR", P_DMR, "DMR_DELTA", DMR_DELTA,
        "DMR_MIN_LEN", DMR_MIN_LEN, "DMR_MIN_CPG", DMR_MIN_CPG, "CORES", CORES)

# ---- 1. Read and filter ----------------------------------------------------
read_cov <- function(f) {
  d <- fread(cmd = paste("zcat", shQuote(f)), header = FALSE,
             select = c(1, 2, 5, 6), col.names = c("chr", "pos", "X", "U"))
  d[, N := X + U][, U := NULL]
  d
}
cov <- lapply(setNames(samples$file, samples$id), function(f) {
  log_msg("reading", basename(f)); read_cov(f)
})

summ <- samples[, .(id, condition)]
summ[, cpgs_ge_min_cov := sapply(cov, function(d) d[N >= MIN_COV, .N])]
summ[, median_cov := sapply(cov, function(d) median(d[N >= MIN_COV, N]))]

# CpGs passing the coverage filters in every library.
keep <- NULL
for (id in samples$id) {
  d  <- cov[[id]]
  hi <- quantile(d$N, HI_PERC / 100)
  k  <- d[N >= MIN_COV & N <= hi, .(chr, pos)]
  keep <- if (is.null(keep)) k else fintersect(keep, k)
}
setkey(keep, chr, pos)
log_msg("CpGs passing coverage filters in all libraries:", nrow(keep))
if (nrow(keep) == 0) stop("no CpGs pass the coverage filters")

# Every library now has the same CpGs in the same (sorted) order, so the
# BSseq object can be built directly instead of via makeBSseqData's merge.
cov <- lapply(cov, function(d) { setkey(d, chr, pos); d[keep, nomatch = NULL] })
stopifnot(all(sapply(cov, nrow) == nrow(keep)))
bs  <- BSseq(chr = keep$chr, pos = keep$pos,
             M   = sapply(cov, `[[`, "X"),
             Cov = sapply(cov, `[[`, "N"),
             sampleNames = samples$id)
rm(cov); invisible(gc())

# ---- 2. Sample QC ----------------------------------------------------------
M <- as.matrix(getCoverage(bs, type = "M"))
C <- as.matrix(getCoverage(bs, type = "Cov"))
pm <- 100 * M / C
colnames(pm) <- samples$id

summ[, cpgs_all_libraries := nrow(pm)]
summ[, mean_pct_meth := round(colMeans(pm), 2)]
summ[, pct_cpgs_over_80 := round(colMeans(pm > 80) * 100, 2)]
summ[, pct_cpgs_under_20 := round(colMeans(pm < 20) * 100, 2)]
fwrite(summ, file.path(OUT_DIR, "sample-summary.tsv"), sep = "\t")

cm <- cor(pm)
fwrite(data.table(sample = rownames(cm), round(cm, 4)),
       file.path(OUT_DIR, "correlation.tsv"), sep = "\t")

# PCA on the 10,000 most variable CpGs.
v   <- matrixStats::rowVars(pm)
top <- order(v, decreasing = TRUE)[seq_len(min(10000, length(v)))]
pca <- prcomp(t(pm[top, , drop = FALSE]), center = TRUE, scale. = FALSE)
ve  <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
fwrite(data.table(sample = samples$id, condition = samples$condition,
                  round(pca$x[, seq_len(min(3, ncol(pca$x))), drop = FALSE], 3)),
       file.path(OUT_DIR, "pca.tsv"), sep = "\t")
fwrite(data.table(pc = paste0("PC", seq_along(ve)), pct_variance = ve),
       file.path(OUT_DIR, "pca-variance.tsv"), sep = "\t")

col <- ifelse(samples$condition == "Vac", "#2a6fdb", "#e0782b")
pdf(file.path(OUT_DIR, "qc-plots.pdf"), width = 7, height = 6)
par(mfrow = c(2, 2))
for (i in seq_along(samples$id)) {
  hist(pm[, i], breaks = 50, col = col[i], border = NA,
       main = samples$id[i], xlab = "% methylation (CpGs in all libraries)")
}
par(mfrow = c(1, 1))
plot(pca$x[, 1], pca$x[, 2], pch = 19, cex = 1.5, col = col,
     xlab = paste0("PC1 (", ve[1], "%)"), ylab = paste0("PC2 (", ve[2], "%)"),
     main = "PCA, 10,000 most variable CpGs")
text(pca$x[, 1], pca$x[, 2], samples$id, pos = 3, cex = 0.8)
plot(hclust(as.dist(1 - cm), method = "average"),
     main = "Clustering, 1 - Pearson correlation", xlab = "", sub = "")
dev.off()
log_msg("QC written")

# ---- 3. Per-CpG test -------------------------------------------------------
# group1 = Vac, group2 = T0, so diff = Vac - T0.
log_msg("DMLtest on", nrow(pm), "CpGs")
dml_test <- DMLtest(bs, group1 = vac, group2 = t0, smoothing = TRUE,
                    ncores = CORES)
fwrite(as.data.table(dml_test), file.path(OUT_DIR, "dmltest.tsv.gz"), sep = "\t")

# ---- 4-5. Calls ------------------------------------------------------------
# Raw per-sample % methylation over a set of regions, and the consistency flag.
cpg_gr <- GRanges(seqnames(bs), ranges(bs))
region_pct <- function(chr, start, end) {
  gr <- GRanges(chr, IRanges(start, end))
  ov <- findOverlaps(gr, cpg_gr)
  q  <- queryHits(ov); s <- subjectHits(ov)
  out <- sapply(samples$id, function(id) {
    as.numeric(tapply(M[s, id], factor(q, seq_along(gr)), sum) /
               tapply(C[s, id], factor(q, seq_along(gr)), sum))
  })
  matrix(round(100 * out, 2), nrow = length(gr),
         dimnames = list(NULL, paste0("pct_", samples$id)))
}
consistent <- function(p) {
  a <- p[, paste0("pct_", vac), drop = FALSE]
  b <- p[, paste0("pct_", t0), drop = FALSE]
  (apply(a, 1, min) > apply(b, 1, max)) | (apply(a, 1, max) < apply(b, 1, min))
}
write_bed <- function(d, prefix, file) {
  bed <- d[, .(chr, start = start - 1L, end,
               name = paste0(prefix, "_", seq_len(.N)),
               score = round(100 * diff_Vac_minus_T0, 1), strand = ".")]
  fwrite(bed, file, sep = "\t", col.names = FALSE)
}

dml <- as.data.table(callDML(dml_test, delta = DELTA, p.threshold = P_DML))
if (nrow(dml)) {
  dml[, `:=`(start = pos, end = pos)]
  dml <- cbind(dml, region_pct(dml$chr, dml$start, dml$end))
  setnames(dml, c("mu1", "mu2", "diff"), c("mu_Vac", "mu_T0", "diff_Vac_minus_T0"))
  dml[, consistent := consistent(as.matrix(.SD)), .SDcols = patterns("^pct_")]
} else {
  dml <- data.table(chr = character(), start = integer(), end = integer(),
                    diff_Vac_minus_T0 = numeric(), consistent = logical())
}
fwrite(dml, file.path(OUT_DIR, "dml.tsv"), sep = "\t")
write_bed(dml, "DML", file.path(OUT_DIR, "dml.bed"))
log_msg("DMLs:", nrow(dml))

dmr <- callDMR(dml_test, delta = DMR_DELTA, p.threshold = P_DMR,
               minlen = DMR_MIN_LEN, minCG = DMR_MIN_CPG, dis.merge = DMR_MERGE)
dmr <- if (is.null(dmr) || !nrow(dmr)) {
  data.table(chr = character(), start = integer(), end = integer(),
             diff_Vac_minus_T0 = numeric(), consistent = logical())
} else {
  d <- as.data.table(dmr)
  d <- cbind(d, region_pct(d$chr, d$start, d$end))
  setnames(d, c("meanMethy1", "meanMethy2", "diff.Methy"),
           c("mean_Vac", "mean_T0", "diff_Vac_minus_T0"))
  d[, consistent := consistent(as.matrix(.SD)), .SDcols = patterns("^pct_")]
  d
}
fwrite(dmr, file.path(OUT_DIR, "dmr.tsv"), sep = "\t")
write_bed(dmr, "DMR", file.path(OUT_DIR, "dmr.bed"))
log_msg("DMRs:", nrow(dmr))

count <- function(d, level, tested) data.table(
  level = level, tested = tested, significant = nrow(d),
  hyper_in_Vac = d[diff_Vac_minus_T0 > 0, .N],
  hypo_in_Vac  = d[diff_Vac_minus_T0 < 0, .N],
  significant_consistent = d[consistent == TRUE, .N])
summary_tab <- rbind(count(dml, "DML", nrow(dml_test)),
                     count(dmr, "DMR", nrow(dml_test)))
fwrite(summary_tab, file.path(OUT_DIR, "summary.tsv"), sep = "\t")
print(summary_tab)

saveRDS(list(samples = samples, bs = bs, dml_test = dml_test),
        file.path(OUT_DIR, "dss-objects.rds"))
log_msg("finished")
