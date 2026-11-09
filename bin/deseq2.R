#!/usr/bin/env Rscript
# Tumour vs normal differential expression with DESeq2.
# Design ~ condition; median-of-ratios normalisation; Wald test; Benjamini-Hochberg.
# Outputs: deseq2_results.csv, significant_genes.csv, summary.txt, volcano.pdf, pca.pdf
suppressPackageStartupMessages(library(DESeq2))

args <- commandArgs(trailingOnly = TRUE)
opt <- list(counts = NULL, samplesheet = NULL, reference = "Normal", padj = 0.05)
i <- 1
while (i <= length(args)) {
  key <- sub("^--", "", args[i]); opt[[key]] <- args[i + 1]; i <- i + 2
}
alpha <- as.numeric(opt$padj)

raw <- read.csv(opt$counts, check.names = FALSE, stringsAsFactors = FALSE)
meta <- read.csv(opt$samplesheet, stringsAsFactors = FALSE, strip.white = TRUE)
stopifnot(!anyDuplicated(meta$sample), all(meta$sample %in% colnames(raw)))

counts <- as.matrix(raw[, meta$sample])
rownames(counts) <- raw$gene_id
meta$condition <- factor(meta$condition)
if (!opt$reference %in% levels(meta$condition))
  stop("Reference level '", opt$reference, "' not found in samplesheet conditions: ",
       paste(levels(meta$condition), collapse = ", "))
meta$condition <- relevel(meta$condition, ref = opt$reference)
if (nlevels(meta$condition) != 2) stop("Exactly two conditions are required")
if (any(table(meta$condition) < 2)) stop("Each condition needs at least 2 samples")
rownames(meta) <- meta$sample

dds <- DESeqDataSetFromMatrix(counts, colData = meta, design = ~ condition)
dds <- DESeq(dds, quiet = TRUE)
test_level <- setdiff(levels(meta$condition), opt$reference)
res <- results(dds, contrast = c("condition", test_level, opt$reference), alpha = alpha)

out <- data.frame(gene_id = rownames(res), gene_name = raw$gene_name[match(rownames(res), raw$gene_id)],
                  as.data.frame(res), check.names = FALSE)
out <- out[order(out$padj, na.last = TRUE), ]
write.csv(out, "deseq2_results.csv", row.names = FALSE)

sig <- out[!is.na(out$padj) & out$padj < alpha, ]
write.csv(sig, "significant_genes.csv", row.names = FALSE)

n_up <- sum(sig$log2FoldChange > 0); n_down <- sum(sig$log2FoldChange < 0)
summary_lines <- c(
  sprintf("contrast\t%s vs %s", test_level, opt$reference),
  sprintf("samples\t%d (%s)", ncol(dds), paste(sprintf("%s=%d", names(table(meta$condition)), table(meta$condition)), collapse = ", ")),
  sprintf("genes_tested\t%d", nrow(out)),
  sprintf("genes_with_padj\t%d", sum(!is.na(out$padj))),
  sprintf("significant_padj_lt_%s\t%d", alpha, nrow(sig)),
  sprintf("up_in_%s\t%d", test_level, n_up),
  sprintf("down_in_%s\t%d", test_level, n_down),
  sprintf("deseq2_version\t%s", as.character(packageVersion("DESeq2")))
)
writeLines(summary_lines, "summary.txt")
cat(summary_lines, sep = "\n")

# Volcano plot (base graphics, no extra packages)
pdf("volcano.pdf", width = 7, height = 6)
ok <- !is.na(out$padj)
y <- -log10(pmax(out$padj[ok], 1e-300))
col <- ifelse(out$padj[ok] < alpha, ifelse(out$log2FoldChange[ok] > 0, "#c0392b", "#2c6fbb"), "grey70")
plot(out$log2FoldChange[ok], y, pch = 20, cex = 0.6, col = col,
     xlab = sprintf("log2 fold change (%s / %s)", test_level, opt$reference), ylab = "-log10 adjusted p",
     main = sprintf("%d significant genes (padj < %s)", nrow(sig), alpha))
abline(h = -log10(alpha), lty = 2, col = "grey40")
top <- head(sig[sig$gene_name != "", ], 10)
if (nrow(top)) text(top$log2FoldChange, -log10(pmax(top$padj, 1e-300)), top$gene_name, cex = 0.6, pos = 3)
dev.off()

# PCA of variance-stabilised counts
vsd <- if (nrow(dds) >= 1000) vst(dds, blind = TRUE) else varianceStabilizingTransformation(dds, blind = TRUE)
mat <- assay(vsd)
mat <- mat[order(apply(mat, 1, var), decreasing = TRUE)[seq_len(min(500, nrow(mat)))], , drop = FALSE]
pc <- prcomp(t(mat))
pv <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
pdf("pca.pdf", width = 6, height = 5.5)
par(mar = c(7, 4, 4, 2))
cols <- c("#2c6fbb", "#c0392b")[as.integer(meta$condition)]
pad <- function(v) { r <- range(v); r + c(-0.15, 0.15) * diff(r) }
plot(pc$x[, 1], pc$x[, 2], pch = 19, col = cols, xlim = pad(pc$x[, 1]), ylim = pad(pc$x[, 2]), xlab = sprintf("PC1 (%s%%)", pv[1]), ylab = sprintf("PC2 (%s%%)", pv[2]),
     main = "PCA, top variable genes")
text(pc$x[, 1], pc$x[, 2], meta$sample, pos = 3, cex = 0.6)
legend("bottom", legend = levels(meta$condition), col = c("#2c6fbb", "#c0392b"), pch = 19, bty = "n", horiz = TRUE, inset = c(0, -0.32), xpd = TRUE)
dev.off()
