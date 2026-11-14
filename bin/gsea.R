#!/usr/bin/env Rscript
# Pre-ranked gene set enrichment (fgsea) on the DESeq2 Wald statistic.
# Genes tested by DESeq2 are ranked by the Wald statistic (sign = direction, size = evidence); a gene name
# that occurs more than once keeps its entry with the largest |stat|.
# Outputs: gsea_results.tsv, gsea_summary.txt, gsea_nes.pdf
suppressPackageStartupMessages(library(fgsea))

args <- commandArgs(trailingOnly = TRUE)
opt <- list(results = NULL, gmt = NULL, padj = 0.05, min_size = 15, max_size = 500, seed = 1)
i <- 1
while (i <= length(args)) {
  key <- sub("^--", "", args[i]); opt[[key]] <- args[i + 1]; i <- i + 2
}
alpha <- as.numeric(opt$padj)

res <- read.csv(opt$results, stringsAsFactors = FALSE)
# Rank only genes that DESeq2 tested (non-NA padj after independent filtering): very low
# counts give many tied statistics and carry no information.
res <- res[!is.na(res$stat) & !is.na(res$padj) & !is.na(res$gene_name) & res$gene_name != "", ]
res <- res[order(-abs(res$stat)), ]
res <- res[!duplicated(res$gene_name), ]
ranks <- setNames(res$stat, res$gene_name)
ranks <- sort(ranks, decreasing = TRUE)

sets <- gmtPathways(opt$gmt)

set.seed(as.integer(opt$seed))
gs <- fgsea(pathways = sets, stats = ranks, minSize = as.integer(opt$min_size),
            maxSize = as.integer(opt$max_size), eps = 0)
gs <- gs[order(gs$padj, -abs(gs$NES)), ]

out <- data.frame(pathway = gs$pathway, size = gs$size, NES = round(gs$NES, 3), ES = round(gs$ES, 3),
                  pval = signif(gs$pval, 4), padj = signif(gs$padj, 4),
                  leading_edge = vapply(gs$leadingEdge, function(x) paste(head(x, 15), collapse = ","), ""))
write.table(out, "gsea_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

sig <- out[!is.na(out$padj) & out$padj < alpha, ]
summary_lines <- c(
  "method\tfgsea (pre-ranked), ranking = DESeq2 Wald statistic",
  sprintf("genes_ranked\t%d", length(ranks)),
  sprintf("gene_sets_in_file\t%d", length(sets)),
  sprintf("gene_sets_tested\t%d (size %s-%s after matching)", nrow(out), opt$min_size, opt$max_size),
  sprintf("significant_padj_lt_%s\t%d", alpha, nrow(sig)),
  sprintf("enriched_up\t%d", sum(sig$NES > 0)),
  sprintf("enriched_down\t%d", sum(sig$NES < 0)),
  sprintf("fgsea_version\t%s", as.character(packageVersion("fgsea")))
)
writeLines(summary_lines, "gsea_summary.txt")
cat(summary_lines, sep = "\n")

# Bar chart of normalised enrichment scores for significant sets (base graphics)
pdf("gsea_nes.pdf", width = 7.5, height = max(3, 0.22 * max(nrow(sig), 1) + 1.8))
if (nrow(sig)) {
  sig <- sig[order(sig$NES), ]
  par(mar = c(4.5, 16, 3, 1))
  barplot(sig$NES, names.arg = sub("^HALLMARK_", "", sig$pathway), horiz = TRUE, las = 1, cex.names = 0.6,
          col = ifelse(sig$NES > 0, "#c0392b", "#2c6fbb"), border = NA, xlab = "Normalised enrichment score",
          main = sprintf("%d gene sets at padj < %s", nrow(sig), alpha))
  abline(v = 0, col = "grey40")
} else {
  plot.new(); title(sprintf("No gene set at padj < %s", alpha))
}
invisible(dev.off())
