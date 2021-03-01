
# BiocManager::install("spatial")  # maybe required
# BiocManager::install("GenomicFeatures")
# BiocManager::install("TxDb.Mmusculus.UCSC.mm39.refGene")
requireNamespace("TxDb.Mmusculus.UCSC.mm39.refGene")

library(dplyr)

tx <- list(
  # Mouse
  mm = merge(
    (
      # https://doi.org/doi:10.18129/B9.bioc.org.Mm.eg.db
      org.Mm.eg.db::org.Mm.egSYMBOL2EG
    ),
    (
      TxDb.Mmusculus.UCSC.mm39.refGene::TxDb.Mmusculus.UCSC.mm39.refGene %>%
        GenomicFeatures::transcriptLengths() %>%
        dplyr::group_by(gene_id) %>%
        dplyr::summarize(avg_tx_len = mean(tx_len))
    ),
    by = "gene_id"
  )
  # data.frame with transcript lengths:
  #     gene_id        symbol avg_tx_len
  # 1 100009600         Zglp1       1010
  # 2 100009609       Vmn2r65       2538
  # 3 100009614       Gm10024        564
  # 4 100009664 F630206G17Rik       2396
  # 5    100012          Oog3       1854
  # 6    100017       Ldlrap1       2671
)
