source("./99_H_extension.R")
require("Repitools")

# Load binding sites for h3k27ac

h3k27ac_binding <- annoDF2GR(h3k27ac_binding.df)

# load paper binding sites
villiers_peaks.promoter_df <- read.table(
    "./data/SRR20653396_peaks_promoter_annotated.txt",
    header = TRUE, sep = "\t"
)
villiers_peaks.promoter <- annoDF2GR(
    villiers_peaks.promoter_df[, -which(names(villiers_peaks.promoter_df) == "width")])

# init overlap columns
villiers_peaks.promoter$olap_hdac <- FALSE
villiers_peaks.promoter$olap_p300 <- FALSE

# find overlaps
olaps_hdac <- findOverlaps(villiers_peaks.promoter, hdac_binding)
olaps_p300 <- findOverlaps(villiers_peaks.promoter, p300_binding)

villiers_peaks.promoter$olap_hdac[queryHits(olaps_hdac)] <- TRUE
villiers_peaks.promoter$olap_p300[queryHits(olaps_p300)] <- TRUE