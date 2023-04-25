source("./99_H_extension.R")

require("Repitools")

# paper genes
paper.peaks <- narrow_to_granges("./data/SRR20653396_peaks.narrowPeak")

# annotate
paper.peaks <- anno_peaks(paper.peaks)
paper.peaks_df <- annoGR2DF(paper.peaks)
write.table(paper.peaks_df,
    file = "./data/SRR20653396_peaks_annotated.txt",
    sep = "\t", row.names = FALSE
)

# plot
pie_from_vec(paper.peaks$annotation_simple, "./00_plots/binding_annotation.png")
