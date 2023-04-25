source("./99_H_extension.R")

library(VennDiagram)
library(ggplot2)
library(readxl)

# load for atlas in g range
tan.peaks <- DataFrame(read_excel("./in_data/bloodbld2020005698-suppl4.xlsx", skip = 1))
tan.peaks <- GRanges(seqnames = tan.peaks$Chr, 
                     ranges = IRanges(start = tan.peaks$Start, end = tan.peaks$End))


# load for paper in Granges
villiers.peaks <- narrow_to_granges( "./data/SRR20653396_peaks.narrowPeak")


peak_overlap <- findOverlaps(tan.peaks, villiers.peaks)

# plot
dev.off(dev.list()["RStudioGD"])
plt <- draw.pairwise.venn(
    area1 = length(tan.peaks),
    area2 = length(villiers.peaks),
    cross.area = length(queryHits(peak_overlap)),
    category = c("Tan", "Me"),
    col = "black", fill = c("blue", "orange"),
    cex = 1.5,
    cat.cex = 1.5
)

ggsave(plt, file = "./00_plots/overlapping_binding.png", device = "png")
