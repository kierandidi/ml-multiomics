source("./99_H_extension.R")
require("Repitools")
require("DiffBind")
require(GenomicRanges)
library(soGGi)
library("ggVennDiagram")

# Load binding sites for HDAC, p300, and h3k27ac_binding
hdac_binding.df <- data.frame(get_atlas_table(8))
p300_binding.df <- data.frame(get_atlas_table(9))
h3k27ac_binding.df <- data.frame(get_atlas_table(10))
hdac_binding <- annoDF2GR(hdac_binding.df)
p300_binding <- annoDF2GR(p300_binding.df)
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
villiers_peaks.promoter$olap_h3k27ac <- FALSE

# find overlaps
olaps_p300 <- findOverlaps(villiers_peaks.promoter, p300_binding)
olaps_h3k27ac <- findOverlaps(villiers_peaks.promoter, h3k27ac_binding)

villiers_peaks.promoter$olap_hdac[queryHits(olaps_hdac)] <- TRUE
villiers_peaks.promoter$olap_p300[queryHits(olaps_p300)] <- TRUE
villiers_peaks.promoter$olap_h3k27ac[queryHits(olaps_h3k27ac)] <- TRUE

# cmoute stats
sum(villiers_peaks.promoter$gene_control & villiers_peaks.promoter$olap_hdac)
sum(villiers_peaks.promoter$gene_control & villiers_peaks.promoter$olap_p300)
sum(villiers_peaks.promoter$gene_control & villiers_peaks.promoter$olap_h3k27ac)
sum(villiers_peaks.promoter$gene_control & villiers_peaks.promoter$olap_p300 & villiers_peaks.promoter$olap_hdac)

# down
sum(villiers_peaks.promoter$gene_control_down & villiers_peaks.promoter$olap_hdac)
sum(villiers_peaks.promoter$gene_control_down & villiers_peaks.promoter$olap_p300)
sum(villiers_peaks.promoter$gene_control_down & villiers_peaks.promoter$olap_p300 & villiers_peaks.promoter$olap_hdac)

dev.off(dev.list()["RStudioGD"])
plt <- draw.triple.venn(
    area1 = sum(villiers_peaks.promoter$gene_control_down),
    area2 = sum(villiers_peaks.promoter$gene_control_down & villiers_peaks.promoter$olap_hdac),
    area3 = sum(villiers_peaks.promoter$gene_control_down & villiers_peaks.promoter$olap_p300),
    n12 = sum(villiers_peaks.promoter$gene_control_down & villiers_peaks.promoter$olap_hdac),
    n13 = sum(villiers_peaks.promoter$gene_control_down & villiers_peaks.promoter$olap_p300),
    n23 = sum(villiers_peaks.promoter$gene_control_down & villiers_peaks.promoter$olap_p300 & villiers_peaks.promoter$olap_hdac),
    n123 = sum(villiers_peaks.promoter$gene_control_down & villiers_peaks.promoter$olap_p300 & villiers_peaks.promoter$olap_hdac),
    col = "black", fill = c("blue", "orange", "red"),
    category = c("PLM::RARA", "PLM::RARA and HDAC", "PLM::RARA and p300"),
    cex = 1.5,
    cat.cex = 1.5
)

# up
sum(villiers_peaks.promoter$gene_control_up & villiers_peaks.promoter$olap_hdac)
sum(villiers_peaks.promoter$gene_control_up & villiers_peaks.promoter$olap_p300)
sum(villiers_peaks.promoter$gene_control_up & villiers_peaks.promoter$olap_p300 & villiers_peaks.promoter$olap_hdac)

dev.off(dev.list()["RStudioGD"])
plt <- draw.triple.venn(
    area1 = sum(villiers_peaks.promoter$gene_control_up),
    area2 = sum(villiers_peaks.promoter$gene_control_up & villiers_peaks.promoter$olap_hdac),
    area3 = sum(villiers_peaks.promoter$gene_control_up & villiers_peaks.promoter$olap_p300),
    n12 = sum(villiers_peaks.promoter$gene_control_up & villiers_peaks.promoter$olap_hdac),
    n13 = sum(villiers_peaks.promoter$gene_control_up & villiers_peaks.promoter$olap_p300),
    n23 = sum(villiers_peaks.promoter$gene_control_up & villiers_peaks.promoter$olap_p300 & villiers_peaks.promoter$olap_hdac),
    n123 = sum(villiers_peaks.promoter$gene_control_up & villiers_peaks.promoter$olap_p300 & villiers_peaks.promoter$olap_hdac),
    col = "black", fill = c("blue", "orange", "red"),
    category = c("PLM::RARA", "PLM::RARA and HDAC", "PLM::RARA and p300"),
    cex = 1.5,
    cat.cex = 1.5
)

# PLOT1: They coolocate
hdac1 <- narrow_to_granges("/Users/maxwuerfek/code/extension/data/SRR8588842_peaks.narrowPeak")
p300 <- narrow_to_granges("/Users/maxwuerfek/code/extension/data/SRR8588843_peaks.narrowPeak")
ext_dba <- dba(sampleSheet = "./samplesheet.csv")

sites <- GRangesList("PML-RARA down regulated"=villiers_peaks.promoter[villiers_peaks.promoter$gene_control_down],
                     "PML-RARA up regulated"=villiers_peaks.promoter[villiers_peaks.promoter$gene_control_up,])
    
profiles <- dba.plotProfile(ext_dba, sites=sites, labels = c("HDAC1", "p300"),
                            samples = list(
    HDAC1=ext_dba$masks$`anti-HDAC1`,
    p300=ext_dba$masks$`anti-P300`
))
dba.plotProfile(profiles)

# PLOT2
pdata_hdac1 <- regionPlot("/Users/maxwuerfek/code/extension/data/SRR8588842_sorted.bam",
                        testRanges = villiers_peaks.promoter,
                        style = "point",
                        format = "bam",
                        distanceAround = 1000)

pdata_p300 <- regionPlot("/Users/maxwuerfek/code/extension/data/SRR8588843_sorted.bam",
                         testRanges = villiers_peaks.promoter,
                         style = "point",
                         format = "bam",
                         distanceAround = 1000)

comb <- c(pdata_hdac1, pdata_p300)

p <- plotRegion(comb, 
                colourBy = "Sample", 
                gts = list(down=villiers_peaks.promoter[villiers_peaks.promoter$gene_control_down,],
                                       up=villiers_peaks.promoter[villiers_peaks.promoter$gene_control_up,]),
                groupBy = "Group"
                )
p

# PLOT3
pdata_h3k27ac <- regionPlot("/Users/maxwuerfek/code/extension/data/SRR8588844_sorted.bam",
                          testRanges = villiers_peaks.promoter,
                          style = "point",
                          format = "bam",
                          distanceAround = 1000)
p <- plotRegion(pdata_h3k27ac, 
                colourBy = "Group", 
                gts = list(down=villiers_peaks.promoter[villiers_peaks.promoter$gene_control_down,],
                           up=villiers_peaks.promoter[villiers_peaks.promoter$gene_control_up,]),
)
p
