source("./99_H_extension.R")

library(readxl)
require("Repitools")

# load the annotated paper binding sites
villiers.peaks_df <- read.table("./data/SRR20653396_peaks_annotated.txt",
    header = TRUE, sep = "\t"
)
villiers.peaks <- annoDF2GR(villiers.peaks_df[, -which(names(villiers.peaks_df) == "width")])


# get atlas differential regulated genes
tan.diff_genes <- read_excel("./in_data/bloodbld2020005698-suppl5.xlsx", skip = 1)
tan.diff_genes.down <- tan.diff_genes[tan.diff_genes$logFC > 0, ]
tan.diff_genes.up <- tan.diff_genes[tan.diff_genes$logFC < 0, ]

# we want to subset our binding sites to those with a promoter
villiers.peaks.promoter <- villiers.peaks[villiers.peaks$annotation_simple == "Promoter", ]

# check whether binding site is controlling differentially expressed gene
villiers.peaks.promoter$gene_control <- villiers.peaks.promoter$geneId %in% tan.diff_genes$GID
villiers.peaks.promoter$gene_control_up <- villiers.peaks.promoter$geneId %in% tan.diff_genes.up$GID
villiers.peaks.promoter$gene_control_down <- villiers.peaks.promoter$geneId %in% tan.diff_genes.down$GID

# save the promoter binding sites
villiers.peaks.promoter_df <- annoGR2DF(villiers.peaks.promoter)
write.table(villiers.peaks.promoter_df,
    file = "./data/SRR20653396_peaks_promoter_annotated.txt",
    sep = "\t", row.names = FALSE
)

# numbers
sum(villiers.peaks.promoter$gene_control)
sum(villiers.peaks.promoter$gene_control_up)
sum(villiers.peaks.promoter$gene_control_down)

length(unique(villiers.peaks.promoter[villiers.peaks.promoter$gene_control, ]$geneId))
length(unique(villiers.peaks.promoter[villiers.peaks.promoter$gene_control_up, ]$geneId))
length(unique(villiers.peaks.promoter[villiers.peaks.promoter$gene_control_down, ]$geneId))


# Make the chart
pipeline <- c(rep("Me NB4 Cut & Run" , 3) , rep("Tan NB4 ChipSeek" , 3) )
condition <- rep(c("Regulated" , "Up-regulated" , "Down-regulated") , 4)
num_pmlrara_regulated_genes <- c(
    length(unique(villiers.peaks.promoter[villiers.peaks.promoter$gene_control, ]$geneId)),
    length(unique(villiers.peaks.promoter[villiers.peaks.promoter$gene_control_up, ]$geneId)),
    length(unique(villiers.peaks.promoter[villiers.peaks.promoter$gene_control_down, ]$geneId)),
    787,
    424,
    363
)
data <- data.frame(pipeline,condition,num_pmlrara_regulated_genes)

# Grouped
ggplot(data, aes(fill=condition, y=num_pmlrara_regulated_genes, x=pipeline)) + 
    geom_bar(position="dodge", stat="identity") +
    geom_text(aes(label=num_pmlrara_regulated_genes), position=position_dodge(width=0.9), vjust=-0.25)