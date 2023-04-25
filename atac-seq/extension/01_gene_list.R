library(readxl)
library(biomaRt)
library(ggplot2)
library(VennDiagram)

# paper genes
villiers.diff_genes <- read_excel(
    "./in_data/41467_2023_36262_MOESM4_ESM.xlsx",
    sheet = "Fig 1 U937-PR9 DEGs",
    skip = 1
)

villiers.diff_genes <- villiers.diff_genes [villiers.diff_genes $adj.P.Val < 0.05, ]
villiers.diff_genes_up <- villiers.diff_genes [villiers.diff_genes $logFC > 0, ]
villiers.diff_genes_down <- villiers.diff_genes [villiers.diff_genes $logFC < 0, ]

# atlas genes
tan.diff_genes <- read_excel("./in_data/bloodbld2020005698-suppl5.xlsx", skip = 1)
tan.diff_genes_up <- tan.diff_genes[tan.diff_genes$logFC > 0, ]
tan.diff_genes_down <- tan.diff_genes[tan.diff_genes$logFC < 0, ]

# matching
ensembl_dataset <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attributes <- c("entrezgene_id", "ensembl_gene_id")

villiers.diff_genes.diff_entrez_ids <- getBM(
    attributes = attributes,
    values = villiers.diff_genes$GeneID,
    mart = ensembl_dataset,
    filters = "entrezgene_id"
)$ensembl_gene_id

villiers.diff_genes_up.up_entrez_ids <- getBM(
    attributes = attributes,
    values = villiers.diff_genes_up$GeneID,
    mart = ensembl_dataset,
    filters = "entrezgene_id"
)$ensembl_gene_id

villiers.diff_genes_down.down_entrez_ids <- getBM(
    attributes = attributes,
    values = villiers.diff_genes_down$GeneID,
    mart = ensembl_dataset,
    filters = "entrezgene_id"
)$ensembl_gene_id

# Overlaps using names
sum(tan.diff_genes$genes %in% villiers.diff_genes$Symbol)
sum(tan.diff_genes_up$genes %in% villiers.diff_genes_up$Symbol)
sum(tan.diff_genes_down$genes %in% villiers.diff_genes_down$Symbol)

# overlaps using entrez ids
sum(tan.diff_genes$GID %in% villiers.diff_genes.diff_entrez_ids)
sum(tan.diff_genes_up$GID %in% villiers.diff_genes_up.up_entrez_ids)
sum(tan.diff_genes_down$GID %in% villiers.diff_genes_down.down_entrez_ids)

# all diff
dev.off(dev.list()["RStudioGD"])
plt <- draw.pairwise.venn(
    area1 = length(tan.diff_genes$genes),
    area2 = length(villiers.diff_genes$Symbol),
    cross.area = sum(tan.diff_genes$GID %in% villiers.diff_genes.diff_entrez_ids),
    category = c("Tan", "Villiers"),
    col = "black", fill = c("blue", "orange"),
    cex = 1.5,
    cat.cex = 1.5
)
ggsave(plt, file = "./00_plots/diff_genes_overlap.png", device = "png")

# all up
dev.off(dev.list()["RStudioGD"])
plt <- draw.pairwise.venn(
    area1 = length(tan.diff_genes.up$genes),
    area2 = length(villiers.diff_genes_up$Symbol),
    cross.area = sum(tan.diff_genes.up$GID %in% villiers.diff_genes_up.up_entrez_ids),
    category = c("Tan", "Villiers"),
    col = "black", fill = c("blue", "orange"),
    cex = 1.5,
    cat.cex = 1.5
)

ggsave(plt, file = "./00_plots/up_genes_overlap.png", device = "png")

# all down
dev.off(dev.list()["RStudioGD"])
plt <- draw.pairwise.venn(
    area1 = length(tan.diff_genes.down$genes),
    area2 = length(villiers.diff_genes_down$Symbol),
    cross.area = sum(tan.diff_genes_down$GID %in% villiers.diff_genes_down.down_entrez_ids),
    category = c("Tan", "Villiers"),
    col = "black", fill = c("blue", "orange"),
    cex = 1.5,
    cat.cex = 1.5
)

ggsave(plt, file = "./00_plots/down_genes_overlap.png",device = "png")
