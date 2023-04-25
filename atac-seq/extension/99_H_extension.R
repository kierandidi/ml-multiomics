library(ChIPseeker)
library(GenomicRanges)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)

pie_from_vec <- function(vec, p, main = "") {
    counts <- table(vec)
    percentages <- round(prop.table(counts) * 100)

    jpeg(filename = p, width = 480, height = 480)
    pie(counts,
        labels = paste0(names(counts), " (", percentages, "%)"),
        main = main
    )
    dev.off()
}

anno_peaks <- function(peaks) {
    EnsHg38 <- EnsDb.Hsapiens.v86
    seqlevelsStyle(EnsHg38) <- "UCSC"

    # annotate
    peaks <- as.GRanges(annotatePeak(peaks, TxDb = EnsHg38))

    # simplify annotation
    peaks$annotation_simple <- sapply(peaks$annotation, function(x) {
        strsplit(x, " (", fixed = TRUE)[[1]][1]
    })
    return(peaks)
}

narrow_to_granges <- function(path) {
    peaks <- read.table(
       path,
        header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = ""
    )

    mcols <- DataFrame(
        signalValue = peaks$V7,
        pValue = peaks$V8,
        qValue = peaks$V9,
        peak = peaks$V10
    )

    peaks <- GRanges(
        seqnames = peaks$V1,
        ranges = IRanges(
            start = peaks$V2,
            end = peaks$V3,
            mcols = mcols
        )
    )
    return(peaks)
}

get_atlas_table <- function(num) {
    df <- data.frame(read_excel(paste("./in_data/bloodbld2020005698-suppl", 
                                      as.character(num), ".xlsx",sep = ""
    ), skip = 1))
    colnames(df) <- tolower(colnames(df))
    return(df)
}
