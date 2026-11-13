#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(GenomicRanges))

# Inputs are GRCh38, 1-based inclusive. Each cache file is read once.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
    stop("Usage: Rscript showcase_ranges.R loci.tsv file_ids.txt cache_dir output.json")
}

loci <- read.delim(args[1], colClasses = c("character", "character", "integer", "integer"))
locus_ranges <- GRanges(loci$chromosome, IRanges(loci$start, loci$end))

overlap_file <- function(file_id) {
    path <- file.path(args[3], paste0(file_id, ".tsv"))
    if (!file.exists(path)) {
        return(list(error = "Segment file missing from cache"))
    }

    tryCatch({
        segments <- read.delim(path, colClasses = "character", fileEncoding = "UTF-8-BOM")
        required <- c("Chromosome", "Start", "End", "Segment_Mean")
        if (!all(required %in% names(segments))) {
            stop("Invalid GDC segment header")
        }
        chromosome <- sub("^chr", "", segments$Chromosome)
        relevant <- chromosome %in% loci$chromosome
        chromosome <- chromosome[relevant]
        starts <- suppressWarnings(as.integer(segments$Start[relevant]))
        ends <- suppressWarnings(as.integer(segments$End[relevant]))
        means <- suppressWarnings(as.numeric(segments$Segment_Mean[relevant]))
        if (anyNA(starts) || anyNA(ends) || any(starts < 1L | starts > ends)) {
            stop("Invalid segment coordinates")
        }

        ranges <- GRanges(chromosome, IRanges(starts, ends))
        hits <- findOverlaps(locus_ranges, ranges, ignore.strand = TRUE)
        evidence <- lapply(seq_len(nrow(loci)), function(index) {
            overlapping <- subjectHits(hits)[queryHits(hits) == index]
            overlapping <- overlapping[order(starts[overlapping], ends[overlapping])]
            if (any(!is.finite(means[overlapping]))) {
                return(list(error = "Nonfinite segment mean"))
            }
            records <- data.frame(
                start = starts[overlapping], end = ends[overlapping], mean = means[overlapping]
            )
            spanning <- starts[overlapping] <= loci$start[index] &
                ends[overlapping] >= loci$end[index]
            resolved <- length(overlapping) == 1L && any(spanning)
            reason <- if (resolved) {
                ""
            } else if (length(overlapping) == 0L) {
                "No overlapping segment"
            } else if (length(overlapping) > 1L) {
                "Multiple overlapping segments"
            } else {
                "Partial locus coverage"
            }
            list(mean = if (resolved) means[overlapping] else NULL,
                 reason = reason, segments = records)
        })
        setNames(evidence, loci$gene)
    }, error = function(error) list(error = conditionMessage(error)))
}

file_ids <- readLines(args[2])
result <- setNames(lapply(file_ids, overlap_file), file_ids)
jsonlite::write_json(result, args[4], auto_unbox = TRUE, digits = NA, dataframe = "rows", null = "null")
