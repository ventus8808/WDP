#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Path to input .rds file")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$input)) {
  stop("Please provide input file with --input")
}

if (!file.exists(opt$input)) {
  stop("File not found: ", opt$input)
}

message("Loading ", opt$input, "...")
data <- readRDS(opt$input)

message("\nObject Size:")
print(object.size(data), units = "auto")

message("\nClass:")
print(class(data))

message("\nLength:")
print(length(data))

if (is.list(data)) {
  message("\nNames (first 10):")
  print(head(names(data), 10))

  if (length(data) > 0) {
    first_item_name <- names(data)[1]
    first_item <- data[[1]]
    message(paste0("\nStructure of first item ('", first_item_name, "'):"))
    str(first_item, max.level = 2)

    if (is.list(first_item) && "draws" %in% names(first_item)) {
       message("\nDimensions of 'draws' in first item:")
       print(dim(first_item$draws))
       message("\nHead of 'draws' in first item:")
       print(head(first_item$draws))
    }
  }
} else {
  message("\nStructure:")
  str(data, max.level = 2)
}
