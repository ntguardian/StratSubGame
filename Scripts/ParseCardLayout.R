#!/usr/bin/Rscript
# ProcessCardLayout.R
# 2022-10-04
# curtis
# Process a card layout CSV file

# argparser: A package for handling command line arguments
if (!suppressPackageStartupMessages(require("argparser"))) {
  install.packages("argparser")
  require("argparser")
}
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))

# EXECUTABLE SCRIPT MAIN FUNCTIONALITY -----------------------------------------

main <- function(type, inputfile) {
  input <- read.csv(inputfile)
  stopifnot(all(c("Sheet", "Card") %in% names(input)))

  for (i in unique(input$Sheet)) {
    filter(input, Sheet == i) %>%
      pull(Card) %>%
      cat(paste0(type, i), ., "\n")
  }
}

# INTERFACE DEFINITION AND COMMAND LINE IMPLEMENTATION -------------------------

if (sys.nframe() == 0) {
  p <- arg_parser("Process a card layout CSV file")
  p <- add_argument(p, "type", type = "character", nargs = 1,
                    help = "The type of card being processed")
  p <- add_argument(p, "inputfile", type = "character", nargs = 1,
                    help = "The input CSV file")

  cl_args <- parse_args(p)
  cl_args <- cl_args[!(names(cl_args) %in% c("help", "opts"))]
  if (any(sapply(cl_args, is.na))) {
    # User did not specify all inputs; print help message
    print(p)
    cat("\n\nNot all needed inputs were given.\n")
    quit()
  }

  do.call(main, cl_args[2:length(cl_args)])
}

