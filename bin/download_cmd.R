#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 1) {
  stop("Exactly one argument must be supplied (CMD study)", call.=FALSE)
} else {
  # default output file
  study = args[1]
}


study_obj <- curatedMetagenomicData(study, dryrun = FALSE)

# write MM to file
writeMM(assay(study_obj[[1]]), "data.mtx")

# write row and col names
write(colnames(study_obj[[1]]), 'cols.txt')
write(rownames(study_obj[[1]]), 'rows.txt')
