library(R.utils)
library(rmarkdown)

argus <- (commandArgs(asValues=TRUE, excludeReserved=TRUE)[-1])
sample_name <- as.character(argus[1])


#run the markdown

rmarkdown::render('../markdown_files/parallelized_reporter3.rmd',output_file=paste0('../patient_reports/',sample_name,'.html'),params=list(sample_id=sample_name))
