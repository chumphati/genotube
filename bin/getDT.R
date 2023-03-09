#!/usr/bin/env Rscript

input <- as.character(commandArgs(trailingOnly=TRUE))
output_name <- paste(strsplit(input, ".", fixed = T)[[1]][1], "refined.txt", sep = ".")
strain_info <- read.table(input, header = T)

DR_type <- "SENSIBLE"
#MDR
if(strain_info[,"isoniazid"] & strain_info[,"rifampicin"]) {
  DR_type <- "MDR"
}
if((strain_info[,"isoniazid"] & strain_info[,"rifampicin"]) & ((strain_info[,"fluoroquinolones"]) | (strain_info[,"amikacin"] | strain_info[,"capreomycin"] | strain_info[,"kanamycin"]))) {
  DR_type <- "PRE-XDR"
}
if(strain_info[,"isoniazid"] & strain_info[,"rifampicin"] & strain_info[,"fluoroquinolones"] & (strain_info[,"amikacin"] | strain_info[,"capreomycin"] | strain_info[,"kanamycin"])) {
  DR_type <- "XDR"
}
output <- as.character(c(strain_info[1], DR_type, strain_info[-1]))

write(output, output_name, sep = "\t", ncolumns = 90)
