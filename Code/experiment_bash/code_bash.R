Sys.setlocale("LC_ALL", "C")
options(stringsAsFactors=F)

### load in arguments
args=(commandArgs(TRUE))
print(args)
if(length(args)==0) {
  stop("Error: No arguments supplied!")
} else {
  n_phenotype = as.numeric(args[[1]])
  n_subject = as.numeric(args[[2]])
}

print(paste("num pheno is", n_phenotype, "; num subject is", n_subject, sep = " "))