# Current array
array.id = as.integer(commandArgs()[length(commandArgs())])

# Create simulation design matrix
num_pheno <- c(4, 6, 8)
num_sub <- c(500, 1000, 1500)
design <- expand.grid(num_sub = num_sub,
                      num_pheno = num_pheno)
n_phenotype = design[array.id, 1]
n_subject = design[array.id, 2]

write.csv(design[array.id,], paste("/mnt/EpsteinFSS/data/sbian/CauchyTest/Code/experiment_bash/result", array.id, ".csv", sep = ""))

print(paste("num pheno is", n_phenotype, "; num subject is", n_subject, sep = " "))