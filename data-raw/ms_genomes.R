x <- "C:/Users/hemst/Documents/ss_file_access/theta4k_1000_10_rho40k.txt"

# genome information
chrl <- 10000000 # chromosome length CURRENTLY MUST BE THE SAME ACROSS ALL CHRs

# read in data
x <- process_ms(x, chrl)
ms_genomes <- x
names(x) <- c("genotypes", "snp_meta")

usethis::use_data(ms_genomes, overwrite = TRUE)
