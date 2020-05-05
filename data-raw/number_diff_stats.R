# names
sn <- read.table("data/diff_names.txt")
sn <- sn[,1]
sn <- as.character(sn)
names_diff_stats <- sn
usethis::use_data(names_diff_stats, overwrite = TRUE)


number_diff_stats <- length(names_diff_stats)

# save
usethis::use_data(number_diff_stats, overwrite = TRUE, internal = T)
