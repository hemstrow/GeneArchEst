# names
sn <- read.table("data/stat_names.txt")
sn <- sn[,1]
sn <- as.character(sn)
names_descriptive_stats <- sn
usethis::use_data(names_descriptive_stats, overwrite = TRUE)


number_descriptive_stats <- length(names_descriptive_stats)

# save
usethis::use_data(number_descriptive_stats, overwrite = TRUE, internal = T)
