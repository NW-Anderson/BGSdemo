library(data.table)
library(dplyr)
library(ggplot2)
setwd("~/Documents/GitHub/BGSdemo/validation/fwdpy")
params <- fread("params.txt")
names(params) <- c("s","N","seed")

setwd("/media/nathan/T7/BGSdemo/genMats")
master <- data.table()
for(file in list.files()){
  par <- params %>% filter(paste(seed,".csv",sep="")==file)
  df <- fread(file) %>% as.matrix()
  sfs <- rowSums(df) %>% table()
  counts <- names(sfs) %>% as.integer()
  obs <- unname(sfs) %>% as.vector()
  frq <- obs / 2 / par$N
  master <- dplyr::bind_rows(master, data.table(s = par$s,
                                                N = par$N,
                                                seed = par$seed,
                                                counts = counts,
                                                obs = obs,
                                                frq =frq))
}

ggplot(master) + 
  geom_line(aes(x = counts, y = obs)) + 
  facet_grid(rows = vars(s),
             cols = vars(N),
             scales = "free")
# this is the deleterious site SFS, need to drop a neutral mutation down