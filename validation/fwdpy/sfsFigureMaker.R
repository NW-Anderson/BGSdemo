library(data.table)
library(dplyr)
library(ggplot2)
setwd("~/Documents/GitHub/BGSdemo/validation/fwdpy")
params <- fread("morereps.txt")
names(params) <- c("s","N","seed")

setwd("/media/nathan/T7/BGSdemo/afsData")
master <- data.table()
count <- 0
totalCount <- length(list.files())
for(file in list.files()){
  count <- count + 1
  if(count %% 100 == 0) print(paste(count, " of ", totalCount))
  par <- params %>% filter(paste(seed,".csv",sep="")==file)
  df <- fread(file) %>% as.matrix()
  numsum <- nrow(df)
  sfs <- colSums(df)
  if(par$s %in% master$s && par$N %in% master$N){
    # break
    master <- master %>% mutate(numberSummed = if_else(par$s == s & par$N == N, 
                                             numberSummed + numsum,
                                             numberSummed),
                      afs = if_else(par$s == s & par$N == N,
                                    afs + sfs[obs+1],
                                    afs))
  }else{
    master <- dplyr::bind_rows(master,
                               data.table(s = par$s,
                                          N = par$N,
                                          seed = par$seed,
                                          numberSummed = numsum,
                                          obs = 0:(length(sfs)-1),
                                          afs = sfs))
  }
}
setwd("/media/nathan/T7/BGSdemo/")
save(master, file = "afs.RData")

# setwd("/media/nathan/T7/BGSdemo/morerepsData")
# master <- data.table()
# count <- 0
# totalCount <- length(list.files())
# for(file in list.files()){
#   count <- count + 1
#   print(paste(count, " of ", totalCount))
#   par <- params %>% filter(paste(seed,".csv",sep="")==file)
#   df <- fread(file) %>% as.matrix()
#   df <- df[,1]
#   master <- dplyr::bind_rows(master,data.table(s = par$s,
#                                                N = par$N,
#                                                seed = par$seed,
#                                                frq = df))
# }
setwd("/media/nathan/T7/BGSdemo/")
# save(master, file = "moreReps.RData")
load(file = "moreReps.RData")

# ggplot(master) + 
#   geom_density(aes(frq, color = as.factor(s))) + 
#   facet_wrap(vars(N)) +
#   scale_y_log10() +
#   scale_x_log10()

dens <- master %>% group_by(N,s) %>% 
  reframe(densx = density(frq, 
                          from =1/(2*unique(N)), 
                          to = 1 - 1/(2*unique(N)))$x,
          densy = density(frq, 
                          from =1/(2*unique(N)), 
                          to = 1 - 1/(2*unique(N)))$y,
          mx = max(densy)) %>%
  mutate(densy = densy / mx) %>% 
  ungroup()

ggplot(dens) + 
  geom_line(aes(x = densx, 
                y = densy, 
                color = as.factor(s))) +
  facet_wrap(vars(N)) +
  scale_y_log10() 





setwd("~/Documents/GitHub/BGSdemo/validation/fwdpy")
params <- fread("params.txt")
names(params) <- c("s","N","seed")


setwd("/media/nathan/T7/BGSdemo/EquilData")
master <- data.table()
for(file in list.files()){
  par <- params %>% filter(paste(seed,".csv",sep="")==file)
  df <- fread(file) %>% as.matrix()
  df <- df[,1]
  master <- dplyr::bind_rows(master,data.table(s = par$s,
                                               N = par$N,
                                               seed = par$seed,
                                               frq = df))
}

# ggplot(master) + 
#   geom_density(aes(frq, color = as.factor(s))) + 
#   facet_wrap(vars(N)) +
#   scale_y_log10() +
#   scale_x_log10()

dens <- master %>% group_by(N,s) %>% 
  reframe(densx = density(frq, 
                          from =1/(2*unique(N)), 
                          to = 1 - 1/(2*unique(N)))$x,
          densy = density(frq, 
                          from =1/(2*unique(N)), 
                          to = 1 - 1/(2*unique(N)))$y,
          mx = max(densy)) %>%
  mutate(densy = densy / mx) %>% 
  ungroup()

ggplot(dens) + 
  geom_line(aes(x = densx, 
                y = densy, 
                color = as.factor(s))) +
  facet_wrap(vars(N)) +
  scale_y_log10() 

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