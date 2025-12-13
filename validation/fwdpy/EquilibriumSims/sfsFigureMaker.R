library(data.table)
library(dplyr)
library(ggplot2)
setwd("~/Documents/GitHub/BGSdemo/validation/fwdpy/EquilibriumSims/humanMaps")
params <- fread("humanmap.txt")
names(params) <- c("s","N", "recMap", "mutMap", "exonMap", "seed")

# setwd("/media/nathan/T7/BGSdemo/humanEquilData")
# master <- data.table()
# count <- 0
# totalCount <- length(list.files())
# done <- c()
# for(file in list.files()){
#   count <- count + 1
#   if(count %% 100 == 0) print(paste(count, " of ", totalCount))
#   par <- params %>% filter(paste(seed,".csv",sep="")==file)
#   df <- fread(file) %>% as.matrix()
#   numsum <- nrow(df)
#   sfs <- colSums(df)
#   foo <- paste(par$s,par$N)
#   if(foo %in% done){
#     # break
#     master <- master %>% mutate(numberSummed = if_else(par$s == s & par$N == N,
#                                              numberSummed + numsum,
#                                              numberSummed),
#                       afs = if_else(par$s == s & par$N == N,
#                                     afs + sfs[obs+1],
#                                     afs))
#   }else{
#     done <- c(done, foo)
#     master <- dplyr::bind_rows(master,
#                                data.table(s = par$s,
#                                           N = par$N,
#                                           seed = par$seed,
#                                           numberSummed = numsum,
#                                           obs = 0:(length(sfs)-1),
#                                           afs = sfs))
#   }
# }
setwd("/media/nathan/T7/BGSdemo/")
# save(master, file = "equilHuman.RData")
load(file = "equilHuman.RData")
df <- master %>% mutate(afs = afs / numberSummed) %>% 
  select(-seed) 
# %>%  filter(obs != 0 & obs != 2 * N)
df %>% mutate(theta = 4 * N * 1e-8 * 1e6) %>% 
  filter(obs == 2 * N)

ggplot(df) + 
  geom_line(aes(x=obs, y = afs, color = as.factor(s))) + 
  facet_wrap(vars(N), scales = "free") + 
  scale_y_log10()

setwd("/media/nathan/T7/BGSdemo/parsedequilHuman")
for(curs in unique(df$s)){
  for(curN in unique(df$N)){
    tmp <- df %>% filter(s == curs,
                         N == curN)
    fwrite(list(tmp$afs),paste(curs,"_",curN,".csv",sep=""))
  }
}

setwd("~/Documents/GitHub/BGSdemo/validation/fwdpy/EquilibriumSims")
params <- fread("morereps.txt")
names(params) <- c("s","N","seed")

# setwd("/media/nathan/T7/BGSdemo/afsData")
# master <- data.table()
# count <- 0
# totalCount <- length(list.files())
# done <- c()
# for(file in list.files()){
#   count <- count + 1
#   if(count %% 100 == 0) print(paste(count, " of ", totalCount))
#   par <- params %>% filter(paste(seed,".csv",sep="")==file)
#   df <- fread(file) %>% as.matrix()
#   numsum <- nrow(df)
#   sfs <- colSums(df)
#   foo <- paste(par$s,par$N)
#   if(foo %in% done){
#     # break
#     master <- master %>% mutate(numberSummed = if_else(par$s == s & par$N == N,
#                                              numberSummed + numsum,
#                                              numberSummed),
#                       afs = if_else(par$s == s & par$N == N,
#                                     afs + sfs[obs+1],
#                                     afs))
#   }else{
#     done <- c(done, foo)
#     master <- dplyr::bind_rows(master,
#                                data.table(s = par$s,
#                                           N = par$N,
#                                           seed = par$seed,
#                                           numberSummed = numsum,
#                                           obs = 0:(length(sfs)-1),
#                                           afs = sfs))
#   }
# }
setwd("/media/nathan/T7/BGSdemo/")
# save(master, file = "afs.RData")
load(file = "afs.RData")
df <- master %>% mutate(afs = afs / numberSummed) %>% 
  select(-seed) 
# %>%  filter(obs != 0 & obs != 2 * N)
df %>% mutate(theta = 4 * N * 1e-8 * 1e6) %>% 
  filter(obs == 2 * N)

ggplot(df) + 
  geom_line(aes(x=obs, y = afs, color = as.factor(s))) + 
  facet_wrap(vars(N), scales = "free") + 
  scale_y_log10()

setwd("/media/nathan/T7/BGSdemo/equilAFS")
for(curs in unique(df$s)){
  for(curN in unique(df$N)){
    tmp <- df %>% filter(s == curs,
                         N == curN)
    fwrite(list(tmp$afs),paste(curs,"_",curN,".csv",sep=""))
  }
}

setwd("/media/nathan/T7/BGSdemo/moreSelData")
master <- data.table()
count <- 0
totalCount <- length(list.files())
done <- c()
for(file in list.files()){
  curSeed = strsplit(file, "_")[[1]][1] %>% as.integer()
  win = strsplit(file, "_")[[1]][2]
  count <- count + 1
  if(count %% 100 == 0) print(paste(count, " of ", totalCount))
  par <- params %>% filter(seed==curSeed)
  df <- fread(file) %>% as.matrix() %>% as.vector() %>% table()
  obsVec <- unname(df) %>% as.vector()
  frqVec <- names(df) %>% as.numeric()
  tmp <- data.table(obs2 = obsVec,
                    frq = frqVec)
  df <- data.table(frq = 0:(2 * par$N) / 2 / par$N,
                   obs = 0) %>%
    merge(tmp, by = c("frq"), all = T) %>%
    mutate(obs2 = if_else(is.na(obs2),
                          0,
                          obs2)) %>%
    mutate(obs = obs + obs2) %>%
    select(-obs2) %>%
    mutate(N = par$N,
           s = par$s,
           wi = win)
  foo <- paste(par$s,par$N,win)
  if(foo %in% done){
    # break
    master <- merge(master, df, by = c("frq", "s", "N", "wi"), all = T) %>%
      mutate(obs.y = if_else(is.na(obs.y),
                             0,
                             obs.y)) %>%
      mutate(obs = obs.x + obs.y) %>%
      select(-c(obs.x, obs.y))
  }else{
    done <- c(done, foo)
    master <- dplyr::bind_rows(master,
                               df)
  }
}

ggplot(master) + 
  geom_line(aes(x = frq,
                y = obs,
                color = as.factor(s))) + 
  facet_grid(rows = vars(N),
             cols = vars(wi)) + 
  scale_y_log10()

setwd("/media/nathan/T7/BGSdemo/parsedMoreSelData")
for(curs in unique(master$s)){
  for(curN in unique(master$N)){
    for(curwi in unique(master$wi)){
      tmp <- master %>% filter(s == curs,
                               N == curN,
                               wi == curwi)
      fwrite(list(tmp$obs),paste(curs,"_",curN,"_",curwi,".csv",sep=""))
    }
  }
}


# setwd("/media/nathan/T7/BGSdemo/equilSelData")
# master <- data.table()
# count <- 0
# totalCount <- length(list.files())
# for(file in list.files()){
#   count <- count + 1
#   print(paste(count, " of ", totalCount))
#   par <- params %>% filter(paste(seed,".csv",sep="")==file)
#   df <- fread(file) %>% as.matrix() %>% as.vector()
#   master <- dplyr::bind_rows(master,data.table(s = par$s,
#                                                N = par$N,
#                                                seed = par$seed,
#                                                frq = df))
# }
setwd("/media/nathan/T7/BGSdemo/")
# save(master, file = "equilSelData.RData")
load(file = "equilSelData.RData")

dens <- master %>% filter(frq != 0,
                          frq != 1) %>%
  group_by(N,s) %>% 
  # summarize(n = n())
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

setwd("/media/nathan/T7/BGSdemo/equilSelData")
master <- data.table()
count <- 0
totalCount <- length(list.files())
done <- c()
for(file in list.files()){
  count <- count + 1
  if(count %% 100 == 0) print(paste(count, " of ", totalCount))
  par <- params %>% filter(paste(seed,".csv",sep="")==file)
  df <- fread(file) %>% as.matrix() %>% as.vector() %>% table()
  obsVec <- unname(df) %>% as.vector()
  frqVec <- names(df) %>% as.numeric()
  tmp <- data.table(obs2 = obsVec,
                    frq = frqVec)
  df <- data.table(frq = 0:(2 * par$N) / 2 / par$N,
                   obs = 0) %>%
    merge(tmp, by = c("frq"), all = T) %>%
    mutate(obs2 = if_else(is.na(obs2),
                          0,
                          obs2)) %>%
    mutate(obs = obs + obs2) %>%
    select(-obs2) %>%
    mutate(N = par$N,
           s = par$s)
  foo <- paste(par$s,par$N)
  if(foo %in% done){
    # break
    master <- merge(master, df, by = c("frq", "s", "N"), all = T) %>%
      mutate(obs.y = if_else(is.na(obs.y),
                             0,
                             obs.y)) %>%
      mutate(obs = obs.x + obs.y) %>%
      select(-c(obs.x, obs.y))
  }else{
    done <- c(done, foo)
    master <- dplyr::bind_rows(master,
                               df)
  }
}

ggplot(master) + 
  geom_line(aes(x = frq,
                y = obs,
                color = as.factor(s))) + 
  facet_wrap(vars(N)) + 
  scale_y_log10()

setwd("/media/nathan/T7/BGSdemo/parsedEquilSelData")
for(curs in unique(master$s)){
  for(curN in unique(master$N)){
    tmp <- master %>% filter(s == curs,
                         N == curN)
    fwrite(list(tmp$obs),paste(curs,"_",curN,".csv",sep=""))
  }
}

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