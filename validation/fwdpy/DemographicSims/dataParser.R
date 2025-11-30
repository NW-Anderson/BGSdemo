library(data.table)
library(dplyr)
library(ggplot2)
setwd("~/Documents/GitHub/BGSdemo/validation/fwdpy/DemographicSims")
params <- fread("bottleneckParams.txt")
names(params) <- c("s","demo","seed")

setwd("/media/nathan/T7/BGSdemo/bottleneckData")
master <- data.table()
count <- 0
totalCount <- length(list.files())
done <- c()
for(file in list.files()){
  count <- count + 1
  if(count %% 100 == 0) print(paste(count, " of ", totalCount))
  par <- params %>% filter(paste(seed,".csv",sep="")==file)
  sfs <- fread(file) %>% as.matrix() %>% as.vector()
  numsum <- 1
  foo <- paste(par$s,par$demo)
  if(foo %in% done){
    # break
    master <- master %>% mutate(numberSummed = if_else(par$s == s & par$demo == demo,
                                             numberSummed + numsum,
                                             numberSummed),
                      afs = if_else(par$s == s & par$demo == demo,
                                    afs + sfs[obs+1],
                                    afs))
  }else{
    done <- c(done, foo)
    master <- dplyr::bind_rows(master,
                               data.table(s = par$s,
                                          demo = par$demo,
                                          seed = par$seed,
                                          numberSummed = numsum,
                                          obs = 0:(length(sfs)-1),
                                          afs = sfs))
  }
}
setwd("/media/nathan/T7/BGSdemo/")
# save(master, file = "bottleneck.RData")
load(file = "bottleneck.RData")
df <- master %>% mutate(afs = afs / numberSummed) %>% 
  select(-seed) 

df %>% mutate(theta = 4 * 1e4 * 1e-8 * 1e6) %>% 
  filter(obs == 2 * 1e4)

ggplot(df) + 
  geom_line(aes(x=obs, y = afs, color = as.factor(s))) + 
  facet_wrap(vars(demo), scales = "free") + 
  scale_y_log10()

setwd("/media/nathan/T7/BGSdemo/parsedbottleneckData")
for(curs in unique(df$s)){
  for(curdemo in unique(df$demo)){
    tmp <- df %>% filter(s == curs,
                         demo == curdemo)
    fwrite(list(tmp$afs),paste(curs,"_",curdemo,".csv",sep=""))
  }
}

setwd("~/Documents/GitHub/BGSdemo/validation/fwdpy/DemographicSims/ooa")
params <- fread("ooaSinglePop.txt")
names(params) <- c("s","demo", "seed")

setwd("/media/nathan/T7/BGSdemo/ooaSinglePopData")
master <- data.table()
count <- 0
totalCount <- length(list.files())
done <- c()
for(file in list.files()){
  count <- count + 1
  if(count %% 100 == 0) print(paste(count, " of ", totalCount))
  par <- params %>% filter(paste(seed,".csv",sep="")==file)
  sfs <- fread(file) %>% as.matrix() %>% as.vector()
  numsum <- 1
  foo <- paste(par$s,par$demo)
  if(foo %in% done){
    # break
    master <- master %>% mutate(numberSummed = if_else(par$s == s & par$demo == demo,
                                                       numberSummed + numsum,
                                                       numberSummed),
                                afs = if_else(par$s == s & par$demo == demo,
                                              afs + sfs[obs+1],
                                              afs))
  }else{
    done <- c(done, foo)
    master <- dplyr::bind_rows(master,
                               data.table(s = par$s,
                                          demo = par$demo,
                                          seed = par$seed,
                                          numberSummed = numsum,
                                          obs = 0:(length(sfs)-1),
                                          afs = sfs))
  }
}
setwd("/media/nathan/T7/BGSdemo/")
# save(master, file = "ooaSinglePop.RData")
load(file = "ooaSinglePop.RData")
df <- master %>% mutate(afs = afs / numberSummed) %>% 
  select(-seed) 

df %>% mutate(theta = 4 * 1e4 * 1e-8 * 1e6) %>% 
  filter(obs == 2 * 1e4)

ggplot(df) + 
  geom_line(aes(x=obs, y = afs, color = as.factor(s))) + 
  facet_wrap(vars(demo), scales = "free") + 
  scale_y_log10()

setwd("/media/nathan/T7/BGSdemo/parsedooaSinglePopData")
for(curs in unique(df$s)){
  for(curdemo in unique(df$demo)){
    tmp <- df %>% filter(s == curs,
                         demo == curdemo)
    fwrite(list(tmp$afs),paste(curs,"_",curdemo,".csv",sep=""))
  }
}

