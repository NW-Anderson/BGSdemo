library(data.table)
library(dplyr)
library(ggplot2)

setwd("~/Documents/GitHub/BGSdemo/validation/fwdpy/EquilibriumSims/humanMaps")
exonMap <- fread("exons_chr22.bed")
names(exonMap) <- c("chr","start","stop")
mutMap <- fread("roulette_tbl_chr22.csv")
exonMutMap <- fread("exonMutMap.csv")
names(exonMutMap) <- c("start", "stop", "rate")

ggplot(exonMap) + 
  geom_point(aes(x = midPoint / 1e6,
                 y = size)) +
  geom_vline(xintercept = 23181599/1e6)
  geom_vline(xintercept = 49885518.61638361/1e6)

exonMap <- exonMap %>% mutate(size = stop - start,
                              midPoint = (start + stop) / 2)
exonMap %>% filter(size == max(exonMap$size))

ggplot(exonMap) + 
  geom_point(aes(x = midPoint,
                 y = size)) +
  geom_vline(xintercept = 49885518.61638361) + 
  coord_cartesian(xlim = c(49885518.61638361 + 1e6/2, 49885518.61638361 - 1e6/2))

tmp <- exonMutMap %>% filter(start >= 46256559 & stop <= 46263323) 
tmp <- tmp %>% mutate(bigU = (stop - start) * rate)
sum(tmp$bigU)

tmp <- exonMutMap %>% filter(start >= 46259941 - 1e6/2 & stop <= 46259941 + 1e6/2) 
tmp <- tmp %>% mutate(bigU = (stop - start) * rate)
sum(tmp$bigU)

exonMap %>% filter(size > 4e3)

ggplot(exonMap) + 
  geom_point(aes(x = midPoint,
                 y = size)) +
  geom_vline(xintercept = 21467417) + 
  coord_cartesian(xlim = c(21467417 + 1e6/2, 21467417 - 1e6/2))

tmp <- exonMutMap %>% filter(start >= 21467417 - 1e6/2 & stop <= 21467417 + 1e6/2) 
tmp <- tmp %>% mutate(bigU = (stop - start) * rate)
sum(tmp$bigU)

mxi = 0
mxbigU = 0
for(i in seq(1e4, 50819000,1e4)){
  tmp <- exonMutMap %>% filter(start >= i - 1e6/2 & stop <= i + 1e6/2) 
  tmp <- tmp %>% mutate(bigU = (stop - start) * rate)
  bigU <- sum(tmp$bigU)
  if(bigU > mxbigU){
    mxbigU = bigU
    mxi = i
  }
}

ggplot(exonMap) + 
  geom_point(aes(x = midPoint,
                 y = size)) +
  geom_vline(xintercept = mxi) + 
  coord_cartesian(xlim = c(mxi + 1e6/2, mxi - 1e6/2))

tmp <- exonMutMap %>% filter(start >= i - 1e6/2 & stop <= i + 1e6/2) 
tmp <- tmp %>% mutate(bigU = (stop - start) * rate)
bigU <- sum(tmp$bigU)

mxi = data.table()
for(i in 0:102 * 5e5){
  df <- data.table(pos = 0,
                   bigU = 0)
  for(j in seq(i-5e5,i,1e4)){
    tmp <- exonMutMap %>% filter(start >= j - 1e6/2 & stop <= j + 1e6/2) 
    tmp <- tmp %>% mutate(bigU = (stop - start) * rate)
    bigU <- sum(tmp$bigU)
    if(bigU > df$bigU){
      df <- data.table(pos = j,
                       bigU = bigU)
    }
  }
  # if(df$pos == 0 & df$bigU == 0) break
  mxi <- bind_rows(mxi, df)
}
mxi <- mxi %>% filter(pos != 0)
pos <- mxi$pos
fwrite(list(pos),"focalPos.csv")
