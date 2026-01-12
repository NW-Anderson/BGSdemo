library(dplyr)
master <- data.frame()
for(s in c(1e-3,5e-3,1e-2)){
  for(n in c(1e3,5e3,1e4)){
    for(rep in 1:100){
      master <- dplyr::bind_rows(master,
                                 data.frame(s,
                                            n))
    }
  }
}
master$seed <- sample(1:1e4,nrow(master))
write.table(master, file = "params.txt",
            row.names = F,col.names = F, quote=F)

library(dplyr)
master <- data.frame()
for(s in c(1e-3,5e-3,1e-2)){
  for(n in c(1e3,5e3,1e4)){
    for(rep in 1:1000){
      master <- dplyr::bind_rows(master,
                                 data.frame(s,
                                            n))
    }
  }
}
master$seed <- sample(1:1e5,nrow(master))
write.table(master, file = "morereps.txt",
            row.names = F,col.names = F, quote=F)

library(dplyr)
master <- data.frame()
for(s in c(1e-3,5e-3,1e-2)){
  for(n in c(1e3,5e3,1e4)){
    for(rep in 1:1000){
      master <- dplyr::bind_rows(master,
                                 data.frame(s,
                                            n,
                                            "YRI_recombination_map_hg38_chr_22.bed",
                                            "roulette_tbl_chr22.csv",
                                            "exons_chr22.bed"))
    }
  }
}
master$seed <- sample(1:1e5,nrow(master))
write.table(master, file = "humanmap.txt",
            row.names = F,col.names = F, quote=F)

library(dplyr)
master <- data.frame()
for(p1 in c(0.25,0.5,0.75)){
  for(n in c(1e3,5e3,1e4)){
    for(rep in 1:1000){
      p2 <- 1 - p1
      master <- dplyr::bind_rows(master,
                                 data.frame(p1,
                                            p2,
                                            n))
    }
  }
}
master$seed <- sample(1:1e4,nrow(master))
write.table(master, file = "dfeparams.txt",
            row.names = F,col.names = F, quote=F)