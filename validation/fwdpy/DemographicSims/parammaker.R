library(dplyr)
master <- data.frame()
for(s in c(1e-3,5e-3,1e-2)){
  for(demo in c("1k.yaml", "5k.yaml")){
    for(rep in 1:1000){
      master <- dplyr::bind_rows(master,
                                 data.frame(s,
                                            demo))
    }
  }
}
master$seed <- sample(1:1e5,nrow(master))
write.table(master, file = "bottleneckParams.txt",
            row.names = F,col.names = F, quote=F)


library(dplyr)
master <- data.frame()
for(s in c(1e-3,5e-3,1e-2)){
  for(demo in c("ooaSinglePop.yaml")){
    for(rep in 1:3000){
      master <- dplyr::bind_rows(master,
                                 data.frame(s,
                                            demo))
    }
  }
}
master$seed <- sample(1:1e5,nrow(master))
write.table(master, file = "ooaSinglePop.txt",
            row.names = F,col.names = F, quote=F)

library(dplyr)
master <- data.frame()
for(s in c(1e-3,5e-3,1e-2)){
  for(demo in c("ooaTwoPop.yaml")){
    for(rep in 1:3000){
      master <- dplyr::bind_rows(master,
                                 data.frame(s,
                                            demo))
    }
  }
}
master$seed <- sample(1:1e5,nrow(master))
write.table(master, file = "ooaTwoPop.txt",
            row.names = F,col.names = F, quote=F)

library(dplyr)
master <- data.frame()
for(s in c(1e-3,5e-3,1e-2)){
  for(demo in c("ooa.yaml")){
    for(rep in 1:3000){
      master <- dplyr::bind_rows(master,
                                 data.frame(s,
                                            demo))
    }
  }
}
master$seed <- sample(1:1e5,nrow(master))
write.table(master, file = "ooa.txt",
            row.names = F,col.names = F, quote=F)