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