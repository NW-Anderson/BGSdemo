library(dplyr)
master <- data.frame()
for(s in c(1e-3,5e-3,1e-2)){
  for(n in c(1e3,5e3,1e4)){
    master <- dplyr::bind_rows(master,
                               data.frame(s,
                                          n))
  }
}
master$tag <- 1:nrow(master)
write.table(master, file = "params.txt",
            row.names = F,col.names = F, quote=F)