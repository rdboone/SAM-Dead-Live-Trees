treetbl <- read.csv("TreeTable.csv")
head(treetbl)
dim(treetbl)
sta.num <- 1:dim(treetbl)[1]
for(i in 1:dim(treetbl)[1]){
  if(treetbl$status[i] == "LIVING"){
    sta.num[i] <- 1
  }else if(treetbl$status[i] == "DEAD"){
    sta.num[i] <- 2
  }
}

unique(sta.num)

treetbl <- cbind(treetbl, "statusCode" = sta.num)
head(treetbl)

write.csv(treetbl, "TreeTable2.csv")
