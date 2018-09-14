#####

merge_counts(file.name, folders) <- function {
  
  temp= c()
  file.list.path = c()
  file.list = c()
  for (i in seq_along(folders)){
    temp = list.files(folders[i], pattern=".genes.results.withGeneName")
    file.list = c(file.list, temp)
    file.list.path = c(file.list.path, paste(folders[i], temp, sep="/"))
  }
  
  for(i in seq_along(file.list)){
    file.data = read.table(file.list.path[i], sep="\t", head=T, quote="", colClasses="character")
    file.sub = data.frame(file.data$gene_id, round(as.numeric(file.data$expected_count), 0))
    
    sample.id = gsub(".genes.results.withGeneName", "", file.list[i])
    
    names(file.sub) = c("gene_id", paste(sample.id, "_counts", sep=""))
    
    if(i == 1){
      file.df = file.sub
    }else{
      file.df = merge(file.df, file.sub, by.x="gene_id", by.y="gene_id")
    }
  }
  
  write.table(file.df, paste(out.dir, file.name, sep=""), sep="\t", row=F, quote=F)
  
}
