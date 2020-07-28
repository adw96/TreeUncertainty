## A first attempt to pull out the terrapene gene tree

write.table(x = NULL, file = "processed_terrapene_trees_all.txt", quote = F, sep = ";",row.names = F, col.names = F)
gene_tree_files <- list.files("Terrapene_raw", pattern = "^T", full.names = F, include.dirs = F)
system("mkdir no_last")
for (my_file in gene_tree_files) {
  system(paste("sed '$ d' Terrapene_raw/",  my_file, " > no_last/", my_file, sep = ""))
  system(paste("sed -e '1,31d' no_last/",  my_file, "> Terrapene_raw/",  my_file, sep = ""))
  
  trees_data <- read.table(paste("Terrapene_raw/", my_file, sep = ""))
  trees_trees <- trees_data[6] #just the trees
  trees_trees <- trees_trees[,1]
  trees_trees <- as.character(trees_trees)
  trees_trees <- gsub("\\[[^\\]]*\\]", "", trees_trees, perl=TRUE)
  trees_trees
  write.table(x = trees_trees, file = "processed_terrapene_trees_all.txt", quote = F, sep = ";",row.names = F, col.names = F, append = T)
}
