create_chromevol_input<-function(nodes=tan, tree=tree){


  for(i in 1:length(tan)){
    sub<-extract.clade(tree, tan[i])
    sub_table<-table[ table$V1 %in% sub$tip.label,]
    sub_table

    write.csv(sub_table, "sub_table.csv")
    tab_chrom(sub, table="sub_table.csv",missing = "X", filename = paste0("subtree",i))

    ##The following file is in this link
    ##https://www.dropbox.com/s/ny1jvkw3diufoyp/PIP_control?dl=0

    files<-readLines("PIP_control")

    write(
      gsub("guide_tree", paste0("subtree",i,"_input.tree"), files[1]), paste0("Control_","subtree",i))
    write(
      gsub("counts", paste0("subtree",i,"_input_table.txt"), files[2]), paste0("Control_","subtree",i),append = TRUE)
    write(
      files[3], paste0("Control_","subtree",i),append = TRUE)
    write(
      gsub("PIP_out", paste0("subtree",i,"_out"), files[4]), paste0("Control_","subtree",i),append = TRUE)
    write(
      files[5], paste0("Control_","subtree",i),append = TRUE)
    write(
      files[6], paste0("Control_","subtree",i),append = TRUE)

    print(i)


  }
  cat("Check your working directory")

}
