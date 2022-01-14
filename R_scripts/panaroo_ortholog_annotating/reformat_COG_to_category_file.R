rm(list = ls(all.names = TRUE))

COG_descrip <- read.table("/data1/gdouglas/db/COG_definitions/cog-20.def.tab",
                          header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "")

COG_to_category <- data.frame(matrix(NA, ncol = 2, nrow = 0))

i = 1
for (row_i in 1:nrow(COG_descrip)) {
 
  COG_id <- COG_descrip[row_i, 1]
  
  COG_categories <- strsplit(COG_descrip[row_i, 2], "")[[1]]
  
  for (COG_category in COG_categories) {
   
    COG_to_category[i, ] <- c(COG_id, COG_category)
  
    i <- i + 1 
  }
}

COG_to_category <- COG_to_category[order(COG_to_category$X2), ]

write.table(x = COG_to_category, file = "/data1/gdouglas/db/COG_definitions/cog-20.to_category.tsv",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
