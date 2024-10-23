#Download gene ontology annotation 
library(data.table)
# Load GAF file
gaf_file <- "goa_human.gaf"  # Ensure you've downloaded the file to this path
gaf <- fread(gaf_file, skip = "!")  # Skip header lines starting with '!'
# Select the relevant columns
gaf_filtered <- gaf[, .(GO_ID = V5, Gene = V3, Aspect = V9)]

# Split the data into different aspects
biological_process <- gaf_filtered[Aspect == "P"]
molecular_function <- gaf_filtered[Aspect == "F"]
cellular_component <- gaf_filtered[Aspect == "C"]


# Install required packages if not already installed
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GSEABase")

library(GSEABase)

# Load the GO annotations
go_obo <- getOBOCollection("go-basic.obo")

GO_non_parsed <- go_obo@.kv

GO_non_parsed_subset <- GO_non_parsed[1:2000, ]

#write GO_non_parsed_subset to a file
write.table(GO_non_parsed,"GO_non_parsed.txt", sep = "\t", row.names = FALSE, col.names = FALSE)

# filter the GO_non_parsed_subset for name
GO_numbers_with_names <- GO_non_parsed[grepl("name", GO_non_parsed$key), ]

# move the values in key column into separate columns
GO_numbers_with_names2 <- GO_numbers_with_names |>
  tidyr::pivot_wider(names_from = key, values_from = value)


#remove first raw
GO_numbers_with_names2 <- GO_numbers_with_names2[-1, ] 
#remove second column
GO_numbers_with_names2 <- GO_numbers_with_names2[, -2]

# View the result
print(GO_numbers_with_names2)

# change column names to "id", "name", "ontology"
GO_numbers_with_names2 <- setNames(GO_numbers_with_names2, c("id", "name", "ontology"))

#write GO_numbers_with_names2 to a file
write.table(GO_numbers_with_names2,"GO_numbers_with_names2.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

#change the colnames of gaf_filtered to match the GO_numbers_with_names2
colnames(gaf_filtered) <- c("id", "gene", "ontology")

# join the gaf_filtered and GO_numbers_with_names2 on id
GO_names_and_genes <- dplyr::left_join(gaf_filtered, GO_numbers_with_names2, by = "id")

#save GO_names_and_genes to a file
write.table(GO_names_and_genes,"GO_names_and_genes.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

