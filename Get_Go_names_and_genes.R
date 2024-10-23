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

#extracting the GO numbers with names https://stackoverflow.com/questions/78532080/getting-namespace-from-a-obo-formatted-text-file

# dat <- read.delim("go-basic.obo", header=F)

# #Getting a data frame format from non-repeated terms
# 
# 
#  GO_numbers_with_names <-  do.call(rbind, 
#           apply(cbind(grep("\\[Term]", dat$V1) + 1, 
#                       grep("\\[Term]", dat$V1) - 1 + 
#                         diff(c(grep("\\[Term]", dat$V1), nrow(dat)))), 1, \(x){
#                           res <- data.table::tstrsplit(dat[x[1]:x[2],], ": ")
#                           id <- grep("^id|name|namespace|def|is_obso", res[[1]])
#                           `length<-`(as.list(res[[2]][id]), 5)})) |> 
#     data.frame() |> 
#     setNames(c("id", "name", "namespace", "definition", "is_obsolete"))



