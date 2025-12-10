library(parallel)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("No arguments supplied. Please provide the taxa file.", call. = FALSE)
}

input_taxa_file <- args[1]

threads <- args[2]

tax_data <- read.table(input_taxa_file, header=FALSE, sep="\t")

colnames(tax_data) <- c("Feature.ID", "Taxon")

rank <- c( "k__" = "kingdom", "p__" = "phylum", "c__" = "class", "o__" = "order", "f__" = "family", "g__" = "genus", "s__" = "species")

tax_dummies <- 3

taxs <- as.data.frame(do.call("rbind", mclapply(tax_data$Taxon, function(aa){
  bb <- unlist(strsplit(aa, split = "\\; "))[1:(7+tax_dummies)]
  bb
}, mc.cores = threads, mc.cleanup = TRUE)), stringsAsFactors = FALSE)

colnames(taxs) <- rank

taxs$header <- tax_data$Feature.ID

rownames(taxs) <- taxs$header

taxs$effective_taxon <- mclapply(1:nrow(taxs), function(aa){
  bb <- taxs[ aa, 1:(7+tax_dummies) ]
  unlist(bb[max(which( !( grepl("__$", bb, perl = TRUE) | is.na(bb) ) ))])
},  mc.cores = threads, mc.cleanup = TRUE)

uniq_tax <- apply(taxs[, 1:(7+tax_dummies)], 2, function(aa){
    bb <- aa[ !( is.na(aa) ) ]
    sort(unique(unlist(bb)))
})

floop <- 1

taxid_proto <- lapply(names(uniq_tax)[ sapply( uniq_tax, function(aa){ !identical( aa, character(0)) })],
    function(aa){
        if( aa == "kingdom"){ floop <<- 1 }
        bb <- uniq_tax[ aa ][[1]]
        cc <- (floop+1):(floop+length(bb))
        floop <<- floop + length(cc)
        list( "name_txt" = bb, "taxid" = cc)
    }
)

names(taxid_proto) <- rank

identify <- function(aa){
    bb_rank <- rank[ gsub("__.*", "__", aa )]
    if( aa %in% taxid_proto[[bb_rank]][[1]] ){
      taxid_proto[[bb_rank]][[2]][ match(aa, taxid_proto[[bb_rank]][[1]] ) ]
    }else{
      print(paste0("taxon ", aa, " wasn't in the taxid map"))
    }
}

identify_parent <- function(aa){
  if( grepl("k__", aa)){
    "1"
  } else {
    bb_rank <- rank[ gsub("__.*", "__", aa) ]
    cc <- taxs[ match(aa, taxs[ , bb_rank]) , (which( names(taxs) == bb_rank) -1) ]
    identify(cc)
    }
}

names_proto <- data.frame(
  "tax_id" = unlist( mclapply(
    unlist( uniq_tax),
    identify, mc.cores = threads, mc.cleanup = TRUE )),
  "name_txt" = unlist( uniq_tax),
  "unique name" = unlist( uniq_tax),
  "name class" = "scientific name",
  stringsAsFactors = FALSE
)

nodes_proto <- data.frame(
  "tax_id" = unlist(mclapply( names_proto$name_txt, identify, mc.cores = threads, mc.cleanup = TRUE )),
  "parent tax_id" = unlist( mclapply( names_proto$name_txt, identify_parent, mc.cores = threads, mc.cleanup = TRUE )),
  "rank" = rank[ unlist( mclapply( names_proto$name_txt, function(aa){ gsub("(^\\w__).*", "\\1", aa) }, mc.cores = threads, mc.cleanup = TRUE ))],
  "embl code" = "",
  "division id" = 0,
  "inherited div flag" = 1,
  "genetic code id" = 11,
  "inherited GC" = 1,
  "mitochondrial genetic code" = 0,
  "inherited MGC" = 1,
  "GenBank hidden" = 0,
  "hidden subtree root" = 0,
  "comments" = ""
)

nodes_proto <- rbind(
  c( "1","1","no rank","","8","0","1","0","0","0","0","0",""),
  nodes_proto
)

names_proto <- rbind(
  c( "1", "all", "all", "synonym"),
  c( "1", "root", "root", "scientific name"),
  names_proto
)

tax_id <- unlist(mclapply(taxs[, "effective_taxon"], identify, mc.cores = threads, mc.cleanup = TRUE))

ref_id <- taxs$header

rename_str <- paste0(ref_id,"|","kraken:taxid","|",tax_id)

kraken_data <- as.data.frame(cbind(ref_id, rename_str))

write.table(nodes_proto, file="nodes.dmp", sep = '\t|\t', col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(names_proto, file="names.dmp", sep = '\t|\t', col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(kraken_data, "seq_id_replacement.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
