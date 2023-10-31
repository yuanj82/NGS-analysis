library(pheatmap)

counts = read.csv(
    'gene.counts', 
    header = T,  
    sep = '\t', 
    row.names = "Geneid", 
    comment.char = '#', 
    check.name = F
)

for (clm in colnames(counts_fpkm)[6:14]) {
    col_fpkm = paste0(clm, "FPKM")
    total = sum(counts_fpkm[clm])
    counts_fpkm[col_fpkm] = (counts_fpkm[clm] * 10^6) / (counts_fpkm$Length * as.numberic(total))
}