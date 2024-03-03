# 数据示例：
# Geneid	Chr	Start	End	Strand	Length	C0_1	C0_2	C0_3	C50_1	C50_2	C50_3	C500_1	C500_2	C500_3
# KYUSt_chr1.5-E1	1	53340	53891	+	552	233	146	195	201	189	264	177	326	243
# KYUSt_chr1.6-E1	1	54648	55664	-	1017	2	4	8	12	6	3	20	47	17
# KYUSt_chr1.7-E1	1	59936	60283	+	348	3	0	6	7	0	0	1	0	0
# KYUSt_chr1.8-E1	1	61782	61841	+	60	0	0	0	0	0	0	2	0	0
# KYUSt_chr1.8-E2	1	62373	62621	+	249	9	21	81	34	17	8	57	12	9
# KYUSt_chr1.9-E10	1	66745	66819	+	75	0	0	0	0	0	0	0	0	0
# KYUSt_chr1.9-E11	1	67360	67443	+	84	0	0	1	0	0	0	0	0	0
# KYUSt_chr1.9-E12	1	69776	69819	+	44	0	0	0	0	0	0	0	0	0
# KYUSt_chr1.9-E13	1	71271	71760	+	490	0	0	0	0	0	0	0	0	0
# KYUSt_chr1.107-E3	1	677255	677949	-	695	4115	26070	15855	3440	11431	5642	6410	5576	11130
# KYUSt_chr1.107-E2	1	678069	679025	-	957	3339	11707	5568	4545	6474	6084	6897	4574	4191
# KYUSt_chr1.107-E1	1	679121	679679	-	559	3645	16809	6030	3822	9043	7565	6074	5887	2603

setwd('./RNA-seq/downstream')
rm(list = ls())  
Sys.setenv(LANGUAGE = "en")

library(DESeq2)

counts = read.csv(
    'gene2.counts', 
    header = T,  
    sep = '\t', 
    row.names = "Geneid", 
    comment.char = '#', 
    check.name = F
)

# 删掉前五列
counts = counts[,-c(1:5)]

# 保留行相加大于10的数据
counts = counts[rowSums(counts)>10, ]

samples = data.frame(
    sampleID = c("C0_1", "C0_2", "C0_3", "C50_1", "C50_2", "C50_3", "C500_1", "C500_2", "C500_3"), 
    sample = c("sample1", "sample1", "sample1", "sample2", "sample2", "sample2", "sample3", "sample3", "sample3")
)

# 按照sampleID更改samples的行名
rownames(samples) = samples$sampleID

# 将因子型数据的默认排序设置为1=sample1, 2=sample2, 3=sample3
samples$sample = factor(samples$sample, levels = c('sample1', 'sample2', 'sample3'))

# 构建 DESeqDataSet 对象
dds = DESeqDataSetFromMatrix(countData = counts, colData = samples, design = ~sample)

# 计算差异倍数并获得 p 值
# parallel = TRUE 可以多线程运行，在数据量较大时建议开启
dds_count <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)

# 查看样本上调还是下调
sampl1_vs_sample2 <- results(dds_count, contrast = c('sample', 'sample1', 'sample2'))

sampl1_vs_sample3 <- results(dds_count, contrast = c('sample', 'sample2', 'sample3'))

result1 <- data.frame(sampl1_vs_sample2, stringsAsFactors = FALSE, check.names = FALSE)
result2 <- data.frame(sampl1_vs_sample3, stringsAsFactors = FALSE, check.names = FALSE)

write.table(result1, 'sampl1_vs_sample2.DESeq2.txt', col.names = NA, sep = '\t', quote = FALSE)
write.table(result2, 'sampl1_vs_sample3.DESeq2.txt', col.names = NA, sep = '\t', quote = FALSE)