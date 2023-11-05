library(DESeq2)
library(pheatmap)
library(clusterProfiler)
library(org.Mm.eg.db)
library(stringr)
library(ggplot2)
library(ggpubr)

counts = read.csv(
    'gene.counts', 
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
    sample = c("sample1", "sample2", "sample3", "sample1", "sample2", "sample3", "sample1", "sample2", "sample3")
)

# 按照sampleID更改samples的行名
rownames(samples) = samples$sampleID

# 将因子型数据的默认排序设置为1=sample1, 2=sample2, 3=sample3
samples$sample = factor(samples$sample, levels = c('sample1', 'sample2', 'sample3'))

# 将数据转为matrix格式
counts = as.matrix(counts[rownames(samples)])

# 从矩阵中抽取DESeq可以运行的数据形式
dds = DESeqDataSetFromMatrix(countData = counts, colData = samples, design = ~sample)

# 使用DESeq做差异分析
dds = DESeq(dds)

# vst方差稳定转换
vsd = vst(dds, blind = F)

# 做PCA
# 使用sample来分组，并返回坐标轴用于下一步作图
plotPCA(vsd, intgroup = c('sample'), returnData = TRUE)

PCAreturn <- "
              PC1        PC2   group  sample   name
C0_1   -16.586637   7.444676 sample1 sample1   C0_1
C0_2   -12.125149 -18.965823 sample2 sample2   C0_2
C0_3   -26.368675  11.700192 sample3 sample3   C0_3
C50_1   26.720069   1.819072 sample1 sample1  C50_1
C50_2  -16.740632  -1.364539 sample2 sample2  C50_2
C50_3    8.106638 -13.586068 sample3 sample3  C50_3
C500_1  16.105152  28.532094 sample1 sample1 C500_1
C500_2   3.091405  -3.977518 sample2 sample2 C500_2
C500_3  17.797829 -11.602086 sample3 sample3 C500_3
"

# 将PCAreturn转换为表格形式
PCAreturn <- read.table(header = TRUE, text = PCAreturn)

# 将sample列改为因子型并且排序
PCAreturn$sample = factor(PCAreturn$sample, levels = c('sample1', 'sample2', 'sample3'))

# ggplot绘图
ggplot(PCAreturn)+
    geom_point(aes(PC1, PC2, color=sample), size = 2)+
    labs(x = "PC1", y = "PC2", face = "bold")+
    theme_light(base_size = 16)


#####################################################################

counts_fpkm = read.csv(
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
    counts_fpkm[col_fpkm] = (counts_fpkm[clm] * 10^6) / (counts_fpkm$Length * as.numeric(total))
}