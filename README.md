# Description

RNA-seq 数据分析所使用的代码示例。

# RNA-seq

RNA-seq 是研究转录组应用最广泛，也是最重要的技术之一，RNA-seq 分析内容包括序列比对、转录本拼接、表达定量、差异分析、融合基因检测、可变剪接、RNA 编辑和突变检测等，具体流程和常用工具如下图所示，通常的分析不一定需要走完全部流程，按需进行，某些步骤可以跳过、简化等。

![](img/rna-seq.png)

RNA-seq 中最常用的分析方法就是找出差异表达基因 (Differential gene expression, DEG)，在实验室中，标准流程就分为三步：

- step1: 构建测序文库，包括提取 RNA, 富集 mRNA 或清除核糖体 RNA, 合成 cDNA, 加上接头
- step2: 在高通量测序平台（通常为 Illumina） 上对文库进行测序，每个样本的测序深度为 10-30M 读长
- step3: 数据分析，具体而言：对测序得到的读长进行比对或组装到转录本上；对覆盖到每个基因区域的读长进行计数；根据统计模型鉴定不同样本间差异表达的基因，（这种分析过程是比较传统的方法）

# Environment

- 上游分析：成都理工大学超算平台 Red Hat 4.8.5-36
- 下游分析：Debian GNU/Linux 12 (bookworm) on Windows 10 x86_64

conda 搭建分析环境：

```bash
git clone https://github.com/yuanj82/RNA-seq.git
conda env create --file env.yml
```

# Read more

访问我的博客：https://yuanj.top 可以获得更多相关内容。

# Special Thanks

欢迎补充代码与文档等，如果代码与文档存在问题欢迎随时补充。
