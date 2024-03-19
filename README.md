# Description

组学数据分析代码示例。

# Environment

- 上游分析：成都理工大学超算平台 Red Hat 4.8.5-36
- 下游分析：Debian GNU/Linux 12 (bookworm) on Windows 10 x86_64

conda 环境搭建：

```bash
git clone https://github.com/yuanj82/NGS-analysis.git
cd NGS-analysis
cp .condarc ~/
conda env create --file env.yml
```

# Repository files

```bash
➜  NGS-analysis git:(main) ✗ tree -a -I '.git' ./ 
./
├── .condarc
├── README.md
├── RNA-seq
│   ├── downstream
│   │   ├── DESeq2.r
│   │   ├── PCA.r
│   │   ├── bubble.r
│   │   ├── evolutionary.tree.r
│   │   ├── fpkm.r
│   │   ├── hs.cluster
│   │   │   ├── go.r
│   │   │   ├── kegg.r
│   │   │   └── symlbol_yo_entrzeid.r
│   │   ├── pheatmap.r
│   │   └── volcano.r
│   └── upstream
│       ├── old.sh
│       ├── single-step
│       │   ├── clean.sh
│       │   ├── example_acc_list.txt
│       │   ├── hista.sh
│       │   ├── samtools.sh
│       │   └── subread.sh
│       └── upstream.sh
├── env.yml
└── py_tools
    └── data_extraction.py

6 directories, 21 files
```

# Read more

- [转录组学分析基础——测序技术](https://yuanj.top/posts/h6g6c2i2/)
- [专题：RNA-Seq 上游分析学习](https://yuanj.top/posts/z6q5z7s5/)
- [RNA Seq 上游分析实践](https://yuanj.top/posts/j9t9v3y1/)
- [详解 ATAC-Seq](https://yuanj.top/posts/x5v6u9k6/)

# Special Thanks

欢迎补充代码与文档等，如果代码与文档存在问题欢迎随时补充。
