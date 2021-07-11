# fgsea包

- 参考资料
- [http://www.bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html](http://www.bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html)
- [https://bioconductor.org/packages/devel/bioc/manuals/fgsea/man/fgsea.pdf](https://bioconductor.org/packages/devel/bioc/manuals/fgsea/man/fgsea.pdf)

```r
rm(list = ls())
options(stringsAsFactors = F)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("fgsea")

library(fgsea)
library(data.table)
library(ggplot2)

data(examplePathways)
data(exampleRanks)
set.seed(42)
```

# fgsea

- 只需要两个文件：
    - 基因排序(Rank list)：把需分析的基因按照表达量或logfc排序
    - 基因集(Gene set): 预设好的一个Gene Set 通路文件，通常来自已进行了功能注释的数据或其他的实验结果。

## examplePathways 通路文件

- The list was obtained by selecting all the pathways from ‘reactome.db‘ package that contain mouse genes

```r
help("examplePathways") 
# 通路文件from 'reactome.db' package that contain mouse genes

length(examplePathways) 
# 1457 pathways
```

## exampleRank 基因数据

- The data were obtained by doing differential expression between Naive and Th1-activated states for GEO dataset GSE14308

```r
head(exampleRanks) 
# 170942    109711     18124     12775     72148     16010 
# -63.33703 -49.74779 -43.63878 -41.51889 -33.26039 -32.77626
class(exampleRanks)
# [1] "numeric"
length(exampleRanks)
# [1] 12000
```

## 富集分析

```r
# pathways：通路文件；
# stats：需要分析的基因+表达量或logfc等;
# nperm：Number of permutations to do. Minimial possible nominal p-value is about
#1/nperm;
# minSize：通路最小基因数;
# maxSize：通路最大基因数;
fgseaRes <- fgsea(pathways = examplePathways, 
                  stats    = exampleRanks,
                  minSize  = 15,
                  nperm=10000,
                  maxSize  = 500)
head(fgseaRes[order(pval), ])
# pathway         pval        padj        ES
# 1:                  5990980_Cell_Cycle 0.0001225490 0.002436671 0.5388497
# 2:         5990979_Cell_Cycle,_Mitotic 0.0001260716 0.002436671 0.5594755
# 3:    5991210_Signaling_by_Rho_GTPases 0.0001326084 0.002436671 0.4238512
# 4:                     5991454_M_Phase 0.0001390241 0.002436671 0.5576247
# 5: 5991023_Metabolism_of_carbohydrates 0.0001404297 0.002436671 0.4944766
# 6:        5991209_RHO_GTPase_Effectors 0.0001406074 0.002436671 0.5248796
# NES nMoreExtreme size                                leadingEdge
# 1: 2.681906            0  369   66336,66977,12442,107995,66442,19361,...
# 2: 2.746538            0  317   66336,66977,12442,107995,66442,12571,...
# 3: 2.011844            0  231 66336,66977,20430,104215,233406,107995,...
# 4: 2.557715            0  173   66336,66977,12442,107995,66442,52276,...
# 5: 2.243426            0  160    11676,21991,15366,58250,12505,20527,...
# 6: 2.376740            0  157 66336,66977,20430,104215,233406,107995,...
# adjPvalue
# 1: significant
# 2: significant
# 3: significant
# 4: significant
# 5: significant
# 6: significant

```

# 条形图

```r
# 删除前面的编号
fgseaRes$p <- substr(fgseaRes$pathway,9,100)
# 分组
fgseaRes$adjPvalue <- ifelse(fgseaRes$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
# 只展示前20
ggplot(fgseaRes[c(1:20),], aes(reorder(p, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=0.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways Enrichment Score from GSEA")
```

![fgsea%E5%8C%85%20ed8deb01f095484b95d1dd7d8ce64b24/Untitled.png](fgsea%E5%8C%85%20ed8deb01f095484b95d1dd7d8ce64b24/Untitled.png)

# plotEnrichment

- 绘制单独的GSEA富集图
- 基本参数

```r
plotEnrichment(pathway, stats, gseaParam = 1, ticksSize = 0.2)
```

- pathway：富集通路；
stats：基因矩阵；
gseaParam：GSEA parameter.
ticksSize：对应的垂直线宽度(default: 0.2)

```r
y <- head(fgseaRes[order(pval), ], 1)$pathway
y
# [1] "5990980_Cell_Cycle"

plotEnrichment(examplePathways[[y]],
               exampleRanks) + labs(title=y)
```

![fgsea%E5%8C%85%20ed8deb01f095484b95d1dd7d8ce64b24/Untitled%201.png](fgsea%E5%8C%85%20ed8deb01f095484b95d1dd7d8ce64b24/Untitled%201.png)

GSEA结果解读：

- 第一步我们需要根据基因的logFC对基因进行排序
- 研究的这个数据集中是否包含我们的目的基因，计算Enrich score的原则就是，从前到后依次检查基因是否是我们当前研究的数据集所包含的，如果包含就加一个正值，如果不包含就加一个负值
- 横坐标表示基因列表的数量
- 黑色的竖线代表的是我们的目的基因，已经被排好序，如果竖线聚集在头部，称为头部效应，如果在尾部，称为尾部效应
- GSEA也可以进行GO和KEGG分析，找到对应的数据集即可

此图：Enrichment score 高, 基因显著富集到左边，绿色的曲线出现高峰，而且，黑色竖线的分布也比较密集，说明该通路受到影响。NES（Normalized Enrichment Score）为正数。

# plotGseaTable

- 基本参数

```r
plotGseaTable(
pathways,
stats,
fgseaRes,
gseaParam = 1,
colwidths = c(5, 3, 0.8, 1.2, 1.2),
render = TRUE
)
```

- pathways：富集通路
stats：基因矩阵
fgseaRes：fgsea结果
gseaParam：Adjusts displayed statistic values, values closer to 0 flatten plots. Default = 1, value of 0.5 is a good choice too. 使gene rank线条更平滑
colwidth：列宽度，为0就不显示
render：If true, the plot is rendered to the current device. Otherwise, the grob is returned.
Default is true.

```r
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
topPathways[1:5]
#[1] "5990980_Cell_Cycle"               "5990979_Cell_Cycle,_Mitotic"     
#[3] "5991210_Signaling_by_Rho_GTPases" "5991454_M_Phase"                 
#[5] "5991209_RHO_GTPase_Effectors"

plotGseaTable(examplePathways[topPathways], exampleRanks, fgseaRes, 
              gseaParam = 0.5)
```

![fgsea%E5%8C%85%20ed8deb01f095484b95d1dd7d8ce64b24/Untitled%202.png](fgsea%E5%8C%85%20ed8deb01f095484b95d1dd7d8ce64b24/Untitled%202.png)

# 个性化分析

- 当有了自己的数据，就需要对感兴趣的基因集进行分析，首先进入GSEA官网进行下载 [https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp](https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp) 。基因集分为八个大类，分别由“H"和“C1”-“C7”开头，每个数据集都有详细的描述。

![fgsea%E5%8C%85%20ed8deb01f095484b95d1dd7d8ce64b24/Untitled%203.png](fgsea%E5%8C%85%20ed8deb01f095484b95d1dd7d8ce64b24/Untitled%203.png)

- 其中C6是关于肿瘤的数据

![fgsea%E5%8C%85%20ed8deb01f095484b95d1dd7d8ce64b24/Untitled%204.png](fgsea%E5%8C%85%20ed8deb01f095484b95d1dd7d8ce64b24/Untitled%204.png)

- c6.all.v7.4.entrez.gmt 通路名称，官网链接，富集的基因

![fgsea%E5%8C%85%20ed8deb01f095484b95d1dd7d8ce64b24/Untitled%205.png](fgsea%E5%8C%85%20ed8deb01f095484b95d1dd7d8ce64b24/Untitled%205.png)

- gmtPathways 就生成对应的pathways
- gene需要提前转换成entrez ID

```r
#fgsea
pathways<-gmtPathways("c6.all.v7.4.entrez.gmt")
fgseaRes <- fgsea(pathways, gene, nperm=1000,minSize=15, maxSize=500)
```

# reactomePathways通路

- 也可以使用reactomePathways通路进行GSEA分析

```r
# BiocManager::install("reactome.db")
library(reactome.db)
pathways <- reactomePathways(names(exampleRanks))
fgsea_reactome <- fgsea(pathways = pathways, 
                        stats = exampleRanks,
                        minSize=15,
                        maxSize=500,
                        nperm=100000)
```