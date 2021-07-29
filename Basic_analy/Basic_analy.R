####------------------加载需要用的程序包-------------------------
p_list = c("FactoMineR", "dplyr", "factoextra", "ggpubr", "pca3d")
for(p in p_list){
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p)
}
# 安装github来源R包
suppressWarnings(suppressMessages(library(devtools)))
if (!requireNamespace("ggbiplot", quietly = TRUE))
  install_github("vqv/ggbiplot")
library("dplyr")
# 加载漂亮热图绘制包pheatmap
library(pheatmap)
library("factoextra")
zoom=1.5 # 控制图片缩放比例
###-------------------数据读取及预处理-------------------------
# 读取数据all-表示所有特征，共7345，filter表示筛选后的特征388个
raw_count <- read.table('genusdata-sort-no-nan.csv',header = TRUE,sep = ',')
count <- raw_count[,-2]
rownames(count) <- count[,1]
count <- count[,-1]
count <- t(count)

####-----------------------------绘制热图------------------------------------
# 读取数据表并标准化百分比
otu=as.data.frame(t(t(count)/colSums(count)*100))
# 通常我们只关注高丰度且显著差异的，按每个OTU的总丰度排序
otu$sum <- rowSums(otu)
# 按每行的和降序排列
otu_order <- otu[order(otu$sum, decreasing=TRUE), ]
# 取丰度前30的OTUs
mat <- otu_order[1:20, -ncol(otu)]
#两种转换
scale_test <- apply(mat, 2, scale)
rownames(scale_test) <- rownames(mat)
# 读取元数据，添加样本注释
annot_data <- data.frame(row.names=colnames(mat), group=raw_count$label)
# 绘制添加样本列注释的图
p1 <- pheatmap(mat=scale_test, treeheight_col=5,
               border_color='grey60',
               cluster_row=F,  cluster_col=F,
               show_rownames=T, show_colnames=F,annotation_col=annot_data,
               clustering_method="complete",cutree_cols=5,fontsize=7)
ggsave(paste0("P1.pdf"), p1, width=89*zoom, height=56*zoom, units="mm")
###----------------------------基于OTU表PCA分析--------------------------------
# 如果变量之间的数据的处于不同数量级或者变量之间的均值/方差相差很大时，建议是进行标准化，常见为：scale(count, center =T, scale. =T)对数据进行标准化，FactoMineR包的PCA()函数和基础包的prcomp()函数自带函数自带标准化参数，如下：
# 基于OTU计算PCA，执行数据标准化
otu.pca <- prcomp(t(count), scale. = TRUE)
get_eigenvalue(otu.pca)[1:3,]

##---------贡献度图-------1-----------
(p2 <- fviz_eig(otu.pca, addlabels = TRUE))
(var <- get_pca_var(otu.pca))

##-------主成分图-单独主成分----------4--------
# 转换原始数据为百分比
norm <- t(t(count)/colSums(count,na=T)) * 10
# 筛选mad(median absolute deviation,中位数偏差的绝对值的中位数,衡量特异波动的方法)值大于0.5的ASV,
#mad.5 <- norm[apply(norm,1,mad)>0.5,]
# 另一种方法：按mad值排序取前N个，如6个波动最大的ASVs
mad.5 <- head(norm[order(apply(norm,1,mad), decreasing=T),],n=5)
# 计算PCA和菌与菌轴的相关性
otu.pca <- prcomp(t(mad.5))
# 绘制观测值PCA主成分分析图，外层()可对保存的图形同时预览
(p3 <- fviz_pca_biplot(otu.pca, col.ind = raw_count$label, palette = "jco", addEllipses = TRUE, label = "var",
                      col.var = "black", repel = TRUE, legend.title = "Group"))
p3 <- p3+theme_bw()+theme(panel.grid = element_blank())
## -----------相关性组合图---------------5----
# 加载ggbiplot并绘制观测值PCA主成分分析图
suppressWarnings(suppressMessages(library("ggbiplot")))
# 绘制观测值PCA图
p<-ggbiplot(otu.pca, obs.scale = 1, var.scale = 1, groups = raw_count$label, ellipse = TRUE,var.axes = T)
ggsave(paste0("sandaintu.pdf"), p, width=89*zoom, height=56*zoom, units="mm")

# 需要制作PCA结果可视化的绘图数据文件
PC1 <- otu.pca$x[,1]
PC2 <- otu.pca$x[,2]
otu.pca.data <- data.frame(rownames(otu.pca$x),PC1,PC2,raw_count$label)
colnames(otu.pca.data) <- c("sample", "PC1", "PC2", "group")
# 这里需要ggpubr包，对组间进行统计检验以及组合图的拼接
library("ggpubr")
# 设置比较组
my_comparisons = list(c('A','B'),c('A','C'),c('A','D'),c('A','E'),
                      c('B','C'),c('B','D'),c('B','E'),
                      c('C','D'),c('C','E'),
                      c('D','E'))

# 绘制y轴为PC1值的分组箱线图
p4 <- ggplot(otu.pca.data, aes(x = group, y= PC1, colour = group)) +
  geom_boxplot() +
  theme(panel.background = element_rect(fill = "white",colour = "white"),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        legend.position="none") +
  xlab("") + ylab("") +
  stat_compare_means(comparisons = my_comparisons,label = "p.signif") # 添加显著性检验
# 绘制y轴为PC2值的分组箱线图
p5 <- ggplot(otu.pca.data, aes(x = group, y= PC2, colour = group)) +
  geom_boxplot(aes()) +
  theme(panel.background = element_rect(fill = "white",colour = "white"),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        legend.position="none") +
  xlab("") + ylab("") +
  coord_flip() + # coord_flip()函数翻转坐标轴
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")
# ggpubr::ggarrange()函数对图进行拼接
(p <- ggarrange(p4, NULL, p3, p5, widths = c(10,4), heights = c(4,8), align = "hv"))
ggsave(paste0("zuhetu.pdf"), p, width=89*zoom, height=56*zoom, units="mm")
ggsave(paste0("p4.pdf"), p4, width=89*zoom, height=56*zoom, units="mm")
ggsave(paste0("p5.pdf"), p5, width=89*zoom, height=56*zoom, units="mm")
