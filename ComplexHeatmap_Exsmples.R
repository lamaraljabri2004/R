rm(list = ls() )  # cleanup the environment
library(ComplexHeatmap)
library(circlize)

expr = readRDS(system.file(package = "ComplexHeatmap", "extdata", "gene_expression.rds"))
mat = as.matrix(expr[, grep("cell", colnames(expr))])
base_mean = rowMeans(mat)
mat_scaled = t(apply(mat, 1, scale))
?apply
type = gsub("s\\d+_", "", colnames(mat))
ha = HeatmapAnnotation(type = type, annotation_name_side = "left")

ht_list = Heatmap(mat_scaled, name = "expression", row_km = 5, 
                  col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
                  top_annotation = ha, 
                  show_column_names = FALSE, row_title = NULL, show_row_dend = FALSE) +
  Heatmap(base_mean, name = "base mean", 
          top_annotation = HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2:6), 
                                                                    height = unit(2, "cm"))),
          width = unit(15, "mm")) +
  rowAnnotation(length = anno_points(expr$length, pch = 16, size = unit(1, "mm"), 
                                     axis_param = list(at = c(0, 2e5, 4e5, 6e5), 
                                                       labels = c("0kb", "200kb", "400kb", "600kb")),
                                     width = unit(2, "cm"))) +
  Heatmap(expr$type, name = "gene type", 
          top_annotation = HeatmapAnnotation(summary = anno_summary(height = unit(2, "cm"))),
          width = unit(15, "mm"))

ht_list = rowAnnotation(block = anno_block(gp = gpar(fill = 2:6, col = NA)), 
                        width = unit(2, "mm")) + ht_list

draw(ht_list)

#Correlations between methylation, expression and other genomic features
res_list = readRDS("data/meth.rds")
type = res_list$type
type=data.frame(type)
mat_meth = res_list$mat_meth
mat_expr = res_list$mat_expr
direction = res_list$direction
cor_pvalue = res_list$cor_pvalue
gene_type = res_list$gene_type
anno_gene = res_list$anno_gene
dist = res_list$dist
anno_enhancer = res_list$anno_enhancer
column_tree = hclust(dist(t(mat_meth)))
column_order = column_tree$order
library(RColorBrewer)
meth_col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
direction_col = c("hyper" = "red", "hypo" = "blue")
expr_col_fun = colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
pvalue_col_fun = colorRamp2(c(0, 2, 4), c("white", "white", "red"))
gene_type_col = structure(brewer.pal(length(unique(gene_type)), "Set3"), 
                          names = unique(gene_type))
anno_gene_col = structure(brewer.pal(length(unique(anno_gene)), "Set1"), 
                          names = unique(anno_gene))
dist_col_fun = colorRamp2(c(0, 10000), c("black", "white"))
enhancer_col_fun = colorRamp2(c(0, 1), c("white", "orange"))


ht_opt(
  legend_title_gp = gpar(fontsize = 8, fontface = "bold"), 
  legend_labels_gp = gpar(fontsize = 8), 
  heatmap_column_names_gp = gpar(fontsize = 8),
  heatmap_column_title_gp = gpar(fontsize = 10),
  heatmap_row_title_gp = gpar(fontsize = 8)
)

ha = HeatmapAnnotation(type = type, 
                       col = list(type = c("Tumor" = "pink", "Control" = "royalblue")),
                       annotation_name_side = "left")
ha2 = HeatmapAnnotation(type = type, 
                        col = list(type = c("Tumor" = "pink", "Control" = "royalblue")), 
                        show_legend = FALSE)

ht_list = Heatmap(mat_meth, name = "methylation", col = meth_col_fun,
                  column_order= column_order,
                  top_annotation = ha, column_title = "Methylation") +
  Heatmap(direction, name = "direction", col = direction_col) +
  Heatmap(mat_expr[, column_tree$order], name = "expression", 
          col = expr_col_fun, 
          column_order = column_order, 
          top_annotation = ha2, column_title = "Expression") +
  Heatmap(cor_pvalue, name = "-log10(cor_p)", col = pvalue_col_fun) +
  Heatmap(gene_type, name = "gene type", col = gene_type_col) +
  Heatmap(anno_gene, name = "anno_gene", col = anno_gene_col) +
  Heatmap(dist, name = "dist_tss", col = dist_col_fun) +
  Heatmap(anno_enhancer, name = "anno_enhancer", col = enhancer_col_fun, 
          cluster_columns = FALSE, column_title = "Enhancer")

draw(ht_list, row_km = 2, row_split = direction,
     column_title = "Comprehensive correspondence between methylation, expression and other genomic features", 
     column_title_gp = gpar(fontsize = 12, fontface = "bold"), 
     merge_legends = TRUE, heatmap_legend_side = "bottom")


#Visualize Methylation Profile with Complex Annotations

library(matrixStats)
library(GenomicRanges)






#Add multiple boxplots for single row
m1 = matrix(sort(rnorm(100)), 10, byrow = TRUE)
m2 = matrix(sort(rnorm(100), decreasing = TRUE), 10, byrow = TRUE)

ht_list = Heatmap(m1, name = "m1") + Heatmap(m2, name = "m2")
m1.df=data.frame(m1)
m2.df=data.frame(m2)

rg = range(c(m1, m2))
rg.df=data.frame(rg)
rg[1] = rg[1] - (rg[2] - rg[1])* 0.02
rg[2] = rg[2] + (rg[2] - rg[1])* 0.02
anno_multiple_boxplot = function(index) {
  nr = length(index)
  pushViewport(viewport(xscale = rg, yscale = c(0.5, nr + 0.5)))
  for(i in seq_along(index)) {
    grid.rect(y = nr-i+1, height = 1, default.units = "native")
    grid.boxplot(m1[ index[i], ], pos = nr-i+1 + 0.2, box_width = 0.3, 
                 gp = gpar(fill = "red"), direction = "horizontal")
    grid.boxplot(m2[ index[i], ], pos = nr-i+1 - 0.2, box_width = 0.3, 
                 gp = gpar(fill = "green"), direction = "horizontal")
  }
  grid.xaxis()
  popViewport()
}

ht_list = ht_list + rowAnnotation(boxplot = anno_multiple_boxplot, width = unit(4, "cm"), 
                                  show_annotation_name = FALSE)
lgd = Legend(labels = c("m1", "m2"), title = "boxplots",
             legend_gp = gpar(fill = c("red", "green")))
draw(ht_list, padding = unit(c(20, 2, 2, 2), "mm"), heatmap_legend_list = list(lgd))






#How to visualize complex heatmaps interactively

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("InteractiveComplexHeatmap")

library(devtools)
install_github("jokergoo/InteractiveComplexHeatmap")


library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
m = matrix(rnorm(100*100), nrow = 100)
ht = Heatmap(m)
ht = draw(ht) # not necessary, but recommended

htShiny(ht)
