# 公共入口可被项目 workflow 调用，也可直接作为 CLI 运行。
main <- function(cli_args = commandArgs(trailingOnly = TRUE)) {
#!/usr/bin/env Rscript

# SCENIC R 独立分析入口
# 输入：主 Seurat 流程输出的注释后 RDS；输出共表达网络、regulon、AUCell 和 RSS。

# 功能：读取 --name value 形式的命令行参数。
arg_value <- function(name, default = NULL) {
  args <- cli_args
  index <- tail(which(args == name), 1L)
  if (!length(index) || index == length(args)) default else args[[index + 1L]]
}

# 功能：解析逗号分隔的 cisTarget ranking 数据库文件名。
split_values <- function(x) trimws(strsplit(x, ",", fixed = TRUE)[[1]])

# 功能：把 pheatmap 同时保存为 PDF 和 600 dpi PNG，并使用 NPG 连续色阶。
save_heatmap <- function(matrix, stem, width = 10, height = 8) {
  colors <- grDevices::colorRampPalette(c("#FFF5F0", "#F39B7F", "#E64B35", "#7A0019"))(101)
  grDevices::pdf(paste0(stem, ".pdf"), width = width, height = height)
  pheatmap::pheatmap(matrix, color = colors, border_color = NA)
  grDevices::dev.off()
  grDevices::png(paste0(stem, ".png"), width = width, height = height, units = "in", res = 600)
  pheatmap::pheatmap(matrix, color = colors, border_color = NA)
  grDevices::dev.off()
}

# 功能：从 Seurat v5 对象提取指定表达层；gene-summarized counts 是官方教程的首选输入。
seurat_expression <- function(obj, assay, layer) {
  obj <- SeuratObject::JoinLayers(obj, assay = assay)
  as.matrix(SeuratObject::LayerData(obj, assay = assay, layer = layer))
}

# 功能：执行 gene filtering、GENIE3、cisTarget pruning、AUCell 和 RSS。
run_scenic <- function(
  input_rds,
  output_dir,
  organism,
  database_dir,
  databases,
  assay = "RNA",
  layer = "counts",
  celltype_col = "cell_type",
  n_cores = 8L,
  min_cell_fraction = 0.01
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  obj <- readRDS(input_rds)
  expr_mat <- seurat_expression(obj, assay, layer)

  old_dir <- getwd()
  on.exit(setwd(old_dir), add = TRUE)
  setwd(output_dir)
  scenic_options <- SCENIC::initializeScenic(
    org = organism,
    dbDir = normalizePath(database_dir, winslash = "/", mustWork = TRUE),
    dbs = databases,
    datasetTitle = "SCENIC",
    nCores = n_cores
  )

  # 官方流程先过滤极低表达基因，减少随机森林网络推断的噪声与计算量。
  genes_kept <- SCENIC::geneFiltering(
    expr_mat,
    scenicOptions = scenic_options,
    minCountsPerGene = 3 * min_cell_fraction * ncol(expr_mat),
    minSamples = min_cell_fraction * ncol(expr_mat)
  )
  expr_filtered <- expr_mat[genes_kept, , drop = FALSE]
  write.csv(
    data.frame(gene = genes_kept),
    file.path(output_dir, "genes_used_for_scenic.csv"),
    row.names = FALSE
  )

  set.seed(1001)
  SCENIC::runCorrelation(expr_filtered, scenic_options)
  SCENIC::runGenie3(expr_filtered, scenic_options)
  SCENIC::runSCENIC_1_coexNetwork2modules(scenic_options)
  SCENIC::runSCENIC_2_createRegulons(scenic_options)
  SCENIC::runSCENIC_3_scoreCells(scenic_options, expr_filtered)

  regulon_auc <- SCENIC::loadInt(scenic_options, "aucell_regulonAUC")
  auc <- t(AUCell::getAUC(regulon_auc))
  write.csv(
    data.frame(cell = rownames(auc), auc, check.names = FALSE),
    file.path(output_dir, "regulon_auc.csv"),
    row.names = FALSE
  )

  # RSS 衡量 regulon 对某一注释群的特异性，不是组间差异显著性检验。
  if (celltype_col %in% colnames(obj[[]])) {
    labels <- as.character(obj[[celltype_col, drop = TRUE]][colnames(AUCell::getAUC(regulon_auc))])
    rss <- SCENIC::calcRSS(AUC = AUCell::getAUC(regulon_auc), cellAnnotation = labels)
    write.csv(
      data.frame(regulon = rownames(rss), rss, check.names = FALSE),
      file.path(output_dir, "regulon_specificity_scores.csv"),
      row.names = FALSE
    )
    variable_rss <- rss[order(matrixStats::rowVars(as.matrix(rss)), decreasing = TRUE), , drop = FALSE]
    save_heatmap(head(variable_rss, 30), file.path(output_dir, "top_regulon_RSS_heatmap"))
  }

  # 把 regulon AUC 作为独立 assay 写回，不污染 RNA 表达层。
  auc_for_seurat <- t(auc[colnames(obj), , drop = FALSE])
  obj[["SCENIC"]] <- SeuratObject::CreateAssayObject(data = auc_for_seurat)
  saveRDS(obj, file.path(object_dir, "scenic_annotated.rds"))

  writeLines(
    c(
      "# SCENIC 结果说明",
      "",
      paste0("- 输入 RDS：", normalizePath(input_rds, winslash = "/", mustWork = FALSE)),
      paste0("- organism：", organism, "；assay/layer：", assay, "/", layer, "。"),
      paste0("- cisTarget 数据库：", paste(databases, collapse = ", "), "。"),
      paste0("- 基因过滤：至少在 ", min_cell_fraction * 100, "% 细胞中表达。"),
      "- GENIE3 建立 TF-target 共表达候选边；RcisTarget 用 motif 富集剪枝；AUCell 计算单细胞 regulon 活性。",
      "- regulon_specificity_scores.csv 用于识别细胞类型特异 regulon。",
      "- SCENIC 是计算推断；motif 支持和共表达不能替代 ChIP、扰动或其他直接实验验证。"
    ),
    file.path(output_dir, "README_results_CN.md")
  )
  invisible(list(object = obj, auc = auc, scenic_options = scenic_options))
}

input_rds <- arg_value("--input")
output_dir <- arg_value("--output")
object_dir <- arg_value("--object-dir", file.path("data/processed", "scenic"))
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
database_dir <- arg_value("--database-dir")
database_files <- arg_value("--databases")
if (any(vapply(list(input_rds, output_dir, database_dir, database_files), is.null, logical(1)))) {
  stop("用法：Rscript scenic_workflow.R --input annotated.rds --output results/scenic --organism hgnc --database-dir cisTarget --databases db1.feather,db2.feather")
}

run_scenic(
  input_rds = input_rds,
  output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE),
  organism = arg_value("--organism", "hgnc"),
  database_dir = database_dir,
  databases = split_values(database_files),
  assay = arg_value("--assay", "RNA"),
  layer = arg_value("--layer", "counts"),
  celltype_col = arg_value("--celltype-col", "cell_type"),
  n_cores = as.integer(arg_value("--cores", "8")),
  min_cell_fraction = as.numeric(arg_value("--min-cell-fraction", "0.01"))
)

}
if (sys.nframe() == 0L) main()
