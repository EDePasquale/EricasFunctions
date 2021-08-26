###########################
#                         #
# Test function calls for #
#    Ericas_Functions     #
#                         #
###########################

qc_and_filter(directories=c("/Volumes/broad_vangalenlab/depasquale/CellRanger/210628_GeneCounts",
                "/Volumes/broad_vangalenlab/depasquale/CellRanger/210422_GeneCounts",
                "/Volumes/broad_vangalenlab/depasquale/CellRanger/210423_GeneCounts"), genome="GRCh38", suffix="QC2")

qc_and_filter(c("/Volumes/broad_vangalenlab/depasquale/CellRanger/210628_GeneCounts"), "GRCh38", "QC3")

qc_and_filter(directories=c("/Volumes/broad_vangalenlab/depasquale/CellRanger/210813_GeneCounts"), min.umis=1000, min.genes=500)

###########################

seurat_transfer_data(reference="~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_0.5_ccrem/TCELL_SUB/Seurat_int_0.5_TCELLSPLIT.rds",
                     datasets=c("~/Documents/Projects/Data_files_temp/Donor/BM191119/Seurat.rds",
                                "~/Documents/Projects/Data_files_temp/Donor/BM191227/Seurat.rds"),
                     results_name="results_test1",
                     ordering_vector=c(12,8,3,7,5,14,13,16,11,10,2,15,9,1,17,4,6),
                     cluster_names=c("CD4NaiveT", "B", "MidEryth", "CD8MemoryT", "Mono", "NK","LateEryth", "EarlyEryth", "CD4MemoryT", "PreB",
                                     "ProB", "HSCProg", "cDC", "ncMono", "Plasma",  "pDC", "CD8NaiveT"),
                     cluster_colors=c("#C0007C", "#2a4bd7", "#f08080", "#e9debb", "#ffee33", "#BE8A66", "#ad2323", "#ff9233", "#ffcdf3", "#9dafff",
                                      "#8126c0", "#000000", "#1d6914", "#29d0d0", "#575757",  "#81c57a", "#708090"))

seurat_transfer_data(reference="~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_0.5_ccrem/TCELL_SUB/Seurat_int_0.5_TCELLSPLIT.rds",
                     datasets=c("~/Documents/Projects/Data_files_temp/Donor/BM191119/Seurat.rds",
                                "~/Documents/Projects/Data_files_temp/Donor/BM191227/Seurat.rds"),
                     results_name="results_test2")


seurat_transfer_data(reference="~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_0.5_ccrem/TCELL_SUB/Seurat_int_0.5_5BM.rds",
                     datasets=c("~/Documents/Projects/Data_files_temp/Donor/BPDCN712/Seurat.rds",
                                "~/Documents/Projects/Data_files_temp/Donor/BPDCN712R/Seurat.rds",
                                "~/Documents/Projects/Data_files_temp/Donor/BPDCN180329/Seurat.rds",
                                "~/Documents/Projects/Data_files_temp/Donor/BPDCN181128/Seurat.rds",
                                "~/Documents/Projects/Data_files_temp/Donor/BPDCN190711/Seurat.rds",
                                "~/Documents/Projects/Data_files_temp/Donor/PT1-DX/Seurat.rds",
                                "~/Documents/Projects/Data_files_temp/Donor/PT1-MRD/Seurat.rds",
                                "~/Documents/Projects/Data_files_temp/Donor/PT5-DX/Seurat.rds",
                                "~/Documents/Projects/Data_files_temp/Donor/PT12-DX/Seurat.rds",
                                "~/Documents/Projects/Data_files_temp/Donor/PT12-REL/Seurat.rds"),
                     results_name="210824_SeuratTD_5BM_10BPDCN",
                     ordering_vector=c(12,8,3,7,5,14,13,16,11,10,2,15,9,1,17,4,6),
                     cluster_names=c("CD4NaiveT", "B", "MidEryth", "CD8MemoryT", "Mono", "NK","LateEryth", "EarlyEryth", "CD4MemoryT", "PreB",
                                     "ProB", "HSCProg", "cDC", "ncMono", "Plasma",  "pDC", "CD8NaiveT"),
                     cluster_colors=c("#C0007C", "#2a4bd7", "#f08080", "#e9debb", "#ffee33", "#BE8A66", "#ad2323", "#ff9233", "#ffcdf3", "#9dafff",
                                      "#8126c0", "#000000", "#1d6914", "#29d0d0", "#575757",  "#81c57a", "#708090"))

###############################

random_forest(reference="~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_0.5_ccrem/TCELL_SUB/Seurat_int_0.5_TCELLSPLIT.rds",
                     datasets=c("~/Documents/Projects/Data_files_temp/Donor/BM191119/Seurat.rds",
                                "~/Documents/Projects/Data_files_temp/Donor/BM191227/Seurat.rds"),
                     results_name="results_test_rf1",
                     ordering_vector=c(12,8,3,7,5,14,13,16,11,10,2,15,9,1,17,4,6),
                     cluster_names=c("CD4NaiveT", "B", "MidEryth", "CD8MemoryT", "Mono", "NK","LateEryth", "EarlyEryth", "CD4MemoryT", "PreB",
                                     "ProB", "HSCProg", "cDC", "ncMono", "Plasma",  "pDC", "CD8NaiveT"),
                     cluster_colors=c("#C0007C", "#2a4bd7", "#f08080", "#e9debb", "#ffee33", "#BE8A66", "#ad2323", "#ff9233", "#ffcdf3", "#9dafff",
                                      "#8126c0", "#000000", "#1d6914", "#29d0d0", "#575757",  "#81c57a", "#708090"))

random_forest(reference="~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_0.5_ccrem/TCELL_SUB/Seurat_int_0.5_TCELLSPLIT.rds",
                     datasets=c("~/Documents/Projects/Data_files_temp/Donor/BM191119/Seurat.rds",
                                "~/Documents/Projects/Data_files_temp/Donor/BM191227/Seurat.rds"),
                     results_name="results_test_rf2")


random_forest(reference="~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_0.5_ccrem/TCELL_SUB/Seurat_int_0.5_5BM.rds",
                     datasets=c("~/Documents/Projects/Data_files_temp/Donor/BPDCN712/Seurat.rds",
                                "~/Documents/Projects/Data_files_temp/Donor/BPDCN712R/Seurat.rds",
                                "~/Documents/Projects/Data_files_temp/Donor/BPDCN180329/Seurat.rds",
                                "~/Documents/Projects/Data_files_temp/Donor/BPDCN181128/Seurat.rds",
                                "~/Documents/Projects/Data_files_temp/Donor/BPDCN190711/Seurat.rds",
                                "~/Documents/Projects/Data_files_temp/Donor/PT1-DX/Seurat.rds",
                                "~/Documents/Projects/Data_files_temp/Donor/PT1-MRD/Seurat.rds",
                                "~/Documents/Projects/Data_files_temp/Donor/PT5-DX/Seurat.rds",
                                "~/Documents/Projects/Data_files_temp/Donor/PT12-DX/Seurat.rds",
                                "~/Documents/Projects/Data_files_temp/Donor/PT12-REL/Seurat.rds"),
                     results_name="210824_RF_5BM_10BPDCN",
                     ordering_vector=c(12,8,3,7,5,14,13,16,11,10,2,15,9,1,17,4,6),
                     cluster_names=c("CD4NaiveT", "B", "MidEryth", "CD8MemoryT", "Mono", "NK","LateEryth", "EarlyEryth", "CD4MemoryT", "PreB",
                                     "ProB", "HSCProg", "cDC", "ncMono", "Plasma",  "pDC", "CD8NaiveT"),
                     cluster_colors=c("#C0007C", "#2a4bd7", "#f08080", "#e9debb", "#ffee33", "#BE8A66", "#ad2323", "#ff9233", "#ffcdf3", "#9dafff",
                                      "#8126c0", "#000000", "#1d6914", "#29d0d0", "#575757",  "#81c57a", "#708090"))

#############################

scpred(reference="~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_0.5_ccrem/TCELL_SUB/Seurat_int_0.5_5BM.rds",
                     datasets=c("~/Documents/Projects/Data_files_temp/Donor/BPDCN712/Seurat.rds",
                                "~/Documents/Projects/Data_files_temp/Donor/BPDCN712R/Seurat.rds",
                                "~/Documents/Projects/Data_files_temp/Donor/BPDCN180329/Seurat.rds",
                                "~/Documents/Projects/Data_files_temp/Donor/BPDCN181128/Seurat.rds",
                                "~/Documents/Projects/Data_files_temp/Donor/BPDCN190711/Seurat.rds",
                                "~/Documents/Projects/Data_files_temp/Donor/PT1-DX/Seurat.rds",
                                "~/Documents/Projects/Data_files_temp/Donor/PT1-MRD/Seurat.rds",
                                "~/Documents/Projects/Data_files_temp/Donor/PT5-DX/Seurat.rds",
                                "~/Documents/Projects/Data_files_temp/Donor/PT12-DX/Seurat.rds",
                                "~/Documents/Projects/Data_files_temp/Donor/PT12-REL/Seurat.rds"),
                     results_name="210824_scPred_5BM_10BPDCN",
                     ordering_vector=c(12,8,3,7,5,14,13,16,11,10,2,15,9,1,17,4,6),
                     cluster_names=c("CD4NaiveT", "B", "MidEryth", "CD8MemoryT", "Mono", "NK","LateEryth", "EarlyEryth", "CD4MemoryT", "PreB",
                                     "ProB", "HSCProg", "cDC", "ncMono", "Plasma",  "pDC", "CD8NaiveT"),
                     cluster_colors=c("#C0007C", "#2a4bd7", "#f08080", "#e9debb", "#ffee33", "#BE8A66", "#ad2323", "#ff9233", "#ffcdf3", "#9dafff",
                                      "#8126c0", "#000000", "#1d6914", "#29d0d0", "#575757",  "#81c57a", "#708090"))

scpred(reference="~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_0.5_ccrem/TCELL_SUB/Seurat_int_0.5_TCELLSPLIT.rds",
              datasets=c("~/Documents/Projects/Data_files_temp/Donor/BM191119/Seurat.rds",
                         "~/Documents/Projects/Data_files_temp/Donor/BM191227/Seurat.rds"),
              results_name="results_test_scpred1",
              ordering_vector=c(12,8,3,7,5,14,13,16,11,10,2,15,9,1,17,4,6),
              cluster_names=c("CD4NaiveT", "B", "MidEryth", "CD8MemoryT", "Mono", "NK","LateEryth", "EarlyEryth", "CD4MemoryT", "PreB",
                              "ProB", "HSCProg", "cDC", "ncMono", "Plasma",  "pDC", "CD8NaiveT"),
              cluster_colors=c("#C0007C", "#2a4bd7", "#f08080", "#e9debb", "#ffee33", "#BE8A66", "#ad2323", "#ff9233", "#ffcdf3", "#9dafff",
                               "#8126c0", "#000000", "#1d6914", "#29d0d0", "#575757",  "#81c57a", "#708090"))

scpred(reference="~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_0.5_ccrem/TCELL_SUB/Seurat_int_0.5_TCELLSPLIT.rds",
              datasets=c("~/Documents/Projects/Data_files_temp/Donor/BM191119/Seurat.rds",
                         "~/Documents/Projects/Data_files_temp/Donor/BM191227/Seurat.rds"),
              results_name="results_test_scpred2")
