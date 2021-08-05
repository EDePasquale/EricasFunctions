###########################
#                         #
# Test function calls for #
#    Ericas_Functions     #
#                         #
###########################

qc_and_filter(c("/Volumes/broad_vangalenlab/depasquale/CellRanger/210628_GeneCounts",
                "/Volumes/broad_vangalenlab/depasquale/CellRanger/210422_GeneCounts",
                "/Volumes/broad_vangalenlab/depasquale/CellRanger/210423_GeneCounts"), "GRCh38", "QC2")

qc_and_filter(c("/Volumes/broad_vangalenlab/depasquale/CellRanger/210628_GeneCounts"), "GRCh38", "QC3")
