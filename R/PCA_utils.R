# Create row/column annotations for distance matrix heatmap.
.create_anno <- function(vars, metadata, colors, anno_type = "row", side = "right") {
    anno <- metadata[, vars, drop = FALSE]

    # Sets color palette to dittoSeq colors instead of random
    out.anno <- NULL
    if (!is.null(anno)) {
        anno.colors <- list()
        i <- 1

        for (n in names(anno)) {
            out <- list()

            for (lev in unique(anno[[n]])) {
                out[[lev]] <- colors[i]
                i <- i + 1
            }

            anno.colors[[n]] <- unlist(out)
        }

        out.anno <- HeatmapAnnotation(df = anno, col = anno.colors, which = anno_type, annotation_name_side = side)
    }
}
