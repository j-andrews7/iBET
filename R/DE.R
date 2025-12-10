#' Create an interactive Shiny app for differential expression results exploration
#'
#' This shiny app is composed of three tabs - an interactive heatmap, MA and volcano plots, and a table of full
#' differential expression results. The interactive heatmap will generate a sub-heatmap for selected rows/columns.
#' This function accepts DE results from any analysis tool (DESeq2, edgeR, limma, etc.) and allows customization
#' of column names for significance, fold change, and abundance metrics.
#'
#' Gene labels can be added to the MAplot and volcano plot by clicking a point. The labels can also be dragged around,
#' though adding labels will reset the positions, so it's recommended to add all labels prior to re-positioning them.
#'
#' @details Note that significance values of 0 will always be pushed above the y-limit of the volcano plot,
#'   as they are infinite values after log transformation.
#'
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @import InteractiveComplexHeatmap
#' @import ComplexHeatmap
#' @import DT
#' @importFrom GetoptLong qq
#' @import circlize
#' @importFrom grid unit grid.newpage grid.text gpar
#' @importFrom SummarizedExperiment colData assay rowData
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics par
#' @importFrom stats quantile loess fitted
#' @importFrom plotly ggplotly plotlyOutput renderPlotly toWebGL plot_ly layout add_annotations config toRGB event_data
#'   add_lines add_markers
#' @import ggplot2
#' @import shinydashboard
#' @import dashboardthemes
#' @importFrom shinyWidgets prettyCheckbox dropdownButton tooltipOptions
#' @importFrom shinycssloaders withSpinner
#' @importFrom shinyjqui jqui_resizable
#' @importFrom shinyjs show useShinyjs hidden
#' @importFrom shinyBS tipify popify
#' @importFrom colourpicker colourInput
#' @importFrom methods is
#' @importFrom htmlwidgets saveWidget
#' @importFrom DESeq2 vst rlog normTransform
#' @importFrom edgeR cpm
#'
#' @param mat A numeric matrix of expression values (e.g., normalized counts, VST-transformed counts, log-CPM),
#'   a \code{\link[DESeq2]{DESeqDataSet}} object, or a \code{\link[edgeR]{DGEList}} object.
#'   Rows are features (genes) and columns are samples. Required.
#' @param res A data frame or named list of data frames containing differential expression results.
#'   Each data frame must contain columns for significance, fold change, and abundance (see sig.col, lfc.col, and abundance.col parameters).
#'   If a named list is provided, users will be able to choose between results in the app. Required.
#' @param metadata A data frame containing sample metadata. Rows should correspond to columns in \code{mat}.
#'   If provided, can be used for heatmap annotations and sample filtering. 
#'   If \code{mat} is a DESeqDataSet or DGEList with sample information, metadata will be extracted automatically if not provided. Optional.
#' @param assay String specifying which transformation/assay to extract when \code{mat} is a DESeqDataSet or DGEList.
#'   For DESeqDataSet: "vst" (default), "rlog", "normTransform", or a specific assay name.
#'   For DGEList: "logCPM" (default), "cpm", or "counts".
#'   Ignored if \code{mat} is already a matrix. Optional.
#' @param annot.by A string or character vector containing the names of metadata columns to be used as heatmap annotations.
#'   Only used if \code{metadata} is provided. Optional.
#' @param sig.col String specifying the column name in \code{res} containing significance values (e.g., "padj", "FDR", "svalue").
#'   Will be auto-detected from common column names if not provided.
#' @param lfc.col String specifying the column name in \code{res} containing log2 fold change values (e.g., "log2FoldChange", "logFC").
#'   Defaults to "log2FoldChange".
#' @param abundance.col String specifying the column name in \code{res} containing abundance/expression values (e.g., "baseMean", "AveExpr").
#'   Defaults to "baseMean".
#' @param h.id String indicating unique ID for interactive heatmaps.
#'   Required if multiple apps are run within the same Rmd file. Defaults to "ht1".
#' @param genesets Optional named list containing genesets that can be interactively highlighted on the plots.
#'   The elements of the list should each be a geneset with gene identifiers matching those used as rownames in \code{mat} and \code{res}.
#' @param swap.rownames String. The column name in \code{res} used to identify features instead of rownames.
#'   Note if this column contains duplicates (e.g. gene symbols), they will be made unique via the addition of ".1", ".2", etc.
#' @param height Number indicating height of app in pixels. Defaults to 800.
#'
#' @return A Shiny app containing interconnected InteractiveComplexHeatmap, MAplot, and volcano plots along with full DE results.
#'
#' @examples
#' \dontrun{
#' # Example with DESeq2 - direct object usage
#' library(DESeq2)
#' dds <- DESeq(dds)
#' res <- results(dds)
#' shinyDE(dds, res = as.data.frame(res))
#'
#' # Example with DESeq2 - multiple comparisons
#' res1 <- results(dds, contrast = c("condition", "treated", "untreated"))
#' res2 <- results(dds, contrast = c("condition", "treated2", "untreated"))
#' res_list <- list("treated vs untreated" = as.data.frame(res1),
#'                  "treated2 vs untreated" = as.data.frame(res2))
#' shinyDE(dds, res = res_list, assay = "vst")
#'
#' # Example with DESeq2 - using rlog transformation
#' shinyDE(dds, res = as.data.frame(res), assay = "rlog")
#'
#' # Example with DESeq2 - manual extraction
#' mat <- assay(vst(dds))
#' shinyDE(mat, res = as.data.frame(res), metadata = colData(dds))
#' 
#' # Example with edgeR - direct object usage
#' library(edgeR)
#' y <- DGEList(counts = counts)
#' y <- calcNormFactors(y)
#' # ... run DE analysis ...
#' res <- topTags(qlf, n = Inf)$table
#' shinyDE(y, res = res, lfc.col = "logFC", abundance.col = "logCPM", sig.col = "FDR")
#' 
#' # Example with edgeR - using CPM without log
#' shinyDE(y, res = res, assay = "cpm", lfc.col = "logFC", abundance.col = "logCPM", sig.col = "FDR")
#' 
#' # Example with manual matrix
#' shinyDE(mat, res = res, lfc.col = "logFC", abundance.col = "logCPM", sig.col = "FDR")
#' }
#'
#' @author Jared Andrews
#' @export
shinyDE <- function(mat, res = NULL, metadata = NULL, annot.by = NULL,
                    sig.col = NULL, lfc.col = "log2FoldChange", abundance.col = "baseMean",
                    assay = NULL, h.id = "ht1", genesets = NULL, swap.rownames = NULL, height = 800) {

    # Store original mat object and feature metadata for swap.rownames handling
    feature_metadata <- NULL
    original_mat <- mat

    # Handle DESeqDataSet objects
    if (is(mat, "DESeqDataSet")) {
        # Extract feature metadata for swap.rownames
        if (!is.null(swap.rownames)) {
            feature_metadata <- as.data.frame(rowData(original_mat))
        }

        if (is.null(metadata)) {
            metadata <- as.data.frame(colData(mat))
        }

        if (is.null(assay)) {
            assay <- "vst"  # default
        }

        # Extract the appropriate assay
        if (assay == "vst") {
            mat <- assay(vst(mat))
        } else if (assay == "rlog") {
            mat <- assay(rlog(mat))
        } else if (assay == "normTransform") {
            mat <- assay(normTransform(mat))
        } else {
            # Try to get a named assay
            mat <- assay(mat, assay)
        }
    }
    
    # Handle DGEList objects
    if (is(mat, "DGEList")) {
        # Extract feature metadata for swap.rownames (edgeR stores in $genes)
        if (!is.null(swap.rownames) && !is.null(original_mat$genes)) {
            feature_metadata <- as.data.frame(original_mat$genes)
        }

        if (is.null(metadata) && !is.null(mat$samples)) {
            metadata <- mat$samples
        }

        # Determine which assay to use
        if (is.null(assay)) {
            assay <- "logCPM"  # default
        }

        # Extract the appropriate values
        if (assay == "logCPM") {
            mat <- cpm(mat, log = TRUE)
        } else if (assay == "cpm") {
            mat <- cpm(mat, log = FALSE)
        } else if (assay == "counts") {
            mat <- mat$counts
        } else {
            stop(paste0("Unknown assay type '", assay, "' for DGEList. Use 'logCPM', 'cpm', or 'counts'."))
        }

        # Set default column names for edgeR if they match DESeq2 defaults
        if (lfc.col == "log2FoldChange") {
            lfc.col <- "logFC"
        }
        if (abundance.col == "baseMean") {
            abundance.col <- "logCPM"
        }
        if (is.null(sig.col)) {
            sig.col <- "FDR"
        }
    }

    # Validate res is provided
    if (is.null(res)) {
        stop("res parameter is required. Please provide differential expression results.")
    }

    # Validate inputs
    if (!is.matrix(mat) && !is.data.frame(mat)) {
        stop("mat must be a matrix, data frame, DESeqDataSet, or DGEList object")
    }

    if (is.data.frame(mat)) {
        mat <- as.matrix(mat)
    }

    # Check if res is individual data frame or named list.
    multi.res <- FALSE
    res.list <- NULL
    if (is(res, "list")) {
        if (!is.null(names(res))) {
            multi.res <- TRUE
            res.list <- res
            res <- res.list[[1]]
        } else {
            stop("Results list elements should be named.")
        }
    }

    # Convert to data frame if needed
    if (!is.data.frame(res)) {
        res <- as.data.frame(res)
    }

    # Handle swap.rownames: create mapping from original IDs to new IDs
    if (!is.null(swap.rownames)) {
        # Try to find the swap column in different sources
        swap_values <- NULL
        swap_source <- NULL

        # 1. Check feature_metadata (from rowData or genes)
        if (!is.null(feature_metadata) && swap.rownames %in% colnames(feature_metadata)) {
            swap_values <- feature_metadata[[swap.rownames]]
            swap_source <- "feature_metadata"
        }
        # 2. Check results dataframe
        else if (swap.rownames %in% colnames(res)) {
            swap_values <- res[[swap.rownames]]
            swap_source <- "results"
        }
        # 3. If res is a list, check all results for the column
        else if (multi.res) {
            for (i in seq_along(res.list)) {
                if (swap.rownames %in% colnames(res.list[[i]])) {
                    swap_values <- res.list[[i]][[swap.rownames]]
                    swap_source <- paste0("results list element ", i)
                    break
                }
            }
        }

        if (is.null(swap_values)) {
            stop(paste0("Column '", swap.rownames, "' not found in feature metadata (rowData/genes) or results dataframe."))
        }

        # Create mapping from original rownames to new identifiers
        original_ids <- rownames(mat)

        # Ensure swap_values matches mat dimensions
        if (length(swap_values) != length(original_ids)) {
            # Try to match by rownames if lengths differ
            if (swap_source == "results" && all(rownames(res) %in% original_ids)) {
                # Reorder swap_values to match mat
                swap_values <- swap_values[match(original_ids, rownames(res))]
            } else {
                warning(paste0("Length mismatch between expression matrix (", length(original_ids), 
                             ") and swap column (", length(swap_values), "). Using first match."))
                min_len <- min(length(original_ids), length(swap_values))
                original_ids <- original_ids[seq_len(min_len)]
                swap_values <- swap_values[seq_len(min_len)]
                mat <- mat[seq_len(min_len), , drop = FALSE]
            }
        }

        # Make identifiers unique if needed
        new_ids <- make.names(swap_values, unique = TRUE)
        rownames(mat) <- new_ids

        # Apply new rownames to res and res.list
        # Create mapping: original rowname -> new rowname
        id_mapping <- setNames(new_ids, original_ids)
        res <- .swap_res_rownames(res, id_mapping, original_ids)

        if (multi.res) {
            res.list <- lapply(res.list, function(r) .swap_res_rownames(r, id_mapping, original_ids))
        }

        message(paste0("Swapped rownames using '", swap.rownames, "' from ", swap_source, 
                       ". Matched ", nrow(res), " features."))
    }

    # Validate required columns exist
    if (!lfc.col %in% colnames(res)) {
        stop(paste0("Column '", lfc.col, "' not found in results. Please specify correct column name with lfc.col parameter."))
    }

    if (!abundance.col %in% colnames(res)) {
        stop(paste0("Column '", abundance.col, "' not found in results. Please specify correct column name with abundance.col parameter."))
    }
    
    # Auto-detect significance column if not provided
    if (is.null(sig.col)) {
        common_sig_cols <- c("padj", "FDR", "adj.P.Val", "svalue", "qvalue", "pvalue", "PValue", "P.Value")
        detected <- intersect(common_sig_cols, colnames(res))
        if (length(detected) > 0) {
            sig.col <- detected[1]
            message(paste0("Auto-detected significance column: ", sig.col))
        } else {
            stop("Could not auto-detect significance column. Please specify with sig.col parameter. Common names include: padj, FDR, adj.P.Val, svalue, qvalue")
        }
    } else {
        if (!sig.col %in% colnames(res)) {
            stop(paste0("Column '", sig.col, "' not found in results. Please specify correct column name with sig.col parameter."))
        }
    }
    
    # Ensure mat and res have matching features
    common_feats <- intersect(rownames(mat), rownames(res))
    if (length(common_feats) == 0) {
        stop("No matching features found between expression matrix and results. Check that rownames match.")
    }
    if (length(common_feats) < nrow(res)) {
        warning(paste0("Only ", length(common_feats), " of ", nrow(res), 
                      " features from results found in expression matrix. Subsetting to common features."))
    }
    mat <- mat[common_feats, , drop = FALSE]
    res <- res[common_feats, , drop = FALSE]
    if (multi.res) {
        res.list <- lapply(res.list, function(r) r[intersect(rownames(r), common_feats), , drop = FALSE])
    }

    # Parameter validation.
    if (!is.null(genesets)) {
        if (is.null(names(genesets))) {
            stop("Genesets list must be named")
        } else if (!is(genesets, "list")) {
            stop("Genesets must be provided as a named list")
        }
    }

    if (!is.null(metadata) && !is.null(annot.by)) {
        if (!all(annot.by %in% names(metadata))) {
            stop("Annotation variables not found in metadata")
        }
    }

    env <- new.env()
    env$lfc.col <- lfc.col
    env$abundance.col <- abundance.col
    env$sig.col <- sig.col

    qa <- quantile(log10(res[[abundance.col]] + 1), 0.99, na.rm = TRUE)
    abundance_col_fun <- circlize::colorRamp2(c(0, qa / 2, qa), c("blue", "white", "red"))

    qa <- quantile(abs(res[[lfc.col]][!is.na(res[[lfc.col]])]), 0.99, na.rm = TRUE)
    log2fc_col_fun <- circlize::colorRamp2(c(-qa, 0, qa), c("blue", "white", "red"))

    sig.term <- sig.col
    coef <- paste0(sig.col, " / ", lfc.col)  # For display purposes

    environment(.make_heatmap) <- env

    body <- mainPanel(
        width = 10,
        tabsetPanel(
            tabPanel(
                "Heatmap",
                fluidRow(
                    column(
                        width = 4,
                        h3("Differential heatmap"),
                        originalHeatmapOutput(h.id, height = 500, width = 400, containment = TRUE)
                    ),
                    column(
                        width = 4,
                        id = "column2",
                        h3("Sub-heatmap"),
                        subHeatmapOutput(h.id, title = NULL, width = 400, containment = TRUE),
                        h3(title = "Output"),
                        HeatmapInfoOutput(h.id, title = NULL, width = 400),
                    ),
                    column(
                        width = 4,
                        h3("Results for the Selected Genes"),
                        div(DT::dataTableOutput("res_table"), style = "font-size:80%")
                    )
                )
            ),
            tabPanel(
                "MA & Volcano Plots",
                fluidRow(
                    column(
                        width = 6,
                        br(),
                        dropdownButton(
                            tags$h3("Plot Settings"),
                            fluidRow(
                                column(
                                    width = 6,
                                    tipify(colourInput("ma.down.color", "Down-reg colour", value = "#0026ff"),
                                        "Color of down-regulated genes.", "right",
                                        options = list(container = "body")
                                    ),
                                    tipify(colourInput("ma.up.color", "Up-reg colour", value = "red"),
                                        "Color of up-regulated genes.", "right",
                                        options = list(container = "body")
                                    ),
                                    tipify(colourInput("ma.insig.color", "Insig colour", value = "black"),
                                        "Color of insignificant genes.", "right",
                                        options = list(container = "body")
                                    ),
                                    tipify(numericInput("ma.sig.opa", label = "Sig opacity:", value = 1, step = 0.05, min = 0),
                                        "Opacity of DE genes.", "right",
                                        options = list(container = "body")
                                    ),
                                    tipify(
                                        prettyCheckbox("ma.loess",
                                            label = "Show LOESS line", TRUE, bigger = TRUE,
                                            animation = "smooth", status = "success",
                                            icon = icon("check"), width = "100%"
                                        ),
                                        "Draw LOESS line over all points (scatterplot smoother).", "right",
                                        options = list(container = "body")
                                    ),
                                    tipify(colourInput("ma.loess.color", "LOESS colour", value = "#fd8f00"),
                                        "Color of LOESS line.", "right",
                                        options = list(container = "body")
                                    )
                                ),
                                column(
                                    width = 6,
                                    numericInput("ma.y", label = "y-axis limits:", value = 5, step = 0.1, min = 0.1),
                                    tipify(numericInput("ma.lab.size", label = "Label size:", value = 10, step = 0.5, min = 1),
                                        "Font size of gene labels.", "right",
                                        options = list(container = "body")
                                    ),
                                    tipify(numericInput("ma.insig.opa", label = "Insig opacity:", value = 1, step = 0.05, min = 0),
                                        "Opacity of non-DE genes.", "right",
                                        options = list(container = "body")
                                    ),
                                    tipify(numericInput("ma.sig.size", label = "Sig pt size:", value = 5, step = 0.1, min = 1),
                                        "Point size of DE genes.", "right",
                                        options = list(container = "body")
                                    ),
                                    tipify(numericInput("ma.insig.size", label = "Insig pt size:", value = 3, step = 0.1, min = 0),
                                        "Point size of non-DE genes.", "right",
                                        options = list(container = "body")
                                    ),
                                    tipify(
                                        prettyCheckbox("ma.loess.hl.genesets",
                                            label = "Show GSets LOESS", TRUE, bigger = TRUE,
                                            animation = "smooth", status = "success",
                                            icon = icon("check"), width = "100%"
                                        ),
                                        "Draw LOESS line over highlighted genesets (scatterplot smoother).", "right",
                                        options = list(container = "body")
                                    ),
                                    tipify(colourInput("ma.loess.hl.genesets.color", "GSets LOESS colour", value = "#23A39D"),
                                        "Color of genesets LOESS line.", "right",
                                        options = list(container = "body")
                                    )
                                )
                            ),
                            splitLayout(
                                tipify(numericInput("ma.loess.span", label = "LOESS span:", value = 0.75, max = 1, step = 0.5, min = 0.01),
                                    "Smoothness of LOESS (higher is smoother).", "right",
                                    options = list(container = "body")
                                ),
                                tipify(numericInput("ma.loess.hl.genesets.span", label = "LOESS genesets span:", value = 0.75, max = 1, step = 0.5, min = 0.01),
                                    "Smoothness of geneset LOESS (higher is smoother).", "right",
                                    options = list(container = "body")
                                ),
                            ),
                            tipify(
                                prettyCheckbox("ma.fcline",
                                    label = "Show MAplot FC Threshold", value = TRUE,
                                    animation = "smooth", status = "success", bigger = TRUE, icon = icon("check")
                                ),
                                "Draw lines for FC threshold.", "right",
                                options = list(container = "body")
                            ),
                            tipify(
                                prettyCheckbox("ma.hl.counts",
                                    label = "Show highlight counts", FALSE, bigger = TRUE,
                                    animation = "smooth", status = "success",
                                    icon = icon("check"), width = "100%"
                                ),
                                "Show count of highlighted genes and genesets.", "right",
                                options = list(container = "body")
                            ),
                            splitLayout(
                                tipify(
                                    prettyCheckbox("ma.webgl",
                                        label = "Use webGL", TRUE, bigger = TRUE,
                                        animation = "smooth", status = "success",
                                        icon = icon("check"), width = "100%"
                                    ),
                                    "Plot with webGL. Faster, but sometimes has visual artifacts.", "right",
                                    options = list(container = "body")
                                ),
                                tipify(
                                    prettyCheckbox("ma.counts",
                                        label = "Show counts", TRUE, bigger = TRUE,
                                        animation = "smooth", status = "success",
                                        icon = icon("check"), width = "100%"
                                    ),
                                    "Show count of DE and total genes.", "right",
                                    options = list(container = "body")
                                )
                            ),
                            splitLayout(
                                tipify(numericInput("ma.webgl.ratio", label = "webGL pixel ratio:", value = 7, step = 0.1, min = 1),
                                    "Controls rasterization resolution when webGL is used. Higher is greater resolution, recommend leaving at default.",
                                    "right",
                                    options = list(container = "body")
                                ),
                                tipify(numericInput("ma.counts.size", label = "Counts size:", value = 8, step = 0.1, min = 0),
                                    "Font size of gene counts.", "right",
                                    options = list(container = "body")
                                )
                            ),
                            circle = FALSE, label = strong("MA-Plot"), status = "danger", size = "lg", icon = icon("gear"),
                            width = "300px", tooltip = tooltipOptions(title = "Click to change plot settings")
                        ),
                        withSpinner(
                            jqui_resizable(
                                plotlyOutput("ma_plot", height = "500px", width = "550px")
                            )
                        ),
                        div(downloadButton("download_plotly_ma", "Download Interactive MAplot"), align = "left", style = "margin-top: 10px;")
                    ),
                    column(
                        width = 6,
                        br(),
                        dropdownButton(
                            tags$h3("Plot Settings"),
                            fluidRow(
                                column(
                                    width = 6,
                                    tipify(colourInput("vol.down.color", "Down-reg colour", value = "#0026ff"),
                                        "Color of down-regulated genes.", "right",
                                        options = list(container = "body")
                                    ),
                                    tipify(colourInput("vol.up.color", "Up-reg colour", value = "red"),
                                        "Color of up-regulated genes.", "right",
                                        options = list(container = "body")
                                    ),
                                    tipify(colourInput("vol.insig.color", "Insig colour", value = "#A6A6A6"),
                                        "Color of insignificant genes.", "right",
                                        options = list(container = "body")
                                    ),
                                    tipify(numericInput("vol.sig.opa", label = "Sig opacity:", value = 1, step = 0.05, min = 0),
                                        "Opacity of DE genes.", "right",
                                        options = list(container = "body")
                                    ),
                                    tipify(numericInput("vol.sig.size", label = "Sig pt size:", value = 5, step = 0.1, min = 0),
                                        "Point size of DE genes.", "right",
                                        options = list(container = "body")
                                    )
                                ),
                                column(
                                    width = 6,
                                    numericInput("vol.x", label = "x-axis limits:", value = 5, step = 0.1, min = 0.1),
                                    numericInput("vol.y",
                                        label = "y-axis limits:",
                                        value = max(-log10(res[[sig.term]][!is.na(res[[sig.term]])])) + 0.1,
                                        step = 0.5, min = 1
                                    ),
                                    tipify(numericInput("vol.lab.size", label = "Label size:", value = 10, step = 0.5, min = 1),
                                        "Font size of gene labels.", "right",
                                        options = list(container = "body")
                                    ),
                                    tipify(numericInput("vol.insig.opa", label = "Insig opacity:", value = 1, step = 0.05, min = 0),
                                        "Opacity of non-DE genes.", "right",
                                        options = list(container = "body")
                                    ),
                                    tipify(numericInput("vol.insig.size", label = "Insig pt size:", value = 3, step = 0.1, min = 0),
                                        "Point size of non-DE genes.", "right",
                                        options = list(container = "body")
                                    ),
                                )
                            ),
                            prettyCheckbox("vol.fcline",
                                label = "Show FC Threshold", value = TRUE,
                                animation = "smooth", status = "success", bigger = TRUE, icon = icon("check")
                            ),
                            tipify(
                                prettyCheckbox("vol.sigline",
                                    label = "Show Signficance Threshold", value = TRUE,
                                    animation = "smooth", status = "success", bigger = TRUE, icon = icon("check")
                                ),
                                "Draw line at significance threshold.", "right",
                                options = list(container = "body")
                            ),
                            tipify(
                                prettyCheckbox("vol.hl.counts",
                                    label = "Show Highlight Counts", value = FALSE,
                                    animation = "smooth", status = "success", bigger = TRUE, icon = icon("check")
                                ),
                                "Show count of highlighted genes and genesets.", "right",
                                options = list(container = "body")
                            ),
                            splitLayout(
                                tipify(
                                    prettyCheckbox("vol.counts",
                                        label = "Show counts", TRUE, bigger = TRUE,
                                        animation = "smooth", status = "success",
                                        icon = icon("check"), width = "100%"
                                    ),
                                    "Show count of DE and total genes.", "right",
                                    options = list(container = "body")
                                ),
                                tipify(
                                    prettyCheckbox("vol.webgl",
                                        label = "Use webGL", TRUE, bigger = TRUE,
                                        animation = "smooth", status = "success",
                                        icon = icon("check"), width = "100%"
                                    ),
                                    "Plot with webGL. Faster, but sometimes has visual artifacts.", "right",
                                    options = list(container = "body")
                                )
                            ),
                            splitLayout(
                                tipify(numericInput("vol.counts.size", label = "Counts size:", value = 8, step = 0.1, min = 0),
                                    "Font size of gene counts.", "right",
                                    options = list(container = "body")
                                ),
                                tipify(numericInput("vol.webgl.ratio", label = "webGL pixel ratio:", value = 7, step = 0.1, min = 1),
                                    "Controls rasterization resolution when webGL is used. Higher is greater resolution, recommend leaving at default.",
                                    "right",
                                    options = list(container = "body")
                                )
                            ),
                            circle = FALSE, label = strong("Volcano Plot"), status = "danger", size = "lg", icon = icon("gear"),
                            width = "300px", tooltip = tooltipOptions(title = "Click to change plot settings")
                        ),
                        withSpinner(jqui_resizable(plotlyOutput("volcano_plot", height = "500px", width = "550px"))),
                        div(downloadButton("download_plotly_volc", "Download Interactive Volcano plot"), align = "left", style = "margin-top: 10px;")
                    )
                )
            ),
            tabPanel(
                "Full DE Results Table",
                div(DT::dataTableOutput("res_table_full"), style = "font-size:80%; margin:10px;")
            ),
            tabPanel(
                "Sample Metadata (Filtering)",
                div(DT::dataTableOutput("metadata"), style = "font-size:80%; margin:10px;")
            )
        )
    )

    # Side bar contains settings for certain cutoffs to select significant genes.
    ui <- dashboardPage(
        dashboardHeader(disable = TRUE),
        dashboardSidebar(disable = TRUE),
        dashboardBody(
            tags$head(
                # Note the wrapping of the string in HTML()
                tags$style(HTML("
          .panel-body {
            padding: 5px;
          }
          .form-group {
            margin-bottom: 5px;
          }
          .well {
            padding: 5px;
            margin-bottom: 10px;
          }
          label {
            font-size: 80%;
          }
          .form-control, .selectize-input{
            padding-bottom: 2px !important;
            padding-top: 2px !important;
            font-size: 10px;
            height: 24px;
          }
        "))
            ),
            useShinyjs(),
            shinyDashboardThemes(
                theme = "onenote"
            ),
            sidebarLayout(
                sidebarPanel(
                    width = 2,
                    tags$label(HTML(qq("Comparison: <code style='font-weight:normal; font-size: 10px;'>@{paste(coef, collapse = ' ')}</code>")),
                        class = "shiny-input-container", style = "font-size:1.2em;"
                    ),
                    hidden(div(id = "mres", tipify(selectInput("res.select", NULL, choices = names(res.list)),
                        "Results to view from provided set.", "right",
                        options = list(container = "body")
                    ))),
                    hr(style = "margin:2px; background-color: #737373;"),
                    bsCollapse(
                        open = "settings",
                        bsCollapsePanel(
                            title = span(icon("plus"), "Plot Settings"), value = "settings", style = "info",
                            tipify(numericInput("sig.thresh", label = "Significance threshold:", value = 0.05, step = 0.001, min = 0.0001),
                                "Significance threshold to consider a gene DE.", "right",
                                options = list(container = "body")
                            ),
                            tipify(
                                selectInput("sig.term",
                                    label = "Significance term:", choices = colnames(res),
                                    selected = sig.col
                                ),
                                "Significance term to use for DE filtering.", "right",
                                options = list(container = "body")
                            ),
                            tipify(numericInput("base_mean", label = paste0("Minimal ", abundance.col, ":"), value = 0, step = 1),
                                paste0("Minimal ", abundance.col, " (abundance) required to consider a feature DE."), "right",
                                options = list(container = "body")
                            ),
                            tipify(numericInput("log2fc", label = paste0("Minimal abs(", lfc.col, "):"), value = 0, step = 0.1, min = 0),
                                paste0(lfc.col, " magnitude threshold required to consider a feature DE."), "right",
                                options = list(container = "body")
                            ),
                            tipify(numericInput("row.km", label = "Row k-means groups:", value = 2, step = 1),
                                "Number of groups to break heatmap into via  k-means clustering on rows.", "right",
                                options = list(container = "body")
                            ),
                            tipify(numericInput("col.km", label = "Column k-means groups:", value = 0, step = 1),
                                "Number of groups to break heatmap into via k-means clustering on columns.", "right",
                                options = list(container = "body")
                            ),
                            tipify(
                                pickerInput("anno.vars", "Annotate Heatmap by:",
                                    choices = if (!is.null(metadata)) c("", names(metadata)) else c(""),
                                    multiple = TRUE, options = list(`live-search` = TRUE, `actions-box` = TRUE),
                                    selected = annot.by
                                ),
                                "Sample metadata columns used for column annotations.", "right",
                                options = list(container = "body")
                            )
                        ),
                        bsCollapsePanel(
                            title = span(icon("plus"), "Highlight Gene(sets)"), value = "genesets", style = "info",
                            textAreaInput("hl.genes", "Highlight Genes:",
                                value = "", rows = 4,
                                placeholder = "Enter space, comma, or newline delimited genes"
                            ),
                            pickerInput("hl.genesets", "Highlight Genesets:",
                                choices = c("", names(genesets)),
                                multiple = TRUE, options = list(
                                    `live-search` = TRUE,
                                    `actions-box` = TRUE
                                )
                            ),
                            fluidRow(
                                column(
                                    6,
                                    tipify(numericInput("hl.genes.opa", label = "Genes opacity:", value = 1, step = 0.05, min = 0),
                                        "Opacity of highlighted genes.", "right",
                                        options = list(container = "body")
                                    ),
                                    tipify(numericInput("hl.genes.size", label = "Genes pt size:", value = 7, step = 0.1, min = 0),
                                        "Point size of highlighted genes.", "right",
                                        options = list(container = "body")
                                    ),
                                    tipify(numericInput("hl.genes.lw", label = "Genes border width:", value = 0.5, step = 0.05, min = 0),
                                        "Width of border for highlighted genes.", "right",
                                        options = list(container = "body")
                                    ),
                                    tipify(colourInput("hl.genes.col", "Genes color:", value = "#FFFF19"),
                                        "Color of genes to highlight.", "right",
                                        options = list(container = "body")
                                    ),
                                    tipify(colourInput("hl.genes.lcol", "Genes border:", value = "#000000"),
                                        "Border color of genes to highlight.", "right",
                                        options = list(container = "body")
                                    )
                                ),
                                column(
                                    6,
                                    tipify(numericInput("hl.genesets.opa", label = "Sets opacity:", value = 1, step = 0.05, min = 0),
                                        "Opacity of highlighted genesets.", "right",
                                        options = list(container = "body")
                                    ),
                                    tipify(numericInput("hl.genesets.size", label = "Sets pt size:", value = 7, step = 0.1, min = 0),
                                        "Point size of highlighted genesets.", "right",
                                        options = list(container = "body")
                                    ),
                                    tipify(numericInput("hl.genesets.lw", label = "Sets border width:", value = 0.5, step = 0.05, min = 0),
                                        "Width of border for highlighted genesets.", "right",
                                        options = list(container = "body")
                                    ),
                                    tipify(colourInput("hl.genesets.col", "Sets color:", value = "#38FFF2"),
                                        "Color of genesets to highlight.", "right",
                                        options = list(container = "body")
                                    ),
                                    tipify(colourInput("hl.genesets.lcol", "Sets border:", value = "#000000"),
                                        "Border color of genesets to highlight.", "right",
                                        options = list(container = "body")
                                    )
                                )
                            )
                        )
                    ),
                    div(actionButton("update", label = "Update Plots"), align = "center")
                ),
                body
            )
        )
    )

    server <- function(input, output, session) {
        if (multi.res) {
            shinyjs::show("mres")
        }

        # Get annotations. If none provided, use design variables.
        anno <- reactiveVal()

        # Used to hold plots for download.
        plot_store <- reactiveValues()

        observeEvent(c(
            input$anno.vars,
            input$res.select,
            input$metadata_rows_all,
            input$update
        ), {
            annos <- metadata

            if (!is.null(annos) && !is.null(input$anno.vars)) {
                annos <- annos[, input$anno.vars, drop = FALSE]
            } else if (is.null(input$anno.vars)) {
                annos <- NULL
            }

            if (!is.null(annos)) {
                l <- sapply(annos, function(x) (is.factor(x) || is.character(x))) | sapply(annos, function(x) length(unique(x)) == 1)
                annos <- annos[, l, drop = FALSE]
            }

            if (!is.null(annos) && ncol(annos) == 0) {
                annos <- NULL
            }

            # Filter samples.
            if (!is.null(input$metadata_rows_all)) {
                annos <- annos[input$metadata_rows_all, , drop = FALSE]
            }

            anno(annos)
        })

        # Keep track of which genes have been clicked
        genes <- reactiveValues(ma = NULL, volc = NULL)

        ress <- reactiveVal({
            res
        })

        # Get the selected results tables.
        observeEvent(input$res.select,
            {
                req(ress, anno)
                ress(res.list[[input$res.select]])

                # Reset selected genes
                genes$ma <- NULL
                genes$volc <- NULL

                # Update significance term choices to new column names.
                updateSelectInput(session, "sig.term",
                    choices = colnames(ress()),
                    selected = sig.col
                )

                if (!is.null(input$metadata_rows_all)) {
                    mat <- mat[, input$metadata_rows_all]
                }

                pdf(NULL)
                ht <- .make_heatmap(mat, ress(), anno(), abundance_col_fun, log2fc_col_fun,
                    sig.term = input$sig.term,
                    lfc.col = lfc.col, abundance.col = abundance.col,
                    sig.thresh = as.numeric(input$sig.thresh), base_mean = input$base_mean, log2fc = input$log2fc,
                    row.km = input$row.km, col.km = input$col.km
                )
                dev.off()

                if (!is.null(ht)) {
                    makeInteractiveComplexHeatmap(input, output, session, ht, h.id,
                        brush_action = .brush_action, click_action = .click_action
                    )
                } else {
                    # The ID for the heatmap plot is encoded as @{heatmap_id}_heatmap.
                    output[[paste0(h.id, "_heatmap")]] <- renderPlot({
                        grid.newpage()
                        grid.text("No row exists after filtering.")
                    })
                }
            },
            ignoreInit = TRUE
        )

        # A self-defined action to respond to heatmap selections to show them in a table.
        .brush_action <- function(df, input, output, session) {
            row_index <- unique(unlist(df$row_index))
            selected <- env$row_index[row_index]

            # Find available columns
            display_cols <- c(abundance.col, lfc.col)
            cnames <- c("stat", "pvalue", "padj", "svalue", "lfcSE")[c("stat", "pvalue", "padj", "svalue", "lfcSE") %in% colnames(ress())]
            if (!is.null(swap.rownames)) {
                cnames <- c(cnames, "ORIGINAL_ROWS")
            }
            all_cols <- c(display_cols, cnames)

            output[["res_table"]] <- DT::renderDataTable({
                DT::formatRound(
                    DT::datatable(as.data.frame(ress()[selected, all_cols]),
                        rownames = TRUE, options = list(
                            lengthMenu = c(5, 10, 25),
                            pageLength = 20
                        )
                    ),
                    columns = seq_along(all_cols), digits = 5
                ) %>%
                    DT::formatStyle(0, target = "row", lineHeight = "40%")
            })
        }

        .click_action <- function(df, input, output, session) {
            row_index <- unique(unlist(df$row_index))
            selected <- env$row_index[row_index]

            output[["res_table"]] <- DT::renderDataTable({
                # Adjust output table columns based on results table.
                display_cols <- c(abundance.col, lfc.col)
                cnames <- c(sig.col)
                if ("lfcSE" %in% colnames(ress())) {
                    display_cols <- c(display_cols, "lfcSE")
                }
                if (!is.null(swap.rownames)) {
                    cnames <- c(cnames, "ORIGINAL_ROWS")
                }
                all_cols <- c(display_cols, cnames)
                
                DT::formatRound(
                    DT::datatable(as.data.frame(ress()[selected, all_cols]),
                        rownames = TRUE, options = list(
                            lengthMenu = c(5, 10, 25),
                            pageLength = 20
                        )
                    ),
                    columns = seq_len(length(display_cols)), digits = 5
                ) %>%
                    DT::formatStyle(0, target = "row", lineHeight = "40%")
            })
        }

        # On click, the key field of the event data contains the gene symbol
        # Add that gene to the set of all "selected" genes
        observeEvent(event_data("plotly_click", source = paste0(h.id, "_ma")), {
            gene <- event_data("plotly_click", source = paste0(h.id, "_ma"))
            gene_old_new <- rbind(genes$ma, gene)
            keep <- gene_old_new[gene_old_new$customdata %in% names(which(table(gene_old_new$customdata) == 1)), ]

            if (nrow(keep) == 0) {
                genes$ma <- NULL
            } else {
                genes$ma <- keep
            }
        })

        observeEvent(event_data("plotly_click", source = paste0(h.id, "_volc")), {
            gene <- event_data("plotly_click", source = paste0(h.id, "_volc"))
            gene_old_new <- rbind(genes$volc, gene)
            keep <- gene_old_new[gene_old_new$customdata %in% names(which(table(gene_old_new$customdata) == 1)), ]

            if (nrow(keep) == 0) {
                genes$volc <- NULL
            } else {
                genes$volc <- keep
            }
        })

        # clear the set of genes when a double-click occurs
        observeEvent(event_data("plotly_doubleclick", source = paste0(h.id, "_ma")), {
            genes$ma <- NULL
        })

        observeEvent(event_data("plotly_doubleclick", source = paste0(h.id, "_volc")), {
            genes$volc <- NULL
        })

        output$ma_plot <- renderPlotly({
            req(genes)
            input$update

            plot_store$ma_plot <- .make_maplot(
                res = ress(),
                ylim = isolate(input$ma.y),
                fc.thresh = isolate(input$log2fc),
                fc.lines = isolate(input$ma.fcline),
                sig.thresh = isolate(input$sig.thresh),
                basemean.thresh = isolate(input$base_mean),
                h.id = h.id,
                sig.term = isolate(input$sig.term),
                lfc.col = lfc.col,
                abundance.col = abundance.col,
                gs = genes$ma,
                up.color = isolate(input$ma.up.color),
                down.color = isolate(input$ma.down.color),
                insig.color = isolate(input$ma.insig.color),
                sig.opacity = isolate(input$ma.sig.opa),
                insig.opacity = isolate(input$ma.insig.opa),
                sig.size = isolate(input$ma.sig.size),
                insig.size = isolate(input$ma.insig.size),
                label.size = isolate(input$ma.lab.size),
                webgl = isolate(input$ma.webgl),
                webgl.ratio = isolate(input$ma.webgl.ratio),
                show.counts = isolate(input$ma.counts),
                show.hl.counts = isolate(input$ma.hl.counts),
                counts.size = isolate(input$ma.counts.size),
                highlight.genesets = isolate(input$hl.genesets),
                highlight.genes = isolate(input$hl.genes),
                genesets = genesets,
                highlight.genes.color = isolate(input$hl.genes.col),
                highlight.genes.size = isolate(input$hl.genes.size),
                highlight.genes.opac = isolate(input$hl.genes.opa),
                highlight.genes.linecolor = isolate(input$hl.genes.lcol),
                highlight.genes.linewidth = isolate(input$hl.genes.lw),
                highlight.genesets.color = isolate(input$hl.genesets.col),
                highlight.genesets.size = isolate(input$hl.genesets.size),
                highlight.genesets.opac = isolate(input$hl.genesets.opa),
                highlight.genesets.linecolor = isolate(input$hl.genesets.lcol),
                highlight.genesets.linewidth = isolate(input$hl.genesets.lw),
                loess = isolate(input$ma.loess),
                loess.color = isolate(input$ma.loess.color),
                loess.span = isolate(input$ma.loess.span),
                loess.genesets = isolate(input$ma.loess.hl.genesets),
                loess.genesets.color = isolate(input$ma.loess.hl.genesets.color),
                loess.genesets.span = isolate(input$ma.loess.hl.genesets.span)
            )


            plot_store$ma_plot
        })

        output$volcano_plot <- renderPlotly({
            req(genes)
            input$update

            plot_store$volcano_plot <- .make_volcano(
                res = ress(),
                xlim = isolate(input$vol.x),
                ylim = isolate(input$vol.y),
                fc.thresh = isolate(input$log2fc),
                fc.lines = isolate(input$vol.fcline),
                sig.thresh = isolate(input$sig.thresh),
                sig.line = isolate(input$vol.sigline),
                basemean.thresh = isolate(input$base_mean),
                h.id = h.id,
                sig.term = isolate(input$sig.term),
                lfc.term = lfc.col,
                abundance.col = abundance.col,
                feat.term = "rows",
                fs = genes$volc,
                up.color = isolate(input$vol.up.color),
                down.color = isolate(input$vol.down.color),
                insig.color = isolate(input$vol.insig.color),
                sig.opacity = isolate(input$vol.sig.opa),
                insig.opacity = isolate(input$vol.insig.opa),
                sig.size = isolate(input$vol.sig.size),
                insig.size = isolate(input$vol.insig.size),
                label.size = isolate(input$vol.lab.size),
                webgl = isolate(input$vol.webgl),
                webgl.ratio = isolate(input$vol.webgl.ratio),
                show.counts = isolate(input$vol.counts),
                show.hl.counts = isolate(input$vol.hl.counts),
                counts.size = isolate(input$vol.counts.size),
                highlight.featsets = isolate(input$hl.genesets),
                highlight.feats = isolate(input$hl.genes),
                featsets = genesets,
                highlight.feats.color = isolate(input$hl.genes.col),
                highlight.feats.size = isolate(input$hl.genes.size),
                highlight.feats.opac = isolate(input$hl.genes.opa),
                highlight.feats.linecolor = isolate(input$hl.genes.lcol),
                highlight.feats.linewidth = isolate(input$hl.genes.lw),
                highlight.featsets.color = isolate(input$hl.genesets.col),
                highlight.featsets.size = isolate(input$hl.genesets.size),
                highlight.featsets.opac = isolate(input$hl.genesets.opa),
                highlight.featsets.linecolor = isolate(input$hl.genesets.lcol),
                highlight.featsets.linewidth = isolate(input$hl.genesets.lw)
            )

            plot_store$volcano_plot
        })

        output[["res_table_full"]] <- DT::renderDT(server = FALSE, {
            req(ress)
            df <- as.data.frame(ress())

            # Collect proper output columns.
            display_cols <- c(abundance.col, lfc.col)
            cnames <- NULL

            for (term in c("lfcSE", "stat", "pvalue", "padj", "svalue", "qvalue", "FDR")) {
                if (term %in% colnames(df)) {
                    cnames <- c(cnames, term)
                }
            }

            # Can't tack on to cnames because of rounding.
            snames <- NULL
            if (!is.null(swap.rownames)) {
                snames <- "ORIGINAL_ROWS"
            }

            all_cols <- c(display_cols, cnames, snames)
            df <- df[, all_cols]

            DT::formatRound(
                DT::datatable(df,
                    rownames = TRUE,
                    filter = "top",
                    extensions = c("Buttons"),
                    options = list(
                        lengthMenu = c(5, 10, 25, 50),
                        pageLength = 25,
                        dom = "Blfrtip",
                        buttons = c("copy", "csv", "excel", "pdf", "print")
                    )
                ),
                columns = c(display_cols, cnames), digits = 5
            ) %>%
                DT::formatStyle(1, target = "row", lineHeight = "40%")
        })

        # Metadata table.
        output$metadata <- DT::renderDT(server = FALSE, {
            if (!is.null(metadata)) {
                df <- as.data.frame(metadata)
                DT::datatable(df,
                    filter = "top",
                    extensions = c("Buttons", "Scroller"),
                    options = list(
                        search = list(regex = TRUE),
                        lengthMenu = list(c(10, 25, 50, -1), c("10", "25", "50", "all")),
                        dom = "Blfrtip",
                        buttons = c("copy", "csv", "excel", "pdf", "print"),
                        scrollX = TRUE,
                        deferRender = TRUE,
                        scrollY = 400,
                        scroller = TRUE
                    )
                )
            } else {
                DT::datatable(data.frame(Message = "No metadata provided"))
            }
        })

        # Only remake heatmap on button click.
        observeEvent(input$update,
            {
                req(ress, anno)

                if (is.null(input$anno.vars)) {
                    anno(NULL)
                }

                if (!is.null(input$metadata_rows_all)) {
                    mat <- mat[, input$metadata_rows_all]
                }

                pdf(NULL)
                ht <- .make_heatmap(mat, ress(), anno(), abundance_col_fun, log2fc_col_fun,
                    sig.term = input$sig.term,
                    lfc.col = lfc.col, abundance.col = abundance.col,
                    sig.thresh = as.numeric(input$sig.thresh), base_mean = input$base_mean, log2fc = input$log2fc,
                    row.km = input$row.km, col.km = input$col.km
                )
                dev.off()

                if (!is.null(ht)) {
                    makeInteractiveComplexHeatmap(input, output, session, ht, h.id,
                        brush_action = .brush_action, click_action = .click_action
                    )
                } else {
                    # The ID for the heatmap plot is encoded as @{heatmap_id}_heatmap.
                    output[[paste0(h.id, "_heatmap")]] <- renderPlot({
                        grid.newpage()
                        grid.text("No row exists after filtering.")
                    })
                }
            },
            ignoreNULL = FALSE
        )

        # Download interactive plots as html.
        output$download_plotly_volc <- downloadHandler(
            filename = function() {
                paste("volcanoplot-", Sys.Date(), ".html", sep = "")
            },
            content = function(file) {
                # export plotly html widget as a temp file to download.
                saveWidget(jqui_resizable(plot_store$volcano_plot),
                    file,
                    selfcontained = TRUE
                )
            }
        )

        output$download_plotly_ma <- downloadHandler(
            filename = function() {
                paste("maplot-", Sys.Date(), ".html", sep = "")
            },
            content = function(file) {
                # export plotly html widget as a temp file to download.
                saveWidget(jqui_resizable(plot_store$ma_plot),
                    file,
                    selfcontained = TRUE
                )
            }
        )
    }

    shinyApp(ui, server, options = list(height = height))
}
