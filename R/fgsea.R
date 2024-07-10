#' Create an interactive Shiny app for visualization & exploration of fgsea analyses
#'
#' This function creates a Shiny app to explore analyis results returned by \code{\link[fgsea]{fgsea}}.
#' It allows the creation of interactive overview and rank plots for specific gene sets. 
#' If RNA-seq data is provided, it also allows for visualization of the leading edge genes in a heatmap.
#' If differential expression results are provided, it allows for incorporation of effect size and significance
#' info into the plots.
#'
#' @param fgsea.res A named list containing \code{\link[fgsea]{fgsea}} results.
#'   Multiple data.frames may be provided, one per element of the list.
#'   Users will be able to swap between them within the app. 
#'   Alternatively, a named list of named lists may be provided, for which the names of the first level of the list
#'   correspond to a given fgsea analysis, and the names of the second level correspond to the fgsea results for a given gene set collection,
#'   e.g. "C5.GOBP".
#' @param ranked.genes A named list containing a named vector of values by which genes should be ranked (as used by fgsea).
#'   If more than one is provided, \code{fgsea.res} must be a named list of named lists,
#'   for which the names of the \code{ranked.genes} list must match the names of the first level of the \code{fgsea.res} list.
#' @param genesets A named list containing gene sets as character vectors of gene identifiers.
#' @param se An optional \code{\link[RangedSummarizedExperiment]{RangedSummarizedExperiment}}-based object containing RNA-seq data.
#'   Used for additional visualizations and providing sample metadata.
#' @param de.results Optional named list of data.frames containing differential expression results. 
#'   If more than one is provided, \code{fgsea.res} must be a named list of named lists,
#'   for which the names of the \code{de.results} list must match the names of the first level of the \code{fgsea.res} list.
#' @param h.id String indicating unique ID for interactive plots.
#'   Required if multiple apps are run within the same Rmd file.
#'
#' @return A Shiny app containing interactive visualizations of \code{\link[fgsea]{fgsea}} results.
#'
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#'
#' @author Jared Andrews
#' @export
#' @examples
#' library(iBET)
#'
#' genes <- list(ESC = d1.genes, plasmid = d2.genes)
#' sgrnas <- list(ESC = d1.sgrnas, plasmid = d2.sgrnas)
#'
#' app <- shinyfgsea(
#'     fgsea.res = fgsea.res, ranked.genes = ranked.genes,
#'     genesets = genesets
#' )
#' if (interactive()) {
#'     shiny::runApp(app)
#' }
shinyfgsea <- function(fgsea.res,
                       ranked.genes,
                       genesets,
                       se = NULL,
                       de.results = NULL,
                       h.id = "sfgsea") {

    # Set initial metadata and dataset choices if input data isn't NULL.
    fgsea.sets <- NULL
    sgrna.choices <- NULL
    geneset.choices <- names(genesets)
    sgrna.gene <- NULL

    if (!is.null(gene.data)) {
        gene.choices <- names(gene.data)
    }

    ui <- navbarPage(
        title = div(a(img(src = "logo/logo.png", height = "50"),
            href = "https://bioconductor.org/packages/IBET"
        ), "IBET"),
        header = list(
            useShinyjs(),
            css,
            tags$head(tags$link(rel = "shortcut icon", href = "logo/logo.png"))
        ),
        ## ---------------Overview tab-----------------
        .create_fgsea_tab_overview(gsea.choices),
        ## ---------------Metadata tab-----------------
        if (!is.null(se)) {
            .create_fgsea_tab_metadata()
        }
    )


    server <- function(input, output, session) {
        ## -------------Reactive Values---------------
        robjects <- reactiveValues(
            fgsea.res = fgsea.res,
            ranked.genes = ranked.genes,
            se = se,
            de.results = de.results,
            genesets = genesets,
            h.id = h.id,
            clicked.swoosh = NULL,
            clicked.ov.le = NULL,
            plot.ov.bar = NULL,
            plot.ov.swoosh = NULL,
            plot.ov.le = NULL
        )

        # Create downloadHander outputs.
        .create_fgsea_dl_outputs(output, robjects)

        ## --------------Disable Inputs-----------------
        # Disable certain inputs if no data is provided.
        .create_fgsea_ui_observers(robjects)

        ## ---------Overview & Summary Tables Tabs-------------
        # Load the gene summaries for easy plotting.
        .create_ov_observers(input, robjects)

        # Summary tables and plots.
        .create_ov_outputs(input, output, robjects)

        # This ensures the rank options are updated even when initially hidden in the collapsible panel.
        # outputOptions(output, "ov.term.options", suspendWhenHidden = FALSE)

        ## ----------------Metadata Tab-----------------
        .create_fgsea_metadata_outputs(input, output, robjects)
    }

    shinyApp(ui, server)
}
