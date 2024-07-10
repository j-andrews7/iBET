#' Create a tabPanel for the metadata tab
#'
#' Create a \code{\link{tabPanel}} with UI elements for the metadata tab.
#'
#' @return
#' A \code{\link{tabPanel}} with UI elements for the metadata tab.
#'
#' @author Jared Andrews
#'
#' @importFrom DT DTOutput
#' @importFrom shiny tabPanel br div
#' @importFrom shinycssloaders withSpinner
#'
#' @rdname INTERNAL_create_tab_metadata
.create_fgsea_tab_metadata <- function() {
    # nocov start
    tabPanel(
        title = "Metadata",
        id = "metadata",
        br(),
        div(withSpinner(DTOutput("fgsea.metadata")), style = "font-size:80%;")
    )
    # nocov end
}