#' @title Interactive Shiny Enrichment Dashboard
#' @description Creates interactive Shiny dashboard for enrichment analysis visualization with Dot Plot, Bar Plot, and Heatmap tabs
#' @param enrichment.obj Enrichment result object (e.g., from clusterProfiler enrichGO/enrichKEGG)
#' @param sim.matrix Reduced enrichment data frame from the rrvgo package using calculatesim.matrix() and reduceSimMatrix() to produce Tree plot
#' @return Shiny application object (interactive dashboard)
#' @importFrom shiny shinyApp dashboardPage dashboardHeader dashboardSidebar dashboardBody
#' @importFrom shiny tabsetPanel tabPanel numericInput selectInput textInput colourInput
#' @importFrom shiny renderPlotly plotlyOutput bsCollapse bsCollapsePanel
#' @importFrom shinydashboard dashboardPage dashboardHeader dashboardSidebar dashboardBody
#' @importFrom shinyBS bsCollapse bsCollapsePanel
#' @importFrom shinyjs useShinyjs
#' @importFrom shinyjqui jqui_resizable
#' @importFrom plotly plotlyOutput
#' @importFrom dashboardthemes shinyDashboardThemes
#' @importFrom bslib bs_theme
#' @importFrom clusterProfiler enrichGO enrichKEGG
#' @author Jacob Martin
#' @export
shinyEnrichment1 <- function(enrichment.obj, sim.matrix = NULL) {
                
    ui <- dashboardPage(
            dashboardHeader(disable = TRUE),
            dashboardSidebar(
            bsCollapse(
                id = "sidebar.collapse", open = "DotPlot.settings",
                # Dot Plot Button UI
                bsCollapsePanel(
                    title = span(icon("plus"), "Dot Plot"),
                    value = "DotPlot.settings", style = "success",
                    selectInput("colour.by", "Colour by:", c("p.adjust", "pvalue", "qvalue")),
                    numericInput("num.categories.dot", "Number of Gene Sets:", value = 10, min = 1, max = 50),
                    textInput("title.dot.plot", "Enter the title for the plot:"),
                    selectInput("x.axis.selection", "X-axis:", c("Count", "GeneRatio")),
                    numericInput("dot.text.size", "Y Axis Text Size: ", value = 10, min = 1, max = 50),
                    colourInput("my.color.dot1", "Select a color", value = "red"),
                    colourInput("my.color.dot2", "Select a color", value = "blue")
                ),
                # Bar Plot Button UI
                bsCollapsePanel(
                    title = span(icon("plus"), "Bar Plot"),
                    value = "BarPlot.settings", style = "warning",
                    numericInput("num.categories.bar", "Categories:", value = 20, min = 1, max = 50),
                    selectInput("colour.by.bar", "Colour by:", c("p.adjust", "pvalue", "qvalue")),
                    textInput("title.bar.plot", "Enter the title for the plot:"),
                    selectInput("x.axis.selectionBar", "X-axis:", c("Count", "GeneRatio")),
                    numericInput("bar.text.size", "Y Axis Text Size: ", value = 10, min = 1, max = 50),
                    colourInput("my.color.bar1", "Select a color", value = "red"),
                    colourInput("my.color.bar2", "Select a color", value = "blue")
                ),
                #Heat PLot
                bsCollapsePanel(
                    title = span(icon("plus"), "HeatMap Plot"),
                    value = "HeatPlot.settings", style = "primary",
                    numericInput("heat.categories", "Categories:", value = 10, min = 1, max = 50),
                    colourInput("colour.heat", "Select a color", value = "blue"),
                    numericInput("text.heat", "Font Size of Y axis: ", value = 10, min = 1, max = 50),
                    textInput("title.heat", "Enter the title for the plot:")
                ),
                bsCollapsePanel(
                    title = span(icon("plus"), "Tree Plot"),
                    value = "TreePlot.settings", style = "primary",
                    uiOutput("tree.controls")
                )

            )
        ),
        dashboardBody(
            tags$head(tags$style(HTML("
            .panel-body { padding: 10px; }
            .form-control { font-size: 10px; height: 24px; }
            label { font-size: 80%; }"
            ))),
            useShinyjs(),
            shinyDashboardThemes(theme = "onenote"),
            tabsetPanel(
                tabPanel(
                    "Dot Plot",
                    jqui_resizable(
                        plotlyOutput("dotplot", height = "100%"),
                            options = list(
                            handles = "all",  
                            minWidth = 300,
                            minHeight = 300,
                            maxWidth = 1200,
                            maxHeight = 800
                        )
                    )
                ),
                tabPanel(
                    "BarPlot",
                    jqui_resizable(
                        plotlyOutput("barPlot", height = "100%"),
                            options = list(
                            handles = "all",  
                            minWidth = 300,
                            minHeight = 300,
                            maxWidth = 1200,
                            maxHeight = 800
                        )
                    )                    
                ),
                tabPanel(
                    "HeatMap",
                    jqui_resizable(
                        plotlyOutput("heatPlot", height = "100%"),
                            options = list(
                            handles = "all",  
                            minWidth = 300,
                            minHeight = 300,
                            maxWidth = 1200,
                            maxHeight = 800
                        )
                    )                    
                ),
                tabPanel(
                    "TreePlot",
                    uiOutput("tree.map.plot")
                )
            )                    
                               
        )

    )    
    server <- function(input, output, session) {
        #Ui Condtion tree.sidebar: 
        output$tree.controls <- renderUI(
            if (is.null(sim.matrix)) {
                tags$div(style = "color:#888; font-size: 90%;",
                    "Tree map unavailable: no similarity matrix provided.")
            } else {
                tagList(
                    textInput("titleTree", "Enter the title for the dot plot:"),
                    selectInput("colour_palette", "RColorBrewerPalette:",
                    c("Blues","BuGn","BuPu","GnBu","Greens","Greys",
                        "Oranges","OrRd","PuBu","PuBuGn","PuRd","Purples",
                        "RdPu","Reds","YlGn","YlGnBu","YlOrBr","YlOrRd","Viridis",
                        "BrBG","PiYG","PRGn","PuOr","RdBu","RdGy","RdYlBu","RdYlGn","Spectral",
                        "Accent","Dark2","Paired","Pastel1","Pastel2","Set1","Set2","Set3"))
                )
            }               
        )
        output$tree.map.plot <- renderUI(
            if (is.null(sim.matrix)) {
                tags$div(style = "color:#888; font-size: 90%;",
                    "Tree map unavailable: no similarity matrix provided.")
            } else {
                plotOutput("TreePlot", height = "700px", width = "100%")
            }
        )
        output$dotplot <- renderPlotly(
            make_dot_plot(enrich = enrichment.obj, num.sets = input$num.categories.dot, colour.by = input$colour.by, title.for.plot = input$title.dot.plot, x.axis = input$x.axis.selection, text.size = input$dot.text.size, colour1 = input$my.color.dot1, colour2 = input$my.color.dot2)
        )
        output$barPlot <- renderPlotly(
            make_bar_plot(enrich = enrichment.obj, num.sets = input$num.categories.bar, colour.by = input$colour.by.bar, title.for.plot = input$title.bar.plot, x.axis = input$x.axis.selectionBar, colour1 = input$my.color.bar1, colour2 = input$my.color.bar2, text.size = input$bar.text.size)
        )
        output$heatPlot <- renderPlotly(
            make_binary_heatmap(enrich = enrichment.obj, num.sets = input$heat.categories, colour = input$colour.heat, size.text = input$text.heat, title = input$title.heat)
        )
        output$TreePlot <- renderPlot(
            treemap::treemap(sim.matrix, index = c("parentTerm", "term"), vSize = "score", type = "index",
                    title = input$titleTree, palette = input$colour_palette, fontcolor.labels = c("#FFFFFFDD", "#00000080"), 
                    bg.labels = 0, border.col = "#00000080")
        )
    }
    shinyApp(ui = ui, server = server)
}
