#' @title Interactive Shiny Enrichment Dashboard
#' @description Creates interactive Shiny dashboard for enrichment analysis visualization with Dot Plot, Bar Plot, and Heatmap tabs
#' @param enrichmentObj Enrichment result object (e.g., from clusterProfiler enrichGO/enrichKEGG)
#' @param simMatrix Reduced enrichment data frame from the rrvgo package using calculateSimMatrix() and reduceSimMatrix() to produce Tree plot
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
shinyEnrichment <- function(enrichmentObj, simMatrix = NULL) {
    sidebar_with_tree <- dashboardSidebar(
            bsCollapse(
                id = "sidebar_collapse", open = "DottPlot.settings",
                # Dot Plot Button UI
                bsCollapsePanel(
                    title = span(icon("plus"), "Dot Plot"),
                    value = "DottPlot.settings", style = "success",
                    selectInput("colour_by", "Colour by:", c("p.adjust", "pvalue", "qvalue")),
                    numericInput("num_categoriesDot", "Number of Gene Sets:", value = 10, min = 1, max = 50),
                    textInput("titleDotPlot", "Enter the title for the plot:"),
                    selectInput("x_axisSelection", "X-axis:", c("Count", "GeneRatio")),
                    numericInput("dot_textSize", "Y Axis Text Size: ", value = 10, min = 1, max = 50),
                    colourInput("my_color_Dot1", "Select a color", value = "red"),
                    colourInput("my_color_Dot2", "Select a color", value = "blue")
                ),
                # Bar Plot Button UI
                bsCollapsePanel(
                    title = span(icon("plus"), "Bar Plot"),
                    value = "BarPlot.settings", style = "warning",
                    numericInput("num_categories_Bar", "Categories:", value = 20, min = 1, max = 50),
                    selectInput("colour_by_Bar", "Colour by:", c("p.adjust", "pvalue", "qvalue")),
                    textInput("titleBarPlot", "Enter the title for the plot:"),
                    selectInput("x_axisSelectionBar", "X-axis:", c("Count", "GeneRatio")),
                    numericInput("bar_textSize", "Y Axis Text Size: ", value = 10, min = 1, max = 50),
                    colourInput("my_color_Bar1", "Select a color", value = "red"),
                    colourInput("my_color_Bar2", "Select a color", value = "blue")
                ),
                #Heat PLot
                bsCollapsePanel(
                    title = span(icon("plus"), "HeatMap Plot"),
                    value = "HeatPlot.settings", style = "primary",
                    numericInput("heat_categories", "Categories:", value = 10, min = 1, max = 50),
                    colourInput("colourHeat", "Select a color", value = "blue"),
                    numericInput("textHeat", "Font Size of Y axis: ", value = 10, min = 1, max = 50),
                    textInput("titleHeat", "Enter the title for the plot:")
                ),
                bsCollapsePanel(
                    title = span(icon("plus"), "Tree Plot"),
                    value = "TreePlot.settings", style = "primary",
                    textInput("titleTree", "Enter the title for the dot plot:"),
                    selectInput("colour_palette", "RColorBrewerPalette: ", 
                    c("Blues", "BuGn", "BuPu", "GnBu", "Greens", "Greys", 
                      "Oranges", "OrRd", "PuBu", "PuBuGn", "PuRd", "Purples", 
                      "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd", "Viridis",
                      "BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral",
                      "Accent", "Dark2", "Paired", "Pastel1", "Pastel2", "Set1", "Set2", "Set3"))
                )
            )
        )
        sidebar_without_tree <- dashboardSidebar(
            bsCollapse(
                id = "sidebar_collapse", open = "DottPlot.settings",
                # Dot Plot Button UI
                bsCollapsePanel(
                    title = span(icon("plus"), "Dot Plot"),
                    value = "DottPlot.settings", style = "success",
                    selectInput("colour_by", "Colour by:", c("p.adjust", "pvalue", "qvalue")),
                    numericInput("num_categoriesDot", "Number of Gene Sets:", value = 10, min = 1, max = 50),
                    textInput("titleDotPlot", "Enter the title for the plot:"),
                    selectInput("x_axisSelection", "X-axis:", c("Count", "GeneRatio")),
                    numericInput("dot_textSize", "Y Axis Text Size: ", value = 10, min = 1, max = 50),
                    colourInput("my_color_Dot1", "Select a color", value = "red"),
                    colourInput("my_color_Dot2", "Select a color", value = "blue")
                ),
                # Bar Plot Button UI
                bsCollapsePanel(
                    title = span(icon("plus"), "Bar Plot"),
                    value = "BarPlot.settings", style = "warning",
                    numericInput("num_categories_Bar", "Categories:", value = 20, min = 1, max = 50),
                    selectInput("colour_by_Bar", "Colour by:", c("p.adjust", "pvalue", "qvalue")),
                    textInput("titleBarPlot", "Enter the title for the plot:"),
                    selectInput("x_axisSelectionBar", "X-axis:", c("Count", "GeneRatio")),
                    numericInput("bar_textSize", "Y Axis Text Size: ", value = 10, min = 1, max = 50),
                    colourInput("my_color_Bar1", "Select a color", value = "red"),
                    colourInput("my_color_Bar2", "Select a color", value = "blue")
                ),
                #Heat PLot
                bsCollapsePanel(
                    title = span(icon("plus"), "HeatMap Plot"),
                    value = "HeatPlot.settings", style = "primary",
                    numericInput("heat_categories", "Categories:", value = 10, min = 1, max = 50),
                    colourInput("colourHeat", "Select a color", value = "blue"),
                    numericInput("textHeat", "Font Size of Y axis: ", value = 10, min = 1, max = 50),
                    textInput("titleHeat", "Enter the title for the plot:")
                )
            )
        )

    tab_without_tree <- dashboardBody(
                    tags$head(tags$style(HTML("
                .panel-body { padding: 10px; }
                .form-control { font-size: 10px; height: 24px; }
                label { font-size: 80%; }
                    "))),
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
                        )       
                    )
                )
    tab_with_tree <- dashboardBody(
                    tags$head(tags$style(HTML("
                .panel-body { padding: 10px; }
                .form-control { font-size: 10px; height: 24px; }
                label { font-size: 80%; }
                "))),
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

                            plotOutput("TreePlot", height = "700px", width = "100%"),

                        )
                    )                    
                               
                )
                
    ui <- dashboardPage(
            dashboardHeader(disable = TRUE),
                if (is.null(simMatrix)) sidebar_without_tree else sidebar_with_tree, # Condtion of simMatrix to determine whether Tree Panel will be visible
                if (is.null(simMatrix)) tab_without_tree else tab_with_tree # Condtion for tree tab

            )    
    server <- function(input, output, session) {
        output$dotplot <- renderPlotly(
            dotPlotEnrichment(enrich = enrichmentObj, numSets = input$num_categoriesDot, colourBy = input$colour_by, titleForPlot = input$titleDotPlot, x_axis = input$x_axisSelection, textSize = input$dot_textSize, colour1 = input$my_color_Dot1, colour2 = input$my_color_Dot2)
        )

        output$barPlot <- renderPlotly(
            barPlotEnrichment(enrich = enrichmentObj, numSets = input$num_categories_Bar, colourBy = input$colour_by_Bar, titleForPlot = input$titleBarPlot, x_axis = input$x_axisSelectionBar, colour1 = input$my_color_Bar1, colour2 = input$my_color_Bar2, textSize = input$bar_textSize)
        )
        output$heatPlot <- renderPlotly(
            heatPlotEnrichment(enrich = enrichmentObj, numSets = input$heat_categories, colour = input$colourHeat, sizeText = input$textHeat, title = input$titleHeat)
        )
        output$TreePlot <- renderPlot(
            treeEnrichment(reduceTerms = simMatrix, title = input$titleTree, colourPalette = input$colour_palette)
        )
    }

    shinyApp(ui = ui, server = server, options = list(launch.browser = TRUE))
}
