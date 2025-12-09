library(devtools)
source("/Volumes/JM/Development/iBET/R/shinyEnrichment_utils.R")
library(shiny)
library(shinydashboard)
library(plotly)
library(shinyBS)
library(shinyjs)
library(dashboardthemes)
library(bslib)
library(shinyjqui)
library(clusterProfiler)
library(ggtree)
library(enrichplot)
library(org.Hs.eg.db)
library(GOSemSim)

shinyEnrichment <- function(enrichmentObj) {
    ui <- dashboardPage(
        dashboardHeader(disable = TRUE),
        dashboardSidebar(
            # Filtering
            bsCollapse(
                id = "sidebar_collapse", open = "general.settings",
                # Dot Plot Button UI
                bsCollapsePanel(
                    title = span(icon("plus"), "Dot Plot"),
                    value = "DottPlot.settings", style = "success",
                    selectInput("colour_by", "Colour by:", c("p.adjust", "pvalue", "qvalue")),
                    numericInput("num_categoriesDot", "Number of Gene Sets:", value = 10, min = 1, max = 50),
                    textInput("titleDotPlot", "Enter the title for the dot plot:"),
                    selectInput("x_axisSelection", "X-axis:", c("Count", "GeneRatio"))
                ),
                # Bar Plot Button UI
                bsCollapsePanel(
                    title = span(icon("plus"), "Bar Plot"),
                    value = "BarPlot.settings", style = "warning",
                    numericInput("num_categories_Bar", "Categories:", value = 20, min = 1, max = 50),
                    selectInput("colour_by_Bar", "Colour by:", c("p.adjust", "pvalue", "qvalue")),
                    textInput("titleBarPlot", "Enter the title for the dot plot:"),
                    selectInput("x_axisSelectionBar", "X-axis:", c("Count", "GeneRatio"))
                ),
                #Heat PLot
                bsCollapsePanel(
                    title = span(icon("plus"), "HeatMap Plot"),
                    value = "HeatPlot.settings", style = "primary",
                    numericInput("heat_categories", "Categories:", value = 10, min = 1, max = 50)
                )
            )
        ),
        dashboardBody(
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
                    "Tree Plot",
                    jqui_resizable(
                        plotOutput("treePlot", height = "100%"),
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
    )    
    server <- function(input, output, session) {
        # Creating Dot Plot Output:
        # DotPlot Output:
        enrich_with_termsim <- pairwise_termsim(enrichmentObj, 
                                       method = "Wang", 
                                       semData = godata("org.Hs.eg.db", ont = "BP"))
        output$dotplot <- renderPlotly(
            dotPlotEnrichment(enrich = enrichmentObj, numSets = input$num_categoriesDot, colourBy = input$colour_by, titleForPlot = input$titleDotPlot, x_axis = input$x_axisSelection)
        )

        output$barPlot <- renderPlotly(
            barPlotEnrichment(enrichObj = enrichmentObj, numSets = input$num_categories_Bar, colourBy = input$colour_by_Bar, titleForPlot = input$titleBarPlot, x_axis = input$x_axisSelectionBar)
        )
        output$heatPlot <- renderPlotly(
            heatPlotEnrichment(enrich = enrichmentObj, numSets = input$heat_categories)
        )
        output$treePlot <- renderPlot(
            treePlotEnrichment(enrich = enrich_with_termsim, numSets = 20)
        )
    }

    shinyApp(ui = ui, server = server, options = list(launch.browser = TRUE))
}
test <- readRDS("/Volumes/JM/Development/iBET/Play/ExampleEnrichment.rds")
app <- shinyEnrichment(enrichmentObj = test)
runApp(app)
