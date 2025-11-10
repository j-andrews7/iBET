# GitHub Copilot Instructions for iBET

## Project Overview
iBET (interactive Bioinformatics Exploratory Tools) is an R package providing drop-in Shiny widgets for common bioinformatics analyses. The package enables interactive exploration of genomics data through visualizations built with PCAtools, DESeq2, and correlation analyses.

## Technology Stack
- **Language**: R (>= 4.0)
- **Framework**: Shiny for interactive web applications
- **Key Dependencies**: 
  - DESeq2, PCAtools for bioinformatics analysis
  - ComplexHeatmap, InteractiveComplexHeatmap for visualizations
  - plotly, ggplot2 for plotting
  - shinydashboard, shinyWidgets for UI components
  - DT for data tables

## Code Style and Conventions

### R Coding Standards
- Use `<-` for assignment, not `=`
- Use `snake_case` for function names and variables
- Use `.` separators in function names where appropriate (e.g., `shinyPCAtools`)
- Indent with 2 spaces, not tabs
- Keep line length under 100 characters when possible
- Use meaningful variable names that describe the data

### Function Documentation
- All exported functions MUST have complete roxygen2 documentation with:
  - `@title` and brief description
  - `@details` for extended explanation
  - `@param` for all parameters with clear descriptions
  - `@return` describing what the function returns
  - `@examples` with working code examples (wrapped in `\dontrun{}` if they require data)
  - `@seealso` for related functions
  - `@author` attribution
  - `@export` tag for exported functions
  - `@importFrom` or `@import` for dependencies

### Import Statements
- Use `@importFrom` for specific functions when only a few are needed
- Use `@import` for packages where many functions are used
- Always declare imports in roxygen2 comments, not in NAMESPACE directly
- Common imports pattern:
  ```r
  #' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
  #' @import DT
  #' @importFrom plotly ggplotly plotlyOutput renderPlotly
  ```

## Shiny Application Structure

### Widget Development Pattern
The package follows a consistent pattern for Shiny widgets:
1. Main exported function (e.g., `shinyPCAtools()`, `shinyDESeq2()`)
2. UI definition using `shinydashboard` components
3. Server logic with reactive expressions
4. Helper utility functions in separate `*_utils.R` files

### UI Guidelines
- Use `dashboardPage`, `dashboardBody` for layout
- Apply custom CSS for styling within `tags$head(tags$style(HTML(...)))`
- Use `shinydashboard` components: `box()`, `tabBox()`, `valueBox()`
- Implement `shinyjqui::jqui_resizable()` for resizable elements
- Use `shinycssloaders::withSpinner()` for loading indicators
- Apply `shinyBS` for tooltips (`tipify()`) and popovers (`popify()`)

### Server Guidelines
- Group related reactive values in `reactiveValues()` objects
- Use `observe()` for side effects
- Use `reactive()` for computed values
- Implement proper error handling with `validate()` and `need()`
- Use `isolate()` to prevent unnecessary reactivity

## Data Handling

### Input Requirements
- Matrix data: features as rows, samples as columns
- Metadata: data.frame with rownames matching matrix column names
- Results data: data.frames or DESeqResults objects
- Gene sets: named lists where names are set names and values are gene vectors

### Data Processing
- Always check for and remove features with zero variance before PCA
- Use `matrixStats::rowVars()` for efficient variance calculations
- Handle missing values explicitly with appropriate warnings
- Validate input dimensions match between matrices and metadata

## Testing and Building

### Package Building
- Use `devtools::document()` to update documentation
- Use `devtools::check()` for R CMD check
- Ensure all examples run without errors (or are wrapped in `\dontrun{}`)
- Version follows semantic versioning (currently 0.0.0.9000 for development)

### Dependencies
- Only add dependencies that are essential
- Prefer Bioconductor packages for bioinformatics functionality
- Use Imports for required packages
- Use Suggests for optional/example packages
- Keep the DESCRIPTION file up to date

## Visualization Standards

### Plotly Integration
- Use `plotly::ggplotly()` to convert ggplot2 objects
- Use `plotly::toWebGL()` for better performance with large datasets
- Implement custom hover text with biologically relevant information
- Make plots interactive with click, hover, and selection events

### Color Schemes
- Use `dittoSeq::dittoColors()` for categorical variables
- Provide color customization options via `colourpicker`
- Ensure colorblind-friendly palettes when possible

### Plot Types
- Scatter plots: PCA biplots, MA plots, volcano plots
- Heatmaps: Use `ComplexHeatmap::Heatmap()` or `pheatmap()`
- Interactive plots: Leverage plotly for all primary visualizations

## File Organization
```
iBET/
├── R/                      # R source code
│   ├── PCAtools.R         # Main PCA widget function
│   ├── PCA_utils.R        # PCA helper functions
│   ├── DESeq2.R           # Main DESeq2 widget function
│   ├── DESeq2_utils.R     # DESeq2 helper functions
│   ├── DECorr.R           # Main DE correlation widget
│   ├── DECorr_utils.R     # DE correlation helpers
│   ├── plot_utils.R       # Shared plotting utilities
│   └── data_utils.R       # Shared data utilities
├── man/                    # Auto-generated documentation
├── vignettes/             # Package vignettes
│   └── iBET.Rmd           # Main package vignette
└── DESCRIPTION            # Package metadata
```

## Common Patterns

### Parameter Validation
```r
if (!is.matrix(mat) && !is.data.frame(mat)) {
  stop("mat must be a matrix or data.frame")
}

if (!all(colnames(mat) %in% rownames(metadata))) {
  stop("All matrix column names must be in metadata rownames")
}
```

### Reactive UI Generation
```r
output$dynamic_ui <- renderUI({
  req(input$some_input)
  # Generate UI based on reactive inputs
})
```

### Safe Data Access
```r
observeEvent(input$button, {
  req(data_available())
  validate(
    need(nrow(dataset()) > 0, "No data available")
  )
  # Proceed with computation
})
```

## Special Considerations

### Bioinformatics Context
- Gene identifiers can be Ensembl IDs, gene symbols, or other formats
- Fold change values are typically log2 transformed
- P-values should be adjusted for multiple testing (FDR, Bonferroni)
- Sample sizes in bioinformatics are often small (n=3-20)

### Performance
- Large datasets (>10,000 features) require optimization
- Use `removeVar` parameter to filter low-variance features
- Consider data.table for large data manipulations
- Use `toWebGL()` for large plotly graphics

### Interactivity Features
- Click-to-label genes on plots (using plotly events)
- Draggable labels with `shinyjqui`
- Dynamic filtering with sliders and checkboxes
- Export functionality for plots and tables

## Error Messages
- Provide clear, actionable error messages
- Include information about expected data formats
- Suggest solutions when validation fails
- Use informative warnings for non-critical issues

## Version Control
- Commit messages should be clear and descriptive
- Update version number in DESCRIPTION for releases
- Update NEWS.md with changes for each version
- Tag releases with semantic version numbers
