.vline <- function(x = 0, color = "red", width = 1, dash = "solid") {
  list(
    type = "line",
    y0 = 0,
    y1 = 1,
    yref = "paper",
    x0 = x,
    x1 = x,
    line = list(color = color, width = width, dash = dash)
  )
}

.hline <- function(y = 0, color = "blue", width = 1, dash = "solid") {
  list(
    type = "line",
    x0 = 0,
    x1 = 1,
    xref = "paper",
    y0 = y,
    y1 = y,
    line = list(color = color, width = width, dash = dash)
  )
}

.fitline <- function(df, color = "black", width = 0.75, dash = "solid") {
  list(
    type = "line",
    line = list(color = color, width = width, dash = dash),
    xref = "x",
    yref = "y",
    y0 = min(df$fv),
    y1 = max(df$fv),
    x0 = df$lfc.x[df$fv == min(df$fv)],
    x1 = df$lfc.x[df$fv == max(df$fv)]
  )


}
