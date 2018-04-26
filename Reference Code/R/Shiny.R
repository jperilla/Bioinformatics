
install.packages('shiny')
library(shiny)

# Shiny app template
ui <- fluidPage()

server <- function(input, output) {}

shinyApp(ui = ui, server = server)

# Shiny app
ui <- fluidPage(
  headerPanel('cdc15 time series'),
  sidebarLayout(
    sidebarPanel(
      selectInput('xcol', 'X Variable', dimnames(data.cdc15cor)[[2]]),
      selectInput('ycol', 'Y Variable', dimnames(data.cdc15cor)[[1]],
                  selected = dimnames(data.cdc15cor)[[2]][2]),
      selectInput('col', 'Point Color', c('red', 'blue', 'green', 'black'))
    ),
    mainPanel(
      plotOutput('plot1')
    )
  )
)


server <- function(input, output) {
  selectedData <- reactive({
    data.cdc15cor[, c(input$xcol, input$ycol)]
  })
  output$plot1 <- renderPlot({plot(data.cdc15cor[, c(input$xcol, input$ycol)], col=input$col)})
}

shinyApp(ui = ui, server = server)