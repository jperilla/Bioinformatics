
data = read.table("c:\\temp\\spellman.txt", header=T, row.names=1)
data[1:4, 1:4]

# show dimensions of data
dim(data)

# show names of samples
dimnames(data)[[2]]

# get only cdc15 sample (23-46)
cdc15Data = data[,23:46]
dimnames(cdc15Data)
cdc15Data[1:4,]

#create pearson's correlation
install.packages('gplots')
library(gplots)
help(cor)
# method = pearson is the default
data.cdc15cor <- cor(cdc15Data, use="pairwise.complete.obs") # only use pairwise observations where there is data on both side
data.cdc15cor

help(legend)
#plot the correlation matrix
help(heatmap)
rc <- rainbow(nrow(data.cdc15cor), start=0, end=.3)
colors <- colorRampPalette(c("black", "blue", "yellow"))(20)
heatmap(x = data.cdc15cor, col = colors, symm = TRUE, 
          main = "Pearson's correlation - cdc15",
          xlab="cdc15 expression over time", ylab="cdc15 expression over time",
          margins = c(7,7), Rowv=NA, Colv=NA)

install.packages('heatmap3')
library(heatmap3)

heatmap3(x = data.cdc15cor, col = colors, symm = TRUE, 
         main = "Pearson's correlation - cdc15",
         xlab="cdc15 expression over time", ylab="cdc15 expression over time",
         margins = c(8,8), Rowv=NA, Colv=NA)


yal002w <- cdc15Data[2,]
yal002w.values <- yal002w[!is.na(yal002w)]
yal002w.mean <- mean(yal002w_values) 
yal002w[is.na(yal002w)] <- yal002w.mean
yal002w

length(yal002w)
plot(range(1:24), range(yal002w), type='n',
     main="YAL002w profile plot", xlab="samples", 
     ylab="intensities")
lines(c(1:24), yal002w)
grid(col='gray')

dimnames(data.cdc15cor)[[1]]

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
        
        
help(renderPlot)
