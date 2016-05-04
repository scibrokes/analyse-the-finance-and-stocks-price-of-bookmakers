require(googleVis)
require(shiny)

shinyUI(pageWithSidebar(
  headerPanel("", windowTitle="Example googleVis with interaction"),
  sidebarPanel(
    selectInput("dataset", "Choose a dataset:", 
                choices = c("pressure", "cars")),
    uiOutput("selectedOut")
  ),
  mainPanel(
    tabsetPanel(
      tabPanel("Main",
               htmlOutput("view"),
               plotOutput("distPlot", width="300px", height="200px")),
      tabPanel("About", includeMarkdown('README.md')
      ))))
)

