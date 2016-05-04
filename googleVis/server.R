require(googleVis)
require(shiny)

shinyServer(function(input, output) {
  datasetInput <- reactive({
    switch(input$dataset, "pressure" = pressure, "cars" = cars)
  })
  output$view <- renderGvis({    
    jscode <- "var sel = chart.getSelection();
    var row = sel[0].row;
    var text = data.getValue(row, 1);               
    $('input#selected').val(text);
    $('input#selected').trigger('change');"    
    gvisScatterChart(data.frame(datasetInput()),
                     options=list(gvis.listener.jscode=jscode,
                                  height=200, width=300))
    
  })
  
  output$distPlot <- renderPlot({
    if (is.null(input$selected))
      return(NULL)
    
    dist <- rnorm(input$selected)
    hist(dist,main=input$selected)
  })
  
  output$selectedOut <- renderUI({
    textInput("selected", "", value="10")
  })
  outputOptions(output, "selectedOut", suspendWhenHidden=FALSE)   
})