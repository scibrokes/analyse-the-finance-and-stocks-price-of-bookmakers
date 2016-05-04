
## Reference : http://stackoverflow.com/questions/17930985/conditional-output-shiny-ui

rm(list = ls())
library(shiny)
runApp(list(
  ui=shinyUI(pageWithSidebar(
    headerPanel("Shiny Example"),
    sidebarPanel(
      wellPanel(
        selectInput(    inputId = "variable1",label = "Select First Variable:", 
                        choices = c("Binary Variable 1" = "binary1",
                                    "Binary Variable 2" = "binary2", 
                                    "Continuous Variable 1" = "cont1",
                                    "Continuous Variable 2" = "cont2"),
                        selected = "Binary Variable 1"
        )
      ),
      wellPanel(
        checkboxInput("bivariate", "Proceed to Bivariate Analysis", FALSE),
        conditionalPanel(
          condition="input.bivariate==true",
          selectInput(inputId = "variable2", 
                      label = "Select Second Variable:",
                      choices = c("Binary Variable 1" = "binary1",
                                  "Binary Variable 2" = "binary2", 
                                  "Continuous Variable 1" = "cont1",
                                  "Continuous Variable 2" = "cont2"),
                      selected = "Binary Variable 2"
          )
        )
      )
    ),
    mainPanel(
      h5("Item Response Rate"),
      verbatimTextOutput("nitem"),
      
      h5(textOutput("caption2")),
      verbatimTextOutput("out2"),
      
      h5(textOutput("caption3")),
      verbatimTextOutput("out3"),
      
      h5(textOutput("caption4")),
      verbatimTextOutput("out4"),
      
      h5(textOutput("caption5")),
      plotOutput("out5")
    )
  )),
  {binary1 <- rbinom(100,1,0.5)
  binary2 <- rbinom(100,1,0.5)
  cont1   <- rnorm(100)
  cont2   <- rnorm(100)
  
  dat <- as.data.frame(cbind(binary1, binary2, cont1, cont2))
  
  dat$binary1 <- as.factor(dat$binary1)
  dat$binary2 <- as.factor(dat$binary2)
  dat$cont1 <- as.numeric(dat$cont1)
  dat$cont2 <- as.numeric(dat$cont2)},
  
  server= shinyServer(function(input, output) {
    inputVar1 <- reactive({
      parse(text=sub(" ","",paste("dat$", input$variable1)))
    })
    
    inputVar2 <- reactive({
      parse(text=sub(" ","",paste("dat$", input$variable2)))
    })
    
    output$nitem <- renderPrint({
      if ( (input$bivariate==FALSE) & (is.factor(eval(inputVar1()))==TRUE) ) {
        n <- sum(table(eval(inputVar1())))
        p <- n/100
        out <- cat(paste(n,gsub(" ","",paste("(",round(as.numeric(p)*100,2),"%",")"))),"\n")
      } else {
        if ( (input$bivariate==FALSE) & (is.numeric(eval(inputVar1()))==TRUE) ) {
          n <- sum(table(eval(inputVar1())))
          p <- n/100
          out <- cat(paste(n,gsub(" ","",paste("(",round(as.numeric(p)*100,2),"%",")"))),"\n")
        } else {
          if ( (input$bivariate==TRUE) & (is.factor(eval(inputVar1()))==TRUE) & (is.factor(eval(inputVar2()))==TRUE) ) {
            n <- sum(table(eval(inputVar1()),eval(inputVar2())))
            p <- n/100
            out <- cat(paste(n,gsub(" ","",paste("(",round(as.numeric(p)*100,2),"%",")"))),"\n")
          } else {
            if ( (input$bivariate==TRUE) & (is.numeric(eval(inputVar1()))==TRUE) & (is.numeric(eval(inputVar2()))==TRUE) ) {
              n <- sum(table(eval(inputVar1()),eval(inputVar2())))
              p <- n/100
              out <- cat(paste(n,gsub(" ","",paste("(",round(as.numeric(p)*100,2),"%",")"))),"\n")
            } else {
              if ( (input$bivariate==TRUE) & (is.factor(eval(inputVar1()))==TRUE) & (is.numeric(eval(inputVar2()))==TRUE) ) {
                n <- sum(table(eval(inputVar1()),eval(inputVar2())))
                p <- n/100
                out <- cat(paste(n,gsub(" ","",paste("(",round(as.numeric(p)*100,2),"%",")"))),"\n")
              } else {
                if ( (input$bivariate==TRUE) & (is.numeric(eval(inputVar1()))==TRUE) & (is.factor(eval(inputVar2()))==TRUE) ) {
                  n <- sum(table(eval(inputVar1()),eval(inputVar2())))
                  p <- n/100
                  out <- cat(paste(n,gsub(" ","",paste("(",round(as.numeric(p)*100,2),"%",")"))),"\n")
                }
              }
            }
          }
        }
      }
      
    })
    
    output$caption2 <- renderText({
      if ( (input$bivariate==FALSE) & (is.factor(eval(inputVar1()))==TRUE) ) {
        caption2 <- "Univariate Table"
      } else {
        if ( (input$bivariate==FALSE) & (is.numeric(eval(inputVar1()))==TRUE) ) {
          caption2 <- "Univariate Summary"
        } else {
          if ( (input$bivariate==TRUE) & (is.factor(eval(inputVar1()))==TRUE) & (is.factor(eval(inputVar2()))==TRUE) ) {
            captions2 <- "Bivariate Table"
          } else {
            if ( (input$bivariate==TRUE) & (is.numeric(eval(inputVar1()))==TRUE) & (is.numeric(eval(inputVar2()))==TRUE) ) {
              caption2 <- "Numeric Summary First Variable"
            } else {
              if ( (input$bivariate==TRUE) & (is.factor(eval(inputVar1()))==TRUE) & (is.numeric(eval(inputVar2()))==TRUE) ) {
                caption2 <- "Numeric Summary By Factor"
              } else {
                if ( (input$bivariate==TRUE) & (is.numeric(eval(inputVar1()))==TRUE) & (is.factor(eval(inputVar2()))==TRUE) ) {
                  caption2 <- "Numeric Summary By Factor"
                }
              }
            }
          }
        }
      }
      
    })
    
    output$out2 <- renderPrint({
      if ( (input$bivariate==FALSE) & (is.factor(eval(inputVar1()))==TRUE) ) {
        table(eval(inputVar1()))
      } else {
        if ( (input$bivariate==FALSE) & (is.numeric(eval(inputVar1()))==TRUE) ) {
          summary(eval(inputVar1()))
        } else {
          if ( (input$bivariate==TRUE) & (is.factor(eval(inputVar1()))==TRUE) & (is.factor(eval(inputVar2()))==TRUE) ) {
            table(eval(inputVar1()), eval(inputVar2()))
          } else {
            if ( (input$bivariate==TRUE) & (is.numeric(eval(inputVar1()))==TRUE) & (is.numeric(eval(inputVar2()))==TRUE) ) {
              summary(eval(inputVar1()))
            } else {
              if ( (input$bivariate==TRUE) & (is.factor(eval(inputVar1()))==TRUE) & (is.numeric(eval(inputVar2()))==TRUE) ) {
                by(eval(inputVar2()), eval(inputVar1()), summary)
              } else {
                if ( (input$bivariate==TRUE) & (is.numeric(eval(inputVar1()))==TRUE) & (is.factor(eval(inputVar2()))==TRUE) ) {
                  by(eval(inputVar1()), eval(inputVar2()), summary)
                }
              }
            }
          }
        }
      }
      
    })
    
    output$caption3 <- renderText({
      if ( (input$bivariate==FALSE) & (is.factor(eval(inputVar1()))==TRUE) ) {
        caption3 <- "Univariate Table of Proportions"
      } else {
        if ( (input$bivariate==FALSE) & (is.numeric(eval(inputVar1()))==TRUE) ) {
          caption3 <- ""
        } else {
          if ( (input$bivariate==TRUE) & (is.factor(eval(inputVar1()))==TRUE) & (is.factor(eval(inputVar2()))==TRUE) ) {
            captions3 <- "Bivariate Table of Row Proportions"
          } else {
            if ( (input$bivariate==TRUE) & (is.numeric(eval(inputVar1()))==TRUE) & (is.numeric(eval(inputVar2()))==TRUE) ) {
              caption3 <- "Numeric Summary Second Variable"
            } else {
              if ( (input$bivariate==TRUE) & (is.factor(eval(inputVar1()))==TRUE) & (is.numeric(eval(inputVar2()))==TRUE) ) {
                caption3 <- "Kruskal Wallis Test"
              } else {
                if ( (input$bivariate==TRUE) & (is.numeric(eval(inputVar1()))==TRUE) & (is.factor(eval(inputVar2()))==TRUE) ) {
                  caption3 <- "Kruskal Wallis Test"
                }
              }
            }
          }
        }
      }
      
    })
    
    
    output$out3 <- renderPrint({
      
      if ( (input$bivariate==FALSE) & (is.factor(eval(inputVar1()))==TRUE) ) {
        prop.table(table(eval(inputVar1())))
      } else {
        
        if ( (input$bivariate==FALSE) & (is.numeric(eval(inputVar1()))==TRUE) ) {
          NULL
        } else {
          
          if ( (input$bivariate==TRUE) & (is.factor(eval(inputVar1()))==TRUE) & (is.factor(eval(inputVar2()))==TRUE) ) {
            prop.table(table(eval(inputVar1()), eval(inputVar2())), margin=1)
          } else {
            
            if ( (input$bivariate==TRUE) & (is.numeric(eval(inputVar1()))==TRUE) & (is.numeric(eval(inputVar2()))==TRUE) ) {
              summary(eval(inputVar2()))
            } else {
              
              if ( (input$bivariate==TRUE) & (is.factor(eval(inputVar1()))==TRUE) & (is.numeric(eval(inputVar2()))==TRUE) ) {
                kruskal.test(eval(inputVar2()) ~ eval(inputVar1()))
              } else { 
                
                if ( (input$bivariate==TRUE) & (is.numeric(eval(inputVar1()))==TRUE) & (is.factor(eval(inputVar2()))==TRUE) ) {
                  kruskal.test(eval(inputVar1()) ~ eval(inputVar2()))
                }
              }
            }
          }
        }
      }
      
    })
    
    output$caption4 <- renderText({
      
      if ( (input$bivariate==FALSE) & (is.factor(eval(inputVar1()))==TRUE) ) {
        caption4 <- ""
      } else {
        
        if ( (input$bivariate==FALSE) & (is.numeric(eval(inputVar1()))==TRUE) ) {
          caption4 <- ""
        } else {
          
          if ( (input$bivariate==TRUE) & (is.factor(eval(inputVar1()))==TRUE) & (is.factor(eval(inputVar2()))==TRUE) ) {
            captions4 <- "Pearsons Chi-Squared Test"
          } else {
            
            if ( (input$bivariate==TRUE) & (is.numeric(eval(inputVar1()))==TRUE) & (is.numeric(eval(inputVar2()))==TRUE) ) {
              caption4 <- "Spearmans Correlation Coefficient"
            } else {
              
              if ( (input$bivariate==TRUE) & (is.factor(eval(inputVar1()))==TRUE) & (is.numeric(eval(inputVar2()))==TRUE) ) {
                caption4 <- ""
              } else { 
                
                if ( (input$bivariate==TRUE) & (is.numeric(eval(inputVar1()))==TRUE) & (is.factor(eval(inputVar2()))==TRUE) ) {
                  caption4 <- ""
                }
              }
            }
          }
        }
      }
      
    })
    
    output$out4 <- renderPrint({
      
      if ( (input$bivariate==FALSE) & (is.factor(eval(inputVar1()))==TRUE) ) {
        NULL
      } else {
        
        if ( (input$bivariate==FALSE) & (is.numeric(eval(inputVar1()))==TRUE) ) {
          NULL
        } else {
          
          if ( (input$bivariate==TRUE) & (is.factor(eval(inputVar1()))==TRUE) & (is.factor(eval(inputVar2()))==TRUE) ) {
            chisq.test(table(eval(inputVar1()), eval(inputVar2())))
          } else {
            
            if ( (input$bivariate==TRUE) & (is.numeric(eval(inputVar1()))==TRUE) & (is.numeric(eval(inputVar2()))==TRUE) ) {
              cor(eval(inputVar1()), eval(inputVar2()), method="spearman", use="pairwise.complete.obs")
            } else {
              
              if ( (input$bivariate==TRUE) & (is.factor(eval(inputVar1()))==TRUE) & (is.numeric(eval(inputVar2()))==TRUE) ) {
                NULL
              } else { 
                
                if ( (input$bivariate==TRUE) & (is.numeric(eval(inputVar1()))==TRUE) & (is.factor(eval(inputVar2()))==TRUE) ) {
                  NULL
                }
              }
            }
          }
        }
      }
      
    })
    
    output$caption5 <- renderText({
      
      if ( (input$bivariate==FALSE) & (is.factor(eval(inputVar1()))==TRUE) ) {
        caption5 <- "Univariate Barplot"
      } else {
        
        if ( (input$bivariate==FALSE) & (is.numeric(eval(inputVar1()))==TRUE) ) {
          caption5 <- "Univariate Histogram"
        } else {
          
          if ( (input$bivariate==TRUE) & (is.factor(eval(inputVar1()))==TRUE) & (is.factor(eval(inputVar2()))==TRUE) ) {
            captions5 <- "Bivariate Barplot"
          } else {
            
            if ( (input$bivariate==TRUE) & (is.numeric(eval(inputVar1()))==TRUE) & (is.numeric(eval(inputVar2()))==TRUE) ) {
              caption5 <- "Bivariate Scatter Graph"
            } else {
              
              if ( (input$bivariate==TRUE) & (is.factor(eval(inputVar1()))==TRUE) & (is.numeric(eval(inputVar2()))==TRUE) ) {
                caption5 <- "Bivariate Boxplot"
              } else { 
                
                if ( (input$bivariate==TRUE) & (is.numeric(eval(inputVar1()))==TRUE) & (is.factor(eval(inputVar2()))==TRUE) ) {
                  caption5 <- "Bivariate Boxplot"
                }
              }
            }
          }
        }
      }
      
    })
    
    output$out5 <- renderPlot({
      
      if ( (input$bivariate==FALSE) & (is.factor(eval(inputVar1()))==TRUE) ) {
        barplot(table(eval(inputVar1())))
      } else {
        
        if ( (input$bivariate==FALSE) & (is.numeric(eval(inputVar1()))==TRUE) ) {
          hist(eval(inputVar1()),main="")
        } else {
          
          if ( (input$bivariate==TRUE) & (is.factor(eval(inputVar1()))==TRUE) & (is.factor(eval(inputVar2()))==TRUE) ) {
            barplot(table(eval(inputVar1()), eval(inputVar2())), beside=TRUE)
          } else {
            
            if ( (input$bivariate==TRUE) & (is.numeric(eval(inputVar1()))==TRUE) & (is.numeric(eval(inputVar2()))==TRUE) ) {
              plot(eval(inputVar1()), eval(inputVar2()), main="")
            } else {
              
              if ( (input$bivariate==TRUE) & (is.factor(eval(inputVar1()))==TRUE) & (is.numeric(eval(inputVar2()))==TRUE) ) {
                boxplot(eval(inputVar2()) ~ eval(inputVar1()))
              } else { 
                
                if ( (input$bivariate==TRUE) & (is.numeric(eval(inputVar1()))==TRUE) & (is.factor(eval(inputVar2()))==TRUE) ) {
                  boxplot(eval(inputVar1()) ~ eval(inputVar2()))
                }
              }
            }
          }
        }
      }
      
    })
    
  })
))
