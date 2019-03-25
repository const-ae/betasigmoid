library(shiny)

ui <- fluidPage(
  titlePanel("BetaSigmoid Distribution Density Plot"),
  sidebarLayout(
    sidebarPanel=  sidebarPanel(
      sliderInput(inputId = "a", label = "a", min= 0, max=1e3, value = 1),
      sliderInput(inputId = "fine_a", label = "fine_a", min= -100, max=0, value = -100),
      sliderInput(inputId = "b", label = "b", min= 0, max=1e3, value = 1),
      sliderInput(inputId = "fine_b", label = "fine_b", min= -100, max=0, value = -100),
      sliderInput(inputId = "infl", label = "infl", min= - 20, max = 20, value = 0),
      sliderInput(inputId = "scale", label = "scale", min= - 20, max = 20, value = 1),
      sliderInput(inputId = "fine_scale", label = "fine_scale", min= -100, max = 0, value = -100)
    ),
     mainPanel = mainPanel(
      plotOutput("distrPlot"),
      plotOutput("logDistrPlot")
     )
  )
)

server <- function(input, output){

  output$distrPlot <- renderPlot({
    a <- input$a
    b <- input$b
    infl <- input$infl
    scale <- input$scale
    if(a == 0){
      a <- exp(input$fine_a)
    }
    if(b == 0){
      b <- exp(input$fine_b)
    }
    if(scale == 0){
      scale <- exp(input$fine_scale)
    }

    dist_mean <- mean_betasigmoid(a, b, infl, scale)
    dist_sd <- sqrt(var_betasigmoid(a, b, infl, scale))

    xg <- seq(dist_mean - 8 * dist_sd, dist_mean + 8 * dist_sd, length.out = 1001)
    plot(xg, dbetasigmoid(xg, a, b, infl, scale), type="l")
  })
  output$logDistrPlot <- renderPlot({
    a <- input$a
    b <- input$b
    infl <- input$infl
    scale <- input$scale
    if(a == 0){
      a <- exp(input$fine_a)
    }
    if(b == 0){
      b <- exp(input$fine_b)
    }
    if(scale == 0){
      scale <- exp(input$fine_scale)
    }

    dist_mean <- mean_betasigmoid(a, b, infl, scale)
    dist_sd <- sqrt(var_betasigmoid(a, b, infl, scale))

    xg <- seq(dist_mean - 8 * dist_sd, dist_mean + 8 * dist_sd, length.out = 1001)
    plot(xg, dbetasigmoid(xg, a, b, infl, scale, log=TRUE), type="l")
  })
}

shinyApp(ui, server)
