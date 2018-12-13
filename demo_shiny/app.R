library(matrixcalc)
library(shiny)
library(coda)
library(stringr)

ui = fluidPage(
  titlePanel("MCMC Financial Variance Reduction"),
  
  sidebarLayout(
    sidebarPanel(
      br(),
      textInput(inputId = "capitalA", label = "Type in your A, seperate with ','", value = "1,0,1,0"),
      br(),
      textInput(inputId = "capitala", label = "Type in your a, seperate with ','", value = "1,0"),
      br(),
      textInput(inputId = "capitala0", label = "Type in your a0", value = "0"),
      br(),
      textInput(inputId = "pencentile", label = "Typy in your Pencentile", value = "0"),
      br(),
      actionButton(inputId = "run", label = "Run")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Delta_S", plotOutput("plot")),
        tabPanel("Probability & Variance Factor Table", tableOutput("table"), h4("It may take a while to show up the table."))
      )
    )))

server = function(input, output) {
  
  active_variance = eventReactive(input$run, {
    
    A = matrix(as.numeric(str_split(input$capitalA, pattern = ",")[[1]]), 2, 2, byrow = T)
    a = matrix(as.numeric(str_split(input$capitala, pattern = ",")[[1]]), 2)
    a0 = as.numeric(input$capitala0)

    Interest <- read.csv("Interest.csv", header = TRUE)
    Bond <- read.csv("Bond.csv", header = TRUE)
    Bond <- Bond[c(625:688),2]
    Interest <- Interest[c(745:808),2]
    delta_Bond <- numeric(length(Bond))
    delta_Interest <- numeric(length(Interest))
    for(i in 2 : length(Bond) ){
      delta_Bond[i-1] <- Bond[i] - Bond[i-1]
      delta_Interest[i-1] <- Interest[i] - Interest[i-1]
    }
   
    Sigma <- cov(data.frame(delta_Bond,delta_Interest))
    Sigma <- as.matrix(Sigma)

    Delta_S <- Multi_Normal(matrix(c(0,0),ncol=1), Sigma)

    x_var <- VAR(a0, a, A, Sigma, Delta_S, percentile)

   variance_general <- General(a0, a, A, Sigma, x_var)
   variance_general <-  c(variance_general[1], variance_general[2]/variance_general[2])
   
   variance_is <- IS(a0, a, A, Sigma, x_var)
   variance_is <-  c(variance_is[1], variance_general[2]/variance_is[2])
   
   variance_cv <- CV(a0, a, A, Sigma, x_var)
   variance_cv <-  c(variance_cv[1], variance_general[2]/variance_cv[2])
   
   variance_ss <- SS(a0, a, A, Sigma, x_var)
   variance_ss <-  c(variance_ss[1], variance_general[2]/variance_ss[2])
   
   variance <-  rbind(variance_general, variance_is, variance_cv, variance_ss)
   colnames(variance) <-  c("Probability", "Variance")
   variance.table <-  cbind(c("General", "Importance Sampling","Control Variates","Stratified Sampling"), variance)
   
   return(variance.table)
  })
  
  active_plot = eventReactive(input$run, {
    A = matrix(as.numeric(str_split(input$capitalA, pattern = ",")[[1]]), 2, 2, byrow = T)
    a = matrix(as.numeric(str_split(input$capitala, pattern = ",")[[1]]), 2)
    a0 = as.numeric(input$capitala0)
    
    Interest <- read.csv("Interest.csv", header = TRUE)
    Bond <- read.csv("Bond.csv", header = TRUE)
    Bond <- Bond[c(625:688),2]
    Interest <- Interest[c(745:808),2]
    delta_Bond <- numeric(length(Bond))
    delta_Interest <- numeric(length(Interest))
    for(i in 2 : length(Bond) ){
      delta_Bond[i-1] <- Bond[i] - Bond[i-1]
      delta_Interest[i-1] <- Interest[i] - Interest[i-1]
    }
    
    Sigma <- cov(data.frame(delta_Bond,delta_Interest))
    Sigma <- as.matrix(Sigma)
    a0 <- 0
    a <- matrix(c(1,1), ncol = 1)
    A <- matrix(c(1,0,0,1), ncol = 2, nrow = 2, byrow = TRUE)
    percentile <- 0.9
    Delta_S <- Multi_Normal(matrix(c(0,0),ncol=1), Sigma)
    gelman.plot(mcmc.list(as.mcmc(Delta_S[,1]), as.mcmc(Delta_S[,2])))
  })
  output$plot = renderPlot({active_plot()})
  output$table = renderTable({active_variance()}, digits = 10)
  }
  

shinyApp(ui = ui, server = server)
