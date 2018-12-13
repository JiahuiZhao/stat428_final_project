library(matrixcalc)
library(shiny)
library(coda)
library(stringr)
library(ggplot2)

ui = fluidPage(
  titlePanel("Variance Reduction for Estimating VAR with MCMC"),
  
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
        tabPanel("Risk Factor", plotOutput("plot")),
        tabPanel("Risk Factor Scatterplot", plotOutput("scatterplot")),
        tabPanel("Probability & Variance Factor Table", tableOutput("table"), h4("It may take a while to show up the table."))
        )
    )))

server = function(input, output) {
  
  active_variance = eventReactive(input$run, {
    
    A = matrix(as.numeric(str_split(input$capitalA, pattern = ",")[[1]]), 2, 2, byrow = T)
    a = matrix(as.numeric(str_split(input$capitala, pattern = ",")[[1]]), 2)
    a0 = as.numeric(input$capitala0)

    set.seed(1234)
    
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
    
    set.seed(1234)

    Delta_S <- Multi_Normal(matrix(c(0,0),ncol=1), Sigma)

    x_var <- VAR(a0, a, A, Sigma, Delta_S, percentile)

   variance_general <- CV(a0, a, A, Sigma, x_var)[3:4]
   variance_general1 <-  c(round(variance_general[1], 4), round(variance_general[2],4),round(variance_general[2]/variance_general[2],4))
   
   variance_is <- IS(a0, a, A, Sigma, x_var)
   variance_is <-  c(round(variance_is[1],4), round(variance_is[2],4),round(variance_general[2]/variance_is[2],4))
   
   variance_cv <- CV(a0, a, A, Sigma, x_var)[1:2]
   variance_cv <-  c(round(variance_cv[1],4), round(variance_cv[2],4),round(variance_general[2]/variance_cv[2],4))
   
   variance_ss <- SS(a0, a, A, Sigma, x_var)
   variance_ss <-  c(round(variance_ss[1],4), round(variance_ss[2],4),round(variance_general[2]/variance_ss[2],4))
   
   variance <-  rbind(variance_general1, variance_cv, variance_is, variance_ss)

   variance.table <-  cbind(c("General" ,"Control Variates", "Importance Sampling","Stratified Sampling"), variance)
   
   colnames(variance.table) <-  c("Variance Deduction Method","Probability", "Variance","Variance Deduction Factor")
   
   return(variance.table)
  })
  
  active_plot = eventReactive(input$run, {
    A = matrix(as.numeric(str_split(input$capitalA, pattern = ",")[[1]]), 2, 2, byrow = T)
    a = matrix(as.numeric(str_split(input$capitala, pattern = ",")[[1]]), 2)
    a0 = as.numeric(input$capitala0)
    
    set.seed(12)
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
    
    set.seed(12)
    Sigma <- cov(data.frame(delta_Bond,delta_Interest))
    Sigma <- as.matrix(Sigma)
    Delta_S <- Multi_Normal(matrix(c(0,0),ncol=1), Sigma)
    gelman.plot(mcmc.list(as.mcmc(Delta_S[,1]), as.mcmc(Delta_S[,2])))
  })
  
  active_scatterplot = eventReactive(input$run, {
    A = matrix(as.numeric(str_split(input$capitalA, pattern = ",")[[1]]), 2, 2, byrow = T)
    a = matrix(as.numeric(str_split(input$capitala, pattern = ",")[[1]]), 2)
    a0 = as.numeric(input$capitala0)
    
    set.seed(1234)
    Interest <- read.csv("Interest.csv", header = TRUE)
    Bond <- read.csv("Bond.csv", header = TRUE)
    Bond <- Bond[c(625:688),2]
    Interest <- Interest[c(745:808),2]
    delta_Bond <- numeric(length(Bond))
    delta_Interest <- numeric(length(Interest))
    for(i in 2 : length(Bond)){
      delta_Bond[i-1] <- Bond[i] - Bond[i-1]
      delta_Interest[i-1] <- Interest[i] - Interest[i-1]
    }
    
    set.seed(1234)
    Sigma <- cov(data.frame(delta_Bond,delta_Interest))
    Sigma <- as.matrix(Sigma)
    Delta_S <- Multi_Normal(matrix(c(0,0),ncol=1), Sigma)
    x = data.frame(Delta_S)
    ggplot(x, aes(x = x[, 1], y= x[, 2]), xlab = "a") + geom_point() + xlab("First Risk Factor") + ylab("Second Risk Factor") + ggtitle("Risk Factors") +
      theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
  })
  
  output$plot = renderPlot({active_plot()})
  output$scatterplot = renderPlot({active_scatterplot()})
  output$table = renderTable({active_variance()})
  }
  

shinyApp(ui = ui, server = server)
