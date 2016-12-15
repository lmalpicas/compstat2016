
library(ggplot2)
library(plotly)
library(shiny)
library(dplyr)

ui <- fluidPage(
  theme = "style.css",
  sidebarPanel(
    h2("Estadística computacional"),
    h3("Otoño 2016"),
    h5("Lorena Malpica Serrano 124841"),
    h3("Índice de tareas:"),
    h5("Tarea1: Método de la función inversa"),
    h5("Tarea2: Integración por Monte Carlo"),
    h5("Tareas 4-6: Regresión lineal con MCMC")
    ), # sidebar panel main 
  mainPanel(
    ###funcion inversa
    tabsetPanel(id="tab", 
                tabPanel("Tarea1",
                         sidebarLayout(
                           sidebarPanel ( 
                             
                             numericInput(inputId = "nsim", label = "Número de simulaciones", min= 10, max= 1000000000000, value = 1000),
                             numericInput(inputId = "lambda", label = "Parámetro lambda", value=10)
                             
                           ), #sidebar panel
                           mainPanel(
                             h3("Simulación x ~ Exp(lambda) por el método de la función inversa."),
                             plotlyOutput(outputId = "plotinv"),
                             sliderInput("nbins", "Número de bins", value=20, min=10, max=100),
                             h4("Prueba Kolmogorov-Smirnoff"),
                             verbatimTextOutput("ks")
                           )#main panel  
                           
                         )#sidebar layout
                         
                         
                         
                ),#tabpanel tarea1
                
                tabPanel(
                  "Tarea 2",
                  
                  sidebarPanel(
                    
                    textInput(inputId = "funcion", label ="Introduce la función a integrar", value = "x^2"),
                    numericInput("lim_a", "Límite inferior", 0),
                    numericInput("lim_b", "Límite superior", 3.141593),
                    strong(h4("Nivel de significancia")),
                    numericInput("nalpha", "Ingresa alpha", value = 0.05,
                                 min = .01, max = .90)
                    
                    
                  ),
                  mainPanel(
                    
                    h3("Integración por método Monte Carlo"),
                    textOutput("rint"),
                    plotOutput("plot1"),
                    h4("Intervalos de confianza para diferentes números de simulaciones"),
                    plotOutput("plot2")
                  )
                ), #tabPanel tarea2
                tabPanel
                ("Tareas 4-6", 
                  h3("Regresión lineal con MCMC"),
                  (fluidRow(column(12, dataTableOutput("table")))),
                  sidebarLayout(
                    sidebarPanel ( 
                      
                      selectInput("x","Escoge tu variable independiente (X):", names(vino), selected = 'pH'),
                      selectInput("y","Escoge tu variable dependiente (Y):", names(vino), selected = 'fixed.acidity'),
                      sliderInput("num_cadenas", label = 'Número de cadenas', value = 1, min = 1, max = 10, step = 1),
                      sliderInput('long_cadenas', label = 'Longitud de las cadenas', min = 500, max = 100000, value =1000, step = 500)
                      
                      
                      
                    ), 
                    mainPanel(
                      
                      plotOutput("scatterpl")
                      
                      
                      
                    )  
                    
                  )
                ), #tab panel tarea 4-6,
                tabPanel("Datos completos tareas 4-6", 
                         (fluidRow(column(12, dataTableOutput("table2"))))
                ) #tab panel solo datos
                  
                
                
    )#tabset panel
  
  )#main panel
  
) #Fluid page 

server <- function(input, output) {
  
  output$plotinv <- renderPlotly({
    set.seed(110104)
    nsim <- input$nsim
    lambda <- input$lambda
    nbins<-input$nbins
    rExp <- function(nsim, lambda){
      return((-1/lambda)*log(1-runif(nsim)))
    }
    dat <- data.frame(Value=rExp(nsim, lambda))
    
    pl<-ggplot(dat, aes(x=Value)) +
      geom_histogram(aes(y=..density..), bins= nbins, fill="skyblue") +
      stat_function(fun = function(x) lambda*exp(-lambda*x),colour = "blue")
    ggplotly(pl)
    
  })
  output$ks <- renderPrint({
    pval <- ks.test(x = Value, 
                    y = "pexp", 
                    rate = input$lambda)
    
    if (pval["p.value"]<0.05){
      paste( "Con p-value",pval, "concluimos es distinta a una distribución exponencial.")
    }else{
      paste( "Con p-value",pval, "concluimos que no es distinta a una distribución exponencial.")
    }
  })
  
#Tarea2   
  
  int_fun <- reactive({
    texto <- paste("aux <- function(x) (",
                   input$lim_b, "-", input$lim_a,")*", 
                   input$funcion)
    eval(parse(text = texto))
    aux
  })
  
  simul <- reactive({
    as.numeric(100000)
  })
  
  int_tbl <- reactive({
    fint <- int_fun()
    SimIntNsims <- function(n.sims){
      set.seed(160911)
      dat.sims <- data.frame(
        us = runif(n.sims, min = input$lim_a, max = input$lim_b)
      )%>% 
        mutate(fus = fint(us))
      estims <- dat.sims %>% 
        summarise(simid = n.sims,
                  muest = abs(mean(fus)), 
                  sdest = sd(fus),
                  centdec = qnorm(input$nalpha/2, lower.tail = F)*sdest/sqrt(n.sims),
                  lowint = muest - centdec,
                  uppint = muest + centdec)
      estims
    }
    step.id <- ifelse(simul() > 1000, 150, 10)
    sims.v <- c(10, 100, 1000, 10000, 100000)
    
    # seq.sims <- seq(10, simul(), by = step.id)
    seq.sims <- sims.v[1:which(sims.v == simul() )]
    
    dat.ints <- lapply( seq.sims, SimIntNsims) %>% 
      rbind_all()
    dat.ints
  })
  
  output$rint<-renderText(
    {
      num.mu <- filter(int_tbl(), simid == max(simid)) %>% .$muest %>% unique
      paste("Resultado de la integral: ", round(num.mu, 2))  
      
    }
  )
  
  
  output$plot1 <- renderPlot({
    num.mu <- filter(int_tbl(), simid == max(simid)) %>% .$muest %>% unique
    lab.mu <- paste("Area estimada: ", round(num.mu, 2)  )
    gg <- ggplot(data.frame(x = c(input$lim_a, input$lim_b)), aes(x)) +
      stat_function(fun = int_fun(), geom = "line", 
                    size = 1.5, color = 'black') +
      stat_function(fun = int_fun(), geom = "ribbon",
                    mapping = aes(ymin = 0, ymax = ..y..),
                    fill = "skyblue", alpha = .3) + 
      
      ylab('f(x)')+ 
      xlab("x") 
    
    print(gg)
  })
  
  
  output$plot2 <- renderPlot({
    gg <- ggplot(int_tbl(), aes(x = simid, y = muest)) + 
      geom_ribbon(aes(ymin = lowint, ymax = uppint), 
                  alpha =.5, fill = "grey") + 
      geom_line(color = 'red') + 
      scale_x_log10() +
      ylab('Valor estimado de la integral') + 
      xlab("Número de simulaciones") 
    print(gg)
    
    
  })
  ##
  
 ##Tareas 4-6 
  output$table <- renderDataTable(head(vino))
  
  
  
  output$table2 <-   renderDataTable(vino)
  
  output$scatterpl<-renderPlot({
    
    ggplot(vino,aes(eval(parse(text=input$x)),eval(parse(text=input$y))))+geom_point()+xlab(input$x)+ylab(input$y)
    
  })
  
  
  ###
  
  
  
  
  
  
  
  
  
  
  
  
}





shinyApp(ui = ui, server = server)