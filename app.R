
library(ggplot2)
library(plotly)
library(shiny)
library(dplyr)
library(DT)
library(invgamma)

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
                    strong(h4("Intervalo de confianza")),
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
                      
                      selectInput("vi","Escoge tu variable independiente (X):", names(vino), selected = 'pH'),
                      selectInput("vd","Escoge tu variable dependiente (Y):", names(vino), selected = 'fixed.acidity'),
                      sliderInput("num_cadenas", label = 'Número de cadenas', value = 1, min = 1, max = 10, step = 1),
                      sliderInput('long_cadenas', label = 'Longitud de las cadenas', min = 500, max = 100000, value =1000, step = 500)
                      
                      
                      
                    ), 
                    mainPanel(
                      
                      plotOutput("scatterpl"),
                      h5("Distribuciones a priori:"),
                      plotOutput("aprioris"),
                      plotOutput('histos'),
                      p('Resultados simulacion'),
                      dataTableOutput('sim'),
                      textOutput('test')
                      
                      
                      
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
 # output$ks <- renderPrint({
    
  #  if (pval["p.value"]<0.05){
   #   paste( "Con p-value",pval, "concluimos es distinta a una distribución exponencial.")
  #  }else{
   #   paste( "Con p-value",pval, "concluimos que no es distinta a una distribución exponencial.")
  #  }
#  })
  
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
    
    ggplot(vino,aes(eval(parse(text=input$vi)),eval(parse(text=input$vd))))+geom_point()+xlab(input$vi)+ylab(input$vd)
    
  })
  df <- reactive({
    f <- vino
    f
  })
  
  data <- reactive({
    f <-vino
    as.matrix(f)
  })
  
  output$data <- renderTable({
    data()
  })
  
  output$vdui <- renderUI({
    selectInput('vd', label = 'Variable dependiente', names(df()))
  })
  
  output$viui <- renderUI({
    selectInput('vi', label = 'Variable independiente', names(df()))
  })
  
  output$varpt <- renderPlot({
    if(!is.null(df))
    {
      y <- df()[[input$vd]]
      x <- df()[[input$vi]]
      plot(x,y)
    }
  })
  
  nreg <- reactive({
    equis<-df()[[input$vi]]
    ye <- df()[[input$vd]]
    sp<-data.frame(equis,ye)
    splm<-lm(ye~equis,data=sp)
    summary_splm<-summary(splm)
    betas<-coefficients(summary_splm)
    list('betas' = betas, 'summary' = summary_splm)
  })
  
  ndist <- reactive({
    x <- seq(-100, 100, length=100)
    dnorm(x,round(nreg()$betas[1,1],digits=2),round(nreg()$betas[1,2],digits=2))
  })
  
  gdist <- reactive({
    x <- seq(-100, 100, length=100)
    dinvgamma(x,13.5,round(25*nreg()$summary$sigma,digits=2))
  })
  

  
 
  
  
  output$aprioris<-renderPlot({
    x <- seq(-100, 100, length=100)
    par(mfrow = c(3,1))
    plot(x, gdist(), type="l", lty=2, xlab="x",ylab="Density", main=paste('A priori eps'))
    plot(x, ndist(), type="l", lty=2, xlab="x",ylab="Density", main=paste('A priori b'))
    plot(x, ndist(), type="l", lty=2, xlab="x",ylab="Density", main=paste('A priori a'))
    })
  
  
  
  likelihood <- function(param){
    b1= param[1]
    b0 = param[2]
    sigma2 = param[3]
    pred = b1*df()[[input$vi]] + b0
    singlelikelihoods = dnorm(df()[[input$vd]], mean = pred, sd = sigma2**.5, log = T)
    sumll = sum(singlelikelihoods)
    return(sumll)
  }
  
  prior <- function(param){
    b1 = param[1]
    b0 = param[2]
    sigma2 = param[3]
    b1prior = dnorm(b1, mean=round(nreg()$betas[1,1],digits=2), sd=round(nreg()$betas[1,2]**.5,digits=2), log = T)
    b0prior = dnorm(b0, mean=round(nreg()$betas[2,1],digits=2), sd=round(nreg()$betas[2,2]**.5,digits=2), log = T)
    sigma2prior = dinvgamma(sigma2,14,round(25*nreg()$summary$sigma,digits=2),log = T)
    return(b1prior+b0prior+sigma2prior)
  }
  
  posterior <- function(param){
    return (likelihood(param) + prior(param))
  }
  
  #Metropolis
  
  proposalfunction <- function(param){
    return(rnorm(3,mean = param, sd= c(0.1,0.5,0.3)))
  }
  
  run_metropolis_MCMC <- function(startvalue, iterations){
    chain <- array(dim = c(iterations+1,3))
    chain[1,] <- startvalue
    for (i in 1:iterations){
      proposal <- proposalfunction(chain[i,])
      
      logprobab =posterior(proposal) - posterior(chain[i,])
      if (log(runif(1)) <= logprobab){
        chain[i+1,] = proposal
      }else{
        chain[i+1,] = chain[i,]
      }
    }
    return(chain)
  }
  
  mcmc <- reactive({
    startvalue = c(rnorm(1,0,1),rnorm(1,0,1),rinvgamma(1,1,1))
    chain = run_metropolis_MCMC(startvalue, input$long_cadenas)
    data.frame(b1=chain[,1],b0=chain[,2],s2=chain[,3])
  })
  
  output$sim <- renderDataTable({
    mcmc()
  })
  
  output$histos<-renderPlot({
    burnIn = input$long_cadenas*.20
    acceptance = 1-mean(duplicated(mcmc()[-(1:burnIn),]))
    par(mfrow = c(2,3))
    hist(mcmc()[-(1:burnIn),1],nclass=30,  main="Posterior of b1", xlab="Parametro" )
    abline(v = mean(mcmc()[-(1:burnIn),1]))
    hist(mcmc()[-(1:burnIn),2],nclass=30, main="Posterior of b0", xlab="Parametro")
    abline(v = mean(mcmc()[-(1:burnIn),2]))
    hist(mcmc()[-(1:burnIn),3],nclass=30, main="Posterior of sigma^2", xlab="Parametro")
    abline(v = mean(mcmc()[-(1:burnIn),3]) )
    plot(mcmc()[-(1:burnIn),1], type = "l", xlab="Iteraciones" , main = "Chain values of b1" )
    plot(mcmc()[-(1:burnIn),2], type = "l", xlab="Iteraciones" , main = "Chain values of b0")
    plot(mcmc()[-(1:burnIn),3], type = "l", xlab="Iteraciones" , main = "Chain values of sigma^2")
  })
  
  
  
  
  
  
  ###
  
  
  
  
  
  
  
  
  
  
  
  
}





shinyApp(ui = ui, server = server)
