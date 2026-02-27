#-------------------------------------------------
# Libraries
#-------------------------------------------------
library(shiny)
library(bslib)
library(ggplot2)
library(thematic)
library(shinycssloaders)
#-------------------------------------------------
#Enable thematic for automatic plot stylin
thematic_shiny(font = "auto")
#-------------------------------------------------
# MATH FUNCTIONS - Unit Sine Kumaraswamy
#-------------------------------------------------
#Probability Distribution Function----
dSKW<-function(y,delta,sigma,log=FALSE)
{
  if (any(delta <= 0)) stop(paste("delta must be positive", "\n", ""))
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", "")) 
  if (any(y <= 0) | any(y >= 1)) stop(paste("y must be between 0 and 1", "\n", ""))
  fy1<-(pi*delta*sigma*y^(delta-1)*(1-y^delta)^sigma*(2*(1-y^delta)^sigma-3)*cos((pi*(1-(1-y^delta)^sigma)*(2-(1-y^delta)^sigma))/4))/(4*(y^delta-1))
  if(log==FALSE) fy<-fy1 else fy<-log(fy1)
  return(fy)
}
#-------------------------------------------------
#Cumulative Distribution Function----
pSKW<-function(q, delta, sigma, lower.tail = TRUE, log.p = FALSE){
  if (any(delta <= 0)) stop(paste("delta must be positive", "\n", ""))
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", ""))
  if (any(q <= 0) | any(q >= 1)) stop(paste("q must be between 0 and 1", "\n", ""))
  cdf1<-sin((pi*(1-(1-q^delta)^sigma)*(2-(1-q^delta)^sigma))/4)
  if(lower.tail==TRUE) cdf<-cdf1 else cdf<- 1-cdf1
  if(log.p==FALSE) cdf<- cdf else cdf<- log(cdf)
  return(cdf)
}
#-------------------------------------------------
#Hazard Function----
hSKW <- function(y, delta, sigma) {
  f <- dSKW(y, delta, sigma)
  F <- pSKW(y, delta, sigma)
  f / (1 - F)
}
#-------------------------------------------------
#Survival Function----
sSKW <- function(y, delta, sigma) {
  1 - pSKW(y, delta, sigma)
}
#-------------------------------------------------
# UI
#-------------------------------------------------
ui <- fluidPage(
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  withMathJax(),
  titlePanel("Unit Sine Kumaraswamy"),
  hr(),
  tabsetPanel(
    tabPanel(title = "SKW Distribution",
             sidebarLayout(
               sidebarPanel(
                 h4("Parameters"),
                 sliderInput(inputId = "delta",
                             label=(helpText(c("$$\\delta$$"))),
                             min = 0.01, max = 20, value = 1.8, step = .2,
                             animate = animationOptions(interval = 2000, loop = T)),
                 sliderInput(inputId = "sigma",
                             label=(helpText(c("$$\\sigma$$"))),
                             min = 0.01, max = 20, value = 0.8, step = .2,
                             animate = animationOptions(interval = 2000, loop = T))),
               mainPanel(
                 fluidRow(
                   column(6, plotOutput("p1")), column(6, plotOutput("p2")),
                   column(6, plotOutput("p3")), column(6, plotOutput("p4"))
                 )
               )
             )
    )
  )
)
#-------------------------------------------------
# SERVER
#-------------------------------------------------
server <- function(input, output) {
  my_plot <- function(df, x, y, title, color) {
    ggplot(df, aes(x = .data[[x]], y = .data[[y]])) +
      geom_line(color = color, linewidth = 1) +
      labs(title = title) +
      theme_minimal()
  }
  sub_df <- reactive({
    y <- seq(0.001, 0.999, length.out = 300)
    dens <- dSKW(y, input$delta, input$sigma)
    cdf <- pSKW(y, input$delta, input$sigma)
    data.frame(y = y, dens = dens, cdf = cdf, surv = 1 - cdf, haz = dens/(1-cdf))
  })
  output$p1 <- renderPlot({ my_plot(sub_df(), "y", "dens", "PDF", "#3498db") })
  output$p2 <- renderPlot({ my_plot(sub_df(), "y", "cdf", "CDF", "#e74c3c") })
  output$p3 <- renderPlot({ my_plot(sub_df(), "y", "surv", "Survival", "#27ae60") })
  output$p4 <- renderPlot({ my_plot(sub_df(), "y", "haz", "Hazard", "#2c3e50") })
}
#-------------------------------------------------
shinyApp(ui, server)
#-------------------------------------------------