library(ggplot2)
library(dplyr)
library(tidyr)
library(optimx)

library(shiny)
library(shinyjs)

theme_set(theme_bw() + theme(text = element_text(size=20)))

ui <- fluidPage(
  
  useShinyjs(),
  withMathJax(),
  
  title = "Threshold Estimation",
  
  h1("Threshold Estimation using the Selker et al. (2019) Method"),
  br(),
  
  tabsetPanel(
    type = "tabs",
    tabPanel(
      p("Introduction", style="font-size:18px"),
      HTML(
        "<div style=\"font-size:18px\">
    <br></br>
    <a href=\"https://doi.org/10.3758/s13428-019-01231-3\">Selker et al. (2019)</a> present a method for describing the locations of an arbitrary number of thresholds on a normal distribution with just two parameters, \\(\\alpha\\) and \\(\\beta\\):
    <ul>
    <li>\\(\\alpha\\) describes the scale of the threshold locations, and</li>
    <li>\\(\\beta\\) describes the degree to which the thresholds are more leftward or rightward (such that \\(\\beta\\) of 0 describes a symmetrical distribution of thresholds)</li>
    </ul>
    <br></br>
    Selker et al. outline a method where threshold \\(\\lambda\\)\\(_c\\) in position \\(c\\), when there are \\(C\\) ordered regions in the distributions, has a location of:
    $$ \\gamma_c = log \\left( \\frac{c/C}{1-c/C} \\right) $$
    $$ \\lambda_c = \\alpha \\gamma_c + \\beta $$
    <br></br>
    In this app, you can generate thresholds using manual locations or Selker et al.'s method, and then see how well the function can describe the locations. This is done by optimising \\(\\alpha\\) and \\(\\beta\\) to minimise MSE.
    <br></br>
    Click the <b>Explore</b> tab above to have a go!
    </div>"
      )
    ),
    
    tabPanel(
      p("Explore", style="font-size:18px"),
      sidebarPanel(
        h4("Thresholds Setup"),
        selectInput("thresh_set_method", "How should thresholds be specified?", choices = c("Manually", "Using Alpha & Beta"), width="100%"),
        numericInput("n_thresh", "Number of Thresholds", value = 4, min = 2, step = 1, width = "100%"),
        div(
          uiOutput("thresh_locs_manual_ui"),
          hidden(div(p("Warning: thresholds are not in order so will be re-sorted", style="color:red"), id="manual_order_warning_div")),
          id="thresh_locs_manual_div"
        ),
        hidden(
          div(
            sliderInput("alpha", "Alpha", min = -5, max = 5, step = 0.01, value = 1),
            sliderInput("beta", "Beta", min = -5, max = 5, step = 0.01, value = 0),
            id="thresh_locs_ab_div"
          )
        ),
        hr(),
        h4("Estimation Setup"),
        selectInput("optim_method", "Optimisation Algorithm", choices = c("Nelder-Mead", "bobyqa"), width="100%")
      ),
      
      mainPanel(
        fluidRow(
          column(12, h2("Actual vs. Estimated Threshold Locations")),
          column(6, plotOutput("latent_plot")),
          column(6, plotOutput("distort_plot"))
        )
      )
    )
    
  )
  
)

server <- function(input, output) {
  
  # show the relevant UI
  observe({
    if (input$thresh_set_method == "Manually") {
      hide("thresh_locs_ab_div")
      show("thresh_locs_manual_div")
    } else {
      hide("thresh_locs_manual_div")
      show("thresh_locs_ab_div")
    }
  })
  
  # sliders for setting the threshold locations
  output$thresh_locs_manual_ui <- renderUI({
    
    req(input$n_thresh)
    
    default_locs <- sort(runif(input$n_thresh, -2.5, 2.5))
    
    lapply(1:input$n_thresh, function(thresh_i) {
      sliderInput(
        inputId = paste("thresh", thresh_i, sep="_"),
        label = sprintf("Threshold %g Location", thresh_i),
        min = -5,
        max = 5,
        step = 0.01,
        value = default_locs[thresh_i],
        width = "100%"
      )
    })
  })
  
  # collect values for manual threshold locations
  manual_thresh_locs <- reactive({
    lapply(1:input$n_thresh, function(thresh_i) {
      input[[paste("thresh", thresh_i, sep="_")]]
    }) %>%
      unlist()
  })
  
  # show a warning if the manual locations aren't in order
  observe({
    req(input$n_thresh, input$thresh_1, manual_thresh_locs())
    
    if (is.unsorted(manual_thresh_locs())) {
      show("manual_order_warning_div")
    } else {
      hide("manual_order_warning_div")
    }
  })
  
  # estimate alpha and beta
  opt_res <- reactive({
    req(input$thresh_1, input$n_thresh, input$thresh_set_method)
    
    # get desired threshold locations
    actual_locs <- if (input$thresh_set_method=="Manually") {
      manual_thresh_locs() %>%
        sort()
    } else {
      selker_thresh(a = input$alpha, b = input$beta, n_resp = input$n_thresh+1)
    }
    
    # estimate a & b
    selker_optim(actual_locs, method = input$optim_method)
  })
  
  # generate plots
  output$latent_plot <- renderPlot({
    req(input$thresh_1, input$n_thresh, input$thresh_set_method)
    plot_latent(opt_res())
  })
  
  output$distort_plot <- renderPlot({
    req(input$thresh_1, input$n_thresh, input$thresh_set_method)
    plot_distort(opt_res())
  })
  
}


shinyApp(ui, server)