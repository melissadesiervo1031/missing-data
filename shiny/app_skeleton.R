# Shiny Application Skeleton
# Heili Lowman
# January 14, 2020

# Shiny application = Packages/Data + User Interface + Server + shinyApp()

#### Setup ####

# Anything that should happen once, and only once, upon launching the app should 
# be placed in the Setup section (e.g., data load-in and processing). All interactive 
# elements will be in the following sections.

# Load packages.
library(tidyverse)
library(shiny)
library(palmerpenguins) # contains sample dataset "penguins"

# Load data.
dat <- penguins

#### User interface ####

# This is the section of the code that builds what the user may interact with.
# "Widgets" are the term for the different items (buttons, slide bars, etc.)
# that the user selects.

# There are also two main ways to set up a UI. Below, I've created a more standard
# layout that displays as an html webpage, but this can also be reconfigured to
# look more like a dashboard. More info at https://rstudio.github.io/shinydashboard/

ui <- fluidPage( # allows the user to scroll through the page

  titlePanel("Penguin Shiny App"), # title
  
  sidebarLayout( # creates a sidebar
    
    sidebarPanel("Widgets go here.", # sidebar title
      
      # widget #1 - dropdown menu
      selectInput(inputId = "species_select", # variable name
        label = "Select a species.", # title for user
        choices = unique(dat$species) # choices displayed
        )
      
      ),
    
    mainPanel("Outputs go here.", # main panel title
      
      # output #1 - plot of penguin counts by island
      p("Penguin Plot:"), # the p() command creates a new labeled paragraph
      plotOutput(outputId = "penguin_plot")
      
      )
  )
)

#### Server ####

# The server is everything that happens "behind the curtain". This is 
# where all the data transformation, plot creation, and simulations might
# take place before being sent back to the user interface to be
# displayed. It is structured grammatically like most R functions, with all
# inputs supplied by variables in the user interface, and all outputs
# created in the server and supplied back to the user interface.

server <- function(input, output){
  
  # create a new dataset based on the user's selection above
  penguin_species <- reactive({ # this `reactive()` tells shiny to monitor the user's choices
    dat %>%
      filter(species == input$species_select)
  })
  
  # create a plot
  output$penguin_plot <- renderPlot({ # there are lots of different `render` functions
    # available depending on what you'd like to create and spit out in the user interface
    ggplot(penguin_species(), aes(x = island)) + # the dataset is presented here as
      # `data()` instead of `data` to indicate that it's a reactive and malleable dataset
      # so shiny also monitors and changes it based on the user's choices
      geom_histogram(stat = "count", aes(fill=island)) +
      theme_minimal() +
      theme(legend.position = "none")
  })
  
}

#### Combine the user interface and server ####
shinyApp(ui = ui, server = server)

# Additional Shiny resources:
# https://shiny.rstudio.com/gallery/widget-gallery.html
# http://shinyapps.dreamrs.fr/shinyWidgets/
# https://shiny.rstudio.com/gallery/#user-showcase

# End of script.
