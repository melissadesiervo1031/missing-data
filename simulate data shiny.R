# Shiny Application Skeleton
# Heili Lowman
# January 14, 2020

# Shiny application = Packages/Data + User Interface + Server + shinyApp()

#### Setup ####

# Anything that should happen once, and only once, upon launching the app should 
# be placed in the Setup section (e.g., data load-in and processing). All interactive 
# elements will be in the following sections.

# Load packages -----------------------------------------------------------

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(plyr)
library(shinystan)
library(faux)
library(bayesplot)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(lubridate)
library(data.table)
library(Amelia)
library(tidyverse)
library(shiny)

# Simulate light data -----------------------------------------------------------

N<-365 #length of data
z<-numeric(N+1) 
t<-numeric(N+1) 
time<-seq(as.POSIXct("2020/1/1"),as.POSIXct("2020/12/30"), by="day") #for light


##light
# From Yard et al. (1995) Ecological Modelling.  Remember your trig?  
# calculate light as umol photon m-2 s-1.
# Arguments are:  
# time = a date and time input (posixct object)
# lat = latitude of field site
# longobs = longitude of field site
# longstd = standard longitude of the field site (NOTE: watch daylight savings time!!!). For PST, longstd is be 120 degrees. But during PDT it is 105 degrees. MST is 105 deg. MDT is 90. 


# convert degrees to radians
radi<-function(degrees){(degrees*pi/180)}

# function to estimate light
lightest<- function (time, lat, longobs, longstd) {
  jday<-yday(time)
  E<- 9.87*sin(radi((720*(jday-81))/365)) - 7.53*cos(radi((360*(jday-81))/365)) - 1.5*sin(radi((360*(jday-81))/365))
  LST<-as.numeric(time-trunc(time))
  ST<-LST+(3.989/1440)*(longstd-longobs)+E/1440
  solardel<- 23.439*sin(radi(360*((283+jday)/365)))
  hourangle<-(0.5-ST)*360
  theta<- acos( sin(radi(solardel)) * sin(radi(lat)) + cos(radi(solardel)) * cos(radi(lat)) * cos(radi(hourangle)) )
  suncos<-ifelse(cos(theta)<0, 0, cos(theta))
  GI<- suncos*2326
  GI	
  
}
light<-lightest(time, 47.8762, -114.03, 105) #Flathead Bio Station just for fun
light.l<-log(light)
light.c<-(light-mean(light))/sd(light)
light.rel<-light/max(light)


# Set range and interval of each variable
sdp<- seq(from=0.001, to=.1, length.out=10)
sdo<-seq(from=0.001, to=.1, length.out=10)
phi<-seq(from=0, to=1, length.out=10)
b0<-seq(from=0, to=3, length.out=10)
b1<-seq(from=0, to=1, length.out=10)
data<-data.frame(sdp, sdo, phi, b0, b1)


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
      
      # widget #1 - phi slider menu
      sliderInput("phi", "Phi:",
                   min=0, max=1, value=0.8, step=0.1,
                  animate=animationOptions(100)),
             
    # widget #2 - sdp slider menu
    sliderInput("sdp", "Sigma_proc:",
                min=0.001, max=.1, value=0.01, step=0.005,
                animate=animationOptions(100)),
    
    # widget #3 - sdo slider menu menu
  #  sliderInput("sdo", "Decimal:",
   #             min=0.001, max=.1, value=0.01, step=0.005,
    #            animate=animationOptions(100)),
   
     # widget #4 - b0 slider menu
    sliderInput("b0", "Intercept:",
                min=0, max=3, value=0, step=0.5,
                animate=animationOptions(100)),
  
      # widget #5 - b1 slider menu
    sliderInput("b1", "Light Beta:",
                min=0, max=1, value=.5, step=0.1,
                animate=animationOptions(100)),
    
  ),
    
  mainPanel(
    
    # Output: Table summarizing the values entered ----
    tableOutput("values"),
    plotOutput(outputId = "time_series_plot"),
    p("Penguin Plot:")
  )
  ),
  )

## TS function
TS<-function(phi, sdp, b0, b1){
  set.seed(4)
  ts<-NA
  ts[1]<-(b0+b1)/(1-phi) #Expected value ts[1] varies with variable shoice
  for (i in 1:364){
    ts[i+1] <- b0+phi*ts[i]+light.rel[i]*b1+rnorm(1, 0, sdp)
  }
  day<-seq(from=1, to=365, by=1)
  dat<-data.frame(ts, day)
  return(list("dat"=dat))
}


#### Server ####

# The server is everything that happens "behind the curtain". This is 
# where all the data transformation, plot creation, and simulations might
# take place before being sent back to the user interface to be
# displayed. It is structured grammatically like most R functions, with all
# inputs supplied by variables in the user interface, and all outputs
# created in the server and supplied back to the user interface.

server <- function(input, output){
  
 # create a new dataset based on the user's selection above
  # this `reactive()` tells shiny to monitor the user's choices
   time_series<- reactive({ # this `reactive()` tells shiny to monitor the user's choices
    sim<-TS(input$phi,input$sdp, input$b0, input$b1)
    sim$dat
    })
  
   
  # create a plot
  output$time_series_plot <- renderPlot({
    ts<-req(time_series())    
    xvar<-ts[,2]
    yvar<-ts[,1]
    fin<-data.frame(xvar, yvar)
    
    ggplot(fin, aes(x = xvar, y=yvar)) + # the dataset is presented here as
      # `data()` instead of `data` to indicate that it's a reactive and malleable dataset
      # so shiny also monitors and changes it based on the user's choices
      geom_point() +
      theme_classic()+
      theme(axis.title.x=element_text(size=18,colour = "black"))+
      theme(axis.title.y=element_text(size=18,colour = "black"))+
      theme(axis.text.y=element_text(size=18,colour = "black"))+
      theme(axis.text.x=element_text(size=18,colour = "black"))+
      theme(legend.position="top")+
      ylab("Data")+
      xlab("Time (days)")+
            theme(legend.position = "none")
  })
  
}

#### Combine the user interface and server ####
#shinyApp(ui = ui, server = server)
if (interactive()) shinyApp(ui, server)
# Additional Shiny resources:
# https://shiny.rstudio.com/gallery/widget-gallery.html
# http://shinyapps.dreamrs.fr/shinyWidgets/
# https://shiny.rstudio.com/gallery/#user-showcase

# End of script.
