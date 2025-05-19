library(shiny)
library(circular)
library(ggplot2)
library(lattice)
# library(extrafont)
library(tidyverse)
library(gridExtra)
library(readxl)
library(reticulate)
library(data.table)
# library(showtext)




source_python("HS_CJPlot.py")
source_python("HS_savePlotSVG.py")

labelBio <- c("A", "B", "C", "D")
labelPhys <- c("0","90","180","270")


meanDistance.func <- function(funVar_df){
  funVar_df_rdcd <- (funVar_df - 0.5) %>% na.omit()
  
  funVar_df_devSqrd <- 
    data.frame(
      x_Sqrd = (funVar_df_rdcd$x - colMeans(funVar_df_rdcd)[1])^2,
      y_Sqrd = (funVar_df_rdcd$y - colMeans(funVar_df_rdcd)[2])^2
    )
  
  funVar_dist <- rowSums(funVar_df_devSqrd) %>% sqrt()
  
  funVar_meanDist <- mean(funVar_dist)
  
  return(funVar_meanDist)
}


# Section: UI -------------------------------------------------------------


## -----UI Side----- ##

ui <- 
  fluidPage(
    tags$head(tags$link(rel = "stylesheet",
                        type = "text/css", href = "styles.css")),
    titlePanel(title = "Plotting circular histograms"),
    sidebarLayout( 
      sidebarPanel(
        fileInput(
          inputId = "file",
          label = "Select two files",
          multiple = T
        ),
        sliderInput(
          inputId = "pointSize",
          label = "Size of Scatter Points",
          min = 10,
          max = 180,
          value = 90,
          step = 10
        ),
        sliderInput(
          inputId = "levelsKDE",
          label = "KDE plot levels",
          min = 0,
          max = 10,
          value = 2,
          step = 1
        ),
        sliderInput(
          inputId = "threshKDE",
          label = "KDE plot threshold",
          min = 0,
          max = 1,
          value = 0.8,
          step = 0.1
        ),
        selectInput(
          inputId = "arrowSize",
          label = "Select Arrow Size",
          choices = c("S", "M", "L")
        ),
        sliderInput(
          inputId = "alphaC",
          label = "Alpha of the contour",
          min = 0,
          max = 1.0,
          value = 0.5,
          step = 0.1
        ),
        sliderInput(
          inputId = "alphaP",
          label = "Alpha of the point",
          min = 0,
          max = 1.0,
          value = 0.5,
          step = 0.1
        ),
        div(
          tags$img(style = "height:50% ; width:50%; border-radius: 50%  ",
                   src = "MonaLisa.jpg"),
          style = "text-align: center;"
          
        )
      ),
      mainPanel(
        uiOutput("tb")
      )
    )
  )


# Section: Server ---------------------------------------------------------


## -----Server Side----- ##

server <-
  function(input, output){
    
    data<-reactive({
      tmp <- 
        lapply(input$file$datapath, 
               FUN =
                 function(fA_filePath)
                 {
                   fread(
                     input = fA_filePath,
                     col.names = c("x","y")
                   )
                 }  
        )
      names(tmp) <- input$file$name
      tmp
      
    })
    
    
    data_download <- function(){
      if(is.null(data())){return()}
      
      fileNames <- names(data()) %>%
        str_remove(pattern = ".csv")
      
      
      cartData1 <- data()[[1]]
      names(cartData1) <- c("X", "Y")
      cartData2 <- data()[[2]]
      names(cartData2) <- c("X", "Y")
      
      SumTable <- cartData1 %>%
        as.data.frame() %>%
        mutate(Type = fileNames[1]) %>%
        bind_rows(cartData2) %>%
        mutate(Type = if_else(is.na(Type), true = fileNames[2], false = Type)) %>%
        filter(!is.na(X), !is.na(Y)) %>%
        mutate(r = sqrt((X-0.5)^2 + (Y-0.5)^2),
               theta = 180*atan2((Y-0.5),(X-0.5))/pi) %>%
        group_by(Type) %>%
        mutate(Size = n()) %>%
        mutate(meanTheta = mean(theta, na.rm = T),
               meanR = mean(r , na.rm = T)) %>%
        mutate(circVar_r1 = sqrt((sum(cospi(theta/180)))^2 + 
                                   (sum(sinpi(theta/180)))^2),
        ) %>%
        mutate(circVar_r1av = circVar_r1/Size,
               circVar_cv1 = 1- circVar_r1av) %>%
        mutate(circSd_usingR = sqrt(-2*log(circVar_r1av))*180/pi)%>%
        mutate(Separator = "||") %>%
        mutate(meanX = mean(X, na.rm = T),
               meanY = mean(Y, na.rm = T),
               meanDistFromMean = sqrt(mean((X-meanX)^2 + (Y-meanY)^2, na.rm = T))
        ) %>%
        select(!c(X,Y,r,theta)) %>%
        distinct()
      
      
      
      summary_table_download1 <<- SumTable
      return(summary_table_download1)
    }
    
    output$filedf <- renderTable({
      if(is.null(data())){return()}
      data()[[1]]
      # class(data())
      # summary_table
    })

# Table Summary -----------------------------------------------------------

    
    output$sum <- renderTable({
      if(is.null(data())){return()}
      fileNames <- names(data()) %>%
        str_remove(pattern = ".csv")
      
      cartData1 <- data()[[1]]
      names(cartData1) <- c("X", "Y")
      cartData2 <- data()[[2]]
      names(cartData2) <- c("X", "Y")
      
      cartData1 %>%
        as.data.frame() %>%
        mutate(Type = fileNames[1]) %>%
        bind_rows(cartData2) %>%
        mutate(Type = if_else(is.na(Type), true = fileNames[2], false = Type)) %>%
        filter(!is.na(X), !is.na(Y)) %>%
        mutate(r = sqrt((X-0.5)^2 + (Y-0.5)^2),
               theta = 180*atan2((Y-0.5),(X-0.5))/pi) -> CircData_comb
      
      CircData_comb %>%
        distinct(Type) %>%
        unlist(use.names = F) -> CircData_comb_types

      CircData_comb %>%
        filter(Type == CircData_comb_types[1]) %>%
        unlist(use.names= F)-> CircData_angles_type1
      # 
      # Circstd_type1 = py_circstd(CircData_angles_type1)
      # 
      # CircData_comb %>%
      #   filter(Type == CircData_comb_types[2]) %>%
      #   unlist(use.names= F)-> CircData_angles_type2
      # 
      # Circstd_type2 = py_circstd(CircData_angles_type2)
      # 
      # data.frame(py_circstd = c(Circstd_type1, Circstd_type2),
      #            Type = CircData_comb_types) -> CircStd_comb_df
      
      
      SumTable <-   
        CircData_comb %>%
        group_by(Type) %>%
        mutate(Size = n()) %>%
        mutate(meanTheta = mean(theta, na.rm = T),
               meanR = mean(r , na.rm = T)) %>%
        mutate(circVar_r1 = sqrt((sum(cospi(theta/180)))^2 + 
                                   (sum(sinpi(theta/180)))^2),
        ) %>%
        mutate(circVar_r1av = circVar_r1/Size,
               circVar_cv1 = 1- circVar_r1av) %>%
        mutate(circSd_usingR = sqrt(-2*log(circVar_r1av))*180/pi)%>%
        mutate(Separator = "||") %>%
        mutate(meanX = mean(X, na.rm = T),
               meanY = mean(Y, na.rm = T),
               meanDistFromMean = sqrt(mean((X-meanX)^2 + (Y-meanY)^2, na.rm = T))
        ) %>%
        select(!c(X,Y,r,theta)) %>%
        distinct()
      
      
      
      summary_table <<- SumTable
        
      
      
      
    })
    
    
    output$downStatsTable_csv  <- downloadHandler(
      # Create the download file name
      filename = function() {
        
        fun_names <-
          names(data()) %>%
          str_remove(pattern = ".csv")
        
        fun_filename <-
          paste(fun_names[1], fun_names[2], sep = "_")%>%
          paste0("Table_",., ".csv")
        
      },
      content = function(file_table) {
        # write.table(data_download(), file = file_table, dec = "." , sep = ";",  fileEncoding = "UTF-16LE")
        write_csv(data_download(), file = file_table)
      },
      contentType = "text/csv"
    )
    
    
    output$table <- renderTable({
      if(is.null(data())){return()}
      cartData <- data()[[1]]
      polarR <- sqrt((cartData$x-0.5)^2 + (cartData$y-0.5)^2)
      polarTheta <- 180*atan2((cartData$y-0.5),(cartData$x-0.5))/pi
      polarData <- data.frame(theta = polarTheta, r = polarR)
      
      polarData %>% na.omit()
    })
    
    
    

    # Polar Plot -----------------------------------------------------

    
    output$plot2 <- 
      renderPlot(res = 72, height ="auto", width = "auto", execOnResize = T , { 
        
        if(is.null(data())){return()}
        else{
          
          cartData1 <- data()[[1]]
          names(cartData1) <- c("X", "Y")
          cartData2 <- data()[[2]]
          names(cartData2) <- c("X", "Y")
          
          cartData12_list <- list()
          cartData12_list[[1]] <- cartData1
          cartData12_list[[2]] <- cartData2
          
          CircPlot <<- plot_CJ(cartData12_list,
                               input$pointSize, 
                               input$arrowSize, 
                               input$alphaC, 
                               input$alphaP,
                               input$levelsKDE,
                               input$threshKDE)
          printPlot(CircPlot)

      }
        })
    output$downPlot2_svg <- downloadHandler(
      filename = function() {
        # fun_filename <- "Trial.svg"
        
        fun_names <-
          names(data()) %>%
          str_remove(pattern = ".csv")
        
        fun_filename <-
          paste(fun_names[1], fun_names[2], sep = "_")%>%
          paste0("Plot_",., ".svg")
        
        },
      
      content = function(file){
        py$FCsvg(CircPlot)$print_svg(file)
      },
      
      contentType = "image/svg"
    )
    
    
    
    output$tb <- renderUI({
      if(is.null(data())){
        tags$iframe(style = "height: 500px; width: 100%",
                    src = "")
      }
      else{
        tags$iframe(style = "height: 10%; width: 100%",
                    src = "cat.png")
        tabsetPanel(
          tabPanel("Data in Cartesian",
                   tableOutput("filedf")
          ),
          
        tabPanel("Data in Polar",
                 tableOutput("table")),
        tabPanel("Summary",
                 tableOutput("sum"),
                 downloadButton(outputId = "downStatsTable_csv",
                                label = "Download Table in CSV"),
        ),
        tabPanel("Polar Plot",
                 plotOutput(width = "auto", height = "600px" , outputId = "plot2"),                                # 
                 downloadButton(outputId = "downPlot2_svg",
                                label = "Download SVG"),
                 # downloadButton(outputId = "downPlot2_png",
                 #                label = "Download PNG")
        )
        )
      }
    })
  }

# options(shiny.host = "172.16.232.222")
# options(shiny.host = "192.168.29.71")
# options(shiny.port = 7775)


# options(shiny.host = "")

shinyApp(ui, server)





# Rough -------------------------------------------------------------------

# cartData <- data()[[1]] %>% as.data.frame() %>% na.omit()
# polarR <- sqrt((cartData$x-0.5)^2 + (cartData$y-0.5)^2)
# polarTheta <- 180*atan2((cartData$y-0.5),(cartData$x-0.5))/pi
# polarData <- data.frame(theta = polarTheta, r = polarR)
# 
# circVar_R <- sqrt((sum(cospi(na.omit(polarTheta)/180)))^2 + (sum(sinpi(na.omit(polarTheta)/180)))^2)
# circVar_Rav <- circVar_R/nrow(polarData)
# 
# circVar_CV <- 1 - circVar_Rav
# 
# cartData_devSqrd <- 
#   data.frame(
#     x_Sqrd = (cartData$x - colMeans(cartData,na.rm = T)[1])^2,
#     y_Sqrd = (cartData$y - colMeans(cartData,na.rm = T)[2])^2
#   )
# 
# cartData_dist <- rowSums(cartData_devSqrd) %>% sqrt()


# data.frame(size = length(polarTheta %>% na.omit()),
#            meanTheta = mean(polarTheta, na.rm = T),
#            meanPolarR = mean(polarR, na.rm=T),
#            circVar = circVar_CV,
#            circVar_Rav = circVar_Rav,
#            circVar_R = circVar_R,
#            
#            Separator = NA,
#            
#            meanX = colMeans(cartData, na.rm = T)[1],
#            meanY = colMeans(cartData, na.rm = T)[2],
#            meanDistanceFromMean = mean(cartData_dist, na.rm = T)
#            
# )
