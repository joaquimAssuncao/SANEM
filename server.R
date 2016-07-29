#################### SAN Generator (SANGe)- shinny server ##########################
# 2015 by Joaquim Assuncao (www.joaquim.pro.br)
# jassuncao@outlook.com
#####################################################

library(shiny)
source("SANGE.r")
source("tools.r")

shinyServer(function(input, output) {
     output$contents <- renderTable({
          
          # input$file1 will be NULL initially. After the user selects
          # and uploads a file, it will be a data frame with 'name',
          # 'size', 'type', and 'datapath' columns. The 'datapath'
          # column will contain the local filenames where the data can
          # be found.
          
          inFile <- input$file1
          
          if (is.null(inFile))
               return(NULL)
          
          usr_data <<- read.csv(inFile$datapath, header=input$header, sep=input$sep, 
                   quote=input$quote)
          
          #Na to 0
          if (input$naTo==TRUE){
               usr_data[is.na(usr_data)] <<- 0
          }
          #Na as previous
          if (input$naAsPrev==TRUE){    
               usr_data <<- na.lomf(usr_data)
              
          }
          

          
          #Normalize
          if (input$norm==TRUE){
               for (c in 1:ncol(usr_data)){
                    if (is.numeric(usr_data[,c])){
                         for (l in 1:nrow(usr_data)){
                              #used to avoid NAs
                              if (!is.na(usr_data[l,c])){
                                   usr_data[l,c] <<- usr_data[l,c]/max(usr_data[!is.na(usr_data[,c]),c])  
                              }   
                         }  
                    }
               }
          }
          
          usr_data
     })
     

     generateSAN <- eventReactive(input$generatesan, {
#           for (i in 1:ncol(usr_data)){
#                nums[i] <- is.numeric(usr_data[,i])
#           }
          selCols  <- as.numeric(strsplit(input$selCols, ",")[[1]]) 
          if (input$intFunc!="A,A,A..."){
               intFunc  <- (strsplit(input$intFunc, ",")[[1]])  
          }else{
               intFunc <- input$intFunc
          }
              
          
          SANEM(usr_data[,selCols], input$sanname, sync= input$syncFact, 
                alphz=input$alphalevel,chron=input$chronS,intFunc=intFunc)
          
     })
     
     #get the quality of the model
     modelMeasures <- eventReactive(input$generateMeasurements, {
          selCols  <- as.numeric(strsplit(input$selCols, ",")[[1]]) 
          tmp <- GetFittingsMeasurements(usr_data[,selCols],1)
          HTML(paste("","mllk: ", tmp[[1]], "avgMatch: ", tmp[[2]], "DTWdist: ",tmp[[3]],sep = '<br/>'))
     })
    
     
     output$downloadData <- downloadHandler(
          filename = function() { 
               paste("usr_data",'.csv', sep='') 
          },
          content = function(file) {
               write.csv(usr_data, file)
          }
     )
     
     
     output$sanCode <- renderText({
          generateSAN()
     })
     output$mMeasures <- renderText({
          modelMeasures()
     })
    
     
    
     
})
