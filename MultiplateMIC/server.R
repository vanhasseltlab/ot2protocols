library(shiny)
library(xlsx)
source("Source.R")


#SERVER MAIN------------
shinyServer(function(input, output) {
  #defining directory-------
  #outputDir_cmdline <- "C:\\Users\\Sebastian\\Desktop\\MSc Leiden 2nd Year\\##LabAst Works\\FinVersion\\OT2_UpstreamExcelProc_MAIN\\Checkpoint - Postpatch\\20201016_AddingPlateMap to Output\\MultiplateMIC_Uploader_v2"
  #outputDir_usrGuide <- "C:\\Users\\Sebastian\\Desktop\\MSc Leiden 2nd Year\\##LabAst Works\\FinVersion\\OT2_UpstreamExcelProc_MAIN\\Checkpoint - Postpatch\\20201016_AddingPlateMap to Output\\MultiplateMIC_Uploader_v2"
  #inputTemplate <- "C:\\Users\\Sebastian\\Desktop\\MSc Leiden 2nd Year\\##LabAst Works\\FinVersion\\OT2_UpstreamExcelProc_MAIN\\Checkpoint - Postpatch\\20201016_AddingPlateMap to Output\\MultiplateMIC_Uploader_v2\\MultiplateMIC_InputTemplate.xlsx"
  outputDir_cmdline <- "/srv/shiny-server/files/Output_CmdList"
  outputDir_usrGuide <- "/srv/shiny-server/files/Output_UsrGuide"
  inputTemplate <- "/srv/shiny-server/MultiplateMIC/MultiplateMIC_InputTemplate.xlsx"
  #outputDir_cmdline <- "/home/sebastian/20201103_6MICtest"
  #outputDir_usrGuide <- "/home/sebastian/20201103_6MICtest"
  #inputTemplate <- "/home/sebastian/20201103_6MICtest"
  
  #Obtain names---------
  new_name <- reactive({
    if(is.null(input$pmid)){pmid <- ''}else{pmid <- input$pmid}
    if(is.null(input$exp_name)){expName <- ''}else{expName <- input$exp_name}
    if(is.null(input$exp_num)){expNum <- ''}else{expNum <- input$exp_num}
    if(is.null(input$f_name)){fName <- ''}else{fName <- input$f_name}
    if(is.null(input$l_name)){lName <- ''}else{lName <- input$l_name}
    
    #paste
    res <- paste("PMID-", pmid, "_EXPID-", 
                 expName, "-", expNum, "_", lName, ".", fName, sep='')
    
    return(res)
  })
  
  #Confirming Upload File-----------
  contents <- reactive({
    infile = input$file
    if(is.null(infile)){return(NULL)}
    
    if(input$do==0){
      dis <- read.xlsx(infile$datapath, 1, rowIndex = c(56:63), colIndex=c(2:13), header=F)
    }else{
      #rename files for safekeeping
      file_name <<- strsplit(infile$name, '.xl')[[1]][1]
      #file.copy(infile$datapath, homeDir)
      #file.rename(paste(homeDir, "0.xlsx", sep='\\'), paste("UserInput_", file_name, ".xlsx", sep=""))
      
      #update table view
      dis <- main(infile$datapath)
      
      #savekeeping output files
      #command line
      cmdLine_name <- paste("CommandList_", new_name(), '.csv', sep='')
      write_dir <- paste(outputDir_cmdline, cmdLine_name, sep='/')
      write.csv(OT2_cmd, write_dir, row.names = FALSE)
      
      #user guide
      usrGuide_name <- paste("RobotHandler_", new_name(), '.xlsx', sep='')
      write_dir <- paste(outputDir_usrGuide, usrGuide_name, sep='/')
      write.xlsx(user_cmd, write_dir, row.names = FALSE, col.names=T)
    }
    return(dis)
  })
  
  output$tab <- renderTable({contents()})
  
  #Enabling download button-------
  output$downloadData <- renderUI({
    req(input$do, contents())
    downloadButton("d_OT2", "Download Robot Commands")
  })
  output$downloadData2 <- renderUI({
    req(input$do, contents())
    downloadButton("guide", "Download Robot Setup Guide")
  })
  
  #Sample File Name--------
  output$tex <- renderText({new_name()})
  
  #Defining download buttons--------
  output$d_OT2 <- downloadHandler(
    filename = function(){paste("CommandList_", new_name(), '.csv', sep='')},
    content = function(file) {
      write.csv(OT2_cmd, file, row.names = FALSE)
    }
  )
  output$guide <- downloadHandler(
    filename = function(){paste("RobotHandler_", new_name(), '.xlsx', sep='')},
    content = function(file) {
      write.xlsx(user_cmd, file, row.names = FALSE, col.names=T)
    }
  )
  
  #download input template
  output$downloadTemplate <- downloadHandler(
    filename = "MultiplateMIC_InputTemplate.xlsx",
    content = function(file) {
      file.copy(inputTemplate, file)
    }
  
  )
})
