library(shiny)

shinyUI(
  pageWithSidebar(
    headerPanel("Multiplate MIC - OT2-Commander"),
    sidebarPanel(
      fileInput("file", "Upload Plate Map", accept=".xlsx"),
      downloadButton("downloadTemplate", label = "Template Input"),
      textInput("pmid", 'Plate Map ID (PMID)'),
      textInput("f_name", 'First Name'),
      textInput("l_name", 'Last Name'),
      textInput("exp_name", 'Experiment Name'),
      textInput("exp_num", 'Experiment Number'),
      textOutput('tex'),
      actionButton("do", "Confirm uploaded file and save"),
      uiOutput('downloadData'),
      uiOutput('downloadData2')
    ),
    mainPanel(
      tableOutput('tab')
    )
))
