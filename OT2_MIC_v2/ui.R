library(shiny)

shinyUI(
    pageWithSidebar(
        headerPanel("MIC Test - OT2-Commander"),
        sidebarPanel(
            fileInput("file", "Upload Plate Map", accept=".xlsx"),
            actionButton("do", "Confirm Uploaded File")
        ),
        mainPanel(
            tableOutput('tab')
        )
    )
)
