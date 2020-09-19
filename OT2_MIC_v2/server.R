library(shiny)
library(xlsx)

StockInf_cr <- function(stock_readout, finAmt){
  #remove NAs
  stock_readout <- stock_readout[colSums(!is.na(stock_readout)) > 0]
  DrugName <- c()
  Conc <- c()
  for(i in c(2:length(stock_readout[1,]))){
    DrugName <- c(DrugName, toString(stock_readout[1,i]))
    Conc <- c(Conc, toString(stock_readout[2,i]))
  }
  stock_info <- cbind.data.frame(DrugName, Conc)
  return(stock_info)
}
PlateInf_cr <- function(filename, finAmt, stockInfo){
  drug_inf <- as.vector(as.matrix(read.xlsx(filename, 1, rowIndex = c(9:16), colIndex=c(2:13), header=F)))
  conc_inf <- as.vector(as.matrix(read.xlsx(filename, 1, rowIndex = c(20:27), colIndex=c(2:13), header=F)))
  solv_inf <- as.vector(as.matrix(read.xlsx(filename, 1, rowIndex = c(32:39), colIndex=c(2:13), header=F)))
  inoc_inf <- as.vector(as.matrix(read.xlsx(filename, 1, rowIndex = c(44:51), colIndex=c(2:13), header=F)))
  
  wells <- c()
  for(i in c(1:8)){
    for(j in c(1:12)){
      wells <- c(wells, paste(LETTERS[i], toString(j), sep=""))
    }
  }
  
  plateinfo <- cbind.data.frame(wells, drug_inf, conc_inf, solv_inf, inoc_inf)
  colnames(plateinfo) <- c('Well', 'Drug', 'Conc', 'Solvent', 'Inoc')
  return(plateinfo)
}
solList_cr <- function(plate_info, finAmt, stockInfo){
  #remove water-only
  sol_list <- plate_info[!(is.na(plate_info$Drug)),]
  
  #get solution IDs
  solIDs <- c()
  for(i in c(1:length(sol_list[,1]))){
    solIDs <- c(solIDs, paste(sol_list[i,2], sol_list[i,3], sol_list[i,4], sep="_"))
  }
  sol_list <- cbind.data.frame(solIDs, sol_list)
  
  #iterate through all drug and solvent combinations; get required amount in well
  new_sol_list <- c()
  solvents <- unique(sol_list$Solvent)
  drugs <- unique(sol_list$Drug)
  for(i in c(1:length(solvents))){
    for(j in c(1:length(drugs))){
      curDat <- subset(sol_list, Drug==drugs[j] & Solvent==solvents[i])
      
      for(m in c(1:length(curDat[,1]))){
        #print(curDat[m,1] %in% new_sol_list[,1])
        if(curDat[m,1] %in% new_sol_list[,1]){
          new_sol_list[(new_sol_list[,1]==curDat[m,1]) ,7] <- 
            new_sol_list[(new_sol_list[,1]==curDat[m,1]),7]+finAmt
        }else{
          nexItem <- cbind(curDat[m,], finAmt)
          new_sol_list <- rbind(new_sol_list,
                                nexItem)
        }
      }
    }
  }
  
  sol_list <- new_sol_list #concatenate for next iteration
  sol_list$finAmt <- sol_list$finAmt+150 #add excess
  new_sol_list <- c()
  #iterate through all drug and solvent combinations; get required amount in well and for dilutions
  for(i in c(1:length(solvents))){
    for(j in c(1:length(drugs))){
      curDat <- subset(sol_list, Drug==drugs[j] & Solvent==solvents[i])
      curDat <- curDat[order(curDat$Conc, decreasing=F),]
      
      req_forLower <- c(0) #always begin from lowest
      
      if(curDat[1,7]>1200){
        dilTube <- c('Dil_15')
      }else{
        dilTube <- c('Dil_1.5')
      }
      
      for(m in c(2:length(curDat[,1]))){
        pass <- curDat[(m-1),4]*curDat[(m-1),7]/curDat[m,4] #amount required for diluting one above
        req_forLower <- c(req_forLower, pass)
        curDat[m, 7] <- curDat[m,7]+pass
        if(curDat[m,7]>1200){
          dilTube <- c(dilTube, 'Dil_15')
        }else{
          dilTube <- c(dilTube, 'Dil_1.5')
        }
      }
      higherV <- c(req_forLower[c(2:length(req_forLower))],
                   curDat[length(curDat[,1]),7]*curDat[length(curDat[,1]),4]/
                     as.numeric(stockInfo$Conc[stockInfo$DrugName==curDat[1,3]]))
      nexItem <- cbind(curDat, req_forLower, dilTube, higherV)
      new_sol_list <- rbind(new_sol_list, nexItem)
    }
  }
  return(new_sol_list)
}
PlateMap_cr <- function(plate_info, finamts, stocks, plates){
  #get solution IDs
  solIDs <- c()
  for(i in c(1:length(plate_info[,1]))){
    solIDs <- c(solIDs, paste(plate_info[i,2], plate_info[i,3], plate_info[i,4], sep="_"))
  }
  plate_info <- cbind.data.frame(solIDs, plate_info)
  return(plate_info)
}
main <- function(fileName){
  #LOAD INPUT-----------
  finAmt <- as.numeric(read.xlsx(fileName, 1, rowIndex=c(5), colIndex = c(3), header=F))
  stockInfo <- StockInf_cr(read.xlsx(fileName, 1, rowIndex=c(1:3), header=F, stringsAsFactors=F), finAmt)
  plateInfo <- PlateInf_cr(fileName, finAmt, stockInfo)
  stockInfo[,2] <- as.numeric(levels(stockInfo[,2]))[stockInfo[,2]]
  cmdList <- c()
  #MAP INITIATION-----------
  #1. Deck Map
  labwares <- c() #create labware positions
  for(i in c(1:12)){
    labwares <- c(labwares, paste("labware_", toString(i), sep=''))
  }
  deckMap <- cbind.data.frame(labwares,
                              c("tip_300", "96-wellA", "96-wellB",
                                "tip_300", "Dil_1.5", "Dil_15",
                                "tip_300", "Sol_50", "Stock",
                                "", "", "TRASH"))
  colnames(deckMap) <- c("position", "item") #create deck map
  deckMap <- deckMap
  
  #2. Stock Map
  stockMap <- c()
  stockCoord <- c(1,1)
  
  for(i in c(1:length(stockInfo[,1]))){
    nexItem <- c(paste(LETTERS[stockCoord[1]], toString(stockCoord[2]), sep=''), toString(stockInfo[i,1]))
    stockMap <- rbind(stockMap, nexItem)
    
    #increase coordinate
    stockCoord[2] <- stockCoord[2]+1
    if(stockCoord[2]>6){
      stockCoord[2] <- 1
      stockCoord[1] <- stockCoord[1]+1
    }
  }
  stockMap <- data.frame(stockMap)
  colnames(stockMap) <- c("Slot", "Fill")
  
  
  #3. Create a list of solutions
  solList <-  solList_cr(plateInfo, finAmt, stockInfo)
  
  #4. Solvent Map
  solvents <- unique(plateInfo$Solvent)
  solventCoord <- c(1,1)
  solventMap <- c()
  for(i in c(1:length(solvents))){
    nexItem <- c(paste(LETTERS[solventCoord[1]], toString(solventCoord[2]), sep=''), toString(solvents[i]))
    solventMap <- rbind(solventMap, nexItem)
    
    #increase coordinate
    solventCoord[2] <- solventCoord[2]+1
    if(solventCoord[2]>3){
      solventCoord[2] <- 1
      solventCoord[1] <- solventCoord[1]+1
    }
   
  }
  solventMap <- data.frame(solventMap)
  colnames(solventMap) <- c("Slot", "Fill")
  
  #------------------
  ############MAIN#################
  cmdList <- c()
  tipID <- 0
  #A. DISTRIBUTE SOLVENTS FOR INITIAL DILUTION--------------
  dilTube_large_coord <- c(1,1)
  dilTube_small_coord <- c(1,1)
  dilTube_large_map <- c()
  dilTube_small_map <- c()
  
  for(i in c(1:length(solList[,1]))){
    #select dilution rack
    target_ware <- toString(deckMap[(deckMap$item==toString(solList$dilTube[i])),1])
    if(solList[i,]$dilTube=="Dil_1.5"){
      coord <- dilTube_small_coord
      map <- dilTube_small_map
      max_row <- 6
    }else{
      coord <- dilTube_large_coord
      map <- dilTube_large_map
      max_row <- 5
    }
    
    #initiate command line
    cmd_line <- c(toString(deckMap[deckMap$item=="Sol_50",1]), #source labware
                  toString(solventMap[solventMap$Fill==toString(solList$Solvent[i]),1]), #source slot
                  target_ware, #target labware
                  paste(LETTERS[coord[1]], toString(coord[2]), sep=''), #target slot
                  toString(solList$finAmt[i]-solList$higherV[i]), "0", "0", "Initial solvent distribution")
    cmdList <- rbind(cmdList, cmd_line)
    
    #update map and coordinate
    map <- rbind(map, c(paste(LETTERS[coord[1]], toString(coord[2]), sep=''), toString(solList$solIDs[i])))
    coord[2] <- coord[2]+1
    if(coord[2]>max_row){
      coord[2] <- 1
      coord[1] <- coord[1]+1
    }
    
    #update back to main tube map/coord
    if(solList[i,]$dilTube=="Dil_1.5"){
      dilTube_small_coord <- coord
      dilTube_small_map <- map 
    }else{
      dilTube_large_coord <- coord
      dilTube_large_map <- map 
    }
  }
  dilTube_large_map <- cbind(dilTube_large_map, replicate(length(dilTube_large_map[,1]), "large"))
  dilTube_small_map <- cbind(dilTube_small_map, replicate(length(dilTube_small_map[,1]), "small"))
  dilTube_map <- rbind.data.frame(dilTube_large_map, dilTube_small_map)
  colnames(dilTube_map) <- c('slot', 'fill', 'TubeType')
  
  #B. DISTRIBUTE HIGHEST STOCK-------------
  #get list of highest concentrations // per-drug and solvent type
  solvents <- unique(solList$Solvent)
  drugs <- unique(solList$Drug)
  for(i in c(1:length(solvents))){
    for(j in c(1:length(drugs))){
      #get drug ID
      nexItem <- subset(solList, Solvent==toString(solvents[i]) & Drug==toString(drugs[j]))
      nexItem <- nexItem[nexItem$Conc==max(nexItem$Conc),]
      
      #get location
      source_location <- c(toString(deckMap[deckMap$item=="Stock",1]), #labware
                           toString(stockMap$Slot[stockMap$Fill==toString(drugs[j])])) #slot
      target_location <- c(toString(dilTube_map$TubeType[dilTube_map$fill==toString(nexItem$solIDs)]), #labware
                           toString(dilTube_map$slot[dilTube_map$fill==toString(nexItem[,1])])) #slot
      if(target_location[1]=="large"){
        target_location[1] <- toString(deckMap[deckMap$item=="Dil_15",1])
      }else{
        target_location[1] <- toString(deckMap[deckMap$item=="Dil_1.5",1])
      }
      #create next command line
      cmd_line <- c(source_location, target_location, nexItem$higherV, nexItem$higherV, 
                    0, "Initial stock dilution")
      cmdList <- rbind(cmdList, cmd_line)
    }
  }
  
  #C. MAIN DILUTION ON RACK-------------
  for(i in c(1:length(solvents))){
    for(j in c(1:length(drugs))){
      #subset current solution list
      cur_solList <- subset(solList, Drug==drugs[j] & Solvent==solvents[i])
      cur_solList <- cur_solList[order(cur_solList$Conc, decreasing=T),]
      #iterate through all rows
      for (m in c(2:length(cur_solList[,1]))){
        nexItem <- cur_solList[m,]
        sourceItem <- cur_solList[(m-1),]
        
        #get location
        source_location <- c(toString(dilTube_map$TubeType[dilTube_map$fill==sourceItem$solIDs]), #labware
                             toString(dilTube_map$slot[dilTube_map$fill==sourceItem$solIDs])) #slot
        target_location <- c(toString(dilTube_map$TubeType[dilTube_map$fill==nexItem$solIDs]), #labware
                             toString(dilTube_map$slot[dilTube_map$fill==nexItem$solIDs])) #slot
        
        if(target_location[1]=="large"){
          target_location[1] <- toString(deckMap[deckMap$item=="Dil_15",1])
        }else{
          target_location[1] <- toString(deckMap[deckMap$item=="Dil_1.5",1])
        }
        
        if(source_location[1]=="large"){
          source_location[1] <- toString(deckMap[deckMap$item=="Dil_15",1])
        }else{
          source_location[1] <- toString(deckMap[deckMap$item=="Dil_1.5",1])
        }
        
        #create command list
        if(nexItem$higherV > 0){
          cmd_line <- c(source_location, target_location, nexItem$higherV, nexItem$higherV,
                        0, "main dilution")
          cmdList <- rbind(cmdList, cmd_line)
        }
      }
    }
  }
  
  #D. DISTRIBUTING TO WELLS-------
  plateMap <- PlateMap_cr(plateInfo, finAmt, stockInfo, plateInfo)
  for(i in c(1:length(dilTube_map[,1]))){
    #get locations
    source_location <- c(toString(dilTube_map$TubeType[i]), toString(dilTube_map$slot[i]))
    target_location <- c(toString(deckMap$position[deckMap$item=='96-wellA']),
                         paste(plateMap$Well[(plateMap$solIDs==toString(dilTube_map$fill[i]))], collapse=', '))
    
    if(source_location[1]=="large"){
      source_location[1] <- toString(deckMap[deckMap$item=="Dil_15",1])
    }else{
      source_location[1] <- toString(deckMap[deckMap$item=="Dil_1.5",1])
    }
    
    #create command line
    tipID <- tipID + 1
    cmd_line <- c(source_location, target_location, finAmt, finAmt,
                  tipID, "final distribution")
    cmdList <- rbind(cmdList, cmd_line)
  }
  
  #E. FILLING EMPTY WELLS---------
  source_location <- c(toString(deckMap$position[deckMap$item=='Sol_50']), #getting solution rack
                       toString(solventMap$Slot[solventMap$Fill=="WATER"])) #getting solution slot
  target_location <- c(toString(deckMap$position[deckMap$item=='96-wellA']), #get 96-well plate location
                       paste(plateMap$Well[grepl("NA", plateMap$solIDs) & grepl("WATER", plateMap$solIDs)], collapse=', '))
  #create command line
  tipID <- tipID + 1
  cmd_line <- c(source_location, target_location, finAmt, finAmt,
                tipID, "filling empty wells")
  cmdList <- rbind(cmdList, cmd_line)
  
  #+1. Calculating Amount of Required Solvents------------
  cmdList <- data.frame(cmdList)
  colnames(cmdList) <- c('SourceWare', 'SourceWell', 'TargetWare', 'TargetWell', 'Amt', 'MixAmt', 'TipID', 'Comment')
  
  solventAmtList <- c()
  stockAmtList <- c()
  
  for(i in c(1:length(cmdList[,1]))){
    if(toString(cmdList$SourceWare[i])==toString(deckMap$position[deckMap$item=='Sol_50'])){
      nexAmt <- length(strsplit(toString(cmdList$TargetWell[i]), split=', ')[[1]])*as.numeric(toString(cmdList$Amt[i]))
      
      #get current source well
      source_well <- toString(cmdList$SourceWell[i])
      #check if current well is already listed in map
      if(source_well %in% solventAmtList[,1]){
        solventAmtList$Amt[solventAmtList[,1]==source_well] <- solventAmtList$Amt[solventAmtList[,1]==source_well]+nexAmt
      }else{
        nex_item <- cbind(solventMap[solventMap$Slot==source_well,], nexAmt)
        colnames(nex_item)[3] <- 'Amt'
        solventAmtList <- rbind(solventAmtList, nex_item)
        colnames(solventAmtList)[3] <- 'Amt' 
      }
    }else if(toString(cmdList$SourceWare[i])==toString(deckMap$position[deckMap$item=='Stock'])){
      
      nexAmt <- length(strsplit(toString(cmdList$TargetWell[i]), split=', ')[[1]])*as.numeric(toString(cmdList$Amt[i]))
      
      #get current source well
      source_well <- toString(cmdList$SourceWell[i])
      #check if current well is already listed in map
      if(source_well %in% stockAmtList[,1]){
        solventAmtList$Amt[solventAmtList[,1]==source_well] <- solventAmtList$Amt[solventAmtList[,1]==source_well]+nexAmt
      }else{
        nex_item <- cbind(stockMap[(stockMap$Slot==source_well),], nexAmt)
        colnames(nex_item)[3] <- 'Amt'
        stockAmtList <- rbind(stockAmtList, nex_item)
        colnames(stockAmtList)[3] <- 'Amt' 
      }
    }
  }
  solventAmtList <- data.frame(solventAmtList)
  stockAmtList <- data.frame(stockAmtList)
  solventAmtList$Amt[solventAmtList$Amt<10000] <- 10000
  stockAmtList$Amt[stockAmtList$Amt<500] <- 500
  solventAmtList$Amt <- ceiling(solventAmtList$Amt/1000)*1000
  stockAmtList$Amt <- ceiling(stockAmtList$Amt/100)*100
  
  ############ OUTPUT #############
  write.table(rbind("", ">Solvent"), "Command List.csv", row.names=F, col.names=F, sep=',', append=F)
  write.table(solventAmtList[,c(2:3)], "Command List.csv", row.names=F, col.names=F, sep=',', append=T)
  write.table(c(">Stock"), "Command List.csv", row.names=F, col.names=F, sep=',', append=T)
  write.table(stockAmtList[,c(2:3)], "Command List.csv", row.names=F, col.names=F, sep=',', append=T)
  write.table(c(">CommandLines"), "Command List.csv", row.names=F, col.names=F, sep=',', append=T)
  write.table(cmdList, "Command List.csv", row.names=F, col.names=F, sep=',', append=T)
  
  return(cmdList)
}

shinyServer(function(input, output, session) {
  #MAIN OPERATION
  #homeDir <- "C:\\Users\\Sebastian\\Documents\\GitHub\\OT2_Controls\\OT2_MIC_v2"
  homeDir <- "/srv/shiny-server/files"
  setwd(homeDir)
  
  #Confirming Upload File
  contents <- reactive({
    infile = input$file
    if(is.null(infile)){return(NULL)}
    
    if(input$do==0){
      dis <- read.xlsx(infile$datapath, 1, rowIndex = c(56:63), colIndex=c(2:13), header=F)
    }else{
      #rename files for safekeeping
      name <- strsplit(input$file$name, '.xl')[[1]][1]
      file.copy(input$file$datapath, getwd())
      file.rename("0.xlsx", paste("UserInput_", name, ".xlsx", sep=""))
      
      #update table view
      dis <- main(input$file$datapath)
      file.rename("Command List.csv", paste("OT2CmdList_", name, ".csv", sep=""))
    }
    return(dis)
  })
  output$tab <- renderTable({contents()})
  
})
