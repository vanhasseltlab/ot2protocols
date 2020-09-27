#EXTERNAL LIBRARIES----------
library(xlsx)

#FUNCTIONS LIBRARY------------
GetStockList <- function(file_name){
  res <- read.xlsx(file_name, 1, 
            endRow=2)
  res <- res[, colSums(is.na(res))==0]
  res <- res[,c(2:length(res))]
  return(res)
}
GetWellVols <- function(file_name){
  res <- read.xlsx(file_name, 1, startRow=5, endRow=6, header=F)
  return(res)
}
GetPlateMap <- function(file_name){
  #read
  res <- read.xlsx(file_name, 1, startRow=56, endRow=64, header=T)
  rownames(res) <- res[,1]
  res <- res[,2:length(res[1,])]
  
  #parse to vector
  map <- c()
  for(row in c(1:8)){
    #subset
    curRow <- unlist(res[row,])
    
    #get info
    well_id <- sapply(c(1:12), function(x) paste(LETTERS[row], toString(x), sep=''))
    curRow <- cbind(well_id, curRow)
    
    #concatenate results
    map <- rbind(map, curRow)
  }
  
  #parse names
  fin_map <- c()
  parsed_names <- sapply(map[,2], function(x) strsplit(x, ' ', fixed=T))
  for(i in c(1:length(parsed_names))){
    if(length(parsed_names[[i]])==4){
      nex_info <- c(paste(parsed_names[[i]][1], parsed_names[[i]][2], parsed_names[[i]][3], sep=' '),
                    parsed_names[[i]][1], #drug name
                    parsed_names[[i]][2], #concentration
                    parsed_names[[i]][3], #solvent
                    parsed_names[[i]][4]) #inoculum
    }else{ #IF FILL
      if(grepl('WATER', map[i,2])){
        nex_info <- c("FILL",
                      "NA",                 #drug name
                      parsed_names[[i]][2], #concentration
                      parsed_names[[i]][3], #solvent
                      "NA")                 #inoculum
      }else{ #IF No inoculum control
        nex_info <- c(paste(parsed_names[[i]][1], parsed_names[[i]][2], parsed_names[[i]][3], sep=' '),
                      parsed_names[[i]][1], #drug name
                      parsed_names[[i]][2], #concentration
                      parsed_names[[i]][3], #solvent
                      "NA")
      }
    }
    fin_map <- rbind(fin_map, nex_info)
    rownames(fin_map) <- c()
  }
  fin_map <- cbind.data.frame(map, fin_map)
  colnames(fin_map) <- c('Well', 'fillID', 'solID', 'DrugType', 'DrugConc', 'Solvent', 'Inoc')
  return(fin_map)
}
CreateSolList <- function(plate_map, total_vol_well, inoc_vol){
  sol_list <- unique(plate_map$solID)
  sol_list <- sol_list[!grepl("FILL", sol_list)]
  
  #get occurence
  occ <- table(plateMap$solID)
  occurences <- as.numeric(occ)
  names(occurences) <- names(occ)
  
  fin_list <- c()
  parsed_names <- sapply(sol_list, function(x) strsplit(x, ' ', fixed=T))
  for(i in c(1:length(parsed_names))){
    nex_info <- c(parsed_names[[i]][1], #drug name
                  parsed_names[[i]][2], #concentration
                  parsed_names[[i]][3], #solvent
                  occurences[sol_list[i]]) #occurence
    
    fin_list <- rbind(fin_list, nex_info)
  }
  rownames(fin_list) <- c()
  
  #final combining
  fin_list <- cbind.data.frame(sol_list, fin_list)
  colnames(fin_list) <- c('SolID', 'DrugType', 'DrugConc', 'Solvent', 'Occurence')
  
  #creating numerics
  fin_list$Occurence <- as.numeric(fin_list$Occurence)
  fin_list$DrugConc <- as.numeric(fin_list$DrugConc)
  
  fin_list <- CalculateDilVolume(fin_list, total_vol_well, inoc_vol)
  return(fin_list)
}
CalculateDilVolume <- function(sol_list, total_vol_well, inoc_vol){
  #calculate initially required amount
  drugSol_well <- total_vol_well - inoc_vol
  solAmt <- sol_list$Occurence * drugSol_well + 150 #adds 10 uL excess
  sol_list <- cbind.data.frame(sol_list, solAmt)
  
  #recalculate amount to adjust with the incoming inoculum
  sol_list$DrugConc <- as.numeric(sol_list$DrugConc) * total_vol_well / (total_vol_well - inoc_vol)
  
  #initiate new list
  new_solList <- c()
  #iterate through all drug and solvent types
  solvents <- unique(sol_list$Solvent)
  drugs <- unique(sol_list$DrugType)
  for(i in c(1:length(solvents))){
    for(j in c(1:length(drugs))){
      #subset
      curList <- subset(sol_list, DrugType==drugs[j] & Solvent==solvents[i])
      
      #perform following actions only if not null
      if(length(curList[,1])>0){
        #order according to concentration
        curList <- curList[order(as.numeric(curList$DrugConc)),]
        
        #check amount needed from above
        needed_from_above <- c()
        for(m in c(1:length(curList[,1]))){
          if(m<length(curList[,1])){ 
            #usual dilution from pre-diluted stock
            amt_needed <- curList$solAmt[m]*curList$DrugConc[m]/curList$DrugConc[m+1]
            
            #check if it is lower than the minimum pipette volume
            if(amt_needed < 30 & amt_needed > 0){
              #set amount from above to 30
              amt_needed <<- amt_needed
              dis <<- curList[m,]
              curList$solAmt[m] <- curList$solAmt[m] * 30 / amt_needed
              amt_needed <- 30 
            }
            
            needed_from_above <- c(needed_from_above, amt_needed)
      
            #add amount to higher concentration
            curList$solAmt[m+1] <- curList$solAmt[m+1] + amt_needed
            
            
          }else{
            #calculate amount for initial dilution
            amt_needed <- curList$solAmt[m]*curList$DrugConc[m]/stockList[curList$DrugType[m]]
            #check if it is lower than the minimum pipette volume
            if(amt_needed < 30 & amt_needed > 0){
              #set amount from above to 30
              amt_needed <<- amt_needed
              dis <<- curList[m,]
              curList$solAmt[m] <- curList$solAmt[m] * 30 / amt_needed
              amt_needed <- 30 
            }
            
            #initial dilution
            needed_from_above <- c(needed_from_above, amt_needed)
          }
        }
      }
      
      nexItem <- cbind(curList, unlist(needed_from_above))
      new_solList <- rbind(new_solList, nexItem)
    }
  }
  
  #renaming
  colnames(new_solList)[length(new_solList[1,])] <- 'AmtHi'
  
  #calculate required solvent amount
  solventAmt <- as.numeric(new_solList$solAmt) - as.numeric(new_solList$AmtHi)
  
  #check required tube size
  reqTube <- replicate(length(new_solList[,1]), "1.5_Eppendorf")
  reqTube[new_solList$solAmt > 1.2*1000] <- "15_Falcon"
  reqTube[new_solList$solAmt > 13*1000] <- "50_Falcon"
  
  #concatenate Info
  new_solList <- cbind.data.frame(new_solList, solventAmt, reqTube)
  
  return(new_solList)
}
CreateDilMap <- function(sol_list, deckMap){
  #initiate map
  small_coords <- c(1, 1)
  large_coords <- c(1, 1)
  spare_coords <- c(1, 1)
  
  small_dilMap <- cbind()
  large_dilMap <- cbind()
  spare_dilMap <- cbind()
  
  #iterate through all items in solution list
  for(i in c(1:length(sol_list[,1]))){
    if(sol_list$solAmt[i]<=1200){
      #check if the main rack is full
      if(small_coords[1]<5){
        #if the main rack is still available
        nexItem <- c(paste(LETTERS[small_coords[1]], toString(small_coords[2]), sep=''),
                     sol_list$SolID[i])
        #concatenate item
        small_dilMap <- rbind(small_dilMap, nexItem)
        
        #update coordinates
        small_coords[2] <- small_coords[2]+1
        if(small_coords[2]>6){
          small_coords[2] <- 1
          small_coords[1] <- small_coords[1] + 1
        }
      }else{
        #if the main rack is full
        nexItem <- c(paste(LETTERS[spare_coords[1]], toString(spare_coords[2]), sep=''),
                     sol_list$SolID[i])
        
        #concatenate item
        spare_dilMap <- rbind(spare_dilMap, nexItem)
        
        #update coordinates
        spare_coords[2] <- spare_coords[2]+1
        if(spare_coords[2]>6){
          spare_coords[2] <- 1
          spare_coords[1] <- spare_coords[1] + 1
        }
      }
    }else{
      #if solution amount is larger
      if(large_coords[1]<4){
        #if the main rack is available
        nexItem <- c(paste(LETTERS[large_coords[1]], toString(large_coords[2]), sep=''),
                     sol_list$SolID[i])
        
        #concatenate item
        large_dilMap <- rbind(large_dilMap, nexItem)
        
        #update coordinates
        large_coords[2] <- large_coords[2]+1
        if(large_coords[2]>5){
          large_coords[2] <- 1
          large_coords[1] <- large_coords[1] + 1
        }
      }else{
        #if the main rack is full
        nexItem <- c(paste(LETTERS[spare_coords[1]], toString(spare_coords[2]), sep=''),
                     sol_list$SolID[i])
        
        #concatenate item
        spare_dilMap <- rbind(spare_dilMap, nexItem)
        
        #update coordinates
        spare_coords[2] <- spare_coords[2]+1
        if(spare_coords[2]>6){
          spare_coords[2] <- 1
          spare_coords[1] <- spare_coords[1] + 1
        }
      }
    }
  }
  
  #get labware locations
  small_dilMap <- cbind(small_dilMap, replicate(length(small_dilMap[,1]),
                                                names(deckMap)[match("1.5_Eppendorf", deckMap)]))
  spare_dilMap <- cbind(spare_dilMap, replicate(length(spare_dilMap[,1]),
                                                names(deckMap)[match("15_Falcon_spare", deckMap)]))
  large_dilMap <- cbind(large_dilMap, replicate(length(large_dilMap[,1]),
                                               names(deckMap)[match("15_Falcon_main", deckMap)]))
  
  #check if racks are empty
  if(length(small_dilMap[,1])==0){
    small_dilMap <- cbind(c('m'), c('m'), c('m'))
  }
  if(length(large_dilMap[,1])==0){
    large_dilMap <- cbind(c('m'), c('m'), c('m'))
  }
  if(length(spare_dilMap[,1])==0){
    spare_dilMap <- cbind(c('m'), c('m'), c('m'))
  }
  
  #re-assign names
  dil_map <- rbind(small_dilMap, large_dilMap, spare_dilMap)
  colnames(dil_map) <- c('Slot', 'Fill', 'Labware')
  dil_map <- data.frame(dil_map)
  dil_map <- dil_map[(dil_map$Slot != 'm'),]
  return(dil_map)
}

#commands
Cmd_InitDist <- function(deck_map, sol_list, solvent_map, dil_map){
  #INITIATE COMMAND LIST
  cmd_list <- c()
  tipID <- 1
  
  for(i in c(1:length(sol_list[,1]))){
    #create next command list
    nexCmd <<- c(names(deck_map)[match('Solvent', deck_map)], #source ware = solvent rack
                solvent_map[solvent_map[,2]==sol_list$Solvent[i],1], #source slot = solvent tube
                dil_map$Labware[dil_map$Fill==sol_list$SolID[i]], #target ware = dilution rack
                dil_map$Slot[dil_map$Fill==sol_list$SolID[i]], #target slot = dilution tube
                sol_list$solventAmt[i], #amount transferred
                "0", tipID, 'Initial solvent distribution')
    dis <<- sol_list$SolventAmt[i]
    #update tip id
    tipID <- tipID + 1
    
    #concatenate command
    cmd_list <- rbind(cmd_list, nexCmd)
  }
  return(cmd_list)
}
Cmd_HiDrug <- function(cmd_list, sol_list, stock_map, deck_map, dil_map){
  #get tip id
  tipID <- as.numeric(cmd_list[length(cmd_list[,1]),7])
  
  #iterate through all drug types
  for(i in c(1:length(stockMap[,1]))){
    #subset
    crSolList <- subset(sol_list, DrugType==stockMap[i,2])
    #iterate through all solvent types
    solvents <- unique(crSolList$Solvent)
    for(j in c(1:length(solvents))){
      #subset highest concentration only
      curSolList <- subset(crSolList, Solvent==solvents[j])
      curSolList <- subset(curSolList, DrugConc==max(curSolList$DrugConc))
      
      #create command
      nexCmd <<- c(names(deck_map)[match("Stock", deck_map)], #source ware = stock rack
                  stock_map[stock_map[,2]==curSolList$DrugType[1],1], #source slot = solvent tube
                  dil_map$Labware[dil_map$Fill==curSolList$SolID[1]], #target ware = dilution rack
                  dil_map$Slot[dil_map$Fill==curSolList$SolID[1]], #target slot = dilution tube
                  curSolList$AmtHi,
                  curSolList$AmtHi, tipID, 'Initial stock dilution')
      
      #update tip ID
      tipID <- tipID + 1
      
      #concatenate command
      cmd_list <- rbind(cmd_list, nexCmd)
    }
  }
  return(cmd_list)
}
Cmd_SerialDil <- function(cmd_list, sol_list, dil_map){
  tipID <- max(as.numeric(cmd_list[,7]), na.rm=T) + 1
  #iterate through all drug and solvents
  drugs <- unique(sol_list$DrugType)
  solvents <- unique(sol_list$Solvent)
  for(i in c(1:length(drugs))){
    for(j in c(1:length(solvents))){
      #subset current list
      cur_SolList <- subset(sol_list, DrugType==drugs[i] & Solvent==solvents[j])
      #check if not null
      if(length(cur_SolList[,1])>0){
        #order list in descending concentration
        cur_SolList <- cur_SolList[order(as.numeric(cur_SolList$DrugConc), decreasing=T), ]
        
        #iterate through the list
        for(m in c(1:(length(cur_SolList[,1])-1))){
          nexCmd <- c(dil_map$Labware[dil_map$Fill==cur_SolList$SolID[m]], #select source labware
                       dil_map$Slot[dil_map$Fill==cur_SolList$SolID[m]], #select source slot
                       dil_map$Labware[dil_map$Fill==cur_SolList$SolID[m+1]],#select target labware
                       dil_map$Slot[dil_map$Fill==cur_SolList$SolID[m+1]],  #select target slot
                       cur_SolList$AmtHi[m+1], #amount to transfer
                       cur_SolList$AmtHi[m+1], tipID, paste('Serially diluting to ', cur_SolList$SolID[m+1]))
          #concatenate command
          if(cur_SolList$AmtHi[m+1]!=0){
            cmd_list <- rbind(cmd_list, nexCmd)
            #update tip ID
            tipID <- tipID+1
          }
        }
      }
    }
  }
  return(cmd_list)
}
Cmd_DrugSolDist <- function(cmd_list, dil_map, plate_map, deck_map, well_info){
  tipID <- max(as.numeric(cmd_list[,7]), na.rm=T)
  transV <- well_info[1,2] - well_info[2,2]
  #iterate through all items in the dilution map
  for(i in c(1:length(dil_map[,1]))){
    target_wells <- plate_map$Well[plate_map$solID==dil_map$Fill[i]]
    nexCmd <- c(dil_map$Labware[i],
                dil_map$Slot[i],
                names(deck_map)[match('96-well', deck_map)],
                paste(target_wells, collapse=', '),
                transV, 0, tipID, paste('Distributing ', dil_map$Fill[i]))
    #concatenate result
    cmd_list <- rbind(cmd_list, nexCmd)
    #update tip ID
    tipID <- tipID + 1
  }
  return(cmd_list)
}
Cmd_Inoculate <- function(plate_map, inoc_map, well_info, cmd_list, deck_map){
  tipID <- max(as.numeric(cmd_list[,7]), na.rm=T) + 1
  #filling blanks with medium
  specialCases <- subset(plate_map, solID != 'FILL' & Inoc == 'NA')
  #iterate through all filled with identical solvent
  solTypes <- unique(specialCases$Solvent)
  for(i in c(1:length(solTypes))){
    curSpec_Case <- subset(specialCases, Solvent==solTypes[i])
    target_wells <- curSpec_Case$Well
    #creating command line
    nexCmd <- c(names(deck_map)[match('Solvent', deck_map)],
                solventMap[solventMap[,2]==curSpec_Case$Solvent[1],1],
                names(deck_map)[match('96-well', deck_map)],
                paste(target_wells, collapse=', '),
                well_info[2,2], 0, tipID, 'Adding blank medium')
    #concatenate command
    cmd_list <- rbind(cmd_list, nexCmd)
    #update tip ID
    tipID <- tipID + 1
  }
  
  #iterate through all inoculum types
  for(i in c(1:length(inoc_map[,1]))){
    #get target wells
    target_wells <- plate_map$Well[plate_map$Inoc==inoc_map[i,2]]
    #create command list
    nexCmd <- c(names(deck_map)[match('Inno', deck_map)],
                inoc_map[i,1],
                names(deck_map)[match('96-well', deck_map)],
                paste(target_wells, collapse=', '),
                well_info[2,2], well_info[2,2], tipID, paste('Inoculating: ', inoc_map[i,2]))
    #concatenate command
    cmd_list <- rbind(cmd_list, nexCmd)
    #update tip ID
    tipID <- tipID + 1
  }
  
  return(cmd_list)
}
Cmd_FillOuter <- function(plate_map, deck_map, solvent_map, well_info, cmd_list){
  tipID <- max(as.numeric(cmd_list[,7]), na.rm=T) + 1
  solvent_map[(solvent_map[,2]=='water' | solvent_map[,2]=='Water'),2] <- 'WATER'
  #select targets
  target_wells <- paste(plate_map$Well[plate_map$solID=='FILL'], collapse=', ')
  #creating command list
  nexCmd <- c(names(deck_map)[match('Solvent', deck_map)],
              solvent_map[solvent_map[,2]=="WATER",1],
              names(deck_map)[match('96-well', deck_map)],
              target_wells,
              well_info[1,2], 0, tipID, 'Filling outer wells with WATER')
  #concatenate result
  cmd_list <- rbind(cmd_list, nexCmd)
  return(cmd_list)
}

#counters
Cal_SolAmt <- function(deck_map, solvent_map, cmd_list){
  req_amt <- c()
  rack_position <- names(deck_map)[match('Solvent', deck_map)]
  for(i in c(1:length(solvent_map[,1]))){
    relCmdList <- subset(cmd_list, SourceLabware==rack_position & SourceSlot==solvent_map[i,1])
    
    #get number of target wells per-row
    n_well <- unlist(sapply(relCmdList$TargetSlot, function(x) length(strsplit(x, split=', ', fixed=T)[[1]])))
    
    #calculate required amount
    solvent_amt <- sum(as.numeric(relCmdList$TransAmt) * n_well)
    
    #concatenate amount
    req_amt <- c(req_amt, solvent_amt)
  }
  
  #put excess of 4 mL
  req_amt <- req_amt + 4000
  #place minimum of 10 mL
  req_amt[req_amt<10000] <- 10000
  #translate to mL
  req_amt <- req_amt/1000
  #round up
  req_amt <- ceiling(req_amt)
  
  #concatenate result
  solvent_map <- cbind(solvent_map, req_amt)
  
  #naming solvent map
  solvent_map <- data.frame(solvent_map)
  names(solvent_map) <- c('Slot', 'Name', 'RequiredAmount')
  
  #Place labware information
  Labware <- replicate(length(solvent_map[,1]), names(deck_map)[match('Solvent', deck_map)])
  Unit <- replicate(length(solvent_map[,1]), "mL")
  Type <- replicate(length(solvent_map[,1]), "50 mL Falcon Tube")
  Category <- replicate(length(solvent_map[,1]), "SOLVENT")
  
  #integrate result
  solvent_map <- cbind.data.frame(Category, Labware, Type, solvent_map, Unit)
  rownames(solvent_map) <- c()
  return(solvent_map)
}
Cal_StockAmt <- function(sol_list, stock_list, stock_map, deck_map){
  #initiate amount list
  amt_list <- replicate(length(stock_list), 0)
  
  #iterate
  drugs <- unique(sol_list$DrugType)
  solvents <- unique(sol_list$Solvent)
  #iterate through all possible combinations of the two
  for(i in c(1:length(drugs))){
    for(j in c(1:length(solvents))){
      #subset
      curList <- subset(sol_list, DrugType==drugs[i] & Solvent==solvents[j])
      #perform if not null
      if(length(curList)>0){
        #get required amount
        cur_reqAmt <- curList$AmtHi[curList$DrugConc==max(curList$DrugConc)]
        
        #add to amount list
        amt_list[names(stock_list)==curList$DrugType[1]] <- amt_list[names(stock_list)==curList$DrugType[1]] + cur_reqAmt
      }
    }
  }
  amt_list <- amt_list
  #make output
  stock_list <- cbind.data.frame(names(stock_list), unlist(stock_list), amt_list)
  colnames(stock_list) <- c('Name', 'Conc', 'RequiredAmount')
  rownames(stock_list) <- c()
  
  #add excess
  stock_list$RequiredAmount <- stock_list$RequiredAmount + 150
  #place minimum
  stock_list$RequiredAmount[stock_list$RequiredAmount<300] <- 300
  #round up
  stock_list$RequiredAmount <- ceiling(stock_list$RequiredAmount/100)*100
  
  #get well locations
  well_loc <- stock_map[stock_map[,2]==stock_list$Name,1]
  stock_map <- cbind.data.frame(well_loc, stock_list)
  colnames(stock_map)[1] <- 'Slot'
  stock_map <- stock_map[,c(1, 2, 4)]
  
  #additional informations
  Labware <- replicate(length(stock_map[,1]), names(deck_map)[match('Stock', deck_map)])
  Unit <- replicate(length(stock_map[,1]), "uL")
  Type <- replicate(length(stock_map[,1]), "1.5 mL Eppendorf Tube")
  Category <- replicate(length(stock_map[,1]), "DRUG STOCK")
  
  #integrate results
  stock_map <- cbind.data.frame(Category, Labware, Type, stock_map, Unit)
  rownames(stock_map) <- c()
  return(stock_map)
}
Cal_InocAmt <- function(inoc_map, cmd_list, deck_map){
  #initiate volume list
  volList <- c()
  #get inoculum position
  inoc_position <- names(deck_map)[match('Inno', deckMap)]
  #iterate through all inoculum types
  for(i in c(1:length(inoc_map[,1]))){
    #get target wells
    selectedRows <- subset(cmd_list, SourceLabware==inoc_position & SourceSlot==inoc_map[i,1])
    #count number of target wells
    n_wells <- unlist(sapply(selectedRows$TargetSlot, function(x) length(strsplit(x, ', ', fixed=T)[[1]])))
    
    #calculate required amount
    req_amt <- sum(n_wells * as.numeric(selectedRows$TransAmt))
    
    #concatenate to list
    volList <- c(volList, req_amt)
  }
  #integrate matrix
  inoc_map <- data.frame(cbind(inoc_map, volList))
  inoc_map[,3] <- as.numeric(inoc_map[,3])
  #get names
  colnames(inoc_map) <- c('Slot', 'Name', 'RequiredAmount')
  
  #add excess 3 mL
  inoc_map$RequiredAmount <- inoc_map$RequiredAmount + 3000
  #place minimum amount
  inoc_map$RequiredAmount[inoc_map$RequiredAmount<5000] <- 5000
  #round up
  inoc_map$RequiredAmount <- ceiling(inoc_map$RequiredAmount/1000) #convert to mL
  
  #additional informations
  Labware <- replicate(length(inoc_map[,1]), names(deck_map)[match('Stock', deck_map)])
  Unit <- replicate(length(inoc_map[,1]), "mL")
  Type <- replicate(length(inoc_map[,1]), "15 mL Falcon Tube")
  Category <- replicate(length(inoc_map[,1]), "INOCULUM")
  
  #integrate results
  inoc_map <- cbind.data.frame(Category, Labware, Type, inoc_map, Unit)
  rownames(inoc_map) <- c()
  return(inoc_map)
}
Cal_DilTubes <- function(dil_map){
  #get number of occurence
  occs <- table(dil_map$Labware)
  
  outputMap <- cbind(names(occs), occs, replicate(length(occs), '15_Falcon'))
  outputMap[outputMap[,1]=='labware_5',3] <- '1.5_Eppendorf'
  outputMap <- data.frame(outputMap)
  colnames(outputMap) <- c('Labware', 'RequiredAmount', 'Name')
  rownames(outputMap) <- c()
  
  #additional informations
  Category <- replicate(length(occs), 'EMPTY TUBES FOR DILUTION')
  Type <- replicate(length(occs), '-')
  Slot <- replicate(length(occs), '-')
  Unit <- replicate(length(occs), 'tubes')
  
  #integrate
  outputMap <- cbind.data.frame(Category, outputMap$Labware, Type, Slot, outputMap$Name,
                                outputMap$RequiredAmount, Unit)
  colnames(outputMap) <- c('Category', 'Labware', 'Type', 'Slot', 'Name', 'RequiredAmount', 'Unit')
  return(outputMap)
}
Cal_DeckAdjustment <- function(cmd_list, deck_map, dil_tubes){
  deck <- matrix(as.character(c(12:1)), ncol=3, byrow=T)
  deck_map <- matrix(deck_map, ncol=3, byrow=T)
  
  fin_deck <- c()
  for(i in c(1:length(deck_map[,1]))){
    fin_deck <- rbind(fin_deck, deck[i,], deck_map[(length(deck_map[,1])-i+1),])
  }
  
  needed <- max(as.numeric(cmd_list$TipID), na.rm=T)
  nbox <- ceiling(needed/96)
  if(nbox<3){
    fin_deck[2,1] <- '(empty)'
    if(nbox<2){
      fin_deck[4,1] <- '(empty)'
    }
  }
  
  #check if tube racks required
  if(!('labware_5' %in% dil_tubes$Location)){
    fin_deck[6,2] <- '(empty)'
  }
  if(!('labware_6' %in% dil_tubes$Location)){
    fin_deck[6,1] <- '(empty)'
  }
  if(!('labware_8' %in% dil_tubes$Location)){
    fin_deck[4,2] <- '(empty)'
  }
  return(fin_deck)
}
#BULK FUNCTION--------
Int_CreateCmdList <- function(deck_map, sol_list, solvent_map, inoc_map,
                              dil_map, stock_map, well_info, plate_map){
  #1. DISTRIBUTE SOLVENT
  cmdlist <- Cmd_InitDist(deck_map, sol_list, solvent_map, dil_map)
  #2. DISTRIBUTE HIGHEST DRUG CONCENTRATION
  cmdlist <- Cmd_HiDrug(cmdlist, sol_list, stock_map, deck_map, dil_map)
  #3. SERIAL DILUTION
  cmdlist <- Cmd_SerialDil(cmdlist, sol_list, dil_map)
  #4. DISTRIBUTING DRUG SOLUTION
  cmdlist <- Cmd_DrugSolDist(cmdlist, dil_map, plate_map, deck_map, well_info)
  #5. FILLING OUTER WELLS
  cmdlist <- Cmd_FillOuter(plate_map, deck_map, solvent_map, well_info, cmdlist)
  #6. INOCULATION
  cmdlist <- Cmd_Inoculate(plate_map, inoc_map, well_info, cmdlist, deck_map)
  
  #NAMING
  cmdlist <- data.frame(cmdlist)
  colnames(cmdlist) <- c('SourceLabware', 'SourceSlot', 'TargetLabware',
                         'TargetSlot', 'TransAmt', 'MixAmt', 'TipID', 'Comment')
  return(cmdlist)
}
#TEST AREA-----------
#file
homeDir <- "C:\\Users\\Sebastian\\Desktop\\MSc Leiden 2nd Year\\##LabAst Works\\FinVersion\\OT2_UpstreamExcelProc_MAIN"
inptFile <- "PlateMapInpt.xlsx"
filename <- paste(homeDir, "\\", inptFile, sep='')

#test
#READ PLATE------
stockList <- GetStockList(filename)
wellInfo <- GetWellVols(filename)
plateMap <- GetPlateMap(filename)

#GET SOLUTION LIST AND DILUTION SCHEME-----------
solList <- CreateSolList(plateMap, wellInfo[1,2], wellInfo[2,2])

#################
#COMMANDS
#################
#-1. LOADING DECK MAP-----------
deckMap <- c('tip', '96-well', 'Inno',
             '15_Falcon_main', '1.5_Eppendorf', 'Solvent',
             'tip', '15_Falcon_spare', 'Stock',
             'tip', '(empty)', 'TRASH')
names(deckMap) <- sapply(c(1:12), function(x) paste('labware', toString(x), sep='_'))

# 0. LOAD LABWARES--------
#solvents
  #initiate map
coords <- c(1, 1)
solventMap <- c()
#iterate through all solvents in platemap
solvents <- c(unique(plateMap$Solvent))
for(i in c(1:length(solvents))){
  nexItem <- c(paste(LETTERS[coords[1]], toString(coords[2]), sep=''), solvents[i])
  #place to map
  solventMap <- rbind(solventMap, nexItem)
  
  #update fill coordinates
  coords[2] <- coords[2]+1
  if(coords[2]>3){
    coords[2] <- 1
    coords[1] <- coords[1] + 1
  }
}

#stock
 #initiate map
coords <- c(1, 1)
stockMap <- c()
#iterate through all items in stockList
for(i in c(1:length(stockList))){
  nexItem <- c(paste(LETTERS[coords[1]], toString(coords[2]), sep=''), names(stockList)[i])
  
  #place to map
  stockMap <- rbind(stockMap, nexItem)
  
  #update fill coordinates
  coords[2] <- coords[2]+1
  if(coords[2]>6){
    coords[2] <- 1
    coords[1] <- coords[1] + 1
  }
}

#assign slots for diluted solutions
dilMap <- CreateDilMap(solList, deckMap)

#inoculum
  #initiate map
coords <- c(1, 1)
inocMap <- c() #assume rack is 15 mL Falcon tube rack
#iterate through all inoculum types in plate
inocs <- unique(plateMap$Inoc)
inocs <- inocs[inocs!='NA']
for(i in c(1:length(inocs))){
  nexItem <- c(paste(LETTERS[coords[1]], toString(coords[2]), sep=''), inocs[i])
  inocMap <- rbind(inocMap, nexItem)
  #update coordinate
  coords[2] <- coords[2]+1
  if(coords[2]>5){
    coords[2] <- 1
    coords[1] <- coords[1] + 1
  }
}

# 1. CREATE COMMAND LIST
cmdList <- Int_CreateCmdList(deckMap, solList, solventMap, inocMap,
                              dilMap, stockMap, wellInfo, plateMap)

# 2. BUNDLING OUTPUT-------
allAmt <- rbind.data.frame(Cal_SolAmt(deckMap, solventMap, cmdList),
                           Cal_StockAmt(solList, stockList, stockMap, deckMap),
                           Cal_InocAmt(inocMap, cmdList, deckMap))

# 3. CALCULATE REQUIRED NUMBER OF ITEMS--------
dilTubes <- Cal_DilTubes(dilMap)

# 4. DECK LAYOUT FOR USER---------
finDeck <- Cal_DeckAdjustment(cmdList, deckMap, dilTubes)

#MAIN FUNCTION---------
main <- function(file_name){
  #READ PLATE------
  stockList <- GetStockList(file_name)
  wellInfo <- GetWellVols(file_name)
  plateMap <- GetPlateMap(file_name)
  
  #GET SOLUTION LIST AND DILUTION SCHEME-----------
  solList <- CreateSolList(plateMap, wellInfo[1,2], wellInfo[2,2])
  
}

