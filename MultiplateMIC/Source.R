#LIBRARIES----------
library(xlsx)
options(stringsAsFactors = F)
#FUNCTIONS--------------
getInput <- function(file_name){
  #read drug stock info
  drugInfo <- t(read.xlsx(file_name, 1, rowIndex=c(1,2,3), header=F, stringsAsFactors=F))
  drugInfo <- drugInfo[!apply(is.na(drugInfo), 1, all),c(1,2)]
  colnames(drugInfo) <- c('DrugName', 'StockConc')
  drugInfo <- drugInfo[c(2:length(drugInfo[,1])),]
  
  if(length(drugInfo)>2){
    drugInfo <- as.data.frame(drugInfo, stringsAsFactors=F)
    drugInfo[,2] <- as.numeric(drugInfo[,2])
  }else{
    drugInfo <- as.data.frame(t(drugInfo), stringsAsFactors=F)
    drugInfo[,2] <- as.numeric(drugInfo[,2])
  }
  
  #read well informations
  wellInfo <- unlist(read.xlsx(file_name, 1, rowIndex=c(5, 6), header=F, colIndex=3, stringsAsFactors=F))
  names(wellInfo) <- c('FinalVol', 'InocVol')
  
  #get number of plates
  nPlate <- as.numeric(read.xlsx(file_name, 1, rowIndex=6, colIndex=10, header=F, stringsAsFactors=F))
  
  #get final plate summary
  plateMap_slot <- c()
  for(i in c(1:12)){
    for(j in c(1:8)){
      plateMap_slot <- c(plateMap_slot, paste(LETTERS[j], toString(i), sep=''))
    }
  }
  plateMap_drugs <- unlist(read.xlsx(file_name, 1, rowIndex=c(10:17), colIndex=c(2:13), header=F, stringsAsFactors=F))
  plateMap_concs <- unlist(read.xlsx(file_name, 1, rowIndex=c(21:28), colIndex=c(2:13), header=F, stringsAsFactors=F))
  plateMap_solvent <- unlist(read.xlsx(file_name, 1, rowIndex=c(33:40), colIndex=c(2:13), header=F, stringsAsFactors=F))
  plateMap_inoc <- unlist(read.xlsx(file_name, 1, rowIndex=c(45:52), colIndex=c(2:13), header=F, stringsAsFactors=F))
  plateMap_solID <- paste(plateMap_drugs, plateMap_solvent, sep='_')
  plateMap_dilConc <- plateMap_concs*wellInfo[1]/(wellInfo[1]-wellInfo[2])
  plateMap <- cbind.data.frame(plateMap_solID, plateMap_slot, plateMap_drugs, plateMap_concs, 
                               plateMap_solvent, plateMap_inoc, plateMap_dilConc, stringsAsFactors=F)
  colnames(plateMap) <- c('solID', 'Slot', 'DrugName', 'Conc', 'Solvent', 'Inoc', 'DilConc')
  
  
  #output
  res <- list(drugInfo, wellInfo, plateMap, nPlate)
  
  return(res)
}
defSolutionsMap <- function(deck_map, plate_map, drug_info){
  #Assigning non-well solution map
  solventRack <- deck_map[grepl('olvent', deck_map[,2]),1]
  stockRack <- deck_map[grepl('tock', deck_map[,2]),1]
  Labware <- c(unlist(replicate(6, solventRack)),
               unlist(replicate(24, stockRack)))
  #get slots
  Slot <- c()
  for(i in c(1,2)){
    for(j in c(1:3)){
      Slot <- c(Slot, paste(LETTERS[i], toString(j), sep=''))
    }
  }
  for(i in c(1:4)){
    for(j in c(1:6)){
      Slot <- c(Slot, paste(LETTERS[i], toString(j), sep=''))
    }
  }
  #create empty fills
  Fill <- replicate(length(Labware), "")
  Conc <- replicate(length(Labware), "")
  Vol <- replicate(length(Labware), "")
  solutionsMap <- cbind.data.frame(Labware, Slot, Fill, Conc, Vol, stringsAsFactors=F)
  
  #place solvents
  solvents <- unique(plate_map$Solvent)
  new_fills <- solutionsMap$Fill[solutionsMap$Labware==toString(solventRack)] 
  new_fills[1:length(solvents)] <- solvents
  solutionsMap$Fill[solutionsMap$Labware==toString(solventRack)] <- new_fills
  
  new_fills[1:length(solvents)] <- 45000
  solutionsMap$Vol[solutionsMap$Labware==toString(solventRack)] <- new_fills
  
  #place stocks
  new_fills <- solutionsMap$Fill[solutionsMap$Labware==toString(stockRack)]
  new_fills[1:length(drug_info$DrugName)] <- drug_info$DrugName
  solutionsMap$Fill[solutionsMap$Labware==toString(stockRack)] <- new_fills
  
  new_fills[1:length(drug_info$DrugName)] <- drug_info$StockConc
  solutionsMap$Conc[solutionsMap$Labware==toString(stockRack)] <- new_fills
  
  new_fills[1:length(drug_info$DrugName)] <- 1200
  solutionsMap$Vol[solutionsMap$Labware==toString(stockRack)] <- new_fills
  
  #make into data frame
  solutionsMap <- data.frame(solutionsMap, stringsAsFactors=F)
  solutionsMap$Vol <- as.numeric(solutionsMap$Vol)
  solutionsMap$Conc <- as.numeric(solutionsMap$Conc)
  
  return(solutionsMap)
}
cmd_PreliminaryDilution <- function(plate_map, deck_map, drug_info, solutions_map){
  solIDs <- unique(subset(plate_map, DrugName!='NA')$solID)
  tipID <- 1
  vol <- 1000 #ul, in dilution tube
  #initiate command list
  cmdList <- c()
  for(i in c(1:length(solIDs))){
    #extracting preliminary info
    cur_drugName <- subset(plate_map, solID==solIDs[i])$DrugName[1]
    cur_stockConc <- subset(drug_info, DrugName==cur_drugName)$StockConc
    cur_targetConc <- max(subset(plate_map, DrugName==cur_drugName)$DilConc)
    labware <- toString(deck_map[grepl('tock', deck_map[,2]),1])
    solventLoc <- c(toString(deck_map[grepl('olvent', deck_map[,2]),1]),
                    solutions_map$Slot[solutions_map$Fill==subset(plate_map, solID==solIDs[i])$Solvent[1]])
    
    #iterate dilution as long as needed
    while(cur_stockConc/cur_targetConc>10){
      #calculate required current required concentration
      current_dil_fac <- min(10, cur_stockConc/(10*cur_targetConc))
      
      #get source location (slot)
      source <- solutions_map$Slot[(solutions_map$Fill==cur_drugName) & (solutions_map$Conc==cur_stockConc)]
      #get target location (slot)
      target <- solutions_map$Slot[(solutions_map$Labware==labware) & (solutions_map$Fill=='')]
      target <- target[1]
      
      #step1: transfer solvent
      solventVol <- vol - vol / current_dil_fac
      nexCmd <- c(solventLoc, labware, target, solventVol, 0, tipID, 'Preliminary Dilution_solvent')
      #concatenate command list, increment tipID
      tipID <- tipID+1
      cmdList <- rbind(cmdList, nexCmd)
      
      #step2: transfer stock
      nexCmd <- c(labware, source, labware, target, 
                  vol / current_dil_fac, vol / current_dil_fac, tipID, 'Preliminary Dilution_drug')
      #concatenate command list, increment tipID
      tipID <- tipID+1
      cmdList <- rbind(cmdList, nexCmd)
      #update solutions list
      emptySlot <- solutions_map$Slot[solutions_map$Labware==labware & solutions_map$Fill==""][1]
      solutions_map$Fill[solutions_map$Labware==labware & solutions_map$Slot==emptySlot] <- cur_drugName
      solutions_map$Conc[solutions_map$Labware==labware & solutions_map$Slot==emptySlot] <- cur_stockConc/current_dil_fac
      solutions_map$Vol[solutions_map$Labware==labware & solutions_map$Slot==emptySlot] <- vol
      
      #update current highest concentration
      cur_stockConc <- cur_stockConc/current_dil_fac
    }
  }
  cmdList <- data.frame(cmdList, stringsAsFactors=F)
  
  if(length(cmdList)==0){
    cmdList <- rbind.data.frame(cmdList, c("", "", "", "", "", "", 0, "skip"))
  }
  
  colnames(cmdList) <- c('SourceLabware', 'SourceSlot', 'TargetLabware', 'TargetSlot',
                         'TransAmt', 'MixAmt', 'tipID', 'Comment')
  cmdList$tipID <- as.numeric(cmdList$tipID)
  return(list(solutions_map, cmdList))
}
cmd_FirstDWDilution <- function(plate_map, deck_map, drug_info, solutions_map, 
                                cmd_list, well_info, serial_step, n_plates){
  #first dilution
  tipID <- max(cmd_list$tipID) + 1
  
  solIDs <- unique(subset(plate_map, DrugName!='NA')$solID)
  finAmt_per_wells <- well_info[1] - well_info[2]
  solvent_labware <- toString(deck_map[grepl('olvent', deck_map[,2]),1])
  stock_labware <- toString(deck_map[grepl('tock', deck_map[,2]),1])
  dwell_labware <- toString(deck_map[grepl('eep', deck_map[,2]),1])
  
  for(i in c(1:length(solIDs))){
    #extract current informations
    cur_drugSet <- subset(plate_map, solID==solIDs[i])
    cur_dilConc <- max(cur_drugSet$DilConc)
    target_wells <- cur_drugSet$Slot[cur_drugSet$DilConc==cur_dilConc]
    
    cur_solMap <- subset(solutions_map, Fill==cur_drugSet$DrugName[1])
    cur_stockConc <- cur_solMap$Conc[cur_solMap$Conc == min(cur_solMap$Conc)]
    
    solvent_slot <- solutions_map$Slot[solutions_map$Fill==cur_drugSet$Solvent[1]]
    
    stock_slot <- cur_solMap$Slot[cur_solMap$Conc == min(cur_solMap$Conc)]
    
    #calculate required volumes
    reqVolume <- (serial_step/(serial_step-1))*((n_plates*finAmt_per_wells)+300)
    reqVolume <- max(reqVolume, 500)
    stockVolume <- reqVolume*cur_dilConc/cur_stockConc
    solventVolume <- reqVolume - stockVolume
    
    #distribute solvents
    nexCmd <- c(solvent_labware, solvent_slot, dwell_labware, paste(target_wells, collapse=', '),
                solventVolume, 0, tipID, 'solvent distribution for first stock dilution in deep wells')
    
    #concatenate command, increment tip ID
    cmd_list <- rbind.data.frame(cmd_list, nexCmd, stringsAsFactors=F)
    tipID <- tipID + 1
    
    #distribute stock solutions; per-well to change tips
    for(j in c(1:length(target_wells))){
      nexCmd <- c(stock_labware, stock_slot, dwell_labware, target_wells[j],
                  stockVolume, stockVolume, tipID, 'First dilution to deep wells')
      #concatenate command, increment tip ID
      cmd_list <- rbind.data.frame(cmd_list, nexCmd, stringsAsFactors=F)
      tipID <- tipID + 1
    }
  }
  
  return(cmd_list)
}
cmd_MainSerialDilution <- function(cmd_list, plate_map, deck_map, solutions_map, well_info, n_plates, col_rat){
  col_rat <<- col_rat
  tipID <- max(as.numeric(cmd_list$tipID)) + 1
  finAmt_per_wells <- well_info[1] - well_info[2]
  
  #calculate amount of solvents needed
  reqVolume <- (serial_step/(serial_step-1))*((n_plates*finAmt_per_wells)+300)
  reqVolume <- max(reqVolume, 500)
  solventVolume <- reqVolume-reqVolume/serial_step
  
  #distributing solvent
  solvents <- unique(plate_map$Solvent)
  for(i in c(1:length(solvents))){
    #obtain target slots
    source_labware <- toString(deck_map[grepl("olvent",deck_map[,2]),1])
    source_slot <- solutions_map$Slot[solutions_map$Fill==solvents[i]]
    target_labware <- toString(deck_map[grepl("eep",deck_map[,2]),1])
    target_slots <- plate_map$Slot[plate_map$Solvent==solvents[i] & 
                                     !(grepl("2", plate_map$Slot)) & 
                                     !(is.na(plate_map$DrugName))]
    
    #give commands
    nexCmd <- c(source_labware, source_slot, target_labware, paste(target_slots, collapse=', '),
                solventVolume, 0, tipID, 'solvent distribution for serial dilution')
    
    #concatenate commands; increment tipID
    cmd_list <- rbind.data.frame(cmd_list, nexCmd, stringsAsFactors=F)
    tipID <- tipID + 1
    
    #add solvent to no drug controls
    target_slots <- plate_map$Slot[plate_map$Solvent==solvents[i] & 
                                     !(grepl("2", plate_map$Slot)) & 
                                     !(is.na(plate_map$DrugName))&
                                     plate_map$Conc==0]
    #give commands
    nexCmd <- c(source_labware, source_slot, target_labware, paste(target_slots, collapse=', '),
                reqVolume/serial_step, 0, tipID, 'solvent distribution for serial dilution')
    
    #concatenate commands; increment tipID
    cmd_list <- rbind.data.frame(cmd_list, nexCmd, stringsAsFactors=F)
    tipID <- tipID + 1
  }
  #serially diluting
  source_labware <- toString(deck_map[grepl("eep",deck_map[,2]),1])
  for(i in c(1:length(col_rat[1,]))){
    if(!Inf %in% col_rat[,i] | TRUE %in% is.na(col_rat[,i])){
      #if zeros not included
      source_slots <- sapply(c(1:8), function(x) paste(LETTERS[x], toString(i+1), sep=''))
      target_slots <- sapply(c(1:8), function(x) paste(LETTERS[x], toString(i+2), sep=''))
      
      #create command
      nexCmd <- c(source_labware, paste(source_slots, collapse=', '), 
                  source_labware, paste(target_slots, collapse=', '),
                  reqVolume/serial_step, reqVolume/serial_step, tipID, 'Serial dilution with MC pipette')
      
      #concatenate command; increment tipID
      cmd_list <- rbind.data.frame(cmd_list, nexCmd, stringsAsFactors=F)
      tipID <- tipID + 8
    }else{
      #if zeros included
      for(j in c(2:7)){
        #manually iterate through all items
        #if current well is not the zero-concentration well
        if(col_rat[j-1,i]!=Inf & col_rat[j-1,i]!=-Inf){
          #get current slots
          source_slots <- paste(LETTERS[j-1], toString(i+1), sep='')
          target_slots <- paste(LETTERS[j-1], toString(i+2), sep='')
          #create command
          nexCmd <- c(source_labware, paste(source_slots, collapse=', '), 
                      source_labware, paste(target_slots, collapse=', '),
                      reqVolume/serial_step, reqVolume/serial_step, tipID, 'Serial dilution with SC pipette')
          
          #concatenate command; increment tipID
          cmd_list <- rbind.data.frame(cmd_list, nexCmd, stringsAsFactors=F)
          tipID <- tipID + 1
        }
      }
    }
  }
  return(cmd_list)
}
cmd_DistributeCol <- function(column_num, targeted_locs, deck_map, cmd_list, finAmt_per_wells){
  #get general information
  tipID <- max(as.numeric(cmd_list$tipID)) + 1
  
  #get source labware (the 96 deep well plate)
  source_labware <- toString(deck_map[grepl('eep',deck_map[,2]),1])
  
  #get slots
  slot <- paste(sapply(c(1:8), function(x) paste(LETTERS[x], toString(column_num), sep='')), collapse=', ')
  
  #iterate through all targeted locations
  for(j in c(1:length(targeted_locs))){
    #create command
    nexCmd <- c(source_labware, slot, targeted_locs[j], slot,
                finAmt_per_wells, finAmt_per_wells, tipID, paste('Distributing column', toString(column_num), 
                                                                 'to', targeted_locs[j]))
    
    #concatenate command; no need to increment tip id (same solution)
    cmd_list <- rbind.data.frame(cmd_list, nexCmd, stringsAsFactors=F)
  }
  return(cmd_list)
}
cmd_Distribute_all <- function(deck_map, n_plate, plate_map, cmd_list, well_info){
  #get location of targeted labwares
  selectedWellLoc <- unlist(deck_map[grepl('96-well', deck_map[,2]),2])
  plate_select <- sapply(selectedWellLoc, 
                         function(x) strsplit(x, '_')[[1]][2] %in% as.character(c(1:n_plate)))
  selectedWellLoc <- unlist(deck_map[deck_map[,2] %in% selectedWellLoc[plate_select],1])
  
  #iterate through all columns in plate
  plateMatrix <- matrix(plate_map$solID, nrow=8)
  for(i in c(1:12)){
    #check if current column is all empty
    if(length(unique(plateMatrix[,i]))!=1){
      #if not, operate
      cmd_list <- cmd_DistributeCol(i, selectedWellLoc, deck_map, cmd_list, well_info[1] - well_info[2])
      cmd_list$tipID[as.numeric(cmd_list$tipID)==max(as.numeric(cmd_list$tipID))] <- max(as.numeric(cmd_list$tipID)) + 7
    }else if(!(grepl("NA",unique(plateMatrix[,i])))){
      #also operate
      cmd_list <- cmd_DistributeCol(i, selectedWellLoc, deck_map, cmd_list, well_info[1] - well_info[2])
      cmd_list$tipID[as.numeric(cmd_list$tipID)==max(as.numeric(cmd_list$tipID))] <- max(as.numeric(cmd_list$tipID)) + 7
    }else{
      #also operate; but put a pre-defined volume
      cmd_list <- cmd_DistributeCol(i, selectedWellLoc, deck_map, cmd_list, well_info[1])
      cmd_list$tipID[as.numeric(cmd_list$tipID)==max(as.numeric(cmd_list$tipID))] <- max(as.numeric(cmd_list$tipID)) + 7
    }
  }
  return(cmd_list)
}
cal_ReqSolutions <- function(solutions_map, plate_map, cmd_list, deck_map){
  #get locations of the input solutions
  filled_locs <- solutions_map[solutions_map$Fill != '',]
  drugs <- unique(plate_map$DrugName[!(is.na(plate_map$DrugName))])
  for(i in c(1:length(drugs))){
    cur_drug <- subset(filled_locs, Fill==drugs[i])
    filled_locs <- rbind.data.frame(filled_locs[filled_locs$Fill!=drugs[i],],
                                    cur_drug[as.numeric(cur_drug$Conc)==max(as.numeric(cur_drug$Conc)),], stringsAsFactors=F)
  }
  #zero the volumes
  filled_locs$Vol <- 0
  
  #iterate through command line
  for(i in c(1:length(cmd_list[,1]))){
    #check if the current line takes from the original input solutions
    if(cmd_list$SourceLabware[i] %in% filled_locs$Labware){
      if(cmd_list$SourceSlot[i] %in% filled_locs$Slot[filled_locs$Labware==cmd_list$SourceLabware[i]]){
        #if true, add the required volume for the current command line
        filled_locs$Vol[filled_locs$Labware==cmd_list$SourceLabware[i] 
                        & filled_locs$Slot==cmd_list$SourceSlot[i]] <- filled_locs$Vol[filled_locs$Labware==cmd_list$SourceLabware[i] 
                                                                                       & filled_locs$Slot==cmd_list$SourceSlot[i]] + 
          as.numeric(cmd_list$TransAmt[i]) * length(strsplit(cmd_list$TargetSlot[i], split=', ')[[1]])
      }
    }
  }
  
  #round up
  filled_locs$Vol[filled_locs$Labware==toString(deck_map[grepl("olvent",deck_map[,2]),1])] <- ceiling(filled_locs$Vol[filled_locs$Labware==toString(deck_map[grepl("olvent",deck_map[,2]),1])]/1000) + 3
  filled_locs$Vol[filled_locs$Labware==toString(deck_map[grepl("tock",deck_map[,2]),1])] <- ceiling(filled_locs$Vol[filled_locs$Labware==toString(deck_map[grepl("tock",deck_map[,2]),1])]/100)*100 + 100
  
  return(filled_locs)
}
cal_amtList_Excess <- function(amt_list, cmd_list, deck_map){
  tubes <- amt_list
  tubes$Vol <- 0
  cmd_list <- cmd_list
  for(i in c(1:length(cmd_list[,1]))){
    #if current tube+slot is sourced
    if(cmd_list$SourceLabware[i] %in% tubes$Labware){
      if(cmd_list$SourceSlot[i] %in% tubes$Slot[tubes$Labware == cmd_list$SourceLabware[i]]){
        #check volume used after transfer
        volUsed_after <- tubes$Vol[tubes$Labware == cmd_list$SourceLabware[i] & tubes$Slot == cmd_list$SourceSlot[i]] + 
          as.numeric(cmd_list$TransAmt[i])*length(strsplit(cmd_list$TargetSlot[i], split=', ')[[1]])
        deltaV <- as.numeric(cmd_list$TransAmt[i])*length(strsplit(cmd_list$TargetSlot[i], split=', ')[[1]])
        
        if((volUsed_after <= 45000 & cmd_list$SourceLabware[i] == 'labware_11') | 
           (volUsed_after <= 1200 & cmd_list$SourceLabware[i] == 'labware_9')){
          #if volume is still acceptable
          #calculate volume
          tubes$Vol[tubes$Labware == cmd_list$SourceLabware[i] & tubes$Slot == cmd_list$SourceSlot[i]] <- tubes$Vol[tubes$Labware == cmd_list$SourceLabware[i] & tubes$Slot == cmd_list$SourceSlot[i]] + 
            as.numeric(cmd_list$TransAmt[i])*length(strsplit(cmd_list$TargetSlot[i], split=', ')[[1]])
        }else{
          #if volume exceeded limit
          #place new tube
          if(cmd_list$SourceLabware[i] == 'labware_9'){
            sol_or_stock <- 'tock'
          }else{
            sol_or_stock <- 'olvent'
          }
          solvent_map <- subset(tubes, Labware==toString(deck_map[grepl(sol_or_stock, deck_map[,2]),1]))
          #get last filled tube
          last_filled <- solvent_map$Slot[length(solvent_map[,1])]
          new_slot <- c(substring(last_filled, 1, 1), substring(last_filled, 2, 2))
          new_slot[2] <- as.numeric(new_slot[2])+1
          if(as.numeric(new_slot[2])>3){
            new_slot[1] <- LETTERS[which(LETTERS==new_slot[1])+1]
            new_slot[2] <- 1
          }
          new_slot <- paste(new_slot, collapse='')
          
          #assign new tube in 'tubes'
          nexDat <- cbind.data.frame(cmd_list$SourceLabware[i], 
                                     new_slot,
                                     tubes$Fill[tubes$Labware==cmd_list$SourceLabware[i] & tubes$Slot==cmd_list$SourceSlot[i]],
                                     tubes$Conc[tubes$Labware==cmd_list$SourceLabware[i] & tubes$Slot==cmd_list$SourceSlot[i]],
                                     deltaV, stringsAsFactors=F)
          colnames(nexDat) <- colnames(tubes)
          tubes <- rbind.data.frame(tubes, nexDat, stringsAsFactors=F)
          #update command lines
          old_slots <- cmd_list[c(i:length(cmd_list[,1])),]
          new_slots <- old_slots
          new_slots$SourceSlot[old_slots$SourceLabware == cmd_list$SourceLabware[i] &
                                 old_slots$SourceSlot == cmd_list$SourceSlot[i]] <- new_slot
          cmd_list$SourceSlot[c(i:length(cmd_list[,1]))] <- new_slots$SourceSlot
        }
      }
    }
  }
  
  #order items
  amt_list <- tubes[order(tubes$Labware),]
  #round up; add excess
  amt_list$Vol[amt_list$Labware==toString(deck_map[grepl('olvent',deck_map[,2]),1])] <- ceiling(amt_list$Vol[amt_list$Labware==toString(deck_map[grepl('olvent',deck_map[,2]),1])]/1000) + 3
  amt_list$Vol[amt_list$Labware==toString(deck_map[grepl('olvent',deck_map[,2]),1])] <- sapply(amt_list$Vol[amt_list$Labware==toString(deck_map[grepl('olvent',deck_map[,2]),1])], function(x) min(x, 49))
  amt_list$Vol[amt_list$Labware==toString(deck_map[grepl('tock',deck_map[,2]),1])] <- ceiling(amt_list$Vol[amt_list$Labware==toString(deck_map[grepl('tock',deck_map[,2]),1])]/100)*100 + 200
  amt_list$Vol[amt_list$Labware==toString(deck_map[grepl('tock',deck_map[,2]),1])] <- sapply(amt_list$Vol[amt_list$Labware==toString(deck_map[grepl('tock',deck_map[,2]),1])], function(x) min(x, 1300))
  
  #return result
  res <- list(cmd_list, amt_list)
  return(res)
}
cal_correctDeck <- function(deck_map, cmd_list, n_plate){
  #get rid of excess tipracks
  n_tipRack <- ceiling(max(as.numeric(cmd_list$tipID)) / 96)
  if(n_tipRack==1){
    deck_map[deck_map[,1]=='labware_10'] <- '(empty)'
  }
  #get rid of excess plates
  reqPlates <- as.character(c(1:n_plate))
  for(i in c(1:length(deck_map[,1]))){
    name <- strsplit(toString(deck_map[i,2]), split='_')[[1]][2]
    if(grepl('96-well', toString(deck_map[i,2])) & !(name %in% reqPlates)){
      deck_map[i,2] <- '(empty)'
    }
  }
  return(deck_map)
}
cal_dilTube <- function(solutions_map, deck_map){
  #checking extra tubes
  extra_tubes <- 0
  drugs <- unique(solutions_map$Fill[solutions_map$Labware==toString(deck_map[grepl('tock',deck_map[,2]),1])])
  drugs <- drugs[drugs!='']
  #iterate through all drug types
  for(i in c(1:length(drugs))){
    cur_solutionsMap <- subset(solutions_map, Fill==drugs[i])
    extra_tubes <- extra_tubes + 
      length(subset(cur_solutionsMap, Conc < max(cur_solutionsMap$Conc))[,1])
  }
  return(extra_tubes)
}
user_AmountList <- function(amt_list, deck_map){
  #get tube types
  tubeTyp <- replicate(length(amt_list[,1]), 'ERROR!')
  tubeTyp[amt_list$Labware==toString(deck_map[grepl('olvent',deck_map[,2]),1])] <- '50 mL Falcon tube'
  tubeTyp[amt_list$Labware==toString(deck_map[grepl('tock',deck_map[,2]),1])] <- '1.5 mL Eppendorf tube'
  
  #get units
  units <- replicate(length(amt_list[,1]), 'uL')
  units[!(grepl('1.5', tubeTyp))] <- 'mL'
  
  
  #remove 'labware'
  amt_list$Labware <- sapply(amt_list$Labware, function(x) substring(x, 9, nchar(x)))
  
  #bind, change names, rearrange
  amt_list <- cbind.data.frame(amt_list, units, tubeTyp, stringsAsFactors=F)
  colnames(amt_list) <- c('Location in OT2', 
                          'Slot',
                          'Solution Name',
                          'Concentration',
                          'Volume in Tube',
                          'Unit',
                          'Item Type' )
  amt_list <- amt_list[,c(1,2,7,3:6)]
  return(amt_list)
}

cmd_FillOuterWells_DeepWell <- function(deck_map, n_plate, cmd_list, 
                                        well_info, solutions_map, plate_map){
  
  #get general information
  tipID <- max(as.numeric(cmd_list$tipID)) + 1
  
  #get slots to fill
  outer_wells <- plate_map[is.na(plate_map$DrugName) & is.na(plate_map$Inoc),]
  
  #get solvent types
  solvents <- unique(outer_wells$Solvent)
  
  #get location of deep well plate
  selectedWellLoc <- toString(deck_map[grepl('eep', deck_map[,2]),1])
  
  #iterate through all solvent types
  for(i in c(1:length(solvents))){
    #get solvent location
    source_labware <- solutions_map$Labware[solutions_map$Fill==solvents[i]]
    source_slot <- solutions_map$Slot[solutions_map$Fill==solvents[i]]
    target_labware <- solutions_map$Labware[solutions_map$Fill==solvents[i]]
    target_slot <- paste(outer_wells$Slot[outer_wells$Solvent==solvents[i]], collapse=', ')
    #get command line
    amt <- well_info[1]*n_plate + 300
    amt <- min(amt, 1800)
    nexCmd <- c(source_labware, source_slot, selectedWellLoc, target_slot,
                amt, 0, tipID, paste('Filling outer wells of plate in', selectedWellLoc))
    
    #concatenate to command list
    cmd_list <- rbind.data.frame(cmd_list, nexCmd, stringsAsFactors=F)
    #increment tip ID
    tipID <- tipID + 1
  }
  return(cmd_list)
  
}
main <- function(file_name){
  error <- c()
  #A. Read Input File
  inputFile <- getInput(file_name)
  drugInfo <- inputFile[[1]]
  wellInfo <- inputFile[[2]]
  plateMap <- inputFile[[3]]
  nPlate <- inputFile[[4]]
  
  #check if dilution factor is equal
  plateMatrix <- matrix(plateMap$DilConc, nrow=8)
  colRat <- c()
  for(col in c(3:11)){
    colRat <- cbind(colRat, plateMatrix[c(2:7),col-1]/plateMatrix[c(2:7),col])
  }
  
  if(length(unique(colRat[colRat!=Inf & colRat != -Inf]))==1){
    serial_step <<- unique(colRat[colRat!=Inf & colRat != -Inf])
  }else{
    serial_step <<- unique(colRat[colRat!=Inf & colRat != -Inf])[1]
    error <- c(error, 'Inequal dilution factor')
  }
  
  #Assigning Deck Map and Slots
  deckMap <- cbind.data.frame(sapply(c(1:12), function(x) paste('labware', x, sep='_')),
                              c('96-well_4', '96-well_5', '96-well_6', '96-well_1',
                                '96-well_2', '96-well_3', 'tip', 'DeepWell', 'Stock_1.5',
                                'tip', 'solvent_50', 'TRASH'), stringsAsFactors=F)
  colnames(deckMap) <- c()
  rownames(deckMap) <- c()
  
  #Get Solutions Map
  solutionsMap <- defSolutionsMap(deckMap, plateMap, drugInfo)
  for_stock <- subset(solutionsMap, Labware=='labware_9' & !is.na(Vol))$Slot
  
  #B. Perform Initial Dilutions
  nex <- cmd_PreliminaryDilution(plateMap, deckMap, drugInfo, solutionsMap)
  cmdList <- nex[[2]]
  solutionsMap <- nex[[1]]
  
  #C. Dilution to Deep Well Plates
  cmdList <- cmd_FirstDWDilution(plateMap, deckMap, drugInfo, solutionsMap, 
                                 cmdList, wellInfo, serial_step, nPlate)
  
  #D. Serial Dilution on Deep Well Plates
  cmdList <- cmd_MainSerialDilution(cmdList, plateMap, deckMap, solutionsMap, wellInfo, nPlate, colRat)
  
  #E. Filling Outer wells in deep well
  cmdList <- cmd_FillOuterWells_DeepWell(deckMap, nPlate, cmdList, 
                                         wellInfo, solutionsMap, plateMap)
  
  #F. Distributing to Plates
  cmdList <- cmd_Distribute_all(deckMap, nPlate, plateMap, cmdList, wellInfo)
  
  #CALCULATIONS----------
  amtList <- cal_ReqSolutions(solutionsMap, plateMap, cmdList, deckMap)
  deckMap <- cal_correctDeck(deckMap, cmdList, nPlate)
  dilTube <- cal_dilTube(solutionsMap, deckMap)
  #calculating excess tubes
  excTubes <- cal_amtList_Excess(amtList, cmdList, deckMap)
  amtList <- excTubes[[2]]
  cmdList <- excTubes[[1]]
  cmdList <- subset(cmdList, Comment!='skip')
  #set output---------
  deck_map <- cbind.data.frame(deckMap, matrix(, nrow=length(deckMap[,1]), ncol=6), stringsAsFactors=F)
  colnames(deck_map) <- colnames(cmdList)
  
  #robot input
  OT2_cmd <- rbind.data.frame(c('>Amount List', replicate(7, "")),
                              cbind.data.frame(amtList, matrix(, nrow=length(amtList[,1]), ncol=3), stringsAsFactors=F),
                              c('>CommandLines', replicate(7, "")), stringsAsFactors=F)
  colnames(OT2_cmd) <- colnames(cmdList)
  
  OT2_cmd <- rbind.data.frame(OT2_cmd, cmdList, c('>DeckMap', replicate(7, "")), deck_map, stringsAsFactors=F)
  rownames(OT2_cmd) <- c()
  OT2_cmd <<- OT2_cmd
  
  #Calculating number of empty tubes required
  for_dil <- subset(solutionsMap, Labware=='labware_9' & !(Slot %in% for_stock) & !is.na(Vol))$Slot
  
  #user output
  amtList <- user_AmountList(amtList, deckMap)
  if(length(for_dil)>0){
    for_dil <- paste(for_dil, collapse=', ')
    dilTubes <- cbind.data.frame('9', for_dil, "1.5 mL Eppendorf tube", "(empty)", "", "", "", stringsAsFactors=F)
    colnames(dilTubes) <- colnames(amtList)
    amtList <- rbind.data.frame(amtList, dilTubes, stringsAsFactors=F)
  }
  
  deck_map <- rbind(as.character(c(10:12)),
                    t(deckMap[c(10:12),2]),
                    as.character(c(7:9)),
                    t(deckMap[c(7:9),2]),
                    as.character(c(4:6)),
                    t(deckMap[c(4:6),2]),
                    as.character(c(1:3)),
                    t(deckMap[c(1:3),2]))
  
  deck_map <- cbind.data.frame(deck_map, matrix(, nrow=length(deck_map[,1]), ncol=4), stringsAsFactors=F)
  colnames(deck_map) <- colnames(amtList)
  user_cmd <- rbind.data.frame(amtList, c('>> OT2 Deck map <<', replicate(6, "")), deck_map, stringsAsFactors=F)
  rownames(user_cmd) <- c()
  user_cmd <<- user_cmd
  return(amtList)
}