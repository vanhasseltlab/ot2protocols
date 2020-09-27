# -*- coding: utf-8 -*-
"""
Created on Sat Sep 26 20:55:41 2020

@author: Sebastian
"""

#FILENAME INPUT--------
path = "C:\\Users\\Sebastian\\Desktop\\MSc Leiden 2nd Year\\##LabAst Works\\FinVersion\\OT2_UpstreamExcelProc_MAIN\\250920200\\ControllerTroubleshooting"
fileName = "CommandList_PlateMapInpt.csv"

#IMPORTING LIBRARIES---------
import csv
import os
import numpy as np
from math import pi
from opentrons import protocol_api

#CUSTOM LIBRARY------
def ReadCSV_Dat(file_name):
    #save all read info into the variable: command_list
    content_list = np.empty(8)
    with open(file_name, 'r') as file:
        cmdCSV = csv.reader(file, delimiter=',')
        for cmdRow in cmdCSV:
            content_list = np.vstack([content_list, cmdRow])
    
    #Find starting point of amount list and command list
    indices = []
    for a in range(len(content_list)):
        if(">" in content_list[a][0]):
            indices.append(a)
    
    #get amount list
    amt_list = content_list[indices[0]+1:indices[1]]
    
    #get command list
    cmd_list = content_list[indices[1]+1:]
            
    return amt_list, cmd_list

def CalTip_Aspirate(solutions_map, cmd_line):
    #get tube type
    tube_loc = [(x[1]==cmd_line[0] and x[2]==cmd_line[1]) for x in solutions_map]
    tube_loc = [i for i, x in enumerate(tube_loc) if x]
    tube_type = solutions_map[tube_loc[0]][4]
    
    #get source amount after aspirated
    src_amt = solutions_map[tube_loc[0]][3]
    
    #if not 96 well plate
    if(tube_type != '96-well'):
        #get dimensions
        if(tube_type=='50 mL Falcon Tube'):
            h_bot = 15.88 #mm
            r = 28.14/2 #mm
            minH = h_bot/2 #mm
        elif(tube_type=='15 mL Falcon Tube'):
            h_bot = 23.36 #mm
            r = 15.62/2 #mm
            minH = h_bot/2 #mm
        else:
            #Tube Dimensions - Eppendorf
            h_bot = 37.8-20 #mm
            r = 8.7/2 #mm
            minH = h_bot/3
        
        #calculate height
        Vmax_bot = pi*r**2*h_bot/3
        
        if(src_amt>Vmax_bot):
            h_tip = h_bot + (src_amt - Vmax_bot)/(pi*r**2)
        else:
            h_tip = ((3*src_amt*h_bot**2)/(pi*r**2))**(1/3)
    
    #if source is a 96-well plate
    else:
        #assign well dimensions
        r = 6.45/2 #mm
        minH = 2 #mm
        
        h_tip = src_amt/(pi*r**2)
    
    #add stab distance; place minimum height into place
    h_tip = max(h_tip-2, minH)
    
    return(h_tip)
 
def CalTip_Dispense(solutions_map, cmd_line, target_well):
    #get tube type
    tube_loc = [(x[1]==cmd_line[2] and x[2]==target_well) for x in solutions_map]
    tube_loc = [i for i, x in enumerate(tube_loc) if x]
    tube_type = solutions_map[tube_loc[0]][4]
    
    #get source amount after dispensed
    src_amt = solutions_map[tube_loc[0]][3]
    
    #if not 96 well plate
    if(tube_type != '96-well'):
        #get dimensions
        if(tube_type=='50 mL Falcon Tube'):
            h_bot = 15.88 #mm
            r = 28.14/2 #mm
            minH = h_bot/2 #mm
        elif(tube_type=='15 mL Falcon Tube'):
            h_bot = 23.36 #mm
            r = 15.62/2 #mm
            minH = h_bot/2 #mm
        else:
            #Tube Dimensions - Eppendorf
            h_bot = 37.8-20 #mm
            r = 8.7/2 #mm
            minH = h_bot/3
        
        #calculate height
        Vmax_bot = pi*r**2*h_bot/3
        
        if(src_amt>Vmax_bot):
            h_tip = h_bot + (src_amt - Vmax_bot)/(pi*r**2)
        else:
            h_tip = ((3*src_amt*h_bot**2)/(pi*r**2))**(1/3)
    
    #if source is a 96-well plate
    else:
        #on top of well
        h_tip = 12.19 #mm
        minH = 12.19 #mm
    
    #add extra distance; place minimum height into place
    h_tip = max(h_tip+2, minH)
    
    return(h_tip)

def Update_Source(solutions_map, cmd_line):
    #get tube location
    tube_loc = [(x[1]==cmd_line[0] and x[2]==cmd_line[1]) for x in solutions_map]
    tube_loc = [i for i, x in enumerate(tube_loc) if x]
    
    #get source amount after dispensed
    solutions_map[tube_loc[0]][3] = solutions_map[tube_loc[0]][3] - float(cmd_line[4])
    
    return(solutions_map)

def Update_Target(solutions_map, cmd_line, target_well, deck_map):
    #get tube location
    tube_loc = [(x[1]==cmd_line[2] and x[2]==target_well) for x in solutions_map]
    tube_loc = [i for i, x in enumerate(tube_loc) if x]
    
    #if tube is not yet registered
    if(len(tube_loc)==0):
        #check target ware type
        ware_type = deck_map[cmd_line[2]]
        
        if('96_wellplate' in ware_type):
            type_target = '96-well'
        elif('1.5ml' in ware_type):
            type_target = '1.5 mL Eppendorf Tube'
        else:
            type_target = '15 mL Falcon Tube'
            
        #generate next item
        regItem = ["NewItem", #name unknown; but it's okay I guess?
                   cmd_line[2], #target labware
                   target_well, #target slot/well
                   float(cmd_line[4]), #initial amount in slot/well
                   type_target] #type of well/slot
        #append
        solutions_map.append(regItem)
    else:
        #get source amount after dispensed
        solutions_map[tube_loc[0]][3] = solutions_map[tube_loc[0]][3] + float(cmd_line[4])
    
    return(solutions_map)
#READ DATA---------
os.chdir(path)
amtList, cmdList = ReadCSV_Dat(fileName)

##############################  SETTINGS  ##############################
dBottom = 4
dTop = 2
aspirateSpeed = 75
dispenseSpeed = 100
deckMap = {
    "labware_1" : 'opentrons_96_tiprack_300ul',
    "labware_2" : 'nest_96_wellplate_200ul_flat',
    "labware_3" : 'opentrons_6_tuberack_falcon_50ml_conical',
    "labware_4" : 'opentrons_15_tuberack_falcon_15ml_conical',
    "labware_5" : 'opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap',
    "labware_6" : 'opentrons_6_tuberack_falcon_50ml_conical',
    "labware_7" : 'opentrons_96_tiprack_300ul',
    "labware_8" : 'opentrons_15_tuberack_falcon_15ml_conical',
    "labware_9" : 'opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap',
    "labware_10": 'opentrons_96_tiprack_300ul',
    "labware_11": '',
    "labware_12": 'TRASH'}
##############################   METADATA   ##############################
metadata = {
    'protocolName': 'OT2_CommandExecuter_vx1',
    'author': 'Sebastian T. Tandar <sebastian.tandar@gmail.com>',
    'description': 'Translating Excel Commands (vx1, 26092020) to Python 2.2 API; height-controlled; Background Structure Check',
    'apiLevel': '2.2'
}
##############################     MAIN     ##############################
def run(protocol: protocol_api.ProtocolContext):
    global cmdList, amtList, deckMap
    
    ############ LOAD LABWARES ############
    tipLocs = []
    for i in range(11):
        #load labware
        labware_name = deckMap["labware_"+str(i+1)]
        if(len(labware_name)>1 and labware_name != 'TRASH'):
            globals()[list(deckMap.keys())[i]] = protocol.load_labware(labware_name, i+1)
        
            #if labware is a tip rack, assign number to tip location(s)
            if('tiprack' in labware_name):
                tipLocs.append(globals()[list(deckMap.keys())[i]])
    
    #load pipettes
    right_pipette = protocol.load_instrument(
        'p300_single', 'right', tip_racks=tipLocs)
    right_pipette.flow_rate.aspirate=aspirateSpeed
    right_pipette.flow_rate.dispense=dispenseSpeed
    
    ########### ASSIGN SOLUTION LOCATIONS ###########
    #initiate map
    solutionsMap = []
    #iterate through all items in amount list
    for item in amtList:
        next_amt = float(item[5])
        
        if(item[0]!='DRUG STOCK'):
            next_amt = next_amt*1000
        next_item = [item[4], item[1], item[3], next_amt, item[2]]
        solutionsMap.append(next_item)
    
    ############ EXECUTE COMMANDS ############
    #iterate through all command lines
    current_tipID = 0 #initiate tip ID
    for i in range(len(cmdList)):
        #subset
        cmdRow = cmdList[i]
        #parse all informations
        source_ware = cmdRow[0]
        source_well = cmdRow[1] #only one source well is allowed
        target_ware = cmdRow[2]
        target_well = cmdRow[3].split(', ')
        transfer_amt = float(cmdRow[4]) #only one transfer amount is allowed
        mix_amt = min(float(cmdRow[5]), 200)
        tipID = int(cmdRow[6]) #each row is performed using a single tip
        
        #pick up tip if needed
        if(tipID != current_tipID):
            right_pipette.pick_up_tip() #pick up tip if tipID changes
            current_tipID = tipID #update tip id
            
        #iterate through all target wells
        for current_target in target_well:
            #update solutions map
            solutionsMap = Update_Source(solutionsMap, cmdRow)
            solutionsMap = Update_Target(solutionsMap, cmdRow, current_target, deckMap)
            
            #calculate aspirate and dispense height
            aspH = CalTip_Aspirate(solutionsMap, cmdRow)
            dspH = CalTip_Dispense(solutionsMap, cmdRow, current_target)
            
            #Main Transfers
            if(mix_amt==0):
                #if no mix
                right_pipette.transfer(transfer_amt,
                                      globals()[source_ware].wells_by_name()[source_well].bottom(aspH),
                                      globals()[target_ware].wells_by_name()[current_target].bottom(dspH),
                                      new_tip='never', blow_out=True, carryover=True)
            else:
                #if mix
                right_pipette.transfer(transfer_amt,
                                      globals()[source_ware].wells_by_name()[source_well].bottom(aspH),
                                      globals()[target_ware].wells_by_name()[current_target].bottom(dspH),
                                      new_tip='never', blow_out=True, mix_before = (3, mix_amt), carryover=True)
            
        #check if tip need to be trashed afterwards
        if(i == len(cmdList)-1):
            #if this is the last operation
            right_pipette.drop_tip()
        elif(int(cmdRow[6]) != int(cmdList[i+1][6])):
            #drop if different tip id is detected
            right_pipette.drop_tip()
        