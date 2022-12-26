#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os as os
import warnings
cmtoEv = 1.24e-4

def folder_scan(folder):
    files = [file.name for file in os.scandir(folder)]
    files.sort()
    for file in files:
        print(file)
    print("\n")
    



def gaussian_db(Gaussian_folder):
    files = [file.name for file in os.scandir(Gaussian_folder)]
    files.sort()


    #ground state files
    gs_files = []
    for file in files:
        if file.split(".")[0][-4::] == "C0S0":
            gs_files.append(file)

    #excited state files
    ex_files = []
    for file in files:
        if file.split(".")[0][-1] != "0":
            ex_files.append(file)



    #constructing database for gaussian --------------------------------------
    missing_ex_files = []

    #constructing db using first file
    mol_serial = gs_files[0].split("_C")[0]
    gs_0 = pd.read_csv(Gaussian_folder+gs_files[0]).iloc[:, 1:7] 
    gs_backup = pd.read_csv(Gaussian_folder+gs_files[0]).iloc[:, 7:9] 

    try:
        ex_0 = pd.read_csv(Gaussian_folder+mol_serial+"_C0S1.csv").iloc[:, 7:9] 
        db_gauss = pd.concat([gs_0, ex_0], axis=1)
    except:
        missing_ex_files.append(mol_serial)
        db_gauss = pd.concat([gs_0, gs_backup], axis=1)

    db_gauss = db_gauss.rename(index = lambda x:mol_serial)

    #adding other files as rows
    for i in range(1, len(gs_files)):
        mol_serial = gs_files[i].split("_C")[0]
        gs_row = pd.read_csv(Gaussian_folder+gs_files[i]).iloc[:, 1:7]
        gs_backup = pd.read_csv(Gaussian_folder+gs_files[i]).iloc[:, 7:9] 

        try:
            ex_row = pd.read_csv(Gaussian_folder+mol_serial+"_C0S1.csv").iloc[:, 7:9] 
            row = pd.concat([gs_row, ex_row], axis=1)
        except:
            missing_ex_files.append(mol_serial)
            row = pd.concat([gs_row, gs_backup], axis=1) 

        row = row.rename(index = lambda x:mol_serial)
        db_gauss = pd.concat([db_gauss, row], ignore_index=False)
    
    return db_gauss

def DA_db(Reorg_folder):
    files = [file.name for file in os.scandir(Reorg_folder)]
    files.sort()

    encountered = []
    reorg_db = pd.DataFrame({"Donor*":[1], "Acceptor*":[1], "Cation":[1], "Anion":[1], "Dyad":[1], "S1-S0":[1], "S2-S0":[1]},index=['start'])


    for file in files:

        mol_serial = file.split("_C")[0]
        transition = file.split("_")[2]+"_"+file.split("_")[3]

        #check if first occurance
        if mol_serial not in encountered:
            reorg_row = pd.DataFrame({"Donor*":[np.float('NaN')], "Acceptor*":[np.float('NaN')], "Cation":[np.float('NaN')], "Anion":[np.float('NaN')], "Dyad":[np.float('NaN')], "S1-S0":[np.float('NaN')], "S2-S0":[np.float('NaN')]},index=[mol_serial])
            reorg_db = pd.concat([reorg_db, reorg_row], ignore_index = False)
            encountered.append(mol_serial)    


        #C0S1_C1S0 (LE on donor to cation)
        if transition == "C0S1_C1S0" or transition == "C1S0_C0S1":
            reorg_db.loc[mol_serial]["Donor*"] = pd.read_csv(Reorg_folder+file)["Reorg Energy"].sum()*cmtoEv
            
        #C0S1_C-1S0 (LE on acceptor to anion)
        if transition == "C0S1_C-1S0" or transition == "C-1S0_C0S1":
            reorg_db.loc[mol_serial]["Acceptor*"] = pd.read_csv(Reorg_folder+file)["Reorg Energy"].sum()*cmtoEv
            
        #C0S0_C1S0 (Cation)
        if transition == "C0S0_C1S0" or transition == "C1S0_C0S0":
            reorg_db.loc[mol_serial]["Cation"] = pd.read_csv(Reorg_folder+file)["Reorg Energy"].sum()*cmtoEv

        #C0S0_C-1S0 (Anion)
        if transition == "C0S0_C-1S0" or transition == "C-1S0_C0S0":
            reorg_db.loc[mol_serial]["Anion"] = pd.read_csv(Reorg_folder+file)["Reorg Energy"].sum()*cmtoEv

        #C0S2_C0S1 (Dyad LE to CT)
        if transition == "C0S2_C0S1" or transition == "C0S1_C0S2":
            reorg_db.loc[mol_serial]["Dyad"] = pd.read_csv(Reorg_folder+file)["Reorg Energy"].sum()*cmtoEv
            
        #C0S1_C0S0 (Dyad S1 to ground)
        if transition == "C0S1_C0S0" or transition == "C0S0_C0S1":
            reorg_db.loc[mol_serial]["S1-S0"] = pd.read_csv(Reorg_folder+file)["Reorg Energy"].sum()*cmtoEv
            
        #C0S2_C0S0 (Dyad S2 to ground)
        if transition == "C0S2_C0S0" or transition == "C0S0_C0S2":
            reorg_db.loc[mol_serial]["S2-S0"] = pd.read_csv(Reorg_folder+file)["Reorg Energy"].sum()*cmtoEv


    reorg_db = reorg_db.iloc[1::,:]
    reorg_db = reorg_db.replace(np.float('NaN'),"-")

        #Other transitions can be added, eg. C0S0_C1S0 (Cation)

    reorg_db.round(decimals = 3)
    return reorg_db


#comparison of transitions
def compare(transitions, Reorg_folder, figures=True):
    
    trans_data_list = {}
    for trans in transitions:
        
        mol_serial = trans.split("_")[0] + "_" + trans.split("_")[1]
        charge_state_1 = trans.split("_")[2]
        charge_state_2 = trans.split("_")[3]
        
        try:
            trans_data_list[trans] = pd.read_csv(Reorg_folder+trans)
        except:
            alt_name = mol_serial+"_"+charge_state_2+"_"+charge_state_1
            trans_data_list[trans] = pd.read_csv(Reorg_folder+alt_name)
                  
    
    #obtain max reorg energy to scale the y axis
    max_energies = []
    for data in trans_data_list.values():
        max_energies.append(max(data["Reorg Energy"]*cmtoEv))

    max_reorg = 1.1*max(max_energies)

    
    if figures:
        fig, axs = plt.subplots(nrows=len(transitions), ncols=1, figsize=(6,2*len(transitions)))
        bar_width = 2000/85
        
        for i in range(0,len(axs)):
            data = trans_data_list[transitions[i]]
            
            axs[i].bar(data["Frequencies"], data["Reorg Energy"]*cmtoEv, bar_width)
            axs[i].set_xlim([0,3000])
            axs[i].set_ylim([0,max_reorg])
            axs[i].set_title(transitions[i])
        
        plt.tight_layout()
        plt.show()
        
    #compare matrix
    
    
    trans_data = trans_data_list[transitions[0]]
    trans_low  = trans_data[trans_data['Frequencies']<500 ]["Reorg Energy"].sum()*cmtoEv
    trans_mid  = trans_data[trans_data['Frequencies']<2000][trans_data['Frequencies']>=500]["Reorg Energy"].sum()*cmtoEv
    trans_high = trans_data[trans_data['Frequencies']>=2000]["Reorg Energy"].sum()*cmtoEv
    
    compare_matrix = pd.DataFrame({"Total": trans_data["Reorg Energy"].sum()*cmtoEv, "0-500": trans_low, "500-2000": trans_mid, "2000+":trans_high}, index=[transitions[0]])
    
    for i in range(1, len(transitions)):
        trans_data = trans_data_list[transitions[i]]
        trans_low  = trans_data[trans_data['Frequencies']<500 ]["Reorg Energy"].sum()*cmtoEv
        trans_mid  = trans_data[trans_data['Frequencies']<2000][trans_data['Frequencies']>=500]["Reorg Energy"].sum()*cmtoEv
        trans_high = trans_data[trans_data['Frequencies']>=2000]["Reorg Energy"].sum()*cmtoEv
    
        compare_row = pd.DataFrame({"Total": trans_data["Reorg Energy"].sum()*cmtoEv, "0-500": trans_low, "500-2000": trans_mid, "2000+":trans_high}, index=[transitions[i]])
        compare_matrix = pd.concat([compare_matrix, compare_row], ignore_index = False)
      
    
    
    print(compare_matrix)
    return compare_matrix
    
    
    
    
    


#comparison of dyad LE to CT vs same transition but in the view of fragments
#comparison of dyad LE to CT vs same transition but in the view of fragments
def Dyad_Fragment_Comparison(donor, acceptor, dyad, starter = "Donor", figures=True, Reorg_folder = "Reorg_Energy/Stage2/", bar_width = 3000/100):
    
    #get dyad data
    try:
        dyad_data = pd.read_csv(Reorg_folder+dyad+"_C0S1_C0S2")
    except:
        dyad_data = pd.read_csv(Reorg_folder+dyad+"_C0S2_C0S1")
    
    if starter == "Donor": #Transition starting from LE on Donor
        
        #get donor data
        try:
            donor_data = pd.read_csv(Reorg_folder+donor+"_C0S1_C1S0")
        except:
            donor_data = pd.read_csv(Reorg_folder+donor+"_C1S0_C0S1")

        #get acceptor data
        try:
            acceptor_data = pd.read_csv(Reorg_folder+acceptor+"_C0S0_C-1S0")
        except:
            acceptor_data = pd.read_csv(Reorg_folder+acceptor+"_C-1S0_C0S0")
            
    if starter == "Acceptor": #Transition starting from LE on Acceptor
        
        #get donor data
        try:
            donor_data = pd.read_csv(Reorg_folder+donor+"_C0S0_C1S0")
        except:
            donor_data = pd.read_csv(Reorg_folder+donor+"_C1S0_C0S0")

        #get acceptor data
        try:
            acceptor_data = pd.read_csv(Reorg_folder+acceptor+"_C0S1_C-1S0")
        except:
            acceptor_data = pd.read_csv(Reorg_folder+acceptor+"_C-1S0_C0S1")
        
        

    
    bins = np.linspace(0,3000,round(len(dyad_data["Frequencies"])*1.5))

    donor_data["bins"] = pd.cut(donor_data["Frequencies"], bins)
    acceptor_data["bins"] = pd.cut(acceptor_data["Frequencies"], bins)
    dyad_data["bins"] = pd.cut(dyad_data["Frequencies"], bins)


    DA_reorg = donor_data.groupby(["bins"]).sum()["Reorg Energy"]*cmtoEv + acceptor_data.groupby(["bins"]).sum()["Reorg Energy"]*cmtoEv
    DA_freqs = []
    for i in range(0,len(bins)-1):
        DA_freqs.append((bins[i] + bins[i+1]) * 0.5)
        
    dyad_reorg = dyad_data.groupby(["bins"]).sum()["Reorg Energy"]*cmtoEv 
    dyad_freqs = []
    for i in range(0,len(bins)-1):
        dyad_freqs.append((bins[i] + bins[i+1]) * 0.5)
    
    
    if figures:
        fig, (ax1, ax2, ax3, ax4) = plt.subplots(4,1, figsize=(10, 10))
        #bar_width = max(dyad_data["Frequencies"])/85

        #Donor
        ax1.bar(donor_data["Frequencies"], donor_data["Reorg Energy"]*cmtoEv, bar_width, color="#0E6BA8")
        ax1.set_xlim([0,3000])
        ax1.set_ylim([0,1.1*max(max(dyad_data["Reorg Energy"]*cmtoEv), max(DA_reorg))])
        donor_title = "Donor: "+donor
        ax1.set_title(donor_title)
        
        #Acceptor
        ax2.bar(acceptor_data["Frequencies"], acceptor_data["Reorg Energy"]*cmtoEv, bar_width, color="red")
        ax2.set_xlim([0,3000])
        ax2.set_ylim([0,1.1*max(max(dyad_data["Reorg Energy"]*cmtoEv), max(DA_reorg))])
        acceptor_title = "Acceptor: "+acceptor
        ax2.set_title(acceptor_title)
 
        #Dyad
        ax3.bar(dyad_data["Frequencies"], dyad_data["Reorg Energy"]*cmtoEv, bar_width, color = "#917C78")
        ax3.set_xlim([0,3000])
        ax3.set_ylim([0,1.1*max(max(dyad_data["Reorg Energy"]*cmtoEv), max(DA_reorg))])
        dyad_title = "Dyad: "+dyad
        ax3.set_title(dyad_title)
        
        #coplotted 
        ax4.bar(dyad_freqs, dyad_reorg, bar_width, label="Dyad", color = "#917C78", alpha=1)
        ax4.bar(DA_freqs, DA_reorg, bar_width, label="Donor + Acceptor", color = "purple", alpha=0.6)
        ax4.set_xlim([0,3000])
        ax4.set_ylim([0,1.1*max(max(dyad_reorg), max(DA_reorg))])
        ax4.set_title("Co Plotted")
        ax4.legend(loc=1)
        
        plt.tight_layout()
        
        plt.show()
    
    dyad_low  = dyad_data[dyad_data['Frequencies']<500 ]["Reorg Energy"].sum()*cmtoEv
    dyad_mid  = dyad_data[dyad_data['Frequencies']<2000][dyad_data['Frequencies']>=500]["Reorg Energy"].sum()*cmtoEv
    dyad_high = dyad_data[dyad_data['Frequencies']>=2000]["Reorg Energy"].sum()*cmtoEv
    
    DA_low  = donor_data[donor_data['Frequencies']<500]["Reorg Energy"].sum()*cmtoEv + acceptor_data[acceptor_data['Frequencies']<500]["Reorg Energy"].sum()*cmtoEv
    DA_mid  = donor_data[donor_data['Frequencies']<2000][donor_data['Frequencies']>=500]["Reorg Energy"].sum()*cmtoEv + acceptor_data[acceptor_data['Frequencies']<2000][acceptor_data['Frequencies']>=500]["Reorg Energy"].sum()*cmtoEv
    DA_high = donor_data[donor_data['Frequencies']>=2000]["Reorg Energy"].sum()*cmtoEv + acceptor_data[acceptor_data['Frequencies']>2000]["Reorg Energy"].sum()*cmtoEv
    
    
    
    compare_matrix = pd.DataFrame({"Total": dyad_reorg.sum(), "0-500": dyad_low, "500-2000": dyad_mid, "2000+":dyad_high}, index=[dyad])
    DA_row = pd.DataFrame({"Total": donor_data["Reorg Energy"].sum()*cmtoEv + acceptor_data["Reorg Energy"].sum()*cmtoEv, "0-500": DA_low, "500-2000": DA_mid, "2000+":DA_high}, index=["DA"])
    
    compare_matrix = pd.concat([compare_matrix, DA_row], ignore_index = False)
    print(compare_matrix)