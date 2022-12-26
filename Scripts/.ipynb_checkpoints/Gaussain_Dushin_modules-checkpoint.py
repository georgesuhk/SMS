#!/usr/bin/env python
# coding: utf-8

# In[2]:

import sys
sys.path.append('../')
sys.path.append('.')

import time, os, pickle
import numpy as np
import pandas as pd

import pybel as pb

import rdkit
from ase import Atoms
from rdkit.Chem import AllChem
from rdkit import rdBase,Chem
from rdkit.Chem import AllChem,Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdMolDescriptors
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

import tempfile                                                                                                                                                                                                   
import os
import subprocess
import time
import cclib

import paramiko
from pathlib import Path

from matplotlib import pyplot as plt



#================================================== Functions ================================================================

def readSmilesFile(fname='gen__2.smi'):
    my_file = Path(fname)
    if my_file.is_file():
        # file exists
        dict_smi_list_file={}
        with open(fname) as out:
                    for i,line in enumerate(out.readlines()):
                        line=line.split()
                        dict_smi_list_file[line[0]]=line[1]

                    print('Number of lines: {}'.format(i))
        #print(dict_smi_list_file.keys())
        ms=[Chem.MolFromSmiles(x) for x in dict_smi_list_file.keys()]
        return ms
    else: 
        print(fname+" does not exist in this folder")
        
        
def generate_XYZgeo_from_mol(Mol_list, mol_list_name, moleculenumber, output_folder, save=True):
    Mol_list[moleculenumber] = Chem.AddHs(Mol_list[moleculenumber]) 
    AllChem.Compute2DCoords(Mol_list[moleculenumber])
    for c in Mol_list[moleculenumber].GetConformers():  #only single conformer is chosen currently, might alter in future
        positions=c.GetPositions()   
    atoms_symbols=np.array([at.GetSymbol() for at in Mol_list[moleculenumber].GetAtoms()])
    atoms=Atoms(atoms_symbols, positions=positions)
    num_atoms = len(atoms)
    types = atoms.get_chemical_symbols()
    all_atoms = zip(types, atoms.get_positions())
    a_str = str(num_atoms) + "\n" + "\n"
    a_str_gauss=""
    for atom in all_atoms:
        a_str += atom[0] + " " + " ".join([str(x) for x in atom[1]]) + "\n"
        a_str_gauss += atom[0] + " " + " ".join([str(x) for x in atom[1]]) + "\n"
        
    #writing xyz coord file
    if save:
        text_file = open(output_folder+mol_list_name+"_"+str(moleculenumber)+".xyz", "wt")
        n = text_file.write(a_str)
        text_file.close()
        
    return a_str_gauss, num_atoms

def get_xyz_optimised(mol_charge_state, output_folder, hostname, username, password, remote_folder):
    gaussian_output = mol_charge_state+".log"
    transfer_from_HPC(hostname, username, password, os.getcwd()+"/", remote_folder, gaussian_output)
    gaussian_data = cclib.io.ccread(gaussian_output)

    #generate xyz file from gaussian log output
    num_atoms = gaussian_data.natom

    coords = gaussian_data.atomcoords
    coords = coords[-1,:,:]


    atomnos = gaussian_data.atomnos

    xyz_str = str(num_atoms) + "\n" + "\n"
    xyz_str_gauss=""

    for i in range(0,len(coords)):
        xyz_str += str(atomnos[i]) + " " + " ".join([str(x) for x in coords[i]]) + "\n"
        xyz_str_gauss += str(atomnos[i]) + " " + " ".join([str(x) for x in coords[i]]) + "\n"

    text_file = open(output_folder+mol_charge_state+".xyz", "wt")
    n = text_file.write(xyz_str)
    text_file.close()
    
    os.remove(gaussian_output)
    print("Obtained opt geometry from optimised log file for: ", mol_charge_state)

    return xyz_str_gauss, num_atoms
    


def generate_Gjf_optfreqStates(Mol_list, mollist_name, Mol_xyz, molnumber, Excited_state, charge, output_folder = os.getcwd()+"/", 
                               functional_basisset='b3lyp/6-311g(d,p)', number_of_procc='32', memory='30GB'):
    # generate gjf for the S0 optimisation
    #molname=Chem.MolToSmiles(Mol_list[molnumber])

    #write the gaussian file
    if charge==0:       
        a_str_gaussian="0 1\n"+Mol_xyz
        gjf_file=mollist_name+"_"+str(molnumber)+"_C"+str(charge)+"S"+str(Excited_state)

    if charge==-1:      
        a_str_gaussian="-1 2\n"+Mol_xyz
        gjf_file=mollist_name+"_"+str(molnumber)+"_C"+str(charge)+"S"+str(Excited_state)
    
    if charge==-2:      
        a_str_gaussian="-2 3\n"+Mol_xyz
        gjf_file=mollist_name+"_"+str(molnumber)+"_C"+str(charge)+"S"+str(Excited_state)
        
    if charge==-3:      
        a_str_gaussian="-3 4\n"+Mol_xyz
        gjf_file=mollist_name+"_"+str(molnumber)+"_C"+str(charge)+"S"+str(Excited_state)

    if charge==1:       
        a_str_gaussian="1 2\n"+Mol_xyz
        gjf_file=mollist_name+"_"+str(molnumber)+"_C"+str(charge)+"S"+str(Excited_state)
        
    if charge==2:      
        a_str_gaussian="2 3\n"+Mol_xyz
        gjf_file=mollist_name+"_"+str(molnumber)+"_C"+str(charge)+"S"+str(Excited_state)

    a_str_gaussian=mollist_name+"_"+str(molnumber)+" S"+str(Excited_state)+" opt\n\n"+a_str_gaussian
    if Excited_state==0:
        a_str_gaussian="#p opt freq NoSymm "+functional_basisset+"\n\n"+a_str_gaussian
    else:
        a_str_gaussian="#p opt freq NoSymm td(NStates=20,Root="+str(Excited_state)+")  "+functional_basisset+"\n\n"+a_str_gaussian
        
    a_str_gaussian="%nprocshared="+number_of_procc+"\n\n"+a_str_gaussian
    a_str_gaussian="%mem="+memory+"\n"+a_str_gaussian
    a_str_gaussian="%chk="+gjf_file+".chk\n"+a_str_gaussian
    a_str_gaussian=a_str_gaussian+"\n"

    #print(a_str_gaussian)
    text_file = open(output_folder+gjf_file+".gjf", "wt")
    n = text_file.write(a_str_gaussian)
    text_file.close()
    return gjf_file


def generate_Gjf_theodore(Mol_list, mollist_name, Mol_xyz, molnumber, Excited_state, charge, output_folder = os.getcwd()+"/", 
                          functional_basisset='b3lyp/6-311g(d,p)', number_of_procc='32', memory='30GB'):
    # generate gjf for running theodore to identify CT states
    # Will only run TDDFT

    #write the gaussian file
    if charge==0:       
        a_str_gaussian="0 1\n"+Mol_xyz
        gjf_file=mollist_name+"_"+str(molnumber)+"_C"+str(charge)+"S"+str(Excited_state)

    if charge==-1:      
        a_str_gaussian="-1 2\n"+Mol_xyz
        gjf_file=mollist_name+"_"+str(molnumber)+"_C"+str(charge)+"S"+str(Excited_state)
    
    if charge==-2:      
        a_str_gaussian="-2 3\n"+Mol_xyz
        gjf_file=mollist_name+"_"+str(molnumber)+"_C"+str(charge)+"S"+str(Excited_state)
        
    if charge==-3:      
        a_str_gaussian="-3 4\n"+Mol_xyz
        gjf_file=mollist_name+"_"+str(molnumber)+"_C"+str(charge)+"S"+str(Excited_state)

    if charge==1:       
        a_str_gaussian="1 2\n"+Mol_xyz
        gjf_file=mollist_name+"_"+str(molnumber)+"_C"+str(charge)+"S"+str(Excited_state)
        
    if charge==2:      
        a_str_gaussian="2 3\n"+Mol_xyz
        gjf_file=mollist_name+"_"+str(molnumber)+"_C"+str(charge)+"S"+str(Excited_state)

    a_str_gaussian=mollist_name+"_"+str(molnumber)+" S"+str(Excited_state)+" opt\n\n"+a_str_gaussian

    a_str_gaussian="#p nosymm Geom=NoCrowd opt pop=full iop(9/40=3) td(NStates=20,Root="+str(Excited_state)+")  "+functional_basisset+"\n\n"+a_str_gaussian
        
    a_str_gaussian="%nprocshared="+number_of_procc+"\n\n"+a_str_gaussian
    a_str_gaussian="%mem="+memory+"\n"+a_str_gaussian
    a_str_gaussian="%chk="+gjf_file+".chk\n"+a_str_gaussian
    a_str_gaussian=a_str_gaussian+"\n"

    #print(a_str_gaussian)
    text_file = open(output_folder+gjf_file+".gjf", "wt")
    n = text_file.write(a_str_gaussian)
    text_file.close()
    return gjf_file

#Reading from a checkpoint file

        
        
        
        
def generate_shfile_gaussian( output_folder= os.getcwd()+"/", walltime='36:00:01',    memory="30gb",    ncpus=32,    timeout="36h",                             molnumber=0,Excited_state=0,mollist_name="",gjf_file=""):
    #write an sh file for the job

    # /rds/general/user/ma11115/home/testssh
    #gjf_file=mollist_name+"_"+str(molnumber)+"_S"+str(Excited_state)
    # Job submission script
    qsub_script="#!/bin/sh\n"
    qsub_script+="#PBS -l walltime="+walltime+"\n"
    qsub_script+="#PBS -l select=1:ncpus="+str(ncpus)+":mem="+memory+":avx=true\n\n"
    qsub_script+="module load gaussian/g16-c01-avx\n"
    qsub_script+="cp $PBS_O_WORKDIR/"+gjf_file+".gjf ./\n"
    qsub_script+="timeout "+timeout+" g16 "+gjf_file+".gjf \n"
    qsub_script+="formchk "+gjf_file+".chk "+gjf_file+".fchk\n"
    qsub_script+="""cp *.log  $PBS_O_WORKDIR
    cp *.chk  $PBS_O_WORKDIR
    cp *.fchk  $PBS_O_WORKDIR

    """
    text_file = open(output_folder+gjf_file+".sh", "wt",newline='\n')
    n = text_file.write(qsub_script)
    text_file.close()
    

def run_job_HPC(hostname = 'login.hpc.ic.ac.uk',username = 'gs920', password="Illusions9", remotefolder="/rds/general/user/gs920/home/PythonTrial", molnumber=0,Excited_state=0,mollist_name="test",gjf_file="test"):
    
    
#run the job in the HPC 
    # Hostname of a 'submission host' in the Grid Engine cluster we want to submit a job to

    #gjf_file=mollist_name+"_"+str(molnumber)+"_S"+str(Excited_state)
    with paramiko.SSHClient() as client:
        client.load_system_host_keys()
        client.set_missing_host_key_policy(paramiko.WarningPolicy)

        # Establish SSH connection
        client.connect(hostname, username=username,password=password)

        # Establish SFTP connection ftp_client=ssh.open_sftp()
        stdin, stdout, stderr = client.exec_command("mkdir -p "+remotefolder)
        with client.open_sftp() as sftp0:
            # Push job submission script to a particular path on the cluster
            sftp0.put(gjf_file+".sh", remotefolder+gjf_file+".sh")
            sftp0.put(gjf_file+".gjf", remotefolder+gjf_file+".gjf")
        # Submit our Grid Engine job by running a remote 'qsub' command over SSH
        stdin, stdout, stderr = client.exec_command("cd "+remotefolder+" \n "+
                                                    "chmod +rwx "+gjf_file+".sh\n"+
                                                    "/opt/pbs/bin/qsub "+gjf_file+".sh"+"\n")
                                                     #"find -type f -name '*.sh.*' -delete")

        # Show the standard output and error of our job
        print("Standard error:")
        print(stderr.read())
        print("Exit status: {}".format(stdout.channel.recv_exit_status()))
        
        
        client.close()

#run qstat
def run_qstat(hostname, username, password):
    #runs qstat on the hpc to check queue progress


    with paramiko.SSHClient() as client:
        client.load_system_host_keys()
        client.set_missing_host_key_policy(paramiko.WarningPolicy)
        client.connect(hostname, username=username,password=password)
        stdin, stdout, stderr = client.exec_command("/opt/pbs/bin/qstat ")
        
        
        q_status = stdout.read()
        print("Standard output:")
        print(q_status)

        print("Standard error:")
        print(stderr.read())

        print("Exit status: {}".format(stdout.channel.recv_exit_status()))
        client.close()
    q_status = str(q_status)
    return q_status

def check_for_completion(hostname, username, password, max_loops = 20, wait_time = 180):
    loops_ran = 0
    while loops_ran < max_loops:
        q_status = run_qstat(hostname, username, password)
        
        if q_status == "b''": #standard output if all queues are finished
            print("\n")
            print("Gaussian completed for all.")
            print("\n")
            return "Gaussian completed"
        
        loops_ran += 1
        time.sleep(wait_time)
    
    return "Timed Out"


def run_dushin(hostname, username, password, remote_folder, remote_bin_path, remote_dushin_folder, mol_list, mol_list_name, mol_num, initial_state, final_state):
    mol_serial= mol_list_name+"_"+str(mol_num)
    print("Dushin operating for: ",mol_serial," ",initial_state," ",final_state)


    with paramiko.SSHClient() as client:
        client.load_system_host_keys()
        client.set_missing_host_key_policy(paramiko.WarningPolicy)
        client.connect(hostname, username=username,password=password)

        stdin, stdout, stderr = client.exec_command("cd "+remote_folder+" \n "+ "export "+remote_bin_path+" \n"+ 
                                                    "find -type f -name '*.sh.*' -delete"+" \n"+
                                                    "cp "+mol_serial+"_"+initial_state+".log "+remote_dushin_folder+" \n"+
                                                    "cp "+mol_serial+"_"+initial_state+".fchk "+remote_dushin_folder+" \n"+
                                                    "cp "+mol_serial+"_"+final_state+".log "+remote_dushin_folder+" \n"+
                                                    "cp "+mol_serial+"_"+final_state+".fchk "+remote_dushin_folder+" \n"+
                                                    "cd "+remote_dushin_folder+" \n"+
                                                    
                                                    "split_log.sh "+mol_serial+"_"+initial_state+".log \n"+                                             
                                                    "mv freq_1_"+mol_serial+"_"+initial_state+".log "\
                                                    +mol_serial+"_"+initial_state+".log" +" \n"+ 
                                                    #"mv "+mol_serial+"_"+initial_state+".fchk "+mol_serial+"_"+initial_state+".fchk" +" \n"+
                                                    
                                                    "split_log.sh "+mol_serial+"_"+final_state+".log \n"+
                                                    "mv freq_1_"+mol_serial+"_"+final_state+".log "\
                                                    + mol_serial+"_"+final_state+".log" +" \n"+
                                                    #"mv "+mol_serial+"_"+final_state+".fchk "+mol_serial+"_"+final_state+".fchk" +" \n"+\
                                                    
                                                    #deleting useless files
                                                    "find -type f -name '*freq_2*' -delete"+" \n"+

                                                    "rundushin.sh "+mol_serial+" "+initial_state+" "+final_state)


        print("Standard output:")
        print(stdout.read())
        print("Standard error:")
        print(stderr.read())
        print("Exit status: {}".format(stdout.channel.recv_exit_status()))
        client.close()

        

def transfer_from_HPC(hostname, username, password, local_folder, remote_folder, filename):
    #copy files to local
    

    #filename = 'dushinbenzene_fam_S0_S1.log'
    with paramiko.SSHClient() as client:
        client.load_system_host_keys()
        client.set_missing_host_key_policy(paramiko.WarningPolicy)
        client.connect(hostname, username=username,password=password)
        sftp_client = client.open_sftp()

        sftp_client.get(remote_folder+filename, local_folder+filename)


        sftp_client.close()
        client.close()
        
        
def reading_gaussian_output(hostname, username, password, remote_folder, gaussian_save_folder, mol_serial, functional_basis):
            
    gaussian_output = mol_serial+".log"
    transfer_from_HPC(hostname, username, password, os.getcwd()+"/", remote_folder, gaussian_output)
    gaussian_data = cclib.io.ccread(gaussian_output)

    #useful information to store
    n_atoms = gaussian_data.natom
    molecular_mass = sum(gaussian_data.atommasses)
    homo_energy = gaussian_data.moenergies[0][gaussian_data.homos[0]]
    lumo_energy = gaussian_data.moenergies[0][gaussian_data.homos[0]+1]
    scf_energy = gaussian_data.scfenergies[-1]
    dipole_moment = [gaussian_data.moments[1]]
    excited_state_energy = float('NaN')
    oscillator_strength = float('NaN')
    
    #Excited state values
    if mol_serial.split("S")[1] == "1":   #for excited state calculations
        excited_energy = [gaussian_data.etenergies]
        oscillator_strength = [gaussian_data.etoscs]
    
    useful_gaussian = pd.DataFrame({"n_atoms":n_atoms, "molecular_mass":molecular_mass, "homo_energy":homo_energy, "lumo_energy":lumo_energy, "scf_energy":scf_energy, "dipole_moment":dipole_moment, "excited_state_energy":excited_state_energy, "oscillator_strength":oscillator_strength})
    
    gaussian_filename = mol_serial+".csv"
    
    #determining which folder to save in

    renamed_fb = functional_basis.split("/")[0] + "_" + functional_basis.split("/")[1].split("(")[0] 
    gaussian_save_folder = gaussian_save_folder+renamed_fb +"/"
    if not os.path.exists(gaussian_save_folder):
        os.makedirs(gaussian_save_folder)

    #storing
    useful_gaussian.to_csv(gaussian_save_folder+gaussian_filename)
    os.remove(gaussian_output)
    print("Gaussian outputs read and stored for {}".format(mol_serial))
            
        
def reading_graphing_dushin_output(dushin_output_folder, dushin_save_folder, image_save_folder, filename, functional_basis):
    #reading dushin output
    cmtoEv = 1.24e-4

    dushin_output = open(dushin_output_folder+filename).readlines()
    split_line_top = ' Displacement: in terms of nc of 1 THEN of nc of 2\n'
    split_line_bot = '\n'
    
    #obtaining chunk of output with useful data
    split_index_top = dushin_output.index(split_line_top)+2
    dushin_output = dushin_output[split_index_top::]

    split_index_bot = dushin_output.index(split_line_bot)
    dushin_output = dushin_output[0:split_index_bot]
    
    if dushin_output[-1] != ' ':    #necessary check as some files have additional line for displ modes
        dushin_output = dushin_output[:-1]

    freq_list = []
    Q_list = []
    lam_list = []

    for line in dushin_output:
        #extracting freq and lamb

        #freq split
        split1 = line.split('freq=')[1::]

        #lamb split
        split2 = [split1[0].split('lam='), split1[1].split('lam=')]


        #freq assign
        freq_list.append(float(split2[0][0].split('Q=')[0]))
        freq_list.append(float(split2[1][0].split('Q=')[0]))


        #q assign
        Q_list.append(float(split2[0][0].split('Q=')[1]))
        Q_list.append(float(split2[1][0].split('Q=')[1]))

        #lam assign
        lam_list.append(float([i for i in split2[0][1].split(' ') if i != ''][0]))
        lam_list.append(float(split2[1][1].split('\n')[0]))
        
    #cleaned file name
    plot_filename = filename.split('dushin')[1].split('.log')[0] 
    
    #storing info -----------------------------------------------------------------------------------------------
    
    #determining which folder to save in

    renamed_fb = functional_basis.split("/")[0] + "_" + functional_basis.split("/")[1].split("(")[0] 
    dushin_save_folder = dushin_save_folder+renamed_fb +"/"
    image_save_folder = image_save_folder+renamed_fb +"/"
    if not os.path.exists(dushin_save_folder):
        os.makedirs(dushin_save_folder)

    if not os.path.exists(image_save_folder):
        os.makedirs(image_save_folder)
    
    
    output_data = pd.DataFrame({"Frequencies":freq_list, "Displacement":Q_list, "Reorg Energy":lam_list})
    output_data.to_csv(dushin_save_folder+plot_filename)

    #plotting
    plt.figure()
    bar_width = max(freq_list)/100
    plt.xlim(0,3500)
    lam_list = [x*cmtoEv for x in lam_list] #conversion to Ev
    
    y_lim_upper = 500
    
    if max(lam_list) > y_lim_upper:
        y_lim_upper = max(lam_list) + 50

    #plt.ylim(0,y_lim_upper)
    plt.bar(freq_list, lam_list, bar_width)
    plt.title("$\lambda$ per normal mode \n"+plot_filename)
    plt.ylabel("Reorganization energy $\lambda$ (eV)")
    plt.xlabel("Normal mode frequency $cm^{-1}$ ")

    plt.savefig(image_save_folder+plot_filename+'_reorg_freq.jpeg')

#==================================================================================================================    
    
#pybel implementation
def pybel_readsmiles(smiles_list, opt_steps):
    pybel_mols = [pb.readstring("smi", x) for x in smiles_list]

    for mol in pybel_mols:
        mol.localopt(steps=opt_steps)
    return pybel_mols
        
def pybel_xyz_fromMol(pybel_mol_list, mol_list_name, mol_num, output_folder):
    mol_serial = mol_list_name+"_"+str(mol_num)
    pybel_mol_list[mol_num].write("xyz", output_folder+mol_serial+".xyz", overwrite = True)

    with open(output_folder+mol_serial+".xyz", 'r') as file:
        split = file.read().split("\n\n")
        pybabelxyz = split[1]
        n_atoms = int(split[0])

    return pybabelxyz, n_atoms

#==================================================================================================================

#theodore functions
def get_fragments(mol_list, mol_list_name, mol_num, substructures, theodore_run_folder):
    mol_serial = mol_list_name+"_"+str(mol_num)
    mol = mol_list[mol_num]
    mol_H = Chem.AddHs(mol)
    
    substructure_1 = substructures[2*mol_num]
    substructure_2 = substructures[2*mol_num+1]
    
    frag1 = list(mol.GetSubstructMatches(substructure_1)[0])
    frag2 = list(mol.GetSubstructMatches(substructure_2)[0])
    others = []
    hydrogens = []
    for atom in mol_H.GetAtoms():
        if (atom.GetIdx() in frag1) == False and (atom.GetIdx() in frag2) == False:
            if atom.GetAtomicNum() == 1:
                hydrogens.append(atom)
            else:
                others.append(atom.GetIdx())

    for atom in hydrogens:
        owner_idx = atom.GetBonds()[0].GetOtherAtom(atom).GetIdx() #the index of the atom the hydrogen is connected to
        if owner_idx in frag1:
            frag1.append(atom.GetIdx())
        elif owner_idx in frag2:
            frag2.append(atom.GetIdx())
        else:
            others.append(atom.GetIdx())

    if others != []:
        fragments = [frag1, frag2, others]
    else:
        fragments = [frag1, frag2]
        
    output_data = pd.DataFrame({"fragments":fragments})
    output_data.to_csv(theodore_run_folder+"fragments/"+mol_serial+"_fragments.csv", index=False)
    
    





