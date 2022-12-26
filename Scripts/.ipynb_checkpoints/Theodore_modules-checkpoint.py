import pandas as pd
import subprocess
import os
import paramiko


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

def generate_theodore_in(mol_serial, mol_serial_state, log_file, xyz_folder, theodore_run_folder):
    fragments = pd.read_csv(theodore_run_folder+"fragments/"+mol_serial+"_fragments.csv")

    combined_fragments = "["
    for i in range(0,len(fragments["fragments"])-1):
        combined_fragments += fragments["fragments"][i]+", "
    combined_fragments += fragments["fragments"][len(fragments["fragments"])-1]+"]"

    theodore_in = "rtype='cclib'\n"
    theodore_in += "rfile="+"'"+log_file+"'"+"\n"
    theodore_in += "coor_file="+"'"+xyz_folder+mol_serial+".xyz"+"'"+"\n"
    theodore_in += "coor_format='xyz'"+"\n"
    theodore_in += "at_lists="+combined_fragments+"\n"
    theodore_in += "Om_formula=2"+"\n"
    theodore_in += "eh_pop=1"+"\n"
    theodore_in += "comp_ntos=True"+"\n"
    theodore_in += "comp_dntos=False"+"\n"
    theodore_in += "jmol_orbitals=True"+"\n"
    theodore_in += "molden_orbitals=False"+"\n"
    theodore_in += "prop_list=['Om', 'POS', 'PR', 'CT', 'COH', 'CTnt', 'PRNTO', 'Z_HE', 'RMSeh']"

    text_file = open(theodore_run_folder+mol_serial_state+".in", "wt",newline='\n')
    n = text_file.write(theodore_in)
    text_file.close()
    
    
def run_theodore(mol_serial_state, theodore_installation_folder, theodore_run_folder): 
    run_theo="#!/bin/bash\n"
    run_theo+='cd '+theodore_installation_folder+" \n"
    run_theo+="source setpaths.bash \n"
    run_theo+='cd '+theodore_run_folder+"\n"
    #run_theo+='conda activate '+environment+"\n"                        #if using conda environment (that is causing issues)
    run_theo+="theodore analyze_tden -f "+mol_serial_state+".in"+" \n"
    run_theo+="printf '\n\n\n\nD\nA\nB\n' | theodore plot_frag_decomp"+" \n"
    run_theo+="mv tden_summ.txt "+"theodore_output/"+mol_serial_state+".txt" +" \n"
    run_theo+="mv ehFrag.txt "+"theodore_output/"+mol_serial_state+"_eh.txt \n"
    run_theo+="mv frag_decomp.png "+"theodore_output/plots/"+mol_serial_state+"_decomp.png \n"

    text_file = open(theodore_run_folder+"theodore_temp.sh", "wt",newline='\n')
    n = text_file.write(run_theo)
    text_file.close()
    
    process = subprocess.run(['cd '+theodore_run_folder+"\n"+
                            'bash theodore_temp.sh'],shell=True, stdout = subprocess.PIPE)


                   
                   
                   
         