#!/bin/bash
cd /home/inferno/Bsc_Project/TheoDORE_3.0/ 
source setpaths.bash 
cd /home/inferno/Bsc_Project/Scripts/Theodore_runs/
conda deactivatetheodore analyze_tden -f phenothiazine-pyrene_0_C0S1.in 
mv tden_summ.txt theodore_output/phenothiazine-pyrene_0_C0S1.txt