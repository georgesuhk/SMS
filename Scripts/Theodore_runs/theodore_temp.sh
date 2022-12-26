#!/bin/bash
cd /home/inferno/Bsc_Project/TheoDORE_3.0/ 
source setpaths.bash 
cd /home/inferno/Bsc_Project/Scripts/Theodore_runs/
theodore analyze_tden -f TT-PDI_6_C0S1.in 
printf '



D
A
B
' | theodore plot_frag_decomp 
mv tden_summ.txt theodore_output/TT-PDI_6_C0S1.txt 
mv ehFrag.txt theodore_output/TT-PDI_6_C0S1_eh.txt 
mv frag_decomp.png theodore_output/plots/TT-PDI_6_C0S1_decomp.png 
