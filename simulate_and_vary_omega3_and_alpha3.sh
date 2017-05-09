n=1; # n is the value of omega 3 increasing from 1 to 3 
until [ $n = "3.25" ];do
  
    #replace . with _ fornaming reasons
    k=${n/./_}
   #echo $k
   #make directories for files if they do not already exists 
   #files have names based on O3_"value of omega 3"_A3_"value of alpha 3"
   #first make directory for no SRV
    if [ ! -d "/home3/sadie/Simulation_redo_2017/Results/1000_Codons_11x11/O3_"$k"_A3_0005" ]; then
    mkdir "/home3/sadie/Simulation_redo_2017/Results/1000_Codons_11x11/O3_"$k"_A3_0001"
    fi
  # then simulate 100 files with no SRV   and with omega3 = n
   sbatch --job-name="Sim_5000_noSRV" -x node[1-9]  --wrap="((echo /home3/sadie/Simulation_redo_2017/LF-Files/1000_Codons_11x11.LF.bf; echo 0.0005; echo 0.333; echo 0.0005; echo 0.3333; echo 0.0005; echo 0.1; echo 0.5; echo 0.75; echo 0.3; echo $n; echo 101; echo /home3/sadie/Simulation_redo_2017/Results/1000_Codons_11x11/O3_"$k"_A3_0001/5000_Codons_O3_"$k"_A3_0001_sim) | /home3/sadie/bin/hyphy/HYPHYMP /home3/sadie/Simulation_redo_2017/Scripts/BUSTED-SRV-sim.bf)"
  
  #then loop thru different alpha3's 
  a=0.5;
  until [ $a = "3.5" ]; do
     z=${a/./_}
       if [ ! -d "/home3/sadie/Simulation_redo_2017/Results/5000_Codons/O3_"$k"_A3_$z" ]; then
      mkdir "/home3/sadie/Simulation_redo_2017/Results/5000_Codons/O3_"$k"_A3_$z"
       fi
       #simulate 100 files with omega3 = n and alpha3 = a
    sbatch --job-name="Sim_$k_$z" -x node[1-9] --wrap="((echo /home3/sadie/Simulation_redo_2017/LF-Files/5000_Codons.0.busted.LF.bf; echo 0.0005; echo 0.333; echo 0.05; echo 0.3333; echo $a; echo 0.1; echo 0.5; echo 0.75; echo 0.3; echo $n; echo 101; echo /home3/sadie/Simulation_redo_2017/Results/5000_Codons/O3_"$k"_A3_"$z"/5000_Codons_O3_"$k"_A3_"$z"_sim) | /home3/sadie/bin/hyphy/HYPHYMP /home3/sadie/Simulation_redo_2017/Scripts/BUSTED-SRV-sim.bf)"
    

    a=$(awk -v "a=$a" 'BEGIN { print a +0.5}')
    done
  
    #echo /home3/sadie/data/SimBUSTEDSRV/simData/Five_seq/5000_Codons/yesSel_yesSRV/O3_"$k"
    n=$(awk -v "n=$n" 'BEGIN { print n +0.25}')
    #echo $n;
done 