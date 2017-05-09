#! /bin/sh
#this is uses the updated BUSTED and BUSTED-SRV bf's from the dev 2.3 branch of hyphy on github
curdir=$1;
#curdir="/home3/sadie/data/SimBUSTEDSRV/BUSTEDSel/*.nex";
#echo "input number of jobs to queue: "
#read jobnum

jobnum=200;
scheduled=0;

#FILES=$1/*[[:digit:]].[[:digit:]][[:digit:]]
#FILES2=$1/*[[:digit:]].[[:digit:]]

#if [ ! -d "/home/swisotsky/data/busted/busted" ]; then
#mkdir "/home/swisotsky/data/busted/busted"
#fi


for k in $curdir/*replicate.[0-9];
do

    echo $k
    bn=$(basename $k)
    sbatch -p short  --job-name="BUSTED"  --wrap="((echo 1; echo $k; echo y; echo 1; echo d) | ~/bin/hyphy/HYPHYMP  ~/bin/hyphy/res/TemplateBatchFiles/SelectionAnalyses/BUSTED.bf )"
    sbatch  --job-name="BUSTED-SRV"  --wrap="((echo 1; echo $k; echo y; echo 1; echo d) | ~/bin/hyphy/HYPHYMP  ~/bin/hyphy/lib/hyphy/TemplateBatchFiles/BUSTED-SRV.bf )" 
done;

for k in $curdir/*replicate.[0-9][0-9];
do

    echo $k
    bn=$(basename $k)
    sbatch -p short  --job-name="BUSTED"  --wrap="((echo 1; echo $k; echo y; echo 1; echo d) | ~/bin/hyphy/HYPHYMP  ~/bin/hyphy/res/TemplateBatchFiles/SelectionAnalyses/BUSTED.bf )"
    sbatch  --job-name="BUSTED-SRV"  --wrap="((echo 1; echo $k; echo y; echo 1; echo d) | ~/bin/hyphy/HYPHYMP  ~/bin/hyphy/lib/hyphy/TemplateBatchFiles/BUSTED-SRV.bf )"
done;

for k in $curdir/*replicate.[0-9][0-9][0-9];
do

        echo $k
    bn=$(basename $k)
    sbatch -p short  --job-name="BUSTED"  --wrap="((echo 1; echo $k; echo y; echo 1; echo d) | ~/bin/hyphy/HYPHYMP  ~/bin/hyphy/res/TemplateBatchFiles/SelectionAnalyses/BUSTED.bf )"
    sbatch  --job-name="BUSTED-SRV"  --wrap="((echo 1; echo $k; echo y; echo 1; echo d) | ~/bin/hyphy/HYPHYMP  ~/bin/hyphy/lib/hyphy/TemplateBatchFiles/BUSTED-SRV.bf )"
done;




