#!/bin/bash

declare -a arr=( "15to30" "30to50" "50to80" "80to120" "120to170" "170to300" "300to470" "470to600" "600to800" "800to1000" "1000to1400" "1400to1800" "1800to2400" "2400to3200" "3200toInf" )
echo -e 'DatasetName\tRegionA\tRegionB\tRegionC\tRegionD\tRegionAWeighted\tRegionBWeighted\tRegionCWeighted\tRegionDWeighted\tRegionAError\tRegionBError\tRegionCError\tRegionDError' > Regions.log
for i in "${arr[@]}"
do
    echo "InputFile /eos/user/c/chensel/Work/Analysis/LongLivedGluinos/MC/2015/Matthias/LLGDVFiles_Spring15/76_v1/QCD_Pt_${i}/QCD_Pt_${i}.root" >> ConfigSync2.txt
    echo "DatasetName QCD_Pt_${i}_TuneCUETP8M1_13TeV_pythia8" >> ConfigSync2.txt
    ./RunAnalysis ConfigSync2.txt
    mv outputHistograms.root ${i}_histos.root
    head -n -2 ConfigSync2.txt > temp.txt ; mv temp.txt ConfigSync2.txt
done
