#!/bin/bash

declare -a arrayV=( "60" )
declare -a arrayH=( "130" "140" "150" "160" "170" )
echo -e "Boundaries\tRegionA\tRegionB\RegionC\tRegionD\RegionAW\tRegionBW\RegionCW\tRegionDW" > Total.tsv
for i in "${arrayV[@]}"
do
  sed -i "s/verticalBoundary.*/verticalBoundary ${i}/g" ConfigSync2.txt
  for j in "${arrayH[@]}"
  do
    sed -i "s/horizontalBoundary.*/horizontalBoundary ${j}/g" ConfigSync2.txt
    echo -e "\n\e[31mNow Running runMultiple with \e[1m  H${j}V${i}  \e[0m"
    ./runMultiple.sh
    mv sumhistos.root sumhistosH${j}V${i}.root
    mv Regions.log RegionsH${j}V${i}.log
    ./regionAnalysis.py -i RegionsH${j}V${i}.log -o Total.tsv -H ${j} -V ${i}
  done
done
