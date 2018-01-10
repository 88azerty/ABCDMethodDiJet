#!/bin/bash

declare -a arrayV=( "45" "50" "55" "60" "65" "70" "75" )
declare -a arrayH=( "150" )
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
  done
done
