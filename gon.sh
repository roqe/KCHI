#!/bin/bash
#usage: ./gon.sh node_number locus_number run_number

echo -e "Node number:\t$1" >> Report
echo -e "Locus number:\t$2" >> Report
echo -e "Running number:\t$3" >> Report

for ((i = 1 ; i <= $3 ; i++))
do	
	let "RM = 2*($RANDOM%4)"
	./rand $1C"$RM"pedigree $2 
	sleep 1
	mv gdata "$i"gdata
	mv psdata "$i"psdata
	mv hsdata "$i"hsdata

	(time ./kchi $1C"$RM"pedigree "$i"gdata "$i"psdata "$i"hsdata >> errorh ) 2>> ot

	let "RS += $RM"
done

echo -| awk -v K=$RS -v J=$3 '{printf "Avg k (loops):\t%.3f\n", K/J}' >> Report
grep "free" errorh | cut -d ':' -f2 | awk '{ errorh += $1; ++n} END {printf "Avg h (free v):\t%.3f\n", errorh/n}' >> Report
grep "correctly phased markers" sum | cut -d ' ' -f4 | awk '{sum += $1; ++n} END {printf "Avg accuracy:\t%.3f\n", sum/n}' >> Report
grep "real" ot | cut -d 'm' -f2 | awk '{ot += $1; ++n} END {printf "Operation time:\t%.3fs\n\n", ot} '>> Report

mkdir $1$2Data
mv *data $1$2Data/ 
mkdir $1$2Info
mv res sum errorh ot $1$2Info/
