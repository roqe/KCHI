#!/bin/bash
#usage: ./gon.sh node_number locus_number run_number

#N=(30 60 100 130 160 200 230 260 300 330 360 400)
#L=(10 30 60 100 130 160 200 230 260 300)

#for ((i = 0; i < ${#N[*]}; i++))
#do
#	for ((j = 0; j < ${#L[*]}; j++))
#	do

#	./gon.sh ${N[$i]} ${L[$j]} 100
	./gon.sh 30 130 100
	./gon.sh 30 200 100
	./gon.sh 60 60 100
	./gon.sh 230 30 100
	./gon.sh 230 230 100
	./gon.sh 260 300 100
	./gon.sh 300 200 100
	./gon.sh 330 30 100
	./gon.sh 330 130 100
	./gon.sh 360 10 100
#	done
#done
