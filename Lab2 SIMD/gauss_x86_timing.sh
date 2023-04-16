# !/bin/sh
timeStr=$(date +%m_%d_%H_%M_%S)
options=(-DSERIAL -DSSE_UNALIGNED -DSSE_ALIGNED -DAVX2_UNALIGNED -DAVX2_ALIGNED)
for option in "${options[@]}"
do
	for step in {1..5}
	do
		echo $option" Step: "$step
		for i in {1..32}
		do 
		    g++ -o2 -march=native -DN=$[64*i] -DBASH_TIMING $option ../gauss.cpp -o ./gauss_x86
		    echo "now calculating: "$[64*i]
		    ./gauss_x86 >> ../timing_result/gauss_x86_$timeStr.xls
		done
		echo ""
	done
done
