# !/bin/sh
options=(-DSERIAL_GROEBNER -DSSE_GROEBNER -DAVX2_GROEBNER )
DATA=("1_130_22_8" "2_254_106_53" "3_562_170_53" "4_1011_539_263" "5_2362_1226_453" "6_3799_2759_1953" "7_8399_6375_4535" "8_23045_18748_14325" "9_37960_29304_14921" "10_43577_39477_54274" "11_85401_5724_756")
for i in {10..10}
do 
	for option in "${options[@]}"
	do
		for step in {1..5}
		do
			echo "${DATA[i]} $option $step"
			g++ -o2 -march=native -DBASH_TIMING -DDATA=\"${DATA[i]}\" -DSTEP=$step $option ./groebner.cpp -o ./groebner_x86
			./groebner_x86 >> ./timing_result/groebner_x86_timing.xls
		done
	done
    
done
