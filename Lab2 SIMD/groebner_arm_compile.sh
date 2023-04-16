# !/bin/sh

options=(-DSERIAL_GROEBNER -DNEON_GROEBNER)
DATA=("1_130_22_8" "2_254_106_53" "3_562_170_53" "4_1011_539_263" "5_2362_1226_453" "6_3799_2759_1953" "7_8399_6375_4535" "8_23045_18748_14325" "9_37960_29304_14921" "10_43577_39477_54274" "11_85401_5724_756")

for i in {0..10}
do 
	for option in "${options[@]}"
	do
		echo "${DATA[i]} $option $step"
		g++ -o2 -std=c++11 -march=native -DBASH_TIMING -DDATA=\"${DATA[i]}\" $option ./groebner.cpp -o ./groebner_arm_compile/$i$option
	done
    
done
