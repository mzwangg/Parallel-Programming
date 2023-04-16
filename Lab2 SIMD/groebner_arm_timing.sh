#PBS -o ./timing_result/groebner_arm_timing.xls
#PBS -e ./timing_result/error.txt
#test .sh
#!/bin/sh
#PBS -N test
#PBS -l nodes=1

options=(-DSERIAL_GROEBNER -DNEON_GROEBNER)
DATA=("1_130_22_8" "2_254_106_53" "3_562_170_53" "4_1011_539_263" "5_2362_1226_453" "6_3799_2759_1953" "7_8399_6375_4535" "8_23045_18748_14325" "9_37960_29304_14921" "10_43577_39477_54274" "11_85401_5724_756")
for i in {0..10}
do 
	for option in "${options[@]}"
	do
		#for step in {1..5}
		#do
		    pssh -h $PBS_NODEFILE mkdir -p /home/ss2113972/SIMD/groebner_arm_compile 1>&2
		    name="${i}${option}"
		    echo name
		    scp master:/home/ss2113972/SIMD/groebner_arm_compile/$name /home/ss2113972/SIMD/groebner_arm_compile
		    pscp -h $PBS_NODEFILE /home/ss2113972/SIMD/groebner_arm_compile/$name /home/ss2113972/SIMD/groebner_arm_compile 1>&2
		    /home/ss2113972/SIMD/groebner_arm_compile/$name
		#done
	done
done
