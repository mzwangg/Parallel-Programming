# !/bin/sh

options=(-DSERIAL -DNEON_ALIGNED_MLS -DNEON_UNALIGNED_MLS -DNEON_ALIGNED_UNMLS -DNEON_UNALIGNED_UNMLS)
for option in "${options[@]}"
do
	echo $option
	for i in {1..32}
	do 
	    g++ -o2 -march=native -DN=$[64*i] -DBASH_TIMING $option gauss.cpp -o ./gauss_arm_compile/"${option}$[64*i]"
	    echo "N=$[64*i]"
	done
	echo ""
done
