#PBS -o ./timing_result/gauss_arm_timing.xls
#PBS -e ./timing_result/error.txt
#test .sh
#!/bin/sh
#PBS -N test
#PBS -l nodes=1

options=(-DSERIAL -DNEON_ALIGNED_MLS -DNEON_UNALIGNED_MLS -DNEON_ALIGNED_UNMLS -DNEON_UNALIGNED_UNMLS)

for option in "${options[@]}"
do
	for step in {1..5}
	do
		for i in {1..32}
		do 
		    pssh -h $PBS_NODEFILE mkdir -p /home/ss2113972/SIMD/gauss_arm_compile 1>&2
		    name="${option}$[64*i]"
		    scp master:/home/ss2113972/SIMD/gauss_arm_compile/$name /home/ss2113972/SIMD/gauss_arm_compile
		    pscp -h $PBS_NODEFILE /home/ss2113972/SIMD/gauss_arm_compile/$name /home/ss2113972/SIMD/gauss_arm_compile 1>&2
		    /home/ss2113972/SIMD/gauss_arm_compile/$name
		done
	done
done

