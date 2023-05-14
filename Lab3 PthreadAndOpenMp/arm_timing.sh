 #PBS -o ./arm_timing.xls
 #PBS -e ./error.txt
 #test .sh
 #!/bin/sh
 #PBS -N test
 #PBS -l nodes=1

 pssh -h $PBS_NODEFILE mkdir -p /home/ss2113972 1>&2
 scp master:/home/ss2113972/gauss /home/ss2113972
 pscp -h $PBS_NODEFILE /home/ss2113972/gauss /home/ss2113972 1>&2
 /home/ss2113972/gauss

 pssh -h $PBS_NODEFILE mkdir -p /home/ss2113972 1>&2
 scp master:/home/ss2113972/groebner /home/ss2113972
 pscp -h $PBS_NODEFILE /home/ss2113972/groebner /home/ss2113972 1>&2
 /home/ss2113972/groebner
