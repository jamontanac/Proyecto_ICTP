rm ./To_run/*
for i in {1..32}
do
echo "#!/bin/bash" >> ./To_run/job_$i.sh
echo "#SBATCH -N $i" >> ./To_run/job_$i.sh
echo "#SBATCH --ntasks-per-node=8" >> ./To_run/job_$i.sh
echo "#SBATCH --time=00:30:00" >> ./To_run/job_$i.sh
echo "" >> ./To_run/job_$i.sh
echo "" >> ./To_run/job_$i.sh

echo "module load intel/mpi" >> ./To_run/job_$i.sh
echo "srun ./jacobi.x 12000 10" >> ./To_run/job_$i.sh
done
