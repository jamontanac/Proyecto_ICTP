nproc=2
for ((i=0; i<$nproc;i++))
do
	echo $i
	mpirun -np 2 ./a.out 10 $i 0 3 3 > matrix_$i.txt
done
python Make_imagen.py

