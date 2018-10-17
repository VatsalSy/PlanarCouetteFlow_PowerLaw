CC99='mpicc -std=c99' qcc -Wall -O2 -D_MPI=1 $file.c -o $file -lm
#mpirun -np 4 ./$file
