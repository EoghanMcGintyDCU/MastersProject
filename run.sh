mpirun -np 6 python -m mpi4py Nanoparticle_2D.py 2 > log_200.out 2>&1
mpirun -np 6 python -m mpi4py Nanoparticle_2D.py 0.2 > log_20.out 2>&1

mpirun -np 6 python -m mpi4py 2_Nanoparticles_2D.py 0.5 0.5 > log_50_50.out 2>&1
mpirun -np 6 python -m mpi4py 2_Nanoparticles_2D.py 0.5 0.25 > log_50_25.out 2>&1
mpirun -np 6 python -m mpi4py 2_Nanoparticles_2D.py 0.5 0.1 > log_50_10.out 2>&1
mpirun -np 6 python -m mpi4py 2_Nanoparticles_2D.py 0.5 0.05 > log_50_5.out 2>&1