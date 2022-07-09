mpirun -np 6 python -m mpi4py 2_Nanoparticles_2D.py 0.25 0.5 > log.out 2>&1
mpirun -np 6 python -m mpi4py 2_Nanoparticles_2D.py 0.25 0.25 > log.out 2>&1
mpirun -np 6 python -m mpi4py 2_Nanoparticles_2D.py 0.25 0.1 > log.out 2>&1
mpirun -np 6 python -m mpi4py 2_Nanoparticles_2D.py 0.25 0.05 > log.out 2>&1