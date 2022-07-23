mpirun -np 6 python -m mpi4py 2_Nanoparticles_2D.py Au materials/Au/R50_G2 -r 50 -g 2 > log.out 2>&1
python Analysis.py materials/Au/R50_G2