#!/bin/bash  

out_dir=materials/Au/Rod/R50_G5/Hz
mkdir -p $out_dir
mpirun -np 6 python -m mpi4py 2_Nanoparticles_2D_DFT.py -f $out_dir -m Au -s rod -r 50 -g 5 > $out_dir/log.out 2>&1

out_dir=materials/Au/Rod/R50_G2/Hz
mkdir -p $out_dir
mpirun -np 6 python -m mpi4py 2_Nanoparticles_2D_DFT.py -f $out_dir -m Au -s rod -r 50 -g 2 > $out_dir/log.out 2>&1

out_dir=materials/Au/Rod/R50_G0.5/Hz
mkdir -p $out_dir
mpirun -np 6 python -m mpi4py 2_Nanoparticles_2D_DFT.py -f $out_dir -m Au -s rod -r 50 -g 0.5 > $out_dir/log.out 2>&1

out_dir=materials/Au/Rod/R20_G5/Hz
mkdir -p $out_dir
mpirun -np 6 python -m mpi4py 2_Nanoparticles_2D_DFT.py -f $out_dir -m Au -s rod -r 20 -g 5 > $out_dir/log.out 2>&1

out_dir=materials/Au/Rod/R20_G2/Hz
mkdir -p $out_dir
mpirun -np 6 python -m mpi4py 2_Nanoparticles_2D_DFT.py -f $out_dir -m Au -s rod -r 20 -g 2 > $out_dir/log.out 2>&1

out_dir=materials/Au/Rod/R20_G0.5/Hz
mkdir -p $out_dir
mpirun -np 6 python -m mpi4py 2_Nanoparticles_2D_DFT.py -f $out_dir -m Au -s rod -r 20 -g 0.5 > $out_dir/log.out 2>&1

# Water
out_dir=materials/Au/Rod/R50_G5/Water/Hz
mkdir -p $out_dir
mpirun -np 6 python -m mpi4py 2_Nanoparticles_2D_DFT.py -f $out_dir -m Au -s rod -r 50 -g 5 -i 1.333 > $out_dir/log.out 2>&1

out_dir=materials/Au/Rod/R20_G2/Water/Hz
mkdir -p $out_dir
mpirun -np 6 python -m mpi4py 2_Nanoparticles_2D_DFT.py -f $out_dir -m Au -s rod -r 20 -g 2 -i 1.333 > $out_dir/log.out 2>&1

# Glass
out_dir=materials/Au/Rod/R50_G5/Glass/Hz
mkdir -p $out_dir
mpirun -np 6 python -m mpi4py 2_Nanoparticles_2D_DFT.py -f $out_dir -m Au -s rod -r 50 -g 5 -i 1.52 > $out_dir/log.out 2>&1

out_dir=materials/Au/Rod/R20_G2/Glass/Hz
mkdir -p $out_dir
mpirun -np 6 python -m mpi4py 2_Nanoparticles_2D_DFT.py -f $out_dir -m Au -s rod -r 20 -g 2 -i 1.52 > $out_dir/log.out 2>&1