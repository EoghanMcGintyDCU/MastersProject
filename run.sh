#!/bin/bash  
out_dir=materials/Au/Sphere/R50_G5/Ex
mkdir -p $out_dir
mpirun -np 6 python -m mpi4py 2_Nanoparticles_2D.py -m Au -f $out_dir -r 50 -g 5 -c Ex > $out_dir/meep_log.out 2>&1
python Analysis.py -f $out_dir > $out_dir/analysis_log.out 2>&1

out_dir=materials/Au/Sphere/R50_G5/Ey
mkdir -p $out_dir
mpirun -np 6 python -m mpi4py 2_Nanoparticles_2D.py -m Au -f $out_dir -r 50 -g 5 > $out_dir/meep_log.out 2>&1
python Analysis.py -f $out_dir > $out_dir/analysis_log.out 2>&1

out_dir=materials/Au/Sphere/R50_G2/Ex
mkdir -p $out_dir
mpirun -np 6 python -m mpi4py 2_Nanoparticles_2D.py -m Au -f $out_dir -r 50 -g 2 -c Ex > $out_dir/meep_log.out 2>&1
python Analysis.py -f $out_dir > $out_dir/analysis_log.out 2>&1

out_dir=materials/Au/Sphere/R50_G2/Ey
mkdir -p $out_dir
mpirun -np 6 python -m mpi4py 2_Nanoparticles_2D.py -m Au -f $out_dir -r 50 -g 2 > $out_dir/meep_log.out 2>&1
python Analysis.py -f $out_dir > $out_dir/analysis_log.out 2>&1

out_dir=materials/Au/Sphere/R20_G5/Ex
mkdir -p $out_dir
mpirun -np 6 python -m mpi4py 2_Nanoparticles_2D.py -m Au -f $out_dir -r 20 -g 5 -c Ex > $out_dir/meep_log.out 2>&1
python Analysis.py -f $out_dir > $out_dir/analysis_log.out 2>&1

out_dir=materials/Au/Sphere/R20_G5/Ey
mkdir -p $out_dir
mpirun -np 6 python -m mpi4py 2_Nanoparticles_2D.py -m Au -f $out_dir -r 20 -g 5 > $out_dir/meep_log.out 2>&1
python Analysis.py -f $out_dir > $out_dir/analysis_log.out 2>&1

out_dir=materials/Au/Sphere/R20_G2/Ex
mkdir -p $out_dir
mpirun -np 6 python -m mpi4py 2_Nanoparticles_2D.py -m Au -f $out_dir -r 20 -g 2 -c Ex > $out_dir/meep_log.out 2>&1
python Analysis.py -f $out_dir > $out_dir/analysis_log.out 2>&1

out_dir=materials/Au/Sphere/R20_G2/Ey
mkdir -p $out_dir
mpirun -np 6 python -m mpi4py 2_Nanoparticles_2D.py -m Au -f $out_dir -r 20 -g 2 > $out_dir/meep_log.out 2>&1
python Analysis.py -f $out_dir > $out_dir/analysis_log.out 2>&1

# Water
out_dir=materials/Au/Sphere/R50_G5/Water/Ey
mkdir -p $out_dir
mpirun -np 6 python -m mpi4py 2_Nanoparticles_2D.py -m Au -f $out_dir -r 50 -g 5 -i 1.333 > $out_dir/meep_log.out 2>&1
python Analysis.py -f $out_dir > $out_dir/analysis_log.out 2>&1

out_dir=materials/Au/Sphere/R50_G2/Water/Ey
mkdir -p $out_dir
mpirun -np 6 python -m mpi4py 2_Nanoparticles_2D.py -m Au -f $out_dir -r 50 -g 2 -i 1.333 > $out_dir/meep_log.out 2>&1
python Analysis.py -f $out_dir > $out_dir/analysis_log.out 2>&1

out_dir=materials/Au/Sphere/R20_G5/Water/Ey
mkdir -p $out_dir
mpirun -np 6 python -m mpi4py 2_Nanoparticles_2D.py -m Au -f $out_dir -r 20 -g 5 -i 1.333 > $out_dir/meep_log.out 2>&1
python Analysis.py -f $out_dir > $out_dir/analysis_log.out 2>&1

out_dir=materials/Au/Sphere/R20_G2/Water/Ey
mkdir -p $out_dir
mpirun -np 6 python -m mpi4py 2_Nanoparticles_2D.py -m Au -f $out_dir -r 20 -g 2 -i 1.333 > $out_dir/meep_log.out 2>&1
python Analysis.py -f $out_dir > $out_dir/analysis_log.out 2>&1

# Glass
out_dir=materials/Au/Sphere/R50_G5/Glass/Ey
mkdir -p $out_dir
mpirun -np 6 python -m mpi4py 2_Nanoparticles_2D.py -m Au -f $out_dir -r 50 -g 5 -i 1.52 > $out_dir/meep_log.out 2>&1
python Analysis.py -f $out_dir > $out_dir/analysis_log.out 2>&1

out_dir=materials/Au/Sphere/R50_G2/Glass/Ey
mkdir -p $out_dir
mpirun -np 6 python -m mpi4py 2_Nanoparticles_2D.py -m Au -f $out_dir -r 50 -g 2 -i 1.52 > $out_dir/meep_log.out 2>&1
python Analysis.py -f $out_dir > $out_dir/analysis_log.out 2>&1

out_dir=materials/Au/Sphere/R20_G5/Glass/Ey
mkdir -p $out_dir
mpirun -np 6 python -m mpi4py 2_Nanoparticles_2D.py -m Au -f $out_dir -r 20 -g 5 -i 1.52 > $out_dir/meep_log.out 2>&1
python Analysis.py -f $out_dir > $out_dir/analysis_log.out 2>&1

out_dir=materials/Au/Sphere/R20_G2/Glass/Ey
mkdir -p $out_dir
mpirun -np 6 python -m mpi4py 2_Nanoparticles_2D.py -m Au -f $out_dir -r 20 -g 2 -i 1.52 > $out_dir/meep_log.out 2>&1
python Analysis.py -f $out_dir > $out_dir/analysis_log.out 2>&1
