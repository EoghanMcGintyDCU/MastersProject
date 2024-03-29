FROM continuumio/miniconda3

RUN conda create -n pmp -c conda-forge pymeep=*=mpi_mpich_*

RUN conda activate pmp

RUN conda install -c conda-forge matplotlib pymiescatt

WORKDIR /app
COPY . .




