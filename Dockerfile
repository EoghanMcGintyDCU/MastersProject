FROM continuumio/miniconda3

RUN conda create -n mp -c conda-forge pymeep pymeep-extras

RUN conda activate mp

RUN conda install -c conda-forge matplotlib

RUN conda install -c anaconda numpy

COPY MaterialDispersion.ipynb /

