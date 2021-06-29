FROM continuumio/miniconda3

# Create specific environment:
RUN conda create -n genomePlotter python=3.6

# Try building the environment:
# ADD my_env.yml /tmp/my_env.yml
# RUN conda env create -f /tmp/my_env.yml

RUN conda config --add channels conda-canary
RUN conda config --add channels anaconda
RUN conda config --add channels bioconda
RUN conda config --set always_yes yes --set changeps1 no
RUN conda update -q conda

RUN conda install --name genomePlotter cairo
RUN conda install --name genomePlotter numpy
RUN conda install --name genomePlotter pandas 
RUN conda install --name genomePlotter pybedtools
RUN conda install --name genomePlotter pytest 
RUN conda install --name genomePlotter coverage
RUN conda install -c conda-forge  --name genomePlotter codacy-coverage

ENV PATH /opt/conda/envs/genomePlotter/bin:$PATH
RUN /bin/bash -c "source activate genomePlotter"
RUN /bin/bash -c "conda env export -n genomePlotter | grep -v prefix > environment.yml"

[![GenomePlotter](https://circleci.com/gh/DSuveges/GenomePlotter.svg?style=svg)](https://app.circleci.com/pipelines/github/DSuveges/GenomePlotter)
