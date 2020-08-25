FROM continuumio/miniconda:4.7.12

LABEL authors="2592722S@student.gla.ac.uk" \
      description="Docker image containing base requirements for the eccDNA-pipeline"

# copy conda environment into docker file storage

COPY /environment.yml .
COPY /main.nf .
COPY /test-datasets .
COPY /conf .
COPY /nextflow.config .
COPY /bin .
COPY testdata .

# run conda environment inside docker container
SHELL ["/bin/bash", "--login", "-c"]

RUN conda env create -f /environment.yml && conda clean -a
RUN echo "source activate myenv" > ~/.bashrc
ENV PATH /opt/conda/envs/env/bin:$PATH
SHELL ["conda", "run", "-n", "myenv", "/bin/bash", "-c"]

# Pull the environment name out of the environment.yml
RUN echo "source activate myenv" > ~/.bashrc
ENV PATH /opt/conda/envs/myenv/bin:$PATH