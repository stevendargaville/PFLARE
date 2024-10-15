ARG BASE_IMAGE=jedbrown/mpich:latest
FROM ${BASE_IMAGE}
ARG PETSC_GIT_BRANCH=release

# Install the python we need
RUN apt-get update \
 && apt-get install wget \
 && apt-get install -y python3.11-venv \
 && apt-get install -y python3.11-dev

# Have to build a venv for python as pip/numpy now complain
ENV VIRTUAL_ENV=/my-venv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# Install python packages (should be in the venv not globally)
RUN wget https://bootstrap.pypa.io/get-pip.py \
 && python3 ./get-pip.py && pip install numpy && pip install Cython

# Install petsc
ENV PETSC_DIR=/build/petsc
ENV PETSC_ARCH=arch-linux-c-opt
RUN git clone --depth=1 --branch=$PETSC_GIT_BRANCH https://gitlab.com/petsc/petsc.git && \
  cd petsc && \
  OPT='-O2 -march=haswell -ffp-contract=fast'; \
  python3 configure \
    --with-cxx-dialect=C++14 \
    --with-debugging=0 COPTFLAGS="$OPT" CXXOPTFLAGS="$OPT" FOPTFLAGS="$OPT" \
    --with-mpi-dir=/usr/local \
    --download-metis \
    --download-parmetis \
    --download-fblaslapack \
    --with-petsc4py \
    && \
  make

# If you want to debug inside the container
# RUN apt-get install -y valgrind 

# Clone PFLARE and run all the tests
RUN git clone https://github.com/stevendargaville/PFLARE.git
RUN cd PFLARE && make && make tests && make python && make tests_python

RUN useradd -ms /bin/bash pflare
USER pflare
WORKDIR /home/PFLARE

LABEL maintainer='Steven Dargaville'
LABEL description='PFLARE'
