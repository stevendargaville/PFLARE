ARG BASE_IMAGE=ubuntu:22.04
FROM ${BASE_IMAGE}

RUN echo Etc/UTC > /etc/timezone && ln -s /usr/share/zoneinfo/Etc/UTC /etc/localtime
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
  autoconf \
  automake \
  bash-completion \
  ca-certificates \
  chrpath \
  cmake \
  curl \
  g++ \
  gcc \
  gfortran \
  git \
  less \
  libblis-serial-dev \
  liblapack-dev \
  libtool \
  locales \
  m4 \
  make \
  openssl \
  pkg-config \
  ripgrep \
  zlib1g-dev \
  python3-virtualenv \
  python3.10-venv \
  python3.10-dev \
  wget \
  && rm -rf /var/lib/apt/lists/*

RUN locale-gen en_US.UTF-8 && update-locale LANG=en_US.UTF-8
ENV LANG=en_US.UTF-8

WORKDIR /build
ARG MPICH_DEVICE=ch3:nemesis
ARG MPICH_VERSION=4.0.2
RUN \
  curl -O https://www.mpich.org/static/downloads/${MPICH_VERSION}/mpich-${MPICH_VERSION}.tar.gz && \
  tar xf mpich-${MPICH_VERSION}.tar.gz && \
  rm mpich-${MPICH_VERSION}.tar.gz && \
  cd mpich-${MPICH_VERSION} && \
  ./configure \
      FFLAGS=-fallow-argument-mismatch \
      FCFLAGS=-fallow-argument-mismatch \
      --disable-wrapper-rpath \
      --with-device=${MPICH_DEVICE} \
      --enable-error-checking=runtime \
      --enable-error-messages=all \
      --enable-g=meminit && \
  make -j$(nproc) && \
  make install && \
  cd /build && \
  rm -rf /build/mpich-${MPICH_VERSION} && \
  ldconfig

# Hide irrelevant MPICH output warnings
ENV HWLOC_HIDE_ERRORS=2
ARG PETSC_GIT_BRANCH=release

# Have to build a venv for python as pip/numpy now complain
ENV VIRTUAL_ENV=/my-venv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# Install python packages (should be in the venv not globally)
RUN wget https://bootstrap.pypa.io/get-pip.py \
 && python3 ./get-pip.py && pip install numpy && pip install Cython && pip install setuptools

# Install petsc
ENV PETSC_DIR=/build/petsc
ENV PETSC_ARCH=arch-linux-c-opt
RUN git clone --depth=1 --branch=$PETSC_GIT_BRANCH https://gitlab.com/petsc/petsc.git && \
  cd petsc && \
  OPT='-O2 -march=haswell -ffp-contract=fast'; \
  python3 configure \
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

WORKDIR /build/PFLARE

LABEL maintainer='Steven Dargaville'
LABEL description='PFLARE'
