ARG BASE_IMAGE=stevendargaville/petsc
FROM ${BASE_IMAGE}

# If you want to debug inside the container
# RUN apt-get install -y valgrind 

# Clone PFLARE and run all the tests
RUN git clone https://github.com/stevendargaville/PFLARE.git
RUN cd PFLARE && make && make tests && make python && make tests_python

WORKDIR /build/PFLARE

LABEL maintainer='Steven Dargaville'
LABEL description='PFLARE'
