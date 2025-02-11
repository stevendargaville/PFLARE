ARG BASE_IMAGE=stevendargaville/petsc

FROM ${BASE_IMAGE}
ARG BRANCH=main

LABEL maintainer="Steven Dargaville"
LABEL description="PFLARE"

# If you want to debug inside the container
# RUN apt-get update && apt-get install -y valgrind 

ENV CHECKOUT_BRANCH=$BRANCH

WORKDIR /build

# Clone PFLARE and run all the tests
RUN set -e; \
    echo "Cloning branch: ${CHECKOUT_BRANCH}" && \
    git clone --branch ${CHECKOUT_BRANCH} https://github.com/stevendargaville/PFLARE.git && \
    cd PFLARE && \
    make && make tests && \
    make python && make tests_python

WORKDIR /build/PFLARE
