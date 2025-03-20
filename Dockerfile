FROM ubuntu

LABEL Description="Sparse2D"

SHELL ["/bin/bash", "-c"]
WORKDIR /workdir

ARG DEBIAN_FRONTEND=noninteractive
ARG CC=gcc-13
ARG CXX=g++-13

RUN apt-get update && \
    apt-get install -y autoconf automake libtool pkg-config libgl1 libglx-mesa0 && \
    apt-get install -y gcc-13 g++-13 && \
    apt-get install -y cmake git wget && \
    apt-get install -y libarmadillo-dev && \
    apt-get install -y libcfitsio-dev && \
    apt-get install -y libfftw3-dev && \
    apt-get install -y libgsl-dev && \
    apt-get install -y libsharp-dev && \
    apt-get install -y libhealpix-cxx-dev && \
    apt-get install -y healpy-data && \
    apt-get install -y python3 python3-pip python3-venv && \
    apt-get clean

RUN cd /home && \
    wget https://github.com/catchorg/Catch2/archive/refs/tags/v3.8.0.tar.gz && \
    tar -xvf v3.8.0.tar.gz && \
    cd Catch2-3.8.0 && \
    cmake -Bbuild -H. -DBUILD_TESTING=OFF && \
    cmake --build build/ --target install

ENV HEALPIX /usr/share/healpy

RUN python3 -m venv /venv

RUN source /venv/bin/activate && \
    pip install --upgrade pip

RUN source /venv/bin/activate && \
    pip install jupyter

COPY . /home

RUN source /venv/bin/activate && \
    cd /home && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make -j8 && \
    make install

ENV PYTHONPATH /home/build/src

RUN echo "source /venv/bin/activate" >> ~/.bashrc

RUN echo -e '#!/bin/bash\njupyter notebook --port=8888 --no-browser --ip=0.0.0.0 --allow-root' > /usr/bin/notebook && chmod +x /usr/bin/notebook