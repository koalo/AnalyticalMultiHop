FROM ubuntu:16.04
ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8

# petsc
RUN apt-get update -y && apt-get install git build-essential gcc g++ gfortran python -y
WORKDIR /opt/
RUN git clone https://github.com/koalo/petsc
WORKDIR petsc
ENV PETSC_DIR /opt/petsc
ENV PETSC_ARCH linux-gnu-c-debug
RUN ./configure --with-cc=gcc --with-cxx=g++ --with-fc=gfortran --download-mpich --download-fblaslapack
RUN make all test

# analyticalmultihop
WORKDIR /opt/
RUN apt-get install libboost-graph-dev libboost-program-options-dev libboost-system-dev libboost-filesystem-dev -y
RUN git clone https://github.com/koalo/AnalyticalMultiHop analyticalmultihop
WORKDIR analyticalmultihop
RUN make

# python environment
RUN apt-get install python3 python3-pip graphviz libgraphviz-dev pkg-config libfreetype6-dev -y
WORKDIR /opt/analyticalmultihop
RUN pip3 install -r ./requirements.txt

CMD ["/bin/bash","-c"]

