FROM debian:stretch

ENV DEBIAN_FRONTEND noninteractive

### RUN set -ex; \

RUN apt-get update -qq; \
    apt-get install -y -qq git \
    apt-utils \
    wget \
    python3-pip \
    ncbi-blast+ \
    ; \
    rm -rf /var/cache/apt/* /var/lib/apt/lists/*;
    
ENV DEBIAN_FRONTEND Teletype

# Install python dependencies
RUN pip3 install -U biopython tabulate cgecore==1.2.2;

# Install kma and download cgefinder.py
RUN mkdir -p /usr/src/kma; \
    wget -O - https://bitbucket.org/genomicepidemiology/resfinder/get/4.0.tar.gz | \
    tar xzf - --strip-components=1 --wildcards  *cge; \
    wget -O - https://bitbucket.org/genomicepidemiology/kma/get/master.tar.gz | \
    tar xzf - --strip-components=1 --wildcards *.c; \
    gcc -O3 -o /usr/local/bin/kma KMA.c -lm -lpthread && \
    gcc -O3 -o /usr/local/bin/kma_index KMA_index.c -lm; \
    mkdir -p /usr/src/CGE; \
    mv cge/cgefinder.py cge/__init__.py  /usr/src/CGE/; \
    rm -rf cge;

COPY ./mlst.py /usr/src/mlst.py

WORKDIR /usr/src/


# Setup .bashrc file for convenience during debugging
RUN echo "alias ls='ls -h --color=tty'\n"\
"alias ll='ls -lrt'\n"\
"alias l='less'\n"\
"alias du='du -hP --max-depth=1'\n"\
"alias cwd='readlink -f .'\n"\
"PATH=$PATH\n">> ~/.bashrc
