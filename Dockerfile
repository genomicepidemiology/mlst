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

# Install kma
RUN mkdir -p /usr/src/kma; \
    wget -O - https://bitbucket.org/genomicepidemiology/kma/get/master.tar.gz | \
    tar xzf - --strip-components=1 --wildcards *.c; \
    gcc -O3 -o /usr/local/bin/kma KMA.c -lm -lpthread && \
    gcc -O3 -o /usr/local/bin/kma_index KMA_index.c -lm; \

# Install cgefinder for kma and blast use
    mkdir -p /usr/src/CGE; \
    git clone --recursive -b 4.0 https://bitbucket.org/genomicepidemiology/resfinder.git; \
    mv resfinder/cge/cgefinder.py resfinder/cge/blaster/ resfinder/cge/__init__.py /usr/src/CGE; \
    rm -rf resfinder;

COPY ./mlst.py /usr/src/mlst.py

WORKDIR /usr/src/

# Execute program when running the container
# ENTRYPOINT ["mlst.py"]
# ENV PATH $PATH:/usr/src/

# Setup .bashrc file for convenience during debugging
RUN echo "alias ls='ls -h --color=tty'\n"\
"alias ll='ls -lrt'\n"\
"alias l='less'\n"\
"alias du='du -hP --max-depth=1'\n"\
"alias cwd='readlink -f .'\n"\
"PATH=$PATH\n">> ~/.bashrc
