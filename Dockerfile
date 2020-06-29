FROM ubuntu:18.04

MAINTAINER Tomer Altman, Altman Analytics LLC

Workdir /root

### Install apt dependencies

RUN DEBIAN_FRONTEND=noninteractive apt-get update
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y cpanminus \
    build-essential \
    libnet-ssleay-perl \
    libcrypt-ssleay-perl \
    git \
    curl \
    zip \
    autoconf

### Perl dependencies:

RUN cpanm install Inline
RUN cpanm install Inline::C
RUN cpanm install LWP::Simple
RUN cpanm install LWP::Protocol::https

### Install VADR:
RUN git clone https://github.com/nawrockie/vadr.git
RUN cd vadr && ./vadr-install.sh linux

### Set up VADR environment variables:
ENV VADRINSTALLDIR="/root/vadr"
ENV VADRSCRIPTSDIR="$VADRINSTALLDIR/vadr"
ENV VADRMODELDIR="$VADRINSTALLDIR/vadr-models"
ENV VADRINFERNALDIR="$VADRINSTALLDIR/infernal/binaries"
ENV VADREASELDIR="$VADRINSTALLDIR/infernal/binaries"
ENV VADRHMMERDIR="$VADRINSTALLDIR/hmmer/binaries"
ENV VADRBIOEASELDIR="$VADRINSTALLDIR/Bio-Easel"
ENV VADRSEQUIPDIR="$VADRINSTALLDIR/sequip"
ENV VADRBLASTDIR="$VADRINSTALLDIR/ncbi-blast/bin"
ENV PERL5LIB="$VADRSCRIPTSDIR:$VADRSEQUIPDIR:$VADRBIOEASELDIR/blib/lib:$VADRBIOEASELDIR/blib/arch:$PERL5LIB"
ENV PATH="$VADRSCRIPTSDIR:$PATH"
