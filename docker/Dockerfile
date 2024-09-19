# Use the base Rocky Linux image
FROM rockylinux/rockylinux:9

# Set some labels
LABEL description="Samtools 1.18"
LABEL version="1.0"
LABEL maintainer="bwbioinfo"
LABEL base_image="rockylinux/rockylinux:9"

# Set environment variables
ENV SAMTOOLS_VERSION=1.21

# Install dependencies
RUN dnf update -y
RUN dnf install -y \
    git \
    gcc \
    make \
    zlib-devel \
    bzip2 \
    bzip2-devel \
    xz \
    xz-devel \
    ncurses \
    ncurses-devel \
    libcurl-devel \
    wget \ 
    zsh


# Set the working directory
WORKDIR /usr/local/src

# Clone Samtools repository
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2
RUN tar -xvjf samtools-${SAMTOOLS_VERSION}.tar.bz2

# Change to the Samtools directory
WORKDIR /usr/local/src/samtools-${SAMTOOLS_VERSION}

# Build and install Samtools
RUN make
RUN make install

# Check the Samtools version
RUN samtools --version

# Add scripts to the PATH
WORKDIR /
COPY scripts/* /usr/local/bin/
RUN chmod a+x /usr/local/bin/*.sh
RUN ln -s /usr/local/bin/samtools_multi_merge_v2.sh /usr/bin/multimerge
RUN ln -s /usr/local/bin/bamstats.sh /usr/bin/bamstats

# Clean up
RUN rm -rf /usr/local/src/samtools
RUN dnf remove -y gcc make zlib-devel bzip2-devel xz-devel ncurses-devel libcurl-devel && \
    dnf clean all

# Set the default command to run Samtools
CMD ["samtools"]

# docker build -t samtools-test -f docker/Dockerfile .