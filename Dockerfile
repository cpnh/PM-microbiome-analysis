FROM rocker/tidyverse:4.4.2
WORKDIR /home

# Install the application dependencies
COPY renv.lock ./
COPY .renvignore ./

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        cmake \
        libbz2-dev \
        libcairo2-dev \
        libcurl4-openssl-dev \
        libfreetype6 \
        libfribidi-dev \
        libgdal-dev \
        libgeos-dev \
        libgmp-dev \
        libgsl-dev \
        libharfbuzz-dev \
        libhdf5-dev \
        libjpeg-dev \
        liblapack-dev \
        liblzma-dev \
        libmagick++-dev \
        libmpfr-dev \
        libopenblas-dev \
        libpng-dev \
        libproj-dev \
        libssh2-1-dev \
        libssl-dev \
        libudunits2-dev \
        libx11-dev \
        libxml2-dev \
        libxslt1-dev \
        libxt-dev \
        nano \
        r-cran-renv \
        wget \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN wget https://github.com/quarto-dev/quarto-cli/releases/download/v1.9.31/quarto-1.9.31-linux-amd64.deb \
    && apt install ./quarto-1.9.31-linux-amd64.deb \
    && rm quarto-1.9.31-linux-amd64.deb

RUN R -e 'install.packages("renv", prompt = FALSE)' && R -e 'renv::init()' && R -e 'renv::restore()'
