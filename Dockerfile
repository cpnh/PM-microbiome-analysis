FROM rocker/tidyverse:4.4.2
WORKDIR /home

# Install the application dependencies
COPY renv.lock ./
COPY .renvignore ./

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        cmake=3.28.3-1build7 \
        libbz2-dev=1.0.8-5.1build0.1 \
        libcairo2-dev=1.18.0-3build1 \
        libcurl4-openssl-dev=8.5.0-2ubuntu10.8 \
        libfreetype6=2.13.2+dfsg-1ubuntu0.1 \
        libfribidi-dev=1.0.13-3build1 \
        libgdal-dev=3.8.4+dfsg-3ubuntu3 \
        libgeos-dev=3.12.1-3build1 \
        libgmp-dev=2:6.3.0+dfsg-2ubuntu6.1 \
        libgsl-dev=2.7.1+dfsg-6ubuntu2 \
        libharfbuzz-dev=8.3.0-2build2 \
        libhdf5-dev=1.10.10+repack-3.1ubuntu4 \
        libjpeg-dev=8c-2ubuntu11 \
        liblapack-dev=3.12.0-3build1.1 \
        liblzma-dev=5.6.1+really5.4.5-1ubuntu0.2 \
        libmagick++-dev=8:6.9.12.98+dfsg1-5.2build2 \
        libmpfr-dev=4.2.1-1build1.1 \
        libopenblas-dev=0.3.26+ds-1ubuntu0.1 \
        libpng-dev=1.6.43-5ubuntu0.5 \
        libproj-dev=9.4.0-1build2 \
        libssh2-1-dev=1.11.0-4.1build2 \
        libssl-dev=3.0.13-0ubuntu3.7 \
        libudunits2-dev=2.2.28-7build1 \
        libx11-dev=2:1.8.7-1build1 \
        libxml2-dev=2.9.14+dfsg-1.3ubuntu3.7 \
        libxslt1-dev=1.1.39-0exp1ubuntu0.24.04.3 \
        libxt-dev=1:1.2.1-1.2build1 \
        nano=7.2-2ubuntu0.1 \
        r-cran-renv=1.0.3-1 \
        wget=1.21.4-1ubuntu4.1 \
        zlib1g-dev=1:1.3.dfsg-3.1ubuntu2.1 \
    && rm -rf /var/lib/apt/lists/*

RUN wget https://github.com/quarto-dev/quarto-cli/releases/download/v1.9.31/quarto-1.9.31-linux-amd64.deb \
    && apt install ./quarto-1.9.31-linux-amd64.deb \
    && rm quarto-1.9.31-linux-amd64.deb

RUN R -e 'install.packages("renv", repos="http://cran.us.r-project.org", prompt = FALSE)' && R -e 'renv::init()' && R -e 'renv::restore()'
